#pragma once

#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/utility/globalindexset.hh>
#include <dune/pdelab/backend/common/tags.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/ordering/orderingbase.hh>

#include <petscis.h>
#include <petscmat.h>

#include <cassert>
#include <petscsystypes.h>
#include <type_traits>
#include <vector>

namespace Dune::PDELab::Petsc {
namespace {
struct MatrixEntryRef {
  MatrixEntryRef(Mat M, PetscInt row, PetscInt col) : M{M}, row{row}, col{col} {}

  MatrixEntryRef operator+=(PetscScalar v)
  {
    PetscFunctionBegin;
    PetscCallAbort(MPI_COMM_SELF, MatSetValuesLocal(M, 1, &row, 1, &col, &v, ADD_VALUES));
    PetscFunctionReturn(*this);
  }

  Mat M;
  PetscInt row;
  PetscInt col;
};
} // namespace

template <class M>
struct MatrixPatternInserter {
  MatrixPatternInserter(M &mat, const PetscScalar &x) : mat{mat}, x{x} {}

  template <class RowIndex, class ColIndex>
  void add_link(const RowIndex &i, const ColIndex &j)
  {
    PetscFunctionBegin;
    PetscCallVoid(MatSetValue(mat, i[0], j[0], x, INSERT_VALUES));
    PetscFunctionReturnVoid();
  }

  M &mat;
  const PetscScalar &x;
};

template <typename GFSV, typename GFSU, typename ET>
class MatrixContainer : public Backend::impl::Wrapper<::Mat> {
public:
  using Container = ::Mat;

private:
  friend Backend::impl::Wrapper<Container>;

public:
  using ElementType = ET;
  using field_type = ET;
  using size_type = std::size_t;
  using index_type = std::size_t;
  using block_type = PetscScalar;

  typedef GFSU TrialGridFunctionSpace;
  typedef GFSV TestGridFunctionSpace;

  using RowIndex = typename GFSV::Ordering::Traits::ContainerIndex;
  using ColIndex = typename GFSU::Ordering::Traits::ContainerIndex;

  template <typename RowCache, typename ColCache>
  using LocalView = UncachedMatrixView<MatrixContainer, RowCache, ColCache>;

  template <typename RowCache, typename ColCache>
  using ConstLocalView = ConstUncachedMatrixView<const MatrixContainer, RowCache, ColCache>;

  using RowOrdering = OrderingBase<typename GFSV::Ordering::Traits::DOFIndex, typename GFSV::Ordering::Traits::ContainerIndex>;
  using ColOrdering = OrderingBase<typename GFSU::Ordering::Traits::DOFIndex, typename GFSU::Ordering::Traits::ContainerIndex>;

  using Pattern = MatrixPatternInserter<Container>;

  template <typename GO>
  explicit MatrixContainer(const GO &go)
  {
    allocate_matrix(go);
  }

  MatrixEntryRef operator()(const RowIndex &ri, const ColIndex &ci)
  {
    PetscFunctionBegin;
    Mat Mlocal;
    PetscCallAbort(MPI_COMM_SELF, MatISGetLocalMat(Mp, &Mlocal));
    MatrixEntryRef ref(Mlocal, ri[0], ci[0]);
    PetscCallAbort(MPI_COMM_SELF, MatISRestoreLocalMat(Mp, &Mlocal));
    PetscFunctionReturn(ref);
  }

  void clear_row(const RowIndex &ri, const ElementType &diagonal_entry)
  {
    PetscFunctionBegin;
    finalize();
    const PetscInt row = static_cast<PetscInt>(ri[0]);
    Mat Mlocal;
    PetscCallVoid(MatISGetLocalMat(Mp, &Mlocal));
    PetscCallVoid(MatZeroRows(Mlocal, 1, &row, diagonal_entry, nullptr, nullptr));
    PetscCallVoid(MatISRestoreLocalMat(Mp, &Mlocal));
    PetscFunctionReturnVoid();
  }

  void flush()
  {
    PetscFunctionBegin;
    PetscCallVoid(MatAssemblyBegin(Mp, MAT_FLUSH_ASSEMBLY));
    PetscCallVoid(MatAssemblyEnd(Mp, MAT_FLUSH_ASSEMBLY));
    PetscFunctionReturnVoid();
  }

  void finalize()
  {
    PetscFunctionBegin;
    PetscBool mat_is_assembled;
    PetscCallVoid(MatAssembled(Mp, &mat_is_assembled));
    if (mat_is_assembled == PETSC_FALSE) {
      PetscCallVoid(MatAssemblyBegin(Mp, MAT_FINAL_ASSEMBLY));
      PetscCallVoid(MatAssemblyEnd(Mp, MAT_FINAL_ASSEMBLY));
    }
    PetscFunctionReturnVoid();
  }

  template <typename GO, typename CC>
  void fix_dirichlet_rows(const GO &go, const CC &cc)
  {
    PetscFunctionBegin;
    using GFS = std::remove_cvref_t<decltype(go.trialGridFunctionSpace())>;
    Dune::PDELab::Backend::Vector<GFS, PetscScalar> dirichlet_mask(go.trialGridFunctionSpace());

    dirichlet_mask = 0;
    set_constrained_dofs(cc, 1., dirichlet_mask);
    dirichlet_mask.assemble();

    Mat Mlocal;
    PetscCallVoid(MatISGetLocalMat(Mp, &Mlocal));

    for (std::size_t i = 0; i < owns_index.size(); ++i) {
      if (owns_index[i])
        continue;

      if (dirichlet_mask[i] > 0) {
        PetscInt ri = i;
        PetscCallVoid(MatZeroRows(Mlocal, 1, &ri, 0, nullptr, nullptr));
      }
    }
    PetscCallVoid(MatISRestoreLocalMat(Mp, &Mlocal));
    PetscFunctionReturnVoid();
  }

  ~MatrixContainer()
  {
    PetscFunctionBegin;
    Mat Mlocal;
    PetscCallVoid(MatISGetLocalMat(Mp, &Mlocal));
    PetscCallVoid(MatDestroy(&Mlocal));
    PetscCallVoid(MatDestroy(&Mp));
    PetscFunctionReturnVoid();
  }

private:
  Container &native() { return Mp; }
  const Container &native() const { return Mp; }

  template <typename GO>
  void allocate_matrix(const GO &go, const ElementType &et = ElementType{0})
  {
    PetscFunctionBegin;
    const auto rank = go.trialGridFunctionSpace().gridView().grid().comm().rank();

    // Create the parallel matrix.
    // First we find out how many "unique" dofs we own. The rule is: a shared dof is owned by
    // the rank with the smallest rank number. This can be determined by creating a vector that
    // has the rank number as its (constant) value and then performing a communication that sets
    // the minimum at shared dofs. If the values at a shared dof still equals our rank number,
    // we are the owner.
    using GFS = std::remove_cvref_t<decltype(go.trialGridFunctionSpace())>;
    using GIVec = Dune::PDELab::Backend::Vector<GFS, PetscScalar>;
    GIVec giv(go.trialGridFunctionSpace());
    giv = rank;
    giv.assemble();
    {
      Dune::PDELab::MinDataHandle mindh(go.trialGridFunctionSpace(), giv);
      go.trialGridFunctionSpace().gridView().communicate(mindh, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
    }
    giv.assemble();

    PetscInt n_local = 0;
    owns_index.resize(giv.size(), false);
    for (std::size_t i = 0; i < giv.size(); ++i) {
      if (giv[i] == rank) {
        n_local += 1;
        owns_index[i] = true;
      }
    }

    std::vector<int> n_local_all(go.trialGridFunctionSpace().gridView().grid().comm().size());
    MPI_Allgather(&n_local, 1, MPI_INT, n_local_all.data(), 1, MPI_INT, MPI_COMM_WORLD);
    auto offset = std::accumulate(n_local_all.begin(), n_local_all.begin() + rank, 0);

    PetscInt cnt = 0;
    for (PetscInt i = 0; i < giv.size(); ++i) {
      if (giv[i] == rank)
        giv[i] = offset + cnt++;
      else
        giv[i] = 0;
    }
    giv.assemble();
    {
      Dune::PDELab::AddDataHandle adddh(go.trialGridFunctionSpace(), giv);
      go.trialGridFunctionSpace().gridView().communicate(adddh, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
    }
    giv.assemble();

    std::vector<PetscInt> gidxs(giv.size());
    for (std::size_t i = 0; i < giv.size(); ++i) {
      gidxs[i] = giv[i];
    }

    // Convert to ISLocalToGlobalMapping
    ISLocalToGlobalMapping ltg;
    PetscCallVoid(ISLocalToGlobalMappingCreate(MPI_COMM_WORLD, 1, gidxs.size(), gidxs.data(), PETSC_COPY_VALUES, &ltg));
    PetscCallVoid(MatCreateIS(MPI_COMM_WORLD, 1, n_local, n_local, PETSC_DECIDE, PETSC_DECIDE, ltg, ltg, &Mp));
    PetscCallVoid(ISLocalToGlobalMappingDestroy(&ltg));

    // Create the processor local matrix
    Mat Mlocal;
    auto rows = go.testGridFunctionSpace().ordering().blockCount();
    auto cols = go.trialGridFunctionSpace().ordering().blockCount();
    PetscCallVoid(MatCreateSeqAIJ(MPI_COMM_SELF, rows, cols, go.matrixBackend().avg_nz_per_row, nullptr, &Mlocal));
    Pattern pattern(Mlocal, et);
    go.fill_pattern(pattern);
    PetscCallVoid(MatAssemblyBegin(Mlocal, MAT_FINAL_ASSEMBLY));
    PetscCallVoid(MatAssemblyEnd(Mlocal, MAT_FINAL_ASSEMBLY));

    PetscCallVoid(MatISSetLocalMat(Mp, Mlocal));
    PetscCallVoid(MatSetUp(Mp));
    PetscFunctionReturnVoid();
  }

  Mat Mp; // The parallel matrix

  std::vector<bool> owns_index;
};

} // namespace Dune::PDELab::Petsc
