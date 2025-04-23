// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <dune/common/timer.hh>
#include <petscmat.h>
#include <petscvec.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <bitset>
#include <iostream>
#include <string>

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/pdelab/backend/petsc.hh>
#include <dune/pdelab/common/partitionviewentityset.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include <dune/pdelab/localoperator/permeability_adapter.hh>

#include <petscksp.h>
#include <petscsys.h>

//===============================================================
//===============================================================
// Solve the Poisson equation
//           - \Delta u = f in \Omega,
//                    u = g on \partial\Omega_D
//  -\nabla u \cdot \nu = j on \partial\Omega_N
//===============================================================
//===============================================================

//===============================================================
// Define parameter functions f,g,j and \partial\Omega_D/N
//===============================================================

template <class GridView, class RF>
class IslandsModelProblem : public Dune::PDELab::ConvectionDiffusionModelProblem<GridView, RF> {
  using BC = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type;

public:
  using Traits = typename Dune::PDELab::ConvectionDiffusionModelProblem<GridView, RF>::Traits;

  typename Traits::PermTensorType A(const typename Traits::ElementType &e, const typename Traits::DomainType &) const
  {
    auto xg = e.geometry().center();

    int ix = std::floor(15.0 * xg[0]);
    int iy = std::floor(15.0 * xg[1]);
    auto x = xg[0];
    auto y = xg[1];

    double kappa = 1.0;

    if (x > 0.3 && x < 0.9 && y > 0.6 - (x - 0.3) / 6 && y < 0.8 - (x - 0.3) / 6) {
      kappa = pow(10, 5.0) * (x + y) * 10.0;
    }

    if (x > 0.1 && x < 0.5 && y > 0.1 + x && y < 0.25 + x) {
      kappa = pow(10, 5.0) * (1.0 + 7.0 * y);
    }

    if (x > 0.5 && x < 0.9 && y > 0.15 - (x - 0.5) * 0.25 && y < 0.35 - (x - 0.5) * 0.25) {
      kappa = pow(10, 5.0) * 2.5;
    }

    if (ix % 2 == 0 && iy % 2 == 0) {
      kappa = pow(10, 5.0) * (1.0 + ix + iy);
    }

    if (GridView::Grid::dimension == 3) {
      const auto radius = 0.02;
      const auto square = [](auto &&v) { return v * v; };

      int divx = 9;
      int divz = 9;
      for (int i = 0; i < divx - 1; ++i) {
        for (int j = 0; j < divz - 1; ++j) {
          typename Traits::DomainType centre{(i + 1.) / divx, (j + 1.) / divz};

          if (square(xg[0] - centre[0]) + square(xg[2] - centre[1]) < square(radius)) {
            kappa = 1e6;
          }
        }
      }
    }

    typename Traits::PermTensorType I;
    for (std::size_t i = 0; i < Traits::dimDomain; i++) {
      for (std::size_t j = 0; j < Traits::dimDomain; j++) {
        I[i][j] = (i == j) ? kappa : 0;
      }
    }

    return I;
  }

  typename Traits::RangeFieldType f(const typename Traits::ElementType &, const typename Traits::DomainType &) const { return 0.0; }

  typename Traits::RangeFieldType g(const typename Traits::ElementType &e, const typename Traits::DomainType &xlocal) const
  {
    auto xglobal = e.geometry().global(xlocal);
    return 1. - xglobal[0];
  }

  BC bctype(const typename Traits::IntersectionType &is, const typename Traits::IntersectionDomainType &x) const
  {
    auto xglobal = is.geometry().global(x);
    if (xglobal[0] < 1e-6) {
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    }
    if (xglobal[0] > 1.0 - 1e-6) {
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    }
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
  }
};

//===============================================================
// Problem setup and solution
//===============================================================

template <typename GV, typename FEM>
void poisson(const GV &gv, const FEM &fem, std::string filename)
{
  PetscFunctionBegin;
  Dune::Timer timer;
  timer.reset();

  using Dune::PDELab::Backend::native;

  using EntitySet = Dune::PDELab::NonOverlappingEntitySet<GV>;
  EntitySet es(gv);

  // make grid function space
  using RF = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
  using Con = Dune::PDELab::ConformingDirichletConstraints;
  using GFS = Dune::PDELab::GridFunctionSpace<EntitySet, FEM, Con, Dune::PDELab::Petsc::VectorBackend>;
  GFS gfs(es, fem);
  gfs.name("solution");

  // make model problem
  using Problem = IslandsModelProblem<EntitySet, RF>;
  Problem problem;

  using CC = typename GFS::template ConstraintsContainer<RF>::Type;
  CC cg;
  cg.clear();
  Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> bctype(es, problem);
  Dune::PDELab::constraints(bctype, gfs, cg);

  // make local operator
  using LocalOperator = Dune::PDELab::ConvectionDiffusionFEM<Problem, FEM>;
  LocalOperator lop(problem);

  using MBE = Dune::PDELab::Petsc::MatrixBackend;
  int avgnz = 27;
  MBE mbe(avgnz);

  // make grid operator
  using GridOperator = Dune::PDELab::GridOperator<GFS, GFS, LocalOperator, MBE, RF, RF, RF, CC, CC>;
  GridOperator gridoperator(gfs, cg, gfs, cg, lop, mbe);

  // make coefficent Vector and initialize it from a function
  // There is some weird shuffling around here - please leave it in,
  // it's there to test the copy constructor and assignment operator of the
  // vector wrapper.
  typedef typename GridOperator::Traits::Domain DV;
  DV x0(gfs, Dune::PDELab::Backend::unattached_container());
  {
    DV x1(gfs);
    DV x2(x1);
    x2 = 0.0;
    x0 = x1;
    x0 = x2;
  }

  // initialize DOFs from Dirichlet extension
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g(es, problem);
  Dune::PDELab::interpolate(g, gfs, x0);
  Dune::PDELab::set_nonconstrained_dofs(cg, 0.0, x0);

  DV dirichlet_mask(gfs);
  dirichlet_mask = 0;
  Dune::PDELab::set_constrained_dofs(cg, 1.0, dirichlet_mask);
  dirichlet_mask.assemble();

  // represent operator as a matrix
  // There is some weird shuffling around here - please leave it in,
  // it's there to test the copy constructor and assignment operator of the
  // matrix wrapper.
  typedef typename GridOperator::Traits::Jacobian M;
  M m(gridoperator);
  {
    // M m1(gridoperator);
    // M m2(m1);
    //   m2 = 0.0;
    //   m = m1;
    //   m = m2;
  }
  // m.fix_dirichlet_rows(gridoperator, cg);

  // Create ghost vec
  PetscInt rstart, rend;
  MatGetOwnershipRange(native(m), &rstart, &rend);
  auto n_local = rend - rstart;
  auto n_ghost = x0.size() - n_local;
  std::cout << "n_local = " << n_local << ", n_ghost = " << n_ghost << "\n";

  using Vec = Dune::PDELab::Backend::Vector<GFS, PetscScalar>;
  const auto rank = gfs.gridView().grid().comm().rank();

  // Assemble the jacobian
  gridoperator.jacobian(x0, m);
  auto elapsed = timer.elapsed();
  if (rank == 0) {
    std::cout << "Jacobian assembly took: " << elapsed << "s" << std::endl;
  }
  timer.reset();

  // Evaluate residual w.r.t initial guess
  using RV = typename GridOperator::Traits::Range;
  RV r(gfs);
  r = 0.0;
  gridoperator.residual(x0, r);
  {
    // Make the residual consistent
    Dune::PDELab::AddDataHandle adddh(gfs, r);
    gfs.gridView().communicate(adddh, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
  }
  r.assemble();
  elapsed = timer.elapsed();
  if (rank == 0) {
    std::cout << "Residual assembly took: " << elapsed << "s" << std::endl;
  }
  timer.reset();

  // Create a PETSc solver
  KSP ksp;
  Mat Ais = native(m);
  Mat A; // The matrix in PETSc's parallel CRS format

  PetscCallVoid(MatConvert(Ais, MATMPIAIJ, MAT_INITIAL_MATRIX, &A));
  PetscCallVoid(MatViewFromOptions(A, nullptr, "-A_mat_view"));

  // Create parallel vectors for the solution and rhs.
  // TODO: This should all be hidden in the backend.
  ::Vec xv, bv;
  PetscCallVoid(MatCreateVecs(A, &xv, &bv));

  // Find out if we are the owner of a shared index (determined by the same rule as in the matrix setup:
  // the rank with the smallest rank number is the owner).
  Vec min_rank(gfs);
  min_rank = rank;
  min_rank.assemble();
  {
    Dune::PDELab::MinDataHandle mindh(gfs, min_rank);
    gfs.gridView().communicate(mindh, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
  }
  PetscScalar *rvdata;
  const PetscScalar *rdata;
  VecGetArrayRead(native(r), &rdata);
  VecGetArrayWrite(bv, &rvdata);
  auto cnt = 0;
  for (std::size_t i = 0; i < min_rank.size(); ++i) {
    if (min_rank[i] == rank)
      rvdata[cnt++] = rdata[i];
  }
  VecRestoreArrayWrite(bv, &rvdata);
  VecRestoreArrayRead(native(r), &rdata);

  // Create the actual solver. Note that the type of Krylov method and the preconditioner can be
  // customised from the command line. To use GMRES preconditioned with PETSc's smoothed aggregation
  // multigrid, for example, add
  //    -ksp_type gmres -pc_type gamg
  // Add the option
  //    -ksp_monitor -ksp_converged_reason
  // to view the convergence history and the reason why the method converged (or why it failed to do so).
  PetscCallVoid(KSPCreate(MPI_COMM_WORLD, &ksp));
  PetscCallVoid(KSPSetFromOptions(ksp));
  PetscCallVoid(KSPSetOperators(ksp, A, A));
  PetscCallVoid(KSPSetUp(ksp));

  // Perform the actual solve
  PetscCallVoid(KSPSolve(ksp, bv, xv));

  // Now copy the solution vector back into the local representation on the GFS
  VecGetArrayRead(xv, &rdata);
  VecGetArrayWrite(native(x0), &rvdata);
  cnt = 0;
  for (std::size_t i = 0; i < min_rank.size(); ++i) {
    if (min_rank[i] == rank)
      rvdata[i] -= rdata[cnt++];
  }
  VecRestoreArrayWrite(native(x0), &rvdata);
  VecRestoreArrayRead(xv, &rdata);

  // We only solved on a non-overlapping set of dofs. We now need to send the solution values
  // from the ranks that own the respective dof to all other ranks. We do this by setting the
  // solution to zero at dofs that we do not own and then adding up the result.
  {
    for (std::size_t i = 0; i < x0.size(); ++i)
      if (min_rank[i] != rank)
        x0[i] = 0;

    Dune::PDELab::AddDataHandle adddh(gfs, x0);
    gfs.gridView().communicate(adddh, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
    x0.assemble();
  }

  // output grid function with VTKWriter
  Dune::VTKWriter vtkwriter(gv);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gfs, x0);
  vtkwriter.addVertexData(x0, "solution");

  // PermeabilityAdapter permdgf(gv, problem); // TODO: Put this in the PDELab namespace
  // using PermVTKDGF = Dune::PDELab::VTKGridFunctionAdapter<decltype(permdgf)>;
  // vtkwriter.addCellData(std::make_shared<PermVTKDGF>(permdgf, "log(K)"));

  vtkwriter.write(filename, Dune::VTK::ascii);

  PetscCallVoid(KSPDestroy(&ksp));
  PetscCallVoid(MatDestroy(&A));
  PetscCallVoid(VecDestroy(&xv));
  PetscCallVoid(VecDestroy(&bv));
  PetscFunctionReturnVoid();
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char **argv)
{
  PetscCall(PetscInitialize(&argc, &argv, nullptr, nullptr));
  {
    try {
      // Maybe initialize Mpi
      const auto &helper = Dune::MPIHelper::instance(argc, argv);
      std::cout << "Rank " << helper.rank() << " / " << helper.size() << "\n";

      // YaspGrid Q1 2D test
      {
        // make grid
        Dune::FieldVector<double, 2> L(1.0);
        std::array<int, 2> N(Dune::filledArray<2, int>(1024));
        Dune::YaspGrid<2> grid(L, N, std::bitset<2>{0ULL}, 0); // Use non-overlapping grid, PETSc backend currently supports only that
        auto gv = grid.leafGridView();
        grid.loadBalance();

        using GridView = Dune::YaspGrid<2>::LeafGridView;
        typedef Dune::PDELab::QkLocalFiniteElementMap<GridView, typename GridView::ctype, double, 1> FEM;
        FEM fem(gv);

        // solve problem
        poisson(gv, fem, "petscbackend_yasp_Q1_2d");
      }

      // // YaspGrid Q2 2D test
      // {
      //   // make grid
      //   Dune::FieldVector<double, 2> L(1.0);
      //   std::array<int, 2> N(Dune::filledArray<2, int>(1));
      //   Dune::YaspGrid<2> grid(L, N);
      //   grid.globalRefine(3);

      //   // get view
      //   typedef Dune::YaspGrid<2>::LeafGridView GV;
      //   const GV gv = grid.leafGridView();

      //   // make finite element map
      //   typedef GV::ctype DF;
      //   typedef Dune::PDELab::QkLocalFiniteElementMap<GV, DF, double, 2> FEM;
      //   FEM fem(gv);

      //   // solve problem
      //   poisson<GV, FEM, Dune::PDELab::ConformingDirichletConstraints, Dune::PDELab::Petsc::MatrixBackend>(gv, fem, "petscbackend_yasp_Q2_2d");
      // }

      // test passed
    }
    catch (Dune::Exception &e) {
      std::cerr << "Dune reported error: " << e << std::endl;
      return 1;
    }
    catch (...) {
      std::cerr << "Unknown exception thrown!" << std::endl;
      return 1;
    }
  }

  PetscCall(PetscFinalize());
}
