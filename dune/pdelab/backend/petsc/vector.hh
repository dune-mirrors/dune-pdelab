#pragma once

#include <dune/pdelab/backend/common/tags.hh>
#include <dune/pdelab/backend/common/uncachedvectorview.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/petsc/descriptors.hh>

#include <petscerror.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include <petscvec.h>

#include <cassert>

namespace Dune::PDELab::Petsc {
template <class ElementType>
struct VecEntryRef {
  using E = ElementType;

  VecEntryRef(Vec v, PetscInt idx, bool &need_assembly) : v{v}, idx{idx}, need_assembly{&need_assembly}
  {
    const PetscScalar *data;
    PetscFunctionBegin;
    if (need_assembly) {
      PetscCallAbort(MPI_COMM_SELF, VecAssemblyBegin(v));
      PetscCallAbort(MPI_COMM_SELF, VecAssemblyEnd(v));
      need_assembly = false;
    }
    PetscCallAbort(MPI_COMM_SELF, VecGetArrayRead(v, &data));
    current_value = data[idx];
    PetscCallAbort(MPI_COMM_SELF, VecRestoreArrayRead(v, &data));
    PetscFunctionReturnVoid();
  }

  VecEntryRef &operator=(PetscScalar x)
  {
    PetscFunctionBegin;
    *need_assembly = true;
    PetscCallAbort(MPI_COMM_SELF, VecSetValue(v, idx, x, INSERT_VALUES));
    current_value = x;
    PetscFunctionReturn(*this);
  }

  VecEntryRef &operator+=(PetscScalar x)
  {
    PetscFunctionBegin;
    *need_assembly = true;
    PetscCallAbort(MPI_COMM_SELF, VecSetValue(v, idx, x, ADD_VALUES));
    current_value += x;
    PetscFunctionReturn(*this);
  }

  VecEntryRef &operator|=(bool b)
  {
    PetscFunctionBegin;
    *this += b;
    PetscFunctionReturn(*this);
  }

  operator E() const { return current_value; }

  bool operator<(const VecEntryRef &other) const { return current_value < other.current_value; }

  PetscScalar current_value{-1};
  Vec v{nullptr};
  PetscInt idx{-1};
  bool *need_assembly{nullptr};
};

template <class Container, class LFSCache>
struct PetscVecView {
  using E = typename Container::E;
  using ElementType = typename Container::E;
  using ContainerIndex = typename Container::ContainerIndex;

  PetscVecView() = default;
  PetscVecView(Container &container) : c(&container)
  {
    PetscFunctionBegin;
    const PetscScalar *vec_data;
    container.assemble();
    PetscCallVoid(VecGetArrayRead(container.vec, &vec_data));
    data.resize(container.size());
    std::copy_n(vec_data, data.size(), data.begin());
    PetscCallVoid(VecRestoreArrayRead(container.vec, &vec_data));

    c->locked = true;
    PetscFunctionReturnVoid();
  }

  ~PetscVecView() { detach(); }

  void bind(const LFSCache &cache) { lfs_cache = &cache; }
  void unbind() { lfs_cache = nullptr; }

  void attach(Container &container)
  {
    c = &container;
    PetscFunctionBegin;
    const PetscScalar *vec_data;
    container.assemble();
    PetscCallVoid(VecGetArrayRead(container.vec, &vec_data));
    data.resize(container.size());
    std::copy_n(vec_data, data.size(), data.begin());
    PetscCallVoid(VecRestoreArrayRead(container.vec, &vec_data));

    c->locked = true;
    PetscFunctionReturnVoid();
  }

  void detach()
  {
    write_back();
    data.clear();
    if (c) {
      c = nullptr;
    }
  }

  const LFSCache &cache() const { return *lfs_cache; }
  std::size_t size() const { return lfs_cache->size(); }

  template <typename LC>
  void write(const LC &local_container)
  {
    assert(lfs_cache);
    assert(c);
    assert(c->vec);
    c->locked = true;
    for (std::size_t i = 0; i < lfs_cache->size(); ++i) {
      data[lfs_cache->containerIndex(i)[0]] = accessBaseContainer(local_container)[i];
    }
  }

  template <typename LC>
  void add(const LC &local_container)
  {
    assert(lfs_cache);
    assert(c);
    assert(c->vec);
    c->locked = true;
    for (std::size_t i = 0; i < lfs_cache->size(); ++i) {
      data[lfs_cache->containerIndex(i)[0]] += accessBaseContainer(local_container)[i];
    }
  }

  template <typename ChildLFS, typename LC>
  void write_sub_container(const ChildLFS &child_lfs, const LC &local_container)
  {
    assert(lfs_cache);
    assert(c);
    assert(c->vec);
    c->locked = true;
    for (std::size_t i = 0; i < child_lfs.size(); ++i) {
      const auto local_index = child_lfs.localIndex(i);
      data[lfs_cache->containerIndex(local_index)[0]] = accessBaseContainer(local_container)[i];
    }
  }

  template <typename ChildLFS, typename LC>
  void add_sub_container(const ChildLFS &child_lfs, const LC &local_container)
  {
    assert(lfs_cache);
    assert(c);
    assert(c->vec);
    c->locked = true;
    for (std::size_t i = 0; i < child_lfs.size(); ++i) {
      const auto local_index = child_lfs.localIndex(i);
      data[lfs_cache->containerIndex(local_index)[0]] += accessBaseContainer(local_container)[i];
    }
  }

  E &operator[](std::size_t i)
  {
    c->locked = true;
    return data[lfs_cache->containerIndex(i)[0]];
  }

  E &operator[](const ContainerIndex &ci)
  {
    c->locked = true;
    return (*this)[ci[0]];
  }

  const E &operator[](std::size_t i) const { return data[lfs_cache->containerIndex(i)[0]]; }
  const E &operator[](const ContainerIndex &ci) const { return (*this)[ci[0]]; }

  void commit() { write_back(); } // TODO: This is called quite often, can we avoid calling write_back here?

  Container &container() { return *c; }

  void write_back()
  {
    if (c && c->vec) {
      c->locked = false;
      PetscScalar *vec_data;
      PetscCallVoid(VecGetArrayWrite(c->vec, &vec_data));
      std::copy(data.begin(), data.end(), vec_data);
      PetscCallVoid(VecRestoreArrayWrite(c->vec, &vec_data));
    }
  }

  const LFSCache *lfs_cache = nullptr;
  Container *c = nullptr;

  std::vector<E> data;
};

template <typename GFS, typename ET>
class VectorContainer : public Backend::impl::Wrapper<::Vec> {
public:
  using Container = ::Vec;

private:
  friend Backend::impl::Wrapper<Container>;

public:
  using GridFunctionSpace = GFS;
  using ElementType = ET; // TODO: It is a bit sketchy to have this as a template parameter because a PETSc vector always holds values of type PetscScalar. Maybe we should warn if ET != PetscScalar?
  using E = ET;
  using size_type = std::size_t;

  using ContainerIndex = typename GFS::Ordering::Traits::ContainerIndex;

  template <class C, class L>
  friend struct PetscVecView;

  template <typename LFSCache>
  using LocalView = PetscVecView<VectorContainer, LFSCache>;

  template <typename LFSCache>
  using ConstLocalView = ConstUncachedVectorView<const VectorContainer, LFSCache>;

  VectorContainer(const VectorContainer &rhs) : gfs{rhs.gfs}
  {
    PetscFunctionBegin;
    PetscCallVoid(VecDuplicate(rhs.vec, &vec));
    PetscCallVoid(VecCopy(rhs.vec, vec));
    PetscFunctionReturnVoid();
  }

  VectorContainer(const GFS &gfs, Backend::attached_container = Backend::attached_container()) : gfs{gfs}
  {
    PetscFunctionBegin;
    PetscCallVoid(VecCreateSeq(MPI_COMM_SELF, gfs.ordering().blockCount(), &vec));
    PetscFunctionReturnVoid();
  }

  VectorContainer(const GFS &gfs, Backend::unattached_container) : gfs{gfs} {}

  VectorContainer(const GFS &gfs, const E &e) : gfs{gfs}
  {
    PetscFunctionBegin;
    PetscCallVoid(VecCreateSeq(MPI_COMM_SELF, gfs.ordering().blockCount(), &vec));
    PetscCallVoid(VecSet(vec, e));
    PetscFunctionReturnVoid();
  }

  void assemble()
  {
    PetscFunctionBegin;
    if (need_assembly) {
      PetscCallVoid(VecAssemblyBegin(vec));
      PetscCallVoid(VecAssemblyEnd(vec));
      need_assembly = false;
    }
    PetscFunctionReturnVoid();
  }

  VectorContainer &operator=(const VectorContainer &rhs)
  {
    assert(!locked);
    PetscFunctionBegin;
    if (this == &rhs)
      PetscFunctionReturn(*this);

    if (vec) {
      PetscCallAbort(MPI_COMM_SELF, VecCopy(rhs.vec, vec));
    }
    else {
      PetscCallAbort(MPI_COMM_SELF, VecDuplicate(rhs.vec, &vec));
      PetscCallAbort(MPI_COMM_SELF, VecCopy(rhs.vec, vec));
    }
    PetscFunctionReturn(*this);
  }

  VectorContainer &operator=(const E &e)
  {
    assert(!locked);
    PetscFunctionBegin;
    assert(vec);
    PetscCallAbort(MPI_COMM_SELF, VecSet(vec, e));
    PetscFunctionReturn(*this);
  }

  E operator[](std::size_t i) const
  {
    PetscFunctionBegin;
    PetscScalar x;
    const PetscScalar *data;
    if (need_assembly) {
      PetscCallAbort(MPI_COMM_SELF, VecAssemblyBegin(vec));
      PetscCallAbort(MPI_COMM_SELF, VecAssemblyEnd(vec));
      need_assembly = false;
    }
    PetscCallAbort(MPI_COMM_SELF, VecGetArrayRead(vec, &data));
    x = data[i];
    PetscCallAbort(MPI_COMM_SELF, VecRestoreArrayRead(vec, &data));
    PetscFunctionReturn(x);
  }

  E operator[](const ContainerIndex &ci) const { return (*this)[ci[0]]; }

  VecEntryRef<E> operator[](std::size_t i)
  {
    assert(!locked);
    return VecEntryRef<E>(vec, i, need_assembly);
  }

  VecEntryRef<E> operator[](const ContainerIndex &ci)
  {
    assert(!locked);
    return (*this)[ci[0]];
  }

  ~VectorContainer()
  {
    PetscFunctionBegin;
    if (vec) {
      PetscCallVoid(VecDestroy(&vec));
    }
    PetscFunctionReturnVoid();
  }

  std::size_t size() const
  {
    PetscFunctionBegin;
    PetscInt s;
    PetscCallAbort(MPI_COMM_SELF, VecGetSize(vec, &s));
    PetscFunctionReturn(s);
  }

  VectorContainer &operator-=(const VectorContainer &rhs)
  {
    assert(!locked);
    PetscFunctionBegin;
    PetscCallAbort(MPI_COMM_SELF, VecAXPY(vec, -1, rhs.vec));
    PetscFunctionReturn(*this);
  }

  VectorContainer &operator+=(const VectorContainer &rhs)
  {
    assert(!locked);
    PetscFunctionBegin;
    PetscCallAbort(MPI_COMM_SELF, VecAXPY(vec, 1, rhs.vec));
    PetscFunctionReturn(*this);
  }

  VectorContainer &operator*=(double a)
  {
    assert(!locked);
    PetscFunctionBegin;
    PetscCallAbort(MPI_COMM_SELF, VecScale(vec, a));
    PetscFunctionReturn(*this);
  }

  void detach() { assert(false && "detach"); }

private:
  bool locked = false; // If we're currently operating on a view of the vector, we "lock" it, to ensure that the view is up-to-date

  Container &native() { return vec; }
  const Container &native() const { return vec; }

  const GFS &gfs;
  Vec vec{nullptr};

  mutable bool need_assembly = false;
};

template <typename GFS, typename E>
struct PetscVectorSelectorHelper {
  using Type = VectorContainer<GFS, E>;
};
} // namespace Dune::PDELab::Petsc

namespace Dune::PDELab::Backend::impl {
template <typename GFS, typename E>
struct BackendVectorSelectorHelper<Dune::PDELab::Petsc::VectorBackend, GFS, E> : public Petsc::PetscVectorSelectorHelper<GFS, E> {};
} // namespace Dune::PDELab::Backend::impl
