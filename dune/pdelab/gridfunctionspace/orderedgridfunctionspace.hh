// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDEREDGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDEREDGRIDFUNCTIONSPACE_HH

#include <dune/pdelab/gridfunctionspace/datahandleprovider.hh>

namespace Dune {
namespace PDELab {

template <class UnorderedGFS>
class OrderedGridFunctionSpace
    : public UnorderedGFS,
      public DataHandleProvider<OrderedGridFunctionSpace<UnorderedGFS>> {
  using ordering_transformation =
      TypeTree::TransformTree<UnorderedGFS, gfs_to_ordering<UnorderedGFS>>;

public:
  using UnorderedGFS::UnorderedGFS;

  using Ordering = typename ordering_transformation::Type;

  //! extract type for storing constraints
  template <typename E> struct ConstraintsContainer {

    using Type =
        ConstraintsTransformation<typename Ordering::Traits::DOFIndex,
                                  typename Ordering::Traits::ContainerIndex, E>;
    //! \brief define Type as the Type of a container of E's TODO: find if every
    //! node has NoConstraints
    // typedef typename std::conditional<
    //     std::is_same<typename UnorderedGFS::Traits::ConstraintsType,
    //     NoConstraints>::value, EmptyTransformation,
    //     ConstraintsTransformation<typename Ordering::Traits::DOFIndex,
    //                               typename Ordering::Traits::ContainerIndex,
    //                               E>>::type Type;
  };

  //------------------------------

  //! Direct access to the DOF ordering.
  const Ordering &ordering() const { return *orderingStorage(); }

  //! Direct access to the DOF ordering.
  Ordering &ordering() { return *orderingStorage(); }

  //! Direct access to the storage of the DOF ordering.
  std::shared_ptr<const Ordering> orderingStorage() const {
    if (!this->isRootSpace()) {
      DUNE_THROW(GridFunctionSpaceHierarchyError,
                 "Ordering can only be obtained for root space in "
                 "GridFunctionSpace tree.");
    }
    if (!_ordering) {
      create_ordering();
      this->update(*_ordering);
    }
    return _ordering;
  }

  //! Direct access to the storage of the DOF ordering.
  std::shared_ptr<Ordering> orderingStorage() {
    if (!this->isRootSpace()) {
      DUNE_THROW(GridFunctionSpaceHierarchyError,
                 "Ordering can only be obtained for root space in "
                 "GridFunctionSpace tree.");
    }
    if (!_ordering) {
      create_ordering();
      this->update(*_ordering);
    }
    return _ordering;
  }

  //! Update the indexing information of the GridFunctionSpace.
  /**
   *
   * \ param force   Set to true if the underlying grid has changed (e.g. due to
   * adaptivity) to force an update of the embedded EntitySet.
   */
  void update(bool force = false) {
    auto entity_set = this->entitySet();
    entity_set.update(force);
    // We bypass the normal access using ordering() here to avoid a double
    // update if the Ordering has not been created yet.
    if (!_ordering)
      create_ordering();
    update(*_ordering);
  }

private:
  template <typename Ordering> void update(Ordering &ordering) const {
    if (!this->isRootSpace()) {
      DUNE_THROW(GridFunctionSpaceHierarchyError,
                 "update() may only be called on the root of the function "
                 "space hierarchy");
    }
    ordering.update();
    TypeTree::applyToTree(
        ordering,
        impl::update_ordering_data<typename UnorderedGFS::Traits::SizeType>(
            ordering));
  }

  // This method here is to avoid a double update of the Ordering when the user
  // calls GFS::update() before GFS::ordering().
  void create_ordering() const {
    _ordering =
        std::make_shared<Ordering>(ordering_transformation::transform(*this));
  }

  mutable std::shared_ptr<Ordering> _ordering;
};

} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDEREDGRIDFUNCTIONSPACE_HH
