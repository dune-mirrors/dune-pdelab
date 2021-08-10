// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDEREDGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDEREDGRIDFUNCTIONSPACE_HH

#include <dune/pdelab/gridfunctionspace/datahandleprovider.hh>

namespace Dune {
namespace PDELab {

template <class GFS, class Tag> class LocalFunctionSpace;
template<class LFS, class C, bool fast> class LFSIndexCache;

namespace Impl {
/**
 * @brief Constraint type accumulator
 * @details This class is meant to accumulate the type of constraints that
 * leaf grid function spaces have. If all leaf nodes have no constraints
 * apply<GFS>() returns true, otherwise false.
 */
class HasNoConstraints {
  template <class N>
  using DynamicDegreeConcept =
      decltype((std::size_t(std::declval<N>().degree()), true));

  template <class N>
  using StaticDegreeConcept =
      decltype((std::integral_constant<std::size_t, N::degree()>{}, true));

  template <class Node, std::enable_if_t<Node::isLeaf, bool> = true>
  auto apply(const Node &, Dune::PriorityTag<4>) {
    // get constraint type from this leaf node
    using Constraints = typename Node::Traits::ConstraintsType;
    // is node constrained?
    return std::is_same_v<Constraints, NoConstraints>;
  }

  template <class Node, DynamicDegreeConcept<Node> = true>
  auto apply(const Node &node, Dune::PriorityTag<2>) {
    // types are equal so is enough to check first child
    return this->apply(node.child(0));
  }

  template <class Node, StaticDegreeConcept<Node> = true>
  auto apply(const Node &node, Dune::PriorityTag<1>) {
    auto sequence = std::make_index_sequence<std::size_t(Node::degree())>();
    return Dune::unpackIntegerSequence(
        [&](auto... indices) {
          // apply this to every child and return conjunction
          // in other words: does every child have empty constrains?
          constexpr bool result =
              (true && ... && decltype(this->apply(node.child(indices))){});
          return std::bool_constant<result>{};
        },
        sequence);
  }

  template <class Node>
  auto apply(const Node &node) {
    return apply(node, Dune::PriorityTag<5>{});
  }

public:
  template <class GFS> static constexpr bool apply() {
    return decltype(HasNoConstraints{}.apply(std::declval<GFS>())){};
  }
};
} // namespace Impl

/**
 * @brief Add ordering and data handler to an unordered grid function space
 * @details This class extends the grid function space interface to have
 * an ordering and a data handler for parallel communication. In particular,
 * this is intened to only be used on the root node.
 *
 * @tparam UnorderedGFS  A grid function space that has no ordering
 *
 * \see Ordering
 */
template <class UnorderedGFS>
class OrderedGridFunctionSpace
    : public UnorderedGFS,
      public DataHandleProvider<OrderedGridFunctionSpace<UnorderedGFS>> {
  using ordering_transformation =
      TypeTree::TransformTree<UnorderedGFS, gfs_to_ordering<UnorderedGFS>>;

public:
  //! Inherit constructor from UnorderedGFS
  using UnorderedGFS::UnorderedGFS;

  //! Ordering tree type
  using Ordering = typename ordering_transformation::Type;

  //! extract type for storing constraints
  template <typename E> struct ConstraintsContainer {
    //! \brief define Type as the Type of a container of E's
    using Type = std::conditional_t<
        Impl::HasNoConstraints::template apply<UnorderedGFS>(),
        EmptyTransformation,
        ConstraintsTransformation<typename Ordering::Traits::DOFIndex,
                                  typename Ordering::Traits::ContainerIndex,
                                  E>>;
  };


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
   * @param force  Set to true if the underlying grid has changed (e.g.
   *      due to adaptivity) to force an update of the embedded EntitySet.
   */
  void update(bool force = false) {
    // every node in the gfs may have different entity sets
    // we only update leaf and root nodes
    auto &entity_set = this->_entity_set;
    if (entity_set)
      entity_set->update(force);
    auto update_leaf_es = impl::update_leaf_entity_set{entity_set, force};
    TypeTree::applyToTree(*this, update_leaf_es);
    // We bypass the normal access using ordering() here to avoid a double
    // update if the Ordering has not been created yet.
    if (!_ordering)
      create_ordering();
    update(*_ordering);
  }

private:
  void update(Ordering &ordering) const {
    if (!this->isRootSpace()) {
      DUNE_THROW(GridFunctionSpaceHierarchyError,
                 "update() may only be called on the root of the function "
                 "space hierarchy");
    }
    ordering.update();
    using SizeType = typename Ordering::Traits::SizeType;
    auto update_visitor = impl::update_ordering_data<SizeType>{ordering};
    TypeTree::applyToTree(ordering, update_visitor);
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
