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



/*
The OrderedGridFunctionSpace intends to demangle the strong relationship that
pdelab currently has between trees of grid function spaces and its ordering.
The problem is that there is not a clear separation of concerns and is very
difficult to reason about one without the other. Moreover, the construction at
compile time and initialization at run time follow different order. This make
the two of them very difficult to debug. In particular, the
OrderedGridFunctionSpace intends to perform all the interactions between the grid
function space and the ordering.

The plan for the future is that the user uses UnorderedXXXGridFunctionSpace
to build a tree of spaces and when this is finished, it pipes the tree to the
OrderedGridFunctionSpace. In other words, we want to deprecate XXXGridFunctionSpace.

*/
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
template <class UnorderedGFS, class AssemblyEntitySet>
class OrderedGridFunctionSpace
    : public UnorderedGFS
    , public impl::GridFunctionSpaceOrderingData<typename UnorderedGFS::Traits::SizeType>
    , public DataHandleProvider<OrderedGridFunctionSpace<UnorderedGFS, AssemblyEntitySet>>
{
  static_assert(isEntitySet<AssemblyEntitySet>{});

  using ordering_transformation =
      TypeTree::TransformTree<OrderedGridFunctionSpace, gfs_to_ordering<OrderedGridFunctionSpace>>;

  using GFSData = impl::GridFunctionSpaceOrderingData<typename UnorderedGFS::Traits::SizeType>;

  // the backwards compatible case requires us to access other gfs data to check initialization
  template<class,class>
  friend class ::Dune::PDELab::OrderedGridFunctionSpace;

public:
  struct Traits : public UnorderedGFS::Traits {
    using EntitySet = AssemblyEntitySet; // overload entity set if already existing
    using GridView = typename EntitySet::GridView;
    using GridViewType /*[[deprecated]]*/ = GridView;
    using Ordering = typename ordering_transformation::Type;
  };

  using OrderingTag [[deprecated]] = typename Traits::OrderingTag;

  // in this case othe order is created and updated at construction time
  OrderedGridFunctionSpace(const UnorderedGFS& ugfs, AssemblyEntitySet entity_set)
    : UnorderedGFS{ugfs}
    , _entity_set{entity_set}
  {
    if constexpr (UnorderedGFS::isLeaf)
      static_assert(std::is_same<typename UnorderedGFS::Traits::EntitySet, AssemblyEntitySet>{});

    // new behavior: we are the only node in the tree that is ordered. But
    // because we allow the whole tree to be ordered in order to be backwards
    // compatible, we need to check that when using this constructor, the whole tree is unordered
    TypeTree::forEachNode(static_cast<UnorderedGFS&>(*this),[&](auto& gfs_node, auto& path){
      if (impl::gfs_data(gfs_node))
        DUNE_THROW(GridFunctionSpaceHierarchyError,"initialized space cannot become part of larger GridFunctionSpace tree");
    });
    _is_root_space = true;
    create_ordering();
    update_ordering();
  }

  // backwards compatible constructor
  // in this case, the ordering is delayed to be constructed and updated until the ordering is called or explicitely constructed
  template<class... Args>
  OrderedGridFunctionSpace(std::in_place_t, Args&&... args)
    : UnorderedGFS(std::forward<Args>(args)...)
  {
    if constexpr (UnorderedGFS::isLeaf)
      static_assert(std::is_same<typename UnorderedGFS::Traits::EntitySet, AssemblyEntitySet>{});

    // the backwards compatible behavior: every node of this tree is also
    // ordered, thus, we need to make sure that no two of them are initialized at the same time
    TypeTree::forEachNode(*this,[&](auto& gfs_node, auto& path){
      if constexpr (std::is_base_of<GFSData,std::decay_t<decltype(gfs_node)>>{}) {
        if (gfs_node._initialized and  gfs_node._is_root_space and not gfs_node.isLeaf)
          DUNE_THROW(GridFunctionSpaceHierarchyError,"initialized space cannot become part of larger GridFunctionSpace tree");
        gfs_node._is_root_space = (path.size() == 0);
      }
    });

    const auto first_entity_set = Impl::first_leaf(static_cast<const UnorderedGFS&>(*this)).entitySet();
     _entity_set.emplace(first_entity_set);
  }


  //! Ordering tree type TODO deprecate in the future
  using Ordering = typename Traits::Ordering;

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

  typename Traits::SizeType size() const {
    if (!_initialized)
      DUNE_THROW(UninitializedGridFunctionSpaceError,
                 "space is not initialized");
    if (!_size_available)
      DUNE_THROW(GridFunctionSpaceHierarchyError,
                 "Size cannot be calculated at this point in the GFS tree.");
    return _size;
  }

  typename Traits::SizeType blockCount() const {
    if (!_initialized)
      DUNE_THROW(UninitializedGridFunctionSpaceError,
                 "space is not initialized");
    if (!_size_available)
      DUNE_THROW(
          GridFunctionSpaceHierarchyError,
          "Block count cannot be calculated at this point in the GFS tree.");
    return _block_count;
  }

  typename Traits::SizeType globalSize() const {
    if (!_initialized)
      DUNE_THROW(UninitializedGridFunctionSpaceError,
                 "space is not initialized");
    return _global_size;
  }

  //! get max dimension of shape function space
  typename Traits::SizeType maxLocalSize() const {
    if (!_initialized)
      DUNE_THROW(UninitializedGridFunctionSpaceError,
                 "space is not initialized");
    return _max_local_size;
  }

  bool isRootSpace() const
  {
    return _is_root_space;
  }

  //! get grid view
  const typename Traits::GridView& gridView () const
  {
    return entitySet().gridView();
  }

  //! get entity set
  typename Traits::EntitySet entitySet () const
  {
    return *_entity_set;
  }

  //! get entity set
  typename Traits::EntitySet entitySet ()
  {
    return *_entity_set;
  }

  //! Direct access to the DOF ordering.
  const Ordering &ordering() const { return *orderingStorage(); }

  //! Direct access to the DOF ordering.
  Ordering &ordering() { return *orderingStorage(); }

  //! Direct access to the storage of the DOF ordering.
  std::shared_ptr<const Ordering> orderingStorage() const { backwards_create_ordering(); return _ordering; }

  //! Direct access to the storage of the DOF ordering.
  std::shared_ptr<Ordering> orderingStorage() { backwards_create_ordering(); return _ordering; }

  //! Update the indexing information of the GridFunctionSpace.
  /**
   * @param force  Set to true if the underlying grid has changed (e.g.
   *      due to adaptivity) to force an update of the embedded EntitySet.
   */
  void update(bool force = false) {
    // every node in the gfs may have different entity sets
    // we only update leaf and root nodes
    auto entity_set = this->_entity_set;
    entity_set->update(force);
    TypeTree::forEachLeafNode(*this, [&](auto& gfs_node, auto& path){
      if (*entity_set != gfs_node.entitySet()) {
        gfs_node.entitySet().update(force);
        entity_set = gfs_node.entitySet();
      }
    });
    update_ordering();
  }

private:

  void update_ordering() {
    check_root_space();
    backwards_create_ordering();

    TypeTree::forEachNode(*_ordering, [&](auto& ordering_node, auto& path){
      // bool is_root = (path.size() == 0);
      if (ordering_node._gfs_data) {
        auto& data = *ordering_node._gfs_data;
        // if (data._initialized && data._is_root_space && !is_root)
        //     DUNE_THROW(GridFunctionSpaceHierarchyError,"former root space is now part of a larger tree");
        data._initialized = true;
        data._global_size = _ordering->size();
        data._max_local_size = _ordering->maxLocalSize();
        data._size_available = ordering_node.update_gfs_data_size(data._size,data._block_count);
      }
    });
  }

  void create_ordering() {
    if (_ordering)
      DUNE_THROW(GridFunctionSpaceError,
                 "Ordering can only be obtained initialized once.");
    _ordering = std::make_shared<Ordering>(ordering_transformation::transform(*this));
    // ordering is in charge of filling out our data
    _ordering->_gfs_data = static_cast<GFSData*>(this);
  }

  void backwards_create_ordering() const {
    // we do this only to be backwards compatible
    if (not _ordering) {
      const_cast<OrderedGridFunctionSpace*>(this)->create_ordering();
      const_cast<OrderedGridFunctionSpace*>(this)->update_ordering();
    }
  }

  void check_root_space() const {
    if (!this->isRootSpace()) {
      DUNE_THROW(GridFunctionSpaceHierarchyError,
                 "update() may only be called on the root of the function "
                 "space hierarchy");
    }
  }

  std::shared_ptr<Ordering> _ordering;
  std::optional<typename Traits::EntitySet> _entity_set;

  using GFSData::_size;
  using GFSData::_block_count;
  using GFSData::_global_size;
  using GFSData::_max_local_size;
  using GFSData::_is_root_space;
  using GFSData::_initialized;
  using GFSData::_size_available;
};

} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDEREDGRIDFUNCTIONSPACE_HH
