// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDEREDGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDEREDGRIDFUNCTIONSPACE_HH

#include <dune/pdelab/gridfunctionspace/datahandleprovider.hh>
#include <dune/pdelab/constraints/common/constraintstransformation.hh>
#include <dune/pdelab/constraints/noconstraints.hh>
#include <dune/pdelab/ordering/lexicographicordering.hh>
#include <dune/pdelab/ordering/entityblockedlocalordering.hh>

namespace Dune {
namespace PDELab {

template <class GFS, class Tag> class LocalFunctionSpace;
template<class LFS, class C, bool fast> class LFSIndexCache;

//! \addtogroup GridFunctionSpace grid function space
//! \ingroup PDELab
//! \{

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
 * an ordering and a data handler for parallel communication.
 *
 * The OrderedGridFunctionSpace intends to demangle the strong relationship that
 * pdelab previously has between trees of grid function spaces and its ordering.
 * The problem is that there was not a clear separation of concerns and is very
 * difficult to reason about one without the other. Moreover, the construction
 * at compile time and initialization at run time follow different order
 * (lazy initialization). This makes the two of them very difficult to debug.
 * The OrderedGridFunctionSpace intends to perform all the interactions between
 * the grid function space and the ordering.
 *
 * This class, in particular, separates the initialization of the ordering in
 * eager and lazy initialization. If you want to understand how grid function
 * spaces interact with orderings, try to do so with the eager initialization
 * where _every_ node un the UnorderedGFS type has no ordering.
 *
 * @tparam UnorderedGFS  A grid function space that may not have an ordering
 *
 * \see Ordering
 */
template <class UnorderedGFS, class AssemblyEntitySet>
class OrderedGridFunctionSpace
    : public UnorderedGFS
    , public DataHandleProvider<OrderedGridFunctionSpace<UnorderedGFS, AssemblyEntitySet>>
{
  static_assert(
      isEntitySet<AssemblyEntitySet>{},
      "OrderedGridFunctionSpace does not accept grid views, use entity sets!");

  using ordering_transformation =
      TypeTree::TransformTree<UnorderedGFS, gfs_to_ordering<UnorderedGFS>>;

  using GFSData = impl::GridFunctionSpaceOrderingData<typename UnorderedGFS::Traits::SizeType>;

  // ensure that we can access other ordered nodes data (lazy ordering only).
  template<class,class>
  friend class ::Dune::PDELab::OrderedGridFunctionSpace;

public:
  struct Traits : public UnorderedGFS::Traits {
    using EntitySet = AssemblyEntitySet;
    using GridView = typename EntitySet::GridView;
    using GridViewType = GridView;
    using Ordering = typename ordering_transformation::Type;
  };

  using OrderingTag [[deprecated("This alias will be removed after "
                                 "PDELab 2.9. Use Traits::OrderingTag.")]] =
      typename Traits::OrderingTag;

  /**
   * @brief Construct a new Eager Ordered Grid Function Space object
   * @details This constructor initializes an ordering object eagerly.
   * Notice that when this constructor is used, this node **cannot** be used to
   * construct bigger grid function space trees. This is because it is only
   * allowed to have one ordering in the whole tree.
   *
   * @param ugfs  Unordered grid function space
   * @param entity_set  Entity set used for assembly algorithms
   */
  OrderedGridFunctionSpace(const UnorderedGFS& ugfs, AssemblyEntitySet entity_set)
    : UnorderedGFS{ugfs}
    , _entity_set{entity_set}
    , _data{std::make_shared<GFSData>()}
  {
    if constexpr (UnorderedGFS::isLeaf)
      static_assert(std::is_same<typename UnorderedGFS::Traits::EntitySet, AssemblyEntitySet>{});

    // Eager behavior: we need to check that no other node is initialized.
    TypeTree::forEachNode(static_cast<UnorderedGFS&>(*this),[&](auto& gfs_node, auto& path){
      if constexpr (HasOrdering<decltype(gfs_node)>)
        if (gfs_node.data()->_initialized)
          DUNE_THROW(GridFunctionSpaceHierarchyError,
            "Initialized space cannot become part of larger GridFunctionSpace tree");
    });
    reset_root_flag();

    // set up ordering eagerly
    create_ordering();
    update_ordering();
  }

  // backwards compatible constructor
  // in this case, the ordering is delayed to be constructed and updated until the ordering is called or explicitely constructed

  /**
   * @brief Construct a new Lazy Ordered Grid Function Space object
   * @details This constructor **does not** initialize an ordering object,
   * instead, this process is delayed until the ordering is required.
   * Notice that when this constructor is used, this node **can** be used to
   * construct bigger grid function space trees. However, only the root node
   * will be allowed to hold an ordering.
   *
   * @tparam Args  Argument types of the UnorderedGFS constructor
   * @param args   Argument types of the UnorderedGFS constructor
   */
  template<class... Args>
  OrderedGridFunctionSpace(std::in_place_t, Args&&... args)
    : UnorderedGFS(std::forward<Args>(args)...)
    , _data{std::make_shared<GFSData>()}
  {
    if constexpr (UnorderedGFS::isLeaf)
      static_assert(std::is_same<typename UnorderedGFS::Traits::EntitySet, AssemblyEntitySet>{});
    reset_root_flag();
    const auto& first_entity_set = Impl::first_leaf(static_cast<const UnorderedGFS&>(*this)).entitySet();
    _entity_set.emplace(first_entity_set);
  }


  //! Ordering tree type
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
    if (!_data->_initialized)
      lazy_create_ordering();
    if (!_data->_size_available)
      DUNE_THROW(GridFunctionSpaceHierarchyError,
                 "Size cannot be calculated at this point in the GFS tree.");
    return _data->_size;
  }

  typename Traits::SizeType blockCount() const {
    if (!_data->_initialized)
      lazy_create_ordering();
    if (!_data->_size_available)
      DUNE_THROW(
          GridFunctionSpaceHierarchyError,
          "Block count cannot be calculated at this point in the GFS tree.");
    return _data->_block_count;
  }

  typename Traits::SizeType globalSize() const {
    if (!_data->_initialized)
      lazy_create_ordering();
    return _data->_global_size;
  }

  //! get max dimension of shape function space
  typename Traits::SizeType maxLocalSize() const {
    if (!_data->_initialized)
      lazy_create_ordering();
    return _data->_max_local_size;
  }

  bool isRootSpace() const
  {
    return _data->_is_root_space;
  }

  //! get grid view
  const typename Traits::GridView& gridView () const
  {
    return entitySet().gridView();
  }

  //! get entity set
  const typename Traits::EntitySet& entitySet () const
  {
    return *_entity_set;
  }

  //! get entity set
  typename Traits::EntitySet& entitySet ()
  {
    return *_entity_set;
  }

  //! Direct access to the DOF ordering.
  const Ordering &ordering() const { return *orderingStorage(); }

  //! Direct access to the DOF ordering.
  Ordering &ordering() { return *orderingStorage(); }

  //! Direct access to the storage of the DOF ordering.
  std::shared_ptr<const Ordering> orderingStorage() const { lazy_create_ordering(); return _ordering; }

  //! Direct access to the storage of the DOF ordering.
  std::shared_ptr<Ordering> orderingStorage() { lazy_create_ordering(); return _ordering; }

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

  std::shared_ptr<GFSData> data() {
    return _data;
  }

private:

  void reset_root_flag() {
    // lazy behavior: every node of this tree may also be ordered, thus, we need
    // to make sure that no two of them are initialized at the same time.
    TypeTree::forEachNode(*this,[&](auto& gfs_node, auto& path){
      if constexpr (HasOrdering<decltype(gfs_node)>) {
        if (gfs_node.data()->_initialized and  gfs_node.data()->_is_root_space and not gfs_node.isLeaf)
          DUNE_THROW(GridFunctionSpaceHierarchyError,
            "Initialized space cannot become part of larger GridFunctionSpace tree");
        gfs_node.data()->_is_root_space = (path.size() == 0);
      }
    });
  }

  void update_ordering() {
    check_root_space();
    lazy_create_ordering();
    _ordering->update();
    TypeTree::forEachNode(*_ordering, [&](auto& ordering_node, auto& path){
      bool is_root = (path.size() == 0);
      if (ordering_node._gfs_data) {
        auto& data = *ordering_node._gfs_data;
        if (data._initialized && data._is_root_space && !is_root)
            DUNE_THROW(GridFunctionSpaceHierarchyError,
              "Former root space is now part of a larger tree");
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
    if (_data->_initialized)
      DUNE_THROW(UninitializedGridFunctionSpaceError,
                 "Space is already initialized by other node");
    _ordering = std::make_shared<Ordering>(ordering_transformation::transform(*this));
    // ordering is in charge of filling out our data
    _ordering->_gfs_data = _data;
  }

  void lazy_create_ordering() const {
    // Since the usage may come from a const object, we need to const_cast to
    // mutate ourselves. Yes, this is as bad as it sounds but we need some
    // backwards compatibility.
    if (not _ordering) {
      const_cast<OrderedGridFunctionSpace*>(this)->create_ordering();
      const_cast<OrderedGridFunctionSpace*>(this)->update_ordering();
    }
  }

  void check_root_space() const {
    if (!isRootSpace()) {
      DUNE_THROW(GridFunctionSpaceHierarchyError,
                 "update() may only be called on the root of the function "
                 "space hierarchy");
    }
  }

  std::shared_ptr<Ordering> _ordering;
  std::optional<typename Traits::EntitySet> _entity_set;
  std::shared_ptr<GFSData> _data;
};

// \}

} // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDEREDGRIDFUNCTIONSPACE_HH
