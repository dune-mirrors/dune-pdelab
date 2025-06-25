#ifndef DUNE_PDELAB_BASIS_TOPOLOGICAL_INTERLEAVING_HH
#define DUNE_PDELAB_BASIS_TOPOLOGICAL_INTERLEAVING_HH

#include <dune/pdelab/common/multiindex.hh>
#include <dune/pdelab/common/utility.hh>

#include <dune/grid/concepts/entity.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/common/typetree/nodeconcepts.hh>
#include <dune/common/typetree/traversal.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/hybridmultiindex.hh>
#include <dune/common/typetree/traversal.hh>

#include <vector>
#include <algorithm>
#include <numeric>
#include <span>
#include <memory>
#include <type_traits>
#include <utility>
#include <version>

#ifndef DUNE_PDELAB_MAX_DOFS_PER_ENTITY
#  define DUNE_PDELAB_MAX_DOFS_PER_ENTITY UCHAR_MAX
#endif

#if DUNE_PDELAB_MAX_DOFS_PER_ENTITY <= UCHAR_MAX
#  define DUNE_PDELAB_ENTITY_OFFSET_TYPE unsigned char
#elif (DUNE_PDELAB_MAX_DOFS_PER_ENTITY > UCHAR_MAX) && (DUNE_PDELAB_MAX_DOFS_PER_ENTITY <= USHRT_MAX)
#  define DUNE_PDELAB_ENTITY_OFFSET_TYPE unsigned short
#elif (DUNE_PDELAB_MAX_DOFS_PER_ENTITY > USHRT_MAX) && (DUNE_PDELAB_MAX_DOFS_PER_ENTITY <= UINT_MAX)
#  define DUNE_PDELAB_ENTITY_OFFSET_TYPE unsigned int
#elif (DUNE_PDELAB_MAX_DOFS_PER_ENTITY > UINT_MAX) && (DUNE_PDELAB_MAX_DOFS_PER_ENTITY <= ULONG_MAX)
#  define DUNE_PDELAB_ENTITY_OFFSET_TYPE unsigned long
#elif (DUNE_PDELAB_MAX_DOFS_PER_ENTITY > ULONG_MAX) && (DUNE_PDELAB_MAX_DOFS_PER_ENTITY <= ULLONG_MAX)
#  define DUNE_PDELAB_ENTITY_OFFSET_TYPE unsigned long long
#else
#  error "DUNE_PDELAB_MAX_DOFS_PER_ENTITY is too large to be represented by a standard unsigned integer type"
#endif


namespace Dune::PDELab {

 /**
 * @ingroup FunctionSpaceBases
 * @brief Represents topological interleaving of local trees for managing degrees of freedom in finite element spaces.
 * @tparam Node The type of the tree node that inherits this class
 *  (@ref ProtoBasisLeaf, @ref ProtoBasisArray, @ref ProtoBasisVector, @ref ProtoBasisTuple), using the Barton–Nackman trick for CRTP.
 *
 * This class is designed to interleave and merge degrees of freedom from product spaces of finite elements.
 * It holds the information of an opaque tree where each child node (called local tree) is associated with a topological entity of the grid view.
 * The local trees are only exposed to the user through their indices via the @ref localTreeIndices and @ref localTreeDegree methods.
 * On the other hand, this class itself forms a tree of topological interleavings that reuses local trees from
 * its children nodes to provide a new interleaving of local trees.
 *
 * The primary purpose of this class is to facilitate the indexing of containers that hold the coefficients of finite
 * element computations. It achieves this by interleaving all local trees of other BasisTopologicalInterleaving instances
 * by their common topological entity, forming a unified local tree structure below each entity of the grid view.
 *
 * Composition of local trees: topological association of 'e' on local trees is always interleaved when BasisTopologicalInterleaving instances
 * are composed in the tree
 *
 * ```
 *                        e
 *                       / \
 *                      /   \
 *    e     e    =>    0     1
 *   / \   / \        / \   / \
 *  0   1 0   1      0   1 0   1
 * ```
 *
 * Lexicographical merging: nodes of underlying children of a local tree are consecutively collapsed into the root node of the local tree.
 *
 * ```
 *       e
 *      / \               e
 *     /   \            /| |\
 *    0     1    =>    | | | |
 *   / \   / \         0 1 2 3
 *  0   1 0   1
 * ```
 *
 * The class supports both fixed-size and variable-size maps:
 * - Fixed-size maps are used when all entities of a given geometry type have the same structure, allowing for more
 *   efficient memory usage and better cache performance.
 * - Variable-size maps are used when each entity may have different sizes, requiring individual storage of offsets per entity.
 *
 * The class provides methods to query the structure of the individual local trees, such as the range of indices for a given entity
 * and the size of blocks at different levels of the hierarchy.
 * It also supports updating the local trees structure based on a new entity set.
 *
 * Key Methods:
 * - @ref localTreeIndices Provides the range of transformed indices for a given path to a node in the tree of topological interleavings.
 * - @ref localTreeDegree Determines the degree of the local tree for a given path to a node of the local trees.
 * - @ref initializeIndices Updates the indices of the interleaving based on a given grid view, using a depth-first algorithm.
 *
 * The merging strategy (i.e., `Node::isIndexBlocked`) influences how the local trees of blocks are merged:
 * - Lexicographical merging: Children nodes are merged in a lexicographical order (@ref Dune::PDELab::BasisFactory::flatByEntity() ).
 * - Blocked merging: The original blocked structure of the tree is preserved (@ref Dune::PDELab::BasisFactory::blockedByEntity() ).
 *
 * @note
 * - An entity in this class is characterized by the tuple (geometry type, entity index) which can be obtained
 *   from Dune::GlobalGeometryTypeIndex::index and the grid view's index set, respectively.
 *
 * - This class does not hold the contents of the coefficients, only its structure. The actual degrees of freedom coefficients are
 *   managed externally, and this class provides the indexing and mapping necessary to access them efficiently.
 *
 * - The multi-indices generated by this map follow a top-to-bottom ordering, meaning that the leftmost indices indicate
 *   a tree path from the upper part of a tree of indices, while the rightmost indices indicate lower parts of the tree.
 *   This ordering is consistent with the Basis in dune-functions and contrary to the GridFunctionSpace.
 *
 * @warning The mappings produced by this class are almost dense. Cases where the map is not dense occur when the blocked
 * structure induces a multi-index whose suffix has a 0-size range on the last index. This happens when there is at
 * least one blocking merging strategy within the tree and the finite elements of the different leaf nodes cover
 * different codimensions. An example is Taylor-Hood elements of degree 1 when all nodes are entity-ordered with
 * a blocked merging strategy: the root node provides the whole local tree suffix to the multi-index, but the child node
 * corresponding to the pressure element has no degrees of freedom in the facets of the elements.
 */
template<class Node>
class BasisTopologicalInterleaving
{
  // Notice that by this point `Node` is an incomplete type because it will only
  // be completed after inheriting from this class (Barton–Nackman trick)
  // Therefore, we can only inquiry its contents after completion, which
  // is why we only deduce types on functions and not in the class scope.

public:
  //! Type to represent offsets and sizes of the interleaving
  using size_type = std::size_t;
private:
  using GridViewIndexType = size_type;
  using EntityDofOffsetType = DUNE_PDELAB_ENTITY_OFFSET_TYPE;

  // key value to identify geometry types that are not mapped by this class
  static constexpr EntityDofOffsetType DOF_UNUSED = std::numeric_limits<EntityDofOffsetType>::max();
  static constexpr std::size_t GT_COUNT = GlobalGeometryTypeIndex::size(Node::GridView::dimension);

  static constexpr std::size_t BF_FIXED_SIZE_POSSIBLE = 0;
  static constexpr std::size_t BF_FIXED_SIZE = 1;
  static constexpr std::size_t BF_CODIM_OFFSET = 2;
  static constexpr std::size_t BF_GT_USED_OFFSET = BF_CODIM_OFFSET+Node::GridView::dimension+1;

  // We share private methods with other nodes of the same kind. This allows
  // to have tree algorithms in one templated class and still keep internals
  // encapsulated to the outside.
  // But with great power comes great responsability:
  //   DO NOT MODIFY OTHER NODE PRIVATE MEMBERS!
  //   Otherwise code becomes -more- unmaintainable and harder to reason about,
  //   instead create a function that describes the performed action
  template<class Node_>
  friend class BasisTopologicalInterleaving;

public:

  //! Constructs an un-initialized topological interleaving object
  BasisTopologicalInterleaving() = default;

  //! Copies the contents of another topological interleaving object (value semantics)
  BasisTopologicalInterleaving(const BasisTopologicalInterleaving&) = default;

  //! Moves the contents of another topological interleaving object
  BasisTopologicalInterleaving(BasisTopologicalInterleaving&&) = default;

  //! Copies the contents of another topological interleaving object (value semantics)
  BasisTopologicalInterleaving& operator=(const BasisTopologicalInterleaving&) = default;

  //! Moves the contents of another topological interleaving object
  BasisTopologicalInterleaving& operator=(BasisTopologicalInterleaving&&) = default;

  //! Maximum size of multi-indices
  [[nodiscard]] static constexpr std::size_t maxMultiIndexSize();

  //! Minimum size of multi-indices
  [[nodiscard]] static constexpr std::size_t minMultiIndexSize();

  //! Whether all local trees associated to entites of the same geometry type are the same at compile time
  [[nodiscard]] static constexpr auto fixedSizePerGeometryTypeStatic();

  //! Whether all local trees associated to entites of the same geometry type are the same
  [[nodiscard]] auto fixedSizePerGeometryType() const noexcept;

  //! Whether a geometry type is mapped to multi-indices within this object
  [[nodiscard]] auto containsGeometryType(size_type gt_index) const noexcept;

  //! Whether a geometry type is mapped to multi-indices within this object
  [[nodiscard]] bool containsGeometryType(const GeometryType& gt) const noexcept;

  //! Whether a codimension is mapped to multi-indices within this object
  [[nodiscard]] auto containsCodim(size_type codim) const noexcept;

  //! @brief Whether DOFs only exist in one codimension
  //! @details Note that when this is true, all the multi-indices can be
  //! constructed from the container index multi-index and the block count
  //! of a given tree (i.e. no need to inspect the other local keys of the
  //! reference local element).
  [[nodiscard]] auto singleCodim() const noexcept;

  //! @brief Whether DOFs are only attached to the volume part of the entities.
  //! @details If true, it should hint to DG and FV function spaces.
  [[nodiscard]] auto disjointCodimClosure() const noexcept;

  //! @brief Maximum number of coefficients that may be associated to local finite elements
  //! @details This is useful to know when finite elements have different sizes per entity
  [[nodiscard]] size_type maxLocalCount() const noexcept;

  //! @brief Print offsets to debug orderings
  void debugOffsets() const;

  //! Counts the number of coefficients in the leaf nodes
  [[nodiscard]] size_type dimension() const noexcept;

  //! Counts the number of coefficients in the leaf nodes of a local tree associated to an entity
  [[nodiscard]] size_type dimension(size_type gt_index, size_type entity_index) const noexcept;

  /**
   * @brief Retrieves a range of container multi-indices for a specified entity and node in the tree of topological interleavings.
   *
   * @details This function provides the range of multi-indices for a given entity and node in the tree
   * of topological interleavings, considering the merging strategy at the current level.
   *
   * @note The multi-indices have a top-to-bottom ordering, meaning that the leftmost indices indicate
   * a local tree path from the upper part of a local tree, while the rightmost indices indicate lower parts of the local tree.
   *
   * @warning The merging strategy at this level affects the result. Therefore, the range obtained from
   * `child(path).localTreeIndices(...)` will differ from `node().localTreeIndices(..., path)`.
   *
   * @tparam TreePath             Type representing the path to a node in the tree of topological interleavings.
   * @param tree_path             Path to a node in the tree of topological interleavings.
   * @param gt_index              Index of the geometry type.
   * @param entity_index          Index of the entity.
   * @return Concept::MultiIndex  Range of multi-indices for every child of the requested node.
   */
  template<class TreePath = HybridMultiIndex<> >
  [[nodiscard]] auto localTreeIndices(
    size_type gt_index,
    size_type entity_index,
    TreePath tree_path = {}) const noexcept;

  /**
   * @brief Retrieves the degree of the local tree for a given path to a node of the local trees.
   *
   * @details This function determines the size induced by a given suffix path to a node of a local tree.
   * The suffix is derived from the multi-indices generated by this object, specifically from the @ref localTreeIndices method.
   *
   * @warning The `local_tree_path` is expected to be reversed with respect to the output of `localTreeIndices`,
   * meaning that the last index in `local_tree_path` corresponds to the topmost index in the local tree path.
   *
   * For example, if a node of a local tree generates the following multi-indices for a given entity:
   * (0,0,0), (0,1,0), ..., (0,10,0), then the suffix (-,-,0) induces a size of 11, while suffixes like (-,0,0),
   * (-,1,0), etc., each induce a size of 1. This is used to appropiately resize containers that can be indexed by the local tree multi-indices.
   *
   *
   * @tparam LocalTreePath    Type representing the suffix of a multi-index.
   * @param local_tree_path   Suffix of a path to a local tree node (reversed order).
   * @param gt_index          Index of the geometry type.
   * @param entity_index      Index of the entity.
   * @return std::size_t      The degree induced by the multi-index suffix.
   */
  template<class LocalTreePath>
  [[nodiscard]] std::size_t localTreeDegree(size_type gt_index,
                                            size_type entity_index,
                                            const LocalTreePath& local_tree_path) const noexcept;

  //! Whether the local trees are merged with their children nodes
  [[nodiscard]] static constexpr auto mergedLocalTrees();


  //! Number of blocks associated to a given entity at this tree level
  [[nodiscard]] auto localTreeDegree(const size_type gt_index, const size_type entity_index) const noexcept;

  //! Number of blocks associated to a given entity (fixed size version) at this tree level
  [[nodiscard]] auto localTreeDegree(std::size_t gt_index) const noexcept;

  /**
   * @brief Updates the indices of the topological interleaving nodes with respect to a new entity set
   * @details Calling this method produces a depth first algorithm:
   * - Before every node: Allocate and set-up node.
   * - Before/After children: Carry partial calculation of the children sizes.
   * - After every node: Calculate final sizes of children and accumulate to obtain offets.
   * @warning This method only needs to be called from the final root node of the tree of topological interleavings.
   */
  void initializeIndices();

  //! @brief Priory compile-time information if a given codimension may contain DOFs wrt a grid entity
  //! @details If false, there is a guarantee that this map will never contain such a codimension.
  //! If true, codimension may or may not be used at run-time.
  template<class EntityCodim>
  [[nodiscard]] static constexpr auto mayContainCodim(EntityCodim entity_codim);

private:
  //! size of multi-indices
  [[nodiscard]] static constexpr std::size_t multiIndexSize(auto f);

  //! Total number of blocks on all entities at this tree level
  [[nodiscard]] auto blockCount() const noexcept;

  //! Compile-time information of maximum codimensions used by this map
  //! @note Less codimenions may be actually used at run-time.
  [[nodiscard]] static auto constexpr maxCodimCount();

  //! Run-time information of maximum codimensions used by this map
  [[nodiscard]] std::size_t codimCount() const noexcept;

  // Whether a codimension is mapped to multi-indices within this object
  [[nodiscard]] auto setContainsCodim(size_type codim) noexcept;

  //! Whether a geometry type is mapped to multi-indices within this object
  [[nodiscard]] auto setContainsGeometryType(size_type gt_index) noexcept;

  //! Setup ordering for fixed size maps
  //! @warning Recursive function, only call this from root node.
  void updateFixedSizeOrderings();

  //! Pre-allocate objects on this node related to the variable size ordering
  void allocateVariableSizeOrdering();

  //! Collect used geometry types on leaf nodes (variable size)
  template<class Entity>
  void collectLeafGeometryTypes(const Entity& entity);

  //! Collect geometry types for all nodes (variable size)
  //! @note Use after every leaf node called `collectLeafGeometryTypes`
  //! @warning Recursive function, only call this from root node.
  void collectGeometryTypes();

  //! Gather sizes for entity sizes on the leaf node from local keys (variable size)
  //! @note Use this after geometry types have been collected
  //! @note This class additionally checks if sizes change between different
  //! entities, thus, allowing us to know if we can compress variable size data
  //! into fixed size data.
  template<class Entity>
  void collectEntitySizes(const Entity& entity, std::vector<size_type>& gt_cache);

  //! Convert size vectors into offset vectors
  //! This class identifies first if fixed size is possible, then accumulates
  //! offsets in the corresponding data structure.
  //! @warning Recursive function, only call this from root node.
  void accumulateGeometryTypeOffsets();

  //! Setup ordering for variable size maps
  //! @note Recursive function, only call this from root node.
  void updateVariableSizeOrderings();

  //! If possible, use shared data from sibling node
  void useDataFromSibling(const auto& other);

  //! If possible, make all children share data among them
  //! @warning Recursive function, only call this from root node.
  void shareChildenData();

  //! Cast to const node implementation (Barton–Nackman trick)
  const Node& node() const noexcept { return static_cast<const Node&>(*this); }

  //! Cast to node implementation (Barton–Nackman trick)
  Node& node() noexcept { return static_cast<Node&>(*this); }

  //! Apply a functor to each child of this node
  void forEachChild(auto&& f) const {
    Dune::Hybrid::forEach(range(this->node().degree()), [&](auto i) {
      f(this->node().child(i), i);
    });
  }

  //! Apply a functor to each child of this node
  void forEachChild(auto&& f) {
    Dune::Hybrid::forEach(range(this->node().degree()), [&](auto i) {
      f(this->node().child(i), i);
    });
  }

  void validateEntityDofOffset(size_type offset) const noexcept {
    assert(offset <= static_cast<size_type>(std::numeric_limits<EntityDofOffsetType>::max()));
  }

  void validateGTOffset(size_type offset) const noexcept {
    assert(offset <= static_cast<size_type>(std::numeric_limits<GridViewIndexType>::max()));
  }

  [[nodiscard]] auto asEntityDofOffset(size_type offset) const noexcept {
    validateEntityDofOffset(offset);
    return static_cast<EntityDofOffsetType>(offset);
  }

  [[nodiscard]] static auto asSizeType(auto value) noexcept {
    return static_cast<size_type>(value);
  }

  /**
   *  Data for this node may be stored in fixed or variable size form:
   *
   * * Fixed Size:
   *    In this case, all entities of a given geometry type look the same.
   *    Thus, we just need to store the offset between nodes for each geometry
   *    type.
   *
   *                               ||    cell   || ... ||   vertex  || ... ||
   *    gt_index                   || *         ||     || *         || ... ||
   *                                  |                   |
   *                                  v                   v
   *    child_index                || 0 | 1 | 2 || ... || 0 | 1 | 2 || ... ||
   *    geometryTypeNodeOffsets    || 2 | 3 | 8 || ... || 1 | 2 | 3 || ... ||
   *
   * * Variable Size:
   *    In this case, every entity has potentially different sizes so we have
   *    to store each offset individully. This information is stored in the
   *    `_variable_dof_offsets`. Since grids have a different indexation for
   *    each geometry type, we store an offset for these geometry types in
   *    `_gt_offsets`       . This works similarly as the
   *    `MultipleCodimMultipleGeomTypeMapper` object in dune-grid works.
   *
   *                            || cell                          || vertex    ||
   *    geometryTypeOffsets     || *                             || *         ||
   *                               |                                |
   *                               v                                v
   *    entity_index            ||   cell0   ||   cell1   || ... ||  vertex0  ||
   *                               |            |                   |
   *                               v            v                   v
   *    child_index             || 0 | 1 | 2 || 0 | 1 | 2 || ... || 0 | 1 | 2 ||
   *    entityOffsets           || 2 | 3 | 8 || 0 | 0 | 0 || ... || 1 | 2 | 3 ||
   *
   */

  auto geometryTypeNodeOffsets() const {
    if constexpr (IsCompileTimeConstant<decltype(node().degree())>::value) {
      const std::size_t sz = GT_COUNT * std::max<std::size_t>(1, Node::degree());
      return std::span<EntityDofOffsetType, sz>(_dof_offsets.get(), sz);
    } else {
      const std::size_t sz = GT_COUNT * std::max<std::size_t>(1, node().degree());
      return std::span<EntityDofOffsetType>(_dof_offsets.get(), sz);
    }
  }

  auto geometryTypeOffsets() const {
    return std::span<GridViewIndexType, GT_COUNT + 1>(_gt_offsets.get(), GT_COUNT + 1);
  }

  auto entityOffsets() const {
    return std::span<EntityDofOffsetType>(
      _dof_offsets.get() + _entity_offsets_base,
      asSizeType(geometryTypeOffsets().back()) * std::max<std::size_t>(node().degree(), 1));
  }

  //! storage (as bits) of fixed size information, used codims, and used geometry types.
  std::vector<bool> _bitfield = std::vector<bool>(BF_GT_USED_OFFSET+GT_COUNT, false);

  // Note that a copy on this class make a shallow copy of these pointers.
  // However, this class is meant to be used with value semantics.
  // Thus, the contents of these arrays should be re-allocated after every initialization stage.
  std::shared_ptr<GridViewIndexType[]> _gt_offsets;
  std::shared_ptr<EntityDofOffsetType[]> _dof_offsets;
  // Base offset of entity data in _dof_offsets. This is non-zero only during
  // variable-size construction while geometryTypeNodeOffsets share the same buffer.
  size_type _entity_offsets_base = 0;

  std::size_t _block_count = 0;
  std::size_t _max_local_coeff_count = 0;
};

template<class Node>
constexpr std::size_t BasisTopologicalInterleaving<Node>::maxMultiIndexSize() {
  return multiIndexSize(Dune::Hybrid::max);
}

template<class Node>
constexpr std::size_t BasisTopologicalInterleaving<Node>::minMultiIndexSize() {
  return multiIndexSize(Dune::Hybrid::min);
}

template<class Node>
constexpr std::size_t BasisTopologicalInterleaving<Node>::multiIndexSize(auto f)
{
  if constexpr (TypeTree::Concept::LeafTreeNode<Node>) {
    return 1;
  } else {
    auto child_depth = [&]() {
      if constexpr (TypeTree::Concept::UniformInnerTreeNode<Node>) {
        return Node::ChildType::multiIndexSize(f);
      } else {
        static_assert(TypeTree::Concept::StaticDegreeInnerTreeNode<Node>);
        return unpackIntegerSequence(
          [f](auto... i) {
            return f(TypeTree::template Child<Node, i>::multiIndexSize(f)...);
          },
          std::make_index_sequence<Node::degree()>{});
      }
    }();
    if constexpr (Node::mergedLocalTrees())
      return child_depth;
    else
      return child_depth + 1;
  }
}

template<class Node>
constexpr auto BasisTopologicalInterleaving<Node>::fixedSizePerGeometryTypeStatic()
{
  if constexpr (TypeTree::Concept::LeafTreeNode<Node>) {
    // base case: query information from finite element map
    return Node::commonSizePerGeometryType().has_value();
  } else if constexpr (TypeTree::Concept::UniformInnerTreeNode<Node>) {
    // all children have the same type so we inherit their fixed size
    return Node::ChildType::fixedSizePerGeometryTypeStatic();
  } else {
    // make a conjunction of all the children types
    static_assert(TypeTree::Concept::StaticDegreeInnerTreeNode<Node>);
    auto unfold_children = [&](auto... i) {
      constexpr bool all_fixed_size =
        (TypeTree::template Child<Node, i>::fixedSizePerGeometryTypeStatic() && ...);
      return std::bool_constant<all_fixed_size>{};
    };
    auto indices = std::make_index_sequence<Node::degree()>{};
    return unpackIntegerSequence(unfold_children, indices);
  }
}

template<class Node>
auto BasisTopologicalInterleaving<Node>::fixedSizePerGeometryType() const noexcept
{
  if constexpr (fixedSizePerGeometryTypeStatic())
    return std::true_type{};
  else
    return bool(_bitfield[BF_FIXED_SIZE]);
}

template<class Node>
auto BasisTopologicalInterleaving<Node>::containsGeometryType(size_type gt_index) const noexcept
{
  return _bitfield[BF_GT_USED_OFFSET + gt_index];
}

template<class Node>
auto BasisTopologicalInterleaving<Node>::setContainsGeometryType(size_type gt_index) noexcept
{
  return _bitfield[BF_GT_USED_OFFSET + gt_index];
}

template<class Node>
bool BasisTopologicalInterleaving<Node>::containsGeometryType(const GeometryType& gt) const noexcept
{
  return containsGeometryType(GlobalGeometryTypeIndex::index(gt));
}

template<class Node>
auto BasisTopologicalInterleaving<Node>::containsCodim(size_type codim) const noexcept
{
  return _bitfield[BF_CODIM_OFFSET + codim];
}

template<class Node>
auto BasisTopologicalInterleaving<Node>::setContainsCodim(size_type codim) noexcept
{
  return _bitfield[BF_CODIM_OFFSET + codim];
}

template<class Node>
std::size_t BasisTopologicalInterleaving<Node>::codimCount() const noexcept
{
  return std::count(_bitfield.begin()+BF_CODIM_OFFSET, _bitfield.begin()+BF_CODIM_OFFSET+Node::GridView::dimension+1, true);
}

template<class Node>
auto BasisTopologicalInterleaving<Node>::singleCodim() const noexcept
{
  if constexpr (maxCodimCount() == 1)
    return std::true_type{};
  else
    return (codimCount() == 1);
}

template<class Node>
auto BasisTopologicalInterleaving<Node>::disjointCodimClosure() const noexcept
{
  if constexpr (maxCodimCount() == 1)
    return std::bool_constant<mayContainCodim(Indices::_0)>{};
  else
    return (codimCount() == 1 and containsCodim(0));
}

template<class Node>
typename BasisTopologicalInterleaving<Node>::size_type
BasisTopologicalInterleaving<Node>::maxLocalCount() const noexcept
{
  return _max_local_coeff_count;
}

template<class Node>
void BasisTopologicalInterleaving<Node>::debugOffsets() const {
  using FEM = typename Node::FiniteElementMap;
  constexpr std::size_t fem_dim = FEM::Traits::FiniteElement::Traits::LocalBasisType::Traits::dimDomain;
  constexpr std::size_t fem_codim = Node::GridView::dimension - fem_dim;

  std::cout << "Fixed Size per Geometry Type: " << bool(fixedSizePerGeometryType()) << "\n";
  std::cout << "Blocked: " << bool(Node::isIndexBlocked) << "\n";
  if (fixedSizePerGeometryType())
    std::cout << "  # Geometry Type, DOF Count\n";
  else
    std::cout << "  # Geometry Type, Entity Index, Acc. Block Count ==\n";

  for (std::size_t codim = fem_codim; codim <= Node::GridView::dimension; ++codim) {
    for (const auto& gt : node().gridView().indexSet().types(codim)) {
      const auto gt_index = GlobalGeometryTypeIndex::index(gt);
      if (fixedSizePerGeometryType()) {
        std::cout << "  " << gt << ", " << localTreeDegree(gt_index) << " \n";
      } else {
        auto gt_begin = std::begin(geometryTypeOffsets()) + gt_index;
        auto size = *(gt_begin+1) - *gt_begin;
        for (std::size_t i = 0; i != size; ++i)
          std::cout << "  " << gt << ", " << i << ", " << entityOffsets()[i + *gt_begin] << "\n";
      }
    }
  }
}

template<class Node>
typename BasisTopologicalInterleaving<Node>::size_type
BasisTopologicalInterleaving<Node>::dimension() const noexcept
{
  size_type count = 0;
  TypeTree::forEachLeafNode(node(), [&](auto& leaf, auto path) {
    count += leaf.blockCount();
  });
  return count;
}

template<class Node>
typename BasisTopologicalInterleaving<Node>::size_type
BasisTopologicalInterleaving<Node>::dimension(size_type gt_index, size_type entity_index) const noexcept
{
  size_type count = 0;
  TypeTree::forEachLeafNode(node(), [&](auto& leaf, auto path) {
    count += leaf.localTreeDegree(gt_index, entity_index);
  });
  return count;
}

template<class Node>
template<class TreePath>
[[nodiscard]] auto BasisTopologicalInterleaving<Node>::localTreeIndices(
  size_type gt_index,
  size_type entity_index,
  TreePath tree_path) const noexcept
{
  // Note: Multi-index is read Outer->Inner
  // only map known geometry indices
  assert(containsGeometryType(gt_index));
  if constexpr (TreePath::size() == 0) {
    // (end of recursion)
    // Simply return iota (in form of a multi-index of size 1) from 0 to the block size of this node.
    return transformedRangeView(range(localTreeDegree(gt_index, entity_index)), [](auto i){ return HybridMultiIndex(i); });
  } else {
    static_assert(not TypeTree::Concept::LeafTreeNode<Node>);
    const auto child = front(tree_path);
    // (continue recursion) get container index of the child node.
    const auto cir = node().child(child).localTreeIndices(gt_index, entity_index, pop_front(tree_path));
    if constexpr (Node::mergedLocalTrees()) {
      size_type offset = 0;
      // lexicopgraphic merging: accumulate front the offest of the (child-1)
      if (child != 0) {
        if (fixedSizePerGeometryType()) {
          offset = asSizeType(geometryTypeNodeOffsets()[gt_index * node().degree() + child - 1]);
        } else {
          const auto index =
            (asSizeType(geometryTypeOffsets()[gt_index]) + entity_index) * node().degree() + child - 1;
          offset = asSizeType(entityOffsets()[index]);
        }
      }
      return transformedRangeView(std::move(cir), [offset](auto ci){ return accumulate_front(ci, offset); });
    } else {
      // blocked merging: simply push front the child index
      return transformedRangeView(std::move(cir), [child](auto ci){ return push_front(ci, child); });
    }
  }
}

template<class Node>
constexpr auto BasisTopologicalInterleaving<Node>::mergedLocalTrees()
{
  // the blocking structure of a local space is given by the tag of its
  // children
  if constexpr (TypeTree::Concept::LeafTreeNode<Node>) {
    return std::false_type{};
  } else if constexpr (TypeTree::Concept::UniformInnerTreeNode<Node>) {
    return std::bool_constant<not Node::ChildType::isIndexBlocked>{};
  } else if constexpr (TypeTree::Concept::StaticDegreeInnerTreeNode<Node>) {
    auto unfold_children = [&](auto... i) {
      constexpr bool any_blocked =
        (TypeTree::template Child<Node, i>::isIndexBlocked || ...);
      constexpr bool all_blocked =
        (TypeTree::template Child<Node, i>::isIndexBlocked && ...);
      static_assert(all_blocked == any_blocked,
                    "All composite nodes topologically interleaved must have the same merging strategy");
      return std::bool_constant<not any_blocked>{};
    };
    auto indices = std::make_index_sequence<Node::degree()>{};
    return unpackIntegerSequence(unfold_children, indices);
  } else {
    static_assert(Dune::AlwaysFalse<Node>{}, "Not known Node Type");
  }
}

template<class Node>
template<class LocalTreePath>
std::size_t
BasisTopologicalInterleaving<Node>::localTreeDegree(size_type gt_index,
                                                         size_type entity_index,
                                                        const LocalTreePath& loca_tree_path) const noexcept
{
  // Note: Multi-index is read Inner->Outer
  // suffix wants the size for this node
  if (loca_tree_path.size() == 0)
    return node().localTreeDegree(gt_index, entity_index);

  if constexpr (TypeTree::Concept::LeafTreeNode<Node>) {
    assert(loca_tree_path.size() == 1);
    return 0; // leaf nodes have no children, so no size
  } else {

    // helper to return from any child with a dynamic child index
    auto childContainerSize = [&](std::size_t child_i, auto next_suffix) -> size_type {
      if constexpr (TypeTree::Concept::UniformInnerTreeNode<Node>) {
        return node().child(child_i).localTreeDegree(gt_index, entity_index, next_suffix);
      } else {
        static_assert(TypeTree::Concept::StaticDegreeInnerTreeNode<Node>);
        // at this point we recoverd the index, but there is no way to
        // propagate its static information outside of this function (i.e. a
        // return type that depends on the child index)
        size_type _size = std::numeric_limits<size_type>::max();
        // make a loop over all nodes and check which one matches the child
        // index
        forEachChild([&](auto& child, auto i) {
          if (i == child_i)
            _size = child.localTreeDegree(gt_index, entity_index, next_suffix);
        });
        return _size;
      }
    };

    // the next index to find out its size
    auto back_index = back(loca_tree_path);
    // task: find child the child node for whom this index corresponds
    if constexpr (Node::mergedLocalTrees()) {
      // here we need to "recover" the child index that describes the
      // back_index (inverse of localTreeIndices operation)
      if (not containsGeometryType(gt_index))
        return 0;
      auto dof_begin = std::begin(fixedSizePerGeometryType() ? geometryTypeNodeOffsets() : entityOffsets());
      if (fixedSizePerGeometryType())
        dof_begin += gt_index * node().degree();
      else
        dof_begin += (asSizeType(geometryTypeOffsets()[gt_index]) + entity_index) * node().degree();
      auto dof_end = dof_begin + node().degree();
      auto dof_it = std::upper_bound(dof_begin, dof_end, asEntityDofOffset(back_index));
      auto next = accumulate_back(loca_tree_path, size_type{ 0 });
      if (dof_it != dof_begin) {
        std::advance(dof_it, -1);
        const auto dof_offset = asSizeType(*dof_it);
        assert(back(loca_tree_path) >= dof_offset);
        next = accumulate_back(next, -dof_offset);
      }
      std::size_t child_index = std::distance(dof_begin, dof_it);
      assert(node().degree() > child_index);
      return childContainerSize(child_index, next);
    } else {
      // easy case, the back_index is exactly the index of the child node
      return childContainerSize(back_index, pop_back(loca_tree_path));
    }
  }
}

template<class Node>
void BasisTopologicalInterleaving<Node>::initializeIndices()
{
  updateFixedSizeOrderings();
  if (not fixedSizePerGeometryTypeStatic())
    updateVariableSizeOrderings();

  shareChildenData();
}

template<class Node>
auto BasisTopologicalInterleaving<Node>::localTreeDegree(const size_type gt_index, const size_type entity_index) const noexcept
{
  if constexpr (TypeTree::Concept::LeafTreeNode<Node>) {
    if constexpr (fixedSizePerGeometryTypeStatic())
      return std::integral_constant<size_type, Node::commonSizePerGeometryType().value()>();
    else if (fixedSizePerGeometryType()) {
      return localTreeDegree(gt_index);
    } else {
      auto gt_offset = asSizeType(geometryTypeOffsets()[gt_index]) + entity_index;
      return gt_offset < entityOffsets().size() ? asSizeType(entityOffsets()[gt_offset]) : 0;
    }
  } else if constexpr (Node::mergedLocalTrees()) {
    if (fixedSizePerGeometryType()) {
      assert(containsGeometryType(gt_index));
      return size_type{ localTreeDegree(gt_index) };
    }

    auto gt_offset = asSizeType(geometryTypeOffsets()[gt_index]) + entity_index;
    const auto degree = node().degree();
    auto gt_i = gt_offset * degree + degree - 1;
    return gt_i < entityOffsets().size() ? asSizeType(entityOffsets()[gt_i]) : 0;
  } else {
    return node().degree();
  }
}

template<class Node>
auto BasisTopologicalInterleaving<Node>::localTreeDegree(std::size_t gt_index) const noexcept
{
  assert(fixedSizePerGeometryType());
  if constexpr (TypeTree::Concept::LeafTreeNode<Node>) {
    if constexpr (fixedSizePerGeometryTypeStatic())
      return std::integral_constant<size_type, Node::commonSizePerGeometryType().value()>();
    else
      return asSizeType(geometryTypeNodeOffsets()[gt_index]);
  } else if constexpr (Node::mergedLocalTrees()) {
    const auto degree = node().degree();
    return asSizeType(geometryTypeNodeOffsets()[gt_index * degree + degree - 1]);
  } else {
    return node().degree();
  }
}

template<class Node>
auto BasisTopologicalInterleaving<Node>::blockCount() const noexcept
{
  if constexpr (Node::isIndexBlocked and not TypeTree::Concept::LeafTreeNode<Node>) {
    return node().degree();
  } else {
    return _block_count;
  }
}


template<class Node>
template<class EntityCodim>
constexpr auto BasisTopologicalInterleaving<Node>::mayContainCodim(EntityCodim entity_codim)
{
  if constexpr (TypeTree::Concept::LeafTreeNode<Node>) {
    constexpr std::size_t entity_dim = Node::GridView::dimension - entity_codim;
    // dimension of the finite element domain (perhaps embedded on sub entities)
    using FEM = typename Node::FiniteElementMap;
    constexpr std::size_t fem_dim = FEM::Traits::FiniteElement::Traits::LocalBasisType::Traits::dimDomain;
    if constexpr (entity_dim > fem_dim) {
      // requested dim is higher than the finite element dimension thus not contained
      return std::false_type{};
    } else if constexpr (requires { FEM::hasDOFs(int{}); }) {
      // in this case we can check if a DOF is included or not
      constexpr std::size_t fem_codim = fem_dim - entity_dim;
      constexpr bool has_dofs = FEM::hasDOFs(static_cast<int>(fem_codim));
      return std::bool_constant<has_dofs>{};
    } else {
      return std::true_type{};
    }
  } else if constexpr (TypeTree::Concept::UniformInnerTreeNode<Node>) {
    return Node::ChildType::mayContainCodim(entity_codim);
  } else if constexpr (TypeTree::Concept::StaticDegreeInnerTreeNode<Node>) {
    auto unfold_children = [&](auto... i) {
      constexpr bool has_dofs =
        (TypeTree::template Child<Node, i>::mayContainCodim(entity_codim) || ...);
      return std::bool_constant<has_dofs>{};
    };
    auto indices = std::make_index_sequence<Node::degree()>{};
    return unpackIntegerSequence(unfold_children, indices);
  } else {
    static_assert(Dune::AlwaysFalse<Node>{}, "Not known Node Type");
  }
}

template<class Node>
auto constexpr BasisTopologicalInterleaving<Node>::maxCodimCount()
{
  auto sequence = std::make_index_sequence<1 + Node::GridView::dimension>{};

  constexpr std::size_t count = Dune::unpackIntegerSequence(
    [&](auto... codim) {
      constexpr std::size_t one{ 1 }, zero{ 0 };
      return ((BasisTopologicalInterleaving::mayContainCodim(codim) ? one : zero) + ...);
    },
    sequence);
  return std::integral_constant<std::size_t, count>{};
}

template<class Node>
void BasisTopologicalInterleaving<Node>::updateFixedSizeOrderings()
{
  const std::size_t dim = Node::GridView::dimension;

  _bitfield[BF_FIXED_SIZE] = fixedSizePerGeometryTypeStatic();
  _max_local_coeff_count = 0;
  _block_count = 0;

  if constexpr (fixedSizePerGeometryTypeStatic()) {
    // reset falgs and offsets
    for (std::size_t codim = 0; codim <= dim; ++codim)
      setContainsCodim(codim) = false;
    for (std::size_t gt = 0; gt < GT_COUNT; ++gt)
      setContainsGeometryType(gt) = false;
    std::size_t sz = GT_COUNT * std::max<std::size_t>(1, node().degree());
    _dof_offsets = std::unique_ptr<EntityDofOffsetType[]>(new EntityDofOffsetType[sz]);
    std::uninitialized_fill_n(_dof_offsets.get(), sz, EntityDofOffsetType{0});
    _gt_offsets.reset();
    _entity_offsets_base = 0;
  }

  // fill out flags and offsets depending on the node type
  if constexpr (TypeTree::Concept::LeafTreeNode<Node>) {

    using FEM = typename Node::FiniteElementMap;
    constexpr std::size_t fem_dim = FEM::Traits::FiniteElement::Traits::LocalBasisType::Traits::dimDomain;
    constexpr std::size_t fem_codim = Node::GridView::dimension - fem_dim;

    if constexpr (fixedSizePerGeometryTypeStatic()) {
      for (std::size_t codim = fem_codim; codim <= Node::GridView::dimension; ++codim) {
        for (const auto& gt : node().gridView().indexSet().types(codim)) {
          size_type size = node().finiteElementMap().size(gt);
          const auto gt_index = GlobalGeometryTypeIndex::index(gt);
          geometryTypeNodeOffsets()[gt_index] = asEntityDofOffset(size);
          setContainsGeometryType(gt_index) = size > 0;
          _block_count += size * node().gridView().size(gt);
          assert(codim == dim - gt.dim());
          setContainsCodim(codim) = containsCodim(codim) or (size > 0);
        }
      }
      _max_local_coeff_count = node().finiteElementMap().maxLocalSize();
    }
  } else {
    forEachChild([&](auto& child, auto i) {
      // first, set up child gt collection
      child.updateFixedSizeOrderings();
      // then, accumulate child results to this node
      if constexpr (fixedSizePerGeometryTypeStatic()) {
        // this node can only be fixed size if child is also fixed size
        assert(child.fixedSizePerGeometryTypeStatic());

        // properties contained in child nodes are also contained here
        for (std::size_t codim = 0; codim <= dim; ++codim)
          setContainsCodim(codim) = containsCodim(codim) or child.containsCodim(codim);

        for (std::size_t gt = 0; gt < GT_COUNT; ++gt)
          setContainsGeometryType(gt) = containsGeometryType(gt) or child.containsGeometryType(gt);

        _max_local_coeff_count += child.maxLocalCount();

        if (not Node::isIndexBlocked)
          _block_count += child.blockCount();

        // get size of child nodes
        for (std::size_t gt = 0; gt < GT_COUNT; ++gt)
          if (child.containsGeometryType(gt))
            geometryTypeNodeOffsets()[gt * node().degree() + i] = asEntityDofOffset(child.localTreeDegree(gt));
      }
    });
  }

  if constexpr (fixedSizePerGeometryTypeStatic()) {
    // finally, convert child gt sizes into offsets
    const auto advance = std::max<std::size_t>(1, node().degree());
    for (std::size_t gt = 0; gt < GT_COUNT; ++gt) {
      size_type carry = 0;
      for (std::size_t i = 0; i < advance; ++i) {
        const auto offset_i = gt * advance + i;
        carry += asSizeType(geometryTypeNodeOffsets()[offset_i]);
        geometryTypeNodeOffsets()[offset_i] = asEntityDofOffset(carry);
      }
    }
  }
}

template<class Node>
void BasisTopologicalInterleaving<Node>::allocateVariableSizeOrdering()
{
  static_assert(not fixedSizePerGeometryTypeStatic());
  // reset flags and offsets
  for (std::size_t codim = 0; codim <= Node::GridView::dimension; ++codim)
    setContainsCodim(codim) = false;
  for (std::size_t gt = 0; gt < GT_COUNT; ++gt)
    setContainsGeometryType(gt) = false;

  std::size_t sz = GT_COUNT * std::max<std::size_t>(1, node().degree());
  _dof_offsets = std::unique_ptr<EntityDofOffsetType[]>(new EntityDofOffsetType[sz]);
  _gt_offsets = std::unique_ptr<GridViewIndexType[]>(new GridViewIndexType[GT_COUNT + 1]);
  _entity_offsets_base = 0;
  std::uninitialized_fill_n(geometryTypeNodeOffsets().begin(), geometryTypeNodeOffsets().size(), DOF_UNUSED);
  std::uninitialized_fill_n(geometryTypeOffsets().begin(), geometryTypeOffsets().size(), GridViewIndexType{0});
}

template<class Node>
template<class Entity>
void BasisTopologicalInterleaving<Node>::collectLeafGeometryTypes(const Entity& entity)
{
  static_assert(not fixedSizePerGeometryTypeStatic());
  static_assert(TypeTree::Concept::LeafTreeNode<Node>);
  static_assert(Dune::Concept::EntityExtended<Entity>);
  assert(not fixedSizePerGeometryType());

  using FEM = typename Node::FiniteElementMap;
  constexpr std::size_t fem_dim = FEM::Traits::FiniteElement::Traits::LocalBasisType::Traits::dimDomain;
  constexpr std::size_t fem_codim = Node::GridView::dimension - fem_dim;

  const FEM& fem = node().finiteElementMap();
  std::size_t max_coeff_count = 0;
  for (const auto& sub_entity : subEntities(entity, Dune::Codim<fem_codim>{})) {
    const auto& finite_element = fem.find(sub_entity);
    using FEM = typename Node::FiniteElementMap;
    using FESwitch = FiniteElementInterfaceSwitch<typename FEM::Traits::FiniteElement>;
    const auto& coeffs = FESwitch::coefficients(finite_element);
    max_coeff_count += coeffs.size();

    const auto& ref_el = referenceElement(sub_entity.geometry());
    for (std::size_t i = 0; i < coeffs.size(); ++i) {
      const auto& key = coeffs.localKey(i);
      auto gt_index = GlobalGeometryTypeIndex::index(ref_el.type(key.subEntity(), key.codim()));
      setContainsGeometryType(gt_index) = true;
      setContainsCodim(fem_codim + key.codim()) = true;
    }
  }
  _max_local_coeff_count = std::max<std::size_t>(_max_local_coeff_count, max_coeff_count);
}

template<class Node>
void BasisTopologicalInterleaving<Node>::collectGeometryTypes()
{
  _bitfield[BF_FIXED_SIZE_POSSIBLE] = true;

  if constexpr (not TypeTree::Concept::LeafTreeNode<Node>) {
    forEachChild([&](auto& child, auto i) {
      child.collectGeometryTypes();

      if constexpr (not fixedSizePerGeometryTypeStatic()) {
        // properties contained in child nodes are also contained here
        for (std::size_t codim = 0; codim <= Node::GridView::dimension; ++codim)
          setContainsCodim(codim) = containsCodim(codim) or child.containsCodim(codim);

        for (std::size_t gt = 0; gt < GT_COUNT; ++gt)
          setContainsGeometryType(gt) = containsGeometryType(gt) or child.containsGeometryType(gt);
      }
    });
  }

  // create offset of indices for contained geometry types
  if constexpr (not fixedSizePerGeometryTypeStatic()) {
    for (std::size_t codim = 0; codim <= Node::GridView::dimension; ++codim) {
      for (const auto& gt : node().gridView().indexSet().types(codim)) {
        const auto gt_index = GlobalGeometryTypeIndex::index(gt);
        if (containsGeometryType(gt_index)) {
          const auto gt_size = node().gridView().indexSet().size(gt);
          validateGTOffset(gt_size);
          geometryTypeOffsets()[gt_index + 1] = static_cast<GridViewIndexType>(gt_size);
        }
      }
    }
    for (std::size_t i = 1; i < geometryTypeOffsets().size(); ++i) {
      const size_type sum = asSizeType(geometryTypeOffsets()[i - 1]) + asSizeType(geometryTypeOffsets()[i]);
      validateGTOffset(sum);
      geometryTypeOffsets()[i] = static_cast<GridViewIndexType>(sum);
    }
    std::size_t sz = asSizeType(geometryTypeOffsets().back()) * std::max<std::size_t>(node().degree(), 1);
    const auto dof_prefix = geometryTypeNodeOffsets().size();
    _dof_offsets = std::unique_ptr<EntityDofOffsetType[]>(new EntityDofOffsetType[dof_prefix + sz]);
    _entity_offsets_base = dof_prefix;
    std::uninitialized_fill_n(entityOffsets().begin(), entityOffsets().size(), EntityDofOffsetType{0});
  }
}

template<class Node>
template<class Entity>
void BasisTopologicalInterleaving<Node>::collectEntitySizes(const Entity& entity, std::vector<size_type>& gt_cache)
{
  static_assert(not fixedSizePerGeometryTypeStatic());
  static_assert(TypeTree::Concept::LeafTreeNode<Node>);
  if (_bitfield[BF_FIXED_SIZE_POSSIBLE])
    std::fill(gt_cache.begin(), gt_cache.end(), asSizeType(DOF_UNUSED));

  using FEM = typename Node::FiniteElementMap;
  constexpr std::size_t fem_dim = FEM::Traits::FiniteElement::Traits::LocalBasisType::Traits::dimDomain;
  constexpr std::size_t fem_codim = Node::GridView::dimension - fem_dim;

  const FEM& fem = node().finiteElementMap();
  for (std::size_t s = 0; s != entity.subEntities(fem_codim); ++s) {
    const auto& sub_entity = [&]{
      if constexpr (fem_codim == 0)
        return entity;
      else
        return entity.template subEntity<fem_codim>(s);
    }();

    if (not node().subDomain().contains(sub_entity)) {
      // if we skip a sub-entity, we discard fixed size optimization
      _bitfield[BF_FIXED_SIZE_POSSIBLE] = false;
      continue;
    }
    const auto& finite_element = fem.find(sub_entity);
    if (finite_element.type() != entity.type())
      DUNE_THROW(RangeError, "Dune::GeometryType of the local finite element and entity do not match!");
    using FEM = typename Node::FiniteElementMap;
    using FESwitch = FiniteElementInterfaceSwitch<typename FEM::Traits::FiniteElement>;
    const auto& coeffs = FESwitch::coefficients(finite_element);

    const auto& ref_el = referenceElement(sub_entity.geometry());
    for (std::size_t dof = 0; dof != coeffs.size(); ++dof) {
      const LocalKey& key = coeffs.localKey(dof);
      const size_type gt_index = GlobalGeometryTypeIndex::index(ref_el.type(key.subEntity(), key.codim()));
      const size_type entity_index = node().gridView().indexSet().subIndex(sub_entity, key.subEntity(), fem_codim + key.codim());
      const size_type index = asSizeType(geometryTypeOffsets()[gt_index]) + entity_index;
      const size_type next = key.index() + 1;
      validateEntityDofOffset(next);
      gt_cache[gt_index] = std::max<size_type>(asSizeType(entityOffsets()[index]), next);
      entityOffsets()[index] = asEntityDofOffset(gt_cache[gt_index]);
    }
  }

  // as long as we haven't discarded possible fixed size, we try to find
  // out if the seen geometry types have fixed size
  if (_bitfield[BF_FIXED_SIZE_POSSIBLE]) {
    for (std::size_t i = 0; i < gt_cache.size(); ++i) {
      // update unused entries
      if (geometryTypeNodeOffsets()[i] == DOF_UNUSED)
        geometryTypeNodeOffsets()[i] = asEntityDofOffset(gt_cache[i]);
      // if cache does not match global gt sizes, we need to discard fixed
      // size computations
      if (asSizeType(geometryTypeNodeOffsets()[i]) != gt_cache[i]) {
        _bitfield[BF_FIXED_SIZE_POSSIBLE] = false;
        break;
      }
    }
  }
}

template<class Node>
void BasisTopologicalInterleaving<Node>::accumulateGeometryTypeOffsets()
{
  _block_count = 0;
  if constexpr (TypeTree::Concept::LeafTreeNode<Node>) {
    // if we did't discard fixed size by this point, we can use the
    // fixed size geometry type sizes
    if (_bitfield[BF_FIXED_SIZE_POSSIBLE]) {
      // _variable_dof_offsets = nullptr; // discard individual entity sizes
      _bitfield[BF_FIXED_SIZE] = true;
    }
    // mask out GT_UNUSED for geometry types that really weren't used
    for (std::size_t codim = 0; codim <= Node::GridView::dimension; ++codim) {
      for (const auto& gt : node().gridView().indexSet().types(codim)) {
        auto& size = geometryTypeNodeOffsets()[GlobalGeometryTypeIndex::index(gt)];
        if (size == DOF_UNUSED)
          size = 0;
        if (_bitfield[BF_FIXED_SIZE])
          _block_count += asSizeType(size) * node().gridView().size(gt);
      }
    }
    if (not _bitfield[BF_FIXED_SIZE])
      _block_count = std::accumulate(std::begin(entityOffsets()), std::end(entityOffsets()), size_type{0},
                                     [](size_type sum, auto value) { return sum + static_cast<size_type>(value); });
  } else {
    // update node properties
    _bitfield[BF_FIXED_SIZE_POSSIBLE] = true;
    forEachChild([&](auto& child, auto i) {
      // update sub-tree
      child.accumulateGeometryTypeOffsets();
      _bitfield[BF_FIXED_SIZE_POSSIBLE] = _bitfield[BF_FIXED_SIZE_POSSIBLE] && child.fixedSizePerGeometryType();
      _max_local_coeff_count += child.maxLocalCount();
    });

    if (_bitfield[BF_FIXED_SIZE_POSSIBLE]) {
      // we need to update gt sizes from children and convert them to offsets
      forEachChild([&](auto& child, auto i) {
        for (std::size_t codim = 0; codim <= Node::GridView::dimension; ++codim) {
          for (const auto& gt : node().gridView().indexSet().types(codim)) {
            const auto gt_index = GlobalGeometryTypeIndex::index(gt);
            const auto gt_offset = gt_index * node().degree() + i;
            const auto cb_count = child.localTreeDegree(gt_index);
            geometryTypeNodeOffsets()[gt_offset] = asEntityDofOffset(cb_count);
            _block_count += cb_count * node().gridView().size(gt);
            if (i != 0) {
              const auto sum = asSizeType(geometryTypeNodeOffsets()[gt_offset - 1]) + cb_count;
              geometryTypeNodeOffsets()[gt_offset] = asEntityDofOffset(sum);
            }
          }
        }
      });

      _bitfield[BF_FIXED_SIZE] = true;
    } else {
      // no fixed size was possible, save offsets for each entity:
      //   for every entity we already have block count,
      //   then, we just need to carry such count to the next entity value
      size_type index = 0;
      for (size_type gt = 0; gt < GT_COUNT; ++gt) {
        if (not containsGeometryType(gt))
          continue;
        assert(geometryTypeOffsets()[gt] <= geometryTypeOffsets()[gt + 1]);
        const size_type entity_count =
          asSizeType(geometryTypeOffsets()[gt + 1]) - asSizeType(geometryTypeOffsets()[gt]);
        for (size_type e_index = 0; e_index < entity_count; ++e_index) {
          size_type carry = 0;
          forEachChild([&](auto& child, auto i) {
            if (child.containsGeometryType(gt))
              carry += child.localTreeDegree(gt, e_index);
            entityOffsets()[index++] = asEntityDofOffset(carry);
          });
          _block_count += Node::isIndexBlocked ? (carry != 0) : carry;
        }
      }
    }
  }

  // Keep prefixed layout only while collecting variable-size geometry data.
  // For regular variable-size usage, entity offsets start at the beginning.
  if (not _bitfield[BF_FIXED_SIZE] and _entity_offsets_base != 0) {
    const auto entity_count = entityOffsets().size();
    auto compact = std::unique_ptr<EntityDofOffsetType[]>(new EntityDofOffsetType[entity_count]);
    std::uninitialized_copy_n(entityOffsets().begin(), entity_count, compact.get());
    _dof_offsets = std::move(compact);
    _entity_offsets_base = 0;
  }
}

template<class Node>
void BasisTopologicalInterleaving<Node>::updateVariableSizeOrderings()
{

  TypeTree::forEachNode(node(), []<class T>(T& node, auto path) {
    if constexpr (not T::fixedSizePerGeometryTypeStatic())
      node.allocateVariableSizeOrdering();
  });

  for (const auto& entity : elements(node().gridView())) {
    TypeTree::forEachLeafNode(node(), [&]<class T>(T& leaf, auto path) {
      if constexpr (not T::fixedSizePerGeometryTypeStatic())
        leaf.collectLeafGeometryTypes(entity);
    });
  }

  collectGeometryTypes();

  { // construct cache outside entity collection to avoid reallocation
    std::vector<size_type> gt_cache(GT_COUNT, 0);

    for (const auto& entity : elements(node().gridView())) {
      TypeTree::forEachLeafNode(node(), [&]<class T>(T& leaf, auto path) {
        if constexpr (not T::fixedSizePerGeometryTypeStatic())
          leaf.collectEntitySizes(entity, gt_cache);
      });
    }
  }
  accumulateGeometryTypeOffsets();
}

template<class Node>
void BasisTopologicalInterleaving<Node>::shareChildenData() {
  // in case of vector space (same underlying function space), enable sharing states in children nodes
  if constexpr (not TypeTree::Concept::LeafTreeNode<Node>) {
    // compare other nodes to the first one
    const auto& ref_node = node().child(Dune::Indices::_0);
    forEachChild([&](auto& child, auto i) {
      child.shareChildenData();
      if (i > 0)
        child.useDataFromSibling(ref_node);
    });
  }
}

template<class Node>
void BasisTopologicalInterleaving<Node>::useDataFromSibling(const auto& other) {
  if (_bitfield[BF_FIXED_SIZE] == other._bitfield[BF_FIXED_SIZE]) {
    if (_bitfield[BF_FIXED_SIZE]
      and std::equal(geometryTypeNodeOffsets().begin(), geometryTypeNodeOffsets().end(), other.geometryTypeNodeOffsets().begin(), other.geometryTypeNodeOffsets().end())) {
      _dof_offsets = other._dof_offsets;
      _entity_offsets_base = other._entity_offsets_base;
    }
    else if (not _bitfield[BF_FIXED_SIZE]
        and std::equal(geometryTypeOffsets().begin(), geometryTypeOffsets().end(), other.geometryTypeOffsets().begin(), other.geometryTypeOffsets().end())
        and std::equal(entityOffsets().begin(), entityOffsets().end(), other.entityOffsets().begin(), other.entityOffsets().end()))
    {
      _gt_offsets = other._gt_offsets;
      _dof_offsets = other._dof_offsets;
      _entity_offsets_base = other._entity_offsets_base;
    }
  }
}

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_BASIS_TOPOLOGICAL_INTERLEAVING_HH
