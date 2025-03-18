#ifndef DUNE_PDELAB_BASIS_ORDERING_TOPOLOGIC_ASSOCIATIVITY_FOREST_COMPOSITE_HH
#define DUNE_PDELAB_BASIS_ORDERING_TOPOLOGIC_ASSOCIATIVITY_FOREST_COMPOSITE_HH

#include <dune/pdelab/basis/topologicassociativityforest/node.hh>

#include <dune/pdelab/common/multiindex.hh>
#include <dune/pdelab/common/tree_traversal.hh>

#include <dune/pdelab/concepts/multiindex.hh>
#include <dune/pdelab/concepts/treenode.hh>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/dynamicpowernode.hh>
#include <dune/typetree/powernode.hh>
#include <dune/typetree/treepath.hh>

#include <vector>
#include <memory>

namespace Dune::PDELab::Impl
{

  /**
   * @class ArrayTopologicAssociativityForest
   * @brief Array-based implementation of a TopologicAssociativityForestNode.
   *
   * This class represents a node in a topological associativity forest that manages an array of child orderings.
   * It uses a fixed-size array to store the child orderings, providing efficient access and management of the
   * degrees of freedom associated with finite element spaces.
   *
   * @tparam MergingStrategy The merging strategy that defines how child indices are merged.
   * @tparam ChildOrdering   The type of child ordering to store. This must satisfy the Concept::TreeNode concept.
   * @tparam degree          The number of child orderings to keep in the array.
   */
  template <class MergingStrategy,
            Concept::TreeNode ChildOrdering,
            std::size_t degree>
  class ArrayTopologicAssociativityForest
      : public TypeTree::PowerNode<ChildOrdering, degree>,
        public TopologicAssociativityForestNode<
            ArrayTopologicAssociativityForest<MergingStrategy, ChildOrdering, degree>,
            MergingStrategy>
  {
    using TreeNode = TypeTree::PowerNode<ChildOrdering, degree>;
    using OrderingNode = TopologicAssociativityForestNode<
        ArrayTopologicAssociativityForest<MergingStrategy, ChildOrdering, degree>,
        MergingStrategy>;

  public:
    /**
     * @brief Constructs an array of entity orderings.
     *
     * @param storage          The storage for the node's children.
     * @param merging_strategy The merging strategy to use for this node.
     */
    ArrayTopologicAssociativityForest(typename TreeNode::NodeStorage &&storage,
                                      const MergingStrategy &merging_strategy)
        : TreeNode{std::move(storage)}, OrderingNode{merging_strategy}
    {
    }

    ArrayTopologicAssociativityForest(const ArrayTopologicAssociativityForest &) = delete;
    ArrayTopologicAssociativityForest(ArrayTopologicAssociativityForest &&) = default;

    ArrayTopologicAssociativityForest &operator=(const ArrayTopologicAssociativityForest &) = delete;
    ArrayTopologicAssociativityForest &operator=(ArrayTopologicAssociativityForest &&) = default;
  };

  /**
 * @class VectorTopologicAssociativityForest
 * @brief Dynamic array-based implementation of a TopologicAssociativityForestNode.
 *
 * This class represents a node in a topological associativity forest that manages a dynamic array of child
 * orderings. It's useful for cases scenarios where the number of children is not known at compile time.
 *
 * @tparam MergingStrategy The merging strategy that defines how child indices are merged.
 * @tparam ChildOrdering   The type of child ordering to store. This must satisfy the Concept::TreeNode concept.
 */
  template <class MergingStrategy, Concept::TreeNode ChildOrdering>
  class VectorTopologicAssociativityForest
      : public TypeTree::DynamicPowerNode<ChildOrdering>,
        public TopologicAssociativityForestNode<
            VectorTopologicAssociativityForest<MergingStrategy, ChildOrdering>,
            MergingStrategy>
  {
    using TreeNode = TypeTree::DynamicPowerNode<ChildOrdering>;
    using OrderingNode =
        TopologicAssociativityForestNode<VectorTopologicAssociativityForest<MergingStrategy, ChildOrdering>,
                                         MergingStrategy>;

  public:

    /**
     * @brief Constructs a dynamic array of entity orderings.
     *
     * @param storage          The storage for the node's children.
     * @param merging_strategy The merging strategy to use for this node.
     */
    VectorTopologicAssociativityForest(typename TreeNode::NodeStorage &&storage,
                                       const MergingStrategy &merging_strategy)
        : TreeNode{std::move(storage)}, OrderingNode{merging_strategy}
    {
      for (std::size_t i = 0; i != this->degree(); ++i)
        assert(this->childStorage(i));
    }

    VectorTopologicAssociativityForest(const VectorTopologicAssociativityForest &) = delete;
    VectorTopologicAssociativityForest(VectorTopologicAssociativityForest &&) = default;

    VectorTopologicAssociativityForest &operator=(const VectorTopologicAssociativityForest &) = delete;
    VectorTopologicAssociativityForest &operator=(VectorTopologicAssociativityForest &&) = default;
  };

  /**
   * @class TupleTopologicAssociativityForest
   * @brief Tuple-based implementation of a TopologicAssociativityForestNode.
   *
   * This class represents a node in a topological associativity forest that manages a tuple of child orderings.
   * It provides a fixed structure for managing a heterogeneous set of child orderings, each potentially of a
   * different type.
   *
   * @tparam MergingStrategy The merging strategy that defines how child indices are merged.
   * @tparam ChildOrdering   Variadic template for the types of child orderings to store. Each must satisfy the
   *                         Concept::TreeNode concept.
   */
  template <class MergingStrategy, Concept::TreeNode... ChildOrdering>
  class TupleTopologicAssociativityForest
      : public TypeTree::CompositeNode<ChildOrdering...>,
        public TopologicAssociativityForestNode<
            TupleTopologicAssociativityForest<MergingStrategy, ChildOrdering...>,
            MergingStrategy>
  {
    using TreeNode = TypeTree::CompositeNode<ChildOrdering...>;
    using OrderingNode =
        TopologicAssociativityForestNode<TupleTopologicAssociativityForest<MergingStrategy, ChildOrdering...>,
                                         MergingStrategy>;

  public:

    /**
     * @brief Constructs a tuple of entity orderings.
     *
     * @param storage          The storage for the node's children.
     * @param merging_strategy The merging strategy to use for this node.
     */
    TupleTopologicAssociativityForest(typename TreeNode::NodeStorage &&storage,
                                      const MergingStrategy &merging_strategy)
        : TreeNode{std::move(storage)}, OrderingNode{merging_strategy}
    {
    }

    TupleTopologicAssociativityForest(const TupleTopologicAssociativityForest &) = delete;
    TupleTopologicAssociativityForest(TupleTopologicAssociativityForest &&) = default;

    TupleTopologicAssociativityForest &operator=(const TupleTopologicAssociativityForest &) = delete;
    TupleTopologicAssociativityForest &operator=(TupleTopologicAssociativityForest &&) = default;
  };

} // namespace Dune::PDELab::Impl

#endif // DUNE_PDELAB_BASIS_ORDERING_TOPOLOGIC_ASSOCIATIVITY_FOREST_COMPOSITE_HH
