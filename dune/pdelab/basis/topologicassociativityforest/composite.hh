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

namespace Dune::PDELab::Impl {

/**
 * @brief Array implementation of an TopologicAssociativityForestNode
 *
 * @tparam MergingStrategy   A merging strategy between the node indices
 * @tparam ChildOrdering     Child ordering type to store
 * @tparam degree            Number of child orderings to keep
 */
template<class MergingStrategy,
         Concept::TreeNode ChildOrdering,
         std::size_t degree>
class ArrayTopologicAssociativityForest
  : public TypeTree::PowerNode<ChildOrdering, degree>
  , public TopologicAssociativityForestNode<
      ArrayTopologicAssociativityForest<MergingStrategy, ChildOrdering, degree>,
      MergingStrategy>
{
  using TreeNode = TypeTree::PowerNode<ChildOrdering, degree>;
  using OrderingNode = TopologicAssociativityForestNode<
    ArrayTopologicAssociativityForest<MergingStrategy, ChildOrdering, degree>,
    MergingStrategy>;

public:
  //! Constructs an array of entity orderings
  ArrayTopologicAssociativityForest(typename TreeNode::NodeStorage&& storage,
                      const MergingStrategy& merging_strategy)
    : TreeNode{ std::move(storage) }
    , OrderingNode{ merging_strategy }
  {
  }

  ArrayTopologicAssociativityForest(const ArrayTopologicAssociativityForest&) = delete;
  ArrayTopologicAssociativityForest(ArrayTopologicAssociativityForest&&) = default;

  ArrayTopologicAssociativityForest& operator=(const ArrayTopologicAssociativityForest&) = delete;
  ArrayTopologicAssociativityForest& operator=(ArrayTopologicAssociativityForest&&) = default;
};

template<class MergingStrategy, Concept::TreeNode ChildOrdering>
class VectorTopologicAssociativityForest
  : public TypeTree::DynamicPowerNode<ChildOrdering>
  , public TopologicAssociativityForestNode<
      VectorTopologicAssociativityForest<MergingStrategy, ChildOrdering>,
      MergingStrategy>
{
  using TreeNode = TypeTree::DynamicPowerNode<ChildOrdering>;
  using OrderingNode =
    TopologicAssociativityForestNode<VectorTopologicAssociativityForest<MergingStrategy, ChildOrdering>,
                       MergingStrategy>;

public:
  VectorTopologicAssociativityForest(typename TreeNode::NodeStorage&& storage,
                       const MergingStrategy& merging_strategy)
    : TreeNode{ std::move(storage) }
    , OrderingNode{ merging_strategy }
  {
    for(std::size_t i = 0; i != this->degree(); ++i)
      assert(this->childStorage(i));
  }

  VectorTopologicAssociativityForest(const VectorTopologicAssociativityForest&) = delete;
  VectorTopologicAssociativityForest(VectorTopologicAssociativityForest&&) = default;

  VectorTopologicAssociativityForest& operator=(const VectorTopologicAssociativityForest&) = delete;
  VectorTopologicAssociativityForest& operator=(VectorTopologicAssociativityForest&&) = default;
};

template<class MergingStrategy, Concept::TreeNode... ChildOrdering>
class TupleTopologicAssociativityForest
  : public TypeTree::CompositeNode<ChildOrdering...>
  , public TopologicAssociativityForestNode<
      TupleTopologicAssociativityForest<MergingStrategy, ChildOrdering...>,
      MergingStrategy>
{
  using TreeNode = TypeTree::CompositeNode<ChildOrdering...>;
  using OrderingNode =
    TopologicAssociativityForestNode<TupleTopologicAssociativityForest<MergingStrategy, ChildOrdering...>,
                       MergingStrategy>;

public:
  TupleTopologicAssociativityForest(typename TreeNode::NodeStorage&& storage,
                      const MergingStrategy& merging_strategy)
    : TreeNode{ std::move(storage) }
    , OrderingNode{ merging_strategy }
  {
  }

  TupleTopologicAssociativityForest(const TupleTopologicAssociativityForest&) = delete;
  TupleTopologicAssociativityForest(TupleTopologicAssociativityForest&&) = default;

  TupleTopologicAssociativityForest& operator=(const TupleTopologicAssociativityForest&) = delete;
  TupleTopologicAssociativityForest& operator=(TupleTopologicAssociativityForest&&) = default;
};

} // namespace Dune::PDELab::Impl

#endif // DUNE_PDELAB_BASIS_ORDERING_TOPOLOGIC_ASSOCIATIVITY_FOREST_COMPOSITE_HH
