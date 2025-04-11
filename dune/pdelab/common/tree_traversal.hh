#ifndef DUNE_PDELAB_COMMON_TREE_TRAVERSAL_HH
#define DUNE_PDELAB_COMMON_TREE_TRAVERSAL_HH

#include <dune/pdelab/concepts/multiindex.hh>
#include <dune/pdelab/concepts/treenode.hh>

#include <dune/typetree/treepath.hh>
#include <dune/typetree/traversal.hh>

#include <dune/common/indices.hh>

#include <utility>
#include <execution>

namespace Dune::PDELab {

/**
 * @brief Traverse each entry of a tree container
 * @details Functors may accept one or two values. The first one being the child
 * node to be evaluated, and the second one the index of such entry
 *
 * @tparam Container          Tree container to be traversed
 * @tparam Callable           Functor to be applied at each child
 */
template<class Container, class Callable>
requires Concept::ParentTreeNode<std::remove_cvref_t<Container>>
constexpr void forEachChild(Container&& container, Callable&& at_value)
{
  auto invoke = [&at_value]<class Value, class Index>(Value&& value, Index index){
    static_assert(std::invocable<Callable&&, Value&&, Index> || std::invocable<Callable&&, Value&&>);
    if constexpr (std::invocable<Callable&&, Value&&, Index>)
      at_value(std::forward<Value>(value), index);
    else
      at_value(std::forward<Value>(value));
  };

  if constexpr (Concept::TupleTreeNode<Container>)
    Dune::unpackIntegerSequence(
      [&](auto... i) { (invoke(container.child(i), i), ...); },
      std::make_index_sequence<std::remove_cvref_t<Container>::degree()>{});
  else
    for (std::size_t i = 0; i != container.degree(); ++i)
      invoke(container.child(i), i);
}

/**
 * @brief Traverse each entry of a tree container
 * @details Functors may accept one or two values. The first one being the entry
 * to be evaluated, and the second one the multi-index of such entry
 *
 * @tparam Container          Recursive tree container to be traversed
 * @tparam PreCall            Functor to be applied at each node
 * @tparam LeafCall           Functor to be applied before each leaf node
 * @tparam PostCall           Functor to be applied after each node
 * @tparam MultiIndex         Multi-index representing the current position of the container
 */
template<Concept::TreeNode Node,
         class PreCall,
         class LeafCall,
         class PostCall,
         Concept::MultiIndex Prefix>
requires Concept::LeafTreeNode<Node>
         && (std::invocable<LeafCall&&, Node&&, const Prefix&> or std::invocable<LeafCall&&, Node&&>)
constexpr auto
forEachNode(Node&& node,
            PreCall&& pre_call,
            LeafCall&& leaf_call,
            PostCall&& post_call,
            Prefix multiindex)
{
  if constexpr (std::invocable<LeafCall&&, Node&&, const Prefix&>)
    leaf_call(std::forward<Node>(node), std::as_const(multiindex));
  else
    leaf_call(std::forward<Node>(node));
}

/**
 * @brief Traverse each entry of a tree container
 * @details Functors may accept one or two values. The first one being the entry
 * to be evaluated, and the second one the multi-index of such entry
 *
 * @tparam Container          Recursive tree container to be traversed
 * @tparam PreCall            Functor to be applied at each node
 * @tparam LeafCall           Functor to be applied before each leaf node
 * @tparam PostCall           Functor to be applied after each node
 * @tparam MultiIndex         Multi-index representing the current position of the container
 */
template<Concept::TreeNode Node,
         class PreCall,
         class LeafCall,
         class PostCall,
         Concept::MultiIndex Prefix>
requires (not Concept::LeafTreeNode<Node>)
         && (std::invocable<PreCall&&, Node&&, const Prefix&> or std::invocable<PreCall&&, Node&&>)
         && (std::invocable<PostCall&&, Node&&, const Prefix&> or std::invocable<PostCall&&, Node&&>)
constexpr auto
forEachNode(Node&& node,
            PreCall&& pre_call,
            LeafCall&& leaf_call,
            PostCall&& post_call,
            Prefix multiindex)
{
  auto invoke = [&multiindex,&node]<class Callable>(Callable&& callable){
    if constexpr (std::invocable<Callable&&, Node&&, const Prefix&>)
      callable(std::forward<Node>(node), std::as_const(multiindex));
    else
      callable(std::forward<Node>(node));
  };

  invoke(std::forward<PreCall>(pre_call));
  Dune::PDELab::forEachChild(
    std::forward<Node>(node),
    [&]<class Child>(Child&& child, auto i) {
      Dune::PDELab::forEachNode(
        std::forward<Child>(child),
        std::forward<PreCall>(pre_call),
        std::forward<LeafCall>(leaf_call),
        std::forward<PostCall>(post_call),
        push_back(multiindex, i)
      );
    });
  invoke(std::forward<PostCall>(post_call));
}

/**
 * @brief Traverse each entry of a tree container
 * @details Functors may accept one or two values. The first one being the entry
 * to be evaluated, and the second one the multi-index of such entry
 *
 * @tparam Container          Recursive tree container to be traversed
 * @tparam Callback           Functor to be applied at each node
 * @tparam MultiIndex         Multi-index representing the current position of the container
 */
template<Concept::TreeNode Node, class Callback>
constexpr auto
forEachNode(Node&& node, Callback&& callback)
{
  Dune::PDELab::forEachNode(
    std::forward<Node>(node),
    std::forward<Callback>(callback),
    std::forward<Callback>(callback),
    TypeTree::NoOp{},
    TypeTree::treePath()
  );
}


/**
 * @brief Traverse each leaf entry of a tree container
 * @details Functors may accept one or two values. The first one being the leaf
 * entry to be evaluated, and the second one the multi-index of such entry
 *
 * @tparam Container          Tree container to be traversed
 * @tparam Callback           Functor to be applied at each leaf node
 */
template<Concept::TreeNode Node, class Callback>
constexpr auto
forEachLeafNode(Node&& node, Callback&& callback)
{
  Dune::PDELab::forEachNode(
    std::forward<Node>(node),
    TypeTree::NoOp{},
    std::forward<Callback>(callback),
    TypeTree::NoOp{},
    TypeTree::treePath()
  );
}

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_COMMON_TREE_TRAVERSAL_HH
