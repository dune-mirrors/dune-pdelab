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
 * @brief Traverse each child of a tree and apply a callable function.
 *
 * This function iterates over each child node of a given tree container and
 * applies a callable function to each child. The callable function can accept
 * either one or two arguments: the child node itself, and optionally, its index.
 *
 * @tparam Container The type of the tree container to be traversed. This must
 *                   satisfy the Concept::ParentTreeNode concept.
 * @tparam Callable  The type of the callable function (functor) to be applied
 *                   to each child node. This can be a lambda, function pointer,
 *                   or any callable object.
 *
 * @param container The tree container whose children will be traversed.
 * @param at_value  The callable function to be applied to each child node.
 *                  This function can accept either one argument (the child node)
 *                  or two arguments (the child node and its index).
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
      [&](auto... i) { (invoke(std::forward<Container>(container).child(i), i), ...); },
      std::make_index_sequence<std::remove_cvref_t<Container>::degree()>{});
  else
    for (std::size_t i = 0; i != container.degree(); ++i)
      invoke(std::forward<Container>(container).child(i), i);
}

/**
 * @brief Traverse each entry of a tree and apply functors at different stages.
 *
 * This function traverses each entry of a recursive tree container and applies
 * specified functors at different stages of the traversal: before each node,
 * before each leaf node, and after each node. The functors can accept either
 * one or two arguments: the entry (node) and optionally its multi-index.
 *
 * @tparam Node       The type of the tree node to be traversed. This must
 *                    satisfy the Concept::TreeNode and Concept::LeafTreeNode
 *                    concepts.
 * @tparam PreCall    The type of the functor to be applied before each node.
 * @tparam LeafCall   The type of the functor to be applied before each leaf node.
 * @tparam PostCall   The type of the functor to be applied after each node.
 * @tparam Prefix     The type representing the multi-index of the current
 *                    position in the tree container. This must satisfy the
 *                    Concept::MultiIndex concept.
 *
 * @param node        The tree node to be traversed.
 * @param pre_call    The functor to be applied before each node.
 * @param leaf_call   The functor to be applied before each leaf node.
 * @param post_call   The functor to be applied after each node.
 * @param multiindex  The multi-index representing the current position in the
 *                    tree container.
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
 * @brief Traverse each entry of a tree and apply functors at different stages.
 *
 * This function traverses each entry of a recursive tree container and applies
 * specified functors at different stages of the traversal: before each node,
 * before each leaf node, and after each node. The functors can accept either
 * one or two arguments: the entry (node) and optionally its multi-index.
 *
 * @tparam Node       The type of the tree node to be traversed. This must
 *                    satisfy the Concept::TreeNode and Concept::LeafTreeNode
 *                    concepts.
 * @tparam PreCall    The type of the functor to be applied before each node.
 * @tparam LeafCall   The type of the functor to be applied before each leaf node.
 * @tparam PostCall   The type of the functor to be applied after each node.
 * @tparam Prefix     The type representing the multi-index of the current
 *                    position in the tree container. This must satisfy the
 *                    Concept::MultiIndex concept.
 *
 * @param node        The tree node to be traversed.
 * @param pre_call    The functor to be applied before each node.
 * @param leaf_call   The functor to be applied before each leaf node.
 * @param post_call   The functor to be applied after each node.
 * @param multiindex  The multi-index representing the current position in the
 *                    tree container.
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
 * @brief Traverse each entry of a tree container and apply a callback function.
 *
 * This function traverses each entry of a recursive tree container and applies
 * a specified callback function to each node. The callback function can accept
 * either one or two arguments: the entry (node) and optionally its multi-index.
 *
 * @tparam Node       The type of the tree node to be traversed. This must
 *                    satisfy the Concept::TreeNode concept.
 * @tparam Callback   The type of the callback function to be applied at each node.
 *                    This can be a lambda, function pointer, or any callable
 *                    object.
 *
 * @param node        The tree node to be traversed.
 * @param callback    The callback function to be applied at each node.
 *                    This function can accept either one argument (the node)
 *                    or two arguments (the node and its multi-index).
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
 * @brief Traverse each leaf entry of a tree container and apply a callback function.
 *
 * This function traverses each leaf entry of a tree container and applies a specified
 * callback function to each leaf node. The callback function can accept either one or
 * two arguments: the leaf entry (node) and optionally its multi-index.
 *
 * @tparam Node       The type of the tree node to be traversed. This must
 *                    satisfy the Concept::TreeNode concept.
 * @tparam Callback   The type of the callback function to be applied at each leaf node.
 *                    This can be a lambda, function pointer, or any callable object.
 *
 * @param node        The tree node to be traversed.
 * @param callback    The callback function to be applied at each leaf node.
 *                    This function can accept either one argument (the leaf node)
 *                    or two arguments (the leaf node and its multi-index).
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
