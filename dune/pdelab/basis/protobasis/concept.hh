#ifndef DUNE_PDELAB_BASIS_PROTOBASIS_CONCEPT_HH
#define DUNE_PDELAB_BASIS_PROTOBASIS_CONCEPT_HH

#include <dune/pdelab/concepts/treenode.hh>
#include <dune/pdelab/common/tree_traversal.hh>

#include <concepts>

namespace Dune::PDELab::Concept {

namespace Impl {

template<class Node>
concept ProtoBasisNode = requires(Node node)
{
  { node.mergingStrategy() } -> std::convertible_to<typename std::remove_cvref_t<Node>::MergingStrategy>;
  { node.name() } -> std::convertible_to<std::string_view>;
  node.name(std::string{});
  node.name(std::string_view{});
  requires std::move_constructible<Node>;
};

void
requireProtoBasisNode(ProtoBasisNode auto&& node);

template<class Leaf>
concept ProtoBasisLeaf = requires(Leaf leaf)
{
  requires ProtoBasisNode<Leaf>;
  { leaf.finiteElementMap() } -> std::convertible_to<typename std::remove_cvref_t<Leaf>::FiniteElementMap>;
  { leaf.constraintsOperator() } -> std::convertible_to<typename std::remove_cvref_t<Leaf>::ConstraintsOperator>;
};

void
requireProtoBasisLeaf(ProtoBasisLeaf auto&& leaf);

template<class PBT>
concept ProtoBasisTree = requires(PBT pre_basis_tree)
{
  requires TreeNode<PBT>;
  Dune::PDELab::forEachNode(
    pre_basis_tree, [](auto&& node, auto tp) {
      Impl::requireProtoBasisNode(node);
      if constexpr (Concept::LeafTreeNode<decltype(node)>)
        Impl::requireProtoBasisLeaf(node);
    });
};

}

template<class Factory>
concept ProtoBasisFactory = requires(Factory factory)
{
  registerProtoBasisFactory(factory);
};

} // namespace Dune::PDELab::Concept::Impl

#endif // DUNE_PDELAB_BASIS_PROTOBASIS_CONCEPT_HH
