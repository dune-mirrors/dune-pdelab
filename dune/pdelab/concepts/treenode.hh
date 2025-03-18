#ifndef DUNE_PDELAB_CONCEPT_TREE_NODE_HH
#define DUNE_PDELAB_CONCEPT_TREE_NODE_HH

#include <dune/common/indices.hh>

#include <concepts>
#include <type_traits>


namespace Dune::PDELab::Concept {

  namespace Impl{

    template<class Node>
    concept DynamicChildAccess = requires(Node node, std::size_t i)
    {
      node.child(i);
    };

    template<class Node, std::size_t I>
    concept StaticChildAccess = requires(Node node, std::integral_constant<std::size_t, I> _i)
    {
      node.child(_i);
    };

    template<class T>
    using asInt = int;

    // this function instantiates and counts all the children (if put directly in gcc-10, it gives a fatal error)
    template<class Node>
    static constexpr std::size_t countChildren() {
      return Dune::unpackIntegerSequence([](auto... i){
          return (asInt<decltype(std::declval<Node>().child(i))>{1} + ...);
        }, std::make_index_sequence<std::remove_cvref_t<Node>::degree()>{});
    }

    template<class I>
    concept IntegralConstant = requires(I i) {
      { I::value } -> std::convertible_to<typename I::value_type>;
      { std::integral_constant<typename I::value_type, I::value>{} }-> std::convertible_to<typename I::value_type>;
    };

    template<class Node>
    concept StaticTreeAccess = requires(Node node)
    {
      { std::remove_cvref_t<Node>::degree() } -> IntegralConstant;
      requires (std::remove_cvref_t<Node>::degree() != 0);
      requires StaticChildAccess<Node, 0>;
      requires StaticChildAccess<Node, (std::remove_cvref_t<Node>::degree()-1)>;
      requires std::remove_cvref_t<Node>::degree() == countChildren<Node>();
    };

  }

  //!@brief Model of a leaf node of a typetree
  template<class Node>
  concept LeafTreeNode = requires
  {
    { std::remove_cvref_t<Node>::degree()} -> Impl::IntegralConstant;
    requires (std::remove_cvref_t<Node>::degree() == 0);
  };

  //!@brief Model of a tuple node of a typetree
  template<class Node>
  concept TupleTreeNode = (not Impl::DynamicChildAccess<Node>) && Impl::StaticTreeAccess<Node>;

  //!@brief Model of a tuple node of a typetree
  template<class Node>
  concept ArrayTreeNode = Impl::DynamicChildAccess<Node> && Impl::StaticTreeAccess<Node>;

  //!@brief Model of a tuple node of a typetree
  template<class Node>
  concept VectorTreeNode = Impl::DynamicChildAccess<Node> && (not Impl::StaticTreeAccess<Node>);

  //!@brief Model of a parent node of a typetree
  template<class Node>
  concept ParentTreeNode = requires
  {
    requires VectorTreeNode<Node> || ArrayTreeNode<Node> || TupleTreeNode<Node>;
  };

  //!@brief Model of a node of a typetree
  template<class Node>
  concept TreeNode = requires
  {
    requires ParentTreeNode<Node> || LeafTreeNode<Node>;
  };

} // namespace Dune::PDELab::Concept

#endif // DUNE_PDELAB_CONCEPT_TREE_NODE_HH
