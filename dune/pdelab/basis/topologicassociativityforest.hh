#ifndef DUNE_PDELAB_BASIS_ORDERING_TOPOLOGIC_ASSOCIATIVITY_FOREST_HH
#define DUNE_PDELAB_BASIS_ORDERING_TOPOLOGIC_ASSOCIATIVITY_FOREST_HH

#include <dune/pdelab/basis/topologicassociativityforest/node.hh>
#include <dune/pdelab/basis/topologicassociativityforest/leaf.hh>
#include <dune/pdelab/basis/topologicassociativityforest/composite.hh>

namespace Dune::PDELab::Impl
{

//! Returns a suitable Entity Ordering tree for a discrete-function-space tree
template<Concept::Impl::ProtoBasisTree ProtoBasis, class GridView>
auto makeTopologicAssociativityForest(const ProtoBasis& proto_basis, const GridView& grid_view)
{
  using MergingStrategy = typename ProtoBasis::MergingStrategy;

  if constexpr (Concept::LeafTreeNode<ProtoBasis>) {
    return LeafTopologicAssociativityForest<ProtoBasis, GridView>(proto_basis, grid_view);
  } else {

    auto makeChildNode = [&]<class Index>(Index i){
      return makeTopologicAssociativityForest(proto_basis.child(i), grid_view);
    };

    if constexpr (Concept::ArrayTreeNode<ProtoBasis>) {
      constexpr std::size_t degree = ProtoBasis::degree();
      using ChildNode = decltype(makeChildNode(0));
      using Node = ArrayTopologicAssociativityForest<MergingStrategy, ChildNode, degree>;
      typename Node::NodeStorage storage;
      for (std::size_t i = 0; i < degree; ++i)
        storage[i] = std::make_unique<ChildNode>(makeChildNode(i));
      return Node(std::move(storage), proto_basis.mergingStrategy());
    } else if constexpr (Concept::VectorTreeNode<ProtoBasis>) {
      std::size_t degree = proto_basis.degree();
      using ChildNode = decltype(makeChildNode(0));
      using Node = VectorTopologicAssociativityForest<MergingStrategy, ChildNode>;
      typename Node::NodeStorage storage(degree);
      for (std::size_t i = 0; i < degree; ++i)
        storage[i] = std::make_unique<Node>(makeChildNode(i));
      return Node(std::move(storage), proto_basis.mergingStrategy());
    } else {
      static_assert(Concept::TupleTreeNode<ProtoBasis>);

      auto unfold_children = [&](auto... i) {
        using Node = TupleTopologicAssociativityForest<MergingStrategy, decltype(makeChildNode(i))...>;
        typename Node::NodeStorage storage{ std::make_unique<decltype(makeChildNode(i))>(makeChildNode(i))... };
        return Node(std::move(storage), proto_basis.mergingStrategy());
      };
      auto range = std::make_index_sequence<ProtoBasis::degree()>{};
      return unpackIntegerSequence(unfold_children, range);
    }
  }
}

template<Concept::Impl::ProtoBasisTree ProtoBasis, class GridView>
using TopologicAssociativityForest = decltype(makeTopologicAssociativityForest(std::declval<ProtoBasis>(), std::declval<GridView>()));

} // namespace Dune::PDELab::Impl

#endif // DUNE_PDELAB_BASIS_ORDERING_TOPOLOGIC_ASSOCIATIVITY_FOREST_HH
