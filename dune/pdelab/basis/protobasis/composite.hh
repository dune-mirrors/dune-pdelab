#ifndef DUNE_PDELAB_BASIS_PROTOBASIS_COMPOSITE_HH
#define DUNE_PDELAB_BASIS_PROTOBASIS_COMPOSITE_HH

#include <dune/pdelab/basis/protobasis/concept.hh>
#include <dune/pdelab/basis/protobasis/node.hh>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/dynamicpowernode.hh>
#include <dune/typetree/powernode.hh>

#include <dune/common/indices.hh>

#include <array>
#include <memory>
#include <tuple>
#include <vector>

namespace Dune::PDELab {

/**
 * @brief Array node of Space
 *
 * @tparam MergingStrategy  Merging strategy
 * @tparam Node             A Space node type
 * @tparam degree           Size of the array
 */
template<class MergingStrategy,
         Concept::Impl::ProtoBasisNode Node,
         std::size_t degree>
class ProtoBasisArray
  : public TypeTree::PowerNode<Node, degree>
  , public ProtoBasisNode<ProtoBasisNodeTraits<MergingStrategy>>
{
private:
  using TreeNode = TypeTree::PowerNode<Node, degree>;
  using BaseNode = ProtoBasisNode<ProtoBasisNodeTraits<MergingStrategy>>;

public:
  ProtoBasisArray(
    const MergingStrategy& merging_strategy,
    const std::array<std::shared_ptr<Node>, degree>& nodes)
    : TreeNode{ nodes }
    , BaseNode{ merging_strategy }
  {
  }

  ProtoBasisArray(const ProtoBasisArray&) = default;
};

/**
 * @brief Make an array node of Space nodes
 *
 * @param merging_strategy  Merging strategy
 * @param nodes             Array of Space
 * @return auto             ProtoBasisArray
 */
template<class MergingStrategy, Concept::Impl::ProtoBasisNode Node, std::size_t degree>
auto
composite(const MergingStrategy& merging_strategy, const std::array<Node, degree>& nodes)
{
  std::array<std::shared_ptr<Node>, degree> storage;
  for (std::size_t i = 0; i < degree; ++i)
    storage[i] = std::make_shared<Node>(nodes[i]);
  using PB = ProtoBasisArray<MergingStrategy, Node, degree>;
  return PB{ merging_strategy, storage };
}

/**
 * @brief Vector node of Space
 *
 * @tparam MergingStrategy  Merging strategy
 * @tparam Node             A Space node type
 */
template<class MergingStrategy, Concept::Impl::ProtoBasisNode Node>
class ProtoBasisVector
  : public TypeTree::DynamicPowerNode<Node>
  , public ProtoBasisNode<ProtoBasisNodeTraits<MergingStrategy>>
{
private:
  using TreeNode = TypeTree::DynamicPowerNode<Node>;
  using BaseNode = ProtoBasisNode<ProtoBasisNodeTraits<MergingStrategy>>;

public:

  ProtoBasisVector(const MergingStrategy& merging_strategy, const std::vector<std::shared_ptr<Node>>& nodes)
    : TreeNode{ nodes }
    , BaseNode{ merging_strategy }
  {
  }

  ProtoBasisVector(const ProtoBasisVector&) = default;
};

/**
 * @brief Make an vector node of Space nodes
 *
 * @param merging_strategy  Merging strategy
 * @param nodes             Vector of Space
 * @return auto             ProtoBasisVector
 */
template<class MergingStrategy, Concept::Impl::ProtoBasisNode Node>
auto
composite(const MergingStrategy& merging_strategy, const std::vector<Node>& nodes)
{
  std::vector<std::shared_ptr<Node>> storage(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i)
    storage[i] = std::make_shared<Node>(nodes[i]);
  using PB = ProtoBasisVector<MergingStrategy, Node>;
  return PB{ merging_strategy, storage };
}

/**
 * @brief Tuple node of Space
 *
 * @tparam MergingStrategy  Merging strategy
 * @tparam Nodes            Space node types
 */
template<class MergingStrategy, Concept::Impl::ProtoBasisNode... Nodes>
class ProtoBasisTuple
  : public TypeTree::CompositeNode<Nodes...>
  , public ProtoBasisNode<ProtoBasisNodeTraits<MergingStrategy>>
{
private:
  using TreeNode = TypeTree::CompositeNode<Nodes...>;
  using BaseNode = ProtoBasisNode<ProtoBasisNodeTraits<MergingStrategy>>;

public:
  ProtoBasisTuple(const MergingStrategy& merging_strategy, const std::tuple<std::shared_ptr<Nodes>...>& nodes)
    : TreeNode{ nodes }
    , BaseNode{ merging_strategy }
  {
  }

  ProtoBasisTuple(const ProtoBasisTuple&) = default;
};

/**
 * @brief Make a tuple node of Space nodes
 *
 * @param merging_strategy  Merging strategy
 * @param nodes             Tuple of Space
 * @return auto             ProtoBasisTuple
 */
template<class MergingStrategy, Concept::Impl::ProtoBasisNode... Nodes>
auto
composite(const MergingStrategy& merging_strategy, const std::tuple<Nodes...>& nodes)
{
  using TypeTuple = std::tuple<Nodes...>;
  using Storage = std::tuple<std::shared_ptr<Nodes>...>;
  auto sequence = std::make_index_sequence<sizeof...(Nodes)>{};
  auto storage = unpackIntegerSequence(
    [&](auto... i) {
      return Storage{ std::make_shared<std::tuple_element_t<i, TypeTuple>>(
        std::get<i>(nodes))... };
    },
    sequence);
  using PB = ProtoBasisTuple<MergingStrategy, Nodes...>;
  return PB{ merging_strategy, storage };
}

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_BASIS_PROTOBASIS_COMPOSITE_HH
