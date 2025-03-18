#ifndef DUNE_PDELAB_BASIS_ORDERING_TOPOLOGIC_ASSOCIATIVITY_FOREST_LEAF_HH
#define DUNE_PDELAB_BASIS_ORDERING_TOPOLOGIC_ASSOCIATIVITY_FOREST_LEAF_HH

#include <dune/pdelab/basis/protobasis/concept.hh>
#include <dune/pdelab/basis/topologicassociativityforest/node.hh>

#include <dune/pdelab/common/multiindex.hh>

#include <dune/typetree/leafnode.hh>
#include <dune/typetree/dynamicpowernode.hh>
#include <dune/typetree/treepath.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>
#include <dune/geometry/referenceelements.hh>

#include <vector>
#include <optional>
#include <numeric>

namespace Dune::PDELab::Impl {

/**
 * @brief Leaf implementation of an TopologicAssociativityForestNode
 * @note This is the first ordering constructed on every ordering tree
 *
 * @tparam ProtoBasisLeaf   Type containing a finite element map and a
 *                        merging strategy
 */
template<Concept::Impl::ProtoBasisLeaf ProtoBasisLeaf>
class LeafTopologicAssociativityForest
  : public TypeTree::LeafNode
  , public TopologicAssociativityForestNode<LeafTopologicAssociativityForest<ProtoBasisLeaf>,
                              typename ProtoBasisLeaf::MergingStrategy>
{
  using TreeNode = TypeTree::LeafNode;
  using OrderingNode =
    TopologicAssociativityForestNode<LeafTopologicAssociativityForest<ProtoBasisLeaf>,
                       typename ProtoBasisLeaf::MergingStrategy>;
  using FE = typename ProtoBasisLeaf::FiniteElementMap::Traits::FiniteElement;
  using ES = typename ProtoBasisLeaf::MergingStrategy::EntitySet;
  static constexpr std::size_t fem_dim = FE::Traits::LocalBasisType::Traits::dimDomain;
  static constexpr std::size_t fem_codim = ES::dimension - fem_dim;
public:
  using ProtoBasis = ProtoBasisLeaf;
  using SizeType = typename OrderingNode::SizeType;

  LeafTopologicAssociativityForest(const ProtoBasis& proto_basis)
    : TreeNode{}
    , OrderingNode{ proto_basis.mergingStrategy() }
    , _proto_basis{ proto_basis }
  {
  }

  LeafTopologicAssociativityForest(const LeafTopologicAssociativityForest&) = delete;
  LeafTopologicAssociativityForest(LeafTopologicAssociativityForest&&) = default;

  LeafTopologicAssociativityForest& operator=(const LeafTopologicAssociativityForest&) = delete;
  LeafTopologicAssociativityForest& operator=(LeafTopologicAssociativityForest&&) = default;

  const ProtoBasis& protoBasis() const { return _proto_basis; }

  ProtoBasis& protoBasis() { return _proto_basis; }

  [[nodiscard]] int maxSubEntities() const { return _max_sub_entities; }
  void setMaxSubEntities(int s) { _max_sub_entities = s; }

  //! gets a common size for all active geometry types if available at compile time
  // this needs the the ProtoBasis::FiniteElementMap::size(...) is static constexpr.
  static constexpr std::optional<std::size_t> commonSizePerGeometryType() {
    using FEM = typename ProtoBasis::FiniteElementMap;
    std::size_t size = 0;
    const auto fem_size = [](auto gt) constexpr {
      if constexpr (requires { { FEM::size(gt) } -> std::convertible_to<std::size_t>; })
        return FEM::size(gt);
      else
        return 0;
    };
    // iterate over all possible geometry types and find out if all share the same size
    for (std::size_t dim = 0; dim <= FEM::dimension ; ++dim) {
      std::size_t gt_size = fem_size(GeometryTypes::none(dim));
      if (gt_size > 0) {
        if (size > 0 and size != gt_size)
          return std::nullopt;
        else
          size = gt_size;
      }
      for (std::size_t topology_id = 0 ; topology_id < (std::size_t{1} << dim) ; ++topology_id) {
        std::size_t gt_size = fem_size(GeometryType(topology_id,dim));
        if (gt_size > 0) {
          if (size > 0 and size != gt_size)
            return std::nullopt;
          else
            size = gt_size;
        }
      }
    }
    if (size == 0)
      return std::nullopt;
    else
      return size;
  }

private:
  ProtoBasis _proto_basis;
  int _max_sub_entities = 0;
};

} // namespace Dune::PDELab::Impl

#endif // DUNE_PDELAB_BASIS_ORDERING_TOPOLOGIC_ASSOCIATIVITY_FOREST_LEAF_HH
