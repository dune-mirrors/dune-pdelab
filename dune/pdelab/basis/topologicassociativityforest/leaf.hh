#ifndef DUNE_PDELAB_BASIS_TOPOLOGIC_ASSOCIATIVITY_FOREST_LEAF_HH
#define DUNE_PDELAB_BASIS_TOPOLOGIC_ASSOCIATIVITY_FOREST_LEAF_HH

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
 * @class LeafTopologicAssociativityForest
 * @brief Leaf implementation of a TopologicAssociativityForestNode.
 *
 * This class represents the leaf node in a topological associativity forest, specifically designed to handle
 * the degrees of freedom associated with finite element spaces. It is the first ordering constructed on every
 * ordering tree and provides the foundation for managing the leaf-level data in the forest structure.
 *
 * @tparam ProtoBasisLeaf The type containing a finite element map and a merging strategy. This type defines
 *                        the characteristics of the finite elements and how they are merged within the forest.
 */
template <Concept::Impl::ProtoBasisLeaf ProtoBasisLeaf, class GridView>
class LeafTopologicAssociativityForest
    : public TypeTree::LeafNode,
      public TopologicAssociativityForestNode<LeafTopologicAssociativityForest<ProtoBasisLeaf,GridView>,
                                              GridView, typename ProtoBasisLeaf::MergingStrategy>
{
  using TreeNode = TypeTree::LeafNode;
  using OrderingNode =
      TopologicAssociativityForestNode<LeafTopologicAssociativityForest<ProtoBasisLeaf, GridView>,
                                       GridView, typename ProtoBasisLeaf::MergingStrategy>;

public:
  using ProtoBasis = ProtoBasisLeaf;
  using SizeType = typename OrderingNode::SizeType;

  //! @brief Constructs a LeafTopologicAssociativityForest with a given proto-basis.
  LeafTopologicAssociativityForest(const ProtoBasis &proto_basis, const GridView& grid_view)
      : TreeNode{}, OrderingNode{grid_view, proto_basis.mergingStrategy()}, _proto_basis{proto_basis}
  {
  }

  LeafTopologicAssociativityForest(const LeafTopologicAssociativityForest &) = delete;
  LeafTopologicAssociativityForest(LeafTopologicAssociativityForest &&) = default;

  LeafTopologicAssociativityForest &operator=(const LeafTopologicAssociativityForest &) = delete;
  LeafTopologicAssociativityForest &operator=(LeafTopologicAssociativityForest &&) = default;

  //! @brief Accesses the proto-basis associated with this leaf node.
  const ProtoBasis &protoBasis() const { return _proto_basis; }

  //! @brief Gets the maximum number of sub-entities.
  [[nodiscard]] int maxSubEntities() const { return _max_sub_entities; }

  //! @brief Sets the maximum number of sub-entities.
  void setMaxSubEntities(int s) { _max_sub_entities = s; }

  /**
   * @brief Gets a common size for all active geometry types if available at compile time.
   *
   * @return An optional size if a common size is available at compile time; otherwise, std::nullopt.
   *
   * @note This requires that ProtoBasis::FiniteElementMap::size(...) is static constexpr.
   */
  static constexpr std::optional<std::size_t> commonSizePerGeometryType();

private:
  ProtoBasis _proto_basis;
  int _max_sub_entities = 0;
};

template<Concept::Impl::ProtoBasisLeaf ProtoBasisLeaf, class GridView>
constexpr std::optional<std::size_t> LeafTopologicAssociativityForest<ProtoBasisLeaf,GridView>::commonSizePerGeometryType() {
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

} // namespace Dune::PDELab::Impl

#endif // DUNE_PDELAB_BASIS_TOPOLOGIC_ASSOCIATIVITY_FOREST_LEAF_HH
