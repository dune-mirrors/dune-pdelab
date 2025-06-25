
#ifndef DUNE_PDELAB_BASIS_LAGRANGE_FACTORY_HH
#define DUNE_PDELAB_BASIS_LAGRANGE_FACTORY_HH

#include <dune/pdelab/basis/protobasis.hh>
#include <dune/pdelab/basis/mergingstrategy.hh>
#include <dune/pdelab/basis/factory.hh>

#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>

#include <dune/common/indices.hh>

namespace Dune::PDELab::BasisFactory {

/**
 * @brief Factory for Qk Lagrange finite element spaces (hypercubes).
 *
 * @tparam k Polynomial order (default 1)
 * @tparam MergingStrategy Merging strategy type (default FlatTopologicalInterleaving)
 * @param strategy Merging strategy (e.g., flatByEntity(), blockedByEntity())
 * @return DefaultPreBasisFactory for Qk Lagrange basis
 */
template<std::size_t k = 1, class MergingStrategy = FlatTopologicalInterleaving>
requires std::same_as<MergingStrategy, FlatTopologicalInterleaving> ||
         std::same_as<MergingStrategy, BlockedTopologicalInterleaving>
auto makeQkProtoBasis(MergingStrategy strategy = {}) {
  return DefaultPreBasisFactory{[strategy]<class GridView>(GridView grid_view) {
    using QkFEM = Dune::PDELab::QkLocalFiniteElementMap<GridView, double, double, k>;
    return makeProtoBasis(QkFEM(grid_view), strategy).protoBasis(grid_view);
  }};
}

/**
 * @brief Factory for Pk Lagrange finite element spaces (simplex).
 *
 * @tparam k Polynomial order (default 1)
 * @tparam MergingStrategy Merging strategy type (default FlatTopologicalInterleaving)
 * @param strategy Merging strategy (e.g., flatByEntity(), blockedByEntity())
 * @return DefaultPreBasisFactory for Pk Lagrange basis
 */
template<std::size_t k = 1, class MergingStrategy = FlatTopologicalInterleaving>
requires std::same_as<MergingStrategy, FlatTopologicalInterleaving> ||
         std::same_as<MergingStrategy, BlockedTopologicalInterleaving>
auto makePkProtoBasis(MergingStrategy strategy = {}) {
  return DefaultPreBasisFactory{[strategy]<class GridView>(GridView grid_view) {
    using PkFEM = Dune::PDELab::PkLocalFiniteElementMap<GridView, double, double, k>;
    return makeProtoBasis(PkFEM(grid_view), strategy).protoBasis(grid_view);
  }};
}

namespace Experimental {

//! Descriptor for Lagrange finite element family
struct LagrangeFiniteElementFamily : FiniteElementFamily {};

//! Global instance of Lagrange finite element family
inline constexpr LagrangeFiniteElementFamily CG;

/**
 * @brief Factory function to create a Lagrange finite element basis for given geometry and degree.
 *
 * @tparam G Geometry type (GeometryType::Id)
 * @tparam Degree Polynomial degree (default 1)
 * @tparam MergingStrategy Merging strategy type (default FlatTopologicalInterleaving)
 * @param family Lagrange finite element family descriptor
 * @param geometry Geometry type
 * @param degree Polynomial degree (default 1)
 * @param strategy Merging strategy (e.g., flatByEntity(), blockedByEntity())
 * @return A ProtoBasis factory for the specified finite element.
 */
template <Dune::GeometryType::Id G, std::size_t Degree = 1, class MergingStrategy = FlatTopologicalInterleaving>
requires std::same_as<MergingStrategy, FlatTopologicalInterleaving> ||
         std::same_as<MergingStrategy, BlockedTopologicalInterleaving>
auto FiniteElement(LagrangeFiniteElementFamily, std::integral_constant<GeometryType::Id, G>, index_constant<Degree> = {}, MergingStrategy strategy = {}) {
  return DefaultPreBasisFactory{[strategy]<class GridView>(GridView grid_view) {
    static_assert(Dune::GeometryType(G).dim() == GridView::dimension, "Geometry dimension and grid dimension do not match");
    if constexpr (Dune::GeometryType(G).isSimplex()) {
      using PkFEM = Dune::PDELab::PkLocalFiniteElementMap<GridView, double, double, Degree>;
      return makeProtoBasis(PkFEM(grid_view), strategy).protoBasis(grid_view);
    } else if constexpr (Dune::GeometryType(G).isCube()) {
      using QkFEM = Dune::PDELab::QkLocalFiniteElementMap<GridView, double, double, Degree>;
      return makeProtoBasis(QkFEM(grid_view), strategy).protoBasis(grid_view);
    } else {
      static_assert(AlwaysFalse<GridView>{}, "Lagrange finite elements are only implemented for simplex and cube geometries");
    }
  }};
}

} // namespace Experimental
} // namespace Dune::PDELab::BasisFactory

#endif // DUNE_PDELAB_BASIS_LAGRANGE_FACTORY_HH
