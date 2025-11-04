#ifndef DUNE_PDELAB_BASIS_FACTORY_HH
#define DUNE_PDELAB_BASIS_FACTORY_HH

#include <dune/common/indices.hh>

#include <dune/pdelab/basis/prebasis.hh>
#include <dune/pdelab/basis/protobasis.hh>
#include <dune/pdelab/basis/mergingstrategy.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>

/**
 * @file
 * @brief Entry point for constructing PDELab bases using factory patterns.
 *
 * This header provides a flexible and composable system for constructing
 * finite element bases in DUNE PDELab, compatible with DUNE Functions.
 * It enables users to build local and global bases using proto-basis and
 * pre-basis factories, and supports composition and power constructions.
 *
 * The main difference to dune-functions is that the factories create bases
 * that can have a topological interleaving of the local trees, i.e.
 * the degrees of freedom are grouped by topological entity.
 *
 * ## Overview
 *
 * The main concept is to use factory objects to build basis trees for a given grid view.
 * These factories can be composed, powered, and customized with merging strategies,
 * allowing users to easily construct complex function spaces.
 *
 * This file is typically used in conjunction with:
 * - DUNE PDELab's proto-basis and pre-basis system
 * - DUNE Functions' global basis system (see `defaultglobalbasis.hh`)
 * - Power and composite basis constructions (see `powerbasis.hh`, `compositebasis.hh`)
 *
 * ## Integration with DUNE Functions
 *
 * The factories here are designed to work seamlessly with DUNE Functions' basis system.
 * For example, you can use these factories with `makeBasis` from DUNE Functions,
 * together with other pre-basis provided by DUNE Functions.
 *
 */

namespace Dune::PDELab::BasisFactory {

using namespace Dune::Indices;

template<class Base>
struct PreBasisFactory;

template<class Callable>
struct DefaultPreBasisFactory;

/**
 * @brief Create a ProtoBasis factory for a given finite element map and merging strategy.
 *
 * It returns a factory that can be used to build PreBasis objects for any grid view.
 *
 * @tparam FiniteElementMap Type of the finite element map.
 * @param finite_element_map The finite element map to use.
 * @return A DefaultPreBasisFactory object.
 *
 */
template<class FiniteElementMap>
auto makeProtoBasis(const FiniteElementMap& finite_element_map, FlatTopologicalInterleaving = {})  {
  return DefaultPreBasisFactory{[finite_element_map]<class GridView, class SubDomain>(const GridView& grid_view, const SubDomain& sub_domain) {
    using ProtoBasis = ProtoBasisLeaf<GridView, FiniteElementMap, false, SubDomain>;
    return ProtoBasis{grid_view, finite_element_map, sub_domain};
  }};
}

/**
 * @brief Create a ProtoBasis factory for a given finite element map and merging strategy.
 *
 * It returns a factory that can be used to build PreBasis objects for any grid view.
 *
 * @tparam FiniteElementMap Type of the finite element map.
 * @param finite_element_map The finite element map to use.
 * @return A DefaultPreBasisFactory object.
 *
 */
template<class FiniteElementMap>
auto makeProtoBasis(const FiniteElementMap& finite_element_map, BlockedTopologicalInterleaving)  {
  return DefaultPreBasisFactory{[finite_element_map]<class GridView, class SubDomain>(const GridView& grid_view, SubDomain sub_domain) {
    using ProtoBasis = ProtoBasisLeaf<GridView, FiniteElementMap, true, SubDomain>;
    return ProtoBasis{grid_view, finite_element_map, sub_domain};
  }};
}

/**
 * @brief Create a power basis factory with compile-time exponent.
 *
 * Constructs a factory for a power basis (vector-valued basis) with `k` copies
 * of the child basis.
 *
 * @tparam k Number of copies (exponent).
 * @tparam Factory Factory for the child proto-basis.
 * @param factory Factory for the child proto-basis.
 * @return A DefaultPreBasisFactory for the power basis.
 *
 */
template<std::size_t k, class Factory>
auto power(Factory&& factory, FlatTopologicalInterleaving) {
  return DefaultPreBasisFactory{[factory = std::forward<Factory>(factory)](const auto& grid_view, const auto& sub_domain){
    auto proto_basis = factory.protoBasis(grid_view, sub_domain);
    using ProtoBasis = decltype(proto_basis);
    std::array<ProtoBasis, k> nodes = unpackIntegerSequence([&](auto... i) {
      return std::array<ProtoBasis, k>{ ((void)i, proto_basis)... };
    }, std::make_index_sequence<k>{});

    return ProtoBasisArray<ProtoBasis, k, false>{ nodes };
  }};
}

/**
 * @brief Create a power basis factory with compile-time exponent.
 *
 * Constructs a factory for a power basis (vector-valued basis) with `k` copies
 * of the child basis.
 *
 * @tparam k Number of copies (exponent).
 * @tparam Factory Factory for the child proto-basis.
 * @param factory Factory for the child proto-basis.
 * @return A DefaultPreBasisFactory for the power basis.
 *
 */
template<std::size_t k, class Factory>
auto power(Factory&& factory, BlockedTopologicalInterleaving) {
  return DefaultPreBasisFactory{[factory = std::forward<Factory>(factory)](const auto& grid_view, const auto& sub_domain){
    auto proto_basis = factory.protoBasis(grid_view, sub_domain);
    using ProtoBasis = decltype(proto_basis);
    std::array<ProtoBasis, k> nodes = unpackIntegerSequence([&](auto... i) {
      return std::array<ProtoBasis, k>{ ((void)i, proto_basis)... };
    }, std::make_index_sequence<k>{});

    return ProtoBasisArray<ProtoBasis, k, true>{ nodes };
  }};
}

/**
 * @brief Create a power basis factory with run-time exponent.
 *
 * Constructs a factory for a power basis (vector-valued basis) with `k` copies,
 * where `k` is determined at run time.
 *
 * @tparam Factory Factory for the child proto-basis.
 * @param factory Factory for the child proto-basis.
 * @param k Number of copies (run-time).
 * @return A DefaultPreBasisFactory for the power basis.
 *
 */
template<class Factory>
auto power(Factory&& factory, std::size_t k, FlatTopologicalInterleaving = {}) {
  return DefaultPreBasisFactory{[factory = std::forward<Factory>(factory), k](const auto& grid_view, const auto& sub_domain){
    auto proto_basis = factory.protoBasis(grid_view, sub_domain);
    using ProtoBasis = decltype(proto_basis);
    std::vector<ProtoBasis> nodes;
    for (std::size_t i = 0; i < k; ++i)
      nodes.emplace_back(proto_basis);
    return ProtoBasisVector<ProtoBasis, false>{ nodes };
  }};
}

/**
 * @brief Create a power basis factory with run-time exponent.
 *
 * Constructs a factory for a power basis (vector-valued basis) with `k` copies,
 * where `k` is determined at run time.
 *
 * @tparam Factory Factory for the child proto-basis.
 * @param factory Factory for the child proto-basis.
 * @param k Number of copies (run-time).
 * @return A DefaultPreBasisFactory for the power basis.
 *
 */
template<class Factory>
auto power(Factory&& factory, std::size_t k, BlockedTopologicalInterleaving){
  return DefaultPreBasisFactory{[factory = std::forward<Factory>(factory), k](const auto& grid_view, const auto& sub_domain){
    auto proto_basis = factory.protoBasis(grid_view, sub_domain);
    using ProtoBasis = decltype(proto_basis);
    std::vector<ProtoBasis> nodes;
    for (std::size_t i = 0; i < k; ++i)
      nodes.emplace_back(proto_basis);
    return ProtoBasisVector<ProtoBasis, true>{ nodes };
  }};
}


/**
 * @brief Create a composite basis factory from several child factories and a merging strategy.
 *
 * This allows you to combine several bases into a single composite basis, e.g. for
 * mixed finite element methods or multi-physics problems.
 *
 * \internal Note that the original composite function in dune-functions is too eager to sink all template arguments.
 * Thus, to provide a proper specialization, we use the first factory as a separate argument.
 * If you come up with a better solution, please refactor this code.
 *
 * @tparam Callable0 Callable of the child proto-basis factory.
 * @tparam Args Additional child factories and the merging strategy.
 * @param factory0 First child factory.
 * @param args Additional child factories and the merging strategy.
 * @return A DefaultPreBasisFactory for the composite basis.
 *
 */
template<class Callable0, class... Args>
  requires (std::convertible_to<std::tuple_element_t<sizeof...(Args)-1, std::tuple<Args...>>, FlatTopologicalInterleaving>)
        || (std::convertible_to<std::tuple_element_t<sizeof...(Args)-1, std::tuple<Args...>>, BlockedTopologicalInterleaving>)
auto composite(DefaultPreBasisFactory<Callable0>&& factory0, Args&&... args)
{
  using ArgsTuple = std::tuple<DefaultPreBasisFactory<Callable0>, Args...>;
  auto argsTuple = std::make_tuple(std::move(factory0), std::forward<Args>(args)...);
  using MergingStrategy = std::tuple_element_t<sizeof...(Args), ArgsTuple>;
  constexpr bool isIndexBlocked = std::is_same_v<MergingStrategy, BlockedTopologicalInterleaving>;

  return DefaultPreBasisFactory{[=](const auto& grid_view, const auto& sub_domain){
    auto sequence = std::make_index_sequence<sizeof...(Args)>{};
    return unpackIntegerSequence(
      [&](auto... i) {
        std::tuple storage = { std::get<i>(argsTuple).protoBasis(grid_view, sub_domain)... };
        return ProtoBasisTuple<isIndexBlocked, decltype(std::get<i>(argsTuple).protoBasis(grid_view, sub_domain))...>{storage};
      },
      sequence);
  }};
}

/**
 * @brief Create a composite basis factory from several child factories and a merging strategy.
 *
 * This allows you to combine several bases into a single composite basis, e.g. for
 * mixed finite element methods or multi-physics problems.
 *
 * \internal Note that the original composite function in dune-functions is too eager to sink all template arguments.
 * Thus, to provide a proper specialization, we use the first factory as a separate argument.
 * If you come up with a better solution, please refactor this code.
 *
 * @tparam Callable0 Callable of the child proto-basis factory.
 * @tparam Args Additional child factories and the merging strategy.
 * @param factory0 First child factory.
 * @param args Additional child factories and the merging strategy.
 * @return A DefaultPreBasisFactory for the composite basis.
 *
 */
template<class Callable0, class... Args>
  requires (std::convertible_to<std::tuple_element_t<sizeof...(Args)-1, std::tuple<Args...>>, FlatTopologicalInterleaving>)
        || (std::convertible_to<std::tuple_element_t<sizeof...(Args)-1, std::tuple<Args...>>, BlockedTopologicalInterleaving>)
auto composite(const DefaultPreBasisFactory<Callable0>& factory0, Args&&... args)
{
  using ArgsTuple = std::tuple<DefaultPreBasisFactory<Callable0>, Args...>;
  auto argsTuple = std::make_tuple(factory0, std::forward<Args>(args)...);
  using MergingStrategy = std::tuple_element_t<sizeof...(Args), ArgsTuple>;
  constexpr bool isIndexBlocked = std::is_same_v<MergingStrategy, BlockedTopologicalInterleaving>;

  return DefaultPreBasisFactory{[=](const auto& grid_view, const auto& sub_domain){
    auto sequence = std::make_index_sequence<sizeof...(Args)>{};
    return unpackIntegerSequence(
      [&](auto... i) {
        std::tuple storage = { std::get<i>(argsTuple).protoBasis(grid_view, sub_domain)... };
        return ProtoBasisTuple<isIndexBlocked, decltype(std::get<i>(argsTuple).protoBasis(grid_view, sub_domain))...>{storage};
      },
      sequence);
  }};
}

/**
 * @brief Restrict a proto-basis factory to a sub-domain.
 *
 * This allows you to create a basis that is only defined on a specific sub-domain
 * of the grid, e.g. for domain decomposition or multi-domain problems.
 *
 * @tparam Factory Type of the child proto-basis factory.
 * @tparam SubDomain Type of the sub-domain (must implement `contains` method).
 * @param factory Factory for the child proto-basis.
 * @param sub_domain The sub-domain to restrict the basis to.
 * @return A DefaultPreBasisFactory for the restricted basis.
 */
template<class Factory, class SubDomain>
auto restrict(PreBasisFactory<Factory>&& factory, const SubDomain& sub_domain) {
  return DefaultPreBasisFactory{[factory = static_cast<Factory&&>(factory), sub_domain](const auto& grid_view, EntireDomain){
    return factory.protoBasis(grid_view, sub_domain);
  }};
}

/**
 * @brief Restrict a proto-basis factory to a sub-domain.
 *
 * This allows you to create a basis that is only defined on a specific sub-domain
 * of the grid, e.g. for domain decomposition or multi-domain problems.
 *
 * @tparam Factory Type of the child proto-basis factory.
 * @tparam SubDomain Type of the sub-domain (must implement `contains` method).
 * @param factory Factory for the child proto-basis.
 * @param sub_domain The sub-domain to restrict the basis to.
 * @return A DefaultPreBasisFactory for the restricted basis.
 */
template<class Factory, class SubDomain>
auto restrict(const PreBasisFactory<Factory>& factory, const SubDomain& sub_domain) {
  return DefaultPreBasisFactory{[factory = static_cast<const Factory&>(factory), sub_domain](const auto& grid_view, EntireDomain){
    return factory.protoBasis(grid_view, sub_domain);
  }};
}


/**
 * @brief Factory for constructing PDELab PreBasis objects
 *
 * This class wraps a ProtoBasis factory and provides methods to construct
 * PreBasis objects for a given grid view. It also allows access to the underlying
 * ProtoBasis for advanced use cases.
 *
 * @tparam Base The derived class implementing the ProtoBasis factory.
 */
template<class Base>
struct PreBasisFactory {

  /**
   * @brief Construct a PreBasis for a given grid view.
   * @param grid_view The grid view to build the basis for.
   * @return A PreBasis object.
   */
  auto operator()(const auto& grid_view) && {
    return PreBasis{impl().protoBasis(grid_view, EntireDomain{})};
  }

  /**
   * @brief Construct a PreBasis for a given grid view.
   * @param grid_view The grid view to build the basis for.
   * @return A PreBasis object.
   */
  auto operator()(const auto& grid_view) const & {
    return PreBasis{impl().protoBasis(grid_view, EntireDomain{})};
  }

  //! @brief Creates a power basis factory from a proto-basis factory and a compile-time degree.
  template <std::size_t k>
  auto operator^(index_constant<k>) && {
    return power<k>(impl(), flatByEntity());
  }

  //! @brief Creates a power basis factory from a proto-basis factory and a compile-time degree.
  template <std::size_t k>
  auto operator^(index_constant<k>) const & {
    return power<k>(impl(), flatByEntity());
  }

  //! @brief Creates a power basis factory from a proto-basis factory and a runtime degree.
  auto operator^(std::size_t k) && {
    return power(impl(), k, flatByEntity());
  }

  //! @brief Creates a power basis factory from a proto-basis factory and a runtime degree.
  auto operator^(std::size_t k) const & {
    return power(impl(), k, flatByEntity());
  }

  //! @brief Restricts the proto-basis to a given sub-domain
  template<class SubDomain>
  auto operator|(const SubDomain& sub_domain) && {
    return restrict(impl(), sub_domain);
  }

  //! @brief Restricts the proto-basis to a given sub-domain
  template<class SubDomain>
  auto operator|(const SubDomain& sub_domain) const & {
    return restrict(impl(), sub_domain);
  }

  //! @brief Creates a composite basis factory from two pre-basis factories.
  template<class Factory>
  auto operator*(Factory&& rightFactory) && {
    return composite(impl(), std::forward<Factory>(rightFactory), flatByEntity());
  }

  //! @brief Creates a composite basis factory from two pre-basis factories.
  template<class Factory>
  auto operator*(Factory&& rightFactory) const & {
    return composite(impl(), std::forward<Factory>(rightFactory), flatByEntity());
  }

  //! @brief Creates a global basis for a given grid view using this pre-basis factory.
  template<class GridView>
  auto operator||(const GridView& grid_view) && {
    return Dune::Functions::BasisFactory::makeBasis(grid_view, impl());
  }

  //! @brief Creates a global basis for a given grid view using this pre-basis factory.
  template<class GridView>
  auto operator||(const GridView& grid_view) const & {
    return Dune::Functions::BasisFactory::makeBasis(grid_view, impl());
  }

protected:
  // Ensure that only derived classes construct this base class
  PreBasisFactory() = default;
  PreBasisFactory(const PreBasisFactory&) = default;
  PreBasisFactory& operator=(const PreBasisFactory&) = default;
  PreBasisFactory(PreBasisFactory&&) = default;
  PreBasisFactory& operator=(PreBasisFactory&&) = default;

private:

  //! Access the underlying implementation
  Base&& impl() && {
    return static_cast<Base&&>(*this);
  }

  //! Access the underlying implementation
  const Base& impl() const & {
    return static_cast<const Base&>(*this);
  }
};

/**
 * @brief Default implementation of a ProtoBasis factory.
 *
 * This class wraps a callable object that constructs ProtoBasis objects
 * for given grid views and sub-domains.
 *
 * @tparam Callable The callable type that constructs ProtoBasis objects.
 */
template<class Callable>
struct DefaultPreBasisFactory : public PreBasisFactory<DefaultPreBasisFactory<Callable>> {

  /**
   * @brief Construct from a Callable.
   * @param factory The proto-basis factory to wrap.
   */
  constexpr DefaultPreBasisFactory(Callable&& factory) : proto_basis_factory_{std::move(factory)} {}

  constexpr DefaultPreBasisFactory(const DefaultPreBasisFactory&) = default;
  constexpr DefaultPreBasisFactory& operator=(const DefaultPreBasisFactory&) = default;

  constexpr DefaultPreBasisFactory(DefaultPreBasisFactory&&) = default;
  constexpr DefaultPreBasisFactory& operator=(DefaultPreBasisFactory&&) = default;

  /**
   * @brief Access the underlying ProtoBasis for a given grid view.
   * @param grid_view The grid view to build the proto-basis for.
   * @param sub_domain The sub-domain to restrict the basis to.
   * @return A ProtoBasis object.
   */
  auto protoBasis(const auto& grid_view, const auto& sub_domain) const {
    return proto_basis_factory_(grid_view, sub_domain);
  }

private:
  Callable proto_basis_factory_;
};

namespace Experimental {

//! integral constant for vertex geometry type
inline constexpr std::integral_constant<GeometryType::Id, GeometryTypes::vertex>        vertex;
//! integral constant for line geometry type
inline constexpr std::integral_constant<GeometryType::Id, GeometryTypes::line>          line;
//! integral constant for triangle geometry type
inline constexpr std::integral_constant<GeometryType::Id, GeometryTypes::triangle>      triangle;
//! integral constant for quadrilateral geometry type
inline constexpr std::integral_constant<GeometryType::Id, GeometryTypes::quadrilateral> quadrilateral;
//! integral constant for tetrahedron geometry type
inline constexpr std::integral_constant<GeometryType::Id, GeometryTypes::tetrahedron>   tetrahedron;
//! integral constant for pyramid geometry type
inline constexpr std::integral_constant<GeometryType::Id, GeometryTypes::pyramid>       pyramid;
//! integral constant for prism geometry type
inline constexpr std::integral_constant<GeometryType::Id, GeometryTypes::prism>         prism;
//! integral constant for hexahedron geometry type
inline constexpr std::integral_constant<GeometryType::Id, GeometryTypes::hexahedron>    hexahedron;

//! Base class for finite element family descriptors
struct FiniteElementFamily {};

/**
 * @brief Factory function to create a vector finite element basis for given family, geometry, and degree.
 *
 * @tparam Family Type of the finite element family (must derive from FiniteElementFamily).
 * @tparam Geometry Type representing the geometry (convertible to GeometryType::Id).
 * @tparam Degree Type representing the polynomial degree (convertible to std::size_t).
 * @tparam LeafMergingStrategy Merging strategy for leaf nodes (default is FlatTopologicalInterleaving).
 * @param family The finite element family descriptor.
 * @param geometry The geometry type.
 * @param degree The polynomial degree.
 * @param leaf_strategy The merging strategy for leaf nodes (default is flat).
 * @return A ProtoBasis factory for the specified finite element.
 */
template <
  std::derived_from<FiniteElementFamily> Family,
  std::convertible_to<GeometryType::Id> Geometry,
  std::convertible_to<std::size_t> Degree,
  class VectorMergingStrategy = FlatTopologicalInterleaving,
  class LeafMergingStrategy = FlatTopologicalInterleaving
>
auto VectorElement(const Family& family, Geometry geometry, Degree degree, VectorMergingStrategy vector_strategy = {}, LeafMergingStrategy leaf_strategy = {}) {
  auto element = FiniteElement(family, geometry, degree, leaf_strategy);
  auto d = index_constant<GeometryType{geometry}.dim()>{};
  return power(element, d, vector_strategy);
}

/**
 * @brief Factory function to create a tensor finite element basis for given family, geometry, and degree.
 *
 * @tparam Family Type of the finite element family (must derive from FiniteElementFamily).
 * @tparam Geometry Type representing the geometry (convertible to GeometryType::Id).
 * @tparam Degree Type representing the polynomial degree (convertible to std::size_t).
 * @tparam TensorMergingStrategy Merging strategy for tensor nodes (default is FlatTopologicalInterleaving).
 * @tparam LeafMergingStrategy Merging strategy for leaf nodes (default is FlatTopologicalInterleaving).
 * @param family The finite element family descriptor.
 * @param geometry The geometry type.
 * @param degree The polynomial degree.
 * @param tensor_strategy The merging strategy for tensor nodes (default is flat).
 * @param leaf_strategy The merging strategy for leaf nodes (default is flat).
 * @return A ProtoBasis factory for the specified finite element.
 */
template <
  std::derived_from<FiniteElementFamily> Family,
  std::convertible_to<GeometryType::Id> Geometry,
  std::convertible_to<std::size_t> Degree,
  class TensorMergingStrategy = FlatTopologicalInterleaving,
  class LeafMergingStrategy = FlatTopologicalInterleaving
>
auto TensorElement(const Family& family, Geometry geometry, Degree degree, TensorMergingStrategy tensor_strategy = {}, LeafMergingStrategy leaf_strategy = {}) {
  auto element = VectorElement(family, geometry, degree, tensor_strategy, leaf_strategy);
  auto d = index_constant<GeometryType{geometry}.dim()>{};
  return power(element, d, tensor_strategy);
}

} // namespace Experimental

} // namespace Dune::PDELab::BasisFactory::Experimental

#endif // DUNE_PDELAB_BASIS_FACTORY_HH
