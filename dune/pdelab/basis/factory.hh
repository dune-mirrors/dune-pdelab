#ifndef DUNE_PDELAB_BASIS_FACTORY_HH
#define DUNE_PDELAB_BASIS_FACTORY_HH

#include <dune/pdelab/basis/prebasis.hh>
#include <dune/pdelab/basis/merging_strategy.hh>

namespace Dune::PDELab {

namespace BasisFactory {

template<class ProtoBasisFactory>
auto preBasis(const ProtoBasisFactory& pre_basis_factory) {
  return [=]<class GridView>(const GridView& grid_view) {
    auto proto_basis = pre_basis_factory(grid_view);
    using ProtoBasis = decltype(proto_basis);
    return EntityInterleavedPreBasis{proto_basis, grid_view};
  };
}

template<class Factory>
struct ProtoBasisFactory : public Factory {

  ProtoBasisFactory(Factory&& factory) : Factory{std::move(factory)} {}

  auto operator()(const auto& grid_view) const {
    return preBasis(static_cast<const Factory&>(*this))(grid_view);
  }

  auto protoBasis(const auto& grid_view) const {
    return static_cast<const Factory&>(*this)(grid_view);
  }

  friend void registerProtoBasisFactory(ProtoBasisFactory);
};

static constexpr auto restrict = [](auto&& pre_basis_factory, auto sub_domaingrid_view_) {
  return [=](const auto& grid_view){
    return pre_basis_factory(sub_domaingrid_view_);
  };
};

template<std::size_t k, Concept::ProtoBasisFactory ChildProtoBasisFactory, bool blocked>
auto power(ChildProtoBasisFactory childPreBasisFactory, const EntityInterleaving<blocked>&){
  return ProtoBasisFactory{[childPreBasisFactory](const auto& grid_view){
    auto proto_basis = childPreBasisFactory.protoBasis(grid_view);
    using ProtoBasis = decltype(proto_basis);
    std::array<std::shared_ptr<ProtoBasis>, k> storage;
    for (std::size_t i = 0; i < k; ++i)
      storage[i] = std::make_shared<ProtoBasis>(proto_basis);
    using PB = Dune::PDELab::ProtoBasisArray<EntityInterleaving<blocked>, ProtoBasis, k>;
    return PB{ EntityInterleaving<blocked>{}, storage };
  }};
}

template<Concept::ProtoBasisFactory ChildProtoBasisFactory, bool blocked>
auto power(ChildProtoBasisFactory&& childPreBasisFactory, std::size_t k, const EntityInterleaving<blocked>&)
{
  return ProtoBasisFactory{[childPreBasisFactory, k](const auto& grid_view){
    auto proto_basis = childPreBasisFactory.protoBasis(grid_view);
    using ProtoBasis = decltype(proto_basis);
    std::vector<std::shared_ptr<ProtoBasis>> storage(k);
    for (std::size_t i = 0; i < k; ++i)
      storage[i] = std::make_shared<ProtoBasis>(proto_basis);
    using PB = Dune::PDELab::ProtoBasisVector<EntityInterleaving<blocked>, ProtoBasis>;
    return PB{ EntityInterleaving<blocked>{}, storage };
  }};
}

template<typename... Args>
  requires (std::convertible_to<std::tuple_element_t<sizeof...(Args)-1, std::tuple<Args...>>, EntityInterleaving<true>>)
        || (std::convertible_to<std::tuple_element_t<sizeof...(Args)-1, std::tuple<Args...>>, EntityInterleaving<false>>)
auto composite(Args&&... args)
{
  using TypeTuple = std::tuple<Args...>;
  using MergingStrategy = std::tuple_element_t<(sizeof...(Args)-1), TypeTuple>;

  return ProtoBasisFactory{[=](const auto& grid_view){
    using Storage = std::tuple<std::shared_ptr<Args>...>;
    auto sequence = std::make_index_sequence<sizeof...(Args)-1>{};
    return unpackIntegerSequence(
      [&](auto... i) {
        auto storage = Storage{ std::make_shared<std::tuple_element_t<i, TypeTuple>>(std::get<i>(args).protoBasis(grid_view))... };
        return Dune::PDELab::ProtoBasisTuple{MergingStrategy{}, storage};
      },
      sequence);
  }};
}

} // namespace BasisFactory

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_BASIS_FACTORY_HH
