#ifndef DUNE_PDELAB_BASIS_MERGING_STRATEGY_HH
#define DUNE_PDELAB_BASIS_MERGING_STRATEGY_HH

namespace Dune::PDELab::BasisFactory {

template<class ES, bool ContainerBlocked>
struct EntityGrouping
{
  using EntitySet = ES;
  static constexpr bool Blocked = ContainerBlocked;

  EntityGrouping(EntitySet entity_set)
    : _entity_set{ std::move(entity_set) }
  {
  }

  EntitySet& entitySet() { return _entity_set; }
  const EntitySet& entitySet() const { return _entity_set; }

  void update(EntitySet entity_set) {
    _entity_set = entity_set;
  }

  friend void registerIndexMergingStrategy(EntityGrouping);

private:
  EntitySet _entity_set;
};

template<class EntitySet>
constexpr auto
flatByEntity(const EntitySet& entity_Set)
{
  return EntityGrouping<EntitySet, false>{ entity_Set };
}

template<class EntitySet>
constexpr auto
blockedByEntity(const EntitySet& entity_Set)
{
  return EntityGrouping<EntitySet, true>{ entity_Set };
}

} // namespace Dune::PDELab::BasisFactory

#endif // DUNE_PDELAB_BASIS_MERGING_STRATEGY_HH
