#ifndef DUNE_PDELAB_BASIS_MERGING_STRATEGY_HH
#define DUNE_PDELAB_BASIS_MERGING_STRATEGY_HH

namespace Dune::PDELab::BasisFactory {

/**
 * @struct EntityInterleaving
 * @brief Represents a grouping of entities with a specified merging strategy.
 *
 * This struct encapsulates an entity set and a boolean flag indicating whether the container is blocked.
 * It is used to manage and update the entity set, providing a way to define how entities are grouped and merged
 * within a topological associativity forest.
 *
 * @tparam ES              The type of the entity set. This defines the collection of entities being managed.
 * @tparam CB              A boolean flag indicating whether the container is blocked. This affects how entities
 *                          are merged and indexed.
 *
 * @details
 * The `EntityInterleaving` struct is designed to work with entity sets in the context of finite element spaces. It
 * provides methods to access and update the entity set, allowing for flexible management of entity groupings.
 */
template<class ES, bool CB>
struct EntityInterleaving
{
  using EntitySet = ES;
  static constexpr bool Blocked = CB;

  /**
   * @brief Constructs an EntityInterleaving with a given entity set.
   *
   * @param entity_set The entity set to be managed by this grouping.
   */
  EntityInterleaving(EntitySet entity_set)
    : _entity_set{ std::move(entity_set) }
  {
  }

  /**
   * @brief Accesses the entity set associated with this grouping.
   *
   * @return A reference to the entity set.
   */
  EntitySet& entitySet() { return _entity_set; }

  /**
   * @brief Accesses the entity set associated with this grouping (const version).
   *
   * @return A const reference to the entity set.
   */
  const EntitySet& entitySet() const { return _entity_set; }

  /**
   * @brief Updates the entity set associated with this strategy.
   *
   * @param entity_set The new entity set to assign to this strategy.
   */
  void update(EntitySet entity_set) {
    _entity_set = entity_set;
  }

private:
  EntitySet _entity_set;
};

template<class ES, bool CB>
void registerIndexMergingStrategy(EntityInterleaving<ES, CB>);

/**
 * @brief Creates a flat (non-blocked) entity interleaving.
 *
 * This function creates a merging strategy that interleaves the topologic entity
 * index of the degrees of freedom and merges them lexicographically.
 *
 * @tparam EntitySet The type of the entity set.
 * @param entity_Set The entity set to be used for the merging strategy.
 * @return An `EntityInterleaving` with a flat merging strategy.
 */
template<class EntitySet>
constexpr auto flatByEntity(const EntitySet& entity_Set)
{
  return EntityInterleaving<EntitySet, false>{ entity_Set };
}

/**
 * @brief Creates a blocked entity grouping.
 *
 * This function creates a merging strategy that interleaves the topologic entity
 * index of the degrees of freedom. This means that every topologic entity will
 * have its own (hiererchical) block of degrees of freedom.
 *
 * @tparam EntitySet The type of the entity set.
 * @param entity_Set The entity set to be used for the merging strategy.
 * @return An `EntityInterleaving` with a blocked merging strategy.
 */
template<class EntitySet>
constexpr auto blockedByEntity(const EntitySet& entity_Set)
{
  return EntityInterleaving<EntitySet, true>{ entity_Set };
}


} // namespace Dune::PDELab::BasisFactory

#endif // DUNE_PDELAB_BASIS_MERGING_STRATEGY_HH
