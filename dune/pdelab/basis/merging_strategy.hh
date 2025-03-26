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
 * @tparam CB              A boolean flag indicating whether the container is blocked. This affects how entities
 *                          are merged and indexed.
 *
 * @details
 * The `EntityInterleaving` struct is designed to work with entity sets in the context of finite element spaces. It
 * provides methods to access and update the entity set, allowing for flexible management of entity groupings.
 */
template<bool CB>
struct EntityInterleaving
{
  static constexpr bool Blocked = CB;
};

template<bool CB>
void registerIndexMergingStrategy(EntityInterleaving<CB>);

/**
 * @brief Creates a flat (non-blocked) entity interleaving.
 *
 * This function creates a merging strategy that interleaves the topologic entity
 * index of the degrees of freedom and merges them lexicographically.
 *
 * @return An `EntityInterleaving` with a flat merging strategy.
 */
constexpr auto flatByEntity()
{
  return EntityInterleaving<false>{};
}

/**
 * @brief Creates a blocked entity grouping.
 *
 * This function creates a merging strategy that interleaves the topologic entity
 * index of the degrees of freedom. This means that every topologic entity will
 * have its own (hiererchical) block of degrees of freedom.
 *
 * @return An `EntityInterleaving` with a blocked merging strategy.
 */
constexpr auto blockedByEntity()
{
  return EntityInterleaving<true>{ };
}


} // namespace Dune::PDELab::BasisFactory

#endif // DUNE_PDELAB_BASIS_MERGING_STRATEGY_HH
