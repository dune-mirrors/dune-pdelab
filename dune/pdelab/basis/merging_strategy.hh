#ifndef DUNE_ASSEMBLER_DISCRETE_FUNCTION_SPACE_MERGING_STRATEGY_HH
#define DUNE_ASSEMBLER_DISCRETE_FUNCTION_SPACE_MERGING_STRATEGY_HH

//! maximum grid world dimension
#ifndef DUNE_ASSEMBLER_MAX_WORLDDIM
#define DUNE_ASSEMBLER_MAX_WORLDDIM 32
#endif

#include <dune/grid/concepts/gridview.hh>

#include <bitset>

namespace Dune::Assembler {
inline namespace Strategy {

template<bool ContainerBlocked>
struct DefaultStrategy
{
  using SizeType = std::size_t;
  using CodimFlag = std::bitset<DUNE_ASSEMBLER_MAX_WORLDDIM>;
  static constexpr bool Blocked = ContainerBlocked;
};

template<bool ContainerBlocked>
void registerIndexMergingStrategy(DefaultStrategy<ContainerBlocked>);

template<bool Blocked>
struct Lexicographic : public DefaultStrategy<Blocked> {};

using FlatLexicographic = Lexicographic<false>;
using BlockedLexicographic = Lexicographic<true>;

constexpr FlatLexicographic
flatLexicographic()
{
  return {};
}
constexpr BlockedLexicographic
blockedLexicographic()
{
  return {};
}


template<Dune::Concept::GridView ES, bool ContainerBlocked>
struct EntityGrouping : public DefaultStrategy<ContainerBlocked>
{
  using EntitySet = ES;

  EntityGrouping(EntitySet entity_set)
    : _entity_set{ std::move(entity_set) }
  {
  }

  EntitySet& entitySet() { return _entity_set; }
  const EntitySet& entitySet() const { return _entity_set; }

private:
  EntitySet _entity_set;
};

template<Dune::Concept::GridView EntitySet>
constexpr auto
flatByEntity(const EntitySet& entity_Set)
{
  return EntityGrouping<EntitySet, false>{ entity_Set };
}

template<Dune::Concept::GridView EntitySet>
constexpr auto
blockedByEntity(const EntitySet& entity_Set)
{
  return EntityGrouping<EntitySet, true>{ entity_Set };
}

} // namespace Strategy

} // namespace Dune::Assembler

#endif // DUNE_ASSEMBLER_DISCRETE_FUNCTION_SPACE_MERGING_STRATEGY_HH
