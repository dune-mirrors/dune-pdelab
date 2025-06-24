#ifndef DUNE_PDELAB_BASIS_MERGING_STRATEGY_HH
#define DUNE_PDELAB_BASIS_MERGING_STRATEGY_HH

namespace Dune::PDELab::BasisFactory {

/**
 * @ingroup FunctionSpaceBases
 * @brief Descriptor of the flat merging strategy for grouping degrees of freedom associated with topological entities.
 * @details This merging strategy instructs the basis objects to interleave the indices associated with the same topological entity.
 * In this form, degrees of freedom associated with the same topological entity are grouped together, and lexicographically merged into a flat index.
 *
 * @see ProtoBasisLeaf, ProtoBasisArray, ProtoBasisVector, ProtoBasisTuple, BasisTopologicalInterleaving
 *
 */
struct FlatTopologicalInterleaving {};

/**
 * @ingroup FunctionSpaceBases
 * @brief Descriptor of the blocked merging strategy for grouping degrees of freedom associated with topological entities.
 * @details This merging strategy instructs the basis objects to interleave the indices associated with the same topological entity.
 * In this form, degrees of freedom associated with the same topological entity are grouped together.
 *
 * @see ProtoBasisLeaf, ProtoBasisArray, ProtoBasisVector, ProtoBasisTuple, BasisTopologicalInterleaving
 *
 */
struct BlockedTopologicalInterleaving {};


#ifndef DOXYGEN
//! Register FlatTopologicalInterleaving as a index merging strategy in dune-functions
void registerIndexMergingStrategy(FlatTopologicalInterleaving);
//! Register BlockedTopologicalInterleaving as a index merging strategy in dune-functions
void registerIndexMergingStrategy(BlockedTopologicalInterleaving);
#endif

/**
 * @brief Creates a flat (non-blocked) entity interleaving descriptor.
 * @ingroup FunctionSpaceBases
 *
 * @return An `FlatTopologicalInterleaving` with a flat merging strategy.
 */
constexpr FlatTopologicalInterleaving flatByEntity()
{
  return {};
}

/**
 * @brief Creates a blocked entity grouping descriptor.
 * @ingroup FunctionSpaceBases
 *
 * @return An `BlockedTopologicalInterleaving` with a blocked merging strategy.
 */
constexpr BlockedTopologicalInterleaving blockedByEntity()
{
  return {};
}

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_BASIS_MERGING_STRATEGY_HH
