#ifndef DUNE_PDELAB_BASIS_PROTOBASIS_HH
#define DUNE_PDELAB_BASIS_PROTOBASIS_HH

#include <dune/pdelab/basis/protobasis/leaf.hh>
#include <dune/pdelab/constraints/noconstraints.hh>

#include <memory>

namespace Dune::PDELab::BasisFactory {

/**
 * @brief Makes a proto basis leaf node
 * A proto-basis is a type descriptor of the tree structure of the final basis.
 *
 * @tparam MergingStrategy            Merging strategy
 * @tparam FiniteElementMapStorage    Shared pointer to a finite element map
 * @tparam ConstraintsOperator                ConstraintsOperator
 * @param merging_strategy            Merging strategy
 * @param fem                         Shared pointer to a finite element map
 * @param con                         ConstraintsOperator
 * @return auto                       The new discrete function space node
 */
template<class MergingStrategy,
         class FiniteElementMap,
         class ConstraintsOperator>
auto
makeProtoBasis(const MergingStrategy& merging_strategy,
              std::shared_ptr<FiniteElementMap> fem,
              const ConstraintsOperator& constraints_op)
{
  return ProtoBasis<MergingStrategy, std::remove_const_t<FiniteElementMap>, ConstraintsOperator>{ merging_strategy, std::move(fem), constraints_op };
}

template<class MergingStrategy,
         class FiniteElementMap>
auto
makeProtoBasis(const MergingStrategy& merging_strategy,
               std::shared_ptr<FiniteElementMap> fem)
{
  return makeProtoBasis(merging_strategy, std::move(fem), NoConstraints{});
}

} // namespace Dune::PDELab::BasisFactory

#endif // DUNE_PDELAB_BASIS_PROTOBASIS_HH
