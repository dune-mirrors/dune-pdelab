#ifndef DUNE_PDELAB_BASIS_PROTOBASIS_LEAF_HH
#define DUNE_PDELAB_BASIS_PROTOBASIS_LEAF_HH

#include <dune/pdelab/basis/protobasis/node.hh>
#include <dune/pdelab/constraints/noconstraints.hh>

#include <dune/typetree/leafnode.hh>

#include <memory>
#include <utility>

namespace Dune::PDELab {

/**
 * @brief Leaf proto-basis node
 *
 * @tparam MS   A merging strategy
 * @tparam FEM  A Finite element map
 * @tparam CON  Constraints operator
 */
template<class MergingStrategy_, class FiniteElementMap_, class ConstraintsOperator_ = NoConstraints>
class ProtoBasis
  : public TypeTree::LeafNode
  , public ProtoBasisNode<MergingStrategy_>
{
  using BaseNode = ProtoBasisNode<MergingStrategy_>;
  using TreeNode = TypeTree::LeafNode;

public:
  using MergingStrategy = MergingStrategy_;
  using FiniteElementMap = FiniteElementMap_;
  using ConstraintsOperator = ConstraintsOperator_;

  /**
   * @brief Construct a new leaf proto-basis object
   * @note Use with value semantics
   *
   * @param merging_strategy  Node merging strategy
   * @param fem               Pointer to a finite element map
   * @param constraints_op    constraints operator
   */
  ProtoBasis(const MergingStrategy& merging_strategy,
            std::shared_ptr<const FiniteElementMap> fem,
            const ConstraintsOperator& constraints_op = {})
    : BaseNode{ merging_strategy }
    , _finite_element_map{ std::move(fem) }
    , _constraints_op{ constraints_op }
  {
  }

  //! Get finite element map
  const FiniteElementMap& finiteElementMap() const
  {
    return *_finite_element_map;
  }

  ConstraintsOperator constraintsOperator() const
  {
    return _constraints_op;
  }

private:
  std::shared_ptr<const FiniteElementMap> _finite_element_map;
  ConstraintsOperator _constraints_op;
};

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_BASIS_PROTOBASIS_LEAF_HH
