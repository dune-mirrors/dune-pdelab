#ifndef DUNE_PDELAB_BASIS_PROTOBASIS_NODE_HH
#define DUNE_PDELAB_BASIS_PROTOBASIS_NODE_HH

#include <string>
#include <utility>

namespace Dune::PDELab {

/**
 * @brief Base node for proto-basis
 *
 * @tparam NodeTraits Traits for this node
 */
template<class MergingStrategy_>
class ProtoBasisNode
{
public:
  using MergingStrategy = MergingStrategy_;

  /**
   * @brief Construct a new proto-basis node object
   *
   * @param merging_strategy  a merging strategy
   */
  ProtoBasisNode(const MergingStrategy merging_strategy)
    : _merging_strategy{ merging_strategy }
  {
  }

  //! Get name of the node
  const std::string& name() const { return _name; }

  //! Set a name for the node
  void name(std::string_view name) { _name = name; }

  //! Get merging strategy
  MergingStrategy& mergingStrategy()
  {
    return _merging_strategy;
  }

  //! Get merging strategy
  const MergingStrategy& mergingStrategy() const
  {
    return _merging_strategy;
  }

private:
  std::string _name;
  MergingStrategy _merging_strategy;
};

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_BASIS_PROTOBASIS_NODE_HH
