// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDOPERATOR_COMMON_MULTIRESIDUAL_HH
#define DUNE_PDELAB_GRIDOPERATOR_COMMON_MULTIRESIDUAL_HH

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

/** \file
 * \brief Residuals for several coefficient vectors in a single grid traversal.
 *
 * Impl::accumulateResiduals(go, xs, r, prepare) adds
 * \f$ \sum_j r(x_j) \f$ to r, where prepare(j) is called before the
 * contributions of xs[j] are evaluated (e.g. to set the time and the weights
 * of the local operator). All contributions share one grid traversal, i.e.
 * the elements, intersections and local function spaces are bound once, and
 * the constraints are applied once. Used by OneStepGridOperator for the
 * constant part of a stage, which collects the previous stages.
 *
 * The overload in this file covers the grid operators built on local
 * assembler engines (GridOperator and FastDGGridOperator on grid function
 * spaces); the GridOperator on dune-functions bases provides its own.
 */

namespace Dune::PDELab::Impl {

  /** \brief Local assembler engine combining residual engines for several solutions
   *
   * Forwards the traversal callbacks to one residual engine per solution;
   * before engine j assembles, prepare(j) is called.
   *
   * \tparam Engine  residual engine of the local assembler (LocalResidualAssemblerEngine)
   * \tparam Prepare callable prepare(std::size_t j)
   */
  template<class Engine, class Prepare>
  class MultiSolutionResidualEngine
  {
  public:
    template<class LocalAssembler, class Solution, class Residual>
    MultiSolutionResidualEngine(const LocalAssembler& localAssembler, const std::vector<const Solution*>& solutions,
                                Residual& residual, Prepare& prepare)
      : prepare_(&prepare)
    {
      // engines refer to their own members, hence they are not copied
      for (const Solution* solution : solutions)
      {
        engines_.push_back(std::make_unique<Engine>(localAssembler));
        engines_.back()->setResidual(residual);
        engines_.back()->setSolution(*solution);
      }
    }

    template<class CU, class CV>
    bool needsConstraintsCaching(const CU& cu, const CV& cv) const
    {
      return engines_.front()->needsConstraintsCaching(cu, cv);
    }

    auto partition() const { return engines_.front()->partition(); }

    template<class EG>
    bool skipEntity(const EG& eg) { return engines_.front()->skipEntity(eg); }

    template<class IG>
    bool skipIntersection(const IG& ig) { return engines_.front()->skipIntersection(ig); }

    bool requireSkeleton() const { return engines_.front()->requireSkeleton(); }
    bool requireSkeletonTwoSided() const { return engines_.front()->requireSkeletonTwoSided(); }
    bool requireUVVolume() const { return engines_.front()->requireUVVolume(); }
    bool requireVVolume() const { return engines_.front()->requireVVolume(); }
    bool requireUVSkeleton() const { return engines_.front()->requireUVSkeleton(); }
    bool requireVSkeleton() const { return engines_.front()->requireVSkeleton(); }
    bool requireUVBoundary() const { return engines_.front()->requireUVBoundary(); }
    bool requireVBoundary() const { return engines_.front()->requireVBoundary(); }
    bool requireUVProcessor() const { return engines_.front()->requireUVProcessor(); }
    bool requireVProcessor() const { return engines_.front()->requireVProcessor(); }
    bool requireUVEnrichedCoupling() const { return engines_.front()->requireUVEnrichedCoupling(); }
    bool requireVEnrichedCoupling() const { return engines_.front()->requireVEnrichedCoupling(); }
    bool requireUVVolumePostSkeleton() const { return engines_.front()->requireUVVolumePostSkeleton(); }
    bool requireVVolumePostSkeleton() const { return engines_.front()->requireVVolumePostSkeleton(); }

    void preAssembly()
    {
      for (auto& engine : engines_)
        engine->preAssembly();
    }

    // all engines add into the same residual; constrain it once
    template<class GFSU, class GFSV>
    void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
    {
      engines_.back()->postAssembly(gfsu, gfsv);
    }

    // binding and loading of coefficients: every engine binds its solution view
#define DUNE_PDELAB_FORWARD_TO_ALL(name)                \
    template<class... Args>                             \
    void name(const Args&... args)                      \
    {                                                   \
      for (auto& engine : engines_)                     \
        engine->name(args...);                          \
    }

    DUNE_PDELAB_FORWARD_TO_ALL(onBindLFSUV)
    DUNE_PDELAB_FORWARD_TO_ALL(onBindLFSV)
    DUNE_PDELAB_FORWARD_TO_ALL(onBindLFSUVInside)
    DUNE_PDELAB_FORWARD_TO_ALL(onBindLFSUVOutside)
    DUNE_PDELAB_FORWARD_TO_ALL(onBindLFSVInside)
    DUNE_PDELAB_FORWARD_TO_ALL(onBindLFSVOutside)
    DUNE_PDELAB_FORWARD_TO_ALL(onBindLFSUVCoupling)
    DUNE_PDELAB_FORWARD_TO_ALL(onBindLFSVCoupling)
    DUNE_PDELAB_FORWARD_TO_ALL(onUnbindLFSUV)
    DUNE_PDELAB_FORWARD_TO_ALL(onUnbindLFSV)
    DUNE_PDELAB_FORWARD_TO_ALL(onUnbindLFSUVInside)
    DUNE_PDELAB_FORWARD_TO_ALL(onUnbindLFSUVOutside)
    DUNE_PDELAB_FORWARD_TO_ALL(onUnbindLFSVInside)
    DUNE_PDELAB_FORWARD_TO_ALL(onUnbindLFSVOutside)
    DUNE_PDELAB_FORWARD_TO_ALL(onUnbindLFSUVCoupling)
    DUNE_PDELAB_FORWARD_TO_ALL(onUnbindLFSVCoupling)
    DUNE_PDELAB_FORWARD_TO_ALL(loadCoefficientsLFSUInside)
    DUNE_PDELAB_FORWARD_TO_ALL(loadCoefficientsLFSUOutside)
    DUNE_PDELAB_FORWARD_TO_ALL(loadCoefficientsLFSUCoupling)

#undef DUNE_PDELAB_FORWARD_TO_ALL

    // assembly: engine j evaluates the local operator after prepare(j)
#define DUNE_PDELAB_FORWARD_PREPARED(name)              \
    template<class... Args>                             \
    void name(const Args&... args)                      \
    {                                                   \
      for (std::size_t j = 0; j < engines_.size(); ++j) \
      {                                                 \
        (*prepare_)(j);                                 \
        engines_[j]->name(args...);                     \
      }                                                 \
    }

    DUNE_PDELAB_FORWARD_PREPARED(assembleUVVolume)
    DUNE_PDELAB_FORWARD_PREPARED(assembleVVolume)
    DUNE_PDELAB_FORWARD_PREPARED(assembleUVSkeleton)
    DUNE_PDELAB_FORWARD_PREPARED(assembleVSkeleton)
    DUNE_PDELAB_FORWARD_PREPARED(assembleUVBoundary)
    DUNE_PDELAB_FORWARD_PREPARED(assembleVBoundary)
    DUNE_PDELAB_FORWARD_PREPARED(assembleUVProcessor)
    DUNE_PDELAB_FORWARD_PREPARED(assembleVProcessor)
    DUNE_PDELAB_FORWARD_PREPARED(assembleUVEnrichedCoupling)
    DUNE_PDELAB_FORWARD_PREPARED(assembleVEnrichedCoupling)
    DUNE_PDELAB_FORWARD_PREPARED(assembleUVVolumePostSkeleton)
    DUNE_PDELAB_FORWARD_PREPARED(assembleVVolumePostSkeleton)

#undef DUNE_PDELAB_FORWARD_PREPARED

  private:
    std::vector<std::unique_ptr<Engine>> engines_;
    Prepare* prepare_;
  };

  /** \brief Add the residuals of several solutions to r in a single grid traversal
   *
   * Grid operators on grid function spaces (GridOperator, FastDGGridOperator).
   * prepare(j) is called before the contributions of xs[j] are evaluated.
   */
  template<class GO, class Prepare>
    requires requires(const GO& go) {
      typename GO::Traits::LocalAssembler::LocalResidualAssemblerEngine;
      go.assembler();
    }
  void accumulateResiduals(const GO& go, const std::vector<const typename GO::Traits::Domain*>& xs,
                           typename GO::Traits::Range& r, Prepare& prepare)
  {
    if (xs.empty())
      return;
    using Engine = typename GO::Traits::LocalAssembler::LocalResidualAssemblerEngine;
    MultiSolutionResidualEngine<Engine, Prepare> engine(go.localAssembler(), xs, r, prepare);
    go.assembler().assemble(engine);
  }

} // namespace Dune::PDELab::Impl

#endif // DUNE_PDELAB_GRIDOPERATOR_COMMON_MULTIRESIDUAL_HH
