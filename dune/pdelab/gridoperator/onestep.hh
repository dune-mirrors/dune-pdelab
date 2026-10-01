// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDOPERATOR_ONESTEP_HH
#define DUNE_PDELAB_GRIDOPERATOR_ONESTEP_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/gridoperator/common/gridoperatorutilities.hh>
#include <dune/pdelab/gridoperator/common/multiresidual.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/instationary/onestepparameter.hh>
#include <dune/pdelab/localoperator/weightedsum.hh>

namespace Dune{
  namespace PDELab{

    template<typename GFSU, typename GFSV, typename LOP,
             typename MB, typename DF, typename RF, typename JF,
             typename CU, typename CV>
    class FastDGGridOperator;

    namespace Impl {

      //! The grid operator of the same kind as GO, but with the local operator LOP
      template<typename GO, typename LOP>
      struct ReplaceLocalOperator;

      template<typename GFSU, typename GFSV, typename LOP0,
               typename MB, typename DF, typename RF, typename JF,
               typename CU, typename CV, typename LOP>
      struct ReplaceLocalOperator<GridOperator<GFSU,GFSV,LOP0,MB,DF,RF,JF,CU,CV>, LOP>
      {
        using type = GridOperator<GFSU,GFSV,LOP,MB,DF,RF,JF,CU,CV>;
      };

      template<typename GFSU, typename GFSV, typename LOP0,
               typename MB, typename DF, typename RF, typename JF,
               typename CU, typename CV, typename LOP>
      struct ReplaceLocalOperator<FastDGGridOperator<GFSU,GFSV,LOP0,MB,DF,RF,JF,CU,CV>, LOP>
      {
        using type = FastDGGridOperator<GFSU,GFSV,LOP,MB,DF,RF,JF,CU,CV>;
      };

    } // namespace Impl

    /** \brief Grid operator for one step methods, combining a spatial and a temporal grid operator
     *
     * Assembles the stage equations of a one step method
     * \f[ \sum_{s\le r} b_{rs} \Delta t\, r_0(x_s,t_s) + a_{rs}\, r_1(x_s,t_s) \f]
     * from the local operators of GO0 (spatial part \f$ r_0 \f$) and GO1
     * (temporal part \f$ r_1 \f$, e.g. the mass term). The two local
     * operators are combined in a WeightedSumLocalOperator, and all
     * assembly is delegated to a grid operator of the same kind as GO0
     * with this sum operator. So the stage contributions of both operators
     * are assembled in a single grid traversal, the constraints are applied
     * once, and the same implementation works for every grid operator
     * (grid function spaces, FastDG, dune-functions bases).
     *
     * GO0 and GO1 provide their local operators, GO0 also the spaces, the
     * constraints and the matrix backend. They are not used for assembly.
     *
     * \tparam GO0      grid operator of the spatial part (zeroth time derivative)
     * \tparam GO1      grid operator of the temporal part (first time derivative)
     * \tparam implicit whether the operator is used by an implicit or an explicit one step method
     */
    template<typename GO0, typename GO1, bool implicit = true>
    class OneStepGridOperator
    {
      using LOP0 = typename GO0::Traits::LocalAssembler::LocalOperator;
      using LOP1 = typename GO1::Traits::LocalAssembler::LocalOperator;

    public:
      //! The type for real numbers, e.g. time
      using Real = typename GO0::Traits::RangeField;

      //! The type of the one step method parameters
      using OneStepParameters = TimeSteppingParameterInterface<Real>;

    private:
      using SumLocalOperator = WeightedSumLocalOperator<Real, LOP0, LOP1>;
      using InnerGridOperator = typename Impl::ReplaceLocalOperator<GO0, SumLocalOperator>::type;

    public:
      //! The grid operator traits
      using Traits = typename InnerGridOperator::Traits;

      using Assembler = typename Traits::Assembler;
      using LocalAssembler = typename Traits::LocalAssembler;

      //! The io types of the operator
      //! @{
      using Domain = typename Traits::Domain;
      using Range = typename Traits::Range;
      using Jacobian = typename Traits::Jacobian;
      //! @}

      template <typename MFT>
      struct MatrixContainer
      {
        typedef Jacobian Type;
      };

      //! Construct from the spatial and the temporal grid operator
      OneStepGridOperator(GO0& go0, GO1& go1)
        : go0_(go0), go1_(go1)
        , lop0_(go0.localAssembler().localOperator())
        , lop1_(go1.localAssembler().localOperator())
        , sum_(lop0_, lop1_)
        , inner_(go0.trialGridFunctionSpace(), go0.localAssembler().trialConstraints(),
                 go0.testGridFunctionSpace(), go0.localAssembler().testConstraints(),
                 sum_, go0.matrixBackend())
        , constResidual_(go0.testGridFunctionSpace())
        , dtMode_(implicit ? MultiplyOperator0ByDT : DoNotAssembleDT)
      {}

      // the sum operator refers to the local operators, the inner grid operator to the sum operator
      OneStepGridOperator(const OneStepGridOperator&) = delete;
      OneStepGridOperator& operator=(const OneStepGridOperator&) = delete;

      //! Determines whether the time step size is multiplied to the
      //! mass term (first order time derivative) or the elliptic term
      //! (zero-th order time derivative).
      void divideMassTermByDeltaT()
      {
        if (not implicit)
          DUNE_THROW(Dune::Exception, "This function should not be called in explicit mode");
        dtMode_ = DivideOperator1ByDT;
      }

      void multiplySpatialTermByDeltaT()
      {
        if (not implicit)
          DUNE_THROW(Dune::Exception, "This function should not be called in explicit mode");
        dtMode_ = MultiplyOperator0ByDT;
      }

      const typename Traits::TrialGridFunctionSpace& trialGridFunctionSpace() const
      {
        return inner_.trialGridFunctionSpace();
      }

      const typename Traits::TestGridFunctionSpace& testGridFunctionSpace() const
      {
        return inner_.testGridFunctionSpace();
      }

      auto globalSizeU() const { return inner_.globalSizeU(); }
      auto globalSizeV() const { return inner_.globalSizeV(); }

      //! The global assembler of the underlying grid operator (grid function spaces only)
      decltype(auto) assembler() const
        requires requires(const InnerGridOperator& go) { go.assembler(); }
      {
        return inner_.assembler();
      }

      //! The local assembler of the underlying grid operator: constraints and the combined local operator
      LocalAssembler& localAssembler() const { return inner_.localAssembler(); }

      //! Fill the sparsity pattern of the jacobian (grid function spaces)
      template<typename Pattern>
        requires requires(const InnerGridOperator& go, Pattern& p) { go.fill_pattern(p); }
      void fill_pattern(Pattern& p) const
      {
        setPatternWeights();
        inner_.fill_pattern(p);
      }

      //! Set up the sparsity pattern of a native matrix (dune-functions bases)
      template<typename NativeMatrix>
        requires requires(const InnerGridOperator& go, NativeMatrix& m) { go.setupMatrixPattern(m); }
      void setupMatrixPattern(NativeMatrix& matrix) const
      {
        setPatternWeights();
        inner_.setupMatrixPattern(matrix);
      }

      //! Assemble the constant part of the residual of stage r from the solutions x[0], ..., x[r-1]
      void preStage(unsigned int stage, const std::vector<Domain*>& x)
      {
        if (not implicit)
          DUNE_THROW(Dune::Exception, "This function should not be called in explicit mode");

        stage_ = stage;
        lop0_.preStage(timeAtStage(stage), stage);
        lop1_.preStage(timeAtStage(stage), stage);

        // all previous stages in one grid traversal
        std::vector<Contribution> contributions;
        for (unsigned int s = 0; s < stage; ++s)
        {
          const Real b = method().b(stage, s);
          const Real a = method().a(stage, s);
          if (nonzero(b) or nonzero(a))
            contributions.push_back({x[s], timeAtStage(s),
                                     nonzero(b) ? b * dtFactor0_ : 0.0, nonzero(a) ? a * dtFactor1_ : 0.0});
        }
        constResidual_ = 0.0;
        accumulate(contributions, constResidual_);
      }

      //! Assemble the residual of the current stage
      void residual(const Domain& x, Range& r) const
      {
        if (not implicit)
          DUNE_THROW(Dune::Exception, "This function should not be called in explicit mode");
        setStageWeights();
        inner_.residual(x, r);
        // the constant part is already constrained
        r += constResidual_;
      }

      //! Assemble the jacobian of the current stage
      void jacobian(const Domain& x, Jacobian& a) const
      {
        if (not implicit)
          DUNE_THROW(Dune::Exception, "This function should not be called in explicit mode");
        setStageWeights();
        inner_.jacobian(x, a);
      }

      //! Apply the jacobian of the current stage to update without assembling it (linear operators)
      void jacobian_apply(const Domain& update, Range& result) const
      {
        setStageWeights();
        inner_.jacobian_apply(update, result);
      }

      //! Apply the jacobian of the current stage at solution to update without assembling it
      void jacobian_apply(const Domain& solution, const Domain& update, Range& result) const
      {
        setStageWeights();
        inner_.jacobian_apply(solution, update, result);
      }

      /** \brief Assemble the stage of an explicit method
       *
       * \param stage the stage r
       * \param x     the solutions x[0], ..., x[r]
       * \param a     minus the jacobian of the temporal part at x[r]
       * \param r1    sum of the temporal parts a_{rs} r_1(x_s) of the previous stages
       * \param r0    sum of the spatial parts b_{rs} r_0(x_s) of the previous stages
       */
      void explicit_jacobian_residual(unsigned int stage, const std::vector<Domain*>& x,
                                      Jacobian& a, Range& r1, Range& r0)
      {
        if (implicit)
          DUNE_THROW(Dune::Exception, "This function should not be called in implicit mode");

        stage_ = stage;
        lop0_.preStage(timeAtStage(stage), stage);
        lop1_.preStage(timeAtStage(stage), stage);

        // spatial and temporal parts of the previous stages, one grid traversal each
        // TODO: with a system assembly (a matrix and two vectors) the stage would
        // need a single traversal instead of three.
        std::vector<Contribution> spatial, temporal;
        for (unsigned int s = 0; s < stage; ++s)
        {
          const Real b = method().b(stage, s);
          const Real a_rs = method().a(stage, s);
          if (nonzero(b))
            spatial.push_back({x[s], timeAtStage(s), b * dtFactor0_, 0.0});
          if (nonzero(a_rs))
            temporal.push_back({x[s], timeAtStage(s), 0.0, a_rs * dtFactor1_});
        }
        r0 = 0.0;
        r1 = 0.0;
        accumulate(spatial, r0);
        accumulate(temporal, r1);

        lop1_.setTime(timeAtStage(stage));
        setWeights(0.0, -1.0);
        inner_.jacobian(*x[stage], a);
      }

      //! Interpolate the constrained values of stage `stage` from f, the others from xold
      template<typename F, typename X>
      void interpolate(unsigned stage, const X& xold, F& f, X& x) const
      {
        const Real time = timeAtStage(stage);
        f.setTime(time);
        lop0_.setTime(time);

        inner_.interpolate(xold, f, x);

        // If there are Dirichlet-constrained DOFs on processor boundaries, we need to make
        // those DOFs consistent again
        if constexpr (requires(const typename Traits::TrialGridFunctionSpace& gfs) { gfs.entitySet(); })
        {
          auto es = trialGridFunctionSpace().entitySet();
          if (es.comm().size() > 1)
          {
            CopyDataHandle<typename Traits::TrialGridFunctionSpace,X> data_handle(trialGridFunctionSpace(),x);
            es.communicate(data_handle,InteriorBorder_All_Interface,ForwardCommunication);
          }
        }
      }

      //! Set the time stepping method
      void setMethod(const OneStepParameters& method)
      {
        method_ = &method;
      }

      //! Parametrize the assembler with a time stepping method, the time and the time step size
      void preStep(const OneStepParameters& method, Real time, Real dt)
      {
        method_ = &method;
        time_ = time;
        dt_ = dt;

        // This switch decides which term will be multiplied with dt
        switch (dtMode_)
        {
        case DivideOperator1ByDT:
          dtFactor0_ = 1.0;
          dtFactor1_ = 1.0 / dt;
          break;
        case MultiplyOperator0ByDT:
          dtFactor0_ = dt;
          dtFactor1_ = 1.0;
          break;
        case DoNotAssembleDT:
          dtFactor0_ = 1.0;
          dtFactor1_ = 1.0;
          break;
        }

        lop0_.preStep(time, dt, method.s());
        lop1_.preStep(time, dt, method.s());
      }

      //! To be called after the step is completed
      void postStep()
      {
        lop0_.postStep();
        lop1_.postStep();
      }

      //! To be called after the stage is completed
      void postStage()
      {
        lop0_.postStage();
        lop1_.postStage();
      }

      //! Time step size suggested by the local operators
      Real suggestTimestep(Real dt) const
      {
        Real suggested_dt = std::min(lop0_.suggestTimestep(dt), lop1_.suggestTimestep(dt));
        if (trialGridFunctionSpace().gridView().comm().size() > 1)
          suggested_dt = trialGridFunctionSpace().gridView().comm().min(suggested_dt);
        return suggested_dt;
      }

      void update()
      {
        go0_.update();
        go1_.update();
        inner_.update();
        constResidual_ = Range(testGridFunctionSpace());
      }

      void make_consistent(Jacobian& a) const
        requires requires(const InnerGridOperator& go, Jacobian& a) { go.make_consistent(a); }
      {
        inner_.make_consistent(a);
      }

      const typename Traits::MatrixBackend& matrixBackend() const
      {
        return inner_.matrixBackend();
      }

      const typename Traits::TrialGridFunctionSpaceConstraints& trialConstraints() const
      {
        return inner_.localAssembler().trialConstraints();
      }

    private:

      enum DTAssemblingMode { DivideOperator1ByDT, MultiplyOperator0ByDT, DoNotAssembleDT };

      // contribution of a previous stage: solution, time, weights of both local operators
      struct Contribution
      {
        const Domain* x;
        Real time;
        Real w0;
        Real w1;
      };

      // r += sum of the residuals of the contributions, in one grid traversal if
      // the inner grid operator supports it
      void accumulate(const std::vector<Contribution>& contributions, Range& r) const
      {
        std::vector<const Domain*> xs;
        for (const auto& contribution : contributions)
          xs.push_back(contribution.x);
        auto prepare = [&](std::size_t j) {
          setTime(contributions[j].time);
          setWeights(contributions[j].w0, contributions[j].w1);
        };
        if constexpr (requires { Impl::accumulateResiduals(inner_, xs, r, prepare); })
          Impl::accumulateResiduals(inner_, xs, r, prepare);
        else
          for (std::size_t j = 0; j < xs.size(); ++j)
          {
            prepare(j);
            inner_.residual(*xs[j], r);
          }
      }

      static bool nonzero(Real w)
      {
        using std::abs;
        return abs(w) > 1e-6;
      }

      const OneStepParameters& method() const
      {
        if (not method_)
          DUNE_THROW(Dune::InvalidStateException, "No time stepping method set, call preStep() or setMethod() first");
        return *method_;
      }

      Real timeAtStage(unsigned int stage) const
      {
        return time_ + method().d(stage) * dt_;
      }

      void setTime(Real time) const
      {
        lop0_.setTime(time);
        lop1_.setTime(time);
      }

      void setWeights(Real w0, Real w1) const
      {
        sum_.setWeight(w0, 0);
        sum_.setWeight(w1, 1);
      }

      // time and weights of the stage equation: b_rr dt r_0(x) + r_1(x)
      void setStageWeights() const
      {
        const Real b = method().b(stage_, stage_);
        setTime(timeAtStage(stage_));
        setWeights(nonzero(b) ? b * dtFactor0_ : 0.0, dtFactor1_);
      }

      // implicit: union of both patterns, explicit: only the temporal part
      void setPatternWeights() const
      {
        setWeights(implicit ? 1.0 : 0.0, 1.0);
      }

      GO0& go0_;
      GO1& go1_;
      LOP0& lop0_;
      LOP1& lop1_;
      mutable SumLocalOperator sum_;
      InnerGridOperator inner_;
      Range constResidual_;

      const OneStepParameters* method_ = nullptr;
      Real time_ = 0.0;
      Real dt_ = 0.0;
      unsigned int stage_ = 0;
      DTAssemblingMode dtMode_;
      Real dtFactor0_ = 1.0;
      Real dtFactor1_ = 1.0;
    };
  }
}
#endif // DUNE_PDELAB_GRIDOPERATOR_ONESTEP_HH
