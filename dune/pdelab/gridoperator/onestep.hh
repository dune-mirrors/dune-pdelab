// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDOPERATOR_ONESTEP_HH
#define DUNE_PDELAB_GRIDOPERATOR_ONESTEP_HH

#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/gridoperator/common/gridoperatorutilities.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep/localassembler.hh>
#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/instationary/onestepparameter.hh>
#include <dune/pdelab/localoperator/weightedsum.hh>

namespace Dune{
  namespace PDELab{

    template<typename GO0, typename GO1, bool implicit = true>
    class OneStepGridOperator
    {
    public:

      //! The sparsity pattern container for the jacobian matrix
      typedef typename GO0::Pattern Pattern;

      //! The global UDG assembler type
      typedef typename GO0::Traits::Assembler Assembler;

      //! The local assembler types of the subordinate grid operators
      //! @{
      typedef typename GO0::Traits::LocalAssembler LocalAssemblerDT0;
      typedef typename GO1::Traits::LocalAssembler LocalAssemblerDT1;
      //! @}

      //! The local assembler type
      typedef OneStepLocalAssembler<OneStepGridOperator,LocalAssemblerDT0,LocalAssemblerDT1> LocalAssembler;

      //! The BorderDOFExchanger
      typedef typename GO0::BorderDOFExchanger BorderDOFExchanger;

      //! The grid operator traits
      typedef Dune::PDELab::GridOperatorTraits
      <typename GO0::Traits::TrialGridFunctionSpace,
       typename GO0::Traits::TestGridFunctionSpace,
       typename GO0::Traits::MatrixBackend,
       typename GO0::Traits::DomainField,
       typename GO0::Traits::RangeField,
       typename GO0::Traits::JacobianField,
       typename GO0::Traits::TrialGridFunctionSpaceConstraints,
       typename GO0::Traits::TestGridFunctionSpaceConstraints,
       Assembler,
       LocalAssembler> Traits;

      //! The io types of the operator
      //! @{
      typedef typename Traits::Domain Domain;
      typedef typename Traits::Range Range;
      typedef typename Traits::Jacobian Jacobian;
      //! @}

      template <typename MFT>
      struct MatrixContainer
      {
        typedef Jacobian Type;
      };

      //! The type for real number e.g. time
      typedef typename LocalAssembler::Real Real;

      //! The type of the one step method parameters
      typedef typename LocalAssembler::OneStepParameters OneStepParameters;

      //! Constructor for non trivial constraints
      OneStepGridOperator(GO0 & go0_, GO1 & go1_)
        : global_assembler(go0_.assembler()),
          go0(go0_), go1(go1_),
          la0(go0_.localAssembler()), la1(go1_.localAssembler()),
          const_residual( go0_.testGridFunctionSpace() ),
          local_assembler(la0,la1, const_residual)
      {
        GO0::setupGridOperators(std::tie(go0_,go1_));
        if(not implicit)
          local_assembler.setDTAssemblingMode(LocalAssembler::DoNotAssembleDT);
      }

      //! Determines whether the time step size is multiplied to the
      //! mass term (first order time derivative) or the elliptic term
      //! (zero-th order time derivative).
      void divideMassTermByDeltaT()
      {
        if(not implicit)
          DUNE_THROW(Dune::Exception,"This function should not be called in explicit mode");
        local_assembler.setDTAssemblingMode(LocalAssembler::DivideOperator1ByDT);
      }
      void multiplySpatialTermByDeltaT()
      {
        if(not implicit)
          DUNE_THROW(Dune::Exception,"This function should not be called in explicit mode");
        local_assembler.setDTAssemblingMode(LocalAssembler::MultiplyOperator0ByDT);
      }

      //! Get the trial grid function space
      const typename Traits::TrialGridFunctionSpace& trialGridFunctionSpace() const
      {
        return global_assembler.trialGridFunctionSpace();
      }

      //! Get the test grid function space
      const typename Traits::TestGridFunctionSpace& testGridFunctionSpace() const
      {
        return global_assembler.testGridFunctionSpace();
      }

      //! Get dimension of space u
      typename Traits::TrialGridFunctionSpace::Traits::SizeType globalSizeU () const
      {
        return trialGridFunctionSpace().globalSize();
      }

      //! Get dimension of space v
      typename Traits::TestGridFunctionSpace::Traits::SizeType globalSizeV () const
      {
        return testGridFunctionSpace().globalSize();
      }

      Assembler & assembler() const { return global_assembler; }

      LocalAssembler & localAssembler() const { return local_assembler; }

      //! Fill pattern of jacobian matrix
      void fill_pattern(Pattern & p) const
      {
        if(implicit)
        {
          typedef typename LocalAssembler::LocalPatternAssemblerEngine PatternEngine;
          PatternEngine & pattern_engine = local_assembler.localPatternAssemblerEngine(p);
          global_assembler.assemble(pattern_engine);
        }
        else
        {
          typedef typename LocalAssembler::LocalExplicitPatternAssemblerEngine PatternEngine;
          PatternEngine & pattern_engine = local_assembler.localExplicitPatternAssemblerEngine(p);
          global_assembler.assemble(pattern_engine);
        }
      }

      //! Assemble constant part of residual
      void preStage(unsigned int stage, const std::vector<Domain*> & x)
      {
        if(not implicit)
          DUNE_THROW(Dune::Exception,"This function should not be called in explicit mode");

        typedef typename LocalAssembler::LocalPreStageAssemblerEngine PreStageEngine;
        local_assembler.setStage(stage);
        PreStageEngine & prestage_engine = local_assembler.localPreStageAssemblerEngine(x);
        global_assembler.assemble(prestage_engine);
      }

      //! Assemble residual
      void residual(const Domain & x, Range & r) const
      {
        if(not implicit)
          DUNE_THROW(Dune::Exception,"This function should not be called in explicit mode");

        typedef typename LocalAssembler::LocalResidualAssemblerEngine ResidualEngine;
        ResidualEngine & residual_engine = local_assembler.localResidualAssemblerEngine(r,x);
        global_assembler.assemble(residual_engine);
      }

      //! Assemble jacobian
      void jacobian(const Domain & x, Jacobian & a) const
      {
        if(not implicit)
          DUNE_THROW(Dune::Exception,"This function should not be called in explicit mode");

        typedef typename LocalAssembler::LocalJacobianAssemblerEngine JacobianEngine;
        JacobianEngine & jacobian_engine = local_assembler.localJacobianAssemblerEngine(a,x);
        global_assembler.assemble(jacobian_engine);
      }

      //! Assemble jacobian and residual simultaneously for explicit treatment
      void explicit_jacobian_residual(unsigned int stage, const std::vector<Domain*> & x,
                                      Jacobian & a, Range & r1, Range & r0)
      {
        if(implicit)
          DUNE_THROW(Dune::Exception,"This function should not be called in implicit mode");

        local_assembler.setStage(stage);

        typedef typename LocalAssembler::LocalExplicitJacobianResidualAssemblerEngine
          ExplicitJacobianResidualEngine;

        ExplicitJacobianResidualEngine & jacobian_residual_engine
          = local_assembler.localExplicitJacobianResidualAssemblerEngine(a,r0,r1,x);
        global_assembler.assemble(jacobian_residual_engine);
      }

      //! Apply jacobian matrix to the vector update without explicitly assembling it
      void jacobian_apply(const Domain & update, Range & result) const
      {
        if ((not la0.localOperator().isLinear) or (not la1.localOperator().isLinear))
          DUNE_THROW(Dune::Exception, "Your trying to use a linear jacobian apply for a non linear problem.");
        global_assembler.assemble(local_assembler.localJacobianApplyAssemblerEngine(update, result));
      }

      //! Apply jacobian matrix to the vector update without explicitly assembling it
      void jacobian_apply(const Domain & solution, const Domain & update, Range & result) const
      {
        if (la0.localOperator().isLinear and la1.localOperator().isLinear)
          DUNE_THROW(Dune::Exception, "Your trying to use a non linear jacobian apply for a linear problem.");
        global_assembler.assemble(local_assembler.localJacobianApplyAssemblerEngine(solution, update, result));
      }

      //! Interpolate constrained values from given function f
      template<typename F, typename X>
      void interpolate (unsigned stage, const X& xold, F& f, X& x) const
      {
        // Set time in boundary value function
        f.setTime(local_assembler.timeAtStage(stage));

        go0.localAssembler().setTime(local_assembler.timeAtStage(stage));

        // Interpolate
        go0.interpolate(xold,f,x);

        // Copy non-constrained dofs from old time step
        Dune::PDELab::copy_nonconstrained_dofs(local_assembler.trialConstraints(),xold,x);

        // If there are Dirichlet-constrained DOFs on processor boundaries, we need to make
        // those DOFs consistent again
        auto es = trialGridFunctionSpace().entitySet();
        if (es.comm().size() > 1)
        {
          CopyDataHandle<typename Traits::TrialGridFunctionSpace,X> data_handle(trialGridFunctionSpace(),x);
          es.communicate(data_handle,InteriorBorder_All_Interface,ForwardCommunication);
        }
      }

      //! set time stepping method
      void setMethod (const TimeSteppingParameterInterface<Real>& method_)
      {
        local_assembler.setMethod(method_);
      }

      //! parametrize assembler with a time-stepping method
      void preStep (const TimeSteppingParameterInterface<Real>& method_, Real time_, Real dt_)
      {
        local_assembler.setMethod(method_);
        local_assembler.preStep(time_,dt_,method_.s());
      }

      //! to be called after step is completed
      void postStep ()
      {
        la0.postStep();
        la1.postStep();
      }

      //! to be called after stage is completed
      void postStage ()
      {
        la0.postStage();
        la1.postStage();
      }

      //! to be called once before each stage
      Real suggestTimestep (Real dt) const
      {
        Real suggested_dt = std::min(la0.suggestTimestep(dt),la1.suggestTimestep(dt));
        if (trialGridFunctionSpace().gridView().comm().size()>1)
          suggested_dt =  trialGridFunctionSpace().gridView().comm().min(suggested_dt);
        return suggested_dt;
      }

      void update()
      {
        go0.update();
        go1.update();
        const_residual = Range(go0.testGridFunctionSpace());
      }

      void make_consistent(Jacobian& a) const
      {
        go0.make_consistent(a);
      }

      const typename Traits::MatrixBackend& matrixBackend() const
      {
        return go0.matrixBackend();
      }

      const typename LocalAssembler::Traits::TrialGridFunctionSpaceConstraints trialConstraints() const
      {
        return local_assembler.trialConstraints();
      }

    private:
      Assembler & global_assembler;
      GO0 & go0;
      GO1 & go1;
      LocalAssemblerDT0 & la0;
      LocalAssemblerDT1 & la1;
      Range const_residual;
      mutable LocalAssembler local_assembler;
    };


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
      class GenericOneStepGridOperator
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
        using InnerGridOperator = typename ReplaceLocalOperator<GO0, SumLocalOperator>::type;

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
        GenericOneStepGridOperator(GO0& go0, GO1& go1)
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
        GenericOneStepGridOperator(const GenericOneStepGridOperator&) = delete;
        GenericOneStepGridOperator& operator=(const GenericOneStepGridOperator&) = delete;

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

          constResidual_ = 0.0;
          // TODO: one grid traversal per previous stage; a local assembler for
          // several coefficient vectors would need only one traversal.
          for (unsigned int s = 0; s < stage; ++s)
          {
            const Real b = method().b(stage, s);
            const Real a = method().a(stage, s);
            if (not nonzero(b) and not nonzero(a))
              continue;
            setTime(timeAtStage(s));
            setWeights(nonzero(b) ? b * dtFactor0_ : 0.0, nonzero(a) ? a * dtFactor1_ : 0.0);
            inner_.residual(*x[s], constResidual_);
          }
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

          r0 = 0.0;
          r1 = 0.0;
          // TODO: up to 2r+1 grid traversals; with a system assembly for a
          // matrix and two vectors this would need a single traversal.
          for (unsigned int s = 0; s < stage; ++s)
          {
            const Real b = method().b(stage, s);
            const Real a_rs = method().a(stage, s);
            setTime(timeAtStage(s));
            if (nonzero(b))
            {
              setWeights(b * dtFactor0_, 0.0);
              inner_.residual(*x[s], r0);
            }
            if (nonzero(a_rs))
            {
              setWeights(0.0, a_rs * dtFactor1_);
              inner_.residual(*x[s], r1);
            }
          }

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

    } // namespace Impl
  }
}
#endif // DUNE_PDELAB_GRIDOPERATOR_ONESTEP_HH
