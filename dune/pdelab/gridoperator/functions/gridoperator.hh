// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDOPERATOR_FUNCTIONS_GRIDOPERATOR_HH
#define DUNE_PDELAB_GRIDOPERATOR_FUNCTIONS_GRIDOPERATOR_HH

#include <memory>
#include <optional>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>

#include <dune/assembler/backends/istlmatrixbackend.hh>
#include <dune/assembler/backends/istlvectorbackend.hh>
#include <dune/assembler/defaultglobalassembler.hh>

#include <dune/pdelab/backend/common/spacetraits.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintstransformation.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridoperator/common/gridoperatorutilities.hh>
#include <dune/pdelab/gridoperator/functions/localoperatorassembler.hh>
#include <dune/pdelab/gridoperator/legacy/gridoperator.hh>

/** \file
 * \brief GridOperator for PDELab local operators on dune-functions bases.
 *
 * Specialization of Dune::PDELab::GridOperator for trial and test spaces
 * given by dune-functions global bases. It offers the interface of the
 * GridOperator on grid function spaces, uses the same PDELab vector and
 * matrix wrappers (Backend::Vector<Basis,E>, Backend::Matrix<...>) and the
 * same constraints containers, so the PDELab solvers work unchanged. The
 * assembly is done by dune-assembler with the local assemblers from
 * functions/localoperatorassembler.hh.
 *
 * \warning This is experimental and sequential only.
 */

namespace Dune::PDELab {

  //! Constraints container for Dirichlet-type constraints on a dune-functions basis
  template<class Basis, class F = double>
  using BasisConstraintsContainer = ConstraintsTransformation<typename Basis::MultiIndex, typename Basis::MultiIndex, F>;

  /** \brief Mark all DOFs of a basis associated with the boundary as Dirichlet constrained
   *
   * \param basis        a global basis or a subspace basis (e.g. the velocity of a Taylor-Hood basis)
   * \param constraints  constraints container indexed by the multi-indices of the root basis
   */
  template<class Basis, class CC>
  void addBoundaryDirichletConstraints(const Basis& basis, CC& constraints)
  {
    Dune::Functions::forEachBoundaryDOF(basis, [&](const auto& index) {
      constraints[index];
    });
  }

  /** \brief Mark the DOFs on the Dirichlet part of the boundary as constrained
   *
   * \param basis        a global basis or a subspace basis
   * \param constraints  constraints container indexed by the multi-indices of the root basis
   * \param bc           PDELab boundary condition, e.g. from makeBoundaryConditionFromCallable()
   *                     or ConvectionDiffusionBoundaryConditionAdapter
   *
   * As for ConformingDirichletConstraints, bc.isDirichlet(intersection, x) is evaluated at
   * the center x of each boundary intersection, and all DOFs of Dirichlet intersections
   * are constrained.
   */
  template<class Basis, class CC, class BC>
  void addBoundaryDirichletConstraints(const Basis& basis, CC& constraints, const BC& bc)
  {
    using Intersection = typename Basis::GridView::Intersection;
    Dune::Functions::forEachBoundaryDOF(basis, [&](auto&& localIndex, const auto& localView, const Intersection& intersection) {
      IntersectionGeometry<Intersection> ig(intersection,
        Impl::Functions::intersectionIndex(basis.gridView(), intersection));
      const auto center = referenceElement(intersection.geometry()).position(0,0);
      if (bc.isDirichlet(ig, center))
        constraints[localView.index(localIndex)];
    });
  }

  namespace Impl::Functions {

    /** \brief Local assembler of the grid operator
     *
     * The assembly itself is done by dune-assembler. As for grid function
     * spaces, this class provides the constraints and their application
     * (derived from LocalAssemblerBase) and access to the local operator.
     */
    template<typename MB, typename CU, typename CV, typename LOP>
    class ConstraintsAssembler
      : public LocalAssemblerBase<MB,CU,CV>
    {
      using Base = LocalAssemblerBase<MB,CU,CV>;

    public:
      using LocalOperator = LOP;

      explicit ConstraintsAssembler(LOP& lop)
        : Base(), lop_(&lop)
      {}

      ConstraintsAssembler(LOP& lop, const CU& cu, const CV& cv)
        : Base(cu, cv), lop_(&lop)
      {}

      using Base::handle_dirichlet_constraints;

      LOP& localOperator() const { return *lop_; }

      static constexpr bool isLinear() { return LOP::isLinear; }

    private:
      LOP* lop_;
    };

  } // namespace Impl::Functions

  /** \brief GridOperator for trial and test spaces given by dune-functions bases
   *
   * Template parameters as for the GridOperator on grid function spaces:
   * \tparam GFSU trial space: dune-functions global basis
   * \tparam GFSV test space: dune-functions global basis
   * \tparam LOP  PDELab local operator
   * \tparam MB   matrix backend, e.g. ISTL::BCRSMatrixBackend<>
   * \tparam DF, RF, JF field types of domain, range and jacobian
   * \tparam CU, CV constraints containers (BasisConstraintsContainer or EmptyTransformation)
   *
   * Not supported (yet): parallel assembly and local
   * operators using doSkeletonTwoSided or *_volume_post_skeleton.
   */
  template<typename GFSU, typename GFSV, typename LOP,
           typename MB, typename DF, typename RF, typename JF,
           typename CU, typename CV>
    requires (Impl::DuneFunctionsBasis<GFSU> and Impl::DuneFunctionsBasis<GFSV>)
  class GridOperator<GFSU,GFSV,LOP,MB,DF,RF,JF,CU,CV>
  {
  public:

    //! The global assembler type (rows: test space, columns: trial space)
    using Assembler = Dune::Assembler::Assembler<GFSV,GFSU>;

    //! Provides the constraints and their application, as in the GridOperator on grid function spaces
    using LocalAssembler = Impl::Functions::ConstraintsAssembler<MB,CU,CV,LOP>;

    using Domain = Dune::PDELab::Backend::Vector<GFSU,DF>;
    using Range = Dune::PDELab::Backend::Vector<GFSV,RF>;
    using Jacobian = Dune::PDELab::Backend::Matrix<MB,Domain,Range,JF>;

    using Traits = Dune::PDELab::GridOperatorTraits<GFSU,GFSV,MB,DF,RF,JF,CU,CV,Assembler,LocalAssembler>;

    template <typename MFT>
    struct MatrixContainer
    {
      typedef typename Traits::Jacobian Type;
    };

    GridOperator(const GFSU& gfsu, const CU& cu, const GFSV& gfsv, const CV& cv, LOP& lop, const MB& mb = MB())
      : gfsu_(&gfsu), gfsv_(&gfsv), lop_(&lop), local_assembler_(lop, cu, cv), backend_(mb)
    {}

    GridOperator(const GFSU& gfsu, const GFSV& gfsv, LOP& lop, const MB& mb = MB())
      : gfsu_(&gfsu), gfsv_(&gfsv), lop_(&lop), local_assembler_(lop), backend_(mb)
    {}

    const GFSU& trialGridFunctionSpace() const { return *gfsu_; }
    const GFSV& testGridFunctionSpace() const { return *gfsv_; }

    std::size_t globalSizeU() const { return gfsu_->dimension(); }
    std::size_t globalSizeV() const { return gfsv_->dimension(); }

    LocalAssembler& localAssembler() const { return local_assembler_; }

    const typename Traits::MatrixBackend& matrixBackend() const { return backend_; }

    //! add the residual r(x) to r and apply the test space constraints
    void residual(const Domain& x, Range& r) const
    {
      auto rBackend = Dune::Assembler::ISTLVectorBackend(Backend::native(r));
      auto localAssembler = Impl::Functions::makeResidualAssembler(*lop_, coefficients(x), *gfsv_, *gfsu_);
      Assembler(*gfsv_, *gfsu_).assembleVectorEntries(localAssembler, rBackend);
      constrain_residual(local_assembler_.testConstraints(), r);
    }

    //! add the jacobian at x to a and replace the constrained rows by unit rows
    void jacobian(const Domain& x, Jacobian& a) const
    {
      auto aBackend = Dune::Assembler::ISTLMatrixBackend(Backend::native(a));
      auto localAssembler = Impl::Functions::makeJacobianAssembler(*lop_, coefficients(x), *gfsv_, *gfsu_);
      Assembler(*gfsv_, *gfsu_).assembleMatrixEntries(localAssembler, aBackend);
      local_assembler_.handle_dirichlet_constraints(*gfsv_, a);
    }

    //! apply the (constant) jacobian of a linear local operator to update and add it to result
    void jacobian_apply(const Domain& update, Range& result) const
    {
      if (not lop_->isLinear)
        DUNE_THROW(Dune::Exception, "Your trying to use a non linear jacobian apply for a linear problem.");
      jacobianApply(update, update, result);
    }

    //! apply the jacobian at solution to update and add it to result
    void jacobian_apply(const Domain& solution, const Domain& update, Range& result) const
    {
      if (lop_->isLinear)
        DUNE_THROW(Dune::Exception, "Your trying to use a non linear jacobian apply for a linear problem.");
      jacobianApply(solution, update, result);
    }

    //! Interpolate the constrained DOFs from f, copy the others from xold (xold and x may be the same vector)
    template<typename F, typename X>
    void interpolate(const X& xold, F& f, X& x) const
    {
      std::optional<X> tmp;
      const X& source = (&xold == &x) ? tmp.emplace(xold) : xold;
      Dune::PDELab::interpolate(f, *gfsu_, x);
      Dune::PDELab::copy_nonconstrained_dofs(local_assembler_.trialConstraints(), source, x);
    }

    //! Only sequential assembly is supported, hence there is nothing to update.
    void update() {}

    /** \brief Set up the sparsity pattern of a native matrix, all entries zero
     *
     * Called by ISTL::BCRSMatrixBackend::buildPattern(), i.e. when constructing a Jacobian from the grid operator.
     */
    template<class NativeMatrix>
    void setupMatrixPattern(NativeMatrix& matrix) const
    {
      auto backend = Dune::Assembler::ISTLMatrixBackend(matrix);
      auto patternBuilder = backend.patternBuilder();
      patternBuilder.resize(*gfsv_, *gfsu_);
      // the coefficients are not accessed while assembling the pattern
      Backend::Native<Domain> unused;
      auto localAssembler = Impl::Functions::makeJacobianAssembler(*lop_, Dune::Functions::istlVectorBackend(std::as_const(unused)), *gfsv_, *gfsu_);
      Assembler(*gfsv_, *gfsu_).assembleMatrixPattern(localAssembler, patternBuilder);
      patternBuilder.setupMatrix();
      backend.setZero();
    }

  private:

    // lightweight accessor x[multiIndex] for the local assemblers
    static auto coefficients(const Domain& x)
    {
      return Dune::Functions::istlVectorBackend(Backend::native(x));
    }

    void jacobianApply(const Domain& solution, const Domain& update, Range& result) const
    {
      auto yBackend = Dune::Assembler::ISTLVectorBackend(Backend::native(result));
      using Coefficients = decltype(coefficients(solution));
      auto localAssembler = Impl::Functions::JacobianApplyLocalAssembler<LOP,Coefficients,Coefficients,GFSV,GFSU>(
        *lop_, coefficients(solution), coefficients(update), *gfsu_);
      Assembler(*gfsv_, *gfsu_).assembleVectorEntries(localAssembler, yBackend);
      constrain_residual(local_assembler_.testConstraints(), result);
    }

    const GFSU* gfsu_;
    const GFSV* gfsv_;
    LOP* lop_;
    mutable LocalAssembler local_assembler_;
    MB backend_;
  };

  namespace Impl {

    /** \brief Add the residuals of several coefficient vectors to r in a single grid traversal
     *
     * GridOperator on dune-functions bases, see gridoperator/common/multiresidual.hh.
     * prepare(j) is called before the contributions of xs[j] are evaluated.
     */
    template<typename GFSU, typename GFSV, typename LOP,
             typename MB, typename DF, typename RF, typename JF,
             typename CU, typename CV, typename Prepare>
      requires (Impl::DuneFunctionsBasis<GFSU> and Impl::DuneFunctionsBasis<GFSV>)
    void accumulateResiduals(const GridOperator<GFSU,GFSV,LOP,MB,DF,RF,JF,CU,CV>& go,
                             const std::vector<const typename GridOperator<GFSU,GFSV,LOP,MB,DF,RF,JF,CU,CV>::Traits::Domain*>& xs,
                             typename GridOperator<GFSU,GFSV,LOP,MB,DF,RF,JF,CU,CV>::Traits::Range& r,
                             Prepare& prepare)
    {
      if (xs.empty())
        return;
      const auto& gfsu = go.trialGridFunctionSpace();
      const auto& gfsv = go.testGridFunctionSpace();
      const auto& lop = go.localAssembler().localOperator();
      using X = decltype(Dune::Functions::istlVectorBackend(Backend::native(*xs.front())));
      using Single = Impl::Functions::ResidualLocalAssembler<LOP, X, GFSV, GFSU>;
      std::vector<Single> copies;
      copies.reserve(xs.size());
      for (const auto* x : xs)
        copies.push_back(Impl::Functions::makeResidualAssembler(lop, Dune::Functions::istlVectorBackend(Backend::native(*x)), gfsv, gfsu));
      Impl::Functions::MultiResidualLocalAssembler<Single, Prepare> localAssembler(std::move(copies), prepare);
      auto rBackend = Dune::Assembler::ISTLVectorBackend(Backend::native(r));
      Dune::Assembler::Assembler<GFSV,GFSU>(gfsv, gfsu).assembleVectorEntries(localAssembler, rBackend);
      constrain_residual(go.localAssembler().testConstraints(), r);
    }

  } // namespace Impl

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_GRIDOPERATOR_FUNCTIONS_GRIDOPERATOR_HH
