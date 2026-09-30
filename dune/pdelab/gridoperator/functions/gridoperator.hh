// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDOPERATOR_FUNCTIONS_GRIDOPERATOR_HH
#define DUNE_PDELAB_GRIDOPERATOR_FUNCTIONS_GRIDOPERATOR_HH

#include <memory>
#include <type_traits>

#include <dune/common/exceptions.hh>

#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>

#include <dune/assembler/backends/istlmatrixbackend.hh>
#include <dune/assembler/backends/istlvectorbackend.hh>
#include <dune/assembler/defaultglobalassembler.hh>

#include <dune/pdelab/backend/common/spacetraits.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintstransformation.hh>
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

  namespace Impl::Functions {

    //! Constraints of the grid operator; as for grid function spaces derived from LocalAssemblerBase
    template<typename MB, typename CU, typename CV>
    class ConstraintsAssembler
      : public LocalAssemblerBase<MB,CU,CV>
    {
      using Base = LocalAssemblerBase<MB,CU,CV>;

    public:
      using Base::Base;
      using Base::handle_dirichlet_constraints;
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
   * Not supported (yet): parallel assembly, interpolate(), and local
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
    using LocalAssembler = Impl::Functions::ConstraintsAssembler<MB,CU,CV>;

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
      : gfsu_(&gfsu), gfsv_(&gfsv), lop_(&lop), local_assembler_(cu, cv), backend_(mb)
    {}

    GridOperator(const GFSU& gfsu, const GFSV& gfsv, LOP& lop, const MB& mb = MB())
      : gfsu_(&gfsu), gfsv_(&gfsv), lop_(&lop), local_assembler_(), backend_(mb)
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

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_GRIDOPERATOR_FUNCTIONS_GRIDOPERATOR_HH
