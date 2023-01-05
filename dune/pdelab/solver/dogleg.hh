// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_SOLVER_DOGLEG_HH
#define DUNE_PDELAB_SOLVER_DOGLEG_HH

#include<dune/pdelab/solver/newton.hh>

namespace Dune::PDELab
{
  /** \brief Newton dogleg Cauchy trust region nonlinear solver
   *
   * Trust region method. Takes Newton step if it is inside trust region,
   * damped Cauchy (steppest descent) step when Newton is outside and
   * Cauchy too, and dogleg step when Newton is outside and Cauchy inside
   * the trust region. Dogleg point is on the line from Cauchy to Newton
   * point and its distance from starting point is radius of trust region.
   *
   * Trust region radius is calculated from estimate how good is
   * quadratic approximation of the problem. The Hessian is approximated
   * by J^TJ, J is Jacobian.
   *
   * Almost completely implemented in Newton's LineSearch step, but
   * setting initial trust region must be made outside the iteration.
   */

  template <typename GridOperator_, typename LinearSolver_>
  class Dogleg
  {
  public:
    //! type of underlying Newton solver
    using Newton = NewtonMethod<GridOperator_,LinearSolver_>;

    //! Type of the grid operator
    using GridOperator = GridOperator_;

    //! Type of the linear solver
    using LinearSolver = LinearSolver_;

    //! Type of the domain (solution)
    using Domain = typename GridOperator::Traits::Domain;

    //! Type of the range (residual)
    using Range = typename GridOperator::Traits::Range;

    //! Type of the Jacobian matrix
    using Jacobian = typename GridOperator::Traits::Jacobian;

    //! Number type
    using Real = typename Dune::FieldTraits<typename Domain::ElementType>::real_type;

    //! Type of results
    using Result = PDESolverResult<Real>;

    //! Type of line search interface
    using LineSearch = LineSearchInterface<Domain>;

    //! Return results
    inline Result& result()
    {
      return _newton.result();
    }

    inline virtual void prepareStep(Domain& solution)
    {
      _newton.prepareStep(solution);
    }

    inline virtual void linearSolve()
    {
      _newton.linearSolve();
    }

    //! Solve the nonlinear problem using solution as initial guess and for storing the result
    inline virtual void apply(Domain& solution)
    {
      _lineSearch.setRadius(solution);
      _newton.apply(solution);
    }

    //! Update _residual and defect in _result
    inline virtual void updateDefect(Domain& solution)
    {
      _newton.updateDefect(solution);
    }

    //! Set how much output you get
    inline void setVerbosityLevel(unsigned int verbosity)
    {
      _newton.setVerbosityLevel(verbosity);
    }

    //! Get verbosity level
    inline unsigned int getVerbosityLevel() const
    {
      return _newton.getVerbosityLevel();
    }

    //! Set reduction Newton needs to achieve
    inline void setReduction(Real reduction)
    {
      _newton.setReduction(reduction);
    }

    //! Get reduction
    inline Real getReduction() const
    {
      return _newton.getReduction;
    }

    //! Set absolute convergence limit
    inline void setAbsoluteLimit(Real absoluteLimit)
    {
      _newton.setAbsoluteLimit(absoluteLimit);
    }

    inline Real getAbsoluteLimit() const
    {
      return _newton.getAbsoluteLimit();
    }

    //! Set whether the jacobian matrix should be kept across calls to apply().
    inline void setKeepMatrix(bool b)
    {
      _newton,setKeepMatrix(b);
    }

    //! Set whether to use the maximum norm for stopping criteria.
    inline void setUseMaxNorm(bool b)
    {
      _newton.setUseMaxNorm(b);
    }

    //! Does the problem have hanging nodes
    inline void setHangingNodeModifications(bool b)
    {
      _newton.setHangingNodeModifications(b);
    }


    //! Return whether the jacobian matrix is kept across calls to apply().
    inline bool keepMatrix() const
    {
      return _newton.keepMatrix();
    }

    //! Discard the stored Jacobian matrix.
    inline void discardMatrix()
    {
      _newton.discardMatrix();
    }

    //! Access stored Jacobian matrix.
    inline std::shared_ptr<Jacobian> getJacobian ()
    {
      return _newton.getJacobian();
    }

    //! Access stored residual
    inline const Range& getResidual () const
    {
      return _newton.getResidual();
    }

    /**\brief Set the minimal reduction in the linear solver
     *
     * \note with minLinearReduction > 0, the linear reduction will be
     * determined as mininum of the minLinearReduction and the linear reduction
     * needed to achieve second order Newton convergence. (As long as you are
     * not using a fixed linear reduction)
     */
    inline void setMinLinearReduction(Real minLinearReduction)
    {
      _newton.setMinLinearReduction(minLinearReduction);
    }

    /** \brief Set wether to use a fixed reduction in the linear solver
     *
     * \note If fixedLinearReduction is true, the linear reduction rate will
     *  always be fixed to minLinearReduction.
     */
    inline void setFixedLinearReduction(bool fixedLinearReduction)
    {
      _newton.setFixedLinearReduction(fixedLinearReduction);
    }

    /** \brief Set a threshold, when the linear operator is reassembled
     *
     * We allow to keep the linear operator over several newton iterations. If
     * the reduction in the newton drops below a given threshold the linear
     * operator is reassembled to ensure convergence.
     */
    inline void setReassembleThreshold(Real reassembleThreshold)
    {
      _newton.setReassembleThreshold(reassembleThreshold);
    }

    /** \brief Interpret a parameter tree as a set of options for the newton solver
     *
     *  Possible parameters:
     *
     *  example configuration:
     *
     *  \code
     *  [newton_parameters]
     *  ReassembleThreshold = 0.1
     *  AbsoluteLimit = 1e-6
     *  Reduction = 1e-4
     *  MinLinearReduction = 1e-3
     *  MaxIterations = 15
     *  LineSearchDampingFactor = 0.7
     *  \endcode
     *
     *  and invocation in the code:
     *  \code
     *  newton.setParameters(param.sub("NewtonParameters"));
     *  \endcode
     *
     *  This can also be used to set single parameters like this
     *
     *  \code
     *  Dune::ParameterTree ptree;
     *  ptree["verbosity"] = "4";
     *  newton.setParameters(ptree);
     *  \endcode
     */
    inline void setParameters(const ParameterTree& parameterTree){
      _newton.setParameters(parameterTree);
    }

    //! Set the termination criterion
    inline void setTerminate(std::shared_ptr<TerminateInterface> terminate)
    {
      _newton.setTerminate(terminate);
    }

    //! Return a pointer to the stored termination criterion
    inline std::shared_ptr<TerminateInterface> getTerminate() const
    {
      return _newton.getTerminate();
    }

    /**\brief Set the line search method is disabled
     *
     * Dogleg uses DoglegLineSearch method.
     */

    //! Return a pointer to the stored line search
    inline std::shared_ptr<LineSearch> getLineSearch() const
    {
      return _newton.getLineSearch();
    }

    /**\brief Output NewtonMethod parameters
     *
     * Setting parameters using ParameterTree is quite error prone.
     * Checking parameters without setting up the debugger can be useful.
     */
    inline void printParameters(const std::string& _name = "NewtonMethod") const
    {
      _newton.printParameters(_name);
    }

    //! Construct Newton using default parameters with default parameters
    /**
       in p
     */
    Dogleg(
      const GridOperator& gridOperator,
      LinearSolver& linearSolver)
      : _newton(gridOperator,linearSolver)
    {
      _lineSearch = std::make_shared<Newton>(_newton);
      _newton.setLineSearch(_lineSearch);
    }

    //! Construct Newton passing a parameter tree
    Dogleg(
      const GridOperator& gridOperator,
      LinearSolver& linearSolver,
      const ParameterTree& parameterTree)
      : _newton(gridOperator,linearSolver)
    {
      _newton.setParameters(parameterTree);
      _lineSearch = std::make_shared<Newton>(_newton);
      _newton.setLineSearch(_lineSearch);
    }

    virtual ~Dogleg() {}

  private:
    Newton _newton;
    std::shared_ptr<DoglegLineSearch<Newton>> _lineSearch;
  };

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_SOLVER_DOGLEG_HH