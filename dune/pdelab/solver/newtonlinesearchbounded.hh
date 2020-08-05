// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_SOLVER_NEWTONLINESEARCH_BOUNDED_HH
#define DUNE_PDELAB_SOLVER_NEWTONLINESEARCH_BOUNDED_HH
#endif

/* \brief BoundedLineSearch
 *
 * This file contains the implementation of the BoundedLineSearchNone
 * and BoundedLineSearchHackbuschReusken. It is included in the middle
 * of newtonlinesearch.hh -because it inherits from LineSearchInterface
 * and must defined for LineSearchStrategy, and both are in that file.
 * For brevity, the implementation is moved to a separate file (this one)
 * and guarded by macro, which also adjusts LineSearchStrategy so it
 * contains BoundedLineSearches.
 *
 * To use this file, simply include it before newtonlinesearch.hh.
 * The "guard" macro on top will unlock parts of its code with
 * BoundedLineSearches and inludes the following code in the middle
 * between LineSearchInterface and LineSearchStrategy.
 */

#ifdef DUNE_PDELAB_SOLVER_NEWTONLINESEARCH_BOUNDED_ENABLED
// Note: it is assumed this file is already in namespace Dune::PDELab
namespace Impl
{
  // Simple struct for calculating the number of unknowns of the system.
  // This works for GFS which are not Powerspace, i.e. have no CHILDREN.
  // Simply putting this std::max(..) formula into code to calculate
  // blocksize does not work.
  template<typename Domain>
  struct BlockSizeExtractor
  {
    std::size_t value = std::max(std::size_t(1),Domain::GridFunctionSpace::CHILDREN);
  };

  /* \brief CorrectSolution
   *
   * 4 classes with single purpose: traverse vector and set specified
   * values into respective bounds. They differentiate based on the type
   * of vector. The correct type is guessed automatically at compile
   * time from GridFunctionSpace (GFS) traits using class SetCorrectSolution.
   * These 4 classes are:
   *   CorrectSolutionSingle: for simple GFS with one unknown.
   *   CorrectSolutionLexicographic: for PowerGFS with Lexicographic OrderingTag
   *   CorrectSolutionEntityBlockedFixedBlocking: for PowerGFS with
   *     EntityBlocked Ordering tag and "fixed" Blocking
   *   CorrectSolutionEntityBlockedNoneBlocking: for PowerGFS with
   *     EntityBlocked OrderingTag and "none" Blocking.
   *
   * Heavily relies on ParameterClass = BoundedLineSearchParametersInterface<Newton>
   * which carries all information about which unknowns are bounded and how.
   */
  template<typename ParameterClass>
  class CorrectSolutionEntityBlockedFixedBlocking
  {
    const ParameterClass& pc;
    using Domain = typename ParameterClass::Domain;
  public:
    CorrectSolutionEntityBlockedFixedBlocking(const ParameterClass& pc_)
      : pc(pc_)
    {}

    void operator() (Domain& solution)
    {
      using Dune::PDELab::Backend::native;
      for (auto& block : native(solution))
      {
        for (std::size_t i=0; i<pc.blockLower.size(); ++i)
          if (block[pc.blockLower[i]] < pc.boundLower[i])
            block[pc.blockLower[i]] = pc.boundLower[i];
        for (std::size_t i=0; i<pc.blockUpper.size(); ++i)
          if (block[pc.blockUpper[i]] > pc.boundUpper[i])
            block[pc.blockUpper[i]] = pc.boundUpper[i];
      }
    }
  };

  template<typename ParameterClass>
  class CorrectSolutionEntityBlockedNoneBlocking
  {
    const ParameterClass& pc;
    using Domain = typename ParameterClass::Domain;
  public:
    CorrectSolutionEntityBlockedNoneBlocking(const ParameterClass& pc_)
      : pc(pc_)
    {}

    void operator() (Domain& solution)
    {
      using Dune::PDELab::Backend::native;
      // EntityBlocked with Blocking "none" is arranged in blocks, but
      // blocks have size 1. Variable i stores the position in the block.
      std::size_t i=0;
      std::size_t low=0; // watches the next lower bound to be corrected
      std::size_t up=0; // watches the next upper bound
      for (auto& block : native(solution))
      {
        // Firstly check whether it is one of the bounded variables.
        if (i==pc.blockLower[low])
        {
          if (block[0]<pc.boundLower[low])
            block[0] = pc.boundLower[low];
          low = (low+1) % pc.blockLower.size();
        }
        if (i==pc.blockUpper[up])
        {
          if (block[0]>pc.boundUpper[up])
            block[0] = pc.boundUpper[up];
          up = (up+1) % pc.blockUpper.size();
        }
        i = (i+1) % pc.blocksize;
      }
    }
  };

  template<typename ParameterClass>
  class CorrectSolutionLexicographic
  {
    const ParameterClass& pc;
    using Domain = typename ParameterClass::Domain;
  public:
    CorrectSolutionLexicographic(const ParameterClass& pc_)
      : pc(pc_)
    {}

    void operator() (Domain& solution)
    {
      using Dune::PDELab::Backend::native;

      std::size_t sol_size = native(solution).size()/pc.blocksize;
      for (std::size_t i=0; i<pc.blockLower.size(); ++i)
        for (std::size_t j=sol_size*pc.blockLower[i]; j<sol_size*(pc.blockLower[i]+1); ++j)
          if (native(solution)[j] < pc.boundLower[i])
            native(solution)[j] = pc.boundLower[i];
      for (std::size_t i=0; i<pc.blockUpper.size(); ++i)
        for (std::size_t j=sol_size*pc.blockUpper[i]; j<sol_size*(pc.blockUpper[i]+1); ++j)
          if (native(solution)[j] > pc.boundUpper[i])
            native(solution)[j] = pc.boundUpper[i];
    }
  };

  template<typename ParameterClass>
  class CorrectSolutionSingle
  {
    const ParameterClass& pc;
    using Domain = typename ParameterClass::Domain;
  public:
    CorrectSolutionSingle(const ParameterClass& pc_)
      : pc(pc_)
    {}

    void operator() (Domain& solution)
    {
      if (pc.blockUpper.size()>0)
        for (auto& v : solution)
          if (v>pc.boundUpper[0])
            v = pc.boundUpper[0];
      if (pc.blockLower.size()>0)
        for (auto& v : solution)
          if (v>pc.boundLower[0])
            v = pc.boundLower[0];
    }
  };

  /* \brief SetCorrectSolution::type
   *
   * This class outputs the correct type for CorrectSolution procedure
   * based on Domain type. 4 types are specialized for different ordering of vector:
   * Single, Lexicographic,and EntityBlocked with None or Fixed blocking.
   * Single is for simple spaces with one unknown, other are for PowerGridFunctionSpace.
   * PowerGFS are specialized based on the vector Ordering -Lexicographic or EntityBlocked.
   * EntityBlocked with Fixed blocking allows further optimization.
   *
   * Other options like bcrs blocking, InterleavedOrderingTag, or CompositeSpaces
   * are not tested/implemented.
   */
  template<typename Domain,typename ParameterClass>
  class SetCorrectSolution
  {
    using currentOrderingTag = typename Domain::GridFunctionSpace::Traits::OrderingTag;
    using currentBlocking = typename Domain::GridFunctionSpace::Traits::Backend;
    static const bool isLexicographic = std::is_same<LexicographicOrderingTag,currentOrderingTag>::value;
    static const bool isEntityBlocked = std::is_same<EntityBlockedOrderingTag,currentOrderingTag>::value;
    static const bool isFixedBlocking = currentBlocking::Traits::blocked;
    using type0 = typename std::conditional<
      isLexicographic,
      CorrectSolutionLexicographic<ParameterClass>,
      CorrectSolutionSingle<ParameterClass>
      >::type;
    using type1 = typename std::conditional<
      isEntityBlocked,
      CorrectSolutionEntityBlockedNoneBlocking<ParameterClass>,
      type0
      >::type;
  public:
    using type = typename std::conditional<
      isFixedBlocking,
      CorrectSolutionEntityBlockedFixedBlocking<ParameterClass>,
      type1
      >::type;
  };

  template<typename Newton>
  class BoundedLineSearchParametersInterface
  {
  public:
    using Domain = typename Newton::Domain;
    using Real = typename Newton::Real;
    using ThisType = BoundedLineSearchParametersInterface<Newton>;

    /* \Brief set restrains
    *
    * Sets bounds for certain values, line search will keep them
    * inside the defined interval. Allows seeting lower bound, upper
    * bound, or both.
    *
    * Bounds are loaded from ParameterTree, which is passed from Newton
    * class via subtree "line_search". To restrain some values, these
    * fields are expected (in subtree line_search passed to Newton class):
    * number = std::size_t, states the number of restrained unknowns,
    *          all following must end with number, starting with 0
    * block0 = std::size_t, the place of restrained unknown in the block,
    *          e.g. block0=2 if we restrain the third unknown of the system
    * block1 = std::size_t, if we restrain two unknowns, etc
    * lowerType0 = bool, true if lower bound is set
    * upperType0 = bool, true if upper bound is set,
    *             default false, at least one of lowerType, upperType must be true
    * lowerBound0 = Real, templated type, lower bound for unknown number block0
    * upperBound1 = Real, templated type, upper bound for unknown number block1
    *
    * Currectly works only for EntityBlockedOrdering.
    */
    void setBoundedParameters(const ParameterTree& parameterTree)
    {
      std::size_t numberOfRestrained = parameterTree.get<std::size_t>("number",0);
      if (numberOfRestrained == 0)
        DUNE_THROW(Exception,"BoundedLineSearch setParameters error: use normal LineSearch if no value is restrained");
      Impl::BlockSizeExtractor<Domain> bse;
      blocksize = bse.value;
      for (std::size_t i=0; i < numberOfRestrained; ++i)
      {
        std::size_t tmpsizet = parameterTree.get<std::size_t>("block"+std::to_string(i),std::size_t(-1));
        if (tmpsizet == std::size_t(-1))
          DUNE_THROW(Exception,"BoundedLineSearch setParameters error: Missing line_search.block"+std::to_string(i)+" field");
        if (tmpsizet >= blocksize)
          DUNE_THROW(Exception,"BoundedLineSearch setParameters error: Restraining unknown "+std::to_string(tmpsizet)+" in a system of size "+std::to_string(blocksize));
        auto tmpboolL = parameterTree.get<bool>("lowerType"+std::to_string(i),false);
        auto tmpboolU = parameterTree.get<bool>("upperType"+std::to_string(i),false);
        if (tmpboolL == false && tmpboolU == false)
          DUNE_THROW(Exception,"BoundedLineSearch setParameters error: Invalid parameters, both lower and upper value are unrestrained in block"+std::to_string(i));
        // bounds may be unused if not restrained (tmpbool==false)
        Real tmpRealL = 0.;
        Real tmpRealU = 0.;
        if (tmpboolL)
          tmpRealL = parameterTree.get<Real>("lowerBound"+std::to_string(i));
        if (tmpboolU)
          tmpRealU = parameterTree.get<Real>("upperBound"+std::to_string(i));
        if (tmpboolL && tmpboolU && tmpRealL >= tmpRealU)
          DUNE_THROW(Exception,"BoundedLineSearch setParameters error: Invalid parameters, block"+std::to_string(i)+" has lower bound greater than upper bound");
        if (tmpboolL)
        {
          blockLower.push_back(tmpsizet);
          boundLower.push_back(tmpRealL);
        }
        if (tmpboolU)
        {
          blockUpper.push_back(tmpsizet);
          boundUpper.push_back(tmpRealU);
        }
      }
    }

  // private:
    std::vector<std::size_t> blockLower; // stores which parts of block are restrained
    std::vector<std::size_t> blockUpper; // stores which parts of block are restrained
    std::size_t blocksize;
    std::vector<Real> boundLower; // remembers lower bounds
    std::vector<Real> boundUpper; // remembers upper bounds
  }; // end BoundedLineSearchParametersInterface */
} // end namespace Impl

//! Class for simply updating the solution without line search,
//  after performing step, solution is changed to fit into given bounds
template <typename Newton>
class BoundedLineSearchNone : public LineSearchInterface<typename Newton::Domain>
{
  using ParameterClass = typename Impl::BoundedLineSearchParametersInterface<Newton>::ThisType;
public:
  using Domain = typename Newton::Domain;
  using Real = typename Newton::Real;

  BoundedLineSearchNone(Newton& newton) : _newton(newton)
    , corSol(param)
  {}

  //! Do line search (in this case just update the solution)
  virtual void lineSearch(Domain& solution, const Domain& correction) override
  {
    solution.axpy(-1., correction);
    corSol(solution);
    _newton.updateDefect(solution);
  }

  virtual void setParameters(const ParameterTree& parameterTree) override
  {
    param.setBoundedParameters(parameterTree);
  }

private:
  Newton& _newton;
  ParameterClass param;
  using CorSol = typename Impl::SetCorrectSolution<Domain,ParameterClass>::type;
  CorSol corSol;
};

/** \brief bounded Hackbusch-Reusken line search
 *
 * Hackbusch-Reusken line search with an extra feature.
 * Allows setting bounds on values which will not be exceeded,
 * after performing line search, solution vector is corrected to stay inside limits.
 *
 * If the parameter line_search_accept_best is set through the setParameters
 * method this line search will simply return the best result even if it did
 * not converge.
 */
template <typename Newton>
class BoundedLineSearchHackbuschReusken : public LineSearchInterface<typename Newton::Domain>
{
  using ParameterClass = Impl::BoundedLineSearchParametersInterface<Newton>;
public:
  using Domain = typename Newton::Domain;
  using Real = typename Newton::Real;

  BoundedLineSearchHackbuschReusken(Newton& newton) : _newton(newton)
    , corSol(param)
  {}

  //! Do line search
  virtual void lineSearch(Domain& solution, const Domain& correction) override
  {
    if ((_newton.result().defect < _newton.getAbsoluteLimit())){
      solution.axpy(-1., correction);
      // param.correctSolution(solution);
      corSol(solution);
      _newton.updateDefect(solution);
      return;
    }

    auto verbosity = _newton.getVerbosityLevel();

    if (verbosity >= 4)
      std::cout << "      Performing line search..." << std::endl;
    Real lambda = 1.;
    Real bestLambda = 0.0;
    Real bestDefect = _newton.result().defect;
    Real previousDefect = _newton.result().defect;
    bool converged = false;

    if (not _previousSolution)
      _previousSolution = std::make_shared<Domain>(solution);
    else
      *_previousSolution = solution;

    for (unsigned int iteration = 0; iteration < _lineSearchMaxIterations; ++iteration){
      if (verbosity >= 4)
        std::cout << "          trying line search damping factor:   "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << lambda
                  << std::endl;

      solution.axpy(-lambda, correction);
      corSol(solution);
      // param.correctSolution(solution);
      _newton.updateDefect(solution);
      if (verbosity >= 4){
        if (not std::isfinite(_newton.result().defect))
          std::cout << "          NaNs detected" << std::endl;
      }       // ignore NaNs and try again with lower lambda

      if (_newton.result().defect <= (1.0 - lambda/4) * previousDefect){
        if (verbosity >= 4)
          std::cout << "          line search converged" << std::endl;
        converged = true;
        break;
      }

      if (_newton.result().defect < bestDefect){
        bestDefect = _newton.result().defect;
        bestLambda = lambda;
      }

      lambda *= _lineSearchDampingFactor;
      solution = *_previousSolution;
    }

    if (not converged){
      if (verbosity >= 4)
        std::cout << "          max line search iterations exceeded" << std::endl;

      if (not _acceptBest){
        solution = *_previousSolution;
        _newton.updateDefect(solution);
        DUNE_THROW(NewtonLineSearchError,
                   "NewtonLineSearch::line_search(): line search failed, "
                   "max iteration count reached, "
                   "defect did not improve enough");
      }
      else{
        if (bestLambda == 0.0){
          solution = *_previousSolution;
          _newton.updateDefect(solution);
          DUNE_THROW(NewtonLineSearchError,
                     "NewtonLineSearch::line_search(): line search failed, "
                     "max iteration count reached, "
                     "defect did not improve in any of the iterations");
        }
        if (bestLambda != lambda){
          solution = *_previousSolution;
          solution.axpy(-bestLambda, correction);
          _newton.updateDefect(solution);
          converged = true;
        }
      }
    }

    if (converged)
      if (verbosity >= 4)
        std::cout << "          line search damping factor:   "
                  << std::setw(12) << std::setprecision(4) << std::scientific
                  << lambda << std::endl;
  }

  /* \brief Set parameters
   *
   * Possible parameters are:
   *
   * - line_search_max_iterations: Maximum number of line search iterations.
   *
   * - line_search_damping_factor: Multiplier to line search parameter after each iteration.
   *
   * - line_search_accept_best: Accept the best line search parameter if
   *   there was any improvement, even if the convergence criterion was not
   *   reached.
   */
  virtual void setParameters(const ParameterTree& parameterTree) override
  {
    _lineSearchMaxIterations = parameterTree.get<unsigned int>("line_search_max_iterations",
                                                               _lineSearchMaxIterations);
    _lineSearchDampingFactor = parameterTree.get<Real>("line_search_damping_factor",
                                                       _lineSearchDampingFactor);
    _acceptBest = parameterTree.get<bool>("line_search_accept_best",
                                          _acceptBest);
    param.setBoundedParameters(parameterTree);
  }

private:
  Newton& _newton;
  std::shared_ptr<Domain> _previousSolution;
  ParameterClass param;
  using CorSol = typename Impl::SetCorrectSolution<Domain,ParameterClass>::type;
  CorSol corSol;

  // Line search parameters
  unsigned int _lineSearchMaxIterations = 10;
  Real _lineSearchDampingFactor = 0.5;
  bool _acceptBest = false;
};
#endif // #ifdef DUNE_PDELAB_SOLVER_NEWTONLINESEARCH_BOUNDED_ENABLED