// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_SOLVER_NEWTONLINESEARCH_BOUNDED_HH
#define DUNE_PDELAB_SOLVER_NEWTONLINESEARCH_BOUNDED_HH
#endif
/* \brief BoundedLineSearch
 *
 * This file contains the implementation of the BoundedLineSearchNone
 * and BoundedLineSearchHackbuschReusken. They inherit from LineSearchInterface
 * and must be defined before LineSearchStrategy. Both of them are in
 * newtonlinesearch.hh, this code must be placed between them.
 *
 * To use this file, simply include it before newtonlinesearch.hh.
 * The "guard" macro on top will unlock code parts dealing with
 * BoundedLineSearches and #inlude the following code in between
 * LineSearchInterface and LineSearchStrategy.
 *
 * Such construction necessitates usage of 3 guard macros:
 * 1. Unlock: is on top of this file. Makes changes to newtonlineseach.hh
 * invisible if unused, and makes this code accessible with #include.
 * 2. Place: is in newtonlinesearch.hh ensuring BoundedLineSearches are
 * loaded after LineSearchInterface and before LineSearchStrategy.
 * 3. Guard: prevents loading BoundedLineSearches for the second time.
 */

#ifdef DUNE_PDELAB_SOLVER_NEWTONLINESEARCH_BOUNDED_ENABLED
#ifndef DUNE_PDELAB_SOLVER_NEWTONLINESEARCH_BOUNDED_READ
#define DUNE_PDELAB_SOLVER_NEWTONLINESEARCH_BOUNDED_READ
// Note: it is assumed this file is included into namespace Dune::PDELab
namespace Impl
{
  // Simple struct for calculating the number of unknowns of the system.
  // This works for GFS which are not Powerspace, i.e. have no CHILDREN.
  // Simply putting this std::max(..) formula into code to calculate
  // systemsize does not work.
  template<typename Domain>
  struct SystemSizeExtractor
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
   * Other options like bcrs blocking, InterleavedOrderingTag, or CompositeSpaces
   * are not tested/implemented.
   *
   * Heavily relies on ParameterClass = BoundedLineSearchParametersInterface<Newton>
   * which carries all information about which unknowns are bounded and how.
   */
  template <typename Newton> class BoundedLineSearchParametersInterface;

  template<typename Newton>
  class CorrectSolutionEntityBlockedFixedBlocking
  {
    using ParameterClass = BoundedLineSearchParametersInterface<Newton>;
    const ParameterClass& pc;
    using Domain = typename Newton::Domain;
  public:
    CorrectSolutionEntityBlockedFixedBlocking(const ParameterClass& pc_)
      : pc(pc_)
    {}

    void operator() (Domain& solution)
    {
      using Dune::PDELab::Backend::native;
      if (pc.placeLower.size()>0 && pc.placeUpper.size()>0)
        for (auto& block : native(solution))
        {
          for (std::size_t i=0; i<pc.placeLower.size(); ++i)
            if (block[pc.placeLower[i]] < pc.boundLower[i])
              block[pc.placeLower[i]] = pc.boundLower[i];
          for (std::size_t i=0; i<pc.placeUpper.size(); ++i)
            if (block[pc.placeUpper[i]] > pc.boundUpper[i])
              block[pc.placeUpper[i]] = pc.boundUpper[i];
        }
        else
        {
          if (pc.placeLower.size()>0)
            for (auto& block : native(solution))
              for (std::size_t i=0; i<pc.placeLower.size(); ++i)
                if (block[pc.placeLower[i]] < pc.boundLower[i])
                  block[pc.placeLower[i]] = pc.boundLower[i];
          if (pc.placeUpper.size()>0)
            for (auto& block : native(solution))
              for (std::size_t i=0; i<pc.placeUpper.size(); ++i)
                if (block[pc.placeUpper[i]] > pc.boundUpper[i])
                  block[pc.placeUpper[i]] = pc.boundUpper[i];
        }
    } // end operator()
  };

  template<typename Newton>
  class CorrectSolutionEntityBlockedNoneBlocking
  {
    using ParameterClass = BoundedLineSearchParametersInterface<Newton>;
    const ParameterClass& pc;
    using Domain = typename Newton::Domain;
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
      std::size_t systemsize = pc.systemsize;

      // Check if both upper and lower bounds are active. Code is copied
      // for cases when only one is active in order to avoid traversing
      // the vector twice or asking frequently for nonemptyness.
      if (pc.placeLower.size()>0 && pc.placeUpper.size()>0)
      {
        std::size_t low=0; // watches the next lower bound to be corrected
        std::size_t up=0; // watches the next upper bound
        std::size_t lowsize = pc.placeLower.size();
        std::size_t upsize = pc.placeUpper.size();
        for (auto& block : native(solution))
        {
          // Firstly check whether it is one of the bounded variables.
          if (i==pc.placeLower[low])
          {
            if (block[0]<pc.boundLower[low])
              block[0] = pc.boundLower[low];
            low = (low+1) % lowsize;
          }
          if (i==pc.placeUpper[up])
          {
            if (block[0]>pc.boundUpper[up])
              block[0] = pc.boundUpper[up];
            up = (up+1) % upsize;
          }
          i = (i+1) % systemsize;
        }
      }
      else
      {
        if (pc.placeLower.size()>0)
        {
          std::size_t low=0; // watches the next lower bound to be corrected
          std::size_t lowsize = pc.placeLower.size();
          for (auto& block : native(solution))
          {
            if (i==pc.placeLower[low])
            {
              if (block[0]<pc.boundLower[low])
                block[0] = pc.boundLower[low];
              low = (low+1) % lowsize;
            }
            i = (i+1) % systemsize;
          }
        }
        if (pc.placeUpper.size()>0)
        {
          std::size_t up=0; // watches the next lower bound to be corrected
          std::size_t upsize = pc.placeUpper.size();
          for (auto& block : native(solution))
          {
            if (i==pc.placeUpper[up])
            {
              if (block[0]>pc.boundUpper[up])
                block[0] = pc.boundUpper[up];
              up = (up+1) % upsize;
            }
            i = (i+1) % systemsize;
          }
        }
      }
    } // end operator()
  };

  template<typename Newton>
  class CorrectSolutionLexicographic
  {
    using ParameterClass = BoundedLineSearchParametersInterface<Newton>;
    const ParameterClass& pc;
    using Domain = typename Newton::Domain;
  public:
    CorrectSolutionLexicographic(const ParameterClass& pc_)
      : pc(pc_)
    {}

    void operator() (Domain& solution)
    {
      using Dune::PDELab::Backend::native;

      std::size_t sol_size = native(solution).size()/pc.systemsize;
      for (std::size_t i=0; i<pc.placeLower.size(); ++i)
        for (std::size_t j=sol_size*pc.placeLower[i]; j<sol_size*(pc.placeLower[i]+1); ++j)
          if (native(solution)[j] < pc.boundLower[i])
            native(solution)[j] = pc.boundLower[i];
      for (std::size_t i=0; i<pc.placeUpper.size(); ++i)
        for (std::size_t j=sol_size*pc.placeUpper[i]; j<sol_size*(pc.placeUpper[i]+1); ++j)
          if (native(solution)[j] > pc.boundUpper[i])
            native(solution)[j] = pc.boundUpper[i];
    }
  };

  template<typename Newton>
  class CorrectSolutionSingle
  {
    using ParameterClass = BoundedLineSearchParametersInterface<Newton>;
    const ParameterClass& pc;
    using Domain = typename Newton::Domain;
  public:
    CorrectSolutionSingle(const ParameterClass& pc_)
      : pc(pc_)
    {}

    void operator() (Domain& solution)
    {
      // avoid traversing the vector twice (or multiple size()>0 checks)
      if (pc.placeUpper.size()>0 && pc.placeLower.size()>0)
        for (auto& v : solution)
        {
          if (v>pc.boundUpper[0])
            v = pc.boundUpper[0];
          if (v<pc.boundLower[0])
            v = pc.boundLower[0];
        }
      else
      {
        if (pc.placeUpper.size()>0)
          for (auto& v : solution)
            if (v>pc.boundUpper[0])
              v = pc.boundUpper[0];
        if (pc.placeLower.size()>0)
          for (auto& v : solution)
            if (v<pc.boundLower[0])
              v = pc.boundLower[0];
      }
    }
  };

  /* \brief SetCorrectSolution::type
   *
   * This class holds the correct type for CorrectSolution procedure.
   * Picks one of above 4 types automatically from Domain type.
   */
  template<typename Newton>
  class SetCorrectSolution
  {
    using Domain = typename Newton::Domain;
    using currentOrderingTag = typename Domain::GridFunctionSpace::Traits::OrderingTag;
    using currentBlocking = typename Domain::GridFunctionSpace::Traits::Backend;
    static const bool isLexicographic = std::is_same<LexicographicOrderingTag,currentOrderingTag>::value;
    static const bool isEntityBlocked = std::is_same<EntityBlockedOrderingTag,currentOrderingTag>::value;
    static const bool isFixedBlocking = currentBlocking::Traits::blocked;
    using type0 = typename std::conditional<
      isLexicographic,
      CorrectSolutionLexicographic<Newton>,
      CorrectSolutionSingle<Newton>
      >::type;
    using type1 = typename std::conditional<
      isEntityBlocked,
      CorrectSolutionEntityBlockedNoneBlocking<Newton>,
      type0
      >::type;
  public:
    using type = typename std::conditional<
      isFixedBlocking,
      CorrectSolutionEntityBlockedFixedBlocking<Newton>,
      type1
      >::type;
  };

  // Forward declaration enabling (template specialized) friendship
  // with BoundedLineSearchParametersInterface.
  template <typename Newton> class BoundedLineSearchNone;
  template <typename Newton> class BoundedLineSearchHackbuschReusken;

  class ExceptionBoundedLineSearch : Dune::Exception {};

  /* \brief BoundedLineSearchParametersInterface
   *
   * Class containing information about bounded unknowns. Stores which
   * are constrained, and their upper or lower bounds.
   *
   * Bounds are set by setBoundedParameters function. Private data are
   * accessible by declaring friendship with CorrectSolution. classes.
   */
  template<typename Newton>
  class BoundedLineSearchParametersInterface
  {
  public:
    friend class CorrectSolutionEntityBlockedFixedBlocking<Newton>;
    friend class CorrectSolutionEntityBlockedNoneBlocking<Newton>;
    friend class CorrectSolutionLexicographic<Newton>;
    friend class CorrectSolutionSingle<Newton>;
    using Domain = typename Newton::Domain;
    using Real = typename Newton::Real;

    /* \brief set restrains
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
    */
    // void setBoundedParameters(const ParameterTree& parameterTree)
    // {
    //   std::size_t numberOfRestrained = parameterTree.get<std::size_t>("number",0);
    //   if (numberOfRestrained == 0)
    //     DUNE_THROW(Exception,"BoundedLineSearch setParameters error: use normal LineSearch if no value is restrained");
    //   Impl::SystemSizeExtractor<Domain> bse;
    //   systemsize = bse.value;
    //   for (std::size_t i=0; i < numberOfRestrained; ++i)
    //   {
    //     std::size_t tmpsizet = parameterTree.get<std::size_t>("block"+std::to_string(i),std::size_t(-1));
    //     if (tmpsizet == std::size_t(-1))
    //       DUNE_THROW(Exception,"BoundedLineSearch setParameters error: Missing line_search.block"+std::to_string(i)+" field");
    //     if (tmpsizet >= systemsize)
    //       DUNE_THROW(Exception,"BoundedLineSearch setParameters error: Restraining unknown "+std::to_string(tmpsizet)+" in a system of size "+std::to_string(systemsize));
    //     auto tmpboolL = parameterTree.get<bool>("lowerType"+std::to_string(i),false);
    //     auto tmpboolU = parameterTree.get<bool>("upperType"+std::to_string(i),false);
    //     if (tmpboolL == false && tmpboolU == false)
    //       DUNE_THROW(Exception,"BoundedLineSearch setParameters error: Invalid parameters, both lower and upper value are unrestrained in block"+std::to_string(i));
    //     // bounds may be unused if not restrained (tmpbool==false)
    //     Real tmpRealL = 0.;
    //     Real tmpRealU = 0.;
    //     if (tmpboolL)
    //       tmpRealL = parameterTree.get<Real>("lowerBound"+std::to_string(i));
    //     if (tmpboolU)
    //       tmpRealU = parameterTree.get<Real>("upperBound"+std::to_string(i));
    //     if (tmpboolL && tmpboolU && tmpRealL >= tmpRealU)
    //       DUNE_THROW(Exception,"BoundedLineSearch setParameters error: Invalid parameters, block"+std::to_string(i)+" has lower bound greater than upper bound");
    //     if (tmpboolL)
    //     {
    //       placeLower.push_back(tmpsizet);
    //       boundLower.push_back(tmpRealL);
    //     }
    //     if (tmpboolU)
    //     {
    //       placeUpper.push_back(tmpsizet);
    //       boundUpper.push_back(tmpRealU);
    //     }
    //   }
    // }

    /* \brief set restrains
    *
    * Sets bounds for certain values, which line search will keep
    * inside the defined interval. Allows seeting lower bound, upper
    * bound, or both.
    *
    * Bounds are loaded from ParameterTree, which is passed from Newton
    * class via subtree "line_search". To restrain some values, these
    * fields are expected (in subtree line_search passed to Newton class):
    * lowerBound* = Real, templated type, lower bound for unknown number *
    * upperBound* = Real, templated type, upper bound for unknown number *
    * numberofrestraints = std:size_t, optional parameter, the total number
    *   of restraints specified (good for avoiding typos)
    *
    * The system size is deduced automatically.
    */
    void setBoundedParameters(const ParameterTree& parameterTree, unsigned int verbosity=0)
    {
      Impl::SystemSizeExtractor<Domain> bse;
      systemsize = bse.value;
      placeLower.clear();
      placeUpper.clear();
      boundLower.clear();
      boundUpper.clear();
      if (systemsize==1)
      {
        bool hasLower = parameterTree.hasKey("lowerBound");
        bool hasLower0 = parameterTree.hasKey("lowerBound0");
        bool hasUpper = parameterTree.hasKey("upperBound");
        bool hasUpper0 = parameterTree.hasKey("upperBound0");
        if (hasLower)
        {
          placeLower.push_back(0);
          Real bound = parameterTree.get<Real>("lowerBound");
          boundLower.push_back(bound);
          if (hasLower0)
          {
            if (bound != parameterTree.get<Real>("lowerBound0"))
              DUNE_THROW(Exception,"BoundedLineSearch parameters error: Two different lower bounds are set!");
          }
        }
        else if (hasLower0)
        {
          placeLower.push_back(0);
          boundLower.push_back(parameterTree.get<Real>("lowerBound0"));
        }
        if (hasUpper)
        {
          placeUpper.push_back(0);
          Real bound = parameterTree.get<Real>("upperBound");
          boundUpper.push_back(bound);
          if (hasUpper0)
          {
            if (bound != parameterTree.get<Real>("upperBound0"))
              DUNE_THROW(Exception,"BoundedLineSearch parameters error: Two different upper bounds are set!");
          }
        }
        else if (hasUpper0)
        {
          placeUpper.push_back(0);
          boundUpper.push_back(parameterTree.get<Real>("upperBound0"));
        }
        if ((hasLower || hasLower0) && (hasUpper || hasUpper0))
          if (boundLower[0] > boundUpper[0])
            DUNE_THROW(Exception,"BoundedLineSearch parameters error: LowerBound higher than upperBound.");

      }
      else // systemsize !=1
      {
        for (std::size_t i=0; i < systemsize; ++i)
        {
          bool hasLower = parameterTree.hasKey("lowerBound"+std::to_string(i));
          bool hasUpper = parameterTree.hasKey("upperBound"+std::to_string(i));
          if (hasLower)
          {
            placeLower.push_back(i);
            boundLower.push_back(parameterTree.get<Real>("lowerBound"+std::to_string(i)));
          }
          if (hasUpper)
          {
            placeUpper.push_back(i);
            boundUpper.push_back(parameterTree.get<Real>("upperBound"+std::to_string(i)));
          }
          if (hasLower && hasUpper)
            if (boundLower[boundLower.size()-1] > boundUpper[boundUpper.size()-1])
              DUNE_THROW(Exception,"BoundedLineSearch parameters error: Unknown "+std::to_string(i)+" has lowerBound higher than upperBound.");
        }
      }
      // optional check
      if (parameterTree.hasKey("numberofrestraints"))
      {
        auto nrestr = parameterTree.get<std::size_t>("numberofrestraints");
        if (nrestr != placeLower.size() + placeUpper.size())
          DUNE_THROW(Exception,"BoundedLineSearch parameters error: The number of restraints is different to the specified number.");
        if (verbosity>=2)
          std::cout << "The system has " << systemsize << " unknowns and "
                    << nrestr << " restraints placed on them." << std::endl;
      }
      if (verbosity >=4)
      {
        if (placeLower.size()==0)
          std::cout << "There are no lower bounds." << std::endl;
        else
        {
          for (std::size_t i=0; i<placeLower.size(); ++i)
          {
            std::cout << "Lower bound of unknown " << placeLower[i] << " is " << boundLower[i] <<std::endl;
          }
        }
        if (placeUpper.size()==0)
          std::cout << "There are no upper bounds." << std::endl;
        else
        {
          for (std::size_t i=0; i<placeUpper.size(); ++i)
          {
            std::cout << "Upper bound of unknown " << placeUpper[i] << " is " << boundUpper[i] <<std::endl;
          }
        }
      }
    }

  private:
    std::vector<std::size_t> placeLower; // stores which parts of block are restrained
    std::vector<std::size_t> placeUpper; // stores which parts of block are restrained
    std::size_t systemsize = 0; // if setBoundedParameters is not used (which definitely should), this leads to zero division
    std::vector<Real> boundLower; // remembers lower bounds
    std::vector<Real> boundUpper; // remembers upper bounds
  }; // BoundedLineSearchParametersInterface
} // namespace Impl

//! Class for simply updating the solution without line search,
//  after performing step, solution is changed to fit into given bounds
template <typename Newton>
class BoundedLineSearchNone : public LineSearchInterface<typename Newton::Domain>
{
  using ParameterClass = typename Impl::BoundedLineSearchParametersInterface<Newton>;
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
    param.setBoundedParameters(parameterTree,_newton.getVerbosityLevel());
  }

private:
  Newton& _newton;
  ParameterClass param;
  using CorSol = typename Impl::SetCorrectSolution<Newton>::type;
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
        DUNE_THROW(LineSearchError,
                   "LineSearch::line_search(): line search failed, "
                   "max iteration count reached, "
                   "defect did not improve enough");
      }
      else{
        if (bestLambda == 0.0){
          solution = *_previousSolution;
          _newton.updateDefect(solution);
          DUNE_THROW(LineSearchError,
                     "LineSearch::line_search(): line search failed, "
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
  using CorSol = typename Impl::SetCorrectSolution<Newton>::type;
  CorSol corSol;

  // Line search parameters
  unsigned int _lineSearchMaxIterations = 10;
  Real _lineSearchDampingFactor = 0.5;
  bool _acceptBest = false;
};
#endif // #ifdef DUNE_PDELAB_SOLVER_NEWTONLINESEARCH_BOUNDED_READ
#endif // #ifdef DUNE_PDELAB_SOLVER_NEWTONLINESEARCH_BOUNDED_ENABLED