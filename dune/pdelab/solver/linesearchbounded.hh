// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_SOLVER_LINESEARCH_BOUNDED_HH
#define DUNE_PDELAB_SOLVER_LINESEARCH_BOUNDED_HH
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

#ifdef DUNE_PDELAB_SOLVER_LINESEARCH_BOUNDED_ENABLED
#ifndef DUNE_PDELAB_SOLVER_LINESEARCH_BOUNDED_READ
#define DUNE_PDELAB_SOLVER_LINESEARCH_BOUNDED_READ
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

  class LineSearchError : public Dune::Exception {};

  /* \brief BoundedLineSearchParametersInterface
   *
   * Class containing information about bounded unknowns. Stores which
   * are constrained, and their upper or lower bounds.
   *
   * Bounds are set by setBoundedParameters function. Get functions
   * provide access to data vectors.
   */
  template<typename Newton>
  class BoundedLineSearchParametersInterface
  {
  public:
    using Domain = typename Newton::Domain;
    using Real = typename Newton::Real;

    /* \brief set restrains
    *
    * Sets bounds for certain values, which line search will keep
    * inside the defined interval. Allows seeting lower bound, upper
    * bound, or both.
    *
    * Bounds are loaded from ParameterTree, which is passed from Newton
    * class via subtree "LineSearch". To restrain some values, these
    * fields are expected (in subtree line_search passed to Newton class):
    * LowerBound* = Real, lower bound for unknown number *
    * UpperBound* = Real, upper bound for unknown number *
    * NumberOfRestraints = std:size_t, optional parameter, the total number
    *   of restraints specified (good for avoiding typos)
    * In case of a problem with only one unknown LowerBound and UpperBound
    * are accepted too, it is the same as if the missing number was 0.
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
        bool hasLower = parameterTree.hasKey("LowerBound");
        bool hasLower0 = parameterTree.hasKey("LowerBound0");
        bool hasUpper = parameterTree.hasKey("UpperBound");
        bool hasUpper0 = parameterTree.hasKey("UpperBound0");
        if (hasLower)
        {
          placeLower.push_back(0);
          Real bound = parameterTree.get<Real>("LowerBound");
          boundLower.push_back(bound);
          if (hasLower0)
          {
            if (bound != parameterTree.get<Real>("LowerBound0"))
              DUNE_THROW(LineSearchError,"BoundedLineSearch parameters error: Two different lower bounds are set!");
          }
        }
        else if (hasLower0)
        {
          placeLower.push_back(0);
          boundLower.push_back(parameterTree.get<Real>("LowerBound0"));
        }
        if (hasUpper)
        {
          placeUpper.push_back(0);
          Real bound = parameterTree.get<Real>("UpperBound");
          boundUpper.push_back(bound);
          if (hasUpper0)
          {
            if (bound != parameterTree.get<Real>("UpperBound0"))
              DUNE_THROW(LineSearchError,"BoundedLineSearch parameters error: Two different upper bounds are set!");
          }
        }
        else if (hasUpper0)
        {
          placeUpper.push_back(0);
          boundUpper.push_back(parameterTree.get<Real>("UpperBound0"));
        }
        if ((hasLower || hasLower0) && (hasUpper || hasUpper0))
          if (boundLower[0] > boundUpper[0])
            DUNE_THROW(LineSearchError,"BoundedLineSearch parameters error: LowerBound higher than UpperBound.");
      }
      else // systemsize !=1
      {
        for (std::size_t i=0; i < systemsize; ++i)
        {
          bool hasLower = parameterTree.hasKey("LowerBound"+std::to_string(i));
          bool hasUpper = parameterTree.hasKey("UpperBound"+std::to_string(i));
          if (hasLower)
          {
            placeLower.push_back(i);
            boundLower.push_back(parameterTree.get<Real>("LowerBound"+std::to_string(i)));
          }
          if (hasUpper)
          {
            placeUpper.push_back(i);
            boundUpper.push_back(parameterTree.get<Real>("UpperBound"+std::to_string(i)));
          }
          if (hasLower && hasUpper)
            if (boundLower[boundLower.size()-1] > boundUpper[boundUpper.size()-1])
              DUNE_THROW(LineSearchError,"BoundedLineSearch parameters error: Unknown "+std::to_string(i)+" has LowerBound higher than UpperBound.");
        }
      }
      // optional check, "NumberOfRestraints" stores the number of restraints.
      // Good for detecting typos, since ParameterTree does not recognize faulty arguments.
      if (parameterTree.hasKey("NumberOfRestraints"))
      {
        auto nrestr = parameterTree.get<std::size_t>("NumberOfRestraints");
        if (nrestr != placeLower.size() + placeUpper.size())
          DUNE_THROW(LineSearchError,"BoundedLineSearch parameters error: The number of restraints is different to the specified number.");
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

    // "get functions" providing access to data
    const std::vector<std::size_t>& getPlaceLower() const
    {
      return placeLower;
    }
    const std::vector<std::size_t>& getPlaceUpper() const
    {
      return placeUpper;
    }
    const std::vector<Real>& getBoundLower() const
    {
      return boundLower;
    }
    const std::vector<Real>& getBoundUpper() const
    {
      return boundUpper;
    }
    std::size_t getSystemSize() const
    {
      return systemsize;
    }

    void printParameters() const
    {
      using std::cout;
      cout << "Number of Bounds on variables: " << placeLower.size()+placeUpper.size() << std::endl;
      if (placeLower.size() > 0)
      {
        cout << "Lower bounds on variable" << (placeLower.size()==1 ? " " : "s ");
        for (const auto& v : placeLower)
          cout << v << ", ";
        cout << (placeLower.size()==1 ? "is " : "are ");
        for (const auto& v : boundLower)
          cout << v << ", ";
        cout << std::endl;
      }
      if (placeUpper.size() > 0)
      {
        cout << "Upper bounds on variable" << (placeUpper.size()==1 ? " " : "s ");
        for (const auto& v : placeUpper)
          cout << v << ", ";
        cout << (placeUpper.size()==1 ? "is " : "are ");
        for (const auto& v : boundUpper)
          cout << v << ", ";
        cout << std::endl;
      }
    }

  private:
    std::vector<std::size_t> placeLower; // stores which parts of block are restrained
    std::vector<std::size_t> placeUpper; // stores which parts of block are restrained
    std::size_t systemsize = 0; // if setBoundedParameters is not used (which definitely should), this leads to zero division
    std::vector<Real> boundLower; // remembers lower bounds
    std::vector<Real> boundUpper; // remembers upper bounds
  }; // BoundedLineSearchParametersInterface

  /* \brief CorrectSolution
   *
   * 4 classes with a single purpose: traverse the vector and set specified
   * values into respective bounds. Classes differ based on the type
   * of the vector. The correct type is guessed automatically at the compile
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
   * are not tested/implemented. Default option is CorrectSolutionSingle,
   * i.e. if neither of other types is recognized, this is the result.
   *
   * We treat separately cases when both upper and lower bounds are present
   * and when only one is - to avoid traversing the vector twice or
   * repeatingly checking the emptyness of the vector storing bounds.
   *
   * Rely on ParameterClass = BoundedLineSearchParametersInterface<Newton>
   * which carries information about which unknowns are bounded and how.
   */

  template<typename Newton>
  class CorrectSolutionEntityBlockedFixedBlocking
  {
  public:
    using ParameterClass = BoundedLineSearchParametersInterface<Newton>;
    using Domain = typename Newton::Domain;

    explicit CorrectSolutionEntityBlockedFixedBlocking(const ParameterClass& pc_)
      : pc(pc_)
    {}

    void operator() (Domain& solution)
    {
      auto& placeLower = pc.getPlaceLower();
      auto& placeUpper = pc.getPlaceUpper();
      auto& boundLower = pc.getBoundLower();
      auto& boundUpper = pc.getBoundUpper();
      using Dune::PDELab::Backend::native;

      // Traverse the vector and set the values inside the bounds.
      if (placeLower.size()>0 && placeUpper.size()>0)
        for (auto& block : native(solution))
        {
          for (std::size_t i=0; i<placeLower.size(); ++i)
            if (block[placeLower[i]] < boundLower[i])
              block[placeLower[i]] = boundLower[i];
          for (std::size_t i=0; i<placeUpper.size(); ++i)
            if (block[placeUpper[i]] > boundUpper[i])
              block[placeUpper[i]] = boundUpper[i];
        }
      else
      {
        if (placeLower.size()>0)
          for (auto& block : native(solution))
            for (std::size_t i=0; i<placeLower.size(); ++i)
              if (block[placeLower[i]] < boundLower[i])
                block[placeLower[i]] = boundLower[i];
        if (placeUpper.size()>0)
          for (auto& block : native(solution))
            for (std::size_t i=0; i<placeUpper.size(); ++i)
              if (block[placeUpper[i]] > boundUpper[i])
                block[placeUpper[i]] = boundUpper[i];
      }
    } // end operator()

  private:
    const ParameterClass& pc;
  };

  template<typename Newton>
  class CorrectSolutionEntityBlockedNoneBlocking
  {
  public:
    using ParameterClass = BoundedLineSearchParametersInterface<Newton>;
    using Domain = typename Newton::Domain;

    explicit CorrectSolutionEntityBlockedNoneBlocking(const ParameterClass& pc_)
      : pc(pc_)
    {}

    void operator() (Domain& solution)
    {
      auto& placeLower = pc.getPlaceLower();
      auto& placeUpper = pc.getPlaceUpper();
      auto& boundLower = pc.getBoundLower();
      auto& boundUpper = pc.getBoundUpper();

      using Dune::PDELab::Backend::native;
      // EntityBlocked with Blocking "none" is arranged in blocks, but
      // blocks have size 1. Variable i stores the position in the block.
      std::size_t i=0;
      std::size_t systemsize = pc.getSystemSize();

      // Check if both upper and lower bounds are active. Code is recycled
      // for cases when only one is active in order to avoid traversing
      // the vector twice or asking frequently for nonemptyness.
      std::size_t lowsize = placeLower.size();
      std::size_t upsize = placeUpper.size();
      if (lowsize>0 && upsize>0)
      {
        std::size_t low=0; // watches the next lower bound to be corrected
        std::size_t up=0; // watches the next upper bound
        for (auto& block : native(solution))
        {
          // Firstly check whether it is one of the bounded variables.
          if (i==placeLower[low])
          {
            if (block[0]<boundLower[low])
              block[0] = boundLower[low];
            low = (low+1) % lowsize;
          }
          if (i==placeUpper[up])
          {
            if (block[0]>boundUpper[up])
              block[0] = boundUpper[up];
            up = (up+1) % upsize;
          }
          i = (i+1) % systemsize;
        }
      }
      else
      {
        if (lowsize>0)
        {
          std::size_t low=0; // watches the next lower bound to be corrected
          for (auto& block : native(solution))
          {
            if (i==placeLower[low])
            {
              if (block[0]<boundLower[low])
                block[0] = boundLower[low];
              low = (low+1) % lowsize;
            }
            i = (i+1) % systemsize;
          }
        }
        if (upsize>0)
        {
          std::size_t up=0; // watches the next lower bound to be corrected
          for (auto& block : native(solution))
          {
            if (i==placeUpper[up])
            {
              if (block[0]>boundUpper[up])
                block[0] = boundUpper[up];
              up = (up+1) % upsize;
            }
            i = (i+1) % systemsize;
          }
        }
      }
    } // end operator()

  private:
    const ParameterClass& pc;
  };

  template<typename Newton>
  class CorrectSolutionLexicographic
  {
  public:
    using ParameterClass = BoundedLineSearchParametersInterface<Newton>;
    using Domain = typename Newton::Domain;

    explicit CorrectSolutionLexicographic(const ParameterClass& pc_)
      : pc(pc_)
    {}

    void operator() (Domain& solution)
    {
      auto& placeLower = pc.getPlaceLower();
      auto& placeUpper = pc.getPlaceUpper();
      auto& boundLower = pc.getBoundLower();
      auto& boundUpper = pc.getBoundUpper();
      using Dune::PDELab::Backend::native;

      std::size_t sol_size = native(solution).size()/pc.getSystemSize();
      for (std::size_t i=0; i<placeLower.size(); ++i)
        for (std::size_t j=sol_size*placeLower[i]; j<sol_size*(placeLower[i]+1); ++j)
          if (native(solution)[j] < boundLower[i])
            native(solution)[j] = boundLower[i];
      for (std::size_t i=0; i<placeUpper.size(); ++i)
        for (std::size_t j=sol_size*placeUpper[i]; j<sol_size*(placeUpper[i]+1); ++j)
          if (native(solution)[j] > boundUpper[i])
            native(solution)[j] = boundUpper[i];
    }

  private:
    const ParameterClass& pc;
  };

  template<typename Newton>
  class CorrectSolutionSingle
  {
  public:
    using ParameterClass = BoundedLineSearchParametersInterface<Newton>;
    using Domain = typename Newton::Domain;

    explicit CorrectSolutionSingle(const ParameterClass& pc_)
      : pc(pc_)
    {}

    void operator() (Domain& solution)
    {
      auto& placeLower = pc.getPlaceLower();
      auto& placeUpper = pc.getPlaceUpper();
      auto& boundLower = pc.getBoundLower();
      auto& boundUpper = pc.getBoundUpper();

      // avoid traversing the vector twice (or multiple size()>0 checks)
      if (placeUpper.size()>0 && placeLower.size()>0)
        for (auto& v : solution)
        {
          if (v>boundUpper[0])
            v = boundUpper[0];
          if (v<boundLower[0])
            v = boundLower[0];
        }
      else
      {
        if (placeUpper.size()>0)
          for (auto& v : solution)
            if (v>boundUpper[0])
              v = boundUpper[0];
        if (placeLower.size()>0)
          for (auto& v : solution)
            if (v<boundLower[0])
              v = boundLower[0];
      }
    }

  private:
    const ParameterClass& pc;
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

  virtual void printParameters() const override
  {
    param.printParameters();
  }

private:
  Newton& _newton;
  using CorSol = typename Impl::SetCorrectSolution<Newton>::type;
  CorSol corSol;
  ParameterClass param;
};

/** \brief bounded Hackbusch-Reusken line search
 *
 * Hackbusch-Reusken line search with an extra feature.
 * Allows setting bounds on values which will not be exceeded,
 * after performing line search, solution vector is corrected to stay inside limits.
 *
 * If the parameter AcceptBest is set through the setParameters
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
   * - MaxIterations: Maximum number of line search iterations.
   *
   * - DampingFactor: Multiplier to line search parameter after each iteration.
   *
   * - AcceptBest: Accept the best line search parameter if
   *   there was any improvement, even if the convergence criterion was not
   *   reached.
   */
  virtual void setParameters(const ParameterTree& parameterTree) override
  {
    _lineSearchMaxIterations = parameterTree.get<unsigned int>("MaxIterations",
                                                               _lineSearchMaxIterations);
    _lineSearchDampingFactor = parameterTree.get<Real>("DampingFactor",
                                                       _lineSearchDampingFactor);
    _acceptBest = parameterTree.get<bool>("AcceptBest",
                                          _acceptBest);
    param.setBoundedParameters(parameterTree);
  }

  virtual void printParameters() const override
  {
    std::cout << "LineSearch.MaxIterations.. " << _lineSearchMaxIterations << std::endl;
    std::cout << "LineSearch.DampingFactor.. " << _lineSearchDampingFactor << std::endl;
    std::cout << "LineSearch.AcceptBest..... " << _acceptBest << std::endl;
    param.printParameters();
  }

private:
  Newton& _newton;
  std::shared_ptr<Domain> _previousSolution;
  using CorSol = typename Impl::SetCorrectSolution<Newton>::type;
  CorSol corSol;

  // Line search parameters
  ParameterClass param;
  unsigned int _lineSearchMaxIterations = 10;
  Real _lineSearchDampingFactor = 0.5;
  bool _acceptBest = false;
};
#endif // #ifdef DUNE_PDELAB_SOLVER_LINESEARCH_BOUNDED_READ
#endif // #ifdef DUNE_PDELAB_SOLVER_LINESEARCH_BOUNDED_ENABLED
