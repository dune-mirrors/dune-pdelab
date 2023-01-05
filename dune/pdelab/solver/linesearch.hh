// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_SOLVER_LINESEARCH_HH
#define DUNE_PDELAB_SOLVER_LINESEARCH_HH

#include <dune/pdelab/solver/newtonerrors.hh>


namespace Dune::PDELab
{

  //! Abstract base class describing the line search interface
  template <typename Domain>
  class LineSearchInterface
  {
  public:
    //! Every abstract base class should have a virtual destructor
    virtual ~LineSearchInterface () {}

    //! Do line search
    virtual void lineSearch(Domain&, const Domain&) = 0;

    //! Set parameters
    virtual void setParameters(const ParameterTree&) = 0;

    //! Print paramters
    virtual void printParameters() const
    {
      std::cout << "NewtonMethod::_lineSearch->printParameters() is not implemented." << std::endl;
    }
  };


  //! Class for simply updating the solution without line search
  template <typename Solver>
  class LineSearchNone : public LineSearchInterface<typename Solver::Domain>
  {
  public:
    using Domain = typename Solver::Domain;
    using Real = typename Solver::Real;

    LineSearchNone(Solver& solver) : _solver(solver) {}

    //! Do line search (in this case just update the solution)
    virtual void lineSearch(Domain& solution, const Domain& correction) override
    {
      solution.axpy(-1.0, correction);
      _solver.updateDefect(solution);
    }

    virtual void setParameters(const ParameterTree&) override {}

    // print line search type
    virtual void printParameters() const override
    {
      std::cout << "LineSearch.Type........... None" << std::endl;
    }

  private:
    Solver& _solver;
  };


  /** \brief Hackbusch-Reusken line search
   *
   * If the parameter line_search_accept_best is set through the setParameters
   * method this line search will simply return the best result even if it did
   * not converge.
   */
  template <typename Solver>
  class LineSearchHackbuschReusken : public LineSearchInterface<typename Solver::Domain>
  {
  public:
    using Domain = typename Solver::Domain;
    using Real = typename Solver::Real;

    LineSearchHackbuschReusken(Solver& solver, bool forceAcceptBest = false) :
      _solver(solver), _forceAcceptBest(forceAcceptBest) {}

    //! Do line search
    virtual void lineSearch(Domain& solution, const Domain& correction) override
    {
      if ((_solver.result().defect < _solver.getAbsoluteLimit())){
        solution.axpy(-1.0, correction);
        _solver.updateDefect(solution);
        return;
      }

      auto verbosity = _solver.getVerbosityLevel();

      if (verbosity >= 4)
        std::cout << "      Performing line search..." << std::endl;
      Real lambda = 1.0;
      Real bestLambda = 0.0;
      Real bestDefect = _solver.result().defect;
      Real previousDefect = _solver.result().defect;
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
        _solver.updateDefect(solution);
        if (verbosity >= 4){
          if (not std::isfinite(_solver.result().defect))
            std::cout << "          NaNs detected" << std::endl;
        }       // ignore NaNs and try again with lower lambda

        if (_solver.result().defect <= (1.0 - lambda/4) * previousDefect){
          if (verbosity >= 4)
            std::cout << "          line search converged" << std::endl;
          converged = true;
          break;
        }

        if (_solver.result().defect < bestDefect){
          bestDefect = _solver.result().defect;
          bestLambda = lambda;
        }

        lambda *= _lineSearchDampingFactor;
        solution = *_previousSolution;
      }

      if (not converged){
        if (verbosity >= 4)
          std::cout << "          max line search iterations exceeded" << std::endl;

        if (not (_acceptBest or _forceAcceptBest)){
          solution = *_previousSolution;
          _solver.updateDefect(solution);
          DUNE_THROW( LineSearchError,
                     "LineSearch::lineSearch(): line search failed, "
                     "max iteration count reached, "
                     "defect did not improve enough");
        }
        else{
          if (bestLambda == 0.0){
            solution = *_previousSolution;
            _solver.updateDefect(solution);
            DUNE_THROW(LineSearchError,
                       "LineSearch::lineSearch(): line search failed, "
                       "max iteration count reached, "
                       "defect did not improve in any of the iterations");
          }
          if (bestLambda != lambda){
            solution = *_previousSolution;
            solution.axpy(-bestLambda, correction);
            _solver.updateDefect(solution);
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
      _acceptBest = parameterTree.get<bool>("AcceptBest", _acceptBest);
    }

    virtual void printParameters() const override
    {
      std::cout << "LineSearch.Type........... Hackbusch-Reusken" << std::endl;
      std::cout << "LineSearch.MaxIterations.. " << _lineSearchMaxIterations << std::endl;
      std::cout << "LineSearch.DampingFactor.. " << _lineSearchDampingFactor << std::endl;
      std::cout << "LineSearch.AcceptBest..... " << (_acceptBest or _forceAcceptBest) << std::endl;
    }

  private:
    Solver& _solver;
    std::shared_ptr<Domain> _previousSolution;

    // Line search parameters
    unsigned int _lineSearchMaxIterations = 10;
    Real _lineSearchDampingFactor = 0.5;
    bool _acceptBest = false;
    bool _forceAcceptBest;
  };

  namespace Impl
  {
    // Default projection for LineSearch with projection step.
    template<typename Domain>
    struct ThrowingProjection
    {
      void operator() (Domain& solution) const
      {
        DUNE_THROW(LineSearchError,"Using projecting LineSearch without user-defined projection.");
      }
    };
  }

  /** \brief Projected line search
   *
   * Projects the solution to the feasible region. The projection is
   * user-defined lambda and passed as a constructor argument.
   *
   * Can be passed to the solver via setLineSearch method.
   */
  template <typename Solver, typename Projection>
  class LineSearchProjectedNone : public LineSearchInterface<typename Solver::Domain>
  {
  public:
    using Domain = typename Solver::Domain;
    using Real = typename Solver::Real;

    LineSearchProjectedNone(Solver& solver, const Projection& projection) : _solver(solver), _projection(projection) {}

    //! Do line search (in this case just update the solution)
    virtual void lineSearch(Domain& solution, const Domain& correction) override
    {
      solution.axpy(-1.0, correction);
      _projection(solution); // project the solution candidate to the feasible region
      _solver.updateDefect(solution);
    }

    virtual void setParameters(const ParameterTree&) override {}

    // print line search type
    virtual void printParameters() const override
    {
      std::cout << "LineSearch.Type........... ProjectedNone" << std::endl;
    }

  private:
    Solver& _solver;
    const Projection& _projection;
  };

  //! Class for mixing Newton and Steepest descent direction
  template <typename Solver>
  class LineSearchWithGradientDescent : public LineSearchInterface<typename Solver::Domain>
  {
  public:
    using Domain = typename Solver::Domain;
    using Real = typename Solver::Real;

    LineSearchWithGradientDescent(Solver& solver) : _solver(solver)
    {}

    //! Update is a mix of Newton direction and Gradient Descent direction
    virtual void lineSearch(Domain& solution, const Domain& correction) override
    {
      if (not _ga)
        _ga = std::make_shared<Domain>(correction);

      try
      {
        uint itnr = _solver.result().iterations;
        if (itnr==0)
          _alphawarp = _alpha;
        else
          _alphawarp = std::min(1., _alphawarp*_arise);
      }
      catch (...)
      {
        _alphawarp = _alpha;
      }
      Backend::native(*(_solver.getJacobian())).mtv(Backend::native(_solver.getResidual()),Backend::native(*_ga)); // J^T f(x), gradient ascent
      solution.axpy(-_alphawarp, correction);
      solution.axpy(-(1.-_alphawarp)*std::sqrt(_alphawarp),*_ga);
      _solver.updateDefect(solution);
    }

    virtual void setParameters(const ParameterTree& parameterTree) override
    {
      _alpha = parameterTree.get<Real>("Alpha", _alpha);
      _arise = parameterTree.get<Real>("AlphaRise", _arise);
    }

    void setAlpha (const Real& alpha)
    {
      _alpha = alpha;
    }

    void setRise (const Real& arise)
    {
      _arise = arise;
    }

    // print line search type
    virtual void printParameters() const override
    {
      std::cout << "LineSearch.Type........... LineSearchWithGradientDescent" << std::endl;
      std::cout << "LineSearch.Alpha.......... " << _alpha << std::endl;
      std::cout << "LineSearch.AlphaRise...... " << _arise << std::endl;
    }

  private:
    Solver& _solver;
    std::shared_ptr<Domain> _ga; // gradient ascent direction
    Real _alpha=0.5; // starting value of mixing
    Real _arise=1.1; // increasing rate, _alphawarp is at maximum 1
    Real _alphawarp = _alpha;  // Mixing parameter, 1=Newton directon, 0 goes to damped gradient descent
  };


  //! Implement Newton Cauchy Dogleg method via line search method
  template<typename Solver>
  class DoglegLineSearch : public LineSearchInterface<typename Solver::Domain>
  {
  public:
    using Domain = typename Solver::Domain;
    using Real = typename Solver::Real;

    DoglegLineSearch(Solver& solver) : _solver(solver)
    {}

    //! Update is a mix of Newton direction and Gradient Descent direction
    virtual void lineSearch(Domain& solution, const Domain& correction) override
    {
      std::shared_ptr<const Domain> pcorr(&correction,[](const Domain*){});
      std::shared_ptr<Domain> _ga; // gradient ascent direction
      std::shared_ptr<Domain> _jv; // temporary, mostly Jacobian*vector

      if (_solver.getVerbosityLevel() >= 4)
        std::cout << "      Performing line search with trust region " << _radius << std::endl;

      const Real& defectold = _solver.result().defect;

      // gradient ascent direction
      if (not _ga)
        _ga = std::make_shared<Domain>(correction);
      Backend::native(*(_solver.getJacobian())).mtv(Backend::native(_solver.getResidual()),Backend::native(*_ga)); // J^T f(x), gradient ascent

      while(true) // reduce radius until rho>nu1
      {
        if (correction.two_norm()<=_radius)
        {
          // Newton direction if Newton guess is inside TrustRegion
          if (_solver.getVerbosityLevel()>=4)
            std::cout << "          Newton direction chosen" << std::endl;
          _pu.set(-1.,pcorr);
        }
        else
        {
          Real ganorm = std::max(_ga->two_norm(),1e-100);
          if (not _jv)
            _jv = std::make_shared<Domain>(*_ga);
          Backend::native(*(_solver.getJacobian())).mtv(Backend::native(*_ga),Backend::native(*_jv)); // J J^T f(x),
          Real gaBga = _jv->two_norm2(); // approximate Hessian by J^TJ
          gaBga = std::max(gaBga,1e-100);
          Real alpha = std::min(_radius/ganorm, ganorm*ganorm/gaBga);
          if (ganorm*alpha>=_radius)
          {
            // damped Cauchy direction if it is big enough
            if (_solver.getVerbosityLevel()>=4)
              std::cout << "          Cauchy direction chosen" << std::endl;
            _pu.set(-alpha,_ga);
          }
          else
          {
            // dogleg direction
            // Point on the line between Newton (N) and Cauchy (C) point at distance TrustRegion (r) from Starting point (S)
            // S-N is correction, S-C is *_ga
            // we search for tau such that \|C-S+tau(N-C)\|=r
            // => quadratic eq. tau^2 NC^2 + 2tau <NC,CS> + CS^2-r^2=0
            // knowing \|C-S\|<r this has one negative root and one in (0,1) that we seek
            // N-C=(S-C)-(S-N) is not constructed.

            Real cs2 = _ga->two_norm2();
            Real sn_sc = Backend::native(correction) * Backend::native(*_ga); // <S-N,S-C>
            Real nc2 = correction.two_norm2()-2.*sn_sc+cs2; // squared norm of N-C
            Real sp = cs2-sn_sc; // <N-C,S-C>
            Real tau = (-sp+std::sqrt(sp*sp-4.*nc2+(cs2-_radius*_radius)))/(2.*nc2);
            if (_solver.getVerbosityLevel()>=4)
              std::cout << "          Dogleg direction chosen" << std::endl;
            _pu.set(tau-1.,_ga); // Cauchy point
            _pu.add(-tau,pcorr); // Dogleg point
          }
        }

        // postprocess, calculate rho
        if (not _jv)
          _jv = std::make_shared<Domain>(solution);
        Backend::native(*_jv)=Backend::native(solution);
        _jv->axpy(-1.,correction);
        _solver.updateDefect(*_jv);
        Real defNewton = _solver.result().defect;
        Backend::native(*(_solver.getJacobian())).mv(Backend::native(correction),Backend::native(*_jv)); // J*p, p is Newton
        Real pJJp = _jv->two_norm2();
        Real resJp = Backend::native(_solver.getResidual()) * Backend::native(*_jv);
        Real rho = (defectold*defectold - defNewton*defNewton)/(resJp-0.5*pJJp);

        // use rho to change radius and decide to continue iterating
        if (rho>=_nu3)
          _radius *= _t2;
        else if (rho<=_nu2)
          _radius *= _t1;
        if (rho>_nu1)
          break;
        else if (_radius<_minRadius)
        {
          _pu.clear();
          DUNE_THROW(LineSearchError,"DoglegLineSearch did not converge, TrustRegion reduced below minimum");
        }
        else
          if (_solver.getVerbosityLevel()>=4)
            std::cout << "        Iteration failed, reducing trust region radius" << std::endl;
      } // end while

      if (_solver.getVerbosityLevel()>=4)
        std::cout << "        Applying chosen direction" << std::endl;
      _pu.apply(solution);
      _solver.updateDefect(solution);
    }

    virtual void setParameters(const ParameterTree& parameterTree) override
    {
      _radius = parameterTree.get<Real>("TrustRegion", _radius);
      _nu1 = parameterTree.get<Real>("Nu1", _nu1);
      _nu2 = parameterTree.get<Real>("Nu2", _nu2);
      _nu3 = parameterTree.get<Real>("Nu3", _nu3);
      _t1 = parameterTree.get<Real>("T1", _t1);
      _t2 = parameterTree.get<Real>("T2", _t2);
    }

    void setTrustRegion(const Real& radius)
    {
      _radius=radius;
    }
    void setTrustRegion(const Domain& solution)
    {
      _radius = _iniRadius*solution.two_norm();
    }

    // print line search type
    virtual void printParameters() const override
    {
      std::cout << "LineSearch.Type........... DoglegLineSearch" << std::endl;
      std::cout << "LineSearch.IniTrustRegion. " << _iniRadius << std::endl;
      std::cout << "LineSearch.MinTrustRegion. " << _minRadius << std::endl;
      std::cout << "LineSearch.Nu1............ " << _nu1 << std::endl;
      std::cout << "LineSearch.Nu2............ " << _nu2 << std::endl;
      std::cout << "LineSearch.Nu3............ " << _nu3 << std::endl;
      std::cout << "LineSearch.T1............. " << _t1 << std::endl;
      std::cout << "LineSearch.T2............. " << _t2 << std::endl;
    }

  private:
    //! Structure holding proposed updates
    // Avoid storing extra result vector, Dogleg makes two steps
    struct ProposedUpdate
    {
      //! set an update, clear
      void set(const Real& coef, const std::shared_ptr<const Domain>& pu)
      {
        clear();
        add(coef,pu);
      }

      //! add another update
      void add(const Real& coef, const std::shared_ptr<const Domain>& pu)
      {
        _pu.push_back(std::pair(coef,pu));
      }

      //! apply proposed updates and clear cache
      void apply(Domain& domain)
      {
        for (const auto& v : _pu)
          domain.axpy(v.first,*(v.second));
        _pu.clear();
      }

      void clear()
      {
        _pu.clear();
      }

    private:
      std::vector<std::pair<Real,std::shared_ptr<const Domain>>> _pu;
    };

    Solver& _solver;
    Real _radius=0.; // trust region _radius, solution-dependent
    Real _iniRadius=0.2;
    Real _minRadius=1e-6;
    Real _nu1=0.001, _nu2=0.25, _nu3=0.75;
    Real _t1=0.25, _t2=2.0;
    ProposedUpdate _pu;
  };

  //! Flags for different line search strategies
  enum class LineSearchStrategy
  {
    noLineSearch,
    hackbuschReusken,
    hackbuschReuskenAcceptBest,
    projectedNoLineSearch,
    lineSearchWithGradientDescent,
    doglegLineSearch
  };

  // we put this into an emty namespace, so that we don't violate the one-definition-rule
  namespace {
    /** \brief Get a LineSearchStrategy from a string identifier
     *
     * \param name Identifier used to pick LineSearchStrategy
     *
     * Possible values for name: "noLineSearch", "hackbuschReusken"
     */
    inline
    LineSearchStrategy lineSearchStrategyFromString (const std::string& name)
    {
      if (name == "noLineSearch")
        return LineSearchStrategy::noLineSearch;
      if (name == "hackbuschReusken")
        return LineSearchStrategy::hackbuschReusken;
      if (name == "hackbuschReuskenAcceptBest")
        return LineSearchStrategy::hackbuschReuskenAcceptBest;
      if (name == "lineSearchWithGradientDescent")
        return LineSearchStrategy::lineSearchWithGradientDescent;
      if (name == "projectedNoLineSearch")
        return LineSearchStrategy::projectedNoLineSearch;
      if (name == "doglegLineSearch")
        return LineSearchStrategy::doglegLineSearch;
      DUNE_THROW(Exception,"Unkown line search strategy: " << name);
    }
  }


  /** \brief fectory function to create an instace of a line-search
   *
   * \tparam Solver A solver
   *
   * \param solver Solver object

   * \param name Identifier to choose line search. Possible values:
   * - "noLineSearch": Return pointer to LineSearchNone
   * - "hackbuschReusken": Return pointer to LineSearchHackbuschReusken
   */
  template <typename Solver>
  std::shared_ptr<LineSearchInterface<typename Solver::Domain>>
  createLineSearch(Solver& solver, LineSearchStrategy strategy)
  {
    if (strategy == LineSearchStrategy::noLineSearch){
      auto lineSearch = std::make_shared<LineSearchNone<Solver>> (solver);
      return lineSearch;
    }
    if (strategy == LineSearchStrategy::hackbuschReusken){
      auto lineSearch = std::make_shared<LineSearchHackbuschReusken<Solver>> (solver);
      return lineSearch;
    }
    if (strategy == LineSearchStrategy::hackbuschReuskenAcceptBest){
      auto lineSearch = std::make_shared<LineSearchHackbuschReusken<Solver>> (solver, true);
      std::cout << "Warning: linesearch hackbuschReuskenAcceptBest is deprecated and will be removed after PDELab 2.7.\n"
                << "         Please use 'hackbuschReusken' and add the parameter 'LineSearchAcceptBest : true'";
      return lineSearch;
    }
    if (strategy == LineSearchStrategy::lineSearchWithGradientDescent){
      auto lineSearch = std::make_shared<LineSearchWithGradientDescent<Solver>> (solver);
      return lineSearch;
    }
    if (strategy == LineSearchStrategy::projectedNoLineSearch){
      // std::cout << "Projected line search requires loading via setLineSearch method. Using LineSearchNone till then." << std::endl;
      // auto lineSearch = std::make_shared<LineSearchNone<Solver>> (solver);
      if (solver.getVerbosityLevel()>=4)
        std::cout << "Setting default projectedNoLineSearch with ThrowingProjection step. "
                  << "Create custom LineSearch and add it via NewtonMethod's setLineSearch method."<< std::endl;
      using Projection = Impl::ThrowingProjection<typename Solver::Domain>;
      Projection throwproj;
      auto lineSearch = std::make_shared<LineSearchProjectedNone<Solver,Projection> > (solver, throwproj);
      return lineSearch;
    }
    if (strategy == LineSearchStrategy::doglegLineSearch){
      auto lineSearch = std::make_shared<DoglegLineSearch<Solver>> (solver);
      return lineSearch;
    }
    DUNE_THROW(Exception,"Unkown line search strategy");
  }

} // namespace Dune::PDELab

#endif
