#ifndef DUNE_PDELAB_TEST_TESTADAPTIVITY_ADAPTIVITY_FUNCTIONS_HH
#define DUNE_PDELAB_TEST_TESTADAPTIVITY_ADAPTIVITY_FUNCTIONS_HH

template<class Grid, class GV>
void adaptivity (Grid& grid, const GV& gv, int startLevel, int maxLevel)
{
  using Dune::PDELab::Backend::native;

  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make basis
  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gv, lagrange<1>());
  typedef decltype(basis) Basis;

  // <<<3>>> assemble constraints on this space
  BCTypeParam bctype; // boundary condition type
  typedef Dune::PDELab::ConstraintsContainer<Basis,Real> CC;
  CC cc;
  Dune::PDELab::constraints( bctype, basis, cc );               // assemble constraints

  // <<<4>>> make DOF vector
  typedef typename Dune::PDELab::Backend::Vector<Basis,Real> U;
  U u(basis,0.0);
  typedef BCExtension<GV,Real> G;                        // boundary value + extension
  G g(gv);
  Dune::PDELab::interpolate(g,basis,u);                  // interpolate coefficient vector

  for (int i = 0; i <= maxLevel - startLevel; i++)
  {
    std::stringstream s;
    s << i;
    std::string iter;
    s >> iter;
    std::cout << "Iteration: " << iter << "\thighest level in grid: " << grid.maxLevel() << std::endl;
    std::cout << "constrained dofs=" << cc.size()
            << " of " << Dune::PDELab::globalSize(basis) << std::endl;

    //  Make grid operator
    typedef Example02LocalOperator<BCTypeParam> LOP;       // operator including boundary
    LOP lop(bctype);
    typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
    MBE mbe(7);
    typedef Dune::PDELab::GridOperator<Basis,Basis,LOP,MBE,Real,Real,Real,CC,CC> GO;
    GO go(basis,cc,basis,cc,lop,mbe);

    // Select a linear solver backend
    typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS;
    LS ls(5000,true);

    // Select linear problem solver
    typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
    SLP slp(go,ls,u,1e-10);

    // Preparation: Define types for the computation of the error estimate eta.
    typedef decltype(makeBasis(gv, lagrange<0>())) P0Basis;
    typedef Dune::PDELab::ExampleErrorEstimator ESTLOP;
    typedef Dune::PDELab::EmptyTransformation NoTrafo;
    typedef typename Dune::PDELab::Backend::Vector<P0Basis,Real> U0;

    // <<<8>>> Solve linear problem.
    slp.apply();

    // <<<9>>> graphical output
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter,basis,u);
    vtkwriter.write("testadaptivity-functions_"+iter,Dune::VTK::appendedraw);

    // <<<10>>> compute estimated error eta
    auto p0basis = makeBasis(gv, lagrange<0>());
    ESTLOP estlop;
    typedef Dune::PDELab::GridOperator<Basis,P0Basis,ESTLOP,MBE,Real,Real,Real,NoTrafo,NoTrafo> ESTGO;
    ESTGO estgo(basis,p0basis,estlop,mbe);
    U0 eta(p0basis,0.0);
    estgo.residual(u,eta);

    for (auto& eta_e : native(eta))
      eta_e = sqrt(eta_e); // eta contains squares

    // Use eta to refine the grid following two different strategies based
    // (1) element fraction
    // (2) error fraction

    double alpha(0.4);       // refinement fraction
    double eta_alpha(0);     // refinement threshold
    double beta(0.0);        // coarsening fraction
    double eta_beta(0);      // coarsening threshold
    int verbose = 0;

    // <<<10>>> Adapt the grid locally...
    // with strategy 1:
    Dune::PDELab::element_fraction( eta, alpha, beta, eta_alpha, eta_beta, verbose );
    // or, alternatively, with strategy 2:
    //Dune::PDELab::error_fraction( eta, alpha, beta, eta_alpha, eta_beta, verbose );

    Dune::PDELab::mark_grid( grid, eta, eta_alpha, 0.0 ,0 , 100, verbose);
    Dune::PDELab::adapt_grid( grid, basis, u, 2 );

    Dune::PDELab::constraints(bctype,basis,cc);
    Dune::PDELab::interpolate(g,basis,u); // this will overwrite the solution !
  }
}

#endif // DUNE_PDELAB_TEST_TESTADAPTIVITY_ADAPTIVITY_FUNCTIONS_HH
