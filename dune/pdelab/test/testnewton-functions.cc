// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include "dune/pdelab.hh"

#include "nonlinearpoissonfem.hh"

// Copy of testnewton.cc using a dune-functions basis instead of a PDELab grid
// function space.

// squared L2 difference of a PDELab grid function and a dune-functions grid function
template<typename GridView, typename F, typename G>
double l2DifferenceSquared(const GridView& gridView, const F& f, const G& g, int order)
{
  double result = 0.0;
  auto localG = localFunction(g);
  for (const auto& element : elements(gridView))
  {
    localG.bind(element);
    const auto geometry = element.geometry();
    for (const auto& qp : Dune::QuadratureRules<double,GridView::dimension>::rule(element.type(), order))
    {
      typename F::Traits::RangeType fValue;
      f.evaluate(element, qp.position(), fValue);
      const auto difference = fValue - localG(qp.position());
      result += difference*difference * qp.weight() * geometry.integrationElement(qp.position());
    }
  }
  return result;
}

template<typename Number>
class NonlinearPoissonProblem
{
  Number eta;
public:
  typedef Number value_type;

  //! Constructor without arg sets nonlinear term to zero
  NonlinearPoissonProblem () : eta(0.0) {}

  //! Constructor takes eta parameter
  NonlinearPoissonProblem (const Number& eta_) : eta(eta_) {}

  //! nonlinearity
  Number q (Number u) const
  {
    return eta*u*u;
  }

  //! derivative of nonlinearity
  Number qprime (Number u) const
  {
    return 2*eta*u;
  }

  //! right hand side
  template<typename E, typename X>
  Number f (const E& e, const X& x) const
  {
    auto global = e.geometry().global(x);
    return -2.0*x.size() + eta*global.two_norm2()*global.two_norm2();
  }

  //! boundary condition type function (true = Dirichlet)
  template<typename I, typename X>
  bool b (const I& i, const X& x) const
  {
    return true;
  }

  //! Dirichlet extension
  template<typename E, typename X>
  Number g (const E& e, const X& x) const
  {
    auto global = e.geometry().global(x);
    return global.two_norm2();
  }

  //! Neumann boundary condition
  template<typename I, typename X>
  Number j (const I& i, const X& x) const
  {
    return 0.0;
  }
};


int main(int argc, char** argv)
{
  try{
    // Maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // Create grid
    const int dim = 2;
    Dune::FieldVector<double,dim> lowerleft(0.0);
    Dune::FieldVector<double,dim> upperright(1.0);
    auto cells = Dune::filledArray<dim,unsigned int>(4);
    using Grid = Dune::YaspGrid<dim>;
    auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerleft, upperright, cells);
    grid -> globalRefine(3);
    using GridView = Grid::LeafGridView;
    GridView gridView = grid -> leafGridView();

    // Basis
    using DomainField = GridView::Grid::ctype;
    using RangeType = double;
    const int degree = 2;
    auto basis = Dune::Functions::BasisFactory::makeBasis(gridView, Dune::Functions::BasisFactory::lagrange<degree>());
    using Basis = decltype(basis);
    using FiniteElementMap = Dune::PDELab::FiniteElementMapFor<typename Basis::LocalView::Tree::FiniteElement>;

    // Solution vector
    using CoefficientVector = Dune::PDELab::Backend::Vector<Basis, DomainField>;
    CoefficientVector coefficientVector(basis);

    // Discrete grid function of solution vetor
    auto discreteGridFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<RangeType>(
      basis, Dune::PDELab::Backend::native(coefficientVector));

    // Local operator (problem is nonlinear and problem class depends on the solution)
    using Problem = NonlinearPoissonProblem<RangeType>;
    Problem problem(2.0);
    using LocalOperator = NonlinearPoissonFEM<Problem, FiniteElementMap>;
    LocalOperator localOperator(problem);

    // Create constraints map
    using ConstraintsContainer = Dune::PDELab::BasisConstraintsContainer<Basis, RangeType>;
    ConstraintsContainer constraintsContainer;
    auto blambda = [&](const auto& i, const auto& x){return problem.b(i,x);};
    auto bctype = Dune::PDELab::makeBoundaryConditionFromCallable(gridView, blambda);
    Dune::PDELab::addBoundaryDirichletConstraints(basis, constraintsContainer, bctype);

    // Print number of DOFs
    std::cout << "basis with " << basis.size() << " dofs generated  "<< std::endl;
    std::cout << "cc with " << constraintsContainer.size() << " dofs generated  "<< std::endl;

    // Grid operator
    using MatrixBackend = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    const int dofestimate = 4 * basis.localView().maxSize();
    MatrixBackend matrixBackend(dofestimate);
    using GridOperator = Dune::PDELab::GridOperator<Basis,
                                                    Basis,
                                                    LocalOperator,
                                                    MatrixBackend,
                                                    DomainField,
                                                    RangeType,
                                                    RangeType,
                                                    ConstraintsContainer,
                                                    ConstraintsContainer>;
    GridOperator gridOperator(basis,
                              constraintsContainer,
                              basis,
                              constraintsContainer,
                              localOperator,
                              matrixBackend);

    // Get grid function from Dirichlet boundary condition
    // (dune-functions interpolates functions of the global coordinate)
    auto glambda = [&](const auto& e, const auto& x){return problem.g(e,x);};
    auto boundaryCondition = Dune::PDELab::makeGridFunctionFromCallable(gridView, glambda);
    Dune::Functions::interpolate(basis, Dune::PDELab::Backend::native(coefficientVector),
                                 [](const auto& x){ return x.two_norm2(); });
    Dune::PDELab::set_nonconstrained_dofs(constraintsContainer, 0.0, coefficientVector);

    // Create Netwon solver
    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_SuperLU;
    LinearSolver linearSolver(false);
    using Solver = Dune::PDELab::NewtonMethod<GridOperator, LinearSolver>;
    Solver solver(gridOperator, linearSolver);
    solver.setReduction(1e-7);

    // Use some nonsense parameters to ensure that setting them explicitly works
    Dune::ParameterTree ptree;
    ptree["VerbosityLevel"] = "4";
    ptree["UseMaxNorm"] = "1";
    ptree["LineSearch.DampingFactor"] = "0.3";
    ptree["Terminate.MaxIterations"] = "0";
    ptree["LineSearch.MaxIterations"] = "0";

    // Set some parameters without loading parameter tree from ini file
    // (just to show that it works)
    solver.setParameters(ptree);

    // Retrieve the terminate interface and set parameters
    auto terminate = solver.getTerminate();
    terminate->setParameters(ptree.sub("Terminate"));

    // Create new terminate with other values and insert it into solver
    auto terminate_new =
      std::make_shared<Dune::PDELab::DefaultTerminate<Solver>>(solver);
    terminate_new->setForceIteration(true);
    terminate_new->setMaxIterations(39);
    solver.setTerminate(terminate_new);

    // Retrieve line search interface and set parameters
    auto line_search = solver.getLineSearch();
    ptree["LineSearch.MaxIterations"] = "20";
    line_search->setParameters(ptree.sub("LineSearch"));

    // Solve PDE
    solver.apply(coefficientVector);

    // Visualization
    using VTKWriter = Dune::SubsamplingVTKWriter<GridView>;
    Dune::RefinementIntervals subint(2);
    VTKWriter vtkwriter(gridView, subint);
    std::string vtkfile("testnewton-functions");
    vtkwriter.addVertexData(discreteGridFunction, Dune::VTK::FieldInfo("numerical_solution", Dune::VTK::FieldInfo::Type::scalar, 1));
    vtkwriter.write(vtkfile, Dune::VTK::ascii);

    // Calculate error
    //
    // Note: The problem is set up in such a way that the Dirichlet boundary
    // condition is the exact solution of the problem.
    const double error = l2DifferenceSquared(gridView, boundaryCondition, discreteGridFunction, 10);
    std::cout << "l2errorsquared: " << error << std::endl;

    // Let the test fail if the error is too large
    bool testfail(false);
    using std::abs;
    using std::isnan;
    if (isnan(error) or abs(error)>1e-7)
      testfail = true;
    return testfail;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
        return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
        return 1;
  }
}
