//===========================================================================
// This is a system test solving the heat equation and comparing against the
// exact solution. The problem is constructed in such a way that the solution
// of the heatequation does not change in time and is given by the dirichlet
// boundary condition.
//
// Copy of testinstationary.cc using a dune-functions basis instead of a
// PDELab grid function space.
// ==========================================================================

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include "dune/pdelab.hh"

template <typename GO>
void testGridOperatorInterface(const GO& go)
{
  typename GO::Traits::Domain u(go.trialGridFunctionSpace(), 0.0);
  typename GO::Traits::Range r(go.trialGridFunctionSpace());
  typename GO::Traits::Jacobian jac(go);
  go.residual(u, r);
  go.jacobian(u, jac);
  go.jacobian_apply(u, r);
}


template <class GridView, class RangeType>
class PoissonProblem
  : public Dune::PDELab::ConvectionDiffusionModelProblem<GridView, RangeType>
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<RangeType>
{
public:
  // Source term
  template<typename Element, typename Coord>
  auto f(const Element& element, const Coord& x) const
  {
    auto global = element.geometry().global(x);
    auto c = (0.5-global[0])*(0.5-global[0]) + (0.5-global[1])*(0.5-global[1]);
    using std::exp;
    auto g = exp(-1.0*c);
    auto f = 4*(1.0-c)*g;
    return f;
  }

  // Boundary condition type
  template<typename Element, typename Coord>
  auto bctype(const Element& element, const Coord& x) const
  {
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  // Dirichlet extension
  template<typename Element, typename Coord>
  RangeType g (const Element& element, const Coord& x) const
  {
    auto global = element.geometry().global(x);
    auto c = (0.5-global[0])*(0.5-global[0]) + (0.5-global[1])*(0.5-global[1]);
    using std::exp;
    auto g = exp(-1.0*c);
    return g;
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

    // Local operator for spatial discretization
    using Problem = PoissonProblem<GridView, RangeType>;
    Problem problem;
    using LocalOperator = Dune::PDELab::ConvectionDiffusionFEM<Problem,FiniteElementMap>;
    LocalOperator localOperator(problem);

    // Local operator for time discretization
    using LocalOperatorTime = Dune::PDELab::L2;
    LocalOperatorTime localOperatorTime;

    // Create constraints map
    using ConstraintsContainer = Dune::PDELab::ConstraintsContainer<Basis, RangeType>;
    ConstraintsContainer constraintsContainer;
    constraintsContainer.clear();
    Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> bctype(gridView, problem);
    Dune::PDELab::constraints(bctype, basis, constraintsContainer);

    // Grid operator for spatial discretization
    using MatrixBackend = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    const int dofestimate = 4 * Dune::PDELab::maxLocalSize(basis);
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

    // Grid operator for time discretization
    using GridOperatorTime = Dune::PDELab::GridOperator<Basis,
                                                        Basis,
                                                        LocalOperatorTime,
                                                        MatrixBackend,
                                                        DomainField,
                                                        RangeType,
                                                        RangeType,
                                                        ConstraintsContainer,
                                                        ConstraintsContainer>;
    GridOperatorTime gridOperatorTime(basis,
                                      constraintsContainer,
                                      basis,
                                      constraintsContainer,
                                      localOperatorTime,
                                      matrixBackend);

    // Combined grid operator
    using InstationaryGridOperator = Dune::PDELab::OneStepGridOperator<GridOperator, GridOperatorTime>;
    InstationaryGridOperator instationaryGridOperator(gridOperator, gridOperatorTime);

    // Solution vector
    using CoefficientVector = Dune::PDELab::Backend::Vector<Basis, DomainField>;
    CoefficientVector coefficientVector(basis);
    using DirichletExtension = Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem>;
    DirichletExtension dirichletExtension(gridView, problem);
    Dune::PDELab::interpolate(dirichletExtension, basis, coefficientVector);
    // Dune::PDELab::set_nonconstrained_dofs(constraintsContainer, 0.0, coefficientVector);

    // Solver
    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_SuperLU;
    LinearSolver linearSolver(false);
    using Solver = Dune::PDELab::NewtonMethod<InstationaryGridOperator, LinearSolver>;
    Solver solver(instationaryGridOperator, linearSolver);

    // Time stepping method
    using TimeSteppingMethod = Dune::PDELab::OneStepThetaParameter<RangeType>;
    TimeSteppingMethod timeSteppingMethod(1.0);
    using OneStepMethod = Dune::PDELab::OneStepMethod<RangeType,
                                                      InstationaryGridOperator,
                                                      Solver,
                                                      CoefficientVector,
                                                      CoefficientVector>;
    OneStepMethod oneStepMethod(timeSteppingMethod, instationaryGridOperator, solver);
    double time = 0.0;

    // Visualization
    using VTKWriter = Dune::SubsamplingVTKWriter<GridView>;
    Dune::RefinementIntervals subint(2);
    VTKWriter vtkwriter(gridView, subint);
    std::string vtkfile("testheat-functions");
    using VTKSW = Dune::VTKSequenceWriter<GridView>;
    VTKSW vtkSequenceWriter(std::make_shared<VTKWriter>(vtkwriter), vtkfile);
    Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter, basis, coefficientVector,
                                         Dune::PDELab::vtk::defaultNameScheme());
    vtkSequenceWriter.write(time, Dune::VTK::appendedraw);

    // Time loop
    double T = 1.0;
    double dt = 0.1;
    while (time<T-1e-8){
      // assemble constraints for new time step
      problem.setTime(time+dt);
      Dune::PDELab::constraints(bctype,basis,constraintsContainer);

      // Do time step
      CoefficientVector newCoefficientVector(coefficientVector);
      oneStepMethod.apply(time,dt,coefficientVector, dirichletExtension, newCoefficientVector);

      // Accept time step
      coefficientVector = newCoefficientVector;
      time+=dt;

      // Output to VTK file
      vtkSequenceWriter.write(time,Dune::VTK::appendedraw);
    }

    // Calculate error
    //
    // Note: The problem is set up in such a way that the Dirichlet boundary
    // condition is the exact solution of the problem.
    using DiscreteGridFunction = Dune::PDELab::DiscreteGridFunction<Basis, CoefficientVector>;
    DiscreteGridFunction discreteGridFunction(basis, coefficientVector);
    using DifferenceSquaredAdapter = Dune::PDELab::DifferenceSquaredAdapter<DirichletExtension, DiscreteGridFunction>;
    DifferenceSquaredAdapter differenceSquaredAdapder(dirichletExtension, discreteGridFunction);
    DifferenceSquaredAdapter::Traits::RangeType error(0.0);
    Dune::PDELab::integrateGridFunction(differenceSquaredAdapder, error, 10);
    std::cout << "l2errorsquared matrix based: " << error << std::endl;

    //=================================
    // Test the grid operator interface
    //=================================
    testGridOperatorInterface(instationaryGridOperator);

    // Let the test fail if the error is too large
    bool testfail(false);
    using std::abs;
    using std::isnan;
    if (isnan(error[0]) or abs(error[0])>1e-7)
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
