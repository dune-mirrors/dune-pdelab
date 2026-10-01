// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/* Stationary lid-driven cavity for the Navier-Stokes equations with
 * Taylor-Hood elements (Q2/Q1), trimmed from the driven cavity example of
 * tutorial08 in dune-pdelab-tutorials: same spaces and boundary conditions,
 * but stationary and with the TaylorHoodNavierStokes local operator of
 * dune-pdelab. One pressure DOF is fixed, so a direct solver can be used.
 *
 * Copy of testnavierstokes.cc using a dune-functions basis instead of a
 * PDELab grid function space.
 *
 * The test checks that Newton converges and compares the kinetic energy of
 * the solution with a reference value.
 */

#include "config.h"

#include <cmath>
#include <iomanip>
#include <iostream>

#include <dune/common/filledarray.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

#include <dune/pdelab.hh>

// Parameters of the driven cavity: no volume force, velocity Dirichlet conditions everywhere
template<typename GV>
class DrivenCavityProblem
{
public:
  static const bool assemble_navier = true;
  static const bool assemble_full_tensor = false;

  using Traits = Dune::PDELab::NavierStokesParameterTraits<GV,double>;

  explicit DrivenCavityProblem(double viscosity)
    : viscosity_(viscosity)
  {}

  void setTime(double) {}

  template<typename Element>
  typename Traits::VelocityRange
  f(const Element&, const typename Traits::Domain&) const
  {
    return typename Traits::VelocityRange(0.0);
  }

  template<typename EG>
  typename Traits::RangeField
  g2(const EG&, const typename Traits::Domain&) const
  {
    return 0.0;
  }

  template<typename EG>
  typename Traits::RangeField
  mu(const EG&, const typename Traits::Domain&) const
  {
    return viscosity_;
  }

  template<typename EG>
  typename Traits::RangeField
  rho(const EG&, const typename Traits::Domain&) const
  {
    return 1.0;
  }

  template<typename IG>
  typename Traits::BoundaryCondition::Type
  bctype(const IG&, const typename Traits::IntersectionDomain&) const
  {
    return Traits::BoundaryCondition::VelocityDirichlet;
  }

  template<typename IG>
  typename Traits::VelocityRange
  j(const IG&, const typename Traits::IntersectionDomain&, const typename Traits::Domain&) const
  {
    return typename Traits::VelocityRange(0.0);
  }

private:
  double viscosity_;
};

// kinetic energy 1/2 int |u|^2 of the discrete velocity
template<typename GV, typename VelocityFunction>
double kineticEnergy(const GV& gridView, const VelocityFunction& velocity)
{
  double energy = 0.0;
  for (const auto& element : elements(gridView))
  {
    const auto geometry = element.geometry();
    for (const auto& qp : Dune::QuadratureRules<double,GV::dimension>::rule(element.type(), 6))
    {
      typename VelocityFunction::Traits::RangeType u;
      velocity.evaluate(element, qp.position(), u);
      energy += 0.5 * u.two_norm2() * qp.weight() * geometry.integrationElement(qp.position());
    }
  }
  return energy;
}


int main(int argc, char** argv)
{
  try{
    // Maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // Create grid
    const int dim = 2;
    Dune::FieldVector<double,dim> lowerleft(0.0);
    Dune::FieldVector<double,dim> upperright(1.0);
    auto cells = Dune::filledArray<dim,unsigned int>(16);
    using Grid = Dune::YaspGrid<dim>;
    auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerleft, upperright, cells);
    using GridView = Grid::LeafGridView;
    GridView gridView = grid->leafGridView();
    using DomainField = GridView::Grid::ctype;
    using RangeType = double;
    const double eps = 1e-7; // accuracy for boundary detection

    // Taylor-Hood basis: velocity Q2, pressure Q1
    using namespace Dune::Functions::BasisFactory;
    using Dune::Indices::_0;
    using Dune::Indices::_1;
    auto basis = makeBasis(gridView,
                           composite(power<dim>(lagrange<2>(), flatInterleaved()),
                                     lagrange<1>(),
                                     flatLexicographic()));
    using Basis = decltype(basis);

    // Constraints: velocity Dirichlet on the whole boundary, no pressure boundary conditions
    auto bux = Dune::PDELab::makeBoundaryConditionFromCallable(gridView, [](const auto& x){ return true; });
    auto buy = Dune::PDELab::makeBoundaryConditionFromCallable(gridView, [](const auto& x){ return true; });
    auto bu = Dune::PDELab::CompositeConstraintsParameters<decltype(bux),decltype(buy)>(bux, buy);
    auto bp = Dune::PDELab::makeBoundaryConditionFromCallable(gridView, [](const auto& x){ return false; });
    auto bconstraints = Dune::PDELab::CompositeConstraintsParameters<decltype(bu),decltype(bp)>(bu, bp);
    using ConstraintsContainer = Dune::PDELab::ConstraintsContainer<Basis, RangeType>;
    ConstraintsContainer constraintsContainer;
    Dune::PDELab::constraints(bconstraints, basis, constraintsContainer);

    // Fix the first pressure DOF of the first element (the pressure is determined up to a constant)
    {
      auto localView = basis.localView();
      localView.bind(*gridView.template begin<0>());
      constraintsContainer[localView.index(localView.tree().child(_1).localIndex(0))];
    }

    // Print number of DOFs
    std::cout << "basis with " << Dune::PDELab::globalSize(basis) << " dofs generated  "<< std::endl;
    std::cout << "cc with " << constraintsContainer.size() << " dofs generated  "<< std::endl;

    // Local operator
    using Problem = DrivenCavityProblem<GridView>;
    Problem problem(1e-2);
    using LocalOperator = Dune::PDELab::TaylorHoodNavierStokes<Problem>;
    LocalOperator localOperator(problem);

    // Grid operator
    using MatrixBackend = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    MatrixBackend matrixBackend(75);
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

    // Solution vector with the lid velocity as Dirichlet values
    using CoefficientVector = Dune::PDELab::Backend::Vector<Basis, DomainField>;
    CoefficientVector coefficientVector(basis, 0.0);
    auto velocityLambda = [&](const auto& x){
      Dune::FieldVector<RangeType,dim> u(0.0);
      if (x[1] > 1.0-eps)
        u[0] = 1.0;
      return u;
    };
    auto velocityFunction = Dune::PDELab::makeGridFunctionFromCallable(gridView, velocityLambda);
    auto pressureFunction = Dune::PDELab::makeGridFunctionFromCallable(gridView, [](const auto& x){ return 0.0; });
    Dune::PDELab::CompositeGridFunction<decltype(velocityFunction),decltype(pressureFunction)>
      boundaryCondition(velocityFunction, pressureFunction);
    Dune::PDELab::interpolate(boundaryCondition, basis, coefficientVector);

    // Create Newton solver
    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_SuperLU;
    LinearSolver linearSolver(false);
    using Solver = Dune::PDELab::NewtonMethod<GridOperator, LinearSolver>;
    Solver solver(gridOperator, linearSolver);
    solver.setReduction(1e-10);
    solver.setVerbosityLevel(2);

    // Solve PDE
    solver.apply(coefficientVector);

    // Visualization
    using VTKWriter = Dune::SubsamplingVTKWriter<GridView>;
    Dune::RefinementIntervals subint(2);
    VTKWriter vtkwriter(gridView, subint);
    std::string vtkfile("testnavierstokes-functions");
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter, basis, coefficientVector,
                                         Dune::PDELab::vtk::defaultNameScheme());
    vtkwriter.write(vtkfile, Dune::VTK::ascii);

    // Kinetic energy of the solution
    auto velocitySubGFS = Dune::Functions::subspaceBasis(basis, _0);
    using VelocitySubGFS = decltype(velocitySubGFS);
    using VelocityDGF = Dune::PDELab::VectorDiscreteGridFunction<VelocitySubGFS, CoefficientVector>;
    VelocityDGF velocity(velocitySubGFS, coefficientVector);
    const double energy = kineticEnergy(gridView, velocity);
    std::cout << "kinetic energy: " << std::setprecision(12) << energy << std::endl;

    // Let the test fail if Newton did not converge or the energy differs from the reference
    const double referenceEnergy = 0.0319937020809;
    bool testfail(false);
    using std::abs;
    if (not solver.result().converged or abs(energy - referenceEnergy) > 1e-8)
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
