#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <iomanip>
#include <iostream>

#include <dune/common/filledarray.hh>
#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/pdelab.hh>

/**
 * \page recipe-functions-navier-stokes Stationary Navier-Stokes with a Taylor-Hood basis
 *
 * This recipe solves the stationary lid-driven cavity problem with Q2/Q1
 * Taylor-Hood elements on a basis from dune-functions.
 * It is a trimmed version of `testnavierstokes-functions.cc`.
 * See \ref pdelab_functions_bases for an overview.
 *
 * The basis is a composite of a power basis for the velocity and a Lagrange
 * basis for the pressure. The velocity components are interleaved
 * (`flatInterleaved`), the two blocks are lexicographic (`flatLexicographic`).
 * This gives one flat `BlockVector<double>`. One level of blocking would work
 * as well, but nested blocking (e.g. `blockedLexicographic` at the top) is not supported.
 * \snippet recipe-functions-navier-stokes.cc Basis
 *
 * Boundary conditions are given per component, with the same parameter classes
 * as for a grid function space. Velocity is Dirichlet everywhere,
 * the pressure is free.
 * \snippet recipe-functions-navier-stokes.cc Constraints
 *
 * Without pressure boundary conditions the pressure is determined up to a constant.
 * We fix one pressure DOF, so a direct solver can be used.
 * `LocalFunctionSpace` and `LFSIndexCache` accept a basis. The code is the same as for a
 * grid function space. Accessing an entry of the container creates it, which
 * marks the DOF as constrained with value zero.
 * \snippet recipe-functions-navier-stokes.cc Pin pressure
 *
 * The local operator and the grid operator are the same as for a grid function space.
 * The basis is used as trial and test space.
 * \snippet recipe-functions-navier-stokes.cc Grid operator
 *
 * The Dirichlet values are interpolated from a composite grid function.
 * The lid moves with velocity 1, everything else is zero.
 * \snippet recipe-functions-navier-stokes.cc Interpolation
 *
 * Newton's method with a direct solver. The result tells us if it converged.
 * \snippet recipe-functions-navier-stokes.cc Newton
 *
 * For the post-processing we use dune-functions. `subspaceBasis` selects the
 * velocity, `makeDiscreteGlobalBasisFunction` wraps the coefficients.
 * The kinetic energy is computed with a local function.
 * \snippet recipe-functions-navier-stokes.cc Post-processing
 *
 * The solution is written with `addSolutionToVTKWriter`. The default name scheme
 * names the components `solution_0` (velocity) and `solution_1` (pressure).
 * \snippet recipe-functions-navier-stokes.cc VTK
 *
 * Full example code: @ref recipe-functions-navier-stokes.cc
 * \example recipe-functions-navier-stokes.cc
 * See explanation at @ref recipe-functions-navier-stokes
 */

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

int main(int argc, char** argv)
{
  try {
    Dune::MPIHelper::instance(argc, argv);

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

    // [Basis]
    using namespace Dune::Functions::BasisFactory;
    using Dune::Indices::_0;
    using Dune::Indices::_1;
    auto basis = makeBasis(gridView,
                           composite(power<dim>(lagrange<2>(), flatInterleaved()),
                                     lagrange<1>(),
                                     flatLexicographic()));
    using Basis = decltype(basis);
    //! [Basis]

    // [Constraints]
    auto bux = Dune::PDELab::makeBoundaryConditionFromCallable(gridView, [](const auto& x){ return true; });
    auto buy = Dune::PDELab::makeBoundaryConditionFromCallable(gridView, [](const auto& x){ return true; });
    auto bu = Dune::PDELab::CompositeConstraintsParameters<decltype(bux),decltype(buy)>(bux, buy);
    auto bp = Dune::PDELab::makeBoundaryConditionFromCallable(gridView, [](const auto& x){ return false; });
    auto bconstraints = Dune::PDELab::CompositeConstraintsParameters<decltype(bu),decltype(bp)>(bu, bp);
    using ConstraintsContainer = Dune::PDELab::ConstraintsContainer<Basis, RangeType>;
    ConstraintsContainer constraintsContainer;
    Dune::PDELab::constraints(bconstraints, basis, constraintsContainer);
    //! [Constraints]

    // [Pin pressure]
    {
      using LFS = Dune::PDELab::LocalFunctionSpace<Basis>;
      LFS lfs(basis);
      Dune::PDELab::LFSIndexCache<LFS> lfsCache(lfs);
      lfs.bind(*gridView.template begin<0>());
      lfsCache.update();
      // first pressure DOF of the first element
      constraintsContainer[lfsCache.containerIndex(lfs.child(_1).localIndex(0))];
    }
    //! [Pin pressure]

    // [Grid operator]
    using Problem = DrivenCavityProblem<GridView>;
    Problem problem(1e-2);
    using LocalOperator = Dune::PDELab::TaylorHoodNavierStokes<Problem>;
    LocalOperator localOperator(problem);

    using MatrixBackend = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    MatrixBackend matrixBackend(75);
    using GridOperator = Dune::PDELab::GridOperator<Basis, Basis,
                                                    LocalOperator,
                                                    MatrixBackend,
                                                    DomainField, RangeType, RangeType,
                                                    ConstraintsContainer,
                                                    ConstraintsContainer>;
    GridOperator gridOperator(basis, constraintsContainer,
                              basis, constraintsContainer,
                              localOperator, matrixBackend);
    //! [Grid operator]

    // [Interpolation]
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
    //! [Interpolation]

    // [Newton]
    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_SuperLU;
    LinearSolver linearSolver(false);
    using Solver = Dune::PDELab::NewtonMethod<GridOperator, LinearSolver>;
    Solver solver(gridOperator, linearSolver);
    solver.setReduction(1e-10);
    solver.setVerbosityLevel(2);
    solver.apply(coefficientVector);
    const bool converged = solver.result().converged;
    //! [Newton]

    // [Post-processing]
    auto velocityBasis = Dune::Functions::subspaceBasis(basis, _0);
    auto velocity = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,dim>>(
      velocityBasis, Dune::PDELab::Backend::native(coefficientVector));
    auto localVelocity = localFunction(velocity);
    double energy = 0.0;
    for (const auto& element : elements(gridView))
    {
      localVelocity.bind(element);
      const auto geometry = element.geometry();
      for (const auto& qp : Dune::QuadratureRules<double,dim>::rule(element.type(), 6))
        energy += 0.5 * localVelocity(qp.position()).two_norm2()
                  * qp.weight() * geometry.integrationElement(qp.position());
    }
    std::cout << "kinetic energy: " << std::setprecision(12) << energy << std::endl;
    //! [Post-processing]

    // [VTK]
    using VTKWriter = Dune::SubsamplingVTKWriter<GridView>;
    VTKWriter vtkwriter(gridView, Dune::RefinementIntervals(2));
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter, basis, coefficientVector,
                                         Dune::PDELab::vtk::defaultNameScheme());
    vtkwriter.write("recipe-functions-navier-stokes", Dune::VTK::ascii);
    //! [VTK]

    // Reference value from testnavierstokes-functions.cc
    const double referenceEnergy = 0.0319937020809;
    using std::abs;
    return (converged and abs(energy - referenceEnergy) < 1e-8) ? 0 : 1;
  }
  catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
