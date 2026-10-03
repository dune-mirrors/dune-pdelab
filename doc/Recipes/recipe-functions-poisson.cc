#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <iostream>

#include <dune/common/filledarray.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/pdelab.hh>

/**
 * \page recipe-functions-poisson Solving the Poisson equation on a dune-functions basis
 *
 * This recipe solves \f$-\Delta u = f\f$ with Dirichlet conditions on a basis from
 * <a href="https://dune-project.org/modules/dune-functions/">dune-functions</a>.
 * The local operator and the solvers are the same as for a grid function space.
 * See \ref pdelab_functions_bases for an overview.
 * The problem is chosen such that the exact solution \f$u = x^2 + y^2\f$ is
 * in the discrete space, so we can check the error at the end.
 *
 * The basis replaces the finite element map and the grid function space:
 * \snippet recipe-functions-poisson.cc Basis
 *
 * Local operators that take a finite element map type get it from the finite
 * element of the basis:
 * \snippet recipe-functions-poisson.cc Local operator
 *
 * The constraints container and the constraints assembly accept a basis.
 * \snippet recipe-functions-poisson.cc Constraints
 *
 * The basis is used as trial and test space of the grid operator.
 * The matrix backend needs the maximal number of entries per row.
 * \snippet recipe-functions-poisson.cc Grid operator
 *
 * The coefficient vector is created from the basis. Interpolating the Dirichlet
 * extension sets the values at constrained DOFs and serves as initial guess.
 * \snippet recipe-functions-poisson.cc Vector
 *
 * Solving works as for a grid function space.
 * \snippet recipe-functions-poisson.cc Solve
 *
 * The point of using a basis is that the result can be used with other
 * dune-functions code. Here, we wrap the coefficients in a grid function
 * and compute the \f$L^2\f$ error with a local function.
 * `Backend::native(x)` gives the underlying ISTL vector.
 * \snippet recipe-functions-poisson.cc Interoperability
 *
 * For the output, `addSolutionToVTKWriter` works for bases.
 * Alternatively, the grid function from above can be passed to the VTK writer
 * with `vtkwriter.addVertexData(u, Dune::VTK::FieldInfo("u", Dune::VTK::FieldInfo::Type::scalar, 1))`.
 * \snippet recipe-functions-poisson.cc VTK
 *
 * Full example code: @ref recipe-functions-poisson.cc
 * \example recipe-functions-poisson.cc
 * See explanation at @ref recipe-functions-poisson
 */

// Poisson problem with exact solution u = x^2 + y^2, i.e. -Laplace(u) = -4
template<typename GV, typename RF>
class PolynomialProblem
  : public Dune::PDELab::ConvectionDiffusionModelProblem<GV,RF>
{
  using Base = Dune::PDELab::ConvectionDiffusionModelProblem<GV,RF>;
public:
  using Traits = typename Base::Traits;

  static RF exact(const typename Traits::DomainType& x)
  {
    return x.two_norm2();
  }

  RF f(const typename Traits::ElementType&, const typename Traits::DomainType&) const
  {
    return -4.0;
  }

  RF g(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return exact(e.geometry().global(x));
  }
};

int main(int argc, char** argv)
{
  try {
    Dune::MPIHelper::instance(argc, argv);

    constexpr int dim = 2;
    using Grid = Dune::YaspGrid<dim>;
    Grid grid({1.0, 1.0}, Dune::filledArray<dim,int>(8));
    using GV = Grid::LeafGridView;
    GV gv = grid.leafGridView();
    using Problem = PolynomialProblem<GV,double>;
    Problem problem;

    // [Basis]
    auto basis = Dune::Functions::BasisFactory::makeBasis(
      gv, Dune::Functions::BasisFactory::lagrange<2>());
    using Basis = decltype(basis);
    //! [Basis]

    // [Local operator]
    using FEM = Dune::PDELab::FiniteElementMapFor<
      typename Basis::LocalView::Tree::FiniteElement>;
    using LOP = Dune::PDELab::ConvectionDiffusionFEM<Problem,FEM>;
    LOP lop(problem);
    //! [Local operator]

    // [Constraints]
    using CC = Dune::PDELab::ConstraintsContainer<Basis,double>;
    CC cc;
    cc.clear();
    Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> bctype(gv,problem);
    Dune::PDELab::constraints(bctype,basis,cc);
    //! [Constraints]

    // [Grid operator]
    using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    MBE mbe(4 * Dune::PDELab::maxLocalSize(basis));
    using GO = Dune::PDELab::GridOperator<Basis,Basis,LOP,MBE,
                                          double,double,double,CC,CC>;
    GO go(basis,cc,basis,cc,lop,mbe);
    //! [Grid operator]

    // [Vector]
    using V = Dune::PDELab::Backend::Vector<Basis,double>;
    V x(basis,0.0);
    Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> g(gv,problem);
    Dune::PDELab::interpolate(g,basis,x);
    //! [Vector]

    // [Solve]
    // The AMG backends are templated on the grid operator.
    using LS = Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO>;
    LS ls(5000,1);
    using SLP = Dune::PDELab::StationaryLinearProblemSolver<GO,LS,V>;
    SLP slp(go,ls,x,1e-10);
    slp.apply();
    //! [Solve]

    // [Interoperability]
    auto u = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(
      basis, Dune::PDELab::Backend::native(x));
    auto localU = localFunction(u);
    double error2 = 0.0;
    for (const auto& element : elements(gv))
    {
      localU.bind(element);
      const auto geometry = element.geometry();
      // order 4 is enough: the integrand is a polynomial of degree 4
      for (const auto& qp : Dune::QuadratureRules<double,dim>::rule(element.type(), 4))
      {
        const double diff = localU(qp.position())
                            - Problem::exact(geometry.global(qp.position()));
        error2 += diff * diff * qp.weight() * geometry.integrationElement(qp.position());
      }
    }
    const double error = std::sqrt(error2);
    std::cout << "L2 error: " << error << std::endl;
    //! [Interoperability]

    // [VTK]
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter,basis,x);
    vtkwriter.write("recipe-functions-poisson",Dune::VTK::ascii);
    //! [VTK]

    // P2 reproduces the exact solution up to the solver tolerance
    return error < 1e-8 ? 0 : 1;
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
