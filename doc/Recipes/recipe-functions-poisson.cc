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
 * \page recipe-functions-poisson Solving the Poisson equation with PDELab
 *
 * This recipe solves the Poisson equation with PDELab. The finite element space
 * is described by a basis from
 * <a href="https://www.dune-project.org/modules/dune-functions/">dune-functions</a>.
 *
 * \section recipe-functions-poisson-problem Problem and weak form
 *
 * We solve
 * \f[
 *   -\Delta u = f \;\mbox{ in } \Omega, \qquad u = g \;\mbox{ on } \partial\Omega .
 * \f]
 * The weak form is: find \f$u \in V\f$ with \f$u = g\f$ on \f$\partial\Omega\f$ such that
 * \f[
 *   \int_\Omega \nabla u \cdot \nabla v \,dx = \int_\Omega f v \,dx
 * \f]
 * for all test functions \f$v\f$ that vanish on \f$\partial\Omega\f$.
 * PDELab works with the residual form
 * \f[
 *   r(u,v) = \int_\Omega \nabla u \cdot \nabla v - f v \,dx .
 * \f]
 * A local operator evaluates \f$r\f$ (and its Jacobian) element by element.
 *
 * The local operator `ConvectionDiffusionFEM` reads its data from a parameter class.
 * It is a general convection-diffusion-reaction operator, which discretizes
 * \f$\nabla\cdot(-A \nabla u + b u) + c u = f\f$ with Dirichlet, Neumann and outflow
 * conditions. Poisson is the special case \f$A = I\f$, \f$b = 0\f$, \f$c = 0\f$.
 * These are the defaults of `ConvectionDiffusionModelProblem`, so we only override what differs:
 * - `A(e,x)`: diffusion tensor (default: identity),
 * - `f(e,x)`: source term,
 * - `g(e,x)`: Dirichlet values,
 * - `bctype(is,x)`: boundary condition type per intersection (default: Dirichlet everywhere).
 *
 * The problem is chosen such that the exact solution \f$u = x^2 + y^2\f$ is in the
 * discrete space. Then \f$f = -4\f$ and we can check the error at the end.
 * \snippet recipe-functions-poisson.cc Problem
 *
 * \section recipe-functions-poisson-basis Discrete space: the basis
 *
 * The finite element space is a dune-functions basis. Here it is the Lagrange basis of
 * order two on the leaf grid view.
 * The basis provides the global numbering of the degrees of freedom (DOFs)
 * and a local finite element for each element.
 * \snippet recipe-functions-poisson.cc Basis
 *
 * More complex spaces (vector-valued spaces, power and composite bases such as Taylor-Hood,
 * index merging strategies and blocking of the coefficient vector) are described in
 * - the dune-functions documentation
 *   (<a href="https://www.dune-project.org/modules/dune-functions/">module page</a>,
 *   Doxygen "Function space bases"),
 * - C. Engwer, C. Gräser, S. Müthing, O. Sander:
 *   <i>Function space bases in the dune-functions module</i>,
 *   <a href="https://arxiv.org/abs/1806.09545">arXiv:1806.09545</a>,
 * - \ref recipe-functions-navier-stokes, which uses a composite basis in PDELab.
 *
 * Dirichlet DOFs are constrained. `constraints` marks all DOFs on boundary
 * intersections where `bctype` says Dirichlet in the `ConstraintsContainer`.
 * \snippet recipe-functions-poisson.cc Constraints
 *
 * \section recipe-functions-poisson-assembly Assembling the local operator
 *
 * The local operator is created from the problem.
 * \remark `FiniteElementMapFor<...>` is only needed because some local operators take a finite
 * element map type as template parameter for historical reasons.
 * It extracts that type from the finite element of the basis.
 * \snippet recipe-functions-poisson.cc Local operator
 *
 * The grid operator loops over the grid and calls the local operator on each element.
 * The assembly itself is done by dune-assembler. The grid operator assembles the global
 * residual vector and the Jacobian matrix.
 * Its template arguments are: trial basis, test basis, local operator, matrix backend,
 * field types of the coefficients, the Jacobian and the residual,
 * and the constraints containers for trial and test space.
 * The matrix backend is told the maximal number of entries per row (an estimate used to
 * preallocate the sparse matrix).
 * \snippet recipe-functions-poisson.cc Grid operator
 *
 * The coefficient vector is created from the basis.
 * Interpolating the Dirichlet values sets the values at the constrained DOFs
 * and serves as initial guess.
 * \snippet recipe-functions-poisson.cc Vector
 *
 * \section recipe-functions-poisson-solve Solving
 *
 * The problem is linear. `StationaryLinearProblemSolver` assembles the residual and the Jacobian
 * and solves the correction equation with a linear solver backend.
 * Here this is CG with an AMG preconditioner (`ISTLBackend_SEQ_CG_AMG_SSOR`).
 * The last argument, `1e-10`, is the required reduction of the linear solver defect.
 * Nonlinear problems use `NewtonMethod` with the same grid operator.
 * \snippet recipe-functions-poisson.cc Solve
 *
 * \section recipe-functions-poisson-result Using the solution
 *
 * `makeDiscreteGlobalBasisFunction` turns the basis and the coefficients into a
 * dune-functions grid function. We evaluate it with a local function and a quadrature rule
 * to compute the \f$L^2\f$ error.
 * `Backend::native(x)` gives the underlying ISTL vector.
 * \snippet recipe-functions-poisson.cc Error
 *
 * `addSolutionToVTKWriter` writes the solution for visualization.
 * \snippet recipe-functions-poisson.cc VTK
 *
 * \note If you know PDELab from grid function spaces: the basis replaces the finite element map
 * and the grid function space. See \ref pdelab_functions_bases for what changes.
 *
 * Full example code: @ref recipe-functions-poisson.cc
 * \example recipe-functions-poisson.cc
 * See explanation at @ref recipe-functions-poisson
 */

// [Problem]
// Poisson problem with exact solution u = x^2 + y^2, i.e. -Laplace(u) = -4.
// A = I, b = 0, c = 0 and Dirichlet conditions everywhere are the defaults.
template<typename GV, typename RF>
class PoissonProblem
  : public Dune::PDELab::ConvectionDiffusionModelProblem<GV,RF>
{
  using Base = Dune::PDELab::ConvectionDiffusionModelProblem<GV,RF>;
public:
  using Traits = typename Base::Traits;

  static RF exact(const typename Traits::DomainType& x)
  {
    return x.two_norm2();
  }

  // source term f
  RF f(const typename Traits::ElementType&, const typename Traits::DomainType&) const
  {
    return -4.0;
  }

  // Dirichlet values g
  RF g(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return exact(e.geometry().global(x));
  }
};
//! [Problem]

int main(int argc, char** argv)
{
  try {
    Dune::MPIHelper::instance(argc, argv);

    constexpr int dim = 2;
    using Grid = Dune::YaspGrid<dim>;
    Grid grid({1.0, 1.0}, Dune::filledArray<dim,int>(8));
    using GV = Grid::LeafGridView;
    GV gv = grid.leafGridView();
    using Problem = PoissonProblem<GV,double>;
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

    // [Error]
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
    //! [Error]

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
