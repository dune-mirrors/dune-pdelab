// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/* Instantiates and applies all sequential ISTL solver backends of
 * backend/istl/seqistlsolverbackend.hh with a grid operator on a PDELab grid
 * function space and on a dune-functions basis.
 *
 * The problem is a Poisson problem with Dirichlet conditions (Q1). Every
 * backend solves it; backends that are expected to converge are checked
 * against a reference solution. Backends that are not meant for this kind of
 * system (Richardson without preconditioner, a few Jacobi loops, the explicit
 * diagonal solver) are only run.
 */

#include "config.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>

#include <dune/common/filledarray.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/istl/solvers.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dune/pdelab.hh>

class Checker
{
public:
  void check(bool condition, const std::string& what)
  {
    std::cout << (condition ? "ok     " : "FAILED ") << what << std::endl;
    failed_ = failed_ or not condition;
  }

  int result() const { return failed_ ? 1 : 0; }

private:
  bool failed_ = false;
};

// Finite element map parameter of the local operator, which only uses its type for the basis cache
template<typename Space>
struct FiniteElementMapType
{
  using GV = std::decay_t<decltype(std::declval<const Space&>().gridView())>;
  using type = Dune::PDELab::QkLocalFiniteElementMap<GV,double,double,1>;
};

template<typename Space>
  requires Dune::PDELab::Impl::DuneFunctionsBasis<Space>
struct FiniteElementMapType<Space>
{
  using type = Dune::PDELab::FiniteElementMapFor<typename Space::LocalView::Tree::FiniteElement>;
};

// Solve the linear problem of go with the backend ls, starting from x0
template<typename GO, typename LS, typename V>
V solveLinearProblem(const GO& go, LS& ls, const V& x0)
{
  V x(x0);
  Dune::PDELab::StationaryLinearProblemSolver<GO,LS,V> solver(go, ls, x, 1e-10, 1e-99, 0);
  solver.apply();
  return x;
}

// Same for matrix-free backends: solve J z = r(x0) with the jacobian applied on the fly
template<typename GO, typename LS, typename V>
V solveMatrixFree(const GO& go, LS& ls, const V& x0)
{
  using W = typename GO::Traits::Range;
  W r(go.testGridFunctionSpace(), 0.0);
  go.residual(x0, r);
  V z(go.trialGridFunctionSpace(), 0.0);
  ls.setLinearizationPoint(x0);
  ls.apply(z, r, 1e-10);
  V x(x0);
  x -= z;
  return x;
}

template<typename Space>
void testBackends(Checker& checker, const std::string& name, const Space& space)
{
  using namespace Dune::PDELab;
  using GV = std::decay_t<decltype(space.gridView())>;
  const auto& gv = space.gridView();

  using Problem = ConvectionDiffusionModelProblem<GV,double>;
  Problem problem;
  ConvectionDiffusionBoundaryConditionAdapter<Problem> bctype(gv, problem);
  using CC = Dune::PDELab::ConstraintsContainer<Space,double>;
  CC cc;
  constraints(bctype, space, cc);

  using FEM = typename FiniteElementMapType<Space>::type;
  using LOP = ConvectionDiffusionFEM<Problem,FEM>;
  LOP lop(problem);
  using MB = ISTL::BCRSMatrixBackend<>;
  using GO = GridOperator<Space,Space,LOP,MB,double,double,double,CC,CC>;
  GO go(space, cc, space, cc, lop, MB(9));

  using V = Backend::Vector<Space,double>;
  V x0(space, 0.0);
  ConvectionDiffusionDirichletExtensionAdapter<Problem> g(gv, problem);
  interpolate(g, space, x0);

  const int maxiter = 5000;
  const int verbose = 0;

  // reference solution
  ISTLBackend_SEQ_CG_ILU0 referenceSolver(maxiter, verbose);
  const V reference = solveLinearProblem(go, referenceSolver, x0);

  auto compare = [&](const std::string& backend, auto& ls, const V& x) {
    V difference(x);
    difference -= reference;
    const double error = difference.infinity_norm() / reference.infinity_norm();
    std::ostringstream message;
    message << name << " " << backend << ": converged " << ls.result().converged
            << ", iterations " << ls.result().iterations
            << ", relative error " << std::scientific << std::setprecision(2) << error;
    checker.check(ls.result().converged and error < 1e-6, message.str());
  };
  auto solveAndCompare = [&](const std::string& backend, auto&& ls) {
    compare(backend, ls, solveLinearProblem(go, ls, x0));
  };
  auto run = [&](const std::string& backend, auto&& ls) {
    solveLinearProblem(go, ls, x0);
    checker.check(true, name + " " + backend + ": instantiated and applied");
  };

  // Krylov solvers with point and incomplete factorization preconditioners
  solveAndCompare("SEQ_BCGS_Jac", ISTLBackend_SEQ_BCGS_Jac(maxiter, verbose));
  solveAndCompare("SEQ_BCGS_SSOR", ISTLBackend_SEQ_BCGS_SSOR(maxiter, verbose));
  solveAndCompare("SEQ_BCGS_ILU0", ISTLBackend_SEQ_BCGS_ILU0(maxiter, verbose));
  solveAndCompare("SEQ_BCGS_ILUn", ISTLBackend_SEQ_BCGS_ILUn(2, 1.0, maxiter, verbose));
  solveAndCompare("SEQ_CG_Jac", ISTLBackend_SEQ_CG_Jac(maxiter, verbose));
  solveAndCompare("SEQ_CG_SSOR", ISTLBackend_SEQ_CG_SSOR(maxiter, verbose));
  solveAndCompare("SEQ_CG_ILU0", ISTLBackend_SEQ_CG_ILU0(maxiter, verbose));
  solveAndCompare("SEQ_CG_ILUn", ISTLBackend_SEQ_CG_ILUn(2, 1.0, maxiter, verbose));
  solveAndCompare("SEQ_MINRES_SSOR", ISTLBackend_SEQ_MINRES_SSOR(maxiter, verbose));
  solveAndCompare("SEQ_GMRES_ILU0", ISTLBackend_SEQ_GMRES_ILU0(200, maxiter, verbose));
  solveAndCompare("SEQ_BCGS_Richardson", ISTLBackend_SEQ_BCGS_Richardson(maxiter, verbose));
  solveAndCompare("SEQ_ILU0<CGSolver>", ISTLBackend_SEQ_ILU0<Dune::CGSolver>(maxiter, verbose));
  solveAndCompare("SEQ_ILUn<BiCGSTABSolver>", ISTLBackend_SEQ_ILUn<Dune::BiCGSTABSolver>(1, 1.0, maxiter, verbose));

  // algebraic multigrid (uses the grid dimension of the space)
  solveAndCompare("SEQ_CG_AMG_SSOR", ISTLBackend_SEQ_CG_AMG_SSOR<GO>(maxiter, verbose));
  solveAndCompare("SEQ_BCGS_AMG_SSOR", ISTLBackend_SEQ_BCGS_AMG_SSOR<GO>(maxiter, verbose));
  solveAndCompare("SEQ_BCGS_AMG_SOR", ISTLBackend_SEQ_BCGS_AMG_SOR<GO>(maxiter, verbose));
  solveAndCompare("SEQ_LS_AMG_SSOR", ISTLBackend_SEQ_LS_AMG_SSOR<GO>(maxiter, verbose));
  solveAndCompare("SEQ_LS_AMG_SOR", ISTLBackend_SEQ_LS_AMG_SOR<GO>(maxiter, verbose));

  // direct solvers
#if HAVE_SUPERLU
  solveAndCompare("SEQ_SuperLU", ISTLBackend_SEQ_SuperLU(verbose));
#endif
#if HAVE_SUITESPARSE_UMFPACK
  solveAndCompare("SEQ_UMFPack", ISTLBackend_SEQ_UMFPack(verbose));
#endif

  // matrix-free: the jacobian is applied on the fly by the grid operator
  {
    ISTLBackend_SEQ_MatrixFree_BCGS_Richardson<GO> ls(go, maxiter, verbose);
    compare("SEQ_MatrixFree_BCGS_Richardson", ls, solveMatrixFree(go, ls, x0));
  }
  {
    ISTLBackend_SEQ_MatrixFree_Richardson<GO, Dune::CGSolver> ls(go, maxiter, verbose);
    compare("SEQ_MatrixFree_Richardson<CGSolver>", ls, solveMatrixFree(go, ls, x0));
  }

  // not meant to solve this system to the reference accuracy: only run them
  run("SEQ_LOOP_Jac", ISTLBackend_SEQ_LOOP_Jac(10, verbose));
  run("SEQ_Richardson<LoopSolver>", ISTLBackend_SEQ_Richardson<Dune::LoopSolver>(10, verbose));
  run("SEQ_ExplicitDiagonal", ISTLBackend_SEQ_ExplicitDiagonal());
}

int main(int argc, char** argv)
{
  try {
    Dune::MPIHelper::instance(argc, argv);

    const int dim = 2;
    using Grid = Dune::YaspGrid<dim>;
    Grid grid(Dune::FieldVector<double,dim>(1.0), Dune::filledArray<dim,int>(16));
    auto gv = grid.leafGridView();
    using GV = decltype(gv);

    Checker checker;

    using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV,double,double,1>;
    FEM fem(gv);
    using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,
                                                Dune::PDELab::ISTL::VectorBackend<>>;
    GFS gfs(gv, fem);
    testBackends(checker, "GridFunctionSpace", gfs);

    auto basis = Dune::Functions::BasisFactory::makeBasis(gv, Dune::Functions::BasisFactory::lagrange<1>());
    testBackends(checker, "dune-functions basis", basis);

    return checker.result();
  }
  catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
}
