// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/* Compares the generic one step grid operator (Impl::GenericOneStepGridOperator,
 * which delegates to a grid operator with a weighted sum of the two local
 * operators) with the engine based OneStepGridOperator on grid function spaces:
 * residual, jacobian and jacobian_apply of every stage of implicit methods,
 * explicit_jacobian_residual of explicit methods, and complete time steps.
 * Continuous Galerkin with GridOperator, DG with GridOperator and with
 * FastDGGridOperator.
 */

#include "config.h"

#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/filledarray.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab.hh>

#include "convectiondiffusionfastdg.hh"
#include "l2fastdg.hh"

// Time dependent convection-diffusion problem with Dirichlet boundary conditions
template<typename GV, typename RF>
class Problem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  static constexpr bool permeabilityIsConstantPerCell()
  {
    return true;
  }

  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType I(0.0);
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      I[i][i] = 1.0;
    return I;
  }

  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType v(0.0);
    v[0] = 1.0 + time_;
    v[1] = 0.5;
    return v;
  }

  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.5;
  }

  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    auto xglobal = e.geometry().global(x);
    return std::sin(1.0 + 3.0*time_) * (1.0 + xglobal[0]*xglobal[1]);
  }

  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    auto xglobal = e.geometry().global(x);
    return std::exp(-time_) * (xglobal[0] + 2.0*xglobal[1]);
  }

  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  void setTime (RF t)
  {
    time_ = t;
  }

private:
  RF time_ = 0.0;
};

template<typename V>
void randomize(V& v, unsigned int seed)
{
  std::mt19937 generator(seed);
  std::uniform_real_distribution<double> distribution(-1.0, 1.0);
  for (auto& block : Dune::PDELab::Backend::native(v))
    for (auto& entry : block)
      entry = distribution(generator);
}

class Checker
{
public:
  explicit Checker(double tolerance) : tolerance_(tolerance) {}

  template<typename V>
  void vectors(const std::string& what, const V& a, const V& b)
  {
    V difference(a);
    difference -= b;
    const double scale = std::max(1.0, a.infinity_norm());
    report(what, difference.infinity_norm() / scale);
  }

  // compare by the number of nonzeros and the application to random vectors
  template<typename M, typename V>
  void matrices(const std::string& what, const M& a, const M& b, const V& prototype)
  {
    using Dune::PDELab::Backend::native;
    if (native(a).nonzeroes() != native(b).nonzeroes())
    {
      std::cout << "FAILED " << what << ": nonzeros " << native(a).nonzeroes()
                << " != " << native(b).nonzeroes() << std::endl;
      failed_ = true;
      return;
    }
    double error = 0.0;
    for (unsigned int seed = 100; seed < 103; ++seed)
    {
      V x(prototype), ya(prototype), yb(prototype);
      randomize(x, seed);
      native(a).mv(native(x), native(ya));
      native(b).mv(native(x), native(yb));
      ya -= yb;
      error = std::max(error, ya.infinity_norm() / std::max(1.0, yb.infinity_norm()));
    }
    report(what, error);
  }

  bool failed() const { return failed_; }

private:
  void report(const std::string& what, double error)
  {
    const bool ok = error <= tolerance_;
    std::cout << (ok ? "ok     " : "FAILED ") << what << ": " << error << std::endl;
    failed_ = failed_ or not ok;
  }

  double tolerance_;
  bool failed_ = false;
};

// residual, jacobian and jacobian_apply of every stage of an implicit method
template<bool withJacobianApply = true, typename IGOOld, typename IGONew, typename GFS>
void compareImplicit(Checker& check, const std::string& name, IGOOld& oldOperator, IGONew& newOperator,
                     const Dune::PDELab::TimeSteppingParameterInterface<double>& method, const GFS& gfs)
{
  using V = typename IGOOld::Traits::Domain;
  using J = typename IGOOld::Traits::Jacobian;

  oldOperator.preStep(method, 0.3, 0.1);
  newOperator.preStep(method, 0.3, 0.1);

  std::vector<V> solutions(method.s(), V(gfs, 0.0));
  std::vector<V*> x;
  for (unsigned int s = 0; s < method.s(); ++s)
  {
    randomize(solutions[s], s + 1);
    x.push_back(&solutions[s]);
  }

  for (unsigned int r = 1; r <= method.s(); ++r)
  {
    const std::string stage = name + " " + method.name() + " stage " + std::to_string(r) + " ";
    std::vector<V*> previous(x.begin(), x.begin() + r);
    oldOperator.preStage(r, previous);
    newOperator.preStage(r, previous);

    V z(gfs);
    randomize(z, 10 + r);

    V rOld(gfs, 0.0), rNew(gfs, 0.0);
    oldOperator.residual(z, rOld);
    newOperator.residual(z, rNew);
    check.vectors(stage + "residual", rOld, rNew);

    J jOld(oldOperator), jNew(newOperator);
    oldOperator.jacobian(z, jOld);
    newOperator.jacobian(z, jNew);
    check.matrices(stage + "jacobian", jOld, jNew, z);

    if constexpr (withJacobianApply)
    {
      V yOld(gfs, 0.0), yNew(gfs, 0.0);
      oldOperator.jacobian_apply(z, yOld);
      newOperator.jacobian_apply(z, yNew);
      check.vectors(stage + "jacobian_apply", yOld, yNew);
    }

    oldOperator.postStage();
    newOperator.postStage();
  }
  oldOperator.postStep();
  newOperator.postStep();
}

// explicit_jacobian_residual of every stage of an explicit method
template<typename IGOOld, typename IGONew, typename GFS, typename CC>
void compareExplicit(Checker& check, const std::string& name, IGOOld& oldOperator, IGONew& newOperator,
                     const Dune::PDELab::TimeSteppingParameterInterface<double>& method, const GFS& gfs,
                     const CC& cc)
{
  using V = typename IGOOld::Traits::Domain;
  using J = typename IGOOld::Traits::Jacobian;

  oldOperator.preStep(method, 0.3, 0.1);
  newOperator.preStep(method, 0.3, 0.1);

  std::vector<V> solutions(method.s() + 1, V(gfs, 0.0));
  std::vector<V*> x;
  for (unsigned int s = 0; s <= method.s(); ++s)
  {
    randomize(solutions[s], 20 + s);
    x.push_back(&solutions[s]);
  }

  for (unsigned int r = 1; r <= method.s(); ++r)
  {
    const std::string stage = name + " " + method.name() + " stage " + std::to_string(r) + " ";
    std::vector<V*> upToStage(x.begin(), x.begin() + r + 1);
    J dOld(oldOperator), dNew(newOperator);
    V alphaOld(gfs, 0.0), alphaNew(gfs, 0.0), betaOld(gfs, 0.0), betaNew(gfs, 0.0);
    oldOperator.explicit_jacobian_residual(r, upToStage, dOld, alphaOld, betaOld);
    newOperator.explicit_jacobian_residual(r, upToStage, dNew, alphaNew, betaNew);
    check.matrices(stage + "D", dOld, dNew, alphaOld);
    // The engine based operator does not constrain the spatial part beta
    // (only the last of the combined grid operators applies the constraints).
    // ExplicitOneStepMethod sets the constrained entries of the combined
    // residual to zero anyway, so compare the unconstrained entries.
    Dune::PDELab::set_constrained_dofs(cc, 0.0, betaOld);
    Dune::PDELab::set_constrained_dofs(cc, 0.0, betaNew);
    check.vectors(stage + "alpha", alphaOld, alphaNew);
    check.vectors(stage + "beta", betaOld, betaNew);
    oldOperator.postStage();
    newOperator.postStage();
  }
  oldOperator.postStep();
  newOperator.postStep();
}

// a few complete time steps of an implicit method
template<typename IGO, typename V, typename G>
V implicitSteps(IGO& igo, const Dune::PDELab::TimeSteppingParameterInterface<double>& method,
                const V& initial, G& g)
{
  using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_SuperLU;
  LinearSolver linearSolver(false);
  using Solver = Dune::PDELab::StationaryLinearProblemSolver<IGO, LinearSolver, V>;
  Solver solver(igo, linearSolver, 1e-12);
  Dune::PDELab::OneStepMethod<double, IGO, Solver, V, V> osm(method, igo, solver);
  osm.setVerbosityLevel(0);
  V x(initial);
  double time = 0.0;
  const double dt = 0.05;
  for (int step = 0; step < 3; ++step)
  {
    V xnew(x);
    osm.apply(time, dt, x, g, xnew);
    x = xnew;
    time += dt;
  }
  return x;
}

// a few complete time steps of an explicit method
template<typename IGO, typename V, typename G>
V explicitSteps(IGO& igo, const Dune::PDELab::TimeSteppingParameterInterface<double>& method,
                const V& initial, G& g)
{
  using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_SuperLU;
  LinearSolver linearSolver(false);
  Dune::PDELab::ExplicitOneStepMethod<double, IGO, LinearSolver, V, V> osm(method, igo, linearSolver);
  osm.setVerbosityLevel(0);
  V x(initial);
  double time = 0.0;
  const double dt = 1e-4;
  for (int step = 0; step < 3; ++step)
  {
    V xnew(x);
    osm.apply(time, dt, x, g, xnew);
    x = xnew;
    time += dt;
  }
  return x;
}

template<typename GV>
void testContinuousGalerkin(Checker& check, const GV& gv)
{
  using RF = double;
  using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV, typename GV::ctype, RF, 2>;
  FEM fem(gv);
  using Constraints = Dune::PDELab::ConformingDirichletConstraints;
  using GFS = Dune::PDELab::GridFunctionSpace<GV, FEM, Constraints, Dune::PDELab::ISTL::VectorBackend<>>;
  GFS gfs(gv, fem);

  using P = Problem<GV,RF>;
  P problem;
  using CC = typename GFS::template ConstraintsContainer<RF>::Type;
  CC cc;
  Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<P> bctype(gv, problem);
  Dune::PDELab::constraints(bctype, gfs, cc);

  using LOP = Dune::PDELab::ConvectionDiffusionFEM<P, FEM>;
  LOP lop(problem);
  using MLOP = Dune::PDELab::L2;
  MLOP mlop;

  using MB = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  MB mb(25);
  using GO0 = Dune::PDELab::GridOperator<GFS,GFS,LOP,MB,RF,RF,RF,CC,CC>;
  using GO1 = Dune::PDELab::GridOperator<GFS,GFS,MLOP,MB,RF,RF,RF,CC,CC>;
  // separate grid operators, since the engine based operator changes the
  // pre- and postprocessing flags of the grid operators it is given
  GO0 go0Old(gfs, cc, gfs, cc, lop, mb), go0New(gfs, cc, gfs, cc, lop, mb);
  GO1 go1Old(gfs, cc, gfs, cc, mlop, mb), go1New(gfs, cc, gfs, cc, mlop, mb);

  using G = Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<P>;
  G g(gv, problem);
  using V = Dune::PDELab::Backend::Vector<GFS, RF>;
  V initial(gfs, 0.0);
  problem.setTime(0.0);
  Dune::PDELab::interpolate(g, gfs, initial);

  {
    using IGOOld = Dune::PDELab::OneStepGridOperator<GO0,GO1>;
    using IGONew = Dune::PDELab::Impl::GenericOneStepGridOperator<GO0,GO1>;
    IGOOld oldOperator(go0Old, go1Old);
    IGONew newOperator(go0New, go1New);

    Dune::PDELab::OneStepThetaParameter<RF> theta(0.5);
    Dune::PDELab::Alexander2Parameter<RF> alexander2;
    Dune::PDELab::FractionalStepParameter<RF> fractionalStep;
    compareImplicit(check, "CG", oldOperator, newOperator, theta, gfs);
    compareImplicit(check, "CG", oldOperator, newOperator, alexander2, gfs);
    compareImplicit(check, "CG", oldOperator, newOperator, fractionalStep, gfs);

    auto xOld = implicitSteps(oldOperator, alexander2, initial, g);
    auto xNew = implicitSteps(newOperator, alexander2, initial, g);
    check.vectors("CG Alexander2 three time steps", xOld, xNew);

    oldOperator.divideMassTermByDeltaT();
    newOperator.divideMassTermByDeltaT();
    compareImplicit(check, "CG divideMassTermByDeltaT", oldOperator, newOperator, alexander2, gfs);
  }

  {
    GO0 go0OldE(gfs, cc, gfs, cc, lop, mb), go0NewE(gfs, cc, gfs, cc, lop, mb);
    GO1 go1OldE(gfs, cc, gfs, cc, mlop, mb), go1NewE(gfs, cc, gfs, cc, mlop, mb);
    using IGOOld = Dune::PDELab::OneStepGridOperator<GO0,GO1,false>;
    using IGONew = Dune::PDELab::Impl::GenericOneStepGridOperator<GO0,GO1,false>;
    IGOOld oldOperator(go0OldE, go1OldE);
    IGONew newOperator(go0NewE, go1NewE);

    Dune::PDELab::ExplicitEulerParameter<RF> explicitEuler;
    Dune::PDELab::HeunParameter<RF> heun;
    Dune::PDELab::RK4Parameter<RF> rk4;
    compareExplicit(check, "CG", oldOperator, newOperator, explicitEuler, gfs, cc);
    compareExplicit(check, "CG", oldOperator, newOperator, heun, gfs, cc);
    compareExplicit(check, "CG", oldOperator, newOperator, rk4, gfs, cc);

    auto xOld = explicitSteps(oldOperator, heun, initial, g);
    auto xNew = explicitSteps(newOperator, heun, initial, g);
    check.vectors("CG Heun three time steps", xOld, xNew);
  }
}

// DG with GridOperator (fast = false) or FastDGGridOperator (fast = true);
// returns the solutions after three Alexander2 steps of both operators
template<bool fast, typename GV>
auto testDG(Checker& check, const GV& gv)
{
  using RF = double;
  const int degree = 1;
  const int dim = GV::dimension;
  using FEM = Dune::PDELab::QkDGLocalFiniteElementMap<typename GV::ctype, RF, degree, dim>;
  FEM fem;
  const int blocksize = Dune::QkStuff::QkSize<degree,dim>::value;
  using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed, blocksize>;
  using GFS = Dune::PDELab::GridFunctionSpace<GV, FEM, Dune::PDELab::NoConstraints, VBE>;
  GFS gfs(gv, fem);

  using P = Problem<GV,RF>;
  P problem;
  using CC = typename GFS::template ConstraintsContainer<RF>::Type;
  CC cc;

  using LOP = std::conditional_t<fast,
                                 Dune::PDELab::ConvectionDiffusionFastDG<P, FEM>,
                                 Dune::PDELab::ConvectionDiffusionDG<P, FEM>>;
  LOP lop(problem, Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,
          Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn, 2.0);
  using MLOP = std::conditional_t<fast, Dune::PDELab::L2FastDG, Dune::PDELab::L2>;
  MLOP mlop(2*degree);

  using MB = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  MB mb(9);
  using GO0 = std::conditional_t<fast,
                                 Dune::PDELab::FastDGGridOperator<GFS,GFS,LOP,MB,RF,RF,RF,CC,CC>,
                                 Dune::PDELab::GridOperator<GFS,GFS,LOP,MB,RF,RF,RF,CC,CC>>;
  using GO1 = std::conditional_t<fast,
                                 Dune::PDELab::FastDGGridOperator<GFS,GFS,MLOP,MB,RF,RF,RF,CC,CC>,
                                 Dune::PDELab::GridOperator<GFS,GFS,MLOP,MB,RF,RF,RF,CC,CC>>;
  GO0 go0Old(gfs, cc, gfs, cc, lop, mb), go0New(gfs, cc, gfs, cc, lop, mb);
  GO1 go1Old(gfs, cc, gfs, cc, mlop, mb), go1New(gfs, cc, gfs, cc, mlop, mb);

  using IGOOld = Dune::PDELab::OneStepGridOperator<GO0,GO1>;
  using IGONew = Dune::PDELab::Impl::GenericOneStepGridOperator<GO0,GO1>;
  IGOOld oldOperator(go0Old, go1Old);
  IGONew newOperator(go0New, go1New);

  const std::string name = fast ? "FastDG" : "DG";
  Dune::PDELab::Alexander2Parameter<RF> alexander2;
  // jacobian_apply of the engine based operator does not compile with FastDG
  // (NumericalJacobianApply writes into the read-only coefficient view)
  compareImplicit<not fast>(check, name, oldOperator, newOperator, alexander2, gfs);

  using G = Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<P>;
  G g(gv, problem);
  using V = Dune::PDELab::Backend::Vector<GFS, RF>;
  V initial(gfs, 0.0);
  problem.setTime(0.0);
  Dune::PDELab::interpolate(g, gfs, initial);
  auto xOld = implicitSteps(oldOperator, alexander2, initial, g);
  auto xNew = implicitSteps(newOperator, alexander2, initial, g);
  check.vectors(name + " Alexander2 three time steps", xOld, xNew);
  // native vectors, the PDELab vectors refer to the local grid function space
  return std::pair(Dune::PDELab::Backend::native(xOld), Dune::PDELab::Backend::native(xNew));
}

int main(int argc, char** argv)
{
  try {
    Dune::MPIHelper::instance(argc, argv);

    const int dim = 2;
    using Grid = Dune::YaspGrid<dim>;
    Grid grid(Dune::FieldVector<double,dim>(1.0), Dune::filledArray<dim,int>(8));
    auto gv = grid.leafGridView();

    Checker check(1e-11);
    testContinuousGalerkin(check, gv);
    auto [dgOld, dgNew] = testDG<false>(check, gv);
    auto [fastOld, fastNew] = testDG<true>(check, gv);
    // the same discretization, assembled by GridOperator and FastDGGridOperator
    check.vectors("DG vs. FastDG (engine based)", dgOld, fastOld);
    check.vectors("DG vs. FastDG (generic)", dgNew, fastNew);

    return check.failed() ? 1 : 0;
  }
  catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
}
