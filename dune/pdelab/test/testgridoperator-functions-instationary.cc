// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/* Compares the OneStepGridOperator on a dune-functions basis with the
 * OneStepGridOperator on a PDELab grid function space spanning the same
 * discrete space: residual, jacobian and jacobian_apply of every stage of an
 * implicit method, explicit_jacobian_residual of an explicit method, and
 * complete time steps. Continuous Galerkin with Dirichlet constraints and
 * DG with skeleton terms, for a time dependent convection-diffusion problem.
 *
 * Multi-stage methods assemble the constant part of a stage (all previous
 * stages) in one grid traversal; this is compared with the stage residual
 * assembled by hand from the spatial and temporal grid operators, also for
 * FastDGGridOperator.
 */

#include "config.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <dune/common/filledarray.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/pdelab.hh>

#include "convectiondiffusionfastdg.hh"
#include "gridoperatorfunctionstesthelpers.hh"
#include "l2fastdg.hh"

using namespace Dune;
using namespace Dune::PDELab::FunctionsBasisTest;

// Time dependent convection-diffusion problem with Dirichlet boundary conditions
template<typename GV, typename RF>
class Problem
{
  typedef PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

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
    return PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
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

// Both one step operators, the grid function space and the basis, and the DOF permutation
template<typename IGO, typename FIGO, typename GFS, typename Basis>
struct Pair
{
  IGO& igo;
  FIGO& figo;
  const GFS& gfs;
  const Basis& basis;
  const std::vector<std::size_t>& perm;
};

template<typename P>
void compareImplicit(Checker& checker, const std::string& name, P p,
                     const PDELab::TimeSteppingParameterInterface<double>& method)
{
  using V = typename std::decay_t<decltype(p.igo)>::Traits::Domain;
  using FV = typename std::decay_t<decltype(p.figo)>::Traits::Domain;
  using J = typename std::decay_t<decltype(p.igo)>::Traits::Jacobian;
  using FJ = typename std::decay_t<decltype(p.figo)>::Traits::Jacobian;

  p.igo.preStep(method, 0.3, 0.1);
  p.figo.preStep(method, 0.3, 0.1);

  std::vector<V> x(method.s(), V(p.gfs, 0.0));
  std::vector<FV> fx(method.s(), FV(p.basis, 0.0));
  std::vector<V*> xp;
  std::vector<FV*> fxp;
  for (unsigned int s = 0; s < method.s(); ++s)
  {
    randomize(PDELab::Backend::native(x[s]), s + 1);
    transfer(PDELab::Backend::native(x[s]), p.perm, PDELab::Backend::native(fx[s]));
    xp.push_back(&x[s]);
    fxp.push_back(&fx[s]);
  }

  for (unsigned int r = 1; r <= method.s(); ++r)
  {
    const std::string stage = name + " " + method.name() + " stage " + std::to_string(r) + " ";
    std::vector<V*> previous(xp.begin(), xp.begin() + r);
    std::vector<FV*> fprevious(fxp.begin(), fxp.begin() + r);
    p.igo.preStage(r, previous);
    p.figo.preStage(r, fprevious);

    V z(p.gfs);
    FV fz(p.basis);
    randomize(PDELab::Backend::native(z), 10 + r);
    transfer(PDELab::Backend::native(z), p.perm, PDELab::Backend::native(fz));

    V res(p.gfs, 0.0);
    FV fres(p.basis, 0.0);
    p.igo.residual(z, res);
    p.figo.residual(fz, fres);
    checker.compareVectors(PDELab::Backend::native(res), PDELab::Backend::native(fres), p.perm, stage + "residual");

    J a(p.igo);
    FJ fa(p.figo);
    p.igo.jacobian(z, a);
    p.figo.jacobian(fz, fa);
    checker.compareMatrices(PDELab::Backend::native(a), PDELab::Backend::native(fa), p.perm, stage + "jacobian");

    V y(p.gfs, 0.0);
    FV fy(p.basis, 0.0);
    p.igo.jacobian_apply(z, y);
    p.figo.jacobian_apply(fz, fy);
    checker.compareVectors(PDELab::Backend::native(y), PDELab::Backend::native(fy), p.perm, stage + "jacobian_apply");

    p.igo.postStage();
    p.figo.postStage();
  }
  p.igo.postStep();
  p.figo.postStep();
}

template<typename P>
void compareExplicit(Checker& checker, const std::string& name, P p,
                     const PDELab::TimeSteppingParameterInterface<double>& method)
{
  using V = typename std::decay_t<decltype(p.igo)>::Traits::Domain;
  using FV = typename std::decay_t<decltype(p.figo)>::Traits::Domain;
  using J = typename std::decay_t<decltype(p.igo)>::Traits::Jacobian;
  using FJ = typename std::decay_t<decltype(p.figo)>::Traits::Jacobian;

  p.igo.preStep(method, 0.3, 0.1);
  p.figo.preStep(method, 0.3, 0.1);

  std::vector<V> x(method.s() + 1, V(p.gfs, 0.0));
  std::vector<FV> fx(method.s() + 1, FV(p.basis, 0.0));
  std::vector<V*> xp;
  std::vector<FV*> fxp;
  for (unsigned int s = 0; s <= method.s(); ++s)
  {
    randomize(PDELab::Backend::native(x[s]), 20 + s);
    transfer(PDELab::Backend::native(x[s]), p.perm, PDELab::Backend::native(fx[s]));
    xp.push_back(&x[s]);
    fxp.push_back(&fx[s]);
  }

  for (unsigned int r = 1; r <= method.s(); ++r)
  {
    const std::string stage = name + " " + method.name() + " stage " + std::to_string(r) + " ";
    std::vector<V*> upToStage(xp.begin(), xp.begin() + r + 1);
    std::vector<FV*> fupToStage(fxp.begin(), fxp.begin() + r + 1);
    J d(p.igo);
    FJ fd(p.figo);
    V alpha(p.gfs, 0.0), beta(p.gfs, 0.0);
    FV falpha(p.basis, 0.0), fbeta(p.basis, 0.0);
    p.igo.explicit_jacobian_residual(r, upToStage, d, alpha, beta);
    p.figo.explicit_jacobian_residual(r, fupToStage, fd, falpha, fbeta);
    checker.compareMatrices(PDELab::Backend::native(d), PDELab::Backend::native(fd), p.perm, stage + "D");
    checker.compareVectors(PDELab::Backend::native(alpha), PDELab::Backend::native(falpha), p.perm, stage + "alpha");
    checker.compareVectors(PDELab::Backend::native(beta), PDELab::Backend::native(fbeta), p.perm, stage + "beta");
    p.igo.postStage();
    p.figo.postStage();
  }
  p.igo.postStep();
  p.figo.postStep();
}

// three implicit time steps of 0.05, starting from x
template<typename IGO, typename V, typename G>
void implicitSteps(IGO& igo, const PDELab::TimeSteppingParameterInterface<double>& method, V& x, G& g)
{
  using LinearSolver = PDELab::ISTLBackend_SEQ_SuperLU;
  LinearSolver linearSolver(false);
  using Solver = PDELab::StationaryLinearProblemSolver<IGO, LinearSolver, V>;
  Solver solver(igo, linearSolver, 1e-12);
  PDELab::OneStepMethod<double, IGO, Solver, V, V> osm(method, igo, solver);
  osm.setVerbosityLevel(0);
  double time = 0.0;
  for (int step = 0; step < 3; ++step)
  {
    V xnew(x);
    osm.apply(time, 0.05, x, g, xnew);
    x = xnew;
    time += 0.05;
  }
}

// The stage residual c_r + b_rr dt A(z,t_r) + M(z,t_r) of every stage, assembled
// by hand from the spatial (go0) and temporal (go1) grid operators, compared
// with the residual of the one step operator after preStage. Checks that the
// constant part of all previous stages is assembled correctly (in one traversal).
template<typename IGO, typename GO0, typename GO1, typename Space>
void compareWithHandAssembly(Checker& checker, const std::string& name, IGO& igo, const GO0& go0, const GO1& go1,
                             const PDELab::TimeSteppingParameterInterface<double>& method, const Space& space)
{
  using V = typename IGO::Traits::Domain;
  auto& lop0 = go0.localAssembler().localOperator();
  auto& lop1 = go1.localAssembler().localOperator();
  const double t = 0.3, dt = 0.1;
  igo.preStep(method, t, dt);

  std::vector<V> x(method.s(), V(space, 0.0));
  std::vector<V*> xp;
  for (unsigned int s = 0; s < method.s(); ++s)
  {
    randomize(PDELab::Backend::native(x[s]), 30 + s);
    xp.push_back(&x[s]);
  }
  V z(space);
  randomize(PDELab::Backend::native(z), 40);
  const auto n = flatValues(PDELab::Backend::native(z)).size();
  std::vector<std::size_t> identity(n);
  for (std::size_t i = 0; i < n; ++i)
    identity[i] = i;

  for (unsigned int r = 1; r <= method.s(); ++r)
  {
    V reference(space, 0.0), part(space, 0.0);
    auto add = [&](const V& u, double time, double w0, double w1) {
      lop0.setTime(time);
      lop1.setTime(time);
      part = 0.0;
      go0.residual(u, part);
      reference.axpy(w0, part);
      part = 0.0;
      go1.residual(u, part);
      reference.axpy(w1, part);
    };
    for (unsigned int s = 0; s < r; ++s)
      add(x[s], t + method.d(s)*dt, method.b(r,s)*dt, method.a(r,s));
    add(z, t + method.d(r)*dt, method.b(r,r)*dt, 1.0);

    std::vector<V*> previous(xp.begin(), xp.begin() + r);
    igo.preStage(r, previous);
    V residual(space, 0.0);
    igo.residual(z, residual);
    auto flat = flatValues(PDELab::Backend::native(residual));
    checker.compareVectors(PDELab::Backend::native(reference), flat, identity,
                           name + " " + method.name() + " stage " + std::to_string(r) + " residual vs. hand assembly");
    igo.postStage();
  }
  igo.postStep();
}

template<typename GV>
void testContinuousGalerkin(Checker& checker, const GV& gv)
{
  std::cout << "=== CG, Q2, Dirichlet constraints" << std::endl;
  using P = Problem<GV,double>;
  P problem;
  PDELab::ConvectionDiffusionBoundaryConditionAdapter<P> bctype(gv, problem);
  using G = PDELab::ConvectionDiffusionDirichletExtensionAdapter<P>;
  G g(gv, problem);
  using MB = PDELab::ISTL::BCRSMatrixBackend<>;
  MB mb(25);

  // grid function space
  using FEM = PDELab::QkLocalFiniteElementMap<GV,double,double,2>;
  FEM fem(gv);
  using GFS = PDELab::GridFunctionSpace<GV,FEM,PDELab::ConformingDirichletConstraints,PDELab::ISTL::VectorBackend<>>;
  GFS gfs(gv, fem);
  using CC = typename GFS::template ConstraintsContainer<double>::Type;
  CC cc;
  PDELab::constraints(bctype, gfs, cc);
  using LOP = PDELab::ConvectionDiffusionFEM<P,FEM>;
  LOP lop(problem);
  PDELab::L2 mlop;
  using GO0 = PDELab::GridOperator<GFS,GFS,LOP,MB,double,double,double,CC,CC>;
  using GO1 = PDELab::GridOperator<GFS,GFS,PDELab::L2,MB,double,double,double,CC,CC>;
  GO0 go0(gfs, cc, gfs, cc, lop, mb);
  GO1 go1(gfs, cc, gfs, cc, mlop, mb);

  // basis
  auto basis = Functions::BasisFactory::makeBasis(gv, Functions::BasisFactory::lagrange<2>());
  using Basis = decltype(basis);
  using FCC = PDELab::BasisConstraintsContainer<Basis,double>;
  FCC fcc;
  PDELab::constraints(bctype, basis, fcc);
  using FLOP = PDELab::ConvectionDiffusionFEM<P,PDELab::FiniteElementMapFor<typename Basis::LocalView::Tree::FiniteElement>>;
  FLOP flop(problem);
  using FGO0 = PDELab::GridOperator<Basis,Basis,FLOP,MB,double,double,double,FCC,FCC>;
  using FGO1 = PDELab::GridOperator<Basis,Basis,PDELab::L2,MB,double,double,double,FCC,FCC>;
  FGO0 fgo0(basis, fcc, basis, fcc, flop, mb);
  FGO1 fgo1(basis, fcc, basis, fcc, mlop, mb);

  // DOF permutation from the coordinates of the Lagrange nodes
  using V = PDELab::Backend::Vector<GFS,double>;
  using FV = PDELab::Backend::Vector<Basis,double>;
  V key1(gfs), key2(gfs);
  PDELab::interpolate(PDELab::makeGridFunctionFromCallable(gv, [](const auto& x) { return x[0]; }), gfs, key1);
  PDELab::interpolate(PDELab::makeGridFunctionFromCallable(gv, [](const auto& x) { return x[1]; }), gfs, key2);
  FV fkey1(basis, 0.0), fkey2(basis, 0.0);
  Functions::interpolate(basis, PDELab::Backend::native(fkey1), [](const auto& x) { return x[0]; });
  Functions::interpolate(basis, PDELab::Backend::native(fkey2), [](const auto& x) { return x[1]; });
  auto perm = matchDofs(flatValues(PDELab::Backend::native(key1)), flatValues(PDELab::Backend::native(key2)),
                        scalarValues(PDELab::Backend::native(fkey1)), scalarValues(PDELab::Backend::native(fkey2)));

  {
    PDELab::OneStepGridOperator<GO0,GO1> igo(go0, go1);
    PDELab::OneStepGridOperator<FGO0,FGO1> figo(fgo0, fgo1);
    Pair<decltype(igo),decltype(figo),GFS,Basis> p{igo, figo, gfs, basis, perm};
    PDELab::OneStepThetaParameter<double> theta(0.5);
    PDELab::Alexander2Parameter<double> alexander2;
    PDELab::Alexander3Parameter<double> alexander3;
    PDELab::FractionalStepParameter<double> fractionalStep;
    compareImplicit(checker, "CG", p, theta);
    compareImplicit(checker, "CG", p, alexander2);
    compareImplicit(checker, "CG", p, alexander3);
    compareImplicit(checker, "CG", p, fractionalStep);
    compareWithHandAssembly(checker, "CG grid function space", igo, go0, go1, alexander3, gfs);
    compareWithHandAssembly(checker, "CG basis", figo, fgo0, fgo1, alexander3, basis);

    V x(gfs, 0.0);
    FV fx(basis, 0.0);
    problem.setTime(0.0);
    PDELab::interpolate(g, gfs, x);
    PDELab::interpolate(g, basis, fx);
    implicitSteps(igo, alexander2, x, g);
    implicitSteps(figo, alexander2, fx, g);
    checker.compareVectors(PDELab::Backend::native(x), PDELab::Backend::native(fx), perm, "CG Alexander2 three time steps");
  }
  {
    PDELab::OneStepGridOperator<GO0,GO1,false> igo(go0, go1);
    PDELab::OneStepGridOperator<FGO0,FGO1,false> figo(fgo0, fgo1);
    Pair<decltype(igo),decltype(figo),GFS,Basis> p{igo, figo, gfs, basis, perm};
    PDELab::HeunParameter<double> heun;
    PDELab::RK4Parameter<double> rk4;
    compareExplicit(checker, "CG", p, heun);
    compareExplicit(checker, "CG", p, rk4);
  }
}

template<typename GV>
void testDG(Checker& checker, const GV& gv)
{
  std::cout << "=== DG, Q1, weak Dirichlet conditions" << std::endl;
  const int dim = GV::dimension;
  using P = Problem<GV,double>;
  P problem;
  using G = PDELab::ConvectionDiffusionDirichletExtensionAdapter<P>;
  G g(gv, problem);
  using MB = PDELab::ISTL::BCRSMatrixBackend<>;
  MB mb(20);
  const auto method = PDELab::ConvectionDiffusionDGMethod::SIPG;
  const auto weights = PDELab::ConvectionDiffusionDGWeights::weightsOn;

  // grid function space
  using FEM = PDELab::QkDGLocalFiniteElementMap<typename GV::ctype,double,1,dim>;
  FEM fem;
  using GFS = PDELab::GridFunctionSpace<GV,FEM,PDELab::NoConstraints,PDELab::ISTL::VectorBackend<>>;
  GFS gfs(gv, fem);
  using LOP = PDELab::ConvectionDiffusionDG<P,FEM>;
  LOP lop(problem, method, weights, 2.0);
  PDELab::L2 mlop;
  using GO0 = PDELab::GridOperator<GFS,GFS,LOP,MB,double,double,double>;
  using GO1 = PDELab::GridOperator<GFS,GFS,PDELab::L2,MB,double,double,double>;
  GO0 go0(gfs, gfs, lop, mb);
  GO1 go1(gfs, gfs, mlop, mb);

  // basis
  auto basis = Functions::BasisFactory::makeBasis(gv, Functions::BasisFactory::lagrangeDG<1>());
  using Basis = decltype(basis);
  using FLOP = PDELab::ConvectionDiffusionDG<P,PDELab::FiniteElementMapFor<typename Basis::LocalView::Tree::FiniteElement>>;
  FLOP flop(problem, method, weights, 2.0);
  using FGO0 = PDELab::GridOperator<Basis,Basis,FLOP,MB,double,double,double>;
  using FGO1 = PDELab::GridOperator<Basis,Basis,PDELab::L2,MB,double,double,double>;
  FGO0 fgo0(basis, basis, flop, mb);
  FGO1 fgo1(basis, basis, mlop, mb);

  // DOF permutation, element by element
  using V = PDELab::Backend::Vector<GFS,double>;
  using FV = PDELab::Backend::Vector<Basis,double>;
  V key1(gfs), key2(gfs);
  PDELab::interpolate(PDELab::makeGridFunctionFromCallable(gv, [](const auto& x) { return x[0]; }), gfs, key1);
  PDELab::interpolate(PDELab::makeGridFunctionFromCallable(gv, [](const auto& x) { return x[1]; }), gfs, key2);
  FV fkey1(basis, 0.0), fkey2(basis, 0.0);
  Functions::interpolate(basis, PDELab::Backend::native(fkey1), [](const auto& x) { return x[0]; });
  Functions::interpolate(basis, PDELab::Backend::native(fkey2), [](const auto& x) { return x[1]; });
  auto perm = matchDofsElementwise(gfs, basis,
                                   flatValues(PDELab::Backend::native(key1)), flatValues(PDELab::Backend::native(key2)),
                                   scalarValues(PDELab::Backend::native(fkey1)), scalarValues(PDELab::Backend::native(fkey2)));

  PDELab::OneStepGridOperator<GO0,GO1> igo(go0, go1);
  PDELab::OneStepGridOperator<FGO0,FGO1> figo(fgo0, fgo1);
  Pair<decltype(igo),decltype(figo),GFS,Basis> p{igo, figo, gfs, basis, perm};
  PDELab::Alexander2Parameter<double> alexander2;
  PDELab::Alexander3Parameter<double> alexander3;
  compareImplicit(checker, "DG", p, alexander2);
  compareImplicit(checker, "DG", p, alexander3);
  compareWithHandAssembly(checker, "DG grid function space", igo, go0, go1, alexander3, gfs);
  compareWithHandAssembly(checker, "DG basis", figo, fgo0, fgo1, alexander3, basis);

  V x(gfs, 0.0);
  FV fx(basis, 0.0);
  problem.setTime(0.0);
  PDELab::interpolate(g, gfs, x);
  PDELab::interpolate(g, basis, fx);
  implicitSteps(igo, alexander2, x, g);
  implicitSteps(figo, alexander2, fx, g);
  checker.compareVectors(PDELab::Backend::native(x), PDELab::Backend::native(fx), perm, "DG Alexander2 three time steps");
}

// FastDG and GridOperator assemble the same DG discretization on a grid function space
template<typename GV>
void testFastDG(Checker& checker, const GV& gv)
{
  std::cout << "=== FastDG vs. GridOperator, DG Q1" << std::endl;
  const int dim = GV::dimension;
  using P = Problem<GV,double>;
  P problem;
  using G = PDELab::ConvectionDiffusionDirichletExtensionAdapter<P>;
  G g(gv, problem);
  using MB = PDELab::ISTL::BCRSMatrixBackend<>;
  MB mb(9);
  using FEM = PDELab::QkDGLocalFiniteElementMap<typename GV::ctype,double,1,dim>;
  FEM fem;
  const int blocksize = QkStuff::QkSize<1,dim>::value;
  using VBE = PDELab::ISTL::VectorBackend<PDELab::ISTL::Blocking::fixed, blocksize>;
  using GFS = PDELab::GridFunctionSpace<GV,FEM,PDELab::NoConstraints,VBE>;
  GFS gfs(gv, fem);
  const auto method = PDELab::ConvectionDiffusionDGMethod::SIPG;
  const auto weights = PDELab::ConvectionDiffusionDGWeights::weightsOn;

  using LOP = PDELab::ConvectionDiffusionDG<P,FEM>;
  LOP lop(problem, method, weights, 2.0);
  PDELab::L2 mlop;
  using FastLOP = PDELab::ConvectionDiffusionFastDG<P,FEM>;
  FastLOP fastLop(problem, method, weights, 2.0);
  PDELab::L2FastDG fastMlop(2);
  using GO0 = PDELab::GridOperator<GFS,GFS,LOP,MB,double,double,double>;
  using GO1 = PDELab::GridOperator<GFS,GFS,PDELab::L2,MB,double,double,double>;
  using FastGO0 = PDELab::FastDGGridOperator<GFS,GFS,FastLOP,MB,double,double,double,PDELab::EmptyTransformation,PDELab::EmptyTransformation>;
  using FastGO1 = PDELab::FastDGGridOperator<GFS,GFS,PDELab::L2FastDG,MB,double,double,double,PDELab::EmptyTransformation,PDELab::EmptyTransformation>;
  GO0 go0(gfs, gfs, lop, mb);
  GO1 go1(gfs, gfs, mlop, mb);
  FastGO0 fastGo0(gfs, gfs, fastLop, mb);
  FastGO1 fastGo1(gfs, gfs, fastMlop, mb);
  PDELab::OneStepGridOperator<GO0,GO1> igo(go0, go1);
  PDELab::OneStepGridOperator<FastGO0,FastGO1> fastIgo(fastGo0, fastGo1);

  PDELab::Alexander3Parameter<double> alexander3;
  compareWithHandAssembly(checker, "FastDG", fastIgo, fastGo0, fastGo1, alexander3, gfs);

  using V = PDELab::Backend::Vector<GFS,double>;
  V x(gfs, 0.0), fastX(gfs, 0.0);
  problem.setTime(0.0);
  PDELab::interpolate(g, gfs, x);
  fastX = x;
  implicitSteps(igo, alexander3, x, g);
  implicitSteps(fastIgo, alexander3, fastX, g);
  const auto values = flatValues(PDELab::Backend::native(x));
  std::vector<std::size_t> identity(values.size());
  for (std::size_t i = 0; i < values.size(); ++i)
    identity[i] = i;
  checker.compareVectors(PDELab::Backend::native(x), flatValues(PDELab::Backend::native(fastX)), identity,
                         "FastDG vs. GridOperator Alexander3 three time steps");
}

int main(int argc, char** argv)
{
  try {
    MPIHelper::instance(argc, argv);

    const int dim = 2;
    using Grid = YaspGrid<dim>;
    Grid grid(FieldVector<double,dim>(1.0), filledArray<dim,int>(8));
    auto gv = grid.leafGridView();

    Checker checker(1e-10);
    testContinuousGalerkin(checker, gv);
    testDG(checker, gv);
    testFastDG(checker, gv);
    return checker.result();
  }
  catch (Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
}
