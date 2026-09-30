// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/* Assemble TaylorHoodNavierStokes with the PDELab GridOperator on a
 * composite(power(Q2),Q1) grid function space and with the GridOperator on
 * the corresponding dune-functions composite(power(lagrange<2>),lagrange<1>)
 * basis, and compare residual,
 * sparsity pattern and jacobian. This exercises the adapter on a
 * hierarchic (composite/power) basis tree.
 */

#include "config.h"

#include <iostream>

#include <dune/common/filledarray.hh>
#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/gridfunctionspace/compositegridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/localoperator/stokesparameter.hh>
#include <dune/pdelab/localoperator/taylorhoodnavierstokes.hh>


#include "gridoperatorfunctionstesthelpers.hh"

// Parameters for TaylorHoodNavierStokes with a volume force and a stress
// boundary condition on x=1, so that volume and boundary terms contribute.
template<typename GV, bool navier, bool tensor>
class StokesProblem
{
public:
  static const bool assemble_navier = navier;
  static const bool assemble_full_tensor = tensor;

  using Traits = Dune::PDELab::NavierStokesParameterTraits<GV,double>;

  void setTime(double) {}

  template<typename Element>
  typename Traits::VelocityRange
  f(const Element& e, const typename Traits::Domain& x) const
  {
    auto xg = e.geometry().global(x);
    return {xg[1], 1.0 - xg[0]};
  }

  template<typename EG>
  typename Traits::RangeField
  mu(const EG&, const typename Traits::Domain&) const
  {
    return 1.5;
  }

  template<typename EG>
  typename Traits::RangeField
  rho(const EG&, const typename Traits::Domain&) const
  {
    return 0.8;
  }

  template<typename EG>
  typename Traits::RangeField
  g2(const EG& eg, const typename Traits::Domain& x) const
  {
    auto xg = eg.geometry().global(x);
    return 0.1 * xg[0];
  }

  template<typename IG>
  typename Traits::BoundaryCondition::Type
  bctype(const IG& ig, const typename Traits::IntersectionDomain& x) const
  {
    auto xg = ig.geometry().global(x);
    if (xg[0] > 1.0-1e-8)
      return Traits::BoundaryCondition::StressNeumann;
    return Traits::BoundaryCondition::VelocityDirichlet;
  }

  template<typename IG>
  typename Traits::VelocityRange
  j(const IG& ig, const typename Traits::IntersectionDomain& x, const typename Traits::Domain& normal) const
  {
    auto xg = ig.geometry().global(x);
    auto result = normal;
    result *= 1.0 + xg[1];
    return result;
  }
};


template<bool navier, bool tensor, class GV>
int testTaylorHood(const GV& gv)
{
  using namespace Dune;
  using namespace Dune::Indices;
  using namespace Dune::PDELab::FunctionsBasisTest;
  constexpr int dim = GV::dimension;
  std::cout << "=== TaylorHoodNavierStokes, navier=" << navier << ", full_tensor=" << tensor << std::endl;

  using Param = StokesProblem<GV,navier,tensor>;
  Param param;
  using LOP = PDELab::TaylorHoodNavierStokes<Param>;
  LOP lop(param);

  // PDELab reference: composite(power<dim>(Q2), Q1), flat lexicographic ordering

  using VBE = PDELab::ISTL::VectorBackend<>;
  using V_FEM = PDELab::QkLocalFiniteElementMap<GV,double,double,2>;
  using P_FEM = PDELab::QkLocalFiniteElementMap<GV,double,double,1>;
  V_FEM vfem(gv);
  P_FEM pfem(gv);
  using V_GFS = PDELab::GridFunctionSpace<GV,V_FEM,PDELab::NoConstraints,VBE>;
  using PV_GFS = PDELab::PowerGridFunctionSpace<V_GFS,dim,VBE,PDELab::LexicographicOrderingTag>;
  using P_GFS = PDELab::GridFunctionSpace<GV,P_FEM,PDELab::NoConstraints,VBE>;
  using GFS = PDELab::CompositeGridFunctionSpace<VBE,PDELab::LexicographicOrderingTag,PV_GFS,P_GFS>;
  V_GFS vgfs(gv,vfem);
  PV_GFS pvgfs(vgfs);
  P_GFS pgfs(gv,pfem);
  GFS gfs(pvgfs,pgfs);

  using MBE = PDELab::ISTL::BCRSMatrixBackend<>;
  using GO = PDELab::GridOperator<GFS,GFS,LOP,MBE,double,double,double,
                                  PDELab::EmptyTransformation,PDELab::EmptyTransformation>;
  GO go(gfs,gfs,lop,MBE(50));

  // dune-functions Taylor-Hood basis with the same tree structure but a different
  // index layout (interleaved velocity components), flat multi-indices.
  // The very same local operator type is used.

  using namespace Functions::BasisFactory;
  auto basis = makeBasis(gv,
    composite(
      power<dim>(lagrange<2>(), flatInterleaved()),
      lagrange<1>(),
      flatLexicographic()));
  using Basis = decltype(basis);
  using FGO = PDELab::GridOperator<Basis,Basis,LOP,MBE,double,double,double>;
  FGO fgo(basis, basis, lop, go.matrixBackend());

  using Vector = typename FGO::Traits::Domain;
  using Matrix = typename FGO::Traits::Jacobian;

  // match DOFs: keys (x+10,y) for u_0, (y+20,x) for u_1 and (x+30,y) for p

  using PVector = typename GO::Traits::Domain;
  using Range = FieldVector<double,dim+1>;
  PVector key1(gfs), key2(gfs);
  PDELab::interpolate(PDELab::makeGridFunctionFromCallable(gv, [](const auto& x) { return Range{x[0]+10, x[1]+20, x[0]+30}; }), gfs, key1);
  PDELab::interpolate(PDELab::makeGridFunctionFromCallable(gv, [](const auto& x) { return Range{x[1], x[0], x[1]}; }), gfs, key2);

  Vector fkey1(basis, 0.0), fkey2(basis, 0.0);
  auto interpolateKeys = [&](Vector& v, auto u0, auto u1, auto p) {
    Functions::interpolate(Functions::subspaceBasis(basis, _0, 0), PDELab::Backend::native(v), u0);
    Functions::interpolate(Functions::subspaceBasis(basis, _0, 1), PDELab::Backend::native(v), u1);
    Functions::interpolate(Functions::subspaceBasis(basis, _1), PDELab::Backend::native(v), p);
  };
  interpolateKeys(fkey1,
    [](const auto& x) { return x[0]+10; }, [](const auto& x) { return x[1]+20; }, [](const auto& x) { return x[0]+30; });
  interpolateKeys(fkey2,
    [](const auto& x) { return x[1]; }, [](const auto& x) { return x[0]; }, [](const auto& x) { return x[1]; });

  auto perm = matchDofs(flatValues(PDELab::Backend::native(key1)), flatValues(PDELab::Backend::native(key2)),
                        scalarValues(PDELab::Backend::native(fkey1)), scalarValues(PDELab::Backend::native(fkey2)));

  Checker checker(1e-10);

  PVector x(gfs);
  randomize(PDELab::Backend::native(x), 42);
  Vector fx(basis);
  transfer(PDELab::Backend::native(x), perm, PDELab::Backend::native(fx));

  // residual

  PVector r(gfs, 0.0);
  go.residual(x, r);
  Vector fr(basis, 0.0);
  fgo.residual(fx, fr);
  checker.compareVectors(PDELab::Backend::native(r), PDELab::Backend::native(fr), perm, "residual");

  // pattern and jacobian (depends on x for Navier-Stokes)

  using PMatrix = typename GO::Traits::Jacobian;
  PMatrix a(go, 0.0);
  go.jacobian(x, a);

  Matrix fa(fgo, 0.0);
  fgo.jacobian(fx, fa);
  checker.compareMatrices(PDELab::Backend::native(a), PDELab::Backend::native(fa), perm, "jacobian");
  checker.check(PDELab::Backend::native(a).nonzeroes() == PDELab::Backend::native(fa).nonzeroes(), "number of pattern entries");

  return checker.result();
}


int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  using Grid = Dune::YaspGrid<2>;
  Grid grid({1.0,1.0}, Dune::filledArray<2,int>(4));
  auto gv = grid.leafGridView();

  int result = 0;
  result |= testTaylorHood<false,false>(gv);
  result |= testTaylorHood<true,true>(gv);
  return result;
}
