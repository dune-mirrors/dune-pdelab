// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/* Assemble ConvectionDiffusionFEM (continuous Lagrange) and
 * ConvectionDiffusionDG (discontinuous Lagrange, with skeleton terms) with
 * the GridOperator on PDELab grid function spaces and with the GridOperator
 * on the corresponding dune-functions bases, and compare residual,
 * jacobian, jacobian application and sparsity pattern.
 */

#include "config.h"

#include <iostream>

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/finiteelementmap/qkdg.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>


#include "gridoperatorfunctionstesthelpers.hh"

// Parameters exercising volume and boundary terms of ConvectionDiffusionFEM:
// anisotropic diffusion, convection, reaction, source, Neumann (x=1) and
// outflow (y=1) boundaries.
template<typename GV, typename RF>
class Problem
{
  using BCType = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type;

public:
  using Traits = Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF>;

  static constexpr bool permeabilityIsConstantPerCell() { return true; }

  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType a(0.0);
    a[0][0] = 1.0; a[1][1] = 2.0; a[0][1] = a[1][0] = 0.25;
    return a;
  }

  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return {0.3, 0.7};
  }

  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.5;
  }

  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    auto xg = e.geometry().global(x);
    return 1.0 + xg[0]*xg[1];
  }

  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    auto xg = is.geometry().global(x);
    if (xg[0] > 1.0-1e-8)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    if (xg[1] > 1.0-1e-8)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow;
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    auto xg = is.geometry().global(x);
    return 1.0 + xg[1];
  }

  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }
};


template<int k, class GV>
int testPoisson(const GV& gv)
{
  using namespace Dune;
  using namespace Dune::PDELab::FunctionsBasisTest;
  std::cout << "=== ConvectionDiffusionFEM, Q" << k << std::endl;

  using Param = Problem<GV,double>;
  Param param;

  // PDELab reference

  using FEM = PDELab::QkLocalFiniteElementMap<GV,double,double,k>;
  FEM fem(gv);
  using GFS = PDELab::GridFunctionSpace<GV,FEM,PDELab::NoConstraints,PDELab::ISTL::VectorBackend<>>;
  GFS gfs(gv,fem);

  using LOP = PDELab::ConvectionDiffusionFEM<Param,FEM>;
  LOP lop(param);
  using MBE = PDELab::ISTL::BCRSMatrixBackend<>;
  using GO = PDELab::GridOperator<GFS,GFS,LOP,MBE,double,double,double,
                                  PDELab::EmptyTransformation,PDELab::EmptyTransformation>;
  GO go(gfs,gfs,lop,MBE(9));

  // dune-functions basis with the same discrete space; the local operator
  // source is reused, only the finite element map parameter differs

  auto basis = Functions::BasisFactory::makeBasis(gv, Functions::BasisFactory::lagrange<k>());
  using Basis = decltype(basis);
  using FunctionsLOP = PDELab::ConvectionDiffusionFEM<Param,
    PDELab::FiniteElementMapFor<typename Basis::LocalView::Tree::FiniteElement>>;
  FunctionsLOP functionsLop(param);
  using FGO = PDELab::GridOperator<Basis,Basis,FunctionsLOP,MBE,double,double,double>;
  FGO fgo(basis, basis, functionsLop, go.matrixBackend());

  using Vector = typename FGO::Traits::Domain;
  using Matrix = typename FGO::Traits::Jacobian;

  // match DOFs via the coordinates of the Lagrange nodes

  using PVector = typename GO::Traits::Domain;
  PVector key1(gfs), key2(gfs);
  PDELab::interpolate(PDELab::makeGridFunctionFromCallable(gv, [](const auto& x) { return x[0]; }), gfs, key1);
  PDELab::interpolate(PDELab::makeGridFunctionFromCallable(gv, [](const auto& x) { return x[1]; }), gfs, key2);
  Vector fkey1(basis, 0.0), fkey2(basis, 0.0);
  Functions::interpolate(basis, PDELab::Backend::native(fkey1), [](const auto& x) { return x[0]; });
  Functions::interpolate(basis, PDELab::Backend::native(fkey2), [](const auto& x) { return x[1]; });
  auto perm = matchDofs(flatValues(PDELab::Backend::native(key1)), flatValues(PDELab::Backend::native(key2)),
                        scalarValues(PDELab::Backend::native(fkey1)), scalarValues(PDELab::Backend::native(fkey2)));

  Checker checker(1e-10);

  // evaluation point and direction

  PVector x(gfs), z(gfs);
  randomize(PDELab::Backend::native(x), 42);
  randomize(PDELab::Backend::native(z), 4711);
  Vector fx(basis), fz(basis);
  transfer(PDELab::Backend::native(x), perm, PDELab::Backend::native(fx));
  transfer(PDELab::Backend::native(z), perm, PDELab::Backend::native(fz));

  // residual

  PVector r(gfs, 0.0);
  go.residual(x, r);
  Vector fr(basis, 0.0);
  fgo.residual(fx, fr);
  checker.compareVectors(PDELab::Backend::native(r), PDELab::Backend::native(fr), perm, "residual");

  // pattern and jacobian

  using PMatrix = typename GO::Traits::Jacobian;
  PMatrix a(go, 0.0);
  go.jacobian(x, a);

  Matrix fa(fgo, 0.0);
  fgo.jacobian(fx, fa);
  checker.compareMatrices(PDELab::Backend::native(a), PDELab::Backend::native(fa), perm, "jacobian");
  checker.check(PDELab::Backend::native(a).nonzeroes() == PDELab::Backend::native(fa).nonzeroes(), "number of pattern entries");

  // jacobian application

  PVector y(gfs, 0.0);
  go.jacobian_apply(z, y);
  Vector fy(basis, 0.0);
  fgo.jacobian_apply(fz, fy);
  checker.compareVectors(PDELab::Backend::native(y), PDELab::Backend::native(fy), perm, "jacobian_apply");

  // consistency: the assembled matrix applied to z. ConvectionDiffusionFEM
  // assembles the jacobian analytically but applies it by finite differences
  // (NumericalJacobianApply, epsilon=1e-7), hence the larger tolerance.
  auto fAz = PDELab::Backend::native(fy);
  PDELab::Backend::native(fa).mv(PDELab::Backend::native(fz), fAz);
  fAz -= PDELab::Backend::native(fy);
  std::cout << "|A z - jacobian_apply(z)|_inf = " << fAz.infinity_norm() << std::endl;
  checker.check(fAz.infinity_norm() < 1e-6 * PDELab::Backend::native(fy).infinity_norm(), "matrix-vector product vs. jacobian_apply");

  return checker.result();
}


template<int k, class GV>
int testDG(const GV& gv)
{
  using namespace Dune;
  using namespace Dune::PDELab::FunctionsBasisTest;
  constexpr int dim = GV::dimension;
  std::cout << "=== ConvectionDiffusionDG, Q" << k << "-DG" << std::endl;

  using Param = Problem<GV,double>;
  Param param;

  // PDELab reference

  using FEM = PDELab::QkDGLocalFiniteElementMap<typename GV::ctype,double,k,dim>;
  FEM fem;
  using GFS = PDELab::GridFunctionSpace<GV,FEM,PDELab::NoConstraints,PDELab::ISTL::VectorBackend<>>;
  GFS gfs(gv,fem);

  using LOP = PDELab::ConvectionDiffusionDG<Param,FEM>;
  LOP lop(param);
  using MBE = PDELab::ISTL::BCRSMatrixBackend<>;
  using GO = PDELab::GridOperator<GFS,GFS,LOP,MBE,double,double,double,
                                  PDELab::EmptyTransformation,PDELab::EmptyTransformation>;
  GO go(gfs,gfs,lop,MBE(5*Dune::power(k+1,dim)));

  // dune-functions DG basis, the same local operator with a different finite element map parameter

  auto basis = Functions::BasisFactory::makeBasis(gv, Functions::BasisFactory::lagrangeDG<k>());
  using Basis = decltype(basis);
  using FunctionsLOP = PDELab::ConvectionDiffusionDG<Param,
    PDELab::FiniteElementMapFor<typename Basis::LocalView::Tree::FiniteElement>>;
  FunctionsLOP functionsLop(param);
  using FGO = PDELab::GridOperator<Basis,Basis,FunctionsLOP,MBE,double,double,double>;
  FGO fgo(basis, basis, functionsLop, go.matrixBackend());

  using Vector = typename FGO::Traits::Domain;
  using Matrix = typename FGO::Traits::Jacobian;

  // match DOFs element by element via the coordinates of the Lagrange nodes

  using PVector = typename GO::Traits::Domain;
  PVector key1(gfs), key2(gfs);
  PDELab::interpolate(PDELab::makeGridFunctionFromCallable(gv, [](const auto& x) { return x[0]; }), gfs, key1);
  PDELab::interpolate(PDELab::makeGridFunctionFromCallable(gv, [](const auto& x) { return x[1]; }), gfs, key2);
  Vector fkey1(basis, 0.0), fkey2(basis, 0.0);
  Functions::interpolate(basis, PDELab::Backend::native(fkey1), [](const auto& x) { return x[0]; });
  Functions::interpolate(basis, PDELab::Backend::native(fkey2), [](const auto& x) { return x[1]; });
  auto perm = matchDofsElementwise(gfs, basis,
                                   flatValues(PDELab::Backend::native(key1)), flatValues(PDELab::Backend::native(key2)),
                                   scalarValues(PDELab::Backend::native(fkey1)), scalarValues(PDELab::Backend::native(fkey2)));

  Checker checker(1e-10);

  PVector x(gfs), z(gfs);
  randomize(PDELab::Backend::native(x), 42);
  randomize(PDELab::Backend::native(z), 4711);
  Vector fx(basis), fz(basis);
  transfer(PDELab::Backend::native(x), perm, PDELab::Backend::native(fx));
  transfer(PDELab::Backend::native(z), perm, PDELab::Backend::native(fz));

  // residual

  PVector r(gfs, 0.0);
  go.residual(x, r);
  Vector fr(basis, 0.0);
  fgo.residual(fx, fr);
  checker.compareVectors(PDELab::Backend::native(r), PDELab::Backend::native(fr), perm, "residual");

  // pattern and jacobian

  using PMatrix = typename GO::Traits::Jacobian;
  PMatrix a(go, 0.0);
  go.jacobian(x, a);

  Matrix fa(fgo, 0.0);
  fgo.jacobian(fx, fa);
  checker.compareMatrices(PDELab::Backend::native(a), PDELab::Backend::native(fa), perm, "jacobian");
  checker.check(PDELab::Backend::native(a).nonzeroes() == PDELab::Backend::native(fa).nonzeroes(), "number of pattern entries");

  // jacobian application (analytic for ConvectionDiffusionDG)

  PVector y(gfs, 0.0);
  go.jacobian_apply(z, y);
  Vector fy(basis, 0.0);
  fgo.jacobian_apply(fz, fy);
  checker.compareVectors(PDELab::Backend::native(y), PDELab::Backend::native(fy), perm, "jacobian_apply");

  auto fAz = PDELab::Backend::native(fy);
  PDELab::Backend::native(fa).mv(PDELab::Backend::native(fz), fAz);
  fAz -= PDELab::Backend::native(fy);
  std::cout << "|A z - jacobian_apply(z)|_inf = " << fAz.infinity_norm() << std::endl;
  checker.check(fAz.infinity_norm() < 1e-10 * PDELab::Backend::native(fy).infinity_norm(), "matrix-vector product vs. jacobian_apply");

  return checker.result();
}



int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  using Grid = Dune::YaspGrid<2>;
  Grid grid({1.0,1.0}, Dune::filledArray<2,int>(8));
  auto gv = grid.leafGridView();

  int result = 0;
  result |= testPoisson<1>(gv);
  result |= testPoisson<2>(gv);
  result |= testDG<1>(gv);
  result |= testDG<2>(gv);
  return result;
}
