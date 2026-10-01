// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/* GridOperator on a dune-functions basis with blocked indices:
 * LinearElasticity on power<2>(lagrange<1>(), blockedInterleaved()), i.e.
 * PDELab wrappers around BlockVector<FieldVector<double,2>> and
 * BCRSMatrix<FieldMatrix<double,2,2>>. Residual and jacobian are compared
 * with the GridOperator on a PDELab power grid function space, and the
 * Dirichlet constraints (unit rows within FieldMatrix blocks) are checked.
 */

#include "config.h"

#include <cmath>
#include <iostream>
#include <vector>

#include <dune/common/filledarray.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/localoperator/linearelasticity.hh>
#include <dune/pdelab/localoperator/linearelasticityparameter.hh>

#include "gridoperatorfunctionstesthelpers.hh"

template<class GV>
class ElasticityProblem
  : public Dune::PDELab::LinearElasticityParameterInterface<
      Dune::PDELab::LinearElasticityParameterTraits<GV,double>, ElasticityProblem<GV>>
{
public:
  using Traits = Dune::PDELab::LinearElasticityParameterTraits<GV,double>;

  void f(const typename Traits::ElementType& e, const typename Traits::DomainType& x,
         typename Traits::RangeType& y) const
  {
    auto xg = e.geometry().global(x);
    y = {1.0, xg[0]*xg[1]};
  }

  template<typename I>
  bool isDirichlet(const I&, const typename Traits::IntersectionDomainType&) const
  {
    return false;
  }

  void g(const typename Traits::ElementType&, const typename Traits::DomainType&,
         typename Traits::RangeType& y) const
  {
    y = 0.0;
  }

  typename Traits::RangeFieldType lambda(const typename Traits::ElementType&, const typename Traits::DomainType&) const
  {
    return 1.2;
  }

  typename Traits::RangeFieldType mu(const typename Traits::ElementType&, const typename Traits::DomainType&) const
  {
    return 0.7;
  }
};

// entry number k of a vector when iterating over all scalar entries
template<class NativeVector>
double& flatEntry(NativeVector& v, std::size_t k)
{
  constexpr std::size_t n = NativeVector::block_type::dimension;
  return v[k/n][k%n];
}

template<class GV>
int testBlocked(const GV& gv)
{
  using namespace Dune;
  using namespace Dune::PDELab::FunctionsBasisTest;
  std::cout << "=== LinearElasticity, power<2>(Q1), blocked indices" << std::endl;

  using Param = ElasticityProblem<GV>;
  Param param;
  using LOP = PDELab::LinearElasticity<Param>;
  LOP lop(param);
  using MBE = PDELab::ISTL::BCRSMatrixBackend<>;

  // PDELab reference: power grid function space, flat lexicographic ordering

  using VBE = PDELab::ISTL::VectorBackend<>;
  using FEM = PDELab::QkLocalFiniteElementMap<GV,double,double,1>;
  FEM fem(gv);
  using LeafGFS = PDELab::GridFunctionSpace<GV,FEM,PDELab::NoConstraints,VBE>;
  using GFS = PDELab::PowerGridFunctionSpace<LeafGFS,2,VBE,PDELab::LexicographicOrderingTag>;
  LeafGFS leafGfs(gv,fem);
  GFS gfs(leafGfs);
  using GO = PDELab::GridOperator<GFS,GFS,LOP,MBE,double,double,double,
                                  PDELab::EmptyTransformation,PDELab::EmptyTransformation>;
  GO go(gfs,gfs,lop,MBE(18));

  // dune-functions basis with blocked indices

  using namespace Functions::BasisFactory;
  auto basis = makeBasis(gv, power<2>(lagrange<1>(), blockedInterleaved()));
  using Basis = decltype(basis);
  using FGO = PDELab::GridOperator<Basis,Basis,LOP,MBE,double,double,double>;
  FGO fgo(basis, basis, lop, MBE(9));

  using Vector = typename FGO::Traits::Domain;
  using Matrix = typename FGO::Traits::Jacobian;
  static_assert(std::is_same_v<PDELab::Backend::Native<Vector>, BlockVector<FieldVector<double,2>>>);
  static_assert(std::is_same_v<PDELab::Backend::Native<Matrix>, BCRSMatrix<FieldMatrix<double,2,2>>>);

  // match DOFs: keys (x + 10 d, y) for component d

  using PVector = typename GO::Traits::Domain;
  using Range2 = FieldVector<double,2>;
  PVector key1(gfs), key2(gfs);
  PDELab::interpolate(PDELab::makeGridFunctionFromCallable(gv, [](const auto& x) { return Range2{x[0], x[0]+10}; }), gfs, key1);
  PDELab::interpolate(PDELab::makeGridFunctionFromCallable(gv, [](const auto& x) { return Range2{x[1], x[1]}; }), gfs, key2);
  Vector fkey1(basis, 0.0), fkey2(basis, 0.0);
  for (int d = 0; d < 2; ++d)
  {
    Functions::interpolate(Functions::subspaceBasis(basis, d), PDELab::Backend::native(fkey1), [d](const auto& x) { return x[0] + 10*d; });
    Functions::interpolate(Functions::subspaceBasis(basis, d), PDELab::Backend::native(fkey2), [](const auto& x) { return x[1]; });
  }
  auto perm = matchDofs(flatValues(PDELab::Backend::native(key1)), flatValues(PDELab::Backend::native(key2)),
                        flatValues(PDELab::Backend::native(fkey1)), flatValues(PDELab::Backend::native(fkey2)));

  Checker checker(1e-10);

  auto transferFlat = [&](const PVector& p, Vector& f) {
    auto values = flatValues(PDELab::Backend::native(p));
    for (std::size_t i = 0; i < values.size(); ++i)
      flatEntry(PDELab::Backend::native(f), perm[i]) = values[i];
  };
  auto compareFlat = [&](const PVector& p, const Vector& f, const std::string& what) {
    auto pv = flatValues(PDELab::Backend::native(p));
    auto fv = flatValues(PDELab::Backend::native(f));
    double maxDiff = 0.0, maxValue = 0.0;
    for (std::size_t i = 0; i < pv.size(); ++i)
    {
      maxDiff = std::max(maxDiff, std::abs(pv[i] - fv[perm[i]]));
      maxValue = std::max(maxValue, std::abs(pv[i]));
    }
    std::cout << what << ": max |pdelab| = " << maxValue << ", max difference = " << maxDiff << std::endl;
    checker.check(maxDiff <= 1e-10 * std::max(1.0, maxValue), what);
  };

  PVector x(gfs);
  randomize(PDELab::Backend::native(x), 42);
  Vector fx(basis, 0.0);
  transferFlat(x, fx);

  // residual
  PVector r(gfs, 0.0);
  go.residual(x, r);
  Vector fr(basis, 0.0);
  fgo.residual(fx, fr);
  compareFlat(r, fr, "residual");

  // jacobian, compared through matrix-vector products and the number of scalar nonzeros
  using PMatrix = typename GO::Traits::Jacobian;
  PMatrix a(go, 0.0);
  go.jacobian(x, a);
  Matrix fa(fgo, 0.0);
  fgo.jacobian(fx, fa);
  std::cout << "jacobian: scalar nonzeros pdelab/functions = " << PDELab::Backend::native(a).nonzeroes()
            << "/" << 4*PDELab::Backend::native(fa).nonzeroes() << std::endl;
  checker.check(PDELab::Backend::native(a).nonzeroes() == 4*PDELab::Backend::native(fa).nonzeroes(), "number of nonzeros");
  for (unsigned int seed : {1u, 2u, 3u})
  {
    PVector z(gfs), y(gfs, 0.0);
    randomize(PDELab::Backend::native(z), seed);
    Vector fz(basis, 0.0), fy(basis, 0.0);
    transferFlat(z, fz);
    PDELab::Backend::native(a).mv(PDELab::Backend::native(z), PDELab::Backend::native(y));
    PDELab::Backend::native(fa).mv(PDELab::Backend::native(fz), PDELab::Backend::native(fy));
    compareFlat(y, fy, "jacobian times random vector");
  }

  // Dirichlet constraints for both components on the boundary x = 0
  using CC = PDELab::ConstraintsContainer<Basis>;
  CC cc;
  Functions::forEachBoundaryDOF(basis, [&](auto&& localIndex, const auto& localView, const auto& intersection) {
    if (intersection.geometry().center()[0] < 1e-8)
      cc[localView.index(localIndex)];
  });
  using CGO = PDELab::GridOperator<Basis,Basis,LOP,MBE,double,double,double,CC,CC>;
  CGO cgo(basis, cc, basis, cc, lop, MBE(9));
  Matrix ca(cgo, 0.0);
  cgo.jacobian(fx, ca);
  Vector cr(basis, 0.0);
  cgo.residual(fx, cr);
  bool unitRows = true;
  for (const auto& [index, row] : cc)
  {
    unitRows = unitRows and cr[index] == 0.0;
    const auto& matrixRow = PDELab::Backend::native(ca)[index[0]];
    for (auto block = matrixRow.begin(); block != matrixRow.end(); ++block)
      for (int j = 0; j < 2; ++j)
      {
        const double expected = (block.index() == index[0] and j == int(index[1])) ? 1.0 : 0.0;
        unitRows = unitRows and (*block)[index[1]][j] == expected;
      }
  }
  std::cout << "constraints: " << cc.size() << " constrained DOFs" << std::endl;
  checker.check(cc.size() > 0 and unitRows, "unit rows and zero residual for constrained DOFs");

  return checker.result();
}


int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  using Grid = Dune::YaspGrid<2>;
  Grid grid({1.0,1.0}, Dune::filledArray<2,int>(6));

  return testBlocked(grid.leafGridView());
}
