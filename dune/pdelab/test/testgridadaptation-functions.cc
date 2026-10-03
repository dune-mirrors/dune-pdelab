// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/* Solution transfer during grid adaptation (adapt_grid, adaptGrid) for
 * dune-functions bases:
 * - polynomial data in the discrete space is transferred exactly when the
 *   grid is refined and when it is coarsened back (interpolation on new
 *   elements, local L2 projection onto vanishing elements), for Lagrange,
 *   DG and Taylor-Hood bases;
 * - the transfer gives the same result for a PDELab grid function space and a
 *   dune-functions basis spanning the same space; both are adapted in one
 *   adaptGrid call.
 */

#include "config.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

#include <dune/pdelab.hh>

#include "gridoperatorfunctionstesthelpers.hh"

using namespace Dune;
using namespace Dune::PDELab::FunctionsBasisTest;

using TestGrid = UGGrid<2>;

std::unique_ptr<TestGrid> makeGrid()
{
  auto grid = StructuredGridFactory<TestGrid>::createSimplexGrid({0.0,0.0}, {1.0,1.0}, {2u,2u});
  grid->globalRefine(2);
  return grid;
}

// refine the elements left of x = 0.5
void markForRefinement(TestGrid& grid)
{
  for (const auto& e : elements(grid.leafGridView()))
    if (e.geometry().center()[0] < 0.5)
      grid.mark(1, e);
}

// coarsen everything that was refined
void markForCoarsening(TestGrid& grid)
{
  for (const auto& e : elements(grid.leafGridView()))
    if (e.level() > 2)
      grid.mark(-1, e);
}

template<class X>
double maxDifference(const X& a, const X& b)
{
  auto av = flatValues(PDELab::Backend::native(a));
  auto bv = flatValues(PDELab::Backend::native(b));
  double d = 0.0;
  for (std::size_t i = 0; i < av.size(); ++i)
    d = std::max(d, std::abs(av[i] - bv[i]));
  return d;
}

// interpolate polynomial data into x (all components); scalar and Taylor-Hood bases
template<class Basis, class X>
void interpolateData(const Basis& basis, X& x, int order)
{
  using Domain = FieldVector<double,2>;
  auto quadratic = [](const Domain& p) { return 1.0 + p[0] - 2.0*p[1] + 0.5*p[0]*p[1]; };
  auto linear = [](const Domain& p) { return 0.5 - p[0] + 3.0*p[1]; };
  using Tree = typename Basis::LocalView::Tree;
  if constexpr (TypeTree::Concept::LeafTreeNode<Tree>)
  {
    if (order >= 2)
      Functions::interpolate(basis, PDELab::Backend::native(x), quadratic);
    else
      Functions::interpolate(basis, PDELab::Backend::native(x), linear);
  }
  else
  {
    using namespace Indices;
    auto velocity = [&](const Domain& p) { return FieldVector<double,2>{quadratic(p), 2.0*quadratic(p) - 1.0}; };
    Functions::interpolate(Functions::subspaceBasis(basis, _0), PDELab::Backend::native(x), velocity);
    Functions::interpolate(Functions::subspaceBasis(basis, _1), PDELab::Backend::native(x), linear);
  }
}

// refine and coarsen back; polynomial data must be preserved exactly
template<class Basis>
bool testExactness(const std::string& name, Basis& basis, TestGrid& grid, int order)
{
  using X = PDELab::Backend::Vector<Basis,double>;
  bool ok = true;
  X x(basis, 0.0);
  interpolateData(basis, x, order);
  const int intorder = 2*order;

  auto check = [&](const std::string& step) {
    X reference(basis, 0.0);
    interpolateData(basis, reference, order);
    const double d = maxDifference(x, reference);
    const bool passed = d < 1e-11;
    std::cout << (passed ? "ok     " : "FAILED ") << name << " " << step
              << ": " << PDELab::globalSize(basis) << " DOFs, max difference " << d << std::endl;
    ok = ok and passed;
  };

  markForRefinement(grid);
  PDELab::adapt_grid(grid, basis, x, intorder);
  check("after refinement");

  markForCoarsening(grid);
  PDELab::adapt_grid(grid, basis, x, intorder);
  check("after coarsening");
  return ok;
}

// same random data on a GFS and a basis, both adapted in one adaptGrid call
template<int k>
bool testGFSvsBasis(TestGrid& grid)
{
  using GV = TestGrid::LeafGridView;
  auto gv = grid.leafGridView();
  using FEM = PDELab::PkLocalFiniteElementMap<GV,double,double,k>;
  FEM fem(gv);
  using GFS = PDELab::GridFunctionSpace<GV,FEM,PDELab::NoConstraints,PDELab::ISTL::VectorBackend<>>;
  GFS gfs(gv, fem);
  auto basis = Functions::BasisFactory::makeBasis(gv, Functions::BasisFactory::lagrange<k>());
  using Basis = decltype(basis);
  using XG = PDELab::Backend::Vector<GFS,double>;
  using XB = PDELab::Backend::Vector<Basis,double>;

  // DOF permutation via the coordinates of the Lagrange nodes
  auto permutation = [&]() {
    XG key1(gfs), key2(gfs);
    PDELab::interpolate(PDELab::makeGridFunctionFromCallable(gfs.gridView(), [](const auto& p) { return p[0]; }), gfs, key1);
    PDELab::interpolate(PDELab::makeGridFunctionFromCallable(gfs.gridView(), [](const auto& p) { return p[1]; }), gfs, key2);
    XB fkey1(basis, 0.0), fkey2(basis, 0.0);
    Functions::interpolate(basis, PDELab::Backend::native(fkey1), [](const auto& p) { return p[0]; });
    Functions::interpolate(basis, PDELab::Backend::native(fkey2), [](const auto& p) { return p[1]; });
    return matchDofs(flatValues(PDELab::Backend::native(key1)), flatValues(PDELab::Backend::native(key2)),
                     scalarValues(PDELab::Backend::native(fkey1)), scalarValues(PDELab::Backend::native(fkey2)));
  };

  XG xg(gfs, 0.0);
  randomize(PDELab::Backend::native(xg), 17);
  XB xb(basis, 0.0);
  transfer(PDELab::Backend::native(xg), permutation(), PDELab::Backend::native(xb));

  const int intorder = 2*k;
  bool ok = true;
  Checker checker(1e-12);
  for (int step = 0; step < 2; ++step)
  {
    if (step == 0) markForRefinement(grid); else markForCoarsening(grid);
    auto packG = PDELab::transferSolutions(gfs, intorder, xg);
    auto packB = PDELab::transferSolutions(basis, intorder, xb);
    PDELab::adaptGrid(grid, packG, packB);
    checker.compareVectors(PDELab::Backend::native(xg), PDELab::Backend::native(xb), permutation(),
                           "GFS vs. basis P" + std::to_string(k) + (step == 0 ? " after refinement" : " after coarsening"));
  }
  return checker.result() == 0;
}

int main(int argc, char** argv)
{
  try {
    MPIHelper::instance(argc, argv);
    using namespace Functions::BasisFactory;
    bool ok = true;

    {
      auto grid = makeGrid();
      auto basis = makeBasis(grid->leafGridView(), lagrange<1>());
      ok = testExactness("lagrange<1>", basis, *grid, 1) and ok;
    }
    {
      auto grid = makeGrid();
      auto basis = makeBasis(grid->leafGridView(), lagrange<2>());
      ok = testExactness("lagrange<2>", basis, *grid, 2) and ok;
    }
    {
      auto grid = makeGrid();
      auto basis = makeBasis(grid->leafGridView(), lagrangeDG<1>());
      ok = testExactness("lagrangeDG<1>", basis, *grid, 1) and ok;
    }
    {
      auto grid = makeGrid();
      auto basis = makeBasis(grid->leafGridView(),
                             composite(power<2>(lagrange<2>(), flatInterleaved()), lagrange<1>(), flatLexicographic()));
      ok = testExactness("Taylor-Hood", basis, *grid, 2) and ok;
    }
    {
      auto grid = makeGrid();
      ok = testGFSvsBasis<1>(*grid) and ok;
    }
    {
      auto grid = makeGrid();
      ok = testGFSvsBasis<2>(*grid) and ok;
    }

    return ok ? 0 : 1;
  }
  catch (Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
}
