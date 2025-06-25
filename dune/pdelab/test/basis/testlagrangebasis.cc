#include "../gridexamples.hh"

#include <dune/pdelab/basis/utility.hh>
#include <dune/pdelab/basis/lagrange.hh>

#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

int main(int argc, char **argv)
{
  Dune::MPIHelper::instance(argc, argv);

  YaspUnitSquare grid;
  grid.globalRefine(1);
  std::cout << "Grid has " << grid.size(0) << " cells, "
            << grid.size(1) << " edges, and "
            << grid.size(2) << " vertices." << std::endl;
  std::cout << "\n--- Lagrange Basis Test Suite ---\n";
  std::cout << "This test covers various Lagrange basis configurations using PDELab and DUNE Functions.\n";

  using GridView = typename YaspUnitSquare::LeafGridView;
  GridView grid_view{grid.leafGridView()};

  using namespace Dune::PDELab::BasisFactory;
  using namespace Dune::Functions::BasisFactory;

  Dune::TestSuite test("Lagrange Test Suite");

  using namespace Dune::PDELab::BasisFactory::Experimental;
  auto Q2 = FiniteElement(CG, quadrilateral, _2);
  auto Q1 = FiniteElement(CG, quadrilateral, _1);
  auto Q1b = FiniteElement(CG, quadrilateral, _1, blockedByEntity());
  auto Q2b = FiniteElement(CG, quadrilateral, _2, blockedByEntity());

  {
    std::cout << "Testing: flatByEntity Q2 Lagrange (QkFEM)" << std::endl;
    auto basis = Q2 || grid_view;
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("flatByEntity_Q2_index_tree.txt", basis);
  }

  {
    std::cout << "Testing: blockedByEntity Q2 Lagrange (QkFEM)" << std::endl;
    auto basis = Q2b || grid_view;
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("blockedByEntity_Q2_index_tree.txt", basis);
  }

  // power nodes
  {
    std::cout << "Testing: power<5> flatByEntity Q2 Lagrange (QkFEM)" << std::endl;
    auto basis = Q2^_5 || grid_view;
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("power5_flatByEntity_Q2_index_tree.txt", basis);
  }

  {
    std::cout << "Testing: power<2> blockedByEntity Q2 Lagrange (QkFEM)" << std::endl;
    auto basis = power<5>(Q2b, blockedByEntity()) || grid_view;
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("power2_blockedByEntity_Q2_index_tree.txt", basis);
  }

  {
    std::cout << "Testing: blockedByEntity-blockedByEntity taylor-hood" << std::endl;
    auto basis = makeBasis(grid_view,
                           composite(
                             power<2>(
                               makeQkProtoBasis<2>(blockedByEntity()),
                               blockedByEntity()),
                             makeQkProtoBasis<1>(blockedByEntity()),
                             blockedByEntity()
                           ));
    test.check(basis.size({}) == basis.gridView().size(0) + basis.gridView().size(1)+basis.gridView().size(2)) << "Expected size " << basis.gridView().size(0) + basis.gridView().size(1)+basis.gridView().size(2) << " but got " << basis.size({});
    for (std::size_t i=0; i != basis.gridView().size(0) + basis.gridView().size(1)+basis.gridView().size(2); ++i) {
      test.check(basis.size({i}) == 2) << "Expected size " << 2 << " but got " << basis.size({i});
      test.check(basis.size({i,0}) == 2) << "Expected size " << 2 << " but got " << basis.size({i,0});
      test.check(basis.size({i,1}) == 1) << "Expected size " << 1 << " but got " << basis.size({i,1});
    }
    Dune::PDELab::writeIndexTree("TaylorHood_blockedByEntity_blockedByEntity_index_tree.txt", basis);
  }

  {
    std::cout << "Testing: blockedByEntity-flatByEntity taylor-hood" << std::endl;
    auto basis = composite(Q2^_2, Q1, blockedByEntity()) || grid_view;
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    test.check(basis.size({}) == basis.gridView().size(0) + basis.gridView().size(1)+basis.gridView().size(2)) << "Expected size " << basis.gridView().size(0) + basis.gridView().size(1)+basis.gridView().size(2) << " but got " << basis.size({});
    for (std::size_t i=0; i != basis.gridView().size(0) + basis.gridView().size(1)+basis.gridView().size(2); ++i)
      test.check(basis.size({i}) == ((i < basis.gridView().size(2)) ? 3 : 2)) << "Expected size " << ((i < basis.gridView().size(2)) ? 3 : 2) << " but got " << basis.size({i});
    Dune::PDELab::writeIndexTree("TaylorHood_blockedByEntity_flatByEntity_index_tree.txt", basis);
  }

  {
    std::cout << "Testing: flatByEntity taylor-hood" << std::endl;
    auto basis = ((Q2^_2) * Q1) || grid_view;
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    test.check(basis.size({}) == 2*grid.size(0) + 2*grid.size(1) + 3*grid.size(2))
      << "Expected size " << (2*grid.size(0) + 2*grid.size(1) + 3*grid.size(2)) << " but got " << basis.size({});
    Dune::PDELab::writeIndexTree("TaylorHood_flatByEntity_index_tree.txt", basis);
  }

  {
    std::cout << "Testing: blockedLexicographic-flatByEntity taylor-hood" << std::endl;
    using namespace Dune::Functions::BasisFactory;
    auto basis = makeBasis(grid_view,
                           composite((Q2^_2), Q1,
                           blockedLexicographic()
                           ));
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    test.check(basis.size({}) == 2) << "Expected size " << 2 << " but got " << basis.size({});
    test.check(basis.size({0}) == 2*grid.size(0) + 2*grid.size(1) + 2*grid.size(2)) << "Expected size " << (2*grid.size(0) + 2*grid.size(1) + 2*grid.size(2)) << " but got " << basis.size({0});
    test.check(basis.size({1}) == grid.size(2)) << "Expected size " << (grid.size(2)) << " but got " << basis.size({1});
    Dune::PDELab::writeIndexTree("TaylorHood_blockedLexicographic_flatByEntity_index_tree.txt", basis);
  }

  {
    std::cout << "Testing: power<5> blockedByEntity Q2 Lagrange (QkFEM), flatByEntity" << std::endl;
    auto basis = Q2b^5 || grid_view;
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("power5_blockedByEntity_flatByEntity_Q2_index_tree.txt", basis);
  }

  {
    std::cout << "Testing: power<5> flatByEntity Q2 Lagrange (QkFEM), blockedByEntity" << std::endl;
    auto basis = power(Q2b, 5, blockedByEntity()) || grid_view;
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("power5_flatByEntity_blockedByEntity_Q2_index_tree.txt", basis);
  }

  // check compatibility with dune-functions
  {
    std::cout << "Testing: power<5> flatByEntity Q2 Lagrange (QkFEM), flatLexicographic" << std::endl;
    using namespace Dune::Functions::BasisFactory;
    auto basis = makeBasis(grid_view, power<5>(Q2, flatLexicographic()));
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("power5_flatByEntity_flatLexicographic_Q2_index_tree.txt", basis);
  }

  {
    std::cout << "Testing: power<2> blockedByEntity Q2 Lagrange (QkFEM), blockedLexicographic" << std::endl;
    using namespace Dune::Functions::BasisFactory;
    auto basis = makeBasis(grid_view, power<2>(Q2b, blockedLexicographic()));
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("power2_blockedByEntity_blockedLexicographic_Q2_index_tree.txt", basis);
  }

  {
    std::cout << "Testing: power<5> blockedByEntity Q2 Lagrange (QkFEM), flatLexicographic" << std::endl;
    auto basis = makeBasis(grid_view, power<5>(Q2b, flatLexicographic()));
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("power5_blockedByEntity_flatLexicographic_Q2_index_tree.txt", basis);
  }

  {
    std::cout << "Testing: power<5> flatByEntity Q2 Lagrange (QkFEM), blockedLexicographic" << std::endl;
    auto basis = makeBasis(grid_view, power<5>(Q2, blockedLexicographic()));
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("power5_flatByEntity_blockedLexicographic_Q2_index_tree.txt", basis);
  }

  // composite nodes
  {
    std::cout << "Testing: composite power<5> blockedByEntity Q1 Lagrange (QkFEM), power<2> blockedByEntity Q1 Lagrange (QkFEM), blockedByEntity" << std::endl;
    auto basis =  composite(
                      power<5>(Q1b, blockedByEntity()),
                      power<2>(Q1b, blockedByEntity()),
                      blockedByEntity()
                  ) || grid_view;
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("composite_power5_power2_blockedByEntity_Q1_index_tree.txt", basis);
  }

  {
    std::cout << "Testing: composite power<5> blockedByEntity Q2 Lagrange (QkFEM), power<2> blockedByEntity Q2 Lagrange (QkFEM), flatByEntity" << std::endl;
    auto basis =  (
                    power<5>(Q2b, blockedByEntity()) * power<2>(Q2b, blockedByEntity())
                  ) || grid_view;
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("composite_power5_power2_blockedByEntity_flatByEntity_Q2_index_tree.txt", basis);
  }

  {
    std::cout << "Testing: composite power<5> flatByEntity Q2 Lagrange (QkFEM), power<2> flatByEntity Q2 Lagrange (QkFEM), flatByEntity" << std::endl;
    auto basis = ((Q2^_5) * (Q2^_2)) || grid_view;
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("composite_power5_power2_flatByEntity_Q2_index_tree.txt", basis);
  }

  {
    std::cout << "Testing: composite power<5> flatByEntity Q2 Lagrange (QkFEM), power<2> blockedByEntity Q2 Lagrange (QkFEM), flatByEntity" << std::endl;
    auto basis = makeBasis(grid_view,
                           composite(
                               power<5>(
                                   Q2,
                                   blockedByEntity()),
                               power<2>(
                                   Q2b,
                                   blockedByEntity()),
                               flatByEntity()));
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("composite_power5_flatByEntity_power2_blockedByEntity_Q2_index_tree.txt", basis);
  }

  // check compatibility with dune-functions
  {
    std::cout << "Testing: composite power<5> flatByEntity Q2 Lagrange (QkFEM), power<2> blockedByEntity Q2 Lagrange (QkFEM), flatLexicographic" << std::endl;
    auto basis = makeBasis(grid_view,
                           composite(
                               power<5>(
                                   Q2,
                                   blockedLexicographic()),
                               power<2>(
                                   Q2b,
                                   blockedByEntity()),
                               flatLexicographic()));
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("composite_power5_flatByEntity_power2_blockedByEntity_flatLexicographic_Q2_index_tree.txt", basis);
  }

  {
    std::cout << "Testing: composite power<5> DUNE Functions lagrange<2>(), blockedLexicographic, power<2> blockedByEntity Q2 Lagrange (QkFEM), blockedByEntity, blockedLexicographic" << std::endl;
    auto basis = makeBasis(grid_view,
                           composite(
                               power<5>(
                                   lagrange<2>(),
                                   blockedLexicographic()),
                               power<2>(
                                   Q2b,
                                   blockedByEntity()),
                               blockedLexicographic()));
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("composite_power5_lagrange2_blockedLexicographic_power2_blockedByEntity_blockedLexicographic_Q2_index_tree.txt", basis);
  }

  {
    std::cout << "Testing: composite power<5> DUNE Functions lagrange<2>(), flatLexicographic, power<2> flatByEntity Q2 Lagrange (QkFEM), flatByEntity, flatLexicographic" << std::endl;
    auto basis = makeBasis(grid_view,
                           composite(
                               power<5>(
                                   lagrange<2>(),
                                   flatLexicographic()),
                               power<2>(
                                   Q2,
                                   flatByEntity()),
                               flatLexicographic()));
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    Dune::PDELab::writeIndexTree("composite_power5_lagrange2_flatLexicographic_power2_flatByEntity_flatLexicographic_Q2_index_tree.txt", basis);
  }

  return test.exit();
}
