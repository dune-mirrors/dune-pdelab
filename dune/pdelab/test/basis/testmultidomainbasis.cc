#include "../gridexamples.hh"

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/functions/backends/istlvectorfactory.hh>
#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/pdelab/basis/utility.hh>
#include <dune/pdelab/basis/lagrange.hh>

#include <numbers>

int main(int argc, char **argv)
{
  Dune::MPIHelper::instance(argc, argv);

  YaspUnitSquare grid;
  grid.globalRefine(2);
  std::cout << "Grid has " << grid.size(0) << " cells, "
            << grid.size(1) << " edges, and "
            << grid.size(2) << " vertices." << std::endl;
  std::cout << "\n--- Lagrange Basis Test Suite ---\n";
  std::cout << "This test covers various Lagrange basis configurations using PDELab and DUNE Functions.\n";

  using GridView = typename YaspUnitSquare::LeafGridView;
  GridView gridView{grid.leafGridView()};
  constexpr auto dim = Dune::index_constant<GridView::dimension>{};

  Dune::TestSuite test("Multi-Domain Basis Test Suite");

  struct SimpleSubDomain
  {
    bool contains(const GridView::template Codim<0>::Entity& entity) const
    {
      return indicator(entity.geometry().center());
    }

    std::function<bool(Dune::FieldVector<double, GridView::dimension>)> indicator;
  };

  // Indicator functions for four square subdomains covering the four quadrants
  // and one overlapping square subdomain in the center
  auto subDomainA = SimpleSubDomain{[](auto x) {
    return (x[0] <= 0.5) and (x[1] <= 0.5);
  }};

  auto subDomainB = SimpleSubDomain{[](auto x) {
    return (x[0] > 0.5) and (x[1] <= 0.5);
  }};

  auto subDomainC = SimpleSubDomain{[](auto x) {
    return (x[0] <= 0.5) and (x[1] > 0.5);
  }};

  using namespace Dune::PDELab::BasisFactory;
  using namespace Dune::PDELab::BasisFactory::Experimental;

  auto Q2 = FiniteElement(CG, quadrilateral, _2);
  auto Q1 = FiniteElement(CG, quadrilateral, _1);

  // Create a basis
  auto basis = composite(
                  Q2 | subDomainA,
                  (Q2^dim | subDomainB) * (Q1 | subDomainC),
                  blockedByEntity()
                ) || gridView;

  // Run basis check
  test.subTest(checkBasis(basis));
  // notice how the middle corner gathers all 4 dofs in one block
  Dune::PDELab::writeIndexTree("multidomain_index_tree.txt", basis);

  auto sinsin = [](auto x) {
    return std::sin(x[0]*2*std::numbers::pi)*std::sin(x[1]*2*std::numbers::pi);
  };

  // fails because dune-functions composite seems to not be able to handle non-uniform vectors yet
  auto c = Dune::Functions::makeISTLVector<double>(basis.preBasis().containerDescriptor());
  auto c_backend = Dune::Functions::istlVectorBackend(c);

  using Dune::Functions::subspaceBasis;
  using Dune::Functions::interpolate;

  interpolate(subspaceBasis(basis, _0), c_backend, sinsin);
  interpolate(subspaceBasis(basis, _1, _0, 0), c_backend, sinsin);
  interpolate(subspaceBasis(basis, _1, _0, 1), c_backend, sinsin);
  interpolate(subspaceBasis(basis, _1, _1), c_backend, sinsin);

  // // Finally we write all grid functions to vtk
  auto c0_gf = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(subspaceBasis(basis, _0), c_backend);
  auto c100_gf = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(subspaceBasis(basis, _1, _0, 0), c_backend);
  auto c101_gf = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(subspaceBasis(basis, _1, _0, 1), c_backend);
  auto c11_gf = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(subspaceBasis(basis, _1,_1), c_backend);

  Dune::SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, Dune::refinementLevels(3));
  vtkWriter.addVertexData(c0_gf, Dune::VTK::FieldInfo("c0", Dune::VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.addVertexData(c100_gf, Dune::VTK::FieldInfo("c100", Dune::VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.addVertexData(c101_gf, Dune::VTK::FieldInfo("c101", Dune::VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.addVertexData(c11_gf, Dune::VTK::FieldInfo("c11", Dune::VTK::FieldInfo::Type::scalar, 1));

  vtkWriter.write("testmultidomainbasis");

  return test.exit();
}
