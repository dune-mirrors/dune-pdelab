#include "../gridexamples.hh"

#include <dune/pdelab/basis/merging_strategy.hh>
#include <dune/pdelab/basis/protobasis.hh>

#include <dune/pdelab/finiteelementmap/qkfem.hh>

#include <dune/pdelab/basis/prebasis.hh>

#include <dune/functions/functionspacebases/concepts.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

template<std::size_t k, typename R=double, class MergingStragetgy>
auto lagrange(const MergingStragetgy& strategy)
{
  return [=]<class GridView>(const GridView& grid_view) {
    using QkFEM = Dune::PDELab::QkLocalFiniteElementMap<GridView, R, R, k>;
    auto qkfem = std::make_shared<QkFEM>(grid_view);
    return protoBasis(strategy, qkfem);
  };
}

int main(int argc, char **argv) {
  Dune::MPIHelper::instance(argc, argv);

  YaspUnitSquare grid;
  grid.globalRefine(3);
  using GridView = typename YaspUnitSquare::LeafGridView;
  GridView grid_view{grid.leafGridView()};

  auto flat_strategy = Dune::PDELab::BasisFactory::flatByEntity();
  auto q1_proto_basis = lagrange<1>(flat_strategy)(grid_view);

  using PreBasis = Dune::PDELab::PreBasis<decltype(q1_proto_basis), GridView>;
  PreBasis pre_basis(q1_proto_basis, grid_view);

  static_assert(Dune::models<Dune::Functions::Concept::PreBasis<GridView>, PreBasis>(), "....");
}
