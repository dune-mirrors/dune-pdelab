#include "../gridexamples.hh"

#include <dune/pdelab/basis/merging_strategy.hh>
#include <dune/pdelab/basis/protobasis.hh>

#include <dune/pdelab/finiteelementmap/qkfem.hh>

#include <dune/pdelab/basis/defaultprebasis.hh>

#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

template<std::size_t k, typename R=double, class MergingStragetgy>
auto lagrange(const MergingStragetgy& strategy)
{
  return [=]<class GridView>(const GridView& grid_view) {
    using QkFEM = Dune::PDELab::QkLocalFiniteElementMap<GridView, R, R, k>;
    return Dune::PDELab::Impl::ProtoBasisLeaf<GridView, QkFEM, MergingStragetgy>(grid_view, QkFEM(grid_view), strategy);
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
  using ProtoBasis = decltype(q1_proto_basis);

  using PreBasis = Dune::PDELab::DefaultPreBasis<ProtoBasis>;
  PreBasis pre_basis(q1_proto_basis);

  static_assert(Dune::models<Dune::Functions::Concept::PreBasis<GridView>, PreBasis>(), "....");

  using Basis = Dune::Functions::DefaultGlobalBasis<PreBasis>;
  Basis basis(q1_proto_basis);

  auto test = checkBasis(basis, EnableContinuityCheck());

  // using ProtoBasisArray = Dune::PDELab::Impl::ProtoBasisArray<ProtoBasis, ProtoBasis, 2>;
  // ProtoBasisArray protoBasisArray(q1_proto_basis, 2, flat_strategy);
}
