#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>

#include <array>
#include <memory>

int
main(int argc, char* argv[])
{
  try
  {
    Dune::MPIHelper::instance(argc, argv);

    // Create a grid
    constexpr int dim = 2;
    using Grid = Dune::UGGrid<dim>;
    Dune::FieldVector<double, dim> lower(0), upper(1);
    std::array<unsigned int, dim> elements = { 10, 10 };
    std::shared_ptr<Grid> grid(
      Dune::StructuredGridFactory<Grid>::createSimplexGrid(
        lower, upper, elements));

    // Entity Set
    constexpr int degree = 2;
    using GridView = typename Grid::LeafGridView;
    using EntitySet = Dune::PDELab::OverlappingEntitySet<GridView>;
    auto es = std::make_shared<EntitySet>(grid->leafGridView());

    // Finite Element Map
    using FEM =
      Dune::PDELab::PkLocalFiniteElementMap<EntitySet, double, double, degree>;
    auto fem = std::make_shared<FEM>(*es);

    // Vector Grid Function Space
    using VectorBackend =
      Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
    using LeafVectorBackend = Dune::PDELab::ISTL::VectorBackend<>;
    using GridFunctionSpace =
      Dune::PDELab::VectorGridFunctionSpace<EntitySet,
                                            FEM,
                                            dim,
                                            VectorBackend,
                                            LeafVectorBackend>;
    auto gfs = std::make_shared<GridFunctionSpace>(*es, fem);

    // Create a DOF vector
    using DOFVector = Dune::PDELab::Backend::Vector<GridFunctionSpace, double>;
    auto dof_vector = std::make_shared<DOFVector>(gfs);

    // Interpolate a function
    auto position = [](const auto& e, const auto& x) {
      return e.geometry().global(x);
    };
    auto gf =
      Dune::PDELab::makeGridFunctionFromCallable(es->gridView(), position);
    Dune::PDELab::interpolate(gf, *gfs, *dof_vector);

    // Create a grid function space gradient
    using GFGradient =
      Dune::PDELab::VectorDiscreteGridFunctionGradient<GridFunctionSpace,
                                                       DOFVector>;
    GFGradient grad(*gfs, *dof_vector); // Breaks here!

    // Try evaluation
    for (auto&& e : Dune::elements(grid->leafGridView()))
    {
      const auto pos = e.geometry().local(e.geometry().center());
      typename GFGradient::Traits::RangeType y;
      grad.evaluate(e, pos, y);
    }

    return 0;
  }
  catch (Dune::Exception& e)
  {
    std::cout << e << std::endl;
    return 1;
  }
  catch (...)
  {
    return 1;
  }
}
