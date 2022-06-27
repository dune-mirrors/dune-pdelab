// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <memory>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/localfunctions/lagrange.hh>
#include <dune/pdelab/finiteelement/qkdglagrange.hh>
#include <dune/pdelab/finiteelement/empty.hh>
#include <dune/pdelab/finiteelementmap/variablefem.hh>

#include <dune/functions/common/signature.hh>

#include <dune/pdelab/ordering/chunkedblockordering.hh>
#include "../function/callableadapter.hh"
#include "../gridfunctionspace/gridfunctionspace.hh"
#include "../backend/istl.hh"
#include "../gridfunctionspace/gridfunctionspaceutilities.hh"
#include "../gridfunctionspace/interpolate.hh"
#include "../ordering/chunkedblockordering.hh"

template<class Backend, class GV, typename OrderingTag, typename GetFEM, typename F, typename R>
int test_variable_fem (const GV& gv, const Backend & backend, const OrderingTag & orderingTag, GetFEM && getFEM, const F && f, const R && r,
        int expected)
{
    // create finite element map with creator FEM
    using Wrapper = typename Dune::Functions::SignatureTraits<GetFEM>::Range;
    using Element = typename GV::template Codim<0>::Entity;
    using FEM = Dune::PDELab::VariableLocalFiniteElementMap<typename Wrapper::LocalBasisTraits, Element>;
    FEM fem(getFEM);

    // create functionspace
    using CON = Dune::PDELab::NoConstraints;
    using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,Backend,OrderingTag>;
    GFS gfs(gv,fem, backend, orderingTag);

    std::cout << "  gfs has " << gfs.ordering().containerSize({}) << " blocks" << std::endl;
    std::cout << "  gfs has " << gfs.ordering().size() << " DOFs" << std::endl;

    // coefficient vector
    using Vec = Dune::PDELab::Backend::Vector<GFS,double>;
    Vec x(gfs);

    // interpolate f
    auto _f = Dune::PDELab::makeGridFunctionFromCallable(gv,f);
    Dune::PDELab::interpolate(_f,gfs,x);

    // paranoia check of size
    const auto & v = Dune::PDELab::Backend::native(x);
    int BCxBS = v.size() * v[0].size();
    std::cout << "  vector " << Dune::className(v) << "\n"
            << "     with " << v.size() << " blocks of size " << v[0].size() << "\n"
            << "     => " << BCxBS << " entries\n"
            << "  expected " << expected << " DOFs\n";

    if (BCxBS != expected)
        return 1;
    return 0;
}


template<class GV, class LeafBackend, class PowerBackend, typename OrderingTag, typename GetFEM>
int test_power_variable_fem (const GV& gv, const LeafBackend & leafbackend, const PowerBackend& powerbackend, const OrderingTag & orderingTag, GetFEM && getFEM, int expected)
{
    // create finite element map with creator FEM
    using Wrapper = typename Dune::Functions::SignatureTraits<GetFEM>::Range;
    using Element = typename GV::template Codim<0>::Entity;
    using FEM = Dune::PDELab::VariableLocalFiniteElementMap<typename Wrapper::LocalBasisTraits, Element>;
    FEM fem(getFEM);

    // create functionspace
    using CON = Dune::PDELab::NoConstraints;
    using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,LeafBackend>;

    GFS gfs(gv,fem, leafbackend);
    using PGFS = Dune::PDELab::PowerGridFunctionSpace<GFS, 2, PowerBackend, OrderingTag>;
    PGFS pgfs(gfs, gfs, powerbackend, orderingTag);

    using MI = typename PGFS::Ordering::Traits::ContainerIndex;
    MI path;
    std::cout << "  gfs() has " << pgfs.ordering().containerSize(path) << " blocks" << std::endl;
    path.push_back(0);
    std::cout << "  gfs(0) has " << pgfs.ordering().containerSize(path) << " blocks" << std::endl;
    path.push_back(0);
    std::cout << "  gfs(0,0) has " << pgfs.ordering().containerSize(path) << " blocks" << std::endl;
    std::cout << "  gfs has " << pgfs.ordering().size() << " DOFs" << std::endl;

    // coefficient vector
    using Vec = Dune::PDELab::Backend::Vector<PGFS,double>;
    Vec x(pgfs);

    const auto & v = Dune::PDELab::Backend::native(x);
    int BCxBS = v.size() * v[0].size() * v[0][0].size();
    std::cout << "  vector " << Dune::className(v) << "\n"
            << "     with blocks of sizes " << v.size() << "x"  << v[0].size() << "x"  << v[0][0].size() << "\n"
            << "     => " << BCxBS << " entries\n"
            << "  expected " << expected << " DOFs\n";

    if (BCxBS != expected)
        return 1;
    return 0;
}

template<class GV, class LBackend, class PBackend, class PPBackend, class POrderingTag, class PPOrderingTag, class GetFEM>
int test_power_power_variable_fem (const GV& gv,
    const LBackend & lbackend,
    const PBackend& pbackend,
    const POrderingTag& potag,
    const PPBackend& ppbackend,
    const PPOrderingTag& ppotag,
    GetFEM && getFEM,
    int expected)
{
    // create finite element map with creator FEM
    using Wrapper = typename Dune::Functions::SignatureTraits<GetFEM>::Range;
    using Element = typename GV::template Codim<0>::Entity;
    using FEM = Dune::PDELab::VariableLocalFiniteElementMap<typename Wrapper::LocalBasisTraits, Element>;
    FEM fem(getFEM);

    // create functionspace
    using CON = Dune::PDELab::NoConstraints;
    using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,LBackend>;

    GFS gfs(gv,fem, lbackend);
    using PGFS = Dune::PDELab::PowerGridFunctionSpace<GFS, 2, PBackend, POrderingTag>;
    PGFS pgfs(gfs, gfs, pbackend, potag);

    using PPGFS = Dune::PDELab::PowerGridFunctionSpace<PGFS, 2, PPBackend, PPOrderingTag>;
    PPGFS ppgfs(pgfs, pgfs, ppbackend, ppotag);

    std::cout << "  gfs has " << ppgfs.ordering().containerSize({}) << " blocks" << std::endl;
    std::cout << "  gfs has " << ppgfs.ordering().size() << " DOFs" << std::endl;

    // coefficient vector
    using Vec = Dune::PDELab::Backend::Vector<PPGFS,double>;
    Vec x(ppgfs);

    const auto & v = Dune::PDELab::Backend::native(x);
    using V = std::decay_t<decltype(v)>;

    if constexpr (Dune::blockLevel<V>() == 2) {
        int BCxBS = v.size() * v[0].size();
        std::cout << "  vector " << Dune::className(v) << "\n"
                << "     with blocks of sizes " << v.size() << "x"  << v[0].size() << "\n"
                << "     => " << BCxBS << " entries\n"
                << "  expected " << expected << " DOFs\n";

        if (BCxBS != expected)
            return 1;
    } else if constexpr (Dune::blockLevel<V>() == 3) {
        int BCxBSxBS = v.size() * v[0].size() * v[0][0].size();
        std::cout << "  vector " << Dune::className(v) << "\n"
                << "     with blocks of sizes " << v.size() << "x"  << v[0].size() << "x"  << v[0][0].size() << "\n"
                << "     => " << BCxBSxBS << " entries\n"
                << "  expected " << expected << " DOFs\n";

        if (BCxBSxBS != expected)
            return 1;
    } else {
        DUNE_THROW(Dune::NotImplemented, "");
    }
    return 0;
}


int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    int result = 0;

    // 2D
    {
        std::cout << "2D tests" << std::endl;
        // need a grid in order to test grid functions
        using Coord = Dune::FieldVector<double,2>;
        Coord L(1.0);
        std::array<int,2> N{8,8};
        Dune::YaspGrid<2> grid(L,N);

        // static const int d = 2;
        using Dune::GeometryTypes::quadrilateral;
        // using dG = Dune::P0LocalFiniteElement<double,double,2>;
        using dG = Dune::QkDGLagrangeLocalFiniteElement<double,double,2,2>;
        using dG1 = Dune::QkDGLagrangeLocalFiniteElement<double,double,1,2>;
        using Q1 = Dune::LagrangeCubeLocalFiniteElement<double,double,2,1>;
        using BasisTraits = typename dG::Traits::LocalBasisType::Traits;
        using Empty = Dune::PDELab::EmptyLocalFiniteElement<BasisTraits>;
        using Element = Dune::YaspGrid<2>::template Codim<0>::Entity;
        using FiniteElementType = typename Dune::PDELab::VariableLocalFiniteElementMapTraits<BasisTraits, Element>::FiniteElementType;

        try
        {
            std::cout << "=== Q1, Blocking::fixed(4), DefaultLeafOrderingTag\n";
            using Backend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,4>;
            using OrderingTag = Dune::PDELab::DefaultLeafOrderingTag;
            result +=
                test_variable_fem<>(grid.leafGridView(),
                    Backend(),
                    OrderingTag(),
                    [](const Element & e) -> FiniteElementType
                    {
                        return FiniteElementType(dG1()); // quadrilateral
                    },
                    [](const Coord & x) { return x[0]; },
                    [](const Coord & x) { return (x[0] < 0.5) ? 0.0 : x[0]; },
                    8*8*4 // expected DOF count
                    );
        }
        catch (const std::exception & e)
        {
            std::cout << "ERROR: " << e.what() << std::endl;
            result++;
        }

        try
        {
            std::cout << "=== Q2/Empty, Blocking::fixed(9), DefaultLeafOrderingTag\n";
            using Backend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,9>;
            using OrderingTag = Dune::PDELab::DefaultLeafOrderingTag;
            result +=
                test_variable_fem<>(grid.leafGridView(),
                    Backend(),
                    OrderingTag(),
                    [](const Element & e) -> FiniteElementType
                    {
                        if (e.geometry().center()[0] < 0.5)
                            return FiniteElementType(Empty(quadrilateral));
                        else
                            return FiniteElementType(dG()); // quadrilateral
                    },
                    [](const Coord & x) { return x[0]; },
                    [](const Coord & x) { return (x[0] < 0.5) ? 0.0 : x[0]; },
                    4*8*9 // expected DOF count
                    );
        }
        catch (const std::exception & e)
        {
            std::cout << "ERROR: " << e.what() << std::endl;
            result++;
        }

        try
        {
            std::cout << "=== Q2/Empty, Blocking::fixed(9), Chunked(9)\n";
            using Backend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,9>;
            using OrderingTag = Dune::PDELab::ordering::Chunked<Dune::PDELab::DefaultLeafOrderingTag>;
            result +=
                test_variable_fem<>(grid.leafGridView(),
                    Backend(),
                    OrderingTag(9),
                    [](const Element & e) -> FiniteElementType
                    {
                        if (e.geometry().center()[0] < 0.5)
                            return FiniteElementType(Empty(quadrilateral));
                        else
                            return FiniteElementType(dG()); // quadrilateral
                    },
                    [](const Coord & x) { return x[0]; },
                    [](const Coord & x) { return (x[0] < 0.5) ? 0.0 : x[0]; },
                    4*8*9
                    );
        }
        catch (const std::exception & e)
        {
            std::cout << "ERROR: " << e.what() << std::endl;
            result++;
        }

        try
        {
            std::cout << "=== Q2x2/Empty, Blocking::fixed(9)/bcrs, DefaultLeafOrderingTag/EntityBlockedOrderingTag\n";
            using LeafBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,9>;
            using PowerBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::bcrs>;
            using PowerOrderingTag = Dune::PDELab::EntityBlockedOrderingTag;
            result +=
                test_power_variable_fem<>(grid.leafGridView(),
                    LeafBackend(),
                    PowerBackend(),
                    PowerOrderingTag(),
                    [](const Element & e) -> FiniteElementType
                    {
                        if (e.geometry().center()[0] < 0.5)
                            return FiniteElementType(Empty(quadrilateral));
                        else
                            return FiniteElementType(dG()); // quadrilateral
                    },
                    /* 4*8*9*2 <- WRONG: since entity blocking happens in the outer node, every entity counts as having two blocks regartheless if one is empty */
                     8*8*9*2
                    );
        }
        catch (const std::exception & e)
        {
            std::cout << "ERROR: " << e.what() << std::endl;
            result++;
        }

        try
        {
            std::cout << "Q2x2, Blocking::fixed(9)/bcrs, DefaultLeafOrderingTag/EntityBlockedOrderingTag\n";
            using LeafBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,9>;
            using PowerBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::bcrs>;
            using PowerOrderingTag = Dune::PDELab::EntityBlockedOrderingTag;
            result +=
                test_power_variable_fem<>(grid.leafGridView(),
                    LeafBackend(),
                    PowerBackend(),
                    PowerOrderingTag(),
                    [](const Element & e) -> FiniteElementType
                    {
                        return FiniteElementType(dG()); // quadrilateral
                    },
                    8*8*9*2
                    );
        }
        catch (const std::exception & e)
        {
            std::cout << "ERROR: " << e.what() << std::endl;
            result++;
        }

        try
        {
            std::cout << "=== Q2x2/Empty, Blocking::none/bcrs, DefaultLeafOrderingTag/LexicographicOrderingTag\n";
            using LeafBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
            using PowerBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::bcrs>;
            using PowerOrderingTag = Dune::PDELab::LexicographicOrderingTag;
            result +=
                test_power_variable_fem<>(grid.leafGridView(),
                    LeafBackend(),
                    PowerBackend(),
                    PowerOrderingTag(),
                    [](const Element & e) -> FiniteElementType
                    {
                        if (e.geometry().center()[0] < 0.5)
                            return FiniteElementType(Empty(quadrilateral));
                        else
                            return FiniteElementType(dG()); // quadrilateral
                    },
                    4*8*9*2
                    );
        }
        catch (const std::exception & e)
        {
            std::cout << "ERROR: " << e.what() << std::endl;
            result++;
        }

        try
        {
            std::cout << "=== Q2DG/Q1Conforming, Blocking::fixed(9), Chunked(DefaultLeafOrderingTag)\n";
            using Backend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,9>;
            using OrderingTag = Dune::PDELab::ordering::Chunked<Dune::PDELab::DefaultLeafOrderingTag>;
            result +=
                test_variable_fem<>(grid.leafGridView(),
                    Backend(),
                    OrderingTag(9), // chunks of size 9
                    [](const Element & e) -> FiniteElementType
                    {
                        if (e.geometry().center()[0] < 0.5)
                            return FiniteElementType(Q1());
                        else
                            return FiniteElementType(dG()); // quadrilateral
                    },
                    [](const Coord & x) { return x[0]; },
                    [](const Coord & x) { return (x[0] < 0.5) ? 0.0 : x[0]; },
                    // left: 5*9*1
                    // right: 4*8*9
                    5*9*1 + 4*8*9 // expected DOF count
                    );
        }
        catch (const std::exception & e)
        {
            std::cout << "ERROR: " << e.what() << std::endl;
            result++;
        }

        try
        {
            std::cout << "=== Q2x2/Empty, Blocking::fixed(9)/none/none, DefaultLeafOrderingTag/EntityBlockedOrderingTag/EntityBlockedOrderingTag\n";
            using LeafBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,9>;
            using PBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
            using POrderingTag = Dune::PDELab::EntityBlockedOrderingTag;
            using PPBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
            using PPOrderingTag = Dune::PDELab::EntityBlockedOrderingTag;
            result +=
                test_power_power_variable_fem<>(grid.leafGridView(),
                    LeafBackend(),
                    PBackend(),
                    POrderingTag(),
                    PPBackend(),
                    PPOrderingTag(),
                    [](const Element & e) -> FiniteElementType
                    {
                        if (e.geometry().center()[0] < 0.5)
                            return FiniteElementType(Empty(quadrilateral));
                        else
                            return FiniteElementType(dG()); // quadrilateral
                    },
                    4*8*9*2*2
                    );
        }
        catch (const std::exception & e)
        {
            std::cout << "ERROR: " << e.what() << std::endl;
            result++;
        }

        try
        {
            std::cout << "=== Q2x2/Empty, Blocking::fixed(9)/none/bcrs, DefaultLeafOrderingTag/EntityBlockedOrderingTag/Chunked(EntityBlockedOrderingTag)\n";
            using LeafBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,9>;
            using PBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
            using POrderingTag = Dune::PDELab::EntityBlockedOrderingTag;
            using PPBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::bcrs>;
            using PPOrderingTag = Dune::PDELab::ordering::Chunked<Dune::PDELab::EntityBlockedOrderingTag>;
            result +=
                test_power_power_variable_fem<>(grid.leafGridView(),
                    LeafBackend(),
                    PBackend(),
                    POrderingTag(),
                    PPBackend(),
                    PPOrderingTag(4*8*2),
                    [](const Element & e) -> FiniteElementType
                    {
                        if (e.geometry().center()[0] < 0.5)
                            return FiniteElementType(Empty(quadrilateral));
                        else
                            return FiniteElementType(dG()); // quadrilateral
                    },
                    4*8*9*2*2
                    );
        }
        catch (const std::exception & e)
        {
            std::cout << "ERROR: " << e.what() << std::endl;
            result++;
        }
    }

    if (result > 0)
        std::cout << result << " variant(s) failed" << std::endl;

    return result; // test passed? ==0: success, >0: error
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
	return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
	return 1;
  }
}
