// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <memory>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/std/make_array.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/localfunctions/lagrange.hh>
#include <dune/pdelab/finiteelement/qkdglagrange.hh>
#include <dune/pdelab/finiteelement/empty.hh>
#include <dune/pdelab/finiteelementmap/variablefem.hh>

#include <dune/functions/common/signature.hh>

#include "../function/callableadapter.hh"
#include "../gridfunctionspace/gridfunctionspace.hh"
#include "../backend/istl.hh"
#include "../gridfunctionspace/gridfunctionspaceutilities.hh"
#include "../gridfunctionspace/interpolate.hh"

template<class Backend, class GV, typename OrderingTag, typename GetFEM, typename F, typename R>
int test_variable_fem (const GV& gv, const Backend & backend, const OrderingTag & orderingTag, GetFEM && getFEM, const F && f, const R && r)
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

    // coefficient vector
    using Vec = Dune::PDELab::Backend::Vector<GFS,double>;
    Vec x(gfs);

    // interpolate f
    auto _f = Dune::PDELab::makeGridFunctionFromCallable(gv,f);
    Dune::PDELab::interpolate(_f,gfs,x);

    const auto & v = Dune::PDELab::Backend::native(x);
    int BCxBS = v.size() * v[0].size();
    int expected = 4*8*9;
    std::cout << "vector " << Dune::className(v) << "\n"
              << "   with " << v.size() << " blocks of size " << v[0].size() << "\n"
              << "   => " << BCxBS << " entries\n"
              << "expected " << expected << " DOFs\n"
              << "gfs has " << gfs.ordering().size() << " DOFs" << std::endl;

    std::cout << "gfs " << gfs.ordering().containerSize({}) << " blocks" << std::endl;

    if (BCxBS != expected)
        return 1;
    return 0;

    // Error:
    // Dune reported error: RangeError [resize:/home/christi/Uni/Dune/functions/dune/functions/backends/istlvectorbackend.hh:204]: Can't resize non-resizable entry v[(  - 64)] of size 9 to size((  - 64))=18446744073709551409

    // compare with r
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

        static const int d = 2;
        using Dune::GeometryTypes::quadrilateral;
        // using dG = Dune::P0LocalFiniteElement<double,double,2>;
        using dG = Dune::QkDGLagrangeLocalFiniteElement<double,double,2,2>;
        using Q1 = Dune::LagrangeCubeLocalFiniteElement<double,double,2,1>;
        using BasisTraits = typename dG::Traits::LocalBasisType::Traits;
        using Empty = Dune::PDELab::EmptyLocalFiniteElement<BasisTraits>;
        using Element = Dune::YaspGrid<2>::template Codim<0>::Entity;
        using FiniteElementType = typename Dune::PDELab::VariableLocalFiniteElementMapTraits<BasisTraits, Element>::FiniteElementType;

        try
        {
            std::cout << "Q2/Empty, Blocking::fixed(9), DefaultLeafOrderingTag\n";
            using Backend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,9>;
            // using OrderingTag = Dune::PDELab::DefaultLeafOrderingTag;
            using OrderingTag = Dune::PDELab::DefaultLeafOrderingTag;
            result +=
                test_variable_fem<>(grid.leafGridView(),
                    Backend(),
                    OrderingTag(),
                    [](const Element & e) -> FiniteElementType
                    {
                        if (e.geometry().center()[0] < 0.5)
                            return FiniteElementType(Empty(quadrilateral));
                        // return FiniteElementType(Q1());
                        else
                            return FiniteElementType(dG()); // quadrilateral
                    },
                    [](const Coord & x) { return x[0]; },
                    [](const Coord & x) { return (x[0] < 0.5) ? 0.0 : x[0]; }
                    // [](const Coord & x) { return x[0]; }
                    );
        }
        catch (const std::exception & e)
        {
            std::cout << "ERROR: " << e.what() << std::endl;
            result++;
        }

        // try
        // {
        //     std::cout << "Q2/Empty, Blocking::fixed";
        //     using Backend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,9>;
        //     result +=
        //         test_variable_fem<>(grid.leafGridView(),
        //             Backend(),
        //             [](const Element & e) -> FiniteElementType
        //             {
        //                 if (e.geometry().center()[0] < 0.5)
        //                     return FiniteElementType(Empty(quadrilateral));
        //                 // return FiniteElementType(Q1());
        //                 else
        //                     return FiniteElementType(dG()); // quadrilateral
        //             },
        //             [](const Coord & x) { return x[0]; },
        //             [](const Coord & x) { return (x[0] < 0.5) ? 0.0 : x[0]; }
        //             // [](const Coord & x) { return x[0]; }
        //             );
        // }
        // catch (const std::exception & e)
        // {
        //     std::cout << "ERROR: " << e.what() << std::endl;
        //     result++;
        // }
    }

    if (result > 0)
        std::cout << result << " variants failed" << std::endl;

	// test passed? ==0: success, >0: error
	return result;
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
