#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <cmath>
#include <algorithm>
#include <memory>
#include <bitset>
#include <array>
#include <string>
#include <sstream>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/pdelab.hh>
#pragma GCC diagnostic pop


// symbolic constants for the dimension of the domain
constexpr int TWODIM { 2 };
constexpr int THREEDIM { 3 };


// Comparison with 1.0 for floating point types
template<typename T>
bool is_numeric_one(const T& x, double eps=1e-8)
{
  return (std::abs(x - 1.0) <= std::max(std::abs(x), 1.0) * eps);
}


std::string generate_filename(const Dune::ParameterTree& ptree)
{
  const int dim { ptree.get<int>("grid.dim") };
  const std::string puType { ptree.get("partition.type", "") };

  if (!(puType == "standard" || puType == "sarkis"))
    DUNE_THROW(Dune::Exception, "Unknown partition of unity type.");

  std::string typeAbbrv { (puType == "standard") ? "std" : "sar" };

  std::stringstream ss;

  ss << "pu_" << typeAbbrv.c_str() << dim << "d"
     << ptree.get<int>("partition.baseDegree") << "b"
     << ptree.get<int>("grid.overlap") << "o"
     << ptree.get<int>("grid.nCellsX") << "x"
     << ptree.get<int>("grid.nCellsY");

  if (dim == TWODIM)
    ss << "c";
  else
    ss << "x" << ptree.get<int>("grid.nCellsZ") << "c";

  ss << ptree.get<int>("grid.nSubdomX") << "x"
     << ptree.get<int>("grid.nSubdomY");

  if (dim == TWODIM)
    ss << "s";
  else
    ss << "x" << ptree.get<int>("grid.nSubdomZ") << "s";

  return ss.str();
}


// The Sarkis partition of unity needs information about the size and structure
// of the grid, so we use this container class for easier passing of arguments
template<typename GVCoord>
class GridParameters2d
{
public:
  const int ovlp;
  const Dune::FieldVector<GVCoord, TWODIM> upperRight;
  const std::bitset<TWODIM> isPeriodic;
  const std::array<int, TWODIM> nCells;
  const std::array<int, TWODIM> subdomLayout;

  GridParameters2d(const Dune::ParameterTree& ptree) :
    ovlp { ptree.get<int>("grid.overlap") },
    upperRight { ptree.get<double>("grid.LX"),
                ptree.get<double>("grid.LY") },
    isPeriodic("00"),
    nCells { ptree.get<int>("grid.nCellsX"),
            ptree.get<int>("grid.nCellsY") },
    subdomLayout { ptree.get<int>("grid.nSubdomX"),
                  ptree.get<int>("grid.nSubdomY") }
  {}
};


// In 3 dimensions, Sarkis partition of unity is not yet supported. We still
// provide the Grid Parameters container for consistency
template<typename GVCoord>
class GridParameters3d
{
public:
  const int ovlp;
  const Dune::FieldVector<GVCoord, THREEDIM> upperRight;
  const std::bitset<THREEDIM> isPeriodic;
  const std::array<int, THREEDIM> nCells;
  const std::array<int, THREEDIM> subdomLayout;

  GridParameters3d(const Dune::ParameterTree& ptree) :
    ovlp { ptree.get<int>("grid.overlap") },
    upperRight { ptree.get<double>("grid.LX"),
                ptree.get<double>("grid.LY"),
                ptree.get<double>("grid.LZ") },
    isPeriodic("000"),
    nCells { ptree.get<int>("grid.nCellsX"),
            ptree.get<int>("grid.nCellsY"),
            ptree.get<int>("grid.nCellsZ") },
    subdomLayout { ptree.get<int>("grid.nSubdomX"),
                  ptree.get<int>("grid.nSubdomY"),
                  ptree.get<int>("grid.nSubdomZ") }
  {}
};


template<class GV, class GFS, class V>
void write_vtk(const std::string& filename, const GV& gv, const GFS& gfs,
               const V& partUnity, const int ord)
{
  using DGF = Dune::PDELab::DiscreteGridFunction<GFS, V>;
  DGF partUnityDGF(gfs, partUnity);

  // prepare the VTKWriter and write to file
  int subsampling { ord };
  using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
  VTKWRITER vtkwriter(gv, Dune::refinementIntervals(subsampling));

  using VTKF = Dune::PDELab::VTKGridFunctionAdapter<DGF>;

  // write partition of unity to file
  std::string outputname { "partition_of_unity" };

  vtkwriter.addVertexData(
    std::shared_ptr<VTKF>(new VTKF(partUnityDGF, outputname)));
  vtkwriter.write(filename, Dune::VTK::appendedraw);

  std::cout << "Partition of unity written to vtk file." << std::endl;
  return;
}


// sanity check: do the partition of unity functions add up to 1 at every
// point? Or more precisely:
// \sum\limits_{j=1}^N R_j^T \Theta_j(v|_{\Omega_j}) = v, \forall v \in V_h
// ? For this, it suffices to verify above relation for the basis functions
template<class GFS, class V>
void perform_sanity_check(const GFS& gfs, V& partUnity)
{
  Dune::PDELab::AddDataHandle<GFS, V> partUnityHandle(gfs, partUnity);
  gfs.gridView().communicate(partUnityHandle, Dune::All_All_Interface,
                             Dune::ForwardCommunication);

  // in theory we could also write the partition of unity after communication
  // to vtk file and check against a reference. However this would
  // unnecessarily clutter the directory with reference vtks containing only
  // ones, so we iterate through the vector and check at runtime.
  for (const auto& entry : partUnity)
  {
    if (!is_numeric_one(entry))
    {
      DUNE_THROW(Dune::Exception,"Sanity check failed! "
                 "Partition of unity does not add up to 1.0.");
    }
  }

  std::cout << "Sanity check succeeded." << std::endl;
  return;
}


/*
 * The sarkisPartitionOfUnity method from geneo/partitionofunity.hh only works
 * in 2 dimensions. In 3 dimensions compilation fails due to invalid type
 * conversion. So, to avoid compiling the sarkis partition of unity in 3
 * dimensions, we split the test routine into two cases: The standard partition
 * of unity in standardTestdriver that works for any dimension, and the
 * sarkisTestdriver that is only called in the 2 dimensional case.
 */


template<class GV, class FEM>
void standardTestdriver(const GV& gv, const FEM& fem, const int ord,
                        const std::string& filename)
{
  using Dune::PDELab::Backend::native;
  using RangeField = double;

  using VBE = Dune::PDELab::ISTL::VectorBackend<>;

  // for partitions of unity, we need constraints on the subdomain boundaries
  using CON = Dune::PDELab::OverlappingConformingDirichletConstraints;

  using GFS = Dune::PDELab::GridFunctionSpace<GV, FEM, CON, VBE>;
  GFS gfs(gv, fem);

  using PartUnityCC =
    typename GFS::template ConstraintsContainer<RangeField>::Type;
  auto cc { PartUnityCC() };

  // exclude constraints at the domain boundary
  Dune::PDELab::NoDirichletConstraintsParameters bctype;

  // enforce constraints
  Dune::PDELab::constraints(bctype, gfs, cc);

  // create partition of unity
  using V = Dune::PDELab::Backend::Vector<GFS, RangeField>;
  auto partUnity { std::make_shared<V>(standardPartitionOfUnity<V>(gfs, cc)) };

  write_vtk(filename, gv, gfs, *partUnity, ord);
  perform_sanity_check(gfs, *partUnity);

  return;
}


template<class GV, class FEM, class GP>
void sarkisTestdriver(const GV& gv, const FEM& fem, const GP& gp,
                      const int ord, const std::string& filename)
{
  using Dune::PDELab::Backend::native;
  using RangeField = double;

  using VBE = Dune::PDELab::ISTL::VectorBackend<>;

  // for partitions of unity, we need constraints on the subdomain boundaries
  using CON = Dune::PDELab::OverlappingConformingDirichletConstraints;

  using GFS = Dune::PDELab::GridFunctionSpace<GV, FEM, CON, VBE>;
  GFS gfs(gv, fem);

  using PartUnityCC =
    typename GFS::template ConstraintsContainer<RangeField>::Type;
  auto cc { PartUnityCC() };

  // exclude constraints at the domain boundary
  Dune::PDELab::NoDirichletConstraintsParameters bctype;

  // enforce constraints
  Dune::PDELab::constraints(bctype, gfs, cc);

  // create an auxiliary local function space for the construction of the
  // Sarkis partition of unity
  using LFS =
    Dune::PDELab::LocalFunctionSpace<GFS, Dune::PDELab::AnySpaceTag>;
  LFS lfs(gfs);

  using V = Dune::PDELab::Backend::Vector<GFS, RangeField>;
  auto partUnity { std::make_shared<V>(sarkisPartitionOfUnity<V>(
                                         gfs, lfs, cc, gp.nCells.at(0),
                                         gp.nCells.at(1), gp.ovlp,
                                         gp.subdomLayout.at(0),
                                         gp.subdomLayout.at(1))) };

  write_vtk(filename, gv, gfs, *partUnity, ord);
  perform_sanity_check(gfs, *partUnity);

  return;
}


int main(int argc, char** argv)
{
  try
  {
    Dune::MPIHelper& mpihelper { Dune::MPIHelper::instance(argc, argv) };

    // open ini file
    Dune::ParameterTree ptree;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree("test-geneo-partitionofunity.ini", ptree);
    ptreeparser.readOptions(argc, argv, ptree);

    // read relevant parameters from ini file
    const int dim { ptree.get<int>("grid.dim") };
    const int ord { ptree.get<int>("partition.baseDegree") };
    const std::string partUnityType { ptree.get("partition.type", "") };

    std::string filename { generate_filename(ptree) };

    if ((partUnityType == "sarkis") && (dim == THREEDIM))
    {
      DUNE_THROW(Dune::Exception, "Sarkis partition of unity is not yet "
                 "supported in 3 dimensions.");
    }

    using Real = double;
    using GVCoord = Real;

    if (dim == TWODIM)
    {
      GridParameters2d<GVCoord> gp(ptree);

      // grid partitioning
      using Partitioner = Dune::YaspFixedSizePartitioner<TWODIM>;
      auto partitioner { std::make_unique<Partitioner>(gp.subdomLayout) };

      // instantiate grid
      using Grid = Dune::YaspGrid<TWODIM>;
      auto grid { std::make_shared<Grid>(
                    gp.upperRight, gp.nCells, gp.isPeriodic, gp.ovlp,
                    Dune::MPIHelper::getCollectiveCommunication(),
                    partitioner.get()) };

      // create GridView to finest level
      using GV = Grid::LevelGridView;
      auto gv { grid->levelGridView(grid->maxLevel()) };

      if (ord == 1)
      {
        Dune::PDELab::QkLocalFiniteElementMap<GV, GVCoord, Real, 1> feMap(gv);

        if (partUnityType == "standard")
          standardTestdriver(gv, feMap, ord, filename);
        else if (partUnityType == "sarkis")
          sarkisTestdriver(gv, feMap, gp, ord, filename);
        else
          DUNE_THROW(Dune::Exception, "Unknown partition of unity type.");
      }
      else if (ord == 2)
      {
        Dune::PDELab::QkLocalFiniteElementMap<GV, GVCoord, Real, 2>  feMap(gv);

        if (partUnityType == "standard")
          standardTestdriver(gv, feMap, ord, filename);
        else if (partUnityType == "sarkis")
          sarkisTestdriver(gv, feMap, gp, ord, filename);
        else
          DUNE_THROW(Dune::Exception, "Unknown partition of unity type.");
      }
      else
        DUNE_THROW(Dune::Exception, "Degree higher than 3 not yet supported.");
    }
    else if (dim == THREEDIM)
    {
      GridParameters3d<GVCoord> gp(ptree);

      // grid partitioning
      using Partitioner = Dune::YaspFixedSizePartitioner<THREEDIM>;
      auto partitioner { std::make_unique<Partitioner>(gp.subdomLayout) };

      // instantiate grid
      using Grid = Dune::YaspGrid<THREEDIM>;
      auto grid { std::make_shared<Grid>(
                    gp.upperRight, gp.nCells, gp.isPeriodic, gp.ovlp,
                    Dune::MPIHelper::getCollectiveCommunication(),
                    partitioner.get()) };

      // create GridView to finest level
      using GV = Grid::LevelGridView;
      auto gv { grid->levelGridView(grid->maxLevel()) };

      if (ord == 1)
      {
        Dune::PDELab::QkLocalFiniteElementMap<GV, GVCoord, Real, 1> feMap(gv);

        standardTestdriver(gv, feMap, ord, filename);
      }
      else if (ord == 2)
      {
        Dune::PDELab::QkLocalFiniteElementMap<GV, GVCoord, Real, 2> feMap(gv);

        standardTestdriver(gv, feMap, ord, filename);
      }
      else
        DUNE_THROW(Dune::Exception, "Degree higher than 2 not yet supported.");
    }

  }
  catch (Dune::Exception& e)
  {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (std::exception& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << "Unknown error!" << std::endl;
    return 1;
  }

  return 0;
}
