// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>
#include <dune/common/test/testsuite.hh>

/** \brief verify consistency of grid communication and cached by checking the global indices

  1. compute communication pattern
  2. save vector with global ideas
  3. cached exchange
  4. compare buffer with local data in "scatter"
 */
template<typename Grid, typename FS>
void checkConsistentCommnunication(const Grid& grid, Dune::TestSuite & test, FS & fs)
{
  using namespace Dune::PDELab;
  using V = typename FS::DOF;
  // setup piece wise constant FEM

  auto comm = grid.comm();

  using GFS = typename FS::GFS;
  using CC = typename FS::CC;
  using LFS = LocalFunctionSpace<GFS>;
  using IndexCache = LFSIndexCache<LFS,CC>;
  GFS & gfs = fs.getGFS();
  LFS lfs(gfs);
  IndexCache cache(lfs,fs.getCC());

  const auto & gid = grid.globalIdSet();
  using ID = typename std::decay_t<decltype(gid)>::IdType;

  using GIDVec = Backend::Vector<GFS,ID>;
  GIDVec v(gfs, ID(99999));

  for (auto && e: elements(grid.leafGridView()))
  {
    lfs.bind(e);
    cache.update();

    for (unsigned int n = 0; n < lfs.size(); n++) {
      auto i = cache.containerIndex(n);
      v[i] = gid.id(e);
    }
  }

  using DOFMapper = ISTL::EntityDOFMapper<GFS>;
  using CachedComm = Dune::ExchangeCommunication< typename DOFMapper::GlobalKeyType >;

  DOFMapper dofmapper(gfs);
  auto all_all_pattern =
    buildCommunicationPatternFromMapper(gfs.gridView(), dofmapper, Dune::All_All_Interface);

  /* check pattern entries */
  auto & torus = grid.torus();
  std::cout << "torus " << torus.dims()[0] << " " << torus.dims()[1] << std::endl;
  for (auto && link : all_all_pattern)
  {
    std::cout << comm.rank()
              << "\ttorus position "
              << torus.coord()[0] << " " << torus.coord()[1]
              << " / "
              << torus.rank_to_coord(link.first)[0]  << " " << torus.rank_to_coord(link.first)[1]  << std::endl;
  }
  // assert(all_all_pattern.size() == 1);
  // comm.barrier();
  // std::cout << "-----------------------\n";
  // for (auto && entry : all_all_pattern.begin()->second.sendIndices)
  //   std::cout << comm.rank() << "\t" << "to   " << all_all_pattern.begin()->first << " ... " << entry << std::endl;
  // comm.barrier();
  // std::cout << "-----------------------\n";
  // for (auto && entry : all_all_pattern.begin()->second.recvIndices)
  //   std::cout << comm.rank() << "\t" << "from " << all_all_pattern.begin()->first << " ... " << entry << std::endl;

  /* and try to communicate global indices */
  CachedComm all_all_comm;
  all_all_comm.init(gfs.gridView().comm());
  all_all_comm.setCommunicationPattern(all_all_pattern);
  all_all_comm.exchange( Backend::native(v),
    [&test](const auto& a, const auto& b) {
      #warning this needs to be improved ...
      using A = std::decay_t<decltype(a)>;
      if constexpr (std::is_same<A,Dune::bigunsignedint<55>>{}) {
        test.check(b == a, "consistend global indices in cached communication")
          << " expected " << b << " got " << a;
      }
      else{
        // this should never be called, but currently the compiler can possibly
        // instantiate this block
        assert(false);
      }
    }
    );
}

/** Parameter class for the stationary convection-diffusion equation of the following form:
 *
 * \f{align*}{
 *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u &=& f \mbox{ in } \Omega,  \ \
 *                                              u &=& g \mbox{ on } \partial\Omega_D (Dirichlet)\ \
 *                (b(x,u) - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N (Flux)\ \
 *                        -(A(x)\nabla u) \cdot n &=& o \mbox{ on } \partial\Omega_O (Outflow)
 * \f}
 * Note:
 *  - This formulation is valid for velocity fields which are non-divergence free.
 *  - Outflow boundary conditions should only be set on the outflow boundary
 *
 * The template parameters are:
 *  - GV a model of a GridView
 *  - RF numeric type to represent results
 */
template<typename GV, typename RF>
class GenericAdvectionProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  GenericAdvectionProblem ()
  {
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 0 : 0;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      v[i] = 0.0;
    v[0] = 1.0;
    v[1] = 1.0/3.0;
    b0 = 0.25;
  }

  //! tensor diffusion constant per cell? return false if you want more than one evaluation of A per cell.
  static constexpr bool permeabilityIsConstantPerCell()
  {
    return true;
  }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! boundary condition type function
  /* return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet for Dirichlet boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann for flux boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow for outflow boundary conditions
   */
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& xlocal) const
  {
    if (is.outerNormal(xlocal)*v<0)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::DomainType x = e.geometry().global(xlocal);
    if (x[1] < v[1]*x[0]+b0)
      return 0.0;
    else
      return 1.0;
  }

  //! flux boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

private:
  typename Traits::PermTensorType I;
  typename Traits::RangeType v;
  typename Traits::RangeFieldType b0;
};

// Solve problem
template <typename Grid, typename FS, typename Problem, Dune::SolverCategory::Category solvertype, int degree>
void solveProblem(const Grid& grid, FS& fs, typename FS::DOF& x, Problem& problem, std::string basename)
{
  // initialize DOF vector it with boundary condition
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
  BCType bctype(grid.leafGridView(),problem);
  typedef typename Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g(grid.leafGridView(),problem);
  Dune::PDELab::interpolate(g,fs.getGFS(),x);

  // assemble constraints
  fs.assembleConstraints(bctype);
  fs.setNonConstrainedDOFS(x,0.0);

  // assembler for finite elemenent problem
  typedef typename Dune::PDELab::ConvectionDiffusionDG<Problem,typename FS::FEM> LOP;
  LOP lop(problem,Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,2.0);
  typedef typename Dune::PDELab::GalerkinGlobalAssembler<FS,LOP,solvertype> ASSEMBLER;
  ASSEMBLER assembler(fs,lop,20);

  // make linear solver and solve problem
  typedef typename Dune::PDELab::ISTLSolverBackend_IterativeDefault<FS,ASSEMBLER,solvertype> SBE;
  SBE sbe(fs,assembler,2,1);
  typedef typename Dune::PDELab::StationaryLinearProblemSolver<typename ASSEMBLER::GO,typename SBE::LS,typename FS::DOF> SLP;
  SLP slp(*assembler,*sbe,x,1e-6);
  slp.apply();

  // output grid to VTK file
  Dune::SubsamplingVTKWriter<typename Grid::LeafGridView> vtkwriter(grid.leafGridView(),Dune::refinementIntervals(degree));
  typename FS::DGF xdgf(fs.getGFS(),x);
  vtkwriter.addVertexData(std::make_shared<typename FS::VTKF>(xdgf,"x_h"));
  vtkwriter.write(basename,Dune::VTK::appendedraw);
}

/** \brief run system test, comparing the solution of a linear transport problem using grid & cached communication */
template<Dune::SolverCategory::Category solvertype, int degree, typename Grid, typename FS>
void checkFullProblem(const Grid& grid, Dune::TestSuite & test, FS & fs)
{
  // DOF vector
  typedef typename FS::DOF X;
  X x(fs.getGFS(),0.0);
  X x2(fs.getGFS(),0.0);

  // make problem parameters
  using NumberType = typename FS::NT;
  typedef GenericAdvectionProblem<typename Grid::LeafGridView,NumberType> Problem;
  Problem problem;

  // solve problem
  std::cout << "Running (parallel) solve with grid communication\n";
  Dune::PDELab::ISTL::ParallelHelper<typename FS::GFS>::useCaches = false;
  solveProblem<Grid,FS,Problem,solvertype,degree>(grid,fs,x,problem,"dg-grid-comm");

  std::cout << "Running (parallel) solve with cached communication\n";
  Dune::PDELab::ISTL::ParallelHelper<typename FS::GFS>::useCaches = true;
  solveProblem<Grid,FS,Problem,solvertype,degree>(grid,fs,x2,problem,"dg-cached-comm");

  // calculate l2 error squared between the two functions
  using DGF = typename FS::DGF;
  DGF xdgf(fs.getGFS(),x);
  DGF xdgf2(fs.getGFS(),x2);
  typedef Dune::PDELab::DifferenceSquaredAdapter<DGF,DGF> DifferenceSquared;
  DifferenceSquared differencesquared(xdgf,xdgf2);
  typename DifferenceSquared::Traits::RangeType l2errorsquared(0.0);
  Dune::PDELab::integrateGridFunction(differencesquared,l2errorsquared,10);

  // global sum of local error
  auto comm = fs.getGFS().gridView().comm();
  l2errorsquared = comm.sum(l2errorsquared);
  test.check(l2errorsquared < 1e-14)
    << "solution with cached comm differs from solution with grid comm: "
    << "L2 error squared = " << l2errorsquared;
}

int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

  // command line args
  int cells=10; if (argc>=2) sscanf(argv[1],"%d",&cells);

  // define parameters
  const unsigned int dim = 2;
  const unsigned int degree = 1;
  const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::cube;
  const Dune::SolverCategory::Category solvertype = Dune::SolverCategory::overlapping;
  typedef double NumberType;

  // make grid
  typedef Dune::YaspGrid<dim> Grid;
  Grid grid({1.0,1.0}, {cells,cells});

  // make a finite element space
  using FS = Dune::PDELab::DGQkSpace<Grid,NumberType,degree,elemtype,solvertype>;
  FS fs(grid.leafGridView());

  // run index test
  Dune::TestSuite test;
  checkConsistentCommnunication(grid,test,fs);

  // run advanced test
  checkFullProblem<solvertype, degree>(grid,test,fs);

  // done
  return test.exit();
}
