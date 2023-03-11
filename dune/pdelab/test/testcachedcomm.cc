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

  if (comm.rank()==0)
    std::cout << "#### checkConsistentCommnunication ... " << std::endl;

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

    auto && fem = lfs.finiteElement();
    for (unsigned int n = 0; n < lfs.size(); n++) {
      // find sub-entity associated with this DOF
      auto && key = fem.localCoefficients().localKey(n);
      // get global id of sub-entity
      auto id = gid.subId(e, key.subEntity(), key.codim());
      // store in the vector
      auto i = cache.containerIndex(n);
      v[i] = id;
    }
  }

  using DOFMapper = ISTL::EntityDOFMapper<GFS>;
  using CachedComm = Dune::ExchangeCommunication< typename DOFMapper::GlobalKeyType >;

  DOFMapper dofmapper(gfs);
  auto [all_all_pattern, skipIndices] =
    buildCommunicationPatternFromMapper(gfs.gridView(), dofmapper, Dune::All_All_Interface);

  /* check pattern entries */
  auto & torus = grid.torus();
  if (comm.rank()==0)
    std::cout << "torus " << torus.dims()[0] << " " << torus.dims()[1] << std::endl;
  for (auto && link : all_all_pattern)
  {
    std::cout << comm.rank()
              << "\ttorus position "
              << torus.coord()[0] << " " << torus.coord()[1]
              << " / "
              << torus.rank_to_coord(link.first)[0]  << " " << torus.rank_to_coord(link.first)[1]  << std::endl;
  }

  /* and try to communicate global indices */
  try {
    CachedComm all_all_comm;
    all_all_comm.init(gfs.gridView().comm());
    all_all_comm.setCommunicationPattern(all_all_pattern);
    all_all_comm.exchange( Backend::native(v),
      [&test,&gfs](const auto& a, const auto& b) {
#warning this needs to be improved ...
        using A = std::decay_t<decltype(a)>;
        if constexpr (std::is_convertible<A,Dune::bigunsignedint<55>>{}) {
          test.check(b == a, "consistent global indices in cached communication")
            << gfs.gridView().comm().rank() << " expected " << b << " got " << a;
        }
        else{
          // this should never be called, but currently the compiler can possibly
          // instantiate this block
          test.check(false, "consistent global indices in cached communication")
            << " expected " << className(b) << " got " << className(a);
          // unsigned int X = 666;
          // throw X;
        }
      }
      );
  }
  catch (unsigned int i)
  {
  }
  comm.barrier();
}

// Poisson problem definition
template<typename GV, typename RF>
class GenericEllipticProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  //! tensor diffusion constant per cell? return false if you want more than one evaluation of A per cell.
  static constexpr bool permeabilityIsConstantPerCell()
  {
    return true;
  }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType I;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 1 : 0;
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType v(0.0);
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
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    auto x = e.geometry().global(xlocal);
    return x[0] * std::sin(5.0*M_PI*x[1]) + std::exp(-((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5)) / 0.02);
  }

  //! Boundary condition type function. Will not be evaluated on periodic boundary, so we simply set Dirichlet.
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    return 0.0;
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
};

// Solve problem
template <typename Grid, typename FS, typename Problem, Dune::SolverCategory::Category solvertype, int degree>
void solveProblem(const Grid& grid, FS& fs, typename FS::DOF& x, Problem& problem, std::string basename, double reduction, int maxiter)
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
  typedef Dune::PDELab::ConvectionDiffusionDG<Problem,typename FS::FEM> LOP;
  LOP lop(problem,Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,2.0);
  typedef Dune::PDELab::GalerkinGlobalAssembler<FS,LOP,solvertype> ASSEMBLER;
  ASSEMBLER assembler(fs,lop,27);

  // make linear solver and solve problem
  using SBE = typename Dune::PDELab::ISTLSolverBackend_IterativeDefault<FS,ASSEMBLER,solvertype>;
  SBE sbe(fs,assembler,maxiter,1);

  typedef typename Dune::PDELab::StationaryLinearProblemSolver<typename ASSEMBLER::GO,typename SBE::LS,typename FS::DOF> SLP;
  SLP slp(*assembler,*sbe,x,reduction);
  slp.apply();

  // output grid to VTK file
  Dune::SubsamplingVTKWriter<typename Grid::LeafGridView> vtkwriter(grid.leafGridView(),Dune::refinementIntervals(degree));
  typename FS::DGF xdgf(fs.getGFS(),x);
  vtkwriter.addVertexData(std::make_shared<typename FS::VTKF>(xdgf,"x_h"));
  vtkwriter.write(basename,Dune::VTK::appendedraw);
}

/** \brief run system test, comparing the solution of a linear transport problem using grid & cached communication */
template<Dune::SolverCategory::Category solvertype, int degree, typename Grid, typename FS, typename Problem>
void checkFullProblem(const Grid& grid, Dune::TestSuite & test, FS & fs, Problem & problem)
{
  auto comm = grid.comm();

  if (comm.rank()==0)
    std::cout << "#### checkFullProblem ... " << std::endl;

  // DOF vector
  typedef typename FS::DOF X;
  X x(fs.getGFS(),0.0);
  X x2(fs.getGFS(),0.0);

  // prescribed reduction
  double reduction = 1e-10;
  int maxiter = 500;

  // solve problem
  if (comm.rank()==0)
    std::cout << "---- Running (parallel) solve with grid communication\n";
  Dune::PDELab::ISTL::ParallelHelper<typename FS::GFS>::useCaches = false;
  solveProblem<Grid,FS,Problem,solvertype,degree>(grid,fs,x,problem,"dg-grid-comm",reduction,maxiter);

  if (comm.rank()==0)
    std::cout << "---- Running (parallel) solve with cached communication\n";
  Dune::PDELab::ISTL::ParallelHelper<typename FS::GFS>::useCaches = true;
  solveProblem<Grid,FS,Problem,solvertype,degree>(grid,fs,x2,problem,"dg-cached-comm",reduction,maxiter);

  // calculate l2 error squared between the two functions
  using DGF = typename FS::DGF;
  DGF xdgf(fs.getGFS(),x);
  DGF xdgf2(fs.getGFS(),x2);
  typedef Dune::PDELab::DifferenceSquaredAdapter<DGF,DGF> DifferenceSquared;
  DifferenceSquared differencesquared(xdgf,xdgf2);
  typename DifferenceSquared::Traits::RangeType l2errorsquared(0.0);
  Dune::PDELab::integrateGridFunction(differencesquared,l2errorsquared,10);

  // global sum of local error
  l2errorsquared = comm.sum(l2errorsquared);
  test.check(l2errorsquared < std::max(reduction,1e-15))
    << "solution with cached comm differs from solution with grid comm: "
    << "L2 error squared = " << l2errorsquared;
}

int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit
  auto & mpi = Dune::MPIHelper::instance(argc,argv);
  int rank = mpi.rank();

  // command line args
  int cells=10; if (argc>=2) sscanf(argv[1],"%d",&cells);

  // define parameters
  const unsigned int dim = 2;
  const Dune::PDELab::MeshType meshtype = Dune::PDELab::MeshType::conforming;
  const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::cube;
  const Dune::SolverCategory::Category solvertype = Dune::SolverCategory::overlapping;
  typedef double NumberType;

  // make grid
  typedef Dune::YaspGrid<dim> Grid;
  Grid grid({1.0,1.0}, {cells,cells});

  // make problem parameters
  typedef GenericEllipticProblem<typename Grid::LeafGridView,NumberType> Problem;
  Problem problem;

  // prepare test suite
  Dune::TestSuite test;

  {
    if (rank==0)
      std::cout << "## Check with DG Space" << std::endl;

    // make a finite element space
    const unsigned int degree = 1;
    using FS = Dune::PDELab::DGQkSpace<Grid,NumberType,degree,elemtype,solvertype>;
    FS fs(grid.leafGridView());

    // run index test
    checkConsistentCommnunication(grid,test,fs);

    // run advanced test
    checkFullProblem<solvertype, degree>(grid,test,fs,problem);
  }

  {
    if (rank==0)
      std::cout << "## Check with Lagrange Space (k=1)" << std::endl;

    // make a finite element space
    const unsigned int degree = 1;
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
    BCType bctype(grid.leafGridView(),problem);
    using FS = Dune::PDELab::CGSpace<Grid,NumberType,degree,BCType,elemtype,meshtype,solvertype>;
    FS fs(grid,bctype);

    // run index test
    checkConsistentCommnunication(grid,test,fs);

    // run advanced test
    checkFullProblem<solvertype, degree>(grid,test,fs,problem);
  }

  {
    if (rank==0)
      std::cout << "## Check with Lagrange Space (k=2)" << std::endl;

    // make a finite element space
    const unsigned int degree = 2;
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
    BCType bctype(grid.leafGridView(),problem);
    using FS = Dune::PDELab::CGSpace<Grid,NumberType,degree,BCType,elemtype,meshtype,solvertype>;
    FS fs(grid,bctype);

    // run index test
    checkConsistentCommnunication(grid,test,fs);

    // run advanced test
    checkFullProblem<solvertype, degree>(grid,test,fs,problem);
  }

  // done
  return test.exit();
}
