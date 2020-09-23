#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab.hh>

/*
 * Defining a Darcy problem with alternating layers of permeability and a high contrast
 */
template<typename GV, typename RF>
class GenericEllipticProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);

    RF perm1 = 1e0;
    RF perm2 = 1e0; // FIXME we want high contrast
    RF layer_thickness = 1.0 / 40.0;

    RF coeff = (int)std::floor(xglobal[0] / layer_thickness) % 2 == 0 ? perm1 : perm2;

    typename Traits::PermTensorType I;
    I[0][0] = coeff;
    I[0][1] = 0;
    I[1][0] = 0;
    I[1][1] = coeff;
    return I;
  }

  //! tensor diffusion constant per cell? return false if you want more than one evaluation of A per cell.
  static constexpr bool permeabilityIsConstantPerCell()
  {
    return true;
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
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 1.0;
  }

  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);
    if (!(xglobal[0]<1E-6 || xglobal[0]>1.0-1E-6))
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    /*if (xglobal[0] > 1.0-1E-6)
      return 1.0;
    else*/
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

template <typename V_>
struct AddGatherScatter
{
    static typename V_::value_type gather(const V_ &a, int i)
    {
        return a[i]; // I am sending my value
    }
    static void scatter(V_ &a, typename V_::value_type v, int i)
    {
        a[i] += v; // add what I receive to my value
    }
};


void driver(std::string basis_type, std::string part_unity_type, Dune::MPIHelper& helper) {

  // ~~~~~~~~~~~~~~~~~~
//  Grid set up
  // ~~~~~~~~~~~~~~~~~~
  // define parameters
  const unsigned int dim = 2;
  const unsigned int degree = 1;
  const std::size_t nonzeros = std::pow(2*degree+1,dim);
  typedef double NumberType;

  typedef Dune::UGGrid<dim> GRID;
  Dune::GridFactory<GRID> factory;
  Dune::GmshReader<GRID>::read(factory, "grids/24x24.msh", true, true);
  std::unique_ptr<GRID> grid (factory.createGrid());

  typedef typename GRID::LeafGridView GV;
  auto gv = grid->leafGridView();

  grid->loadBalance();

  // ~~~~~~~~~~~~~~~~~~
//  Type definitions
  // ~~~~~~~~~~~~~~~~~~
  const int components = 1;
  using K = double;
  using Vector = Dune::BlockVector<Dune::FieldVector<K,components>>;
  using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<K,components,components>>;

  typedef Dune::BlockVector<Dune::FieldVector<K, 1>> CoarseVector;
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<K, 1, 1>> CoarseMatrix;

  using ESExcluder = Dune::PDELab::EntitySetExcluder<Vector, GV>;
  auto ghost_excluder = std::make_shared<Dune::PDELab::EntitySetGhostExcluder<Vector, GV>>();


  using ES = Dune::PDELab::ExcluderEntitySet<GV,Dune::Partitions::All, ESExcluder>;
  ES es(gv, ghost_excluder);

  // make problem parameters
  typedef GenericEllipticProblem<ES,NumberType> Problem;
  Problem problem;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
  BCType bctype(es,problem);

  using Dune::PDELab::Backend::native;

  // ~~~~~~~~~~~~~~~~~~
//  FE fine Space definition
  // ~~~~~~~~~~~~~~~~~~
  typedef typename ES::Grid::ctype DF;
  // instantiate finite element maps
  typedef Dune::PDELab::QkLocalFiniteElementMap<ES,DF,double,1> FEM;
  FEM fem(es);
  // function space with no constraints on processor boundaries, needed for the GenEO eigenproblem
  typedef Dune::PDELab::GridFunctionSpace<ES,FEM,
                                          Dune::PDELab::ConformingDirichletConstraints,
                                          Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,1>
                                          > GFS;
  GFS gfs(es,fem);
  int verbose = 0;
  if (gfs.gridView().comm().rank()==0) verbose = 2;
  // make a degree of freedom vector on fine grid and initialize it with interpolation of Dirichlet condition
  typedef Dune::PDELab::Backend::Vector<GFS,NumberType> V;
  V x(gfs,0.0);
  // Extract domain boundary constraints from problem definition, apply trace to solution vector
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g(es,problem);
  Dune::PDELab::interpolate(g,gfs,x);
  // Set up constraints containers with boundary constraints, but without processor constraints
  typedef typename GFS::template ConstraintsContainer<NumberType>::Type CC;
  auto cc = CC();
  // assemble constraints
  Dune::PDELab::constraints(bctype,gfs,cc);
  // set initial guess
  V x0(gfs,0.0);
  Dune::PDELab::copy_nonconstrained_dofs(cc,x0,x);
  // LocalOperator for given problem
  typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,FEM> LOP;
  LOP lop(problem);
  // LocalOperator wrapper zeroing out subdomains' interiors in order to set up overlap matrix
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  // Construct GridOperators from LocalOperators
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,NumberType,NumberType,NumberType,CC,CC> GO;
  auto go = GO(gfs,cc,gfs,cc,lop,MBE(nonzeros));
  // Assemble fine grid matrix defined without processor constraints
  typedef typename GO::Jacobian M;
  M A(go);
  go.jacobian(x,A);
  // set up and assemble right hand side w.r.t. l(v)-a(u_g,v)
  V d(gfs,0.0);
  go.residual(x,d);


  // ~~~~~~~~~~~~~~~~~~
//  Solving process begin here: First some parameters
  // ~~~~~~~~~~~~~~~~~~
  double eigenvalue_threshold = -1;
  const int algebraic_overlap = 1;
  int nev = 30;
  int nev_arpack = nev;
  // double shift = 0.001;

  // ~~~~~~~~~~~~~~~~~~
  // Create offline hierarchy
  // ~~~~~~~~~~~~~~~~~~
  std::string path_to_storage = "Offline/";


  // ~~~~~~~~~~~~~~~~~~
//  ADAPTER
  // ~~~~~~~~~~~~~~~~~~
  Dune::NonoverlappingOverlapAdapter<GV, Vector, Matrix> adapter(gv, native(A), nonzeros, algebraic_overlap);
  using Attribute = Dune::EPISAttribute;
  Dune::AllSet<Attribute> allAttribute;
  auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
  allinterface->build(*adapter.getRemoteIndices(), allAttribute, allAttribute); // all to all communication
  // build up buffered communicator allowing communication over a dof vector
  auto communicator = std::shared_ptr<Dune::BufferedCommunicator>(new Dune::BufferedCommunicator());
  communicator->build<Vector>(*allinterface);

  auto geneo_matrices = setupGenEOMatrices(go, adapter, A);
  std::shared_ptr<Matrix> A_extended = std::get<0>(geneo_matrices);
  std::shared_ptr<Matrix> A_overlap_extended = std::get<1>(geneo_matrices);
  std::shared_ptr<Vector> part_unity = std::get<2>(geneo_matrices);

  auto rank = adapter.gridView().comm().rank();

  if (adapter.gridView().comm().rank()==0)
    std::cout << A.N() << std::endl;

  if (adapter.gridView().comm().rank()==1)
    std::cout << A.N() << std::endl;

  // ~~~~~~~~~~~~~~~~~~
  // Save Partition of unity
  // ~~~~~~~~~~~~~~~~~~

  std::cout << part_unity->size() << std::endl;
  std::string filename_PoU = path_to_storage + std::to_string(rank) + "_PoU.mm";
  Dune::storeMatrixMarket(*part_unity, filename_PoU, 15);

  // ~~~~~~~~~~~~~~~~~~
  // Save new2old / old2new index to be able to reproduce extension/reduction online
  // ~~~~~~~~~~~~~~~~~~

  auto new2old_localindex = adapter.get_new2old_localindex();
  auto old2new_localindex = adapter.get_old2new_localindex();

  Dune::BlockVector<Dune::FieldVector<int, 1>> n2o(new2old_localindex.size()), o2n(old2new_localindex.size());
  for (int i=0; i<new2old_localindex.size();i++) { n2o[i]=new2old_localindex[i]; }
  for (int i=0; i<old2new_localindex.size();i++) { o2n[i]=old2new_localindex[i]; }

  std::string filename_n2o = path_to_storage + std::to_string(rank) + "_new2old.mm";
  Dune::storeMatrixMarket(n2o, filename_n2o, 15);
  std::string filename_o2n = path_to_storage + std::to_string(rank) + "_old2new.mm";
  Dune::storeMatrixMarket(o2n, filename_o2n, 15);

  // ~~~~~~~~~~~~~~~~~~
  // Save neighbor ranks vector
  // ~~~~~~~~~~~~~~~~~~

  std::vector<int> neighbor_ranks = adapter.findNeighboringRanks();

  Dune::BlockVector<Dune::FieldVector<int, 1>> NR(neighbor_ranks.size());
  for (int i=0; i<neighbor_ranks.size();i++) { NR[i]=neighbor_ranks[i]; }

  std::string filename_NR = path_to_storage + std::to_string(rank) + "_neighborRanks.mm";
  Dune::storeMatrixMarket(NR, filename_NR, 15);

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Subdomain basis computations
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> subdomainbasis;
  subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GO, Matrix, Vector>>(adapter, A_extended, A_overlap_extended, part_unity, eigenvalue_threshold, nev, nev_arpack);

  // ~~~~~~~~~~~~~~~~~~
  //  Particular solution
  // ~~~~~~~~~~~~~~~~~~
  // TODO: compute the particular solution just on certain subdomains not all

  // Extend the load vector d (right hand side) of the fine to space to its virtual overlap
  Vector b(adapter.getExtendedSize());

  //~ Add a solution to another rhs here fill with ones or
  // for(auto it = b.begin(); it!=b.end(); ++it){
  //   b[it.index()] += 1.0;
    // std::srand(std::time(0));
    // b[it.index()] = -1.0 + 2.0* (std::rand()+0.0) / (RAND_MAX + 1.0);
    // if(adapter.gridView().comm().rank()==0) std::cout<< b[it.index()] << std::endl;
  // }

  // Use of the real rhs
  adapter.extendVector(native(d), b);

  communicator->forward<AddGatherScatter<Vector>>(b,b); // make function known in other subdomains

  // ~~~~~~~~~~~~~~~~~~
  // Set up boundary condition for the computation of the particular solution
  // ~~~~~~~~~~~~~~~~~~
  const int block_size = Vector::block_type::dimension;
  Matrix A_dirichlet = *A_extended;
  // Apply Dirichlet conditions to matrix on processor boundaries, inferred from partition of unity
  for (auto rIt=A_dirichlet.begin(); rIt!=A_dirichlet.end(); ++rIt){
      for(int block_i = 0; block_i < block_size; block_i++){
          if ((*part_unity)[rIt.index()][block_i] == .0){
              for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt){
                  for(int block_j = 0; block_j < block_size; block_j++){
                      (*cIt)[block_i][block_j] = (rIt.index() == cIt.index() && block_i == block_j) ? 1.0 : 0.0;
                  }
              }
          }
      }
  }
  Vector b_cpy(b);
  for (auto rIt=A_dirichlet.begin(); rIt!=A_dirichlet.end(); ++rIt) {
      for(int block_i = 0; block_i < block_size; block_i++){
          bool isDirichlet = true;
          for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt){
              for(int block_j = 0; block_j < block_size; block_j++){
                  if ((rIt.index() != cIt.index() || block_i!=block_j) && (*cIt)[block_i][block_j] != 0.0){
                      isDirichlet = false;
                      break;
                  }
              }
              if(!isDirichlet) break;
          }
          if (isDirichlet){
              b_cpy[rIt.index()] = .0;
              b[rIt.index()] = .0;
          }
      }
  }

  // ~~~~~~~~~~~~~~~~~~
  // Solve the local fine space with ovlp system
  // ~~~~~~~~~~~~~~~~~~
  Vector ui(adapter.getExtendedSize());
  Dune::UMFPack<Matrix> subdomain_solver(A_dirichlet, false);
  Dune::InverseOperatorResult result1;
  subdomain_solver.apply(ui,b_cpy,result1);
  subdomainbasis->append(ui);

  // save the subdomain basis for the online go (here to keep the particular solution)
  // subdomainbasis->to_file(basename, rank);

  // ~~~~~~~~~~~~~~~~~~
  // Save local basis sizes
  // ~~~~~~~~~~~~~~~~~~
  // Communicate local coarse space dimensions
  Dune::BlockVector<Dune::FieldVector<int, 1>> local_basis_sizes(adapter.gridView().comm().size());
  int local_size = subdomainbasis->basis_size();
  adapter.gridView().comm().allgather(&local_size, 1, local_basis_sizes.data());
  if (rank==0) {
    std::string filename_ls = path_to_storage + "localBasisSizes.mm";
    Dune::storeMatrixMarket(local_basis_sizes, filename_ls, 15);
  }

  // ~~~~~~~~~~~~~~~~~~
  // Save local basis
  // ~~~~~~~~~~~~~~~~~~
  for (int basis_index = 0; basis_index < subdomainbasis->basis_size(); basis_index++) {
    std::string filename_EV = path_to_storage + std::to_string(rank) + "_EV_" + std::to_string(basis_index) + ".mm";
    Dune::storeMatrixMarket(*subdomainbasis->get_basis_vector(basis_index), filename_EV, 15);
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Coarse space construction and saving
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  auto coarse_space = std::make_shared<Dune::PDELab::NonoverlappingSubdomainProjectedCoarseSpace<GV, Matrix, Vector>>(adapter, gv, *A_extended, subdomainbasis, verbose);
  // Get the coarse matrix from the coarse space
  std::shared_ptr<CoarseMatrix> AH = coarse_space->get_coarse_system();
  // Save the coarse space
  std::string filename = path_to_storage + "OfflineAH.mm";
  Dune::storeMatrixMarket(*AH, filename, 15);

  // // Save the sizes of local basis
  // coarse_space->local_basis_sizes_to_file(path_to_storage + "local_basis_sizes");

  // ~~~~~~~~~~~~~~~~~~
//  Coarse space right hand side and saving
  // ~~~~~~~~~~~~~~~~~~
  // Use the correct rhs to solve the final system
  // adapter.extendVector(native(d), b); // if comment: this could already been declared and filled for the particular solution
  // Initializate a load vector (right hand side) in the coarse space : coarse_d = RH * b
  CoarseVector coarse_d(coarse_space->basis_size(), coarse_space->basis_size());
  coarse_space->restrict(b,coarse_d);
  // Save the coarse rhs
  std::string filename_cb = path_to_storage + "OfflineCoarseb.mm";
  Dune::storeMatrixMarket(coarse_d, filename_cb, 15);

  // ~~~~~~~~~~~~~~~~~~
// Solve the coarse space system  ::TODO:: MAKE THAT OPTIONAL for a normal offline go
  // ~~~~~~~~~~~~~~~~~~
  // Objective :  find x = RH^T * AH^-1 * RH * b
  // Use of UMFPack to solve the problem [AH * coarse_v = RH * b] instead of inversing AH :: objective is to have [coarse_v = AH^-1 *  RH * b]
  Dune::UMFPack<CoarseMatrix> coarse_solver(*AH, false); // Is there something better than UMFPACK?
  CoarseVector coarse_v(coarse_space->basis_size(), coarse_space->basis_size());
  Dune::InverseOperatorResult result;
  coarse_solver.apply(coarse_v,coarse_d,result);

  // ~~~~~~~~~~~~~~~~~~
// Prolongate the solution in order to have v_fine_vovlp on the fine space : v_fine_vovlp = RH^T * coarse_v = RH^T * AH^-1 * RH * b
  // ~~~~~~~~~~~~~~~~~~
  Vector v_fine_vovlp(adapter.getExtendedSize());
  coarse_space->prolongate(coarse_v, v_fine_vovlp);
  communicator->forward<AddGatherScatter<Vector>>(v_fine_vovlp, v_fine_vovlp);
  V v(gfs, 0.0);
  adapter.restrictVector(v_fine_vovlp, v);

  native(x) -= v;

  // Write solution to VTK
  Dune::VTKWriter<GV> vtkwriter(gv);
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  typedef Dune::PDELab::VTKGridFunctionAdapter<DGF> ADAPT;
  auto adapt = std::make_shared<ADAPT>(xdgf,"solution");
  vtkwriter.addVertexData(adapt);
  vtkwriter.write("nonovlptestgeneo_basis_" + basis_type + "_part_unity_" + part_unity_type);

  // ~~~~~~~~~~~~~~~~~~
// Visualise all the basis in vtk format
  // ~~~~~~~~~~~~~~~~~~
  // for (int basis_index = 0; basis_index < subdomainbasis->basis_size(); basis_index++) {
  //     Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,Dune::refinementLevels(0));
  //     V vect(gfs, 0.0);
  //     adapter.restrictVector(native(*subdomainbasis->get_basis_vector(basis_index)), native(vect));

  //     int rank = adapter.gridView().comm().rank();
  //     std::string filename = "BasisVector_"+std::to_string(basis_index);

  //     Dune::PDELab::vtk::DefaultFunctionNameGenerator fieldname;
  //     fieldname.prefix("EV");

  //     Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,vect,fieldname);
  //     vtkwriter.write(filename,Dune::VTK::ascii);
  // }
}


int main(int argc, char **argv)
{
  using Dune::PDELab::Backend::native;

  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc,argv);

  driver("geneo", "standard", helper);

  return 0;
}
