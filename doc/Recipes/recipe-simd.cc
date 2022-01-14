#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/simd/loop.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab.hh>

/**
 * \page recipe-simd Solving a linear system with multiple right-hand sides
 *
 * Here we show how to solve a linear system of equations originating from a PDE
 * with multiple right-hand sides. I.e. Multiple linear systems with the same
 * linear operator. For that we use the SIMD datatype LoopSIMD from dune-common.
 *
 * This recipe builds up on the @ref recipe-linear-system-solution-istl.
 *
 * We define a RHS type in the `Problem` class te be a SIMD type.
 * \snippet recipe-simd.cc SIMD define RHSType
 *
 * The functions for the
 * source term \f$f\f$ as well as of the boundary conditions \f$g,j,o,\f$
 * return a SIMD vector.
 * \snippet recipe-simd.cc SIMD Source term
 * \snippet recipe-simd.cc SIMD Dirichlet boundary
 * \snippet recipe-simd.cc SIMD Neumann boundary
 * \snippet recipe-simd.cc SIMD Outflow boundary
 *
 * The `GridOperator` now takes the SIMD datatype as Domain and Range type.
 * \snippet recipe-simd.cc SIMD Grid operator
 *
 * Use a `VTKGridFunctionAdapter` to specify which lane should be written in the vtk file.
 * \snippet recipe-simd.cc Solution output
 *
 * Full example code: @ref recipe-simd.cc
 * \example recipe-simd.cc
 * See explanation at @ref recipe-simd
 */



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
class GenericEllipticProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;
  // [SIMD define RHSType]
  typedef Dune::LoopSIMD<RF, 8> RHSFieldType;
  //! [SIMD define RHSType]

  constexpr static double PI = acos(-1.);

  //! tensor diffusion constant per cell? return false if you want more than one evaluation of A per cell.
  static constexpr bool permeabilityIsConstantPerCell()
  {
    return true;
  }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::PermTensorType I;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 1 : 0;
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::RangeType v(0.0);
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    return 0.0;
  }

  // [SIMD Source term]
  //typename Traits::RangeFieldType
  RHSFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    using namespace Dune::Simd;
    RHSFieldType fx = 0.;
    typename Traits::DomainType x = e.geometry().global(xlocal);
    if(x[0] > 0.5){
      for(size_t i=0;i<lanes(fx); ++i){
        lane(i,fx) = -100*sin(2*PI*x[0]*(i+1));
      }
    }
    return fx;
  }
  //! [SIMD Source term]

  //! boundary condition type function
  /* return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet for Dirichlet boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann for flux boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow for outflow boundary conditions
   */
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& xlocal) const
  {
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  // [SIMD Dirichlet boundary]
  //! Dirichlet boundary condition value
  RHSFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    using std::sin;
    using namespace Dune::Simd;
    typename Traits::DomainType x = e.geometry().global(xlocal);
    RHSFieldType result = 0.;
    if (x[0] < 0.5)
      for(size_t l = 0; l<lanes(result); ++l)
        lane(l, result) = sin(PI*l*x[1]);
    return result;
  }
  //! [SIMD Dirichlet boundary]

  // [SIMD Neumann boundary]
  RHSFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& xlocal) const
  {
    return 0.0;
  }
  //! [SIMD Neumann boundary]

  // [SIMD Outflow boundary]
  RHSFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& xlocal) const
  {
    return 0.0;
  }
  //! [SIMD Outflow boundary]
};


int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

  // define parameters
  typedef double NumberType;

  // need a grid in order to test grid functions
  constexpr unsigned int dim = 2;
  constexpr unsigned int degree = 1;
  constexpr std::size_t nonzeros = Dune::power(2*degree+1,dim);

  Dune::FieldVector<NumberType,dim> L(1.0);
  std::array<int,dim> N(Dune::filledArray<dim,int>(50));

  typedef Dune::YaspGrid<dim> Grid;
  Grid grid(L,N);
  const auto& gv = grid.leafGridView();


  // make grid
  typedef Dune::YaspGrid<dim> GM;


  // make problem parameters
  typedef GenericEllipticProblem<typename GM::LeafGridView,NumberType> Problem;
  Problem problem;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
  BCType bctype(gv,problem);


  typedef typename GM::ctype DF;
  typedef Dune::PDELab::QkLocalFiniteElementMap<GM::LeafGridView,DF,NumberType,degree> FEM;
  FEM fem(gv);

  typedef Dune::PDELab::GridFunctionSpace<GM::LeafGridView,FEM,
  Dune::PDELab::ConformingDirichletConstraints,
  Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,1>> GFS;
  GFS gfs(gv,fem);
  gfs.name("solution");

  typedef typename GFS::template ConstraintsContainer<NumberType>::Type CC;
  CC cc;
  Dune::PDELab::constraints(bctype,gfs,cc);

  // local operator for finite elemenent problem
  typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,FEM> LOP;
  LOP lop(problem);

  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  // [SIMD Grid operator]
  typedef typename Problem::RHSFieldType DomainFieldType;
  typedef typename Problem::RHSFieldType RangeFieldType;
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,DomainFieldType,RangeFieldType,NumberType,CC,CC> GO;
  auto go = GO(gfs,cc,gfs,cc,lop,MBE(nonzeros));
  //! [SIMD Grid operator]

  // [Make degree of freedom vector]
  typedef GO::Traits::Domain X;
  X x(gfs,0.0);
  //! [Make degree of freedom vector]
  // [Set it to match boundary conditions]
  using namespace std::placeholders;
  auto g = std::bind(&Problem::g, problem, _1, _2);
  auto ggf = Dune::PDELab::makeGridFunctionFromCallable(gv, g);
  Dune::PDELab::interpolate(ggf, gfs, x);
  //! [Set it to match boundary conditions]

  //! [Create solver]
  typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS;
  LS ls;
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,X> SLP;
  SLP slp(go,ls,x,1e-10);
  slp.apply(); // here all the work is done!
  // [Create solver]

  // [Solution output]
  Dune::VTKWriter<std::decay_t<decltype(gv)>> vtkwriter(gv);
  using DGF = Dune::PDELab::DiscreteGridFunction<GFS, X>;
  DGF dgf_x(gfs, x);
  for(size_t l = 0; l < Dune::Simd::lanes<typename X::field_type>(); ++l){
    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF>>(dgf_x, "solution"+std::to_string(l), l));
  }
  vtkwriter.write("recipe-simd",Dune::VTK::appendedraw);
  //! [Solution output]
}
