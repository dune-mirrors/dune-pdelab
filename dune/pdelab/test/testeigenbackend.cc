// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <string>

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/operators.hh>

#include <dune/pdelab.hh>

#include "gridexamples.hh"

//===============================================================
//===============================================================
// Solve the Poisson equation
//           - \Delta u = f in \Omega,
//                    u = g on \partial\Omega_D
//  -\nabla u \cdot \nu = j on \partial\Omega_N
//===============================================================
//===============================================================

//===============================================================
// Define parameter functions f,g,j and \partial\Omega_D/N
//===============================================================

template<typename GV, typename RF>
class PoissonModelProblem
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
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    const auto& xglobal = e.geometry().global(x);
    if (xglobal[0]>0.25 && xglobal[0]<0.375 && xglobal[1]>0.25 && xglobal[1]<0.375)
      return 50.0;
    else
      return 0.0;
  }

  //! boundary condition type function
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);

    if (xglobal[1]<1E-6 || xglobal[1]>1.0-1E-6)
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
      }
    if (xglobal[0]>1.0-1E-6 && xglobal[1]>0.5+1E-6)
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
      }
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    xglobal -= 0.5;
    return exp(-xglobal.two_norm2());
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);

    if (xglobal[0] > 1.0 - 1E-6 && xglobal[1] > 0.5 + 1E-6) {
      return -5.0;
    } else {
      return 0.0;
    }
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }
};

//===============================================================
// Problem setup and solution
//===============================================================

// generate a P1 function and output it
template<typename GV, typename FEM, typename CON>
void poisson (const GV& gv, const FEM& fem, std::string filename, int q)
{
  // constants and types
  typedef typename FEM::Traits::FiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<
    GV,
    FEM,
    CON,
    Dune::PDELab::Eigen::VectorBackend
    > GFS;
  GFS gfs(gv,fem);
  gfs.name("solution");

  // make model problem
  typedef PoissonModelProblem<GV,R> Problem;
  Problem problem;

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<R>::Type C;
  C cg;
  cg.clear();
  Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> bctype(gv,problem);
  Dune::PDELab::constraints(bctype,gfs,cg);

  // make local operator
  typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,FEM> LOP;
  LOP lop(problem);

  // matrix backend
  typedef Dune::PDELab::Eigen::MatrixBackend<> MBE;
  int avgnz = 27;
  MBE mbe(avgnz);

  // make grid operator
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,
                                     MBE,
                                     double,double,double,
                                     C,C> GridOperator;
  GridOperator gridoperator(gfs,cg,gfs,cg,lop,mbe);

  Dune::Timer timer;

  // make coefficent Vector and initialize it from a function
  // There is some weird shuffling around here - please leave it in,
  // it's there to test the copy constructor and assignment operator of the
  // matrix wrapper
  typedef typename GridOperator::Traits::Domain DV;
  DV x0(gfs,Dune::PDELab::Backend::unattached_container());
  {
    DV x1(gfs);
    DV x2(x1);
    x2 = 0.0;
    x0 = x1;
    x0 = x2;
  }

  // initialize DOFs from Dirichlet extension
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g(gv,problem);
  Dune::PDELab::interpolate(g,gfs,x0);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x0);

  // represent operator as a matrix
  // There is some weird shuffling around here - please leave it in,
  // it's there to test the copy constructor and assignment operator of the
  // matrix wrapper
  typedef typename GridOperator::Traits::Jacobian M;
  M m;
  {
    M m1(gridoperator);
    M m2(m1);
    m2 = 0.0;
    m = m1;
    m = m2;
  }
  gridoperator.jacobian(x0,m);
  std::cout << "two-norm: " << x0.two_norm() << std::endl
            << "one-norm: " << x0.one_norm() << std::endl
            << "inf-norm: " << x0.infinity_norm() << std::endl
            << "dot-prod: " << x0.dot(x0) << std::endl;

  // evaluate residual w.r.t initial guess
  typedef typename GridOperator::Traits::Range RV;
  RV r(gfs);
  r = 0.0;
  gridoperator.residual(x0,r);

  // test the container interface
  {
    // make ISTL solver
    Dune::MatrixAdapter<M,DV,RV> opa(m);
    // only Richardson - everything else would need ISTL support
    Dune::Richardson<DV,RV> richardson(1.0);
  }

  DV x(gfs,0.0);

  // test the solver backend
  {
    M m(gridoperator);
    // make EIGEN solver
    Dune::PDELab::EigenBackend_BiCGSTAB_Diagonal solver;
    x = 0.0;
    gridoperator.jacobian(x0,m);
    // solve the jacobian system
    r *= -1.0; // need -residual
    solver.apply(m,x,r,1e-10);
    x += x0;
    std::cout << "Time manually: " << timer.elapsed() << std::endl;
  }

  // test the solver backend as part of a pdelab solver
  {
    Dune::Timer timer;
    std::cout << "StationaryLinearProblemSolver" << std::endl;
    typedef Dune::PDELab::EigenBackend_BiCGSTAB_Diagonal LS;
    LS linearSolver(5000);
    x = 0.0;
    Dune::PDELab::StationaryLinearProblemSolver<GridOperator,LS,DV> slp(gridoperator,linearSolver,x,1e-10);
    slp.apply();
    x += x0;
    std::cout << "Time SLPSolver: " << timer.elapsed() << std::endl;
  }

  // test the CG solver backend as part of a pdelab solver
  {
    Dune::Timer timer;
    std::cout << "StationaryLinearProblemSolver" << std::endl;
    typedef Dune::PDELab::EigenBackend_CG_Diagonal_Up LS;
    LS linearSolver(5000);
    x = 0.0;
    Dune::PDELab::StationaryLinearProblemSolver<GridOperator,LS,DV> slp(gridoperator,linearSolver,x,1e-10);
    slp.apply();
    x += x0;
    std::cout << "Time SLPSolver: " << timer.elapsed() << std::endl;
  }

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,x);
  vtkwriter.write(filename,Dune::VTK::ascii);
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // YaspGrid Q1 2D test
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      std::array<int,2> N(Dune::filledArray<2,int>(1));
      Dune::YaspGrid<2> grid(L,N);
      grid.globalRefine(5);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
      FEM fem(gv);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"eigenbackend_yasp_Q1_2d",2);
    }

    // YaspGrid Q2 2D test
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      std::array<int,2> N(Dune::filledArray<2,int>(1));
      Dune::YaspGrid<2> grid(L,N);
      grid.globalRefine(5);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,2> FEM;
      FEM fem(gv);

      // and again with the space matrix
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"eigenbackend_yasp_Q2_2d",2);
    }

    // test passed
    return 0;
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
