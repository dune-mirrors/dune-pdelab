// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#define DISABLE_LOCAL_INTERFACE

#include <iostream>
#include <string>

#include <unistd.h>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/smartpointer.hh>

#include <dune/grid/albertagrid/agrid.hh>

#ifdef HAVE_ALUGRID
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#endif

#ifdef HAVE_ALBERTA
#include <dune/grid/albertagrid/dgfparser.hh>
#endif

#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>

#include "../backend/istlmatrixbackend.hh"
#include "../backend/istlsolverbackend.hh"
#include "../backend/istlvectorbackend.hh"
#include "../common/function.hh"
#include "../common/geometrywrapper.hh"
#include "../finiteelementmap/conformingconstraints.hh"
#include "../finiteelementmap/edges02dfem.hh"
#include "../finiteelementmap/edges03dfem.hh"
#include "../gridfunctionspace/constraints.hh"
#include "../gridfunctionspace/gridfunctionspace.hh"
#include "../gridfunctionspace/gridfunctionspaceutilities.hh"
#include "../gridfunctionspace/interpolate.hh"
#include "../gridoperatorspace/gridoperatorspace.hh"
#include "../localoperator/electrodynamic.hh"

#include "divergence-probe.hh"
#include "electricenergy-probe.hh"
#include "fmt.hh"
#include "globalerror-probe.hh"
#include "gnuplotgraph.hh"
#include "gridexamples.hh"
#include "l2error-probe.hh"
#include "l2norm-probe.hh"
#include "physicalconstants.hh"
#include "probe.hh"
#include "resonatorsolution.hh"
#include "vtk-probe.hh"

//===============================================================
//===============================================================
// Solve the simple electrodynamic wave equation
//    rot(1/mu * rot E) + eps * \partial_t^2 E = 0 in \Omega, 
//                                       n x E = 0 on \partial\Omega_D
//===============================================================
//===============================================================

//
//  CONFIGURATION
//

// default limit for the total convergence
// (alberta with the triangulated unit cube uses a limit derived from this
// since albertas refinement algorithm seem to produce bad results)
const double conv_limit = 0.85;

// stop refining after the grid has more than this many elements (that means
// in 3D that the fine grid may have up to 8 times as may elements)
const unsigned maxelements = 1<<11;

// multiplier for the stepsize obtained from the FDTD criterion
const double stepadjust = 1.0/4;

// how long to run the simulation
const double duration = 10/c0;

// whether to run the check for all levels of refinement.  If false, will run
// the check only for the first and the last level.  Running the check for all
// levels is useful when debugging.
const bool do_all_levels = true;

// Order of quadrature rules to use
const unsigned int quadrature_order = 3;

// where to place a probe to measure the E-field
const std::string probe_location = ".33333333333333333333 .2 .14285714285714285714";

// whether to use the exact solution directly as reference solution or to
// interpolate it into a discrete grid function first
const bool use_interpolated_reference = false;

const double excitation_sigma = .5/c0;
const double excitation_t0 = 2/c0;
const double excitation_amp[] = {1, 0, 0};

//
//  CODE
//

//===============================================================
// Define parameter function mu
//===============================================================

// function for defining the source term
template<typename GV, typename RF, unsigned dimRange = 1>
class ConstFunc
  : public Dune::PDELab::AnalyticGridFunctionBase<
      Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dimRange>,
      ConstFunc<GV,RF,dimRange>
    >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dimRange> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<
    Traits,
    ConstFunc<GV,RF, dimRange> > BaseT;

  ConstFunc (const GV& gv, const typename Traits::RangeType& val_ = 1)
    : BaseT(gv)
    , val(val_)
  {}

  inline void
  evaluateGlobal (const typename Traits::DomainType& x, 
                  typename Traits::RangeType& y) const
  {
    y = val;
  }

private:
  typename Traits::RangeType val;
};

// boundary grid function selecting boundary conditions 
template<typename GV>
class B
  : public Dune::PDELab::BoundaryGridFunctionBase<
      Dune::PDELab::BoundaryGridFunctionTraits<
        GV,
        int,1,Dune::FieldVector<int,1>
      >,
      B<GV>
    >
{
  const GV& gv;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B<GV> > BaseT;

  B (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void
  evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig, 
            const typename Traits::DomainType& x,
            typename Traits::RangeType& y) const
  {
    y = 1; // Dirichlet everywhere
  }

  //! get a reference to the GridView
  inline const GV& getGridView () const
  {
    return gv;
  }
};

// function for initialization
template<typename GV, typename RF>
class Init
  : public Dune::PDELab::AnalyticGridFunctionBase<
      Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,3>,
      Init<GV,RF>
    >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,3> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,Init<GV,RF> > BaseT;

  Init (const GV& gv, typename Traits::DomainFieldType radius = .1)
    : BaseT(gv), center(0.5), radius2(radius*radius)
  { }

  inline void
  evaluateGlobal (const typename Traits::DomainType& x, 
                  typename Traits::RangeType& y) const
  {
    y = 0;
    if((x-center).two_norm2() < radius2)
      y[0] = -1;
  }

private:
  typename Traits::DomainType center;
  typename Traits::DomainFieldType radius2;
};

//! calculate length of smallest edge.  If no edge was found, return +inf
template<typename GV>
typename GV::Grid::ctype smallestEdge(const GV& gv)
{
  typedef typename GV::Grid::ctype ct;
  static const int dim = GV::dimension;
  static const int dimw = GV::dimensionworld;

  ct min = std::numeric_limits<ct>::infinity();
  typedef Dune::FieldVector<ct, dimw> CV;
  typedef typename GV::template Codim<0>::Iterator Iterator;
  const Iterator end = gv.template end<0>();

  for(Iterator it = gv.template begin<0>(); it!=end; ++it) {
    const Dune::GenericReferenceElement<ct, dim> &refElem
      = Dune::GenericReferenceElements<ct, dim>::general(it->type());

    for(int i = 0; i < refElem.size(dim-1); ++i) {
      CV diff = it->geometry().corner(refElem.subEntity(i,dim-1, 0,dim));
      diff   -= it->geometry().corner(refElem.subEntity(i,dim-1, 1,dim));
      ct distance = diff.two_norm();
      if(distance < min)
        min = distance;
    }
  }
  return min;
}

// find an edge near the given position
/*
 * \tparam GFS The GridFunctionSpace
 *
 * \param [in]     gfs The GridFunctionSpace.
 * \param [in,out] pos On entry: the position where to find an edge.  On exit:
 *                 the center of the edge actually found.
 * \returns The index of the edge found within the GridFunctionSpace.
 */
template<typename GFS>
unsigned edge_dof_near_pos
(const GFS& gfs,
 Dune::FieldVector<typename GFS::Traits::GridViewType::ctype,
                   GFS::Traits::GridViewType::dimensionworld>& pos)
{
  typedef typename GFS::Traits::GridViewType GV;
  typedef typename GFS::LocalFunctionSpace LFS;
  typedef typename GV::ctype DF;
  static const unsigned dimD = GV::dimensionworld;
  typedef Dune::FieldVector<DF, dimD> D;
  typedef typename GV::template Codim<0>::EntityPointer EP;
  typedef Dune::GenericReferenceElements<DF, dimD> REs;
  typedef Dune::GenericReferenceElement<DF, dimD> RE;

  Dune::HierarchicSearch<typename GV::Grid, typename GV::IndexSet>
    hsearch(gfs.gridview().grid(), gfs.gridview().indexSet());
  EP ep = hsearch.findEntity(pos);
  const typename GV::template Codim<0>::Geometry& geo = ep->geometry();
  const RE& re = REs::general(geo.type());

  D bestpos = geo.global(re.position(0,dimD-1));
  unsigned best_edge = 0;
  DF bestdist2 = (pos-bestpos).two_norm2();

  for(unsigned i = 1; i < re.size(dimD-1); ++i) {
    D thispos = geo.global(re.position(i,dimD-1));
    DF thisdist2 = (pos-thispos).two_norm2();
    if(thisdist2 < bestdist2) {
      bestpos = thispos;
      bestdist2 = thisdist2;
      best_edge = i;
    }
  }

  pos = bestpos;
  LFS lfs(gfs);
  lfs.bind(*ep);
  return lfs.globalIndex(best_edge);
}

template<typename RF, unsigned dimR>
class Excitation {
  const double t0;
  const double sigma;
  const Dune::FieldVector<RF, dimR> amp;

  static Dune::FieldVector<RF, dimR> default_amp() {
    Dune::FieldVector<RF, dimR> amp;
    std::copy(&excitation_amp[0], &excitation_amp[dimR], amp.begin());
    return amp;
  }

public:
  Excitation(const Dune::FieldVector<RF, dimR>& amp_ = default_amp(),
             double sigma_ = excitation_sigma,
             double t0_ = excitation_t0)
    : t0(t0_), sigma(sigma_), amp(amp_)
  {}

  Dune::FieldVector<RF, dimR>
  operator()(double t) const {
    t -= t0;
    Dune::FieldVector<RF, dimR> result = amp;

    result *= (-t/std::sqrt(2*pi)/sigma/sigma/sigma*
               std::exp(-t*t/2/sigma/sigma));
    return result;
  }
};

//===============================================================
// Problem setup and solution 
//===============================================================

// generate a P1 function and output it
template<typename GV, typename FEM, typename CON, typename ReferenceFactory,
         typename FProbe, typename CProbe> 
void electrodynamic (const GV& gv, const FEM& fem, unsigned integrationOrder,
                       const ReferenceFactory &referenceFactory,
                       double Delta_t, unsigned steps,
                       const std::string filename, FProbe &fprobe, CProbe &cprobe)
{
  // constants and types
  typedef typename GV::ctype DF;
  static const unsigned dimD = GV::dimensionworld;
  typedef Dune::FieldVector<DF, dimD> D;
  typedef typename FEM::Traits::LocalFiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType RangeField;

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,
    Dune::PDELab::ISTLVectorBackend<1> > GFS; 
  GFS gfs(gv,fem);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<RangeField>::Type C;
  C cg;
  cg.clear();
  typedef B<GV> BType;
  BType b(gv);
  Dune::PDELab::constraints(b,gfs,cg);

  // make coefficent Vector and initialize it from a function
  typedef typename GFS::template VectorContainer<RangeField>::Type V;
  typedef Dune::PDELab::DiscreteGridFunctionGlobal<GFS,V> DGF;
  typedef Dune::PDELab::DiscreteGridFunctionGlobalCurl<GFS,V> CDGF;

  //initial solution for time-step -1
  Dune::SmartPointer<V> xprev(new V(gfs));
  *xprev = 0.0;
  Dune::PDELab::interpolateGlobal(*referenceFactory.function(gv, -Delta_t),gfs,*xprev);
  // actually, it is an error if any constrained dof is != 0, but set it here nevertheless
  Dune::PDELab::set_constrained_dofs(cg,0.0,*xprev);

  //initial solution for time-step 0
  Dune::SmartPointer<V> xcur(new V(gfs));
  *xcur = 0.0;
  Dune::PDELab::interpolateGlobal(*referenceFactory.function(gv, 0),gfs,*xcur);
  // actually, it is an error if any constrained dof is != 0, but set it here nevertheless
  Dune::PDELab::set_constrained_dofs(cg,0.0,*xcur);

  Dune::SmartPointer<V> xnext(0);

  // we're using dirichlet 0 everywhere, simply leave everything as 0
  V affineShift(gfs);
  affineShift = 0.0;
  //Dune::PDELab::interpolateGlobal(g,gfs,affineShift);
  //Dune::PDELab::set_nonconstrained_dofs(cg,0.0,affineShift);

  // make grid function operator
  typedef ConstFunc<GV,RangeField> MuType;
  MuType mu(gv,mu0);
  typedef ConstFunc<GV,RangeField> EpsType;
  EpsType eps(gv,eps0);
  typedef ConstFunc<GV,RangeField,GV::dimensionworld> DtJType;
  DtJType dtJ(gv,typename DtJType::Traits::RangeType(0));
  typedef Dune::PDELab::Electrodynamic<EpsType,MuType,DtJType,V> LOP;
  D dtJPos(0.5);
  LOP lop(eps, mu, dtJ, edge_dof_near_pos(gfs, dtJPos), Delta_t,
          integrationOrder);
  std::cout << "Setting dipole at " << fmt(dtJPos) << std::endl;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,
    LOP,C,C,Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > GOS;
  GOS gos(gfs,cg,gfs,cg,lop);
  Excitation<DF, dimD> excitation;

  // represent operator as a matrix
  typedef typename GOS::template MatrixContainer<RangeField>::Type M;
  M m(gos);
  m = 0;
  // just use some values here
  lop.setEprev(*xprev);
  lop.setEcur(*xcur);
  lop.setDtJPCur(excitation(0));
  gos.jacobian(affineShift,m);
//   Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);
  Dune::Richardson<V,V> prec(1.0);
//   Dune::SeqILU0<M,V,V> prec(m,1.0);
  Dune::MatrixAdapter<M,V,V> op(m);
  Dune::CGSolver<V> solver(op,prec,1E-10,1000,0);

//   typedef Dune::SuperLU<typename M::BaseT> Solver;
//   Solver solver(m);

  if(use_interpolated_reference) {
    Dune::SmartPointer<V> ref(new V(gfs));
    *ref = 0.0;
    Dune::PDELab::interpolateGlobal(*referenceFactory.function(gv, -Delta_t),
                                    gfs,*ref);
    fprobe.measureExact(DGF(gfs, *xprev),
                        DGF(gfs, *ref),
                        -Delta_t);
  } else
    fprobe.measureExact(DGF(gfs, *xprev),
                        *referenceFactory.function(gv, -Delta_t),
                        -Delta_t);
  cprobe.measure(CDGF(gfs, *xprev), -Delta_t);
//  std::cout << "u[-1]\n" << *xprev << std::endl;
  if(use_interpolated_reference) {
    Dune::SmartPointer<V> ref(new V(gfs));
    *ref = 0.0;
    Dune::PDELab::interpolateGlobal(*referenceFactory.function(gv, 0),
                                    gfs,*ref);
    fprobe.measureExact(DGF(gfs, *xcur),
                        DGF(gfs, *ref),
                        0);
  } else
    fprobe.measureExact(DGF(gfs, *xcur),
                        *referenceFactory.function(gv, 0),
                        0);
  cprobe.measure(CDGF(gfs, *xcur), 0);
//   std::cout << "u[0]\n" << *xcur << std::endl;

  std::cout << "Number of steps " << steps << std::endl;
  for(unsigned step = 1; step <= steps; ++step) {
    std::cout << "Doing step " << step << std::endl;


    lop.setEprev(*xprev);
    lop.setEcur(*xcur);
    lop.setDtJPCur(excitation(Delta_t*(step-1)));

    // evaluate residual w.r.t initial guess
    V rhs(gfs);
    rhs = 0.0;
    gos.residual(affineShift,rhs);
    rhs *= -1;

//     std::cout << "RHS\n" << rhs << std::endl;

    Dune::InverseOperatorResult stat;

    xnext = new V(gfs);
    *xnext = *xcur; // copy values to have a better start value for the solver
    *xnext -= affineShift;
    solver.apply(*xnext,rhs,stat);
    std::cout << "Solver results:" << std::endl
              << "  iterations: " << stat.iterations << std::endl
              << "  reduction:  " << stat.reduction << std::endl
              << "  converged:  " << stat.converged << std::endl
              << "  conv_rate:  " << stat.conv_rate << std::endl
              << "  elapsed:    " << stat.elapsed << std::endl;
    if(!stat.converged)
      DUNE_THROW(Dune::Exception, "Solver did not converge");

    *xnext += affineShift;

    if(use_interpolated_reference) {
      Dune::SmartPointer<V> ref(new V(gfs));
      *ref = 0.0;
      Dune::PDELab::interpolateGlobal(*referenceFactory.function
                                      (gv, Delta_t*step),
                                      gfs, *ref);
      fprobe.measureExact(DGF(gfs, *xnext),
                          DGF(gfs, *ref),
                          Delta_t*step);
    } else
      fprobe.measureExact(DGF(gfs, *xnext),
                          *referenceFactory.function(gv, Delta_t*step),
                          Delta_t*step);
    cprobe.measure(CDGF(gfs, *xnext), Delta_t*step);
//     std::cout << "u[" << step << "]\n" << *xnext << std::endl;

    xprev = xcur;
    xcur = xnext;
  }
  
  if(use_interpolated_reference) {
    Dune::SmartPointer<V> ref(new V(gfs));
    *ref = 0.0;
    Dune::PDELab::interpolateGlobal(*referenceFactory.function
                                    (gv, Delta_t*steps),
                                    gfs, *ref);
    fprobe.measureFinalExact(DGF(gfs, *xcur),
                             DGF(gfs, *ref),
                             Delta_t*steps);
  } else
    fprobe.measureFinalExact(DGF(gfs, *xcur),
                             *referenceFactory.function(gv, Delta_t*steps),
                             Delta_t*steps);
  cprobe.measureFinal(CDGF(gfs, *xcur), Delta_t*steps);
}

template<typename P>
struct GridPtrTraits;
template<typename G>
struct GridPtrTraits<Dune::SmartPointer<G> > {
  typedef G Grid;
};
template<typename G>
struct GridPtrTraits<Dune::GridPtr<G> > {
  typedef G Grid;
};

template<typename GV, typename RF, unsigned dimDomain = GV::dimension>
struct FEMTraits;
template<typename GV, typename RF>
struct FEMTraits<GV, RF, 2>
{
  typedef Dune::PDELab::EdgeS02DLocalFiniteElementMap<GV, RF> FEM;
};
template<typename GV, typename RF>
struct FEMTraits<GV, RF, 3>
{
  typedef Dune::PDELab::EdgeS03DLocalFiniteElementMap<GV, RF> FEM;
};

template<typename GV, typename LPF, typename CLPF, typename ELPF>
void testLevel(const GV &gv, unsigned level, const std::string &prefix,
               LPF &lpf, CLPF &clpf, ELPF &elpf, double &error, double &mean_h)
{
  typedef typename FEMTraits<GV, double>::FEM FEM;
  std::ostringstream levelprefix;
  levelprefix << prefix << ".level" << level;

  typedef typename LPF::template Traits<GV>::Probe Probe;
  Dune::SmartPointer<Probe> probe(lpf.getProbe(gv, level));

  typedef typename CLPF::template Traits<GV>::Probe CProbe;
  Dune::SmartPointer<CProbe> cprobe(clpf.getProbe(gv, level));

  typedef typename ELPF::template Traits<GV>::Probe EProbe;
  Dune::SmartPointer<EProbe> eprobe(elpf.getProbe(gv, level));

  typedef Dune::PDELab::ProbePair<Probe, EProbe> FProbe;
  Dune::SmartPointer<FProbe> fprobe(new FProbe(probe, eprobe));

  std::cout << "electrodynamic level " << level << std::endl;
  // time step
  double Delta_t = stepadjust*smallestEdge(gv)/std::sqrt(double(GV::dimension))/c0;
  unsigned steps = duration/Delta_t;
  if(Delta_t*steps < duration) {
    //adjust such that duration is hit exactly
    ++steps;
    Delta_t = duration/steps;
  }
  electrodynamic
    <GV,FEM,Dune::PDELab::OverlappingConformingDirichletConstraints,ResonatorSolutionFactory<GV,double>,FProbe,CProbe>
    (gv, FEM(gv), quadrature_order,
     ResonatorSolutionFactory<GV,double>(),
     Delta_t, steps, levelprefix.str(), *fprobe, *cprobe);
  mean_h = std::pow(1/double(gv.size(0)), 1/double(GV::dimension));
  error = eprobe->get_error();
}
template<typename GV, typename LPF, typename CLPF, typename ELPF>
void testLevel(const GV &gv, unsigned level, const std::string &prefix,
               LPF &lpf, CLPF &clpf, ELPF &elpf)
{
  double error, mean_h;
  testLevel(gv, level, prefix, lpf, clpf, elpf, error, mean_h);
}


template<typename Grid, typename GPF, typename CGPF, typename EGPF>
void test(Grid &grid, int &result, GPF &gpf, CGPF &cgpf, EGPF &egpf,
          double conv_limit, std::string name = "")
{
  typedef typename Grid::LeafGridView GV;

  if(name == "") name = grid.name();

  typedef typename GPF::template Traits<Grid>::LevelProbeFactory LPF;
  Dune::SmartPointer<LPF> lpf(gpf.levelProbeFactory(grid, name));

  typedef typename CGPF::template Traits<Grid>::LevelProbeFactory CLPF;
  Dune::SmartPointer<CLPF> clpf(cgpf.levelProbeFactory(grid, name));

  typedef typename EGPF::template Traits<Grid>::LevelProbeFactory ELPF;
  Dune::SmartPointer<ELPF> elpf(egpf.levelProbeFactory(grid, name));

  std::cout << std::endl
            << "Testing Electrodynamic problem with EdgeS03D and " << name << std::endl;

  std::string filename = "electrodynamic-" + name;

  unsigned level = 0;
  double error0, mean_h0;
  testLevel(grid.leafView(), level, filename, *lpf, *clpf, *elpf, error0, mean_h0);

  while(true) {
    ++level;
    grid.globalRefine(1);
    if(unsigned(grid.leafView().size(0)) >= maxelements) break;

    if(do_all_levels)
      testLevel(grid.leafView(), level, filename, *lpf, *clpf, *elpf);
  }
  
  double errorf, mean_hf;
  testLevel(grid.leafView(), level, filename, *lpf, *clpf, *elpf, errorf, mean_hf);

  double total_convergence = std::log(errorf/error0)/std::log(mean_hf/mean_h0);
  std::cout << "electrodynamic total convergence: "
            << std::scientific << total_convergence << std::endl;

  if(result != 1)
    result = 0;

  if(total_convergence < conv_limit) {
    std::cout << "Error: electrodynamic total convergence < " << conv_limit << std::endl;
    result = 1;
  }
}


//===============================================================
// Main program with grid setup
//===============================================================

template<typename GPF, typename CGPF, typename EGPF>
void testAll(int &result, GPF &gpf, CGPF & cgpf, EGPF &egpf) {
#ifdef HAVE_ALBERTA
#if (ALBERTA_DIM != 3)
#error ALBERTA_DIM is not set to 3 -- please check the Makefile.am
#endif
//   test(*UnitTetrahedronMaker         <Dune::AlbertaGrid<3, 3>    >::create(),
//        result, graph, conv_limit,    "alberta-tetrahedron");
  // test(*KuhnTriangulatedUnitCubeMaker<Dune::AlbertaGrid<3, 3>    >::create(),
  //      result, gpf, cgpf, egpf, .7*conv_limit, "alberta-triangulated-cube-6");
//   {
//     Dune::GridPtr<Dune::AlbertaGrid<3, 3> > gridptr("grids/brick.dgf");
//     test(*gridptr,
//          result, graph, conv_limit,    "alu-triangulated-brick-6");
//   }
#endif

#ifdef HAVE_ALUGRID
//   test(*UnitTetrahedronMaker         <Dune::ALUSimplexGrid<3, 3> >::create(),
//        result, graph, conv_limit,    "alu-tetrahedron");
  // test(*KuhnTriangulatedUnitCubeMaker<Dune::ALUSimplexGrid<3, 3> >::create(),
  //      result, gpf, cgpf, egpf, conv_limit,    "alu-triangulated-cube-6");
#endif // HAVE_ALUGRID

#ifdef HAVE_UG
//   test(*UnitTetrahedronMaker         <Dune::UGGrid<3>            >::create(),
//        result, graph, conv_limit,    "ug-tetrahedron");
  test(*KuhnTriangulatedUnitCubeMaker<Dune::UGGrid<3>            >::create(),
       result, gpf, cgpf, egpf, conv_limit,    "ug-triangulated-cube-6");
  // test(*TriangulatedUnitSquareMaker<Dune::UGGrid<2>            >::create(),
  //      result, gpf, cgpf, egpf, conv_limit,    "ug-triangulated-square");
#endif // HAVE_UG
}

int main(int argc, char** argv)
{
  Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);
  if(mpiHelper.rank() == 0)
    std::cout << "Number of processes: " << mpiHelper.size() << std::endl;
  std::cout << "Rank: " << mpiHelper.rank() << std::endl;
//   if(mpiHelper.rank() == 1) {
//     int i = 0;
//     std::cout << "PID " << getpid() << " ready for attach\n", getpid();
//     while (0 == i)
//       sleep(5);
//   }

  try{
    // 77 is special and means "test was skipped".  Return that if non of the
    // supported grids were available
    int result = 77;

    // Initialize probe_location_fv from probe_location
    typedef Dune::PDELab::GnuplotGridProbeFactory<double> PointPF;
    Dune::SmartPointer<PointPF>
      pointPF(new PointPF("electrodynamic-probe", probe_location));

    typedef VTKGridProbeFactory VTKOutput;
    Dune::SmartPointer<VTKOutput>
      vtkOutput(new VTKOutput("electrodynamic"));

    typedef GlobalErrorGridProbeFactory GlobalError;
    Dune::SmartPointer<GlobalError>
      globalError(new GlobalError(quadrature_order, "electrodynamic-globalerror", "electrodynamic-globalevolution"));

    typedef L2ErrorGridProbeFactory L2Error;
    Dune::SmartPointer<L2Error>
      l2Error(new L2Error("electrodynamic-l2error", quadrature_order));

    typedef L2ErrorEvolutionGridProbeFactory L2Evolution;
    Dune::SmartPointer<L2Evolution>
      l2Evolution(new L2Evolution("electrodynamic-l2evolution", quadrature_order));

    typedef ElectricEnergyGridProbeFactory EnergyEvolution;
    Dune::SmartPointer<EnergyEvolution>
      energyEvolution(new EnergyEvolution("electrodynamic-electricenergy-evolution", quadrature_order));

    typedef DivergenceGridProbeFactory DivergenceEvolution;
    Dune::SmartPointer<DivergenceEvolution>
      divergenceEvolution
      (new DivergenceEvolution
       ("electrodynamic-divergence-evolution", quadrature_order));

    typedef L2NormGridProbeFactory L2NormEvolution;
    Dune::SmartPointer<L2NormEvolution>
      curlEvolution
      (new L2NormEvolution
       ("electrodynamic-curl-evolution", "Curl", quadrature_order));

    testAll(result,
            *Dune::PDELab::makeGridProbeFactoryList(pointPF,
                                                    vtkOutput,
                                                    globalError,
                                                    l2Evolution,
                                                    energyEvolution,
                                                    divergenceEvolution),
            *Dune::PDELab::makeGridProbeFactoryList(curlEvolution),
            *l2Error);

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
