// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/* Solve problems with the PDELab solvers (StationaryLinearProblemSolver,
 * NewtonMethod) and linear solver backends using the GridOperator on
 * dune-functions bases. The exact solutions are contained in the discrete
 * spaces, so the discrete solutions must coincide with their interpolants:
 *
 * - linear Poisson problem, Q2, ConvectionDiffusionFEM, ISTLBackend_SEQ_UMFPack
 * - nonlinear reaction-diffusion problem, Q2, numerical jacobian, ISTLSolver
 * - Navier-Stokes, Taylor-Hood, TaylorHoodNavierStokes, Dirichlet
 *   constraints on the velocity subspace and one pressure DOF
 */

#include "config.h"

#include <cmath>
#include <iostream>
#include <string>

#include <dune/common/filledarray.hh>
#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/backend/istl/istlsolver.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/numericaljacobian.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/stokesparameter.hh>
#include <dune/pdelab/localoperator/taylorhoodnavierstokes.hh>
#include <dune/pdelab/solver/newton.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;

// initial guess: exact values on constrained DOFs, zero elsewhere
template<class V, class CC>
V initialGuess(const V& exact, const CC& constraints)
{
  V x(exact.gridFunctionSpace(), 0.0);
  for (const auto& [index, row] : constraints)
    x[index] = exact[index];
  return x;
}

template<class V>
int check(const V& x, const V& exact, double tolerance, const std::string& what)
{
  auto error = Dune::PDELab::Backend::native(x);
  error -= Dune::PDELab::Backend::native(exact);
  std::cout << what << ": max nodal error " << error.infinity_norm() << std::endl;
  if (error.infinity_norm() > tolerance)
  {
    std::cerr << "FAILED: " << what << std::endl;
    return 1;
  }
  return 0;
}


// -Delta u = f, u = x^2 + y^2 on the whole boundary
template<typename GV, typename RF>
class QuadraticPoissonProblem
  : public Dune::PDELab::ConvectionDiffusionModelProblem<GV,RF>
{
public:
  using Traits = Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF>;

  typename Traits::RangeFieldType
  f (const typename Traits::ElementType&, const typename Traits::DomainType&) const
  {
    return -4.0;
  }

  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    auto xg = e.geometry().global(x);
    return xg.two_norm2();
  }
};

template<class GV>
int testLinearPoisson(const GV& gv)
{
  using namespace Dune;
  std::cout << "=== StationaryLinearProblemSolver, Poisson, Q2" << std::endl;

  auto basis = Functions::BasisFactory::makeBasis(gv, Functions::BasisFactory::lagrange<2>());
  using Basis = decltype(basis);
  using Param = QuadraticPoissonProblem<GV,double>;
  Param param;
  using LOP = PDELab::ConvectionDiffusionFEM<Param,
    PDELab::FiniteElementMapFor<typename Basis::LocalView::Tree::FiniteElement>>;
  LOP lop(param);

  using CC = PDELab::BasisConstraintsContainer<Basis>;
  CC cc;
  PDELab::addBoundaryDirichletConstraints(basis, cc);

  using GO = PDELab::GridOperator<Basis,Basis,LOP,MBE,double,double,double,CC,CC>;
  GO go(basis, cc, basis, cc, lop, MBE(25));

  using V = typename GO::Traits::Domain;
  V exact(basis, 0.0);
  Functions::interpolate(basis, PDELab::Backend::native(exact), [](const auto& x) { return x.two_norm2(); });
  auto x = initialGuess(exact, cc);

#if HAVE_SUITESPARSE_UMFPACK
  PDELab::ISTLBackend_SEQ_UMFPack ls(0);
#else
  PDELab::ISTLBackend_SEQ_BCGS_ILU0 ls(5000, 0);
#endif
  PDELab::StationaryLinearProblemSolver<GO,decltype(ls),V> slp(go, ls, x, 1e-12, 1e-99, 0);
  slp.apply();

  return check(x, exact, 1e-10, "linear Poisson");
}


// -Delta u + u^3 = f, u = x^2 + y^2 on the whole boundary; the jacobian is computed numerically
class NonlinearReactionDiffusion
  : public Dune::PDELab::NumericalJacobianVolume<NonlinearReactionDiffusion>
  , public Dune::PDELab::FullVolumePattern
  , public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  enum { doPatternVolume = true };
  enum { doAlphaVolume = true };
  enum { isLinear = false };

  template<class EG, class LFSU, class X, class LFSV, class R>
  void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    using FE = typename LFSU::Traits::FiniteElementType;
    using RangeType = typename FE::Traits::LocalBasisType::Traits::RangeType;
    using JacobianType = typename FE::Traits::LocalBasisType::Traits::JacobianType;
    constexpr int dim = EG::Entity::dimension;

    const auto& localBasis = lfsu.finiteElement().localBasis();
    const auto geo = eg.geometry();
    std::vector<RangeType> phi(lfsu.size());
    std::vector<JacobianType> referenceGradients(lfsu.size());
    std::vector<Dune::FieldVector<double,dim>> gradients(lfsu.size());

    const auto rule = Dune::QuadratureRules<double,dim>::rule(geo.type(), 2*localBasis.order()+2);
    for (const auto& ip : rule)
    {
      localBasis.evaluateFunction(ip.position(), phi);
      localBasis.evaluateJacobian(ip.position(), referenceGradients);
      const auto jit = geo.jacobianInverseTransposed(ip.position());
      double u = 0.0;
      Dune::FieldVector<double,dim> gradu(0.0);
      for (std::size_t i = 0; i < lfsu.size(); ++i)
      {
        jit.mv(referenceGradients[i][0], gradients[i]);
        u += x(lfsu,i) * phi[i];
        gradu.axpy(x(lfsu,i), gradients[i]);
      }
      const auto xg = geo.global(ip.position());
      const double exact = xg.two_norm2();
      const double f = -4.0 + exact*exact*exact;
      const double factor = ip.weight() * geo.integrationElement(ip.position());
      for (std::size_t i = 0; i < lfsv.size(); ++i)
        r.accumulate(lfsv, i, (gradu * gradients[i] + (u*u*u - f) * phi[i]) * factor);
    }
  }
};

template<class GV>
int testNonlinear(const GV& gv)
{
  using namespace Dune;
  std::cout << "=== NewtonMethod, nonlinear reaction-diffusion, Q2" << std::endl;

  auto basis = Functions::BasisFactory::makeBasis(gv, Functions::BasisFactory::lagrange<2>());
  using Basis = decltype(basis);
  NonlinearReactionDiffusion lop;

  using CC = PDELab::BasisConstraintsContainer<Basis>;
  CC cc;
  PDELab::addBoundaryDirichletConstraints(basis, cc);

  using GO = PDELab::GridOperator<Basis,Basis,NonlinearReactionDiffusion,MBE,double,double,double,CC,CC>;
  GO go(basis, cc, basis, cc, lop, MBE(25));

  using V = typename GO::Traits::Domain;
  using M = typename GO::Traits::Jacobian;
  V exact(basis, 0.0);
  Functions::interpolate(basis, PDELab::Backend::native(exact), [](const auto& x) { return x.two_norm2(); });
  auto x = initialGuess(exact, cc);

  // a dune-istl inverse operator through the generic ISTLSolver adapter
#if HAVE_SUITESPARSE_UMFPACK
  auto ls = PDELab::makeUMFPackSolver<PDELab::Backend::Native<M>,PDELab::Backend::Native<V>>();
#else
  auto ls = PDELab::makeBiCGSTABILU0Solver<PDELab::Backend::Native<M>,PDELab::Backend::Native<V>>();
#endif
  PDELab::NewtonMethod<GO,decltype(ls)> newton(go, ls);
  newton.setReduction(1e-12);
  newton.setVerbosityLevel(1);
  newton.apply(x);

  int result = check(x, exact, 1e-8, "nonlinear reaction-diffusion");
  if (newton.result().iterations < 2)
  {
    std::cerr << "FAILED: Newton needed only " << newton.result().iterations << " iterations" << std::endl;
    result = 1;
  }
  return result;
}


// Navier-Stokes with exact solution u = (x^2, -2xy), p = x + y
template<typename GV>
class NavierStokesProblem
{
public:
  static const bool assemble_navier = true;
  static const bool assemble_full_tensor = false;

  using Traits = Dune::PDELab::NavierStokesParameterTraits<GV,double>;

  static constexpr double viscosity = 1.5;
  static constexpr double density = 0.8;

  void setTime(double) {}

  // f = -mu Delta u + grad p + rho (u.grad) u
  template<typename Element>
  typename Traits::VelocityRange
  f(const Element& e, const typename Traits::Domain& x) const
  {
    auto xg = e.geometry().global(x);
    return {-2.0*viscosity + 1.0 + density*2.0*xg[0]*xg[0]*xg[0],
            1.0 + density*2.0*xg[0]*xg[0]*xg[1]};
  }

  template<typename EG>
  typename Traits::RangeField
  g2(const EG&, const typename Traits::Domain&) const
  {
    return 0.0;
  }

  template<typename EG>
  typename Traits::RangeField
  mu(const EG&, const typename Traits::Domain&) const { return viscosity; }

  template<typename EG>
  typename Traits::RangeField
  rho(const EG&, const typename Traits::Domain&) const { return density; }

  template<typename IG>
  typename Traits::BoundaryCondition::Type
  bctype(const IG&, const typename Traits::IntersectionDomain&) const
  {
    return Traits::BoundaryCondition::VelocityDirichlet;
  }

  template<typename IG>
  typename Traits::VelocityRange
  j(const IG&, const typename Traits::IntersectionDomain&, const typename Traits::Domain&) const
  {
    return typename Traits::VelocityRange(0.0);
  }
};

template<class GV>
int testNavierStokes(const GV& gv)
{
  using namespace Dune;
  using namespace Dune::Indices;
  constexpr int dim = GV::dimension;
  std::cout << "=== NewtonMethod, Navier-Stokes, Taylor-Hood" << std::endl;

  using namespace Functions::BasisFactory;
  auto basis = makeBasis(gv,
    composite(
      power<dim>(lagrange<2>(), flatInterleaved()),
      lagrange<1>(),
      flatLexicographic()));
  using Basis = decltype(basis);

  using Param = NavierStokesProblem<GV>;
  Param param;
  using LOP = PDELab::TaylorHoodNavierStokes<Param>;
  LOP lop(param);

  auto velocityBasis = Functions::subspaceBasis(basis, _0);
  auto pressureBasis = Functions::subspaceBasis(basis, _1);

  // velocity Dirichlet values on the whole boundary; the pressure is fixed in one DOF
  using CC = PDELab::BasisConstraintsContainer<Basis>;
  CC cc;
  PDELab::addBoundaryDirichletConstraints(velocityBasis, cc);
  {
    auto localView = pressureBasis.localView();
    localView.bind(*gv.template begin<0>());
    cc[localView.index(localView.tree().localIndex(0))];
  }

  using GO = PDELab::GridOperator<Basis,Basis,LOP,MBE,double,double,double,CC,CC>;
  GO go(basis, cc, basis, cc, lop, MBE(50));

  using V = typename GO::Traits::Domain;
  V exact(basis, 0.0);
  auto& nativeExact = PDELab::Backend::native(exact);
  Functions::interpolate(Functions::subspaceBasis(basis, _0, 0), nativeExact, [](const auto& x) { return x[0]*x[0]; });
  Functions::interpolate(Functions::subspaceBasis(basis, _0, 1), nativeExact, [](const auto& x) { return -2.0*x[0]*x[1]; });
  Functions::interpolate(pressureBasis, nativeExact, [](const auto& x) { return x[0] + x[1]; });
  auto x = initialGuess(exact, cc);

#if HAVE_SUITESPARSE_UMFPACK
  PDELab::ISTLBackend_SEQ_UMFPack ls(0);
#else
  PDELab::ISTLBackend_SEQ_BCGS_ILU0 ls(5000, 0);
#endif
  PDELab::NewtonMethod<GO,decltype(ls)> newton(go, ls);
  newton.setReduction(1e-12);
  newton.setVerbosityLevel(1);
  newton.apply(x);

  int result = check(x, exact, 1e-8, "Navier-Stokes Taylor-Hood");
  if (newton.result().iterations < 2)
  {
    std::cerr << "FAILED: Newton needed only " << newton.result().iterations << " iterations" << std::endl;
    result = 1;
  }
  return result;
}


int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  using Grid = Dune::YaspGrid<2>;
  Grid grid({1.0,1.0}, Dune::filledArray<2,int>(8));
  auto gv = grid.leafGridView();

  int result = 0;
  result |= testLinearPoisson(gv);
  result |= testNonlinear(gv);
  result |= testNavierStokes(gv);
  return result;
}
