// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_ISTLSOLVER_HH
#define DUNE_PDELAB_BACKEND_ISTL_ISTLSOLVER_HH

#include <functional>
#include <memory>

#include <dune/common/ftraits.hh>

#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solver.hh>
#include <dune/istl/solvers.hh>
#if HAVE_SUITESPARSE_UMFPACK
#include <dune/istl/umfpack.hh>
#endif

#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/solver.hh>

/** \file
 * \brief PDELab linear solver interface for arbitrary dune-istl inverse operators.
 *
 * NewtonMethod and StationaryLinearProblemSolver expect a linear solver
 * with apply(A,z,r,reduction), norm(r) and result(). ISTLSolver provides
 * this for any dune-istl InverseOperator created by a user-given factory,
 * working on PDELab vector/matrix wrappers (via Backend::native) as well as
 * on plain ISTL containers.
 *
 * \warning This is experimental and sequential only.
 */

namespace Dune::PDELab {

  /** \brief PDELab linear solver backend wrapping a dune-istl InverseOperator
   *
   * \tparam M native ISTL matrix type
   * \tparam V native ISTL vector type
   */
  template<class M, class V>
  class ISTLSolver
    : public Dune::PDELab::LinearResultStorage
  {
  public:
    using Matrix = M;
    using Vector = V;
    using InverseOperator = Dune::InverseOperator<Vector,Vector>;
    using Real = typename Dune::FieldTraits<typename Dune::FieldTraits<Vector>::field_type>::real_type;

    //! creates the solver for a given matrix; the matrix outlives the solver
    using Factory = std::function<std::shared_ptr<InverseOperator>(const Matrix&)>;

    explicit ISTLSolver(Factory factory)
      : factory_(std::move(factory))
    {}

    //! solve A z = r up to the given relative reduction; r is not modified
    template<class MatrixOrWrapper, class VectorOrWrapper, class RhsOrWrapper>
    void apply(MatrixOrWrapper& a, VectorOrWrapper& z, const RhsOrWrapper& r, Real reduction)
    {
      const Matrix& matrix = Backend::native(a);
      if (not solver_ or not reuse_ or &matrix != matrix_)
      {
        solver_ = factory_(matrix);
        matrix_ = &matrix;
      }
      Vector rhs = Backend::native(r);
      Dune::InverseOperatorResult statistics;
      solver_->apply(Backend::native(z), rhs, reduction, statistics);
      res.converged = statistics.converged;
      res.iterations = statistics.iterations;
      res.elapsed = statistics.elapsed;
      res.reduction = statistics.reduction;
      res.conv_rate = statistics.conv_rate;
    }

    //! euclidean norm of a residual
    template<class VectorOrWrapper>
    Real norm(const VectorOrWrapper& v) const
    {
      return Backend::native(v).two_norm();
    }

    //! reuse the solver (e.g. a factorization) for the same matrix; used by NewtonMethod
    void setReuse(bool reuse) { reuse_ = reuse; }
    bool getReuse() const { return reuse_; }

  private:
    Factory factory_;
    std::shared_ptr<InverseOperator> solver_;
    const Matrix* matrix_ = nullptr;
    bool reuse_ = false;
  };

#if HAVE_SUITESPARSE_UMFPACK
  //! ISTLSolver using the UMFPack direct solver
  template<class Matrix, class Vector>
  auto makeUMFPackSolver(int verbosity = 0)
  {
    return ISTLSolver<Matrix,Vector>([verbosity](const Matrix& a) {
      return std::make_shared<Dune::UMFPack<Matrix>>(a, verbosity);
    });
  }
#endif

  //! ISTLSolver using BiCGSTAB preconditioned with ILU(0)
  template<class Matrix, class Vector>
  auto makeBiCGSTABILU0Solver(int maxIterations = 1000, int verbosity = 0)
  {
    return ISTLSolver<Matrix,Vector>([=](const Matrix& a) {
      struct Solver : Dune::InverseOperator<Vector,Vector>
      {
        Solver(const Matrix& a, int maxIterations, int verbosity)
          : op(a), ilu(a, 1.0), bicgstab(op, ilu, 1e-10, maxIterations, verbosity)
        {}
        void apply(Vector& x, Vector& b, Dune::InverseOperatorResult& r) override { bicgstab.apply(x, b, r); }
        void apply(Vector& x, Vector& b, double reduction, Dune::InverseOperatorResult& r) override { bicgstab.apply(x, b, reduction, r); }
        Dune::SolverCategory::Category category() const override { return Dune::SolverCategory::sequential; }
        Dune::MatrixAdapter<Matrix,Vector,Vector> op;
        Dune::SeqILU<Matrix,Vector,Vector> ilu;
        Dune::BiCGSTABSolver<Vector> bicgstab;
      };
      return std::make_shared<Solver>(a, maxIterations, verbosity);
    });
  }

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_BACKEND_ISTL_ISTLSOLVER_HH
