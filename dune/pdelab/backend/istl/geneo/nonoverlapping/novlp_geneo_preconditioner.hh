#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_NOVLP_GENEO_PRECONDITIONER_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_NOVLP_GENEO_PRECONDITIONER_HH

#include <dune/pdelab/backend/istl/geneo/nonoverlapping/geneobasis.hh>

// #include <dune/composites/Driver/Solvers/zembasis_virtual_overlap.hh>

namespace Dune {
  namespace PDELab {

    /**
     * @brief Black-box GenEO preconditioner working on nonoverlapping matrices
     */
    template <typename GO, typename Matrix, typename Vector>
    class NonoverlappingGenEOPreconditioner : public Dune::Preconditioner<Vector,Vector> {

      using GFS = typename GO::Traits::TrialGridFunctionSpace;
      using GV = typename GFS::Traits::GridView;


    public:

      // define the category
      virtual Dune::SolverCategory::Category category() const override
      {
        return prec->category();
      }


      /*!
      *          \brief Prepare the preconditioner.
      *
      *          \copydoc Preconditioner::pre(X&,Y&)
      */
      virtual void pre (Vector& x, Vector& b) override
      {
        prec->pre(x, b);
      }

      /*!
      *          \brief Apply the precondioner.
      *
      *          \copydoc Preconditioner::apply(X&,const Y&)
      */
      virtual void apply (Vector& v, const Vector& d) override
      {
        prec->apply(v, d);
      }
      /*!
      *          \brief Clean up.
      *
      *          \copydoc Preconditioner::post(X&)
      */
      virtual void post (Vector& x) override
      {
        prec->post(x);
      }

      NonoverlappingGenEOPreconditioner(const GO& go, const typename GO::Jacobian& A, int algebraic_overlap, int avg_nonzeros, const double eigenvalue_threshold, int& nev,
                          int nev_arpack = -1, const double shift = 0.001, int verbose = 0, std::string preconditionerType = "fullyAdditive")
      {
        using Dune::PDELab::Backend::native;

        const GV& gv = go.trialGridFunctionSpace().gridView();

        Dune::NonoverlappingOverlapAdapter<GV, Vector, Matrix> adapter(gv, native(A), avg_nonzeros, algebraic_overlap);

        auto geneo_matrices = setupGenEOMatrices(go, adapter, A);
        std::shared_ptr<Matrix> A_extended = std::get<0>(geneo_matrices);
        std::shared_ptr<Matrix> A_overlap_extended = std::get<1>(geneo_matrices);
        std::shared_ptr<Vector> part_unity = std::get<2>(geneo_matrices);
        std::vector<int> IntBndDofs = std::get<3>(geneo_matrices);

        auto subdomainbasis = std::make_shared<Dune::PDELab::NonoverlappingGenEOBasis<GO, Matrix, Vector>>(adapter, A_extended, A_overlap_extended, part_unity, eigenvalue_threshold, nev, nev_arpack, shift, true);

        /* Test preconditionner using zeros energy modes to build the coarse space*/
        // auto subdomainbasis = std::make_shared<Dune::PDELab::ZEMBasis_virtual_overlap<GO,Vector,Matrix,3,3> >(go, adapter, part_unity);

        // std::cout << adapter.gridView().comm().rank() << " before the coarse space construction." << std::endl;

        // auto coarse_space = std::make_shared<Dune::PDELab::NonoverlappingSubdomainProjectedCoarseSpaceHeldByRank0<GV, Matrix, Vector>>(adapter, gv, *A_extended, subdomainbasis, verbose);
        auto coarse_space = std::make_shared<Dune::PDELab::NonoverlappingSubdomainProjectedCoarseSpace<GV, Matrix, Vector>>(adapter, gv, *A_extended, subdomainbasis, verbose);
        if(verbose)
          std::cout << adapter.gridView().comm().rank() << " have built the coarse space." << std::endl;

        if (preconditionerType == "restrictedHybrid") {
          prec = std::make_shared<Dune::PDELab::ISTL::RestrictedHybridTwoLevelSchwarz<GV, Matrix, Vector>>(adapter, A_extended, avg_nonzeros, *part_unity, IntBndDofs, coarse_space, true, verbose);
        } else {
          prec = std::make_shared<Dune::PDELab::ISTL::NonoverlappingTwoLevelOverlappingAdditiveSchwarz<GV, Matrix, Vector>>(adapter, A_extended, *part_unity, coarse_space, true, verbose);
        }
        if(verbose)
          std::cout << adapter.gridView().comm().rank() << " prec initialized." << std::endl;

      }

    private:

      std::shared_ptr<Dune::Preconditioner<Vector,Vector>> prec = nullptr;

    };

  }
}

#endif
