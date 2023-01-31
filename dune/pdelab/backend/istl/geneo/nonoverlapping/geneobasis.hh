#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_NonoverlappingGenEOBasis_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_NonoverlappingGenEOBasis_HH

#include <algorithm>
#include <functional>

#include <dune/pdelab/backend/istl/geneo/subdomainbasis.hh>
#include <dune/pdelab/backend/istl/geneo/arpackpp_geneo.hh>

#if HAVE_ARPACKPP

namespace Dune {
  namespace PDELab {

    /*!
     * \brief Implementation of the GenEO coarse basis
     * See Spillane et al., 2014, 'Abstract robust coarse spaces for systems of PDEs via generalized eigenproblems in the overlaps'
     * This coarse space is based on generalized eigenpoblems defined on the full stiffness matrix of a subdomain and one assembled only
     * on the area where this subdomain overlaps with others.
     */
    template<class GO, class Matrix, class Vector>
    class NonoverlappingGenEOBasis : public SubdomainBasis<Vector>
    {
      using GFS = typename GO::Traits::TrialGridFunctionSpace;
      using GridView = typename GFS::Traits::GridView;
      using GV = GridView;
      using V = typename GO::Domain;
      using M = typename GO::Jacobian;

    public:


      /**
       * @brief Constructor
       *
       * @param adapter Adapter for extending matrices by virtual overlap
       * @param A_extended Neumann matrix extended by virtual overlap
       * @param A_ovlp_extended Neumann matrix restricted to overlap region and extended by virtual overlap
       * @param part_unity Partition of unity extended by virtual overlap
       * @param eigenvalue_threshold Threshold up to which eigenvalues should be included in the GenEO space; no threshold is applied if value negative
       * @param nev Number of eigenvectors to be included in coarse space; returns chosen value in case of thresholding
       * @param nev_arpack Number of eigenvectors to be computed by ARPACK, possibly a larger value for stability
       * @param shift Shift value for shift invert mode solution of eigenproblem
       * @param add_part_unity Whether to add the partition of unity to the coarse space
       */
      NonoverlappingGenEOBasis(NonoverlappingOverlapAdapter<GridView, Vector, Matrix>& adapter, std::shared_ptr<Matrix> A_extended, std::shared_ptr<Matrix> A_ovlp_extended, std::shared_ptr<Vector> part_unity, const double eigenvalue_threshold,
                    int& nev, int nev_arpack = -1, const double shift = 0.001, const bool add_part_unity = false, const int verbose = 0) {

        if (nev_arpack == -1)
          nev_arpack = std::max(nev, 2);
        if (nev_arpack < nev)
          DUNE_THROW(Dune::Exception,"nev_arpack is less then nev!");

        // X * A_0 * X
        Matrix ovlp_mat(*A_ovlp_extended);
        using Dune::PDELab::Backend::native;
        const int block_size = Vector::block_type::dimension;
        for (auto row_iter = native(ovlp_mat).begin(); row_iter != native(ovlp_mat).end(); row_iter++) {
          for (auto col_iter = row_iter->begin(); col_iter != row_iter->end(); col_iter++) {
            for (int i = 0; i < block_size; i++) {
              for (int j = 0; j < block_size; j++) {
                (*col_iter)[i][j] *= native(*part_unity)[row_iter.index()][i] * native(*part_unity)[col_iter.index()][j];
              }
            }
          }
        }

        // Setup Arpack for solving generalized eigenproblem
        if(verbose)
          std::cout << "ARPACK setup...";
        ArpackGeneo::ArPackPlusPlus_Algorithms<Matrix, Vector> arpack(*A_extended, 1000, verbose);
        if(verbose)
          std::cout << " done" << std::endl;
        double eps = 0.0;

        std::vector<double> eigenvalues(nev_arpack,0.0);
        std::vector<Vector> eigenvectors(nev_arpack,Vector(adapter.getExtendedSize()));


        if(verbose)
          std::cout << "ARPACK solve...";
        arpack.computeGenNonSymMinMagnitude(ovlp_mat, eps, eigenvectors, eigenvalues, shift);
        if(verbose)
          std::cout << " done" << std::endl;

        // Dune::BlockVector<Dune::FieldVector<double, 1>> export_eigenvalues(nev_arpack);
        // // Dune::BlockVector<Dune::FieldVector<double, 1>> error(nev_arpack);

        // for(int g =0; g<nev_arpack; g++){
        //   // Vector left(A_extended->N()), right(A_extended->N());
        //   // A_extended->mv(eigenvectors[g],left);
        //   // ovlp_mat.mv(eigenvectors[g],right);
        //   // right *= eigenvalues[g];

        //   // double L2 = 0.0;
        //   // for (int i=0; i< A_extended->N(); i++) {
        //   //   double local_dist=0.0;
        //   //   local_dist += (right[i][0] - left[i][0])*(right[i][0] - left[i][0]);
        //   //   local_dist += (right[i][1] - left[i][1])*(right[i][1] - left[i][1]);
        //   //   local_dist += (right[i][2] - left[i][2])*(right[i][2] - left[i][2]);
        //   //   L2 += local_dist;
        //   // }
        //   // L2 = sqrt(L2);
        //   // std::cout << "Subdomain : " << adapter.gridView().comm().rank() << " ==> L2 error norm EV " << g << " : " << L2 << std::endl;
        //   // std::cout << "EValue: " << eigenvalues[g] << std::endl;
        //   export_eigenvalues[g] = eigenvalues[g];
        //   // error[g] = L2;
        // }

        // std::string filename_eigenvalues = "Offline/" + std::to_string(adapter.gridView().comm().rank()) + "_" + std::to_string(nev_arpack) + "_GenEO-eigenvalues.mm";
        // Dune::storeMatrixMarket(export_eigenvalues, filename_eigenvalues, 15);
        // std::string filename_error = "Offline/" + std::to_string(adapter.gridView().comm().rank()) + "_" + std::to_string(nev_arpack) + "_GenEO-error.mm";
        // Dune::storeMatrixMarket(error, filename_error, 15);


        // Count eigenvectors below threshold
        int cnt = -1;
        if (eigenvalue_threshold >= 0) {
          for (int i = 0; i < nev; i++) {
            if (eigenvalues[i] > eigenvalue_threshold) {
              cnt = i;
              break;
            }
          }
          if (verbose > 0)
            std::cout << "Process " << adapter.gridView().comm().rank() << " picked " << cnt << " eigenvectors" << std::endl;
          if (cnt == -1)
            cnt = nev;
            // DUNE_THROW(Dune::Exception,"No eigenvalue above threshold - not enough eigenvalues computed!");
        } else {
          cnt = nev;
        }

        // std::cout << adapter.gridView().comm().rank() << " ici." << std::endl;

        // Write results to basis
        this->local_basis.resize(cnt);
        for (int base_id = 0; base_id < cnt; base_id++) {
          this->local_basis[base_id] = std::make_shared<Vector>(adapter.getExtendedSize());
          *this->local_basis[base_id] = eigenvectors[base_id];
          for (int i = 0; i < part_unity->N(); i++) {
            for (int j = 0; j < block_size; j++) {
              (*(this->local_basis[base_id]))[i][j] *= (*part_unity)[i][j];
            }
          }
        }

        // std::cout << adapter.gridView().comm().rank() << " wrote results to basis." << std::endl;

        // Normalize basis vectors
        for (auto& v : this->local_basis) {
          *v *= 1.0 / v->two_norm2();
        }

        // Optionally add partition of unity to eigenvectors
        // Only if there is no near-zero eigenvalue (that usually already corresponds to a partition of unity!)
        if (add_part_unity && eigenvalues[0] > 1E-10) {
          this->local_basis.insert (this->local_basis.begin(), std::make_shared<Vector>(*part_unity));
          this->local_basis.pop_back();
        }

        if(verbose)
          std::cout << adapter.gridView().comm().rank() << " have added PoU in the basis." << std::endl;
      }
    };


  }
}

#endif

#endif
