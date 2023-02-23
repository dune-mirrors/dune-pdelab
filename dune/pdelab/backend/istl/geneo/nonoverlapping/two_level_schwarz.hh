#ifndef NONOVERLAPPING_TWO_LEVEL_SCHWARZ_HH
#define NONOVERLAPPING_TWO_LEVEL_SCHWARZ_HH

#if HAVE_SUITESPARSE_UMFPACK

#include <dune/common/timer.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/umfpack.hh>
#include <dune/pdelab/gridfunctionspace/dunefunctionsgridfunctionspace.hh>

#include <dune/pdelab/backend/istl/geneo/coarsespace.hh>
#include "overlaptools.hh"

#include <dune/common/parallel/communicator.hh>

namespace Dune {



  namespace PDELab {
    namespace ISTL {

      /*!
      * \brief Two level overlapping Schwarz preconditioner with arbitrary coarse space.
      */
      template<typename GridView, typename Matrix, typename Vector>
      class NonoverlappingTwoLevelOverlappingAdditiveSchwarz
       : public Dune::Preconditioner<Vector,Vector>
      {

        template<typename V>
        struct AddGatherScatter
        {
          static typename V::value_type gather(const V& a, int i)
          {
            return a[i]; // I am sending my value
          }
          static void scatter (V& a, typename V::value_type v, int i)
          {
            a[i]+=v; // add what I receive to my value
          }
        };

      public:
        // define the category
        virtual Dune::SolverCategory::Category category() const
        {
          return Dune::SolverCategory::nonoverlapping;
        }

        NonoverlappingTwoLevelOverlappingAdditiveSchwarz (const NonoverlappingOverlapAdapter<GridView, Vector, Matrix>& adapter, std::shared_ptr<Matrix> A, const Vector& part_unity, std::shared_ptr<CoarseSpace<Vector> > coarse_space, bool coarse_space_active = true, int verbosity = 0)
        : verbosity_(verbosity),
          coarse_space_active_(coarse_space_active),
          adapter_(adapter),
          A_(A),
          coarse_space_(coarse_space),
          coarse_solver_ (*coarse_space_->get_coarse_system()),
          coarse_defect_(coarse_space_->basis_size(), coarse_space_->basis_size()),
          prolongated_(adapter_.getExtendedSize())
        {
          const int block_size = Vector::block_type::dimension;

          // Apply Dirichlet conditions to matrix on processor boundaries, inferred from partition of unity
          for (auto rIt=A_->begin(); rIt!=A_->end(); ++rIt){
            for(int block_i = 0; block_i < block_size; block_i++){
              if (part_unity[rIt.index()][block_i] == .0) {
                for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
                {
                  for(int block_j = 0; block_j < block_size; block_j++){
                    (*cIt)[block_i][block_j] = (rIt.index() == cIt.index() && block_i == block_j) ? 1.0 : 0.0;
                  }
                }
              }
            }
          }

          solverf_ = std::make_shared<Dune::UMFPack<Matrix>>(*A_,false);
        }

        /*!
          \brief Prepare the preconditioner.

          \copydoc Preconditioner::pre(X&,Y&)
        */
        virtual void pre (Vector& x, Vector& b)
        { }

        /*!
          \brief Apply the precondioner.

          \copydoc Preconditioner::apply(X&,const Y&)
        */
        virtual void apply (Vector& v, const Vector& d)
        {

          using Attribute = EPISAttribute;
          Dune::AllSet<Attribute> allAttribute;
          auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
          allinterface->build(*adapter_.getRemoteIndices(),allAttribute,allAttribute); // all to all communication

          // build up buffered communicator allowing communication over a dof vector
          auto communicator = std::shared_ptr<Dune::BufferedCommunicator>(new Dune::BufferedCommunicator());
          communicator->build<Vector>(*allinterface);

          if (verbosity_ > 2) Dune::printvector(std::cout, d, "defect (local)", "", 1, 10, 17);

          // first the subdomain solves
          Vector b(adapter_.getExtendedSize()), correction(adapter_.getExtendedSize()); // need copy, since solver overwrites right hand side
          adapter_.extendVector(d, b);
          if (verbosity_ > 2) Dune::printvector(std::cout, b, "defect (extended)", "", 1, 10, 17);
          communicator->forward<AddGatherScatter<Vector>>(b,b); // make function known in other subdomains
          if (verbosity_ > 2) Dune::printvector(std::cout, b, "defect (distributed)", "", 1, 10, 17);


          const int block_size = Vector::block_type::dimension;
          Vector b_cpy(b);// FIXME: Avoid this
          for (auto rIt=A_->begin(); rIt!=A_->end(); ++rIt) {
              for(int block_i = 0; block_i < block_size; block_i++){
                  bool isDirichlet = true;
                  for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
                  {
                      for(int block_j = 0; block_j < block_size; block_j++){
                          if ((rIt.index() != cIt.index() || block_i!=block_j) && (*cIt)[block_i][block_j] != 0.0) {
                              isDirichlet = false;
                              break;
                          }
                      }
                      if(!isDirichlet) break;
                  }
                  if (isDirichlet) {
                      b_cpy[rIt.index()] = .0;
                      b[rIt.index()] = .0;
                  }
              }
          }
          if (verbosity_ > 2) Dune::printvector(std::cout, b_cpy, "defect (Dirichlet applied) ", "", 1, 10, 17);

          Dune::InverseOperatorResult result;
          solverf_->apply(correction,b_cpy,result);

          if (verbosity_ > 2) Dune::printvector(std::cout, correction, "correction (1lvl) ", "", 1, 10, 17);

          if (!coarse_space_active_) {

            communicator->forward<AddGatherScatter<Vector>>(correction,correction); // make function known in other subdomains

            adapter_.restrictVector(correction, v);

          } else {

            Dune::Timer timer_coarse_solve;

            coarse_space_->restrict (b, coarse_defect_);

            if (verbosity_ > 2) Dune::printvector(std::cout, coarse_defect_, "coarse_defect_ ", "", 1, 10, 17);

            typename CoarseSpace<Vector>::COARSE_V v0(coarse_space_->basis_size());
            coarse_solver_.apply(v0, coarse_defect_, result);
            if (verbosity_ > 2) Dune::printvector(std::cout, v0, "v0 ", "");
            coarse_space_->prolongate(v0, prolongated_);

            if (verbosity_ > 2) Dune::printvector(std::cout, prolongated_, "prolongated_ ", "", 1, 10, 17);

            correction += prolongated_;
            if (verbosity_ > 2) Dune::printvector(std::cout, correction, "correction ", "", 1, 10, 17);

            communicator->forward<AddGatherScatter<Vector>>(correction,correction); // make function known in other subdomains
            if (verbosity_ > 2) Dune::printvector(std::cout, correction, "correction (sum) ", "", 1, 10, 17);

            adapter_.restrictVector(correction, v);
            if (verbosity_ > 2) Dune::printvector(std::cout, v, "correction (restricted) ", "", 1, 10, 17);

            coarse_time_ += timer_coarse_solve.elapsed();
          }
          apply_calls_++;
        }

        /*!
          \brief Clean up.

          \copydoc Preconditioner::post(X&)
        */
        virtual void post (Vector& x) {
          if (verbosity_ > 0) std::cout << "Coarse time CT=" << coarse_time_ << std::endl;
          if (verbosity_ > 0) std::cout << "Coarse time per apply CTA=" << coarse_time_ / apply_calls_ << std::endl;
        }

      private:
        int verbosity_;
        bool coarse_space_active_;

        NonoverlappingOverlapAdapter<GridView, Vector, Matrix> adapter_;

        std::shared_ptr<Matrix> A_ = nullptr;
        std::shared_ptr<Dune::UMFPack<Matrix>> solverf_ = nullptr;

        double coarse_time_ = 0.0;
        int apply_calls_ = 0;

        std::shared_ptr<CoarseSpace<Vector> > coarse_space_;
        Dune::UMFPack<typename CoarseSpace<Vector>::COARSE_M> coarse_solver_;

        typename CoarseSpace<Vector>::COARSE_V coarse_defect_;
        Vector prolongated_;

      };


      /*!
      * \brief Restricted, hybrid (i.e. parallel in local solves, coarse solve apllied muliplicatively),
       overlapping Two level Schwarz preconditioner with arbitrary coarse space for use on a nonoverlapping grid.
      */
      template<typename GridView, typename Matrix, typename Vector>
      class RestrictedHybridTwoLevelSchwarz
       : public Dune::Preconditioner<Vector,Vector>
      {

        template<typename V>
        struct AddGatherScatter
        {
          static typename V::value_type gather(const V& a, int i)
          {
            return a[i]; // I am sending my value
          }
          static void scatter (V& a, typename V::value_type v, int i)
          {
            a[i]+=v; // add what I receive to my value
          }
        };

      public:
        // define the category
        virtual Dune::SolverCategory::Category category() const
        {
          return Dune::SolverCategory::nonoverlapping;
        }

        RestrictedHybridTwoLevelSchwarz (const NonoverlappingOverlapAdapter<GridView, Vector, Matrix>& adapter, std::shared_ptr<Matrix> A, int avg_nonzeros, const Vector& part_unity,std::vector<int>& intBndDofs, std::shared_ptr<CoarseSpace<Vector> > coarse_space, bool coarse_space_active = true, int verbosity = 0)
        : verbosity_(verbosity),
          coarse_space_active_(coarse_space_active),
          adapter_(adapter),
          A_(A),
          coarse_space_(coarse_space),
          coarse_solver_ (*coarse_space_->get_coarse_system()),
          coarse_defect_(coarse_space_->basis_size(), coarse_space_->basis_size()),
          prolongated_(adapter_.getExtendedSize()),
          part_unity_(part_unity),
          intBndDofs_(intBndDofs)
        {
          using Dune::PDELab::Backend::native;
          const int block_size = Vector::block_type::dimension;

          // make a copy of A because we need it with and without Dirichlet bc on processor boundaries
          A_Neumann_ = std::shared_ptr<Matrix>(new Matrix(adapter_.getExtendedSize(),adapter_.getExtendedSize(),avg_nonzeros,0.1,Matrix::implicit));
          for (auto rIt=A_->begin(); rIt!=A_->end(); ++rIt) {// loop over entries of A
            for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt) {
              auto i = rIt.index();
              auto j = cIt.index();
              A_Neumann_->entry(i,j) = *cIt;
            }
          }
          A_Neumann_->compress();
          if (verbosity_> 2) {std::cout << "size of A_Neumann_ = " << A_Neumann_->N() << std::endl;}

          // Apply Dirichlet conditions to matrix A_ on processor boundaries
          int j = 0;
          for (auto rIt=A_->begin(); rIt!=A_->end(); ++rIt){
            if (rIt.index() == intBndDofs_[j]){
              if (verbosity_> 2) {std::cout << "j = " << j << ". intBndDof index = " << intBndDofs_[j] << std::endl;}
              j++;
              for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt){
                for(int block_i = 0; block_i < block_size; block_i++){
                  for(int block_j = 0; block_j < block_size; block_j++){
                    (*cIt)[block_i][block_j] = (rIt.index() == cIt.index() && block_i == block_j) ? 1.0 : 0.0;
                  }
                }
              }
            }
          }
          if (verbosity_> 2) {std::cout << "j = " << j << ". intBndDof.size() = " << intBndDofs_.size() << ". (Should be equal.)" << std::endl;}

          solverf_ = std::make_shared<Dune::UMFPack<Matrix>>(*A_,false);

          // // is A_Neumann_ actually different from A_?
          // std::cout << "Frobenius norm of A_ = " << native(*A_).frobenius_norm() << std::endl;
          // std::cout << "Frobenius norm of A_Neumann_ = " << native(*A_Neumann_).frobenius_norm() << std::endl;
        }

        /*!
          \brief Prepare the preconditioner.

          \copydoc Preconditioner::pre(X&,Y&)
        */
        virtual void pre (Vector& x, Vector& b)
        { }

        /*!
          \brief Apply the precondioner.

          \copydoc Preconditioner::apply(X&,const Y&)
        */
        virtual void apply (Vector& v, const Vector& d_const)
        {
          using Dune::PDELab::Backend::native;
          using Attribute = EPISAttribute;
          Dune::AllSet<Attribute> allAttribute;
          auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
          allinterface->build(*adapter_.getRemoteIndices(),allAttribute,allAttribute); // all to all communication

          // build up buffered communicator allowing communication over a dof vector
          auto communicator = std::shared_ptr<Dune::BufferedCommunicator>(new Dune::BufferedCommunicator());
          communicator->build<Vector>(*allinterface);

          Dune::Timer timer_local_solve;
          // ------------------------------------- subdomain solves --------------------------------------
          // first, prepare rhs. Since d_const is constant, we have to make a copy.
          Vector d(d_const);
          if (verbosity_ > 2) Dune::printvector(std::cout, d, "defect (local)", "", 1, 10, 17);

          // extend nonoverlapping vector d to overlapping vector b.
          Vector b(adapter_.getExtendedSize()), correction(adapter_.getExtendedSize());
          adapter_.extendVector(d, b);
          // std::cout << "Maximum entry in b = " << b.infinity_norm() << std::endl;

          if (verbosity_ > 2) Dune::printvector(std::cout, b, "defect (extended)", "", 1, 10, 17);
          // make function known in other subdomains. Here it is important that we masked the dof on the nonoverlapping processor boundaries earlier
          communicator->forward<AddGatherScatter<Vector>>(b,b);
          if (verbosity_ > 2) Dune::printvector(std::cout, b, "defect (distributed)", "", 1, 10, 17);

          // write zeros into b at processor boudnaries, for the Dirichlet BCs there.
          const int block_size = Vector::block_type::dimension;
          Vector b_cpy(b);// need copy, since solver overwrites right hand side
          for (int j = 0; j<intBndDofs_.size() ; ++j) {
            for (int i = 0; i < block_size; ++i){
              b_cpy[intBndDofs_[j]][i] = .0;
            }
          }

          if (verbosity_ > 2) Dune::printvector(std::cout, b_cpy, "defect (Dirichlet applied) ", "", 1, 10, 17);

          // solve the local system
          Dune::InverseOperatorResult result;
          solverf_->apply(correction,b_cpy,result);

          // multiply local correction with partition of unity
          for (int i = 0; i < correction.size(); ++i) {
            for (int j = 0; j<correction[i].size(); ++j) {
              correction[i][j]*=part_unity_[i][j];
            }
          }
          local_time_ += timer_local_solve.elapsed();

          communicator->forward<AddGatherScatter<Vector>>(correction,correction); // add all local corrections up

          if (verbosity_ > 2) Dune::printvector(std::cout, correction, "correction (1lvl) ", "", 1, 10, 17);

          // -------------------------------------- Coarse Solve --------------------------------------
          if (!coarse_space_active_) {
            // if coarse space is inactive, do nothing.
            // just bring correction into nonoverlapping format by restricting
            adapter_.restrictVector(correction, v);
          } else {
            // account for local corrections in rhs b
            native(*A_Neumann_).mmv(native(correction), native(b)); // b = b - A_ * correction

            Dune::Timer timer_coarse_solve;

            // project rhs b into coarse space
            coarse_space_->restrict(b, coarse_defect_);
            if (verbosity_ > 2) Dune::printvector(std::cout, coarse_defect_, "coarse_defect_ ", "", 1, 10, 17);

            // solve the coarse system
            typename CoarseSpace<Vector>::COARSE_V v0(coarse_space_->basis_size());
            coarse_solver_.apply(v0, coarse_defect_, result);
            if (verbosity_ > 2) Dune::printvector(std::cout, v0, "v0 ", "");
            // iclusion of coarse vector into fine space
            coarse_space_->prolongate(v0, prolongated_);
            // Note that the function prolongate only uses the coefficient of each processes' own eigenfunctions.
            // Thus we have to communicate and add the results of the neighbors.
            communicator->forward<AddGatherScatter<Vector>>(prolongated_,prolongated_);

            if (verbosity_ > 2) Dune::printvector(std::cout, prolongated_, "prolongated_ ", "", 1, 10, 17);

            // add coarse correction to local corrections
            correction += prolongated_;
            adapter_.restrictVector(correction, v);

            if (verbosity_ > 2) Dune::printvector(std::cout, v, "correction (restricted) ", "", 1, 10, 17);

            coarse_time_ += timer_coarse_solve.elapsed();
          }
          apply_calls_++;
        }

        /*!
          \brief Clean up.

          \copydoc Preconditioner::post(X&)
        */
        virtual void post (Vector& x) {
          if (verbosity_ > 0) std::cout << "Local time =" << local_time_ << std::endl;
          if (verbosity_ > 0) std::cout << "Local time per apply =" << local_time_ / apply_calls_ << std::endl;
          if (verbosity_ > 0) std::cout << "Coarse time CT=" << coarse_time_ << std::endl;
          if (verbosity_ > 0) std::cout << "Coarse time per apply CTA=" << coarse_time_ / apply_calls_ << std::endl;
        }

      private:
        int verbosity_;
        bool coarse_space_active_;

        NonoverlappingOverlapAdapter<GridView, Vector, Matrix> adapter_;

        std::shared_ptr<Matrix> A_ = nullptr; // stiffness matrix. Will be overwritten in constructor with homogeneous Dirichlet BC on processor boundaries.
        std::shared_ptr<Matrix> A_Neumann_ = nullptr; // explicit BCRSMatrix(); // stiffness matrix. Will not be overwritten in constructor.
        std::shared_ptr<Dune::UMFPack<Matrix>> solverf_ = nullptr;
        Vector part_unity_ = nullptr;
        // Vector mask_ = nullptr;
        std::vector<int> intBndDofs_;

        double coarse_time_ = 0.0;
        double local_time_ = 0.0;
        int apply_calls_ = 0;

        std::shared_ptr<CoarseSpace<Vector> > coarse_space_;
        Dune::UMFPack<typename CoarseSpace<Vector>::COARSE_M> coarse_solver_;

        typename CoarseSpace<Vector>::COARSE_V coarse_defect_;
        Vector prolongated_;

      };
    }
  }
}
#endif

#endif
