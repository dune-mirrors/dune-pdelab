#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_NOVLP_GENEO_MATRIX_SETUP_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_NOVLP_GENEO_MATRIX_SETUP_HH

#include "overlaptools.hh"
#include "overlaptools_serendipity.hh"

namespace Dune {
  namespace PDELab {

    template <typename Vector, typename GV>
    class EntitySetExcluder {
    public:

      typedef typename GV::template Codim<0>::Entity Entity;
      typedef typename GV::template Codim<1>::Entity Entity1;
      typedef typename GV::template Codim<2>::Entity Entity2;

      bool includeEntity(const Entity1& entity) const {
        return true;
      }
      bool includeEntity(const Entity2& entity) const {
        return true;
      }

      virtual bool includeEntity(const Entity& entity) const {
        return true;
      }

    };

    /**
     * @brief Excludes elements where a given partition of unity is nonzero
     */
    template <typename Vector, typename GV, typename GFS>
    class EntitySetPartUnityExcluder : public EntitySetExcluder<Vector, GV> {
    public:

      typedef typename GV::template Codim<0>::Entity Entity;

      EntitySetPartUnityExcluder(const GFS& gfs, std::shared_ptr<Vector> partUnity) : lfs(gfs), partUnity_(partUnity) {}

      bool includeEntity(const Entity& entity) const override {

        if (entity.partitionType() == Dune::PartitionType::GhostEntity)
          return false;

        typedef Dune::PDELab::LFSIndexCache<LFS,Dune::PDELab::EmptyTransformation> LFSCache;
        LFSCache lfs_cache(lfs);
        lfs.bind( entity );
        lfs_cache.update();

        for (std::size_t i = 0; i < lfs_cache.size(); i++)
        {
          if ((*partUnity_)[lfs_cache.containerIndex(i).back()][0] > 0.0 &&
            (*partUnity_)[lfs_cache.containerIndex(i).back()][0] < 1.0) {
            return true;
            }
        }
        return false;
      }

    private:
      typedef Dune::PDELab::LocalFunctionSpace<GFS, Dune::PDELab::TrialSpaceTag> LFS;
      mutable LFS lfs;
      std::shared_ptr<Vector> partUnity_ = nullptr;
    };

    /**
     * @brief Excludes ghost elements
     */
    template <typename Vector, typename GV>
    class EntitySetGhostExcluder : public EntitySetExcluder<Vector, GV> {
    public:

      typedef typename GV::template Codim<0>::Entity Entity;

      bool includeEntity(const Entity& entity) const override {
        return entity.partitionType() != Dune::PartitionType::GhostEntity;
      }
    };


    /**
     * @brief Create GenEO matrices with virtual overlap from nonoverlapping Neumann matrix
     *
     * @param go Grid operator representing weak form of PDE
     * @param adapter Virtual overlap adapter
     * @param A Nonoverlapping matrix assembled with Neumann boundaries
     * @return A tuple containing the GenEO matrices extended by virtual overlap and the corresponding partition of unity
     */
    template<class GO, class Matrix, class Vector>
    std::tuple<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>, std::shared_ptr<Vector>, std::vector<int>> setupGenEOMatrices(const GO& go, NonoverlappingOverlapAdapter<typename GO::Traits::TrialGridFunctionSpace::Traits::GridView, Vector, Matrix>& adapter, const typename GO::Jacobian& A, bool need_enforced_dirichlet_on_PoU=1, int PoU_restriction=0) {

      using Dune::PDELab::Backend::native;
      Dune::Timer timer_detailed(true);

      using GFS = typename GO::Traits::TrialGridFunctionSpace;
      using GV = typename GFS::Traits::GridView;

      const GFS& gfs = go.trialGridFunctionSpace();
      const GV& gv = gfs.gridView();

      std::shared_ptr<Matrix> A_extended = adapter.extendMatrix(native(A));

      std::shared_ptr<Vector> part_unity = Dune::makePartitionOfUnity<GV, Matrix, Vector>(adapter, *A_extended); // The one used for A_extended setup


      auto PoU_and_IntBndDofs = Dune::makePartitionOfUnityRestricted<GV, Matrix, Vector>(adapter, *A_extended, PoU_restriction); // The one used later one, potentially restricted on boundaries to have an oversampled subdomain
      std::shared_ptr<Vector> part_unity_R = part_unity;
      // std::shared_ptr<Vector> part_unity_R = std::get<0>(PoU_and_IntBndDofs);
      std::vector<int> IntBndDofs = std::get<1>(PoU_and_IntBndDofs);


      using Attribute = Dune::EPISAttribute;
      Dune::AllSet<Attribute> allAttribute;
      auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
      allinterface->build(*adapter.getRemoteIndices(),allAttribute,allAttribute); // all to all communication
      auto communicator = std::shared_ptr<DuneWithRank::BufferedCommunicator>(new DuneWithRank::BufferedCommunicator());
      communicator->build<Vector>(*allinterface);

      Dune::PDELab::MultiVectorBundle<GV, Vector, Matrix> remotePartUnities(adapter);
      remotePartUnities.localVector_ = part_unity;
      communicator->forward<Dune::PDELab::MultiGatherScatter<Dune::PDELab::MultiVectorBundle<GV, Vector, Matrix>>>(remotePartUnities,remotePartUnities); // make function known in other subdomains


      // Assemble fine grid matrix defined only on overlap region
      auto part_unity_restricted = std::make_shared<Vector>(native(A).N());
      adapter.restrictVector(*part_unity, *part_unity_restricted);

      auto es_pou_excluder = std::make_shared<EntitySetPartUnityExcluder<Vector, GV, GFS>> (gfs, part_unity_restricted);
      gfs.entitySet().setExcluder(es_pou_excluder);


      using V = typename GO::Domain;
      using M = typename GO::Jacobian;

      V x(gfs,0.0); // NOTE: We assume linear problems, so simply set x to zero here!
      M A_ovlp(go);
      go.jacobian(x,A_ovlp);

      M newmat(go);
      // Provide neighbors with matrices assembled exclusively on respective overlap area
      std::pair<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>> extended_matrices = adapter.lambdaMultiExtendMatrix(native(A), native(A_ovlp), [&](int i){
        std::shared_ptr<Vector> neighbor_part_unity = remotePartUnities.getVectorForRank(i);

        adapter.restrictVector(*neighbor_part_unity, *part_unity_restricted);

        std::shared_ptr<EntitySetExcluder<Vector, GV>> es_local_pou_excluder = std::make_shared<EntitySetPartUnityExcluder<Vector, GV, GFS>> (gfs, part_unity_restricted);
        gfs.entitySet().setExcluder(es_local_pou_excluder);

        for (auto rIt=native(newmat).begin(); rIt!=native(newmat).end(); ++rIt)
          for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt) {
            *cIt = 0.0;
          }

        go.jacobian(x,newmat);

        return stackobject_to_shared_ptr(native(newmat));
      });



      // double time_A_extended_and_A_ovlp = timer_detailed.elapsed(); timer_detailed.reset();
      // std::cout << "A_extended_and_A_ovlp creation: " << time_A_extended_and_A_ovlp << std::endl;

      if (need_enforced_dirichlet_on_PoU==true){
        // Enforce problem's Dirichlet condition on PoU
        const int block_size = Vector::block_type::dimension;
        for (auto rIt=A_extended->begin(); rIt!=A_extended->end(); ++rIt) {
          for(int block_i = 0; block_i < block_size; block_i++){ //loop over block
            bool isDirichlet = true;
            for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
            {
              for(int block_j = 0; block_j<block_size; block_j++){
                if ((rIt.index() != cIt.index() || block_i!=block_j) && (*cIt)[block_i][block_j] != 0.0) {
                  isDirichlet = false;
                  break;
                }
              }
              if(!isDirichlet) break;
            }
            if (isDirichlet) {
              (*part_unity)[rIt.index()][block_i] = .0;
            }
          }
        }
      }

      using ESExcluder = Dune::PDELab::EntitySetExcluder<Vector, GV>;
      gfs.entitySet().setExcluder(std::make_shared<ESExcluder>());

      return std::make_tuple(extended_matrices.first, extended_matrices.second, part_unity_R, IntBndDofs);
    }


    /**
     * @brief Create GenEO matrices with virtual overlap from nonoverlapping Neumann matrix
     *
     * @param go Grid operator representing weak form of PDE
     * @param adapter Virtual overlap adapter
     * @param A Nonoverlapping matrix assembled with Neumann boundaries
     * @return A tuple containing the GenEO matrices extended by virtual overlap and the corresponding partition of unity
     */
    template<class GO, class Matrix, class Vector>
    std::tuple<std::shared_ptr<Matrix>, std::shared_ptr<Vector>, std::vector<int>> setupGenEOMatrices2(const GO& go, NonoverlappingOverlapAdapter<typename GO::Traits::TrialGridFunctionSpace::Traits::GridView, Vector, Matrix>& adapter, const typename GO::Jacobian& A, bool need_enforced_dirichlet_on_PoU=1, int PoU_restriction=0) {

      using Dune::PDELab::Backend::native;
      Dune::Timer timer_detailed(true);

      using GFS = typename GO::Traits::TrialGridFunctionSpace;
      using GV = typename GFS::Traits::GridView;

      const GFS& gfs = go.trialGridFunctionSpace();
      const GV& gv = gfs.gridView();

      std::shared_ptr<Matrix> A_extended = adapter.extendMatrix(native(A));

      std::shared_ptr<Vector> part_unity = Dune::makePartitionOfUnity<GV, Matrix, Vector>(adapter, *A_extended); // The one used for A_extended setup


      auto PoU_and_IntBndDofs = Dune::makePartitionOfUnityRestricted<GV, Matrix, Vector>(adapter, *A_extended, PoU_restriction); // The one used later one, potentially restricted on boundaries to have an oversampled subdomain
      // std::shared_ptr<Vector> part_unity_R = part_unity;
      std::shared_ptr<Vector> part_unity_R = std::get<0>(PoU_and_IntBndDofs);
      std::vector<int> IntBndDofs = std::get<1>(PoU_and_IntBndDofs);


      using Attribute = Dune::EPISAttribute;
      Dune::AllSet<Attribute> allAttribute;
      auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
      allinterface->build(*adapter.getRemoteIndices(),allAttribute,allAttribute); // all to all communication
      auto communicator = std::shared_ptr<DuneWithRank::BufferedCommunicator>(new DuneWithRank::BufferedCommunicator());
      communicator->build<Vector>(*allinterface);

      Dune::PDELab::MultiVectorBundle<GV, Vector, Matrix> remotePartUnities(adapter);
      remotePartUnities.localVector_ = part_unity;
      communicator->forward<Dune::PDELab::MultiGatherScatter<Dune::PDELab::MultiVectorBundle<GV, Vector, Matrix>>>(remotePartUnities,remotePartUnities); // make function known in other subdomains


      // Assemble fine grid matrix defined only on overlap region
      auto part_unity_restricted = std::make_shared<Vector>(native(A).N());
      adapter.restrictVector(*part_unity, *part_unity_restricted);

      auto es_pou_excluder = std::make_shared<EntitySetPartUnityExcluder<Vector, GV, GFS>> (gfs, part_unity_restricted);
      gfs.entitySet().setExcluder(es_pou_excluder);


      using V = typename GO::Domain;
      using M = typename GO::Jacobian;

      V x(gfs,0.0); // NOTE: We assume linear problems, so simply set x to zero here!
      M A_ovlp(go);
      go.jacobian(x,A_ovlp);

      M newmat(go);
      // Provide neighbors with matrices assembled exclusively on respective overlap area
      std::pair<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>> extended_matrices = adapter.lambdaMultiExtendMatrix(native(A), native(A_ovlp), [&](int i){
        std::shared_ptr<Vector> neighbor_part_unity = remotePartUnities.getVectorForRank(i);

        adapter.restrictVector(*neighbor_part_unity, *part_unity_restricted);

        std::shared_ptr<EntitySetExcluder<Vector, GV>> es_local_pou_excluder = std::make_shared<EntitySetPartUnityExcluder<Vector, GV, GFS>> (gfs, part_unity_restricted);
        gfs.entitySet().setExcluder(es_local_pou_excluder);

        for (auto rIt=native(newmat).begin(); rIt!=native(newmat).end(); ++rIt)
          for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt) {
            *cIt = 0.0;
          }

        go.jacobian(x,newmat);

        return stackobject_to_shared_ptr(native(newmat));
      });

      if (need_enforced_dirichlet_on_PoU==true){
        // Enforce problem's Dirichlet condition on PoU
        const int block_size = Vector::block_type::dimension;
        for (auto rIt=A_extended->begin(); rIt!=A_extended->end(); ++rIt) {
          for(int block_i = 0; block_i < block_size; block_i++){ //loop over block
            bool isDirichlet = true;
            for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
            {
              for(int block_j = 0; block_j<block_size; block_j++){
                if ((rIt.index() != cIt.index() || block_i!=block_j) && (*cIt)[block_i][block_j] != 0.0) {
                  isDirichlet = false;
                  break;
                }
              }
              if(!isDirichlet) break;
            }
            if (isDirichlet) {
              (*part_unity)[rIt.index()][block_i] = .0;
            }
          }
        }
      }

      using ESExcluder = Dune::PDELab::EntitySetExcluder<Vector, GV>;
      gfs.entitySet().setExcluder(std::make_shared<ESExcluder>());

      // return std::make_tuple(extended_matrices.first, part_unity, IntBndDofs);
      return std::make_tuple(extended_matrices.first, part_unity_R, IntBndDofs);
    }



    /**
     * @brief Create GenEO matrices with virtual overlap from nonoverlapping Neumann matrix
     *
     * @param go Grid operator representing weak form of PDE
     * @param adapter Virtual overlap adapter
     * @param A Nonoverlapping matrix assembled with Neumann boundaries
     * @return A tuple containing the GenEO matrices extended by virtual overlap and the corresponding partition of unity
     */
    template<typename MODEL, class GO, typename GV, class Matrix, class Vector>
    std::tuple<std::shared_ptr<Matrix>, std::shared_ptr<Vector>, std::vector<int>> setupGenEOMatricesNL(
        MODEL& myModel, const GO& go, NonoverlappingOverlapAdapter<GV, Vector,
        Matrix>& adapter, const typename GO::Jacobian& A, const typename GO::Domain& u,
        bool need_enforced_dirichlet_on_PoU=1, int PoU_restriction=0
      ) {

      using Dune::PDELab::Backend::native;
      Dune::Timer timer_detailed(true);

      using GFS = typename GO::Traits::TrialGridFunctionSpace;
      // using GV = typename GFS::Traits::GridView;

      const GFS& gfs = go.trialGridFunctionSpace();
      const GV& gv = gfs.gridView();

      std::shared_ptr<Matrix> A_extended = adapter.extendMatrix(native(A));

      std::shared_ptr<Vector> part_unity = Dune::makePartitionOfUnity<GV, Matrix, Vector>(adapter, *A_extended); // The one used for A_extended setup


      auto PoU_and_IntBndDofs = Dune::makePartitionOfUnityRestricted<GV, Matrix, Vector>(adapter, *A_extended, PoU_restriction); // The one used later one, potentially restricted on boundaries to have an oversampled subdomain
      std::shared_ptr<Vector> part_unity_R = std::get<0>(PoU_and_IntBndDofs);
      std::vector<int> IntBndDofs = std::get<1>(PoU_and_IntBndDofs);


      using Attribute = Dune::EPISAttribute;
      Dune::AllSet<Attribute> allAttribute;
      auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
      allinterface->build(*adapter.getRemoteIndices(),allAttribute,allAttribute); // all to all communication
      auto communicator = std::shared_ptr<DuneWithRank::BufferedCommunicator>(new DuneWithRank::BufferedCommunicator());
      communicator->build<Vector>(*allinterface);

      Dune::PDELab::MultiVectorBundle<GV, Vector, Matrix> remotePartUnities(adapter);
      remotePartUnities.localVector_ = part_unity;
      communicator->forward<Dune::PDELab::MultiGatherScatter<Dune::PDELab::MultiVectorBundle<GV, Vector, Matrix>>>(remotePartUnities,remotePartUnities); // make function known in other subdomains


      // Assemble fine grid matrix defined only on overlap region
      auto part_unity_restricted = std::make_shared<Vector>(native(A).N());
      adapter.restrictVector(*part_unity, *part_unity_restricted);

      auto es_pou_excluder = std::make_shared<EntitySetPartUnityExcluder<Vector, GV, GFS>> (gfs, part_unity_restricted);
      gfs.entitySet().setExcluder(es_pou_excluder);

      using V = typename GO::Domain;
      using M = typename GO::Jacobian;

      // V x(gfs,0.0); // NOTE: We assume linear problems, so simply set x to zero here!
      // if ( gv.comm().rank()==1){
      // int iter = 0;
      // for(auto val : native(u)){
      //   std::cout << gv.comm().rank() << " u(" << iter << ") : " << val[0] << ", " << val[1] << ", " << val[2] << std::endl;
      //   iter++;
      // }}
      M A_ovlp(go);
      go.jacobian(u, A_ovlp);

      M newmat(go);
      auto part_unity_restricted_ = std::make_shared<Vector>(native(A).N());
      // Provide neighbors with matrices assembled exclusively on respective overlap area
      std::pair<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>> extended_matrices = adapter.lambdaMultiExtendMatrix(native(A), native(A_ovlp), [&](int i){
        std::shared_ptr<Vector> neighbor_part_unity = remotePartUnities.getVectorForRank(i);

        adapter.restrictVector(*neighbor_part_unity, *part_unity_restricted_);

        std::shared_ptr<EntitySetExcluder<Vector, GV>> es_local_pou_excluder = std::make_shared<EntitySetPartUnityExcluder<Vector, GV, GFS>> (gfs, part_unity_restricted_);
        gfs.entitySet().setExcluder(es_local_pou_excluder);

        for (auto rIt=native(newmat).begin(); rIt!=native(newmat).end(); ++rIt)
          for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt) {
            *cIt = 0.0;
          }

        go.jacobian(u,newmat);
        myModel.CorrectDirichletStiffness(gv, gfs, newmat, *part_unity_R);

        return stackobject_to_shared_ptr(native(newmat));
      });

      if (need_enforced_dirichlet_on_PoU==true){
        // Enforce problem's Dirichlet condition on PoU
        const int block_size = Vector::block_type::dimension;
        for (auto rIt=A_extended->begin(); rIt!=A_extended->end(); ++rIt) {
          for(int block_i = 0; block_i < block_size; block_i++){ //loop over block
            bool isDirichlet = true;
            for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
            {
              for(int block_j = 0; block_j<block_size; block_j++){
                if ((rIt.index() != cIt.index() || block_i!=block_j) && (*cIt)[block_i][block_j] != 0.0) {
                  isDirichlet = false;
                  break;
                }
              }
              if(!isDirichlet) break;
            }
            if (isDirichlet) {
              (*part_unity)[rIt.index()][block_i] = .0;
            }
          }
        }
      }

      using ESExcluder = Dune::PDELab::EntitySetExcluder<Vector, GV>;
      gfs.entitySet().setExcluder(std::make_shared<ESExcluder>());

      // return std::make_tuple(extended_matrices.first, part_unity, IntBndDofs);
      return std::make_tuple(extended_matrices.first, part_unity_R, IntBndDofs);
      // return std::make_tuple(A_extended, part_unity_R, IntBndDofs);
    }




    /**
     * @brief Create GenEO matrices with virtual overlap from nonoverlapping Neumann matrix
     *
     * @param go Grid operator representing weak form of PDE
     * @param adapter Virtual overlap adapter
     * @param A Nonoverlapping matrix assembled with Neumann boundaries
     * @return A tuple containing the GenEO matrices extended by virtual overlap and the corresponding partition of unity
     */
    template<class GO, typename GFS, typename GV, class Matrix, class Vector>
    std::tuple<std::shared_ptr<Matrix>, std::shared_ptr<Vector>, std::vector<int>> setupGenEOMatricesNLSerendipity(
              const GO& go, NonoverlappingOverlapAdapterSerendipity<GV, GFS, Vector, Matrix>& adapter,
              const typename GO::Jacobian& A, const typename GO::Domain& u, bool need_enforced_dirichlet_on_PoU=1,
              int PoU_restriction=0) {

      using Dune::PDELab::Backend::native;
      Dune::Timer timer_detailed(true);

      const GFS& gfs = go.trialGridFunctionSpace();
      const GV& gv = gfs.gridView();

      std::shared_ptr<Matrix> A_extended = adapter.extendMatrix(native(A));

      std::shared_ptr<Vector> part_unity = Dune::makePartitionOfUnity<GV, GFS, Matrix, Vector>(adapter, *A_extended); // The one used for A_extended setup
      auto PoU_and_IntBndDofs = Dune::makePartitionOfUnityRestricted<GV, GFS, Matrix, Vector>(adapter, *A_extended, PoU_restriction); // The one used later one, potentially restricted on boundaries to have an oversampled subdomain

      std::shared_ptr<Vector> part_unity_R = std::get<0>(PoU_and_IntBndDofs);
      std::vector<int> IntBndDofs = std::get<1>(PoU_and_IntBndDofs);


      using Attribute = Dune::EPISAttribute;
      Dune::AllSet<Attribute> allAttribute;
      auto allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
      allinterface->build(*adapter.getRemoteIndices(),allAttribute,allAttribute); // all to all communication
      auto communicator = std::shared_ptr<DuneWithRank::BufferedCommunicator>(new DuneWithRank::BufferedCommunicator());
      communicator->build<Vector>(*allinterface);

      Dune::PDELab::MultiVectorBundleSerendipity<GV, GFS, Vector, Matrix> remotePartUnities(adapter);
      remotePartUnities.localVector_ = part_unity;
      communicator->forward<Dune::PDELab::MultiGatherScatter<Dune::PDELab::MultiVectorBundleSerendipity<GV, GFS, Vector, Matrix>>>(remotePartUnities,remotePartUnities); // make function known in other subdomains


      // Assemble fine grid matrix defined only on overlap region
      auto part_unity_restricted = std::make_shared<Vector>(native(A).N());
      adapter.restrictVector(*part_unity, *part_unity_restricted);

      auto es_pou_excluder = std::make_shared<EntitySetPartUnityExcluder<Vector, GV, GFS>> (gfs, part_unity_restricted);
      gfs.entitySet().setExcluder(es_pou_excluder);

      using V = typename GO::Domain;
      using M = typename GO::Jacobian;

      // V x(gfs,0.0); // NOTE: We assume linear problems, so simply set x to zero here!
      // if ( gv.comm().rank()==1){
      // int iter = 0;
      // for(auto val : native(u)){
      //   std::cout << gv.comm().rank() << " u(" << iter << ") : " << val[0] << ", " << val[1] << ", " << val[2] << std::endl;
      //   iter++;
      // }}
      M A_ovlp(go);
      go.jacobian(u, A_ovlp);

      M newmat(go);
      // Provide neighbors with matrices assembled exclusively on respective overlap area
      std::pair<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>> extended_matrices = adapter.lambdaMultiExtendMatrix(native(A), native(A_ovlp), [&](int i){
        std::shared_ptr<Vector> neighbor_part_unity = remotePartUnities.getVectorForRank(i);

        adapter.restrictVector(*neighbor_part_unity, *part_unity_restricted);

        std::shared_ptr<EntitySetExcluder<Vector, GV>> es_local_pou_excluder = std::make_shared<EntitySetPartUnityExcluder<Vector, GV, GFS>> (gfs, part_unity_restricted);
        gfs.entitySet().setExcluder(es_local_pou_excluder);

        for (auto rIt=native(newmat).begin(); rIt!=native(newmat).end(); ++rIt)
          for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt) {
            *cIt = 0.0;
          }

        go.jacobian(u,newmat);

        return stackobject_to_shared_ptr(native(newmat));
      });

      if (need_enforced_dirichlet_on_PoU==true){
        // Enforce problem's Dirichlet condition on PoU
        const int block_size = Vector::block_type::dimension;
        for (auto rIt=A_extended->begin(); rIt!=A_extended->end(); ++rIt) {
          for(int block_i = 0; block_i < block_size; block_i++){ //loop over block
            bool isDirichlet = true;
            for (auto cIt=rIt->begin(); cIt!=rIt->end(); ++cIt)
            {
              for(int block_j = 0; block_j<block_size; block_j++){
                if ((rIt.index() != cIt.index() || block_i!=block_j) && (*cIt)[block_i][block_j] != 0.0) {
                  isDirichlet = false;
                  break;
                }
              }
              if(!isDirichlet) break;
            }
            if (isDirichlet) {
              (*part_unity)[rIt.index()][block_i] = .0;
            }
          }
        }
      }

      using ESExcluder = Dune::PDELab::EntitySetExcluder<Vector, GV>;
      gfs.entitySet().setExcluder(std::make_shared<ESExcluder>());

      // return std::make_tuple(extended_matrices.first, part_unity, IntBndDofs);
      return std::make_tuple(extended_matrices.first, part_unity_R, IntBndDofs);
      // return std::make_tuple(A_extended, part_unity_R, IntBndDofs);
    }
  }
}
#endif