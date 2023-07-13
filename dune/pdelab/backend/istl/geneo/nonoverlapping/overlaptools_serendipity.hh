// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_EXTEND_OVERLAP_SERENDIPITY_TOOLS_HH
#define DUNE_PDELAB_EXTEND_OVERLAP_SERENDIPITY_TOOLS_HH

// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/enumset.hh>
#include<dune/common/parallel/indexset.hh>
#include<dune/common/parallel/plocalindex.hh>
#include<dune/common/parallel/interface.hh>
#include<dune/common/parallel/remoteindices.hh>
#include<dune/common/parallel/communicator.hh>
#include<dune/common/parallel/variablesizecommunicator.hh>

#include "communicator_with_rank.hh"
#include "variablesizecommunicator_with_rank.hh"
#include "overlaptools.hh"

namespace Dune {

  // Attribute type for communication and transfering partition type from the grid
  // enum EPISAttribute {interior=0,border=1,overlap=2,front=3,ghost=4};

  /** A class extending a given index set by adding overlap with respect to a matrix graph
   *
   * \tparam CollectiveCommunication     Type collective communication
   * \tparam GlobalId                    Type representing globally unique ids
   * \tparam Matrix                      ISTL sparse matrix type
   */
  template<typename CollectiveCommunication, typename GlobalId, typename Matrix>
  class ExtendedParallelIndexSetSerendipity
  {
  public:
    // make enum type for attributes public
    using AttributedLocalIndex = Dune::ParallelLocalIndex<EPISAttribute>;
    using ParallelIndexSet = Dune::ParallelIndexSet<GlobalId,AttributedLocalIndex,256>;
    using RemoteIndices = Dune::RemoteIndices<ParallelIndexSet>;
    using LocalIndex = typename Matrix::size_type;

  private:
    // some types we need
    using FieldType = typename Matrix::field_type;

    // local data
    const CollectiveCommunication& comm; // collective communication object
    int rank;     // my rank
    int p;        // number of subdomains in total
    int overlapsize;// overlap in terms of graph distance
    size_t N_orig; // number of degrees of freedom in the input
    size_t N_prec; // number of degrees of freedom of vectors used in preconditioner
    std::vector<LocalIndex> new2old_localindex; // maps new local index to a local index in the original input
    std::vector<LocalIndex> old2new_localindex; // maps local index from input set to its new value or the value N_orig if it is not contained
    std::shared_ptr<ParallelIndexSet> pis; // newly built up index set
    std::shared_ptr<RemoteIndices> si; // shared index information

    // data handle that constructs a set of new indices on each rank
    // - the parallel index set is not modified
    // - this works on the input matrix, so we must care about ghosts
    // - after completion the new indices can be read off
    class ExtenderDataHandle
    {
      const int rank;
      const ParallelIndexSet& pis;
      const Matrix& A;
      const std::vector<LocalIndex>& n2o;
      const std::vector<LocalIndex>& o2n;
      std::vector<GlobalId> globalid;
      std::set<GlobalId> newindices;

    public:
      // export type of items to be sent around
      typedef GlobalId DataType;

      // constructor
      ExtenderDataHandle (int rank_,
                          const ParallelIndexSet& pis_,
                          const Matrix& A_,
                          const std::vector<LocalIndex>& n2o_,
                          const std::vector<LocalIndex>& o2n_)
        : rank(rank_), pis(pis_), A(A_), n2o(n2o_), o2n(o2n_), globalid(pis.size())
      {
        // build a map from local index to global id
        auto endit = pis.end();
        for (auto it=pis.begin(); it!=endit; ++it)
          globalid[it->local().local()] = it->global();
      }

      // enable variable size
      bool fixedsize()
      {
        return false; // we have variable size
      }

      // return size for given (local) index
      std::size_t size (int i)
      {
        // since local indices are always appended at the end (i<n2o.size())
        // indicates that this local index was already in the initial set (including the ghosts)
        // then we have part of the graph and send information about our part of the row
        std::size_t count=1;
        if ((unsigned)i<n2o.size()) // i is a new local index since it comes from the parallel index set
          {
            auto I = n2o[i]; // I is the corresponding index in the original local index set of A
            auto cIt = A[I].begin();
            auto cEndIt = A[I].end();
            for (; cIt!=cEndIt; ++cIt)
              if (cIt.index()<o2n.size()) {
                if (o2n[cIt.index()]<n2o.size()) // this indicates a valid index
                  count++;
              } else std::cout << "BING THIS SHOULD NOT HAPPEN" << std::endl;
          }
        return count;
      }

      // gather to buffer
      template<class B>
      void gather (B& buffer, int i)
      {
        buffer.write(GlobalId()); // write dummy since we want message length >= 1
        if ((unsigned)i<n2o.size()) // i is a new local index since it comes from the parallel index set
          {
            auto I = n2o[i]; // I is the corresponding index in the original local index set of A
            auto cIt = A[I].begin();
            auto cEndIt = A[I].end();
            for (; cIt!=cEndIt; ++cIt)
              if (o2n[cIt.index()]<n2o.size()) // this indicates a valid index
                {
                  buffer.write(globalid[o2n[cIt.index()]]);
                  //std::cout << rank << ": gather " << globalid[i] << "," << globalid[o2n[cIt.index()]] << std::endl;
                }
          }
      }

      template<class B>
      void scatter(B& buffer, int i, int size)
      {
        DataType x;
        buffer.read(x); // read dummy
        // now read size items from buffer
        for (int k=1; k<size; k++)
          {
            buffer.read(x); // here we receive a global index
            if (!pis.exists(x))
              {
                newindices.insert(x);
                // std::cout << rank << ": scatter " << globalid[i] << " -> " << x << std::endl;
              }
          }
      }

      // give out the result
      std::vector<GlobalId> getNewIndices ()
      {
        // convert set to std::vector
        std::vector<GlobalId> temp;
        for (auto it=newindices.begin(); it!=newindices.end(); ++it)
          temp.push_back(*it);
        return temp;
      }
    };

  public:

    //! \brief Constructor.
    ExtendedParallelIndexSetSerendipity (const CollectiveCommunication& comm_,    // collective communication object
                              const std::vector<int>& neighbors,       // ranks with whom we possibly ever share degrees of freedom
                              const Matrix& A,                         // matrix assembled on interior elements and with dirichlet rows replaced
                              const std::vector<EPISAttribute>& partitiontype,   // vector giving partitiontype for each degree of freedom
                              const std::vector<GlobalId>& globalid,   // vector mapping local index to globally unique id
                              int overlap_,                            // layers of overlap to add
                              bool verbose=false)                      // be verbose if true

      : comm(comm_), rank(comm.rank()), p(comm.size()), overlapsize(overlap_), N_orig(A.N()), old2new_localindex(A.N())
    {
      // compute dist of every dof to the boundary (or border dofs) up to a certain depth
      // this is needed later to make all indices public that are potentially in the overlap

      std::vector<int> dist(N_orig); // Number of DoF here different from partitiontype.size --> number of node

      for (size_t i=0; i<N_orig; i++)
        dist[i] = (partitiontype[i]==EPISAttribute::border) ? 0 : (overlapsize+1);

      for (int round=0; round<overlapsize; round++)
        for (size_t i=0; i<N_orig; i++)
          {
            auto cIt = A[i].begin();
            auto cEndIt = A[i].end();
            for (; cIt!=cEndIt; ++cIt)
              dist[i] = std::min(dist[i],dist[cIt.index()]+1);
          }

      // now we can set up a parallel index which is initialized as follows:
      // - it contains only interior and border dofs, ghost dofs are skipped
      // - a new local index is built up so that these are consecutive
      // - we need a map to transfer between the the new local index without ghosts and the original local index including ghosts
      //   -- new2old_localindex takes a new local index and gives an old local index
      //   -- old2new_localindex does the reverse map; this is only needed in setup phase when copying the input matrix
      // - only indices which have dist <= overlapsize are public and may be shared with indices
      pis = std::shared_ptr<ParallelIndexSet>(new ParallelIndexSet());
      pis->beginResize();
      for (size_t i=0; i<N_orig; i++)
        {
          if (partitiontype[i]==EPISAttribute::ghost)
            {
              old2new_localindex[i] = N_orig; // this index will be out of bounds to indicate that it was a ghost
              continue; // skip ghosts
            }
          EPISAttribute attr=partitiontype[i];
          bool pub = (partitiontype[i]==EPISAttribute::overlap)||(partitiontype[i]==EPISAttribute::front)||(dist[i]<=overlapsize);
          new2old_localindex.push_back(i);
          old2new_localindex[i] = new2old_localindex.size()-1;
          pis->add(globalid[i],AttributedLocalIndex(new2old_localindex.size()-1,attr,pub));
        }
      pis->endResize();
      if (verbose)
        std::cout << "rank["<<rank<<"] newsize="<<new2old_localindex.size()<<" oldsize="<<N_orig<<std::endl;

      // now let us extend the index set
      // - each round adds one layer of overlap
      for (int round=0; round<overlapsize; round++)
        {
          // build remote indices information
          RemoteIndices tempsi(*pis,*pis,comm,neighbors);
          // comm.barrier();
          // std::cout << rank << ": before building remote indices" << std::endl;
          // comm.barrier();
          tempsi.template rebuild<false>(); // find shared indices

          // build communication interface
          Dune::AllSet<EPISAttribute> allAttribute;
          Dune::Interface tempinterface;
          tempinterface.build(tempsi,allAttribute,allAttribute);

          // find new indices
          ExtenderDataHandle extdh(rank,*pis,A,new2old_localindex,old2new_localindex);
          DuneWithRank::VariableSizeCommunicator<> varcommunicator(tempinterface);
          // comm.barrier();
          // std::cout << rank << ": before finding new indices" << std::endl;
          // comm.barrier();
          varcommunicator.forward(extdh);
          // comm.barrier();
          // std::cout << rank << ": after finding new indices" << std::endl;
          // comm.barrier();
          auto newindices = extdh.getNewIndices();
          // std::cout << "rank " << rank << " has " << newindices.size() << " new indices" << std::endl;

          // and lets enlarge the index set ...
          LocalIndex localindex = pis->size(); // append new local indices at the end
          pis->beginResize();
          for (auto it=newindices.begin(); it!=newindices.end(); ++it)
            {
              //std::cout << "Proc [" << helper.rank() << "] id=" << globalidset.id(v) << " index=" << indexset.index(v) << std::endl;
              pis->add(*it,AttributedLocalIndex(localindex++,overlap,true));
            }
          pis->endResize();

          if (verbose)
            for (int i=0; i<p; i++)
              {
                comm.barrier();
                if (rank==i) {
                  std::cout <<rank<<": round " << round << " : additional indices="<<newindices.size()<<" new size now="<<pis->size()<<std::endl;
                }
              }
        }
      N_prec = pis->size(); // this now the new size after adding overlap
      //std::cout <<rank<<": input size=" << N_orig <<" new size="<<pis->size()<<std::endl;

      // now find out about shared dofs
      si = std::shared_ptr<RemoteIndices>(new RemoteIndices(*pis,*pis,comm,neighbors));
      si->template rebuild<false>(); // exchanges all sharing information
    }

    //! get number of indices in input index set taken from A.N()
    size_t originalSize () const
    {
      return N_orig;
    }

    //! get number of indices after adding overlap
    size_t extendedSize () const
    {
      return N_prec;
    }

    int overlapSize () const
    {
      return overlapsize;
    }

    //! get the ParallelIndexSet
    std::shared_ptr<ParallelIndexSet> parallelIndexSet () const
    {
      return pis;
    }

    //! get the RemoteIndices
    std::shared_ptr<RemoteIndices> remoteIndices () const
    {
      return si;
    }

    //! returns a map mapping a local index from the extended index set to the original index set
    // This only works for those degrees of freedom that are in both sets
    const std::vector<LocalIndex>& extendedToOriginalLocalIndex () const
    {
      return new2old_localindex;
    }

    //! returns a map mapping a local index from the original index set to the extended index set
    //  This map works on all indices in the original set. If the index is not present in the extended set
    //  (i.e. it corresponded to a ghost dof) then originalSize() is returned by the map.
    const std::vector<LocalIndex>& originalToExtendedLocalIndex () const
    {
      return old2new_localindex;
    }

    const CollectiveCommunication& collectiveCommunication () const
    {
      return comm;
    }
  };


// namespace Dune {
  template<typename GridView, typename GFS, typename Vector, typename Matrix>
  class NonoverlappingOverlapAdapterSerendipity {

    using GlobalId = typename GridView::Grid::GlobalIdSet::IdType;
    using CollectiveCommunication = typename GridView::CollectiveCommunication;
    using EPIS = Dune::ExtendedParallelIndexSetSerendipity<CollectiveCommunication,GlobalId,Matrix>;
    using LocalIndex = typename Vector::size_type;
    using FieldType = typename Vector::field_type;
    using ScalarVector = Dune::BlockVector<Dune::FieldVector<FieldType,1>>;

    using Attribute = EPISAttribute;
    using AttributedLocalIndex = Dune::ParallelLocalIndex<Attribute>;
    using ParallelIndexSet = Dune::ParallelIndexSet<GlobalId,AttributedLocalIndex,256>;
    using RemoteIndices = Dune::RemoteIndices<ParallelIndexSet>;

  public:
    NonoverlappingOverlapAdapterSerendipity(const GridView& gv,       // ranks with whom we share degrees of freedom
                                 const GFS& gfs,          // Matrix used to determine matrix graph
                                 const Matrix& A,          // Matrix used to determine matrix graph
                                 int avg,                  // average number of matrix entries
                                 int overlap)
    : gv_(gv),
    gfs_(gfs),
    //A_(A),
    avg_(avg),
    globalid_(buildGlobalIdVector(gv, gfs)), // TODO: Avoid copy?
    overlap_(overlap),
    epis_(gv.comm(),Dune::PDELab::findNeighboringRanks(gv,overlap_),
          A, buildPartitionTypeVector(gv, gfs),globalid_,overlap_,false)
    {
      new2old_localindex_ = epis_.extendedToOriginalLocalIndex(); // TODO: avoid copy?
    }

    size_t getExtendedSize() const {
      return epis_.parallelIndexSet()->size();
    }
    std::shared_ptr<RemoteIndices> getRemoteIndices() const {
      return epis_.remoteIndices();
    }
    const EPIS& getEpis() const {
      return epis_;
    }
    std::vector<int> findNeighboringRanks() {
      return Dune::PDELab::findNeighboringRanks(gv_,overlap_);
    }

    std::shared_ptr<Matrix> extendMatrix(const Matrix& A) {
      return Dune::copyAndExtendMatrix(epis_,A,avg_,globalid_);
    }

    std::shared_ptr<Matrix> extendMatrix(const Matrix& A_local, const Matrix& A) {
      return Dune::copyAndExtendMatrix(epis_,A_local,A,avg_,globalid_);
    }

    std::pair<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>> multiExtendMatrix(const Matrix& A_local, const Matrix& A2_local, std::map<int,std::shared_ptr<Matrix>> A) {
      return Dune::multiCopyAndExtendMatrix(epis_,A_local,A2_local,A,avg_,globalid_);
    }

    std::pair<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>> lambdaMultiExtendMatrix(const Matrix& A_local, const Matrix& A2_local, std::function<std::shared_ptr<Matrix>(int)> MatrixLambda) {
      return Dune::lambdaMultiCopyAndExtendMatrix(epis_,A_local,A2_local,MatrixLambda,avg_,globalid_);
    }

    const GridView& gridView() const {
      return gv_;
    }

    void extendVector(const Vector& restricted, Vector& extended) const {
      extended = 0.0;
      for (typename Vector::size_type i=0; i<new2old_localindex_.size(); i++)
        extended[i] = restricted[new2old_localindex_[i]];
    }
    void restrictVector(const Vector& extended, Vector& restricted) const {
      for (typename ScalarVector::size_type i=0; i<new2old_localindex_.size(); i++)
        restricted[new2old_localindex_[i]] = extended[i];
    }

    std::vector<LocalIndex> get_old2new_localindex() {
      return epis_.originalToExtendedLocalIndex();
    }

    std::vector<LocalIndex> get_new2old_localindex() {
      return new2old_localindex_;
    }

    std::vector<typename GridView::Grid::GlobalIdSet::IdType> get_globalid() {
      // auto& indexset = gv_.indexSet();
      // auto& globalidset = gv_.grid().globalIdSet();
      // using GlobalID = typename GridView::Grid::GlobalIdSet::IdType;
      // std::vector<GlobalID> globalid(indexset.size(GridView::Grid::dimension));
      // for (const auto& v : vertices(gv_,Dune::Partitions::all))
      //   globalid[indexset.index(v)] = globalidset.id(v);
      // return globalid;

      return buildGlobalIdVector(gv_, gfs_);
    }

  private:

    std::vector<EPISAttribute> buildPartitionTypeVector(const GridView& gv, const GFS& gfs) {

      auto& indexset = gv.indexSet();
      std::vector<Dune::EPISAttribute> partitiontype(indexset.size(GridView::Grid::dimension)+\
                                                     indexset.size(GridView::Grid::dimension-1)); // codim3: vertices + codim2: edges


      using ES = typename GFS::Traits::EntitySet;
      ES es = gfs.entitySet();
      using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
      typedef typename LFS::template Child<0>::Type LFSE;
      LFS lfs(gfs);
      LFSE lfs_e = lfs.template child<0>();

      for (const auto& e : elements(gv,Dune::Partitions::all)){
        lfs_e.bind(e);
        auto& coeffs = lfs_e.finiteElement().localCoefficients();
        for (std::size_t i = 0; i < coeffs.size(); ++i) {
          auto uniquesubindex = es.indexSet().uniqueSubIndex(e, coeffs.localKey(i).subEntity(), coeffs.localKey(i).codim());
          // auto subId = globalidset.subId(e, coeffs.localKey(i).subEntity(), coeffs.localKey(i).codim());

          // auto dof = e.subEntity(coeffs.localKey(i).subEntity(), coeffs.localKey(i).codim());
          // auto dof = e.subEntities(coeffs.localKey(i).subEntity(), coeffs.localKey(i).codim());
          static const int coeff_codim = coeffs.localKey(i).codim();


          if (coeffs.localKey(i).codim() == 3) { // vertices
            auto dof_ptype = e.template subEntity<3>(coeffs.localKey(i).subEntity()).partitionType();
            if (dof_ptype==Dune::InteriorEntity) partitiontype[uniquesubindex] = Dune::EPISAttribute::interior;
            if (dof_ptype==Dune::BorderEntity) partitiontype[uniquesubindex] = Dune::EPISAttribute::border;
            if (dof_ptype==Dune::OverlapEntity) partitiontype[uniquesubindex] = Dune::EPISAttribute::overlap;
            if (dof_ptype==Dune::GhostEntity) partitiontype[uniquesubindex] = Dune::EPISAttribute::ghost;
          } else if (coeffs.localKey(i).codim() == 2) { // edges
            auto dof_ptype = e.template subEntity<2>(coeffs.localKey(i).subEntity()).partitionType();
            if (dof_ptype==Dune::InteriorEntity) partitiontype[uniquesubindex] = Dune::EPISAttribute::interior;
            if (dof_ptype==Dune::BorderEntity) partitiontype[uniquesubindex] = Dune::EPISAttribute::border;
            if (dof_ptype==Dune::OverlapEntity) partitiontype[uniquesubindex] = Dune::EPISAttribute::overlap;
            if (dof_ptype==Dune::GhostEntity) partitiontype[uniquesubindex] = Dune::EPISAttribute::ghost;
          }

        }
      }
      return partitiontype;
    }

    std::vector<typename GridView::Grid::GlobalIdSet::IdType> buildGlobalIdVector(const GridView& gv, const GFS& gfs) {
      auto& indexset = gv.indexSet();
      auto& globalidset = gv.grid().globalIdSet();
      using GlobalID = typename GridView::Grid::GlobalIdSet::IdType;

      std::vector<GlobalID> globalid(indexset.size(GridView::Grid::dimension)+\
                                     indexset.size(GridView::Grid::dimension-1)); // codim3: vertices + codim2: edges

      using ES = typename GFS::Traits::EntitySet;
      ES es = gfs.entitySet();
      using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
      typedef typename LFS::template Child<0>::Type LFSE;
      LFS lfs(gfs);
      LFSE lfs_e = lfs.template child<0>();

      for (const auto& e : elements(gv,Dune::Partitions::all)){
        lfs_e.bind(e);
        auto& coeffs = lfs_e.finiteElement().localCoefficients();
        for (std::size_t i = 0; i < coeffs.size(); ++i) {
          auto uniquesubindex = es.indexSet().uniqueSubIndex(e, coeffs.localKey(i).subEntity(), coeffs.localKey(i).codim());
          auto subId = globalidset.subId(e, coeffs.localKey(i).subEntity(), coeffs.localKey(i).codim());
          globalid[uniquesubindex] = subId;
        }
      }

      // if(gv.comm().rank()==0){
      //   for (int i=0; i<globalid.size(); i++) {
      //     std::cout << "globalid["<<i<<"] = " <<  globalid[i] << std::endl;
      //   }
      // }

      return globalid;
    }

    const GridView& gv_;
    const GFS& gfs_;
    const int avg_;
    const std::vector<GlobalId> globalid_;
    const int overlap_;
    EPIS epis_;
  public:
    std::vector<LocalIndex> new2old_localindex_;
  };

  template<typename GridView, typename GFS, typename Matrix, typename Vector>
  std::shared_ptr<Vector> makePartitionOfUnity(NonoverlappingOverlapAdapterSerendipity<GridView, GFS, Vector, Matrix>& adapter, const Matrix& A) {
    using GlobalId = typename GridView::Grid::GlobalIdSet::IdType;
    using CollectiveCommunication = typename GridView::CollectiveCommunication;
    using EPIS = Dune::ExtendedParallelIndexSetSerendipity<CollectiveCommunication,GlobalId,Matrix>;
    return Dune::makePartitionOfUnity<EPIS, Matrix, Vector>(adapter.getEpis(), A, adapter.getEpis().overlapSize());
  }

  template<typename GridView, typename GFS, typename Matrix, typename Vector>
  std::tuple<std::shared_ptr<Vector>, std::vector<int>> makePartitionOfUnityRestricted(NonoverlappingOverlapAdapterSerendipity<GridView, GFS, Vector, Matrix>& adapter, const Matrix& A, int extra=0) {
    using GlobalId = typename GridView::Grid::GlobalIdSet::IdType;
    using CollectiveCommunication = typename GridView::CollectiveCommunication;
    using EPIS = Dune::ExtendedParallelIndexSetSerendipity<CollectiveCommunication,GlobalId,Matrix>;
    return Dune::makePartitionOfUnityRestricted<EPIS, Matrix, Vector>(adapter.getEpis(), A, adapter.getEpis().overlapSize(), extra);
  }
}


#endif // DUNE_PDELAB_EXTEND_OVERLAP_TOOLS_HH
