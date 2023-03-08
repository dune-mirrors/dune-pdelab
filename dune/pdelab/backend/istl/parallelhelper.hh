// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_PARALLELHELPER_HH
#define DUNE_PDELAB_BACKEND_ISTL_PARALLELHELPER_HH

#include <limits>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/typetraits.hh>

#if HAVE_DUNE_UGGRID && PDELAB_SEQUENTIAL_UG
// We need the UGGrid declaration for the assertion
#include <dune/grid/uggrid.hh>
#endif

//- dune-grid includes
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/utility/entitycommhelper.hh>

#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/access.hh>

#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/istl/vector.hh>
#include <dune/pdelab/backend/istl/utility.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>

#include <dune/pdelab/common/cachedcommmanager.hh>

namespace Dune {
  namespace PDELab {
    namespace ISTL {

      /** --PatternBuilder
       *
       *  BlockMapperInterface:
       *    bool contains( codim ) returns true if dofs are attached to entities with this codim
       *
       *    size_t numEntityDofs( entity ) -> return number of dofs located at this entity
       *
       *    void mapEntityDofs( entity, std::vector< GlobalKey >& indices ) which fills a vector with global keys (vector indices) of the dofs
       *
       */
      template< class Communication, class BlockMapper, InterfaceType CommInterface >
      class PatternBuilder :
        public CommDataHandleIF< PatternBuilder< Communication, BlockMapper, CommInterface >,
                                 typename BlockMapper::GlobalKeyType
                                 /* we currently assume that the rank can be sent in the same format as the indices */
                                 >
      {
      public:
        typedef Communication  CommunicationType;
        typedef BlockMapper    BlockMapperType;

        typedef typename BlockMapperType::GlobalKeyType Index;

        using Pattern = typename ExchangeCommunication<Index>::CommunicationPattern;

      protected:
        const CommunicationType& comm_;
        const BlockMapperType &blockMapper_;

        const int myRank_;
        const int mySize_;

        Pattern &pattern_;

        std::vector<Index> &skip_;

      public:
        PatternBuilder( const CommunicationType& comm,
          const BlockMapperType& blockMapper,
          Pattern &pattern,
          std::vector<Index> &skip)
          : comm_( comm ),
            blockMapper_( blockMapper ),
            myRank_( comm.rank() ),
            mySize_( comm.size() ),
            pattern_( pattern ),
            skip_( skip )
        {}

      protected:
        void sendBackSendMaps()
        {
          typedef MPIPack BufferType;
          typedef MPIFuture< BufferType > FutureType;

          std::map< int, FutureType > recvFutures, sendFutures;

          for( const auto& link : pattern_ )
          {
            auto dest = link.first;
            auto& sendIndexMap_ = link.second.sendIndices;
            BufferType buffer( comm_ );
            buffer << sendIndexMap_; // [ dest ];
            //sendFutures.insert( std::make_pair( dest, comm_.isend(std::move(buffer), dest, 123) ) );
            comm_.send(buffer, dest, 124);
          }

          for( auto& link : pattern_ )
          {
            auto source = link.first;
            //recvFutures.insert( std::make_pair( source, comm_.irecv( BufferType(comm_), source, 123 ) ) );
            auto buffer = comm_.rrecv( BufferType(comm_), source, 124 );
            auto& indices = link.second.sendIndices;
            buffer >> indices;
          }
        }

      public:
        //! destructor
        ~PatternBuilder()
        {
          sendBackSendMaps();
        }

        //! returns true if combination is contained
        bool contains( int dim, int codim ) const
        {
          return blockMapper_.contains( codim );
          // _communication_descriptor.contains(_gfs,dim,codim);
        }

        //! return whether we have a fixed size
        bool fixedSize( int dim, int codim ) const
        {
          return false;
          // _communication_descriptor.fixedSize(_gfs,dim,codim);
        }

        template<typename Idx,
                 std::enable_if_t<IsNumber<Idx>::value, int> = 0>
        static Idx int2index(int i)
        {
          return Idx(i);
        }

#warning TODO better SFINAE
        template<typename Idx,
                 std::enable_if_t<not IsNumber<Idx>::value, int> = 0>
        static Idx int2index(int i)
        {
          Idx idx;
          idx[0] = i;
          return idx;
        }

        template<typename Idx,
                 std::enable_if_t<IsNumber<Idx>::value, int> = 0>
        static int index2int(Idx i)
        {
          return int(i);
        }

#warning TODO better SFINAE
        template<typename Idx,
                 std::enable_if_t<not IsNumber<Idx>::value, int> = 0>
        static int index2int(Idx idx)
        {
          int i;
          return idx[0];
        }

        //! read buffer and apply operation
        template< class MessageBuffer, class Entity >
        void gather( MessageBuffer &buffer, const Entity &entity ) const
        {
          // check whether we are a sending entity
          const auto myPartitionType = entity.partitionType();
          const bool send = EntityCommHelper< CommInterface >::send( myPartitionType );

          // if we send data then send rank and dofs
          if( send )
          {
            // send rank for linkage
            Index rank = int2index<Index>(myRank_);
            buffer.write( rank );

            const unsigned int numDofs = blockMapper_.numEntityDofs( entity );

            std::vector< Index > indices( numDofs );

            // copy all global keys
            blockMapper_.obtainEntityDofs( entity, indices );

            // write global keys to message buffer
            for( unsigned int i = 0; i < numDofs; ++i )
              buffer.write( indices[ i ] );
          }
        }

        //! read buffer and apply operation
        template< class MessageBuffer, class Entity >
        void scatter( MessageBuffer &buffer, const Entity &entity, const size_t dataSize )
        {
          // if data size > 0 then other side is sender
          if( dataSize > 0 )
          {
            // read rank of other side
            Index rankIndex;
            buffer.read( rankIndex );
            int rank = index2int(rankIndex);
            assert( (rank >= 0) && (rank < mySize_) );

            // check whether we are a sending entity
            const auto myPartitionType = entity.partitionType();
            const bool receive = EntityCommHelper< CommInterface >::receive( myPartitionType );

            // insert rank of link into set of links
            // linkStorage_.insert( rank );
            pattern_[rank].rank = rank; // Dune::MPIComm::Link<Index>(rank); // TODO Index typ

            // read indices from stream
            std::vector< Index > indices( dataSize - 1 );
            for(size_t i=0; i<dataSize-1; ++i)
              buffer.read( indices[i] );

            // if we are a receiving entity
            if( receive )
            {
              //////////////////////////////////////////////////////////
              //
              // Problem here: sending and receiving order might differ
              // Solution: sort all dofs after receiving order and send
              // senders dofs back at the end
              //
              //////////////////////////////////////////////////////////

              // if data has been send and we are receive entity
              // then insert indices into send map of rank
              insert( pattern_[rank].sendIndices, indices );
              // insert( sendIndexMap_[rank], indices );

              // resize vector
              indices.resize( blockMapper_.numEntityDofs( entity ) );

              // build local mapping for receiving of dofs
              // copy all global keys
              blockMapper_.obtainEntityDofs( entity, indices );

              insert( pattern_[rank].recvIndices, indices );
              // insert( recvIndexMap_[rank], indices );

              //////////////////////////////////////////////////////////
              // construct disjoint set of DOFs
              //////////////////////////////////////////////////////////
              switch(myPartitionType)
              {
                case InteriorEntity:
                  // we own all interior DOFs, so no need to skip
                  break;
                case BorderEntity:
                  // Assign DOFs to processor with minimum rank
                  // -> skip if our rank is larger than remote rank
                  if (myRank_ > rank) break;
                case OverlapEntity:
                case FrontEntity:
                case GhostEntity:
                  // add all indices belonging to current entity to skip list
                  insert( skip_, indices );
              }

            }
          }
        }

        template <class Vector>
        void insert(Vector& idxMap, const Vector& indices)
        {
          {
            const size_t size = indices.size();
            size_t count = idxMap.size();
            // std::cout << myRank_ << "\tInsert2 " << count << std::endl;

            // reserve memory
            idxMap.resize( count + size );
            assert( idxMap.size() == (count + size) );

            // copy indices to index vector
            for( size_t i = 0; i < size; ++i, ++count )
            {
              #warning generalize to work with arbitrary nested spaces
              // WORKS ONLY FOR SCALAR INDICES: assert( indices[ i ] >= 0 );
              idxMap[ count ] = indices[ i ];
            }
          }
        }

        //! return local dof size to be communicated
        template< class Entity >
        size_t size( const Entity &entity ) const
        {
          const PartitionType myPartitionType = entity.partitionType();
          const bool send = EntityCommHelper< CommInterface >::send( myPartitionType );
          return (send) ? (
            blockMapper_.numEntityDofs( entity ) // dofs
            +1 // rank
            ) : 0;
        }
      };

      /** expected Mapper interface:
       *
       *  - true if any DOFs are associated with entites of codim
       *    `Mapper.contains( codim );`
       *  - number of DOFs associated with entity
       *    `Mapper.numEntityDofs( entity );`
       *  - write multi-indices of DOFs associated with entity into indices
       *    `Mapper.obtainEntityDofs( entity, indices );`
       *
       */
      template< class GridView, class BlockMapper >
      std::tuple<
                 typename ExchangeCommunication<typename BlockMapper::GlobalKeyType>::CommunicationPattern,
                 std::vector<typename BlockMapper::GlobalKeyType>>
      buildCommunicationPatternFromMapper( const GridView& gv, const BlockMapper& blockMapper, const InterfaceType interface )
      {
        typename ExchangeCommunication<typename BlockMapper::GlobalKeyType>::CommunicationPattern pattern;
        std::vector<typename BlockMapper::GlobalKeyType> skip;
        // return Hybrid::switchCases(
        //   std::make_index_sequence<Size::value>{},
        //   {InteriorBorder_All_Interface},
        //   mi[mi.size()-i-1],
        //   [&](auto ii) {
        //     PatternBuilder< CommunicationType, BlockMapper, InteriorBorder_All_Interface >
        //       handle( gv.comm(), blockMapper, pattern);
        //     // make one all to all communication to build up communication pattern
        //     gv.communicate( handle, InteriorBorder_All_Interface , ForwardCommunication );
        //   });

        // TODO replace by switchCases(const Cases& cases, const Value& value, Branches&& branches)
        typedef typename GridView::Communication CommunicationType;
        if( interface == InteriorBorder_All_Interface )
        {
          PatternBuilder< CommunicationType, BlockMapper, InteriorBorder_All_Interface >
            handle( gv.comm(), blockMapper, pattern, skip);
          // make one all to all communication to build up communication pattern
          gv.communicate( handle, InteriorBorder_All_Interface , ForwardCommunication );
        }
        else if( interface == InteriorBorder_InteriorBorder_Interface )
        {
          PatternBuilder< CommunicationType, BlockMapper, InteriorBorder_InteriorBorder_Interface >
            handle( gv.comm(), blockMapper, pattern, skip);
          // make one all to all communication to build up communication pattern
          gv.communicate( handle, InteriorBorder_InteriorBorder_Interface , ForwardCommunication );
        }
        else if( interface == All_All_Interface )
        {
          PatternBuilder< CommunicationType, BlockMapper, All_All_Interface >
            handle( gv.comm(), blockMapper, pattern, skip);
          // make one all to all communication to build up communication pattern
          gv.communicate( handle, All_All_Interface , ForwardCommunication );
        }
        else
          DUNE_THROW( NotImplemented, "CommunicationPattern for the given interface has not been implemented, yet." );
        return std::make_tuple(pattern,skip);
      }

      //! \addtogroup Backend
      //! \ingroup PDELab
      //! \{

#define USEGFS 1

      // Helper class to map DOFs per (sub-)entity
      template<typename GFS>
      struct EntityDOFMapper
      {
        using IndexSet = typename GFS::Traits::GridView::IndexSet;

        using Ordering = typename GFS::Ordering;
        using GlobalKeyType = typename Ordering::Traits::ContainerIndex;

        const GFS& gfs_;
        const IndexSet& indexSet_;
#if not USEGFS
        const int localSize_;
#endif

        // temporary storage
        // container to store local DOFIndices (needing during calculation)
        using DOFIndex = typename Ordering::Traits::DOFIndex;
        static const std::size_t leaf_count = TypeTree::TreeInfo<Ordering>::leafCount;
        using Offsets = std::array<std::size_t,leaf_count + 1>;
        mutable std::vector<DOFIndex> dofIndices_;
        mutable Offsets offsets_;

        EntityDOFMapper( const GFS& gfs ) :
          gfs_(gfs), indexSet_(gfs_.gridView().indexSet())
#if not USEGFS
          , localSize_(1) /* hard coded for our initial example */
#endif
        {
        }

        bool contains( const int codim ) const {
#if USEGFS
          return gfs_.dataHandleContains(codim);
#else
          return codim == 0;
#endif
        }

        template <class Entity>
        unsigned int numEntityDofs(const Entity& entity) const {
          // we should add an option to not communicate indices of the
          // leaf entities (if the blocking is appropriate, e.g. dG
          // spaces or PowerGFS)
#if USEGFS
          return gfs_.dataHandleSize(entity);
#else
          assert(Entity::codimension == 0);
          return localSize_;
#endif
        }

        // we use blocked vectors and only store the cell indices
        template <class Entity, class Vector> // Vector = std::vector< GlobalKeyType >
        void obtainEntityDofs( const Entity& entity, Vector& containerIndices ) const
        {
          assert( containerIndices.size() == numEntityDofs( entity ) ); // numEntityDofs
#if USEGFS
          assert(gfs_.dataHandleContains(Entity::codimension));

          std::fill(offsets_.begin(),offsets_.end(),0);
          gfs_.dataHandleIndices(
            entity,
            containerIndices,
            dofIndices_,
            offsets_.begin(),
            std::integral_constant<bool,false>());

          // -- reverse indices
          // PDELab uses (inner -> outer) ordering [for technical reasons],
          // whereas dune-istl expects (outer -> inner) [which is more intuitive].
          for (auto it = containerIndices.end()-numEntityDofs(entity);
               it < containerIndices.end(); ++it)
            reverseMultiIndex(*it);
#else
          if (localSize_ == 1) ////// SCALAR
          {
            containerIndices[ 0 ] = indexSet_.index( entity );
          }
          else /////// BLOCKED
          {
            for (unsigned int i=0; i<localSize_; i++) {
              containerIndices[i].set(indexSet_.index( entity ));
              vector[i].push_back(i);
            }
          }
#endif
        }

      private:
        template<typename MI>
        void reverseMultiIndex(MI & mi) const
        {
          auto S = mi.size();
          for (unsigned int i=0; i<S/2; i++)
            std::swap(mi[i], mi[S-i-1]);
        }
      };

      //========================================================
      // A parallel helper class providing a nonoverlapping
      // decomposition of all degrees of freedom
      //========================================================

      template<typename GFS>
      class ParallelHelper
      {

        //! Type for storing rank values.
        typedef int RankIndex;

        //! Type used to store owner rank values of all DOFs.
        using RankVector = Dune::PDELab::Backend::Vector<GFS,RankIndex>;

        //! ContainerIndex of the underlying GridFunctionSpace.
        typedef typename GFS::Ordering::Traits::ContainerIndex ContainerIndex;

        // ...
        using DOFMapper = EntityDOFMapper<GFS>;
        using CachedComm = Dune::ExchangeCommunication< typename DOFMapper::GlobalKeyType >;
        using SkipList = std::vector< typename DOFMapper::GlobalKeyType >;

      public:

        inline static bool useCaches = false;

        ParallelHelper (const GFS& gfs, int verbose = 1)
          : _gfs(gfs)
          , _rank(gfs.gridView().comm().rank())
          , _rank_partition(gfs,_rank)
          , _verbose(verbose)
        {

          // Let's try to be clever and reduce the communication overhead by picking the smallest
          // possible communication interface depending on the overlap structure of the GFS.
          // FIXME: Switch to simple comparison as soon as dune-grid:1b3e83ec0 is reliably available!
          if (gfs.entitySet().partitions().value == Partitions::interiorBorder.value)
            {
              // The GFS only spans the interior and border partitions, so we can skip sending or
              // receiving anything else.
              _interiorBorder_all_interface = InteriorBorder_InteriorBorder_Interface;
              _all_all_interface = InteriorBorder_InteriorBorder_Interface;
            }
          else
            {
              // In general, we have to transmit more.
              _interiorBorder_all_interface = InteriorBorder_All_Interface;
              _all_all_interface = All_All_Interface;
            }

          if (useCaches)
          {
            DOFMapper dofmapper(gfs);
            auto [all_all_pattern, skip] =
              buildCommunicationPatternFromMapper(_gfs.gridView(), dofmapper, All_All_Interface);
            _all_all_comm.init(_gfs.gridView().comm());
            _all_all_comm.setCommunicationPattern(all_all_pattern);
            auto [interiorBorder_all_pattern, _ignore] =
              buildCommunicationPatternFromMapper(_gfs.gridView(), dofmapper, InteriorBorder_All_Interface);
            _interiorBorder_all_comm.init(_gfs.gridView().comm());
            _interiorBorder_all_comm.setCommunicationPattern(interiorBorder_all_pattern);
            //
            std::sort(skip.begin(), skip.end());
            std::unique(skip.begin(), skip.end());
            std::swap(_skip_indices, skip);
          }

          if (_gfs.gridView().comm().size()>1)
            {
              // create disjoint DOF partitioning
              //            GFSDataHandle<GFS,RankVector,DisjointPartitioningGatherScatter<RankIndex> >
              //  ibdh(_gfs,_rank_partition,DisjointPartitioningGatherScatter<RankIndex>(_rank));
              DisjointPartitioningDataHandle<GFS,RankVector> pdh(_gfs,_rank_partition);
              _gfs.gridView().communicate(pdh,_interiorBorder_all_interface,Dune::ForwardCommunication);
            }

          // Generate list of neighbors' ranks
          if(useCaches) {
            auto& pattern = _all_all_comm.pattern();
            for (const auto & [rank,link] : pattern)
              _neighbor_ranks.push_back(rank);
          }
          else {
            std::set<RankIndex> rank_set;
            for (RankIndex rank : _rank_partition)
              if (rank != _rank)
                rank_set.insert(rank);

            for (RankIndex rank : rank_set)
              _neighbor_ranks.push_back(rank);
          }
        }

        const CachedComm& allToAllCommunication() const
        {
          return _all_all_comm;
        }

        const CachedComm& interiorBorderToAllCommunication() const
        {
          return _interiorBorder_all_comm;
        }

        //! Returns a sorted list of the ranks of all neighboring processes
        const std::vector<RankIndex>& getNeighborRanks() const
        {
          return _neighbor_ranks;
        }

        //! Mask out all DOFs not owned by the current process with 0.
        template<typename X>
        void maskForeignDOFs(X& x) const
        {
          using Backend::native;
          if (useCaches)
          {
            ::Dune::ISTL::forEachIndex (
              [](auto & val, auto & mi) {
                val = 0.0; // set skipped entries to zero
              },
              _skip_indices, native(x)
              );
          }
          {
            // dispatch to implementation.
            maskForeignDOFs(ISTL::container_tag(native(x)),native(x),native(_rank_partition));
          }
        }

      private:

        // Implementation for block vector; recursively masks blocks.
        template<typename X, typename Mask>
        void maskForeignDOFs(ISTL::tags::block_vector, X& x, const Mask& mask) const
        {
          typename Mask::const_iterator mask_it = mask.begin();
          for (typename X::iterator it = x.begin(),
                 end_it = x.end();
               it != end_it;
               ++it, ++mask_it)
            maskForeignDOFs(ISTL::container_tag(*it),*it,*mask_it);
        }

        // Implementation for field vector, iterates over entries and masks them individually.
        template<typename X, typename Mask>
        void maskForeignDOFs(ISTL::tags::field_vector, X& x, const Mask& mask) const
        {
          typename Mask::const_iterator mask_it = mask.begin();
          for (typename X::iterator it = x.begin(),
                 end_it = x.end();
               it != end_it;
               ++it, ++mask_it)
            *it = (*mask_it == _rank ? *it : typename X::field_type(0));
        }

      public:

        //! Tests whether the given index is owned by this process.
        bool owned(const ContainerIndex& i) const
        {
          return _rank_partition[i] == _rank;
        }

        //! Calculates the (rank-local) dot product of x and y on the disjoint partition defined by the helper.
        template<typename X, typename Y>
        typename PromotionTraits<
          typename X::field_type,
          typename Y::field_type
          >::PromotedType
        disjointDot(const X& x, const Y& y) const
        {
          using Backend::native;
          if (useCaches)
          {
            auto sp = native(x).dot(native(y));
            decltype(sp) skip = 0.0;
            ::Dune::ISTL::forEachIndex (
              [&skip](auto & v, auto & w, auto & mi) { skip += dot(v,w); },
              _skip_indices, native(x), native(y)
              );
            return sp-skip;
          }
          else
          {
            return disjointDot(ISTL::container_tag(native(x)),
                               native(x),
                               native(y),
                               native(_rank_partition)
                               );
          }
        }

      private:

        // Implementation for BlockVector, collects the result of recursively
        // invoking the algorithm on the vector blocks.
        template<typename X, typename Y, typename Mask>
        typename PromotionTraits<
        typename X::field_type,
        typename Y::field_type
        >::PromotedType
        disjointDot(ISTL::tags::block_vector, const X& x, const Y& y, const Mask& mask) const
        {
          typedef typename PromotionTraits<
            typename X::field_type,
            typename Y::field_type
            >::PromotedType result_type;

          result_type r(0);

          typename Y::const_iterator y_it = y.begin();
          typename Mask::const_iterator mask_it = mask.begin();
          for (typename X::const_iterator x_it = x.begin(),
                 end_it = x.end();
               x_it != end_it;
               ++x_it, ++y_it,  ++mask_it)
            r += disjointDot(ISTL::container_tag(*x_it),*x_it,*y_it,*mask_it);

          return r;
        }

        // Implementation for FieldVector, iterates over the entries and calls Dune::dot() for DOFs
        // associated with the current rank.
        template<typename X, typename Y, typename Mask>
        typename PromotionTraits<
          typename X::field_type,
          typename Y::field_type
          >::PromotedType
        disjointDot(ISTL::tags::field_vector, const X& x, const Y& y, const Mask& mask) const
        {
          typedef typename PromotionTraits<
            typename X::field_type,
            typename Y::field_type
            >::PromotedType result_type;

          result_type r(0);

          typename Y::const_iterator y_it = y.begin();
          typename Mask::const_iterator mask_it = mask.begin();
          for (typename X::const_iterator x_it = x.begin(),
                 end_it = x.end();
               x_it != end_it;
               ++x_it, ++y_it, ++mask_it)
            r += (*mask_it == _rank ? Dune::dot(*x_it,*y_it) : result_type(0));

          return r;
        }

      public:

        //! Returns the MPI rank of this process.
        RankIndex rank() const
        {
          return _rank;
        }

#if HAVE_MPI

        //! Makes the matrix consistent and creates the parallel information for AMG.
        /**
         * This function accomplishes two things:
         *
         * 1. Makes the matrix consistent w.r.t. to the disjoint partitioning of the DOF space,
         *    i.e. aggregates matrix entries for border entries from neighboring ranks.
         *
         * 2. Sets up the parallel communication information for AMG.
         *
         * \warning  This function silenty assumes that the matrix only has a single level
         *           of blocking and will not work correctly otherwise. Also note that AMG
         *           will only work correctly for P1 discretisations.
         *
         * \param m  The PDELab matrix container.
         * \param c  The parallel information object providing index set, interfaces and
         *           communicators.
         */
        template<typename MatrixType, typename Comm>
        void createIndexSetAndProjectForAMG(MatrixType& m, Comm& c);

      private:

        // Checks whether a matrix block is owned by the current process. Used for the AMG
        // construction and thus assumes a single level of blocking and blocks with ownership
        // restricted to a single DOF.
        bool owned_for_amg(std::size_t i) const
        {
          return Backend::native(_rank_partition)[i][0] == _rank;
        }

#endif // HAVE_MPI

      private:

        const GFS& _gfs;
        const RankIndex _rank;
        RankVector _rank_partition; // vector to identify unique decomposition
        std::vector<RankIndex> _neighbor_ranks; // list of neighbors' ranks
        int _verbose; //verbosity

        //! The actual communication interface used when algorithm requires InteriorBorder_All_Interface.
        InterfaceType _interiorBorder_all_interface;

        //! The actual communication interface used when algorithm requires All_All_Interface.
        InterfaceType _all_all_interface;

        //! cached communicators
        CachedComm _all_all_comm;
        CachedComm _interiorBorder_all_comm;

        //! list of DOF indices not owned by the local rank
        SkipList _skip_indices;
      };

#if HAVE_MPI

      template<typename GFS>
      template<typename M, typename C>
      void ParallelHelper<GFS>::createIndexSetAndProjectForAMG(M& m, C& c)
      {

        using Backend::native;

        const bool is_bcrs_matrix =
          std::is_same<
            typename ISTL::tags::container<
              Backend::Native<M>
              >::type::base_tag,
          ISTL::tags::bcrs_matrix
          >::value;

        const bool block_type_is_field_matrix =
          std::is_same<
            typename ISTL::tags::container<
              typename Backend::Native<M>::block_type
              >::type::base_tag,
          ISTL::tags::field_matrix
          >::value;

        // We assume M to be a BCRSMatrix in the following, so better check for that
        static_assert(is_bcrs_matrix && block_type_is_field_matrix, "matrix structure not compatible with AMG");

        // ********************************************************************************
        // In the following, the code will always assume that all DOFs stored in a single
        // block of the BCRSMatrix are attached to the same entity and can be handled
        // identically. For that reason, the code often restricts itself to inspecting the
        // first entry of the blocks in the diverse BlockVectors.
        // ********************************************************************************

        typedef typename GFS::Traits::GridViewType GV;
        typedef typename RankVector::size_type size_type;
        const GV& gv = _gfs.gridView();

        // Do we need to communicate at all?
        const bool need_communication = _gfs.gridView().comm().size() > 1;

        // First find out which dofs we share with other processors
        using BoolVector = Backend::Vector<GFS,bool>;
        BoolVector sharedDOF(_gfs, false);

        if (need_communication)
          {
            SharedDOFDataHandle<GFS,BoolVector> data_handle(_gfs,sharedDOF,false);
            _gfs.gridView().communicate(data_handle,_all_all_interface,Dune::ForwardCommunication);
          }

        // Count shared dofs that we own
        typedef typename C::ParallelIndexSet::GlobalIndex GlobalIndex;
        GlobalIndex count = 0;

        for (size_type i = 0; i < sharedDOF.N(); ++i)
          if (owned_for_amg(i) && native(sharedDOF)[i][0])
            ++count;

        dverb << gv.comm().rank() << ": shared block count is " << count.touint() << std::endl;

        // Communicate per-rank count of owned and shared DOFs to all processes.
        std::vector<GlobalIndex> counts(_gfs.gridView().comm().size());
        _gfs.gridView().comm().allgather(&count, 1, &(counts[0]));

        // Compute start index start_p = \sum_{i=0}^{i<p} counts_i
        GlobalIndex start = std::accumulate(counts.begin(),counts.begin() + _rank,GlobalIndex(0));

        using GIVector = Dune::PDELab::Backend::Vector<GFS,GlobalIndex>;
        GIVector scalarIndices(_gfs, std::numeric_limits<GlobalIndex>::max());

        for (size_type i = 0; i < sharedDOF.N(); ++i)
          if (owned_for_amg(i) && native(sharedDOF)[i][0])
            {
              native(scalarIndices)[i][0] = start;
              ++start;
            }

        // Publish global indices for the shared DOFS to other processors.
        if (need_communication)
          {
            MinDataHandle<GFS,GIVector> data_handle(_gfs,scalarIndices);
            _gfs.gridView().communicate(data_handle,_interiorBorder_all_interface,Dune::ForwardCommunication);
          }

        // Setup the index set
        c.indexSet().beginResize();
        for (size_type i=0; i<scalarIndices.N(); ++i)
          {
            Dune::OwnerOverlapCopyAttributeSet::AttributeSet attr;
            if(native(scalarIndices)[i][0] != std::numeric_limits<GlobalIndex>::max())
              {
                // global index exist in index set
                if (owned_for_amg(i))
                  {
                    // This dof is managed by us.
                    attr = Dune::OwnerOverlapCopyAttributeSet::owner;
                  }
                else
                  {
                    attr = Dune::OwnerOverlapCopyAttributeSet::copy;
                  }
                c.indexSet().add(native(scalarIndices)[i][0], typename C::ParallelIndexSet::LocalIndex(i,attr));
              }
          }
        c.indexSet().endResize();

        // Compute neighbors using communication
        std::set<int> neighbors;

        if (need_communication)
          {
            GFSNeighborDataHandle<GFS,int> data_handle(_gfs,_rank,neighbors);
            _gfs.gridView().communicate(data_handle,_all_all_interface,Dune::ForwardCommunication);
          }

        c.remoteIndices().setNeighbours(neighbors);
        c.remoteIndices().template rebuild<false>();
      }

#endif // HAVE_MPI

      template<int s, bool isFakeMPIHelper>
      struct CommSelector
      {
        typedef Dune::Amg::SequentialInformation type;
      };


#if HAVE_MPI

      // Need MPI for OwnerOverlapCopyCommunication
      template<int s>
      struct CommSelector<s,false>
      {
        typedef OwnerOverlapCopyCommunication<bigunsignedint<s>,int> type;
      };

#endif // HAVE_MPI

      template<typename T>
      void assertParallelUG(T comm)
      {}

#if HAVE_DUNE_UGGRID && PDELAB_SEQUENTIAL_UG
      template<int dim>
      void assertParallelUG(Dune::CollectiveCommunication<Dune::UGGrid<dim> > comm)
      {
        static_assert(Dune::AlwaysFalse<Dune::UGGrid<dim> >::value, "Using sequential UG in parallel environment");
      };
#endif
      //! \} group Backend

    } // namespace ISTL
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_PARALLELHELPER_HH
