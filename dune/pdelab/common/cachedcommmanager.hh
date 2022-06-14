#ifndef DUNE_CACHED_COMMUNICATION_MANAGER_HH
#define DUNE_CACHED_COMMUNICATION_MANAGER_HH

#include <cassert>
#include <cstddef>

//- system includes
#include <iostream>
#include <map>
#include <queue>
#include <memory>
#include <vector>

//- dune-common includes
#include <dune/common/math.hh>
#include <dune/common/timer.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/mpicommunication.hh>
#include <dune/common/parallel/mpipack.hh>
#include <dune/common/parallel/mpifuture.hh>

//- dune-grid includes
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/utility/entitycommhelper.hh>

/*
  Issues Grid-Comm Interface

  * grid knows source rank, but doesn't tell the user
  * difficult to collect data in subentities without knowing the entity
  *

 */

namespace Dune
{

    namespace MPIComm
    {

        template<typename IndexType>
        struct Link
        {
            using Indices = std::vector< IndexType >;
            int rank;
            Indices recvIndices;
            Indices sendIndices;

            Link () : rank(-1) {
                assert(sendIndices.size() == 0);
            }
            Link (int r) : rank(r) {}
        };

        // todo, add rank, sort by rank in Link
        template<typename IndexType>
        using CommunicationPattern = std::map<int, Link<IndexType>>; // rank -> Link

    } // end namespace MPIComm

    /** @addtogroup Communication Communication
        @{
     **/

    /** \brief Handle "halo"-exchange between neighboring processes
     *
     *
     * CommunicationPattern is a convenience class to build up a map
     * of all dofs of entities to be exchanged during a communication procedure.
     * This speeds up the communication procedure, because no grid traversal is
     * necessary anymore to exchange data. This class is singleton for different
     * discrete function spaces, depending on the BlockMapper.
     */
    template< class IndexType >
    class CommunicationPattern
    {
    public:
      //! type of block mapper of discrete function space (may be the same for
      //! different space (i.e. various DG spaces)
        // typedef BlockMapper BlockMapperType;
      using GlobalKeyType = IndexType; // typedef typename BlockMapperType :: GlobalKeyType

    protected:
      typedef Dune::MPIHelper::MPICommunicator MPICommunicatorType;
      typedef Dune::Communication< MPICommunicatorType > CommunicationType;

    protected:

      MPIComm::CommunicationPattern<Index> pattern_;

      // communicator
      std::unique_ptr< CommunicationType > comm_;

      // exchange time
      mutable double exchangeTime_;
      // setup time
      mutable double buildTime_;

      //! know grid sequence number
      int sequence_;

      mutable int nonBlockingObjects_ ;

    protected:
      template< class Communication, class BlockMapper, InterfaceType CommInterface >
      class PatternBuilder;

      /////////////////////////////////////////////////////////////////
      //  begin NonBlockingCommunication
      /////////////////////////////////////////////////////////////////

      class NonBlockingCommunication
      {
        typedef CommunicationPattern < IndexType > CommunicationPatternType;

        // create an unique tag for the communication
        DUNE_EXPORT static int getMessageTag()
        {
          enum { initial = 665 };
          static int tagCounter = initial ;
          ++ tagCounter;
          int messageTag = tagCounter ;

          // avoid overflow
          if( messageTag < 0 )
          {
            messageTag = initial ;
            tagCounter = initial ;
          }
          return messageTag;
        }

      public:
        NonBlockingCommunication( const CommunicationPatternType& dependencyCache )
          : dependencyCache_( dependencyCache ),
            exchangeTime_( 0.0 ),
            tag_( getMessageTag() )
        {
          // notify dependency cache of open communication
          dependencyCache_.attachComm();
        }

        // copy constructor
        NonBlockingCommunication( const NonBlockingCommunication& other )
          : dependencyCache_( other.dependencyCache_ ),
            exchangeTime_( 0.0 ),
            tag_( other.tag_ )
        {
          // notify dependency cache of open communication
          dependencyCache_.attachComm();
        }

        ~NonBlockingCommunication()
        {
          // notify dependency cache that comm is finished
          dependencyCache_.detachComm() ;
        }

        template < class Data >
        void send( const Data& data )
        {
          sendFutures_.clear();
          recvFutures_.clear();

          // take time
          Dune::Timer sendTimer ;

          const auto& comm = dependencyCache_.comm();

          auto& pattern = dependencyCache_.pattern();
          // for( const auto& dest : linkStorage )
          for( auto& link : pattern )
          {
            int dest = link.first;
            BufferType buffer( comm );
            pack( dest, buffer, data);
            sendFutures_.insert( std::make_pair( dest, comm.isend(std::move(buffer), dest, tag_) ) );
            //comm.send(buffer, dest, tag_);
          }

          // store time needed for sending
          exchangeTime_ = sendTimer.elapsed();
        }

        //! receive data for discrete function and given operation
        template < class Data, class Operation >
        double receive( Data& data, const Operation& operation )
        {
          // take time
          Dune::Timer recvTimer ;

          typedef typename Data::value_type value_type;
          typedef typename value_type::field_type DofType;

          static const int blockSize = value_type::dimension; // data[0].size();
          const size_t dataSize = static_cast< std::size_t >( blockSize * sizeof( DofType ) );

          const auto& comm = dependencyCache_.comm();
          for( const auto& link : dependencyCache_.pattern() )
          {
            int source = link.first;
            const size_t bufSize = dependencyCache_.recvBufferSize( source ) * dataSize;
            recvFutures_.insert( std::make_pair( source, comm.irecv( BufferType(comm, bufSize), source, tag_ ) ) );
            //auto buffer = comm.rrecv( BufferType(comm), source, tag_ );
            //unpack( source, buffer, data, operation );
          }

          for( auto& [rank, future] : recvFutures_ )
          {
            auto buffer = future.get();
            unpack( rank, buffer, data, operation );
          }

          // wait for all sends to finish
          for( auto& [rank, future] : sendFutures_ )
          {
            if( !future.ready() )
              future.wait();
          }

          // store time needed for sending
          exchangeTime_ += recvTimer.elapsed();
          return exchangeTime_;
        }

        //! receive method with default operation
        template < class Data >
        double receive( Data& data )
        {
          // TODO
          // get type of default operation
          auto op = [](const double& a, double& b) { b = a; };
          return receive( data, op );
        }

      protected:
        template <class Buffer, class Data>
        void pack( const int rank, Buffer& buffer, const Data& data )
        {
          // write data of discrete function to message buffer
          dependencyCache_.writeBuffer( rank, buffer, data );
        }

        template <class Buffer, class Data, class Operation>
        void unpack( const int rank, Buffer& buffer,
                     Data& data, const Operation& operation )
        {
          // read data of discrete function from message buffer
          dependencyCache_.readBuffer( rank, buffer, data, operation );
        }

      protected:
        const CommunicationPatternType& dependencyCache_;
        typedef MPIPack BufferType;
        typedef MPIFuture< BufferType > FutureType;

        std::map< int, FutureType > recvFutures_;
        std::map< int, FutureType > sendFutures_;

        double exchangeTime_ ;
        const int tag_;
      };

    public:
      typedef NonBlockingCommunication NonBlockingCommunicationType;

      //! return object for non-blocking communication
      template <class Space>
      NonBlockingCommunicationType nonBlockingCommunication( const Space& space )
      {
        // create non-blocking communication object
        return NonBlockingCommunicationType( space, *this );
      }
      /////////////////////////////////////////////////////////////////
      //  end NonBlockingCommunication
      /////////////////////////////////////////////////////////////////

      //! constructor taking communicator object
      CommunicationPattern()
      : pattern_(),
        comm_(),
        exchangeTime_( 0.0 ),
        buildTime_( 0.0 ),
        sequence_( -1 ),
        nonBlockingObjects_( 0 )
      {
      }

      template <class Communication>
      void init( const Communication& comm )
      {
        if( ! comm_ )
        {
          comm_.reset( new CommunicationType( comm ) );
        }
      }

      auto& pattern() { return pattern_; }

      const auto& pattern() const { return pattern_; }

      operator bool() const
      {
        return pattern_.size() > 0;
      }

      // no copying
      CommunicationPattern( const CommunicationPattern & ) = delete;

      // notify for open non-blocking communications
      void attachComm() const
      {
        ++nonBlockingObjects_;
      }

      // notify for finished non-blocking communication
      void detachComm() const
      {
        --nonBlockingObjects_;
        assert( nonBlockingObjects_ >= 0 );
      }

      bool noOpenCommunications() const
      {
        return true ;
      }

      size_t recvBufferSize( const int rank ) const
      {
        auto it = pattern_.find( rank );
        const auto &indexMap = it->second.recvIndices;
        return indexMap.size();
      }

    public:

      void setCommunicationPattern (Dune::MPIComm::CommunicationPattern<GlobalKeyType> /* TODO */)
      {
      }

      /** \brief Rebuild underlying exchange dof mapping.
       *  \note: Different spaces may have the same exchange dof mapping!
       */
      template <class GridView, class BlockMapperType>
      inline void rebuild( const GridView& gridView,
                           const BlockMapperType& blockMapper,
                           const InterfaceType interface,
                           const bool force = true )
      {
        const auto& comm = gridView.comm();

        // only in parallel we have to do something
        if( comm.size() <= 1 ) return;

        // make sure all non-blocking communications have been finished by now
        assert( noOpenCommunications() );

        // check whether grid has changed.
        if( force )
        {
          // create communicator
          init( comm );

          // take timer needed for rebuild
          Dune::Timer buildTime;

          // rebuild maps holding exchange dof information
          buildMaps( gridView, blockMapper, interface );

          // store time needed
          buildTime_ = buildTime.elapsed();
        }

        // for (int r = 0; r < gridView.comm().size(); r++)
        // {
        //     if (r == gridView.comm().rank())
        //         for (auto && link : pattern_)
        //         {
        //             int rank = link.first;
        //             const auto & Old = sendIndexMap_[rank];
        //             const auto & New = link.second.sendIndices;
        //             assert(Old.size() == New.size());
        //             for (int i=0; i<Old.size(); i++)
        //                 std::cout << gridView.comm().rank() << "/" << rank << " ::::::: "
        //                           << Old[i] << " <-> " << New[i]
        //                           << std::endl;
        //         }
        //     gridView.comm().barrier();
        // }
      }

      //! exchange data of discrete function
      template< class Data, class Operation >
      inline void exchange( Data &data, const Operation& operation ) const;

      //! return reference to communication object
      inline CommunicationType &comm()
      {
        assert( comm_ );
        return *comm_;
      }

      //! return reference to communication object
      inline const CommunicationType &comm() const
      {
        assert( comm_ );
        return *comm_;
      }

    protected:
      // build linkage and index maps
      template < class GridView, class BlockMapper >
      inline void buildMaps( const GridView& gv, const BlockMapper& blockMapper, const InterfaceType interface );

      // check consistency of maps
      inline void checkConsistency();

      template< class GridView, class Comm, class BlockMapper, InterfaceType CI >
      inline void buildMaps( const GridView& gv, PatternBuilder< Comm, BlockMapper, CI > &handle );

    protected:
      // serialize data of DataImp& vector to object stream
      // --writeBuffer
      template< class Buffer, class Data >
      inline void writeBuffer( const int dest,
                               Buffer &buffer,
                               const Data &data ) const
      {
        // auto it = sendIndexMap_.find( dest );
        // const auto &indexMap = it->second;
        auto it = pattern_.find( dest );
        const auto &indexMap = it->second.sendIndices;
        const int size = indexMap.size();

        //typedef typename Data :: DofType DofType;
        typedef typename Data::value_type value_type;
        typedef typename value_type :: field_type DofType;

        {
          //static const int blockSize = Data::blockSize;
          static const int blockSize = value_type::dimension;
          buffer.enlarge( size * blockSize * sizeof( DofType ) );
          for( int i = 0; i < size; ++i )
          {
            const auto &block = data[ indexMap[ i ] ];
            for( int k=0; k<blockSize; ++k )
              buffer << block[ k ];
          }
        }
      }

      // deserialize data from object stream to DataImp& data vector
      // --readBuffer
      template< class Buffer, class Data, class Operation >
      inline void readBuffer( const int source,
                              Buffer& buffer,
                              Data &data,
                              const Operation& operation ) const
      {
        static_assert( ! std::is_pointer< Operation > :: value,
                       "CommunicationPattern::readBuffer: Operation needs to be a reference!");

        // get index map of rank belonging to source
        // auto it = recvIndexMap_.find( source );
        // const auto &indexMap = it->second;
        auto it = pattern_.find( source );
        const auto &indexMap = it->second.recvIndices;

        const int size = indexMap.size();

        typedef typename Data::value_type value_type;
        typedef typename value_type :: field_type DofType;

        {
          static const int blockSize = value_type::dimension;
          assert( static_cast< std::size_t >( size * blockSize * sizeof( DofType ) ) <= static_cast< std::size_t >( (buffer.size()-buffer.tell()) ) );
          for( int i = 0; i < size; ++i )
          {
            auto &&block = data[ indexMap[ i ] ];
            for( int k=0; k<blockSize; ++k )
            {
              DofType value;
              buffer >> value;
              // apply operation, i.e. COPY, ADD, etc.
              operation( value, block[ k ] );
            }
          }
        }
      }
    };

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
    template< class IndexType >
    template< class Communication, class BlockMapper, InterfaceType CommInterface >
    class CommunicationPattern< IndexType > :: PatternBuilder
    : public CommDataHandleIF
    < PatternBuilder< Communication, BlockMapper, CommInterface >,
      IndexType /* we currently assume that the rank can be sent in the same format as the indices */ >
    {
    public:
      typedef Communication  CommunicationType;
      typedef BlockMapper    BlockMapperType;

      typedef typename BlockMapperType::GlobalKeyType GlobalKeyType;

      typedef GlobalKeyType DataType;
      using Pattern = MPIComm::CommunicationPattern<DataType>;

    protected:
      const CommunicationType& comm_;
      const BlockMapperType &blockMapper_;

      const GlobalKeyType myRank_;
      const GlobalKeyType mySize_;

      Pattern &pattern_;
      // LinkStorageType &linkStorage_;
      // IndexVectorMapType &sendIndexMap_;
      // IndexVectorMapType &recvIndexMap_;


    public:
      PatternBuilder( const CommunicationType& comm,
                      const BlockMapperType& blockMapper,
                      Pattern &pattern )
      : comm_( comm ),
        blockMapper_( blockMapper ),
        myRank_( comm.rank() ),
        mySize_( comm.size() ),
        pattern_( pattern )
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
          // auto& indices = sendIndexMap_[ source ];
          buffer >> indices;
        }

        /*
        for( auto& [rank, future] : recvFutures )
        {
          auto buffer = future.get();
          auto& indices = recvIndexMap_[ rank ];
          buffer >> indices;
        }
        */

        /*
        // wait for all sends to finish
        for( auto& [rank, future] : sendFutures )
        {
          if( future.valid() )
          {
            if( !future.ready() )
              future.wait();
          }
        }
        */

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
      }

      //! return whether we have a fixed size
      bool fixedSize( int dim, int codim ) const
      {
        return false;
      }

      //! read buffer and apply operation
      template< class MessageBuffer, class Entity >
      void gather( MessageBuffer &buffer, const Entity &entity ) const
      {
        // check whether we are a sending entity
        const auto myPartitionType = entity.partitionType();
        const bool send = EntityCommHelper< CommInterface > :: send( myPartitionType );

        // if we send data then send rank and dofs
        if( send )
        {
          // send rank for linkage
          buffer.write( myRank_ );

          const unsigned int numDofs = blockMapper_.numEntityDofs( entity );

          std::vector< GlobalKeyType > indices( numDofs );

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
          GlobalKeyType rank;
          buffer.read( rank );
          assert( (rank >= 0) && (rank < mySize_) );

          // check whether we are a sending entity
          const auto myPartitionType = entity.partitionType();
          const bool receive = EntityCommHelper< CommInterface > :: receive( myPartitionType );

          // insert rank of link into set of links
          // linkStorage_.insert( rank );
          pattern_[rank].rank = rank; // Dune::MPIComm::Link<GlobalKeyType>(rank); // TODO Index typ

          // read indices from stream
          std::vector< GlobalKeyType > indices( dataSize - 1 );
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
            assert( indices[ i ] >= 0 );
            idxMap[ count ] = indices[ i ];
          }
        }
      }

      //! return local dof size to be communicated
      template< class Entity >
      size_t size( const Entity &entity ) const
      {
        const PartitionType myPartitionType = entity.partitionType();
        const bool send = EntityCommHelper< CommInterface > :: send( myPartitionType );
        return (send) ? (
            blockMapper_.numEntityDofs( entity ) // dofs
            +1 // rank
            ) : 0;
      }
    };



    template< class IndexType >
    template< class GridView, class BlockMapper >
    inline void CommunicationPattern< IndexType > :: buildMaps( const GridView& gv, const BlockMapper& blockMapper, const InterfaceType interface )
    {
        // TODO replace by switchCases(const Cases& cases, const Value& value, Branches&& branches)
      typedef typename GridView::CollectiveCommunication CommunicationType;
      if( interface == InteriorBorder_All_Interface )
      {
        PatternBuilder< CommunicationType, BlockMapper, InteriorBorder_All_Interface >
          handle( gv.comm(), blockMapper, pattern_);
        buildMaps( gv, handle );
      }
      else if( interface == InteriorBorder_InteriorBorder_Interface )
      {
        PatternBuilder< CommunicationType, BlockMapper, InteriorBorder_InteriorBorder_Interface >
          handle( gv.comm(), blockMapper, pattern_);
        buildMaps( gv, handle );
      }
      else if( interface == All_All_Interface )
      {
        PatternBuilder< CommunicationType, BlockMapper, All_All_Interface >
          handle( gv.comm(), blockMapper, pattern_);
        buildMaps( gv, handle );
      }
      else
        DUNE_THROW( NotImplemented, "CommunicationPattern for the given interface has not been implemented, yet." );
#ifndef NDEBUG
      // checks that sizes of index maps are equal on sending and receiving proc
      checkConsistency();
#endif
    }


    template< class IndexType >
    template< class GridView, class Comm, class BlockMapper, InterfaceType CI >
    inline void CommunicationPattern< IndexType >
    :: buildMaps( const GridView& gv, PatternBuilder< Comm, BlockMapper, CI > &handle )
    {
      pattern_.clear();

      // make one all to all communication to build up communication pattern
      gv.communicate( handle, All_All_Interface , ForwardCommunication );
    }

    template< class IndexType >
    inline void CommunicationPattern< IndexType > :: checkConsistency()
    {
    }

    template< class IndexType >
    template< class Vector, class Operation >
    inline void CommunicationPattern< IndexType >
    :: exchange( Vector &vector, const Operation& operation ) const
    {
      // create non-blocking communication object
      NonBlockingCommunicationType nbc(*this);

      // perform send operation
      nbc.send(vector);

      // store time for send and receive of data
      exchangeTime_ = nbc.receive(vector, operation);
    }
    //@}

} // namespace Dune

#endif // #ifndef DUNE_FEM_CACHED_COMMUNICATION_MANAGER_HHbuildMap
