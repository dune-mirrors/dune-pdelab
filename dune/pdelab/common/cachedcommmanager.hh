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

//- necessary dune-istl
#include <dune/istl/access.hh>
#include <dune/istl/foreach.hh>

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

    // todo (1): add rank, sort by rank in Link
    // todo (2): make this a lightweight object
    template<typename IndexType>
    using CommunicationPattern = std::map<int, Link<IndexType>>; // rank -> Link

  } // end namespace MPIComm

  /** @addtogroup Communication Communication
      @{
  **/

    template< class Communication, class BlockMapper, InterfaceType CommInterface >
    class PatternBuilder;

    /** \brief Handle "halo"-exchange between neighboring processes
     *
     *
     * CommunicationPattern is a convenience class to build up a map
     * of all dofs of entities to be exchanged during a communication procedure.
     * This speeds up the communication procedure, because no grid traversal is
     * necessary anymore to exchange data. This class is singleton for different
     * discrete function spaces, depending on the BlockMapper.
     */
    template<typename IndexType>
    class ExchangeCommunication
    {
    public:
      //! type of block mapper of discrete function space (may be the same for
      //! different space (i.e. various DG spaces)
        // typedef BlockMapper BlockMapperType;
      using Index = IndexType; // typedef typename BlockMapperType :: GlobalKeyType

    protected:
      typedef Dune::MPIHelper::MPICommunicator MPICommunicatorType;
      typedef Dune::Communication< MPICommunicatorType > CommunicationType;

    protected:

      MPIComm::CommunicationPattern<Index> pattern_;

      // communicator
      std::unique_ptr< CommunicationType > comm_;

      // exchange time
      mutable double exchangeTime_;

      mutable int nonBlockingObjects_ ;

    protected:
      /////////////////////////////////////////////////////////////////
      //  begin NonBlockingCommunication
      /////////////////////////////////////////////////////////////////

      class NonBlockingCommunication
      {
        typedef ExchangeCommunication < IndexType > ExchangeCommunicationType;

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
        NonBlockingCommunication( const ExchangeCommunicationType& dependencyCache )
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
            #warning needs update for inhomogeneous containers
            const size_t bufSize = dependencyCache_.recvBufferSize( source ) * dataSize;
            recvFutures_.insert( std::make_pair( source, comm.irecv( BufferType(comm, bufSize), source, tag_ ) ) );
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
        const ExchangeCommunicationType& dependencyCache_;
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
      ExchangeCommunication()
      : pattern_(),
        comm_(),
        exchangeTime_( 0.0 ),
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
      ExchangeCommunication( const ExchangeCommunication & ) = delete;

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

      void setCommunicationPattern (const Dune::MPIComm::CommunicationPattern<Index> & pattern)
      {
        pattern_ = pattern;
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

    protected:

      template< typename Indices, typename Data >
      std::size_t calcSize(const Indices & indices,
                           const Data &data ) const
      {
        using namespace Dune::ISTL;
        if constexpr (false) // (uniform_access<Data>())
        {
          std::size_t blocksize = 0;
          applyAtIndex(indices[0], data,
            [&](auto && block, const auto & index)
            {
              blocksize = sizeof(block);
            });
          return blocksize * indices.size();
        }
        else
        {
          std::size_t size = 0;
          forEachIndex(indices, data,
            [&](auto && block, const auto & index)
            {
              // serialize the block
              flatVectorForEach(block,
                [&](auto && entry, const auto & index){
                  size += sizeof(entry);
                });
            }
            );
          return size;
        }
      }

      // serialize data of DataImp& vector to object stream
      // --writeBuffer
      template< class Buffer, class Data >
      inline void writeBuffer( const int dest,
                               Buffer &buffer,
                               const Data &data ) const
      {
        auto it = pattern_.find( dest );
        const auto &indices = it->second.sendIndices;

        // static const int blockSize = Data::blockSize;
        const int size = calcSize(indices, data);
        buffer.enlarge( size );

        using namespace Dune::ISTL;
        auto writeEntry = [&](const auto & entry, const auto & index){ buffer << entry; };
        auto writeBlock = [&](const auto & block, const auto & index)
        {
          // serialize the block
          flatVectorForEach(block, writeEntry);
        };

        // for each index
        forEachIndex(indices, data, writeBlock);
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
                       "ExchangeCommunication::readBuffer: Operation needs to be a reference!");

        // get index map of rank belonging to source
        // auto it = recvIndexMap_.find( source );
        // const auto &indexMap = it->second;
        auto it = pattern_.find( source );
        const auto &indices = it->second.recvIndices;

        using namespace Dune::ISTL;
        // apply operation, i.e. COPY, ADD, etc.
        // depending on the block type (scalar or vector),
        // we need to specialize the operation
        forEachIndex(indices, data,
          [&](auto && block, const auto & index)
          {
            auto value = block; // not only correct type, but also correct size
            buffer >> value;
            // must work if branch is true
            #warning should this nesting should be the responsibility of the user?
            operation(value, block);
          }
          );
      }
    };

    template< class IndexType >
    template< class Vector, class Operation >
    inline void ExchangeCommunication< IndexType >
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
