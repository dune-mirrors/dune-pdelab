// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTL_BCRSMATRIX_HH
#define DUNE_PDELAB_BACKEND_ISTL_BCRSMATRIX_HH

#include <dune/common/typetraits.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/version.hh>
#if DUNE_VERSION_GTE(ISTL,2,8)
#include <dune/istl/blocklevel.hh>
#endif
#include <dune/common/typetraits.hh>

#include <dune/assembler/backends/istlmatrixbackend.hh>

#include <dune/pdelab/backend/common/spacetraits.hh>
#include <dune/pdelab/backend/common/tags.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/common/aliasedmatrixview.hh>
#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>

namespace Dune {
  namespace PDELab {

    namespace ISTL {

#ifndef DOXYGEN
      // pattern used for grid function spaces; bases set up the matrix via dune-assembler
      template<typename C, typename GFSV, typename GFSU>
      struct wrapper_pattern_type
      {
        typedef typename ISTL::build_pattern_type<C,GFSV,GFSU,typename GFSV::Ordering::ContainerAllocationTag>::type type;
      };

      template<typename C, typename GFSV, typename GFSU>
        requires Dune::PDELab::Impl::DuneFunctionsBasis<GFSV>
      struct wrapper_pattern_type<C,GFSV,GFSU>
      {
        typedef void type;
      };

      // Clear a row of a matrix addressed by an outer-to-inner multi-index
      // (dune-functions bases) and write the diagonal entry if it exists.
      // clearBlock clears the whole block row (as clear_row_block).
      template<typename C, typename MI, typename E>
      void clear_matrix_row_multiindex(C& c, const MI& ri, const E& diagonal_entry, bool clearBlock)
      {
        using Block = typename C::block_type;
        auto& row = c[ri[0]];
        auto diagonal = row.find(ri[0]);
        if constexpr (Dune::IsNumber<Block>::value)
        {
          assert(ri.size() == 1);
          for (auto& entry : row)
            entry = 0.0;
          if (diagonal != row.end())
            *diagonal = diagonal_entry;
        }
        else
        {
          static_assert(Dune::IsNumber<typename Block::field_type>::value and requires { Block::rows; },
            "clear_row is only implemented for BCRSMatrix with scalar or FieldMatrix blocks");
          assert(ri.size() == 2);
          for (auto& block : row)
          {
            if (clearBlock)
              block = 0.0;
            else
              block[ri[1]] = 0.0;
          }
          if (diagonal != row.end())
          {
            if (clearBlock)
              for (int i = 0; i < Block::rows; ++i)
                (*diagonal)[i][i] = diagonal_entry;
            else
              (*diagonal)[ri[1]][ri[1]] = diagonal_entry;
          }
        }
      }
#endif // DOXYGEN

      /** \brief PDELab wrapper around a native ISTL BCRS matrix
       *
       * \tparam GFSV test space (rows): PDELab grid function space or dune-functions basis
       * \tparam GFSU trial space (columns): PDELab grid function space or dune-functions basis
       *
       * Element access is delegated to the dune-assembler ISTL matrix
       * backend, see Impl::SpaceTraits for the index conventions.
       */
      template<typename GFSV, typename GFSU, typename C, typename Stats>
      class BCRSMatrix
        : public Backend::impl::Wrapper<C>
      {

        friend Backend::impl::Wrapper<C>;

        using RowTraits = Dune::PDELab::Impl::SpaceTraits<GFSV>;
        using ColTraits = Dune::PDELab::Impl::SpaceTraits<GFSU>;

      public:

        typedef typename C::field_type ElementType;
        typedef ElementType E;
        typedef C Container;
        typedef typename C::field_type field_type;
        typedef typename C::block_type block_type;
        typedef typename C::size_type size_type;

        typedef GFSU TrialGridFunctionSpace;
        typedef GFSV TestGridFunctionSpace;

        typedef typename RowTraits::ContainerIndex RowIndex;
        typedef typename ColTraits::ContainerIndex ColIndex;

        typedef typename ISTL::wrapper_pattern_type<C,GFSV,GFSU>::type Pattern;

        typedef Stats PatternStatistics;

        using value_type = E;

#ifndef DOXYGEN

        // some trickery to avoid exposing average users to the fact that there might
        // be multiple statistics objects
        typedef typename std::conditional<
#if DUNE_VERSION_LT(DUNE_ISTL,2,8)
          (C::blocklevel > 2),
#else
          (blockLevel<C>() > 2),
#endif
          std::vector<PatternStatistics>,
          PatternStatistics
          >::type StatisticsReturnType;

#endif // DOXYGEN

        template<typename RowCache, typename ColCache>
        using LocalView = UncachedMatrixView<BCRSMatrix,RowCache,ColCache>;

        template<typename RowCache, typename ColCache>
        using ConstLocalView = ConstUncachedMatrixView<const BCRSMatrix,RowCache,ColCache>;

        template<typename RowCache, typename ColCache>
        using AliasedLocalView = AliasedMatrixView<BCRSMatrix,RowCache,ColCache>;

        template<typename RowCache, typename ColCache>
        using ConstAliasedLocalView = ConstAliasedMatrixView<const BCRSMatrix,RowCache,ColCache>;

        template<typename GO>
        explicit BCRSMatrix (const GO& go)
          : _container(std::make_shared<Container>())
        {
          _stats = go.matrixBackend().buildPattern(go,*this);
        }

        /** \brief Construct matrix container using an externally given matrix as storage
         *
         * \tparam GO GridOperator type used to assemble into the matrix
         *
         * \param go GridOperator object used to assemble into the matrix
         * \param container ISTL matrix type that stores the actual data
         *
         * This BCRSMatrix constructor will reassemble the matrix occupation pattern.
         */
        template<typename GO>
        BCRSMatrix (const GO& go, Container& container)
          : _container(Dune::stackobject_to_shared_ptr(container))
        {
          _stats = go.matrixBackend().buildPattern(go,*this);
        }

        template<typename GO>
        BCRSMatrix (const GO& go, const E& e)
          : _container(std::make_shared<Container>())
        {
          _stats = go.matrixBackend().buildPattern(go,*this);
          (*_container) = e;
        }

        //! Creates an BCRSMatrix without allocating an underlying ISTL matrix.
        explicit BCRSMatrix (Backend::unattached_container = Backend::unattached_container())
        {}

        //! Creates an BCRSMatrix with an empty underlying ISTL matrix.
        explicit BCRSMatrix (Backend::attached_container)
          : _container(std::make_shared<Container>())
        {}

        BCRSMatrix(const BCRSMatrix& rhs)
          : _container(std::make_shared<Container>(*(rhs._container)))
        {}

        BCRSMatrix& operator=(const BCRSMatrix& rhs)
        {
          if (this == &rhs)
            return *this;
          _stats.clear();
          if (attached())
            {
              (*_container) = (*(rhs._container));
            }
          else
            {
              _container = std::make_shared<Container>(*(rhs._container));
            }
          return *this;
        }

        //! Returns pattern statistics for all contained BCRSMatrix objects.
        const StatisticsReturnType& patternStatistics() const
        {
#if DUNE_VERSION_LT(DUNE_ISTL,2,8)
          return patternStatistics(std::integral_constant<bool,(C::blocklevel > 2)>());
#else
          return patternStatistics(std::integral_constant<bool,(blockLevel<C>() > 2)>());
#endif
        }

#ifndef DOXYGEN

      private:

        const PatternStatistics& patternStatistics(std::false_type multiple) const
        {
          if (_stats.empty())
            DUNE_THROW(InvalidStateException,"no pattern statistics available");
          return _stats[0];
        }

        const std::vector<PatternStatistics>& patternStatistics(std::true_type multiple) const
        {
          if (_stats.empty())
            DUNE_THROW(InvalidStateException,"no pattern statistics available");
          return _stats;
        }

      public:

#endif

        void detach()
        {
          _container.reset();
          _stats.clear();
        }

        void attach(std::shared_ptr<Container> container)
        {
          _container = container;
        }

        bool attached() const
        {
          return bool(_container);
        }

        const std::shared_ptr<Container>& storage() const
        {
          return _container;
        }

        size_type N() const
        {
          return _container->N();
        }

        size_type M() const
        {
          return _container->M();
        }

        BCRSMatrix& operator= (const E& e)
        {
          (*_container) = e;
          return *this;
        }

        BCRSMatrix& operator*= (const E& e)
        {
          (*_container) *= e;
          return *this;
        }

        E& operator()(const RowIndex& ri, const ColIndex& ci)
        {
          return Dune::Assembler::ISTLMatrixBackend<Container>(*_container)(
            RowTraits::upstreamIndex(ri, blockLevel<C>()), ColTraits::upstreamIndex(ci, blockLevel<C>()));
        }

        const E& operator()(const RowIndex& ri, const ColIndex& ci) const
        {
          const Dune::Assembler::ISTLMatrixBackend<Container> backend(*_container);
          return backend(RowTraits::upstreamIndex(ri, blockLevel<C>()), ColTraits::upstreamIndex(ci, blockLevel<C>()));
        }

      private:

        const Container& native() const
        {
          return *_container;
        }

        Container& native()
        {
          return *_container;
        }

      public:

        template<typename RowCache, typename ColCache>
        value_type* data(const RowCache& row_cache, const ColCache& col_cache)
        {
          return &((*this)(row_cache.containerIndex(0),col_cache.containerIndex(0)));
        }

        template<typename RowCache, typename ColCache>
        const value_type* data(const RowCache& row_cache, const ColCache& col_cache) const
        {
          return &((*this)(row_cache.containerIndex(0),col_cache.containerIndex(0)));
        }

        void flush()
        {}

        void finalize()
        {}

        void clear_row(const RowIndex& ri, const E& diagonal_entry)
        {
          if constexpr (Dune::PDELab::Impl::DuneFunctionsBasis<GFSV>)
            ISTL::clear_matrix_row_multiindex(*_container,ri,diagonal_entry,false);
          else
          {
            ISTL::clear_matrix_row(ISTL::container_tag(*_container),*_container,ri,ri.size()-1);
            ISTL::write_matrix_element_if_exists(diagonal_entry,ISTL::container_tag(*_container),*_container,ri,ri,ri.size()-1,ri.size()-1);
          }
        }

        void clear_row_block(const RowIndex& ri, const E& diagonal_entry)
        {
          if constexpr (Dune::PDELab::Impl::DuneFunctionsBasis<GFSV>)
            ISTL::clear_matrix_row_multiindex(*_container,ri,diagonal_entry,true);
          else
          {
            ISTL::clear_matrix_row_block(ISTL::container_tag(*_container),*_container,ri,ri.size()-1);
            ISTL::write_matrix_element_if_exists_to_block(diagonal_entry,ISTL::container_tag(*_container),*_container,ri,ri,ri.size()-1,ri.size()-1);
          }
        }

      private:

        std::shared_ptr<Container> _container;
        std::vector<PatternStatistics> _stats;

      };

    } // namespace ISTL

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_BCRSMATRIX_HH
