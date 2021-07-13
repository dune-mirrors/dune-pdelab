// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTL_BCRSMATRIXBACKEND_HH
#define DUNE_PDELAB_BACKEND_ISTL_BCRSMATRIXBACKEND_HH

#include <dune/pdelab/backend/istl/bcrsmatrix.hh>
#include <dune/pdelab/backend/istl/bcrspattern.hh>
#include <dune/pdelab/backend/istl/patternstatistics.hh>

namespace Dune {
  namespace PDELab {
    namespace ISTL {

      // ********************************************************************************
      // infrastructure for deducing the pattern type from row and column orderings
      // ********************************************************************************

      namespace {

        template<class M, class GFSV, class GFSU, class = void>
        struct build_bcrs_pattern_type
        {
          static_assert(Dune::blockLevel<M>() > 1, "There is not enough blocking to create a pattern in M!");

          // // descend into blocks
          static_assert(Dune::AlwaysFalse<M>{}, "No pattern supported for type M");
          // using BlockPattern = typename build_bcrs_pattern_type<typename M::block_type>::type;
          // using type = DensePattern<BlockPattern>;
        };

        template<class B, class A, class GFSV, class GFSU>
        struct build_bcrs_pattern_type<Dune::BCRSMatrix<B,A>, GFSV, GFSU, std::enable_if_t<(Dune::blockLevel<B>() == 1)>>
        {
          using RowSizeType = typename GFSV::Ordering::Traits::SizeType;
          using RowSizePrefix = typename GFSV::Ordering::Traits::SizePrefix;
          using RowSizeProvider = SizeProviderAdapter<RowSizeType, RowSizePrefix>;
          using ColSizeType = typename GFSU::Ordering::Traits::SizeType;
          using ColSizePrefix = typename GFSU::Ordering::Traits::SizePrefix;
          using ColSizeProvider = SizeProviderAdapter<ColSizeType, ColSizePrefix>;
          using type = SparcityPattern<RowSizeProvider,ColSizeProvider>;
        };

        template<class B, class A, class GFSV, class GFSU>
        struct build_bcrs_pattern_type<Dune::BCRSMatrix<B,A>, GFSV, GFSU, std::enable_if_t<(Dune::blockLevel<B>() > 1)>>
        {
          using BlockPattern = typename build_bcrs_pattern_type<B,GFSV,GFSU>::type;
          using type = BlockSparsityPattern<BlockPattern>;
        };


        // // leaf BCRSMatrix
        template<class Pattern, class B, class A, typename StatsVector>
        void allocate_bcrs_matrix(Pattern& p,
                             typename Pattern::RowSizePrefix prefix_v,
                             typename Pattern::ColSizePrefix prefix_u,
                             Dune::BCRSMatrix<B,A>& c,
                             StatsVector& stats)
        {
          using size_type = typename Pattern::size_type;
          static_assert(Pattern::ColSizeProvider::size_prefix_order == MultiIndexOrder::Outer2Inner);
          static_assert(Pattern::ColSizeProvider::size_prefix_order == MultiIndexOrder::Outer2Inner);
          auto row_size = p.rowSizeProvider().size(prefix_v);
          auto col_size = p.colSizeProvider().size(prefix_u);
          c.setSize(row_size, col_size, 0);
          c.setBuildMode(Dune::BCRSMatrix<B, A>::random);

          auto row_sizes = p.sizes();
          size_type nnz = 0;
          size_type longest_row = 0;

          for (size_type i = 0; i < c.N(); ++i) {
            nnz += row_sizes[i];
            longest_row = std::max(longest_row, row_sizes[i]);
            c.setrowsize(i, row_sizes[i]);
          }
          c.endrowsizes();

          //   stats.push_back(typename
          //   StatsVector::value_type(nnz,longest_row,p.overflowCount(),p.entriesPerRow(),row_size));

          for (size_type i = 0; i < c.N(); ++i){
            for (auto it = p.begin(i) ; it != p.end(i); ++it)
              std::cout << *it << std::endl;
            c.setIndices(i,p.begin(i),p.end(i));
          }

          if constexpr (std::is_void<typename Pattern::SubPattern>{}) {
            // free temporary index storage in pattern before allocating data array in matrix
            p.clear();
            // allocate data array
            c.endindices();
          } else {
            // allocate data array
            c.endindices();
            // add sub patterns on sub matrices
            typename Pattern::RowSizePrefix row_sub_prefix = prefix_v;
            row_sub_prefix.push_back(0);
            typename Pattern::RowSizePrefix col_sub_prefix = prefix_u;
            col_sub_prefix.push_back(0);

            for (auto row_it = c.begin(); row_it != c.end(); ++row_it) {
              auto row = row_it.index();
              row_sub_prefix.back() = row;

              auto pattern_it = p.begin(row);
              // pattern iterator is not necessarely order when has overflow, thus we iterate over
              while (pattern_it != p.end(row)) {
                auto col = *pattern_it;
                col_sub_prefix.back() = col;
                assert(c.exists(row,col));

                // TODO check if pattern matches another already created matrix and copy matrix to reuse pattern
                auto& sub_container = (*row_it)[col];
                allocate_bcrs_matrix(pattern_it.pattern(), row_sub_prefix, col_sub_prefix, sub_container, stats);
                ++pattern_it;
              }
              assert(pattern_it == p.end(row));
            }
            p.clear();
          }
        }

      } // anonymous namespace



      //! Backend using (possibly nested) ISTL BCRSMatrices.
      /**
       * BCRSMatrixBackend is a matrix backend descriptor for ISTL matrices. It expects that
       * both the ansatz and the test function space use ISTL vectors and automatically deduces
       * the correct matrix type from those two vector backends.
       *
       * The backend uses an accelerated pattern construction scheme, which requires the average
       * number of non-zero entries per matrix row as a priori information. In constrast to the
       * older construction scheme, the improved version never requires more memory than the matrix
       * does after pattern construction and runs a lot faster, as long as it is provided with a
       * reasonable estimate for the number of non-zero entries per row.
       *
       */
      template<typename EntriesPerRow = std::size_t>
      struct BCRSMatrixBackend
      {

        //! The size type of the BCRSMatrix.
        typedef std::size_t size_type;

        //! The type of the object holding the statistics generated during pattern construction.
        typedef PatternStatistics<size_type> Statistics;

        //! The type of the pattern object passed to the GridOperator for pattern construction.
        template<typename Matrix, typename GFSV, typename GFSU>
        using Pattern = typename build_bcrs_pattern_type<
          typename Matrix::Container,
          GFSV,
          GFSU
          >::type;

        template<typename VV, typename VU, typename E>
        struct MatrixHelper
        {
          typedef BCRSMatrix<
            typename VV::GridFunctionSpace,
            typename VU::GridFunctionSpace,
            typename build_matrix_type<
              E,
              typename VV::Container,
              typename VU::Container
              >::type,
            Statistics
            > type;
        };

        //! Builds the matrix pattern associated with grid_operator and initializes matrix with it.
        /**
         * \returns  a vector with statistics object for all leaf BCRSMatrices in row-major order.
         */
        template<typename GridOperator, typename Matrix>
        std::vector<Statistics> buildPattern(const GridOperator& grid_operator, Matrix& matrix) const
        {
          using P = Pattern<
            Matrix,
            typename GridOperator::Traits::TestGridFunctionSpace,
            typename GridOperator::Traits::TrialGridFunctionSpace
            >;

         auto row_size_provider =
              SizeProviderAdapter{grid_operator.testGridFunctionSpace().orderingStorage()};
          auto col_size_provider =
              SizeProviderAdapter{grid_operator.trialGridFunctionSpace().orderingStorage()};
          P pattern(row_size_provider, {}, col_size_provider, {},
                                  _entries_per_row);
          grid_operator.fill_pattern(pattern);
          std::vector<Statistics> stats;
          allocate_bcrs_matrix(pattern, {}, {}, Backend::native(matrix), stats);
          return stats;
        }

        //! Constructs a BCRSMatrixBackend.
        /**
         * TODO: Document and flesh out the way this should work for nested matrices (use a nested array as entries_per_row).
         *
         * \param entries_per_row  The average number of nonzero entries per row in matrices created with this backend.
         */
        BCRSMatrixBackend(const EntriesPerRow& entries_per_row)
          : _entries_per_row(entries_per_row)
        {}

      private:

        EntriesPerRow _entries_per_row;

      };

    } // namespace ISTL
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_BCRSMATRIXBACKEND_HH
