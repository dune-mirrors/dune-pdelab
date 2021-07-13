// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTL_BCRSPATTERN_HH
#define DUNE_PDELAB_BACKEND_ISTL_BCRSPATTERN_HH

#include <utility>
#include <vector>
#include <algorithm>
#include <set>
#include <map>

#include <dune/common/iteratorfacades.hh>

#include <dune/pdelab/backend/common/tags.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>

namespace Dune {
  namespace PDELab {
    namespace ISTL {

      //! Pattern builder for generic BCRS-like sparse matrices.
      /**
       * BCRSPattern is a pattern builder for unstructured sparse matrices
       * for operators mapping from a vector that conforms to RowSizeProvider to a vector
       * that conforms to ColSizeProvider.
       *
       * BCRSPattern has much better runtime performance and requires far less memory
       * than the older pattern constructon method in PDELab. By letting the user specify
       * the average number of nonzeroes per row, it is possible to use a more efficient
       * array-based storage scheme for the majority of the pattern entries, only using
       * expensive map-like lookups for those entries that exceed that average.
       *
       * BCRSPattern requires a recent version of the BCRSMatrix with support for row-wise
       * setting of column indices and split allocation of column index and data arrays.
       *
       * Note that unlike the implicit construction mode of the BCRSMatrix itself, this
       * pattern builder will neither throw an exception if the number of nonzeroes was set
       * too low nor retain excess memory if it was set too high after the pattern construction
       * is complete. Performance will degrade if the user-provided estimate is too far away
       * from the real value.
       */
      template<class RowSizeProvider_, class ColSizeProvider_>
      class SparcityPattern
      {

      public:
        //! SparcityPattern cannot contain nested subpatterns. This entry is only required for TMP purposes.
        using SubPattern = void;

        using RowSizeProvider = RowSizeProvider_;
        using ColSizeProvider = ColSizeProvider_;

        static_assert(RowSizeProvider::size_prefix_order == MultiIndexOrder::Outer2Inner);
        static_assert(ColSizeProvider::size_prefix_order == MultiIndexOrder::Outer2Inner);

        using RowSizePrefix = typename RowSizeProvider::SizePrefix;
        using ColSizePrefix = typename ColSizeProvider::SizePrefix;

        //! size type used by SparcityPattern.
        using size_type = std::size_t;


      private:

        //! Marker value indicating an empty array entry.
        static const size_type empty = ~static_cast<size_type>(0);

        typedef typename std::vector<size_type>::iterator IndicesIterator;
        typedef typename std::set<std::pair<size_type,size_type> >::iterator OverflowIterator;

        typedef typename std::vector<size_type>::const_iterator ConstIndicesIterator;
        typedef typename std::set<std::pair<size_type,size_type> >::const_iterator ConstOverflowIterator;

      public:

        //! Add a link between the row indicated by ri and the column indicated by ci.
        template<typename RI, typename CI>
        void add_link(const RI& ri, const CI& ci)
        {
          // extract block indices for current level
          size_type i = ri.back();
          size_type j = ci.back();

          IndicesIterator start = _indices.begin();
          IndicesIterator begin = start + _entries_per_row*i;
          IndicesIterator end = start + _entries_per_row*(i+1);

          // Looking up column indices requires a special comparison iterator,
          // as we want to either return the position of the actual index if it
          // has already been inserted or of the first empty matrix entry that we
          // can use to store the index (Only use for sequential search algorithms).
          auto padded_col_criterion = [j](auto k){return k == empty || k == j;};

          // Does the entry (i,j) already exist?
          IndicesIterator it = std::find_if(begin,end,padded_col_criterion);
          if (it != end)
            {
              // Yes, just write out j. This does the right thing regardless of whether
              // it points at the correct column value or at an empty entry.
              *it = j;
            }
          else
            {
              // The row is already full -> spill into map
              _overflow.insert(std::make_pair(i,j));
            }
        }

#ifndef DOXYGEN

        // implementation detail for nested patterns
        template<typename RI, typename CI>
        void recursive_add_entry(const RI& ri, const CI& ci)
        {
          add_link(ri,ci);
        }

#endif //DOXYGEN

        //! Stream the sizes of all rows into the output iterator rit.
        template<typename I>
        void sizes(I rit) const
        {
          ConstIndicesIterator it = _indices.begin();
          ConstIndicesIterator end = _indices.begin() + _entries_per_row;
          ConstOverflowIterator oit = _overflow.begin();
          ConstOverflowIterator oend = _overflow.end();
          for (size_type i = 0; i < _rows; ++i, ++rit, end+=_entries_per_row)
            {
              size_type s = 0;
              // count non-empty column entries, break when first empty one is found.
              for (; it != end; ++it)
                {
                  if (*it == empty)
                    break;
                  ++s;
                }
              it = end;
              // add overflow entries
              for (; oit != oend && oit->first == i; ++oit)
                ++s;
              *rit = s;
            }
        }

        //! Returns a vector with the size of all rows in the pattern.
        std::vector<size_type> sizes() const
        {
          std::vector<size_type> r(_rows);
          sizes(r.begin());
          return r;
        }

        //! Iterator over all column indices for a given row, unique but in arbitrary order.
        struct iterator
          : public ForwardIteratorFacade<iterator, const size_type>
        {

#ifndef DOXYGEN

          const size_type& dereference() const
          {
            if (_in_overflow)
              return _oit->second;
            else
              return *_it;
          }

          void increment()
          {
            if (_in_overflow)
              {
                if (++_oit == _oend or _oit->first != _row)
                  {
                    // we have exhausted the row, invalidate iterator
                    _at_end = true;
                  }
              }
            else
              {
                if (_it != _end)
                  ++_it;
                if (_it == _end || *_it == empty)
                  {
                    _in_overflow = true;
                    // we have exhausted the row, invalidate iterator
                    if (_oit == _oend || _oit->first > _row)
                      _at_end = true;
                  }
              }
          }

          bool equals(const iterator& other) const
          {
            if (_row != other._row)
              return false;
            if (_at_end || other._at_end)
              return _at_end && other._at_end;
            if (_in_overflow)
              return _oit == other._oit;
            else
              return _it == other._it;
          }

          iterator(const SparcityPattern& p, size_type row, bool at_end)
            : _row(row)
            , _in_overflow(false)
            , _at_end(at_end)
            , _it(p._indices.begin() + row * p._entries_per_row)
            , _end(p._indices.begin() + (row+1) * p._entries_per_row)
            , _oit(p._overflow.lower_bound(std::make_pair(row,0)))
            , _oend(p._overflow.end())
          {
            // catch corner case with completely empty row
            if ((!_at_end) && (_it == _end || *_it == empty))
              {
                _in_overflow = true;
                _at_end = _oit == _oend || _oit->first != _row;
              }
          }

          size_type _row;
          bool _in_overflow;
          bool _at_end;
          typename std::vector<size_type>::const_iterator _it;
          typename std::vector<size_type>::const_iterator _end;
          typename std::set<std::pair<size_type,size_type> >::const_iterator _oit;
          const typename std::set<std::pair<size_type,size_type> >::const_iterator _oend;

#endif // DOXYGEN

        };

        //! Returns an iterator to the first column index of row i.
        iterator begin(size_type i) const
        {
          return iterator(*this,i,false);
        }

        //! Returns an iterator past the last column index of row i.
        iterator end(size_type i) const
        {
          return iterator(*this,i,true);
        }

        //! Constructs a SparcityPattern for the given pair of orderings and reserves space for the provided average number of entries per row.
        /**
         * \param row_size_provider    Ordering describing the row structure
         * \param col_size_provider    Ordering describing the column structure
         * \param entries_per_row An estimate of the average number of entries per row.
         */
        SparcityPattern(const RowSizeProvider row_size_provider,
                    typename RowSizeProvider::SizePrefix row_prefix,
                    const ColSizeProvider col_size_provider,
                    typename ColSizeProvider::SizePrefix col_prefix,
                    size_type entries_per_row)
          : _row_size_provider(row_size_provider)
          , _col_size_provider(col_size_provider)
          , _rows(_row_size_provider.size(row_prefix))
          , _entries_per_row(entries_per_row)
          , _indices(_rows*entries_per_row,size_type(empty))
        {}

        //! Discard all internal data.
        /**
         * The purpose of this method is to release all internal memory before calling
         * BCRSMatrix::endindices(). That way, the matrix creation process never consumes
         * substantially more memory as required by the matrix after construction, as the
         * second copy of the column indices is about as large as the data array.
         */
        void clear()
        {
          _indices = std::vector<size_type>();
          _overflow = std::set<std::pair<size_type,size_type> >();
        }

        size_type entriesPerRow() const
        {
          return _entries_per_row;
        }

        size_type overflowCount() const
        {
          return _overflow.size();
        }

        const RowSizeProvider& rowSizeProvider() const
        {
          return _row_size_provider;
        }

        const ColSizeProvider& colSizeProvider() const
        {
          return _col_size_provider;
        }

      private:

        const RowSizeProvider _row_size_provider;
        const ColSizeProvider _col_size_provider;

        const size_type _rows, _entries_per_row;

        std::vector<size_type> _indices;
        std::set<std::pair<size_type,size_type> > _overflow;

      };


      //! Pattern builder for nested hierarchies of generic BCRS-like sparse matrices.
      /**
       * NestedPattern contains a dense set of subpatterns for each matrix block. Those
       * blocks can be nested (i.e. be NestedPatterns again) or BCRSPattern instances.
       */

      template<class Pattern>
      class BlockSparsityPattern
      {
      public:
        //! The pattern type used for each block.
        using SubPattern = Pattern;

        using RowSizeProvider = typename Pattern::RowSizeProvider;
        using ColSizeProvider = typename Pattern::ColSizeProvider;

        static_assert(RowSizeProvider::size_prefix_order == MultiIndexOrder::Outer2Inner);
        static_assert(ColSizeProvider::size_prefix_order == MultiIndexOrder::Outer2Inner);

        using RowSizePrefix = typename RowSizeProvider::SizePrefix;
        using ColSizePrefix = typename ColSizeProvider::SizePrefix;


        //! size type used by BlockSparsityPattern.
        typedef typename SubPattern::size_type size_type;


        //! Marker value indicating an empty array entry.
        static constexpr size_type empty = ~static_cast<size_type>(0);

        //! Add a link between the row indicated by ri and the column indicated by ci.
        /**
         * This method just forwards the call to the relevant block as indicated by the
         * tail members of ri and ci.
         */
        // template<typename RI, typename CI>
        // void add_link(const RI& ri, const CI& ci)
        // {
        //   // recursive_add_entry(ri.view(),ci.view());
        // }

        //! Add a link between the row indicated by ri and the column indicated by ci.
        template<class RowIndex, class ColIndex>
        void add_link(RowIndex row_index, ColIndex col_index)
        {
          // in case of empty indices, no link is needed ???
          if(row_index.size() == 0 or col_index.size() == 0)
            return;

          // extract block indices for current level
          size_type i = row_index.back();
          size_type j = col_index.back();

          auto start = _indices.begin();
          auto begin = start + _entries_per_row*i;
          auto end = start + _entries_per_row*(i+1);

          // Looking up column indices requires a special comparison iterator,
          // as we want to either return the position of the actual index if it
          // has already been inserted or of the first empty matrix entry that we
          // can use to store the index (Only use for sequential search algorithms).
          auto padded_col_criterion = [j](auto k){return k.first == empty || k.first == j;};

          // Does the entry (i,j) already exist?
          auto it = std::find_if(begin,end,padded_col_criterion);
          std::shared_ptr<SubPattern> sub_pattern;
          if (it != end) {
            RowSizePrefix sub_row_prefix = _row_prefix;
            ColSizePrefix sub_col_prefix = _col_prefix;
            sub_row_prefix.push_back(i);
            sub_col_prefix.push_back(j);
            if (not it->second)
              it->second = std::make_shared<SubPattern>(
                  _row_size_provider, sub_row_prefix, _col_size_provider,
                  sub_col_prefix, _entries_per_row);
            // Yes, just write out j. This does the right thing regardless of
            // whether it points at the correct column value or at an empty
            // entry.
            it->first = j;
            sub_pattern = it->second;
          } else {
            DUNE_THROW(NotImplemented, "Foo");
            // The row is already full -> spill into map
            // _overflow.insert(std::make_pair(i,j));
          }
          row_index.pop_back();
          col_index.pop_back();
          sub_pattern->add_link(row_index,col_index);
        }

        //! Stream the sizes of all rows into the output iterator rit.
        template<typename OutIt>
        void sizes(OutIt out) const
        {
          auto idx_it = std::begin(_indices);
          auto idx_end = std::begin(_indices) + _entries_per_row;
          auto ovf_it = std::begin(_overflow);
          auto ovf_end = std::end(_overflow);
          for (size_type i = 0; i < _rows; ++i) {
            size_type row_size = 0;
            // count non-empty column entries
            while (idx_it != idx_end and idx_it->first != empty)
              ++idx_it, ++row_size;
            // add overflow entries
            while (ovf_it != ovf_end and ovf_it->first.first == i)
              ++ovf_it, ++row_size;
            // write row size into out iterator and advance to next row
            *out = row_size;
            ++out;
            // advance index range to next row
            idx_it = idx_end;
            idx_end += _entries_per_row;
          }
        }


        //! Returns a vector with the size of all rows in the pattern.
        std::vector<size_type> sizes() const
        {
          std::vector<size_type> row_sizes(_rows);
          this->sizes(std::begin(row_sizes));
          return row_sizes;
        }


        //! Iterator over all column indices for a given row, unique but in arbitrary order.
        struct iterator
          : public ForwardIteratorFacade<iterator, const size_type>
        {

#ifndef DOXYGEN

          const size_type& dereference() const
          {
            if (_in_overflow)
              return _oit->first.second;
            else
              return _it->first;
          }

          SubPattern& pattern()
          {
            if (_in_overflow)
              return *_oit->second;
            else
              return *_it->second;
          }

          void increment()
          {
            if (_in_overflow)
              {
                if (++_oit == _oend || _oit->first.first != _row)
                  {
                    // we have exhausted the row, invalidate iterator
                    _at_end = true;
                  }
              }
            else
              {
                if (_it != _end)
                  ++_it;
                if (_it == _end || _it->first == empty)
                  {
                    _in_overflow = true;
                    // we have exhausted the row, invalidate iterator
                    if (_oit == _oend || _oit->first.first > _row)
                      _at_end = true;
                  }
              }
          }

          bool equals(const iterator& other) const
          {
            if (_row != other._row)
              return false;
            if (_at_end || other._at_end)
              return _at_end && other._at_end;
            if (_in_overflow)
              return _oit == other._oit;
            else
              return _it == other._it;
          }

          iterator(const BlockSparsityPattern& p, size_type row, bool at_end)
            : _row(row)
            , _in_overflow(false)
            , _at_end(at_end)
            , _it(p._indices.begin() + row * p._entries_per_row)
            , _end(p._indices.begin() + (row+1) * p._entries_per_row)
            , _oit(p._overflow.lower_bound(std::make_pair(row,0)))
            , _oend(p._overflow.end())
          {
            // catch corner case with completely empty row
            if ((!_at_end) && (_it == _end || _it->first == empty))
              {
                _in_overflow = true;
                _at_end = _oit == _oend || _oit->first.first != _row;
              }
          }

          size_type _row;
          bool _in_overflow;
          bool _at_end;
          typename  std::vector<std::pair<size_type,std::shared_ptr<SubPattern> > >::const_iterator _it;
          typename  std::vector<std::pair<size_type,std::shared_ptr<SubPattern> > >::const_iterator _end;
          typename std::map<std::pair<size_type,size_type>,std::shared_ptr<SubPattern> >::const_iterator _oit;
          const typename std::map<std::pair<size_type,size_type>,std::shared_ptr<SubPattern> >::const_iterator _oend;

#endif // DOXYGEN

        };

        //! Returns an iterator to the first column index of row i.
        iterator begin(size_type i) const
        {
          return iterator(*this,i,false);
        }

        //! Returns an iterator past the last column index of row i.
        iterator end(size_type i) const
        {
          return iterator(*this,i,true);
        }

#ifndef DOXYGEN

        // template<typename RI, typename CI>
        // void recursive_add_entry(const RI& ri, const CI& ci)
        // {
        //   _sub_patterns[ri.back() * _cols + ci.back()].recursive_add_entry(ri.back_popped(),ci.back_popped());
        // }

#endif // DOXYGEN

        BlockSparsityPattern(const RowSizeProvider& row_size_provider,
                      RowSizePrefix row_prefix,
                      const ColSizeProvider& col_size_provider,
                      ColSizePrefix col_prefix,
                      const std::size_t& entries_per_row)
          : _row_size_provider(row_size_provider)
          , _col_size_provider(col_size_provider)
          , _row_prefix(row_prefix)
          , _col_prefix(col_prefix)
          , _rows(row_size_provider.size(_row_prefix))
          , _cols(col_size_provider.size(_col_prefix))
          , _entries_per_row(entries_per_row)
          , _indices(_rows*entries_per_row,std::make_pair(empty,nullptr))
        {

        }

        // NestedPattern(const RowSizeProvider& row_size_provider,
        //               typename RowSizeProvider::SizePrefix row_prefix,
        //               const ColSizeProvider& col_size_provider,
        //               typename ColSizeProvider::SizePrefix col_prefix,
        //               const size_type& entries_per_row)
        //   : _rows(row_size_provider.size(row_prefix))
        //   , _cols(col_size_provider.size(col_prefix))
        //   , _indices(_rows*entries_per_row,std::make_pair(empty,nullptr))
        // {

        // }

        // //! Returns the subpattern associated with block (i,j).
        // SubPattern& subPattern(size_type i, size_type j)
        // {
        //   return _sub_patterns[i * _cols + j];
        // }

        void clear()
        {
          _indices.clear();
          _overflow.clear();
        }

        const RowSizeProvider& rowSizeProvider() const
        {
          return _row_size_provider;
        }

        const ColSizeProvider& colSizeProvider() const
        {
          return _col_size_provider;
        }

      private:

        const RowSizeProvider _row_size_provider;
        const ColSizeProvider _col_size_provider;

        const RowSizePrefix _row_prefix;
        const ColSizePrefix _col_prefix;
        const size_type _rows, _cols, _entries_per_row;

        std::vector<std::pair<size_type,std::shared_ptr<SubPattern> > > _indices;
        std::map<std::pair<size_type,size_type>,std::shared_ptr<SubPattern> > _overflow;
      };


    } // namespace ISTL
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_BCRSPATTERN_HH
