// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_ORDERING_INTERLEAVEDORDERING_HH
#define DUNE_PDELAB_ORDERING_INTERLEAVEDORDERING_HH

#include <array>
#include <string>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/powernode.hh>

#include <dune/pdelab/common/exceptions.hh>

#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/ordering/orderingbase.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

    namespace interleaved_ordering {

      //! Interface for merging index spaces
      template<typename DI, typename CI, typename Node>
      class Base
        : public OrderingBase<DI,CI>
      {

        typedef OrderingBase<DI,CI> BaseT;

      public:

        typedef typename OrderingBase<DI,CI>::Traits Traits;

        typedef InterleavedOrderingTag OrderingTag;

        static const bool consume_tree_index = true;

        //! Construct ordering object
        /**
         * In general, an ordering object is not properly setup after
         * construction.  This must be done by a seperate call to update()
         * after all the children have been properly set up.
         */
        Base(Node& node, bool container_blocked, const OrderingTag& ordering_tag, typename BaseT::GFSData* gfs_data)
          : BaseT(node,container_blocked,ordering_tag.offsets(),gfs_data,nullptr)
        {
          // This check looks a little weird, but there is always one offset more than
          // there are blocks (the first offsets is 0, and the last one is the "offset
          // beyond the end" to encode the size of the final child).
          if (node.CHILDREN != ordering_tag.offsets().size() - 1)
            DUNE_THROW(OrderingStructureError,
                       "Invalid block structure for InterleavedOrdering: "
                       << node.CHILDREN << " children, but "
                       << (ordering_tag.offsets().size() - 1) << " block sizes.");
        }

        template<typename ItIn, typename ItOut>
        void map_lfs_indices(const ItIn begin, const ItIn end, ItOut out) const
        {
          typedef typename Traits::SizeType size_type;
          if (this->_container_blocked)
            {
              for (ItIn in = begin; in != end; ++in, ++out)
                {
                  size_type child_index = in->treeIndex().back();
                  size_type child_block_offset = this->_child_block_merge_offsets[child_index];
                  size_type child_block_size = this->_child_block_merge_offsets[child_index + 1] - child_block_offset;
                  size_type index = out->back();
                  size_type block_index = index / child_block_size;
                  size_type offset = index % child_block_size;
                  out->back() = child_block_offset + offset;
                  out->push_back(block_index);
                }
            }
          else
            {
              for (ItIn in = begin; in != end; ++in, ++out)
                {
                  size_type child_index = in->treeIndex().back();
                  size_type child_block_offset = this->_child_block_merge_offsets[child_index];
                  size_type child_block_size = this->_child_block_merge_offsets[child_index + 1] - child_block_offset;
                  size_type block_size = this->_child_block_merge_offsets.back();
                  size_type index = out->back();
                  size_type block_index = index / child_block_size;
                  size_type offset = index % child_block_size;
                  out->back() = block_index * block_size + child_block_offset + offset;
                }
            }
        }

        template<typename CIOutIterator, typename DIOutIterator = DummyDOFIndexIterator>
        typename Traits::SizeType
        extract_entity_indices(const typename Traits::DOFIndex::EntityIndex& ei,
                               typename Traits::SizeType child_index,
                               CIOutIterator ci_out, const CIOutIterator ci_end) const
        {
          typedef typename Traits::SizeType size_type;
          if (this->_container_blocked)
            {
              for (; ci_out != ci_end; ++ci_out)
                {
                  size_type child_block_offset = this->_child_block_merge_offsets[child_index];
                  size_type child_block_size = this->_child_block_merge_offsets[child_index + 1] - child_block_offset;
                  size_type index = ci_out->back();
                  size_type block_index = index / child_block_size;
                  size_type offset =index % child_block_size;
                  ci_out->back() = child_block_offset + offset;
                  ci_out->push_back(block_index);
                }
            }
          else
            {
              for (; ci_out != ci_end; ++ci_out)
                {
                  size_type child_block_offset = this->_child_block_merge_offsets[child_index];
                  size_type child_block_size = this->_child_block_merge_offsets[child_index + 1] - child_block_offset;
                  size_type block_size = this->_child_block_merge_offsets.back();
                  size_type index = ci_out->back();
                  size_type block_index = index / child_block_size;
                  size_type offset =index % child_block_size;
                  ci_out->back() = block_index * block_size + child_block_offset + offset;
                }
            }

          // The return value is not used for non-leaf orderings.
          return 0;
        }

      };

    } // namespace interleaved_ordering


    template<typename DI, typename CI, typename Child, std::size_t k>
    class PowerInterleavedOrdering
      : public TypeTree::PowerNode<Child, k>
      , public interleaved_ordering::Base<DI,
                                          CI,
                                          PowerInterleavedOrdering<DI,CI,Child,k>
                                          >
    {
      typedef TypeTree::PowerNode<Child, k> Node;

      typedef interleaved_ordering::Base<DI,
                                         CI,
                                         PowerInterleavedOrdering<DI,CI,Child,k>
                                         > Base;

    public:

      //! Construct ordering object
      /**
       * In general, an ordering object is not properly setup after
       * construction.  This must be done by a seperate call to update() after
       * all the children have been properly set up.
       *
       * \note This constructor must be present for ordering objects not at
       *       the leaf of the tree.
       */
      PowerInterleavedOrdering(bool container_blocked, const InterleavedOrderingTag& ordering_tag, const typename Node::NodeStorage& children, typename Base::GFSData* gfs_data)
        : Node(children)
        , Base(*this,container_blocked,ordering_tag,gfs_data)
      {}

      void update()
      {
        for (std::size_t i = 0; i < k; ++i)
          {
            this->child(i).update();
          }
        Base::update();
      }

      std::string name() const { return "PowerInterleavedOrdering"; }
    };


    template<typename GFS, typename Transformation>
    struct power_gfs_to_interleaved_ordering_descriptor
    {

      static const bool recursive = true;

      template<typename TC>
      struct result
      {

        typedef PowerInterleavedOrdering<
          typename Transformation::DOFIndex,
          typename Transformation::ContainerIndex,
          TC,
          TypeTree::StaticDegree<GFS>::value
          > type;

        typedef std::shared_ptr<type> storage_type;

      };

      template<typename TC>
      static typename result<TC>::type transform(const GFS& gfs, const Transformation& t, const std::array<std::shared_ptr<TC>,TypeTree::StaticDegree<GFS>::value>& children)
      {
        return typename result<TC>::type(gfs.backend().blocked(gfs),gfs.orderingTag(),children,const_cast<GFS*>(&gfs));
      }

      template<typename TC>
      static typename result<TC>::storage_type transform_storage(std::shared_ptr<const GFS> gfs, const Transformation& t, const std::array<std::shared_ptr<TC>,TypeTree::StaticDegree<GFS>::value>& children)
      {
        return std::make_shared<typename result<TC>::type>(gfs->backend().blocked(*gfs),gfs->orderingTag(),children,const_cast<GFS*>(gfs.get()));
      }

    };

    template<typename GFS, typename Transformation>
    power_gfs_to_interleaved_ordering_descriptor<GFS,Transformation>
    register_power_gfs_to_ordering_descriptor(GFS*,Transformation*,InterleavedOrderingTag*);



    template<typename DI, typename CI, typename... Children>
    class CompositeInterleavedOrdering :
      public TypeTree::CompositeNode<Children...>,
      public interleaved_ordering::Base<DI,
                                        CI,
                                        CompositeInterleavedOrdering<
                                          DI,
                                          CI,
                                          Children...
                                          >
                                        >
    {
      typedef TypeTree::CompositeNode<Children...> Node;

      typedef interleaved_ordering::Base<
        DI,
        CI,
        CompositeInterleavedOrdering<
          DI,
          CI,
          Children...
          >
        > Base;

    public:
      //! Construct ordering object
      /**
       * In general, an ordering object is not properly setup after
       * construction.  This must be done by a seperate call to update() after
       * all the children have been properly set up.
       *
       * \note This constructor must be present for ordering objects not at
       *       the leaf of the tree.
       */
      CompositeInterleavedOrdering(bool backend_blocked, const InterleavedOrderingTag& ordering_tag, typename Base::GFSData* gfs_data, std::shared_ptr<Children>... children)
        : Node(children...)
        , Base(*this,backend_blocked,ordering_tag,gfs_data)
      { }

      std::string name() const { return "CompositeInterleavedOrdering"; }

      void update()
      {
        TypeTree::applyToTree(*this,ordering::update_direct_children());
        Base::update();
      }
    };

    template<typename GFS, typename Transformation>
    struct composite_gfs_to_interleaved_ordering_descriptor
    {

      static const bool recursive = true;

      template<typename... TC>
      struct result
      {

        typedef CompositeInterleavedOrdering<
          typename Transformation::DOFIndex,
          typename Transformation::ContainerIndex,
          TC...
          > type;

        typedef std::shared_ptr<type> storage_type;

      };

      template<typename... TC>
      static typename result<TC...>::type transform(const GFS& gfs, const Transformation& t, std::shared_ptr<TC>... children)
      {
        return typename result<TC...>::type(gfs.backend().blocked(gfs),gfs.orderingTag(),const_cast<GFS*>(&gfs),children...);
      }

      template<typename... TC>
      static typename result<TC...>::storage_type transform_storage(std::shared_ptr<const GFS> gfs, const Transformation& t, std::shared_ptr<TC>... children)
      {
        return std::make_shared<typename result<TC...>::type>(gfs->backend().blocked(*gfs),gfs.orderingTag(),const_cast<GFS*>(gfs.get()),children...);
      }

    };

    template<typename GFS, typename Transformation>
    composite_gfs_to_interleaved_ordering_descriptor<GFS,Transformation>
    register_composite_gfs_to_ordering_descriptor(GFS*,Transformation*,InterleavedOrderingTag*);

   //! \} group Ordering
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_INTERLEAVEDORDERING_HH
