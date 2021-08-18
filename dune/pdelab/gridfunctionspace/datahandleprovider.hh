// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_DATAHANDLEPROVIDER_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_DATAHANDLEPROVIDER_HH

#include <vector>
#include <stack>

#include <dune/common/typetraits.hh>
#include <dune/common/reservedvector.hh>
#include <dune/typetree/visitor.hh>

#include <dune/pdelab/ordering/utility.hh>

namespace Dune {
  namespace PDELab {

    namespace Impl {

      template<class OrderingNode>
      using GridViewOrderingConcept = decltype((std::declval<OrderingNode>().localOrdering(),true));

      template<class Entity, class EntityIndex>
      struct entity_index_visitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<class OrderingNode, class TreePath>
        void pre(const OrderingNode& ordering, const TreePath& tp)
        {
          if constexpr (Std::is_detected<GridViewOrderingConcept,OrderingNode>{}) {
            OrderingNode::Traits::DOFIndexAccessor::GeometryIndex::store(
              _entity_index,
              _entity.type(),
              ordering.entitySet().gridView().indexSet().index(_entity)
            );
          }
        }

        entity_index_visitor(const Entity& entity) : _entity{entity} {}

        EntityIndex _entity_index;
      private:
        const Entity& _entity;
      };

      template<class Entity, class EntityIndex>
      struct get_size_for_entity
        : public entity_index_visitor<Entity, EntityIndex>
      {
        using Base = entity_index_visitor<Entity, EntityIndex>;

        template<typename Ordering, typename TreePath>
        void leaf(const Ordering& ordering, const TreePath& tp)
        {
          _size += ordering.size(this->_entity_index);
        }

        get_size_for_entity(const Entity& entity)
          : Base(entity)
          , _size{0}
        {}

        std::size_t size() const
        {
          return _size;
        }

      private:

        std::size_t _size;
      };


      template<class Entity, class EntityIndex, class OffsetIterator>
      struct get_leaf_offsets_for_entity
        : public entity_index_visitor<Entity, EntityIndex>
      {
        using Base = entity_index_visitor<Entity, EntityIndex>;

        template<typename Ordering, typename TreePath>
        void leaf(const Ordering& ordering, const TreePath& tp)
        {
          *(++_oit) = ordering.size(this->_entity_index);
        }

        get_leaf_offsets_for_entity(const Entity& entity, OffsetIterator oit)
          : Base(entity)
          , _oit(oit)
        {}

        //! Export current position of offset iterator - required for MultiDomain support
        OffsetIterator offsetIterator() const
        {
          return _oit;
        }

      private:

        OffsetIterator _oit;
      };


      template<class Entity, class DOFIndex, class ContainerIndex, std::size_t tree_depth, bool map_dof_indices = false>
      struct indices_for_entity
        : public entity_index_visitor<Entity, typename DOFIndex::EntityIndex>
      {
        using Base = entity_index_visitor<Entity, typename DOFIndex::EntityIndex>;

        typedef std::size_t size_type;
        typedef typename std::vector<ContainerIndex>::iterator CIIterator;
        typedef typename std::conditional<
          map_dof_indices,
          typename std::vector<DOFIndex>::iterator,
          DummyDOFIndexIterator
          >::type DIIterator;


        template<typename Ordering, typename Child, typename TreePath, typename ChildIndex>
        void beforeChild(const Ordering& ordering, const Child& child, const TreePath& tp, const ChildIndex& childIndex)
        {
          _stack.push(std::make_pair(_ci_it,_di_it));
        }

        template<typename Ordering, typename TreePath>
        void leaf(const Ordering& ordering, const TreePath& tp)
        {
          static_assert(TreePath::size()>=1, "Cannot call method 'indices_for_entity::leaf' for empty tree path");
          size_type size = ordering.extract_entity_indices(this->_entity_index,
                                                           tp.back(),
                                                           _ci_it,
                                                           _ci_end,
                                                           _di_it);

          _ci_end += size;
          _ci_it = _ci_end;
          _di_end += size;
          _di_it = _di_end;
        }

        template<typename Ordering, typename Child, typename TreePath, typename ChildIndex>
        void afterChild(const Ordering& ordering, const Child& child, const TreePath& tp, const ChildIndex& childIndex)
        {
          // pop
          ordering.extract_entity_indices(this->_entity_index,
                                          childIndex,
                                          _stack.top().first,
                                          _ci_end);

          if (Ordering::consume_tree_index)
            for (DIIterator it = _stack.top().second;
                 it != _di_end;
                 ++it)
              it->treeIndex().push_back(childIndex);

          _stack.pop();
        }


        indices_for_entity(const Entity& entity,
                           CIIterator ci_begin,
                           DIIterator di_begin = DIIterator())
          : Base(entity)
          , _ci_it(ci_begin)
          , _ci_end(ci_begin)
          , _di_it(di_begin)
          , _di_end(di_begin)
        {}


        // Exposed for multidomain support
        CIIterator ci_end() const
        {
          return _ci_end;
        }

        // Exposed for multidomain support
        DIIterator di_end() const
        {
          return _di_end;
        }

      private:

        CIIterator _ci_it;
        CIIterator _ci_end;
        DIIterator _di_it;
        DIIterator _di_end;

        std::stack<
          std::pair<
            CIIterator,
            DIIterator
            >,
          ReservedVector<
            std::pair<
              CIIterator,
              DIIterator
              >,
            tree_depth
            >
          > _stack;
      };

    } // anonymous namespace


    template<typename GFS>
    class DataHandleProvider
    {

    public:

      typedef std::size_t size_type;

      //------------------------------
      // generic data handle interface
      //------------------------------

      //! returns true if data for this codim should be communicated
      bool dataHandleContains (int codim) const
      {
        return gfs().ordering().contains(codim);
      }

      //! returns true if size per entity of given dim and codim is a constant
      bool dataHandleFixedSize (int codim) const
      {
        return gfs().ordering().fixedSize(codim);
      }

      //! Returns true if the sizes of the leaf orderings in this tree should be sent as part of the communcation.
      /**
       * The MultiDomain extensions require knowledge about the size of the individual
       * orderings, which might belong to separate subdomains. Otherwise it is possible
       * to have size mismatches for entities with codim > 0 if there are protruding edges
       * in the parallel mesh partitioning.
       *
       * By default, this method will always return false. It must be overridden for cases
       * where the data actually needs to be sent.
       *
       * This flag also modifies the behavior of the generic data handles, which will automatically
       * send, receive and process the additional information. Note that if sendLeafSizes() returns
       * true, the underlying DataHandleIF of the grid will always use the data type char to be able
       * to send different types of data, which will automatically be marshalled to / from a byte stream.
       */
      constexpr bool sendLeafSizes() const
      {
        using OrderingTag = typename GFS::Traits::OrderingTag;
        return not (std::is_same_v<OrderingTag, EntityBlockedOrderingTag> or std::is_same_v<OrderingTag, DefaultLeafOrderingTag>);
      }

      /*! how many objects of type DataType have to be sent for a given entity

        Note: Only the sender side needs to know this size.
      */
      template<typename Entity>
      size_type dataHandleSize (const Entity& e) const
      {
        typedef typename GFS::Ordering Ordering;

        typedef typename Ordering::Traits::DOFIndex::EntityIndex EntityIndex;
        EntityIndex ei;

        Ordering::Traits::DOFIndexAccessor::GeometryIndex::store(
          ei,
          e.type(),
          gfs().gridView().indexSet().index(e)
        );

        Impl::get_size_for_entity<Entity,EntityIndex> get_size(e);
        TypeTree::applyToTree(gfs().ordering(),get_size);

        return get_size.size();
      }

      template<typename V, typename EntityIndex>
      void setup_dof_indices(V& v, size_type n, const EntityIndex& ei, std::integral_constant<bool,true>) const
      {
        v.resize(n);
        for (typename V::iterator it = v.begin(),
               endit = v.end();
             it != endit;
             ++it)
          {
            it->treeIndex().clear();
            it->entityIndex() = ei;
          }
      }

      template<typename V, typename EntityIndex>
      void setup_dof_indices(V& v, size_type n, const EntityIndex& ei, std::integral_constant<bool,false>) const
      {}

      template<typename V>
      typename V::iterator dof_indices_begin(V& v, std::integral_constant<bool,true>) const
      {
        return v.begin();
      }

      template<typename V>
      DummyDOFIndexIterator dof_indices_begin(V& v, std::integral_constant<bool,false>) const
      {
        return DummyDOFIndexIterator();
      }

      //! return vector of global indices associated with the given entity
      template<typename Entity, typename ContainerIndex, typename DOFIndex, typename OffsetIterator, bool map_dof_indices>
      void dataHandleIndices (const Entity& e,
                              std::vector<ContainerIndex>& container_indices,
                              std::vector<DOFIndex>& dof_indices,
                              OffsetIterator oit,
                              std::integral_constant<bool,map_dof_indices> map_dof_indices_value
                              ) const
      {
        typedef typename GFS::Ordering Ordering;

        static_assert((std::is_same<ContainerIndex,typename Ordering::Traits::ContainerIndex>::value),
                      "dataHandleContainerIndices() called with invalid ContainerIndex type.");

        typedef typename Ordering::Traits::DOFIndex::EntityIndex EntityIndex;

        Impl::get_leaf_offsets_for_entity<Entity,EntityIndex,OffsetIterator> get_offsets{e,oit};
        TypeTree::applyToTree(gfs().ordering(),get_offsets);
        OffsetIterator end_oit = oit + (TypeTree::TreeInfo<Ordering>::leafCount + 1);

        // convert sizes to offsets - last entry contains total size
        std::partial_sum(oit,end_oit,oit);
        size_type size = *(oit + TypeTree::TreeInfo<Ordering>::leafCount);

        container_indices.resize(size);
        // Clear index state
        for (auto& ci : container_indices)
          ci.clear();

        // setup_dof_indices
        if constexpr (map_dof_indices_value)
        {
          dof_indices.resize(size);
          for (auto& di : dof_indices)
              di.clear();
        }

        Impl::indices_for_entity<
          Entity,
          DOFIndex,
          ContainerIndex,
          TypeTree::TreeInfo<Ordering>::depth,
          map_dof_indices
          > extract_indices(e,container_indices.begin(),dof_indices_begin(dof_indices,map_dof_indices_value));
        TypeTree::applyToTree(gfs().ordering(),extract_indices);

      }

    protected:

      const GFS& gfs() const
      {
        return static_cast<const GFS&>(*this);
      }

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_DATAHANDLEPROVIDER_HH
