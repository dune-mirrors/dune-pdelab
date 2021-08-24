// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACEBASE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACEBASE_HH

#include <optional>

#include <dune/typetree/visitor.hh>
#include <dune/typetree/traversal.hh>

#include <dune/pdelab/common/exceptions.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

#ifndef DOXYGEN

    // forward declaration for friend declaration
    template<typename, typename>
    class GridFunctionSpaceBase;

    template<typename, typename>
    class OrderedGridFunctionSpace;

    namespace impl {

      // helper class with minimal dependencies. Orderings keep a pointer to this structure and populate it
      // during their update procedure.

      template<typename size_type>
      class GridFunctionSpaceOrderingData
      {

        template<typename,typename>
        friend class ::Dune::PDELab::GridFunctionSpaceBase;

        template<class,class>
        friend class ::Dune::PDELab::OrderedGridFunctionSpace;

        GridFunctionSpaceOrderingData()
          : _size(0)
          , _block_count(0)
          , _global_size(0)
          , _max_local_size(0)
          , _is_root_space(true)
          , _initialized(false)
          , _size_available(true)
        {}

        size_type _size;
        size_type _block_count;
        size_type _global_size;
        size_type _max_local_size;
        bool _is_root_space;
        bool _initialized;
        bool _size_available;

      };

      template<class GFS>
      auto gfs_data(GFS& gfs) {
        using SizeType = typename GFS::Traits::SizeType;
        using GFSData = impl::GridFunctionSpaceOrderingData<SizeType>;
        if constexpr (std::is_base_of<GFSData,GFS>{})
          return static_cast<GFSData*>(&gfs);
        else
          return (GFSData*)nullptr;
      }

      //! Checks that every leaf node has the same entity set
      template<class EntitySet>
      struct common_entity_set
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {
        template<typename T, typename TreePath>
        void leaf(T&& t, TreePath treePath) {
          if (not _entity_set)
            _entity_set = t.entitySet();
          else if (*_entity_set != t.entitySet())
            DUNE_THROW(
                GridFunctionSpaceHierarchyError,
                "Use same entity sets for every space that is entity blocked! "
                "A reason for getting this error is creating GridFunctionSpaces with "
                "a grid view in the constructor. To solve this, create an entity set"
                "(e.g. AllEntitySet<GV>) and use one instance to construct all of your GridFunctionSpaces."
);
        }

        std::optional<EntitySet> _entity_set;
      };

    } // namespace impl

#endif // DOXYGEN


    template<class B, class O>
    struct GridFunctionSpaceNodeTraits
    {
      //! vector backend
      [[deprecated]] typedef B BackendType;

      typedef B Backend;

      //! short cut for size type exported by Backend
      typedef typename B::size_type SizeType;

      //! tag describing the ordering.
      /**
       * The tag type may contain additional constants and typedefs to
       * control the behavior of the created ordering.
       */
      typedef O OrderingTag;

    };

    template<typename GFSTraits>
    class GridFunctionSpaceNode
    {
    public:

      typedef GFSTraits Traits;

      template<typename Backend_, typename OrderingTag_>
      GridFunctionSpaceNode(Backend_&& backend, OrderingTag_&& ordering_tag)
        : _backend(std::forward<Backend_>(backend))
        , _ordering_tag(std::forward<OrderingTag_>(ordering_tag))
      {}

      const std::string& name() const
      {
        return _name;
      }

      void name(const std::string& name)
      {
        _name = name;
      }

      typename Traits::Backend& backend()
      {
        return _backend;
      }

      const typename Traits::Backend& backend() const
      {
        return _backend;
      }

      typename Traits::OrderingTag& orderingTag()
      {
        return _ordering_tag;
      }

      const typename Traits::OrderingTag& orderingTag() const
      {
        return _ordering_tag;
      }


    private:

      std::string _name;
      typename Traits::Backend _backend;
      typename Traits::OrderingTag _ordering_tag;
    };


    template<typename GFS, typename GFSTraits>
    class [[deprecated]] GridFunctionSpaceBase
      : public GridFunctionSpaceNode<GFSTraits>
      , public impl::GridFunctionSpaceOrderingData<typename GFSTraits::SizeType>
    {
      using Base = GridFunctionSpaceNode<GFSTraits>;
      GFS& gfs()
      {
        return static_cast<GFS&>(*this);
      }

      const GFS& gfs() const
      {
        return static_cast<const GFS&>(*this);
      }
    public:
      using Base::Base;
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACEBASE_HH
