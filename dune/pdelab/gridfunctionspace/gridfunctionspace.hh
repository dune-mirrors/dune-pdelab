// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACE_HH

#include <cstddef>
#include <map>
#include <ostream>
#include <set>
#include <vector>
#include <memory>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/common/version.hh>
#include <dune/typetree/typetree.hh>

// This alias should be removed after a PDELab 2.7 release.
#if DUNE_VERSION_LT_REV(DUNE_TYPETREE,2,7,1)
namespace Dune {
  namespace TypeTree {
    template<std::size_t... i>
    using StaticTreePath = TreePath<i...>;
  }
}
#endif

#include <dune/pdelab/common/partitionviewentityset.hh>
#include <dune/pdelab/backend/interface.hh>

// we just want the descriptors here, so we temporarily switch off the warning for
// directly including ISTL backend headers
#define _DUNE_PDELAB_SUPPRESS_ISTL_HH_WARNING
#include <dune/pdelab/backend/istl/descriptors.hh>
#undef _DUNE_PDELAB_SUPPRESS_ISTL_HH_WARNING

#include <dune/pdelab/constraints/noconstraints.hh>
#include <dune/pdelab/gridfunctionspace/compositegridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspacebase.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/orderedgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/ordering/gridviewordering.hh>
#include <dune/pdelab/ordering/lexicographicordering.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

#ifndef DOXYGEN

    namespace impl {

      // Helper structs to avoid compilation failures in the
      // backwards compatibility mode where users stuff a
      // GridView into a GridFunctionSpace.
      // In that case, we cannot extract the GridView type from
      // the GridView itself, so we use a std::conditional in the
      // Traits class to pick either one of the following structs
      // and then use the correct class to do the lookup.

      struct _lazy_identity
      {
        template<typename T>
        struct evaluate
        {
          using type = T;
        };
      };

      struct _lazy_extract_gridview
      {
        template<typename T>
        struct evaluate
        {
          using type = typename T::GridView;
        };
      };

      // Returns a GridView, regardless of whether GV_or_ES is a GridView or an EntitySet
      template<typename GV_or_ES>
      using GridView = typename std::conditional<
        isEntitySet<GV_or_ES>::value,
        impl::_lazy_extract_gridview,
        impl::_lazy_identity
        >::type::template evaluate<GV_or_ES>::type;


      // Returns an EntitySet, regardless of whether GV_or_ES is a GridView or an EntitySet
      template<typename GV_or_ES>
      using EntitySet = typename std::conditional<
        isEntitySet<GV_or_ES>::value,
        GV_or_ES,
        AllEntitySet<GV_or_ES>
        >::type;

    }

#endif // DOXYGEN

    //=======================================
    // grid function space : single component case
    //=======================================

    //! collect types exported by a leaf grid function space
    /**
     * This is based on a global FiniteElementMap
     */
    template<typename G, typename L, typename C, typename B, typename O>
    struct GridFunctionSpaceLeafTraits : public GridFunctionSpaceNodeTraits<B,O>
    {
      //! True if this grid function space is composed of others.
      [[deprecated]] static const bool isComposite = false;

      //! the grid view where grid function is defined upon
      using GridView = impl::GridView<G>;

      //! the entity set of this function space.
      using EntitySet = impl::EntitySet<G>;

      using GridViewType = GridView;

      typedef L FiniteElementMapType;

      //! finite element map
      typedef L FiniteElementMap;

      typedef typename L::Traits::FiniteElementType FiniteElementType;

      //! finite element
      typedef typename L::Traits::FiniteElementType FiniteElement;

      //! type representing constraints
      typedef C /*[[deprecated]]*/ ConstraintsType;

      typedef C Constraints;
    };

    /** \brief An unordered grid function space.
     *  \details Representation of an function space in a grid.
     *
     *  - This class should be prefered over an (ordered) GridFunctionSpace when
     *    building a tree of spaces.
     *  - This class is not able to provide an \ref Ordering by itself. Thus,
     *    it's not suitable for assembling local operators. If that's necessary,
     *    use OrderedGridFunctionSpace or GridFunctionSpace.
     *
     *  \tparam ES   Entity Set type (See PartitionViewEntitySet). Respresents
     *               the sub set of entities where the finite element map has
     *               support.
     *  \tparam FEM  Type implementing FiniteElementMapInterface. A map from
     *               entity to local finite element.
     *  \tparam CE   Type for constraints assembler
     *  \tparam B    Backend type
     *  \tparam O    Ordering tag
     *
     *  \see Ordering
     *  \see GridFunctionSpace
     *  \see PartitionViewEntitySet
     */
    template<typename ES, typename FEM, typename CE=NoConstraints,
             typename B=ISTL::VectorBackend<>, typename O=DefaultLeafOrderingTag>
    class UnorderedGridFunctionSpace
      : public TypeTree::LeafNode
      , public GridFunctionSpaceNode<GridFunctionSpaceLeafTraits<ES,FEM,CE,B,O>>
      , public GridFunctionOutputParameters
    {
      static_assert(isEntitySet<ES>{});

      template<typename,typename>
      friend class GridFunctionSpaceBase;

    public:
      //! export Traits class
      typedef GridFunctionSpaceLeafTraits<ES,FEM,CE,B,O> Traits;

    private:

      typedef GridFunctionSpaceNode<Traits> BaseT;

    public:

      typedef LeafGridFunctionSpaceTag ImplementationTag;


      // ****************************************************************************************************
      // Construct from EntitySet
      // ****************************************************************************************************

      /**
       * @copybrief UnorderedGridFunctionSpace
       * @warning The entity set internals will be modified according to the
       *          finite element map used codimensions
       *
       * @param entitySet     Copy of an entity set
       * @param fem           Finite Element Map pointer
       * @param ce            Constraints Assembler pointer
       * @param backend       Vector backend
       * @param ordering_tag  Ordering tag
       */
      UnorderedGridFunctionSpace(typename Traits::EntitySet entitySet,
                        const std::shared_ptr<const FEM> &fem,
                        const std::shared_ptr<const CE> &ce,
                        const B &backend = B(),
                        const O &ordering_tag = O())
          : BaseT(backend, ordering_tag)
          , pfem(fem)
          , _pce(ce)
          , _entity_set{entitySet}
      {}

      /**
       * @copybrief UnorderedGridFunctionSpace
       * @warning The entity set internals will be modified according to the
       *          finite element map used codimensions
       *
       * @param entitySet     Copy of an entity set
       * @param fem           Finite Element Map pointer
       * @param backend       Vector backend
       * @param ordering_tag  Ordering tag
       */
      UnorderedGridFunctionSpace(typename Traits::EntitySet entitySet,
                        const std::shared_ptr<const FEM> &fem,
                        const B &backend = B(),
                        const O &ordering_tag = O())
          : BaseT(backend, ordering_tag)
          , pfem(fem)
          , _pce(std::make_shared<CE>())
          , _entity_set{entitySet}
      {}

      //! get finite element map
      const FEM& finiteElementMap () const
      {
        return *pfem;
      }

      //! get finite element map
      std::shared_ptr<const FEM> finiteElementMapStorage () const
      {
        return pfem;
      }

      //! return constraints engine
      const typename Traits::ConstraintsType& constraints () const
      {
        return *_pce;
      }

      //! return storage of constraints engine
      std::shared_ptr<const CE> constraintsStorage () const
      {
        return _pce;
      }

      //! get grid view
      const typename Traits::GridView& gridView () const
      {
        return _entity_set.gridView();
      }

      //! get entity set
      const typename Traits::EntitySet& entitySet () const
      {
        return _entity_set;
      }

      //! get entity set
      typename Traits::EntitySet& entitySet ()
      {
        return _entity_set;
      }

  private:

    std::shared_ptr<FEM const> pfem;
    std::shared_ptr<CE const> _pce;
    typename Traits::EntitySet _entity_set;
  };


    /** \brief An ordered grid function space.
     *  \details Representation of an function space in a grid.
     *
     *  - Constructor isthe same as for UnorderedGridFunctionSpace.
     *  - This class should be prefered over an UnorderedGridFunctionSpace when
     *    only one space is required (i.e. no space trees).
     *  - This class is able to provide an \ref Ordering by itself. Thus,
     *    it is suitable for assembling local operators.
     *
     *  \tparam ES   Entity Set type (See PartitionViewEntitySet). Respresents
     *               the sub set of entities where the finite element map has
     *               support.
     *  \tparam FEM  Type implementing FiniteElementMapInterface. A map from
     *               entity to local finite element.
     *  \tparam CE   Type for constraints assembler
     *  \tparam B    Backend type
     *  \tparam O    Ordering tag
     *
     *  \see Ordering
     *  \see GridFunctionSpace
     *  \see PartitionViewEntitySet
     */
    template<typename ES, typename FEM, typename CE=NoConstraints,
             typename B=ISTL::VectorBackend<>, typename O=DefaultLeafOrderingTag>
    class GridFunctionSpace
      : public OrderedGridFunctionSpace<UnorderedGridFunctionSpace<impl::EntitySet<ES>,FEM,CE,B,O>,impl::EntitySet<ES>>
    {
      using UnorderedGFS = UnorderedGridFunctionSpace<impl::EntitySet<ES>,FEM,CE,B,O>;
      using OrderedGFS = OrderedGridFunctionSpace<UnorderedGFS,impl::EntitySet<ES>>;

    public:

      using Traits = typename OrderedGFS::Traits;

      typedef typename ES::Traits::template Codim<0>::Entity Element;
      typedef typename ES::Traits::template Codim<0>::Iterator ElementIterator;

      using OrderingTag [[deprecated]] = O;

      GridFunctionSpace(typename Traits::EntitySet entitySet,
                        const std::shared_ptr<const FEM> &fem,
                        const std::shared_ptr<const CE> &ce,
                        const B &backend = B(),
                        const O &ordering_tag = O())
        : OrderedGFS(std::in_place, std::move(entitySet), fem, ce, backend, ordering_tag)
      {}

      GridFunctionSpace(typename Traits::EntitySet entitySet,
                        const std::shared_ptr<const FEM> &fem,
                        const B &backend = B(),
                        const O &ordering_tag = O())
        : OrderedGFS(std::in_place, std::move(entitySet), fem, backend, ordering_tag)
      {}

      // ****************************************************************************************************
      // Construct from GridView
      // ****************************************************************************************************

       //! constructor
      GridFunctionSpace(const typename Traits::GridView& gridview, const FEM& fem, const CE& ce, const B& backend = B(), const O& ordering_tag = O())
#if DUNE_PDELAB_WARN_ON_GRIDVIEW_BASED_GFS
        DUNE_DEPRECATED_MSG("GridFunctionSpaces now internally use an EntitySet instead of a GridView, please replace the template parameter and the first constructor parameter by an EntitySet").
#endif
        : OrderedGFS(std::in_place, typename Traits::EntitySet{gridview}, stackobject_to_shared_ptr(fem), stackobject_to_shared_ptr(ce), backend, ordering_tag)
      {}

       //! constructor
      GridFunctionSpace(const typename Traits::GridView& gridview, const std::shared_ptr<const FEM>& fem, const std::shared_ptr<const CE>& ce, const B& backend = B(), const O& ordering_tag = O())
#if DUNE_PDELAB_WARN_ON_GRIDVIEW_BASED_GFS
        DUNE_DEPRECATED_MSG("GridFunctionSpaces now internally use an EntitySet instead of a GridView, please replace the template parameter and the first constructor parameter by an EntitySet").
#endif
        : OrderedGFS(std::in_place, typename Traits::EntitySet{gridview}, fem, ce, backend, ordering_tag)
      {}

       //! constructor
      GridFunctionSpace(const typename Traits::GridView& gridview, const FEM& fem, const B& backend = B(), const O& ordering_tag = O())
#if DUNE_PDELAB_WARN_ON_GRIDVIEW_BASED_GFS
        DUNE_DEPRECATED_MSG("GridFunctionSpaces now internally use an EntitySet instead of a GridView, please replace the template parameter and the first constructor parameter by an EntitySet").
#endif
        : OrderedGFS(std::in_place, typename Traits::EntitySet{gridview}, stackobject_to_shared_ptr(fem), backend, ordering_tag)
      {}

       //! constructor
      GridFunctionSpace(const typename Traits::GridView& gridview, const std::shared_ptr<const FEM>& fem, const B& backend = B(), const O& ordering_tag = O())
#if DUNE_PDELAB_WARN_ON_GRIDVIEW_BASED_GFS
        DUNE_DEPRECATED_MSG("GridFunctionSpaces now internally use an EntitySet instead of a GridView, please replace the template parameter and the first constructor parameter by an EntitySet").
#endif
        : OrderedGFS(std::in_place, typename Traits::EntitySet{gridview}, fem, backend, ordering_tag)
      {}


      /**
       * @brief Construct a new Grid Function Space object
       * @warning The entity set internals will be modified according to the
       *          finite element map used codimensions
       * @warning This version of the constructor takes a reference on the fem
       *          and the ce. Therefore, these objects shall live longer than
       *          the grid function space!
       *
       * @param entitySet     Copy of an entity set
       * @param fem           Finite Element Map
       * @param ce            Constraints Assembler
       * @param backend       Vector backend
       * @param ordering_tag  Ordering tag
       */
      GridFunctionSpace(typename Traits::EntitySet entitySet, const FEM &fem,
                        const CE &ce, const B &backend = B(),
                        const O &ordering_tag = O())
          : OrderedGFS(std::in_place, std::move(entitySet), stackobject_to_shared_ptr(fem), stackobject_to_shared_ptr(ce), backend, ordering_tag)
      {}

      /**
       * @copybrief UnorderedGridFunctionSpace
       * @warning The entity set internals will be modified according to the
       *          finite element map used codimensions
       * @warning This version of the constructor takes a reference on the fem.
       *          Therefore, these objects shall live longer than the grid
       *          function space!
       *
       * @param entitySet     Copy of an entity set
       * @param fem           Finite Element Map
       * @param backend       Vector backend
       * @param ordering_tag  Ordering tag
       */
      GridFunctionSpace(typename Traits::EntitySet entitySet, const FEM &fem,
                        const B &backend = B(),
                        const O &ordering_tag = O())
          : OrderedGFS(std::in_place, std::move(entitySet), stackobject_to_shared_ptr(fem), backend, ordering_tag)
      {}
    };

    //! \}

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDFUNCTIONSPACE_HH
