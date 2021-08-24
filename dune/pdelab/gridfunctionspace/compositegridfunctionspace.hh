// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEGRIDFUNCTIONSPACE_HH

#include <memory>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/utility.hh>

#include <dune/pdelab/common/utility.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspacebase.hh>
#include <dune/pdelab/gridfunctionspace/powercompositegridfunctionspacebase.hh>
#include <dune/pdelab/gridfunctionspace/orderedgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>

namespace Dune {
  namespace PDELab {

    //=======================================
    // composite grid function space
    //=======================================

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

    /** \brief base class for tuples of grid function spaces
        base class that holds implementation of the methods
        this is the default version with lexicographic ordering
        \tparam Mapper is the ordering parameter. Use e.g.
        \link GridFunctionSpaceLexicographicMapper GridFunctionSpaceLexicographicMapper \endlink
        or \link  GridFunctionSpaceComponentBlockwiseMapper  GridFunctionSpaceComponentBlockwiseMapper \endlink
        or \link  GridFunctionSpaceBlockwiseMapper  GridFunctionSpaceBlockwiseMapper \endlink
        or \link  GridFunctionSpaceDynamicBlockwiseMapper  GridFunctionSpaceDynamicBlockwiseMapper \endlink
        \tparam Ti are all grid function spaces
    */
    template<typename Backend,
             typename OrderingTag,
             typename... Children>
    class UnorderedCompositeGridFunctionSpace
      : public TypeTree::CompositeNode<Children...>
      , public GridFunctionSpaceNode<GridFunctionSpaceNodeTraits<Backend,OrderingTag>>
    {
      typedef TypeTree::CompositeNode<Children...> NodeT;

      typedef GridFunctionSpaceNode<GridFunctionSpaceNodeTraits<Backend,OrderingTag>> ImplementationBase;

    public:
      typedef CompositeGridFunctionSpaceTag ImplementationTag;

      typedef typename ImplementationBase::Traits Traits;

      //! Initialize the CompositeNode with a copy of the passed-in storage
      //! type.
      UnorderedCompositeGridFunctionSpace(
          const std::tuple<std::shared_ptr<Children>...> &children,
          const Backend &backend = Backend(),
          const OrderingTag ordering_tag = OrderingTag())
          : NodeT(children), ImplementationBase(backend, ordering_tag) {}
    };

    //! \copydoc UnorderedCompositeGridFunctionSpace
    template<class Backend,
             class OrderingTag,
             class... Children>
    class CompositeGridFunctionSpace
      : public OrderedGridFunctionSpace<UnorderedCompositeGridFunctionSpace<Backend,OrderingTag,Children...>,typename Impl::FirstLeaf<Children...>::Traits::EntitySet>
    {
      using Base = OrderedGridFunctionSpace<UnorderedCompositeGridFunctionSpace<Backend,OrderingTag,Children...>, typename Impl::FirstLeaf<Children...>::Traits::EntitySet>;

    public:

      // ********************************************************************************
      // constructors for stack-constructed children passed in by reference
      // ********************************************************************************

      CompositeGridFunctionSpace(const Backend& backend, Children&... children)
        : Base(std::in_place, std::make_tuple(stackobject_to_shared_ptr(children)...), backend)
      { }

      CompositeGridFunctionSpace(const OrderingTag& ordering_tag, Children&... children)
        : Base(std::in_place, std::make_tuple(stackobject_to_shared_ptr(children)...), Backend{}, ordering_tag)
      { }

      CompositeGridFunctionSpace(const Backend& backend, const OrderingTag& ordering_tag, Children&... children)
        : Base(std::in_place, std::make_tuple(stackobject_to_shared_ptr(children)...), backend, ordering_tag)
      { }

      CompositeGridFunctionSpace(Children&... children)
        : Base(std::in_place, std::make_tuple(stackobject_to_shared_ptr(children)...))
      { }

      // ********************************************************************************
      // constructors for heap-constructed children passed in as shared_ptrs
      // ********************************************************************************

      CompositeGridFunctionSpace(const Backend& backend, const std::shared_ptr<Children>&... children)
        : Base(std::in_place, std::make_tuple(children...), backend)
      { }

      CompositeGridFunctionSpace(const OrderingTag& ordering_tag, const std::shared_ptr<Children>&... children)
        : Base(std::in_place, std::make_tuple(children...), Backend{}, ordering_tag)
      { }

      CompositeGridFunctionSpace(const Backend& backend, const OrderingTag& ordering_tag, const std::shared_ptr<Children>&... children)
        : Base(std::in_place, std::make_tuple(children...), Backend{}, ordering_tag)
      { }

      CompositeGridFunctionSpace(const std::shared_ptr<Children>&... children)
        : Base(std::in_place, std::make_tuple(children...))
      { }
    };


    //! \}

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEGRIDFUNCTIONSPACE_HH
