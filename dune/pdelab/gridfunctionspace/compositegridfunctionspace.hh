// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEGRIDFUNCTIONSPACE_HH

#include <memory>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/utility.hh>

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
      , public PowerCompositeGridFunctionSpaceBase<
          UnorderedCompositeGridFunctionSpace<
            Backend,
            OrderingTag,
            Children...>,
          typename TypeTree::Child<TypeTree::CompositeNode<Children...>,0>::Traits::EntitySet,
          Backend,
          OrderingTag,
          sizeof...(Children)
        >
    {
      typedef TypeTree::CompositeNode<Children...> NodeT;

      typedef PowerCompositeGridFunctionSpaceBase<
        UnorderedCompositeGridFunctionSpace,
        typename TypeTree::Child<NodeT,0>::Traits::EntitySet,
        Backend,
        OrderingTag,
        sizeof...(Children)> ImplementationBase;

      friend class PowerCompositeGridFunctionSpaceBase<
        UnorderedCompositeGridFunctionSpace,
        typename TypeTree::Child<NodeT,0>::Traits::EntitySet,
        Backend,
        OrderingTag,
        sizeof...(Children)>;

      template<typename,typename>
      friend class GridFunctionSpaceBase;

    public:
      typedef CompositeGridFunctionSpaceTag ImplementationTag;

      typedef typename ImplementationBase::Traits Traits;

      // ********************************************************************************
      // constructors for stack-constructed children passed in by reference
      // ********************************************************************************

      UnorderedCompositeGridFunctionSpace(const Backend& backend, Children&... children)
        : NodeT(TypeTree::assertGridViewType<typename NodeT::template Child<0>::Type>(children)...)
        , ImplementationBase(backend,OrderingTag())
      { }

      UnorderedCompositeGridFunctionSpace(const OrderingTag& ordering_tag, Children&... children)
        : NodeT(TypeTree::assertGridViewType<typename NodeT::template Child<0>::Type>(children)...)
        , ImplementationBase(Backend(),ordering_tag)
      { }

      UnorderedCompositeGridFunctionSpace(const Backend& backend, const OrderingTag& ordering_tag, Children&... children)
        : NodeT(TypeTree::assertGridViewType<typename NodeT::template Child<0>::Type>(children)...)
        , ImplementationBase(backend,ordering_tag)
      { }

      UnorderedCompositeGridFunctionSpace(Children&... children)
        : NodeT(TypeTree::assertGridViewType<typename NodeT::template Child<0>::Type>(children)...)
        , ImplementationBase(Backend(),OrderingTag())
      { }

      // ********************************************************************************
      // constructors for heap-constructed children passed in as shared_ptrs
      // ********************************************************************************

      UnorderedCompositeGridFunctionSpace(const Backend& backend, std::shared_ptr<Children>... children)
        : NodeT(children...)
        , ImplementationBase(backend,OrderingTag())
      { }

      UnorderedCompositeGridFunctionSpace(const OrderingTag& ordering_tag, std::shared_ptr<Children>... children)
        : NodeT(children...)
        , ImplementationBase(Backend(),ordering_tag)
      { }

      UnorderedCompositeGridFunctionSpace(const Backend& backend, const OrderingTag& ordering_tag, std::shared_ptr<Children>... children)
        : NodeT(children...)
        , ImplementationBase(backend,ordering_tag)
      { }

      UnorderedCompositeGridFunctionSpace(std::shared_ptr<Children>... children)
        : NodeT(children...)
        , ImplementationBase(Backend(),OrderingTag())
      { }

    };


    template<typename Backend,
             typename OrderingTag,
             typename... Children>
    class CompositeGridFunctionSpace
      : public OrderedGridFunctionSpace<UnorderedCompositeGridFunctionSpace<Backend,OrderingTag,Children...>>
    {
      using Base = OrderedGridFunctionSpace<UnorderedCompositeGridFunctionSpace<Backend,OrderingTag,Children...>>;
    public:
      using Base::Base;
    };


    //! \}

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_COMPOSITEGRIDFUNCTIONSPACE_HH
