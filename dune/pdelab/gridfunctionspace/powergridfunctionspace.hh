// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_POWERGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_POWERGRIDFUNCTIONSPACE_HH

#include <cstddef>
#include <memory>

#include <dune/typetree/powernode.hh>

#include <dune/pdelab/gridfunctionspace/powercompositegridfunctionspacebase.hh>
#include <dune/pdelab/gridfunctionspace/orderedgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>

namespace Dune {
  namespace PDELab {

    //=======================================
    // power grid function space
    //=======================================

    /** \brief base class for tuples of grid function spaces
        product of identical grid function spaces
        base class that holds implementation of the methods

        PGFS(T,k) = {T}^k

        \tparam T the underlying are all grid function spaces
        \tparam k power factor
        \tparam Mapper is the ordering parameter. Use e.g.
        \link GridFunctionSpaceLexicographicMapper GridFunctionSpaceLexicographicMapper \endlink
        or \link  GridFunctionSpaceComponentBlockwiseMapper  GridFunctionSpaceComponentBlockwiseMapper \endlink
        or \link  GridFunctionSpaceBlockwiseMapper  GridFunctionSpaceBlockwiseMapper \endlink
        or \link  GridFunctionSpaceDynamicBlockwiseMapper  GridFunctionSpaceDynamicBlockwiseMapper \endlink
    */
    template<typename T, std::size_t k,
             typename Backend,
             typename OrderingTag = LexicographicOrderingTag>
    class UnorderedPowerGridFunctionSpace
      : public TypeTree::PowerNode<T,k>
      , public PowerCompositeGridFunctionSpaceBase<
          UnorderedPowerGridFunctionSpace<T, k, Backend, OrderingTag>,
          typename T::Traits::EntitySet,
          Backend,
          OrderingTag,
          k>
    {

    public:

      typedef PowerGridFunctionSpaceTag ImplementationTag;

      typedef TypeTree::PowerNode<T,k> BaseT;

    private:

      typedef PowerCompositeGridFunctionSpaceBase<
        UnorderedPowerGridFunctionSpace,
        typename T::Traits::EntitySet,
        Backend,
        OrderingTag,
        k
        > ImplementationBase;

      friend class PowerCompositeGridFunctionSpaceBase<
        UnorderedPowerGridFunctionSpace,
        typename T::Traits::EntitySet,
        Backend,
        OrderingTag,
        k>;

      template<typename,typename>
      friend class GridFunctionSpaceBase;

    public:

      //! export traits class
      typedef typename ImplementationBase::Traits Traits;

      /**
       * @brief Construct a new Power Grid Function Space object
       *
       * @param container     array with pointers to child spaces
       * @param backend       backend object
       * @param ordering_tag  ordering tag object
       */
      UnorderedPowerGridFunctionSpace(const std::array<std::shared_ptr<T>,k>& container, const Backend& backend = Backend(), const OrderingTag ordering_tag = OrderingTag())
        : BaseT(container)
        , ImplementationBase(backend,ordering_tag)
      {}

      UnorderedPowerGridFunctionSpace(T& c, const Backend& backend = Backend(), const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c)
        , ImplementationBase(backend,ordering_tag)
      {}

      UnorderedPowerGridFunctionSpace (T& c0,
                              T& c1,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1)
        , ImplementationBase(backend,ordering_tag)
      {}

      UnorderedPowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1,c2)
        , ImplementationBase(backend,ordering_tag)
      {}

      UnorderedPowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1,c2,c3)
        , ImplementationBase(backend,ordering_tag)
      {}

      UnorderedPowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1,c2,c3,c4)
        , ImplementationBase(backend,ordering_tag)
      {}

      UnorderedPowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1,c2,c3,c4,c5)
        , ImplementationBase(backend,ordering_tag)
      {}

      UnorderedPowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1,c2,c3,c4,c5,c6)
        , ImplementationBase(backend,ordering_tag)
      {}

      UnorderedPowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7)
        , ImplementationBase(backend,ordering_tag)
      {}

      UnorderedPowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              T& c8,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8)
        , ImplementationBase(backend,ordering_tag)
      {}

      UnorderedPowerGridFunctionSpace (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              T& c8,
                              T& c9,
                              const Backend& backend = Backend(),
                              const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
        , ImplementationBase(backend,ordering_tag)
      {}

      template<typename Child0, typename... Children>
      UnorderedPowerGridFunctionSpace(std::shared_ptr<Child0> child0, std::shared_ptr<Children>... children)
        : BaseT(child0, children...)
        , ImplementationBase(Backend(),OrderingTag())
      {}

    };

    template<typename T, std::size_t k,
             typename Backend,
             typename OrderingTag = LexicographicOrderingTag>
    class PowerGridFunctionSpace
      : public OrderedGridFunctionSpace<UnorderedPowerGridFunctionSpace<T,k,Backend,OrderingTag>>
    {
      using Base = OrderedGridFunctionSpace<UnorderedPowerGridFunctionSpace<T,k,Backend,OrderingTag>>;
    public:
      using Base::Base;
    };


  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_POWERGRIDFUNCTIONSPACE_HH
