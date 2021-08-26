// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_VECTORGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_VECTORGRIDFUNCTIONSPACE_HH

#include <algorithm>
#include <cstddef>
#include <memory>

#include <dune/common/shared_ptr.hh>

#include <dune/typetree/powernode.hh>

#include <dune/pdelab/gridfunctionspace/orderedgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>

namespace Dune {
  namespace PDELab {

    //=======================================
    // vector grid function space
    //=======================================

    /** \brief tensorproduct space representing a vector valued function space

        In its structure this space is very similar to a
        PowerGridFunctionSpace:

        VGFS(FEM,k) = PGFS(GFS(FEM),k) = {GFS(FEM)}^k

        Stating explicitly that a space is a VectorGridFunctionSpace
        mainly changes the way the data is interpreted. One can
        immediatelly create a discrete function as a member of a
        VectorGridFunctionSpace and visualize it via VTK. In this case
        the output data is automatically tagged as vector valued data,
        allowing for a better visualization.

        \tparam GV               Type implementing GridView
        \tparam FEM              Type implementing FiniteElementMapInterface
        \tparam k                Physical dimension of the space
        \tparam Backend          Linear algebra backend type at the level of the tensorproduct construction (the same backend one might pass to a PowerGridFunctionSpace)
        \tparam LeafBackend      Linear algebra backend type at the level of the underlying scalar space (GFS(FEM))
        \tparam Constraints      Type for constraints assembler
        \tparam OrderingTag      ordering of DOFs at the level of the tensorproduct construction (usually on will choose either \link LexicographicOrderingTag or \link EntityBlockedOrderingTag)
        \tparam LeafOrderingTag  ordering of DOFs at the level of the underlying scalar space (default: DefaultLeafOrderingTag)
    */
    template<typename ES,
             typename FEM,
             std::size_t k,
             typename Backend,
             typename LeafBackend,
             typename Constraints = NoConstraints,
             typename OrderingTag = LexicographicOrderingTag,
             typename LeafOrderingTag = DefaultLeafOrderingTag>
    class UnorderedVectorGridFunctionSpace
      : public TypeTree::PowerNode<UnorderedGridFunctionSpace<
                                     ES,
                                     FEM,
                                     Constraints,
                                     LeafBackend,
                                     LeafOrderingTag
                                     >,
                                   k>
       , public GridFunctionSpaceNode<GridFunctionSpaceNodeTraits<Backend,OrderingTag>>
      , public GridFunctionOutputParameters
    {

      static_assert(isEntitySet<ES>{});

      typedef UnorderedGridFunctionSpace<
        ES,
        FEM,
        Constraints,
        LeafBackend,
        LeafOrderingTag
        > LeafGFS;

    public:

      typedef VectorGridFunctionSpaceTag ImplementationTag;

      typedef TypeTree::PowerNode<LeafGFS,k> BaseT;

      typedef GridFunctionSpaceNode<GridFunctionSpaceNodeTraits<Backend,OrderingTag>> ImplementationBase;

    protected:

      // Preconstruct children - it is important that the children are set before entering the constructor
      // of ImplementationBase!
      static typename BaseT::NodeStorage create_components(const typename LeafGFS::Traits::EntitySet& es,
                                                           std::shared_ptr<const FEM> fem_ptr,
                                                           const LeafBackend& leaf_backend,
                                                           const LeafOrderingTag& leaf_ordering_tag)
      {
        typename BaseT::NodeStorage r;
        auto leaf_gfs = std::make_shared<LeafGFS>(es,fem_ptr,leaf_backend,leaf_ordering_tag);
        for (std::size_t i = 0; i < k; ++i)
          r[i] = leaf_gfs;
        return r;
      }

    public:

      UnorderedVectorGridFunctionSpace(const typename LeafGFS::Traits::EntitySet& es, std::shared_ptr<const FEM> fem,
                                    const Backend& backend = Backend(), const LeafBackend& leaf_backend = LeafBackend(),
                                    const OrderingTag& ordering_tag = OrderingTag(), const LeafOrderingTag& leaf_ordering_tag = LeafOrderingTag())
              : BaseT(create_components(es, fem, leaf_backend, leaf_ordering_tag))
              , ImplementationBase(backend,ordering_tag)
      {}

      using ImplementationBase::name;

      void name(const std::string& name)
      {
        ImplementationBase::name(name);
        for (std::size_t i = 0; i < k; ++i)
          {
            std::stringstream ns;
            ns << name << "_" << i;
            this->child(i).name(ns.str());
          }
      }
    };

    //! \copydoc UnorderedVectorGridFunctionSpace
    template<typename ES,
             typename FEM,
             std::size_t k,
             typename Backend,
             typename LeafBackend,
             typename Constraints = NoConstraints,
             typename OrderingTag = LexicographicOrderingTag,
             typename LeafOrderingTag = DefaultLeafOrderingTag>
    class VectorGridFunctionSpace
      : public OrderedGridFunctionSpace<UnorderedVectorGridFunctionSpace<impl::EntitySet<ES>,FEM,k,Backend,LeafBackend,Constraints,OrderingTag,LeafOrderingTag>,impl::EntitySet<ES>>
    {
      using Base = OrderedGridFunctionSpace<UnorderedVectorGridFunctionSpace<impl::EntitySet<ES>,FEM,k,Backend,LeafBackend,Constraints,OrderingTag,LeafOrderingTag>,impl::EntitySet<ES>>;
      using LeafGFS = typename Base::ChildType;
    public:

      struct Traits : public Base::Traits {
        enum {
          //! \brief True if this grid function space is composed of others.
          isComposite
          [[deprecated("This enum will be removed after PDELab 2.9.")]] = 1,
          //! \brief number of child spaces
          noChilds
          [[deprecated("This enum will be removed after PDELab 2.9.")]] = k
        };

        [[deprecated(
            "This enum will be removed after PDELab 2.9. Use degree() from the "
            "TypeTree base class")]] const static std::size_t CHILDREN = k;

        //! \brief mapper
        using MapperType [[deprecated("This enum will be removed after PDELab "
                                      "2.9. Use OrderingTag instead.")]] =
            typename Base::Traits::OrderingTag;
      };

      VectorGridFunctionSpace(const typename LeafGFS::Traits::GridView& gv, const FEM& fem,
                              const Backend& backend = Backend(), const LeafBackend& leaf_backend = LeafBackend(),
                              const OrderingTag& ordering_tag = OrderingTag(), const LeafOrderingTag& leaf_ordering_tag = LeafOrderingTag())
        : Base(std::in_place, typename LeafGFS::Traits::EntitySet{gv}, stackobject_to_shared_ptr(fem), backend, leaf_backend, ordering_tag, leaf_ordering_tag)
      {}

      VectorGridFunctionSpace(typename LeafGFS::Traits::EntitySet es, const FEM& fem,
                              const Backend& backend = Backend(), const LeafBackend& leaf_backend = LeafBackend(),
                              const OrderingTag& ordering_tag = OrderingTag(), const LeafOrderingTag& leaf_ordering_tag = LeafOrderingTag())
        : Base(std::in_place, es, stackobject_to_shared_ptr(fem), backend, leaf_backend, ordering_tag, leaf_ordering_tag)
      {}

      VectorGridFunctionSpace(const typename LeafGFS::Traits::GridView& gv, std::shared_ptr<const FEM> fem,
                              const Backend& backend = Backend(), const LeafBackend& leaf_backend = LeafBackend(),
                              const OrderingTag& ordering_tag = OrderingTag(), const LeafOrderingTag& leaf_ordering_tag = LeafOrderingTag())
        : Base(std::in_place, typename LeafGFS::Traits::EntitySet{gv}, fem, backend, leaf_backend, ordering_tag, leaf_ordering_tag)
      {}

      VectorGridFunctionSpace(typename LeafGFS::Traits::EntitySet es, std::shared_ptr<const FEM> fem,
                                    const Backend& backend = Backend(), const LeafBackend& leaf_backend = LeafBackend(),
                                    const OrderingTag& ordering_tag = OrderingTag(), const LeafOrderingTag& leaf_ordering_tag = LeafOrderingTag())
        : Base(std::in_place, es, fem, backend, leaf_backend, ordering_tag, leaf_ordering_tag)
      {}
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_VECTORGRIDFUNCTIONSPACE_HH
