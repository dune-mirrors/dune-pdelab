// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_ORDERING_TRANSFORMATIONS_HH
#define DUNE_PDELAB_ORDERING_TRANSFORMATIONS_HH

#include <cstddef>
#include <algorithm>

#include <dune/typetree/traversal.hh>
#include <dune/typetree/accumulate_static.hh>

#include <dune/pdelab/common/typetraits.hh>
#include <dune/pdelab/common/multiindex.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>

/**
 * \file
 * \brief Define and register ordering related transformations.
 * This header defines the two transformations gfs_to_ordering and gfs_to_local_ordering.
 * As both the PowerGFS and the CompositeGFS transformations use a single descriptor that
 * is specialized for the individual orderings, the prototype of those descriptors are also
 * declared in this file and are registered with the TypeTree transformation system.
 */

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

#ifndef DOXYGEN

  /**
   * @brief Integral constant visitor to extract max blocking depth
   * @details This visitor is intended to be used with an unordered grid
   *          function space tree where each node gives a vector backend
   *          containing the information whether such node should be blocked or
   *          not (i.e. `Node::Traits::Backend::Traits::max_blocking_depth`).
   */
  struct extract_max_container_depth
      : public TypeTree::Experimental::DefaultHybridVisitor,
        public TypeTree::StaticTraversal,
        public TypeTree::VisitDirectChildren
  {
    template <class Node, class TreePath, class U>
    auto leaf(const Node &, TreePath, U) const {
      // is leaf node blocked?
      return index_constant<Node::Traits::Backend::Traits::max_blocking_depth>{};
    }

    template <class Node, class Child, class TreePath, class ChildIndex,
              class CarryValue>
    auto beforeChild(const Node&, const Child& child, TreePath, ChildIndex,
                     CarryValue) const {
      using namespace TypeTree::Experimental;
      // extract max container for depth for child i
      using ChildDepth = decltype(hybridApplyToTree(child, *this, Indices::_0));
      // return max value between child and carried value (previous children)
      return index_constant<std::max(CarryValue::value, ChildDepth::value)>{};
    }

    template <class Node, class TreePath, class ChildrenCount>
    auto post(const Node &, TreePath, ChildrenCount) const {
      // After all children accumulated the max of their depths, we add the
      // blocking value of this node
      const std::size_t blocked = Node::Traits::Backend::Traits::max_blocking_depth;
      return index_constant<ChildrenCount::value + blocked>{};
    }
  };

    //! GridFunctionSpace to Ordering transformation descriptor
    template<typename RootGFS>
    struct gfs_to_ordering
    {
      // extract max blocking depth (notice that this is an unevaluated context,
      // and the resulting type know the depth at compile-time)
      static const std::size_t ci_depth = decltype(
        TypeTree::Experimental::hybridApplyToTree(
            std::declval<RootGFS>(),
            extract_max_container_depth{},
            Indices::_0)
          )::value + 1;

      typedef typename gfs_to_lfs<RootGFS>::DOFIndex DOFIndex;
      typedef MultiIndex<std::size_t,ci_depth> ContainerIndex;
    };

    //! GridFunctionSpace to LocalOrdering transformation descriptor
    template<typename GlobalTransformation>
    struct gfs_to_local_ordering
    {
      typedef typename GlobalTransformation::DOFIndex DOFIndex;
      typedef typename GlobalTransformation::ContainerIndex ContainerIndex;
    };


    // Declare PowerGFS to ordering descriptor and register transformation

    template<typename GFS, typename Transformation, typename OrderingTag>
    struct power_gfs_to_ordering_descriptor
      : public TypeTree::meta_function
    {
      typedef decltype(
        register_power_gfs_to_ordering_descriptor(
          TypeTree::declptr<GFS>(),
          TypeTree::declptr<Transformation>(),
          TypeTree::declptr<OrderingTag>()
          )
        ) type;
    };

    template<typename GridFunctionSpace, typename Params>
    power_gfs_to_ordering_descriptor<
      GridFunctionSpace,
      gfs_to_ordering<Params>,
      typename GridFunctionSpace::OrderingTag
      >
    registerNodeTransformation(GridFunctionSpace*, gfs_to_ordering<Params>*, PowerGridFunctionSpaceTag*);


    // Declare LeafGFS to ordering descriptor and register transformation

    template<typename GFS, typename Transformation, typename OrderingTag>
    struct leaf_gfs_to_ordering_descriptor
      : public TypeTree::meta_function
    {
      typedef decltype(
        register_leaf_gfs_to_ordering_descriptor(
          TypeTree::declptr<GFS>(),
          TypeTree::declptr<Transformation>(),
          TypeTree::declptr<OrderingTag>()
          )
        ) type;
    };

    template<typename GridFunctionSpace, typename Params>
    leaf_gfs_to_ordering_descriptor<
      GridFunctionSpace,
      gfs_to_ordering<Params>,
      typename GridFunctionSpace::Traits::OrderingTag
      >
    registerNodeTransformation(GridFunctionSpace*, gfs_to_ordering<Params>*, LeafGridFunctionSpaceTag*);


    // Declare CompositeGFS to ordering descriptor and register transformation

    template<typename GFS, typename Transformation, typename OrderingTag>
    struct composite_gfs_to_ordering_descriptor
      : public TypeTree::meta_function
    {
      typedef decltype(
        register_composite_gfs_to_ordering_descriptor(
          TypeTree::declptr<GFS>(),
          TypeTree::declptr<Transformation>(),
          TypeTree::declptr<OrderingTag>()
          )
        ) type;
    };

    template<typename GridFunctionSpace, typename Params>
    composite_gfs_to_ordering_descriptor<
      GridFunctionSpace,
      gfs_to_ordering<Params>,
      typename GridFunctionSpace::OrderingTag
      >
    registerNodeTransformation(GridFunctionSpace*, gfs_to_ordering<Params>*, CompositeGridFunctionSpaceTag*);


    // Declare PowerGFS to local ordering descriptor and register transformation

   template<typename GFS, typename Transformation, typename OrderingTag>
   struct power_gfs_to_local_ordering_descriptor;

    template<typename GFS, typename Params>
    power_gfs_to_local_ordering_descriptor<
      GFS,
      gfs_to_local_ordering<Params>,
      typename GFS::OrderingTag
      >
    registerNodeTransformation(GFS*, gfs_to_local_ordering<Params>*, PowerGridFunctionSpaceTag*);


    // Declare LeafGFS to local ordering descriptor and register transformation

    template<typename GFS, typename Transformation, typename OrderingTag>
    struct leaf_gfs_to_local_ordering_descriptor;

    template<typename GFS, typename Params>
    leaf_gfs_to_local_ordering_descriptor<
      GFS,
      gfs_to_local_ordering<Params>,
      typename GFS::Traits::OrderingTag
      >
    registerNodeTransformation(GFS*, gfs_to_local_ordering<Params>*, LeafGridFunctionSpaceTag*);


    // Declare CompositeGFS to ordering descriptor and register transformation

    template<typename GFS, typename Transformation, typename OrderingTag>
    struct composite_gfs_to_local_ordering_descriptor;

    template<typename GFS, typename Params>
    composite_gfs_to_local_ordering_descriptor<
      GFS,
      gfs_to_local_ordering<Params>,
      typename GFS::OrderingTag
      >
    registerNodeTransformation(GFS*, gfs_to_local_ordering<Params>*, CompositeGridFunctionSpaceTag*);


#endif // DOXYGEN

    //! \} group ordering

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_TRANSFORMATIONS_HH
