// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDEREDGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDEREDGRIDFUNCTIONSPACE_HH

#include <dune/pdelab/gridfunctionspace/datahandleprovider.hh>

#include <dune/typetree/visitor.hh>

namespace Dune {
  namespace PDELab {

    namespace Impl {
    /**
     * @brief Constraint type accumulator
     * @details This class is meant to accumulate the type of constraints that
     * leaf grid function spaces have. If all leaf nodes have no constraints, an
     * hybrid accumulate on the tree will return true, otherwise false.
     */
      struct HasNoConstraintsVisitor
          : public TypeTree::Experimental::DefaultHybridVisitor,
            public TypeTree::StaticTraversal,
            public TypeTree::VisitTree
      {

        template <class LeafGFS, class TreePath, class HasNoConstraints>
        auto leaf(const LeafGFS &, TreePath, HasNoConstraints) const {
          // get constraint type from this leaf node
          using Constraints = typename LeafGFS::Traits::ConstraintsType;
          // is node constrained?
          const bool node_has_no_constraints =
              std::is_same_v<Constraints, NoConstraints>;
          // return conjunction on carried value and this node value
          return std::bool_constant < HasNoConstraints{} &&
                node_has_no_constraints > {};
        }
      };
    } // namespace Impl

    /**
     * @brief Add ordering and data handler to an unordered grid function space
     * @details This class extends the grid function space interface to have
     * an ordering and a data handler for parallel communication. In particular,
     * this is intened to only be used on the root node.
     *
     * @tparam UnorderedGFS  A grid function space that has no ordering
     *
     * \see Ordering
     */
    template <class UnorderedGFS>
    class OrderedGridFunctionSpace
        : public UnorderedGFS,
          public DataHandleProvider<OrderedGridFunctionSpace<UnorderedGFS>>
    {
      using ordering_transformation =
          TypeTree::TransformTree<UnorderedGFS, gfs_to_ordering<UnorderedGFS>>;

    public:
      //! Inherit constructor from UnorderedGFS
      using UnorderedGFS::UnorderedGFS;

      //! Ordering tree type
      using Ordering = typename ordering_transformation::Type;

      //! extract type for storing constraints
      template <typename E>
      class ConstraintsContainer {

        //! evaluation of HasNoConstraintsVisitor on the UnorderedGFS
        static constexpr bool has_no_constraints =
            decltype(
              TypeTree::Experimental::hybridApplyToTree(
                std::declval<UnorderedGFS>(),
                Impl::HasNoConstraintsVisitor{},
                std::true_type{})
              ){};

      public:
        //! \brief define Type as the Type of a container of E's
        using Type = std::conditional_t<
            has_no_constraints,
            EmptyTransformation,
            ConstraintsTransformation<typename Ordering::Traits::DOFIndex,
                                      typename Ordering::Traits::ContainerIndex,
                                      E>>;
      };

      //------------------------------

      //! Direct access to the DOF ordering.
      const Ordering &ordering() const { return *orderingStorage(); }

      //! Direct access to the DOF ordering.
      Ordering &ordering() { return *orderingStorage(); }

      //! Direct access to the storage of the DOF ordering.
      std::shared_ptr<const Ordering> orderingStorage() const {
        if (!this->isRootSpace()) {
          DUNE_THROW(GridFunctionSpaceHierarchyError,
                    "Ordering can only be obtained for root space in "
                    "GridFunctionSpace tree.");
        }
        if (!_ordering) {
          create_ordering();
          this->update(*_ordering);
        }
        return _ordering;
      }

      //! Direct access to the storage of the DOF ordering.
      std::shared_ptr<Ordering> orderingStorage() {
        if (!this->isRootSpace()) {
          DUNE_THROW(GridFunctionSpaceHierarchyError,
                    "Ordering can only be obtained for root space in "
                    "GridFunctionSpace tree.");
        }
        if (!_ordering) {
          create_ordering();
          this->update(*_ordering);
        }
        return _ordering;
      }

        //! Update the indexing information of the GridFunctionSpace.
        /**
         * @param force  Set to true if the underlying grid has changed (e.g.
         *      due to adaptivity) to force an update of the embedded EntitySet.
         */
        void update(bool force = false)
        {
          // every node in the gfs may have different entity sets
          // we only update leaf and root nodes
          auto& entity_set = this->_entity_set;
          if (entity_set)
            entity_set->update(force);
          auto update_leaf_es = impl::update_leaf_entity_set{entity_set, force};
          TypeTree::applyToTree(*this, update_leaf_es);
          // We bypass the normal access using ordering() here to avoid a double
          // update if the Ordering has not been created yet.
          if (!_ordering)
            create_ordering();
          update(*_ordering);
        }

    private:
      void update(Ordering &ordering) const {
        if (!this->isRootSpace()) {
          DUNE_THROW(GridFunctionSpaceHierarchyError,
                    "update() may only be called on the root of the function "
                    "space hierarchy");
        }
        ordering.update();
        TypeTree::applyToTree(
            ordering,
            impl::update_ordering_data<typename UnorderedGFS::Traits::SizeType>(
                ordering));
      }

      // This method here is to avoid a double update of the Ordering when the user
      // calls GFS::update() before GFS::ordering().
      void create_ordering() const {
        _ordering =
            std::make_shared<Ordering>(ordering_transformation::transform(*this));
      }

      mutable std::shared_ptr<Ordering> _ordering;
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDEREDGRIDFUNCTIONSPACE_HH
