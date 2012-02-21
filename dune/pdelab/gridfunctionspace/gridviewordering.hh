// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDVIEWORDERING_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_GRIDVIEWORDERING_HH

#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/gridfunctionspace/orderingutility.hh>
#include <dune/pdelab/gridfunctionspace/localorderingdynamicbase.hh>

namespace Dune {
  namespace PDELab {



    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //! Dummy ordering for leaf gridfunctionspaces
    template<typename LocalOrdering>
    class LeafGridViewOrdering
      : public TypeTree::VariadicCompositeNode<LocalOrdering>
      , public VirtualOrderingBase<typename LocalOrdering::Traits::DOFIndex,typename LocalOrdering::Traits::ContainerIndex>
    {
    public:
      typedef typename LocalOrdering::Traits Traits;

    private:

      typedef typename Traits::GridView GV;

      typedef TypeTree::VariadicCompositeNode<LocalOrdering> NodeT;

      GV _gv;

    public:

      LocalOrdering& localOrdering()
      {
        return this->template child<0>();
      }

      const LocalOrdering& localOrdering() const
      {
        return this->template child<0>();
      }


      LeafGridViewOrdering(const typename NodeT::NodeStorage& localOrdering)
        : NodeT(localOrdering)
        , _gv(this->template child<0>().gridView())
      {}

      virtual void map_index_dynamic(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {
        map_index(di,ci);
      }

      void map_index(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {
        const typename Traits::SizeType geometry_type_index = di.entityIndex()[0];
        const typename Traits::SizeType entity_index = di.entityIndex()[1];
        assert (di.treeIndex().size() == 1);
        ci.push_back(di.treeIndex().back());
        if (_backend_blocked)
          {
            ci.push_back(localOrdering()._gt_entity_offsets[geometry_type_index] + entity_index);
          }
        else if (localOrdering()._fixed_size)
          {
            ci.back() += _gt_dof_offsets[geometry_type_index] + entity_index * localOrdering()._gt_dof_offsets[geometry_type_index];
          }
        else
          {
            ci.back() += localOrdering()._entity_dof_offsets[localOrdering()._gt_entity_offsets[geometry_type_index] + entity_index];
          }
      }

      //! update internal data structures
      /**
       * In general this method must be called after initialization and every
       * time the structure of the dof-vector of one of gfs's children
       * changes.  For this particular ordering however this method does
       * nothing.
       */
      void _recursive_update()
      {
        LocalOrdering& lo = localOrdering();
        lo.update_a_priori_fixed_size();

        const std::size_t dim = GV::dimension;

        typedef typename Traits::SizeType size_type;
        typedef std::vector<GeometryType> GTVector;
        GTVector geom_types;

        for (size_type cc = 0; cc <= dim; ++cc)
          {
            const GTVector& per_codim_geom_types = _gv.indexSet().geometryTypes(cc);
            std::copy(per_codim_geom_types.begin(),per_codim_geom_types.end(),std::back_inserter(geom_types));
          }

        if (lo._fixed_size)
          {
            lo.update_fixed_size(geom_types);

            _gt_dof_offsets.resize(GlobalGeometryTypeIndex::size(dim) + 1);

            const GTVector::const_iterator end_it = geom_types.end();
            for (GTVector::const_iterator it = geom_types.begin(); it != end_it; ++it)
              {
                const size_type gt_index = GlobalGeometryTypeIndex::index(*it);
                _gt_dof_offsets[gt_index + 1] = lo.size(gt_index,0) * _gv.indexSet().size(*it);
              }
            std::partial_sum(_gt_dof_offsets.begin(),_gt_dof_offsets.end(),_gt_dof_offsets.begin());
          }
        else
          {

          }
      }

      //! dofs are blocked per entity/intersection on the leafs
      bool blocked() const { return true; }

      //! \brief whether all entites of the same geometry type/all
      //!        intersections have the same number of dofs
      /**
       * On the leaf this is realized by iterating over the grid during update
       * an checking.
       *
       * \note Even if fixedSize()==true the number of dofs may still vary
       *       between entities od different geometry type or between entities
       *       and intersections.
       */
      bool fixedSize() const { return localOrdering()._fixed_size; }

      //! number of indices in this ordering
      typename Traits::SizeType size() const { return _size; }

      //! \brief maximum number of dofs attached to any given element and all
      //!        of its subentities and intersections
      /**
       * This is generally not an exact maximum and may be bigger than the
       * actual maximum.  There is however one special case: it is guaranteed
       * to be the exact maximum for fixedSize()==true.
       */
      typename Traits::SizeType maxLocalSize() const { return _max_local_size; }

#if 0
      //! \brief number of indices attached to a given entity (of arbitrary
      //!        codimension)
      /**
       * \note If the grid does not support a given entity type, it may still
       *       be possible to get this information using entitySize(const
       *       Element &e, std::size_t codim, std::size_t subentity).
       */
      template<class Entity>
      typename Traits::SizeType entitySize(const Entity &e) const { return gfs.entitySize(e); }
      //! number of indices attached to a given subentity of an element
      /**
       * This method determines the number of indices attached to a subentity
       * of the given codim 0 entity.  If the grid (and the ordering) directly
       * supports entities of the given codimension, this is equivalent to
       * calling entitySize((*e.subEntity<codim>(subentity)).
       */
      template<class Element>
      typename Traits::SizeType entitySize(const Element &e, std::size_t codim,
                          std::size_t subentity) const
      { return gfs.entitySize(e, codim, subentity); }
      //! number of indices attached to a given intersection
      template<class Intersection>
      typename Traits::SizeType intersectionSize(const Intersection &i) const
      { return gfs.intersectionSize(i); }

      //! \brief offset of the block of dofs attached to a given entity (of
      //!        arbitrary codimension)
      /**
       * \note If the grid does not support a given entity type, it may still
       *       be possible to get this information using entityOffset(const
       *       Element &e, std::size_t codim, std::size_t subentity).
       *
       * \throw NotImplemented        If this EntityType is not supported by
       *                              the ordering.
       * \throw InvalidStateException If blocked()==false.
       */
      template<class Entity>
      typename Traits::SizeType entityOffset(const Entity &e) const
      { return gfs.entityOffset(e); }
      //! \brief offset of the blocks of dofs attached to a given subentity of
      //!        an element
      /**
       * This method determines the starting offset of the block of dofs
       * attached to a subentity of the given codim 0 entity.  If the grid
       * (and the ordering) directly support entities of the given
       * codimension, this is equivalent to calling
       * entityOffset(*e.subEntity<codim>(subentity)).
       */
      template<class Element>
      typename Traits::SizeType entityOffset(const Element &e, std::size_t codim,
                            std::size_t subentity) const
      { return gfs.entityOffset(e, codim, subentity); }
      //! offset of the block of dofs attached to a given intersection
      template<class Intersection>
      typename Traits::SizeType intersectionOffset(const Intersection &i) const
      { return gfs.intersectionOffset(i); }

#endif // 0

    private:

      const bool _backend_blocked;
      std::size_t _size;
      std::size_t _max_local_size;
      std::vector<typename Traits::SizeType> _gt_dof_offsets;

    };


    template<typename GridFunctionSpace, typename Params>
    Dune::PDELab::TypeTree::GenericLeafNodeTransformation<
      GridFunctionSpace,
      gfs_to_ordering<Params>,
      LeafLocalFunctionSpaceNode<GridFunctionSpace,typename gfs_to_lfs<Params>::MultiIndex>
      >
    lookupNodeTransformation(GridFunctionSpace* gfs, gfs_to_ordering<Params>* t, LeafGridFunctionSpaceTag tag);


    struct collect_a_priori_fixed_size
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp) const
      {
        node._fixed_size = node.gridFunctionSpace().alwaysFixedSize();
      }

      template<typename Node, typename TreePath>
      void pre(Node& node, TreePath tp) const
      {
        node._fixed_size = true;
      }

      template<typename Node, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(Node& node, const Child& child, TreePath tp, ChildIndex childIndex) const
      {
        node._fixed_size = node._fixed_size && child._fixed_size;
      }

    };


    template<typename GV>
    struct update_fixed_size
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      typedef std::vector<Dune::GeometryType> GTVector;

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp) const
      {
        if (node._fixed_size)
          {
            typedef typename Node::Traits::SizeType size_type;
            const size_type dim = GV::dimension;
            node._codim_used.assign(dim,false);
            node._gt_used.assign(Dune::GlobalGeometryTypeIndex::size(dim),false);
            node._gt_dof_offsets.assign(Dune::GlobalGeometryTypeIndex::size(dim),0);
            for (GTVector::const_iterator it = geom_types.begin(); it != geom_types.end(); ++it)
              {
                size_type size = node.gridFunctionSpace().size(*it);
                node._gt_dof_offsets[Dune::GlobalGeometryTypeIndex::index(*it)] = size;
                node._gt_used[Dune::GlobalGeometryTypeIndex::index(*it)] = size > 0;
                node._codim_used[cc] |= size > 0;
              }
          }
      }

      template<typename Node, typename TreePath>
      void pre(Node& node, TreePath tp) const
      {
        if (node._fixed_size)
          {
            typedef typename Node::Traits::SizeType size_type;
            const size_type dim = GV::dimension;
            node._codim_used.assign(dim,false);
            node._gt_used.assign(Dune::GlobalGeometryTypeIndex::size(dim),false);
            node._gt_dof_offsets.assign(Dune::GlobalGeometryTypeIndex::size(dim) * Node::CHILDREN,0);
          }
      }

      template<typename Node, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(Node& node, const Child& child, TreePath tp, ChildIndex childIndex) const
      {
        if (node._fixed_size)
          {
            std::transform(node._codim_used.begin(),
                           node._codim_used.end(),
                           child._codim_used.begin(),
                           node._codim_used.begin(),
                           std::logical_or<bool>());
            std::transform(node._gt_used.begin(),
                           node._gt_used.end(),
                           child._gt_used.begin(),
                           node._gt_used.begin(),
                           std::logical_or<bool>());

            typedef typename Node::Traits::SizeType size_type;

            for (size_type gt = 0; gt < Dune::GlobalGeometryTypeIndex::size(GV::dimension); ++gt)
              node._gt_dof_offsets[gt * Node::CHILDREN + childIndex] = child._gt_dof_offsets[gt * child._child_count + child._child_count - 1];
          }
      }

      template<typename Node, typename TreePath>
      void post(Node& node, TreePath tp) const
      {
        if (node._fixed_size)
          {
            typedef std::vector<typename Node::Traits::SizeType>::iterator iterator;

            iterator next_gt_it = node._gt_dof_offsets.begin() + Node::CHILDREN;
            const iterator end_it = node._gt_dof_offsets.end();

            for (iterator it = node._gt_dof_offsets.begin();
                 it != end_it;
                 it += Node::CHILDREN, next_gt_it += Node::CHILDREN)
              std::partial_sum(it,next_gt_it,it);
          }
      }

      update_fixed_size(const GV gv_, const GTVector& geom_types_)
        : gv(gv_)
        , geom_types(geom_types_)
      {}

      GV gv;
      const GTVector& geom_types;

    };


    struct pre_collect_used_geometry_types
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp) const
      {
        if (!node._fixed_size)
          {
            node._codim_used.assign(dim,false);
            node._gt_used.assign(Dune::GlobalGeometryTypeIndex::size(dim),false);
            node._gt_dof_offsets.assign(Dune::GlobalGeometryTypeIndex::size(dim) * node._child_count,0);
            node._gt_entity_offsets.assign(Dune::GlobalGeometryTypeIndex::size(dim) + 1,0);
          }
      }

      template<typename Node, typename TreePath>
      void pre(Node& node, TreePath tp) const
      {
        leaf(node,tp);
      }

      pre_collect_used_geometry_types(std::size_t dimension)
        : dim(dimension)
      {}

      const std::size_t dim;

    };


    template<typename Cell>
    struct collect_used_geometry_types_from_cell
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp) const
      {
        if (!node._fixed_size)
          {
            node.gridFunctionSpace().usedCodimsAndGeometryTypes(cell,node._codim_used,node._gt_used);
          }
      }

      collect_used_geometry_types_from_cell(const Cell& cell_)
        : cell(cell_)
      {}

      const Cell& cell;

    };


    template<typename GV>
    struct post_collect_used_geometry_types
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      typedef std::vector<Dune::GeometryType> GTVector;


      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp) const
      {
        if (!node._fixed_size)
          {
            typedef typename Node::Traits::SizeType size_type;

            for (GTVector::const_iterator it = geom_types.begin(); it != geom_types.end(); ++it)
              {
                if (node._gt_used[Dune::GlobalGeometryTypeIndex::index(*it)])
                  node._gt_entity_offsets[Dune::GlobalGeometryTypeIndex::index(*it) + 1] = gv.indexSet().size(*it);
              }

            std::partial_sum(node._gt_entity_offsets.begin(),node._gt_entity_offsets.end(),node._gt_entity_offsets.begin());
            node._entity_dof_offsets.assign(node._gt_entity_offsets.back() * node._child_count,0);
          }
      }

      template<typename Node, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(Node& node, const Child& child, TreePath tp, ChildIndex childIndex) const
      {
        if (!node._fixed_size)
          {
            std::transform(node._codim_used.begin(),
                           node._codim_used.end(),
                           child._codim_used.begin(),
                           node._codim_used.begin(),
                           std::logical_or<bool>());
            std::transform(node._gt_used.begin(),
                           node._gt_used.end(),
                           child._gt_used.begin(),
                           node._gt_used.begin(),
                           std::logical_or<bool>());
          }
      }

      template<typename Node, typename TreePath>
      void post(Node& node, TreePath tp) const
      {
        leaf(node,tp);
      }

      post_collect_used_geometry_types(const GV& gv_, const GTVector& geom_types_)
        : gv(gv_)
        , geom_types(geom_types_)
      {}

      GV gv;
      const GTVector& geom_types;

    };


    template<typename GV>
    struct extract_per_entity_sizes_from_cell
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {

      static const std::size_t dim = GV::dimension;
      typedef typename GV::template Codim<0>::Entity Cell;
      typedef std::size_t size_type;

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp) const
      {
        if (!node._fixed_size)
          {

            if (node._fixed_size_possible)
              std::fill(gt_sizes.begin(),gt_sizes.end(),0);

            ... coeffs = ...;

            for (CoeffIterator it = coeffs.begin(); it != coeffs.end(); ++it)
              {
                Dune::GeometryType gt = ref_el.type(it->subEntity(),it->codim());
                const size_type geometry_type_index = Dune::GlobalGeometryTypeIndex::index(gt);

                const size_type entity_index = gv.indexSet().subIndex(*cell,it->subEntity(),it->codim());
                const size_type index = node._gt_entity_offsets[geometry_type_index] + entity_index;
                gt_sizes[geometry_type_index] = node._entity_dof_offsets[index] = std::max(node._entity_dof_offsets[index],it->index() + 1);
              }

            îf (node._fixed_size_possible)
              {
                for (size_type i = 0; i < gt_sizes.size(); ++i)
                  if (gt_sizes[i] > 0)
                    {
                      if (node._gt_dof_offsets[i + 1] == 0)
                        node._gt_dof_offsets[i + 1] = gt_sizes[i];
                      else if (node._gt_dof_offsets[i + 1] != gt_sizes[i])
                        {
                          node._fixed_size_possible = false;
                          break;
                        }
                    }
              }

          }
      }

      extract_per_entity_sizes_from_cell(const GV& gv_)
        : gv(gv_)
        , cell(nullptr)
        , ref_el(Dune::GenericReferenceElements<typename Entity::ctype,dim>::general(cell.type()))
        , gt_sizes(Dune::GlobalGeometryTypeIndex::size(dim),0)
      {}

      GV gv;
      const Cell* cell;
      const Dune::GenericReferenceElement<typename Entity::ctype,dim>& ref_el;
      std::vector<size_type> gt_sizes;

    };


    template<typename GV>
    struct post_extract_per_entity_sizes
      : public TypeTree::TreeVisitor
      , public TypeTree::DynamicTraversal
    {


      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath tp) const
      {
        if (!node._fixed_size)
          {
            if (node._fixed_size_possible)
              {
                std::swap(node._entity_dof_offsets,typename Node::SizeVector());
                node._fixed_size = true;
              }
          }
      }

      template<typename Node, typename TreePath>
      void pre(Node& node, TreePath tp) const
      {
        if (!node._fixed_size)
          node._fixed_size_possible = true;
      }


      template<typename Node, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(Node& node, const Child& child, TreePath tp, ChildIndex childIndex) const
      {
        if (!node._fixed_size)
          node._fixed_size_possible = node._fixed_size_possible && child._fixed_size;
      }


      template<typename Node, typename TreePath>
      void post(Node& node, TreePath tp) const
      {
        if (!node._fixed_size)
          {
            if (node._fixed_size_possible)
              {
                for (size_type gt = 0; gt < Dune::GlobalGeometryTypeIndex::size(GV::dimension); ++gt)
                  {
                    size_type carry = 0;
                    for (size_type child_index = 0; child_index < Node::CHILDREN; ++child_index)
                      node._gt_dof_offsets[gt * Node::CHILDREN + childIndex] = (carry += child._gt_dof_offsets[gt * child._child_count + child._child_count - 1]);
                  }
                node._fixed_size = true;
              }
            else
              {
                typedef typename Node::Traits::SizeType size_type;

                for (GTVector::const_iterator it = geom_types.begin(); it != geom_types.end(); ++it)
                  {
                    if (node._gt_used[Dune::GlobalGeometryTypeIndex::index(*it)])
                      node._gt_entity_offsets[Dune::GlobalGeometryTypeIndex::index(*it) + 1] = gv.indexSet().size(*it);
                  }

                std::partial_sum(node._gt_entity_offsets.begin(),node._gt_entity_offsets.end(),node._gt_entity_offsets.begin());
                node._entity_dof_offsets.assign(node._gt_entity_offsets.back() * node._child_count,0);

                for (size_type geometry_type_index = 0; geometry_type_index < Dune::GlobalGeometryTypeIndex::size(dim); ++geometry_type_index)
                  {
                    if (!node._gt_used[geometry_type_index])
                      continue;
                    for (size_type entity_index = 0; entity_index < node._gt_entity_offsets[geometry_type_index+1]; ++entity_index)
                      {
                        size_type carry = 0;
                        for (size_type child_index = 0; child_index < Node::CHILDREN; ++child_index)
                          node._entity_dof_offsets[node._gt_entity_offsets[geometry_type_index] + entity_index] = (carry += node._children[child_index].size(geometry_type_index,entity_index));
                      }
                  }
              }
          }
      }

      post_collect_used_geometry_types(const GV& gv_, const GTVector& geom_types_)
        : gv(gv_)
        , geom_types(geom_types_)
      {}

      GV gv;
      const GTVector& geom_types;

    };


    template<typename GV, typename LocalOrdering>
    class NonLeafGridViewOrdering
      : public TypeTree::VariadicCompositeNode<LocalOrdering>
      , public VirtualOrderingBase
    {
    public:
      typedef typename LocalOrdering::Traits::SizeType SizeType;

    private:

      typedef TypeTree::VariadicCompositeNode<LocalOrdering> NodeT;

      GV _gv;

    public:
      //! Construct ordering object
      /**
       * In general, an ordering object is not properly setup after
       * construction.  This must be done by a seperate call to update().
       * This particular ordering however can be used right away.
       */
      NonLeafGridViewOrdering(const GV& gv, const NodeT::NodeStorage& local_ordering)
        : NodeT(local_ordering)
        , _gv(gv)
      {}

      virtual void map_index_dynamic(const MultiIndex& mi, ContainerIndex& ci) const
      {
        map_index(mi,ci);
      }

      void map_index(const MultiIndex& mi, ContainerIndex& ci) const
      {
        const SizeType geometry_type_index = ...;
        const SizeType entity_index = ...;
        const SizeType child_index = mi.back();
        localOrdering().map_local_index(geometry_type_index,entity_index,mi,ci);
        if (_backend_blocked)
          {
            ci.push_back(_gt_entity_offsets[geometry_type_index] + entity_index);
          }
        else if (_fixed_size)
          {
            ci.back() += _gt_dof_offsets[geometry_type_index] + entity_index * _gt_dof_size[geometry_type_index];
          }
        else
          {
            ci.back() += _entity_dof_offsets[_gt_entity_offsets[geometry_type_index] + entity_index];
          }
      }

      void update()
      {

        GTVector geom_types;

        for (size_type cc = 0; cc <= dim; ++cc)
          {
            const GTVector& per_codim_geom_types = _gv.indexSet().geometryTypes(cc);
            std::copy(per_codim_geom_types.begin(),per_codim_geom_types.end(),std::back_inserter(geom_types));
          }

        // Do we already know that we have fixed per-GeometryType sizes?
        TypeTree::applyToTree(localOrdering(),collect_a_priori_fixed_size());
        _fixed_size = localOrdering().fixedSize();

        typedef std::vector<Dune::GeometryType> GTVector;
        const Traits::SizeType gt_index_count = Dune::GlobalGeometryTypeIndex::size(GV::dimension);

        if (_fixed_size)
          {
            // collect used GeometryTypes
            TypeTree::applyToTree(localOrdering(),update_fixed_size<GV>(_gv,geom_types));

            _gt_dof_offsets.resize(gt_index_count + 1);

            const GTVector::const_iterator end_it = geom_types.end();
            for (GTVector::const_iterator it = geom_types.begin(); it != end_it; ++it)
              {
                const Traits::SizeType gt_index = Dune::GlobalGeometryTypeIndex::index(*it);
                _gt_dof_offsets[gt_index + 1] = localOrdering().size(gt_index,0) * _gv.indexSet().size(*it);
              }
            std::partial_sum(_gt_dof_offsets.begin(),_gt_dof_offsets.end(),_gt_dof_offsets.begin());
          }
        else
          {
            TypeTree::applyToTree(localOrdering(),pre_collect_used_geometry_types(GV::dimension));
            const CellIterator end_it = _gv.template end<0>();
            for (CellIterator it = _gv.template begin<0>(); it != end_it; ++it)
              {
                TypeTree::applyToTree(localOrdering(),collect_used_geometry_types_from_cell<Cell>(*it));
              }
            TypeTree::applyToTree(localOrdering(),post_collect_used_geometry_types<GV>(_gv,geom_types));
            // allocate

            TypeTree::applyToTree(localOrdering(),pre_extract_per_entity_sizes<GV>(_gv));
            for (CellIterator it = _gv.template begin<0>(); it != end_it; ++it)
              {
                TypeTree::applyToTree(localOrdering(),extract_per_entity_sizes_from_cell<Cell>(*it));
              }
            TypeTree::applyToTree(localOrdering(),post_extract_per_entity_sizes<GV>(_gv,geom_types));

            if (localOrdering().fixedSize())
              {
                _fixed_size = true;
                _gt_dof_offsets.resize(gt_index_count + 1);

                for (GTVector::const_iterator it = geom_types.begin(); it != geom_types.end(); ++it)
                  {
                    const Traits::SizeType gt_index = Dune::GlobalGeometryTypeIndex::index(*it);
                    _gt_dof_offsets[gt_index + 1] = localOrdering().size(gt_index,0) * _gv.indexSet().size(*it);
                  }

                std::partial_sum(_gt_dof_offsets.begin(),_gt_dof_offsets.end(),_gt_dof_offsets.begin());
              }
            else
              {
                _gt_entity_offsets.assign(gt_index_count + 1,0);

                for (GTVector::const_iterator it = geom_types.begin(); it != geom_types.end(); ++it)
                  {
                    if (!localOrdering().contains(*it))
                      continue;
                    const Traits::SizeType gt_index = Dune::GlobalGeometryTypeIndex::index(*it);
                    _gt_entity_offsets[gt_index + 1] = _gv.indexSet().size(*it);
                  }

                std::partial_sum(_gt_entity_offsets.begin(),_gt_entity_offsets.end(),_gt_entity_offsets.begin());
                _entity_dof_offsets.assign(_gt_entity_offsets.back()+1,0);
                Traits::SizeType carry = 0;
                for (GTVector::const_iterator it = geom_types.begin(); it != geom_types.end(); ++it)
                  {
                    if (!localOrdering().contains(*it))
                      continue;
                    const Traits::SizeType gt_index = Dune::GlobalGeometryTypeIndex::index(*it);
                    Traits::SizeType entity_pos = _gt_entity_offsets[gt_index] + 1;
                    const Traits::SizeType entity_count = _gt_entity_offsets[gt_index + 1] - entity_pos + 1;
                    for (Traits::SizeType entity_index = 0; entity_index < entity_count; ++entity_index, ++entity_pos)
                      _entity_dof_offsets[entity_pos] = (carry += localOrdering.size(gt_index,entity_index));
                  }

              }
          }
      }
    };


   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_LEAFORDERING_HH
