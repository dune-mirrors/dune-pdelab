// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_LEAFGRIDVIEWORDERING_HH
#define DUNE_PDELAB_ORDERING_LEAFGRIDVIEWORDERING_HH

#include <dune/pdelab/ordering/directleaflocalordering.hh>
#include <dune/pdelab/ordering/leaforderingbase.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

    //! Gridview ordering for leaf spaces
    template<typename LocalOrdering>
    class LeafGridViewOrdering
      : public LeafOrderingBase<LocalOrdering>
    {
    public:
      typedef typename LocalOrdering::Traits Traits;

    private:

      using ES = typename Traits::EntitySet;

      typedef LeafOrderingBase<LocalOrdering> BaseT;
      typedef typename BaseT::NodeT NodeT;

    public:

      static std::string name() {return "LeafGridViewOrdering";}

      LeafGridViewOrdering(const typename NodeT::NodeStorage& local_ordering, typename BaseT::GFSData* gfs_data)
        : BaseT(local_ordering, true, gfs_data)
        , _es(this->template child<0>().entitySet())
      {}

#ifndef DOXYGEN

// we need to override the default copy / move ctor to fix the delegate pointer, but that is
// hardly interesting to our users...

      LeafGridViewOrdering(const LeafGridViewOrdering& r)
        : BaseT(r)
        , _es(r._es)
      {}

      LeafGridViewOrdering(LeafGridViewOrdering&& r)
        : BaseT(std::move(r))
        , _es(r._es)
      {}

#endif // DOXYGEN

      virtual ~LeafGridViewOrdering() override = default;


      using BaseT::size;

      const ES& entitySet() const {
        return _es;
      }

      /**
       * @brief Gives the size for a given suffix
       * @param suffix  MultiIndex with a partial path to a container
       * @return Traits::SizeType  The size required for such a path.
       */
      typename Traits::SizeType size(typename Traits::ContainerIndex suffix) const
      {
        // notice that this algorithm is the same as for GridViewOrdering
        // the only difference is that here we can borrow some offsets from the
        // local ordering

        using size_type = typename Traits::SizeType;
        if (suffix.size() == Traits::ContainerIndex::max_depth)
          return 0; // all indices in suffix were consumed, no more sizes to provide
        if (suffix.size() == 0) // suffix wants the size of this depth
          return _block_count;

        // we first have to figure out the entity index
        typename Traits::DOFIndex::EntityIndex entity_index;

        // the next index to find out its size
        auto back_index = suffix.back();

        // borrow offsets & fixed_size from local ordering
        auto& _gt_entity_offsets = this->localOrdering()._gt_entity_offsets;
        auto& _entity_dof_offsets = this->localOrdering()._entity_dof_offsets;
        bool _fixed_size = this->localOrdering()._fixed_size;
        // we just need to make the inverse computation of the mapIndex funtion to find the entity index
        if (this->localOrdering()._container_blocked)
          suffix.pop_back();

        auto dof_begin = _fixed_size ? _gt_dof_offsets.begin() : _entity_dof_offsets.begin();
        auto dof_end = _fixed_size ? _gt_dof_offsets.end() : _entity_dof_offsets.end();
        auto dof_it = std::prev(std::upper_bound(dof_begin, dof_end, back_index));
        size_type dof_dist = std::distance(dof_begin, dof_it);
        if (_fixed_size) {
          // On fixed size, entity index is not used down the tree. Set max to trigger segfault if this does not hold.
          Traits::DOFIndexAccessor::GeometryIndex::store(entity_index,dof_dist,~size_type{0});
        } else {
          auto gt_begin = _gt_entity_offsets.begin();
          auto gt_end = _gt_entity_offsets.end();
          auto gt_it = std::prev(std::upper_bound(gt_begin, gt_end, dof_dist));
          size_type gt = std::distance(gt_begin, gt_it);
          assert(dof_dist >= *gt_it);
          size_type ei = dof_dist - *gt_it;
          Traits::DOFIndexAccessor::GeometryIndex::store(entity_index,gt,ei);
        }

        // then, the local ordering knows the size for a given entity.
        return this->localOrdering().size(suffix, entity_index);
      }


      virtual void update() override
      {
        LocalOrdering& lo = this->localOrdering();
        lo.update_a_priori_fixed_size();

        const std::size_t dim = ES::dimension;

        typename ES::CodimMask codims;
        codims.set(0); // always need cells
        lo.collect_used_codims(codims);

        for (typename ES::dim_type codim = 0; codim <= ES::dimension; ++codim)
          if (codims.test(codim))
            _es.addCodim(codim);

        _es.update();

        typedef typename Traits::SizeType size_type;
        auto geom_types = _es.indexSet().types();

        if (lo._fixed_size)
          {
            lo.update_fixed_size(geom_types.begin(),geom_types.end());
          }
        else
          {
            lo.pre_collect_used_geometry_types_from_cell();

            for (const auto& element : elements(_es))
              {
                lo.collect_used_geometry_types_from_cell(element);
              }

            lo.allocate_entity_offset_vector(geom_types.begin(),geom_types.end());

            for (const auto& element : elements(_es))
              {
                lo.extract_per_entity_sizes_from_cell(element);
              }

            // FIXME: handling of blocked containers!
            lo.finalize_non_fixed_size_update();
          }

        // we need to re-test here, as the local ordering could have detected a fixed size ordering
        // and switched its implementation
        if (lo._fixed_size)
          {
            _gt_dof_offsets.assign(GlobalGeometryTypeIndex::size(dim) + 1,0);
            _size = 0;

            for (auto gt : geom_types)
              {
                const size_type gt_index = GlobalGeometryTypeIndex::index(gt);
                size_type gt_size = lo.size(gt_index,0);
                size_type entity_count = _es.indexSet().size(gt);
                _size += gt_size * entity_count;
                if (this->localOrdering()._container_blocked)
                  gt_size = gt_size > 0;
                _gt_dof_offsets[gt_index + 1] = gt_size * entity_count;
              }
            std::partial_sum(_gt_dof_offsets.begin(),_gt_dof_offsets.end(),_gt_dof_offsets.begin());
            _block_count = _gt_dof_offsets.back();
            _codim_fixed_size.set();
          }
        else
          {
            _size = lo._entity_dof_offsets.back();
            if (lo._container_blocked)
              _block_count = lo._gt_entity_offsets.back();
            else
              _block_count = _size;
            _codim_fixed_size.reset();
          }

        _fixed_size = lo._fixed_size;
        _max_local_size = lo.maxLocalSize();

        _codim_used = lo._codim_used;
        _codim_fixed_size = lo._codim_fixed_size;

      }

      using BaseT::fixedSize;

    private:

      using BaseT::_max_local_size;
      using BaseT::_size;
      using BaseT::_block_count;
      using BaseT::_fixed_size;
      using BaseT::_codim_used;
      using BaseT::_codim_fixed_size;
      using BaseT::_gt_dof_offsets;

      typename Traits::EntitySet _es;
    };


    template<typename GFS, typename Transformation>
    struct direct_leaf_gfs_to_gridview_ordering_descriptor
    {

      static const bool recursive = false;

      typedef DirectLeafLocalOrdering<typename GFS::Traits::OrderingTag,
                                      typename GFS::Traits::FiniteElementMap,
                                      typename GFS::Traits::EntitySet,
                                      typename Transformation::DOFIndex,
                                      typename Transformation::ContainerIndex
                                      > LocalOrdering;

      typedef LeafGridViewOrdering<LocalOrdering> GridViewOrdering;

      typedef GridViewOrdering transformed_type;
      typedef std::shared_ptr<transformed_type> transformed_storage_type;

      static transformed_type transform(const GFS& gfs, const Transformation& t)
      {
        auto local_ordering = std::make_shared<LocalOrdering>(gfs.finiteElementMapStorage(),gfs.entitySet(), gfs.backend().blocked(gfs));
        return transformed_type(make_tuple(local_ordering),const_cast<GFS*>(&gfs));
      }

      static transformed_storage_type transform_storage(std::shared_ptr<const GFS> gfs, const Transformation& t)
      {
        auto local_ordering = std::make_shared<LocalOrdering>(gfs->finiteElementMapStorage(),gfs->entitySet(),gfs->backend().blocked(*gfs));
        return std::make_shared<transformed_type>(make_tuple(local_ordering),const_cast<GFS*>(gfs.get()));
      }

    };


    template<typename GFS, typename Transformation, typename Params>
    direct_leaf_gfs_to_gridview_ordering_descriptor<GFS,Transformation>
    register_leaf_gfs_to_ordering_descriptor(GFS*,Transformation*,LeafOrderingTag<Params>*);

   //! \} group Ordering
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_LEAFGRIDVIEWORDERING_HH
