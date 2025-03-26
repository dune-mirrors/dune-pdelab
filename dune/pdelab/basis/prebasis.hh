#ifndef DUNE_PDELAB_BASIS_PREBASIS_HH
#define DUNE_PDELAB_BASIS_PREBASIS_HH

#include <dune/pdelab/basis/topologicassociativityforest.hh>

#include <dune/pdelab/basis/nodes.hh>

namespace Dune::PDELab {

template<class ProtoBasisTree, class GV>
class PreBasis {
  using TopologicAssociativityForest = Impl::TopologicAssociativityForest<ProtoBasisTree, GV>;
  using Element = GV::Traits::template Codim<0>::Entity;
public:
  using size_type = typename TopologicAssociativityForest::SizeType;
  static constexpr size_type maxMultiIndexSize = TopologicAssociativityForest::maxContainerDepth();
  static constexpr size_type minMultiIndexSize = TopologicAssociativityForest::maxContainerDepth();
  static constexpr size_type multiIndexBufferSize = TopologicAssociativityForest::maxContainerDepth();
  using SizePrefix = Dune::PDELab::MultiIndex<size_type,maxMultiIndexSize>;
  using Node = LocalBasisViewTree<Element, ProtoBasisTree>;
  using GridView = GV;

  PreBasis(const ProtoBasisTree& proto_basis, GridView grid_view)
    : taf_{Impl::makeTopologicAssociativityForest(proto_basis, grid_view)}
  {}

  GridView gridView() const {
    return taf_.entitySet();
  }

  Node makeNode() const {
    return makeLocalBasisView<Element>(taf_.protoBasis());
  }

  size_type size() const {
    return dimension();
  }

  size_type size(auto multi_index) const {
    // Note: Multi-index is read Outer->Inner
    return containerSize(reverse(multi_index));
  }

  [[nodiscard]] size_type dimension() const noexcept {
    if (taf_.fixedSize())
      return _gt_dof_offsets->back();
    else
      return _entity_dof_offsets->back();
  }

  size_type maxNodeSize() const {
    return taf_.maxLocalCount();
  }

  template<typename It>
  It indices(const Node& node, It it) const {

    std::array<std::size_t,Node::Element::dimension+2> _codim_offsets;
    std::vector<std::array<size_type,2>> _index_cache; // TODO: move to the node

    // find the sub entities with topological associativity
    if (taf_.disjointCodimClosure()) {
      // DG/FV branch
      size_type gt_index = GlobalGeometryTypeIndex::index(node.element().type());
      size_type entity_index = gridView().indexSet().index(node.element());
      _index_cache.assign(1, {gt_index, entity_index});
    } else {
      std::fill(std::begin(_codim_offsets), std::end(_codim_offsets), 0);
      Hybrid::forEach(range(index_constant<Node::Element::dimension+1>{}), [&](auto codim) {
        if (codim == 0 or (TopologicAssociativityForest::mayContainCodim(codim) and taf_.containsCodim(codim))) {
          std::size_t sub_entities = node.element().subEntities(codim);
          _index_cache.resize(_index_cache.size() + sub_entities);
          std::size_t codim_offset = _codim_offsets[codim];
          for (std::size_t s = 0; s != sub_entities; ++s) {
            const auto& sub_entity = node.element().template subEntity<codim>(s);
            size_type gt_index = GlobalGeometryTypeIndex::index(sub_entity.type());
            size_type entity_index = gridView().indexSet().index(sub_entity);
            _index_cache[codim_offset + s] = {gt_index, entity_index};
          }
        }
        _codim_offsets[codim+1] = _index_cache.size();
      });
    }

    PDELab::forEachLeafNode(node, [&]<class LeafNode>(LeafNode& leaf_node, const auto& path) {

      const auto& fem = taf_.protoBasis(path).finiteElementMap();
      leaf_node.bindFiniteElement(fem.find(leaf_node.element()));
      const auto& fe = leaf_node.finiteElement();
      assert(fe.type() == leaf_node.element().type());

      if (leaf_node.disjointCodimClosure()) {
        // DG/FV branch
        const auto [gt_index, entity_index] = _index_cache[0];
        assert(fe.size() == taf_.blockCount(gt_index, entity_index));
        // In this case, we know that indices are contiguous and can be
        //   reconstructed from the local index. Thus, there is no need to
        //   fill the whole vector.
        auto cir = containerIndexRange(path, gt_index, entity_index);
        // notice that the other indices are garbage!
        // but we need them for the size to be correct
        if constexpr (not std::decay_t<LeafNode>::optimizeFastDG()) {
          // This the observed behavior for non-optimizaed cases:
          const auto& coeffs = fe.localCoefficients();
          assert(coeffs.localKey(0).index() == 0 && "Indices in DG/FV elements shall be ordered");
          for (std::size_t dof = 0; dof != fe.size(); ++dof) {
            assert(coeffs.localKey(dof).index() == dof && "Indices in DG/FV elements shall be ordered");
            *it = cir[dof];
            ++it;
          }
        }
      } else {
        const auto& coeffs = fe.localCoefficients();
        // for each local dof, we need to figure out its assigned container index
        for (std::size_t dof = 0; dof < coeffs.size(); ++dof) {
          const auto& key = coeffs.localKey(dof);
          assert(taf_.containsCodim(key.codim()));
          const auto [gt_index, entity_index] = [&]{
            auto offset = _codim_offsets[key.codim()];
            offset += key.subEntity();
            assert(offset < _index_cache.size());
            return _index_cache[offset];
          }();
          auto cir = containerIndexRange(path, gt_index, entity_index);
          assert(key.index() < taf_.blockCount(gt_index, entity_index));
          *it = cir[key.index()];
          ++it;
        }
      }
    });
    return it;
  }

  void update(GridView grid_view) {
    taf_.update(grid_view);
  }

  void initializeIndices()
  {
    taf_.initializeIndices();
    if (taf_.fixedSize())
      updateFixedSizeOrdering();
    else
      updateVariableSizeOrdering();
  }

  // template<Entity>
  // auto containerIndexRange(auto path, const Entity& entity) const noexcept;

private:

  static constexpr auto containerBlocked() {
    return std::integral_constant<bool,TopologicAssociativityForest::MergingStrategy::Blocked>{};
  }

  auto containerIndexRange(auto path, const size_type& gt_index, const size_type& e_index) const noexcept {
    // Note: Multi-index is read Outer->Inner
    // get local container suffix
    assert(taf_.containsGeometry(gt_index));
    const auto lcir = taf_.containerIndexRange(path, gt_index, e_index);
    auto offset = blockOffset(gt_index, e_index);

    if constexpr (containerBlocked())
      return transformedRangeView(lcir, [](auto mi){ return push_front(mi, offset); });
    else
      return transformedRangeView(lcir, [](auto mi){ return accumulate_front(mi, offset); });
  }

  auto getGeometryEntityIndex(size_type i) const {
    // we just need to make the inverse computation of the mapIndex funtion to find the entity index
    size_type gt_index, entity_index;
    auto dof_begin = std::begin(taf_.fixedSize() ? _gt_dof_offsets : _entity_dof_offsets);
    auto dof_end = std::end(taf_.fixedSize() ? _gt_dof_offsets : _entity_dof_offsets);
    auto dof_it = std::prev(std::upper_bound(dof_begin, dof_end, i));
    auto dof_dist = std::distance(dof_begin, dof_it);

    if (taf_.fixedSize()) {
      gt_index = dof_dist;
      entity_index = i - *dof_it;
      entity_index /= taf_.blockCount(gt_index);
    } else {
      auto gt_begin = std::begin(_gt_entity_offsets);
      auto gt_end = std::end(_gt_entity_offsets);
      auto gt_it = std::prev(std::upper_bound(gt_begin, gt_end, dof_dist));
      gt_index = std::distance(gt_begin, gt_it);
      entity_index = dof_dist - *gt_it;
    }
    return std::make_pair(gt_index, entity_index);
  }


  size_type blockOffset(size_type geometry_index, size_type entity_index) const {
    if (taf_.fixedSize()) {
      const auto& gt_offset = _gt_dof_offsets[geometry_index];
      if constexpr (containerBlocked())
        return gt_offset + entity_index;
      else
        return gt_offset + entity_index * taf_.blockCount(geometry_index);
    } else {
      const auto& gt_offset = _gt_entity_offsets[geometry_index];
      return _entity_dof_offsets[gt_offset + entity_index];
    }
  }


  template<Concept::MultiIndex ContainerSuffix>
  std::size_t containerSize(ContainerSuffix suffix) const noexcept
  {
    // Note: Multi-index is read Inner->Outer
    if (suffix.size() == 0)
      return dimension();
    else {
      // the next index to find out its size
      auto back_index = suffix.back();

      // we first have to figure out the entity index
      const auto [gt_index, entity_index] = getGeometryEntityIndex(back_index);

      // then, the entity ordering knows the size for a given entity.
      if constexpr (containerBlocked())
        return taf_.containerSize(pop_back(suffix), gt_index, entity_index);
      else {
        // remove offset introduced by the entity_index
        suffix.back() -= blockOffset(gt_index, entity_index);
        return taf_.containerSize(suffix, gt_index, entity_index);
      }
    }
  }

  void updateFixedSizeOrdering() {
    assert(taf_.fixedSize());

    const size_type gt_count = GlobalGeometryTypeIndex::size(GridView::dimension);
    std::tie(_gt_dof_offsets_storage, _gt_dof_offsets) = Impl::make_span_storage(gt_count + 1, size_type{0});

    for (std::size_t codim = 0; codim <= GridView::dimension; ++codim) {
      for (const auto& gt : gridView().indexSet().types(codim)) {
        const size_type gt_index = GlobalGeometryTypeIndex::index(gt);
        if (not taf_.containsGeometry(gt_index))
          continue;
        size_type gt_size = taf_.blockCount(gt_index);
        const size_type gt_entity_count = gridView().indexSet().size(gt);
        if (containerBlocked())
          gt_size = (gt_size > 0);
        _gt_dof_offsets[gt_index + 1] = gt_size * gt_entity_count;
      }
    }

    std::partial_sum(std::begin(_gt_dof_offsets), std::end(_gt_dof_offsets), std::begin(_gt_dof_offsets));
  }


  void updateVariableSizeOrdering() {
    // in this case, we need to count the number of _used_ sub blocks
    assert(not taf_.fixedSize());
    const size_type gt_count = GlobalGeometryTypeIndex::size(GridView::dimension);
    std::tie(_gt_entity_offsets_storage, _gt_entity_offsets) = Impl::make_span_storage(gt_count + 1, size_type{0});

    for (std::size_t codim = 0; codim <= GridView::dimension; ++codim) {
      for (const auto& gt : gridView().indexSet().types(codim)) {
        if (!taf_.containsGeometry(gt))
          continue;
        const size_type gt_index = GlobalGeometryTypeIndex::index(gt);
        _gt_entity_offsets[gt_index + 1] = gridView().indexSet().size(gt);
      }
    }

    std::partial_sum(std::begin(_gt_entity_offsets), std::end(_gt_entity_offsets), std::begin(_gt_entity_offsets));
    std::tie(_entity_dof_offsets_storage, _entity_dof_offsets) = Impl::make_span_storage(_gt_entity_offsets->back()+1, size_type{0});

    size_type carry_block = 0;
    size_type index = 0;
    for (size_type gt_index = 0; gt_index < gt_count; ++gt_index) {
      if (not taf_.containsGeometry(gt_index))
        continue;
      const size_type entity_count = _gt_entity_offsets[gt_index + 1] - _gt_entity_offsets[gt_index];
      for (size_type entity_index = 0; entity_index < entity_count; ++entity_index) {
        const size_type size = taf_.blockCount(gt_index, entity_index);
        carry_block += (containerBlocked() ? (size > 0) : size);
        _entity_dof_offsets[++index] = carry_block;
      }
    }
  }

  TopologicAssociativityForest taf_;

  std::span<size_type> _gt_dof_offsets;
  std::span<size_type> _gt_entity_offsets;
  std::span<size_type> _entity_dof_offsets;

  std::unique_ptr<size_type[]> _gt_dof_offsets_storage;
  std::unique_ptr<size_type[]> _gt_entity_offsets_storage;
  std::unique_ptr<size_type[]> _entity_dof_offsets_storage;
};

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_BASIS_PREBASIS_HH
