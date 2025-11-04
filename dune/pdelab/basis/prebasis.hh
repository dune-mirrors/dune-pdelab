#ifndef DUNE_PDELAB_BASIS_PREBASIS_HH
#define DUNE_PDELAB_BASIS_PREBASIS_HH

#include <dune/pdelab/basis/localview.hh>
#include <dune/pdelab/basis/utility.hh>

#include <dune/functions/functionspacebases/containerdescriptors.hh>

namespace Dune::PDELab
{

  // Forward declaration for PreBasis
  template <class ProtoBasis>
  class PreBasis;

  /**
   * @brief Node type for PreBasis, representing a local basis view tree bound to a grid entity.
   *
   * This class extends LocalBasisViewTree and adds caching and offset management
   * for efficient binding and traversal of local basis functions on a grid entity.
   *
   * @tparam ProtoBasis The proto-basis type used to construct the local basis view tree.
   */
  template <class ProtoBasis>
  class PreBasisNode : public LocalBasisViewTree<ProtoBasis>
  {
    using GridView = typename ProtoBasis::GridView;

  public:
    using size_type = typename ProtoBasis::size_type;
    using Element = typename GridView::template Codim<0>::Entity;

    /**
     * @brief Constructs a PreBasisNode from a local basis view tree and proto-basis.
     * @param super The local basis view tree to use as the base.
     * @param proto_basis The proto-basis used for construction.
     */
    PreBasisNode(LocalBasisViewTree<ProtoBasis> &&super, const ProtoBasis &proto_basis);

    /**
     * @brief Binds the node to a specific grid entity, updating indices and offsets.
     * @param entity The entity to bind to.
     */
    void bind(const Element &entity) noexcept;


    /**
     * @brief Returns a range of entity indices for all basis functions with topological association on sub-entities of the bound entity.
     * @return A span of pairs representing (geometry type index, entity index).
     */
    std::span<const std::array<size_type, 2>> subEntityIndices() const noexcept {
      return std::span{index_cache_}.subspan(0, codim_offsets_.back());
    }

  private:
    friend PreBasis<ProtoBasis>;

    std::vector<std::array<size_type, 2>> index_cache_;             ///< Cache for entity indices
    std::array<std::size_t, Element::dimension + 2> codim_offsets_; ///< Offsets for codimensions
    const ProtoBasis &_proto_basis;                                 ///< Reference to the proto-basis
  };

  /**
   * @ingroup FunctionSpaceBases
   * @brief Pre-basis for a given proto-basis.
   *
   * Implements the interface of a pre-basis using the details from a proto-basis tree
   * (@ref ProtoBasisLeaf, @ref ProtoBasisArray, @ref ProtoBasisVector, or @ref ProtoBasisTuple).
   * Pre-basis trees can be composed and transformed into global bases using dune-functions.
   * They manage merging strategies and ordering for local-to-global index mappings.
   *
   * In essence, proto-basis trees control the merging strategy concerning ordering with respect to individual entities of the grid view,
   * whereas pre-basis trees control the merging merging concerning generic transformations of the ordering.
   *
   * @tparam ProtoBasis_ The proto-basis type providing finite element mapping and indexation.
   */
  template <class ProtoBasis_>
  class PreBasis
  {
    //! Maximum count of geometry types
    static constexpr std::size_t GT_COUNT = GlobalGeometryTypeIndex::size(ProtoBasis_::GridView::dimension);

  public:
    //! Proto-basis type
    using ProtoBasis = ProtoBasis_;
    //! Type to represent offsets and sizes
    using size_type = typename ProtoBasis::size_type;
    //! Maximum size of multi-indices
    static constexpr size_type maxMultiIndexSize = ProtoBasis::maxMultiIndexSize() + (ProtoBasis::isIndexBlocked ? 1 : 0);
    //! Minimum size of multi-indices
    static constexpr size_type minMultiIndexSize = ProtoBasis::minMultiIndexSize() + (ProtoBasis::isIndexBlocked ? 1 : 0);
    //! Minimum size of multi-indices during their construction
    static constexpr size_type multiIndexBufferSize = maxMultiIndexSize;
    //! Multi-index type for container index prefixes
    using SizePrefix = Dune::PDELab::MultiIndex<size_type, maxMultiIndexSize>;
    //! Grid view type
    using GridView = typename ProtoBasis::GridView;

    using Node = PreBasisNode<ProtoBasis>;

    /**
     * @brief Constructs a PreBasis from a proto-basis.
     * @param proto_basis The proto-basis to use.
     */
    PreBasis(const ProtoBasis &proto_basis);

    /**
     * @brief Returns the associated grid view.
     */
    const GridView &gridView() const;

    /**
     * @brief Creates a new PreBasisNode for this pre-basis.
     */
    Node makeNode() const;

    /**
     * @brief Returns the total size (number of degrees of freedom).
     */
    size_type size() const;

    /**
     * @brief Returns the size for a given multi-index prefix.
     * @param size_prefix The prefix to query.
     */
    size_type size(auto size_prefix) const;

    /**
     * @brief Returns the dimension of the underlying proto-basis.
     */
    [[nodiscard]] size_type dimension() const noexcept;

    /**
     * @brief Returns the maximum size of a node in the tree.
     */
    size_type maxNodeSize() const;

    /**
     * @brief Fills container indices for a node using an iterator.
     * @param node The node to query.
     * @param it Output iterator for indices.
     */
    template <typename It>
    It indices(const Node &node, It it) const;

    /**
     * @brief Updates the grid view for all leaf nodes.
     * @param grid_view The new grid view.
     */
    void update(const GridView& grid_view);

    /**
     * @brief Initializes index offsets for the pre-basis.
     */
    void initializeIndices();

    /**
     * @brief Returns a container descriptor for the entire pre-basis.
     *
     * The returned descriptor describes how the pre-basis coefficients can be organized in memory.
     */
    auto containerDescriptor() const {

      const auto& entity_it = gridView().template begin<0>();
      size_type gt_index0{0}, e_index0{0};
      if (entity_it != gridView().template end<0>()) {
        gt_index0 = GlobalGeometryTypeIndex::index(entity_it->type());
        e_index0 = _proto_basis.gridView().indexSet().index(*entity_it);
      }
      auto descriptor0 = protoBasis().containerDescriptor(gt_index0, e_index0);

      using ProtoBasisDescriptor = std::remove_cvref_t<decltype(descriptor0)>;
      using namespace Dune::PDELab::ContainerDescriptors;
      if constexpr (std::same_as<ProtoBasisDescriptor, Unknown>)
        return Unknown{};
      else if constexpr (ProtoBasis::isIndexBlocked) {
        if constexpr (IsCompileTimeUniform<ProtoBasisDescriptor>)
          return UniformVector<ProtoBasisDescriptor>(size(), descriptor0);
        else {
          Vector<ProtoBasisDescriptor> vector;
          for (auto&& element : elements(gridView())) {
            size_type gt_index = GlobalGeometryTypeIndex::index(element.type());
            size_type e_index = _proto_basis.gridView().indexSet().index(element);
            auto desc = protoBasis().containerDescriptor(gt_index, e_index);
            vector.emplace_back(std::move(desc));
          }
          return vector;
        }
      } else {
        using BlockType = BlockType<ProtoBasisDescriptor>;
        if constexpr (std::same_as<BlockType, Unknown>)
          return Unknown{};
        else if constexpr (IsCompileTimeUniform<ProtoBasisDescriptor>)
          return UniformVector<BlockType>(size(), descriptor0[0]);
        else {
          Vector<BlockType> vector;
          for (auto&& element : elements(gridView())) {
            size_type gt_index = GlobalGeometryTypeIndex::index(element.type());
            size_type e_index = _proto_basis.gridView().indexSet().index(element);
            auto desc = protoBasis().containerDescriptor(gt_index, e_index);
            Hybrid::forEach(range(Hybrid::size(desc)), [&](auto i){
              vector.emplace_back(std::move(desc[i]));
            });
          }
          return vector;
        }
      }
    }

    /**
     * @brief Returns the underlying proto-basis.
     */
    const ProtoBasis &protoBasis() const;

  private:
    /// Returns whether the pre-basis is merged.
    static constexpr auto merged();

    auto containerIndexRange(auto path, const size_type &gt_index, const size_type &e_index) const noexcept;

    auto getGeometryEntityIndex(size_type i) const;

    size_type blockOffset(size_type geometry_index, size_type entity_index) const;

    template <class ContainerSuffix>
    size_type treeDegree(ContainerSuffix suffix) const noexcept;

    void updateFixedSizeOrdering();

    void updateVariableSizeOrdering();

    /**
     * @brief Returns a span of geometry type offsets.
     */
    auto geometryTypeOffsets() const
    {
      return std::span<size_type, GT_COUNT + 1>(_offsets.get(), GT_COUNT + 1);
    }

    /**
     * @brief Returns a span of entity offsets.
     */
    auto entityOffsets() const
    {
      return std::span<size_type>(geometryTypeOffsets().data() + geometryTypeOffsets().size(), geometryTypeOffsets().back() + 1);
    }

    // Note that a copy on this class make a shallow copy of these pointers.
    // However, this class is meant to be used with value semantics.
    // Thus, the contents of this array should be re-allocated after every initialization stage.
    std::shared_ptr<size_type[]> _offsets; ///< Storage for offsets
    ProtoBasis _proto_basis;               ///< The underlying proto-basis
  };

  template <class ProtoBasis>
  PreBasisNode<ProtoBasis>::PreBasisNode(LocalBasisViewTree<ProtoBasis> &&super, const ProtoBasis &proto_basis)
      : LocalBasisViewTree<ProtoBasis>{std::move(super)}, _proto_basis{proto_basis}
  {
    // prepare index cache for binding
    // we need to be able to store all sub-entities indices in the cache
    std::size_t max_entities = 0;
    for (auto gt : _proto_basis.gridView().indexSet().types(0))
    {
      auto ref_el = referenceElement<double, Element::dimension>(gt);
      std::size_t entities = 0;
      for (int codim = 0; codim <= Element::dimension; ++codim)
        entities += ref_el.size(codim);
      max_entities = std::max(max_entities, entities);
    }
    index_cache_.resize(max_entities);
  }

  // Binds the node to a grid entity, updating index caches and offsets
  template <class ProtoBasis>
  void PreBasisNode<ProtoBasis>::bind(const Element &entity) noexcept
  {
    // find the sub entities with topological associativity
    if (_proto_basis.disjointCodimClosure())
    {
      // DG/FV branch
      const size_type gt_index = GlobalGeometryTypeIndex::index(entity.type());
      const size_type entity_index = _proto_basis.gridView().indexSet().index(entity);
      index_cache_[0] = {gt_index, entity_index};
    }
    else
    {
      std::fill(std::begin(codim_offsets_), std::end(codim_offsets_), 0);
      std::size_t codim_offset = 0;
      Hybrid::forEach(range(index_constant<Element::dimension + 1>{}), [&](auto codim) {
        std::size_t sub_entities = 0;
        if (ProtoBasis::mayContainCodim(codim) and _proto_basis.containsCodim(codim)) {
          sub_entities = entity.subEntities(codim);
          for (std::size_t s = 0; s != sub_entities; ++s) {
            const auto& sub_entity = entity.template subEntity<codim>(s);
            const size_type gt_index = GlobalGeometryTypeIndex::index(sub_entity.type());
            const size_type entity_index = _proto_basis.gridView().indexSet().index(sub_entity);
            index_cache_[codim_offset + s] = {gt_index, entity_index};
          }
        }
        codim_offset = codim_offsets_[codim+1] = (codim_offset + sub_entities); });
    }

    // update nodes offsets, sizes, element views, and finite elements
    using Accessor = Impl::LocalBasisViewAccessor;
    size_type offset = this->offset();
    TypeTree::forEachNode(*this,
      /* pre-inner node visit */
      [&]<class Node>(Node &node, const auto &path) {
        Accessor::setOffset(node, offset);
      },
      /* leaf node visit */
      [&]<class Node>(Node &node, const auto &path) {
        Accessor::setElementView(node, &entity);
        Accessor::setOffset(node, offset);
        if (TypeTree::child(_proto_basis, path).subDomain().contains(entity)) {
          Accessor::setFiniteElement(node, TypeTree::child(_proto_basis, path).finiteElementMap().find(node.element()));
          Accessor::setSize(node, node.finiteElement().size());
        } else {
          Accessor::setFiniteElement(node, nullptr);
          Accessor::setSize(node, 0);
        }
        offset += node.size();
      },
      /* post-inner node visit */
      [&]<class Node>(Node &node, const auto &path) {
        Accessor::setSize(node, offset - Accessor::offset(node));
      });
  }

  // PreBasis constructor
  template <class ProtoBasis_>
  PreBasis<ProtoBasis_>::PreBasis(const ProtoBasis &proto_basis)
      : _proto_basis{proto_basis}
  {
  }

  // Returns the grid view
  template <class ProtoBasis_>
  const typename PreBasis<ProtoBasis_>::GridView &PreBasis<ProtoBasis_>::gridView() const
  {
    return _proto_basis.gridView();
  }

  // Creates a new PreBasisNode
  template <class ProtoBasis_>
  typename PreBasis<ProtoBasis_>::Node PreBasis<ProtoBasis_>::makeNode() const
  {
    return Node{localBasisViewTree(_proto_basis), _proto_basis};
  }

  // Returns the total size (number of degrees of freedom)
  template <class ProtoBasis_>
  typename PreBasis<ProtoBasis_>::size_type PreBasis<ProtoBasis_>::size() const
  {
    return treeDegree(SizePrefix{});
  }

  // Returns the size for a given multi-index prefix
  template <class ProtoBasis_>
  typename PreBasis<ProtoBasis_>::size_type PreBasis<ProtoBasis_>::size(auto size_prefix) const
  {
    // Note: Multi-index is read Outer->Inner
    SizePrefix size_sufix;
    size_sufix.resize(size_prefix.size());
    std::reverse_copy(size_prefix.begin(), size_prefix.end(), size_sufix.begin());
    return treeDegree(size_sufix);
  }

  // Returns the dimension of the proto-basis
  template <class ProtoBasis_>
  [[nodiscard]] typename PreBasis<ProtoBasis_>::size_type PreBasis<ProtoBasis_>::dimension() const noexcept
  {
    return _proto_basis.dimension();
  }

  // Returns the maximum node size
  template <class ProtoBasis_>
  typename PreBasis<ProtoBasis_>::size_type PreBasis<ProtoBasis_>::maxNodeSize() const
  {
    return _proto_basis.maxLocalCount();
  }

  // Fills container indices for a node using an iterator
  template <class ProtoBasis_>
  template <typename It>
  It PreBasis<ProtoBasis_>::indices(const Node &node, It it) const
  {
    TypeTree::forEachLeafNode(node, [&]<class LeafNode>(const LeafNode &leaf_node, const auto &path) {
        if (leaf_node.size() == 0)
          return; // skip empty nodes
        const auto& fe = leaf_node.finiteElement();
        assert(fe.type() == leaf_node.element().type());

        if (TypeTree::child(_proto_basis, path).disjointCodimClosure()) {
          // DG/FV branch
          const auto [gt_index, entity_index] = node.index_cache_.front();

          // In this case, we know that indices are contiguous and can be
          //   reconstructed from the local index. Thus, there is no need to
          //   fill the whole vector.
          auto cir = containerIndexRange(path, gt_index, entity_index);
          assert(fe.size() == cir.size());
          [[maybe_unused]] const auto& coeffs = fe.localCoefficients();
          assert(coeffs.localKey(0).index() == 0 && "Indices in DG/FV elements shall be ordered");
          for (std::size_t dof = 0; dof != fe.size(); ++dof) {
            assert(coeffs.localKey(dof).index() == dof && "Indices in DG/FV elements shall be ordered");
            auto ci = cir[dof];
            it->resize(ci.size());
            for (std::size_t i = 0; i != ci.size(); ++i)
              (*it)[i] = ci[i];
            ++it;
          }
        } else {
          const auto& coeffs = fe.localCoefficients();
          // for each local dof, we need to figure out its assigned container index
          for (std::size_t dof = 0; dof < coeffs.size(); ++dof) {
            const auto& key = coeffs.localKey(dof);
            assert(_proto_basis.containsCodim(key.codim()));
            const auto [gt_index, entity_index] = [&]{
              auto offset = node.codim_offsets_[key.codim()];
              offset += key.subEntity();
              assert(offset < node.index_cache_.size());
              assert(offset < node.codim_offsets_[key.codim()+1]);
              return node.index_cache_[offset];
            }();
            auto cir = containerIndexRange(path, gt_index, entity_index);
            assert(key.index() < cir.size());
            auto ci = cir[key.index()];
            it->resize(ci.size());
            for (std::size_t i = 0; i != ci.size(); ++i)
              (*it)[i] = ci[i];
            ++it;
          }
        } });
    return it;
  }

  // Updates the grid view for all leaf nodes
  template <class ProtoBasis_>
  void PreBasis<ProtoBasis_>::update(const GridView& grid_view)
  {
    TypeTree::forEachLeafNode(_proto_basis, [&grid_view](auto &leaf, auto path)
                    { leaf.update(grid_view); });
    _offsets = nullptr;
  }

  // Initializes index offsets for the pre-basis
  template <class ProtoBasis_>
  void PreBasis<ProtoBasis_>::initializeIndices()
  {
    _proto_basis.initializeIndices();

    _offsets = std::unique_ptr<size_type[]>(new size_type[GT_COUNT + 1]);
    std::uninitialized_fill_n(_offsets.get(), GT_COUNT + 1, size_type{0});

    if (_proto_basis.fixedSizePerGeometryType())
      updateFixedSizeOrdering();
    else
      updateVariableSizeOrdering();
  }

  // Returns the underlying proto-basis
  template <class ProtoBasis_>
  const typename PreBasis<ProtoBasis_>::ProtoBasis &PreBasis<ProtoBasis_>::protoBasis() const
  {
    return _proto_basis;
  }

  // Returns whether local trees are merged
  template <class ProtoBasis_>
  constexpr auto PreBasis<ProtoBasis_>::merged()
  {
    return std::bool_constant<not ProtoBasis::isIndexBlocked>{};
  }

  // Returns the container index range for a given path, geometry type, and entity index
  template <class ProtoBasis_>
  auto PreBasis<ProtoBasis_>::containerIndexRange(auto path, const size_type &gt_index, const size_type &e_index) const noexcept
  {
    // Note: Multi-index is read Outer->Inner
    // get local container suffix
    assert(_proto_basis.containsGeometryType(gt_index));
    const auto lcir = _proto_basis.localTreeIndices(gt_index, e_index, path);
    auto offset = blockOffset(gt_index, e_index);

    if constexpr (merged())
      return transformedRangeView(std::move(lcir), [offset](auto mi)
                                  { return accumulate_front(mi, offset); });
    else
      return transformedRangeView(std::move(lcir), [offset](auto mi)
                                  { return push_front(mi, offset); });
  }

  // Computes geometry and entity indices from a global index
  template <class ProtoBasis_>
  auto PreBasis<ProtoBasis_>::getGeometryEntityIndex(size_type i) const
  {
    // we just need to make the inverse computation of the mapIndex funtion to find the entity index
    size_type gt_index, entity_index;
    auto dof_begin = std::begin(_proto_basis.fixedSizePerGeometryType() ? geometryTypeOffsets() : entityOffsets());
    auto dof_end = std::end(_proto_basis.fixedSizePerGeometryType() ? geometryTypeOffsets() : entityOffsets());
    auto dof_it = std::prev(std::upper_bound(dof_begin, dof_end, i));
    auto dof_dist = std::distance(dof_begin, dof_it);

    if (_proto_basis.fixedSizePerGeometryType())
    {
      gt_index = dof_dist;
      entity_index = i - *dof_it;
      entity_index /= _proto_basis.localTreeDegree(gt_index);
    }
    else
    {
      auto gt_begin = std::begin(geometryTypeOffsets());
      auto gt_end = std::end(geometryTypeOffsets());
      auto gt_it = std::prev(std::upper_bound(gt_begin, gt_end, dof_dist));
      gt_index = std::distance(gt_begin, gt_it);
      entity_index = dof_dist - *gt_it;
    }
    return std::make_pair(gt_index, entity_index);
  }

  // Computes the block offset for a given geometry and entity index
  template <class ProtoBasis_>
  typename PreBasis<ProtoBasis_>::size_type PreBasis<ProtoBasis_>::blockOffset(size_type geometry_index, size_type entity_index) const
  {
    if (_proto_basis.fixedSizePerGeometryType())
    {
      const auto &gt_offset = geometryTypeOffsets()[geometry_index];
      if constexpr (merged())
        return gt_offset + entity_index * _proto_basis.localTreeDegree(geometry_index);
      else
        return gt_offset + entity_index;
    }
    else
    {
      const auto &gt_offset = geometryTypeOffsets()[geometry_index];
      return entityOffsets()[gt_offset + entity_index];
    }
  }

  // Returns the degree (size) of the tree for a given container suffix
  template <class ProtoBasis_>
  template <class ContainerSuffix>
  typename PreBasis<ProtoBasis_>::size_type PreBasis<ProtoBasis_>::treeDegree(ContainerSuffix suffix) const noexcept
  {
    // Note: Multi-index is read Inner->Outer
    if (suffix.size() == 0)
    {
      if (_proto_basis.fixedSizePerGeometryType())
        return geometryTypeOffsets().back();
      else
        return entityOffsets().back();
    }
    else
    {
      // the next index to find out its size
      auto back_index = suffix.back();

      // we first have to figure out the entity index
      const auto [gt_index, entity_index] = getGeometryEntityIndex(back_index);

      // then, the entity ordering knows the size for a given entity.
      if constexpr (merged())
      {
        // remove offset introduced by the entity_index
        suffix.back() -= blockOffset(gt_index, entity_index);
        return _proto_basis.localTreeDegree(gt_index, entity_index, suffix);
      }
      else
      {
        return _proto_basis.localTreeDegree(gt_index, entity_index, pop_back(suffix));
      }
    }
  }

  // Updates offsets for fixed-size ordering
  template <class ProtoBasis_>
  void PreBasis<ProtoBasis_>::updateFixedSizeOrdering()
  {
    assert(_proto_basis.fixedSizePerGeometryType());

    for (std::size_t codim = 0; codim <= GridView::dimension; ++codim)
    {
      for (const auto &gt : gridView().indexSet().types(codim))
      {
        const size_type gt_index = GlobalGeometryTypeIndex::index(gt);
        if (not _proto_basis.containsGeometryType(gt_index))
          continue;
        size_type gt_size = _proto_basis.localTreeDegree(gt_index);
        const size_type gt_entity_count = gridView().indexSet().size(gt);
        if (not merged())
          gt_size = (gt_size > 0);
        geometryTypeOffsets()[gt_index + 1] = gt_size * gt_entity_count;
      }
    }

    std::partial_sum(std::begin(geometryTypeOffsets()), std::end(geometryTypeOffsets()), std::begin(geometryTypeOffsets()));
  }

  // Updates offsets for variable-size ordering
  template <class ProtoBasis_>
  void PreBasis<ProtoBasis_>::updateVariableSizeOrdering()
  {
    // in this case, we need to count the number of _used_ sub blocks
    assert(not _proto_basis.fixedSizePerGeometryType());

    for (std::size_t codim = 0; codim <= GridView::dimension; ++codim)
    {
      for (const auto &gt : gridView().indexSet().types(codim))
      {
        if (!_proto_basis.containsGeometryType(gt))
          continue;
        const size_type gt_index = GlobalGeometryTypeIndex::index(gt);
        geometryTypeOffsets()[gt_index + 1] = gridView().indexSet().size(gt);
      }
    }

    std::partial_sum(std::begin(geometryTypeOffsets()), std::end(geometryTypeOffsets()), std::begin(geometryTypeOffsets()));
    auto ptr = std::unique_ptr<size_type[]>(new size_type[geometryTypeOffsets().size() + geometryTypeOffsets().back() + 1]);
    std::copy_n(std::begin(geometryTypeOffsets()), geometryTypeOffsets().size(), ptr.get());
    _offsets = std::move(ptr);
    std::uninitialized_fill_n(entityOffsets().begin(), entityOffsets().size(), size_type{0});

    size_type carry_block = 0;
    size_type index = 0;
    for (size_type gt_index = 0; gt_index < GT_COUNT; ++gt_index)
    {
      if (not _proto_basis.containsGeometryType(gt_index))
        continue;
      const size_type entity_count = geometryTypeOffsets()[gt_index + 1] - geometryTypeOffsets()[gt_index];
      for (size_type entity_index = 0; entity_index < entity_count; ++entity_index)
      {
        const size_type dimension = _proto_basis.dimension(gt_index, entity_index);
        if (dimension != 0) {
          const size_type block_size = _proto_basis.localTreeDegree(gt_index, entity_index);
          carry_block += (merged() ? block_size : (block_size > 0));
        }
        entityOffsets()[++index] = carry_block;
      }
    }
    assert(entityOffsets().back() == carry_block);
  }

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_BASIS_PREBASIS_HH
