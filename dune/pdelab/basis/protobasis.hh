#ifndef DUNE_PDELAB_BASIS_PROTOBASIS_HH
#define DUNE_PDELAB_BASIS_PROTOBASIS_HH

#include <array>
#include <string>
#include <utility>
#include <memory>
#include <tuple>
#include <vector>

#include <dune/pdelab/basis/topologicalinterleaving.hh>
#include <dune/pdelab/basis/utility.hh>

namespace Dune::PDELab
{

  //! \ingroup FunctionSpaceBases
  //! \{

  //! \brief Represents the entire computational domain.
  //! This class can be used as a sub-domain that includes all entities of the grid view.
  //! It provides a method to check if an entity belongs to the domain (always true).
  class EntireDomain
  {
  public:
    template <class Entity>
    constexpr static auto contains(const Entity&)
    {
      return std::true_type();
    }
  };

  /**
   * @brief Leaf proto-basis node
   *
   * @tparam GridView_              Type of the grid view
   * @tparam FiniteElementMap_      Type of the finite element map @ref FiniteElementMap
   * @tparam IndexBlocked           A boolean flag indicating whether the container is blocked or lexicographically merged into a flat index (see @ref BasisTopologicalInterleaving).
   *
   * @details This is the fundamental class to build trees of bases.
   * It provides a map to the local finite elements and the indexation of their degrees of freedom through the @ref BasisTopologicalInterleaving object.
   * This object could be composed with @ref ProtoBasisArray, @ref ProtoBasisVector, or @ref ProtoBasisTuple to form a tree of proto-bases
   * that interleave (or groups) the degrees of freedom by entity, meaning that all degrees of freedom associated with the same topological entity
   * are ordered together. Or it can used be to form a @ref PreBasis, which can be used to transform multi-indices even further and expose them as
   * a final basis object with the dune-functions interface.
   */
  template <class GridView_, class FiniteElementMap_, bool IndexBlocked = false, class SubDomain_ = EntireDomain>
  class ProtoBasisLeaf
      : public BasisTopologicalInterleaving<ProtoBasisLeaf<GridView_, FiniteElementMap_, IndexBlocked, SubDomain_>>
  {
    using Base = BasisTopologicalInterleaving<ProtoBasisLeaf<GridView_, FiniteElementMap_, IndexBlocked, SubDomain_>>;
  public:
    //! Whether the indices of this node are blocked or lexicographically merged into a flat index (see @ref BasisTopologicalInterleaving).
    static constexpr bool isIndexBlocked = IndexBlocked;

    //! Type to represent offsets and sizes
    using typename Base::size_type;
    //! Type of the grid view
    using GridView = GridView_;
    //! Type of the finite element map
    using FiniteElementMap = FiniteElementMap_;
    //! Type of the sub-domain
    using SubDomain = SubDomain_;

    /**
     * @brief Construct a new leaf proto-basis object
     *
     * @param grid_view           Grid view associated with this leaf node
     * @param finite_element_map  Pointer to a finite element map
     * @param sub_domain          Sub-domain associated with this leaf node (default is entire domain)
     */
    ProtoBasisLeaf(GridView grid_view,
                   const FiniteElementMap& finite_element_map,
                   SubDomain sub_domain = {})
        : Base()
        , _grid_view{grid_view}
        , _finite_element_map{finite_element_map}
        , _sub_domain{sub_domain}
    {
    }

    //! Get a reference to the finite element map
    const FiniteElementMap &finiteElementMap() const
    {
      return _finite_element_map;
    }

    //! @brief Gets the grid view associated with this leaf node
    const GridView &gridView() const
    {
      return _grid_view;
    }

    //! @brief Gets the sub-domain associated with this leaf node
    const SubDomain& subDomain() const
    {
      return _sub_domain;
    }

    //! @brief Updates the grid view
    void update(GridView grid_view)
    {
      _grid_view = grid_view;
    }

    //! Returns the number of children nodes
    static constexpr auto degree()
    {
      return std::integral_constant<std::size_t, 0>{};
    }

    //! @brief Returns a container descriptor for a given geometry type and entity index
    auto containerDescriptor(size_type gt_index, size_type e_index) const
    {
      // common size for all geometry types
      constexpr auto gt_common_size = commonSizePerGeometryType();
      using namespace Dune::PDELab::ContainerDescriptors;

      if constexpr (gt_common_size)
        return FlatArray<gt_common_size.value()>();
      else
        return FlatVector(this->localTreeDegree(gt_index, e_index));
    }

    /**
     * @brief Gets a common size for all active geometry types if available at compile time.
     *
     * @return An optional size if a common size is available at compile time; otherwise, std::nullopt.
     *
     * @note This requires that ProtoBasis::FiniteElementMap::size(...) is static constexpr.
     */
    static constexpr std::optional<std::size_t> commonSizePerGeometryType()
    {
      if constexpr (not std::same_as<SubDomain, EntireDomain>)
        return std::nullopt;
      std::size_t size = 0;
      const auto fem_size = [](auto gt) constexpr
      {
        if constexpr (requires { { FiniteElementMap::size(gt) } -> std::convertible_to<std::size_t>; })
          return FiniteElementMap::size(gt);
        else
          return 0;
      };
      // iterate over all possible geometry types and find out if all share the same size
      for (std::size_t dim = 0; dim <= FiniteElementMap::dimension; ++dim)
      {
        std::size_t gt_size = fem_size(GeometryTypes::none(dim));
        if (gt_size > 0)
        {
          if (size > 0 and size != gt_size)
            return std::nullopt;
          else
            size = gt_size;
        }
        for (std::size_t topology_id = 0; topology_id < (std::size_t{1} << dim); ++topology_id)
        {
          std::size_t gt_size = fem_size(GeometryType(topology_id, dim));
          if (gt_size > 0)
          {
            if (size > 0 and size != gt_size)
              return std::nullopt;
            else
              size = gt_size;
          }
        }
      }
      if (size == 0)
        return std::nullopt;
      else
        return size;
    }

  private:
    GridView _grid_view;
    FiniteElementMap _finite_element_map;
    SubDomain _sub_domain;
  };


  /**
   * @brief Array of proto-basis nodes
   *
   * @tparam ProtoBasisChild  A proto-basis node (@ref ProtoBasisLeaf, @ref ProtoBasisArray, @ref ProtoBasisVector, or @ref ProtoBasisTuple)
   * @tparam degree           Number of children proto-bases
   * @tparam IndexBlocked     A boolean flag indicating whether the container is blocked or lexicographically merged into a flat index (see @ref BasisTopologicalInterleaving).
  *
   * @details This is the fundamental class to build trees of bases.
   * It provides the indexation of their degrees of freedom through the @ref BasisTopologicalInterleaving object.
   * This object could be composed with @ref ProtoBasisArray, @ref ProtoBasisVector, or @ref ProtoBasisTuple to form a tree of proto-bases
   * that interleave (or groups) the degrees of freedom by entity, meaning that all degrees of freedom associated with the same topological entity
   * are ordered together. Or it can used be to form a @ref PreBasis, which can be used to transform multi-indices even further and expose them as
   * a final basis object with the dune-functions interface.
   */
  template <class ProtoBasisChild,
            std::size_t degree_,
            bool IndexBlocked = false>
  class ProtoBasisArray
    : public BasisTopologicalInterleaving<ProtoBasisArray<ProtoBasisChild, degree_, IndexBlocked>>
  {
    using Base = BasisTopologicalInterleaving<ProtoBasisArray<ProtoBasisChild, degree_, IndexBlocked>>;
  public:
    //! Whether the indices of this node are blocked or lexicographically merged into a flat index (see @ref BasisTopologicalInterleaving).
    static constexpr bool isIndexBlocked = IndexBlocked;

    //! Type to represent offsets and sizes
    using typename Base::size_type;
    //! Type of the grid view
    using GridView = typename ProtoBasisChild::GridView;
    //! Type of the children nodes
    using ChildType = ProtoBasisChild;

    //! @brief Constructs a array of proto-basis
    //! @pre All grid views in the child nodes should be the same
    ProtoBasisArray(
        const std::array<ProtoBasisChild, degree_> &nodes)
        : Base()
        , _nodes{nodes}
    {
    }

    //! Copies the contents of another proto-basis object (value semantics)
    ProtoBasisArray(const ProtoBasisArray &) = default;

    //! Copies the contents of another proto-basis object (value semantics)
    ProtoBasisArray &operator=(const ProtoBasisArray &) = default;

    //! Moves the contents from another proto-basis object
    ProtoBasisArray(ProtoBasisArray &&) = default;

    //! Moves the contents from another proto-basis object
    ProtoBasisArray &operator=(ProtoBasisArray &&) = default;

    //! Returns the grid view associated with this node.
    const GridView& gridView() const { return this->child(0).gridView(); }

    //! Returns the number of children nodes
    static constexpr auto degree()
    {
      return std::integral_constant<std::size_t, degree_>{};
    }

    //! Returns the i-th proto-basis child
    const ProtoBasisChild &child(std::size_t i) const
    {
      return _nodes[i];
    }

    //! Returns the i-th proto-basis child
    ProtoBasisChild &child(std::size_t i)
    {
      return _nodes[i];
    }

    //! @brief Returns a container descriptor for a given geometry type and entity index
    auto containerDescriptor(size_type gt_index, std::size_t e_index) const
    {
      using namespace Dune::PDELab::ContainerDescriptors;
      auto descriptor0 = child(0).containerDescriptor(gt_index, e_index);
      using ChildDescriptor = std::remove_cvref_t<decltype(descriptor0)>;
      if constexpr (std::same_as<ChildDescriptor, Unknown>)
        return Unknown{};
      else if constexpr (Base::mergedLocalTrees()) {
        using Block = BlockType<ChildDescriptor>;
        if constexpr (requires{ { ChildDescriptor::size() } -> std::convertible_to<std::size_t>; } ) {
          if constexpr (IsCompileTimeUniform<ChildDescriptor>)
            return UniformArray<Block, ChildDescriptor::size() * degree()>(descriptor0[0]);
          else {
            Array<Block, ChildDescriptor::size() * degree()> array;
            for (std::size_t i = 0; i < degree(); ++i) {
              auto desc = child(i).containerDescriptor(gt_index, e_index);
              for (std::size_t j = 0; j < desc.size(); ++j)
                array[i * desc.size() + j] = std::move(desc[j]);
            }
            return array;
          }
        } else {
          if constexpr (IsCompileTimeUniform<ChildDescriptor>)
            return UniformVector<Block>(descriptor0.size() * degree(), descriptor0[0]);
          else {
            Vector<Block> vector;
            for (std::size_t i = 0; i < degree(); ++i) {
              auto desc = child(i).containerDescriptor(gt_index, e_index);
              for (std::size_t j = 0; j < desc.size(); ++j)
                vector.emplace_back(std::move(desc[j]));
            }
            return vector;
          }
        }
      } else {
        if constexpr (IsCompileTimeUniform<ChildDescriptor> )
          return UniformArray<ChildDescriptor, degree()>(descriptor0);
        else {
          Array<ChildDescriptor, degree()> array;
          for (std::size_t i = 0; i < degree(); ++i)
            array[i] = child(i).containerDescriptor(gt_index, e_index);
          return array;
        }
      }
    }

  private:
    std::array<ProtoBasisChild, degree_> _nodes;
  };

  /**
   * @brief Vector of proto-basis nodes
   *
   * @tparam ProtoBasisChild  A proto-basis node (@ref ProtoBasisLeaf, @ref ProtoBasisArray, @ref ProtoBasisVector, or @ref ProtoBasisTuple)
   * @tparam IndexBlocked     A boolean flag indicating whether the container is blocked or lexicographically merged into a flat index (see @ref BlockedTopologicalInterleaving and @ref FlatTopologicalInterleaving).
  *
   * @details This is the fundamental class to build trees of bases.
   * It provides the indexation of their degrees of freedom through the @ref BasisTopologicalInterleaving object.
   * This object could be composed with @ref ProtoBasisArray, @ref ProtoBasisVector, or @ref ProtoBasisTuple to form a tree of proto-bases
   * that interleave (or groups) the degrees of freedom by entity, meaning that all degrees of freedom associated with the same topological entity
   * are ordered together. Or it can used be to form a @ref PreBasis, which can be used to transform multi-indices even further and expose them as
   * a final basis object with the dune-functions interface.
   */
  template <class ProtoBasisChild, bool IndexBlocked = false>
  class ProtoBasisVector
    : public BasisTopologicalInterleaving<ProtoBasisVector<ProtoBasisChild, IndexBlocked>>
  {
    using Base = BasisTopologicalInterleaving<ProtoBasisVector<ProtoBasisChild, IndexBlocked>>;
  public:
    //! Whether the indices of this node are blocked or lexicographically merged into a flat index (see @ref BasisTopologicalInterleaving).
    static constexpr bool isIndexBlocked = IndexBlocked;

    //! Type to represent offsets and sizes
    using typename Base::size_type;
    //! Type of the grid view
    using GridView = typename ProtoBasisChild::GridView;
    //! Type of the children nodes
    using ChildType = ProtoBasisChild;

    //! @brief Constructs a vector of proto-basis
    //! @pre All grid views in the child nodes should be the same
    ProtoBasisVector(const std::vector<ProtoBasisChild> &nodes)
        : Base()
        , _nodes{nodes}
    {
    }

    //! Copies the contents of another proto-basis object (value semantics)
    ProtoBasisVector(const ProtoBasisVector &) = default;

    //! Copies the contents of another proto-basis object (value semantics)
    ProtoBasisVector &operator=(const ProtoBasisVector &) = default;

    //! Moves the contents from another proto-basis object
    ProtoBasisVector(ProtoBasisVector &&) = default;

    //! Moves the contents from another proto-basis object
    ProtoBasisVector &operator=(ProtoBasisVector &&) = default;

    //! Returns the grid view associated with this node
    const GridView& gridView() const { return this->child(0).gridView(); }

    //! Returns the number of children nodes
    std::size_t degree() const
    {
      return _nodes.size();
    }

    //! Returns the i-th proto-basis child
    const ProtoBasisChild &child(std::size_t i) const
    {
      return _nodes[i];
    }

    //! Returns the i-th proto-basis child
    ProtoBasisChild &child(std::size_t i)
    {
      return _nodes[i];
    }

    //! @brief Returns a container descriptor for a given geometry type and entity index
    auto containerDescriptor(size_type gt_index, std::size_t e_index) const
    {
      using namespace Dune::PDELab::ContainerDescriptors;
      auto descriptor0 = child(0).containerDescriptor(gt_index, e_index);
      using ChildDescriptor = std::remove_cvref_t<decltype(descriptor0)>;
      if constexpr (std::same_as<ChildDescriptor, Unknown>)
        return Unknown{};
      else if constexpr (Base::mergedLocalTrees()) {
        using Block = BlockType<ChildDescriptor>;
        if constexpr (IsCompileTimeUniform<ChildDescriptor>)
          return UniformVector<Block>(descriptor0.size() * degree(), descriptor0[0]);
        else {
          Vector<Block> vector;
          for (std::size_t i = 0; i < degree(); ++i) {
            auto desc = child(i).containerDescriptor(gt_index, e_index);
            for (std::size_t j = 0; j < desc.size(); ++j)
              vector.emplace_back(std::move(desc[j]));
          }
          return vector;
        }
      } else {
        if constexpr (IsCompileTimeUniform<ChildDescriptor>)
          return UniformVector<ChildDescriptor>(degree(), descriptor0);
        else {
          Vector<ChildDescriptor> vector;
          for (std::size_t i = 0; i < degree(); ++i)
            vector.emplace_back(child(i).containerDescriptor(gt_index, e_index));
          return vector;
        }
      }
    }

  private:
    std::vector<ProtoBasisChild> _nodes;
  };

  /**
   * @brief Tuple of proto-basis nodes
   *
   * @tparam IndexBlocked        A boolean flag indicating whether the container is blocked or lexicographically merged into a flat index (see @ref FlatTopologicalInterleaving and @ref BlockedTopologicalInterleaving).
   * @tparam ProtoBasisChildren  A variadic set of proto-basis nodes (@ref ProtoBasisLeaf, @ref ProtoBasisArray, @ref ProtoBasisVector, or @ref ProtoBasisTuple)
  *
   * @details This is the fundamental class to build trees of bases.
   * It provides the indexation of their degrees of freedom through the @ref BasisTopologicalInterleaving object.
   * This object could be composed with @ref ProtoBasisArray, @ref ProtoBasisVector, or @ref ProtoBasisTuple to form a tree of proto-bases
   * that interleave (or groups) the degrees of freedom by entity, meaning that all degrees of freedom associated with the same topological entity
   * are ordered together. Or it can used be to form a @ref PreBasis, which can be used to transform multi-indices even further and expose them as
   * a final basis object with the dune-functions interface.
   */
  template <bool IndexBlocked, class... ProtoBasisChildren>
  class ProtoBasisTuple
    : public BasisTopologicalInterleaving<ProtoBasisTuple<IndexBlocked, ProtoBasisChildren...>>
  {
    using Base = BasisTopologicalInterleaving<ProtoBasisTuple<IndexBlocked, ProtoBasisChildren...>>;
  public:

    //! Whether the indices of this node are blocked or lexicographically merged into a flat index (see @ref BasisTopologicalInterleaving).
    static constexpr bool isIndexBlocked = IndexBlocked;

    //! Type to represent offsets and sizes
    using typename Base::size_type;
    //! Type of the grid view
    using GridView = std::common_type_t<typename ProtoBasisChildren::GridView...>;

    //! @brief Constructs a tuple of proto-basis
    //! @pre All grid views in the child nodes should be the same
    ProtoBasisTuple(const std::tuple<ProtoBasisChildren...> &nodes)
        : Base()
        , _nodes{nodes}
    {
    }

    //! Copies the contents of another proto-basis object (value semantics)
    ProtoBasisTuple(const ProtoBasisTuple &) = default;

    //! Copies the contents of another proto-basis object (value semantics)
    ProtoBasisTuple &operator=(const ProtoBasisTuple &) = default;

    //! Moves the contents from another proto-basis object
    ProtoBasisTuple(ProtoBasisTuple &&) = default;

    //! Moves the contents from another proto-basis object
    ProtoBasisTuple &operator=(ProtoBasisTuple &&) = default;

    //! Returns the grid view associated with this node
    const GridView& gridView() const { return this->child(index_constant<0>()).gridView(); }

    //! Returns the number of children nodes
    static constexpr auto degree()
    {
      return std::integral_constant<std::size_t, sizeof...(ProtoBasisChildren)>{};
    }

    //! Returns the i-th proto-basis child
    template <std::size_t i>
    const auto &child(index_constant<i> = {}) const
    {
      return std::get<i>(_nodes);
    }

    //! Returns the i-th proto-basis child
    template <std::size_t i>
    auto &child(index_constant<i> = {})
    {
      return std::get<i>(_nodes);
    }

    //! @brief Returns a container descriptor for a given geometry type and entity index
    auto containerDescriptor(size_type gt_index, std::size_t e_index) const
    {
      using namespace Dune::PDELab::ContainerDescriptors;
      auto childDescriptors = std::apply([=](const auto&... child_i){
        return std::tuple<std::remove_cvref_t<decltype(child_i.containerDescriptor(gt_index, e_index))>...>{
          child_i.containerDescriptor(gt_index, e_index)...
        };
      }, _nodes);

      if constexpr (std::apply([&]<class... Desc>(const Desc&... descriptors) {
                               return (std::same_as<Desc, Unknown> || ...);
                             },
                             childDescriptors))
        return Unknown{};
      else if constexpr (Base::mergedLocalTrees()) {
        return std::apply(
          [&]<class... Descriptors>(const Descriptors&... descriptors) {
            if constexpr (requires { std::common_type_t<BlockType<Descriptors>...>(); }) {
              using Block = std::common_type_t<BlockType<Descriptors>...>;
              if constexpr (requires{ { (Descriptors::size() + ...) } -> std::convertible_to<std::size_t>; } ) {
                if constexpr (IsCompileTimeUniform<Block>)
                  return UniformArray<Block, (Descriptors::size() + ...)>(std::get<0>(childDescriptors)[Indices::_0]);
                else {
                  Array<Block, (Descriptors::size() + ...)> array;
                  Hybrid::forEach([&](auto i){
                    auto& desc = std::get<i>(childDescriptors);
                    for (std::size_t j = 0; j < desc.size(); ++j)
                      array[i * desc.size() + j] = std::move(desc[j]);
                  }, range(degree()));
                  return array;
                }
              } else {
                if constexpr (IsCompileTimeUniform<Block>)
                  return UniformVector<Block>((descriptors.size() + ...), std::get<0>(childDescriptors)[Indices::_0]);
                else {
                  Vector<Block> vector;
                  Hybrid::forEach([&](auto i){
                    auto& desc = std::get<i>(childDescriptors);
                    for (std::size_t j = 0; j < desc.size(); ++j)
                      vector.emplace_back(std::move(desc[j]));
                  }, range(degree()));
                  return vector;
                }
              }
            } else {
              return Unknown();
            }
          },
          childDescriptors);
      } else {
        return TupleVector{childDescriptors};
      }
    }

  private:
    std::tuple<ProtoBasisChildren...> _nodes;
  };

  //! \}
} // namespace Dune::PDELab

#endif // DUNE_PDELAB_BASIS_PROTOBASIS_HH
