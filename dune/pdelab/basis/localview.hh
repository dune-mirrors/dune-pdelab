#ifndef DUNE_PDELAB_BASIS_LOCALVIEW_HH
#define DUNE_PDELAB_BASIS_LOCALVIEW_HH

#include <dune/pdelab/basis/protobasis.hh>

#include <dune/functions/functionspacebases/nodes.hh>

#include <dune/common/typetree/nodeconcepts.hh>

namespace Dune::PDELab
{

  namespace Impl
  {
    /**
     * @brief Internal accessor for protected members of BasisNodeMixin.
     *
     * This struct provides static methods to access and modify the protected members of
     * BasisNodeMixin. It is intended as a workaround for integration with Dune::Functions,
     * which requires LocalView nodes to be based on their own classes, preventing direct
     * node re-implementation. The accessor enables setting and getting internal state such
     * as offset, size, tree index, element view, and finite element, facilitating
     * compatibility with different visitor patterns used in pdelab and Dune::Functions.
     *
     * @note This is an internal implementation detail and should be used with care.
     */
    struct LocalBasisViewAccessor
    {
      /// Set the offset of a node
      static void setOffset(auto &node, auto offset)
      {
        node.setOffset(offset);
      }

      /// Set the size of a node
      static void setSize(auto &node, auto size)
      {
        node.setSize(size);
      }

      /// Set the tree index of a node
      static void setTreeIndex(auto &node, auto treeIndex)
      {
        node.setTreeIndex(treeIndex);
      }

      /// Set the entity view pointer of a node
      static void setElementView(auto &node, auto element_view)
      {
        node.setElementView(element_view);
      }

      /// Set the finite element pointer or value of a node
      static void setFiniteElement(auto &node, auto &&fe)
      {
        node.setFiniteElement(std::forward<decltype(fe)>(fe));
      }

      /// Get the offset of a node
      static auto offset(auto &node)
      {
        return node.offset();
      }
    };
  }

  /**
   * @brief Leaf node representing a local finite element basis on a single entity
   * @ingroup FunctionSpaceBases
   *
   * This class provides an ordering to the space induced by the discrete-function-space,
   * restricted to a particular entity.
   *
   * @tparam FE Type of the finite element
   * @tparam E Type of the grid entity
   *
   * @note ContainerIndex usually differs between different nodes of the same ordering.
   * In particular, it is a hybrid multi-index where known parent indices for this node
   * may be pre-filled at compile-time.
   */
  template <class FE, class E>
  class LocalBasisViewLeaf : public Functions::LeafBasisNode
  {
    friend struct Impl::LocalBasisViewAccessor;

  public:
    using FiniteElement = FE;
    using Element = E;
    using size_type = std::size_t;

    /**
     * @brief Default constructor. Creates an unbound node.
     */
    LocalBasisViewLeaf()
        : _entity_view{nullptr}, _fe_view{nullptr}
    {
    }

    LocalBasisViewLeaf(const LocalBasisViewLeaf &) = delete;
    LocalBasisViewLeaf(LocalBasisViewLeaf &&) = default;
    LocalBasisViewLeaf &operator=(const LocalBasisViewLeaf &) = delete;
    LocalBasisViewLeaf &operator=(LocalBasisViewLeaf &&) = default;

    /**
     * @brief Unbinds the node from its finite element and entity.
     */
    void unbind() noexcept
    {
      _fe_view = nullptr;
      _entity_view = nullptr;
    }

    /**
     * @brief Returns a reference to the bound finite element.
     * @warning If the size of the finite element is 0, this operation is undefined.
     */
    [[nodiscard]] const FiniteElement &finiteElement() const noexcept
    {
      assert(_fe_view);
      return *_fe_view;
    }

    /**
     * @brief Returns a reference to the bound entity.
     * @warning If the size of the entity is 0, this operation is undefined.
     */
    [[nodiscard]] const Element &element() const noexcept
    {
      assert(boundElement() && "Entity is not bound: local function is not bound or it has no support on the bound entity");
      return *_entity_view;
    }

    /**
     * @brief Checks if the node is bound to an entity.
     * @return True if bound, false otherwise.
     */
    [[nodiscard]] bool boundElement() const noexcept
    {
      return _entity_view != nullptr;
    }

  private:
    /**
     * @brief Binds a finite element to the node.
     * Stores rvalues, references lvalues.
     * @tparam FE_ Type of finite element
     * @param finite_element The finite element to bind
     */
    void setFiniteElement(FiniteElement &&finite_element)
    {
      static_assert(std::move_constructible<FiniteElement>);
      static_assert(std::is_move_assignable_v<FiniteElement>);

      if (_fe_store) [[likely]]
        (*_fe_store) = std::move(finite_element);
      else
        _fe_store = std::make_unique<FiniteElement>(std::move(finite_element));
      _fe_view = _fe_store.get();
    }

    void setFiniteElement(const FiniteElement &finite_element) noexcept
    {
      _fe_view = &finite_element;
    }

    void setFiniteElement(std::nullptr_t) noexcept
    {
      _fe_view = nullptr;
    }

    /**
     * @brief Binds a view on the entity.
     * @param entity_view Pointer to the entity
     */
    void setElementView(const Element *entity_view) noexcept
    {
      _entity_view = entity_view;
    }

    std::unique_ptr<FiniteElement> _fe_store;    ///< Storage for owned finite element
    Element const *_entity_view;                 ///< Pointer to bound entity
    FiniteElement const *_fe_view;               ///< Pointer to bound finite element
  };

  /**
   * @brief Array node for compile-time power local basis structures.
   * @tparam LocalBasisViewChild Child node type
   * @tparam degree Number of children (degree)
   */
  template <typename LocalBasisViewChild, std::size_t degree>
  struct LocalBasisViewArray : public Functions::PowerBasisNode<LocalBasisViewChild, degree>
  {
    friend struct Impl::LocalBasisViewAccessor;
    using Functions::PowerBasisNode<LocalBasisViewChild, degree>::PowerBasisNode;
  };

  /**
   * @brief Vector node for run-time power local basis structures.
   * @tparam LocalBasisViewChild Child node type
   */
  template <typename LocalBasisViewChild>
  struct LocalBasisViewVector : public Functions::DynamicPowerBasisNode<LocalBasisViewChild>
  {
    friend struct Impl::LocalBasisViewAccessor;
    using Functions::DynamicPowerBasisNode<LocalBasisViewChild>::DynamicPowerBasisNode;
  };

  /**
   * @brief Tuple node for composite local basis structures.
   * @tparam LocalBasisViewChildren Types of child nodes
   */
  template <typename... LocalBasisViewChildren>
  class LocalBasisViewTuple : public Functions::CompositeBasisNode<LocalBasisViewChildren...>
  {
    friend struct Impl::LocalBasisViewAccessor;
    using Functions::CompositeBasisNode<LocalBasisViewChildren...>::CompositeBasisNode;
  };

  /**
   * @brief Constructs a local basis view leaf node for a given leaf proto-basis.
   *
   * @tparam ProtoBasis Type of the proto basis
   * @param proto_basis The proto basis to construct the tree from
   * @return The constructed local basis view tree
   */
  template<TypeTree::Concept::LeafTreeNode ProtoBasis>
  auto localBasisViewTree(const ProtoBasis &proto_basis)
  {
    using FiniteElement = typename ProtoBasis::FiniteElementMap::Traits::FiniteElement;
    using Element = typename ProtoBasis::GridView::template Codim<0>::Entity;
    return LocalBasisViewLeaf<FiniteElement, Element>();
  }

  /**
   * @brief Constructs a local basis view array, vector, or tuple node for a given proto-basis inner node.
   *
   * @tparam ProtoBasis Type of the proto basis
   * @param proto_basis The proto basis to construct the tree from
   * @return The constructed local basis view tree
   */
  template<TypeTree::Concept::InnerTreeNode ProtoBasis>
  auto localBasisViewTree(const ProtoBasis &proto_basis)
  {
    if constexpr (TypeTree::Concept::UniformInnerTreeNode<ProtoBasis>) {

      auto degree = proto_basis.degree();
      using ChildNode = decltype(localBasisViewTree(proto_basis.child(0)));
      if constexpr (TypeTree::Concept::StaticDegreeInnerTreeNode<ProtoBasis>) {
        LocalBasisViewArray<ChildNode, ProtoBasis::degree()> node{};
        for (std::size_t i = 0; i < degree; ++i)
          node.setChild(i, localBasisViewTree(proto_basis.child(i)));
        return node;
      } else {
        LocalBasisViewVector<ChildNode> node(degree);
        for (std::size_t i = 0; i < degree; ++i)
          node.setChild(i, localBasisViewTree(proto_basis.child(i)));
        return node;
      }
    } else {
      auto range = std::make_index_sequence<ProtoBasis::degree()>{};
      auto makeChildNode = [&]<class Index>(Index i) {
        return localBasisViewTree(proto_basis.child(i));
      };
      // instantiate LocalBasisViewTuple with child local view types.
      auto node = unpackIntegerSequence([&](auto... i) {
        return LocalBasisViewTuple<decltype(makeChildNode(i))...>();
      }, range);
      // fill children of the tuple node.
      Hybrid::forEach(range, [&](auto i) {
        node.setChild(makeChildNode(i), i);
      });
      return node;
    }
  }

  /**
   * @brief Alias for the type of a local basis view tree for a given ProtoBasis.
   * @tparam ProtoBasis Type of the proto basis
   */
  template <TypeTree::Concept::TreeNode ProtoBasis>
  using LocalBasisViewTree = decltype(localBasisViewTree(std::declval<ProtoBasis>()));

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_BASIS_LOCALVIEW_HH
