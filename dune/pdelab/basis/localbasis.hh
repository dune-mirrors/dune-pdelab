#ifndef DUNE_PDELAB_BASIS_NODE_HH
#define DUNE_PDELAB_BASIS_NODE_HH

#include <dune/pdelab/basis/protobasis.hh>

#include <dune/functions/functionspacebases/nodes.hh>

namespace Dune::PDELab {

  /**
   * @brief A realization of a local space node
   * This class provides an ordering to the space induced by the
   * discrete-function-space but restricted to a particular entity. Countrary to
   * LeafEntityOrdering, this class is only instantiated after the whole
   * ordering tree has been constructed.
   *
   * @note This class fullfils the Concept::LeafLocalBasisViewLeaf, thus, is the
   * facade class for the user.
   * @note ContainerIndex usually differs between different nodes of the same
   * ordering. In particular, it is an hybrid multi-index where known parent
   * indices for this node may be pre-filled at compile-time.
   *
   * @tparam ContainerIndex  The type to store container indices
   */
  template<class FE, class E>
  class LeafLocalBasisView
    : public Functions::LeafBasisNode
  {
    // static constexpr std::size_t fem_dim = FE::Traits::LocalBasisType::Traits::dimDomain;
    // static constexpr std::size_t fem_codim = E::dimension - fem_dim;
  public:
    using FiniteElement = FE;
    // using EntitySet = ES;
    using Element = E;
    using size_type = std::size_t;
    // using MultiIndex = ContainerIndex;
    // using Path = ViewPath;


    //! Local space traits (compatiblility layer for PDELab)
    struct Traits
    {
      using FiniteElementType /*[[deprecated]]*/ = FE;
      using FiniteElement /*[[deprecated]]*/ = FE;
      using SizeType /*[[deprecated]]*/ = size_type;
    };

    //! Class constructing the local space
    LeafLocalBasisView()
      : _entity_view{ nullptr }
      , _fe_view{ nullptr }
    {
    }

    LeafLocalBasisView(const LeafLocalBasisView&) = delete;
    LeafLocalBasisView(LeafLocalBasisView&&) = default;

    LeafLocalBasisView& operator=(const LeafLocalBasisView&) = delete;
    LeafLocalBasisView& operator=(LeafLocalBasisView&&) = default;


    //! Binds a finite element: rvalues are stored, lvalues are referenced
    template<class FE_>
    void bindFiniteElement(FE_&& finite_element) noexcept
    {
      static_assert(std::same_as<std::decay_t<FE_>, FiniteElement>);
      if constexpr (std::is_rvalue_reference_v<FE_&&>) {
        static_assert(std::move_constructible<FiniteElement>);
        static_assert(std::is_move_assignable_v<FiniteElement>);

        if (_fe_store) [[likely]]
          (*_fe_store) = std::move(finite_element);
        else
          _fe_store = std::make_unique<FiniteElement>(std::move(finite_element));
        _fe_view = _fe_store.get();
      } else {
        _fe_view = &finite_element;
      }

      Functions::LeafBasisNode::setSize(finiteElement().size());
    }

    //! Binds a view on the entity. Internally, we keep a reference the object
    void bindElement(const Element* entity) noexcept {
      _entity_view = entity;
    }

    void unbind() noexcept {
      _fe_view = nullptr;
      _entity_view = nullptr;
    }

    //! Returns a view on the local finite element
    [[nodiscard]] const FiniteElement& finiteElement() const noexcept
    {
      assert(_fe_view);
      return *_fe_view;
    }

    //! Returns a view on the bound entity
    [[nodiscard]] const Element& element() const noexcept
    {
      assert(boundElement() && "Entity is not bound: local function is not bound or it has no support on the bound entity");
      return *_entity_view;
    }

    [[nodiscard]] bool boundElement() const noexcept {
      return _entity_view != nullptr;
    }

    // // returns a local view path: join(orderingViewPath(), subEntityPath())
    // [[nodiscard]] Path path() const noexcept {
    //   return _view_path;
    // }

    // [[nodiscard]] auto orderingViewPath() const noexcept {
    //   if constexpr (fem_dim != E::dimension)
    //     return pop_back(_view_path);
    //   else
    //     return _view_path;
    // }

    // [[nodiscard]] auto subEntityPath() const noexcept {
    //   if constexpr (fem_dim != E::dimension)
    //     return back(_view_path);
    //   else
    //     return TypeTree::treePath();
    // }


  private:
    std::unique_ptr<FiniteElement> _fe_store;
    Element const* _entity_view;
    FiniteElement const* _fe_view;
    // Path _view_path;
  };




//! Returns a suitable Entity Ordering tree for a discrete-function-space tree
template<class Element, Concept::Impl::ProtoBasisTree ProtoBasis>
auto makeLocalBasisView(const ProtoBasis& proto_basis)
{
  if constexpr (Concept::LeafTreeNode<ProtoBasis>) {
    using FiniteElement = typename ProtoBasis::FiniteElementMap::Traits::FiniteElement;
    return LeafLocalBasisView<FiniteElement, Element>();
  } else {

    auto makeChildNode = [&]<class Index>(Index i){
      return makeLocalBasisView<Element>(proto_basis.child(i));
    };

    if constexpr (Concept::ArrayTreeNode<ProtoBasis>) {
      constexpr std::size_t degree = ProtoBasis::degree();
      using ChildNode = decltype(makeChildNode(0));
      using Node = Functions::PowerBasisNode<ChildNode, degree>;
      typename Node::NodeStorage storage;
      for (std::size_t i = 0; i < degree; ++i)
        storage[i] = std::make_unique<ChildNode>(makeChildNode(i));
      return Node(std::move(storage));
    } else if constexpr (Concept::VectorTreeNode<ProtoBasis>) {
      std::size_t degree = proto_basis.degree();
      using ChildNode = decltype(makeChildNode(0));
      using Node = Functions::DynamicPowerBasisNode<ChildNode>;
      typename Node::NodeStorage storage(degree);
      for (std::size_t i = 0; i < degree; ++i)
        storage[i] = std::make_unique<Node>(makeChildNode(i));
      return Node(std::move(storage));
    } else {
      static_assert(Concept::TupleTreeNode<ProtoBasis>);

      auto unfold_children = [&](auto... i) {
        using Node = Functions::CompositeBasisNode<decltype(makeChildNode(i))...>;
        typename Node::NodeStorage storage{ std::make_unique<decltype(makeChildNode(i))>(makeChildNode(i))... };
        return Node(std::move(storage));
      };
      auto range = std::make_index_sequence<ProtoBasis::degree()>{};
      return unpackIntegerSequence(unfold_children, range);
    }
  }
}

template<class Element, Concept::Impl::ProtoBasisTree ProtoBasis>
using LocalBasisViewTree = decltype(makeLocalBasisView<Element>(std::declval<ProtoBasis>()));


} // namespace Dune::PDELab

#endif // DUNE_PDELAB_BASIS_NODE_HH
