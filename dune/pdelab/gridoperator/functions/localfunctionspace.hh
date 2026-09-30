// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDOPERATOR_FUNCTIONS_LOCALFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDOPERATOR_FUNCTIONS_LOCALFUNCTIONSPACE_HH

#include <array>
#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/indices.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/common/typetree/childaccess.hh>
#include <dune/common/typetree/nodeconcepts.hh>

#include <dune/functions/functionspacebases/globalvaluedlocalfiniteelement.hh>

/** \file
 * \brief Adapter presenting a dune-functions basis tree as a PDELab local function space.
 *
 * PDELab local operators only use a small part of the local function space
 * interface: size(), localIndex(), finiteElement() on leaves, child access,
 * Traits::FiniteElementType and Traits::SizeType. dune-functions basis nodes
 * already provide size(), localIndex() (offset into the element-local numbering
 * of the LocalView, which is exactly PDELab's notion) and the tree structure.
 * The adapter adds the remaining PDELab-specific parts and mirrors the tree.
 *
 * Used by the GridOperator specialization for dune-functions bases.
 *
 * \warning This is experimental.
 */

namespace Dune::PDELab {

  /** \brief Common base of the local function space adapter nodes
   *
   * Makes Dune::PDELab an associated namespace of the adapter nodes, so that
   * unqualified calls of PDELab helpers with a local function space argument
   * in local operators (e.g. makeJacobianContainer(lfsu)) are found by
   * argument-dependent lookup, as for the PDELab local function spaces.
   */
  struct FunctionsBasisLocalFunctionSpaceTag {};

} // namespace Dune::PDELab

namespace Dune::PDELab::Impl::Functions {

  // Existing local operators call child(lfs,i) unqualified and rely on
  // argument-dependent lookup finding TypeTree::child (PDELab local function
  // spaces derive from TypeTree nodes). Make it visible for the adapter nodes.
  using Dune::TypeTree::child;


    template<class FE>
    struct IsGlobalValuedFiniteElement : std::false_type {};

    template<class Transformator, class LocalValuedLFE, class Element>
    struct IsGlobalValuedFiniteElement<
      Dune::Functions::Impl::GlobalValuedLocalFiniteElement<Transformator,LocalValuedLFE,Element>>
      : std::true_type {};

    template<class Node>
    concept StaticPowerNode =
      TypeTree::Concept::StaticDegreeInnerTreeNode<Node> && TypeTree::Concept::UniformInnerTreeNode<Node>;

    template<class Node>
    concept CompositeNode =
      TypeTree::Concept::StaticDegreeInnerTreeNode<Node> && not TypeTree::Concept::UniformInnerTreeNode<Node>;

    template<class Node>
    concept DynamicPowerNode =
      TypeTree::Concept::UniformInnerTreeNode<Node> && not TypeTree::Concept::StaticDegreeInnerTreeNode<Node>;


  /** \brief Grid-level information exposed via lfs.gridFunctionSpace()
   *
   * A few local operators (e.g. the block and point diagonal wrappers) query
   * the grid view from the local function space. This provides the minimal
   * interface for them.
   */
  template<class B>
  class SpaceInfo
  {
  public:
    using Basis = B;
    using GridView = typename B::GridView;
    using EntitySet = GridView;

    explicit SpaceInfo(const Basis& basis)
      : basis_(&basis)
    {}

    const Basis& basis() const { return *basis_; }
    const GridView& gridView() const { return basis_->gridView(); }
    const EntitySet& entitySet() const { return basis_->gridView(); }

  private:
    const Basis* basis_;
  };

  /** \brief Node of the local function space adapter tree
   *
   * \tparam Space SpaceInfo of the global basis
   * \tparam Node  dune-functions basis node type wrapped by this adapter node
   *
   * The adapter holds a pointer to the corresponding node in the
   * LocalView's tree. Binding the LocalView updates that node in place,
   * so the adapter never needs to be rebound itself.
   */
  template<class Space, class Node>
  class LocalFunctionSpaceNode;

  //! Common part of all adapter nodes
  template<class Space, class Node>
  class LocalFunctionSpaceNodeBase
    : public FunctionsBasisLocalFunctionSpaceTag
  {
  public:
    using size_type = std::size_t;
    using BasisNode = Node;

    struct Traits
    {
      using SizeType = std::size_t;
      using GridFunctionSpace = Space;
      using GridView = typename Space::GridView;
      using GridViewType = GridView;
      using EntitySet = typename Space::EntitySet;
    };

    LocalFunctionSpaceNodeBase(const Space& space, const Node& node)
      : space_(&space), node_(&node)
    {}

    //! number of DOFs in this subtree for the currently bound element
    size_type size() const { return node_->size(); }

    //! map index in this subtree to the index in the element-local vector of the root
    size_type localIndex(size_type i) const { return node_->localIndex(i); }

    const Space& gridFunctionSpace() const { return *space_; }

    //! the wrapped dune-functions basis node
    const Node& basisNode() const { return *node_; }

  protected:
    const Space* space_;
    const Node* node_;
  };

  //! Leaf node, provides the finite element
  template<class Space, class Node>
    requires TypeTree::Concept::LeafTreeNode<Node>
  class LocalFunctionSpaceNode<Space,Node>
    : public LocalFunctionSpaceNodeBase<Space,Node>
  {
    using Base = LocalFunctionSpaceNodeBase<Space,Node>;

    using FE = std::decay_t<decltype(std::declval<const Node&>().finiteElement())>;

  public:
    static_assert(not IsGlobalValuedFiniteElement<FE>::value,
      "The dune-functions basis uses global-valued finite elements (e.g. Raviart-Thomas, "
      "Nedelec). Their shape functions are already transformed to the element, while PDELab "
      "local operators expect local-valued finite elements.");

    struct Traits : Base::Traits
    {
      using FiniteElementType = FE;
      using FiniteElement = FE;
    };

    LocalFunctionSpaceNode(const Space& space, const Node& node)
      : Base(space, node)
    {}

    static constexpr auto degree() { return index_constant<0>{}; }

    const FE& finiteElement() const { return this->node_->finiteElement(); }
  };

  //! Inner node with a compile-time number of identical children (PowerBasisNode)
  template<class Space, class Node>
    requires StaticPowerNode<Node>
  class LocalFunctionSpaceNode<Space,Node>
    : public LocalFunctionSpaceNodeBase<Space,Node>
  {
    using Base = LocalFunctionSpaceNodeBase<Space,Node>;
    static constexpr std::size_t n = Node::degree();
    using ChildNode = std::decay_t<decltype(std::declval<const Node&>().child(0))>;
    using Child = LocalFunctionSpaceNode<Space,ChildNode>;

  public:
    LocalFunctionSpaceNode(const Space& space, const Node& node)
      : Base(space, node)
      , children_(makeChildren(space, node, std::make_index_sequence<n>{}))
    {}

    static constexpr auto degree() { return index_constant<n>{}; }

    template<class Index>
      requires std::is_convertible_v<Index,std::size_t>
    const Child& child(Index i) const { return children_[i]; }

    //! old-style child access, as used by some existing local operators
    template<std::size_t i>
    const Child& child() const { return children_[i]; }

  private:
    template<std::size_t... i>
    static std::array<Child,n> makeChildren(const Space& space, const Node& node, std::index_sequence<i...>)
    {
      return {{ Child(space, node.child(i))... }};
    }

    std::array<Child,n> children_;
  };

  //! Inner node with a runtime number of identical children (DynamicPowerBasisNode)
  template<class Space, class Node>
    requires DynamicPowerNode<Node>
  class LocalFunctionSpaceNode<Space,Node>
    : public LocalFunctionSpaceNodeBase<Space,Node>
  {
    using Base = LocalFunctionSpaceNodeBase<Space,Node>;
    using ChildNode = std::decay_t<decltype(std::declval<const Node&>().child(0))>;
    using Child = LocalFunctionSpaceNode<Space,ChildNode>;

  public:
    LocalFunctionSpaceNode(const Space& space, const Node& node)
      : Base(space, node)
    {
      children_.reserve(node.degree());
      for (std::size_t i = 0; i < node.degree(); ++i)
        children_.emplace_back(space, node.child(i));
    }

    std::size_t degree() const { return children_.size(); }

    template<class Index>
      requires std::is_convertible_v<Index,std::size_t>
    const Child& child(Index i) const { return children_[i]; }

  private:
    std::vector<Child> children_;
  };

  //! Inner node with a compile-time number of heterogeneous children (CompositeBasisNode)
  template<class Space, class Node>
    requires CompositeNode<Node>
  class LocalFunctionSpaceNode<Space,Node>
    : public LocalFunctionSpaceNodeBase<Space,Node>
  {
    using Base = LocalFunctionSpaceNodeBase<Space,Node>;
    static constexpr std::size_t n = Node::degree();

    template<std::size_t i>
    using ChildNode = std::decay_t<decltype(std::declval<const Node&>().child(index_constant<i>{}))>;

    template<class Seq> struct ChildrenTuple;
    template<std::size_t... i>
    struct ChildrenTuple<std::index_sequence<i...>>
    {
      using type = std::tuple<LocalFunctionSpaceNode<Space,ChildNode<i>>...>;
    };
    using Children = typename ChildrenTuple<std::make_index_sequence<n>>::type;

  public:
    LocalFunctionSpaceNode(const Space& space, const Node& node)
      : Base(space, node)
      , children_(makeChildren(space, node, std::make_index_sequence<n>{}))
    {}

    static constexpr auto degree() { return index_constant<n>{}; }

    template<std::size_t i>
    const auto& child(index_constant<i>) const { return std::get<i>(children_); }

    //! old-style child access, as used by some existing local operators
    template<std::size_t i>
    const auto& child() const { return std::get<i>(children_); }

  private:
    template<std::size_t... i>
    static Children makeChildren(const Space& space, const Node& node, std::index_sequence<i...>)
    {
      return Children(LocalFunctionSpaceNode<Space,ChildNode<i>>(space, node.child(index_constant<i>{}))...);
    }

    Children children_;
  };


  /** \brief Local function space on top of an existing dune-functions LocalView
   *
   * Does not own the LocalView: whoever owns it (e.g. the dune-assembler
   * global assembler) binds it to elements. The adapter tree (tree()) is
   * what is handed to PDELab local operators; its localIndex() refers to
   * the element-local numbering of the LocalView, so index(localIndex(i))
   * is the global multi-index.
   *
   * The adapter tree stores pointers into the LocalView's tree and into
   * this object, hence the class is neither copyable nor movable.
   */
  template<class LV>
  class LocalFunctionSpaceView
  {
  public:
    using LocalView = LV;
    using Basis = typename LocalView::GlobalBasis;
    using Space = SpaceInfo<Basis>;
    using Element = typename LocalView::Element;
    using Tree = LocalFunctionSpaceNode<Space, typename LocalView::Tree>;
    using MultiIndex = typename LocalView::MultiIndex;
    using size_type = std::size_t;

    explicit LocalFunctionSpaceView(const LocalView& localView)
      : localView_(&localView)
      , space_(localView.globalBasis())
      , tree_(space_, localView.tree())
    {}

    LocalFunctionSpaceView(const LocalFunctionSpaceView&) = delete;
    LocalFunctionSpaceView& operator=(const LocalFunctionSpaceView&) = delete;

    //! adapter tree to be passed to PDELab local operators
    const Tree& tree() const { return tree_; }

    //! number of DOFs on the bound element (size of the element-local vector)
    size_type size() const { return localView_->size(); }
    size_type maxSize() const { return localView_->maxSize(); }

    //! global multi-index of the i-th entry of the element-local vector
    const MultiIndex& index(size_type i) const { return localView_->index(i); }

    const LocalView& localView() const { return *localView_; }
    const Basis& basis() const { return space_.basis(); }
    const Element& element() const { return localView_->element(); }

  private:
    const LocalView* localView_;
    Space space_;
    Tree tree_;
  };


  /** \brief Local function space owning a LocalView of a dune-functions global basis
   *
   * Like LocalFunctionSpaceView, but owns and binds the LocalView itself.
   */
  template<class B>
  class LocalFunctionSpace
  {
  public:
    using Basis = B;
    using LocalView = typename Basis::LocalView;
    using View = LocalFunctionSpaceView<LocalView>;
    using Element = typename LocalView::Element;
    using Tree = typename View::Tree;
    using MultiIndex = typename LocalView::MultiIndex;
    using size_type = std::size_t;

    explicit LocalFunctionSpace(const Basis& basis)
      : localView_(basis.localView())
      , view_(localView_)
    {}

    // The view refers to the own LocalView, so copies have to rebuild it.
    LocalFunctionSpace(const LocalFunctionSpace& other)
      : localView_(other.localView_)
      , view_(localView_)
    {}

    LocalFunctionSpace& operator=(const LocalFunctionSpace&) = delete;

    void bind(const Element& element) { localView_.bind(element); }
    void unbind() { localView_.unbind(); }

    //! the non-owning view on the own LocalView
    const View& view() const { return view_; }

    //! adapter tree to be passed to PDELab local operators
    const Tree& tree() const { return view_.tree(); }

    //! number of DOFs on the bound element (size of the element-local vector)
    size_type size() const { return localView_.size(); }
    size_type maxSize() const { return localView_.maxSize(); }

    //! global multi-index of the i-th entry of the element-local vector
    const MultiIndex& index(size_type i) const { return localView_.index(i); }

    const LocalView& localView() const { return localView_; }
    const Basis& basis() const { return view_.basis(); }
    const Element& element() const { return localView_.element(); }

  private:
    LocalView localView_;
    View view_;
  };

} // namespace Dune::PDELab::Impl::Functions

namespace Dune::PDELab {

  /** \brief Stand-in for the FiniteElementMap template parameter of some local operators
   *
   * Some PDELab local operators (e.g. ConvectionDiffusionFEM) are templated on the
   * FiniteElementMap only to extract the finite element type, e.g. for a
   * LocalBasisCache. This provides the required traits for a given finite element.
   */
  template<class FE>
  struct FiniteElementMapFor
  {
    struct Traits
    {
      using FiniteElementType = FE;
      using FiniteElement = FE;
    };
  };

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_GRIDOPERATOR_FUNCTIONS_LOCALFUNCTIONSPACE_HH
