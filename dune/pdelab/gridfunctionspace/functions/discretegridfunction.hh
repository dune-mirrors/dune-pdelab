// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONS_DISCRETEGRIDFUNCTION_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONS_DISCRETEGRIDFUNCTION_HH

#include <cstddef>
#include <concepts>
#include <memory>
#include <optional>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/typetree/childaccess.hh>
#include <dune/common/typetree/nodeconcepts.hh>

#include <dune/pdelab/backend/common/spacetraits.hh>
#include <dune/pdelab/common/function.hh>
// primary templates; this header is included at its end
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>

/** \file
 * \brief DiscreteGridFunction and VectorDiscreteGridFunction for dune-functions bases.
 *
 * Specializations of Dune::PDELab::DiscreteGridFunction and
 * Dune::PDELab::VectorDiscreteGridFunction, so that code using the PDELab
 * function interface (e.g. DifferenceSquaredAdapter, integrateGridFunction,
 * VTKGridFunctionAdapter) works with coefficient vectors of dune-functions
 * bases. As GridFunctionSubSpace for grid function spaces, a subspace basis
 * (Dune::Functions::subspaceBasis()) selects a part of the coefficient vector
 * of the root basis.
 *
 * \warning This is experimental.
 */

namespace Dune::PDELab {

  namespace Impl::Functions {

    //! A dune-functions global basis or subspace basis
    template<class B>
    concept GlobalBasisLike = requires(const B& basis)
    {
      typename B::GridView;
      typename B::LocalView;
      typename B::MultiIndex;
      { basis.gridView() } -> std::convertible_to<const typename B::GridView&>;
      { basis.localView() } -> std::convertible_to<typename B::LocalView>;
      basis.dimension();
    };

    template<typename Basis>
    using DiscreteGridFunctionTraits = GridFunctionTraits<
      typename Basis::GridView,
      typename Basis::LocalView::Tree::FiniteElement::Traits::LocalBasisType::Traits::RangeFieldType,
      Basis::LocalView::Tree::FiniteElement::Traits::LocalBasisType::Traits::dimRange,
      typename Basis::LocalView::Tree::FiniteElement::Traits::LocalBasisType::Traits::RangeType>;

  } // namespace Impl::Functions

  /** \brief Discrete grid function given by a dune-functions basis and a coefficient vector
   *
   * \tparam B dune-functions global basis or subspace basis with a single leaf (e.g. a Lagrange basis)
   * \tparam X coefficient vector of the root basis, e.g. Backend::Vector<RootBasis,...>
   */
  template<typename B, typename X>
    requires Impl::Functions::GlobalBasisLike<B>
  class DiscreteGridFunction<B,X>
    : public GridFunctionBase<Impl::Functions::DiscreteGridFunctionTraits<B>, DiscreteGridFunction<B,X>>
  {
    static_assert(Dune::TypeTree::Concept::LeafTreeNode<typename B::LocalView::Tree>,
      "DiscreteGridFunction is only implemented for dune-functions bases with a single leaf, use a subspace basis");

    using BaseT = GridFunctionBase<Impl::Functions::DiscreteGridFunctionTraits<B>, DiscreteGridFunction<B,X>>;
    using Element = typename B::LocalView::Element;

  public:
    using Traits = typename BaseT::Traits;

    DiscreteGridFunction(const B& basis, const X& x)
      : basis_(stackobject_to_shared_ptr(basis))
      , x_(stackobject_to_shared_ptr(x))
      , localView_(basis.localView())
    {}

    DiscreteGridFunction(std::shared_ptr<const B> basis, std::shared_ptr<const X> x)
      : basis_(basis)
      , x_(x)
      , localView_(basis->localView())
    {}

    inline void evaluate(const typename Traits::ElementType& e,
                         const typename Traits::DomainType& x,
                         typename Traits::RangeType& y) const
    {
      if (not boundElement_ or not (*boundElement_ == e))
      {
        localView_.bind(e);
        boundElement_ = e;
      }
      const auto& node = localView_.tree();
      node.finiteElement().localBasis().evaluateFunction(x, values_);
      y = 0;
      for (std::size_t i = 0; i < node.size(); ++i)
        y.axpy((*x_)[localView_.index(node.localIndex(i))], values_[i]);
    }

    inline const typename Traits::GridViewType& getGridView() const
    {
      return basis_->gridView();
    }

  private:
    std::shared_ptr<const B> basis_;
    std::shared_ptr<const X> x_;
    mutable typename B::LocalView localView_;
    mutable std::optional<Element> boundElement_;
    mutable std::vector<typename Traits::RangeType> values_;
  };

  namespace Impl {

    //! For bases, the default number of components is the degree of the (power) tree
    template<typename B>
      requires Functions::GlobalBasisLike<B>
    struct VectorDiscreteGridFunctionDimension<B>
    {
      static constexpr std::size_t value = std::size_t(B::LocalView::Tree::degree());
    };

  } // namespace Impl

  /** \brief Vector-valued discrete grid function given by a dune-functions (subspace) basis and a coefficient vector
   *
   * The basis tree is a power node of scalar leaves, e.g. the velocity
   * subspaceBasis(basis, _0) of a Taylor-Hood basis; component k of the
   * function is given by the child remap[k].
   *
   * \tparam B    dune-functions global basis or subspace basis with a power tree
   * \tparam X    coefficient vector of the root basis, e.g. Backend::Vector<RootBasis,...>
   * \tparam dimR number of components
   */
  template<typename B, typename X, std::size_t dimR>
    requires Impl::Functions::GlobalBasisLike<B>
  class VectorDiscreteGridFunction<B,X,dimR>
    : public GridFunctionBase<
        GridFunctionTraits<
          typename B::GridView,
          typename Dune::TypeTree::Child<typename B::LocalView::Tree,0>::FiniteElement::Traits::LocalBasisType::Traits::RangeFieldType,
          dimR,
          Dune::FieldVector<typename Dune::TypeTree::Child<typename B::LocalView::Tree,0>::FiniteElement::Traits::LocalBasisType::Traits::RangeFieldType, dimR>>,
        VectorDiscreteGridFunction<B,X,dimR>>
  {
    using Tree = typename B::LocalView::Tree;
    using ChildType = Dune::TypeTree::Child<Tree,0>;
    static_assert(Dune::TypeTree::Concept::LeafTreeNode<ChildType>,
      "VectorDiscreteGridFunction is only implemented for dune-functions bases with a power tree of leaves");

    using RF = typename ChildType::FiniteElement::Traits::LocalBasisType::Traits::RangeFieldType;
    using RT = typename ChildType::FiniteElement::Traits::LocalBasisType::Traits::RangeType;
    using BaseT = GridFunctionBase<GridFunctionTraits<typename B::GridView, RF, dimR, Dune::FieldVector<RF,dimR>>,
                                   VectorDiscreteGridFunction<B,X,dimR>>;
    using Element = typename B::LocalView::Element;

  public:
    using Traits = typename BaseT::Traits;

    /** \brief construct
     *
     * \param basis dune-functions (subspace) basis
     * \param x     coefficient vector
     * \param start number of the first child of the basis tree to use
     */
    VectorDiscreteGridFunction(const B& basis, const X& x, std::size_t start = 0)
      : VectorDiscreteGridFunction(stackobject_to_shared_ptr(basis), stackobject_to_shared_ptr(x), start)
    {}

    VectorDiscreteGridFunction(std::shared_ptr<const B> basis, std::shared_ptr<const X> x, std::size_t start = 0)
      : basis_(basis)
      , x_(x)
      , localView_(basis->localView())
    {
      for (std::size_t k = 0; k < dimR; ++k)
        remap_[k] = k + start;
    }

    inline void evaluate(const typename Traits::ElementType& e,
                         const typename Traits::DomainType& x,
                         typename Traits::RangeType& y) const
    {
      if (not boundElement_ or not (*boundElement_ == e))
      {
        localView_.bind(e);
        boundElement_ = e;
      }
      const auto& tree = localView_.tree();
      for (std::size_t k = 0; k < dimR; ++k)
      {
        const auto& node = tree.child(remap_[k]);
        node.finiteElement().localBasis().evaluateFunction(x, values_);
        y[k] = 0.0;
        for (std::size_t i = 0; i < node.size(); ++i)
          y[k] += (*x_)[localView_.index(node.localIndex(i))] * values_[i][0];
      }
    }

    inline const typename Traits::GridViewType& getGridView() const
    {
      return basis_->gridView();
    }

  private:
    std::shared_ptr<const B> basis_;
    std::shared_ptr<const X> x_;
    std::size_t remap_[dimR];
    mutable typename B::LocalView localView_;
    mutable std::optional<Element> boundElement_;
    mutable std::vector<RT> values_;
  };

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONS_DISCRETEGRIDFUNCTION_HH
