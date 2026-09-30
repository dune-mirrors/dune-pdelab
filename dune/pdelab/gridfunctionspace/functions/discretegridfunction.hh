// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONS_DISCRETEGRIDFUNCTION_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONS_DISCRETEGRIDFUNCTION_HH

#include <memory>
#include <optional>
#include <vector>

#include <dune/common/shared_ptr.hh>
#include <dune/common/typetree/nodeconcepts.hh>

#include <dune/pdelab/backend/common/spacetraits.hh>
#include <dune/pdelab/common/function.hh>

/** \file
 * \brief DiscreteGridFunction for coefficient vectors of scalar dune-functions bases.
 *
 * Specialization of Dune::PDELab::DiscreteGridFunction, so that code using
 * the PDELab function interface (e.g. DifferenceSquaredAdapter,
 * integrateGridFunction, VTKGridFunctionAdapter) works with coefficient
 * vectors of dune-functions bases.
 *
 * \warning This is experimental; only bases with a single leaf are supported.
 */

namespace Dune::PDELab {

  template<typename T, typename X>
  class DiscreteGridFunction;

  namespace Impl::Functions {

    template<typename Basis>
    using DiscreteGridFunctionTraits = GridFunctionTraits<
      typename Basis::GridView,
      typename Basis::LocalView::Tree::FiniteElement::Traits::LocalBasisType::Traits::RangeFieldType,
      Basis::LocalView::Tree::FiniteElement::Traits::LocalBasisType::Traits::dimRange,
      typename Basis::LocalView::Tree::FiniteElement::Traits::LocalBasisType::Traits::RangeType>;

  } // namespace Impl::Functions

  /** \brief Discrete grid function given by a dune-functions basis and a coefficient vector
   *
   * \tparam B dune-functions global basis with a single leaf (e.g. a Lagrange basis)
   * \tparam X coefficient vector, e.g. Backend::Vector<B,...>
   */
  template<typename B, typename X>
    requires Impl::DuneFunctionsBasis<B>
  class DiscreteGridFunction<B,X>
    : public GridFunctionBase<Impl::Functions::DiscreteGridFunctionTraits<B>, DiscreteGridFunction<B,X>>
  {
    static_assert(Dune::TypeTree::Concept::LeafTreeNode<typename B::LocalView::Tree>,
      "DiscreteGridFunction is only implemented for dune-functions bases with a single leaf");

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

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONS_DISCRETEGRIDFUNCTION_HH
