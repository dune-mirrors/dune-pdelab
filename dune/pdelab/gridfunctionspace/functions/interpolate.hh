// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONS_INTERPOLATE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONS_INTERPOLATE_HH

#include <cstddef>
#include <type_traits>
#include <vector>

#include <dune/common/ftraits.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/typetree/nodeconcepts.hh>
#include <dune/common/typetree/traversal.hh>

#include <dune/typetree/typetree.hh>

#include <dune/pdelab/backend/common/spacetraits.hh>
// needs dune-typetree, which function/localfunctionhelper.hh does not include itself
#include <dune/pdelab/function/localfunction.hh>

/** \file
 * \brief Interpolation of PDELab grid functions into dune-functions bases.
 *
 * Overload of Dune::PDELab::interpolate() for dune-functions bases, with the
 * same semantics as for grid function spaces: a scalar function is
 * interpolated into every leaf, a vector-valued function leaf by leaf (one
 * component per leaf), and a tree of functions (e.g. CompositeGridFunction)
 * is matched with the tree of the basis.
 *
 * \warning This is experimental.
 */

namespace Dune::PDELab {

  namespace Impl::Functions {

    template<class Node, class LocalView, class Coefficients, class LocalFunction>
    void interpolateLeaf(const Node& node, const LocalView& localView, Coefficients& x, const LocalFunction& f)
    {
      using E = std::decay_t<decltype(x[localView.index(0)])>;
      std::vector<E> xl(node.size());
      node.finiteElement().localInterpolation().interpolate(f, xl);
      for (std::size_t i = 0; i < node.size(); ++i)
        x[localView.index(node.localIndex(i))] = xl[i];
    }

    template<class LF, class Node, class LocalView, class Coefficients>
    void interpolateTree(const LF& lf, const Node& node, const LocalView& localView, Coefficients& x)
    {
      if constexpr (LF::isLeaf)
      {
        using Domain = typename LocalView::Element::Geometry::LocalCoordinate;
        using Range = std::decay_t<decltype(lf(std::declval<Domain>()))>;
        if constexpr (Dune::TypeTree::Concept::LeafTreeNode<Node>)
          interpolateLeaf(node, localView, x, lf);
        else if constexpr (std::is_convertible_v<Range, typename Dune::FieldTraits<Range>::field_type>)
        {
          // scalar function: interpolate into every leaf
          Dune::TypeTree::forEachLeafNode(node, [&](const auto& leaf, auto&&) {
            interpolateLeaf(leaf, localView, x, lf);
          });
        }
        else
        {
          // vector-valued function: one component per leaf
          std::size_t index = 0;
          Dune::TypeTree::forEachLeafNode(node, [&](const auto& leaf, auto&&) {
            auto component = [&](const Domain& xi) { return lf(xi)[index]; };
            interpolateLeaf(leaf, localView, x, component);
            ++index;
          });
        }
      }
      else
      {
        // tree of functions: match the children
        Dune::Hybrid::forEach(Dune::range(node.degree()), [&](auto i) {
          interpolateTree(lf.child(i), node.child(i), localView, x);
        });
      }
    }

  } // namespace Impl::Functions

  /** \brief Interpolate a PDELab grid function (or tree of grid functions) into a dune-functions basis
   *
   * \param f     PDELab grid function, callable adapter or tree of functions
   * \param basis dune-functions global basis
   * \param xg    coefficient vector, e.g. Backend::Vector<Basis,...>
   */
  template<typename F, typename Basis, typename XG>
    requires Impl::DuneFunctionsBasis<Basis>
  void interpolate(const F& f, const Basis& basis, XG& xg)
  {
    auto lf = makeLocalFunctionTree(f, basis.gridView());
    auto localView = basis.localView();
    for (const auto& element : elements(basis.gridView()))
    {
      lf.bind(element);
      localView.bind(element);
      Impl::Functions::interpolateTree(lf, localView.tree(), localView, xg);
      lf.unbind();
    }
  }

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONS_INTERPOLATE_HH
