// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_CONSTRAINTS_COMMON_FUNCTIONS_CONSTRAINTS_HH
#define DUNE_PDELAB_CONSTRAINTS_COMMON_FUNCTIONS_CONSTRAINTS_HH

#include <cstddef>
#include <iostream>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/typetree/traversal.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/functions/functionspacebases/subentitydofs.hh>

#include <dune/pdelab/backend/common/spacetraits.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintstransformation.hh>

/** \file
 * \brief Dirichlet constraints from PDELab constraints parameters for dune-functions bases.
 *
 * Overloads of Dune::PDELab::constraints() for dune-functions bases, with the
 * same semantics as for grid function spaces using
 * ConformingDirichletConstraints: the parameter tree (a leaf, or a
 * CompositeConstraintsParameters / PowerConstraintsParameters tree) is
 * matched with the tree of the basis, a leaf parameter is used for all
 * leaves of the basis below it, and all DOFs of the boundary intersections
 * for which p.isDirichlet(intersection, x) holds at the intersection center x
 * are constrained.
 *
 * \warning This is experimental. Other types of constraints (hanging nodes,
 * processor boundaries, grid functions as parameters) are not supported.
 */

namespace Dune::PDELab {

  namespace Impl::Functions {

    // Call f(p, leaf) for each leaf of the basis tree together with the matching parameter leaf
    template<class P, class Node, class F>
    void forEachLeafParameter(const P& p, const Node& node, F&& f)
    {
      if constexpr (P::isLeaf)
        Dune::TypeTree::forEachLeafNode(node, [&](const auto& leaf, auto&&) {
          f(p, leaf);
        });
      else
      {
        constexpr std::size_t degree = P::degree();
        if (node.degree() != degree)
          DUNE_THROW(Dune::RangeError, "The constraints parameter tree does not match the tree of the basis");
        Dune::Hybrid::forEach(Dune::range(std::integral_constant<std::size_t, degree>{}), [&](auto i) {
          forEachLeafParameter(p.child(i), node.child(i), f);
        });
      }
    }

    template<class CG>
    void printConstraints(const CG& cg)
    {
      std::cout << "constraints:" << std::endl;
      std::cout << cg.size() << " constrained degrees of freedom" << std::endl;
      for (const auto& col : cg)
      {
        std::cout << col.first << ": "; // col index
        for (const auto& row : col.second)
          std::cout << "(" << row.first << "," << row.second << ") "; // row index , value
        std::cout << std::endl;
      }
    }

  } // namespace Impl::Functions

  /** \brief Construct the Dirichlet constraints of a dune-functions basis from a constraints parameter tree
   *
   * \param p       constraints parameters, e.g. from makeBoundaryConditionFromCallable(), or a tree of them
   * \param basis   dune-functions global basis
   * \param cg      constraints container, e.g. ConstraintsContainer<Basis>
   * \param verbose print the constraints at the end
   */
  template<typename P, typename Basis, typename CG>
    requires Impl::DuneFunctionsBasis<Basis>
  void constraints(const P& p, const Basis& basis, CG& cg, const bool verbose = false)
  {
    if constexpr (not std::is_same_v<CG, EmptyTransformation>)
    {
      cg.clear();
      const auto& gridView = basis.gridView();
      auto localView = basis.localView();
      auto seDOFs = Dune::Functions::subEntityDOFs(basis);
      using Intersection = typename Basis::GridView::Intersection;
      for (const auto& element : elements(gridView))
      {
        if (not element.hasBoundaryIntersections())
          continue;
        localView.bind(element);
        unsigned int intersectionIndex = 0;
        for (const auto& intersection : intersections(gridView, element))
        {
          if (intersection.boundary())
          {
            IntersectionGeometry<Intersection> ig(intersection, intersectionIndex);
            const auto center = referenceElement(intersection.geometry()).position(0,0);
            seDOFs.bind(localView, intersection);
            Impl::Functions::forEachLeafParameter(p, localView.tree(), [&](const auto& pLeaf, const auto& leaf) {
              if (pLeaf.isDirichlet(ig, center))
                for (std::size_t i = 0; i < leaf.size(); ++i)
                  if (seDOFs.contains(leaf.localIndex(i)))
                    cg[localView.index(leaf.localIndex(i))];
            });
          }
          ++intersectionIndex;
        }
      }
      if (verbose)
        Impl::Functions::printConstraints(cg);
    }
  }

  /** \brief Construct the constraints of a dune-functions basis without constraints parameters
   *
   * A dune-functions basis has no constraints of its own (as e.g. hanging
   * nodes), hence nothing is added; as for grid function spaces, cg is not cleared.
   */
  template<typename Basis, typename CG>
    requires Impl::DuneFunctionsBasis<Basis>
  void constraints(const Basis& basis, CG& cg, const bool verbose = false)
  {
    if constexpr (not std::is_same_v<CG, EmptyTransformation>)
      if (verbose)
        Impl::Functions::printConstraints(cg);
  }

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_CONSTRAINTS_COMMON_FUNCTIONS_CONSTRAINTS_HH
