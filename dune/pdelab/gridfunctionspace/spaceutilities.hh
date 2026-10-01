// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_SPACEUTILITIES_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_SPACEUTILITIES_HH

#include <cstddef>

#include <dune/pdelab/backend/common/spacetraits.hh>
#include <dune/pdelab/constraints/common/constraintstransformation.hh>

/** \file
 * \brief Uniform access to properties of PDELab grid function spaces and dune-functions bases.
 *
 * Code written against these helpers works with both kinds of function
 * spaces, e.g.
 * \code
 * using CC = Dune::PDELab::ConstraintsContainer<Space, double>;
 * MatrixBackend mb(4 * Dune::PDELab::maxLocalSize(space));
 * std::cout << Dune::PDELab::globalSize(space) << " DOFs" << std::endl;
 * \endcode
 */

namespace Dune::PDELab {

  namespace Impl {

    template<class Space, class E>
    struct ConstraintsContainerFor;

    // grid function spaces: the nested type, EmptyTransformation for spaces without constraints
    template<class Space, class E>
      requires requires { typename Space::template ConstraintsContainer<E>::Type; }
    struct ConstraintsContainerFor<Space,E>
    {
      using type = typename Space::template ConstraintsContainer<E>::Type;
    };

    // dune-functions bases: constraints indexed by the multi-indices of the basis
    template<class Space, class E>
      requires DuneFunctionsBasis<Space>
    struct ConstraintsContainerFor<Space,E>
    {
      using type = ConstraintsTransformation<typename Space::MultiIndex, typename Space::MultiIndex, E>;
    };

  } // namespace Impl

  /** \brief Type of the constraints container of a grid function space or a dune-functions basis
   *
   * For a grid function space this is `typename Space::template ConstraintsContainer<E>::Type`
   * (EmptyTransformation for spaces without constraints).
   *
   * \tparam Space PDELab grid function space or dune-functions global basis
   * \tparam E     field type of the constraints coefficients
   */
  template<class Space, class E = double>
  using ConstraintsContainer = typename Impl::ConstraintsContainerFor<Space,E>::type;

  /** \brief Maximal number of degrees of freedom on an element
   *
   * \c space.maxLocalSize() for grid function spaces, the maximal size of a local view for bases.
   */
  template<class Space>
  std::size_t maxLocalSize(const Space& space)
  {
    if constexpr (Impl::DuneFunctionsBasis<Space>)
      return space.localView().maxSize();
    else
      return space.maxLocalSize();
  }

  /** \brief Total number of degrees of freedom
   *
   * \c space.globalSize() for grid function spaces, \c space.dimension() for bases.
   *
   * \note For blocked bases \c space.size() is the number of top-level blocks, not the number of DOFs.
   */
  template<class Space>
  std::size_t globalSize(const Space& space)
  {
    if constexpr (Impl::DuneFunctionsBasis<Space>)
      return space.dimension();
    else
      return space.globalSize();
  }

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_SPACEUTILITIES_HH
