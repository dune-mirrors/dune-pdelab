// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONS_LOCALFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONS_LOCALFUNCTIONSPACE_HH

#include <cstddef>

#include <dune/pdelab/backend/common/spacetraits.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspacetags.hh>
#include <dune/pdelab/gridoperator/functions/localfunctionspace.hh>

/** \file
 * \brief LocalFunctionSpace and LFSIndexCache for dune-functions bases.
 *
 * Specializations of Dune::PDELab::LocalFunctionSpace and
 * Dune::PDELab::LFSIndexCache, so that code written for grid function spaces
 * (bind, child access, localIndex(), containerIndex()) also works with
 * dune-functions bases, e.g. to access single DOFs:
 * \code
 * LocalFunctionSpace<Basis> lfs(basis);
 * LFSIndexCache<decltype(lfs)> cache(lfs);
 * lfs.bind(element);
 * cache.update();
 * constraints[cache.containerIndex(lfs.child(_1).localIndex(0))];
 * \endcode
 * The container indices are the multi-indices of the basis, i.e. the keys of
 * ConstraintsContainer<Basis> and the indices of Backend::Vector<Basis>.
 */

namespace Dune::PDELab {

  template<typename GFS, typename TAG>
  class LocalFunctionSpace;

  namespace Impl::Functions {

    /** \brief Local function space of a dune-functions basis with the interface of a PDELab LFS
     *
     * Owns the local view and behaves like the root node of a PDELab local
     * function space tree: child access, degree(), size(), localIndex(i) and,
     * for leaf bases, finiteElement().
     */
    template<class B>
    class BasisLocalFunctionSpace
      : public LocalFunctionSpace<B>
    {
      using Base = LocalFunctionSpace<B>;

    public:
      using Base::Base;

      static constexpr auto degree() { return Base::Tree::degree(); }

      template<class... Indices>
      decltype(auto) child(Indices... indices) const
        requires requires (const typename Base::Tree& t) { t.child(indices...); }
      {
        return this->tree().child(indices...);
      }

      std::size_t localIndex(std::size_t i) const { return this->tree().localIndex(i); }

      decltype(auto) finiteElement() const
        requires requires (const typename Base::Tree& t) { t.finiteElement(); }
      {
        return this->tree().finiteElement();
      }

      //! the basis, in the role of the grid function space
      const B& gridFunctionSpace() const { return this->basis(); }
    };

  } // namespace Impl::Functions

  //! LocalFunctionSpace of a dune-functions basis
  template<typename B, typename TAG>
    requires Impl::DuneFunctionsBasis<B>
  class LocalFunctionSpace<B,TAG>
    : public Impl::Functions::BasisLocalFunctionSpace<B>
  {
    using Impl::Functions::BasisLocalFunctionSpace<B>::BasisLocalFunctionSpace;
  };

  //! LocalFunctionSpace of a dune-functions basis (resolves the ambiguity with the AnySpaceTag specialization)
  template<typename B>
    requires Impl::DuneFunctionsBasis<B>
  class LocalFunctionSpace<B,AnySpaceTag>
    : public Impl::Functions::BasisLocalFunctionSpace<B>
  {
    using Impl::Functions::BasisLocalFunctionSpace<B>::BasisLocalFunctionSpace;
  };

  template<typename LFS, typename C, bool fast>
  class LFSIndexCache;

  /** \brief LFSIndexCache of a dune-functions basis
   *
   * The container index of a local DOF is the multi-index of the basis. The
   * constraints container is not used; constraints are accessed directly with
   * the container indices.
   */
  template<typename B, typename TAG, typename C, bool fast>
    requires Impl::DuneFunctionsBasis<B>
  class LFSIndexCache<LocalFunctionSpace<B,TAG>,C,fast>
  {
  public:
    using LocalFunctionSpace = Dune::PDELab::LocalFunctionSpace<B,TAG>;
    using ContainerIndex = typename B::MultiIndex;
    using size_type = std::size_t;

    explicit LFSIndexCache(const LocalFunctionSpace& lfs)
      : lfs_(&lfs)
    {}

    LFSIndexCache(const LocalFunctionSpace& lfs, const C&, bool = true)
      : lfs_(&lfs)
    {}

    //! nothing to do, the indices are computed when the local function space is bound
    void update() {}

    const ContainerIndex& containerIndex(size_type i) const { return lfs_->index(i); }

    size_type size() const { return lfs_->size(); }

    const LocalFunctionSpace& localFunctionSpace() const { return *lfs_; }

  private:
    const LocalFunctionSpace* lfs_;
  };

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_FUNCTIONS_LOCALFUNCTIONSPACE_HH
