// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_COMMON_SPACETRAITS_HH
#define DUNE_PDELAB_BACKEND_COMMON_SPACETRAITS_HH

#include <algorithm>
#include <cstddef>
#include <utility>

#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/backends/istlvectorfactory.hh>
#include <dune/functions/functionspacebases/containerdescriptors.hh>

#include <dune/pdelab/ordering/utility.hh>

/** \file
 * \brief Properties of the function spaces the linear algebra backends are built for.
 *
 * The PDELab vector and matrix wrappers can be built either for a PDELab
 * grid function space or for a dune-functions global basis. The access to
 * the nested native containers is delegated to the dune-functions /
 * dune-assembler backends, which use multi-indices ordered from the
 * outermost to the innermost container level. PDELab container indices are
 * ordered the other way round, so they are presented through a reversed view.
 */

namespace Dune::PDELab::Impl {

  //! dune-functions global bases, as opposed to PDELab grid function spaces
  template<class B>
  concept DuneFunctionsBasis = requires(const B& b)
  {
    typename B::PreBasis;
    typename B::LocalView;
    typename B::MultiIndex;
    typename B::SizePrefix;
    b.preBasis();
    b.gridView();
    b.localView();
    b.dimension();
  };

  /** \brief View on a PDELab container index (inner to outer) as an outer to inner multi-index
   *
   * Only stores a pointer to the container index. PDELab omits the index of
   * blocks of size one; these innermost indices are provided as zero, and the
   * view is padded to at least minSize entries. The dune-functions vector
   * access pads by itself, the dune-assembler matrix access needs indices of
   * full length (the block level of the matrix).
   */
  template<class CI>
  class ReversedMultiIndex
  {
  public:
    explicit ReversedMultiIndex(const CI& ci, std::size_t minSize = 0)
      : ci_(&ci), size_(std::max<std::size_t>(ci.size(), minSize))
    {}

    std::size_t size() const
    {
      return size_;
    }

    std::size_t operator[](std::size_t i) const
    {
      return i < ci_->size() ? (*ci_)[ci_->size()-1-i] : 0;
    }

  private:
    const CI* ci_;
    std::size_t size_;
  };

  //! Properties of a PDELab grid function space
  template<class GFS>
  struct SpaceTraits
  {
    //! index of an entry of the vector, inner to outer
    using ContainerIndex = typename GFS::Ordering::Traits::ContainerIndex;

    //! the index as seen by the dune-functions / dune-assembler backends (outer to inner), padded to minSize
    static ReversedMultiIndex<ContainerIndex> upstreamIndex(const ContainerIndex& ci, std::size_t minSize = 0)
    {
      return ReversedMultiIndex<ContainerIndex>(ci, minSize);
    }

    //! resize a native vector for the space
    template<class Vector>
    static void resize(const GFS& gfs, Vector& v)
    {
      auto sizeProvider = SizeProviderAdapter{gfs.orderingStorage()};
      static_assert(decltype(sizeProvider)::ContainerIndexOrder == MultiIndexOrder::Outer2Inner);
      Dune::Functions::istlVectorBackend(v).resize(sizeProvider);
    }
  };

  //! Properties of a dune-functions global basis
  template<class B>
    requires DuneFunctionsBasis<B>
  struct SpaceTraits<B>
  {
    //! index of an entry of the vector, the global multi-index of the basis (outer to inner)
    using ContainerIndex = typename B::MultiIndex;

    static const ContainerIndex& upstreamIndex(const ContainerIndex& ci, std::size_t = 0)
    {
      return ci;
    }

    template<class Vector>
    static void resize(const B& basis, Vector& v)
    {
      Dune::Functions::istlVectorBackend(v).resize(basis);
    }

    //! native ISTL vector with the block structure of the basis
    template<class E>
    using NativeVector = decltype(Dune::Functions::makeISTLVector<E>(
      Dune::Functions::containerDescriptor(std::declval<const B&>().preBasis())));
  };

} // namespace Dune::PDELab::Impl

#endif // DUNE_PDELAB_BACKEND_COMMON_SPACETRAITS_HH
