// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_VARIABLEFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_VARIABLEFEM_HH

#include <functional>
#include <memory>

#include <dune/geometry/type.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include "finiteelementmap.hh"

namespace Dune {
  namespace PDELab {

    template<class T>
    class VirtualLocalFiniteElement
    {
      using Interface = LocalFiniteElementVirtualInterface<T>;
      template<typename LFE>
      using Wrapper = LocalFiniteElementVirtualImp<LFE>;
    public:
      using LocalBasisTraits = T;
      using Traits =
        LocalFiniteElementTraits<
          LocalBasisVirtualInterface<LocalBasisTraits>,
          LocalCoefficientsVirtualInterface,
          LocalInterpolationVirtualInterface<
            typename LocalBasisTraits::DomainType,
            typename LocalBasisTraits::RangeType> >;

      VirtualLocalFiniteElement(const VirtualLocalFiniteElement & other) : _lfe(other._lfe)
      {}

      VirtualLocalFiniteElement(const Interface & iface) : _lfe(iface.clone)
      {}

      template<class LFE>
      VirtualLocalFiniteElement(LocalFiniteElementVirtualImp<LFE> && imp) : _lfe(std::move(imp))
      {}

      template<class LFE>
      VirtualLocalFiniteElement(LFE && imp) : _lfe(std::make_shared<Wrapper<LFE>>(imp))
      {}

      //! \copydoc LocalFiniteElementVirtualInterface::localBasis
      const typename Traits::LocalBasisType& localBasis () const
      {
        return _lfe->localBasis();
      }

      //! \copydoc LocalFiniteElementVirtualInterface::localCoefficients
      const typename Traits::LocalCoefficientsType& localCoefficients () const
      {
        return _lfe->localCoefficients();
      }

      //! \copydoc LocalFiniteElementVirtualInterface::localInterpolation
      const typename Traits::LocalInterpolationType& localInterpolation () const
      {
        return _lfe->localInterpolation();
      }

      //! \copydoc LocalFiniteElementVirtualInterface::size
      unsigned int size () const
      {
        return _lfe->size();
      }

      //! \copydoc LocalFiniteElementVirtualInterface::type
      const GeometryType type () const
      {
        return _lfe->type();
      }

    private:
      std::shared_ptr<Interface> _lfe;
    };

    // template<typename LFE, typename... Args>
    // auto makeWrappedLocalFiniteElement(Args... && args)
    // {
    // }

#warning pass GetFEM as only parameter to VariableLocalFiniteElementMap
#warning use Signature of GetFEM to deduce BasisTraits and Element in VariableLocalFiniteElementMapTraits
    template<class BasisTraits, class Element>
    struct VariableLocalFiniteElementMapTraits
    {
      //! Type of finite element from local functions
      using FiniteElementType = VirtualLocalFiniteElement<BasisTraits>;
      using Signature = FiniteElementType(Element);
      using DefaultGetFEM = std::function<Signature>;
    };

    //! FiniteElementMap which provides ialLocalFiniteElement instances, depending on the local polynomial degree
    //! \ingroup FiniteElementMap
    template<class BasisTraits,
             class Element,
             class GetFEM = typename VariableLocalFiniteElementMapTraits<BasisTraits,Element>::DefaultGetFEM>
    class VariableLocalFiniteElementMap
    {
      using T = VariableLocalFiniteElementMapTraits<BasisTraits,Element>;
      //! Type of finite element from local functions
      using FiniteElementType = typename T::FiniteElementType;
    public:
      using Traits = FiniteElementMapTraits<FiniteElementType>;

      //! The dimension of the finite elements returned by this map.
      static constexpr int dimension = BasisTraits::dimDomain;

      /** Construct a VariableLocalFiniteElementMap with user-defined LocalFiniteElement provided by getFEM */
      VariableLocalFiniteElementMap (GetFEM && getFEM) : _getFEM(getFEM) {}

      //! \brief get local basis functions for entity
      typename Traits::FiniteElementType find (const Element& e) const
      {
        return _getFEM(e);
      }

      static constexpr bool fixedSize()
      {
        // return true;
        return false;
      }

      static constexpr bool hasDOFs(int codim)
      {
        return true; // we never know, so let's say we might have DOFs at that codimension...
        // return codim == 0;
      }

      std::size_t size(GeometryType gt) const
      {
        // return gt.isQuadrilateral() ? 9 : 0;
        DUNE_THROW(VariableElementSize,"VariableLocalFiniteElementMap can contain elements of variable order.");
      }

      std::size_t maxLocalSize() const
      {
        // return 9;
        DUNE_THROW(VariableElementSize,"VariableLocalFiniteElementMap can contain elements of variable order.");
      }

    private:
      GetFEM _getFEM;
    };



  }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_VARIABLEMONOMFEM_HH
