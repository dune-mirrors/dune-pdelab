// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_FINITEELEMENT_EMPTY_HH
#define DUNE_PDELAB_FINITEELEMENT_EMPTY_HH

#include <dune/geometry/type.hh>

namespace Dune {
  namespace PDELab {

    /**
       @defgroup EmptyLocalFiniteElement Empty Finite Element

       \brief This is the implementation of a local finite element which has
       zero degrees of freedom.

       This element can thus only represent a trivial solution, which is constant zero.

       It may be used to \e deactivate elements (e.g. for multi-domain or unfitted finite element approaches).

       @{
     */
    template<class T>
    class EmptyLocalBasis
    {
    public:

      typedef T Traits;

      EmptyLocalBasis(){}

      //! \brief number of shape functions
      unsigned int size () const
      {
        return 0;
      }

      //! \brief Evaluate all shape functions
      inline void evaluateFunction (const typename Traits::DomainType& in,
                                    std::vector<typename Traits::RangeType>& out) const
      {
        out.clear();
      }

      //! return given derivative of all components
      template<unsigned int k>
      inline void evaluate (const std::array<int,k>& directions,
                            const typename Traits::DomainType& in,
                            std::vector<typename Traits::RangeType>& out) const
      {
        out.clear();
      }

      //! \brief Evaluate Jacobian of all shape functions
      inline void
      evaluateJacobian (const typename Traits::DomainType& in,         // position
                        std::vector<typename Traits::JacobianType>& out) const      // return value
      {
        out.resize(0);
      }

      //! \brief Polynomial order of the shape functions
      unsigned int order () const
      {
        return 0;
      }

      //! \brief Evaluate partial derivative of all shape functions
      void partial(const std::array<unsigned int, Traits::dimDomain>& order,
                          const typename Traits::DomainType& in,
                          std::vector<typename Traits::RangeType>& out) const
      {
        out.clear();
      }

    };

    class EmptyLocalCoefficients
    {
      LocalKey key;

    public:
      //! \brief Standard constructor
      EmptyLocalCoefficients ()
        :  key(LocalKey(0,0,0))
      {}

      //! number of coefficients
      std::size_t size () const
      {
        return 0;
      }

      //! get i'th index
      const LocalKey& localKey (std::size_t i) const
      {
        return key;
      }

    };

    class EmptyLocalInterpolation
    {
    public:

      EmptyLocalInterpolation ()
      { }

      //! determine coefficients interpolating a given function
      template<typename F, typename C>
      void interpolate (const F& f, std::vector<C>& out) const
      {
        out.clear();
      }
    };

    template<class BT>
    class EmptyLocalFiniteElement
    {

    public:

      typedef LocalFiniteElementTraits<
          EmptyLocalBasis<BT>,
          EmptyLocalCoefficients,
          EmptyLocalInterpolation
          > Traits;

      EmptyLocalFiniteElement (const Dune::GeometryType & _gtype) :
        gtype(_gtype)
      {}

      const typename Traits::LocalBasisType& localBasis () const
      {
        return basis;
      }

      const typename Traits::LocalCoefficientsType& localCoefficients () const
      {
        return coefficients;
      }

      const typename Traits::LocalInterpolationType& localInterpolation () const
      {
        return interpolation;
      }

      GeometryType type () const
      {
        return gtype;
      }

      std::size_t size() const
      {
        return basis.size();
      }

    private:
      const Dune::GeometryType gtype;

      EmptyLocalBasis<BT> basis;
      EmptyLocalCoefficients coefficients;
      EmptyLocalInterpolation interpolation;
    };
    /**
        @}
     */

  }
}

#endif // DUNE_PDELAB_FINITEELEMENT_EMPTY_HH
