// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_FUNCTION_DISCRETEGRIDVIEWFUNCTION_HH
#define DUNE_PDELAB_FUNCTION_DISCRETEGRIDVIEWFUNCTION_HH

#include <cstdlib>
#include <vector>
#include <memory>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/common/jacobiantocurl.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>

namespace std {
  /** \brief specialization of common_type for FieldVector and scalar */
  template<typename T, int N, typename R2>
  struct common_type<Dune::FieldVector<T,N>, R2>
  {
    using type = Dune::FieldVector<typename std::common_type<T,R2>::type,N>;
  };
  /** \brief specialization of common_type for two FieldVectors */
  template<typename T, int N, typename R2>
  struct common_type<Dune::FieldVector<T,N>, Dune::FieldVector<R2,N>>
  {
    using type = Dune::FieldVector<typename std::common_type<T,R2>::type,N>;
  };
}

namespace Dune {
namespace PDELab {

template<typename Signature, typename E, template<class> class D, int B,
         int diffOrder>
struct DiscreteGridViewFunctionTraits :
  Functions::Imp::GridFunctionTraits<
    typename DiscreteGridViewFunctionTraits<Signature,E,D,B,diffOrder-1>::DerivativeSignature
    ,E,D,B>
{};

template<typename Signature, typename E, template<class> class D, int B>
struct DiscreteGridViewFunctionTraits<Signature,E,D,B,0> :
  Functions::Imp::GridFunctionTraits<Signature,E,D,B>
{};

template<typename GFS, typename V, int diffOrder = 0>
class DiscreteGridViewFunction
{
public:
  using GridView = typename GFS::Traits::GridView;
  using EntitySet = Functions::GridViewEntitySet<GridView, 0>;

  using Domain = typename EntitySet::GlobalCoordinate;
  using LocalBasisTraits = typename GFS::Traits::FiniteElementMap::Traits::FiniteElement::Traits::LocalBasisType::Traits;
  using LocalBasisRange = typename LocalBasisTraits::RangeType;
  using VectorRange = typename V::ElementType;
  using ElementaryRange = typename std::common_type<LocalBasisRange, VectorRange>::type;

  using LocalDomain = typename EntitySet::LocalCoordinate;
  using Element = typename EntitySet::Element;

  using Traits =
    DiscreteGridViewFunctionTraits<ElementaryRange(Domain), EntitySet,
                                   Functions::DefaultDerivativeTraits, 16, diffOrder>;

  using Range = typename Traits::Range; // this is actually either the Range, the Jacobian or Hessian

  using Basis = GFS;
  using GridFunctionSpace = GFS;
  using Vector = V;

  class LocalFunction
  {
    using LFS = LocalFunctionSpace<GridFunctionSpace>;
    using LFSCache = LFSIndexCache<LFS>;
    using XView = typename Vector::template ConstLocalView<LFSCache>;

  public:

    using GlobalFunction = DiscreteGridViewFunction;
    using Domain = LocalDomain;
    using Range = GlobalFunction::Range;
    using Element = GlobalFunction::Element;
    using size_type = std::size_t;

    LocalFunction(const shared_ptr<const GridFunctionSpace> gfs, const shared_ptr<const Vector> v)
      : pgfs_(gfs)
      , v_(v)
      , lfs_(*pgfs_)
      , lfs_cache_(lfs_)
      , x_view_(*v_)
      , xl_(pgfs_->maxLocalSize())
      , yb_(pgfs_->maxLocalSize())
      , element_(nullptr)
    {}

    /**
     * \brief Bind LocalFunction to grid element.
     *
     * You must call this method before evaluate()
     * and after changes to the coefficient vector.
     */
    void bind(const Element& element)
    {
      element_ = &element;
      lfs_.bind(element);
      lfs_cache_.update();
      x_view_.bind(lfs_cache_);
      x_view_.read(xl_);
      x_view_.unbind();
    }

    void unbind()
    {
      element_ = nullptr;
    }

    const Element& localContext() const
    {
#ifndef NDEBUG
      if (!element_)
        DUNE_THROW(InvalidStateException,"can't call localContext on unbound DiscreteGridViewFunction::LocalFunction");
#endif
      return *element_;
    }

    friend typename DiscreteGridViewFunction<GFS,V,diffOrder+1>::LocalFunction derivative(const LocalFunction& t)
    {
      typename DiscreteGridViewFunction<GFS,V,diffOrder+1>::LocalFunction
        diff(t.pgfs_, t.v_);
      // TODO: do we really want this?
      if (t.element_) diff.bind(*t.element_);
      return diff;
    }

    /**
     * \brief Evaluate LocalFunction at bound element.
     *
     * The result of this method is undefined if you did
     * not call bind() beforehand or changed the coefficient
     * vector after the last call to bind(). In the latter case
     * you have to call bind() again in order to make operator()
     * usable.
     */
    Range
    operator()(const Domain& coord)
    {
      Range r;
      evaluate<LocalBasisTraits::diffOrder>(coord, r);
      return r;
    };

  private:
    using ElementaryJacobian =
      DiscreteGridViewFunctionTraits<ElementaryRange(Domain), EntitySet,
                                     Functions::DefaultDerivativeTraits, 16, 1>;
    using ElementaryHessian =
      DiscreteGridViewFunctionTraits<ElementaryRange(Domain), EntitySet,
                                     Functions::DefaultDerivativeTraits, 16, 2>;

    template<int maxDiffOrder, typename... T>
    void evaluate(const T&...) const
    {
      if (diffOrder > 2) DUNE_THROW(NotImplemented,
        "Derivatives are only implemented up to degree 2");
      if (diffOrder > maxDiffOrder) DUNE_THROW(NotImplemented,
        "Derivative of degree " << diffOrder << "is not provided by the local basis");
      DUNE_THROW(Exception, "unexpected error");
    };

    template<int maxDiffOrder>
    void evaluate(const Domain& coord,
      ElementaryRange& r) const
    {
      auto& basis = lfs_.finiteElement().localBasis();
      basis.evaluateFunction(coord,yb_);
      for (size_type i = 0; i < yb_.size(); ++i)
      {
        r.axpy(xl_[i],yb_[i]);
      }
    }

    template<int maxDiffOrder>
    typename std::enable_if<maxDiffOrder >= 1>::type
    evaluate(const Domain& coord, ElementaryJacobian & r) const
    {
      // get Jacobian of geometry
      const typename Element::Geometry::JacobianInverseTransposed
        JgeoIT = element_->geometry().jacobianInverseTransposed(coord);

      // get local Jacobians/gradients of the shape functions
      lfs_.finiteElement().localBasis().evaluateJacobian(coord,yb_);

      Range gradphi;
      r = 0;
      for(std::size_t i = 0; i < yb_.size(); ++i) {
        assert(gradphi.size() == yb_[i].size());
        for(std::size_t j = 0; j < gradphi.size(); ++j) {
          // compute global gradient of shape function i
          // graphi += {J^{-1}}^T * yb_i0
          JgeoIT.mv(yb_[i][j], gradphi[j]);

          // sum up global gradients, weighting them with the appropriate coeff
          // r \in R^{1,dim}
          // r_0 += xl_i * grad \phi
          r[j].axpy(xl_[i], gradphi[j]);
        }
      }
      return r;
    }

    template<int maxDiffOrder>
    typename std::enable_if<maxDiffOrder >= 2>::type
    evaluate(const Domain& coord, ElementaryHessian& r) const
    {
      // TODO: we currently require affine geometries.
      if (! element_->geometry().affine())
        DUNE_THROW(NotImplemented, "Due to missing features in the Geometry interface, "
          "the computation of higher derivatives (>=2) works only for affine transformations.");
      // get Jacobian of geometry
      const typename Element::Geometry::JacobianInverseTransposed
        JgeoIT = element_->geometry().jacobianInverseTransposed(coord);

      // TODO: we currently only implement the hessian...
      //       a proper implementation will require TMP magic.
      static const unsigned int dim = GridView::dimensionworld;
      // static_assert(
      //   isHessian<Range>::value,
      //   "We currently only higher order derivative we support is the Hessian of scalar functions");

      // get local hessian of the shape functions
      array<std::size_t, dim> directions;
      for(std::size_t i = 0; i < dim; ++i) {
        for(std::size_t j = i; j < dim; ++j) {
          directions[0] = 0;
          directions[1] = 0;
          directions[i]++;
          directions[j]++;
          lfs_.finiteElement().localBasis().evaluate(directions,coord,yb_);
          assert( yb_.size() == 1); // TODO: we only implement the hessian of scalar functions
          for(std::size_t n = 0; n < yb_.size(); ++n) {
            // sum up derivatives, weighting them with the appropriate coeff
            r[i][j] += xl_[i] * yb_[j];
          }
          // use symmetry of the hessian
          if (i != j) r[i][j] = r[j][i];
        }
      }
      // transform back to global coordinates
      for(std::size_t i = 0; i < dim; ++i)
        for(std::size_t j = i; j < dim; ++j)
          r[i][j] *= JgeoIT[i][j] * JgeoIT[i][j];

      return r;
    }

  protected:

    const shared_ptr<const GridFunctionSpace> pgfs_;
    const shared_ptr<const Vector> v_;
    LFS lfs_;
    LFSCache lfs_cache_;
    XView x_view_;
    mutable std::vector<ElementaryRange> xl_;
    mutable std::vector<ElementaryRange> yb_;
    const Element* element_;
  };

  DiscreteGridViewFunction(const GridFunctionSpace& gfs, const Vector& v)
    : pgfs_(stackobject_to_shared_ptr(gfs)),v_(stackobject_to_shared_ptr(v))
  {}

  DiscreteGridViewFunction(std::shared_ptr<const GridFunctionSpace> pgfs, std::shared_ptr<const Vector> v)
    : pgfs_(pgfs),v_(v)
  {}

  // this is part of the interface in dune-functions
  const Basis& basis() const
  {
    return *pgfs_;
  }
  const GridFunctionSpace& gridFunctionSpace() const
  {
    return *pgfs_;
  }

  const V& dofs() const
  {
    return *v_;
  }

  // TODO: Implement this using hierarchic search
  Range operator() (const Domain& x) const
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  friend DiscreteGridViewFunction<GFS,V,diffOrder+1> derivative(const DiscreteGridViewFunction& t)
  {
    return DiscreteGridViewFunction<GFS,V,diffOrder+1>(t.pgfs_, t.v_);
  }

  friend LocalFunction localFunction(const DiscreteGridViewFunction& t)
  {
    return LocalFunction(t.pgfs_, t.v_);
  }

  /**
   * \brief Get associated EntitySet
   */
  EntitySet entitySet() const
  {
    return pgfs_->gridView();
  }

private:

  const shared_ptr<const GridFunctionSpace> pgfs_;
  const shared_ptr<const Vector> v_;

};

} // end of namespace Dune::PDELab
} // end of namespace Dune

#endif // DUNE_PDELAB_FUNCTION_DISCRETEGRIDVIEWFUNCTION_HH
