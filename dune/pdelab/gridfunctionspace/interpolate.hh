// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_INTERPOLATE_HH
#define DUNE_PDELAB_INTERPOLATE_HH

#include <vector>
#include <utility>

#include <dune/common/exceptions.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/functions/common/functionfromcallable.hh>

#include <dune/typetree/typetree.hh>
#include <dune/typetree/pairtraversal.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/function/localfunction.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    // Backend for standard local interpolation
    struct InterpolateBackendStandard
    {
      template<typename FE, typename ElemFunction, typename XL>
      void interpolate(const FE &fe, const ElemFunction &elemFunction,
                       XL &xl) const
      {
        FiniteElementInterfaceSwitch<FE>::interpolation(fe).
          interpolate(elemFunction,xl);
      }
    };

    namespace {

      template<typename IB, typename LF, typename XG>
      struct InterpolateLeafFromScalarVisitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename LFS, typename TreePath>
        void leaf(const LFS& lfs, TreePath treePath) const
        {
          std::vector<typename XG::ElementType> xl(lfs.size());

          // call interpolate for the basis
          ib.interpolate(lfs.finiteElement(), lf, xl);

          // write coefficients into local vector
          xg.write_sub_container(lfs,xl);
        }

        InterpolateLeafFromScalarVisitor(const IB& ib_, const LF& lf_, XG& xg_)
          : ib(ib_)
          , lf(lf_)
          , xg(xg_)
        {}

        const IB& ib;
        const LF& lf;
        XG& xg;

      };


      template<typename IB, typename LF, typename XG>
      struct InterpolateLeafFromVectorVisitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {
        using Domain = typename Functions::SignatureTraits<LF>::Domain;
        using Range = typename Functions::SignatureTraits<LF>::Range;
        using RangeField = typename FieldTraits<Range>::field_type;

        template<typename LFS, typename TreePath>
        void leaf(const LFS& lfs, TreePath treePath) const
        {
          std::vector<typename XG::ElementType> xl(lfs.size());

          // call interpolate for the basis
          auto f = [&](const Domain& x) -> RangeField { return lf(x)[index]; };

          using LocalFunction = typename Dune::Functions::FunctionFromCallable<RangeField(Domain), decltype(f), TypeTree::EmptyNode>;
          LocalFunction fnkt(f);

          ib.interpolate(lfs.finiteElement(), fnkt, xl);

          // write coefficients into local vector
          xg.write_sub_container(lfs,xl);

          // increment index
          assert(index == treePath.back());
          index++;
        }

        InterpolateLeafFromVectorVisitor(const IB& ib_, const LF& lf_, XG& xg_)
          : ib(ib_)
          , lf(lf_)
          , index(0)
          , xg(xg_)
        {}

        const IB& ib;
        const LF& lf;
        std::size_t index;
        XG& xg;

      };


      template<typename IB, typename E, typename XG>
      struct InterpolateVisitor
        : public TypeTree::TreePairVisitor
        , public TypeTree::DynamicTraversal
      {
        template<typename F, typename LFS, typename TreePath>
        typename enable_if<F::isLeaf && LFS::isLeaf>::type
        leaf(const F& f, const LFS& lfs, TreePath treePath) const
        {
          std::vector<typename XG::ElementType> xl(lfs.size());

          // call interpolate for the basis
          using Domain = typename Functions::SignatureTraits<F>::Domain;
          using Range = typename Functions::SignatureTraits<F>::Range;
          using LocalFunction = typename Dune::Functions::FunctionFromCallable<Range(Domain), F, TypeTree::EmptyNode>;
          LocalFunction lf(f);
          ib.interpolate(lfs.finiteElement(), lf, xl);

          // write coefficients into local vector
          xg.write_sub_container(lfs,xl);
        }

        // interpolate PowerLFS from scalar-valued function
        template<typename F, typename LFS, typename TreePath,
                 typename Range = typename Functions::SignatureTraits<F>::Range>
        typename enable_if<F::isLeaf &&
                           std::is_convertible<Range, typename FieldTraits< Range >::field_type>::value &&
                           (!LFS::isLeaf)>::type
        leaf(const F& f, const LFS& lfs, TreePath treePath) const
        {
          static_assert((TypeTree::TreeInfo<LFS>::depth == 2),
                        "Automatic interpolation of vector-valued function " \
                        "is restricted to trees of depth 1");

          // call interpolate for the basis
          using Domain = typename Functions::SignatureTraits<F>::Domain;
          using LocalFunction = typename Dune::Functions::FunctionFromCallable<Range(Domain), F, TypeTree::EmptyNode>;
          LocalFunction lf(f);

          TypeTree::applyToTree(lfs,InterpolateLeafFromScalarVisitor<IB,LocalFunction,XG>(ib,lf,xg));

        }

        // interpolate PowerLFS from vector-valued function
        template<typename F, typename LFS, typename TreePath,
                 typename Range = typename Functions::SignatureTraits<F>::Range>
        typename enable_if<F::isLeaf &&
                          (!std::is_convertible<Range, typename FieldTraits< Range >::field_type>::value) &&
                          (!LFS::isLeaf)>::type
        leaf(const F& f, const LFS& lfs, TreePath treePath) const
        {
          static_assert((TypeTree::TreeInfo<LFS>::depth == 2),
                        "Automatic interpolation of vector-valued function " \
                        "is restricted to trees of depth 1");
          static_assert((TypeTree::TreeInfo<LFS>::leafCount == Range::dimension),
                        "Number of leaf nodes and dimension of range type " \
                        "must match for automatic interpolation of "    \
                        "vector-valued function");

          TypeTree::applyToTree(lfs,InterpolateLeafFromVectorVisitor<IB,F,XG>(ib,f,xg));
        }

        InterpolateVisitor(IB ib_, const E& e_, XG& xg_)
          : ib(ib_)
          , e(e_)
          , xg(xg_)
        {}

      private:
        IB ib;
        const E& e;
        XG& xg;
      };

    } // anonymous namespace

    //! interpolation from a given grid function
    /**
     * \code
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
     * \endcode
     * \param f   Function to interpolate from.
     * \param gfs GridFunctionSpace to use for interpoaltion.
     * \param xg  Global vector of dofs to interpolate into.
     *
     * \note \c xg needs to be initialized to the correct size, but there is
     *       no need to initialize its contents.
     */
    template<typename F, typename GFS, typename XG>
    void interpolate (const F& f, const GFS& gfs, XG& xg)
    {
      // get some types
      typedef typename GFS::Traits::GridViewType GV;
      typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
      typedef typename GV::Traits::template Codim<0>::Entity Element;

      // make local function
      auto lf = makeLocalFunctionTree(f, gfs.gridView());

      // make local function space
      typedef LocalFunctionSpace<GFS> LFS;
      LFS lfs(gfs);
      typedef LFSIndexCache<LFS> LFSCache;
      LFSCache lfs_cache(lfs);
      typedef typename XG::template LocalView<LFSCache> XView;

      XView x_view(xg);

      // loop once over the grid
      for (ElementIterator it = gfs.gridView().template begin<0>();
           it!=gfs.gridView().template end<0>(); ++it)
        {
          // bind local function space to element
          lf.bind(*it);
          lfs.bind(*it);
          lfs_cache.update();
          x_view.bind(lfs_cache);

          // call interpolate
          TypeTree::applyToTreePair(lf,lfs,InterpolateVisitor<InterpolateBackendStandard,Element,XView>(InterpolateBackendStandard(),*it,x_view));

          x_view.unbind();
          lf.unbind();
        }

      x_view.detach();
    }

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
