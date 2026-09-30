// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDOPERATOR_FUNCTIONS_LOCALOPERATORASSEMBLER_HH
#define DUNE_PDELAB_GRIDOPERATOR_FUNCTIONS_LOCALOPERATORASSEMBLER_HH

#include <cstddef>
#include <optional>
#include <type_traits>
#include <utility>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspacetags.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridoperator/common/localmatrix.hh>
#include <dune/pdelab/localoperator/callswitch.hh>

#include <dune/pdelab/gridoperator/functions/localfunctionspace.hh>

/** \file
 * \brief Local assemblers for dune-assembler wrapping PDELab local operators.
 *
 * The classes in this file model the local assembler concepts of
 * dune-assembler (Dune::Assembler::Concepts) and forward to a PDELab local
 * operator via LocalOperatorApply, i.e. with the flags and calling
 * conventions of the PDELab GridOperator:
 *
 * - ResidualLocalAssembler: vector assembler for the residual r(x)
 * - JacobianLocalAssembler: matrix assembler for the jacobian at x and its
 *   sparsity pattern
 * - JacobianApplyLocalAssembler: vector assembler for the application of
 *   the jacobian at x to z
 *
 * Use them with Dune::Assembler::Assembler (rows: test basis, columns:
 * trial basis), e.g.
 * \code
 * auto assembler = Dune::Assembler::Assembler(basis);
 * auto residual = makeResidualAssembler(lop, Functions::istlVectorBackend(x), basis);
 * assembler.assembleVectorEntries(residual, rBackend);
 * \endcode
 *
 * Coefficient vectors are passed as lightweight accessors providing
 * operator[] for the global multi-indices of the trial basis (e.g.
 * Functions::istlVectorBackend(x)); they are stored by value.
 *
 * Intersection-based methods are only provided if the local operator has
 * the corresponding terms, so dune-assembler skips intersections otherwise.
 *
 * Not supported (rejected at compile or run time):
 * - doSkeletonTwoSided and the *_volume_post_skeleton methods, which have no
 *   counterpart in the dune-assembler traversal
 * - periodic intersections (boundary() and neighbor()), which dune-assembler
 *   does not visit. Note that YaspGrid does not use this convention: its
 *   periodic neighbors are overlap copies (plain interior intersections) or,
 *   without overlap, intersections that are neither boundary nor neighbor,
 *   which are ignored as processor boundaries, as in the PDELab assembler.
 *
 * Local operators are called through a shared reference. PDELab local
 * operators often contain mutable caches, so use a sequential executor.
 *
 * \warning This is experimental.
 */

namespace Dune::PDELab::Impl::Functions {


    template<class LocalVector, class LFS, class Coefficients>
    void gather(LocalVector& xl, const LFS& lfs, const Coefficients& x)
    {
      xl.resize(lfs.size());
      for (std::size_t i = 0; i < lfs.size(); ++i)
        xl.base()[i] = x[lfs.index(i)];
    }

    template<class T, class Tag, class W, class Target>
    void addTo(const LocalVector<T,Tag,W>& rl, Target& v)
    {
      for (std::size_t i = 0; i < rl.size(); ++i)
        v[i] += rl.base()[i];
    }

    template<class T, class W, class Target>
    void addTo(const LocalMatrix<T,W>& al, Target& m)
    {
      for (std::size_t i = 0; i < al.nrows(); ++i)
        for (std::size_t j = 0; j < al.ncols(); ++j)
          m[i][j] += al.getEntry(i,j);
    }

    template<class Target>
    void addTo(const LocalSparsityPattern& pl, Target& p)
    {
      for (const auto& link : pl)
        p.add(link.i(), link.j());
    }

    template<class Coefficients, class MultiIndex>
    using FieldOf = std::decay_t<decltype(std::declval<const Coefficients&>()[std::declval<const MultiIndex&>()])>;

    // Position of the intersection within the intersection range of its
    // inside element, as used by PDELab's IntersectionGeometry.
    template<class GridView, class Intersection>
    unsigned int intersectionIndex(const GridView& gridView, const Intersection& intersection)
    {
      unsigned int index = 0;
      for (const auto& other : intersections(gridView, intersection.inside()))
      {
        if (other == intersection)
          return index;
        ++index;
      }
      DUNE_THROW(Dune::InvalidStateException, "Intersection not found in the intersection range of its inside element");
    }

    template<class GridView, class Element>
    void checkNoPeriodicIntersections(const GridView& gridView, const Element& element)
    {
      for (const auto& intersection : intersections(gridView, element))
        if (intersection.boundary() and intersection.neighbor())
          DUNE_THROW(Dune::NotImplemented, "Periodic intersections (boundary() and neighbor()) are not visited "
                     "by dune-assembler, so the skeleton and boundary terms of the local operator would be incomplete");
    }

    //! Common part of the local assemblers: local operator, weight and supported features
    template<class LOP, class TestBasis, class TrialBasis, class Field>
    class LocalOperatorAssemblerBase
    {
      static_assert(not LOP::doSkeletonTwoSided,
        "Local operators with doSkeletonTwoSided are not supported by the dune-assembler traversal");
      static_assert(not (LOP::doPatternVolumePostSkeleton or LOP::doAlphaVolumePostSkeleton or LOP::doLambdaVolumePostSkeleton),
        "Local operators with *_volume_post_skeleton methods are not supported by the dune-assembler traversal");

    protected:
      using GridView = typename TestBasis::GridView;
      using Element = typename GridView::template Codim<0>::Entity;
      using Intersection = typename GridView::Intersection;
      using TestLocalView = typename TestBasis::LocalView;
      using TrialLocalView = typename TrialBasis::LocalView;
      using TestSpace = LocalFunctionSpaceView<TestLocalView>;
      using TrialSpace = LocalFunctionSpaceView<TrialLocalView>;
      using EG = ElementGeometry<Element>;
      using IG = IntersectionGeometry<Intersection>;

    public:
      using LocalOperator = LOP;
      using WeightType = Field;

      explicit LocalOperatorAssemblerBase(const LOP& lop)
        : lop_(&lop)
      {}

      //! weight of all contributions, as in the PDELab local assemblers
      void setWeight(WeightType weight) { weight_ = weight; }
      WeightType weight() const { return weight_; }

      const LOP& localOperator() const { return *lop_; }

    protected:
      // bound state is not copied; dune-assembler copies local assemblers before binding
      LocalOperatorAssemblerBase(const LocalOperatorAssemblerBase& other)
        : lop_(other.lop_), weight_(other.weight_)
      {}

      template<class Space>
      void bindElementCommon(const Space& lfsv, const Element& element, bool needsIntersections)
      {
        EG eg(element);
        skip_ = LocalOperatorApply::skipEntity(*lop_, eg);
        if (needsIntersections and not skip_)
          checkNoPeriodicIntersections(lfsv.basis().gridView(), element);
      }

      bool skip() const { return skip_; }

      template<class Space>
      std::optional<IG> intersectionGeometry(const Space& lfsv, const Intersection& intersection) const
      {
        IG ig(intersection, intersectionIndex(lfsv.basis().gridView(), intersection));
        if (LocalOperatorApply::skipIntersection(*lop_, ig))
          return std::nullopt;
        return ig;
      }

      const LOP* lop_;
      WeightType weight_ = 1.0;
      bool skip_ = false;
    };

    /* Vector assemblers get the local views of the test basis only. The
     * trial space either shares the test local view (same basis object) or
     * uses an own LocalView of the trial basis.
     */
    template<class TestBasis, class TrialBasis>
    class TrialSpaceProvider
    {
      using TestLocalView = typename TestBasis::LocalView;
      using TrialLocalView = typename TrialBasis::LocalView;
      using TrialSpace = LocalFunctionSpaceView<TrialLocalView>;
      using Element = typename TestLocalView::Element;

    public:
      explicit TrialSpaceProvider(const TrialBasis& trialBasis) : trialBasis_(&trialBasis) {}
      TrialSpaceProvider(const TrialSpaceProvider& other) : trialBasis_(other.trialBasis_) {}

      void bindLocalView(const TestLocalView& testLocalView)
      {
        if constexpr (std::is_same_v<TestLocalView,TrialLocalView>)
          if (static_cast<const void*>(&testLocalView.globalBasis()) == static_cast<const void*>(trialBasis_))
          {
            shared_.emplace(testLocalView);
            return;
          }
        own_.emplace(*trialBasis_);
      }

      void bind(const Element& element)
      {
        if (own_)
          own_->bind(element);
      }

      const TrialSpace& space() const { return own_ ? own_->view() : *shared_; }

    private:
      const TrialBasis* trialBasis_;
      std::optional<TrialSpace> shared_;
      std::optional<LocalFunctionSpace<TrialBasis>> own_;
    };



  /** \brief dune-assembler vector assembler for the residual of a PDELab local operator
   *
   * \tparam LOP        PDELab local operator
   * \tparam X          accessor for the coefficients of the trial function (x[multiIndex])
   * \tparam TestBasis  dune-functions basis of the test space (rows)
   * \tparam TrialBasis dune-functions basis of the trial space
   */
  template<class LOP, class X, class TestBasis, class TrialBasis = TestBasis>
  class ResidualLocalAssembler
    : public LocalOperatorAssemblerBase<LOP, TestBasis, TrialBasis,
                                              FieldOf<X, typename TrialBasis::LocalView::MultiIndex>>
  {
    using Field = FieldOf<X, typename TrialBasis::LocalView::MultiIndex>;
    using Base = LocalOperatorAssemblerBase<LOP, TestBasis, TrialBasis, Field>;
    using typename Base::Element;
    using typename Base::Intersection;
    using typename Base::TestLocalView;
    using typename Base::TestSpace;
    using typename Base::EG;

    static constexpr bool skeleton = LOP::doAlphaSkeleton or LOP::doLambdaSkeleton;
    static constexpr bool boundary = LOP::doAlphaBoundary or LOP::doLambdaBoundary;

  public:
    ResidualLocalAssembler(const LOP& lop, const X& x, const TrialBasis& trialBasis)
      : Base(lop), x_(x), trial_(trialBasis), trialOutside_(trialBasis)
    {}

    ResidualLocalAssembler(const ResidualLocalAssembler& other)
      : Base(other), x_(other.x_), trial_(other.trial_), trialOutside_(other.trialOutside_)
    {}

    void bindLocalView(const TestLocalView& localView)
    {
      lfsv_.emplace(localView);
      trial_.bindLocalView(localView);
    }

    void bindOutsideLocalView(const TestLocalView& localView)
      requires skeleton
    {
      lfsvn_.emplace(localView);
      trialOutside_.bindLocalView(localView);
    }

    void bindElement(const Element& element)
    {
      this->bindElementCommon(*lfsv_, element, skeleton or boundary);
      if (this->skip())
        return;
      trial_.bind(element);
      gather(xl_, trial_.space(), x_);
    }

    template<class TargetVector>
    void assembleElementVector(const Element& element, TargetVector& localVector)
    {
      if (this->skip())
        return;
      const auto& lfsu = trial_.space();
      rl_.assign(lfsv_->size(), 0.0);
      auto rv = rl_.weightedAccumulationView(this->weight());
      EG eg(element);
      LocalOperatorApply::alphaVolume(*this->lop_, eg, lfsu.tree(), xl_, lfsv_->tree(), rv);
      LocalOperatorApply::lambdaVolume(*this->lop_, eg, lfsv_->tree(), rv);
      addTo(rl_, localVector);
    }

    template<class TargetVector>
    void assembleBoundaryIntersectionVector(const Intersection& intersection, TargetVector& localVector)
      requires boundary
    {
      if (this->skip())
        return;
      auto ig = this->intersectionGeometry(*lfsv_, intersection);
      if (not ig)
        return;
      const auto& lfsu = trial_.space();
      rl_.assign(lfsv_->size(), 0.0);
      auto rv = rl_.weightedAccumulationView(this->weight());
      LocalOperatorApply::alphaBoundary(*this->lop_, *ig, lfsu.tree(), xl_, lfsv_->tree(), rv);
      LocalOperatorApply::lambdaBoundary(*this->lop_, *ig, lfsv_->tree(), rv);
      addTo(rl_, localVector);
    }

    template<class TargetVectors>
    void assembleInteriorIntersectionVector(const Intersection& intersection, TargetVectors& localVectors)
      requires skeleton
    {
      if (this->skip())
        return;
      auto ig = this->intersectionGeometry(*lfsv_, intersection);
      if (not ig)
        return;
      const auto& lfsu = trial_.space();
      trialOutside_.bind(intersection.outside());
      const auto& lfsun = trialOutside_.space();
      gather(xn_, lfsun, x_);
      rl_.assign(lfsv_->size(), 0.0);
      rn_.assign(lfsvn_->size(), 0.0);
      auto rv = rl_.weightedAccumulationView(this->weight());
      auto rnv = rn_.weightedAccumulationView(this->weight());
      LocalOperatorApply::alphaSkeleton(*this->lop_, *ig, lfsu.tree(), xl_, lfsv_->tree(),
                                        lfsun.tree(), xn_, lfsvn_->tree(), rv, rnv);
      LocalOperatorApply::lambdaSkeleton(*this->lop_, *ig, lfsv_->tree(), lfsvn_->tree(), rv, rnv);
      addTo(rl_, localVectors[0]);
      addTo(rn_, localVectors[1]);
    }

  private:
    X x_;
    TrialSpaceProvider<TestBasis,TrialBasis> trial_, trialOutside_;
    std::optional<TestSpace> lfsv_, lfsvn_;
    LocalVector<Field,TrialSpaceTag> xl_, xn_;
    LocalVector<Field,TestSpaceTag> rl_, rn_;
  };


  /** \brief dune-assembler matrix assembler for the jacobian of a PDELab local operator
   *
   * Provides the jacobian at x and the sparsity pattern of the local operator.
   *
   * \tparam LOP        PDELab local operator
   * \tparam X          accessor for the coefficients of the trial function (x[multiIndex])
   * \tparam TestBasis  dune-functions basis of the test space (rows)
   * \tparam TrialBasis dune-functions basis of the trial space (columns)
   */
  template<class LOP, class X, class TestBasis, class TrialBasis = TestBasis>
  class JacobianLocalAssembler
    : public LocalOperatorAssemblerBase<LOP, TestBasis, TrialBasis,
                                              FieldOf<X, typename TrialBasis::LocalView::MultiIndex>>
  {
    using Field = FieldOf<X, typename TrialBasis::LocalView::MultiIndex>;
    using Base = LocalOperatorAssemblerBase<LOP, TestBasis, TrialBasis, Field>;
    using typename Base::Element;
    using typename Base::Intersection;
    using typename Base::TestLocalView;
    using typename Base::TrialLocalView;
    using typename Base::TestSpace;
    using typename Base::TrialSpace;
    using typename Base::EG;

    static constexpr bool skeleton = LOP::doAlphaSkeleton or LOP::doPatternSkeleton;
    static constexpr bool boundary = LOP::doAlphaBoundary or LOP::doPatternBoundary;

  public:
    JacobianLocalAssembler(const LOP& lop, const X& x)
      : Base(lop), x_(x)
    {}

    JacobianLocalAssembler(const JacobianLocalAssembler& other)
      : Base(other), x_(other.x_)
    {}

    void bindLocalViews(const TestLocalView& testLocalView, const TrialLocalView& trialLocalView)
    {
      lfsv_.emplace(testLocalView);
      lfsu_.emplace(trialLocalView);
    }

    void bindOutsideLocalViews(const TestLocalView& testLocalView, const TrialLocalView& trialLocalView)
      requires skeleton
    {
      lfsvn_.emplace(testLocalView);
      lfsun_.emplace(trialLocalView);
    }

    void bindElement(const Element& element)
    {
      this->bindElementCommon(*lfsv_, element, skeleton or boundary);
      gathered_ = false;
    }

    template<class TargetMatrix>
    void assembleElementMatrix(const Element& element, TargetMatrix& localMatrix)
    {
      if (this->skip())
        return;
      gatherInside();
      a_ss_.assign(lfsv_->size(), lfsu_->size(), 0.0);
      auto av = a_ss_.weightedAccumulationView(this->weight());
      EG eg(element);
      LocalOperatorApply::jacobianVolume(*this->lop_, eg, lfsu_->tree(), xl_, lfsv_->tree(), av);
      addTo(a_ss_, localMatrix);
    }

    template<class TargetPattern>
    void assembleElementMatrixPattern(const Element&, TargetPattern& localPattern)
    {
      if (this->skip())
        return;
      p_ss_.clear();
      LocalOperatorApply::patternVolume(*this->lop_, lfsu_->tree(), lfsv_->tree(), p_ss_);
      addTo(p_ss_, localPattern);
    }

    template<class TargetMatrix>
    void assembleBoundaryIntersectionMatrix(const Intersection& intersection, TargetMatrix& localMatrix)
      requires boundary
    {
      if (this->skip())
        return;
      auto ig = this->intersectionGeometry(*lfsv_, intersection);
      if (not ig)
        return;
      gatherInside();
      a_ss_.assign(lfsv_->size(), lfsu_->size(), 0.0);
      auto av = a_ss_.weightedAccumulationView(this->weight());
      LocalOperatorApply::jacobianBoundary(*this->lop_, *ig, lfsu_->tree(), xl_, lfsv_->tree(), av);
      addTo(a_ss_, localMatrix);
    }

    template<class TargetPattern>
    void assembleBoundaryIntersectionMatrixPattern(const Intersection& intersection, TargetPattern& localPattern)
      requires boundary
    {
      if (this->skip())
        return;
      if (not this->intersectionGeometry(*lfsv_, intersection))
        return;
      p_ss_.clear();
      LocalOperatorApply::patternBoundary(*this->lop_, lfsu_->tree(), lfsv_->tree(), p_ss_);
      addTo(p_ss_, localPattern);
    }

    template<class TargetMatrices>
    void assembleInteriorIntersectionMatrix(const Intersection& intersection, TargetMatrices& localMatrices)
      requires skeleton
    {
      if (this->skip())
        return;
      auto ig = this->intersectionGeometry(*lfsv_, intersection);
      if (not ig)
        return;
      gatherInside();
      gather(xn_, *lfsun_, x_);
      a_ss_.assign(lfsv_->size(), lfsu_->size(), 0.0);
      a_sn_.assign(lfsv_->size(), lfsun_->size(), 0.0);
      a_ns_.assign(lfsvn_->size(), lfsu_->size(), 0.0);
      a_nn_.assign(lfsvn_->size(), lfsun_->size(), 0.0);
      auto av_ss = a_ss_.weightedAccumulationView(this->weight());
      auto av_sn = a_sn_.weightedAccumulationView(this->weight());
      auto av_ns = a_ns_.weightedAccumulationView(this->weight());
      auto av_nn = a_nn_.weightedAccumulationView(this->weight());
      LocalOperatorApply::jacobianSkeleton(*this->lop_, *ig, lfsu_->tree(), xl_, lfsv_->tree(),
                                           lfsun_->tree(), xn_, lfsvn_->tree(),
                                           av_ss, av_sn, av_ns, av_nn);
      addTo(a_ss_, localMatrices[0][0]);
      addTo(a_sn_, localMatrices[0][1]);
      addTo(a_ns_, localMatrices[1][0]);
      addTo(a_nn_, localMatrices[1][1]);
    }

    template<class TargetPatterns>
    void assembleInteriorIntersectionMatrixPattern(const Intersection& intersection, TargetPatterns& localPatterns)
      requires skeleton
    {
      if (this->skip())
        return;
      if (not this->intersectionGeometry(*lfsv_, intersection))
        return;
      p_sn_.clear();
      p_ns_.clear();
      p_nn_.clear();
      LocalOperatorApply::patternSkeleton(*this->lop_, lfsu_->tree(), lfsv_->tree(),
                                          lfsun_->tree(), lfsvn_->tree(), p_sn_, p_ns_);
      // PDELab's pattern_skeleton only provides the couplings between inside
      // and outside; the outside-outside block is covered by the volume pattern
      // of the outside element. dune-assembler scatters each block of the
      // skeleton matrices with its own local pattern, so provide it here.
      // (The inside-inside block shares the local pattern of the element.)
      LocalOperatorApply::patternVolume(*this->lop_, lfsun_->tree(), lfsvn_->tree(), p_nn_);
      addTo(p_sn_, localPatterns[0][1]);
      addTo(p_ns_, localPatterns[1][0]);
      addTo(p_nn_, localPatterns[1][1]);
    }

  private:
    // the coefficients are only needed for matrix entries, not for the pattern
    void gatherInside()
    {
      if (not gathered_)
        gather(xl_, *lfsu_, x_);
      gathered_ = true;
    }

    X x_;
    bool gathered_ = false;
    std::optional<TestSpace> lfsv_, lfsvn_;
    std::optional<TrialSpace> lfsu_, lfsun_;
    LocalVector<Field,TrialSpaceTag> xl_, xn_;
    LocalMatrix<Field> a_ss_, a_sn_, a_ns_, a_nn_;
    LocalSparsityPattern p_ss_, p_sn_, p_ns_, p_nn_;
  };


  /** \brief dune-assembler vector assembler applying the jacobian of a PDELab local operator
   *
   * Assembles the application of the jacobian at x to z. For linear local
   * operators x is ignored.
   *
   * \tparam LOP        PDELab local operator
   * \tparam X          accessor for the linearization point (x[multiIndex])
   * \tparam Z          accessor for the vector the jacobian is applied to
   * \tparam TestBasis  dune-functions basis of the test space (rows)
   * \tparam TrialBasis dune-functions basis of the trial space
   */
  template<class LOP, class X, class Z, class TestBasis, class TrialBasis = TestBasis>
  class JacobianApplyLocalAssembler
    : public LocalOperatorAssemblerBase<LOP, TestBasis, TrialBasis,
                                              FieldOf<X, typename TrialBasis::LocalView::MultiIndex>>
  {
    using Field = FieldOf<X, typename TrialBasis::LocalView::MultiIndex>;
    using Base = LocalOperatorAssemblerBase<LOP, TestBasis, TrialBasis, Field>;
    using typename Base::Element;
    using typename Base::Intersection;
    using typename Base::TestLocalView;
    using typename Base::TestSpace;
    using typename Base::EG;

    static constexpr bool skeleton = LOP::doAlphaSkeleton;
    static constexpr bool boundary = LOP::doAlphaBoundary;

  public:
    JacobianApplyLocalAssembler(const LOP& lop, const X& x, const Z& z, const TrialBasis& trialBasis)
      : Base(lop), x_(x), z_(z), trial_(trialBasis), trialOutside_(trialBasis)
    {}

    JacobianApplyLocalAssembler(const JacobianApplyLocalAssembler& other)
      : Base(other), x_(other.x_), z_(other.z_), trial_(other.trial_), trialOutside_(other.trialOutside_)
    {}

    void bindLocalView(const TestLocalView& localView)
    {
      lfsv_.emplace(localView);
      trial_.bindLocalView(localView);
    }

    void bindOutsideLocalView(const TestLocalView& localView)
      requires skeleton
    {
      lfsvn_.emplace(localView);
      trialOutside_.bindLocalView(localView);
    }

    void bindElement(const Element& element)
    {
      this->bindElementCommon(*lfsv_, element, skeleton or boundary);
      if (this->skip())
        return;
      trial_.bind(element);
      gather(xl_, trial_.space(), x_);
      gather(zl_, trial_.space(), z_);
    }

    template<class TargetVector>
    void assembleElementVector(const Element& element, TargetVector& localVector)
    {
      if (this->skip())
        return;
      const auto& lfsu = trial_.space();
      yl_.assign(lfsv_->size(), 0.0);
      auto yv = yl_.weightedAccumulationView(this->weight());
      EG eg(element);
      LocalOperatorApply::jacobianApplyVolume(*this->lop_, eg, lfsu.tree(), xl_, zl_, lfsv_->tree(), yv);
      addTo(yl_, localVector);
    }

    template<class TargetVector>
    void assembleBoundaryIntersectionVector(const Intersection& intersection, TargetVector& localVector)
      requires boundary
    {
      if (this->skip())
        return;
      auto ig = this->intersectionGeometry(*lfsv_, intersection);
      if (not ig)
        return;
      const auto& lfsu = trial_.space();
      yl_.assign(lfsv_->size(), 0.0);
      auto yv = yl_.weightedAccumulationView(this->weight());
      LocalOperatorApply::jacobianApplyBoundary(*this->lop_, *ig, lfsu.tree(), xl_, zl_, lfsv_->tree(), yv);
      addTo(yl_, localVector);
    }

    template<class TargetVectors>
    void assembleInteriorIntersectionVector(const Intersection& intersection, TargetVectors& localVectors)
      requires skeleton
    {
      if (this->skip())
        return;
      auto ig = this->intersectionGeometry(*lfsv_, intersection);
      if (not ig)
        return;
      const auto& lfsu = trial_.space();
      trialOutside_.bind(intersection.outside());
      const auto& lfsun = trialOutside_.space();
      gather(xn_, lfsun, x_);
      gather(zn_, lfsun, z_);
      yl_.assign(lfsv_->size(), 0.0);
      yn_.assign(lfsvn_->size(), 0.0);
      auto yv = yl_.weightedAccumulationView(this->weight());
      auto ynv = yn_.weightedAccumulationView(this->weight());
      LocalOperatorApply::jacobianApplySkeleton(*this->lop_, *ig, lfsu.tree(), xl_, zl_, lfsv_->tree(),
                                                lfsun.tree(), xn_, zn_, lfsvn_->tree(), yv, ynv);
      addTo(yl_, localVectors[0]);
      addTo(yn_, localVectors[1]);
    }

  private:
    X x_;
    Z z_;
    TrialSpaceProvider<TestBasis,TrialBasis> trial_, trialOutside_;
    std::optional<TestSpace> lfsv_, lfsvn_;
    LocalVector<Field,TrialSpaceTag> xl_, xn_, zl_, zn_;
    LocalVector<Field,TestSpaceTag> yl_, yn_;
  };


  //! local assembler for the residual r(x), test and trial space given by the same basis
  template<class LOP, class X, class Basis>
  auto makeResidualAssembler(const LOP& lop, const X& x, const Basis& basis)
  {
    return ResidualLocalAssembler<LOP,X,Basis>(lop, x, basis);
  }

  //! local assembler for the residual r(x) with different test and trial bases
  template<class LOP, class X, class TestBasis, class TrialBasis>
  auto makeResidualAssembler(const LOP& lop, const X& x, const TestBasis&, const TrialBasis& trialBasis)
  {
    return ResidualLocalAssembler<LOP,X,TestBasis,TrialBasis>(lop, x, trialBasis);
  }

  //! local assembler for the jacobian at x and its pattern, test and trial space given by the same basis
  template<class LOP, class X, class Basis>
  auto makeJacobianAssembler(const LOP& lop, const X& x, const Basis&)
  {
    return JacobianLocalAssembler<LOP,X,Basis>(lop, x);
  }

  //! local assembler for the jacobian at x and its pattern with different test and trial bases
  template<class LOP, class X, class TestBasis, class TrialBasis>
  auto makeJacobianAssembler(const LOP& lop, const X& x, const TestBasis&, const TrialBasis&)
  {
    return JacobianLocalAssembler<LOP,X,TestBasis,TrialBasis>(lop, x);
  }

  //! local assembler applying the jacobian at x to z, test and trial space given by the same basis
  template<class LOP, class X, class Z, class Basis>
  auto makeJacobianApplyAssembler(const LOP& lop, const X& x, const Z& z, const Basis& basis)
  {
    return JacobianApplyLocalAssembler<LOP,X,Z,Basis>(lop, x, z, basis);
  }

  //! local assembler applying the (constant) jacobian of a linear local operator to z
  template<class LOP, class Z, class Basis>
    requires (bool(LOP::isLinear))
  auto makeJacobianApplyAssembler(const LOP& lop, const Z& z, const Basis& basis)
  {
    return JacobianApplyLocalAssembler<LOP,Z,Z,Basis>(lop, z, z, basis);
  }

} // namespace Dune::PDELab::Impl::Functions

#endif // DUNE_PDELAB_GRIDOPERATOR_FUNCTIONS_LOCALOPERATORASSEMBLER_HH
