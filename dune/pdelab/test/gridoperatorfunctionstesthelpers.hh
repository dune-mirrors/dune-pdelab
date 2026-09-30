// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_TEST_GRIDOPERATORFUNCTIONSTESTHELPERS_HH
#define DUNE_PDELAB_TEST_GRIDOPERATORFUNCTIONSTESTHELPERS_HH

/* Helpers for comparing assembly results of the PDELab GridOperator with
 * the GridOperator on a dune-functions basis spanning the
 * same discrete space.
 *
 * The DOF numberings of both frameworks differ. For Lagrange-type spaces
 * every DOF is identified by the values of two suitably chosen interpolated
 * functions ("keys"), which gives a permutation between the numberings.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

namespace Dune::PDELab::FunctionsBasisTest {

  //! entries of a PDELab ISTL vector (native BlockVector<FieldVector<K,1>>), flattened
  template<class NativeVector>
  std::vector<double> flatValues(const NativeVector& v)
  {
    std::vector<double> result;
    for (const auto& block : v)
      for (const auto& entry : block)
        result.push_back(entry);
    return result;
  }

  //! entries of a flat dune-functions compatible vector (BlockVector<K>)
  template<class Vector>
  std::vector<double> scalarValues(const Vector& v)
  {
    return std::vector<double>(v.begin(), v.end());
  }

  /** \brief permutation perm with perm[pdelabIndex] = functionsIndex
   *
   * Keys are rounded to a fixed precision; all keys must be pairwise distinct.
   */
  inline std::vector<std::size_t> matchDofs(const std::vector<double>& pdelabKey1, const std::vector<double>& pdelabKey2,
                                            const std::vector<double>& functionsKey1, const std::vector<double>& functionsKey2)
  {
    if (pdelabKey1.size() != functionsKey1.size())
      DUNE_THROW(Dune::Exception, "Spaces have different dimensions: " << pdelabKey1.size() << " vs. " << functionsKey1.size());

    auto key = [](double a, double b) {
      return std::make_pair(std::llround(a*1e8), std::llround(b*1e8));
    };

    std::map<std::pair<long long,long long>,std::size_t> functionsIndex;
    for (std::size_t i = 0; i < functionsKey1.size(); ++i)
      if (not functionsIndex.emplace(key(functionsKey1[i],functionsKey2[i]), i).second)
        DUNE_THROW(Dune::Exception, "DOF keys of the dune-functions basis are not unique");

    std::vector<std::size_t> perm(pdelabKey1.size());
    for (std::size_t i = 0; i < pdelabKey1.size(); ++i)
    {
      auto it = functionsIndex.find(key(pdelabKey1[i],pdelabKey2[i]));
      if (it == functionsIndex.end())
        DUNE_THROW(Dune::Exception, "PDELab DOF " << i << " has no counterpart in the dune-functions basis");
      perm[i] = it->second;
    }
    return perm;
  }

  /** \brief permutation perm with perm[pdelabIndex] = functionsIndex, matched per element
   *
   * For discontinuous spaces the keys are only unique within an element.
   * DOFs are matched element by element: the local DOFs of the PDELab space
   * (via LFSIndexCache) and of the dune-functions basis (via its LocalView)
   * are paired by their keys. Both spaces must use flat indices.
   */
  template<class GFS, class Basis>
  std::vector<std::size_t> matchDofsElementwise(const GFS& gfs, const Basis& basis,
                                                const std::vector<double>& pdelabKey1, const std::vector<double>& pdelabKey2,
                                                const std::vector<double>& functionsKey1, const std::vector<double>& functionsKey2)
  {
    if (pdelabKey1.size() != functionsKey1.size())
      DUNE_THROW(Dune::Exception, "Spaces have different dimensions: " << pdelabKey1.size() << " vs. " << functionsKey1.size());

    auto key = [](double a, double b) {
      return std::make_pair(std::llround(a*1e8), std::llround(b*1e8));
    };

    using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
    LFS lfs(gfs);
    Dune::PDELab::LFSIndexCache<LFS> cache(lfs);
    auto localView = basis.localView();

    const std::size_t unset = pdelabKey1.size();
    std::vector<std::size_t> perm(pdelabKey1.size(), unset);
    for (const auto& element : elements(basis.gridView()))
    {
      lfs.bind(element);
      cache.update();
      localView.bind(element);

      std::map<std::pair<long long,long long>,std::size_t> functionsIndex;
      for (std::size_t i = 0; i < localView.size(); ++i)
      {
        auto fi = localView.index(i)[0];
        if (not functionsIndex.emplace(key(functionsKey1[fi],functionsKey2[fi]), fi).second)
          DUNE_THROW(Dune::Exception, "DOF keys of the dune-functions basis are not unique within an element");
      }

      for (std::size_t i = 0; i < cache.size(); ++i)
      {
        auto pi = cache.containerIndex(i)[0];
        auto it = functionsIndex.find(key(pdelabKey1[pi],pdelabKey2[pi]));
        if (it == functionsIndex.end())
          DUNE_THROW(Dune::Exception, "PDELab DOF " << pi << " has no counterpart in the dune-functions basis");
        if (perm[pi] != unset and perm[pi] != it->second)
          DUNE_THROW(Dune::Exception, "PDELab DOF " << pi << " is matched inconsistently");
        perm[pi] = it->second;
      }
    }
    return perm;
  }

  //! fill a PDELab native vector with random values
  template<class NativeVector>
  void randomize(NativeVector& v, unsigned int seed)
  {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-1.0,1.0);
    for (auto& block : v)
      for (auto& entry : block)
        entry = dist(gen);
  }

  //! transfer a PDELab native vector to a flat dune-functions vector using perm
  template<class NativeVector, class Vector>
  void transfer(const NativeVector& pdelab, const std::vector<std::size_t>& perm, Vector& functions)
  {
    auto values = flatValues(pdelab);
    functions.resize(values.size());
    for (std::size_t i = 0; i < values.size(); ++i)
      functions[perm[i]] = values[i];
  }

  class Checker
  {
  public:
    explicit Checker(double tolerance) : tolerance_(tolerance) {}

    void check(bool condition, const std::string& what)
    {
      if (not condition)
      {
        std::cerr << "FAILED: " << what << std::endl;
        ++failures_;
      }
    }

    //! compare PDELab vector with flat dune-functions vector entrywise
    template<class NativeVector, class Vector>
    void compareVectors(const NativeVector& pdelab, const Vector& functions,
                        const std::vector<std::size_t>& perm, const std::string& what)
    {
      auto p = flatValues(pdelab);
      double maxDiff = 0.0, maxValue = 0.0;
      for (std::size_t i = 0; i < p.size(); ++i)
      {
        maxDiff = std::max(maxDiff, std::abs(p[i] - functions[perm[i]]));
        maxValue = std::max(maxValue, std::abs(p[i]));
      }
      std::cout << what << ": max |pdelab| = " << maxValue << ", max difference = " << maxDiff << std::endl;
      check(maxDiff <= tolerance_ * std::max(1.0, maxValue), what);
    }

    /** \brief compare PDELab matrix (BCRS of 1x1 blocks) with a scalar BCRS matrix entrywise
     *
     * Every nonzero entry of the PDELab matrix must match the corresponding
     * entry of the dune-functions matrix, and the Frobenius norms must agree
     * (so the dune-functions matrix contains no additional nonzeros).
     */
    template<class NativeMatrix, class Matrix>
    void compareMatrices(const NativeMatrix& pdelab, const Matrix& functions,
                         const std::vector<std::size_t>& perm, const std::string& what)
    {
      double maxDiff = 0.0, maxValue = 0.0, normP = 0.0, normF = 0.0;
      std::size_t missing = 0;
      for (auto row = pdelab.begin(); row != pdelab.end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
        {
          const double value = (*col)[0][0];
          normP += value*value;
          maxValue = std::max(maxValue, std::abs(value));
          const auto& functionsRow = functions[perm[row.index()]];
          auto it = functionsRow.find(perm[col.index()]);
          if (it == functionsRow.end())
          {
            if (value != 0.0)
              ++missing;
            continue;
          }
          maxDiff = std::max(maxDiff, std::abs(value - *it));
        }
      for (auto row = functions.begin(); row != functions.end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
          normF += (*col) * (*col);
      normP = std::sqrt(normP);
      normF = std::sqrt(normF);
      std::cout << what << ": nonzeros pdelab/functions = " << pdelab.nonzeroes() << "/" << functions.nonzeroes()
                << ", |A|_F = " << normP << "/" << normF
                << ", max difference = " << maxDiff << std::endl;
      check(missing == 0, what + " (entries missing from the pattern)");
      check(maxDiff <= tolerance_ * std::max(1.0, maxValue), what + " (entries)");
      check(std::abs(normP - normF) <= tolerance_ * std::max(1.0, normP), what + " (Frobenius norm)");
    }

    int result() const { return failures_ == 0 ? 0 : 1; }

  private:
    double tolerance_;
    int failures_ = 0;
  };

} // namespace Dune::PDELab::FunctionsBasisTest

#endif // DUNE_PDELAB_TEST_GRIDOPERATORFUNCTIONSTESTHELPERS_HH
