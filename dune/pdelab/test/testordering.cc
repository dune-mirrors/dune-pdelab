// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/filledarray.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#include <dune/grid/common/scsgmapper.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/pdelab.hh>
#include <dune/pdelab/common/tree2graphviz.hh>

#include "gridexamples.hh"

// hard coded exceptions for sparse patterns
static auto is_ordering_dense = [](auto& test_case, auto& suffix) {
  // When blocking local orderings it can happen that patterns are not
  // dense anymore: consider a composite node with (p1, p0) that is entity
  // blocked (e.g. BlockVector<FieldVector<FEMBlock,2>>). The sub-blocks
  // (FieldVector<FEMBlock,2>) for p1 will be empty in codim 0 entities and
  // likewise for p0 sub-blocks for codim 1 entities.

  // Thus, we can only be sure that a
  // pattern has to be dense in leaf and root nodes. In that case, we
  // conform ourselves to know that the sampled pattern is not out of
  // bounds
  // The test cases marked as not dense are hard coded after manual inspection.

  if (test_case >= 108 or test_case <= 112) {
    // GFS: Dune::PDELab::OrderedGridFunctionSpace<Dune::PDELab::UnorderedCompositeGridFunctionSpace<Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)2>, Dune::PDELab::EntityBlockedOrderingTag, Dune::PDELab::UnorderedGridFunctionSpace<Dune::PDELab::PartitionViewEntitySet<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, Dune::PartitionSet<31u> >, Dune::PDELab::VariableMonomLocalFiniteElementMap<Dune::SingleCodimSingleGeomTypeMapper<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, 0>, float, double, 2, 6>, Dune::PDELab::NoConstraints, Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)1>, Dune::PDELab::LeafOrderingTag<Dune::PDELab::EmptyParams> >, Dune::PDELab::UnorderedGridFunctionSpace<Dune::PDELab::PartitionViewEntitySet<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, Dune::PartitionSet<31u> >, Dune::PDELab::QkLocalFiniteElementMap<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, float, double, 2ul>, Dune::PDELab::NoConstraints, Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)1>, Dune::PDELab::LeafOrderingTag<Dune::PDELab::EmptyParams> > > >
    // GFS: Dune::PDELab::OrderedGridFunctionSpace<Dune::PDELab::UnorderedPowerGridFunctionSpace<Dune::PDELab::UnorderedCompositeGridFunctionSpace<Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)2>, Dune::PDELab::EntityBlockedOrderingTag, Dune::PDELab::UnorderedGridFunctionSpace<Dune::PDELab::PartitionViewEntitySet<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, Dune::PartitionSet<31u> >, Dune::PDELab::VariableMonomLocalFiniteElementMap<Dune::SingleCodimSingleGeomTypeMapper<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, 0>, float, double, 2, 6>, Dune::PDELab::NoConstraints, Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)1>, Dune::PDELab::LeafOrderingTag<Dune::PDELab::EmptyParams> >, Dune::PDELab::UnorderedGridFunctionSpace<Dune::PDELab::PartitionViewEntitySet<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, Dune::PartitionSet<31u> >, Dune::PDELab::QkLocalFiniteElementMap<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, float, double, 2ul>, Dune::PDELab::NoConstraints, Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)1>, Dune::PDELab::LeafOrderingTag<Dune::PDELab::EmptyParams> > >, 3ul, Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)0>, Dune::PDELab::LexicographicOrderingTag> >
    // GFS: Dune::PDELab::OrderedGridFunctionSpace<Dune::PDELab::UnorderedPowerGridFunctionSpace<Dune::PDELab::UnorderedCompositeGridFunctionSpace<Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)2>, Dune::PDELab::EntityBlockedOrderingTag, Dune::PDELab::UnorderedGridFunctionSpace<Dune::PDELab::PartitionViewEntitySet<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, Dune::PartitionSet<31u> >, Dune::PDELab::VariableMonomLocalFiniteElementMap<Dune::SingleCodimSingleGeomTypeMapper<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, 0>, float, double, 2, 6>, Dune::PDELab::NoConstraints, Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)1>, Dune::PDELab::LeafOrderingTag<Dune::PDELab::EmptyParams> >, Dune::PDELab::UnorderedGridFunctionSpace<Dune::PDELab::PartitionViewEntitySet<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, Dune::PartitionSet<31u> >, Dune::PDELab::QkLocalFiniteElementMap<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, float, double, 2ul>, Dune::PDELab::NoConstraints, Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)1>, Dune::PDELab::LeafOrderingTag<Dune::PDELab::EmptyParams> > >, 3ul, Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)2>, Dune::PDELab::LexicographicOrderingTag> >
    // GFS: Dune::PDELab::OrderedGridFunctionSpace<Dune::PDELab::UnorderedPowerGridFunctionSpace<Dune::PDELab::UnorderedCompositeGridFunctionSpace<Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)2>, Dune::PDELab::EntityBlockedOrderingTag, Dune::PDELab::UnorderedGridFunctionSpace<Dune::PDELab::PartitionViewEntitySet<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, Dune::PartitionSet<31u> >, Dune::PDELab::VariableMonomLocalFiniteElementMap<Dune::SingleCodimSingleGeomTypeMapper<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, 0>, float, double, 2, 6>, Dune::PDELab::NoConstraints, Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)1>, Dune::PDELab::LeafOrderingTag<Dune::PDELab::EmptyParams> >, Dune::PDELab::UnorderedGridFunctionSpace<Dune::PDELab::PartitionViewEntitySet<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, Dune::PartitionSet<31u> >, Dune::PDELab::QkLocalFiniteElementMap<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, float, double, 2ul>, Dune::PDELab::NoConstraints, Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)1>, Dune::PDELab::LeafOrderingTag<Dune::PDELab::EmptyParams> > >, 3ul, Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)0>, Dune::PDELab::EntityBlockedOrderingTag> >
    // GFS: Dune::PDELab::OrderedGridFunctionSpace<Dune::PDELab::UnorderedPowerGridFunctionSpace<Dune::PDELab::UnorderedCompositeGridFunctionSpace<Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)2>, Dune::PDELab::EntityBlockedOrderingTag, Dune::PDELab::UnorderedGridFunctionSpace<Dune::PDELab::PartitionViewEntitySet<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, Dune::PartitionSet<31u> >, Dune::PDELab::VariableMonomLocalFiniteElementMap<Dune::SingleCodimSingleGeomTypeMapper<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, 0>, float, double, 2, 6>, Dune::PDELab::NoConstraints, Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)1>, Dune::PDELab::LeafOrderingTag<Dune::PDELab::EmptyParams> >, Dune::PDELab::UnorderedGridFunctionSpace<Dune::PDELab::PartitionViewEntitySet<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, Dune::PartitionSet<31u> >, Dune::PDELab::QkLocalFiniteElementMap<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> >, float, double, 2ul>, Dune::PDELab::NoConstraints, Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)1>, Dune::PDELab::LeafOrderingTag<Dune::PDELab::EmptyParams> > >, 3ul, Dune::PDELab::ISTL::BackendOptions<(Dune::PDELab::ISTL::Blocking)2>, Dune::PDELab::EntityBlockedOrderingTag> >
    std::size_t pos = (test_case == 110 or test_case == 112) ? 2 : 1;
    if (suffix.size() == pos and (suffix.front() % 25) < 21){
      // reason: monomial has no DOFs in codim != 0 and there exist blocking in composite node with q2
      // in cases 110 and 112, power of 3 indices are merged lexicographically for 24 entities
      return false;
    }
  }

  return true;
};

static std::size_t test_case = 0;
template <class GFS> void check_ordering(const GFS &gfs, bool debug = false) {
  std::cout << "Test Case : " << test_case << std::endl << std::endl;
  std::cout << "GFS: " << Dune::className<GFS>() << std::endl << std::endl;
  Dune::PDELab::writeTreeToGraphviz(std::cout, gfs);
  std::cout << std::endl;

  std::cout << "Ordering: " << Dune::className<typename GFS::Ordering>()
            << std::endl
            << std::endl;
  Dune::PDELab::writeTreeToGraphviz(std::cout, gfs.ordering());
  std::cout << std::endl;

  using V = Dune::PDELab::Backend::Vector<GFS, double>;
  std::cout << "Vector type: " << Dune::className<typename V::Container>()
            << std::endl
            << std::endl;

  using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
  LFS lfs{gfs};
  Dune::PDELab::LFSIndexCache<LFS> lfs_cache{lfs};
  Dune::PDELab::EntityIndexCache<GFS,true> ei_cache{gfs};

  using DOFIndex = typename GFS::Ordering::Traits::DOFIndex;
  using ContainerIndex = typename GFS::Ordering::Traits::ContainerIndex;
  using SizeType = typename GFS::Ordering::Traits::SizeType;

  // for each suffix we store:
  //  * max seen sub-block index
  //  * a counter for the suffix
  // after accessing every container index, these should match size info from
  // ordering
  std::unordered_map<ContainerIndex, std::set<SizeType>> gfs_suffix;

  std::unordered_map<DOFIndex, ContainerIndex> lfs_map;
  std::unordered_map<DOFIndex, ContainerIndex> ec_map;

  for (const auto &element : elements(gfs.entitySet())) {
    // lfs computes dof indices
    lfs.bind(element);
    // cache computes container indices
    lfs_cache.update();

    lfs_map.clear();

    // loop over local DOFs in the lfs.
    for (unsigned ldof = 0; ldof < lfs.size(); ++ldof) {
      // get global DOF from local DOF index
      auto gdof = lfs.dofIndex(ldof);
      // map global DOF index to a container index (from ordering)
      auto ci = gfs.ordering().mapIndex(gdof);
      // map local DOF index to a container index (from cache)
      auto ci_cache_ldof = lfs_cache.containerIndex(ldof);
      // map dof index to a container index (from cache)
      auto ci_cache_gdof = lfs_cache.containerIndex(gdof);
      // the different mappings shall give the same result
      if (ci != ci_cache_ldof or ci != ci_cache_gdof)
        DUNE_THROW(Dune::RangeError, "Container index mappings do not match:"
          << ci << " " << ci_cache_ldof << " " << ci_cache_gdof);

      if (debug)
        std::cout << "DOFIndex: " << gdof << "   <-->    "
                  << "ContainerIndex: " << ci << std::endl;

      if (lfs_map.find(gdof) != lfs_map.end())
        DUNE_THROW(Dune::RangeError,
          "Local function spaces shall contain every DOF index only once!");

      lfs_map[gdof] = ci;

      // loop all possible suffixes for the current container index
      auto suffix = ci;
      auto ci_begin = ci.begin();

      while (suffix.size() != 0) {
        // get outer container index in the suffix
        auto block_index = suffix.front();
        // reduce current suffix by one and copy contents from CI
        suffix.resize(suffix.size() - 1);
        std::copy_backward(++ci_begin, ci.end(), suffix.end());
        // update max and count checks
        gfs_suffix[suffix].insert(block_index);
        if (debug) {
          std::cout << "  Suffix: " << suffix << std::endl;
          std::cout << "    Current size: " << gfs_suffix[suffix].size()
                    << std::endl;
          std::cout << "    Current range: [" << *gfs_suffix[suffix].begin()
                    << ", " << *gfs_suffix[suffix].rbegin() << "]" << std::endl;
          std::cout << "    Size: " << gfs.ordering().size(suffix) << std::endl;
        }
      }
    }

    // loop all codimensions
    auto codims = std::make_index_sequence<GFS::Traits::EntitySet::dimension+1>{};
    Dune::Hybrid::forEach(codims,[&](auto codim){
    // loop over all sub-entities...
      for (std::size_t se = 0; se < element.subEntities(codim); ++se) {
        const auto& sub_element = element.template subEntity<codim>(se);
        ei_cache.update(sub_element);
        for (unsigned ldof = 0; ldof < ei_cache.size(); ++ldof) {
          const auto& gdof = ei_cache.dofIndex(ldof);
          if (lfs_map.find(gdof) == lfs_map.end())
            DUNE_THROW(Dune::RangeError,
              "Ordering does not map the same DOFs as the EntityCache!");
          const auto& ci = ei_cache.containerIndex(ldof);
          if (lfs_map[gdof] != ci)
            DUNE_THROW(Dune::RangeError,
              "Ordering does not map the same Container Index as the EntityCache!" << std::endl
              << "DOFIndex:       " << gdof << std::endl
              << "EntityCache CI: " << ci << std::endl
              << "LFS CI:         " << lfs_map[gdof]);
        }
      }
    });
  }

  [[maybe_unused]] bool failed = false;
  for (const auto &[suffix, set] : gfs_suffix) {
    bool dense_pattern = true; // whether suffix sampled a dense pattern
    assert(not set.empty());
    auto front = *set.begin();
    auto back = *set.rbegin();
    std::cout << "Suffix: " << suffix << "   --->   Range: [" << front << ", "
              << (back + 1) << ")" << std::endl;
    if (((back + 1) != set.size()) or (front != 0)) {
      dense_pattern = false;
      std::cout << "  CI suffix '" << suffix
                << "' has no dense pattern: " << std::endl;
      for (auto i : set)
        std::cout << "    * " << i << std::endl;
    }
    auto size = gfs.ordering().size(suffix);
    if ((back + 1) != size) {
      dense_pattern = false;
      std::cout << "  CI suffix has size '" << size
                << "' but it does not match size of sampled pattern '"
                << set.size() << "'" << std::endl;
    }

    if (not dense_pattern) {
      if (is_ordering_dense(test_case, suffix)) {
        failed |= true;
      } else {
        // if ordering is not dense, we at least check that we are not out of bounds
        failed |= (back + 1) > size;
      }
    }
  }

  if (failed)
    DUNE_THROW(Dune::RangeError,
               "CI vs CI Sufix size information does not match");
  using V = Dune::PDELab::Backend::Vector<GFS,double>;
  V x(gfs); // dune-functions needs allow dynamic vectors to be 0-sized
  x = 0.0; // dune-common DenseVector needs to fix this overload on field_type!
  std::cout << "**************************************" << std::endl;
  std::cout << std::endl;
  ++test_case;
}

template <bool t0 = true, bool t1 = true, bool t2 = true, bool t3 = true,
          class UGFS>
void test_power_ordering(std::shared_ptr<UGFS> ugfs) {
  using namespace Dune::PDELab;

  check_ordering(OrderedGridFunctionSpace{ugfs});

  if constexpr (t0) {
    auto pugfs = std::make_shared<UnorderedPowerGridFunctionSpace<
        UGFS, 3, ISTL::BackendOptions<>, LexicographicOrderingTag>>(ugfs, ugfs,
                                                                    ugfs);
    check_ordering(OrderedGridFunctionSpace{pugfs});
  }

  if constexpr (t1) {
    auto pugfs = std::make_shared<UnorderedPowerGridFunctionSpace<
        UGFS, 3, ISTL::BackendOptions<ISTL::Blocking::Static>,
        LexicographicOrderingTag>>(ugfs, ugfs, ugfs);
    check_ordering(OrderedGridFunctionSpace{pugfs});
  }

  if constexpr (t2) {
    auto pugfs = std::make_shared<UnorderedPowerGridFunctionSpace<
        UGFS, 3, ISTL::BackendOptions<>, EntityBlockedOrderingTag>>(ugfs, ugfs,
                                                                    ugfs);
    check_ordering(OrderedGridFunctionSpace{pugfs});
  }

  if constexpr (t3) {
    auto pugfs = std::make_shared<UnorderedPowerGridFunctionSpace<
        UGFS, 3, ISTL::BackendOptions<ISTL::Blocking::Static>,
        EntityBlockedOrderingTag>>(ugfs, ugfs, ugfs);
    check_ordering(OrderedGridFunctionSpace{pugfs});
  }
};

template <bool blockedLeaf = true, class UGFS1, class UGFS2>
void test_power_composite_ordering(std::shared_ptr<UGFS1> ugfs1,
                                   std::shared_ptr<UGFS2> ugfs2) {
  using namespace Dune::PDELab;

  test_power_ordering<true, blockedLeaf, true, blockedLeaf>(ugfs1);
  if constexpr (not std::is_same<UGFS1, UGFS2>{})
    test_power_ordering<true, blockedLeaf, true, blockedLeaf>(ugfs2);
  else if (ugfs1 != ugfs2)
    test_power_ordering<true, blockedLeaf, true, blockedLeaf>(ugfs2);

  {
    auto cugfs = std::make_shared<UnorderedCompositeGridFunctionSpace<
        ISTL::BackendOptions<>, LexicographicOrderingTag, UGFS1, UGFS2>>(ugfs1,
                                                                         ugfs2);
    test_power_ordering<true, blockedLeaf, false, false>(cugfs);
  }

  if constexpr (blockedLeaf) {
    using CUGFS = UnorderedCompositeGridFunctionSpace<
        ISTL::BackendOptions<ISTL::Blocking::Static>, LexicographicOrderingTag,
        UGFS1, UGFS2>;
    auto cugfs = std::make_shared<CUGFS>(ugfs1, ugfs2);
    test_power_ordering<true, true, false, false>(cugfs);
  }

  {
    auto cugfs = std::make_shared<UnorderedCompositeGridFunctionSpace<
        ISTL::BackendOptions<>, EntityBlockedOrderingTag, UGFS1, UGFS2>>(ugfs1,
                                                                         ugfs2);
    test_power_ordering<true, blockedLeaf, true, blockedLeaf>(cugfs);
  }

  if constexpr (blockedLeaf) {
    auto cugfs = std::make_shared<UnorderedCompositeGridFunctionSpace<
        ISTL::BackendOptions<ISTL::Blocking::Static>, EntityBlockedOrderingTag,
        UGFS1, UGFS2>>(ugfs1, ugfs2);
    test_power_ordering<true, true, true, true>(cugfs);
  }
};

int main(int argc, char **argv) {
  // Maybe initialize Mpi
  Dune::MPIHelper::instance(argc, argv);

  using namespace Dune::PDELab;
  using Grid = YaspUnitSquare;
  using GV = typename Grid::LeafGridView;
  Grid grid;
  grid.globalRefine(1);
  auto gv = grid.leafGridView();

  using ES = AllEntitySet<typename Grid::LeafGridView>;
  ES es{gv};

  using Q2FEM = Dune::PDELab::QkLocalFiniteElementMap<GV, float, double, 2>;
  Q2FEM q2fem(gv);

  auto mgfs = std::make_shared<UnorderedGridFunctionSpace<
      ES, Q2FEM, NoConstraints, ISTL::BackendOptions<>>>(es, q2fem);
  test_power_composite_ordering<false>(mgfs, mgfs);

  auto bgfs = std::make_shared<UnorderedGridFunctionSpace<
      ES, Q2FEM, NoConstraints, ISTL::BackendOptions<ISTL::Blocking::Dynamic>>>(
      es, q2fem);
  test_power_composite_ordering<true>(bgfs, bgfs);

  typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> CellMapper;
  CellMapper cellmapper(gv);
  typedef Dune::PDELab::VariableMonomLocalFiniteElementMap<
      CellMapper, float, double, GV::dimension>
      MonomFEM;
  MonomFEM monomfem(cellmapper, Dune::GeometryTypes::quadrilateral, 2);

  // test _fixed_size_possible branch
  {
    auto monomgfs = std::make_shared<UnorderedGridFunctionSpace<
        ES, MonomFEM, NoConstraints, ISTL::BackendOptions<>>>(es, monomfem);
    test_power_composite_ordering<false>(monomgfs, monomgfs);

    auto bmonomgfs = std::make_shared<UnorderedGridFunctionSpace<
        ES, MonomFEM, NoConstraints,
        ISTL::BackendOptions<ISTL::Blocking::Dynamic>>>(es, monomfem);
    test_power_composite_ordering<true>(bmonomgfs, bmonomgfs);
  }

  // force _fixed_size to be false
  std::size_t p = 0;
  for (const auto &e : elements(gv))
    monomfem.setOrder(e, p++ % 2);

  for (auto i : {0, 1}) {
    auto &fe = monomfem.getFEM(i);
    std::cout << "Degree: " << i << std::endl;
    for (std::size_t i = 0; i < fe.size(); ++i)
      std::cout << fe.localCoefficients().localKey(i) << std::endl;
  }

  {
    auto monomgfs = std::make_shared<UnorderedGridFunctionSpace<
        ES, MonomFEM, NoConstraints, ISTL::BackendOptions<>>>(es, monomfem);
    test_power_composite_ordering<false>(monomgfs, monomgfs);

    auto bmonomgfs = std::make_shared<UnorderedGridFunctionSpace<
        ES, MonomFEM, NoConstraints,
        ISTL::BackendOptions<ISTL::Blocking::Dynamic>>>(es, monomfem);
    test_power_composite_ordering<true>(bmonomgfs, bmonomgfs);
  }

  {
    auto bmonomgfs = std::make_shared<UnorderedGridFunctionSpace<
        ES, MonomFEM, NoConstraints,
        ISTL::BackendOptions<ISTL::Blocking::Dynamic>>>(es, monomfem);
    auto bgfs = std::make_shared<UnorderedGridFunctionSpace<
        ES, Q2FEM, NoConstraints,
        ISTL::BackendOptions<ISTL::Blocking::Dynamic>>>(es, q2fem);
    test_power_composite_ordering<true>(bmonomgfs, bgfs);
  }

  // TODO: use Q1 on composite node for force empty blocks
  // TODO: use empty fems to force empty orderings

  // test passed
  return 0;
}
