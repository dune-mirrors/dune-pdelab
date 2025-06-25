#include "../gridexamples.hh"

#include <dune/pdelab/basis/protobasis.hh>

#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/variablemonomfem.hh>

#include <dune/pdelab/constraints/noconstraints.hh>

#include <dune/typetree/treepath.hh>

#include <dune/grid/common/scsgmapper.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

template<class Q1TopologicalInterleaving>
void testQ1TopologicalInterleaving(Dune::TestSuite& testSuite, Q1TopologicalInterleaving& entity_ordering) {
  static_assert(Q1TopologicalInterleaving::fixedSizePerGeometryTypeStatic());
  testSuite.check(entity_ordering.fixedSizePerGeometryType());

  static_assert(Q1TopologicalInterleaving::maxMultiIndexSize() == 1);
  testSuite.check(entity_ordering.containsGeometryType(Dune::GeometryTypes::vertex));
  testSuite.check(not entity_ordering.containsGeometryType(Dune::GeometryTypes::line));
  testSuite.check(not entity_ordering.containsGeometryType(Dune::GeometryTypes::quadrilateral));
  testSuite.check(not entity_ordering.containsGeometryType(Dune::GeometryTypes::triangle));
  testSuite.check(not entity_ordering.containsCodim(0));
  testSuite.check(not entity_ordering.containsCodim(1));
  testSuite.check(entity_ordering.containsCodim(2));
  testSuite.check(entity_ordering.singleCodim());
  testSuite.check(not entity_ordering.disjointCodimClosure());
  testSuite.check(entity_ordering.maxLocalCount() == 4);

  auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::vertex);
  testSuite.check(entity_ordering.localTreeIndices(gt_index, 0).size() == 1);

  for(const auto& vertex : vertices(entity_ordering.gridView())) {
    auto entity_index = entity_ordering.gridView().indexSet().index(vertex);
    testSuite.check(entity_ordering.localTreeIndices(gt_index, entity_index).size() == 1);
    auto tpath = Dune::TypeTree::treePath();
    auto cir = entity_ordering.localTreeIndices(gt_index, entity_index, tpath);
    testSuite.check(cir[0].size() == 1);
    testSuite.check(cir[0][0] == 0);
    testSuite.check(entity_ordering.localTreeDegree(gt_index, entity_index, cir[0]) == 0);
  }

  testSuite.check(entity_ordering.dimension() == entity_ordering.gridView().indexSet().size(2));
}

void testQ1TopologicalInterleaving() {
  Dune::TestSuite testSuite("Q1 Topologic Associativity Forest");
  YaspUnitSquare grid;
  grid.globalRefine(3);
  using GridView = typename YaspUnitSquare::LeafGridView;
  GridView grid_view{grid.leafGridView()};
  {
    constexpr auto blocked = false;
    using QkFEM = Dune::PDELab::QkLocalFiniteElementMap<GridView, double, double, 1>;
    using ProtoBasis = Dune::PDELab::ProtoBasisLeaf<GridView, QkFEM, blocked>;
    ProtoBasis q1_proto_basis(grid_view, QkFEM{grid_view});
    ProtoBasis& entity_ordering = q1_proto_basis;
    entity_ordering.initializeIndices();
    testQ1TopologicalInterleaving(testSuite, entity_ordering);

    using ProtoBasisVector = Dune::PDELab::ProtoBasisVector<ProtoBasis, blocked>;
    ProtoBasisVector p_entity_ordering{std::vector{entity_ordering,entity_ordering}};
    p_entity_ordering.initializeIndices();
    testQ1TopologicalInterleaving(testSuite, p_entity_ordering.child(0));

    testSuite.check(p_entity_ordering.degree() == 2);
    static_assert(not ProtoBasisVector::isIndexBlocked);
    static_assert(ProtoBasisVector::fixedSizePerGeometryTypeStatic());
    testSuite.check(p_entity_ordering.fixedSizePerGeometryType());

    static_assert(ProtoBasisVector::maxMultiIndexSize() == 1);
    testSuite.check(p_entity_ordering.containsGeometryType(Dune::GeometryTypes::vertex));
    testSuite.check(not p_entity_ordering.containsGeometryType(Dune::GeometryTypes::line));
    testSuite.check(not p_entity_ordering.containsGeometryType(Dune::GeometryTypes::quadrilateral));
    testSuite.check(not p_entity_ordering.containsGeometryType(Dune::GeometryTypes::triangle));
    testSuite.check(not p_entity_ordering.containsCodim(0));
    testSuite.check(not p_entity_ordering.containsCodim(1));
    testSuite.check(p_entity_ordering.containsCodim(2));
    testSuite.check(p_entity_ordering.singleCodim());
    testSuite.check(not p_entity_ordering.disjointCodimClosure());
    testSuite.check(p_entity_ordering.maxLocalCount() == 2*4);

    auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::vertex);
    testSuite.check(p_entity_ordering.localTreeIndices(gt_index, 0).size() == 2);

    for(const auto& vertex : vertices(p_entity_ordering.gridView())) {
      auto entity_index = p_entity_ordering.gridView().indexSet().index(vertex);
      testSuite.check(p_entity_ordering.localTreeIndices(gt_index, entity_index).size() == 2);
      auto tpath0 = Dune::TypeTree::treePath(0);

      auto cir0 = p_entity_ordering.localTreeIndices(gt_index, entity_index, tpath0);
      testSuite.check(cir0[0].size() == 1);
      testSuite.check(cir0[0][0] == 0);
      testSuite.check(p_entity_ordering.localTreeDegree(gt_index, entity_index, cir0[0]) == 0);

      auto tpath1 = Dune::TypeTree::treePath(1);
      auto cir1 = p_entity_ordering.localTreeIndices(gt_index, entity_index, tpath1);
      testSuite.check(cir1[0].size() == 1);
      testSuite.check(cir1[0][0] == 1);
      testSuite.check(p_entity_ordering.localTreeDegree(gt_index, entity_index, cir1[0]) == 0);
    }

    testSuite.check(p_entity_ordering.dimension() == 2*p_entity_ordering.gridView().indexSet().size(2));
  }

  {
    constexpr auto blocked = true;
    using QkFEM = Dune::PDELab::QkLocalFiniteElementMap<GridView, double, double, 1>;
    using ProtoBasis = Dune::PDELab::ProtoBasisLeaf<GridView, QkFEM, blocked>;
    ProtoBasis q1_proto_basis(grid_view, QkFEM{grid_view});
    ProtoBasis& entity_ordering = q1_proto_basis;
    entity_ordering.initializeIndices();
    testQ1TopologicalInterleaving(testSuite, entity_ordering);

    using ProtoBasisArray = Dune::PDELab::ProtoBasisArray<ProtoBasis, 2, blocked>;

    ProtoBasisArray p_entity_ordering{std::array{entity_ordering,entity_ordering}};
    p_entity_ordering.initializeIndices();
    testQ1TopologicalInterleaving(testSuite, p_entity_ordering.child(0));

    static_assert(ProtoBasisArray::degree() == 2);
    static_assert(ProtoBasisArray::isIndexBlocked);
    static_assert(ProtoBasisArray::fixedSizePerGeometryTypeStatic());
    testSuite.check(p_entity_ordering.fixedSizePerGeometryType());

    static_assert(ProtoBasisArray::maxMultiIndexSize() == 2);
    testSuite.check(p_entity_ordering.containsGeometryType(Dune::GeometryTypes::vertex));
    testSuite.check(not p_entity_ordering.containsGeometryType(Dune::GeometryTypes::line));
    testSuite.check(not p_entity_ordering.containsGeometryType(Dune::GeometryTypes::quadrilateral));
    testSuite.check(not p_entity_ordering.containsGeometryType(Dune::GeometryTypes::triangle));
    testSuite.check(not p_entity_ordering.containsCodim(0));
    testSuite.check(not p_entity_ordering.containsCodim(1));
    testSuite.check(p_entity_ordering.containsCodim(2));
    testSuite.check(p_entity_ordering.singleCodim());
    testSuite.check(not p_entity_ordering.disjointCodimClosure());
    testSuite.check(p_entity_ordering.maxLocalCount() == 2*4);

    auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::vertex);

    auto gt_count = p_entity_ordering.localTreeIndices(gt_index, 0).size();
    testSuite.check(gt_count == 2);
    static_assert(std::same_as<decltype(gt_count), std::integral_constant<std::size_t,2>>);

    for(const auto& vertex : vertices(p_entity_ordering.gridView())) {
      auto entity_index = p_entity_ordering.gridView().indexSet().index(vertex);
      testSuite.check(Dune::Hybrid::equal_to(p_entity_ordering.localTreeIndices(gt_index, entity_index).size(), Dune::Indices::_2));

      using namespace Dune::Indices;
      auto tpath0 = Dune::TypeTree::treePath(_0);
      auto cir0 = p_entity_ordering.localTreeIndices(gt_index, entity_index, tpath0);
      testSuite.check(cir0.size() == 1);
      testSuite.check(cir0[0].size() == 2);
      testSuite.check(cir0[0][0] == 0);
      testSuite.check(cir0[0][1] == 0);
      testSuite.check(p_entity_ordering.localTreeDegree(gt_index, entity_index, reverse(cir0[0])) == 0);
      testSuite.check(p_entity_ordering.localTreeDegree(gt_index, entity_index, pop_front(reverse(cir0[0]))) == 1);

      auto tpath1 = Dune::TypeTree::treePath(1);
      auto cir1 = p_entity_ordering.localTreeIndices(gt_index, entity_index, tpath1);
      testSuite.check(cir1.size() == 1);
      testSuite.check(cir1[0].size() == 2);
      testSuite.check(cir1[0][0] == 1);
      testSuite.check(cir1[0][1] == 0);
      testSuite.check(p_entity_ordering.localTreeDegree(gt_index, entity_index, reverse(cir1[0])) == 0);
      testSuite.check(p_entity_ordering.localTreeDegree(gt_index, entity_index, pop_front(reverse(cir1[0]))) == 1);

      testSuite.check(p_entity_ordering.localTreeDegree(gt_index, entity_index) == 2);
    }

    testSuite.check(p_entity_ordering.dimension() == 2*p_entity_ordering.gridView().indexSet().size(2));
  }
}

template<class Q2TopologicalInterleaving>
void testQ2TopologicalInterleaving(Dune::TestSuite& testSuite, Q2TopologicalInterleaving& entity_ordering) {
  static_assert(not Q2TopologicalInterleaving::isIndexBlocked);
  static_assert(Q2TopologicalInterleaving::fixedSizePerGeometryTypeStatic());
  testSuite.check(entity_ordering.fixedSizePerGeometryType());

  static_assert(Q2TopologicalInterleaving::maxMultiIndexSize() == 1);
  testSuite.check(entity_ordering.containsGeometryType(Dune::GeometryTypes::vertex));
  testSuite.check(entity_ordering.containsGeometryType(Dune::GeometryTypes::line));
  testSuite.check(entity_ordering.containsGeometryType(Dune::GeometryTypes::quadrilateral));
  testSuite.check(not entity_ordering.containsGeometryType(Dune::GeometryTypes::triangle));
  testSuite.check(entity_ordering.containsCodim(0));
  testSuite.check(entity_ordering.containsCodim(1));
  testSuite.check(entity_ordering.containsCodim(2));
  testSuite.check(not entity_ordering.singleCodim());
  testSuite.check(not entity_ordering.disjointCodimClosure());
  testSuite.check(entity_ordering.maxLocalCount() == 9);

  auto gt_index_vertex = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::vertex);
  testSuite.check(entity_ordering.localTreeIndices(gt_index_vertex,0).size() == 1);

  for(const auto& vertex : vertices(entity_ordering.gridView())) {
    auto entity_index = entity_ordering.gridView().indexSet().index(vertex);
    testSuite.check(entity_ordering.localTreeIndices(gt_index_vertex, entity_index).size() == 1);
    auto tpath = Dune::TypeTree::treePath();
    auto cir = entity_ordering.localTreeIndices(gt_index_vertex, entity_index, tpath);
    testSuite.check(cir.size() == 1);
    testSuite.check(cir[0].size() == 1);
    testSuite.check(cir[0][0] == 0);
    testSuite.check(entity_ordering.localTreeDegree(gt_index_vertex, entity_index, cir[0]) == 0);
  }

  auto gt_index_edge = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::line);
  testSuite.check(entity_ordering.localTreeIndices(gt_index_edge, 0).size() == 1);

  for(const auto& edge : edges(entity_ordering.gridView())) {
    auto entity_index = entity_ordering.gridView().indexSet().index(edge);
    testSuite.check(entity_ordering.localTreeIndices(gt_index_edge, entity_index).size() == 1);
    auto tpath = Dune::TypeTree::treePath();
    auto cir = entity_ordering.localTreeIndices(gt_index_edge, entity_index, tpath);
    testSuite.check(cir.size() == 1);
    testSuite.check(cir[0].size() == 1);
    testSuite.check(cir[0][0] == 0);
    testSuite.check(entity_ordering.localTreeDegree(gt_index_edge, entity_index, cir[0]) == 0);
  }

  auto gt_index_quad = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::quadrilateral);
  testSuite.check(entity_ordering.localTreeIndices(gt_index_quad, 0).size() == 1);

  for(const auto& edge : edges(entity_ordering.gridView())) {
    auto entity_index = entity_ordering.gridView().indexSet().index(edge);
    testSuite.check(entity_ordering.localTreeIndices(gt_index_quad, entity_index).size() == 1);
    auto tpath = Dune::TypeTree::treePath();
    auto cir = entity_ordering.localTreeIndices(gt_index_edge, entity_index, tpath);
    testSuite.check(cir.size() == 1);
    testSuite.check(cir[0].size() == 1);
    testSuite.check(cir[0][0] == 0);
    testSuite.check(entity_ordering.localTreeDegree(gt_index_edge, entity_index, cir[0]) == 0);
  }

  testSuite.check(entity_ordering.dimension() ==
    entity_ordering.gridView().indexSet().size(2) +
    entity_ordering.gridView().indexSet().size(1) +
    entity_ordering.gridView().indexSet().size(0)
  );
}

void testQ2TopologicalInterleaving() {
  Dune::TestSuite testSuite("Q2 Topologic Associativity Forest");

  YaspUnitSquare grid;
  grid.globalRefine(3);
  using GridView = typename YaspUnitSquare::LeafGridView;
  GridView grid_view{grid.leafGridView()};

  constexpr auto blocked = false;
  using QkFEM = Dune::PDELab::QkLocalFiniteElementMap<GridView, double, double, 2>;
  using ProtoBasis = Dune::PDELab::ProtoBasisLeaf<GridView, QkFEM, blocked>;
  ProtoBasis q2_pb(grid_view, QkFEM{grid_view});

  q2_pb.initializeIndices();
  testQ2TopologicalInterleaving(testSuite, q2_pb);
}


template<class FixedMonomialTopologicalInterleaving>
void testFixedMonomialTopologicalInterleaving(Dune::TestSuite& testSuite, FixedMonomialTopologicalInterleaving& entity_ordering) {
  static_assert(not FixedMonomialTopologicalInterleaving::isIndexBlocked);
  static_assert(not FixedMonomialTopologicalInterleaving::fixedSizePerGeometryTypeStatic());
  testSuite.check(entity_ordering.fixedSizePerGeometryType());

  static_assert(FixedMonomialTopologicalInterleaving::maxMultiIndexSize() == 1);
  testSuite.check(not entity_ordering.containsGeometryType(Dune::GeometryTypes::vertex));
  testSuite.check(not entity_ordering.containsGeometryType(Dune::GeometryTypes::line));
  testSuite.check(entity_ordering.containsGeometryType(Dune::GeometryTypes::quadrilateral));
  testSuite.check(not entity_ordering.containsGeometryType(Dune::GeometryTypes::triangle));
  testSuite.check(entity_ordering.containsCodim(0));
  testSuite.check(not entity_ordering.containsCodim(1));
  testSuite.check(not entity_ordering.containsCodim(2));
  testSuite.check(entity_ordering.singleCodim());
  testSuite.check(entity_ordering.disjointCodimClosure());
  testSuite.check(entity_ordering.maxLocalCount() == 6);

  auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::quadrilateral);
  testSuite.check(entity_ordering.localTreeIndices(gt_index, 0).size() == 6);

  for(const auto& element : elements(entity_ordering.gridView())) {
    auto entity_index = entity_ordering.gridView().indexSet().index(element);
    testSuite.check(entity_ordering.localTreeIndices(gt_index, entity_index).size() == 6);
    auto tpath = Dune::TypeTree::treePath();
    auto cir = entity_ordering.localTreeIndices(gt_index, entity_index, tpath);
    testSuite.check(cir.size() == 6);
    testSuite.check(cir[0].size() == 1);
    testSuite.check(cir[0][0] == 0);
    testSuite.check(entity_ordering.localTreeDegree(gt_index, entity_index, cir[0]) == 0);
  }

  testSuite.check(entity_ordering.dimension() ==
    entity_ordering.gridView().indexSet().size(0) * entity_ordering.localTreeIndices(gt_index, 0).size()
  );
}


void testFixedMonomialTopologicalInterleaving() {
  Dune::TestSuite testSuite("Fixed Size Monomial Topologic Associativity Forest");

  YaspUnitSquare grid;
  grid.globalRefine(3);
  using GridView = typename YaspUnitSquare::LeafGridView;
  GridView grid_view{grid.leafGridView()};

  using CellMapper = Dune::SingleCodimSingleGeomTypeMapper<GridView, 0> ;
  CellMapper cellmapper(grid_view);
  using MonomFEM = Dune::PDELab::VariableMonomLocalFiniteElementMap<CellMapper, float, double, GridView::dimension>;
  auto monomfem = MonomFEM(cellmapper, Dune::GeometryTypes::quadrilateral, 2);

  constexpr bool merged = false;
  using ProtoBasis = Dune::PDELab::ProtoBasisLeaf<GridView, MonomFEM, merged>;
  ProtoBasis monom_pb(grid_view, monomfem);
  auto& entity_ordering = monom_pb;

  entity_ordering.initializeIndices();

  testFixedMonomialTopologicalInterleaving(testSuite, entity_ordering);

  using ProtoBasisArray = Dune::PDELab::ProtoBasisArray<ProtoBasis, 2, merged>;

  ProtoBasisArray p_entity_ordering{std::array{entity_ordering,entity_ordering}};
  p_entity_ordering.initializeIndices();
  testFixedMonomialTopologicalInterleaving(testSuite, p_entity_ordering.child(0));

  static_assert(not ProtoBasisArray::isIndexBlocked);
  static_assert(not p_entity_ordering.fixedSizePerGeometryTypeStatic());
  testSuite.check(p_entity_ordering.fixedSizePerGeometryType());

  static_assert(p_entity_ordering.maxMultiIndexSize() == 1);
  testSuite.check(not p_entity_ordering.containsGeometryType(Dune::GeometryTypes::vertex));
  testSuite.check(not p_entity_ordering.containsGeometryType(Dune::GeometryTypes::line));
  testSuite.check(p_entity_ordering.containsGeometryType(Dune::GeometryTypes::quadrilateral));
  testSuite.check(not p_entity_ordering.containsGeometryType(Dune::GeometryTypes::triangle));
  testSuite.check(p_entity_ordering.containsCodim(0));
  testSuite.check(not p_entity_ordering.containsCodim(1));
  testSuite.check(not p_entity_ordering.containsCodim(2));
  testSuite.check(p_entity_ordering.singleCodim());
  testSuite.check(p_entity_ordering.disjointCodimClosure());
  testSuite.check(p_entity_ordering.maxLocalCount() == 2*6);

  auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::quadrilateral);
  testSuite.check(p_entity_ordering.localTreeIndices(gt_index, 0).size() == 2*6);

  for(const auto& element : elements(p_entity_ordering.gridView())) {
    auto entity_index = p_entity_ordering.gridView().indexSet().index(element);
    testSuite.check(p_entity_ordering.localTreeIndices(gt_index, entity_index).size() == 2*6);

    auto tpath0 = Dune::TypeTree::treePath(0);
    auto cir0 = p_entity_ordering.localTreeIndices(gt_index, entity_index, tpath0);
    testSuite.check(cir0.size() == 6);
    testSuite.check(cir0[0].size() == 1);
    testSuite.check(cir0[0][0] == 0);
    testSuite.check(p_entity_ordering.localTreeDegree(gt_index, entity_index, cir0[0]) == 0);

    auto tpath1 = Dune::TypeTree::treePath(1);
    auto cir1 = p_entity_ordering.localTreeIndices(gt_index, entity_index, tpath1);
    testSuite.check(cir1.size() == 6);
    testSuite.check(cir1[0].size() == 1);
    testSuite.check(cir1[0][0] == 6);
    testSuite.check(p_entity_ordering.localTreeDegree(gt_index, entity_index, cir1[0]) == 0);
  }

  testSuite.check(p_entity_ordering.dimension() ==
    p_entity_ordering.gridView().indexSet().size(0) * p_entity_ordering.localTreeIndices(gt_index, 0).size()
  );
}

template<class VariableMonomialTopologicalInterleaving>
void testVariableMonomialTopologicalInterleaving(Dune::TestSuite& testSuite, VariableMonomialTopologicalInterleaving& entity_ordering) {
  static_assert(not VariableMonomialTopologicalInterleaving::isIndexBlocked);
  static_assert(not VariableMonomialTopologicalInterleaving::fixedSizePerGeometryTypeStatic());
  testSuite.check(not entity_ordering.fixedSizePerGeometryType());

  static_assert(VariableMonomialTopologicalInterleaving::maxMultiIndexSize() == 1);
  testSuite.check(not entity_ordering.containsGeometryType(Dune::GeometryTypes::vertex));
  testSuite.check(not entity_ordering.containsGeometryType(Dune::GeometryTypes::line));
  testSuite.check(entity_ordering.containsGeometryType(Dune::GeometryTypes::quadrilateral));
  testSuite.check(not entity_ordering.containsGeometryType(Dune::GeometryTypes::triangle));
  testSuite.check(entity_ordering.containsCodim(0));
  testSuite.check(not entity_ordering.containsCodim(1));
  testSuite.check(not entity_ordering.containsCodim(2));
  testSuite.check(entity_ordering.singleCodim());
  testSuite.check(entity_ordering.disjointCodimClosure());
  testSuite.check(entity_ordering.maxLocalCount() == 6);

  auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::quadrilateral);

  std::size_t order = 0;
  for(const auto& element : elements(entity_ordering.gridView())) {
    auto& fe = entity_ordering.finiteElementMap().getFEM(order++ % 3);
    auto entity_index = entity_ordering.gridView().indexSet().index(element);
    testSuite.check(entity_ordering.localTreeIndices(gt_index, entity_index).size() == fe.size());
    auto tpath = Dune::TypeTree::treePath();
    auto cir = entity_ordering.localTreeIndices(gt_index, entity_index, tpath);
    testSuite.check(cir.size() == fe.size());
    testSuite.check(cir[0].size() == 1);
    testSuite.check(cir[0][0] == 0);
    testSuite.check(entity_ordering.localTreeDegree(gt_index, entity_index, cir[0]) == 0);
  }

  testSuite.check(entity_ordering.gridView().indexSet().size(0) == 4);

  testSuite.check(entity_ordering.dimension() ==
    1 + 3 + 6 + 1
  );
}


void testVariableMonomialTopologicalInterleaving() {
  Dune::TestSuite testSuite("Fixed Size Monomial Topologic Associativity Forest");

  Dune::FieldVector<double,2> L(1.0);
  std::array<int,2> N;
  std::fill(begin(N), end(N), 2); // for this test we require exactly 4 entities
  Dune::YaspGrid<2> grid{L,N};
  using GridView = typename Dune::YaspGrid<2>::LeafGridView;
  GridView grid_view{grid.leafGridView()};

  using CellMapper = Dune::SingleCodimSingleGeomTypeMapper<GridView, 0> ;
  CellMapper cellmapper(grid_view);
  using MonomFEM = Dune::PDELab::VariableMonomLocalFiniteElementMap<CellMapper, float, double, GridView::dimension>;
  auto monomfem = MonomFEM(cellmapper, Dune::GeometryTypes::quadrilateral, 2);

  std::size_t order = 0;
  for (const auto &e : elements(grid_view))
    monomfem.setOrder(e, order++ % 3);

  constexpr bool merged = false;
  using ProtoBasis = Dune::PDELab::ProtoBasisLeaf<GridView, MonomFEM, merged>;
  ProtoBasis monom_pb(grid_view, monomfem);
  auto& entity_ordering = monom_pb;

  entity_ordering.initializeIndices();

  testVariableMonomialTopologicalInterleaving(testSuite, entity_ordering);

  using namespace Dune::Indices;
  using ProtoBasisTuple = Dune::PDELab::ProtoBasisTuple<merged, ProtoBasis, ProtoBasis>;

  ProtoBasisTuple p_entity_ordering{std::tuple{entity_ordering,entity_ordering}};
  p_entity_ordering.initializeIndices();
  testVariableMonomialTopologicalInterleaving(testSuite, p_entity_ordering.child(_0));
  testVariableMonomialTopologicalInterleaving(testSuite, p_entity_ordering.child(_1));

  static_assert(not ProtoBasisTuple::isIndexBlocked);
  static_assert(not ProtoBasisTuple::fixedSizePerGeometryTypeStatic());
  testSuite.check(not p_entity_ordering.fixedSizePerGeometryType());

  static_assert(ProtoBasisTuple::maxMultiIndexSize() == 1);
  testSuite.check(not p_entity_ordering.containsGeometryType(Dune::GeometryTypes::vertex));
  testSuite.check(not p_entity_ordering.containsGeometryType(Dune::GeometryTypes::line));
  testSuite.check(p_entity_ordering.containsGeometryType(Dune::GeometryTypes::quadrilateral));
  testSuite.check(not p_entity_ordering.containsGeometryType(Dune::GeometryTypes::triangle));
  testSuite.check(p_entity_ordering.containsCodim(0));
  testSuite.check(not p_entity_ordering.containsCodim(1));
  testSuite.check(not p_entity_ordering.containsCodim(2));
  testSuite.check(p_entity_ordering.singleCodim());
  testSuite.check(p_entity_ordering.disjointCodimClosure());
  testSuite.check(p_entity_ordering.maxLocalCount() == 2*6);

  auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::quadrilateral);

  order = 0;
  for(const auto& element : elements(p_entity_ordering.gridView())) {
    auto& fe = p_entity_ordering.child(_0).finiteElementMap().getFEM(order++ % 3);
    auto entity_index = p_entity_ordering.gridView().indexSet().index(element);
    testSuite.check(p_entity_ordering.localTreeIndices(gt_index, entity_index).size() == 2*fe.size());

    testSuite.check(p_entity_ordering.child(_0).localTreeIndices(gt_index, entity_index).size() == fe.size());
    auto tpath0 = Dune::TypeTree::treePath(_0);
    auto cir0 = p_entity_ordering.localTreeIndices(gt_index, entity_index, tpath0);
    testSuite.check(cir0.size() == fe.size());
    testSuite.check(cir0[0].size() == 1);
    testSuite.check(cir0[0][0] == 0);
    testSuite.check(p_entity_ordering.localTreeDegree(gt_index, entity_index, cir0[0]) == 0);

    testSuite.check(p_entity_ordering.child(_1).localTreeIndices(gt_index, entity_index).size() == fe.size());
    auto tpath1 = Dune::TypeTree::treePath(_1);
    auto cir1 = p_entity_ordering.localTreeIndices(gt_index, entity_index, tpath1);
    testSuite.check(cir1.size() == fe.size());
    testSuite.check(cir1[0].size() == 1);
    testSuite.check(cir1[0][0] == fe.size());
    testSuite.check(p_entity_ordering.localTreeDegree(gt_index, entity_index, cir1[0]) == 0);
  }

  testSuite.check(p_entity_ordering.dimension() ==
    2*(1 + 3 + 6 + 1)
  );
}

int main(int argc, char **argv) {
  Dune::MPIHelper::instance(argc, argv);
  testQ1TopologicalInterleaving();
  testQ2TopologicalInterleaving();
  testFixedMonomialTopologicalInterleaving();
  testVariableMonomialTopologicalInterleaving();
}
