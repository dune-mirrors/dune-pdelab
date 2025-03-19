#include "../gridexamples.hh"

#include <dune/pdelab/basis/topologicassociativityforest/node.hh>
#include <dune/pdelab/basis/topologicassociativityforest/leaf.hh>
#include <dune/pdelab/basis/topologicassociativityforest/composite.hh>

#include <dune/pdelab/basis/merging_strategy.hh>
#include <dune/pdelab/basis/protobasis.hh>

#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/variablemonomfem.hh>

#include <dune/pdelab/constraints/noconstraints.hh>

#include <dune/typetree/treepath.hh>

#include <dune/grid/common/scsgmapper.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

template<class MergingStragetgy, std::size_t k>
auto makeQkProtoBasis(const MergingStragetgy& strategy, Dune::index_constant<k>) {
  using Q1FEM = Dune::PDELab::QkLocalFiniteElementMap<typename MergingStragetgy::EntitySet, double, double, k>;
  auto qkfem = std::make_shared<Q1FEM>(strategy.entitySet());
  return makeProtoBasis(strategy, qkfem);
}

template<class Q1TopologicAssociativityForest>
void testQ1TopologicAssociativityForest(Dune::TestSuite& testSuite, Q1TopologicAssociativityForest& entity_ordering) {
  static_assert(not Q1TopologicAssociativityForest::containerBlocked());
  static_assert(Q1TopologicAssociativityForest::fixedSizePerGeometryTypeStatic());
  testSuite.check(entity_ordering.fixedSizePerGeometryType());

  static_assert(Q1TopologicAssociativityForest::maxContainerDepth() == 1);
  testSuite.check(entity_ordering.containsGeometry(Dune::GeometryTypes::vertex));
  testSuite.check(not entity_ordering.containsGeometry(Dune::GeometryTypes::line));
  testSuite.check(not entity_ordering.containsGeometry(Dune::GeometryTypes::quadrilateral));
  testSuite.check(not entity_ordering.containsGeometry(Dune::GeometryTypes::triangle));
  testSuite.check(not entity_ordering.containsCodim(0));
  testSuite.check(not entity_ordering.containsCodim(1));
  testSuite.check(entity_ordering.containsCodim(2));
  testSuite.check(entity_ordering.singleCodim());
  testSuite.check(not entity_ordering.disjointCodimClosure());
  testSuite.check(entity_ordering.maxLocalCount() == 4);

  auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::vertex);
  testSuite.check(entity_ordering.blockCount(gt_index) == 1);

  for(const auto& vertex : vertices(entity_ordering.entitySet())) {
    auto entity_index = entity_ordering.entitySet().indexSet().index(vertex);
    testSuite.check(entity_ordering.blockCount(gt_index, entity_index) == 1);
    auto tpath = Dune::TypeTree::treePath();
    auto cir = entity_ordering.containerIndexRange(tpath, gt_index, entity_index);
    testSuite.check(cir[0].size() == 1);
    testSuite.check(cir[0][0] == 0);
    testSuite.check(entity_ordering.containerSize(cir[0], gt_index, entity_index) == 0);
  }

  testSuite.check(entity_ordering.blockCount() == entity_ordering.dimension());
  testSuite.check(entity_ordering.dimension() == entity_ordering.entitySet().indexSet().size(2));
}

void testQ1TopologicAssociativityForest() {
  Dune::TestSuite testSuite("Q1 Topologic Associativity Forest");
  YaspUnitSquare grid;
  grid.globalRefine(3);
  using GridView = typename YaspUnitSquare::LeafGridView;
  GridView entity_set{grid.leafGridView()};
  {
    auto flat_strategy = Dune::PDELab::BasisFactory::flatByEntity(entity_set);
    auto q1_proto_basis = makeQkProtoBasis(flat_strategy, Dune::Indices::_1);
    using ProtoBasis = decltype(q1_proto_basis);
    using TopologicAssociativityForest = Dune::PDELab::Impl::LeafTopologicAssociativityForest<ProtoBasis>;
    auto entity_ordering = std::make_shared<TopologicAssociativityForest>(q1_proto_basis);
    entity_ordering->update();
    testQ1TopologicAssociativityForest(testSuite, *entity_ordering);

    using VectorTopologicAssociativityForest = Dune::PDELab::Impl::VectorTopologicAssociativityForest<decltype(flat_strategy), TopologicAssociativityForest>;

    VectorTopologicAssociativityForest p_entity_ordering{std::vector{entity_ordering,entity_ordering}, flat_strategy};
    p_entity_ordering.update();
    testQ1TopologicAssociativityForest(testSuite, p_entity_ordering.child(0));
    testSuite.check(p_entity_ordering.childStorage(0) == p_entity_ordering.childStorage(1));

    testSuite.check(p_entity_ordering.degree() == 2);
    static_assert(not VectorTopologicAssociativityForest::containerBlocked());
    static_assert(VectorTopologicAssociativityForest::fixedSizePerGeometryTypeStatic());
    testSuite.check(p_entity_ordering.fixedSizePerGeometryType());

    static_assert(VectorTopologicAssociativityForest::maxContainerDepth() == 1);
    testSuite.check(p_entity_ordering.containsGeometry(Dune::GeometryTypes::vertex));
    testSuite.check(not p_entity_ordering.containsGeometry(Dune::GeometryTypes::line));
    testSuite.check(not p_entity_ordering.containsGeometry(Dune::GeometryTypes::quadrilateral));
    testSuite.check(not p_entity_ordering.containsGeometry(Dune::GeometryTypes::triangle));
    testSuite.check(not p_entity_ordering.containsCodim(0));
    testSuite.check(not p_entity_ordering.containsCodim(1));
    testSuite.check(p_entity_ordering.containsCodim(2));
    testSuite.check(p_entity_ordering.singleCodim());
    testSuite.check(not p_entity_ordering.disjointCodimClosure());
    testSuite.check(p_entity_ordering.maxLocalCount() == 2*4);

    auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::vertex);
    testSuite.check(p_entity_ordering.blockCount(gt_index) == 2);

    for(const auto& vertex : vertices(p_entity_ordering.entitySet())) {
      auto entity_index = p_entity_ordering.entitySet().indexSet().index(vertex);
      testSuite.check(p_entity_ordering.blockCount(gt_index, entity_index) == 2);
      auto tpath0 = Dune::TypeTree::treePath(0);

      auto cir0 = p_entity_ordering.containerIndexRange(tpath0, gt_index, entity_index);
      testSuite.check(cir0[0].size() == 1);
      testSuite.check(cir0[0][0] == 0);
      testSuite.check(p_entity_ordering.containerSize(cir0[0], gt_index, entity_index) == 0);

      auto tpath1 = Dune::TypeTree::treePath(1);
      auto cir1 = p_entity_ordering.containerIndexRange(tpath1, gt_index, entity_index);
      testSuite.check(cir1[0].size() == 1);
      testSuite.check(cir1[0][0] == 1);
      testSuite.check(p_entity_ordering.containerSize(cir1[0], gt_index, entity_index) == 0);
    }

    testSuite.check(p_entity_ordering.blockCount() == 2*p_entity_ordering.entitySet().indexSet().size(2));
    testSuite.check(p_entity_ordering.dimension() == p_entity_ordering.blockCount());
  }

  {
    auto blocked_strategy = Dune::PDELab::BasisFactory::blockedByEntity(entity_set);
    auto q1_proto_basis = makeQkProtoBasis(blocked_strategy, Dune::Indices::_1);
    using ProtoBasis = decltype(q1_proto_basis);
    using TopologicAssociativityForest = Dune::PDELab::Impl::LeafTopologicAssociativityForest<ProtoBasis>;
    auto entity_ordering = std::make_shared<TopologicAssociativityForest>(q1_proto_basis);
    entity_ordering->update();
    testQ1TopologicAssociativityForest(testSuite, *entity_ordering);

    using ArrayTopologicAssociativityForest = Dune::PDELab::Impl::ArrayTopologicAssociativityForest<decltype(blocked_strategy), TopologicAssociativityForest, 2>;

    ArrayTopologicAssociativityForest p_entity_ordering{std::array{entity_ordering,entity_ordering}, blocked_strategy};
    p_entity_ordering.update();
    testQ1TopologicAssociativityForest(testSuite, p_entity_ordering.child(0));

    static_assert(ArrayTopologicAssociativityForest::degree() == 2);
    static_assert(ArrayTopologicAssociativityForest::containerBlocked());
    static_assert(ArrayTopologicAssociativityForest::fixedSizePerGeometryTypeStatic());
    testSuite.check(p_entity_ordering.fixedSizePerGeometryType());

    static_assert(ArrayTopologicAssociativityForest::maxContainerDepth() == 2);
    testSuite.check(p_entity_ordering.containsGeometry(Dune::GeometryTypes::vertex));
    testSuite.check(not p_entity_ordering.containsGeometry(Dune::GeometryTypes::line));
    testSuite.check(not p_entity_ordering.containsGeometry(Dune::GeometryTypes::quadrilateral));
    testSuite.check(not p_entity_ordering.containsGeometry(Dune::GeometryTypes::triangle));
    testSuite.check(not p_entity_ordering.containsCodim(0));
    testSuite.check(not p_entity_ordering.containsCodim(1));
    testSuite.check(p_entity_ordering.containsCodim(2));
    testSuite.check(p_entity_ordering.singleCodim());
    testSuite.check(not p_entity_ordering.disjointCodimClosure());
    testSuite.check(p_entity_ordering.maxLocalCount() == 2*4);

    auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::vertex);

    auto gt_count = p_entity_ordering.blockCount(gt_index);
    testSuite.check(gt_count == 2);
    static_assert(std::same_as<decltype(gt_count), std::integral_constant<std::size_t,2>>);

    for(const auto& vertex : vertices(p_entity_ordering.entitySet())) {
      auto entity_index = p_entity_ordering.entitySet().indexSet().index(vertex);
      testSuite.check(Dune::Hybrid::equal_to(p_entity_ordering.blockCount(gt_index, entity_index), Dune::Indices::_2));

      using namespace Dune::Indices;
      auto tpath0 = Dune::TypeTree::treePath(_0);
      auto cir0 = p_entity_ordering.containerIndexRange(tpath0, gt_index, entity_index);
      testSuite.check(cir0.size() == 1);
      testSuite.check(cir0[0].size() == 2);
      testSuite.check(cir0[0][0] == 0);
      testSuite.check(cir0[0][1] == 0);
      testSuite.check(p_entity_ordering.containerSize(reverse(cir0[0]), gt_index, entity_index) == 0);
      testSuite.check(p_entity_ordering.containerSize(pop_front(reverse(cir0[0])), gt_index, entity_index) == 1);

      auto tpath1 = Dune::TypeTree::treePath(1);
      auto cir1 = p_entity_ordering.containerIndexRange(tpath1, gt_index, entity_index);
      testSuite.check(cir1.size() == 1);
      testSuite.check(cir1[0].size() == 2);
      testSuite.check(cir1[0][0] == 1);
      testSuite.check(cir1[0][1] == 0);
      testSuite.check(p_entity_ordering.containerSize(reverse(cir1[0]), gt_index, entity_index) == 0);
      testSuite.check(p_entity_ordering.containerSize(pop_front(reverse(cir1[0])), gt_index, entity_index) == 1);

      testSuite.check(p_entity_ordering.containerSize(Dune::TypeTree::treePath(), gt_index, entity_index) == 2);
    }

    auto count = p_entity_ordering.blockCount();
    testSuite.check(count == 2);
    static_assert(std::same_as<decltype(count), std::integral_constant<std::size_t,2>>);
    testSuite.check(p_entity_ordering.dimension() == 2*p_entity_ordering.entitySet().indexSet().size(2));
  }
}

template<class Q2TopologicAssociativityForest>
void testQ2TopologicAssociativityForest(Dune::TestSuite& testSuite, Q2TopologicAssociativityForest& entity_ordering) {
  static_assert(not Q2TopologicAssociativityForest::containerBlocked());
  static_assert(Q2TopologicAssociativityForest::fixedSizePerGeometryTypeStatic());
  testSuite.check(entity_ordering.fixedSizePerGeometryType());

  static_assert(Q2TopologicAssociativityForest::maxContainerDepth() == 1);
  testSuite.check(entity_ordering.containsGeometry(Dune::GeometryTypes::vertex));
  testSuite.check(entity_ordering.containsGeometry(Dune::GeometryTypes::line));
  testSuite.check(entity_ordering.containsGeometry(Dune::GeometryTypes::quadrilateral));
  testSuite.check(not entity_ordering.containsGeometry(Dune::GeometryTypes::triangle));
  testSuite.check(entity_ordering.containsCodim(0));
  testSuite.check(entity_ordering.containsCodim(1));
  testSuite.check(entity_ordering.containsCodim(2));
  testSuite.check(not entity_ordering.singleCodim());
  testSuite.check(not entity_ordering.disjointCodimClosure());
  testSuite.check(entity_ordering.maxLocalCount() == 9);

  auto gt_index_vertex = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::vertex);
  testSuite.check(entity_ordering.blockCount(gt_index_vertex) == 1);

  for(const auto& vertex : vertices(entity_ordering.entitySet())) {
    auto entity_index = entity_ordering.entitySet().indexSet().index(vertex);
    testSuite.check(entity_ordering.blockCount(gt_index_vertex, entity_index) == 1);
    auto tpath = Dune::TypeTree::treePath();
    auto cir = entity_ordering.containerIndexRange(tpath, gt_index_vertex, entity_index);
    testSuite.check(cir.size() == 1);
    testSuite.check(cir[0].size() == 1);
    testSuite.check(cir[0][0] == 0);
    testSuite.check(entity_ordering.containerSize(cir[0], gt_index_vertex, entity_index) == 0);
  }

  auto gt_index_edge = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::line);
  testSuite.check(entity_ordering.blockCount(gt_index_edge) == 1);

  for(const auto& edge : edges(entity_ordering.entitySet())) {
    auto entity_index = entity_ordering.entitySet().indexSet().index(edge);
    testSuite.check(entity_ordering.blockCount(gt_index_edge, entity_index) == 1);
    auto tpath = Dune::TypeTree::treePath();
    auto cir = entity_ordering.containerIndexRange(tpath, gt_index_edge, entity_index);
    testSuite.check(cir.size() == 1);
    testSuite.check(cir[0].size() == 1);
    testSuite.check(cir[0][0] == 0);
    testSuite.check(entity_ordering.containerSize(cir[0], gt_index_edge, entity_index) == 0);
  }

  auto gt_index_quad = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::quadrilateral);
  testSuite.check(entity_ordering.blockCount(gt_index_quad) == 1);

  for(const auto& edge : edges(entity_ordering.entitySet())) {
    auto entity_index = entity_ordering.entitySet().indexSet().index(edge);
    testSuite.check(entity_ordering.blockCount(gt_index_quad, entity_index) == 1);
    auto tpath = Dune::TypeTree::treePath();
    auto cir = entity_ordering.containerIndexRange(tpath, gt_index_edge, entity_index);
    testSuite.check(cir.size() == 1);
    testSuite.check(cir[0].size() == 1);
    testSuite.check(cir[0][0] == 0);
    testSuite.check(entity_ordering.containerSize(cir[0], gt_index_edge, entity_index) == 0);
  }

  testSuite.check(entity_ordering.dimension() ==
    entity_ordering.entitySet().indexSet().size(2) +
    entity_ordering.entitySet().indexSet().size(1) +
    entity_ordering.entitySet().indexSet().size(0)
  );
  testSuite.check(entity_ordering.dimension() == entity_ordering.blockCount());
}

void testQ2TopologicAssociativityForest() {
  Dune::TestSuite testSuite("Q2 Topologic Associativity Forest");

  YaspUnitSquare grid;
  grid.globalRefine(3);
  using GridView = typename YaspUnitSquare::LeafGridView;
  GridView entity_set{grid.leafGridView()};

  auto flat_strategy = Dune::PDELab::BasisFactory::flatByEntity(entity_set);
  auto q2_pb = makeQkProtoBasis(flat_strategy, Dune::Indices::_2);
  Dune::PDELab::Impl::LeafTopologicAssociativityForest entity_ordering{q2_pb};

  entity_ordering.update();
  testQ2TopologicAssociativityForest(testSuite, entity_ordering);
}


template<class FixedMonomialTopologicAssociativityForest>
void testFixedMonomialTopologicAssociativityForest(Dune::TestSuite& testSuite, FixedMonomialTopologicAssociativityForest& entity_ordering) {
  static_assert(not FixedMonomialTopologicAssociativityForest::containerBlocked());
  static_assert(not FixedMonomialTopologicAssociativityForest::fixedSizePerGeometryTypeStatic());
  testSuite.check(entity_ordering.fixedSizePerGeometryType());

  static_assert(FixedMonomialTopologicAssociativityForest::maxContainerDepth() == 1);
  testSuite.check(not entity_ordering.containsGeometry(Dune::GeometryTypes::vertex));
  testSuite.check(not entity_ordering.containsGeometry(Dune::GeometryTypes::line));
  testSuite.check(entity_ordering.containsGeometry(Dune::GeometryTypes::quadrilateral));
  testSuite.check(not entity_ordering.containsGeometry(Dune::GeometryTypes::triangle));
  testSuite.check(entity_ordering.containsCodim(0));
  testSuite.check(not entity_ordering.containsCodim(1));
  testSuite.check(not entity_ordering.containsCodim(2));
  testSuite.check(entity_ordering.singleCodim());
  testSuite.check(entity_ordering.disjointCodimClosure());
  testSuite.check(entity_ordering.maxLocalCount() == 6);

  auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::quadrilateral);
  testSuite.check(entity_ordering.blockCount(gt_index) == 6);

  for(const auto& element : elements(entity_ordering.entitySet())) {
    auto entity_index = entity_ordering.entitySet().indexSet().index(element);
    testSuite.check(entity_ordering.blockCount(gt_index, entity_index) == 6);
    auto tpath = Dune::TypeTree::treePath();
    auto cir = entity_ordering.containerIndexRange(tpath, gt_index, entity_index);
    testSuite.check(cir.size() == 6);
    testSuite.check(cir[0].size() == 1);
    testSuite.check(cir[0][0] == 0);
    testSuite.check(entity_ordering.containerSize(cir[0], gt_index, entity_index) == 0);
  }

  testSuite.check(entity_ordering.dimension() ==
    entity_ordering.entitySet().indexSet().size(0) * entity_ordering.blockCount(gt_index)
  );
  testSuite.check(entity_ordering.dimension() == entity_ordering.blockCount());
}


void testFixedMonomialTopologicAssociativityForest() {
  Dune::TestSuite testSuite("Fixed Size Monomial Topologic Associativity Forest");

  YaspUnitSquare grid;
  grid.globalRefine(3);
  using GridView = typename YaspUnitSquare::LeafGridView;
  GridView entity_set{grid.leafGridView()};

  using CellMapper = Dune::SingleCodimSingleGeomTypeMapper<GridView, 0> ;
  CellMapper cellmapper(entity_set);
  using MonomFEM = Dune::PDELab::VariableMonomLocalFiniteElementMap<CellMapper, float, double, GridView::dimension>;
  auto monomfem = std::make_shared<MonomFEM>(cellmapper, Dune::GeometryTypes::quadrilateral, 2);

  auto strategy = Dune::PDELab::BasisFactory::flatByEntity(entity_set);
  auto monom_pb = makeProtoBasis(strategy, monomfem);

  using ProtoBasis = decltype(monom_pb);
  using TopologicAssociativityForest = Dune::PDELab::Impl::LeafTopologicAssociativityForest<ProtoBasis>;
  auto entity_ordering = std::make_shared<TopologicAssociativityForest>(monom_pb);

  entity_ordering->update();

  testFixedMonomialTopologicAssociativityForest(testSuite, *entity_ordering);

  using ArrayTopologicAssociativityForest = Dune::PDELab::Impl::ArrayTopologicAssociativityForest<decltype(strategy), TopologicAssociativityForest, 2>;

  ArrayTopologicAssociativityForest p_entity_ordering{std::array{entity_ordering,entity_ordering}, strategy};
  p_entity_ordering.update();
  testFixedMonomialTopologicAssociativityForest(testSuite, p_entity_ordering.child(0));

  static_assert(not p_entity_ordering.containerBlocked());
  static_assert(not p_entity_ordering.fixedSizePerGeometryTypeStatic());
  testSuite.check(p_entity_ordering.fixedSizePerGeometryType());

  static_assert(p_entity_ordering.maxContainerDepth() == 1);
  testSuite.check(not p_entity_ordering.containsGeometry(Dune::GeometryTypes::vertex));
  testSuite.check(not p_entity_ordering.containsGeometry(Dune::GeometryTypes::line));
  testSuite.check(p_entity_ordering.containsGeometry(Dune::GeometryTypes::quadrilateral));
  testSuite.check(not p_entity_ordering.containsGeometry(Dune::GeometryTypes::triangle));
  testSuite.check(p_entity_ordering.containsCodim(0));
  testSuite.check(not p_entity_ordering.containsCodim(1));
  testSuite.check(not p_entity_ordering.containsCodim(2));
  testSuite.check(p_entity_ordering.singleCodim());
  testSuite.check(p_entity_ordering.disjointCodimClosure());
  testSuite.check(p_entity_ordering.maxLocalCount() == 2*6);

  auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::quadrilateral);
  testSuite.check(p_entity_ordering.blockCount(gt_index) == 2*6);

  for(const auto& element : elements(p_entity_ordering.entitySet())) {
    auto entity_index = p_entity_ordering.entitySet().indexSet().index(element);
    testSuite.check(p_entity_ordering.blockCount(gt_index, entity_index) == 2*6);

    auto tpath0 = Dune::TypeTree::treePath(0);
    auto cir0 = p_entity_ordering.containerIndexRange(tpath0, gt_index, entity_index);
    testSuite.check(cir0.size() == 6);
    testSuite.check(cir0[0].size() == 1);
    testSuite.check(cir0[0][0] == 0);
    testSuite.check(p_entity_ordering.containerSize(cir0[0], gt_index, entity_index) == 0);

    auto tpath1 = Dune::TypeTree::treePath(1);
    auto cir1 = p_entity_ordering.containerIndexRange(tpath1, gt_index, entity_index);
    testSuite.check(cir1.size() == 6);
    testSuite.check(cir1[0].size() == 1);
    testSuite.check(cir1[0][0] == 6);
    testSuite.check(p_entity_ordering.containerSize(cir1[0], gt_index, entity_index) == 0);
  }

  testSuite.check(p_entity_ordering.blockCount() ==
    p_entity_ordering.entitySet().indexSet().size(0) * p_entity_ordering.blockCount(gt_index)
  );
  testSuite.check(p_entity_ordering.dimension() == p_entity_ordering.blockCount());
}

template<class VariableMonomialTopologicAssociativityForest>
void testVariableMonomialTopologicAssociativityForest(Dune::TestSuite& testSuite, VariableMonomialTopologicAssociativityForest& entity_ordering) {
  static_assert(not VariableMonomialTopologicAssociativityForest::containerBlocked());
  static_assert(not VariableMonomialTopologicAssociativityForest::fixedSizePerGeometryTypeStatic());
  testSuite.check(not entity_ordering.fixedSizePerGeometryType());

  static_assert(VariableMonomialTopologicAssociativityForest::maxContainerDepth() == 1);
  testSuite.check(not entity_ordering.containsGeometry(Dune::GeometryTypes::vertex));
  testSuite.check(not entity_ordering.containsGeometry(Dune::GeometryTypes::line));
  testSuite.check(entity_ordering.containsGeometry(Dune::GeometryTypes::quadrilateral));
  testSuite.check(not entity_ordering.containsGeometry(Dune::GeometryTypes::triangle));
  testSuite.check(entity_ordering.containsCodim(0));
  testSuite.check(not entity_ordering.containsCodim(1));
  testSuite.check(not entity_ordering.containsCodim(2));
  testSuite.check(entity_ordering.singleCodim());
  testSuite.check(entity_ordering.disjointCodimClosure());
  testSuite.check(entity_ordering.maxLocalCount() == 6);

  auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::quadrilateral);

  std::size_t order = 0;
  for(const auto& element : elements(entity_ordering.entitySet())) {
    auto& fe = entity_ordering.protoBasis().finiteElementMap().getFEM(order++ % 3);
    auto entity_index = entity_ordering.entitySet().indexSet().index(element);
    testSuite.check(entity_ordering.blockCount(gt_index, entity_index) == fe.size());
    auto tpath = Dune::TypeTree::treePath();
    auto cir = entity_ordering.containerIndexRange(tpath, gt_index, entity_index);
    testSuite.check(cir.size() == fe.size());
    testSuite.check(cir[0].size() == 1);
    testSuite.check(cir[0][0] == 0);
    testSuite.check(entity_ordering.containerSize(cir[0], gt_index, entity_index) == 0);
  }

  testSuite.check(entity_ordering.entitySet().indexSet().size(0) == 4);

  testSuite.check(entity_ordering.dimension() ==
    1 + 3 + 6 + 1
  );
  testSuite.check(entity_ordering.dimension() == entity_ordering.blockCount());
}


void testVariableMonomialTopologicAssociativityForest() {
  Dune::TestSuite testSuite("Fixed Size Monomial Topologic Associativity Forest");

  Dune::FieldVector<double,2> L(1.0);
  std::array<int,2> N;
  std::fill(begin(N), end(N), 2); // for this test we require exactly 4 entities
  Dune::YaspGrid<2> grid{L,N};
  using GridView = typename Dune::YaspGrid<2>::LeafGridView;
  GridView entity_set{grid.leafGridView()};

  using CellMapper = Dune::SingleCodimSingleGeomTypeMapper<GridView, 0> ;
  CellMapper cellmapper(entity_set);
  using MonomFEM = Dune::PDELab::VariableMonomLocalFiniteElementMap<CellMapper, float, double, GridView::dimension>;
  auto monomfem = std::make_shared<MonomFEM>(cellmapper, Dune::GeometryTypes::quadrilateral, 2);

  std::size_t order = 0;
  for (const auto &e : elements(entity_set))
    monomfem->setOrder(e, order++ % 3);

  auto strategy = Dune::PDELab::BasisFactory::flatByEntity(entity_set);
  auto monom_pb = makeProtoBasis(strategy, monomfem);

  using ProtoBasis = decltype(monom_pb);
  using TopologicAssociativityForest = Dune::PDELab::Impl::LeafTopologicAssociativityForest<ProtoBasis>;
  auto entity_ordering = std::make_shared<TopologicAssociativityForest>(monom_pb);

  entity_ordering->update();

  testVariableMonomialTopologicAssociativityForest(testSuite, *entity_ordering);

  using namespace Dune::Indices;
  using TupleTopologicAssociativityForest = Dune::PDELab::Impl::TupleTopologicAssociativityForest<decltype(strategy), TopologicAssociativityForest, TopologicAssociativityForest>;

  TupleTopologicAssociativityForest p_entity_ordering{std::tuple{entity_ordering,entity_ordering}, strategy};
  p_entity_ordering.update();
  testVariableMonomialTopologicAssociativityForest(testSuite, p_entity_ordering.child(_0));
  testVariableMonomialTopologicAssociativityForest(testSuite, p_entity_ordering.child(_1));

  static_assert(not TupleTopologicAssociativityForest::containerBlocked());
  static_assert(not TupleTopologicAssociativityForest::fixedSizePerGeometryTypeStatic());
  testSuite.check(not p_entity_ordering.fixedSizePerGeometryType());

  static_assert(TupleTopologicAssociativityForest::maxContainerDepth() == 1);
  testSuite.check(not p_entity_ordering.containsGeometry(Dune::GeometryTypes::vertex));
  testSuite.check(not p_entity_ordering.containsGeometry(Dune::GeometryTypes::line));
  testSuite.check(p_entity_ordering.containsGeometry(Dune::GeometryTypes::quadrilateral));
  testSuite.check(not p_entity_ordering.containsGeometry(Dune::GeometryTypes::triangle));
  testSuite.check(p_entity_ordering.containsCodim(0));
  testSuite.check(not p_entity_ordering.containsCodim(1));
  testSuite.check(not p_entity_ordering.containsCodim(2));
  testSuite.check(p_entity_ordering.singleCodim());
  testSuite.check(p_entity_ordering.disjointCodimClosure());
  testSuite.check(p_entity_ordering.maxLocalCount() == 2*6);

  auto gt_index = Dune::GlobalGeometryTypeIndex::index(Dune::GeometryTypes::quadrilateral);

  order = 0;
  for(const auto& element : elements(p_entity_ordering.entitySet())) {
    auto& fe = p_entity_ordering.child(_0).protoBasis().finiteElementMap().getFEM(order++ % 3);
    auto entity_index = p_entity_ordering.entitySet().indexSet().index(element);
    testSuite.check(p_entity_ordering.blockCount(gt_index, entity_index) == 2*fe.size());

    testSuite.check(p_entity_ordering.child(_0).blockCount(gt_index, entity_index) == fe.size());
    auto tpath0 = Dune::TypeTree::treePath(_0);
    auto cir0 = p_entity_ordering.containerIndexRange(tpath0, gt_index, entity_index);
    testSuite.check(cir0.size() == fe.size());
    testSuite.check(cir0[0].size() == 1);
    testSuite.check(cir0[0][0] == 0);
    testSuite.check(p_entity_ordering.containerSize(cir0[0], gt_index, entity_index) == 0);

    testSuite.check(p_entity_ordering.child(_1).blockCount(gt_index, entity_index) == fe.size());
    auto tpath1 = Dune::TypeTree::treePath(_1);
    auto cir1 = p_entity_ordering.containerIndexRange(tpath1, gt_index, entity_index);
    testSuite.check(cir1.size() == fe.size());
    testSuite.check(cir1[0].size() == 1);
    testSuite.check(cir1[0][0] == fe.size());
    testSuite.check(p_entity_ordering.containerSize(cir1[0], gt_index, entity_index) == 0);
  }

  testSuite.check(p_entity_ordering.dimension() ==
    2*(1 + 3 + 6 + 1)
  );
  testSuite.check(p_entity_ordering.blockCount() == p_entity_ordering.dimension());
}

int main(int argc, char **argv) {
  Dune::MPIHelper::instance(argc, argv);
  testQ1TopologicAssociativityForest();
  testQ2TopologicAssociativityForest();
  testFixedMonomialTopologicAssociativityForest();
  testVariableMonomialTopologicAssociativityForest();
}
