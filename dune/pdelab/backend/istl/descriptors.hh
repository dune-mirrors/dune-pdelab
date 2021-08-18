// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_DESCRIPTORS_HH
#define DUNE_PDELAB_BACKEND_ISTL_DESCRIPTORS_HH

#include <dune/pdelab/finiteelementmap/utility.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/istl/forwarddeclarations.hh>
#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/backend/istl/utility.hh>
#include <dune/common/tuplevector.hh>
#include <cstddef>

namespace Dune {
  namespace PDELab {

#ifndef DOXYGEN
      template<typename T>
      constexpr bool deactivate_standard_blocking_for_ordering(const T&)
      {
        return false;
      }
#endif

    namespace ISTL {

      //! The type of blocking employed at this node in the function space tree.
      enum class Blocking
      {
        //! No blocking at this level.
        none,
        Merged = none,
        //! Creates one macro block for each child space, each block is a BlockVector / BCRS matrix.
        bcrs,
        Dynamic = bcrs,
        //! Create fixed size blocks that each group together a fixed number of DOFs from each child space.
        /**
         * Creates a block structure with fixed size child blocks that each group together a fixed number
         * of DOFs from each child space. Typically used for entity-wise blocking of DOFs across child spaces
         * for e.g. velocity in flow problems or concentrations in multi-component transport.
         *
         * \note This type of blocking cannot be nested due to limitations in ISTL.
         */
        fixed,
        Static = fixed,
      };

      //! Tag describing an ISTL BlockVector backend.
      struct vector_backend_tag {};

      template<Blocking blocking = Blocking::none, std::size_t block_size_ = 0>
      struct VectorBackend
      {

        using tag = vector_backend_tag;

        using size_type = std::size_t;

        static const size_type blockSize = block_size_;

        struct Traits
        {
          static const Blocking block_type = blocking;
          static const size_type block_size = block_size_;

          static const bool blocked = blocking != Blocking::none;

          static const size_type max_blocking_depth = blocked ? 1 : 0;
        };

        template<typename GFS>
        bool blocked(const GFS& gfs) const
        {
          if (deactivate_standard_blocking_for_ordering(gfs.orderingTag()))
            return false;
          // We have to make an exception for static blocking and block_size == 1:
          // In that case, the ISTL backends expect the redundant index information
          // at that level to be elided, and keeping it in here will break vector
          // and matrix accesses.
          // To work around that problem, we override the user and just turn off
          // blocking internally.
          // A bock size of 0 also needs special handling, as it is actually a marker for
          // automatic block size deduction
          return Traits::blocked && (blocking != Blocking::fixed || !GFS::isLeaf || block_size_ > 1 || block_size_ == 0);
        }

      };


      template<Blocking NodeBlocking = Blocking::Merged>
      struct BackendOptions {
        using size_type = std::size_t;

        template<typename GFS>
        auto blocked(const GFS& gfs) const {
          return blocked<GFS>();
        }

        template<typename GFS>
        static constexpr auto blocked() {
          return std::bool_constant<(NodeBlocking != Blocking::Merged)>{};
        }
      };

      namespace Impl {
        template<class N>
        using StaticDegreeConcept = decltype((std::integral_constant<std::size_t, N::degree()>{}, true));

        template<class N>
        using StaticSizeConcept = decltype((std::integral_constant<std::size_t, N::size()>{}, true));

        template<class... T>
        using CommonTypeConcept = decltype((std::declval<std::common_type_t<T...>>(), true));

        template<class T>
        using IndexableConcept = decltype((std::declval<T>()[0], true));

        // This function is found by ADL
        template <bool isLocalOrdering, class GFS, class Field, Blocking NodeBlocking>
        static auto registerVectorContainer(BackendOptions<NodeBlocking>) {
          using OrderingTag = typename GFS::OrderingTag;
          const bool GridViewBlocking = (std::is_same<OrderingTag,DefaultLeafOrderingTag>{} or std::is_same<OrderingTag,EntityBlockedOrderingTag>{});
          if constexpr ((not isLocalOrdering) and GridViewBlocking)
          {
            auto block = registerVectorContainer<true, GFS,Field>(BackendOptions<NodeBlocking>{});
            return Dune::BlockVector<decltype(block)>{};
            // return block;
          } else {
            if constexpr (NodeBlocking == Blocking::Static) {
              if constexpr (GFS::isLeaf) {
                  using FEM = typename GFS::Traits::FiniteElementMap;
                  static_assert(Std::is_detected<FiniteElementMapBlockSize, FEM>{},
                                "Blocking::Static is only allowed in basis with "
                                "compile-time known finite lement map sizes");
                  static_assert(isLocalOrdering == true);
                  const std::size_t size = finiteElementMapBlockSize<FEM>();
                  static_assert(size > 0, "Static blocking is not possible for this finite element map."
                    "Try specify Dynamic blocking, or to create your own FiniteElementMap with contexpr"
                    "`size(GeometryType)`");
                  using Block = Dune::FieldVector<Field, size>;
                  return Block{};
              } else if constexpr (GFS::isPower) {
                using Child = typename GFS::ChildType;
                auto block = registerVectorContainer<isLocalOrdering,Child, Field>(
                    typename Child::Traits::Backend{});
                using Block = decltype(block);
                static_assert(
                  not std::is_same<Block,Field>{},
                  "Blocking::Static needs that sub-blocks are container (different than Blocking::Merged)");
                static_assert(
                    Std::is_detected<Impl::StaticDegreeConcept, GFS>{},
                    "Blocking::Static is not possible in dynamic power nodes");
                return Dune::FieldVector<Block, GFS::degree()>{};
              } else if constexpr (GFS::isComposite) {
                const auto sequence =
                    std::make_index_sequence<std::size_t(GFS::degree())>();
                auto child_container = [&](auto i) {
                  using Child = TypeTree::Child<GFS, i>;
                  return registerVectorContainer<isLocalOrdering,Child, Field>(
                      typename Child::Traits::Backend{});
                };
                return Dune::unpackIntegerSequence(
                    [&](auto... indices) {
                      if constexpr (Std::is_detected<Impl::CommonTypeConcept,
                                                      decltype(child_container(
                                                          indices))...>{}) {
                        using Block = std::common_type_t<decltype(child_container(
                            indices))...>;
                        static_assert(
                          not std::is_same<Block,Field>{},
                          "Blocking::Static needs that sub-blocks are container (different than Blocking::Merged)");
                        return Dune::FieldVector<Block, sizeof...(indices)>{};
                      } else {
                        return Dune::makeTupleVector(child_container(indices)...);
                      }
                    },
                    sequence);
              }
            } else if constexpr (NodeBlocking == Blocking::Dynamic) {
              if constexpr (GFS::isLeaf) {
                static_assert(isLocalOrdering == true);
                return Dune::BlockVector<Field>{};
              } else if constexpr (GFS::isPower) {
                using Child = typename GFS::ChildType;
                auto block = registerVectorContainer<isLocalOrdering,Child, Field>(
                    typename Child::Traits::Backend{});
                using Block = decltype(block);
                return Dune::BlockVector<Block>{};
              } else if constexpr (GFS::isComposite) {
                const auto sequence =
                    std::make_index_sequence<std::size_t(GFS::degree())>();

                auto child_container = [&](auto i) {
                  using Child = TypeTree::Child<GFS, i>;
                  return registerVectorContainer<isLocalOrdering,Child, Field>(
                      typename Child::Traits::Backend{});
                };
                auto block = Dune::unpackIntegerSequence(
                    [&](auto... indices) {
                      static_assert(
                          Std::is_detected<Impl::CommonTypeConcept,
                                            decltype(child_container(
                                                indices))...>{},
                          "Blocking::Merged is only possible in composite nodes "
                          "with children able to create common containers");
                      return std::common_type_t<decltype(child_container(
                          indices))...>{};
                    },
                    sequence);
                using Block = decltype(block);
                return Dune::BlockVector<Block>{};
              }
            } else if constexpr (NodeBlocking == Blocking::Merged) {
              if constexpr (GFS::isLeaf) {
                static_assert(isLocalOrdering == true);
                return Field{};
              } else if constexpr (GFS::isPower) {
                using Child = typename GFS::ChildType;
                auto block = registerVectorContainer<isLocalOrdering,Child, Field>(
                    typename Child::Traits::Backend{});
                using Block = decltype(block);
                if constexpr (Std::is_detected<Impl::StaticDegreeConcept,
                                                GFS>{} and
                              Std::is_detected<Impl::StaticSizeConcept,
                                                Block>{}) {
                  // this is the only case where we can ensure static sized field
                  const std::size_t size = Block::size() * GFS::degree();
                  using block_type = typename Block::block_type;
                  return Dune::FieldVector<block_type, size>{};
                } else if constexpr (Std::is_detected<IndexableConcept,Block>{}) {
                  using block_type = std::decay_t<decltype(std::declval<Block>()[0])>;
                  return Dune::BlockVector<block_type>{};
                } else {
                  return Block{};;
                }
              } else if constexpr (GFS::isComposite) {
                const auto sequence =
                    std::make_index_sequence<std::size_t(GFS::degree())>();

                auto child_container = [&](auto i) {
                  using Child = TypeTree::Child<GFS, i>;
                  return registerVectorContainer<isLocalOrdering,Child, Field>(
                      typename Child::Traits::Backend{});
                };
                auto block = Dune::unpackIntegerSequence(
                    [&](auto... indices) {
                      static_assert(
                          Std::is_detected<Impl::CommonTypeConcept,
                                            decltype(child_container(
                                                indices))...>{},
                          "Blocking::Merged is only possible in composite nodes "
                          "with children able to create common containers");
                      return std::common_type_t<decltype(child_container(
                          indices))...>{};
                    },
                    sequence);
                using Block = decltype(block);
                if constexpr (Std::is_detected<Impl::StaticDegreeConcept,
                                                GFS>{} and
                              Std::is_detected<Impl::StaticSizeConcept,
                                                Block>{}) {
                  // this is the only case where we can ensure static sized field
                  const std::size_t size = Block::size() * GFS::degree();
                  using block_type = typename Block::block_type;
                  return Dune::FieldVector<block_type, size>{};
                } else if constexpr (Std::is_detected<IndexableConcept,Block>{}) {
                  using block_type = std::decay_t<decltype(std::declval<Block>()[0])>;
                  return Dune::BlockVector<block_type>{};
                } else {
                  return Block{};
                }
              }
            } else {
              static_assert(AlwaysFalse<GFS>{}, "Not known blocking");
            }
          }
        }
      }

      template <class GFS, class Field, Blocking NodeBlocking>
      static auto registerVectorContainer(BackendOptions<NodeBlocking> opts) {
        return Impl::registerVectorContainer<false,GFS,Field>(opts);
      }

    } // namespace ISTL

  } // namespace PDELab
} // namespace Dune


#endif // DUNE_PDELAB_BACKEND_ISTL_DESCRIPTORS_HH
