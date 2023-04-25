#ifndef DUNE_PDELAB_BASIS_BACKEND_INDEX_TREE_HH
#define DUNE_PDELAB_BASIS_BACKEND_INDEX_TREE_HH

#include <dune/functions/functionspacebases/indextree.hh>

#include <vector>
#include <array>
#include <tuple>
#include <concepts>
#include <type_traits>

namespace Dune::PDELab {

  struct IndexTreeBackend {

  private:
    template<class T>
    static constexpr bool recursiveStaticUniform() {
      if constexpr (T::isTypeUniform and requires { index_constant<T::size()>{}; })
        return (T::size() == 0) or recursiveStaticUniform<std::decay_t<decltype(std::declval<T>()[Indices::_0])>>();
      else
        return false;
    }

  public:
    template<class T>
    static auto makeVector() {
      if constexpr (recursiveStaticUniform<T>())
        return Dune::Functions::UniformIndexTree<T>{};
      else
        return Dune::Functions::TypeUniformIndexTree<T>{};
    }

    template<class T, std::size_t k>
    static auto makeArray() {
      if constexpr (T::isUniform)
        return Dune::Functions::StaticUniformIndexTree<T, k>{};
      else
        return Dune::Functions::StaticTypeUniformIndexTree<T, k>{};
    }

    template<class... T>
    static auto makeTuple() {
      return Dune::Functions::StaticNonUniformIndexTree<T...>{};
    }

    template<class>
    static auto makeField() {
      return Dune::Functions::EmptyIndexTree{};
    }

    template<class C>
    struct block_type {
      static_assert(C::isTypeUniform);
      using type = std::decay_t<decltype(std::declval<C>()[Indices::_0])>;
    };

    template<class T0, class... T>
    struct block_type<Dune::Functions::StaticNonUniformIndexTree<T0, T...>> {
      static_assert((std::common_with<T0,T> && ...), "Tuple arguments do not have a common block type");
      using type = std::common_type_t<T0,T...>;
    };
  };

} // namespace Dune::PDELab::inline Experimental

#endif // DUNE_PDELAB_BASIS_BACKEND_INDEX_TREE_HH
