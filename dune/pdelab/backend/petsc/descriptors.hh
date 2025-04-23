#pragma once

#include <petscsys.h>

namespace Dune::PDELab::Petsc {
template <typename GFS, typename E>
class VectorContainer;

template <typename GFSV, typename GFSU, typename ET>
class MatrixContainer;

template <class M>
struct MatrixPatternInserter;

struct VectorBackend {
  using size_type = std::size_t;

  struct Traits {
    static const size_type max_blocking_depth = 0;
  };

  template <typename GFS>
  bool blocked(const GFS &gfs) const
  {
    return false;
  }
};

struct MatrixBackend {
  using size_type = std::size_t;
  size_type avg_nz_per_row;

  MatrixBackend(size_type avg_nz_per_row) : avg_nz_per_row{avg_nz_per_row} {}

  template <typename Matrix, typename GFSV, typename GFSU>
  using Pattern = MatrixPatternInserter<typename Matrix::Container>;

  template <typename VV, typename VU, typename E>
  struct MatrixHelper {
    using type = MatrixContainer<typename VV::GridFunctionSpace, typename VU::GridFunctionSpace, E>;
  };
};
} // namespace Dune::PDELab::Petsc
