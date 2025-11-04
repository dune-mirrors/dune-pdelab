#ifndef DUNE_PDELAB_BASIS_LOCALLOCK_HH
#define DUNE_PDELAB_BASIS_LOCALLOCK_HH

#include <dune/pdelab/common/entitymutexes.hh>

#include <dune/geometry/typeindex.hh>
#include <dune/common/hybridutilities.hh>

#include <algorithm>
#include <climits>
#include <limits>
#include <numeric>
#include <vector>
#include <array>



namespace Dune::PDELab {

  /** Local handle for managing locks on basis functions.
   * \details This class provides a lockable handle for basis functions,
   * allowing to lock and unlock basis functions associated with sub-entities
   * of the grid. The handle is bound to a specific basis lock and local basis view.
   *
   * \tparam BasisLock The BasisLock type managing the locks.
   */
  template <class GridView, class Locks>
  class LocalBasisLockHandle
  {
    static constexpr std::size_t unused_index = std::numeric_limits<std::size_t>::max();
    static constexpr std::array<std::size_t, 2> unused_slot = {unused_index, unused_index};

  public:
    //! Constructs a LocalBasisLockHandle for a given EntityMutexes.
    LocalBasisLockHandle(EntityMutexes<GridView, Locks> &entity_mutexes)
      : _entity_mutexes{std::addressof(entity_mutexes)}
    {
      // determine maximum number of subentities per node
      std::size_t max_subentities = 0;
      for (GeometryType gt : _entity_mutexes->gridView().indexSet().types(0))
      {
        std::size_t gt_max_subentities = 0;
        const auto &ref_el = referenceElement<double, GridView::dimension>(gt);
        for (int codim = 0; codim <= GridView::dimension; ++codim)
          gt_max_subentities += ref_el.size(codim);
        max_subentities = std::max<std::size_t>(max_subentities, gt_max_subentities);
      }

      // allocate storage for indices
      _indices.resize(max_subentities, unused_slot);
    }

    LocalBasisLockHandle(const LocalBasisLockHandle &) = default;
    LocalBasisLockHandle(LocalBasisLockHandle &&) = default;

    LocalBasisLockHandle &operator=(const LocalBasisLockHandle &) = default;
    LocalBasisLockHandle &operator=(LocalBasisLockHandle &&) = default;

    //! Binds the handle to a specific local basis view, preparing for locking.
    void bind(const auto &view) noexcept
    {
      unbind(); // clear previous indices
      // get subentity indices and insert them (assume gridview is the same for each node)
      bindNode(view.tree());
    }

    //! Locks all basis functions associated with the bound local basis view.
    void lock() noexcept
    {
      // spin lock until we acquire all of the mutexes
      while (not try_lock()) {};
    }

    //! Tries to lock all basis functions associated with the bound local basis view.
    [[nodiscard]] bool try_lock() noexcept
    {
      for (std::size_t i = 0; i != _index_count; ++i)
      {
        assert(i < _indices.size() and _indices[i][0] != unused_index);
        // try to lock every index on the vector
        if (not _entity_mutexes->mutex(_indices[i][0], _indices[i][1]).try_lock())
        {
          // entity was already locked, we have to roll back
          unlock(i);
          // ...and inform that we could not adquire the lock
          return false;
        }
      }
      // if all entities were locked, we succeded
      return true;
    }

    //! Unlocks all basis functions associated with the bound local basis view.
    void unlock() noexcept
    {
      unlock(_index_count);
    }

  private:
    // unlock the first lock_count locks in reverse order
    void unlock(std::size_t lock_count) noexcept
    {
      assert(lock_count <= _index_count);
      for (std::size_t j = lock_count; j != 0; --j)
        _entity_mutexes->mutex(_indices[j - 1][0], _indices[j - 1][1]).unlock();
    }

    // bind the node, collecting all sub-entity indices
    void bindNode(const auto &node) noexcept
      requires requires { node.subEntityIndices(); }
    {
      const auto &subentity_indices = node.subEntityIndices();
      if (_index_count == 0) {
        // buffer is empty, just copy all indices
        std::copy(subentity_indices.begin(), subentity_indices.end(), _indices.begin());
        _index_count = subentity_indices.size();
      } else {
        // buffer is partially filled, insert only new indices
        for (auto index : subentity_indices) {
          auto it = std::find_if(_indices.begin(), _indices.end(), [&](const auto &slot)
                                  { return (slot == index) | (slot == unused_slot); });
          assert(it != _indices.end());
          *it = index;
          _index_count = std::max<std::size_t>(_index_count, std::distance(_indices.begin(), it) + 1);
        }
      }
    }

    // bind the node, recursively visiting children
    template <class Node>
    void bindNode(const Node &node) noexcept
      requires(!requires { node.subEntityIndices(); })
    {
      auto degree = node.degree();
      Hybrid::forEach(range(node.degree()), [&](const auto i)
                      { bindNode(node.child(i)); });
    }

    //! Unbinds the handle from the current local basis view.
    void unbind() noexcept
    {
      if (std::exchange(_index_count, 0) == 0)
        return; // already unbound, nothing to do
      std::fill_n(_indices.begin(), _indices.size(), unused_slot);
    }

    EntityMutexes<GridView, Locks> * _entity_mutexes;
    std::vector<std::array<std::size_t, 2>> _indices;
    std::size_t _index_count = 0;
  };

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_BASIS_LOCALLOCK_HH
