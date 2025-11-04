
#ifndef DUNE_PDELAB_COMMON_ENTITYMUTEXES_HH
#define DUNE_PDELAB_COMMON_ENTITYMUTEXES_HH

#include <dune/pdelab/common/spinlocks.hh>

#include <dune/geometry/typeindex.hh>

#include <array>
#include <numeric>
#include <concepts>

namespace Dune::PDELab
{

  // locks is a range of lockable objects

  /**
   * @brief Mutexes for entities of a grid view
   * @details Provides mutexes for each entity of a given grid view.
   * Each entity mutex can be accessed by the entity, or its geometry type index and entity index.
   * The geometry type index can be obtained using `GlobalGeometryTypeIndex::index(entity.type())`
   * and the entity index using `gridView.indexSet().index(entity)`.
   *
   * The type of the lock container can be customized with a template parameter whose type returned by
   * an indexation with the bracket operator `operator[]` on the lock container fulfills the Lockable
   * named requirement and is constructible by the number of locks (e.g. `std::vector<std::mutex>`).
   *
   * @tparam GV Type of the grid view
   * @tparam Locks_ Type of the lock container (default SpinLocks)
   */
  template <class GV, class Locks_ = SpinLocks>
    requires std::move_constructible<Locks_> && std::constructible_from<Locks_, std::size_t> && requires(Locks_ locks, std::size_t i) {
      { locks[i].lock() };
      { locks[i].unlock() };
      { locks[i].try_lock() } -> std::convertible_to<bool>;
    }
  class EntityMutexes {
  public:

    //! Type of the grid view
    using GridView = GV;
    //! Type of the locks container
    using Locks = Locks_;

    /**
     * @brief Constructs entity mutexes for a given grid view
     * @param grid_view The grid view whose entities will be managed
     */
    EntityMutexes(GridView const &grid_view)
    {
      update(grid_view);
    }

    //! Accesses the mutex for a given geometry type index and entity index
    decltype(auto) mutex(std::size_t gt_index, std::size_t entity_index) noexcept
    {
      return _locks[_offsets[gt_index] + entity_index];
    }

    //! Accesses the mutex for a given entity
    template <class Entity>
    decltype(auto) mutex(const Entity &entity)
    {
      return mutex(GlobalGeometryTypeIndex::index(entity.type()), _grid_view->indexSet().index(entity));
    }

    //! Updates the entity mutexes for a new grid view
    void update(const GridView &grid_view)
    {
      _grid_view = &grid_view;
      const auto &index_set = _grid_view->indexSet();
      std::fill(std::begin(_offsets), std::end(_offsets), 0);

      for (std::size_t codim = 0; codim <= GridView::dimension; ++codim)
        for (const auto &gt : index_set.types(codim))
          _offsets[GlobalGeometryTypeIndex::index(gt) + 1] = index_set.size(gt);

      std::partial_sum(std::begin(_offsets), std::end(_offsets), std::begin(_offsets));
      _locks = Locks(_offsets.back());
    }

    //! Accesses the associated grid view
    GridView const &gridView() const noexcept
    {
      return *_grid_view;
    }

  private:
    GridView const *_grid_view;
    std::array<std::size_t, GlobalGeometryTypeIndex::size(GridView::dimension) + 1> _offsets;
    Locks _locks;
  };

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_COMMON_ENTITYMUTEXES_HH
