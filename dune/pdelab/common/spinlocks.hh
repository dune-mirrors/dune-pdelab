#ifndef DUNE_PDELAB_COMMON_SPINLOCKS_HH
#define DUNE_PDELAB_COMMON_SPINLOCKS_HH

#include <vector>
#include <atomic>
#include <climits>
#include <cassert>
#include <utility>
#include <cstddef>

#include <iostream>

namespace Dune::PDELab {

/**
 * @brief Vector of bit spin locks
 * Each bit in the vector is a different lock.
 *
 * @warning Notice that locks don't respect the cache line sizes, so
 * heterogeneous thread access of neighboring locks will result in false sharing.
 */
class SpinLocks
{
    using T = unsigned char;
  static constexpr std::size_t T_BIT =  CHAR_BIT*sizeof(T);

public:

  /**
   * @brief Handle of an spin lock on a bit of a given memory address
   * @details Manupulates the acquisition and release of the lock and fullfils
   * the standard Lockable requirement by modifying a bit of a given memory address
   * @note The lifetime of the lock handle shall not outlive the SpinLocks object
   */
  class LockHandle {
    static_assert(std::atomic_ref<T>::is_always_lock_free);

    LockHandle(T& chunk, std::size_t bit)
      : _lock{static_cast<T>(1 << bit%T_BIT)}
      , _chunk_ptr{std::addressof(chunk)}
    {}

    friend class SpinLocks;

  public:
    // Copy constructor
    LockHandle(const LockHandle&) = default;

    // Copy assignment
    LockHandle& operator=(const LockHandle&) = default;

    // Move constructor
    LockHandle(LockHandle&& other) {
      _lock = std::exchange(other._lock, 0);
      _chunk_ptr = std::exchange(other._chunk_ptr, nullptr);
    }

    // Move assignment
    LockHandle& operator=(LockHandle&& other) noexcept {
      _lock = std::exchange(other._lock, 0);
      _chunk_ptr = std::exchange(other._chunk_ptr, nullptr);
      return *this;
    }

    /**
     * @brief Acquire ownership of the lock
     * If locked, this lock will spin until unlocked by another thread
     */
    void lock() {
      if (try_lock())
        return;
      do {
        wait();
      } while (not try_lock());
    }

    //! Try to acquire ownership of the lock. Return true in success
      [[nodiscard]] bool try_lock() noexcept {
      return not (std::atomic_ref{*_chunk_ptr}.fetch_or(_lock, std::memory_order_acquire) & _lock);
    }

    //! Relesase ownership of the lock
      void unlock() noexcept {
      std::atomic_ref{*_chunk_ptr}.fetch_and(~_lock, std::memory_order_release);
    }

    //! Compare two lock handles for equality
    [[nodiscard]] bool operator==(const LockHandle& other) const {
      return (_chunk_ptr == other._chunk_ptr) and (_lock == other._lock);
    }

  private:

    //! Wait until lock could be acquired
    void wait() const noexcept {
      do {
        pause();
      } while (locked());
    }

    //! Check if the lock is currently held by some thread
    [[nodiscard]] bool locked() const noexcept {
      return std::atomic_ref{*_chunk_ptr}.load(std::memory_order_relaxed) & _lock;
    }

    //! CPU pause instruction to reduce power consumption while spinning
    inline void pause() const noexcept {
#     if defined(__x86_64__) || defined(__i386__)
#       if defined(__GNUC__) || defined(__GNUG__)
          __builtin_ia32_pause();
#       elif defined(__clang__)
          __asm__ __volatile__ ("pause");
#       endif
#     endif
    }

    T _lock;
    T* _chunk_ptr;
  };

  //! Constructs a vector with a given size
  SpinLocks(std::size_t size = 0) {
    resize(size);
  }

  //! Deleted copy constructor
  SpinLocks(const SpinLocks&) = delete;
  //! Default move constructor
  SpinLocks(SpinLocks&&) = default;

  //! Deleted copy assignment
  SpinLocks& operator=(const SpinLocks&) = delete;
  //! Default move assignment
  SpinLocks& operator=(SpinLocks&&) = default;

  //! Get the size of the vector
  std::size_t size() const
  {
    return _size;
  }

  //! Access the i-th lock
  LockHandle operator[](std::size_t i) {
    return {_data[i/T_BIT], i - T_BIT*(i/T_BIT)};
  }

  //! Sets new size of the lock vector. All previous locks are dumped
  void resize(std::size_t new_size) {
    _size = new_size;
    std::size_t chunks = (_size+T_BIT-1)/T_BIT;
    _data.assign(chunks, T{} ^ T{});
  }

private:
  std::vector<T> _data;
  std::size_t _size;
};

} // namespace Dune::PDELab

#endif // DUNE_PDELAB_COMMON_SPINLOCKS_HH
