#include "lbm_parallel.hpp"
#include "lbm.hpp"
#include "lbm_c.hpp"
#include <array>
#include <cstddef>
#include <iostream>
#include <pthread.h>
#include <thread>

constexpr std::size_t THREAD_COUNT = 8;

std::barrier<> sync_barrier{THREAD_COUNT + 1};

std::atomic<bool> exit_cond(false);

std::array<std::thread, THREAD_COUNT> threads;

void LBMParallel::init_threads() {
  LBM::initialize();

  std::size_t rows_per_thread =
      (LBM_CONSTANTS::HEIGHT + THREAD_COUNT - 1) / THREAD_COUNT;

  for (std::size_t i{0}; i < THREAD_COUNT; i++) {

    threads[i] = std::thread([rows_per_thread, i, &sync_barrier, &exit_cond]() {
      cpu_set_t cpuset;
      CPU_ZERO(&cpuset);
      CPU_SET(i, &cpuset);
      pthread_setaffinity_np(pthread_self(), sizeof(cpuset), &cpuset);

      bool k = 0;
      while (!exit_cond) {
        // std::cout << "DSA" << "\n";
        LBM::update(k, rows_per_thread * i,
                    std::min(static_cast<std::size_t>(LBM_CONSTANTS::HEIGHT),
                             (i + 1) * rows_per_thread),
                    sync_barrier);
        k ^= 1;
      }
    });
  }
}

void LBMParallel::run_iteration(bool is_last) {
  sync_barrier.arrive_and_wait();

  sync_barrier.arrive_and_wait();

  if (!is_last) {
    sync_barrier.arrive_and_wait();
  }
}
void LBMParallel::stop_threads() {
  exit_cond.exchange(true, std::memory_order_seq_cst);

  sync_barrier.arrive_and_wait();

  for (auto &x : threads) {
    x.join();
  }
}
