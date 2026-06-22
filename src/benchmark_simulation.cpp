#include "benchmark_simulation.hpp"
#include "lbm.hpp"
#include "lbm_parallel.hpp"

#include <atomic>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <ratio>

constexpr int32_t TIMESTEPS = 100'000;

void BenchmarkRun::run() {
  std::cout << std::format("Simulation started\n");

  LBMParallel::init_threads();

  auto tic = std::chrono::high_resolution_clock::now();

  for (int32_t i{0}; i < TIMESTEPS; i++) {
    LBMParallel::run_iteration(i == TIMESTEPS - 1);
  }
  LBMParallel::stop_threads();

  auto tac = std::chrono::high_resolution_clock::now();

  std::cout << std::format(
      "CPU time passed: {}ms\n",
      std::chrono::duration<float, std::milli>(tac - tic).count());
}
