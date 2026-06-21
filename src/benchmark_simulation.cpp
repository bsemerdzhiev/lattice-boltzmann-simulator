#include "benchmark_simulation.hpp"
#include "lbm.hpp"

#include <chrono>
#include <cstdint>
#include <iostream>
#include <ratio>

constexpr int32_t TIMESTEPS = 150;

void BenchmarkRun::run() {
  std::cout << std::format("Simulation started\n");
  LBM::initialize();

  auto tic = std::chrono::high_resolution_clock::now();

  bool k = 0;
  for (int32_t i{0}; i < TIMESTEPS; i++) {
    LBM::update(k);
    k ^= 1;
  }

  auto tac = std::chrono::high_resolution_clock::now();

  std::cout << std::format(
      "CPU time passed: {}ms\n",
      std::chrono::duration<float, std::milli>(tac - tic).count());
}
