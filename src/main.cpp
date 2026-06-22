#include "benchmark_simulation.hpp"
#include "render_simulation.hpp"
#include <cstdint>

constexpr bool RENDER = false;

int32_t main() {
  if constexpr (RENDER) {
    RenderSimulation::render();
  } else {
    BenchmarkRun::run();
  }

  return 0;
}
