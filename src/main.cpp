#include "benchmark_simulation.hpp"
#include "render_simulation.hpp"
#include <cstdint>
#include <string>

constexpr bool RENDER = false;

int32_t main(int32_t argc, char *argv[]) {
  if (argc > 1 && std::string(argv[1]) == "render") {
    RenderSimulation::render();
  } else {
    BenchmarkRun::run();
  }

  return 0;
}
