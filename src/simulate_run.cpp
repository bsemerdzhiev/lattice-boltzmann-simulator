#include "simulate_run.hpp"
#include "GLFW/glfw3.h"
#include "lbm.hpp"

int32_t k = 1;

void SimulateRun::run() {
  LBM::update(k);
  k ^= 1;
}
