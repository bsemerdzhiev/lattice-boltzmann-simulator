#include "simulate_run.hpp"
#include "lbm.hpp"

void SimulateRun::run() {
  LBM::initialize();

  int32_t k = 1;
  while (true) {
    LBM::update(k);
    k ^= 1;
  }
}
