#pragma once

#include "lbm_c.hpp"
#include "math_util.hpp"
#include <array>

struct Cell {
  double density;
  Vect<double> velocity;

  std::array<double, LBM_CONSTANTS::LATTICE_COUNT> pdf[2];
  std::array<double, LBM_CONSTANTS::LATTICE_COUNT> pdf_eq[2];
};
