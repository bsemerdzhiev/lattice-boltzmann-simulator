#pragma once

#include "lbm_c.hpp"
#include "math_util.hpp"
#include <array>

struct Cell {
  float density;
  Vect<float> velocity;

  std::array<float, LATTICE_COUNT> pdf[2];
  std::array<float, LATTICE_COUNT> pdf_eq[2];
};
