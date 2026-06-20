#pragma once

#include "lbm_c.hpp"
#include "math_util.hpp"
#include <array>

struct Cell {
  float density;
  Vect velocity;

  std::array<float, LATTICE_COUNT> pdf;
};
