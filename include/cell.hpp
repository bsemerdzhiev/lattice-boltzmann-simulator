#pragma once

#include "lbm_c.hpp"
#include "math_util.hpp"
#include <array>

template <typename T>
using Lattice =
    std::array<std::array<T, LBM_CONSTANTS::WIDTH>, LBM_CONSTANTS::HEIGHT>;

namespace Cell {
extern Lattice<float> density;
extern Lattice<Vect<float>> velocity;

extern Lattice<bool> blockade;

extern Lattice<std::array<float, LBM_CONSTANTS::LATTICE_COUNT>> pdf[2];
extern Lattice<std::array<float, LBM_CONSTANTS::LATTICE_COUNT>> pdf_eq[2];
}; // namespace Cell
