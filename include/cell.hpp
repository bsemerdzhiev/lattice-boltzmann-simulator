#pragma once

#include "lbm_c.hpp"
#include <array>

template <typename T>
using Lattice =
    std::array<std::array<T, LBM_CONSTANTS::WIDTH>, LBM_CONSTANTS::HEIGHT>;

template <typename T>
using PdfLattice = std::array<std::array<std::array<T, LBM_CONSTANTS::WIDTH>,
                                         LBM_CONSTANTS::LATTICE_COUNT>,
                              LBM_CONSTANTS::HEIGHT>;

namespace Cell {
extern Lattice<float> density;
extern Lattice<float> velocity_x;
extern Lattice<float> velocity_y;

extern Lattice<int32_t> blockade;

extern PdfLattice<float> pdf[2];
extern PdfLattice<float> pdf_eq[2];
}; // namespace Cell
