#pragma once

#include "cell.hpp"
#include "lbm_c.hpp"
#include "math_util.hpp"

#include <array>
namespace LBM {

constexpr std::array<double, LBM_CONSTANTS::LATTICE_COUNT> WEIGHTS = {
    1.0 / 9.0,  1.0 / 36.0, 1.0 / 9.0,  1.0 / 36.0, 1.0 / 9.0,
    1.0 / 36.0, 1.0 / 9.0,  1.0 / 36.0, 4.0 / 9.0};

constexpr std::array<Vect<double>, LBM_CONSTANTS::LATTICE_COUNT> DISC_VELOCITY =
    {
        Vect{1.0, 0.0},   Vect{+1.0, -1.0}, Vect{0.0, -1.0},
        Vect{-1.0, -1.0}, Vect{-1.0, 0.0},  Vect{-1.0, 1.0},
        Vect{0.0, 1.0},   Vect{+1.0, 1.0},  Vect{0.0, 0.0},
};

constexpr std::array<Vect<int32_t>, LBM_CONSTANTS::LATTICE_COUNT>
    DISC_VELOCITY_INT = {
        Vect{1, 0},  Vect{+1, -1}, Vect{0, -1}, Vect{-1, -1}, Vect{-1, 0},
        Vect{-1, 1}, Vect{0, 1},   Vect{+1, 1}, Vect{0, 0},
};

template <typename T>
using Lattice =
    std::array<std::array<T, LBM_CONSTANTS::WIDTH>, LBM_CONSTANTS::HEIGHT>;

extern Lattice<Cell> lattice;
extern Lattice<bool> blockade;

void initialize();

void update(bool turn);
} // namespace LBM
