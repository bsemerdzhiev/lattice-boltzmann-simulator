#pragma once

#include "cell.hpp"
#include "lbm_c.hpp"
#include "math_util.hpp"

#include <array>
namespace LBM {
constexpr std::array<float, LATTICE_COUNT> WEIGHTS = {
    1.f / 9.f,  1.f / 36.f, 1.f / 9.f,  1.f / 36.f, 1.f / 9.f,
    1.f / 36.f, 1.f / 9.f,  1.f / 36.f, 4.f / 9.f};

constexpr std::array<Vect<float>, LATTICE_COUNT> DISC_VELOCITY = {
    Vect{1.f, 0.f},   Vect{+1.f, -1.f}, Vect{0.f, -1.f},
    Vect{-1.f, -1.f}, Vect{-1.f, 0.f},  Vect{-1.f, 1.f},
    Vect{0.f, 1.f},   Vect{+1.f, 1.f},  Vect{0.f, 0.f},
};

constexpr std::array<Vect<int32_t>, LATTICE_COUNT> DISC_VELOCITY_INT = {
    Vect{1, 0},  Vect{+1, -1}, Vect{0, -1}, Vect{-1, -1}, Vect{-1, 0},
    Vect{-1, 1}, Vect{0, 1},   Vect{+1, 1}, Vect{0, 0},
};

extern std::array<std::array<Cell, HEIGHT>, WIDTH> lattice;

void initialize();

void update(bool turn);
} // namespace LBM
