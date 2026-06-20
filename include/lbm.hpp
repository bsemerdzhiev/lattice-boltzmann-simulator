#pragma once

#include "cell.hpp"
#include "lbm_c.hpp"
#include "math_util.hpp"

#include <array>
#include <utility>
namespace LBM {
constexpr std::array<float, LATTICE_COUNT> WEIGHTS = {
    1.f / 9.f,  1.f / 36.f, 1.f / 9.f,  1.f / 36.f, 1.f / 9.f,
    1.f / 36.f, 1.f / 9.f,  1.f / 36.f, 4.f / 9.f};

constexpr std::array<Vect, LATTICE_COUNT> DISC_VELOCITY = {
    Vect{1.f, 0.f},   Vect{+1.f, -1.f}, Vect{0.f, -1.f},
    Vect{-1.f, -1.f}, Vect{-1.f, 0.f},  Vect{-1.f, 1.f},
    Vect{0.f, 1.f},   Vect{+1.f, 1.f},  Vect{0.f, 0.f},
};

extern std::array<std::array<Cell, HEIGHT>, WIDTH> lattice;

void initialize();

void update();
} // namespace LBM
