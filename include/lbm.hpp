#pragma once

#include "lbm_c.hpp"
#include "math_util.hpp"

#include <array>
#include <barrier>
namespace LBM {

constexpr std::array<float, LBM_CONSTANTS::LATTICE_COUNT> WEIGHTS = {
    1.0 / 9.0,  1.0 / 36.0, 1.0 / 9.0,  1.0 / 36.0, 1.0 / 9.0,
    1.0 / 36.0, 1.0 / 9.0,  1.0 / 36.0, 4.0 / 9.0};

constexpr std::array<Vect<float>, LBM_CONSTANTS::LATTICE_COUNT> DISC_VELOCITY =
    {
        Vect{1.f, 0.f},   Vect{+1.f, -1.f}, Vect{0.f, -1.f},
        Vect{-1.f, -1.f}, Vect{-1.f, 0.f},  Vect{-1.f, 1.f},
        Vect{0.f, 1.f},   Vect{+1.f, 1.f},  Vect{0.f, 0.f},
};

constexpr std::array<Vect<int32_t>, LBM_CONSTANTS::LATTICE_COUNT>
    DISC_VELOCITY_INT = {
        Vect{1, 0},  Vect{+1, -1}, Vect{0, -1}, Vect{-1, -1}, Vect{-1, 0},
        Vect{-1, 1}, Vect{0, 1},   Vect{+1, 1}, Vect{0, 0},
};

void initialize();

void update(bool k, std::size_t from, std::size_t to,
            std::barrier<> &sync_barrier);
} // namespace LBM
