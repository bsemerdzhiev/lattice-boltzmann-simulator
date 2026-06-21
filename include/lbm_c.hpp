#pragma once

#include "math_util.hpp"
#include <cstdint>
namespace LBM_CONSTANTS {
constexpr int32_t HEIGHT = 50;
constexpr int32_t WIDTH = 300;
constexpr int32_t TOTAL_CELLS = HEIGHT * WIDTH;

constexpr int32_t LATTICE_COUNT = 9;

constexpr float delta_x = 1.0f;
constexpr float delta_t = 1.0f;

constexpr Vect<float> HORIZONTAL_VELOCITY{0.04f, 0.f};

constexpr float CYLINDER_RADIUS = HEIGHT / 9.f;

constexpr float REYNOLDS_NUMBER = 80.f;

constexpr float KINEMATIC_VISCOSITY =
    (HORIZONTAL_VELOCITY.x * CYLINDER_RADIUS) / REYNOLDS_NUMBER;

constexpr float RELAXATION_OMEGA = 1.f / (3.f * KINEMATIC_VISCOSITY + 0.5f);

}; // namespace LBM_CONSTANTS
