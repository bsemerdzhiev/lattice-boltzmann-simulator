#pragma once

#include "math_util.hpp"
#include <cstdint>
namespace LBM_CONSTANTS {
constexpr int32_t HEIGHT = 50;
constexpr int32_t WIDTH = 300;

constexpr int32_t LATTICE_COUNT = 9;

constexpr double delta_x = 1.0;
constexpr double delta_t = 1.0;

constexpr Vect<double> HORIZONTAL_VELOCITY{0.040, 0.0};

constexpr double CYLINDER_RADIUS = HEIGHT / 9.0;

constexpr double REYNOLDS_NUMBER = 80.0;

constexpr double KINEMATIC_VISCOSITY =
    (HORIZONTAL_VELOCITY.x * CYLINDER_RADIUS) / REYNOLDS_NUMBER;

constexpr double RELAXATION_OMEGA = 1.0 / (3.0 * KINEMATIC_VISCOSITY + 0.50);

}; // namespace LBM_CONSTANTS
