#pragma once

#include <cstdint>
constexpr int32_t HEIGHT = 100;
constexpr int32_t WIDTH = 100;

constexpr int32_t LATTICE_COUNT = 9;

constexpr float delta_x = 1.f;
constexpr float delta_t = 1.f;
constexpr float THAO = 0.01f;

constexpr float lattice_speed = 1.f;
constexpr float lattice_speed_squred = lattice_speed * lattice_speed;
constexpr float lattice_speed_quadrupelled =
    lattice_speed_squred * lattice_speed_squred;
