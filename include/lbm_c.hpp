#pragma once

#include <cstdint>
constexpr int32_t HEIGHT = 300;
constexpr int32_t WIDTH = 600;

constexpr int32_t LATTICE_COUNT = 9;

constexpr float delta_x = 1.f;
constexpr float delta_t = 1.f;
constexpr float THAO = 1.f;
constexpr float lattice_speed = 1.f;

constexpr float C_s_sq = 1.f / 3.f;

constexpr float U_lid = 0.515f;

constexpr float lattice_speed_squred = lattice_speed * lattice_speed;
constexpr float lattice_speed_quadrupelled =
    lattice_speed_squred * lattice_speed_squred;
