#include "lbm.hpp"
#include "cell.hpp"
#include "lbm_c.hpp"
#include "math_util.hpp"
#include <cstddef>
#include <iostream>

void LBM::initialize() {
  for (std::size_t i{0}; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j++) {
      Cell::density[i][j] = 1.0f;
      Cell::velocity[i][j] = LBM_CONSTANTS::HORIZONTAL_VELOCITY;

      for (std::size_t z = 0; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        Cell::pdf[0][i][j][z] = WEIGHTS[z];
      }
    }
  }

  // compute the equilibrium velocities
  for (std::size_t i{0}; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j++) {
      if (Cell::blockade[i][j]) {
        continue;
      }
      // calculate the equilibrium distribution
      for (std::size_t z{0}; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        float velocity_q = (DISC_VELOCITY[z].dot_product(Cell::velocity[i][j]));

        float linear_term = (3.f * velocity_q);

        float quadratic_term = linear_term * linear_term / 2.f;

        float last_term =
            (3.f * Cell::velocity[i][j].dot_product(Cell::velocity[i][j])) /
            (2.f);

        Cell::pdf[0][i][j][z] =
            WEIGHTS[z] * Cell::density[i][j] *
            (1.f + linear_term + quadratic_term - last_term);
      }
    }
  }

  Vect<float> circle_coord{LBM_CONSTANTS::WIDTH / 5.f,
                           LBM_CONSTANTS::HEIGHT / 2.f};

  for (std::size_t i{0}; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j++) {
      if (circle_coord.euclid_dist(Vect<float>{1.f * j, 1.f * i}) <
          LBM_CONSTANTS::CYLINDER_RADIUS) {
        Cell::blockade[i][j] = 1;
      }
    }
  }
}

void LBM::update(bool k) {
  // set the pdf of the next step to zero
  for (std::size_t i{0}; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j++) {
      for (std::size_t z{0}; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        Cell::pdf[k ^ 1][i][j][z] = 0.f;
        Cell::pdf_eq[k ^ 1][i][j][z] = 0.f;
      }
    }
  }

  // right boundary outlet
  for (std::size_t i{0}; i < LBM_CONSTANTS::HEIGHT; i++) {
    // for the outlet, perform interpolation from the results of the
    // previous cell
    Cell::pdf[k][i][LBM_CONSTANTS::WIDTH - 1][3] =
        Cell::pdf[k][i][LBM_CONSTANTS::WIDTH - 2][3];
    Cell::pdf[k][i][LBM_CONSTANTS::WIDTH - 1][4] =
        Cell::pdf[k][i][LBM_CONSTANTS::WIDTH - 2][4];
    Cell::pdf[k][i][LBM_CONSTANTS::WIDTH - 1][5] =
        Cell::pdf[k][i][LBM_CONSTANTS::WIDTH - 2][5];
  }

  // compute the macroscopic density and velocity
  for (std::size_t i{0}; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j++) {
      if (Cell::blockade[i][j]) {
        continue;
      }
      Cell::density[i][j] = 0.f;
      Cell::velocity[i][j] = Vect{0.f, 0.f};
      for (std::size_t z{0}; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        Cell::density[i][j] += Cell::pdf[k][i][j][z];
        Cell::velocity[i][j] =
            Cell::velocity[i][j] + DISC_VELOCITY[z] * Cell::pdf[k][i][j][z];
      }
      Cell::velocity[i][j] = Cell::velocity[i][j] / Cell::density[i][j];
    }
  }

  // prescribe the Zou/He scheme
  for (std::size_t i{0}; i < LBM_CONSTANTS::HEIGHT; i++) {
    float f2 = Cell::pdf[k][i][0][2];
    float f3 = Cell::pdf[k][i][0][3];
    float f4 = Cell::pdf[k][i][0][4];
    float f5 = Cell::pdf[k][i][0][5];
    float f6 = Cell::pdf[k][i][0][6];
    float f8 = Cell::pdf[k][i][0][8];

    // skip settings the velocity for
    // the first and last elements of the first column
    if ((i > 0) & (i < LBM_CONSTANTS::HEIGHT - 1)) {
      Cell::velocity[i][0] = LBM_CONSTANTS::HORIZONTAL_VELOCITY;
    }

    float cell_density =
        (f2 + f6 + f8 + 2.f * (f3 + f4 + f5)) / (1.f - Cell::velocity[i][0].x);

    Cell::density[i][0] = cell_density;
  }

  // compute the equilibrium velocities
  for (std::size_t i{0}; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j++) {
      if (Cell::blockade[i][j]) {
        continue;
      }
      // calculate the equilibrium distribution
      for (std::size_t z{0}; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        float velocity_q = (DISC_VELOCITY[z].dot_product(Cell::velocity[i][j]));

        float linear_term = (3.f * velocity_q);

        float quadratic_term = linear_term * linear_term / 2.f;

        float last_term =
            (3.f * Cell::velocity[i][j].dot_product(Cell::velocity[i][j])) /
            (2.f);

        Cell::pdf_eq[k][i][j][z] =
            WEIGHTS[z] * Cell::density[i][j] *
            (1.f + linear_term + quadratic_term - last_term);
      }
    }
  }

  // apply Zou He boundary condition for the inlet flow
  // https://arxiv.org/abs/comp-gas/9611001
  for (std::size_t i{0}; i < LBM_CONSTANTS::HEIGHT; i++) {
    Cell::pdf[k][i][0][0] = Cell::pdf_eq[k][i][0][0];
    Cell::pdf[k][i][0][1] = Cell::pdf_eq[k][i][0][1];
    Cell::pdf[k][i][0][7] = Cell::pdf_eq[k][i][0][7];
  }

  // collide
  for (std::size_t i{0}; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j++) {
      if (Cell::blockade[i][j]) {
        continue;
      }
      for (std::size_t z{0}; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        // calculate collisions
        float relaxation_term =
            (Cell::pdf[k][i][j][z] - Cell::pdf_eq[k][i][j][z]) *
            LBM_CONSTANTS::RELAXATION_OMEGA;

        Vect<int32_t> new_pos = Vect{int32_t(j + LBM_CONSTANTS::WIDTH),
                                     int32_t(i + LBM_CONSTANTS::HEIGHT)};

        new_pos = new_pos + DISC_VELOCITY_INT[z];
        new_pos.y = new_pos.y % LBM_CONSTANTS::HEIGHT;
        new_pos.x = new_pos.x % LBM_CONSTANTS::WIDTH;

        if (Cell::blockade[new_pos.y][new_pos.x]) {
          // (out of bounds)
          // apply the half-way bounce back technique by reversing the
          // discrete velocity direction
          Cell::pdf[k ^ 1][i][j][(z + 4) % 8] +=
              Cell::pdf[k][i][j][z] - relaxation_term;
        } else {
          Cell::pdf[k ^ 1][new_pos.y][new_pos.x][z] +=
              Cell::pdf[k][i][j][z] - relaxation_term;
        }
      }
    }
  }
}
