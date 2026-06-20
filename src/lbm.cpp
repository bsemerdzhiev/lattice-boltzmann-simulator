#include "lbm.hpp"
#include "lbm_c.hpp"
#include "math_util.hpp"
#include <cstddef>

std::array<std::array<Cell, WIDTH>, HEIGHT> LBM::lattice;

void LBM::initialize() {
  for (std::size_t i = 0; i < HEIGHT; i++) {
    for (std::size_t j = 0; j < WIDTH; j++) {
      lattice[i][j].density = 1.0f;
      lattice[i][j].velocity = Vect{0.f, 0.f};

      for (std::size_t z = 0; z < LATTICE_COUNT; z++) {
        lattice[i][j].pdf[0][z] = WEIGHTS[z];
      }
    }
  }
}

void LBM::update(bool k) {
  // set the pdf of the next step to zero
  for (std::size_t i{0}; i < HEIGHT; i++) {
    for (std::size_t j{0}; j < WIDTH; j++) {
      for (std::size_t z{0}; z < LATTICE_COUNT; z++) {
        lattice[i][j].pdf[k ^ 1][z] = 0.f;
      }
    }
  }

  // compute the macroscopic density and velocity
  for (std::size_t i = 0; i < HEIGHT; i++) {
    for (std::size_t j = 0; j < WIDTH; j++) {
      // calculate the density and velocity for the current cell
      lattice[i][j].density = 0.f;
      lattice[i][j].velocity = Vect{0.f, 0.f};
      for (std::size_t z = 0; z < LATTICE_COUNT; z++) {
        lattice[i][j].density += lattice[i][j].pdf[k][z];
        lattice[i][j].velocity =
            lattice[i][j].velocity + DISC_VELOCITY[z] * lattice[i][j].pdf[k][z];
      }
      lattice[i][j].velocity = lattice[i][j].velocity / lattice[i][j].density;

      // calculate the equilibrium distribution
      for (std::size_t z = 0; z < LATTICE_COUNT; z++) {
        float velocity_q =
            (DISC_VELOCITY[z].dot_product(lattice[i][j].velocity));

        float linear_term = (3.f * velocity_q) / lattice_speed_squred;

        float quadratic_term = linear_term * linear_term / 2.f;

        float last_term =
            (3.f * lattice[i][j].velocity.dot_product(lattice[i][j].velocity)) /
            (2.f * lattice_speed_squred);

        lattice[i][j].pdf_eq[k][z] =
            WEIGHTS[z] * lattice[i][j].density *
            (1.f + linear_term + quadratic_term - last_term);
      }

      // calculate collisions
      for (std::size_t z = 0; z < LATTICE_COUNT; z++) {
        float relaxation_term =
            (lattice[i][j].pdf[k][z] - lattice[i][j].pdf_eq[k][z]) / THAO;

        Vect<int32_t> new_coord =
            Vect{static_cast<int32_t>(j), static_cast<int32_t>(i)} +
            DISC_VELOCITY_INT[z];

        if (new_coord.x < 0 || new_coord.y < 0 || new_coord.y >= HEIGHT ||
            new_coord.x >= WIDTH) {
          // (out of bounds)
          // apply bounce back technique by reversing the
          // discrete velocity direction
          lattice[i][j].pdf[k ^ 1][(z + 4) % 8] +=
              lattice[i][j].pdf[k][z] - relaxation_term;
        } else {
          lattice[new_coord.y][new_coord.x].pdf[k ^ 1][z] +=
              lattice[i][j].pdf[k][z] - relaxation_term;
        }
      }
    }
  }
}
