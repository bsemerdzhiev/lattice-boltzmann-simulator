#include "lbm.hpp"
#include "lbm_c.hpp"
#include <cstddef>

std::array<std::array<Cell, HEIGHT>, WIDTH> lattice;

void LBM::initialize() {
  for (std::size_t i = 0; i < HEIGHT; i++) {
    for (std::size_t j = 0; j < WIDTH; j++) {
      lattice[i][j].density = 1.0f;
      lattice[i][j].velocity = Vect{0.f, 0.f};

      for (std::size_t z = 0; z < LATTICE_COUNT; z++) {
        lattice[i][j].pdf[z] = WEIGHTS[z];
      }
    }
  }
}

void LBM::update() {
  // compute the macroscopic density and velocity

  for (std::size_t i = 0; i < HEIGHT; i++) {
    for (std::size_t j = 0; j < WIDTH; j++) {
      // calculate the density and velocity for the current cell
      lattice[i][j].density = 0.f;
      lattice[i][j].velocity = Vect{0.f, 0.f};
      for (std::size_t z = 0; z < LATTICE_COUNT; z++) {
        lattice[i][j].density += lattice[i][j].pdf[z];
        lattice[i][j].velocity =
            lattice[i][j].velocity + DISC_VELOCITY[z] * lattice[i][j].pdf[z];
      }
      lattice[i][j].velocity = lattice[i][j].velocity / lattice[i][j].density;

      // calculate the equilibrium distribution
      for (std::size_t z = 0; z < LATTICE_COUNT; z++) {
        float velocity_q =
            (DISC_VELOCITY[z].dot_product(lattice[i][j].velocity));

        float linear_term = (3 * velocity_q) / lattice_speed_squred;

        float quadratic_term = linear_term * linear_term / 2.f;

        float last_term =
            (3.f * lattice[i][j].velocity.dot_product(lattice[i][j].velocity)) /
            (2.f * lattice_speed_squred);

        lattice[i][j].pdf[z] = WEIGHTS[z] * lattice[i][j].density *
                               (1.f + linear_term + quadratic_term - last_term);
      }
    }
  }
}
