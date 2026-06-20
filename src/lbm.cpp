#include "lbm.hpp"
#include "lbm_c.hpp"
#include "math_util.hpp"
#include <cstddef>
#include <format>
#include <iostream>

LBM::Lattice<Cell> LBM::lattice;
LBM::Lattice<bool> LBM::blockade;
LBM::Lattice<Vect<float>> LBM::forces;

void LBM::initialize() {
  for (std::size_t i = 0; i < HEIGHT; i++) {
    for (std::size_t j = 0; j < WIDTH; j++) {
      lattice[i][j].density = 1.0f;
      lattice[i][j].velocity = Vect{0.f, 0.f};

      for (std::size_t z = 0; z < LATTICE_COUNT; z++) {
        lattice[i][j].pdf[0][z] = WEIGHTS[z];
        lattice[i][j].pdf[1][z] = WEIGHTS[z];
      }
    }
  }

  // initialize the blockade array
  for (std::size_t i = 0; i < WIDTH; i++) {
    blockade[0][i] = true;
    blockade[HEIGHT - 1][i] = true;
  }
  for (std::size_t i = 0; i < HEIGHT; i++) {
    blockade[i][0] = true;
    blockade[i][WIDTH - 1] = true;

    if (i >= HEIGHT / 4.f && i <= HEIGHT * 3.f / 4.f)
      forces[i][1] = Vect{1.f, 0.f};
  }

  Vect<float> circle_coord{WIDTH / 10.f, HEIGHT / 2.f};

  for (std::size_t i = 0; i < HEIGHT; i++) {
    for (std::size_t j = 0; j < WIDTH; j++) {
      if (circle_coord.euclid_dist(Vect<float>{1.f * j, 1.f * i}) < 40) {
        blockade[i][j] = 1;
      }
    }
    // std::cout << "\n";
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
      if (blockade[i][j])
        continue;

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

        float body_force_term =
            WEIGHTS[z] *
            (DISC_VELOCITY[z].dot_product(forces[i][j]) * delta_t) / C_s_sq;

        Vect<int32_t> new_coord =
            Vect{static_cast<int32_t>(j), static_cast<int32_t>(i)} +
            DISC_VELOCITY_INT[z];

        if (blockade[new_coord.y][new_coord.x]) {
          // (out of bounds)
          // apply bounce back technique by reversing the
          // discrete velocity direction

          lattice[i][j].pdf[k ^ 1][(z + 4) % 8] +=
              lattice[i][j].pdf[k][z] - relaxation_term;
        } else {
          lattice[new_coord.y][new_coord.x].pdf[k ^ 1][z] +=
              lattice[i][j].pdf[k][z] - relaxation_term + body_force_term;
        }
      }
    }
  }
}
