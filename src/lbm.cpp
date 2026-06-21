#include "lbm.hpp"
#include "lbm_c.hpp"
#include "math_util.hpp"
#include <cstddef>

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

  Vect<float> circle_coord{100, HEIGHT / 2.f};

  for (std::size_t i = 0; i < HEIGHT; i++) {
    for (std::size_t j = 0; j < WIDTH; j++) {
      if (circle_coord.euclid_dist(Vect<float>{1.f * j, 1.f * i}) < 1000) {
        blockade[i][j] = 1;
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

        Vect<int32_t> new_coord =
            Vect{static_cast<int32_t>(j), static_cast<int32_t>(i)} +
            DISC_VELOCITY_INT[z];

        if (new_coord.x < 0 || new_coord.y < 0 || new_coord.x >= WIDTH ||
            new_coord.y >= HEIGHT) {
          continue;
        }

        if (blockade[new_coord.y][new_coord.x]) {
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
  for (std::size_t i = 0; i < HEIGHT; i++) {
    Vect<float> velocity = wanted_velocity;

    float f2 = lattice[i][0].pdf[k ^ 1][2];
    float f3 = lattice[i][0].pdf[k ^ 1][3];
    float f4 = lattice[i][0].pdf[k ^ 1][4];
    float f5 = lattice[i][0].pdf[k ^ 1][5];
    float f6 = lattice[i][0].pdf[k ^ 1][6];
    float f8 = lattice[i][0].pdf[k ^ 1][8];

    float cell_density =
        (f2 + f6 + f8 + 2.f * (f3 + f4 + f5)) / (1.f - velocity.x);

    lattice[i][0].pdf[k ^ 1][0] = f4 + (2.f * cell_density * velocity.x) / 3.f;
    lattice[i][0].pdf[k ^ 1][7] =
        f3 + (cell_density * velocity.x) / 6.f + (f2 - f6) / 2.f;
    lattice[i][0].pdf[k ^ 1][1] =
        f5 + (cell_density * velocity.x) / 6.f + (f6 - f2) / 2.f;

    lattice[i][WIDTH - 1].pdf[k ^ 1][3] = lattice[i][WIDTH - 2].pdf[k ^ 1][3];
    lattice[i][WIDTH - 1].pdf[k ^ 1][4] = lattice[i][WIDTH - 2].pdf[k ^ 1][4];
    lattice[i][WIDTH - 1].pdf[k ^ 1][5] = lattice[i][WIDTH - 2].pdf[k ^ 1][5];
  }
}
