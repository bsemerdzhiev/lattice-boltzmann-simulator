#include "lbm.hpp"
#include "lbm_c.hpp"
#include "math_util.hpp"
#include <cstddef>

LBM::Lattice<Cell> LBM::lattice;
LBM::Lattice<bool> LBM::blockade;

void LBM::initialize() {
  for (std::size_t i = 0; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j = 0; j < LBM_CONSTANTS::WIDTH; j++) {
      lattice[i][j].density = 1.00;
      lattice[i][j].velocity = LBM_CONSTANTS::HORIZONTAL_VELOCITY;

      for (std::size_t z = 0; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        lattice[i][j].pdf[0][z] = WEIGHTS[z];
      }
    }
  }

  // compute the equilibrium velocities
  for (std::size_t i = 0; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j = 0; j < LBM_CONSTANTS::WIDTH; j++) {
      if (blockade[i][j]) {
        continue;
      }
      // calculate the equilibrium distribution
      for (std::size_t z = 0; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        double velocity_q =
            (DISC_VELOCITY[z].dot_product(lattice[i][j].velocity));

        double linear_term = (3.0 * velocity_q);

        double quadratic_term = linear_term * linear_term / 2.0;

        double last_term =
            (3.0 * lattice[i][j].velocity.dot_product(lattice[i][j].velocity)) /
            (2.0);

        lattice[i][j].pdf[0][z] =
            WEIGHTS[z] * lattice[i][j].density *
            (1.0 + linear_term + quadratic_term - last_term);
      }
    }
  }

  Vect<double> circle_coord{LBM_CONSTANTS::WIDTH / 5.0,
                            LBM_CONSTANTS::HEIGHT / 2.0};

  for (std::size_t i = 0; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j = 0; j < LBM_CONSTANTS::WIDTH; j++) {
      if (circle_coord.euclid_dist(Vect<double>{1.0 * j, 1.0 * i}) <
          LBM_CONSTANTS::CYLINDER_RADIUS) {
        blockade[i][j] = 1;
      }
    }
  }
}

void LBM::update(bool k) {
  // set the pdf of the next step to zero
  for (std::size_t i{0}; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j++) {
      for (std::size_t z{0}; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        lattice[i][j].pdf[k ^ 1][z] = 0.0;
        lattice[i][j].pdf_eq[k ^ 1][z] = 0.0;
      }
    }
  }

  // right boundary outlet
  for (std::size_t i = 0; i < LBM_CONSTANTS::HEIGHT; i++) {
    // for the outlet, perform interpolation from the results of the
    // previous cell
    lattice[i][LBM_CONSTANTS::WIDTH - 1].pdf[k][3] =
        lattice[i][LBM_CONSTANTS::WIDTH - 2].pdf[k][3];
    lattice[i][LBM_CONSTANTS::WIDTH - 1].pdf[k][4] =
        lattice[i][LBM_CONSTANTS::WIDTH - 2].pdf[k][4];
    lattice[i][LBM_CONSTANTS::WIDTH - 1].pdf[k][5] =
        lattice[i][LBM_CONSTANTS::WIDTH - 2].pdf[k][5];
  }

  // compute the macroscopic density and velocity
  for (std::size_t i = 0; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j = 0; j < LBM_CONSTANTS::WIDTH; j++) {
      if (blockade[i][j]) {
        continue;
      }
      lattice[i][j].density = 0.0;
      lattice[i][j].velocity = Vect{0.0, 0.0};
      for (std::size_t z = 0; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        lattice[i][j].density += lattice[i][j].pdf[k][z];
        lattice[i][j].velocity =
            lattice[i][j].velocity + DISC_VELOCITY[z] * lattice[i][j].pdf[k][z];
      }
      lattice[i][j].velocity = lattice[i][j].velocity / lattice[i][j].density;
    }
  }

  // prescribe the Zou/He scheme
  for (std::size_t i{0}; i < LBM_CONSTANTS::HEIGHT; i++) {
    double f2 = lattice[i][0].pdf[k][2];
    double f3 = lattice[i][0].pdf[k][3];
    double f4 = lattice[i][0].pdf[k][4];
    double f5 = lattice[i][0].pdf[k][5];
    double f6 = lattice[i][0].pdf[k][6];
    double f8 = lattice[i][0].pdf[k][8];

    // skip settings the velocity for
    // the first and last elements of the first column
    if (i > 0 && i < LBM_CONSTANTS::HEIGHT - 1) {
      lattice[i][0].velocity = LBM_CONSTANTS::HORIZONTAL_VELOCITY;
    }

    double cell_density = (f2 + f6 + f8 + 2.0 * (f3 + f4 + f5)) /
                          (1.0 - lattice[i][0].velocity.x);

    lattice[i][0].density = cell_density;
  }

  // compute the equilibrium velocities
  for (std::size_t i = 0; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j = 0; j < LBM_CONSTANTS::WIDTH; j++) {
      if (blockade[i][j]) {
        continue;
      }
      // calculate the equilibrium distribution
      for (std::size_t z = 0; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        double velocity_q =
            (DISC_VELOCITY[z].dot_product(lattice[i][j].velocity));

        double linear_term = (3.0 * velocity_q);

        double quadratic_term = linear_term * linear_term / 2.0;

        double last_term =
            (3.0 * lattice[i][j].velocity.dot_product(lattice[i][j].velocity)) /
            (2.0);

        lattice[i][j].pdf_eq[k][z] =
            WEIGHTS[z] * lattice[i][j].density *
            (1.0 + linear_term + quadratic_term - last_term);
      }
    }
  }

  // apply Zou He boundary condition for the inlet flow
  // https://arxiv.org/abs/comp-gas/9611001
  for (std::size_t i = 0; i < LBM_CONSTANTS::HEIGHT; i++) {
    lattice[i][0].pdf[k][0] = lattice[i][0].pdf_eq[k][0];
    lattice[i][0].pdf[k][1] = lattice[i][0].pdf_eq[k][1];
    lattice[i][0].pdf[k][7] = lattice[i][0].pdf_eq[k][7];
  }

  // collide
  for (std::size_t i{0}; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j++) {
      if (blockade[i][j]) {
        continue;
      }
      for (std::size_t z = 0; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        // calculate collisions
        double relaxation_term =
            (lattice[i][j].pdf[k][z] - lattice[i][j].pdf_eq[k][z]) *
            LBM_CONSTANTS::RELAXATION_OMEGA;

        Vect<int32_t> new_pos = Vect{int32_t(j + LBM_CONSTANTS::WIDTH),
                                     int32_t(i + LBM_CONSTANTS::HEIGHT)};

        new_pos = new_pos + DISC_VELOCITY_INT[z];
        new_pos.y = new_pos.y % LBM_CONSTANTS::HEIGHT;
        new_pos.x = new_pos.x % LBM_CONSTANTS::WIDTH;

        if (blockade[new_pos.y][new_pos.x]) {
          // (out of bounds)
          // apply the half-way bounce back technique by reversing the
          // discrete velocity direction
          lattice[i][j].pdf[k ^ 1][(z + 4) % 8] +=
              lattice[i][j].pdf[k][z] - relaxation_term;
        } else {
          lattice[new_pos.y][new_pos.x].pdf[k ^ 1][z] +=
              lattice[i][j].pdf[k][z] - relaxation_term;
        }
      }
    }
  }
}
