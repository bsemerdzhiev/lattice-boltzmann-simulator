#include "lbm.hpp"
#include "cell.hpp"
#include "lbm_c.hpp"
#include "math_util.hpp"
#include <cstddef>
#include <immintrin.h>

void LBM::initialize() {
  for (std::size_t i{0}; i < LBM_CONSTANTS::HEIGHT; i++) {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j++) {

      Cell::density[i][j] = 1.0f;
      Cell::velocity_x[i][j] = LBM_CONSTANTS::HORIZONTAL_VELOCITY.x;
      Cell::velocity_y[i][j] = LBM_CONSTANTS::HORIZONTAL_VELOCITY.y;

      for (std::size_t z = 0; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        Cell::pdf[0][i][z][j] = WEIGHTS[z];
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
        float velocity_q = (DISC_VELOCITY[z].dot_product(
            Vect<float>{Cell::velocity_x[i][j], Cell::velocity_y[i][j]}));

        float linear_term = (3.f * velocity_q);

        float quadratic_term = linear_term * linear_term / 2.f;

        float last_term =
            (3.f * Vect<float>{Cell::velocity_x[i][j], Cell::velocity_y[i][j]}
                       .dot_product(Vect<float>{Cell::velocity_x[i][j],
                                                Cell::velocity_y[i][j]})) /
            (2.f);

        Cell::pdf[0][i][z][j] =
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

void LBM::update(bool k, std::size_t from, std::size_t to,
                 std::barrier<> &sync_barrier) {
  sync_barrier.arrive_and_wait();
  __m256 half = _mm256_set1_ps(0.5f);
  __m256 zero = _mm256_setzero_ps();
  __m256 one = _mm256_set1_ps(1.f);
  __m256 two = _mm256_set1_ps(2.f);
  __m256 three = _mm256_set1_ps(3.f);
  __m256 relaxation_constant = _mm256_set1_ps(LBM_CONSTANTS::RELAXATION_OMEGA);

  // set the pdf of the next step to zero
  for (std::size_t i{from}; i < to; i++) {
    for (std::size_t z{0}; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
      for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j += 8) {
        _mm256_store_ps(&Cell::pdf[k ^ 1][i][z][j], zero);
        _mm256_store_ps(&Cell::pdf_eq[k ^ 1][i][z][j], zero);
      }
    }
  }

  // right boundary outlet
  for (std::size_t i{from}; i < to; i++) {
    // for the outlet, perform interpolation from the results of the
    // previous cell
    Cell::pdf[k][i][3][LBM_CONSTANTS::WIDTH - 1] =
        Cell::pdf[k][i][3][LBM_CONSTANTS::WIDTH - 2];
    Cell::pdf[k][i][4][LBM_CONSTANTS::WIDTH - 1] =
        Cell::pdf[k][i][4][LBM_CONSTANTS::WIDTH - 2];
    Cell::pdf[k][i][5][LBM_CONSTANTS::WIDTH - 1] =
        Cell::pdf[k][i][5][LBM_CONSTANTS::WIDTH - 2];
  }

  // compute the macroscopic density and velocity
  for (std::size_t i{from}; i < to; i++) {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j += 8) {

      __m256 density_arr = zero;
      __m256 velocity_x = zero;
      __m256 velocity_y = zero;

      for (std::size_t z{0}; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        __m256 res = _mm256_load_ps(&Cell::pdf[k][i][z][j]);
        density_arr = _mm256_add_ps(density_arr, res);

        __m256 disc_arr_x = _mm256_set1_ps(DISC_VELOCITY[z].x);
        __m256 disc_arr_y = _mm256_set1_ps(DISC_VELOCITY[z].y);

        disc_arr_x = _mm256_mul_ps(disc_arr_x, res);
        disc_arr_y = _mm256_mul_ps(disc_arr_y, res);

        velocity_x = _mm256_add_ps(velocity_x, disc_arr_x);
        velocity_y = _mm256_add_ps(velocity_y, disc_arr_y);
      }
      velocity_x = _mm256_div_ps(velocity_x, density_arr);
      velocity_y = _mm256_div_ps(velocity_y, density_arr);
      _mm256_store_ps(&Cell::velocity_x[i][j], velocity_x);
      _mm256_store_ps(&Cell::velocity_y[i][j], velocity_y);

      _mm256_store_ps(&Cell::density[i][j], density_arr);
    }
  }

  // prescribe the Zou/He scheme
  for (std::size_t i{from}; i < to; i++) {
    float f2 = Cell::pdf[k][i][2][0];
    float f3 = Cell::pdf[k][i][3][0];
    float f4 = Cell::pdf[k][i][4][0];
    float f5 = Cell::pdf[k][i][5][0];
    float f6 = Cell::pdf[k][i][6][0];
    float f8 = Cell::pdf[k][i][8][0];

    // skip settings the velocity for
    // the first and last elements of the first column
    if ((i > 0) & (i < LBM_CONSTANTS::HEIGHT - 1)) {
      Cell::velocity_x[i][0] = LBM_CONSTANTS::HORIZONTAL_VELOCITY.x;
      Cell::velocity_y[i][0] = LBM_CONSTANTS::HORIZONTAL_VELOCITY.y;
    }

    float cell_density =
        (f2 + f6 + f8 + 2.f * (f3 + f4 + f5)) / (1.f - Cell::velocity_x[i][0]);

    Cell::density[i][0] = cell_density;
  }

  // compute the equilibrium velocities
  for (std::size_t i{from}; i < to; i++) {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j += 8) {
      __m256 velocity_x = _mm256_load_ps(&Cell::velocity_x[i][j]);
      __m256 velocity_y = _mm256_load_ps(&Cell::velocity_y[i][j]);

      __m256 density = _mm256_load_ps(&Cell::density[i][j]);

      // calculate the equilibrium distribution
      for (std::size_t z{0}; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        __m256 disc_velocity_x = _mm256_set1_ps(DISC_VELOCITY[z].x);
        __m256 disc_velocity_y = _mm256_set1_ps(DISC_VELOCITY[z].y);

        __m256 velocity_q_x = _mm256_mul_ps(disc_velocity_x, velocity_x);

        __m256 dot_product_result =
            _mm256_fmadd_ps(disc_velocity_y, velocity_y, velocity_q_x);

        __m256 linear_term = _mm256_mul_ps(dot_product_result, three);

        __m256 quadratic_term = _mm256_mul_ps(linear_term, linear_term);
        quadratic_term = _mm256_mul_ps(quadratic_term, half);

        __m256 magnitude = _mm256_mul_ps(velocity_x, velocity_x);
        magnitude = _mm256_fmadd_ps(velocity_y, velocity_y, magnitude);

        __m256 last_term = _mm256_mul_ps(three, magnitude);
        last_term = _mm256_mul_ps(last_term, half);

        __m256 weights_vect = _mm256_set1_ps(WEIGHTS[z]);

        __m256 final_result = _mm256_mul_ps(weights_vect, density);

        __m256 internal_summation = _mm256_add_ps(linear_term, quadratic_term);
        internal_summation = _mm256_sub_ps(internal_summation, last_term);
        internal_summation = _mm256_add_ps(internal_summation, one);

        final_result = _mm256_mul_ps(final_result, internal_summation);

        _mm256_store_ps(&Cell::pdf_eq[k][i][z][j], final_result);
      }
    }
  }

  // apply Zou He boundary condition for the inlet flow
  // https://arxiv.org/abs/comp-gas/9611001
  for (std::size_t i{from}; i < to; i++) {
    Cell::pdf[k][i][0][0] = Cell::pdf_eq[k][i][0][0];
    Cell::pdf[k][i][1][0] = Cell::pdf_eq[k][i][1][0];
    Cell::pdf[k][i][7][0] = Cell::pdf_eq[k][i][7][0];
  }

  __m256i width_vec = _mm256_set1_epi32(LBM_CONSTANTS::WIDTH);
  __m256i height_vec = _mm256_set1_epi32(LBM_CONSTANTS::HEIGHT);

  // collide (Push stage)
  for (std::size_t i{from}; i < to; i++) {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j += 8) {

      for (std::size_t z{0}; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        // calculate collisions
        // float relaxation_term =
        __m256 cur_pdf = _mm256_load_ps(&Cell::pdf[k][i][z][j]);
        __m256 cur_pdf_eq = _mm256_load_ps(&Cell::pdf_eq[k][i][z][j]);
        __m256 relaxation_term = _mm256_sub_ps(cur_pdf, cur_pdf_eq);
        relaxation_term = _mm256_mul_ps(relaxation_term, relaxation_constant);

        __m256 final_result = _mm256_sub_ps(cur_pdf, relaxation_term);

        _mm256_store_ps(&Cell::pdf[k][i][z][j], final_result);
      }
    }
  }
  sync_barrier.arrive_and_wait();
  for (std::size_t i{from}; i < to; i++) {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j += 8) {
      __m256i pos_x =
          _mm256_add_epi32(_mm256_set1_epi32(j + LBM_CONSTANTS::WIDTH),
                           _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7));
      __m256i pos_y = _mm256_set1_epi32(i + LBM_CONSTANTS::HEIGHT);

      for (std::size_t z{0}; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        __m256i disc_velocity_x = _mm256_set1_epi32(DISC_VELOCITY_INT[z].x);
        __m256i disc_velocity_y = _mm256_set1_epi32(DISC_VELOCITY_INT[z].y);

        disc_velocity_x = _mm256_sub_epi32(pos_x, disc_velocity_x);

        __m256i mask = _mm256_cmpgt_epi32(
            disc_velocity_x, _mm256_set1_epi32(LBM_CONSTANTS::WIDTH - 1));
        disc_velocity_x = _mm256_sub_epi32(disc_velocity_x,
                                           _mm256_and_si256(mask, width_vec));
        mask = _mm256_cmpgt_epi32(disc_velocity_x,
                                  _mm256_set1_epi32(LBM_CONSTANTS::WIDTH - 1));
        disc_velocity_x = _mm256_sub_epi32(disc_velocity_x,
                                           _mm256_and_si256(mask, width_vec));

        disc_velocity_y = _mm256_sub_epi32(pos_y, disc_velocity_y);
        mask = _mm256_cmpgt_epi32(disc_velocity_y,
                                  _mm256_set1_epi32(LBM_CONSTANTS::HEIGHT - 1));
        disc_velocity_y = _mm256_sub_epi32(disc_velocity_y,
                                           _mm256_and_si256(mask, height_vec));
        mask = _mm256_cmpgt_epi32(disc_velocity_y,
                                  _mm256_set1_epi32(LBM_CONSTANTS::HEIGHT - 1));
        disc_velocity_y = _mm256_sub_epi32(disc_velocity_y,
                                           _mm256_and_si256(mask, height_vec));

        __m256i new_ind = _mm256_mullo_epi32(disc_velocity_y, width_vec);
        new_ind = _mm256_add_epi32(new_ind, disc_velocity_x);

        __m256i blockade_vals =
            _mm256_i32gather_epi32(&Cell::blockade[0][0], new_ind, 4);

        __m256 bounce_back = _mm256_load_ps(&Cell::pdf[k][i][(z + 4) % 8][j]);

        __m256i base_offset = _mm256_set1_epi32(
            k * LBM_CONSTANTS::HEIGHT * LBM_CONSTANTS::LATTICE_COUNT *
                LBM_CONSTANTS::WIDTH +
            z * LBM_CONSTANTS::WIDTH);

        __m256i source_idx = _mm256_add_epi32(
            _mm256_mullo_epi32(disc_velocity_y,
                               _mm256_set1_epi32(LBM_CONSTANTS::LATTICE_COUNT *
                                                 LBM_CONSTANTS::WIDTH)),

            _mm256_add_epi32(base_offset, disc_velocity_x));
        __m256 pulled_pdf =
            _mm256_i32gather_ps(&Cell::pdf[0][0][0][0], source_idx, 4);
        __m256 blockade_mask = _mm256_castsi256_ps(_mm256_andnot_si256(
            _mm256_cmpeq_epi32(blockade_vals, _mm256_setzero_si256()),
            _mm256_set1_epi32(0xFFFFFFFF)));

        __m256 result =
            _mm256_blendv_ps(pulled_pdf, bounce_back, blockade_mask);

        __m256i is_current_cell_blocked =
            _mm256_load_epi32(&Cell::blockade[i][j]);
        __m256 blocked_mask = _mm256_castsi256_ps(_mm256_cmpeq_epi32(
            is_current_cell_blocked, _mm256_setzero_si256()));

        __m256 existing = _mm256_load_ps(&Cell::pdf[k ^ 1][i][z][j]);

        __m256 final = _mm256_blendv_ps(existing, result, blocked_mask);

        _mm256_store_ps(&Cell::pdf[k ^ 1][i][z][j], final);
      }
    }
  }

  sync_barrier.arrive_and_wait();
}
