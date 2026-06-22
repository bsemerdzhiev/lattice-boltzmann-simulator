#include "lbm.hpp"
#include "cell.hpp"
#include "lbm_c.hpp"
#include "math_util.hpp"
#include <cstddef>
#include <immintrin.h>

constexpr bool KARMAN_VORTEX = true;

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
  // for (std::size_t i{0}; i < LBM_CONSTANTS::WIDTH; i++) {
  //   Cell::blockade[0][i] = 1;
  //   Cell::blockade[LBM_CONSTANTS::HEIGHT - 1][i] = 1;
  // }

  if constexpr (KARMAN_VORTEX) {
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
  } else {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j++) {
      Cell::blockade[0][j] = 1;
      Cell::blockade[LBM_CONSTANTS::HEIGHT - 1][j] = 1;
    }
  }
}

void LBM::update(bool k, std::size_t from, std::size_t to,
                 std::barrier<> &sync_barrier) {
  sync_barrier.arrive_and_wait();
  __m512 half = _mm512_set1_ps(0.5f);
  __m512 zero = _mm512_setzero_ps();
  __m512 one = _mm512_set1_ps(1.f);
  __m512 two = _mm512_set1_ps(2.f);
  __m512 three = _mm512_set1_ps(3.f);
  __m512 relaxation_constant = _mm512_set1_ps(LBM_CONSTANTS::RELAXATION_OMEGA);

  // set the pdf of the next step to zero
  for (std::size_t i{from}; i < to; i++) {
    for (std::size_t z{0}; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
      for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j += 16) {
        _mm512_store_ps(&Cell::pdf[k ^ 1][i][z][j], zero);
        _mm512_store_ps(&Cell::pdf_eq[k ^ 1][i][z][j], zero);
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
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j += 16) {

      __m512 density_arr = zero;
      __m512 velocity_x = zero;
      __m512 velocity_y = zero;

      for (std::size_t z{0}; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        __m512 res = _mm512_load_ps(&Cell::pdf[k][i][z][j]);
        density_arr = _mm512_add_ps(density_arr, res);

        __m512 disc_arr_x = _mm512_set1_ps(DISC_VELOCITY[z].x);
        __m512 disc_arr_y = _mm512_set1_ps(DISC_VELOCITY[z].y);

        disc_arr_x = _mm512_mul_ps(disc_arr_x, res);
        disc_arr_y = _mm512_mul_ps(disc_arr_y, res);

        velocity_x = _mm512_add_ps(velocity_x, disc_arr_x);
        velocity_y = _mm512_add_ps(velocity_y, disc_arr_y);
      }
      velocity_x = _mm512_div_ps(velocity_x, density_arr);
      velocity_y = _mm512_div_ps(velocity_y, density_arr);

      _mm512_store_ps(&Cell::velocity_x[i][j], velocity_x);
      _mm512_store_ps(&Cell::velocity_y[i][j], velocity_y);

      _mm512_store_ps(&Cell::density[i][j], density_arr);
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
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j += 16) {
      __m512 velocity_x = _mm512_load_ps(&Cell::velocity_x[i][j]);
      __m512 velocity_y = _mm512_load_ps(&Cell::velocity_y[i][j]);

      __m512 density = _mm512_load_ps(&Cell::density[i][j]);

      // calculate the equilibrium distribution
      for (std::size_t z{0}; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        __m512 disc_velocity_x = _mm512_set1_ps(DISC_VELOCITY[z].x);
        __m512 disc_velocity_y = _mm512_set1_ps(DISC_VELOCITY[z].y);

        __m512 velocity_q_x = _mm512_mul_ps(disc_velocity_x, velocity_x);

        __m512 dot_product_result =
            _mm512_fmadd_ps(disc_velocity_y, velocity_y, velocity_q_x);

        __m512 linear_term = _mm512_mul_ps(dot_product_result, three);

        __m512 quadratic_term = _mm512_mul_ps(linear_term, linear_term);
        quadratic_term = _mm512_mul_ps(quadratic_term, half);

        __m512 magnitude = _mm512_mul_ps(velocity_x, velocity_x);
        magnitude = _mm512_fmadd_ps(velocity_y, velocity_y, magnitude);

        __m512 last_term = _mm512_mul_ps(three, magnitude);
        last_term = _mm512_mul_ps(last_term, half);

        __m512 weights_vect = _mm512_set1_ps(WEIGHTS[z]);

        __m512 final_result = _mm512_mul_ps(weights_vect, density);

        __m512 internal_summation = _mm512_add_ps(linear_term, quadratic_term);
        internal_summation = _mm512_sub_ps(internal_summation, last_term);
        internal_summation = _mm512_add_ps(internal_summation, one);

        final_result = _mm512_mul_ps(final_result, internal_summation);

        _mm512_store_ps(&Cell::pdf_eq[k][i][z][j], final_result);
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

  __m512i width_vec = _mm512_set1_epi32(LBM_CONSTANTS::WIDTH);
  __m512i height_vec = _mm512_set1_epi32(LBM_CONSTANTS::HEIGHT);

  // collide (Push stage)
  for (std::size_t i{from}; i < to; i++) {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j += 16) {

      for (std::size_t z{0}; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        // calculate collisions
        // float relaxation_term =
        __m512 cur_pdf = _mm512_load_ps(&Cell::pdf[k][i][z][j]);
        __m512 cur_pdf_eq = _mm512_load_ps(&Cell::pdf_eq[k][i][z][j]);
        __m512 relaxation_term = _mm512_sub_ps(cur_pdf, cur_pdf_eq);
        relaxation_term = _mm512_mul_ps(relaxation_term, relaxation_constant);

        __m512 final_result = _mm512_sub_ps(cur_pdf, relaxation_term);

        _mm512_store_ps(&Cell::pdf[k][i][z][j], final_result);
      }
    }
  }
  sync_barrier.arrive_and_wait();
  for (std::size_t i{from}; i < to; i++) {
    for (std::size_t j{0}; j < LBM_CONSTANTS::WIDTH; j += 16) {
      __m512i pos_x =
          _mm512_add_epi32(_mm512_set1_epi32(j + LBM_CONSTANTS::WIDTH),
                           _mm512_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                             11, 12, 13, 14, 15));
      __m512i pos_y = _mm512_set1_epi32(i + LBM_CONSTANTS::HEIGHT);

      for (std::size_t z{0}; z < LBM_CONSTANTS::LATTICE_COUNT; z++) {
        __m512i disc_velocity_x = _mm512_set1_epi32(DISC_VELOCITY_INT[z].x);
        __m512i disc_velocity_y = _mm512_set1_epi32(DISC_VELOCITY_INT[z].y);

        disc_velocity_x = _mm512_sub_epi32(pos_x, disc_velocity_x);

        __mmask16 mask = _mm512_cmpgt_epi32_mask(
            disc_velocity_x, _mm512_set1_epi32(LBM_CONSTANTS::WIDTH - 1));

        disc_velocity_x = _mm512_mask_sub_epi32(disc_velocity_x, mask,
                                                disc_velocity_x, width_vec);

        mask = _mm512_cmpgt_epi32_mask(
            disc_velocity_x, _mm512_set1_epi32(LBM_CONSTANTS::WIDTH - 1));

        disc_velocity_x = _mm512_mask_sub_epi32(disc_velocity_x, mask,
                                                disc_velocity_x, width_vec);

        disc_velocity_y = _mm512_sub_epi32(pos_y, disc_velocity_y);

        mask = _mm512_cmpgt_epi32_mask(
            disc_velocity_y, _mm512_set1_epi32(LBM_CONSTANTS::HEIGHT - 1));

        disc_velocity_y = _mm512_mask_sub_epi32(disc_velocity_y, mask,
                                                disc_velocity_y, height_vec);
        mask = _mm512_cmpgt_epi32_mask(
            disc_velocity_y, _mm512_set1_epi32(LBM_CONSTANTS::HEIGHT - 1));

        disc_velocity_y = _mm512_mask_sub_epi32(disc_velocity_y, mask,
                                                disc_velocity_y, height_vec);

        __m512i new_ind = _mm512_mullo_epi32(disc_velocity_y, width_vec);
        new_ind = _mm512_add_epi32(new_ind, disc_velocity_x);

        __m512i blockade_vals =
            _mm512_i32gather_epi32(new_ind, &Cell::blockade[0][0], 4);

        __m512 bounce_back = _mm512_load_ps(&Cell::pdf[k][i][(z + 4) % 8][j]);

        __m512i base_offset = _mm512_set1_epi32(
            k * LBM_CONSTANTS::HEIGHT * LBM_CONSTANTS::LATTICE_COUNT *
                LBM_CONSTANTS::WIDTH +
            z * LBM_CONSTANTS::WIDTH);

        __m512i source_idx = _mm512_add_epi32(
            _mm512_mullo_epi32(disc_velocity_y,
                               _mm512_set1_epi32(LBM_CONSTANTS::LATTICE_COUNT *
                                                 LBM_CONSTANTS::WIDTH)),

            _mm512_add_epi32(base_offset, disc_velocity_x));
        __m512 pulled_pdf =
            _mm512_i32gather_ps(source_idx, &Cell::pdf[0][0][0][0], 4);

        __mmask16 blockade_mask =
            _mm512_cmpeq_epi32_mask(blockade_vals, _mm512_setzero_si512());

        __m512 result =
            _mm512_mask_blend_ps(blockade_mask, bounce_back, pulled_pdf);

        __m512i is_current_cell_blocked =
            _mm512_load_epi32(&Cell::blockade[i][j]);

        __mmask16 blocked_mask = _mm512_cmpeq_epi32_mask(
            is_current_cell_blocked, _mm512_setzero_si512());

        __m512 existing = _mm512_load_ps(&Cell::pdf[k ^ 1][i][z][j]);

        __m512 final = _mm512_mask_blend_ps(blocked_mask, existing, result);

        _mm512_store_ps(&Cell::pdf[k ^ 1][i][z][j], final);
      }
    }
  }

  sync_barrier.arrive_and_wait();
}
