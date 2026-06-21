#include <cell.hpp>

Lattice<float> Cell::density;
Lattice<Vect<float>> Cell::velocity;

Lattice<bool> Cell::blockade;

Lattice<std::array<float, LBM_CONSTANTS::LATTICE_COUNT>> Cell::pdf[2];
Lattice<std::array<float, LBM_CONSTANTS::LATTICE_COUNT>> Cell::pdf_eq[2];
