#include <cell.hpp>

Lattice<float> Cell::density;
// Lattice<Vect<float>> velocity;

Lattice<float> Cell::velocity_x;
Lattice<float> Cell::velocity_y;

Lattice<int32_t> Cell::blockade;

PdfLattice<float> Cell::pdf[2];
PdfLattice<float> Cell::pdf_eq[2];
