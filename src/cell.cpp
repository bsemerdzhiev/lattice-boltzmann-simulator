#include <cell.hpp>

alignas(32) Lattice<float> Cell::density;
// Lattice<Vect<float>> velocity;

alignas(32) Lattice<float> Cell::velocity_x;
alignas(32) Lattice<float> Cell::velocity_y;

alignas(32) Lattice<int32_t> Cell::blockade;

alignas(32) PdfLattice<float> Cell::pdf[2];
alignas(32) PdfLattice<float> Cell::pdf_eq[2];
