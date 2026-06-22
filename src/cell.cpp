#include <cell.hpp>

alignas(64) Lattice<float> Cell::density;

alignas(64) Lattice<float> Cell::velocity_x;
alignas(64) Lattice<float> Cell::velocity_y;

alignas(64) Lattice<int32_t> Cell::blockade;

alignas(64) PdfLattice<float> Cell::pdf[2];
alignas(64) PdfLattice<float> Cell::pdf_eq[2];
