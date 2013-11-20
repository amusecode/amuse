#include "aton_fortran.h"

extern "C" int aton_cell_index4_(int *i, int *j, int *k, int *component) {
  int ncellx, ncelly, ncellz, nbnd;
  aton_get_grid_size_(&ncellx, &ncelly, &ncellz, &nbnd);
  return  (*i+nbnd) +
    (*j+nbnd)*(ncellx+2*nbnd) +
    (*k+nbnd)*(ncellx+2*nbnd)*(ncelly+2*nbnd) +
    *component*(ncellx+2*nbnd)*(ncelly+2*nbnd)*(ncellz+2*nbnd);
}

extern "C" int aton_cell_index_(int *i, int *j, int *k) {
  int component = 0;
  return aton_cell_index4_(i, j, k, &component);
}
