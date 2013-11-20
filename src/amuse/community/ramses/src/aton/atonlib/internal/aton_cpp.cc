#include <algorithm>
#include <cstdlib>

#include "aton_cpp.h"
#include "aton_fortran.h"
#include "gpu.h"

namespace aton {

void State::Init(double E, double nH, double T, double xHII,
		 double photon_source) {
  std::fill(this->E, this->E+size, E);
  std::fill(this->nH, this->nH+size, nH);
  std::fill(this->T, this->T+size, T);
  std::fill(this->xHII, this->xHII+size, xHII);
  std::fill(this->photon_source, this->photon_source+size, photon_source);
  std::fill(this->F, this->F+3*size, 0.0);
}

int init_cuda_device(bool allow_gpu_overload) {
  int overload = allow_gpu_overload;
  return aton_init_gpu_(&overload);
}

void gpu_allocate() {
  int num_sources = 0;
  aton_gpu_malloc_(&num_sources);
}

void get_grid_size(int *ncellx, int *ncelly, int *ncellz, int *nbound) {
  aton_get_grid_size_(ncellx, ncelly, ncellz, nbound);
}

int cell_index4(int i, int j, int k, int component) {
  return aton_cell_index4_(&i, &j, &k, &component);
}
  
void cpu_allocate(State* state) {
  int ncellx, ncelly, ncellz, nbnd;
  aton_get_grid_size_(&ncellx, &ncelly, &ncellz, &nbnd);
  int n = (ncellx + 2*nbnd)*(ncelly + 2*nbnd)*(ncellz + 2*nbnd);

  state->size = n;
  state->E = new double[n];
  state->nH = new double[n];
  state->T = new double[n];
  state->xHII = new double[n];
  state->photon_source = new double[n];
  state->F = new double[3*n];
}

void cpu_to_gpu(State state) {
  int zero_sources = 0;
  aton_cpu_to_gpu_full_(state.E, state.F, state.xHII, state.T, state.nH,
		        NULL, NULL, &zero_sources,
		        state.photon_source);
}

void gpu_to_cpu(State state) {
  int zero_sources = 0;
  aton_gpu_to_cpu_full_(state.E, state.F, state.xHII, state.T, state.nH,
		        NULL, NULL, &zero_sources);
}


void gpu_transport(double c_light, double dx, double dt) {
  gpu_rad_transport(c_light, dx, dt);
}

void gpu_add_sources(double dx, double dt) {
  gpu_rad_add_sources(dx, dt, 0);
}

void gpu_cooling(double c_light, double dx, double dt,
		 double aexp, double hubblet, double fudgecool) {
  gpu_rad_cooling(c_light, dx, dt, aexp, hubblet, fudgecool);
}

}
