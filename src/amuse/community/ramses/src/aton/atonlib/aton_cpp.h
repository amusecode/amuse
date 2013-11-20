// Public C++ interface for ATON.

// WARNING:
// This interface is not finalized.
// Clearly it needs to be improved a lot.
// Incompatible changes are likely to be made.

// The GPU functions operate on global GPU arrays.
// Therefore, there can only be a single GPU ATON instance per process.
// However, there is also a CPU implementation which doesn't use global arrays.
// It's possible to have several CPU ATON instances per process.

#ifndef ATON_H
#define ATON_H

namespace aton {

  struct State {
    int size;
    
    // NOTE: We use Fortran memory layout for the arrays.

    // Scalar quantities:
    double* E;    // Photon number density [photons / m^3]
    double* nH;   // Hydrogen gas number density [atoms / m^3]
    double* T;    // Temperature [Kelvin]
    double* xHII; // Ionization fraction [dimensionless]
    double* photon_source; // photon source field [photons / m^3 / s]

    // 3-vector quantities:
    double* F; // Photon flux [photons / m^2 / s]

    // Point sources
    // TODO: These should be deleted. Rather use the photon_source field. 
    int point_source_count;
    const int* point_source_pos;
    const double* point_source;

    // Set the arrays to these constant values.
    void Init(double E, double nH, double T, double xHII, double photon_source);
  };

  int init_cuda_device(bool allow_gpu_overload);
  void gpu_allocate();

  void get_grid_size(int *ncellx, int *ncelly, int *ncellz, int *nbound);
  int get_boundary_buffer_size();
  int cell_index4(int i, int j, int k, int component);

  void cpu_allocate(State* state);

  void cpu_to_gpu(State state);
  void gpu_to_cpu(State state);
  void cpu_to_gpu_boundary_values(const double *values);
  void gpu_to_cpu_boundary_values(double *values);

  void gpu_transport(double c_light, double dx, double dt);
  void gpu_add_sources(double dx, double dt);
  void gpu_cooling(double c_light, double dx, double dt,
		   double aexp, double hubblet, double fudgecool);

  void cpu_transport(State state, double c_light, double dx, double dt);
  void cpu_add_sources(State state, double dx, double dt);
  void cpu_cooling(State state, double c_light, double dx, double dt,
		   double aexp, double hubblet, double fudgecool);
  
  // Returns false if the state is inconsistent or invalid.
  // Also prints an error.
  // FIXME: Currently this always returns true.
  bool validate(State state, double c_light);
}

#endif
