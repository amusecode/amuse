// Public Fortran interface to ATON.

// The library operates with global arrays on the GPU.
// Therefore, there can only be a single ATON instance per process.
// Also, the GPU grid size is fixed at compile time.
// You need to edit internal/gpu.h to change the size.
//
// ATON always uses 'physical' quantities, never comoving quantities.
// TODO: Explain the arrays we use.


// Initialize the CUDA device.
//
// If hosts have several CPU cores and GPU devices then multiple MPI processes
// can run on the same host. This function ensures that each process on the host
// accesses a different GPU device.
//
// allow_gpu_overload:
//   If this flag is true and there are more MPI processes than GPU devices on
//   the host then an error is raised.
//
// Returns 1 on success, 0 on failure.
//
extern "C" int aton_init_gpu_(const int* allow_gpu_overload);

// The GPU grid size is fixed at compile time.
// This function returns the compiled-in dimensions.
extern "C" void aton_get_grid_size_(int *ncellx, int *ncelly, int *ncellz, int *nbound);

// We always use Fortran array layout. These functions return the array index for scalar and vector fields.
extern "C" int aton_cell_index_(int *i, int *j, int *k);
extern "C" int aton_cell_index4_(int *i, int *j, int *k, int *component);

// Allocate the GPU global arrays.
extern "C" void aton_gpu_malloc_(int *nsrc);

// Copy the full arrays between the CPU and GPU.
// The arguments refer to CPU arrays in both cases.
//
// The array dimensions should match aton_get_grid_size.
//
// rad_N: photon number density [photons / m^3], scalar
// rad_F: photon flux [photons / m^2 / s], vector
// xion: Ionization fraction [dimensionless], scalar
// temperature: Temperature [Kelvin], scalar
// nH: Gas number density [atoms / m^3], scalar
//
// source_N: Emission rates of the point sources [photons / m^3 / s]
// source_pos: Grid positions of the point sources.
// source_count: Number of point sources.
//
// photon_source: photon source field [photons / m^3 / s], scalar
//
extern "C" void aton_gpu_to_cpu_full_(
    double *rad_N, double *rad_F, double *xion, double *T, double *nH,
    double *source_N, int *source_pos, int *source_count);
// FIXME: The arguments should be const.
extern "C" void aton_cpu_to_gpu_full_(
    double *rad_N, double *rad_F, double *xion, double *T, double *nH,
    double *source_N, int *source_pos, int *source_count,
    double *photon_source);

// Pack the boundary cells into buffers and transfer them between the CPU and GPU.
// This is faster than transferring the full arrays.
// The size of boundary_values is 4*(max(NCELLX, max(NCELLY, NCELLZ))).
// TODO: Provide a function to get the buffer size.
extern "C" void aton_gpu_to_cpu_boundary_(double *boundary_values);
extern "C" void aton_cpu_to_gpu_boundary_(const double *boundary_values);

// Dump the given CPU arrays to a file for debugging.
extern "C" void aton_debug_dump_(double *e, double *f, double *x, double *temp, double *src0, int *src0pos, double *dens, int *nsource, double *photon_source, double *time, int *isnap);

// Run the full ATON step.
// The function operates on the ATON global arrays.
// All the arguments are pointers so that it can be called from Fortran.
//
// c: speed of light [m/s]
// dx: grid spacing [m]
// dt: time step [s]
// source_count: number of point sources
// fudgecool: factor by which to reduce the cooling time step if the initial step is unacceptable
// aexp: cosmological scale factor (aexp=1 today)
// hubblet: hubble time 1/H(aexp) [s]
//
extern "C" void aton_gpu_step_(
	const double* c, const double* dx, const double* dt,
	const int* source_count, const double* fudgecool,
	const double* aexp, const double* hubblet);

// Check that the quantities make sense.
//
// For example, the following checks are included:
//  nH > 0
//  T > 0
//  0 <= x <= 1
//  F > cN
// If any cell is invalid, an error message is printed.
//
// For debugging, I recommend calling this function a lot. It can help to find the code that is generating problems.
//
// label: a label to include in the error message (useful for debugging)
// c: speed of light [m/s]
//
extern "C" void aton_validate_(
    const int *label, const double *c,
    const double *cpu_e, const double* cpu_d, const double* cpu_t,
    const double* cpu_x, const double *cpu_f);
