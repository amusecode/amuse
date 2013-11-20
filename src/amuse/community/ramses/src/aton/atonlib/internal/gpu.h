// GPU definitions.

// Unfortunately the grid size is fixed at compile time.
// It is defined by the NCELLX, NCELLY and NCELLZ macros.
// These are always passed on the compiler command line
// so that they can be configured in the Makefile.

// Number of boundary cells.
// This must be at least 1.
// Larger values (e.g. 4) can improve performance due to memory alignment.
#define NBOUND 1  // This must match nbnd in coupling.f90
#define NBOUND2 (2*NBOUND)


void gpu_rad_transport(double c_light, double dx, double dt);
void gpu_rad_add_sources(double dx, double dt, int nsource);
void gpu_rad_cooling(double c_light, double dx, double dt,
                     double aexp, double hubblet, double fudgecool);


// GPU global arrays
// They are defined in gpu_memory.cc.
// See aton_cpp.h or aton.f90 for descriptions.
// FIXME: These need more descriptive names.

extern double *cuegy, *cuegy_new;
extern double *cuflx, *cuflx_new;

extern double *cudedd;
extern double *cusrc0;
extern int  *cusrc0pos;

extern double *cutemperature;
extern double *cuxion;
extern double *cudensity;

extern double *cu_photon_source;

// Boundary values of cuegy and cuflx for efficient transfer at radiation
// substeps.
// length: 6 * 4 * max(ncellx,ncelly,ncellz)^2
// contents: [surface][{e,fx,fy,fz}][i][j]
extern double *cu_boundary_values;
