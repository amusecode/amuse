#include "src/proto.hpp"
#include "src/tags.hpp"
#include <gsl/gsl_rng.h>


typedef struct {
    double mass;                                        /// mass
    double x, y, z;                                     /// position
    double vx, vy, vz;                                  /// velocity
} dynamics_state;

typedef struct {
    double mass;                                        /// mass
    double x, y, z;                                     /// position
    double vx, vy, vz;                                  /// velocity
    double u;                                           /// entropy
#ifdef MORRIS97VISC
    double alpha, dalphadt;				///viscosity
#endif
} sph_state;

//void   gadgetmg2::begrun(void);
//double gadgetmg2::second(void);

int found_particle(int index_of_the_particle, int *local_index);
void update_particle_map(void);

/*void gadgetmg2::hydro_state_at_point(double pos[3], double vel[3], double *h_out,
  double *ngb_out, double *dhsml_out, double *rho_out, double *rhov_out,
  double *rhov2_out, double *rhoe_out);
*/
