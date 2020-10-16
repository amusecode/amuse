/* ################################################################################## */
/* ###                                                                            ### */
/* ###                                 Gadgetmp2                                  ### */
/* ###                                                                            ### */
/* ###   Original: Gadget2 in the version used in Amuse                           ### */
/* ###   Author: Gadget2 and Amuse contributors                                   ### */
/* ###                                                                            ### */
/* ###   Modified: July 2020                                                      ### */
/* ###   Author: Thomas Schano                                                    ### */
/* ###                                                                            ### */
/* ###   Changes are intended to enable precise calculations in                   ### */
/* ###   non periodic small domain simulations in which comoving parts            ### */
/* ###   are simulated in std precision                                           ### */
/* ###                                                                            ### */
/* ################################################################################## */
#include "src/tags.hpp"
#include <gsl/gsl_rng.h>
#include "src/proto.hpp"


typedef struct {
    my_float mass;                                        /// mass
    my_float x, y, z;                                     /// position
    my_float vx, vy, vz;                                  /// velocity
    my_float radius;                                      /// radius
} dynamics_state;

typedef struct {
    my_float mass;                                        /// mass
    my_float x, y, z;                                     /// position
    my_float vx, vy, vz;                                  /// velocity
    my_float u;                                           /// entropy
#ifdef MORRIS97VISC
    my_float alpha, dalphadt;				///viscosity
#endif
} sph_state;

//void   gadgetmp2::begrun(void);
//double gadgetmp2::second(void);

int found_particle(int index_of_the_particle, int *local_index);
void update_particle_map(void);

/*void gadgetmp2::hydro_state_at_point(double pos[3], double vel[3], double *h_out,
  double *ngb_out, double *dhsml_out, double *rho_out, double *rhov_out,
  double *rhov2_out, double *rhoe_out);
*/
