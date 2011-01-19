extern "C" {
    #include "src/allvars.h"
    #include "src/proto.h"
    #include "src/tags.h"
    #include <gsl/gsl_rng.h>
}

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
} sph_state;


extern "C" void   begrun(void);
extern "C" double second(void);

int find_particle(int index_of_the_particle, struct particle_data **Pfound, int *at_task);
int found_particle(int index_of_the_particle, int *local_index);
int find_sph_particle(int index_of_the_particle, struct particle_data **Pfound, 
    struct sph_particle_data **SphPfound, int *at_task);
void update_particle_map(void);

extern "C" void hydro_state_at_point(FLOAT pos[3], FLOAT vel[3], FLOAT *h_out,
  FLOAT *ngb_out, FLOAT *dhsml_out, FLOAT *rho_out, FLOAT *rhov_out,
  FLOAT *rhov2_out, FLOAT *rhoe_out);


