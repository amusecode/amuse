extern "C" {
    #include "src/allvars.h"
    #include "src/proto.h"
    #include "src/tags.h"
    #include <gsl/gsl_rng.h>
}

typedef struct {
    int id;                                             /// identifier
    double mass;                                        /// mass
    double x, y, z;                                     /// position
    double vx, vy, vz;                                  /// velocity
} dynamics_state;

typedef struct {
    int id;                                             /// identifier
    double mass;                                        /// mass
    double x, y, z;                                     /// position
    double vx, vy, vz;                                  /// velocity
    double u;                                           /// entropy
} sph_state;


extern "C" void   begrun(void);
extern "C" double second(void);

int find_particle(int index_of_the_particle, struct particle_data *Pfound);
int find_sph_particle(int index_of_the_particle, struct sph_particle_data *SphPfound);




