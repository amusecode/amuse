#include<fstream>
#include<vector>

#include "src/stdinc.h"
#include "src/vec.h"
#include "src/nbody_particle.h"

typedef double  real;

typedef struct {
    int id;                                             /// identifier
    double mass;                                        /// mass
    double radius;                                      /// radius
    double x, y, z;                                     /// position
    double vx, vy, vz;                                  /// velocity
} dynamics_state;

#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << "\n"

typedef nbody_particle real_particle;
typedef nbody_system real_system;
typedef nbody_VF_ptr real_VF_ptr;
typedef nbody_RF_ptr real_RF_ptr;
typedef nbody_RRF_ptr real_RRF_ptr;

extern "C" double cpusec();

int  pgetopt(int argc, char ** argv,  char * optstr);
void pskipopt();

//#define USE_VEC

typedef real_system BHTC_SYSTEM;


int get_identity_from_index(int i);

int get_index_from_identity(int id);

int find_colliding_primary();

int find_colliding_secondary(int id);
