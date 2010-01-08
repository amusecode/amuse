
#include "smallN2.h"
#include <vector>

typedef double real;

// This is the same class as dynamics_state in the gravity modules.
// Using it here with a different name because we don't quite know how
// to handle multiple instances of the same object, or how to return
// compound objects using swig...

typedef struct {
    int id;
    double mass;
    double radius;
    double x, y, z;
    double vx, vy, vz;
} mult_cm_state;

void set_rmax_system(real rm);
real get_rmax_system();

void report_multiples(int level = 0);
void set_name_ranges(int nstart=100, int delta_n=100);

real semi_from_period(real m1, real m2, real period);
vec random_unit_vec();

// Stored multiples:

int add_binary(int id1, int id2, real mass1, real mass2,
	       real period, real eccentricity);
int remove_multiple(int id);

real get_energy(int id);
real get_total_energy();
real get_top_level_energy_error();
real get_mass(int id);
real get_radius(int id);
real get_n(int id);
bool is_multiple(int id);

// Set up an interaction:

void clear_multiple();
void add_to_interaction(int i, real m, real x[3], real v[3]);
void add_to_interaction(int i, real m, real x, real y, real z,
			real vx, real vy, real vz);

// Get the results of an interaction:

int get_status(int i);
void get_particle_result(int k, int *id, real *mass, real *x, real *y, real *z, real *vx, real *vy, real *vz);
void get_particle_original(int k, int *id, real *mass, real *x, real *y, real *z, real *vx, real *vy, real *vz);
int is_new_particle(int k);

// Integrate a multipe interaction:

int integrate_multiple(real *end_time, int verbose = 1, real eps2 = 0);
