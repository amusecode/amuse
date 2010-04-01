#include<fstream>
#include<vector>

#include "integrator.h"


typedef double  real;


// Interface functions for any MUSE gravity module.
// The structure of this file follows/defines the content of gravity.py,
// also in this directory.  (We probably should come up with a way to
// construct one file from the other...)

//-------------------------------------------------------------------------
//
/// Perform any necessary global initialization, and set the system
/// time to the specified value. Use re-initialize to make changes
/// in the middle of a run without changing the time, eg. when adding
/// a particle (see gravity.h for more details

int initialize_particles(double time);
int reinitialize_particles();

//-------------------------------------------------------------------------
//
/// Perform any necessary global cleanup, possibly prior to a new start.

int setup_module();
int cleanup_module();


//-------------------------------------------------------------------------
//
/// Structure defining particle dynamics data.

//  Use components to avoid possible SWIG problems with vectors.

typedef struct {
    int id;                                             /// identifier
    double mass;                                        /// mass
    double radius;                                      /// radius
    double x, y, z;                                     /// position
    double vx, vy, vz;                                  /// velocity
} dynamics_state;


//-------------------------------------------------------------------------
//
/// Add a particle described by d to the dynamical system.  Print an
/// error message and return without action if the particle with the
/// specified ID already exists.

/// Return the number of particles in the system after the particle
/// has been added.


int add_particle(dynamics_state d);


//-------------------------------------------------------------------------
//
/// Set the parameters of the particle described by d in the dynamical
/// system.  Return an error message and do nothing if the particle
/// does not exist.

int set_particle(int id, dynamics_state d);

//-------------------------------------------------------------------------

/// Set the mass of particle id.  Return 0 iff successful.

int set_mass(int id, double m);

//-------------------------------------------------------------------------
//
/// Set the radius of particle id.  Return 0 iff successful.

int set_radius(int id, double r);

//-------------------------------------------------------------------------
//
/// Set the position vector of particle id.  Return 0 iff successful.

int set_pos(int id, double pos[]);

//-------------------------------------------------------------------------
//
/// Set the velocity vector of particle id.  Return 0 iff successful.

int set_vel(int id, double pos[]);

//-------------------------------------------------------------------------
//
/// Delete the particle described by id from the dynamical system.

/// Return the number of particles in the system after the particle
/// has been removed.

int remove_particle(int id);

//-------------------------------------------------------------------------
//
/// Return the number of particles in the dynamical system.

int get_number();

//-------------------------------------------------------------------------
//
/// Return the current dynamical time scale of the dynamical system.

double get_dynamical_time_scale();

//-------------------------------------------------------------------------
//
/// Return the current time in the dynamical system.

double get_time_dynamics();

//-------------------------------------------------------------------------
//
/// Return the next time (t+dt) of the dynamical system.

double get_time_step();

//-------------------------------------------------------------------------
//
/// Perform any actions needed before the current time step.

int initialize_time_step();

//-------------------------------------------------------------------------
//
/// Perform any actions needed after the current time step.

int finalize_time_step();

//-------------------------------------------------------------------------
//
/// Evolve the dynamical system to the specified time.

/// On return, the dynamical system time will be greater than or equal
/// to t_end.  Return the total number of dynamical steps taken for
/// this segment of the integration.
int evolve(double t_end);

//-------------------------------------------------------------------------
//
/// Return the identifier of the more massive member of a colliding pair.

/// "Collision" means that the separation is less than the sum of the
/// radii.  Return the id of the primary if found, -1 otherwise.

int find_colliding_primary();

//-------------------------------------------------------------------------
//
/// Return the identifier of the other member of a colliding pair.

//* The primary is id; return -1 if no secondary found. */

int find_colliding_secondary(int id);

//-------------------------------------------------------------------------
//
/// Return the current dynamical state of particle id.

dynamics_state get_state(int id);

//-------------------------------------------------------------------------
//
/// Return the mass of particle id, or -1 if id doesn't exist.

double get_mass(int id);

//-------------------------------------------------------------------------
//
/// Return the radius of particle id, or -1 if id doesn't exist.

double get_radius(int id);

//-------------------------------------------------------------------------
//
/// Return the total kinetic energy of the system.

double get_kinetic_energy();

//-------------------------------------------------------------------------
//
/// Return the total potential energy of the system.

double get_potential_energy();

//-------------------------------------------------------------------------
//
/// Return the id of the next escaper found, -1 if none exists.
int get_escaper();

#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << "\n"

int get_index_from_identity(int id);
