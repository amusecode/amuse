#include<fstream>
#include<vector>

#include "src/stdinc.h"
#include "src/vec.h"
#include "src/nbody_particle.h"

typedef double  real;


// Interface functions for any MUSE gravity module.
// The structure of this file follows/defines the content of gravity.py,
// also in this directory.  (We probably should come up with a way to
// construct one file from the other...)

//-------------------------------------------------------------------------
//
/// Perform any necessary global initialization, before loading
/// particles into the system.

int setup_module();

//-------------------------------------------------------------------------
//
/// Perform any necessary global cleanup at end, possibly prior
/// to a new start.

int cleanup_module();

//-------------------------------------------------------------------------
//
/// Initialize all particles in the system (calculate accelerations,
/// jerks, time steps, etc.) and set the system time.  It is assumed
/// that the particles are already synchronized.  If not, this call
/// will (probably incorrectly) force the appearance of synchronization,
/// but no actual dynamical update will be performed.

int initialize_particles(double time = 0);

//-------------------------------------------------------------------------
//
/// Recompute the accelerations, jerks, etc. of all particles in the
/// system whose base times are equal to the current system time
/// (i.e. confine attention only to particles which have just been
/// advanced or added).  Recompute particle time steps (if relevant)
/// only if they aren't already set.  Use after adding new particles
/// to an existing system.

int reinitialize_particles();

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
int add_particle(int id, double mass, double radius, double x, double y, double z, double vx, double vy, double vz);

    
    

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

double get_time();

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
/// to t_end.  By default, particles will not be synchronized at the
/// specified time; rather, all particles will have
///
///	t0 <= t_end < t0 + dt
///
/// so the system is ready to take a new step on return without any
/// reinitialization.  To force synchronization, set the optional second
/// argument to 1.

/// OLD: Return the total number of dynamical steps taken for this
/// segment of the integration.  

/// NEW: Return the id of the more massive star involved in a collision.

int evolve(double t_end, int synchronize = 0);

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


void get_state(int id, int * id_out,  double * mass, double * radius, double * x, double * y, double * z, double * vx, double * vy, double * vz);

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

typedef nbody_particle real_particle;
typedef nbody_system real_system;
typedef nbody_VF_ptr real_VF_ptr;
typedef nbody_RF_ptr real_RF_ptr;
typedef nbody_RRF_ptr real_RRF_ptr;

extern "C" double cpusec();
int  pgetopt(int argc, char ** argv,  char * optstr);
void pskipopt();

// dynamics_state:  Structure describing the state of a particle in the
//		    dynamics integration.

//#define USE_VEC

typedef real_system BHTC_SYSTEM;

int get_index_from_identity(int id);
