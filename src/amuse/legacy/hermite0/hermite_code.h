
// Interface functions for the MUSE dynamics module.
//
// This file differs in content from ../gravity.h only in the
// declarations of get_colliding_primary/secondary below, to
// illustrate (in muse_dynamics.py) how to overcome minor
// differences in specifications to comply with the interface.

//-------------------------------------------------------------------------
//
/// Perform any necessary global initialization, before loading
/// particles into the system.

int setup_module(bool reeval = false, bool test_mode = false);

//-------------------------------------------------------------------------
//
/// Perform any necessary global cleanup at end, possibly prior
/// to a new start.

int cleanup_module();

//-------------------------------------------------------------------------
//
/// Recompute the accelerations, jerks, etc. of all particles in the
/// system whose base times are equal to the current system time
/// (i.e. confine attention only to particles which have just been
/// advanced or added).  Recompute the time step only if it isn't
/// already set.

int reinitialize_particles();


//-------------------------------------------------------------------------
//
/// Add a particle described by d to the dynamical system.  Print an
/// error message and return without action if the particle with the
/// specified ID already exists.

/// Return the number of particles in the system after the particle
/// has been added.

int new_particle(int *id, double mass, double radius, double x, double y, double z, double vx, double vy, double vz);

/// Delete the particle described by id from the dynamical system.

/// Return the number of particles in the system after the particle
/// has been removed.

//int remove_particle(int id);
//
int delete_particle(int id);

//-------------------------------------------------------------------------

/// Return the current dynamical state of particle id.

int get_state(int id, double *_mass, double *_radius, double *x, double *y, double *z, double *vx, double *vy, double *vz);
//-------------------------------------------------------------------------

int set_state(int id, double _mass, double _radius, double x, double y, double z, double vx,double vy, double vz);
//-------------------------------------------------------------------------

/// Return the mass of particle id, or -1 if id doesn't exist.

int get_mass(int id, double *mass);

//-------------------------------------------------------------------------

/// Set the mass of particle id.  Return 0 if successful.

int set_mass(int id, double m);

/// Return the radius of particle id, or -1 if id doesn't exist.

int get_radius(int id, double *_radius);

//-------------------------------------------------------------------------

/// Set the radius of particle id.  Return 0 iff successful.

int set_radius(int id, double r);

//-------------------------------------------------------------------------

/// Get the position of particle id. Return 0 if successful.

int get_position(int id, double *x, double *y, double *z);

//-------------------------------------------------------------------------

/// Set the position vector of particle id.  Return 0 if successful.
int set_position(int id, double x, double y, double z);
//-------------------------------------------------------------------------


int get_velocity(int id, double *vx, double *vy, double *vz);

/// Set the velocity vector of particle id.  Return 0 iff successful.
int set_velocity(int id, double vx, double vy, double vz);
//-------------------------------------------------------------------------

int get_acceleration(int id, double *ax, double *ay, double *az);
//-------------------------------------------------------------------------

int set_acceleration(int id, double ax, double ay, double az);
//-------------------------------------------------------------------------

///pot @ point
///BE AWARE if point is @ particle and if smoothing is naught 
///it is a singularity 
int get_potential(double x, double y, double z,  double *V);
//-------------------------------------------------------------------------

/*
Do we want a int get_potential @ particle (excluding the particle)??

int get_potential_at_particle(int id, double *V);

*/

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

int evolve(double t_end);

//-------------------------------------------------------------------------
//
/// Initialize all particles in the system (calculate accelerations,
/// jerks, time steps, etc.) and set the system time.  It is assumed
/// that all particles are already synchronized.  If not, this call
/// will (probably incorrectly) force the appearance of synchronization,
/// but no actual dynamical update will be performed.

int initialize_particles(double time = 0);

//-------------------------------------------------------------------------

///
int initialize_code();
//-------------------------------------------------------------------------


int get_eps2(double *eps2);

int set_eps2(double eps2);

/// Return the total kinetic energy of the system.

int get_kinetic_energy(double *kinetic_energy);

//-------------------------------------------------------------------------
/// Return the total potential energy of the system.

int get_potential_energy(double *potential_energy);

//-------------------------------------------------------------------------

int get_indices_of_colliding_particles(int *id1, int *id2);
/// Return the current time in the dynamical system.

int get_time(double *time);

//-------------------------------------------------------------------------

int get_total_mass(double *_mass);
int get_center_of_mass_position(double *x, double *y, double *z);
int get_center_of_mass_velocity(double *vx, double *vy,double *vz);
int get_total_radius(double *radius);
int get_gravity_at_point(double eps, double x, double y, double z,  double *forcex, double *forcey, double *forcez);
/// Return the number of particles in the dynamical system.
int get_number_of_particles(int *no_parts);
//int get_number();
///////////////
/// weet nog niet ########################
int get_index_of_first_particle(int *index_of_first_particle);
int get_index_of_next_particle(int id, int *index_of_the_next_particle);
//
/// Set the parameters of the particle described by d in the dynamical
/// system.  Return an error message and do nothing if the particle
/// does not exist.

//int set_particle(int id, double mass, double radius, double x, double y, double z, double vx, double vy, double vz);
    

//-------------------------------------------------------------------------

//-------------------------------------------------------------------------


//-------------------------------------------------------------------------

/// Return the current dynamical time scale of the dynamical system.

int get_dynamical_time_scale(double *ts);

//-------------------------------------------------------------------------


/// Return the time step of the dynamical system.

int get_time_step(double *time_step);

//-------------------------------------------------------------------------

/// Perform any actions needed before the current time step.

int initialize_time_step();

//-------------------------------------------------------------------------

/// Perform any actions needed after the current time step.

int finalize_time_step();

/// Return the total number of time steps taken.

int get_n_steps();

//-------------------------------------------------------------------------

/// Return the identifier of the more massive member of a colliding pair.

/// "Collision" means that the separation is less than the sum of the
/// radii.  Return the id of the primary if found, -1 otherwise.

////**** Note that the next two functions are "misnamed" relative to ****
////**** the definitions in gravity.h, to illustrate how internal    ****
////**** names can be reworked to create python functions.           ****
////****                                               (Steve, 1/07) ****

int get_colliding_primary();

//-------------------------------------------------------------------------

/// Return the identifier of the other member of a colliding pair.

//* The primary is id; return -1 if no secondary found. */

int get_colliding_secondary(int id1);

//-------------------------------------------------------------------------


int get_indices_of_colliding_particles(int *id1, int *id2);

/// Return the id of the next escaper found, -1 if none exists.

int get_escaper();

int get_potential_at_point(double eps, double x, double y, double z, double * phi);

