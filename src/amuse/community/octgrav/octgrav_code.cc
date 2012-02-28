#include <iostream>
#include "octgrav_code.h"
#include "worker_code.h"
#include <vector>
#include <algorithm>
#include <string>
#include <map>
//AMUSE STOPPING CONDITIONS SUPPORT
#include <stopcond.h>
#include <time.h>

using namespace std;

// global static parameters
// N-body data:

// Control parameters:

// Module data:

vector<float4> bodies_pos(0);    // (w = mass)
vector<float4> bodies_grav(0);   // (w = potential)
vector<float4> bodies_vel(0);    // (w serves no purpose currently)
vector<int>    starids(0);       // list of identifiers
vector<float>  radii(0);         // list of radii

double timestep = 0.01;
double total_mass = 0.0;
bool initialized = false;
int verbose_mode = 0;

int id_counter = 0; // for each new particle is this is
                    // being incremented to give unique id for new particle.

int get_index_from_identity(int id)
{
  for (int i = 0; i < n_bodies; i++)
    {
      if (id == starids[i])
        {
          return i;
        }
    }
  return -1;
}

// Interface functions:
int new_particle(int *id, double mass, double x, double y, double z, double vx, double vy, double vz, double radius)
{
  // Add the new particle and reinitialize immediately.

  bodies_pos.resize(n_bodies+1);
  bodies_pos[n_bodies].x = x;
  bodies_pos[n_bodies].y = y;
  bodies_pos[n_bodies].z = z;
  bodies_pos[n_bodies].w = mass;

  bodies_vel.resize(n_bodies+1);
  bodies_vel[n_bodies].x = vx;
  bodies_vel[n_bodies].y = vy;
  bodies_vel[n_bodies].z = vz;
  bodies_vel[n_bodies].w = 0;

  bodies_grav.resize(n_bodies+1);

  n_bodies++;
  id_counter++;

  starids.push_back(id_counter);
  radii.push_back(radius);

  total_mass += mass;

  *id = id_counter;

  return 0;
}

int delete_particle(int id)
{
    int i = get_index_from_identity(id);

    if (i >= 0 && i < n_bodies)
      {
        bodies_pos.erase(bodies_pos.begin()+i);
        bodies_vel.erase(bodies_vel.begin()+i);
        bodies_grav.erase(bodies_grav.begin()+i);
        starids.erase(starids.begin()+i);
        radii.erase(radii.begin()+i);
        n_bodies--;
        return 0;
      }
    else
      {
        return -1;
      }
}

static void create_treecode_system();

void print_tree_counters();

// Horrible hack to extract colliders from the tree traversals
// prompted by integrate() below.

///TODO: implement at some point?
///static int id_primary = -1, id_secondary = -1;

///extern int id_collision_1, id_collision_2;   // determined in BHTC/BHtree.C
///extern double r_collision_12;


int evolve_model(double t_end)
{
  // advance from the current time to t_end
  /* 1. Create equal-sized time steps by adjusting timestep */

  double delta_t = t_end - t_now;
  int nsteps;
  //  dtime = timestep;
  nsteps = (int) (delta_t/timestep)+1;
  double dtime = delta_t/nsteps;

  /* 2. Calculate energies */
  ///calculate the inital energies

  double E_kin = calcEkin(bodies_pos, bodies_vel);
  double E_pot = calcEpot(bodies_pos, bodies_grav);

  /* 4. Integrate as long as necessary */
  static octgrav system;
  system.set_softening(eps);
  system.set_opening_angle(theta);

  //AMUSE STOPPING CONDITIONS
  int is_number_of_steps_detection_enabled;
  int number_of_steps_innerloop = 0;
  int max_number_of_steps;
  int error;
  int number_of_steps_detection;
  
  error = is_stopping_condition_enabled(NUMBER_OF_STEPS_DETECTION, 
					&is_number_of_steps_detection_enabled);
  get_stopping_condition_number_of_steps_parameter(&max_number_of_steps);    

  //stop @ nsteps or max_number_of_steps whichever is smallest

  number_of_steps_detection = 0; //false
  if (is_number_of_steps_detection_enabled) {
      if (max_number_of_steps<=nsteps) {
	  nsteps = max_number_of_steps;
	  number_of_steps_detection =1;
      }
  }

  fprintf(stdout, "eps:%f theta:%f ", eps, theta);
  fflush(stdout);

  for (int i = 0; i < nsteps ; i++) {
      leapfrog(dtime, bodies_pos, bodies_vel, bodies_grav, system);
  }

  if (number_of_steps_detection) {
      int stopping_index  = next_index_for_stopping_condition();
      set_stopping_condition_info(stopping_index, NUMBER_OF_STEPS_DETECTION);
  }

  return 0;
}

int get_number_of_particles(int *number_of_particles)
{
  *number_of_particles = n_bodies;
  return 0;
}

int remove_particle(int id)
{
    if (!initialized) return -1;

    int i = get_index_from_identity(id);

    if (i >= 0 && i < n_bodies)
      {
        bodies_pos.erase(bodies_pos.begin()+i);
        bodies_vel.erase(bodies_vel.begin()+i);
        bodies_grav.erase(bodies_grav.begin()+i);
        starids.erase(starids.begin()+i);
        radii.erase(radii.begin()+i);
      }

    return n_bodies;
}

int get_potential(int id, double * value)
{
  int i = get_index_from_identity(id);
  if (i >= 0 && i < n_bodies)
    {
      *value = bodies_grav[i].w;
      return 0;
    }
    else
    {
        *value = 0.0;
        return -1;
    }
}

int get_potential_energy(double *potential_energy)
{
  *potential_energy = calcEpot(bodies_pos, bodies_grav);
  return 0;
}

int get_kinetic_energy(double *kinetic_energy)
{
  *kinetic_energy = calcEkin(bodies_pos, bodies_vel);
  return 0;
}


int get_gravity_at_points(int n, double *x, double *y, double *z,  double *forcex, double *forcey, double *forcez, double *pot)
{
  int i;
  if(n_bodies <=0 ) return -1;
  bodies_pos.resize(n_bodies+n);
  bodies_vel.resize(n_bodies+n);
  bodies_grav.resize(n_bodies+n);
  for(i=0;i<n;i++){
    bodies_pos[n_bodies+i].x = *(x+i);
    bodies_pos[n_bodies+i].y = *(y+i);
    bodies_pos[n_bodies+i].z = *(z+i);
    bodies_pos[n_bodies+i].w = 0.;
  }

  {
    octgrav system;
    system.set_softening(eps);
    system.set_opening_angle(theta);
    system.evaluate_gravity(bodies_pos, bodies_grav);
  }
  for(i=0;i<n;i++)
  {
    if(forcex != NULL) *(forcex+i)=bodies_grav[n_bodies+i].x;
    if(forcey != NULL) *(forcey+i)=bodies_grav[n_bodies+i].y;
    if(forcez != NULL) *(forcez+i)=bodies_grav[n_bodies+i].z;
    if(pot != NULL)*(pot+i)=bodies_grav[n_bodies+i].w;
  }
  bodies_pos.resize(n_bodies);
  bodies_vel.resize(n_bodies);
  bodies_grav.resize(n_bodies);
  return 0;
}

int get_gravity_at_point(double *eps, double *x, double *y, double *z,  double *forcex, double *forcey, double *forcez,int n)
{
  return get_gravity_at_points(n, x,y,z, forcex,forcey,forcez, NULL);
}

int get_potential_at_point(double *eps, double *x, double *y, double *z, double *pot, int n)
{
  return get_gravity_at_points(n, x,y,z, NULL, NULL,NULL,pot);
}


double get_dynamical_time_scale()
{
    if (!initialized) return 0;

    double mtot = total_mass;
    double ekin;
    double epot;
    get_kinetic_energy(&ekin);
    get_potential_energy(&epot);
    //double etot = ekin + epot;

    double tdyn = (-0.5*mtot*mtot/epot) / sqrt(2*ekin/mtot);
    return tdyn;
}

double get_time_dynamics()
{
  return t_now;
}

int get_state(int id, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *radius)
{
  int i = get_index_from_identity(id);
  if (i >= 0 && i < n_bodies)
    {
      *mass = bodies_pos[i].w;
      *radius = radii[i];

      *x = bodies_pos[i].x;
      *y = bodies_pos[i].y;
      *z = bodies_pos[i].z;

      *vx = bodies_vel[i].x;;
      *vy = bodies_vel[i].y;
      *vz = bodies_vel[i].z;
      return 0;
    }
  else
    {
      return -1;
    }
}

int get_time_step(double *_timestep)
{
  *_timestep = timestep;
  return 0;
}

int set_time_step(double _timestep)
{
  timestep = _timestep;
  return 0;
}

int get_total_radius(double *radius)
{
  return -2;
}

int get_total_mass(double *mass)
{
  *mass = 0;

  for (int i=0; i < n_bodies; i++)
    {
      *mass += bodies_pos[i].w;
    }
  return 0;
}

int get_center_of_mass_position(double *x, double *y, double *z)
{
  double M;
  double mass;

  *x = 0; *y=0; *z=0;

  get_total_mass(&M);

  for(int i = 0; i<n_bodies; i++)
    {
      mass = bodies_pos[i].w;

      *x += mass*bodies_pos[i].x;
      *y += mass*bodies_pos[i].y;
      *z += mass*bodies_pos[i].z;
    }

  *x /= M;
  *y /= M;
  *y /= M;

  return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz)
{
  double M;
  double mass;

  *vx = 0; *vy=0; *vz=0;

  get_total_mass(&M);

  for(int i = 0; i < n_bodies; i++)
    {
      mass = bodies_pos[i].w;

      *vx += mass*bodies_vel[i].x;
      *vy += mass*bodies_vel[i].y;
      *vz += mass*bodies_vel[i].z;
    }

  *vx /= M;
  *vy /= M;
  *vy /= M;

  return 0;
}

int get_acceleration(int id, double * ax, double * ay, double * az)
{
  return -2;
}

void clear_all()
{
  bodies_pos.clear();    //(w = mass)
  bodies_grav.clear();   //(w = potential)
  bodies_vel.clear();    //(w serves no purpose currently)
  starids.clear();  //list of identifiers

  n_bodies = 0;
  total_mass = 0.0;
}

int get_eps2(double *epsilon_squared)
{
  *epsilon_squared = eps*eps;
  return 0;
}

int get_time(double *time)
{
  *time = t_now;
  return 0;
}

int initialize_particles()
{
  //clear_all();

  //t_now = t;

  /* Do an initial gravity calculation */

  cerr << "initialize()..." << endl;

  ///Load the system (this currently has to be done before every force integration)
  octgrav system;
  system.set_softening(eps);
  system.set_opening_angle(theta);

  cerr << "...(2)... n_bodies = " << bodies_pos.size() << endl;

  if(verbose_mode == 1) {
    fprintf(stderr,"stats for particle 0: %f,%f,%f,%f,%f,%f.\n", bodies_pos[0].x, bodies_pos[0].y, bodies_pos[0].z, bodies_vel[0].x, bodies_vel[0].y,bodies_vel[0].z);
    fprintf(stderr,"stats for particle 32767: %f,%f,%f,%f,%f,%f.\n", bodies_pos[32767].x, bodies_pos[32767].y, bodies_pos[32767].z, bodies_vel[32767].x, bodies_vel[32767].y, bodies_vel[32767].z);
  }

  ///Single force integration
  fprintf(stderr, "Performing preliminary gravity evaluation.\n");
  system.evaluate_gravity(bodies_pos, bodies_grav);

  if(verbose_mode == 1) {
    fprintf(stderr,"stats for particle 0: %f,%f,%f,%f,%f,%f.\n", bodies_pos[0].x, bodies_pos[0].y, bodies_pos[0].z, bodies_vel[0].x, bodies_vel[0].y,bodies_vel[0].z);
    fprintf(stderr,"stats for particle 32767: %f,%f,%f,%f,%f,%f.\n", bodies_pos[32767].x, bodies_pos[32767].y, bodies_pos[32767].z, bodies_vel[32767].x, bodies_vel[32767].y,bodies_vel[32767].z);
  }

  /* 2. Calculate energies */
  ///calculate the initial energies
  double E_kin = calcEkin(bodies_pos,bodies_vel);
  double E_pot = calcEpot(bodies_pos,bodies_grav);

  E_init = E_kin + E_pot;
  //fprintf(stdout, "t: %lg  E_kin: %lg E_pot: %lg E_init: %lg\n", t_now, E_kin,E_pot,E_init);

  //fprintf(stdout, "t_now: %lg, t_end: %lg\n", t_now, t_end);

  initialized = true;

  return 0;
}

// reinitialize not needed(?) but here for completness of interface
int reinitialize_particles()
{
  return 0;
}

// setting and cleaning up module
// new interface definitions, setup_module added, finalize was renamed
// to cleanup                      (by SH 1-/NOV/08)
int setup_module()
{
  return 0;
}

int cleanup_module()
{
  clear_all();
  initialized = false;
  return 0;
}

int initialize_time_step()
{
  return 0;
}

int finalize_time_step()
{
  return 0;
}

int set_particle(int id, dynamics_state d)
{
  int i = get_index_from_identity(id);
  if (i >= 0 && i < n_bodies) {

    bodies_pos[i].x = d.x;
    bodies_pos[i].y = d.y;
    bodies_pos[i].z = d.z;
    bodies_pos[i].w = d.mass;

    bodies_vel[i].x = d.vx;
    bodies_vel[i].y = d.vy;
    bodies_vel[i].z = d.vz;
    bodies_vel[i].w = 0;

    radii[i] = d.radius;
    return 0;
  }
  else {
    cerr << "set_particle: " << id << " doesn't exist.  Use add_particle." << endl;
    return -1;
  }
}

int set_pos(int id, double x[])
{
  int i = get_index_from_identity(id);
  if (i >= 0 && i < n_bodies) {
    bodies_pos[i].x = x[0];
    bodies_pos[i].y = x[1];
    bodies_pos[i].z = x[2];
    return 0;
  }
  return -1;
}

int set_position(int id, double x, double y, double z)
{
  int i = get_index_from_identity(id);
  if (i >= 0 && i < n_bodies)
    {
      bodies_pos[i].x = x;
      bodies_pos[i].y = y;
      bodies_pos[i].z = z;
      return 0;
    }
  return -1;
}

int set_vel(int id, double v[])
{
  int i = get_index_from_identity(id);
  if (i >= 0 && i < n_bodies) {
    bodies_vel[i].x = v[0];
    bodies_vel[i].y = v[1];
    bodies_vel[i].z = v[2];
    return 0;
  } else
  return -1;
}

int set_velocity(int id, double vx, double vy, double vz)
{
  int i = get_index_from_identity(id);
  if (i >= 0 && i < n_bodies)
    {
      bodies_vel[i].x = vx;
      bodies_vel[i].y = vy;
      bodies_vel[i].z = vz;
      return 0;
    } else
    return -1;
}

int get_mass(int id, double *mass)
{
  int i = get_index_from_identity(id);
  if (i >= 0 && i < n_bodies)
    {
      *mass = bodies_pos[i].w;
      return 0;
    }
  return -1;
}

int get_radius(int id, double *radius)
{
  int i = get_index_from_identity(id);
  if (i >= 0 && i < n_bodies)
    {
      *radius = radii[i];
      return 0;
    }
  return -1;
}

int get_position(int id, double *x, double *y, double *z)
{
  int i = get_index_from_identity(id);
  if (i >= 0 && i < n_bodies)
    {
      *x = bodies_pos[i].x;
      *y = bodies_pos[i].y;
      *z = bodies_pos[i].z;
      return 0;
    }
  return -1;
}

int get_velocity(int id, double *vx, double *vy, double *vz)
{
  int i = get_index_from_identity(id);
  if (i >= 0 && i < n_bodies)
    {
      *vx = bodies_vel[i].x;
      *vy = bodies_vel[i].y;
      *vz = bodies_vel[i].z;
      return 0;
    }
  return -1;
}

int find_colliding_primary()
{
  return -1;
}

int get_index_of_first_particle(int * index_of_the_particle)
{
  return -2;
}

int get_index_of_next_particle(int index_of_the_particle, int *index_of_the_next_particle)
{
  return -2;
}


int find_colliding_secondary(int id1)
{
  return -1;
}


int get_escaper()
{
  return -1;
}

int set_radius(int id, double radius)
{
  int i = get_index_from_identity(id);
  if (i >= 0 && i < n_bodies)
    {
      radii[i] = radius;
      return 0;
    }
  return -1;
}

int set_mass(int id, double mass)
{
  int i = get_index_from_identity(id);
  if (i >= 0 && i < n_bodies)
    {
      bodies_pos[i].w = mass;
      return 0;
    }
  return -1;
}

int set_state(int id, double mass, double x, double y, double z, double vx, double vy, double vz, double radius)
{
  int i = get_index_from_identity(id);
  if (i >= 0 && i < n_bodies)
    {
      bodies_pos[i].w = mass;
      radii[i] = radius;

      bodies_pos[i].x = x;
      bodies_pos[i].y = y;
      bodies_pos[i].z = z;

      bodies_vel[i].x = vx;
      bodies_vel[i].y = vy;
      bodies_vel[i].z = vz;
    }

  return 0;
}

int set_eps2(double epsilon_squared)
{
  eps = sqrt(epsilon_squared);
  return 0;
}

int set_acceleration(int index_of_the_particle, double ax, double ay, double az)
{
  return -2;
}

int synchronize_model()
{
  return 0;
}

int recommit_particles()
{
  //tree setup
  return 0;
}

int set_theta_for_tree(double theta_for_tree)
{
  theta = theta_for_tree;
  return 0;
}

int get_theta_for_tree(double *theta_for_tree)
{
  *theta_for_tree = theta;
  return 0;
}

int commit_parameters()
{
  return 0;
}

int recommit_parameters()
{
  return commit_parameters();
}

int commit_particles()
{
  //tree setup 
  initialize_particles();
  return 0;
}

int cleanup_code()
{
  clear_all();
  initialized = false;
  return 0;
}

int initialize_code()
{
    set_support_for_condition(NUMBER_OF_STEPS_DETECTION);
    return 0;
}
