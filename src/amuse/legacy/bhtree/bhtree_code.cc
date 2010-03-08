#include <iostream>
#include "bhtree_code.h"
#include "worker_code.h"
#include "local.h"
#include <vector>
#include <algorithm>
#include <string>
#include <map>

using namespace std;

const bool debug = false;

// global static parameters
// N-body data:

// Control parameters and accessors:

const real DT_DIA = 1;
real dt_dia = DT_DIA;                // time interval between diagnostics output
void set_dt_dia(real d)                        {dt_dia = d;}

//bool x_flag = false;

//int nsteps = 0;                // number of integration time steps completed
//real einit = 0;                // initial total energy of the system
//bool init_flag = false;
//real t_dia = 0;

// Module data:

bool initialized = false;

BHTC_SYSTEM bhtcs;
vector<dynamics_state> ds;        // for initialization only

// External parameters to be used inside the treecode.
// Their definitions are extended to python class modules in the file
// parameters.h; the accessor functions are for the driver program.

const real TIMESTEP = 0.015625;
real timestep =        TIMESTEP;
void set_timestep(real dt)                {timestep = dt;}

const real EPS2_FOR_GRAVITY = 0.125;
real eps2_for_gravity = EPS2_FOR_GRAVITY;
void set_eps2_for_gravity(real eps2)        {eps2_for_gravity = eps2;}

const real THETA_FOR_TREE = 0.75;
real theta_for_tree = THETA_FOR_TREE;
void set_theta_for_tree(real theta)        {theta_for_tree = theta;}

const int USE_SELF_GRAVITY = 1;
int use_self_gravity = USE_SELF_GRAVITY;
void set_use_self_gravity(int use)        {use_self_gravity = use;}

const int NCRIT_FOR_TREE = 1024;
int ncrit_for_tree = NCRIT_FOR_TREE;
void set_ncrit_for_tree(int ncrit)        {ncrit_for_tree = ncrit;}

// Interface functions:

int _new_particle(int *id, dynamics_state d)
{
//    cerr << d.id << " " << d.x << " " << d.y << " " << d.z << " "
//         << d.radius << endl;


    if (!initialized) 
      {
        // Defer reinitialization; save in ds.
        ds.push_back(d);
        *id = ds.size();
        //return ds.size();
        return 0;
      } 
    else 
      {
        //cello, I think this is obsolete...the code is responsible for
        //new particle id, 

        
        //int i = get_index_from_identity(d.id);
        //if (i >= 0 && i < bhtcs.n) 
        //  {
        //    // Particle already exists.  Do nothing.
        //    /*
        //    cerr << "add_particle: " << d.id
        //         << " already exists.  Use set_particle." << endl;
        //    */
        //    return -1;
        //  } 
        //else 
        //  { 
        

        // Add the new particle and reinitialize immediately.
        // cerr << "Adding particle " << d.id << endl;
        
        int n1 = bhtcs.n + 1;
        nbody_particle *np = bhtcs.get_particle_pointer();
        nbody_particle *np1 = new nbody_particle[n1];
        // Copy the system.
        for (int i = 0; i < bhtcs.n; i++) 
          {
            np1[i] = np[i];
          }
        // Add the particle.
        vec v;
        //cello
        //np1[n1-1].set_index(d.id);
        np1[n1-1].set_index(n1);
        np1[n1-1].set_mass(d.mass);
        np1[n1-1].set_radius(d.radius);
        v[0] = d.x;
        v[1] = d.y;
        v[2] = d.z;
        np1[n1-1].set_pos(v);
        v[0] = d.vx;
        v[1] = d.vy;
        v[2] = d.vz;
        np1[n1-1].set_vel(v);
        
        bhtcs.mass += d.mass;
        bhtcs.n = n1;
        bhtcs.set_nsize(bhtcs.n);
        bhtcs.set_particle_pointer(np1);
        delete [] np;
        
        real pos_scale = 1;
        real vel_scale = 1;
        bhtcs.apply_vf(&real_particle::scale_pos, pos_scale);
        bhtcs.apply_vf(&real_particle::scale_vel, vel_scale);
        //return bhtcs.n;
        *id = bhtcs.n;
        return 0;
        // PRC(pos_scale); PRL(vel_scale);
      }

    
    //return bhtcs.n;
}

static void create_treecode_system()
{
    nbody_particle *np = new nbody_particle[ds.size()];
    bhtcs.n = 0;
    bhtcs.mass = 0;

    for (unsigned int i=0; i<ds.size(); i++) {
        vec v;
        np[i].set_index(ds[i].id);
        np[i].set_mass(ds[i].mass);
        np[i].set_radius(ds[i].radius);
        v[0] = ds[i].x;
        v[1] = ds[i].y;
        v[2] = ds[i].z;
        np[i].set_pos(v);
        v[0] = ds[i].vx;
        v[1] = ds[i].vy;
        v[2] = ds[i].vz;
        np[i].set_vel(v);
        bhtcs.n++;
        bhtcs.mass += ds[i].mass;
    }
    // cerr << "create_treecode_system: "; PRC(bhtcs.n); PRL(bhtcs.mass);

    bhtcs.set_nsize(bhtcs.n);
    bhtcs.set_particle_pointer(np);

    real pos_scale = 1;
    real vel_scale = 1;
    bhtcs.apply_vf(&real_particle::scale_pos, pos_scale);
    bhtcs.apply_vf(&real_particle::scale_vel, vel_scale);
    // PRC(pos_scale); PRL(vel_scale);
}

void print_tree_counters();

// Horrible hack to extract colliders from the tree traversals
// prompted by integrate() below.

static int id_primary = -1, id_secondary = -1;

extern int id_collision_1, id_collision_2;        // determined in BHTC/BHtree.C
extern real r_collision_12;

int get_mass(int id, double *mass)
{
    int i = get_index_from_identity(id);
    
    if (i >= 0 && i < bhtcs.n) 
      {
        nbody_particle *np = bhtcs.get_particle_pointer();
        *mass = np[i].get_mass(); 
        return 0;
      }

    return -1;
}

int get_time(double *time)
{
  *time = bhtcs.time;
  return 0;
}

int set_mass(int id, double m)
{
    int i = get_index_from_identity(id);
    if (i >= 0 && i < bhtcs.n) {
        nbody_particle *np = bhtcs.get_particle_pointer();
        np[i].set_mass(m);
        return 0;
    } else
        return -1;
}

int evolve(real t_end)                // default sync = 0
{

  // Advance from the current time to t_end.

    real dt = bhtcs.timestep;
    real delta_t = t_end - bhtcs.time;
    int nsteps = (int)(delta_t/dt+0.1);
    //    dt = t_end/nsteps;  //<--- bug, should be the following line
    dt = delta_t/nsteps;

    //     PRC(nsteps); PRL(dt);
    //     PRL(bhtcs.n);
    //     PRL(bhtcs.time);
    //     PRL(bhtcs.timestep);
    //     PRL(bhtcs.eps2_for_gravity);
 
    //     PRL(bhtcs.use_self_gravity);
    //     PRL(bhtcs.get_particle_pointer());
    //     for (int j = 0; j < bhtcs.n; j++)
    //         (bhtcs.get_particle_pointer()+j)->dump();

    bhtcs.calculate_gravity();
    real E0 = bhtcs.energy();
    real KE0 = bhtcs.kinetic_energy();
    real PE0 = E0 - KE0;
    real Q0 = KE0/PE0;
    if (debug) 
      {
        PRC(KE0); PRC(PE0); PRC(E0); PRL(Q0);
      }

    for (int i = 0; i < nsteps; i++)
      {

        // bhtcs.time = (i+1)*dt;        // <--- bug part 2
        bhtcs.time += dt;                // this is how it should be

        id_collision_1 = id_collision_2 = id_primary = id_secondary = -1;

        bhtcs.integrate(dt);                // advance the entire system by time dt

        if (id_collision_1 >= 0) 
          {
            id_primary = id_collision_1;
            id_secondary = id_collision_2;        // for pick-up below
            if (debug) 
              {
                cerr << "collision: ";
                PRC(id_primary); PRC(id_secondary); PRL(r_collision_12);
              }
            break;
          }
#if 0
        real KE = bhtcs.kinetic_energy();
        real E = bhtcs.energy();
        real dE = E-E0;
        real Q = KE/(KE-E);
        real t = bhtcs.time;
        PRC(t); PRC(KE); PRC(E); PRC(dE); PRL(Q);
        bhtcs.calculate_cmterms();
        cerr << endl;
#endif

      }

    if (id_primary >= 0) 
      {
        if (debug)
          { 
            print_tree_counters();
          }
        bhtcs.calculate_cmterms();
        real KE = bhtcs.kinetic_energy();
        real E = bhtcs.energy();
        real dE = E-E0;
        real Q = KE/(KE-E);
        if (debug) 
          {
            PRL(bhtcs.time); PRC(KE); PRC(E); PRC(dE); PRL(Q);
            cerr << "CPU sec = " << cpusec() << endl;
          }
      }
    //id_primary
    return 0;
}

int get_index_of_first_particle(int * index_of_the_particle)
{
  *index_of_the_particle =  get_identity_from_index(0);
  return 0;
}

int get_total_radius(double *radius)
{
  return -2;
}

int testme(int *id)
{
  *id = 3;
  return 0;
}

int new_particle(int *id, double mass, double radius, double x, double y, double z, double vx, double vy, double vz) 
{
  dynamics_state state;
  //state.id = id; do this if _new_particle if needed?
  
  state.mass = mass;
  state.radius = radius;
  state.x = x;
  state.y = y;
  state.z = z;
  state.vx = vx;
  state.vy = vy;
  state.vz = vz;
  return _new_particle(id, state);
}

int get_total_mass(double *mass)
{
  return -2;
}

int reinitialize_particles()
{
  return 0;
}

int set_eps2(double epsilon_squared)
{
  return -2;
}

int cleanup_module()
{
  ds.clear();
  initialized = false;
  return 0;
}

int get_eps2(double *epsilon_squared)
{
  return -2;
}

int get_index_of_next_particle(int index_of_the_particle, int *index_of_the_next_particle)
{
  int i = get_index_from_identity(index_of_the_particle);
  if (i==-1)
    {
      //we have an error, so make sure one doesn't get the prev answer
      *index_of_the_next_particle = 0;
      return i;
    }
  else
    {
      int j = get_identity_from_index(i+1);
      if (j>0)
        {
          *index_of_the_next_particle = j;
          return 0;
        }
      else
        {
          *index_of_the_next_particle = 0;
          return 1;
        }
    }
}

int delete_particle(int id)
{
  //used to be "remove_particle"
    if (!initialized) return -1;

    //cerr << "Deleting particle " << id << endl;

    int i = get_index_from_identity(id);
    if (i >= 0 && i < bhtcs.n) 
      {

        // Remove this particle from the tree N-body system
        // by swapping it with the last particle.

        nbody_particle *np = bhtcs.get_particle_pointer();

        if (i < bhtcs.n-1)
          {
            np[i] = np[bhtcs.n-1];
          }
        bhtcs.n--;

        bhtcs.set_nsize(bhtcs.n);

        real pos_scale = 1;
        real vel_scale = 1;
        bhtcs.apply_vf(&real_particle::scale_pos, pos_scale);
        bhtcs.apply_vf(&real_particle::scale_vel, vel_scale);
      }

    //return bhtcs.n;
    return 0;
}

int get_potential(double x, double y, double z, double *V)
{
  return -2;
}

int set_state(int index_of_the_particle, double mass, double radius, double x, double y, double z, double vx, double vy, double vz)
{
  return -2;
}

int get_state(int id, double *mass, double *radius, double *x, double *y, double *z, double *vx, double *vy, double *vz) 
{
    int i = get_index_from_identity(id);
    if (i >= 0 && i < bhtcs.n) 
      {
        nbody_particle *np = bhtcs.get_particle_pointer();
        //*id_out = np[i].get_index();
        *mass = np[i].get_mass();
        *radius = np[i].get_radius();
        vec v = np[i].get_pos();
        *x = v[0];
        *y = v[1];
        *z = v[2];
        v = np[i].get_vel();
        *vx = v[0];
        *vy = v[1];
        *vz = v[2];
        return 0;
      } 
    else 
      {
        return -1;
      }
}

int get_kinetic_energy(double *kinetic_energy)
{
    if (!initialized)
      {
        // some error??
        return 0;
      }
    else
      {
        *kinetic_energy = bhtcs.kinetic_energy();
        return 0;
      }
}

//get number of particles?
int get_number_of_particles(int *number_of_particles)
{
  *number_of_particles = bhtcs.n;
  return 0; 
}

int set_acceleration(int index_of_the_particle, double ax, double ay, double az)
{
  return -2;
}

int get_indices_of_colliding_particles(int *index_of_particle1, int *index_of_particle2)
{
  *index_of_particle1 = find_colliding_primary();
  *index_of_particle2 = find_colliding_secondary(*index_of_particle1);
  return 0;
}

int get_center_of_mass_position(double *x, double *y, double *z)
{
  return -2;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz)
{
  return -2;
}

int get_radius(int id, double *radius)
{
    int i = get_index_from_identity(id);
    if (i >= 0 && i < bhtcs.n) {
        nbody_particle *np = bhtcs.get_particle_pointer();
        *radius = np[i].get_radius();
        return 0;
    }

    return -1;
}

int set_radius(int id, double radius)
{
    int i = get_index_from_identity(id);
    if (i >= 0 && i < bhtcs.n) 
      {
        nbody_particle *np = bhtcs.get_particle_pointer();
        np[i].set_radius(radius);
        return 0;
      }
    else
      {
        return -1;
      }
}

int initialize_code()
{
  return -2;
}

int get_identity_from_index(int i)
{
  if (i >= 0 && i < bhtcs.n)
    {
      nbody_particle *np = bhtcs.get_particle_pointer();
      return np[i].get_index();
    }
  else
    {
      return -1;
    }
}

int get_index_from_identity(int id)
{
  //CHECK IF IN RANGE TO GIVE CORRECT ERROR
  nbody_particle *np = bhtcs.get_particle_pointer();
  for (int i = 0; i < bhtcs.n; i++) 
    {
      if (id == np[i].get_index())
        return i;
    }
  return -1;
}

int initialize_particles()
{
    bhtcs.time = 0.0;
    bhtcs.setup_tree();

    return 0;
}

int get_potential_energy(double *potential_energy)
{
    if (!initialized) 
      {
        return 0;
      }

    bhtcs.calculate_gravity();
    *potential_energy=bhtcs.energy() - bhtcs.kinetic_energy();
    return 0;
}

int get_gravity_at_point(double eps, double x, double y, double z,  double *forcex, double *forcey, double *forcez)
{
    vec p = 0;
    p[0] = x;
    p[1] = y;
    p[2] = z;
    
    
    vec acc = bhtcs.calculate_gravity_at_point(p, bhtcs.eps2_for_gravity, bhtcs.theta_for_tree * bhtcs.theta_for_tree);
    
    *forcex = acc[0];
    *forcey = acc[1];
    *forcez = acc[2];
    
    return 0;
}

int get_potential_at_point(double eps, double x, double y, double z, double * phi)
{
    vec p = 0;
    p[0] = x;
    p[1] = y;
    p[2] = z;
    
    
    *phi  = bhtcs.calculate_potential_at_point(p, bhtcs.eps2_for_gravity, bhtcs.theta_for_tree * bhtcs.theta_for_tree);
    cerr << "phi : "<< *phi << endl;
    return 0;
}

int get_velocity(int id, double *vx, double *vy, double *vz)
{
  return -2;
}

int setup_module()
{
  bhtcs.timestep = timestep;
  bhtcs.eps2_for_gravity = eps2_for_gravity;
  bhtcs.use_self_gravity = use_self_gravity;
  bhtcs.theta_for_tree = theta_for_tree;
  bhtcs.ncrit_for_tree = ncrit_for_tree;

  create_treecode_system();        // note that ds is never used again after this
  initialized = true;

  return 0;
}

int get_position(int id, double *x, double *y, double *z)
{
    int i = get_index_from_identity(id);
    
    if (i >= 0 && i < bhtcs.n) 
    {
        nbody_particle *np = bhtcs.get_particle_pointer();
        vec v = np[i].get_pos();
        *x = v[0];
        *y = v[1];
        *z = v[2];
        return 0;
    }

    return -1;
}

int set_position(int id, double x, double y, double z)
{
    int i = get_index_from_identity(id);
    if (i >= 0 && i < bhtcs.n) 
    {
        nbody_particle *np = bhtcs.get_particle_pointer();
        vec pos(x, y, z);
        np[i].set_pos(pos);
        return 0;
    } 
    else
    {
        return -1;
    }
}

int get_acceleration(int id, double * ax, double * ay, double * az)
{
  return -2;
}

int set_velocity(int id, double vx, double vy, double vz)
{
    int i = get_index_from_identity(id);
    if (i >= 0 && i < bhtcs.n) 
      {
        nbody_particle *np = bhtcs.get_particle_pointer();
        vec vel(vx, vy, vz);
        np[i].set_vel(vel);
        return 0;
      } 
    else
      {
        return -1;
      }
}

double get_dynamical_time_scale()
{
    // PRL(initialized);
    if (!initialized) return 0;

    bhtcs.calculate_gravity();
    real mtot = bhtcs.mass;
    real ekin = bhtcs.kinetic_energy();
    real etot = bhtcs.energy();
    real epot = etot - ekin;

    real tdyn = (-0.5*mtot*mtot/epot) / sqrt(2*ekin/mtot);
    return tdyn;
}

double get_time_step()
{
    //return bhtcs.time + bhtcs.timestep;
    return bhtcs.timestep;
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
    if (i >= 0 && i < bhtcs.n) {
        nbody_particle *np = bhtcs.get_particle_pointer();
        np[i].set_mass(d.mass);
        np[i].set_radius(d.radius);
        np[i].set_pos(vec(d.x, d.y, d.z));
        np[i].set_vel(vec(d.vx, d.vy, d.vz));
        return 0;

    } else {

        cerr << "set_particle: " << id
             << " doesn't exist.  Use add_particle." << endl;
        return -1;
    }
}


int find_colliding_primary()
{
  id_primary = id_secondary = -1;        // defined earlier...

    int n = bhtcs.n;
    nbody_particle *np = bhtcs.get_particle_pointer();
    
    for (int i = 0; i < n-1; i++) {        // very inefficient! -- better to
                                        // flag collisions using the tree
        for (int j = i+1; j < n; j++) {
            real r2 = square(np[i].get_pos() - np[j].get_pos());
            real rsum = np[i].get_radius() + np[j].get_radius();
            // PRC(i); PRC(j); PRC(r2); PRL(rsum);
            if (r2 <= rsum*rsum) {
                if (np[i].get_mass() >= np[j].get_mass()) {
                    id_primary = np[i].get_index();
                    id_secondary = np[j].get_index();
                } else {
                    id_primary = np[j].get_index();
                    id_secondary = np[i].get_index();
                }
                return id_primary;
            }
        }
    }
    return -1;
}

int find_colliding_secondary(int id1)
{
    if (id1 >= 0 && id1 == id_primary && id_secondary >= 0)
        return id_secondary;
    else
        return -1;
}

int get_escaper()
{
    return -1;                                // not implemented yet
}

