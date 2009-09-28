#include <iostream>
#include "muse_dynamics.h"
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
real dt_dia = DT_DIA;		// time interval between diagnostics output
void set_dt_dia(real d)			{dt_dia = d;}

//bool x_flag = false;

//int nsteps = 0;		// number of integration time steps completed
//real einit = 0;		// initial total energy of the system
//bool init_flag = false;
//real t_dia = 0;

// Module data:

bool initialized = false;

BHTC_SYSTEM bhtcs;
vector<dynamics_state> ds;	// for initialization only

// External parameters to be used inside the treecode.
// Their definitions are extended to python class modules in the file
// parameters.h; the accessor functions are for the driver program.

const real TIMESTEP = 0.015625;
real timestep =	TIMESTEP;
void set_timestep(real dt)		{timestep = dt;}

const real EPS2_FOR_GRAVITY = 0.125;
real eps2_for_gravity = EPS2_FOR_GRAVITY;
void set_eps2_for_gravity(real eps2)	{eps2_for_gravity = eps2;}

const real THETA_FOR_TREE = 0.75;
real theta_for_tree = THETA_FOR_TREE;
void set_theta_for_tree(real theta)	{theta_for_tree = theta;}

const int USE_SELF_GRAVITY = 1;
int use_self_gravity = USE_SELF_GRAVITY;
void set_use_self_gravity(int use)	{use_self_gravity = use;}

const int NCRIT_FOR_TREE = 1024;
int ncrit_for_tree = NCRIT_FOR_TREE;
void set_ncrit_for_tree(int ncrit)	{ncrit_for_tree = ncrit;}

// Interface functions:

int _add_particle(dynamics_state d)
{
//    cerr << d.id << " " << d.x << " " << d.y << " " << d.z << " "
//	 << d.radius << endl;

    if (!initialized) {

	// Defer reinitialization; save in ds.

	ds.push_back(d);
	return ds.size();

    } else {

	int i = get_index_from_identity(d.id);
	if (i >= 0 && i < bhtcs.n) {

	    // Particle already exists.  Do nothing.

	    cerr << "add_particle: " << d.id
		 << " already exists.  Use set_particle." << endl;

	} else {

	    // Add the new particle and reinitialize immediately.

	    // cerr << "Adding particle " << d.id << endl;
	
	    int n1 = bhtcs.n + 1;
	    nbody_particle *np = bhtcs.get_particle_pointer();
	    nbody_particle *np1 = new nbody_particle[n1];

	    // Copy the system.

	    for (int i = 0; i < bhtcs.n; i++) np1[i] = np[i];

	    // Add the particle.

	    vec v;
	    np1[n1-1].set_index(d.id);
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

	    // PRC(pos_scale); PRL(vel_scale);
	}
    }

    return bhtcs.n;
}

 int add_particle(int id, double mass, double radius, double x, double y, double z, double vx, double vy, double vz) {
        dynamics_state state;
        state.id = id;
        state.mass = mass;
        state.radius = radius;
        state.x = x;
        state.y = y;
        state.z = z;
        state.vx = vx;
        state.vy = vy;
        state.vz = vz;
        return _add_particle(state);
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

extern int id_collision_1, id_collision_2;	// determined in BHTC/BHtree.C
extern real r_collision_12;

int evolve(real t_end, int sync)		// default sync = 0
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
    // 	(bhtcs.get_particle_pointer()+j)->dump();

    bhtcs.calculate_gravity();
    real E0 = bhtcs.energy();
    real KE0 = bhtcs.kinetic_energy();
    real PE0 = E0 - KE0;
    real Q0 = KE0/PE0;
    if (debug) {
      PRC(KE0); PRC(PE0); PRC(E0); PRL(Q0);
    }

    for (int i = 0; i < nsteps; i++){

	// bhtcs.time = (i+1)*dt;	// <--- bug part 2
	bhtcs.time += dt;		// this is how it should be

	id_collision_1 = id_collision_2 = id_primary = id_secondary = -1;

	bhtcs.integrate(dt);		// advance the entire system by time dt

	if (id_collision_1 >= 0) {
	    id_primary = id_collision_1;
	    id_secondary = id_collision_2;	// for pick-up below
	    if (debug) {
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

    if (id_primary >= 0) {
      if (debug) print_tree_counters();
	bhtcs.calculate_cmterms();
	real KE = bhtcs.kinetic_energy();
	real E = bhtcs.energy();
	real dE = E-E0;
	real Q = KE/(KE-E);
	if (debug) {
	  PRL(bhtcs.time); PRC(KE); PRC(E); PRC(dE); PRL(Q);
	  cerr << "CPU sec = " << cpusec() << endl;
	}
    }

    return id_primary;
}

int get_number()
{
    return bhtcs.n;
}

int get_index_from_identity(int id)
{
    nbody_particle *np = bhtcs.get_particle_pointer();
    for (int i = 0; i < bhtcs.n; i++) {
	if (id == np[i].get_index())
	    return i;
    }
  return -1;
}

int remove_particle(int id)
{
    if (!initialized) return -1;

    cerr << "Deleting particle " << id << endl;

    int i = get_index_from_identity(id);
    if (i >= 0 && i < bhtcs.n) {

	// Remove this particle from the tree N-body system
	// by swapping it with the last particle.

	nbody_particle *np = bhtcs.get_particle_pointer();

	if (i < bhtcs.n-1) np[i] = np[bhtcs.n-1];
	bhtcs.n--;

	bhtcs.set_nsize(bhtcs.n);

	real pos_scale = 1;
	real vel_scale = 1;
	bhtcs.apply_vf(&real_particle::scale_pos, pos_scale);
	bhtcs.apply_vf(&real_particle::scale_vel, vel_scale);
    }

    return bhtcs.n;
}

double get_potential_energy()
{
    if (!initialized) return 0;

    bhtcs.calculate_gravity();
    return bhtcs.energy() - bhtcs.kinetic_energy();
}

double get_kinetic_energy()
{
    if (!initialized) return 0;

    return bhtcs.kinetic_energy();
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

double get_time()
{
    return bhtcs.time;
}

void get_state(int id, int * id_out,  double * mass, double * radius, double * x, double * y, double * z, double * vx, double * vy, double * vz) {
    int i = get_index_from_identity(id);
    if (i >= 0 && i < bhtcs.n) {
        nbody_particle *np = bhtcs.get_particle_pointer();
        *id_out = np[i].get_index();
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
    } else {
        *id_out = -1;
    }
}

double get_time_step()
{
    //return bhtcs.time + bhtcs.timestep;
    return bhtcs.timestep;
}

int setup_module()
{
  bhtcs.timestep = timestep;
  bhtcs.eps2_for_gravity = eps2_for_gravity;
  bhtcs.use_self_gravity = use_self_gravity;
  bhtcs.theta_for_tree = theta_for_tree;
  bhtcs.ncrit_for_tree = ncrit_for_tree;

  create_treecode_system();	// note that ds is never used again after this
  initialized = true;

  return 0;
}

int cleanup_module()
{
  ds.clear();
  initialized = false;
  return 0;
}

int initialize_particles(real t)
{
  bhtcs.time = t;
//   PRL(bhtcs.eps2_for_gravity);
//   for (int id = 1; id < 10; id++) {
//     dynamics_state s = get_state(id);
//     real m = s.mass;
//     vec pos = vec(s.x, s.y, s.z);
//     PRC(id); PRC(m); PRL(pos);
//   }
  return 0;
}

int reinitialize_particles()
{
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

int set_radius(int id, double r)
{
    int i = get_index_from_identity(id);
    if (i >= 0 && i < bhtcs.n) {
	nbody_particle *np = bhtcs.get_particle_pointer();
	np[i].set_radius(r);
	return 0;
    } else
	return -1;
}

int set_pos(int id, double x[])
{
    int i = get_index_from_identity(id);
    if (i >= 0 && i < bhtcs.n) {
	nbody_particle *np = bhtcs.get_particle_pointer();
	vec pos(x[0], x[1], x[2]);
	np[i].set_pos(pos);
	return 0;
    } else
	return -1;
}

int set_vel(int id, double v[])
{
    int i = get_index_from_identity(id);
    if (i >= 0 && i < bhtcs.n) {
	nbody_particle *np = bhtcs.get_particle_pointer();
	vec vel(v[0], v[1], v[2]);
	np[i].set_vel(vel);
	return 0;
    } else
	return -1;
}

double get_mass(int id)
{
    int i = get_index_from_identity(id);
    if (i >= 0 && i < bhtcs.n) {
	nbody_particle *np = bhtcs.get_particle_pointer();
	return np[i].get_mass();
    }

    return -1;
}

double get_radius(int id)
{
    int i = get_index_from_identity(id);
    if (i >= 0 && i < bhtcs.n) {
	nbody_particle *np = bhtcs.get_particle_pointer();
	return np[i].get_radius();
    }

    return -1;
}

int find_colliding_primary()
{
  id_primary = id_secondary = -1;	// defined earlier...

    int n = bhtcs.n;
    nbody_particle *np = bhtcs.get_particle_pointer();
    
    for (int i = 0; i < n-1; i++) {	// very inefficient! -- better to
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
    return -1;				// not implemented yet
}

