#include <cmath>
#include <map>

extern "C" {
#include "src/universal_kepler_solver.h"
}

extern "C" {
#include "src/universal_variable_kepler.h"
}

#include "interface.h"

using namespace std;


typedef struct {
    double mass;                                        /// mass
    double x, y, z;                                     /// position
    double vx, vy, vz;                                  /// velocity
} dynamics_state;

static double simtime;
static map<int, dynamics_state*> kmap;
static int kcounter=1;

static double central_mass=0.;
static double central_radius=0.;
static double central_x=0.,central_y=0.,central_z=0.;
static double central_vx=0.,central_vy=0.,central_vz=0.;
static int central_id=-1;

static double eps2=0;

static int method=0;

int initialize_code()
{
    return 0;
}

int commit_parameters()
{
    return 0;
}

int commit_particles()
{
    return 0;
}

int recommit_particles()
{
    return 0;
}

int synchronize_model()
{
    return 0;
}

int recommit_parameters()
{
    return 0;
}

int cleanup_code()
{
    return 0;
}

int set_acceleration(int particle_identifier, double ax, double ay, double az) {
    return -2; // Not implemented
}

int get_acceleration(int particle_identifier, double *ax, double *ay, double *az) {
    return -2; // Not implemented
}
int get_potential(int id,double *p){
    return -2; // Not implemented
}

int set_central_mass(int id, double m)
{
  if(id!=0) return -1;
  central_mass=m;
  return 0;
}
int get_central_mass(int id, double *m)
{ 
  if(id!=0) return -1;
  *m=central_mass;
  return 0;
}

int set_central_radius(int id, double m)
{
  if(id!=0) return -1;
  central_radius=m;
  return 0;
}
int get_central_radius(int id, double *m)
{ 
  if(id!=0) return -1;
  *m=central_radius;
  return 0;
}


int set_central_pos(int id,double x,double y,double z)
{ 
  if(id!=0) return -1;
  central_x=x;
  central_y=y;
  central_z=z;
  return 0;
}
int get_central_pos(int id,double *x,double *y,double *z)
{
  if(id!=0) return -1;
  *x=central_x;
  *y=central_y;
  *z=central_z;
  return 0;
}

int set_central_vel(int id, double vx,double vy,double vz)
{ 
  if(id!=0) return -1;
  central_vx=vx;
  central_vy=vy;
  central_vz=vz;
  return 0;
}
int get_central_vel(int id,double *vx,double *vy,double *vz)
{ 
  if(id!=0) return -1;
  *vx=central_vx;
  *vy=central_vy;
  *vz=central_vz;
  return 0;
}

int new_central_particle(int *id, double mass,
			double x, double y, double z,
			double vx, double vy, double vz,double r)
{
    if(central_id!=-1) return -1;
    central_id=0;
    *id=0;
    central_mass=mass;
    central_x=x;
    central_y=y;
    central_z=z;
    central_vx=vx;
    central_vy=vy;
    central_vz=vz;
    central_radius=r;
    return 0;
}

int delete_central_particle(int id)
{
    if(id!=0) return -1;
    central_id=-1;
    return 0;
}

int set_eps2(double m){ eps2=m;return 0;}
int get_eps2(double *m){ *m=eps2;return 0;}

int set_method(int m){ method=m;return 0;}
int get_method(int *m){ *m=method;return 0;}

int get_time_step(double *m){ *m=0.;return 0;}

int new_particle(int *id, double mass,
			double x, double y, double z,
			double vx, double vy, double vz,double r)
{
    dynamics_state *k = new dynamics_state;
    k->mass=mass;
    k->x=x;
    k->y=y;
    k->z=z;
    k->vx=vx;
    k->vy=vy;
    k->vz=vz;
    *id=kcounter++;
   	kmap.insert(kmap.end(), std::pair<int, dynamics_state*>( *id, k));
    return 0;
}

int delete_particle(int id)
{
    map<int, dynamics_state*>::iterator iter = kmap.find(id);
    
    if (iter != kmap.end()){
        delete (*iter).second;
        kmap.erase(iter);
    } else {
      return -3;
    }
    return 0;
}

int get_state(int id, double *m, double *x, double *y, double *z,
                 double *vx, double *vy, double *vz,double *rad)
{
    map<int, dynamics_state*>::iterator iter = kmap.find(id);
    
    if (iter != kmap.end()){
        dynamics_state *k=iter->second;
        *m=k->mass;
        *x = k->x;
        *y = k->y;
        *z = k->z;
        *vx = k->vx;
        *vy = k->vy;
        *vz = k->vz;
        *rad = 0.;
    } else {
      return -3;
    }
    return 0;
}

int set_state(int id, double m, double x, double y, double z,
                 double vx, double vy, double vz,double r)
{
    map<int, dynamics_state*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
        dynamics_state *k=iter->second;
        k->mass=m;
        k->x=x;
        k->y=y;
        k->z=z;
        k->vx=vx;
        k->vy=vy;
        k->vz=vz;
    } else {
      return -3;
    }
    return 0;
}

int set_mass(int id, double m)
{
    map<int, dynamics_state*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
        dynamics_state *k=iter->second;
        k->mass=m;
    } else {
      return -3;
    }
    return 0;
}

int set_radius(int id, double r)
{
    map<int, dynamics_state*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
    } else {
      return -3;
    }
    return 0;
}

int get_radius(int id, double *r)
{
    map<int, dynamics_state*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
      *r=0;
    } else {
      return -3;
    }
    return 0;
}

int set_position(int id, double x, double y, double z)
{
    map<int, dynamics_state*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
        dynamics_state *k=iter->second;
        k->x=x;
        k->y=y;
        k->z=z;
    } else {
      return -3;
    }
    return 0;
}

int set_velocity(int id,double vx, double vy, double vz)
{
    map<int, dynamics_state*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
        dynamics_state *k=iter->second;
        k->vx=vx;
        k->vy=vy;
        k->vz=vz;
    } else {
      return -3;
    }
    return 0;
}



int evolve_model(double tend)
{
    map<int, dynamics_state*>::iterator iter;
    int *keys= new int[kmap.size()];
    int i=0;
    for(iter=kmap.begin(); iter!=kmap.end(); iter++)
    {
      keys[i]=iter->first;
      i++;
    }

    double new_x=central_x;
    double new_y=central_y;
    double new_z=central_z;

    new_x+=central_vx*(tend-simtime);
    new_y+=central_vy*(tend-simtime);
    new_z+=central_vz*(tend-simtime);

#pragma omp parallel for
    for(i=0;i<kmap.size();i++)
    {
        dynamics_state *d=kmap[keys[i]];
     
        double pos0[3], vel0[3];
        double pos[3],vel[3];
        double dt, mu;
        
        pos0[0] = d->x;
        pos0[1] = d->y;
        pos0[2] = d->z;
        
        vel0[0] = d->vx;
        vel0[1] = d->vy;
        vel0[2] = d->vz;
        
        dt = tend-simtime;
        mu = central_mass + d->mass;
        
        if(method)
             universal_variable_kepler_solver(dt, mu, pos0, vel0, pos, vel);
        else
        universal_kepler_solver(dt, mu, eps2, pos0[0], pos0[1], pos0[2], vel0[0], vel0[1], vel0[2], 
          &pos[0], &pos[1], &pos[2], &vel[0], &vel[1], &vel[2]);
        
        d->x = pos[0]+new_x;
        d->y = pos[1]+new_y;
        d->z = pos[2]+new_z;
        d->vx = vel[0]+central_vx;
        d->vy = vel[1]+central_vy;
        d->vz = vel[2]+central_vz;
    }

    central_x=new_x;
    central_y=new_y;
    central_z=new_z;

    simtime=tend;
    
    delete[] keys; 

    return 0;
}

int get_mass(int id, double *mass)
{
    map<int, dynamics_state*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
       dynamics_state *k=iter->second;
       *mass=k->mass;
    } else {
      return -3;
    }
    return 0;
}

int get_time(double * t)
{
    *t = simtime;
    return 0;
}

int get_position(int id,double * x, double * y, double * z)
{
    map<int, dynamics_state*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
       dynamics_state *k=iter->second;
       *x = k->x;
       *y = k->y;
       *z = k->z;
    } else {
      return -3;
    }
    return 0;
}

int get_velocity(int id, double * vx, double * vy, double * vz)
{
    map<int, dynamics_state*>::iterator iter = kmap.find(id);
    if (iter != kmap.end()){
       dynamics_state *k=iter->second;
       *vx = k->vx;
       *vy = k->vy;
       *vz = k->vz;
    } else {
      return -3;
    }
    return 0;
}


int get_index_of_first_particle(int *particle_identifier) {
    return -2; // Not implemented
}
int get_index_of_next_particle(int particle_identifier, int *next_particle_identifier) {
    return -2; // Not implemented
}

int get_indices_of_colliding_particles(int *index_of_particle1, int *index_of_particle2) {
    return -2; // Not implemented
}

// simulation property getters:
int get_total_mass(double *total_mass) {
    return -2; // Not implemented
}
int get_total_radius(double *total_radius) {
    return -2; // Not implemented
}

int set_begin_time(double input) {
    simtime=input;
    return 0;
}

int get_begin_time(double * output) {
    return -2; // Not implemented
}

int get_center_of_mass_position(double *x, double *y, double *z){
    return -2; // Not implemented
}
int get_center_of_mass_velocity(double *vx, double *vy, double *vz) {
    return -2; // Not implemented
}
int get_kinetic_energy(double *kinetic_energy) {
    return -2; // Not implemented
}
int get_potential_energy(double *potential_energy) {
    return -2; // Not implemented
}
int get_number_of_particles(int *number_of_particles) {
    *number_of_particles = kmap.size();
    return 0;
}


