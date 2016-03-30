////////////////////////////////////////////////////////
// Includes
////////////////////////////////////////////////////////
#include "Communicator.h"

#include <iostream>
using namespace std;

#include <fstream>
#include <cstdlib>
#include <map>
#include <vector>

#include "Particle.h"
#include "Particles.h"
#include "Sakura.h"
#include "Diagnostics.h"

// AMUSE STOPPING CONDITIONS SUPPORT
#include <stopcond.h>
#include "worker_code.h"

////////////////////////////////////////////////////////
// Declare global variables
////////////////////////////////////////////////////////
int particle_id_counter = 0;
vector<double> data, data_radius;

double t_begin;
double dt;
double t;

/*------------- MPI data ---------------------*/
int mpi_rank = 0;
int mpi_size = 1;
/*------------- MPI data ---------------------*/

////////////////////////////////////////////////////////
// Declare global objects 
////////////////////////////////////////////////////////
string out_directory;
std::map<int, int> local_index_map;

Sakura *sakura = NULL;
Communicator *communicator = NULL;

////////////////////////////////////////////////////////
// Amuse interface functions
////////////////////////////////////////////////////////
int initialize_code() {
#ifndef NOMPI
    int error = 0;
    error = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if(error) {
        cerr << error << endl;
        return -1;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if(error) {
        cerr << error << endl;
        return -1;
    }
#else
    mpi_rank = 0;
    mpi_size = 1;
#endif

    communicator = new Communicator();
    communicator->start_mpi();

    sakura = new Sakura();

    particle_id_counter = 0;
    data.clear();
    data_radius.clear();

    t_begin = 0;
    dt = 1e-3;
    t = t_begin;

#ifndef NOMPI
    mpi_setup_stopping_conditions();
#endif

    return 0;
}

int new_particle_float64(int *particle_identifier, double mass, 
        double x, double y, double z, double vx, double vy, double vz, double radius) {

    data.push_back(mass);
    data.push_back(x);
    data.push_back(y);
    data.push_back(z);
    data.push_back(vx);
    data.push_back(vy);
    data.push_back(vz);

    data_radius.push_back(radius);

    *particle_identifier = particle_id_counter;
    particle_id_counter++;
    return 0;
}
int commit_particles() {
    Particles particles;
    particles.set_t(t_begin);
    particles.set_N(data.size()/7);
    particles.set_data(data);
    sakura->set_particles(particles);

    int numStar = data.size()/7;
    communicator->divide_work(numStar);

    return 0;
}

int set_t_begin(double tb) {
    t_begin = tb;
    return 0;
}
int get_t_begin(double *tb) {
    *tb = t_begin;
    return 0;
}

int set_dt(double ts) {
    dt = ts;
    sakura->set_dt(dt);
    return 0;
}
int get_dt(double *ts) {
    *ts = dt;
    return 0;
}

int set_t(double tt) {
    t = tt;
    return 0;
}
int get_t(double *tt) {
    *tt = t;
    return 0;
}

int commit_parameters() {
    return 0;
}
int recommit_parameters() {
    return commit_parameters();
}

int get_mass(int id, double*mass) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  *mass = data[id*7+0];
  return 0;
} 
int set_mass(int id, double mass) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  data[id*7+0] = mass;
  return 0;
}
int get_position(int id, double* x, double* y, double* z) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  *x = data[id*7+1];
  *y = data[id*7+2];
  *z = data[id*7+3];
  return 0;
}
int set_position(int id, double x, double y, double z) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  data[id*7+1] = x;
  data[id*7+2] = y;
  data[id*7+3] = z;
  return 0;
}
int get_velocity(int id, double* vx, double* vy, double* vz) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  *vx = data[id*7+4];
  *vy = data[id*7+5];
  *vz = data[id*7+6];
  return 0;
}
int set_velocity(int id, double vx, double vy, double vz) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  data[id*7+4] = vx;
  data[id*7+5] = vy;
  data[id*7+6] = vz;
  return 0;
}
int get_state(int id, double* m, double* x, double* y, double* z, double* vx, double* vy, double* vz, double* radius) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  *radius = data_radius[id];
  *m = data[id*7+0];
  *x = data[id*7+1];
  *y = data[id*7+2];
  *z = data[id*7+3];
  *vx = data[id*7+4];
  *vy = data[id*7+5];
  *vz = data[id*7+6];
  return 0;
}
int set_state(int id, double m, double x, double y, double z, double vx, double vy, double vz, double radius) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  data_radius[id] = radius;
  data[id*7+0] = m;
  data[id*7+1] = x;
  data[id*7+2] = y;
  data[id*7+3] = z;
  data[id*7+4] = vx;
  data[id*7+5] = vy;
  data[id*7+6] = vz;
  return 0;
}
int get_radius(int id, double* radius){ 
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  *radius = data_radius[id];
  return 0;
}
int set_radius(int id, double radius) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  data_radius[id] = radius;
  return 0;
}

int evolve_model(double t_end) {
    sakura->set_t(t);
    sakura->update_particles(data);    
    sakura->set_dt(dt);

    sakura->evolve(t_end, *communicator);

    t = t_end;
    data = sakura->get_data();

    return 0;
}
int synchronize_model() {
    return 0;
}
int cleanup_code() {
    delete sakura;
    communicator->stop_mpi();
    delete communicator;
    particle_id_counter = 0;
    return 0;
}

int delete_particle(int id) {
  return -2;
}
int recommit_particles() {
  return -2;
}

int set_begin_time(double input) {
    t_begin = input;
    return 0;
}

int get_begin_time(double * output) {
    *output = t_begin;
    return 0;
}

int get_time(double* time){
    *time = t;
    return 0;
}

int set_eps2(double eps2) {
    return 0;
}
int get_eps2(double *eps2) {
    return 0;
}

int set_sakura_output_directory(char *output_directory){
    out_directory = std::string(output_directory);
    if(out_directory.length() > 0){
        if(*out_directory.rbegin() != '/'){
            out_directory.append("/");
        }
    }
    return 0;
}
int get_sakura_output_directory(char **output_directory){
    *output_directory = (char*) out_directory.c_str();
    return 0;
}
int get_time_step(double* dtt){
    *dtt = dt;
    return 0;
}

int get_potential(int id, double* pot){return -2;}
int get_gravity_at_point(double m, double x, double y, double z, double* rx, double* ry, double* rz){return -2;}
int get_number_of_particles(int* N){return -2;}
int get_potential_at_point(double m, double x, double y, double z, double* p){return -2;}
int get_total_radius(double* R){return -2;}
int get_index_of_first_particle(int* id){return -2;}
int get_index_of_next_particle(int id, int* idnext){return -2;}

int get_total_mass(double* M){ 
  Diagnostics diag;
  *M = diag.get_mass(data);
  return 0;
}
int get_kinetic_energy(double* ek) {
  Diagnostics diag;
  *ek = diag.get_kinetic_energy(data);
  return 0;
}
int get_potential_energy(double* ep) {
  Diagnostics diag;
  *ep = diag.get_potential_energy(data);
  return 0;
}
int get_center_of_mass_position(double* x , double* y, double* z){ 
  Diagnostics diag;
  vector<double> rcm = diag.get_rcm(data);
  *x = rcm[0];
  *y = rcm[1];
  *z = rcm[2];
  return 0;
}
int get_center_of_mass_velocity(double* vx, double* vy, double* vz){ 
  Diagnostics diag;
  vector<double> vcm = diag.get_vcm(data);
  *vx = vcm[0];
  *vy = vcm[1];
  *vz = vcm[2];
  return 0;
}

int get_acceleration(int id, double* ax, double* ay, double* az){return -2;}
int set_acceleration(int id, double ax, double ay, double az){return -2;}

