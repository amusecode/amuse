#include "src/include/octree.h"



#include <string>
#include <iostream>
#include <fstream>
#include "src/include/octree.h"


octree *bonsai;
bool initialized = false;
string bonsai_source_directory = "./";

long long my_dev::base_mem::currentMemUsage;
long long my_dev::base_mem::maxMemUsage;

vector<float4> bodies_pos;
vector<float4> bodies_vel;
vector<float4> bodies_grav(0);   // (w = potential)
vector<int>    starids(0);       // list of identifiers
vector<float>  radii(0);         // list of radii

int id_counter          = 0;
int n_bodies            = 0;
double total_mass       = 0;
double t_now		= 0.0;

bool curStateOnHost	= false;

double timestep;

std::string logFileName = "bonsaiLog.txt";
std::ofstream logFile;

//Helper functions

int getCurrentStateToHost()
{
  if(!curStateOnHost)
  {
     //Retrieve the current state from the device
     bonsai->desort_bodies(bonsai->localTree);

     bonsai->localTree.bodies_pos.d2h();
     bonsai->localTree.bodies_vel.d2h();
     bonsai->localTree.bodies_ids.d2h();
     bonsai->localTree.bodies_acc1.d2h();

     curStateOnHost = true;
  }
  return 0;
}

//End helper functions

//Initialise the code
int initialize_code()
{
  int devID = 0;
  float theta = 0.75;
  float eps = 0.05;
  std::string snapshotFile = "test";
  int snapshotIter = -1;
  float timeStep =  1.0 / 64;
  float tEnd = 10.0;

  logFile.open(logFileName.c_str());
 //Creat the octree class and set the properties
  bonsai = new octree(devID, theta, eps, snapshotFile, snapshotIter, timeStep, tEnd);
  bonsai->set_src_directory(bonsai_source_directory);
  bonsai->set_context(logFile, false); //Do logging to file and enable timing (false = enabled)
  bonsai->load_kernels();

  initialized = true;

  return 0;
}


int echo(int input)
{
	initialize_code();
	return input;
}

int set_src_directory(char * src_dir)
{
    bonsai_source_directory.assign(src_dir);
    return 0;
}


// Interface functions:
int new_particle(int *id, double mass, double radius, double x, double y, double z, double vx, double vy, double vz)
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

int delete_particle(int index_of_the_particle)
{
  //In order to delete a particle we have to reinitialize
  //the particle arrays, first copy the properties to vectors
  getCurrentStateToHost();

  //Resize should not be needed
  bodies_pos.resize(n_bodies);
  bodies_vel.resize(n_bodies);
  bodies_grav.resize(n_bodies);
  starids.resize(n_bodies);
  radii.resize(n_bodies);

  //From code into std::vectors
  memcpy(&bodies_pos[0], &bonsai->localTree.bodies_pos[0], sizeof(float4)*n_bodies);
  memcpy(&bodies_vel[0], &bonsai->localTree.bodies_vel[0], sizeof(float4)*n_bodies);
  memcpy(&bodies_grav[0], &bonsai->localTree.bodies_acc1[0], sizeof(float4)*n_bodies);
  memcpy(&starids[0], &bonsai->localTree.bodies_ids[0], sizeof(float4)*n_bodies);

  int i = index_of_the_particle;

  if (i >= 0 && i < n_bodies)
  {
      total_mass -= bodies_pos[i].w;
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


int commit_particles()
{
  assert(initialized == true);

  bonsai->localTree.setN(n_bodies);
  bonsai->allocateParticleMemory(bonsai->localTree);

  //Load data into the host buffers
  for(int i=0; i < n_bodies; i++)
  {
    bonsai->localTree.bodies_pos[i] = bodies_pos[i];
    bonsai->localTree.bodies_vel[i] = bodies_vel[i];
    bonsai->localTree.bodies_ids[i] = starids[i];

    bonsai->localTree.bodies_Ppos[i] = bodies_pos[i];
    bonsai->localTree.bodies_Pvel[i] = bodies_vel[i];
  }

  //Copy the particles to the device
  bonsai->localTree.bodies_pos.h2d();
  bonsai->localTree.bodies_vel.h2d();
  bonsai->localTree.bodies_Ppos.h2d();
  bonsai->localTree.bodies_Pvel.h2d();
  bonsai->localTree.bodies_ids.h2d();

  //Build a tree-structure for initial initialization
  bonsai->sort_bodies(bonsai->localTree);
  bonsai->build(bonsai->localTree);
  bonsai->allocateTreePropMemory(bonsai->localTree);
  bonsai->compute_properties(bonsai->localTree);

  return 0;
}


int commit_parameters()
{
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return 0;
}


int get_state(int id, double *mass, double *radius, double *x, double *y, double *z, double *vx, double *vy, double *vz)
{
  assert(initialized == true);

  getCurrentStateToHost();

  int i = id;
  if (i >= 0 && i < n_bodies)
  {
    *mass   = bonsai->localTree.bodies_pos[i].w;
    *radius = radii[i];

    *x   = bonsai->localTree.bodies_pos[i].x;
    *y   = bonsai->localTree.bodies_pos[i].y;
    *z   = bonsai->localTree.bodies_pos[i].z;

    *vx   = bonsai->localTree.bodies_vel[i].x;
    *vy   = bonsai->localTree.bodies_vel[i].y;
    *vz   = bonsai->localTree.bodies_vel[i].z;
    return 0;
  }
  else
  {
    return -1;
  }
}



int evolve_model(double t_end)
{
  // advance from the current time to t_end
  /* 1. Create equal-sized time steps by adjusting timestep */

//  double delta_t = t_end - t_now;
 // int nsteps;
  //  dtime = timestep;
 // nsteps = (int) (delta_t/timestep)+1;
  //double dtime = delta_t/nsteps;


  curStateOnHost = false;

  /* 2. Calculate energies */
  ///calculate the inital energies

//  double E_kin = calcEkin(bodies_pos, bodies_vel);
//  double E_pot = calcEpot(bodies_pos, bodies_grav);

  /* 4. Integrate as long as necessary */
//  system.set_softening(eps);
//  system.set_opening_angle(theta);


  bonsai->setTEnd(t_end);
  bonsai->iterate();

  t_now = t_end;

  return 0;
}

//Parameters

int get_time_step(double *_timestep)
{
  *_timestep = bonsai->getDt();
  return 0;
}

int set_time_step(double _timestep)
{
  bonsai->setDt(timestep);
  return 0;
}

int get_eps2(double *epsilon_squared)
{
  double  eps = bonsai->getEps();
  *epsilon_squared = eps*eps;
  return 0;
}


int set_eps2(double epsilon_squared)
{
  double eps = sqrt(epsilon_squared);
  bonsai->setEps(eps);
  return 0;
}

int set_time(double time_now)
{
  bonsai->setTime(time_now);
  t_now = time_now;
  return 0;
}

int get_time(double *time)
{
  *time = bonsai->getTime();
  t_now = *time;
  return 0;
}


int set_theta_for_tree(double theta_for_tree)
{
  bonsai->setTheta(theta_for_tree);
  return 0;
}

int get_theta_for_tree(double *theta_for_tree)
{
  *theta_for_tree = bonsai->getTheta();
  return 0;
}
//End parameters


int recommit_parameters(){
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return 0;
}

int synchronize_model(){
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return 0;
}

int recommit_particles(){
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return 0;
}


//Particle property functions
int get_mass(int index_of_the_particle, double * mass){
  getCurrentStateToHost();
  *mass = bonsai->localTree.bodies_pos[index_of_the_particle].w;
  return 0;
}

int set_mass(int index_of_the_particle, double mass){
  getCurrentStateToHost();
  bonsai->localTree.bodies_pos[index_of_the_particle].w = mass;
  return 0;
}

int get_radius(int index_of_the_particle, double * radius){
  *radius = radii[index_of_the_particle];
  return 0;
}

int set_radius(int index_of_the_particle, double radius){
  radii[index_of_the_particle] = radius;
  return 0;
}

int get_potential(int index_of_the_particle, double * potential){
  getCurrentStateToHost();
  *potential = bonsai->localTree.bodies_acc1[index_of_the_particle].w;
  return 0;
}

int set_acceleration(int index_of_the_particle, double ax, double ay,
  double az){
  getCurrentStateToHost();
  bonsai->localTree.bodies_acc1[index_of_the_particle].x = ax;
  bonsai->localTree.bodies_acc1[index_of_the_particle].x = ay;
  bonsai->localTree.bodies_acc1[index_of_the_particle].x = az;
  return 0;
}

int get_acceleration(int index_of_the_particle, double * ax, double * ay,
  double * az){
  getCurrentStateToHost();
  *ax = bonsai->localTree.bodies_acc1[index_of_the_particle].x;
  *ay = bonsai->localTree.bodies_acc1[index_of_the_particle].y;
  *az = bonsai->localTree.bodies_acc1[index_of_the_particle].z;
  return 0;
}

int get_velocity(int index_of_the_particle, double * vx, double * vy,
  double * vz){
  getCurrentStateToHost();
  *vx = bonsai->localTree.bodies_vel[index_of_the_particle].x;
  *vy = bonsai->localTree.bodies_vel[index_of_the_particle].y;
  *vz = bonsai->localTree.bodies_vel[index_of_the_particle].z;
  return 0;
}

int get_position(int index_of_the_particle, double * x, double * y,
  double * z){
  getCurrentStateToHost();
  *x = bonsai->localTree.bodies_pos[index_of_the_particle].x;
  *y = bonsai->localTree.bodies_pos[index_of_the_particle].y;
  *z = bonsai->localTree.bodies_pos[index_of_the_particle].z;
  return 0;
}

int set_position(int index_of_the_particle, double x, double y, double z){
  getCurrentStateToHost();
  bonsai->localTree.bodies_pos[index_of_the_particle].x = x;
  bonsai->localTree.bodies_pos[index_of_the_particle].x = y;
  bonsai->localTree.bodies_pos[index_of_the_particle].x = z;
  return 0;
}

int set_velocity(int index_of_the_particle, double vx, double vy,
  double vz){
  getCurrentStateToHost();
  bonsai->localTree.bodies_vel[index_of_the_particle].x = vx;
  bonsai->localTree.bodies_vel[index_of_the_particle].x = vy;
  bonsai->localTree.bodies_vel[index_of_the_particle].x = vz;
  return 0;
}

int set_state(int index_of_the_particle, double mass, double radius,
  double x, double y, double z, double vx, double vy, double vz){
  set_mass(index_of_the_particle, mass);
  set_radius(index_of_the_particle, radius);
  set_position(index_of_the_particle, x, y, z);
  set_velocity(index_of_the_particle, x, y, z);
  return 0;
}

//End particle properties

//Simulation system properties, stats, etc

int get_total_radius(double * radius){
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return 0;
}

int get_total_mass(double * mass){
  *mass = total_mass;
  return 0;
}

int get_number_of_particles(int *number){
 *number = n_bodies;
 return 0;
}


//End Simulation system properties


int get_kinetic_energy(double * kinetic_energy){
  *kinetic_energy = bonsai->getKin();
  return 0;
}

int get_potential_energy(double * potential_energy){
  *potential_energy = bonsai->getPot();
  return 0;
}


int get_index_of_first_particle(int * index_of_the_particle){
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return 0;
}
int get_index_of_next_particle(int index_of_the_particle,
                int * index_of_the_next_particle){
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return 0;
}

int get_potential_at_point(double eps, double x, double y, double z,
  double * phi){
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return 0;
}




int get_center_of_mass_position(double * x, double * y, double * z){
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return 0;
}


int get_gravity_at_point(double eps, double x, double y, double z,
  double * forcex, double * forcey, double * forcez){
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return 0;
}




int cleanup_code(){
  delete bonsai;
  bonsai = NULL;
  return 0;
}


