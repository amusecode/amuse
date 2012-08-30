#include "src/include/octree.h"


#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include "src/include/octree.h"

#ifdef _AMUSE_STOPPING_CONDITIONS_
  // AMUSE STOPPING CONDITIONS SUPPORT
  #include <stopcond.h>
#endif  


octree *bonsai;
bool initialized = false;
string bonsai_source_directory = "./";
bool maxlevels_exceeded = false;

long long my_dev::base_mem::currentMemUsage;
long long my_dev::base_mem::maxMemUsage;

vector<float4> bodies_pos;
vector<float4> bodies_vel;
vector<double2> bodies_time;
vector<float4> bodies_grav(0);   // (w = potential)
vector<int>    starids(0);       // list of identifiers
vector<float>  radii(0);         // list of radii

int id_counter          = 0;
int n_bodies            = 0;
double total_mass       = 0;
double t_now		= 0.0;
double timestep_h	= 0.0;
static double begin_time = 0;

bool curStateOnHost	= false;


std::map<int, int> idToIndex;
std::map<int, int>::iterator iTIIter;

std::string logFileName = "bonsaiLog.txt";
std::ofstream logFile;

//Helper functions

int getCurrentStateToHost()
{
  if(!curStateOnHost)
  {
     //Retrieve the current state from the device
  //   bonsai->desort_bodies(bonsai->localTree);

     bonsai->localTree.bodies_pos.d2h();
     bonsai->localTree.bodies_vel.d2h();
     bonsai->localTree.bodies_ids.d2h();
     bonsai->localTree.bodies_acc1.d2h();
     bonsai->localTree.bodies_time.d2h();

     //From code into std::vectors
     memcpy(&bodies_pos[0], &bonsai->localTree.bodies_pos[0], sizeof(float4)*n_bodies);
     memcpy(&bodies_vel[0], &bonsai->localTree.bodies_vel[0], sizeof(float4)*n_bodies);
     memcpy(&bodies_grav[0], &bonsai->localTree.bodies_acc1[0], sizeof(float4)*n_bodies);
     memcpy(&bodies_time[0], &bonsai->localTree.bodies_time[0], sizeof(double2)*n_bodies);
     memcpy(&starids[0], &bonsai->localTree.bodies_ids[0], sizeof(int)*n_bodies);

	
     //Build the map to index
     idToIndex.clear();
     for(int i=0; i < n_bodies; i++)
     {
        idToIndex[starids[i]] = i;
     }

     curStateOnHost = true;
  }
  return 0;
}


int getIdxFromId(int id)
{
  iTIIter = idToIndex.find(id);
  if(iTIIter != idToIndex.end())
  {
    return (*iTIIter).second;
  }
  else
  {
    fprintf(stderr, "Request for unknown particle-id: %d \n", id);
    return -1;
  }
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
  char **argv = NULL;
  bonsai = new octree(argv,devID, theta, eps, snapshotFile, snapshotIter, timeStep, tEnd);
  bonsai->set_src_directory(bonsai_source_directory);
  bonsai->set_context(logFile, false); //Do logging to file and enable timing (false = enabled)
  bonsai->load_kernels();

  initialized = true;
  curStateOnHost = true;
  
  #ifdef _AMUSE_STOPPING_CONDITIONS_
    // AMUSE STOPPING CONDITIONS SUPPORT
    set_support_for_condition(COLLISION_DETECTION);  
  #endif
  

  return 0;
}

int set_src_directory(char * src_dir)
{
    bonsai_source_directory.assign(src_dir);
    return 0;
}


// Interface functions:
int new_particle(int *id, double mass, double x, double y, double z, double vx, double vy, double vz, double radius)
{
  *id = id_counter;

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
  //bodies_vel[n_bodies].w = 0;
  bodies_vel[n_bodies].w = radius; //Store radius in 'w' component for easy access stopping conditions

  bodies_time.resize(n_bodies+1);
  bodies_time[n_bodies].x = t_now;

  double tempTimeStep =  bonsai->getDt();

  //TODO change this into probably t_now to make sure it will do a step to get proper dt in block step mode
  bodies_time[n_bodies].y = t_now + tempTimeStep;


  bodies_grav.resize(n_bodies+1);

  n_bodies++;

  starids.push_back(id_counter);
  radii.push_back(radius);

  id_counter++;
  total_mass += mass;

  return 0;
}

int delete_particle(int id)
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
  bodies_time.resize(n_bodies);

  //From code into std::vectors
  //memcpy(&bodies_pos[0], &bonsai->localTree.bodies_pos[0], sizeof(float4)*n_bodies);
  //memcpy(&bodies_vel[0], &bonsai->localTree.bodies_vel[0], sizeof(float4)*n_bodies);
  //memcpy(&bodies_grav[0], &bonsai->localTree.bodies_acc1[0], sizeof(float4)*n_bodies);
  //memcpy(&starids[0], &bonsai->localTree.bodies_ids[0], sizeof(int)*n_bodies);
  
  int index = getIdxFromId(id);
  if(index < 0) return -3;  

 
  if(id == starids[index])
  {
    total_mass -= bodies_pos[index].w;
    bodies_pos.erase(bodies_pos.begin()+index);
    bodies_vel.erase(bodies_vel.begin()+index);
    bodies_grav.erase(bodies_grav.begin()+index);
    bodies_time.erase(bodies_time.begin()+index);
    starids.erase(starids.begin()+index);
    radii.erase(radii.begin()+index);
    n_bodies--;

    //Have to rebuild the idx to id since it has become
    //invalid since we removed something
    idToIndex.clear();
    for(int i=0; i < n_bodies; i++) {
      idToIndex[starids[i]] = i;
    }

    return 0;
  }     
  else
  {
    //The id array is out of sync??
    return -1;
  }

}






int commit_particles()
{
  assert(initialized == true);
  maxlevels_exceeded = false;
  
  idToIndex.clear();

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
  bonsai->sort_bodies(bonsai->localTree, true);
  bonsai->build(bonsai->localTree);
  if (maxlevels_exceeded) return -4;
  bonsai->allocateTreePropMemory(bonsai->localTree);
  bonsai->compute_properties(bonsai->localTree);


  curStateOnHost = false;

  return 0;
}


int commit_parameters()
{
    bonsai->setTime(begin_time);
    t_now = begin_time;
    return 0;
}


int get_state(int id, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *radius)
{
  assert(initialized == true);

  getCurrentStateToHost();

  int i = getIdxFromId(id);
  if(i < 0) return -3;


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
    return -3;
  }
}



int evolve_model(double t_end)
{
  maxlevels_exceeded = false;
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
  if (maxlevels_exceeded) return -4;

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
  bonsai->setDt(_timestep);
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
  //return 0;
 

  //Recommit also copies the acc and time from host to device
  //while those are not needed/available in the commit_particles

  assert(initialized == true);
  maxlevels_exceeded = false;
  
  idToIndex.clear();

  bonsai->localTree.setN(n_bodies);
  bonsai->allocateParticleMemory(bonsai->localTree);

  //Load data into the host buffers
  for(int i=0; i < n_bodies; i++)
  {
    bonsai->localTree.bodies_pos[i] = bodies_pos[i];
    bonsai->localTree.bodies_vel[i] = bodies_vel[i];
    bonsai->localTree.bodies_ids[i] = starids[i];
    bonsai->localTree.bodies_time[i] = bodies_time[i];
    bonsai->localTree.bodies_ids[i] = starids[i];
    bonsai->localTree.bodies_acc1[i] = bodies_grav[i];

    bonsai->localTree.bodies_Ppos[i] = bodies_pos[i];
    bonsai->localTree.bodies_Pvel[i] = bodies_vel[i];
   
    //Force the time to t_now to make sure all particles
    //will take one step after the particles have been
    //recommitted. 
    bonsai->localTree.bodies_time[i].x = t_now; 
    bonsai->localTree.bodies_time[i].y = t_now; 
  }

  //Copy the particles to the device
  bonsai->localTree.bodies_pos.h2d();
  bonsai->localTree.bodies_vel.h2d();
  bonsai->localTree.bodies_Ppos.h2d();
  bonsai->localTree.bodies_Pvel.h2d();
  bonsai->localTree.bodies_ids.h2d();
  bonsai->localTree.bodies_acc1.h2d();
  bonsai->localTree.bodies_time.h2d();


  //Build a tree-structure for initial initialization
  bonsai->sort_bodies(bonsai->localTree, true);
  bonsai->build(bonsai->localTree);
  if (maxlevels_exceeded) return -4;
  bonsai->allocateTreePropMemory(bonsai->localTree);
  bonsai->compute_properties(bonsai->localTree);


  curStateOnHost = false;

  return 0;
}


//Particle property functions
int get_mass(int id, double * mass){
  getCurrentStateToHost();
  int index_of_the_particle = getIdxFromId(id);  
  if(index_of_the_particle < 0)     return -3;
  
  *mass = bonsai->localTree.bodies_pos[index_of_the_particle].w;
  return 0;
}

int set_mass(int *index, double *mass, int length){
  getCurrentStateToHost();
  
  for (int i = 0; i < length; i++)
  { 
    int index_of_the_particle = getIdxFromId(index[i]);  
    if(index_of_the_particle < 0)    return -3;

    bonsai->localTree.bodies_pos[index_of_the_particle].w = mass[i];
  }
  
  bonsai->localTree.bodies_pos.h2d();
  
  return 0;
}

int get_radius(int id, double * radius){
  int index_of_the_particle = getIdxFromId(id);    
  if(index_of_the_particle < 0)     return -3;
  
  *radius = radii[index_of_the_particle];
  return 0;
}

int set_radius(int id, double radius){
  int index_of_the_particle = getIdxFromId(id);  
  if(index_of_the_particle < 0)     return -3;
  
  radii[index_of_the_particle] = radius;
  return 0;
}

int get_potential(int id, double * potential){
  getCurrentStateToHost();
  int index_of_the_particle = getIdxFromId(id);  
  if(index_of_the_particle < 0)     return -3;
  
  *potential = bonsai->localTree.bodies_acc1[index_of_the_particle].w;
  return 0;
}

int set_acceleration(int *id, double *ax, double *ay, double *az, int length){
  getCurrentStateToHost();
  
  for(int i=0; i < length; i++)
  {
    int index_of_the_particle = getIdxFromId(id[i]);  
    if(index_of_the_particle < 0)    return -3;
    
    bonsai->localTree.bodies_acc1[index_of_the_particle].x = ax[i];
    bonsai->localTree.bodies_acc1[index_of_the_particle].y = ay[i];
    bonsai->localTree.bodies_acc1[index_of_the_particle].z = az[i];
  }
  
  bonsai->localTree.bodies_acc1.h2d();
  
  return 0;
}

int get_acceleration(int id, double * ax, double * ay, double * az){
  getCurrentStateToHost();
  int index_of_the_particle = getIdxFromId(id);  
  if(index_of_the_particle < 0)     return -3;
  
  *ax = bonsai->localTree.bodies_acc1[index_of_the_particle].x;
  *ay = bonsai->localTree.bodies_acc1[index_of_the_particle].y;
  *az = bonsai->localTree.bodies_acc1[index_of_the_particle].z;
  return 0;
}

int get_velocity(int id, double * vx, double * vy, double * vz){
  getCurrentStateToHost();
  int index_of_the_particle = getIdxFromId(id);  
  if(index_of_the_particle < 0)     return -3;
  
  *vx = bonsai->localTree.bodies_vel[index_of_the_particle].x;
  *vy = bonsai->localTree.bodies_vel[index_of_the_particle].y;
  *vz = bonsai->localTree.bodies_vel[index_of_the_particle].z;
  return 0;
}

int get_position(int id, double * x, double * y, double * z){
  getCurrentStateToHost();
  int index_of_the_particle = getIdxFromId(id);  
  if(index_of_the_particle < 0)    return -3;
  
  *x = bonsai->localTree.bodies_pos[index_of_the_particle].x;
  *y = bonsai->localTree.bodies_pos[index_of_the_particle].y;
  *z = bonsai->localTree.bodies_pos[index_of_the_particle].z;
  return 0;
}

int set_position(int *id, double *x, double *y, double *z, int length){
  getCurrentStateToHost();
    
  for(int i=0; i < length; i++)
  {
    int index_of_the_particle = getIdxFromId(id[i]);  
    if(index_of_the_particle < 0)    return -3;
    
    bonsai->localTree.bodies_pos[index_of_the_particle].x = x[i];
    bonsai->localTree.bodies_pos[index_of_the_particle].y = y[i];
    bonsai->localTree.bodies_pos[index_of_the_particle].z = z[i];
  }
  
  bonsai->localTree.bodies_pos.h2d();
  
  return 0;
}

int set_velocity(int *id, double *vx, double *vy, double *vz, int length){
  getCurrentStateToHost();
  
  for(int i=0; i < length; i++)
  {
    int index_of_the_particle = getIdxFromId(id[i]);  
    if(index_of_the_particle < 0)    return -3;
  
    bonsai->localTree.bodies_vel[index_of_the_particle].x = vx[i];
    bonsai->localTree.bodies_vel[index_of_the_particle].y = vy[i];
    bonsai->localTree.bodies_vel[index_of_the_particle].z = vz[i];
  }

  bonsai->localTree.bodies_vel.h2d();
  return 0;
}


int set_state(int *index, double *mass, double *x, double *y, double *z, 
              double *vx, double *vy, double *vz, double * radius, int length){
  getCurrentStateToHost();
  
  for (int i = 0; i < length; i++)    
  {
    int index_of_the_particle = getIdxFromId(index[i]);  
    if(index_of_the_particle < 0)    return -3;
    
    bonsai->localTree.bodies_pos[index_of_the_particle].x = x[i];
    bonsai->localTree.bodies_pos[index_of_the_particle].y = y[i];
    bonsai->localTree.bodies_pos[index_of_the_particle].z = z[i];
    bonsai->localTree.bodies_pos[index_of_the_particle].w = mass[i];
    
    bonsai->localTree.bodies_vel[index_of_the_particle].x = vx[i];
    bonsai->localTree.bodies_vel[index_of_the_particle].y = vy[i];
    bonsai->localTree.bodies_vel[index_of_the_particle].z = vz[i];            
    
    radii[index_of_the_particle] = radius[i];
  }
  
  
  bonsai->localTree.bodies_pos.h2d();    
  bonsai->localTree.bodies_vel.h2d();
  
  return 0;
}

//End particle properties

//Simulation system properties, stats, etc

int get_total_radius(double * radius){
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return -2;
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
  return -2;
}
int get_index_of_next_particle(int index_of_the_particle,
                int * index_of_the_next_particle){
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return -2;
}

int get_potential_at_point(double eps, double x, double y, double z,
  double * phi){
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return -2;
}


int get_center_of_mass_position(double * x, double * y, double * z){
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return -2;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return -2;
}


int get_gravity_at_point(double eps, double x, double y, double z,
  double * forcex, double * forcey, double * forcez){
  fprintf(stderr,"NOT IMPLEMENTED: %s:%d \n", __FILE__, __LINE__);
  return -2;
}

int cleanup_code(){
  if (initialized) {
    getCurrentStateToHost();
  }
  
  bodies_pos.clear();
  bodies_vel.clear();
  bodies_grav.clear();
  starids.clear();
  radii.clear();
  bodies_time.clear();
  
  total_mass = 0;
  n_bodies = 0;
  id_counter = 0;
  idToIndex.clear();

  //TODO waarom de test op n? 
  if (initialized && bonsai->localTree.n > 0) delete bonsai;
  bonsai = NULL;
  initialized = false;
  return 0;
}


int set_begin_time(double input) {
    begin_time = input;
    return 0;
}

int get_begin_time(double * output) {
    *output = begin_time;
    return 0;
}
