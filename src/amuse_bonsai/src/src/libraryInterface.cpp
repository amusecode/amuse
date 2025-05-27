#include "octree.h"

void octree::setEps(float eps)
{
    eps2  = eps*eps;  
}

float octree::getEps()
{
  return sqrt(eps2);
}

void octree::setDt(float dt)
{
  timeStep = dt;
}

float octree::getDt()
{
  return timeStep;
}


void octree::setEta(float _eta)
{
  eta = _eta;
}

float octree::getEta()
{
  return eta;
}


void octree::setTheta(float thetax)
{
    theta = thetax; 
    inv_theta   = 1.0/theta;    
}

float octree::getTheta()
{
  return theta;
}

void octree::setTEnd(float tEndx)
{
    tEnd  = tEndx;  
}

float octree::getTEnd()
{
  return tEnd;
}

void octree::setTime(float t_now)
{
   t_current = t_now;
}

float octree::getTime()
{
  return t_current;
}

float octree::getPot()
{
  return Epot;
}

float octree::getKin()
{
  return Ekin;
}

//   void calcGravityOnParticles(real4 *bodyPositions, real4 *bodyVelocities, int *bodyIDs);
// 


#if 0

#include <string>
#include <iostream>
#include <fstream>
#include "src/include/octree.h"



octree *bonsai;
bool initialized = false;
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

  std::string logFileName = "bonsaiLog.txt";
  std::ofstream logFile(logFileName.c_str());
 //Creat the octree class and set the properties
  octree *tree = new octree(devID, theta, eps, snapshotFile, snapshotIter, timeStep, tEnd);  
  tree->set_context(logFile, false); //Do logging to file and enable timing (false = enabled)
  tree->load_kernels();

  initialized = true;

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


int commit_particles()
{
  assert(initialized == true);
  //tree setup
  initialize_particles();

  bonsai->localTree.setN(n_bodies);
  bonsai->allocateParticleMemory(tree->localTree);

  //TODO remove this when its in the bonsai code itself
  //Should put in function allocSupportMem
  int nblock_reduce = NBLOCK_REDUCE;
  bonsai->tnext.ccalloc(nblock_reduce,false);

  //Load data onto the device
  for(uint i=0; i < bodyPositions.size(); i++)
  {
    tree->localTree.bodies_pos[i] = bodyPositions[i];
    tree->localTree.bodies_vel[i] = bodyVelocities[i];
    tree->localTree.bodies_ids[i] = bodyIDs[i];

    tree->localTree.bodies_Ppos[i] = bodyPositions[i];
    tree->localTree.bodies_Pvel[i] = bodyVelocities[i];
  }

   tree->localTree.bodies_pos.h2d();
   tree->localTree.bodies_vel.h2d();
   tree->localTree.bodies_Ppos.h2d();
   tree->localTree.bodies_Pvel.h2d();
   tree->localTree.bodies_ids.h2d();



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


int get_potential_energy(double *potential_energy)
{
//TODO  *potential_energy = calcEpot(bodies_pos,bodies_grav);
  return 0;
}

int get_kinetic_energy(double *kinetic_energy)
{
//TODO  *kinetic_energy = calcEkin(bodies_pos,bodies_vel);
  return 0;
}


int get_state(int id, double *mass, double *radius, double *x, double *y, double *z, double *vx, double *vy, double *vz)
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


int get_eps2(double *epsilon_squared)
{
  *epsilon_squared = eps*eps;
  return 0;
}


int set_eps2(double epsilon_squared)
{
  eps = sqrt(epsilon_squared);
  return 0;
}


int get_time(double *time)
{
  *time = t_now;
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





int echo(int input)
{
	initialize_code();
	return input;
}
#endif
