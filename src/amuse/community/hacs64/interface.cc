#include "interface.h"

// A stub of this file is machine generated, but the content is
// hand-coded.  SAVE A COPY (here interface.cc.1) to avoid accidental
// overwriting!

#include "src/hacs64.h"

std::vector<particle> node::ptcl;
std::vector<node>     node::node_heap;
std::vector<std::pair<node*, node*> > node::pair_list;


#ifndef __MACOSX_
#define __LINUX__
#endif

#ifdef __MACOSX__
#include <Accelerate/Accelerate.h>
#include <xmmintrin.h>
inline void fpe_catch() {
	_mm_setcsr( _MM_MASK_MASK &~
			(_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO) );
}
#elif defined __LINUX__
#include <fenv.h>
void fpe_catch(void) 
{
	/* Enable some exceptions. At startup all exceptions are masked. */
	feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
#else
crap
void fpe_catch(void) {}
#endif


static hacs64::Nbody *nbody_ptr = NULL;

inline int get_id_from_idx(const int index_of_the_particle)
{
  if (nbody_ptr->index2id_map.find(index_of_the_particle) == nbody_ptr->index2id_map.end())  
    return -1;
  else
    return nbody_ptr->index2id_map[index_of_the_particle];
}

int initialize_code()
{
  assert(nbody_ptr == NULL); 
  fpe_catch();
  nbody_ptr = new hacs64::Nbody;

  return 0;
}
int cleanup_code()
{
  assert(nbody_ptr != NULL);
  delete nbody_ptr;
  nbody_ptr = NULL;
  return 0;
}

/******************/

int set_nmax(int nmax)
{
  assert(nbody_ptr != NULL);
  if (nbody_ptr->nmax > 0) return -1;

  nbody_ptr->nmax = nmax;
  return 0;
}
int get_nmax(int * nmax)
{
  assert(nbody_ptr != NULL);
  *nmax = nbody_ptr->nmax;
  return 0;
}

/******************/

int set_dtmax(double dtmax)
{
  assert(nbody_ptr != NULL);
  if (nbody_ptr->dtmax > 0.0) return -1;
  nbody_ptr->dtmax = dtmax;
  return 0;
}
int get_dtmax(double * dtmax)
{
  assert(nbody_ptr != NULL);
  *dtmax = nbody_ptr->dtmax;
  return 0;
}

/******************/

int set_eps2(double epsilon_squared)
{
  assert(nbody_ptr != NULL);
  nbody_ptr->eps2 = epsilon_squared;
  return 0;
}
int get_eps2(double * epsilon_squared)
{
  assert(nbody_ptr != NULL);
  *epsilon_squared = nbody_ptr->eps2;
  return 0;
}

/*****************/

int set_h2max(double h2max)
{
  assert(nbody_ptr != NULL);
  nbody_ptr->h2max = h2max;
  return 0;
}
int get_h2max(double * h2max)
{
  assert(nbody_ptr != NULL);
  *h2max = nbody_ptr->h2max;
  return 0;
}

/*****************/

int set_eta_reg(double eta_reg)
{
  assert(nbody_ptr != NULL);
  nbody_ptr->eta_reg = eta_reg;
  return 0;
}
int get_eta_reg(double *eta_reg)
{
  assert(nbody_ptr != NULL);
  *eta_reg = nbody_ptr->eta_reg;
  return 0;
}

/******************/

int set_eta_irr(double eta_irr)
{
  assert(nbody_ptr != NULL);
  nbody_ptr->eta_irr = eta_irr;
  return 0;
}
int get_eta_irr(double *eta_irr)
{
  assert(nbody_ptr != NULL);
  *eta_irr = nbody_ptr->eta_irr;
  return 0;
}

/****************/

int commit_parameters()
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->irr_ptr == NULL);
  assert(nbody_ptr->reg_ptr == NULL);
  nbody_ptr->commit_parameters();

  return 0;
}
int recommit_parameters()
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  nbody_ptr->recommit_parameters();
  return 0;
}

/****************/
/****************/
/****************/

int new_particle(
    int *index_of_the_particle,
    double mass, double radius, 
    double x, double y, double z,
    double vx, double vy, double vz)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  *index_of_the_particle = nbody_ptr->cyclical_idx++;
  nbody_ptr->ptcl2add.push_back(hacs64::Particle(mass, radius, dvec3(x,y,z), dvec3(vx,vy,vz), *index_of_the_particle));
  return 0;
}
int delete_particle(int index_of_the_particle)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  const int id = get_id_from_idx(index_of_the_particle);
  if (id == -1) return -1;
  nbody_ptr->ptcl2remove.push_back(id);
  return 0;
}

/****************/

int set_state(int index_of_the_particle,
    double mass, double radius, 
    double x, double y, double z,
    double vx, double vy, double vz)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  const int id = get_id_from_idx(index_of_the_particle);
  if (id == -1) return -1;
  nbody_ptr->ptcl[id] = hacs64::Particle(mass, radius, dvec3(x,y,z), dvec3(vx,vy,vz), index_of_the_particle);
#if 0
  nbody_ptr->ptcl2modify.push_back(std::make_pair(index_of_the_particle, hacs64::Particle::ALL));
#endif
  return 0;
}
int get_state(int index_of_the_particle,
    double * mass, double * radius, 
    double * x, double * y, double * z,
    double * vx, double * vy, double * vz)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  const int id = get_id_from_idx(index_of_the_particle);
  if (id == -1) return -1;
  const hacs64::Particle &pi = nbody_ptr->ptcl[id];
  assert(pi.id == index_of_the_particle);
  *mass   = pi.mass;
  *radius = pi.radius;
  *x      = pi.pos.x;
  *y      = pi.pos.y;
  *z      = pi.pos.z;
  *vx     = pi.vel.x;
  *vy     = pi.vel.y;
  *vz     = pi.vel.z;
  return 0;
}

/****************/

int set_mass(int index_of_the_particle, double mass)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  const int id = get_id_from_idx(index_of_the_particle);
  if (id == -1) return -1;
  nbody_ptr->ptcl[id].mass = mass;
#if 0
  nbody_ptr->ptcl2modify.push_back(std::make_pair(index_of_the_particle, hacs64::Particle::MASS));
#endif
  return 0;
}
int get_mass(int index_of_the_particle, double * mass)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  const int id = get_id_from_idx(index_of_the_particle);
  if (id == -1) return -1;
  *mass = nbody_ptr->ptcl[id].mass;
  return 0;
}

/****************/

int set_radius(int index_of_the_particle, double radius)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  const int id = get_id_from_idx(index_of_the_particle);
  if (id == -1) return -1;
  nbody_ptr->ptcl[id].radius = radius;
#if 0
  nbody_ptr->ptcl2modify.push_back(std::make_pair(index_of_the_particle, hacs64::Particle::RADIUS));
#endif
  return 0;
}
int get_radius(int index_of_the_particle, double * radius)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  const int id = get_id_from_idx(index_of_the_particle);
  if (id == -1) return -1;
  *radius = nbody_ptr->ptcl[id].radius;
  return 0;
}

/****************/

int set_position(int index_of_the_particle,
    double x, double y, double z)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  const int id = get_id_from_idx(index_of_the_particle);
  if (id == -1) return -1;
  nbody_ptr->ptcl[id].pos = dvec3(x,y,z);
#if 0
  nbody_ptr->ptcl2modify.push_back(std::make_pair(index_of_the_particle, hacs64::Particle::POS));
#endif
  return 0;
}
int get_position(int index_of_the_particle,
    double * x, double * y, double * z)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  const int id = get_id_from_idx(index_of_the_particle);
  if (id == -1) return -1;
  *x = nbody_ptr->ptcl[id].pos.x;
  *y = nbody_ptr->ptcl[id].pos.y;
  *z = nbody_ptr->ptcl[id].pos.z;
  return 0;
}

/****************/

int set_velocity(int index_of_the_particle,
    double vx, double vy, double vz)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  const int id = get_id_from_idx(index_of_the_particle);
  if (id == -1) return -1;
  nbody_ptr->ptcl[id].vel = dvec3(vx, vy, vz);
#if 0
  nbody_ptr->ptcl2modify.push_back(std::make_pair(index_of_the_particle, hacs64::Particle::VEL));
#endif
  return 0;
}
int get_velocity(int index_of_the_particle,
    double * vx, double * vy, double * vz)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  const int id = get_id_from_idx(index_of_the_particle);
  if (id == -1) return -1;
  *vx = nbody_ptr->ptcl[id].vel.x;
  *vy = nbody_ptr->ptcl[id].vel.y;
  *vz = nbody_ptr->ptcl[id].vel.z;
  return 0;
}

/****************/

int set_acceleration(int index_of_the_particle,
    double ax, double ay, double az)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  return -1;
}
int get_acceleration(int index_of_the_particle,
    double * ax, double * ay, double * az)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  const int id = get_id_from_idx(index_of_the_particle);
  if (id == -1) return -1;
  *ax = nbody_ptr->ptcl[id].ftot.acc.x;
  *ay = nbody_ptr->ptcl[id].ftot.acc.y;
  *az = nbody_ptr->ptcl[id].ftot.acc.z;
  return 0;
}

/****************/

int get_potential(int index_of_the_particle, double * pot)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  const int id = get_id_from_idx(index_of_the_particle);
  if (id == -1) return -1;
  *pot = nbody_ptr->ptcl[id].pot;
  return 0;
}

/****************/

int commit_particles()
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  nbody_ptr->commit_particles();
  return 0;
}
int recommit_particles()
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  nbody_ptr->recommit_particles();
  return 0;
}

/****************/

int evolve_model(double time)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  nbody_ptr->evolve_model(time);

  return 0;	
}

/****************/

int synchronize_model()
{
  // Synchronize all particles at the current system time.  The
  // default is not to reinitialize the scheduler, as this will be
  // handled later, in recommit_particles().

  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  nbody_ptr->__synchmodel();
  return 0;
}

int predict_model()
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  nbody_ptr->__predictmodel();
  return 0;
}

/****************/

int get_time(double * sys_time)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  *sys_time = nbody_ptr->t_global;
  return 0;
}

int get_time_step(double * time_step)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  *time_step = nbody_ptr->dt_global;	
  return 0;
}

int get_index_of_first_particle(int * index_of_the_particle)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  *index_of_the_particle = nbody_ptr->ptcl[0].id;
  return 0;
}
int get_index_of_next_particle(int id, int *index_of_the_next_particle)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  return -1;
}

#if 0
int get_indices_of_colliding_particles(int * index_of_particle1, 
    int * index_of_particle2)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  return -1;
}
#endif

int get_number_of_particles(int * number_of_particles)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  *number_of_particles = nbody_ptr->ptcl.size();
  return 0;
}

int get_total_mass(double * mass)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  *mass = nbody_ptr->get_total_mass();
  return 0;
}

int get_potential_energy(double * potential_energy)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  *potential_energy = nbody_ptr->get_epot();
  return 0;
}

int get_kinetic_energy(double * kinetic_energy)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  *kinetic_energy = nbody_ptr->get_ekin();
  return 0;
}

int get_center_of_mass_position(double * x, double * y, double * z)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  const dvec3 pos = nbody_ptr->get_com_pos();
  *x = pos.x;
  *y = pos.y;
  *z = pos.z;
  return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  const dvec3 vel = nbody_ptr->get_com_vel();
  *vx = vel.x;
  *vy = vel.y;
  *vz = vel.z;
  return 0;
}

int get_total_radius(double * radius)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  *radius = nbody_ptr->get_total_radius();
  return 0;
}

int get_potential_at_point(double eps,
    double x, double y, double z, 
    double * phi)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  return -1;
}

int get_gravity_at_point(double eps, double x, double y, double z, 
    double * forcex, double * forcey, double * forcez)
{
  assert(nbody_ptr != NULL);
  assert(nbody_ptr->is_sane());
  return -1;
}

