#include "interface.h"
#include <stopcond.h>
#include <iostream>

// A stub of this file is machine generated, but the content is
// hand-coded.  SAVE A COPY (here interface.cc.1) to avoid accidental
// overwriting!

#include "src/localassert.h"
#include "src/hacs64.h"

std::vector<particle> node::ptcl;
std::vector<node>     node::node_heap;
std::vector<std::pair<node*, node*> > node::pair_list;
std::vector< std::pair<int, int> > UpdatedPtcl;

inline double SQR(const double x) {return x*x;}


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
static bool is_commited = false;

inline int get_id_from_idx(const int index_of_the_particle)
{
  if (nbody_ptr->index2id_map.find(index_of_the_particle) == nbody_ptr->index2id_map.end())  
  {
    fprintf(stderr, " attempted to get non-existend id= %d [ %u ]\n", 
		    index_of_the_particle, nbody_ptr->cyclical_idx);
    return -1;
  }
  else
    return nbody_ptr->index2id_map[index_of_the_particle];
}

int handle_assert_failed(assert_failed & ex)
{
    std::cerr << std::endl;
    std::cerr << ex;
    std::cerr << std::endl;
    return -10;
}

int initialize_code()
{
    try
    {
        assert(nbody_ptr == NULL); 
        fpe_catch();
        nbody_ptr = new hacs64::Nbody;

        // set default parameters to sane values
        is_commited = false;
        
        nbody_ptr->nmax = 10000; // will be enough for small calculations but will fail for bigger
        nbody_ptr->dtmax = 1;    // dtmax will influence timestep scaling, for now take 1 nbody time unit
    // AMUSE STOPPING CONDITIONS SUPPORT
        set_support_for_condition(COLLISION_DETECTION);
#if 0
        set_support_for_condition(PAIR_DETECTION);
        set_support_for_condition(TIMEOUT_DETECTION);
        set_support_for_condition(NUMBER_OF_STEPS_DETECTION);
        set_support_for_condition(OUT_OF_BOX_DETECTION);
#endif
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
  return 0;
}

int cleanup_code()
{
    try
    {
        if(nbody_ptr == NULL){
        return 0;
        }
        delete nbody_ptr;
        nbody_ptr = NULL;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
   }    
}

/******************/

int set_nmax(int nmax)
{
    try
    {
        assert(nbody_ptr != NULL);
        if (is_commited) return -1;

        nbody_ptr->nmax = nmax;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}
int get_nmax(int * nmax)
{
    try
    {
        assert(nbody_ptr != NULL);
        *nmax = nbody_ptr->nmax;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/******************/

int set_dtmax(double dtmax)
{
    try
    {
        assert(nbody_ptr != NULL);
        if (is_commited) return -1;
        nbody_ptr->dtmax = dtmax;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}
int get_dtmax(double * dtmax)
{
    try
    {
        assert(nbody_ptr != NULL);
        *dtmax = nbody_ptr->dtmax;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/******************/

int set_eps2(double epsilon_squared)
{
    try
    {
        assert(nbody_ptr != NULL);
        nbody_ptr->eps2 = epsilon_squared;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}
int get_eps2(double * epsilon_squared)
{
    try
    {
        assert(nbody_ptr != NULL);
        *epsilon_squared = nbody_ptr->eps2;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/*****************/

int set_h2max(double h2max)
{
    try
    {
        assert(nbody_ptr != NULL);
        nbody_ptr->h2max = h2max;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}
int get_h2max(double * h2max)
{
    try
    {
        assert(nbody_ptr != NULL);
        *h2max = nbody_ptr->h2max;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/*****************/

int set_eta_reg(double eta_reg)
{
    try
    {
        assert(nbody_ptr != NULL);
        nbody_ptr->eta_reg = eta_reg;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}
int get_eta_reg(double *eta_reg)
{
    try
    {
        assert(nbody_ptr != NULL);
        *eta_reg = nbody_ptr->eta_reg;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/******************/

int set_eta_irr(double eta_irr)
{
    try
    {
        assert(nbody_ptr != NULL);
        nbody_ptr->eta_irr = eta_irr;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}
int get_eta_irr(double *eta_irr)
{
    try
    {
        assert(nbody_ptr != NULL);
        *eta_irr = nbody_ptr->eta_irr;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/****************/

int commit_parameters()
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->irr_ptr == NULL);
        assert(nbody_ptr->reg_ptr == NULL);
        assert(!is_commited);
        nbody_ptr->commit_parameters();
        is_commited = true;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}
int recommit_parameters()
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        nbody_ptr->recommit_parameters();
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/****************/
/****************/
/****************/

int new_particle(
    int *index_of_the_particle,
    double mass, 
    double x, double y, double z,
    double vx, double vy, double vz,
    double radius)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        nbody_ptr->cyclical_idx &= 0x7FFFFFFF;
        
        *index_of_the_particle = nbody_ptr->cyclical_idx++;
#if 0
        fprintf(stderr , "--new-particle-added= %d %d \n",
		  *index_of_the_particle, index_to_set);
#endif
        nbody_ptr->ptcl2add.push_back(hacs64::Particle(mass, radius, dvec3(x,y,z), dvec3(vx,vy,vz), *index_of_the_particle));
        UpdatedPtcl.push_back(std::make_pair(*index_of_the_particle, 2));
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}
int delete_particle(int index_of_the_particle)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        const int id = get_id_from_idx(index_of_the_particle);
        if (id == -1) return -1;
        nbody_ptr->ptcl[id].id = -1;
        //  nbody_ptr->ptcl2remove.push_back(id);
        UpdatedPtcl.push_back(std::make_pair(index_of_the_particle, 1));
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/****************/

int set_state(int index_of_the_particle,
    double mass, 
    double x, double y, double z,
    double vx, double vy, double vz,
    double radius)
{
    try
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
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}
int get_state(int index_of_the_particle,
    double * mass,
    double * x, double * y, double * z,
    double * vx, double * vy, double * vz,
    double * radius)
{
    try
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
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/****************/

int set_mass(int index_of_the_particle, double mass)
{
    try
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
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}
int get_mass(int index_of_the_particle, double * mass)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        const int id = get_id_from_idx(index_of_the_particle);
        if (id == -1) return -1;
        *mass = nbody_ptr->ptcl[id].mass;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    } 
}

/****************/

int set_radius(int index_of_the_particle, double radius)
{
    try
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
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}
int get_radius(int index_of_the_particle, double * radius)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        const int id = get_id_from_idx(index_of_the_particle);
        if (id == -1) return -1;
        *radius = nbody_ptr->ptcl[id].radius;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/****************/

int set_position(int index_of_the_particle,
    double x, double y, double z)
{
    try
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
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}
int get_position(int index_of_the_particle,
    double * x, double * y, double * z)
{
    try
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
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/****************/

int set_velocity(int index_of_the_particle,
    double vx, double vy, double vz)
{
    try
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
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}
int get_velocity(int index_of_the_particle,
    double * vx, double * vy, double * vz)
{
    try
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
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/****************/

int set_acceleration(int index_of_the_particle,
    double ax, double ay, double az)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        return -1;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}
int get_acceleration(int index_of_the_particle,
    double * ax, double * ay, double * az)
{
    try
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
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/****************/

int get_potential(int index_of_the_particle, double * pot)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        const int id = get_id_from_idx(index_of_the_particle);
        if (id == -1) return -1;
        *pot = nbody_ptr->ptcl[id].pot;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/****************/

int commit_particles()
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        nbody_ptr->commit_particles();
        std::cerr<<"HALLO!"<<std::endl;
        nbody_ptr->get_epot();
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}
int recommit_particles()
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        nbody_ptr->recommit_particles();
        nbody_ptr->get_epot();
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/****************/

int synchronize_model()
{
    try
    {
        // Synchronize all particles at the current system time.  The
        // default is not to reinitialize the scheduler, as this will be
        // handled later, in recommit_particles().

        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        UpdatedPtcl.clear();
        nbody_ptr->__synchmodel();
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

int predict_model()
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        nbody_ptr->__predictmodel();
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/****************/

int get_time(double * sys_time)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        *sys_time = nbody_ptr->t_global;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

int get_time_step(double * time_step)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        *time_step = nbody_ptr->dt_global;	
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

int get_index_of_first_particle(int * index_of_the_particle)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        *index_of_the_particle = nbody_ptr->ptcl[0].id;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

int get_index_of_next_particle(int index, int *index_of_the_next_particle)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        const int id = get_id_from_idx(index);
        if (id == -1) return -1;
        if (id+1 < nbody_ptr->ptcl.size())
        {
            *index_of_the_next_particle = nbody_ptr->ptcl[id+1].id;
            return 0;
        }
        else
        {
            *index_of_the_next_particle = -1;
            return 1;
        }
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

#if 0
int get_indices_of_colliding_particles(int * index_of_particle1, 
    int * index_of_particle2)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        return -1;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}
#endif

int get_number_of_particles(int * number_of_particles)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        *number_of_particles = nbody_ptr->ptcl.size();
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

int get_total_mass(double * mass)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        *mass = nbody_ptr->get_total_mass();
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

int get_potential_energy(double * potential_energy)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        *potential_energy = nbody_ptr->get_epot();
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

int get_kinetic_energy(double * kinetic_energy)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        *kinetic_energy = nbody_ptr->get_ekin();
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

int get_center_of_mass_position(double * x, double * y, double * z)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        const dvec3 pos = nbody_ptr->get_com_pos();
        *x = pos.x;
        *y = pos.y;
        *z = pos.z;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        const dvec3 vel = nbody_ptr->get_com_vel();
        *vx = vel.x;
        *vy = vel.y;
        *vz = vel.z;
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

int get_total_radius(double * radius)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        *radius = nbody_ptr->get_total_radius();
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

int get_potential_at_point(double eps,
    double x, double y, double z, 
    double * phi)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        return -1;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

int get_gravity_at_point(double eps, double x, double y, double z, 
    double * forcex, double * forcey, double * forcez)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        return -1;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/****************/

int get_number_of_particles_updated(int * value)
{
    try
    {
        *value = UpdatedPtcl.size();
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

int get_id_of_updated_particle(int index, int * index_of_particle, int * status)
{
    try
    {
        if (index < 0 || index >= (int)UpdatedPtcl.size())
        return -1;

        *index_of_particle = UpdatedPtcl[index].first;  // id
        *status            = UpdatedPtcl[index].second; // status  1 - remove, 2- add

        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

/****************/
/****************/
/****************/

int evolve_model(double time)
{
    try
    {
        assert(nbody_ptr != NULL);
        assert(nbody_ptr->is_sane());
        UpdatedPtcl.clear();
        nbody_ptr->evolve_model(time);
        return 0;
    }
    catch(assert_failed ex)
    {
        return handle_assert_failed(ex);
    }    
}

namespace hacs64
{
  void Nbody::evolve_model(const double t_next)
  {
    assert(ptcl2add.empty());
    assert(ptcl2remove.empty());

    reset_stopping_conditions();    

    int is_collision_detection_enabled;
    int is_pair_detection_enabled;
    is_stopping_condition_enabled(COLLISION_DETECTION, 
        &is_collision_detection_enabled);
    is_stopping_condition_enabled(PAIR_DETECTION, 
        &is_pair_detection_enabled);

#if 0
    is_collision_detection_enabled = false;
#endif

    while (t_global < t_next)
    {
      iterate();
#if 0
      if (iteration%33 == 0) break;
#endif

      const int nirr = irr_list.size();
      int   i_close = -1;
      float r2close = HUGE;

      int   i_coll = -1;

      for (int ix = 0; ix < nirr; ix++)
      {
        const int i = irr_list[ix];
        if (ptcl[i].jr2 < r2close)
        {
          i_close = ptcl[i].jnb;
          r2close = ptcl[i].jr2;
        }
        if (ptcl[i].jr2 < SQR(ptcl[i].radius + ptcl[ptcl[i].jnb].radius))
        {
          i_coll = ptcl[i].jnb;
        }
      }
      assert(r2close > 0.0);

      const float r2crit = SQR(2.0/(ptcl.size()));

      if (is_collision_detection_enabled) 
      {
        if (ptcl[i_close].jr2 < r2crit)
        {
          int stopping_index = next_index_for_stopping_condition();
          set_stopping_condition_info(stopping_index, COLLISION_DETECTION);
          set_stopping_condition_particle_index(stopping_index, 0, ptcl[i_close].id);
          set_stopping_condition_particle_index(stopping_index, 1, ptcl[ptcl[i_close].jnb].id);
          break;
        }
      }


      if (is_pair_detection_enabled)
      {
        assert(false);
      }
    }

  }
}





