#include "src/types.h"
#include "interface.h"
#include "src/evolve.h"


int highest_particle_index = 0;
int highest_external_particle_index = 0;
ParticlesMap particlesMap;
External_ParticlesMap external_particlesMap;

double relative_tolerance = 1.0e-16;
double absolute_tolerance_eccentricity_vectors = 1.0e-14;
bool include_quadrupole_order_terms = true;
bool include_octupole_order_binary_pair_terms = true;
bool include_octupole_order_binary_triplet_terms = false;
bool include_hexadecupole_order_binary_pair_terms = false;
bool include_dotriacontupole_order_binary_pair_terms = false;
int orbital_phases_random_seed = 0;

/*******************
/* basic interface *
 ******************/
 
int new_particle(int * index_of_the_particle, bool is_binary)
{

    *index_of_the_particle = highest_particle_index;
    Particle * p = new Particle(highest_particle_index, is_binary);
    particlesMap[highest_particle_index] = p;

    highest_particle_index++;

    return 0;
}
int delete_particle(int index_of_the_particle)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    particlesMap.erase(index_of_the_particle);
  
    return 0;
}

int set_children(int index_of_the_particle, int child1, int child2)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    p->child1 = child1;
    p->child2 = child2;
    
    return 0;
}
int get_children(int index_of_the_particle, int *child1, int *child2)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *child1 = p->child1;
    *child2 = p->child2;
    
    return 0;
}

int set_mass(int index_of_the_particle, double value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    p->mass = value;
    
    return 0;
}
int get_mass(int index_of_the_particle, double *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->mass;
    
    return 0;
}

int set_mass_dot_external(int index_of_the_particle, double value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    p->mass_dot_external = value;
    
    return 0;
}
int get_mass_dot_external(int index_of_the_particle, double *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->mass_dot_external;
    
    return 0;
}

int set_radius(int index_of_the_particle, double value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    p->radius = value;
    
    return 0;
}
int get_radius(int index_of_the_particle, double *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->radius;
    
    return 0;
}

int set_radius_dot_external(int index_of_the_particle, double value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    p->radius_dot_external = value;
    
    return 0;
}
int get_radius_dot_external(int index_of_the_particle, double *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->radius_dot_external;
    
    return 0;
}

int set_radius_ddot_external(int index_of_the_particle, double value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    p->radius_ddot_external = value;
    
    return 0;
}
int get_radius_ddot_external(int index_of_the_particle, double *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->radius_ddot_external;
    
    return 0;
}

int get_level(int index_of_the_particle, int *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->level;
    
    return 0;
}

int set_stellar_type(int index_of_the_particle, int value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    p->stellar_type = value;
    
    return 0;
}
int get_stellar_type(int index_of_the_particle, int *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->stellar_type;
    
    return 0;
}

int set_true_anomaly(int index_of_the_particle, double value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    p->true_anomaly = value;
    
    return 0;
}
int get_true_anomaly(int index_of_the_particle, double *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->true_anomaly;
    
    return 0;
}
int set_sample_orbital_phases_randomly(int index_of_the_particle, bool value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    p->sample_orbital_phases_randomly = value;
    
    return 0;
}
int get_sample_orbital_phases_randomly(int index_of_the_particle, bool *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->sample_orbital_phases_randomly;
    
    return 0;
}


/*******************************
 * instantaneous perturbations *
 * ****************************/

int set_instantaneous_perturbation_delta_mass(int index_of_the_particle, double value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
    
    Particle * p = particlesMap[index_of_the_particle];
    p->instantaneous_perturbation_delta_mass = value;
    
    return 0;
}
int get_instantaneous_perturbation_delta_mass(int index_of_the_particle, double *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
    
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->instantaneous_perturbation_delta_mass;
    
    return 0;
}

int set_instantaneous_perturbation_delta_position(int index_of_the_particle, double x, double y, double z)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
    
    Particle * p = particlesMap[index_of_the_particle];
    p->instantaneous_perturbation_delta_position_x = x;
    p->instantaneous_perturbation_delta_position_y = y;
    p->instantaneous_perturbation_delta_position_z = z;
    
    return 0;
}
int get_instantaneous_perturbation_delta_position(int index_of_the_particle, double *x, double *y, double *z)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
    
    Particle * p = particlesMap[index_of_the_particle];
    *x = p->instantaneous_perturbation_delta_position_x;
    *y = p->instantaneous_perturbation_delta_position_y;
    *z = p->instantaneous_perturbation_delta_position_z;
    
    return 0;
}

int set_instantaneous_perturbation_delta_velocity(int index_of_the_particle, double x, double y, double z)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
    
    Particle * p = particlesMap[index_of_the_particle];
    p->instantaneous_perturbation_delta_velocity_x = x;
    p->instantaneous_perturbation_delta_velocity_y = y;
    p->instantaneous_perturbation_delta_velocity_z = z;
    
    return 0;
}
int get_instantaneous_perturbation_delta_velocity(int index_of_the_particle, double *x, double *y, double *z)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
    
    Particle * p = particlesMap[index_of_the_particle];
    *x = p->instantaneous_perturbation_delta_velocity_x;
    *y = p->instantaneous_perturbation_delta_velocity_y;
    *z = p->instantaneous_perturbation_delta_velocity_z;
    
    return 0;
}


/************
 * external *
 * *********/
 
int new_external_particle(int *index_of_the_external_particle, double mass)
{
    
    *index_of_the_external_particle = highest_external_particle_index;
    External_Particle *f = new External_Particle(highest_external_particle_index);
    external_particlesMap[highest_external_particle_index] = f;

    highest_external_particle_index++;
    f->mass = mass;

    return 0;
}
int delete_external_particle(int index_of_the_external_particle)
{
    if (index_of_the_external_particle > highest_external_particle_index)
    {
        return -1;
    }
  
    external_particlesMap.erase(index_of_the_external_particle);
  
    return 0;
}

int set_external_mass(int index_of_the_external_particle, double value)
{
    //printf("set_external_mass\n");
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }

    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    f->mass = value;
    
    return 0;
}
int get_external_mass(int index_of_the_external_particle, double *value)
{
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }
  
    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    *value = f->mass;
    
    return 0;
}

int set_external_path(int index_of_the_external_particle, int value)
{
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }

    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    f->path = value;
    
    return 0;
}
int get_external_path(int index_of_the_external_particle, int *value)
{
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }
  
    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    *value = f->path;
    
    return 0;
}

int set_external_mode(int index_of_the_external_particle, int value)
{
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }

    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    f->mode = value;
    
    return 0;
}
int get_external_mode(int index_of_the_external_particle, int *value)
{
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }
  
    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    *value = f->mode;
    
    return 0;
}

int set_external_t_ref(int index_of_the_external_particle, double value)
{
    //printf("set_external_t_ref\n");
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }

    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    f->t_ref = value;
    
    return 0;
}
int get_external_t_ref(int index_of_the_external_particle, double *value)
{
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }
  
    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    *value = f->t_ref;
    
    return 0;
}

int set_external_t_passed(int index_of_the_external_particle, double value)
{
    //printf("set_external_t_ref\n");
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }

    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    f->t_passed = value;
    
    return 0;
}
int get_external_t_passed(int index_of_the_external_particle, double *value)
{
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }
  
    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    *value = f->t_passed;
    
    return 0;
}

int set_external_r0_vectors(int index_of_the_external_particle, double vec_x, double vec_y, double vec_z)
{
    //printf("set_external_r0_vectors\n");
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }
    
    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    f->r0_vec_x = vec_x;
    f->r0_vec_y = vec_y;
    f->r0_vec_z = vec_z;
    
    return 0;
}
int get_external_r0_vectors(int index_of_the_external_particle, double *vec_x, double *vec_y, double *vec_z)
{
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }
  
    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    *vec_x = f->r0_vec_x;
    *vec_y = f->r0_vec_y;
    *vec_z = f->r0_vec_z;
    
    return 0;
}

int set_external_rdot_vectors(int index_of_the_external_particle, double rdot_vec_x, double rdot_vec_y, double rdot_vec_z)
{
    //printf("set_external_rdot_vectors\n");
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }
  
    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    f->rdot_vec_x = rdot_vec_x;
    f->rdot_vec_y = rdot_vec_y;
    f->rdot_vec_z = rdot_vec_z;
    
    return 0;
}
int get_external_rdot_vectors(int index_of_the_external_particle, double *rdot_vec_x, double *rdot_vec_y, double *rdot_vec_z)
{
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }
  
    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    *rdot_vec_x = f->rdot_vec_x;
    *rdot_vec_y = f->rdot_vec_y;
    *rdot_vec_z = f->rdot_vec_z;
    
    return 0;
}


int set_external_periapse_distance(int index_of_the_external_particle, double value)
{
    //printf("set_external_t_ref\n");
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }

    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    f->periapse_distance = value;
    
    return 0;
}
int get_external_periapse_distance(int index_of_the_external_particle, double *value)
{
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }
  
    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    *value = f->periapse_distance;
    
    return 0;
}

int set_external_eccentricity(int index_of_the_external_particle, double value)
{
    //printf("set_external_t_ref\n");
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }

    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    f->eccentricity = value;
    
    return 0;
}
int get_external_eccentricity(int index_of_the_external_particle, double *value)
{
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }
  
    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    *value = f->eccentricity;
    
    return 0;
}

int set_external_e_hat_vectors(int index_of_the_external_particle, double x, double y, double z)
{
    //printf("set_external_rdot_vectors\n");
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }
  
    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    f->e_hat_vec_x = x;
    f->e_hat_vec_y = y;
    f->e_hat_vec_z = z;
    
    return 0;
}
int get_external_e_hat_vectors(int index_of_the_external_particle, double *x, double *y, double *z)
{
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }
  
    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    *x = f->e_hat_vec_x;
    *y = f->e_hat_vec_y;
    *z = f->e_hat_vec_z;
    
    return 0;
}

int set_external_h_hat_vectors(int index_of_the_external_particle, double x, double y, double z)
{
    //printf("set_external_rdot_vectors\n");
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }
  
    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    f->h_hat_vec_x = x;
    f->h_hat_vec_y = y;
    f->h_hat_vec_z = z;
    
    return 0;
}
int get_external_h_hat_vectors(int index_of_the_external_particle, double *x, double *y, double *z)
{
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }
  
    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    *x = f->h_hat_vec_x;
    *y = f->h_hat_vec_y;
    *z = f->h_hat_vec_z;
    
    return 0;
}


int get_external_r_vectors(int index_of_the_external_particle, double *r_vec_x, double *r_vec_y, double *r_vec_z)
{
    if (index_of_the_external_particle > highest_external_particle_index)
    {
      return -1;
    }
  
    External_Particle *f = external_particlesMap[index_of_the_external_particle];
    *r_vec_x = f->r_vec_x;
    *r_vec_y = f->r_vec_y;
    *r_vec_z = f->r_vec_z;
    
    return 0;
}







/****************
/* spin vectors *
 ****************/

int set_spin_vector(int index_of_the_particle, double spin_vec_x, double spin_vec_y, double spin_vec_z)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->spin_vec_x = spin_vec_x;
    p->spin_vec_y = spin_vec_y;
    p->spin_vec_z = spin_vec_z;
    
    return 0;
}
int get_spin_vector(int index_of_the_particle, double *spin_vec_x, double *spin_vec_y, double *spin_vec_z)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *spin_vec_x = p->spin_vec_x;
    *spin_vec_y = p->spin_vec_y;
    *spin_vec_z = p->spin_vec_z;
    
    return 0;
}

int set_spin_vector_dot_external(int index_of_the_particle, double spin_vec_x_dot_external, double spin_vec_y_dot_external, double spin_vec_z_dot_external)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->spin_vec_x_dot_external = spin_vec_x_dot_external;
    p->spin_vec_y_dot_external = spin_vec_y_dot_external;
    p->spin_vec_z_dot_external = spin_vec_z_dot_external;
    
    return 0;
}
int get_spin_vector_dot_external(int index_of_the_particle, double *spin_vec_x_dot_external, double *spin_vec_y_dot_external, double *spin_vec_z_dot_external)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *spin_vec_x_dot_external = p->spin_vec_x_dot_external;
    *spin_vec_y_dot_external = p->spin_vec_y_dot_external;
    *spin_vec_z_dot_external = p->spin_vec_z_dot_external;
    
    return 0;
}

/****************************
/* orbital vectors/elements *
 ****************************/

int set_orbital_vectors(int index_of_the_particle, double e_vec_x, double e_vec_y, double e_vec_z, \
    double h_vec_x, double h_vec_y, double h_vec_z)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->e_vec_x = e_vec_x;
    p->e_vec_y = e_vec_y;
    p->e_vec_z = e_vec_z;
    p->h_vec_x = h_vec_x;
    p->h_vec_y = h_vec_y;
    p->h_vec_z = h_vec_z;
    
    return 0;
}
int get_orbital_vectors(int index_of_the_particle, double *e_vec_x, double *e_vec_y, double *e_vec_z, \
    double *h_vec_x, double *h_vec_y, double *h_vec_z)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *e_vec_x = p->e_vec_x;
    *e_vec_y = p->e_vec_y;
    *e_vec_z = p->e_vec_z;
    *h_vec_x = p->h_vec_x;
    *h_vec_y = p->h_vec_y;
    *h_vec_z = p->h_vec_z;

    return 0;
}

int set_orbital_vectors_dot_external(int index_of_the_particle, double e_vec_x_dot_external, double e_vec_y_dot_external, double e_vec_z_dot_external, \
    double h_vec_x_dot_external, double h_vec_y_dot_external, double h_vec_z_dot_external)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->e_vec_x_dot_external = e_vec_x_dot_external;
    p->e_vec_y_dot_external = e_vec_y_dot_external;
    p->e_vec_z_dot_external = e_vec_z_dot_external;
    p->h_vec_x_dot_external = h_vec_x_dot_external;
    p->h_vec_y_dot_external = h_vec_y_dot_external;
    p->h_vec_z_dot_external = h_vec_z_dot_external;
    
    return 0;
}
int get_orbital_vectors_dot_external(int index_of_the_particle, double *e_vec_x_dot_external, double *e_vec_y_dot_external, double *e_vec_z_dot_external, \
    double *h_vec_x_dot_external, double *h_vec_y_dot_external, double *h_vec_z_dot_external)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *e_vec_x_dot_external = p->e_vec_x_dot_external;
    *e_vec_y_dot_external = p->e_vec_y_dot_external;
    *e_vec_z_dot_external = p->e_vec_z_dot_external;
    *h_vec_x_dot_external = p->h_vec_x_dot_external;
    *h_vec_y_dot_external = p->h_vec_y_dot_external;
    *h_vec_z_dot_external = p->h_vec_z_dot_external;
    
    return 0;
}

int set_orbital_elements(int index_of_the_particle, double semimajor_axis, double eccentricity, \
    double inclination, double argument_of_pericenter, double longitude_of_ascending_node)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];    
    
    if (p->is_binary == false)
    {
        return 0;
    }

    /* determine masses in all binaries */
    int N_bodies, N_binaries, N_root_finding;
    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding);
    set_binary_masses_from_body_masses(&particlesMap);
    
    compute_orbital_vectors_from_orbital_elements(p->child1_mass, p->child2_mass, semimajor_axis, eccentricity, \
        inclination, argument_of_pericenter, longitude_of_ascending_node, \
        &(p->e_vec_x), &(p->e_vec_y), &(p->e_vec_z), &(p->h_vec_x), &(p->h_vec_y), &(p->h_vec_z) );
    
    return 0;
}
int get_orbital_elements(int index_of_the_particle, double *semimajor_axis, double *eccentricity, \
    double *inclination, double *argument_of_pericenter, double *longitude_of_ascending_node)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];    
    
    if (p->is_binary == false)
    {
        return 0;
    }

    double h_tot_vec[3];
    compute_h_tot_vector(&particlesMap,h_tot_vec);

    /* determine masses in all binaries */
    int N_bodies, N_binaries, N_root_finding;
    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding);
    set_binary_masses_from_body_masses(&particlesMap);
    
    compute_orbital_elements_from_orbital_vectors(p->child1_mass, p->child2_mass, h_tot_vec, \
        p->e_vec_x,p->e_vec_y,p->e_vec_z,p->h_vec_x,p->h_vec_y,p->h_vec_z,
        semimajor_axis, eccentricity, inclination, argument_of_pericenter, longitude_of_ascending_node);
    return 0;
}



int set_position_vector(int index_of_the_particle, double x, double y, double z)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->position_x = x;
    p->position_y = y;
    p->position_z = z;
   
    return 0;
}
int get_position_vector(int index_of_the_particle, double *x, double *y, double *z)
{
    //printf("get_position_vector\n");
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    set_positions_and_velocities(&particlesMap);
    
    Particle * p = particlesMap[index_of_the_particle];
    *x = p->position_x;
    *y = p->position_y;
    *z = p->position_z;
    
    return 0;
}

int set_velocity_vector(int index_of_the_particle, double x, double y, double z)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->velocity_x = x;
    p->velocity_y = y;
    p->velocity_z = z;
   
    return 0;
}
int get_velocity_vector(int index_of_the_particle, double *x, double *y, double *z)
{
    if (index_of_the_particle > highest_particle_index)
    {
      return -1;
    }

    set_positions_and_velocities(&particlesMap);
    
    Particle * p = particlesMap[index_of_the_particle];
    *x = p->velocity_x;
    *y = p->velocity_y;
    *z = p->velocity_z;
    
    return 0;
}

/************
/* PN terms *
 ************/

int set_include_pairwise_1PN_terms(int index_of_the_particle, bool include_pairwise_1PN_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    
    p->include_pairwise_1PN_terms = include_pairwise_1PN_terms;
        
    return 0;
}
int get_include_pairwise_1PN_terms(int index_of_the_particle, bool *include_pairwise_1PN_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    
    *include_pairwise_1PN_terms = p->include_pairwise_1PN_terms;
        
    return 0;
}

int set_include_pairwise_25PN_terms(int index_of_the_particle, bool include_pairwise_25PN_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    
    p->include_pairwise_25PN_terms = include_pairwise_25PN_terms;
        
    return 0;
}
int get_include_pairwise_25PN_terms(int index_of_the_particle, bool *include_pairwise_25PN_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    
    *include_pairwise_25PN_terms = p->include_pairwise_25PN_terms;
        
    return 0;
}


/*********
/* tides *
 *********/
int set_include_tidal_friction_terms(int index_of_the_particle, bool include_tidal_friction_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    
    p->include_tidal_friction_terms = include_tidal_friction_terms;
        
    return 0;
}
int get_include_tidal_friction_terms(int index_of_the_particle, bool *include_tidal_friction_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    
    *include_tidal_friction_terms = p->include_tidal_friction_terms;
        
    return 0;
}

int set_tides_method(int index_of_the_particle, int tides_method)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    
    p->tides_method = tides_method;
        
    return 0;
}
int get_tides_method(int index_of_the_particle, int *tides_method)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    
    *tides_method = p->tides_method;
        
    return 0;
}

int set_include_tidal_bulges_precession_terms(int index_of_the_particle, bool include_tidal_bulges_precession_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    
    p->include_tidal_bulges_precession_terms = include_tidal_bulges_precession_terms;
        
    return 0;
}
int get_include_tidal_bulges_precession_terms(int index_of_the_particle, bool *include_tidal_bulges_precession_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    
    *include_tidal_bulges_precession_terms = p->include_tidal_bulges_precession_terms;
        
    return 0;
}

int set_include_rotation_precession_terms(int index_of_the_particle, bool include_rotation_precession_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    
    p->include_rotation_precession_terms = include_rotation_precession_terms;
        
    return 0;
}

int get_include_rotation_precession_terms(int index_of_the_particle, bool *include_rotation_precession_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    
    *include_rotation_precession_terms = p->include_rotation_precession_terms;
        
    return 0;
}

int set_minimum_eccentricity_for_tidal_precession(int index_of_the_particle, double minimum_eccentricity_for_tidal_precession)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    
    p->minimum_eccentricity_for_tidal_precession = minimum_eccentricity_for_tidal_precession;
        
    return 0;
}
int get_minimum_eccentricity_for_tidal_precession(int index_of_the_particle, double *minimum_eccentricity_for_tidal_precession)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    
    *minimum_eccentricity_for_tidal_precession = p->minimum_eccentricity_for_tidal_precession;
        
    return 0;
}

int set_tides_apsidal_motion_constant(int index_of_the_particle, double value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];

    p->tides_apsidal_motion_constant = value;
        
    return 0;
}
int get_tides_apsidal_motion_constant(int index_of_the_particle, double *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    
    *value = p->tides_apsidal_motion_constant;
    
    return 0;
}

int set_tides_gyration_radius(int index_of_the_particle, double value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];

    p->tides_gyration_radius = value;
        
    return 0;
}
int get_tides_gyration_radius(int index_of_the_particle, double *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    
    *value = p->tides_gyration_radius;
    
    return 0;
}

int set_tides_viscous_time_scale(int index_of_the_particle, double value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];

    p->tides_viscous_time_scale = value;
        
    return 0;
}
int get_tides_viscous_time_scale(int index_of_the_particle, double *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    
    *value = p->tides_viscous_time_scale;
    
    return 0;
}

int set_tides_viscous_time_scale_prescription(int index_of_the_particle, int value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];

    p->tides_viscous_time_scale_prescription = value;
        
    return 0;
}
int get_tides_viscous_time_scale_prescription(int index_of_the_particle, int *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    
    *value = p->tides_viscous_time_scale_prescription;
    
    return 0;
}

int set_convective_envelope_mass(int index_of_the_particle, double value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];

    p->convective_envelope_mass = value;
        
    return 0;
}
int get_convective_envelope_mass(int index_of_the_particle, double *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    
    *value = p->convective_envelope_mass;
    
    return 0;
}

int set_convective_envelope_radius(int index_of_the_particle, double value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];

    p->convective_envelope_radius = value;
        
    return 0;
}
int get_convective_envelope_radius(int index_of_the_particle, double *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    
    *value = p->convective_envelope_radius;
    
    return 0;
}

int set_luminosity(int index_of_the_particle, double value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];

    p->luminosity = value;
        
    return 0;
}
int get_luminosity(int index_of_the_particle, double *value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    
    *value = p->luminosity;
    
    return 0;
}

/****************
/* root finding *
 ****************/

/* secular breakdown*/
int set_check_for_secular_breakdown(int index_of_the_particle, bool value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->check_for_secular_breakdown = value;

    return 0;
}
int get_check_for_secular_breakdown(int index_of_the_particle, bool* value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->check_for_secular_breakdown;

    return 0;
}

/* dynamical instablity*/
int set_check_for_dynamical_instability(int index_of_the_particle, bool value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->check_for_dynamical_instability = value;

    return 0;
}
int get_check_for_dynamical_instability(int index_of_the_particle, bool* value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->check_for_dynamical_instability;

    return 0;
}

int set_dynamical_instability_criterion(int index_of_the_particle, int value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->dynamical_instability_criterion = value;

    return 0;
}
int get_dynamical_instability_criterion(int index_of_the_particle, int* value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->dynamical_instability_criterion;

    return 0;
}


int set_dynamical_instability_central_particle(int index_of_the_particle, int value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->dynamical_instability_central_particle = value;

    return 0;
}
int get_dynamical_instability_central_particle(int index_of_the_particle, int* value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->dynamical_instability_central_particle;

    return 0;
}

int set_dynamical_instability_K_parameter(int index_of_the_particle, double value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->dynamical_instability_K_parameter = value;

    return 0;
}
int get_dynamical_instability_K_parameter(int index_of_the_particle, double* value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->dynamical_instability_K_parameter;

    return 0;
}

/* physical collision / orbit crossing*/
int set_check_for_physical_collision_or_orbit_crossing(int index_of_the_particle, bool value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->check_for_physical_collision_or_orbit_crossing = value;

    return 0;
}
int get_check_for_physical_collision_or_orbit_crossing(int index_of_the_particle, bool* value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->check_for_physical_collision_or_orbit_crossing;

    return 0;
}

/* minimum periapse distance reached */
int set_check_for_minimum_periapse_distance(int index_of_the_particle, bool value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->check_for_minimum_periapse_distance = value;

    return 0;
}
int get_check_for_minimum_periapse_distance(int index_of_the_particle, bool* value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->check_for_minimum_periapse_distance;

    return 0;
}
int set_check_for_minimum_periapse_distance_value(int index_of_the_particle, double value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->check_for_minimum_periapse_distance_value = value;

    return 0;
}
int get_check_for_minimum_periapse_distance_value(int index_of_the_particle, double* value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->check_for_minimum_periapse_distance_value;

    return 0;
}

/* RLOF at pericentre */
int set_check_for_RLOF_at_pericentre(int index_of_the_particle, bool value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->check_for_RLOF_at_pericentre = value;

    return 0;
}
int get_check_for_RLOF_at_pericentre(int index_of_the_particle, bool* value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->check_for_RLOF_at_pericentre;

    return 0;
}

int set_check_for_RLOF_at_pericentre_use_sepinsky_fit(int index_of_the_particle, bool value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->check_for_RLOF_at_pericentre_use_sepinsky_fit = value;

    return 0;
}
int get_check_for_RLOF_at_pericentre_use_sepinsky_fit(int index_of_the_particle, bool* value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->check_for_RLOF_at_pericentre_use_sepinsky_fit;

    return 0;
}


/* retrieve root finding state */
int set_root_finding_state(int index_of_the_particle, bool secular_breakdown_has_occurred, bool dynamical_instability_has_occurred, bool physical_collision_or_orbit_crossing_has_occurred, bool minimum_periapse_distance_has_occurred, bool RLOF_at_pericentre_has_occurred)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];

    p->secular_breakdown_has_occurred = secular_breakdown_has_occurred;
    p->dynamical_instability_has_occurred = dynamical_instability_has_occurred;
    p->physical_collision_or_orbit_crossing_has_occurred = physical_collision_or_orbit_crossing_has_occurred;
    p->minimum_periapse_distance_has_occurred = minimum_periapse_distance_has_occurred;
    p->RLOF_at_pericentre_has_occurred = RLOF_at_pericentre_has_occurred;
    
    return 0;
}
int get_root_finding_state(int index_of_the_particle, bool *secular_breakdown_has_occurred, bool *dynamical_instability_has_occurred, bool *physical_collision_or_orbit_crossing_has_occurred, bool* minimum_periapse_distance_has_occurred, bool *RLOF_at_pericentre_has_occurred)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];

    *secular_breakdown_has_occurred = p->secular_breakdown_has_occurred;
    *dynamical_instability_has_occurred = p->dynamical_instability_has_occurred;
    *physical_collision_or_orbit_crossing_has_occurred = p->physical_collision_or_orbit_crossing_has_occurred;
    *minimum_periapse_distance_has_occurred = p->minimum_periapse_distance_has_occurred;
    *RLOF_at_pericentre_has_occurred = p->RLOF_at_pericentre_has_occurred;
    
    return 0;
}


/********************
/* evolve interface *
 ********************/

int evolve_interface(double start_time, double time_step, double *output_time, double *hamiltonian, int *flag, int *error_code)
{
    int result = evolve(&particlesMap, &external_particlesMap, start_time, time_step, output_time, hamiltonian, flag, error_code);
    
    return result;
}

/* set levels and masses */
int determine_binary_parents_levels_and_masses_interface()
{
    //printf("determine_binary_parents_levels_and_masses_interface\n");
    int N_bodies, N_binaries, N_root_finding;
    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding);
    set_binary_masses_from_body_masses(&particlesMap);
    
    return 0;
}

int apply_external_perturbation_assuming_integrated_orbits_interface()
{
    //printf("apply_external_perturbation_assuming_integrated_orbits_interface\n");
    apply_external_perturbation_assuming_integrated_orbits(&particlesMap, &external_particlesMap);

    return 0;
}

int apply_user_specified_instantaneous_perturbation_interface()
{
    //printf("apply_user_specified_instantaneous_perturbation\n");
    apply_user_specified_instantaneous_perturbation(&particlesMap);
    
    return 0;
}

int set_positions_and_velocities_interface()
{
    set_positions_and_velocities(&particlesMap);
}

/**********************************************
/* orbital element/vector conversion routines *
 **********************************************/
void compute_h_tot_vector(ParticlesMap* particlesMap, double h_tot_vec[3])
{
    for (int i=0; i<3; i++)
    {
        h_tot_vec[i] = 0.0;
    }
    
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == 1)
        {
            h_tot_vec[0] += p->h_vec_x;
            h_tot_vec[1] += p->h_vec_y;
            h_tot_vec[2] += p->h_vec_z;
        }
    }
//    printf("compute_h_tot_vector %g %g %g\n",h_tot_vec[0],h_tot_vec[1],h_tot_vec[2]);
}
    
int compute_orbital_vectors_from_orbital_elements(double child1_mass, double child2_mass, double semimajor_axis, double eccentricity, double inclination, double argument_of_pericenter,double longitude_of_ascending_node, double *e_vec_x, double *e_vec_y, double *e_vec_z, double *h_vec_x, double *h_vec_y, double *h_vec_z)
{
    double cos_INCL = cos(inclination);
    double sin_INCL = sin(inclination);
    double cos_AP = cos(argument_of_pericenter);
    double sin_AP = sin(argument_of_pericenter);
    double cos_LAN = cos(longitude_of_ascending_node);
    double sin_LAN = sin(longitude_of_ascending_node);
           
    double h = (child1_mass*child2_mass*sqrt(CONST_G*semimajor_axis/(child1_mass+child2_mass)))*sqrt(1.0 - eccentricity*eccentricity);

    *e_vec_x = eccentricity*(cos_LAN*cos_AP - sin_LAN*sin_AP*cos_INCL);
    *e_vec_y = eccentricity*(sin_LAN*cos_AP + cos_LAN*sin_AP*cos_INCL);
    *e_vec_z = eccentricity*(sin_AP*sin_INCL);
    
    *h_vec_x = h*sin_LAN*sin_INCL;
    *h_vec_y = -h*cos_LAN*sin_INCL;
    *h_vec_z = h*cos_INCL;

    return 0;
}

int compute_orbital_elements_from_orbital_vectors(double child1_mass, double child2_mass, double h_tot_vec[3], double e_vec_x, double e_vec_y, double e_vec_z, double h_vec_x, double h_vec_y, double h_vec_z, double *semimajor_axis, double *eccentricity, double *inclination, double *argument_of_pericenter,double *longitude_of_ascending_node)
{
    double e_vec[3] = {e_vec_x,e_vec_y,e_vec_z};
    double h_vec[3] = {h_vec_x,h_vec_y,h_vec_z};
    double eccentricity_squared = norm3_squared(e_vec);
    *eccentricity = sqrt(eccentricity_squared);
    double h_squared = norm3_squared(h_vec);
    *semimajor_axis = h_squared*(child1_mass+child2_mass)/( CONST_G*child1_mass*child1_mass*child2_mass*child2_mass*(1.0 - eccentricity_squared) );
    double h = sqrt(h_squared);
    
//    double x_vec[3] = {1.0,0.0,0.0};
//    double y_vec[3] = {0.0,1.0,0.0};
//    double z_vec[3] = {0.0,0.0,1.0};

    double h_tot = norm3(h_tot_vec);
//    printf("h_tot %g x %g y %g z %g\n",h_tot,h_tot_vec[0],h_tot_vec[1],h_tot_vec[2]);
    double x_vec[3], y_vec[3], z_vec[3];
    for (int i=0; i<3; i++)
    {
        z_vec[i] = h_tot_vec[i]/h_tot;
    }

//    printf("test %g %g %g\n",z_vec[0],z_vec[1],z_vec[2]);
    z_vec[0] = 0.0;
    z_vec[1] = 0.0;
    z_vec[2] = 1.0;

    /* the above assumes that the total angular momentum vector does not change (i.e. no SNe effects etc.) */
    
    double f = 1.0/sqrt( z_vec[0]*z_vec[0] + z_vec[2]*z_vec[2] );
    x_vec[0] = z_vec[2]*f;
    x_vec[1] = 0.0;
    x_vec[2] = -z_vec[0]*f;
    cross3(z_vec,x_vec,y_vec);

    double cos_INCL = dot3(h_vec,z_vec)/h;

    double LAN_vec[3],LAN_vec_unit[3];
    cross3(z_vec,h_vec,LAN_vec);
    double LAN_vec_norm = norm3(LAN_vec);

    double e_vec_unit[3],h_vec_unit[3];

    for (int i=0; i<3; i++)
    {
        LAN_vec_unit[i] = LAN_vec[i]/LAN_vec_norm;
        e_vec_unit[i] = e_vec[i]/(*eccentricity);
        h_vec_unit[i] = h_vec[i]/h;
    }

    double sin_LAN = dot3(LAN_vec_unit,y_vec);
    double cos_LAN = dot3(LAN_vec_unit,x_vec);

    double e_vec_unit_cross_h_vec_unit[3];
    cross3(e_vec_unit,h_vec_unit,e_vec_unit_cross_h_vec_unit);
    double sin_AP = dot3(LAN_vec_unit,e_vec_unit_cross_h_vec_unit);
    double cos_AP = dot3(LAN_vec_unit,e_vec_unit);

    *inclination = acos(cos_INCL);
    *argument_of_pericenter = atan2(sin_AP,cos_AP);
    *longitude_of_ascending_node = atan2(sin_LAN,cos_LAN);
    
    return 0;
}

int get_inclination_relative_to_parent(int index_of_the_particle, double *inclination_relative_to_parent)
{

    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    if (p->is_binary == 0)
    {
        *inclination_relative_to_parent = 0.0;
        return 0;
    }
    if (p->parent == -1)
    {
        *inclination_relative_to_parent = 0.0;
        return 0;
    }

    Particle *parent = particlesMap[p->parent];
    

    double h1_vec[3] = {p->h_vec_x,p->h_vec_y,p->h_vec_z};
    double h2_vec[3] = {parent->h_vec_x,parent->h_vec_y,parent->h_vec_z};

    double h1 = norm3(h1_vec);
    double h2 = norm3(h2_vec);
    
    *inclination_relative_to_parent = acos( dot3(h1_vec,h2_vec)/(h1*h2) );
    
    return 0;
}


void compute_eccentric_anomaly_from_mean_anomaly(double mean_anomaly, double eccentricity, double *cos_eccentric_anomaly, double *sin_eccentric_anomaly)
{
    double eccentric_anomaly;
    double eccentric_anomaly_next = mean_anomaly;
    double epsilon = 1e-10;
    double error = 2.0*epsilon;
    int j = 0;
    while (error > epsilon || j < 15)
    {
        j += 1;
        eccentric_anomaly = eccentric_anomaly_next;
        eccentric_anomaly_next = eccentric_anomaly - (eccentric_anomaly - eccentricity*sin(eccentric_anomaly) - mean_anomaly)/(1.0 - eccentricity*cos(eccentric_anomaly));
        error = fabs(eccentric_anomaly_next - eccentric_anomaly);
    }
    *cos_eccentric_anomaly = cos(eccentric_anomaly);
    *sin_eccentric_anomaly = sin(eccentric_anomaly);
}

void compute_true_anomaly_from_eccentric_anomaly(double cos_eccentric_anomaly, double sin_eccentric_anomaly, double eccentricity, double *cos_true_anomaly, double *sin_true_anomaly)
{
    *cos_true_anomaly = (cos_eccentric_anomaly - eccentricity)/(1.0 - eccentricity*cos_eccentric_anomaly);
    *sin_true_anomaly = sqrt(1.0 - eccentricity*eccentricity)*sin_eccentric_anomaly/(1.0 - eccentricity*cos_eccentric_anomaly);
}

double compute_true_anomaly_from_mean_anomaly(double mean_anomaly, double eccentricity)
{
    double cos_eccentric_anomaly,sin_eccentric_anomaly;
    double cos_true_anomaly,sin_true_anomaly;
    
    compute_eccentric_anomaly_from_mean_anomaly(mean_anomaly,eccentricity,&cos_eccentric_anomaly,&sin_eccentric_anomaly);
    compute_true_anomaly_from_eccentric_anomaly(cos_eccentric_anomaly,sin_eccentric_anomaly,eccentricity,&cos_true_anomaly,&sin_true_anomaly);
    double true_anomaly = atan2(sin_true_anomaly,cos_true_anomaly);

    return true_anomaly;
}

double sample_random_true_anomaly(double eccentricity,int seed)
{
    srand(seed);
    double x = ((double) rand() / (RAND_MAX));
    double mean_anomaly = (2.0*x - 1.0)*M_PI;
    double true_anomaly = compute_true_anomaly_from_mean_anomaly(mean_anomaly,eccentricity);

    return true_anomaly;
}

void from_orbital_vectors_to_cartesian(double child1_mass, double child2_mass, double e_vec[3], double h_vec[3], double true_anomaly, double r[3], double v[3])
{
    double total_mass = child1_mass + child2_mass;
    
    double e = norm3(e_vec);
    double h = norm3(h_vec);

    double e_vec_unit[3],q_vec_unit[3],q_vec[3];
    cross3(h_vec,e_vec,q_vec);
    double q = norm3(q_vec);

    int i;
    for (i=0; i<3; i++)
    {        
        e_vec_unit[i] = e_vec[i]/e;
        q_vec_unit[i] = q_vec[i]/q;        
    }
    
    double e_p2 = e*e;
    double j_p2 = 1.0 - e_p2;
   
    double a = h*h*total_mass/( CONST_G*child1_mass*child2_mass*child1_mass*child2_mass*j_p2 );

    double cos_f = cos(true_anomaly);
    double sin_f = sin(true_anomaly);
    
    double r_norm = a*j_p2/(1.0 + e*cos_f);
    double v_norm = sqrt( CONST_G*total_mass/(a*j_p2) );
    
    for (i=0; i<3; i++)
    {
        r[i] = r_norm*( cos_f*e_vec_unit[i] + sin_f*q_vec_unit[i]);
        v[i] = v_norm*( -sin_f*e_vec_unit[i] + (e + cos_f)*q_vec_unit[i] );
    }
}

void from_cartesian_to_orbital_vectors(double child1_mass, double child2_mass, double r[3], double v[3], double e_vec[3], double h_vec[3])
{
    double total_mass = child1_mass + child2_mass;
       
    double v_dot_v = dot3(v,v);
    double r_dot_v = dot3(r,v);
    double r_norm = norm3(r);
    for (int i=0; i<3; i++)
    {
        e_vec[i] = (r[i]*v_dot_v - v[i]*r_dot_v)/(CONST_G*total_mass) - r[i]/r_norm;
    }

    double mu = child1_mass*child2_mass/total_mass;
    cross3(r,v,h_vec);
    for (int i=0; i<3; i++)
    {
        h_vec[i] *= mu;
    }
}


int get_de_dt(int index_of_the_particle, double *de_dt)
{

    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    if (p->is_binary == 0)
    {
        *de_dt = 0.0;
        return 0;
    }

    *de_dt = dot3(p->e_vec_unit,p->de_vec_dt);

    return 0;
}



void get_position_and_velocity_vectors_from_particle(Particle *p, double r[3], double v[3])
{
    r[0] = p->position_x;
    r[1] = p->position_y;
    r[2] = p->position_z;
    v[0] = p->velocity_x;
    v[1] = p->velocity_y;
    v[2] = p->velocity_z;
}
void set_position_and_velocity_vectors_in_particle(Particle *p,  double r[3], double v[3])
{
    p->position_x = r[0];
    p->position_y = r[1];
    p->position_z = r[2];
    p->velocity_x = v[0];
    p->velocity_y = v[1];
    p->velocity_z = v[2];
}
void get_e_and_h_vectors_from_particle(Particle *p, double e_vec[3], double h_vec[3])
{
    e_vec[0] = p->e_vec_x;
    e_vec[1] = p->e_vec_y;
    e_vec[2] = p->e_vec_z;
    h_vec[0] = p->h_vec_x;    
    h_vec[1] = p->h_vec_y;    
    h_vec[2] = p->h_vec_z;    
}
void set_e_and_h_vectors_in_particle(Particle *p, double e_vec[3], double h_vec[3])
{
    p->e_vec_x = e_vec[0];
    p->e_vec_y = e_vec[1];
    p->e_vec_z = e_vec[2];
    p->h_vec_x = h_vec[0];
    p->h_vec_y = h_vec[1];
    p->h_vec_z = h_vec[2];
}


/************************
/* interface parameters *
 ************************/
 
int get_relative_tolerance(double *value)
{
    *value = relative_tolerance;
    return 0;
}
int set_relative_tolerance(double value)
{
    relative_tolerance = value;
    return 0;
}
int get_absolute_tolerance_eccentricity_vectors(double *value)
{
    *value = absolute_tolerance_eccentricity_vectors;
    return 0;
}
int set_absolute_tolerance_eccentricity_vectors(double value)
{
    absolute_tolerance_eccentricity_vectors = value;
    return 0;
}

int get_include_quadrupole_order_terms(bool *value){
    *value = include_quadrupole_order_terms;
    return 0;
}
int set_include_quadrupole_order_terms(bool value){
    include_quadrupole_order_terms = value;
    return 0;
}

int get_include_octupole_order_binary_pair_terms(bool *value){
    *value = include_octupole_order_binary_pair_terms;
    return 0;
}
int set_include_octupole_order_binary_pair_terms(bool value){
    include_octupole_order_binary_pair_terms = value;
    return 0;
}

int get_include_octupole_order_binary_triplet_terms(bool *value){
    *value = include_octupole_order_binary_triplet_terms;
    return 0;
}
int set_include_octupole_order_binary_triplet_terms(bool value){
    include_octupole_order_binary_triplet_terms = value;
    return 0;
}

int get_include_hexadecupole_order_binary_pair_terms(bool *value){
    *value = include_hexadecupole_order_binary_pair_terms;
    return 0;
}
int set_include_hexadecupole_order_binary_pair_terms(bool value){
    include_hexadecupole_order_binary_pair_terms = value;
    return 0;
}

int get_include_dotriacontupole_order_binary_pair_terms(bool *value){
    *value = include_dotriacontupole_order_binary_pair_terms;
    return 0;
}
int set_include_dotriacontupole_order_binary_pair_terms(bool value){
    include_dotriacontupole_order_binary_pair_terms = value;
    return 0;
}

int get_orbital_phases_random_seed(int *value)
{
    *value = orbital_phases_random_seed;
    return 0;
}
int set_orbital_phases_random_seed(int value)
{
    orbital_phases_random_seed = value;
    return 0;
}
