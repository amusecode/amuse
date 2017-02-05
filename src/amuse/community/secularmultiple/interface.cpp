#include "src/types.h"
#include "interface.h"
#include "src/evolve.h"


int highest_particle_index = 0;
ParticlesMap particlesMap;

double relative_tolerance = 1.0e-16;
double absolute_tolerance_eccentricity_vectors = 1.0e-14;
double absolute_tolerance_angular_momentum_vectors = 1.0e-2;
double absolute_tolerance_spin_vectors = 1.0e4;
bool compute_orbital_elements_with_respect_to_total_angular_momentum_vector = false;
bool include_quadrupole_order_terms = true;
bool include_octupole_order_binary_pair_terms = true;
bool include_octupole_order_binary_triplet_terms = false;
bool include_hexadecupole_order_binary_pair_terms = false;
bool include_dotriacontupole_order_binary_pair_terms = false;

/*******************
/* basic interface *
 ******************/
 
int new_particle(int * index_of_the_particle, int is_binary)
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
    determine_binary_parents_levels_and_masses(&particlesMap, &N_bodies, &N_binaries, &N_root_finding);
    
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
    determine_binary_parents_levels_and_masses(&particlesMap, &N_bodies, &N_binaries,&N_root_finding);
    
    compute_orbital_elements_from_orbital_vectors(p->child1_mass, p->child2_mass, h_tot_vec, \
        p->e_vec_x,p->e_vec_y,p->e_vec_z,p->h_vec_x,p->h_vec_y,p->h_vec_z,
        semimajor_axis, eccentricity, inclination, argument_of_pericenter, longitude_of_ascending_node);
    
    return 0;
}


/************
/* PN terms *
 ************/

int set_include_pairwise_1PN_terms(int index_of_the_particle, int include_pairwise_1PN_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    
    p->include_pairwise_1PN_terms = include_pairwise_1PN_terms;
        
    return 0;
}
int get_include_pairwise_1PN_terms(int index_of_the_particle, int *include_pairwise_1PN_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    
    *include_pairwise_1PN_terms = p->include_pairwise_1PN_terms;
        
    return 0;
}

int set_include_pairwise_25PN_terms(int index_of_the_particle, int include_pairwise_25PN_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    
    p->include_pairwise_25PN_terms = include_pairwise_25PN_terms;
        
    return 0;
}
int get_include_pairwise_25PN_terms(int index_of_the_particle, int *include_pairwise_25PN_terms)
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
int set_include_tidal_friction_terms(int index_of_the_particle, int include_tidal_friction_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    
    p->include_tidal_friction_terms = include_tidal_friction_terms;
        
    return 0;
}
int get_include_tidal_friction_terms(int index_of_the_particle, int *include_tidal_friction_terms)
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

int set_include_tidal_bulges_precession_terms(int index_of_the_particle, int include_tidal_bulges_precession_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    
    p->include_tidal_bulges_precession_terms = include_tidal_bulges_precession_terms;
        
    return 0;
}
int get_include_tidal_bulges_precession_terms(int index_of_the_particle, int *include_tidal_bulges_precession_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    
    *include_tidal_bulges_precession_terms = p->include_tidal_bulges_precession_terms;
        
    return 0;
}

int set_include_rotation_precession_terms(int index_of_the_particle, int include_rotation_precession_terms)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index_of_the_particle];
    
    p->include_rotation_precession_terms = include_rotation_precession_terms;
        
    return 0;
}

int get_include_rotation_precession_terms(int index_of_the_particle, int *include_rotation_precession_terms)
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


/****************
/* root finding *
 ****************/

/* secular breakdown*/
int set_check_for_secular_breakdown(int index_of_the_particle, int value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->check_for_secular_breakdown = value;

    return 0;
}
int get_check_for_secular_breakdown(int index_of_the_particle, int* value)
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
int set_check_for_dynamical_instability(int index_of_the_particle, int value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->check_for_dynamical_instability = value;

    return 0;
}
int get_check_for_dynamical_instability(int index_of_the_particle, int* value)
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
int set_check_for_physical_collision_or_orbit_crossing(int index_of_the_particle, int value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->check_for_physical_collision_or_orbit_crossing = value;

    return 0;
}
int get_check_for_physical_collision_or_orbit_crossing(int index_of_the_particle, int* value)
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
int set_check_for_minimum_periapse_distance(int index_of_the_particle, int value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->check_for_minimum_periapse_distance = value;

    return 0;
}
int get_check_for_minimum_periapse_distance(int index_of_the_particle, int* value)
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
int set_check_for_RLOF_at_pericentre(int index_of_the_particle, int value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->check_for_RLOF_at_pericentre = value;

    return 0;
}
int get_check_for_RLOF_at_pericentre(int index_of_the_particle, int* value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    *value = p->check_for_RLOF_at_pericentre;

    return 0;
}

int set_check_for_RLOF_at_pericentre_use_sepinsky_fit(int index_of_the_particle, int value)
{
    if (index_of_the_particle > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index_of_the_particle];
    p->check_for_RLOF_at_pericentre_use_sepinsky_fit = value;

    return 0;
}
int get_check_for_RLOF_at_pericentre_use_sepinsky_fit(int index_of_the_particle, int* value)
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
int set_root_finding_state(int index_of_the_particle, int secular_breakdown_has_occurred, int dynamical_instability_has_occurred, int physical_collision_or_orbit_crossing_has_occurred, int minimum_periapse_distance_has_occurred, int RLOF_at_pericentre_has_occurred)
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
int get_root_finding_state(int index_of_the_particle, int *secular_breakdown_has_occurred, int *dynamical_instability_has_occurred, int *physical_collision_or_orbit_crossing_has_occurred, int* minimum_periapse_distance_has_occurred, int *RLOF_at_pericentre_has_occurred)
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
    int result = evolve(&particlesMap, start_time, time_step, output_time, hamiltonian, flag, error_code);
    
    return result;
}

/* set levels and masses */
int determine_binary_parents_levels_and_masses_interface()
{
    
    int N_bodies, N_binaries, N_root_finding;
    determine_binary_parents_levels_and_masses(&particlesMap,&N_bodies,&N_binaries,&N_root_finding);
    
    return 0;
}

/**********************************************
/* orbital element/vector conversion routines *
 **********************************************/
int compute_h_tot_vector(ParticlesMap* particlesMap, double h_tot_vec[3])
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
    
    double h_tot = norm3(h_tot_vec);

    /* The reference coordinate frame is given by x_vec[3],y_vec[3] and z_vec[3] */
    double x_vec[3], y_vec[3], z_vec[3];

    z_vec[0] = 0.0;
    z_vec[1] = 0.0;
    z_vec[2] = 1.0;

    if (compute_orbital_elements_with_respect_to_total_angular_momentum_vector == true)
    {
        for (int i=0; i<3; i++)
        {
            z_vec[i] = h_tot_vec[i]/h_tot;
        }
    }
    
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

int get_absolute_tolerance_angular_momentum_vectors(double *value)
{
    *value = absolute_tolerance_angular_momentum_vectors;
    return 0;
}
int set_absolute_tolerance_angular_momentum_vectors(double value)
{
    absolute_tolerance_angular_momentum_vectors = value;
    return 0;
}

int get_absolute_tolerance_spin_vectors(double *value)
{
    *value = absolute_tolerance_spin_vectors;
    return 0;
}
int set_absolute_tolerance_spin_vectors(double value)
{
    absolute_tolerance_spin_vectors = value;
    return 0;
}

int get_compute_orbital_elements_with_respect_to_total_angular_momentum_vector(int *value){
    *value = compute_orbital_elements_with_respect_to_total_angular_momentum_vector ? 1 : 0;
    return 0;
}
int set_compute_orbital_elements_with_respect_to_total_angular_momentum_vector(int value){
    compute_orbital_elements_with_respect_to_total_angular_momentum_vector = value == 1;
    return 0;
}

int get_include_quadrupole_order_terms(int *value){
    *value = include_quadrupole_order_terms ? 1 : 0;
    return 0;
}
int set_include_quadrupole_order_terms(int value){
    include_quadrupole_order_terms = value == 1;
    return 0;
}

int get_include_octupole_order_binary_pair_terms(int *value){
    *value = include_octupole_order_binary_pair_terms ? 1 : 0;
    return 0;
}
int set_include_octupole_order_binary_pair_terms(int value){
    include_octupole_order_binary_pair_terms = value == 1;
    return 0;
}

int get_include_octupole_order_binary_triplet_terms(int *value){
    *value = include_octupole_order_binary_triplet_terms ? 1 : 0;
    return 0;
}
int set_include_octupole_order_binary_triplet_terms(int value){
    include_octupole_order_binary_triplet_terms = value == 1;
    return 0;
}

int get_include_hexadecupole_order_binary_pair_terms(int *value){
    *value = include_hexadecupole_order_binary_pair_terms ? 1 : 0;
    return 0;
}
int set_include_hexadecupole_order_binary_pair_terms(int value){
    include_hexadecupole_order_binary_pair_terms = value == 1;
    return 0;
}

int get_include_dotriacontupole_order_binary_pair_terms(int *value){
    *value = include_dotriacontupole_order_binary_pair_terms ? 1 : 0;
    return 0;
}
int set_include_dotriacontupole_order_binary_pair_terms(int value){
    include_dotriacontupole_order_binary_pair_terms = value == 1;
    return 0;
}
