#include "src/types.h"

/*******************
/* basic interface *
 ******************/
int new_particle(int * index_of_the_particle, int is_binary);
int delete_particle(int index_of_the_particle);

int set_children(int index_of_the_particle, int child1, int child2);
int get_children(int index_of_the_particle, int *child1, int *child2);

int set_mass(int index_of_the_particle, double value);
int get_mass(int index_of_the_particle, double *value);

int set_radius(int index_of_the_particle, double value);
int get_radius(int index_of_the_particle, double *value);

int get_level(int index_of_the_particle, int *level);


/****************
/* spin vectors *
 ****************/
int set_spin_vector(int index_of_the_particle, double spin_vec_x, double spin_vec_y, double spin_vec_z);
int get_spin_vector(int index_of_the_particle, double *spin_vec_x, double *spin_vec_y, double *spin_vec_z);


/****************************
/* orbital vectors/elements *
 ****************************/
int set_orbital_vectors(int index_of_the_particle, double e_vec_x, double e_vec_y, double e_vec_z, \
    double h_vec_x, double h_vec_y, double h_vec_z);
int get_orbital_vectors(int index_of_the_particle, double *e_vec_x, double *e_vec_y, double *e_vec_z, \
    double *h_vec_x, double *h_vec_y, double *h_vec_z);
    
int set_orbital_elements(int index_of_the_particle, double semimajor_axis, double eccentricity, \
    double inclination, double argument_of_pericenter, double longitude_of_ascending_node);
int get_orbital_elements(int index_of_the_particle, double *semimajor_axis, double *eccentricity, \
    double *inclination, double *argument_of_pericenter, double *longitude_of_ascending_node);
    
    
/************
/* PN terms *
 ************/
int set_include_pairwise_1PN_terms(int index_of_the_particle, int include_pairwise_1PN_terms);
int get_include_pairwise_1PN_terms(int index_of_the_particle, int *include_pairwise_1PN_terms);
int set_include_pairwise_25PN_terms(int index_of_the_particle, int include_pairwise_25PN_terms);
int get_include_pairwise_25PN_terms(int index_of_the_particle, int *include_pairwise_25PN_terms);

/*********
/* tides *
 *********/
int set_include_tidal_friction_terms(int index_of_the_particle, int include_tidal_friction_terms);
int get_include_tidal_friction_terms(int index_of_the_particle, int *include_tidal_friction_terms);
int set_tides_method(int index_of_the_particle, int tides_method);
int get_tides_method(int index_of_the_particle, int *tides_method);
int set_include_tidal_bulges_precession_terms(int index_of_the_particle, int include_tidal_bulges_precession_terms);
int get_include_tidal_bulges_precession_terms(int index_of_the_particle, int *include_tidal_bulges_precession_terms);
int set_include_rotation_precession_terms(int index_of_the_particle, int include_rotation_precession_terms);
int get_include_rotation_precession_terms(int index_of_the_particle, int *include_rotation_precession_terms);
int set_minimum_eccentricity_for_tidal_precession(int index_of_the_particle, double minimum_eccentricity_for_tidal_precession);
int get_minimum_eccentricity_for_tidal_precession(int index_of_the_particle, double *minimum_eccentricity_for_tidal_precession);
int set_tides_apsidal_motion_constant(int index_of_the_particle, double value);
int get_tides_apsidal_motion_constant(int index_of_the_particle, double *value);
int set_tides_gyration_radius(int index_of_the_particle, double value);
int get_tides_gyration_radius(int index_of_the_particle, double *value);
int set_tides_viscous_time_scale(int index_of_the_particle, double value);
int get_tides_viscous_time_scale(int index_of_the_particle, double *value);


/****************
/* root finding *
 ****************/
int set_check_for_secular_breakdown(int index_of_the_particle, int value);
int get_check_for_secular_breakdown(int index_of_the_particle, int* value);

int set_check_for_dynamical_instability(int index_of_the_particle, int value);
int get_check_for_dynamical_instability(int index_of_the_particle, int* value);
int set_dynamical_instability_criterion(int index_of_the_particle, int value);
int get_dynamical_instability_criterion(int index_of_the_particle, int* value);
int set_dynamical_instability_central_particle(int index_of_the_particle, int value);
int get_dynamical_instability_central_particle(int index_of_the_particle, int* value);
int set_dynamical_instability_K_parameter(int index_of_the_particle, double value);
int get_dynamical_instability_K_parameter(int index_of_the_particle, double* value);

int set_check_for_physical_collision_or_orbit_crossing(int index_of_the_particle, int value);
int get_check_for_physical_collision_or_orbit_crossing(int index_of_the_particle, int* value);

int set_check_for_minimum_periapse_distance(int index_of_the_particle, int value);
int get_check_for_minimum_periapse_distance(int index_of_the_particle, int* value);
int set_check_for_minimum_periapse_distance_value(int index_of_the_particle, double value);
int get_check_for_minimum_periapse_distance_value(int index_of_the_particle, double* value);

int set_check_for_RLOF_at_pericentre(int index_of_the_particle, int value);
int get_check_for_RLOF_at_pericentre(int index_of_the_particle, int* value);
int set_check_for_RLOF_at_pericentre_use_sepinsky_fit(int index_of_the_particle, int value);
int get_check_for_RLOF_at_pericentre_use_sepinsky_fit(int index_of_the_particle, int* value);

int set_root_finding_state(int index_of_the_particle, int secular_breakdown_has_occurred, int dynamical_instability_has_occurred, int physical_collision_or_orbit_crossing_has_occurred, int minimum_periapse_distance_has_occurred, int RLOF_at_pericentre_has_occurred);
int get_root_finding_state(int index_of_the_particle, int *secular_breakdown_has_occurred, int *dynamical_instability_has_occurred, int *physical_collision_or_orbit_crossing_has_occurred, int* minimum_periapse_distance_has_occurred, int *RLOF_at_pericentre_has_occurred);


/********************
/* evolve interface *
 ********************/
int evolve_interface(double start_time, double time_step, double *output_time, double *hamiltonian, int *flag, int *error_code);
int determine_binary_parents_levels_and_masses_interface();


/**********************************************
/* orbital element/vector conversion routines *
 **********************************************/
int compute_h_tot_vector(ParticlesMap* particlesMap, double h_tot_vector[3]);
int compute_orbital_vectors_from_orbital_elements(double child1_mass, double child2_mass, double semimajor_axis, double eccentricity, double inclination, double argument_of_pericenter,double longitude_of_ascending_node, double *e_vec_x, double *e_vec_y, double *e_vec_z, double *h_vec_x, double *h_vec_y, double *h_vec_z);
int compute_orbital_elements_from_orbital_vectors(double child1_mass, double child2_mass, double h_tot_vec[3], double e_vec_x, double e_vec_y, double e_vec_z, double h_vec_x, double h_vec_y, double h_vec_z, double *semimajor_axis, double *eccentricity, double *inclination, double *argument_of_pericenter,double *longitude_of_ascending_node);
int get_inclination_relative_to_parent(int index_of_the_particle, double *inclination_relative_to_parent);
int get_de_dt(int index_of_the_particle, double *de_dt);

/************************
/* interface parameters *
 ************************/
extern double relative_tolerance;
extern double absolute_tolerance_eccentricity_vectors;
extern double absolute_tolerance_angular_momentum_vectors;
extern double absolute_tolerance_spin_vectors;
extern bool compute_orbital_elements_with_respect_to_total_angular_momentum_vector;
extern bool include_quadrupole_order_terms;
extern bool include_octupole_order_binary_pair_terms;
extern bool include_octupole_order_binary_triplet_terms;
extern bool include_hexadecupole_order_binary_pair_terms;
extern bool include_dotriacontupole_order_binary_pair_terms;

int get_relative_tolerance(double *value);
int set_relative_tolerance(double value);

int get_absolute_tolerance_eccentricity_vectors(double *value);
int set_absolute_tolerance_eccentricity_vectors(double value);

int get_absolute_tolerance_angular_momentum_vectors(double *value);
int set_absolute_tolerance_angular_momentum_vectors(double value);

int get_absolute_tolerance_spin_vectors(double *value);
int set_absolute_tolerance_spin_vectors(double value);

int get_compute_orbital_elements_with_respect_to_total_angular_momentum_vector(int *value);
int set_compute_orbital_elements_with_respect_to_total_angular_momentum_vector(int value);

int get_include_quadrupole_order_terms(int *value);
int set_include_quadrupole_order_terms(int value);

int get_include_octupole_order_binary_pair_terms(int *value);
int set_include_octupole_order_binary_pair_terms(int value);

int get_include_octupole_order_binary_triplet_terms(int *value);
int set_include_octupole_order_binary_triplet_terms(int value);

int get_include_hexadecupole_order_binary_pair_terms(int *value);
int set_include_hexadecupole_order_binary_pair_terms(int value);

int get_include_dotriacontupole_order_binary_pair_terms(int *value);
int set_include_dotriacontupole_order_binary_pair_terms(int value);
