#include "src/types.h"

/*******************
/* basic interface *
 ******************/
int new_particle(int * index_of_the_particle, bool is_binary);
int delete_particle(int index_of_the_particle);

int set_children(int index_of_the_particle, int child1, int child2);
int get_children(int index_of_the_particle, int *child1, int *child2);

int set_mass(int index_of_the_particle, double value);
int get_mass(int index_of_the_particle, double *value);

int set_mass_dot_external(int index_of_the_particle, double value);
int get_mass_dot_external(int index_of_the_particle, double *value);

int set_radius(int index_of_the_particle, double value);
int get_radius(int index_of_the_particle, double *value);

int set_radius_dot_external(int index_of_the_particle, double value);
int get_radius_dot_external(int index_of_the_particle, double *value);

int set_radius_ddot_external(int index_of_the_particle, double value);
int get_radius_ddot_external(int index_of_the_particle, double *value);

int get_level(int index_of_the_particle, int *level);

int set_stellar_type(int index_of_the_particle, int value);
int get_stellar_type(int index_of_the_particle, int *stellar_type);

int set_true_anomaly(int index_of_the_particle, double value);
int get_true_anomaly(int index_of_the_particle, double *value);

int set_sample_orbital_phases_randomly(int index_of_the_particle, bool value);
int get_sample_orbital_phases_randomly(int index_of_the_particle, bool *value);


/*******************************
 * instantaneous perturbations *
 * ****************************/

int set_instantaneous_perturbation_delta_mass(int index_of_the_particle, double value);
int get_instantaneous_perturbation_delta_mass(int index_of_the_particle, double *value);

int set_instantaneous_perturbation_delta_position(int index_of_the_particle, double x, double y, double z);
int get_instantaneous_perturbation_delta_position(int index_of_the_particle, double *x, double *y, double *z);

int set_instantaneous_perturbation_delta_velocity(int index_of_the_particle, double x, double y, double z);
int get_instantaneous_perturbation_delta_velocity(int index_of_the_particle, double *x, double *y, double *z);


/************
 * external *
 * *********/
 
int new_external_particle(int * index_of_the_particle, double mass);
int delete_external_particle(int index_of_the_particle);

int set_external_mass(int index_of_the_particle, double value);
int get_external_mass(int index_of_the_particle, double *value);

int set_external_path(int index_of_the_external_particle, int value);
int get_external_path(int index_of_the_external_particle, int *value);

int set_external_mode(int index_of_the_external_particle, int value);
int get_external_mode(int index_of_the_external_particle, int *value);

int set_external_t_ref(int index_of_the_particle, double value);
int get_external_t_ref(int index_of_the_particle, double *value);

int set_external_t_passed(int index_of_the_external_particle, double value);
int get_external_t_passed(int index_of_the_external_particle, double *value);

int set_external_r0_vectors(int index_of_the_particle, double vec_x, double vec_y, double vec_z);
int get_external_r0_vectors(int index_of_the_particle, double *vec_x, double *vec_y, double *vec_z);

int set_external_rdot_vectors(int index_of_the_external_particle, double rdot_vec_x, double rdot_vec_y, double rdot_vec_z);
int get_external_rdot_vectors(int index_of_the_external_particle, double *rdot_vec_x, double *rdot_vec_y, double *rdot_vec_z);

int set_external_periapse_distance(int index_of_the_external_particle, double value);
int get_external_periapse_distance(int index_of_the_external_particle, double *value);

int set_external_eccentricity(int index_of_the_external_particle, double value);
int get_external_eccentricity(int index_of_the_external_particle, double *value);

int set_external_e_hat_vectors(int index_of_the_external_particle, double x, double y, double z);
int get_external_e_hat_vectors(int index_of_the_external_particle, double *x, double *y, double *z);

int set_external_h_hat_vectors(int index_of_the_external_particle, double x, double y, double z);
int get_external_h_hat_vectors(int index_of_the_external_particle, double *x, double *y, double *z);

int get_external_r_vectors(int index_of_the_external_particle, double *r_vec_x, double *r_vec_y, double *r_vec_z);

/****************
/* spin vectors *
 ****************/
int set_spin_vector(int index_of_the_particle, double spin_vec_x, double spin_vec_y, double spin_vec_z);
int get_spin_vector(int index_of_the_particle, double *spin_vec_x, double *spin_vec_y, double *spin_vec_z);

int set_spin_vector_dot_external(int index_of_the_particle, double spin_vec_x_dot_external, double spin_vec_y_dot_external, double spin_vec_z_dot_external);
int get_spin_vector_dot_external(int index_of_the_particle, double *spin_vec_x_dot_external, double *spin_vec_y_dot_external, double *spin_vec_z_dot_external);

/****************************
/* orbital vectors/elements *
 ****************************/
int set_orbital_vectors(int index_of_the_particle, double e_vec_x, double e_vec_y, double e_vec_z, \
    double h_vec_x, double h_vec_y, double h_vec_z);
int get_orbital_vectors(int index_of_the_particle, double *e_vec_x, double *e_vec_y, double *e_vec_z, \
    double *h_vec_x, double *h_vec_y, double *h_vec_z);
    
int set_orbital_vectors_dot_external(int index_of_the_particle, double e_vec_x_dot_external, double e_vec_y_dot_external, double e_vec_z_dot_external, \
    double h_vec_x_dot_external, double h_vec_y_dot_external, double h_vec_z_dot_external);
int get_orbital_vectors_dot_external(int index_of_the_particle, double *e_vec_x_dot_external, double *e_vec_y_dot_external, double *e_vec_z_dot_external, \
    double *h_vec_x_dot_external, double *h_vec_y_dot_external, double *h_vec_z_dot_external);

int set_orbital_elements(int index_of_the_particle, double semimajor_axis, double eccentricity, \
    double inclination, double argument_of_pericenter, double longitude_of_ascending_node);
int get_orbital_elements(int index_of_the_particle, double *semimajor_axis, double *eccentricity, \
    double *inclination, double *argument_of_pericenter, double *longitude_of_ascending_node);


int set_position_vector(int index_of_the_particle, double x, double y, double z);
int get_position_vector(int index_of_the_particle, double *x, double *y, double *z);

int set_velocity_vector(int index_of_the_particle, double x, double y, double z);
int get_velocity_vector(int index_of_the_particle, double *x, double *y, double *z);

    
/************
/* PN terms *
 ************/
int set_include_pairwise_1PN_terms(int index_of_the_particle, bool include_pairwise_1PN_terms);
int get_include_pairwise_1PN_terms(int index_of_the_particle, bool *include_pairwise_1PN_terms);
int set_include_pairwise_25PN_terms(int index_of_the_particle, bool include_pairwise_25PN_terms);
int get_include_pairwise_25PN_terms(int index_of_the_particle, bool *include_pairwise_25PN_terms);

/*********
/* tides *
 *********/
int set_include_tidal_friction_terms(int index_of_the_particle, bool include_tidal_friction_terms);
int get_include_tidal_friction_terms(int index_of_the_particle, bool *include_tidal_friction_terms);
int set_tides_method(int index_of_the_particle, int tides_method);
int get_tides_method(int index_of_the_particle, int *tides_method);
int set_include_tidal_bulges_precession_terms(int index_of_the_particle, bool include_tidal_bulges_precession_terms);
int get_include_tidal_bulges_precession_terms(int index_of_the_particle, bool *include_tidal_bulges_precession_terms);
int set_include_rotation_precession_terms(int index_of_the_particle, bool include_rotation_precession_terms);
int get_include_rotation_precession_terms(int index_of_the_particle, bool *include_rotation_precession_terms);
int set_minimum_eccentricity_for_tidal_precession(int index_of_the_particle, double minimum_eccentricity_for_tidal_precession);
int get_minimum_eccentricity_for_tidal_precession(int index_of_the_particle, double *minimum_eccentricity_for_tidal_precession);
int set_tides_apsidal_motion_constant(int index_of_the_particle, double value);
int get_tides_apsidal_motion_constant(int index_of_the_particle, double *value);
int set_tides_gyration_radius(int index_of_the_particle, double value);
int get_tides_gyration_radius(int index_of_the_particle, double *value);
int set_tides_viscous_time_scale(int index_of_the_particle, double value);
int get_tides_viscous_time_scale(int index_of_the_particle, double *value);
int set_tides_viscous_time_scale_prescription(int index_of_the_particle, int value);
int get_tides_viscous_time_scale_prescription(int index_of_the_particle, int *value);
int set_convective_envelope_mass(int index_of_the_particle, double value);
int get_convective_envelope_mass(int index_of_the_particle, double *value);
int set_convective_envelope_radius(int index_of_the_particle, double value);
int get_convective_envelope_radius(int index_of_the_particle, double *value);
int set_luminosity(int index_of_the_particle, double luminosity);
int get_luminosity(int index_of_the_particle, double *luminosity);


/****************
/* root finding *
 ****************/

int set_check_for_secular_breakdown(int index_of_the_particle, bool value);
int get_check_for_secular_breakdown(int index_of_the_particle, bool* value);

int set_check_for_dynamical_instability(int index_of_the_particle, bool value);
int get_check_for_dynamical_instability(int index_of_the_particle, bool* value);
int set_dynamical_instability_criterion(int index_of_the_particle, int value);
int get_dynamical_instability_criterion(int index_of_the_particle, int* value);
int set_dynamical_instability_central_particle(int index_of_the_particle, int value);
int get_dynamical_instability_central_particle(int index_of_the_particle, int* value);
int set_dynamical_instability_K_parameter(int index_of_the_particle, double value);
int get_dynamical_instability_K_parameter(int index_of_the_particle, double* value);

int set_check_for_physical_collision_or_orbit_crossing(int index_of_the_particle, bool value);
int get_check_for_physical_collision_or_orbit_crossing(int index_of_the_particle, bool* value);

int set_check_for_minimum_periapse_distance(int index_of_the_particle, bool value);
int get_check_for_minimum_periapse_distance(int index_of_the_particle, bool* value);
int set_check_for_minimum_periapse_distance_value(int index_of_the_particle, double value);
int get_check_for_minimum_periapse_distance_value(int index_of_the_particle, double* value);

int set_check_for_RLOF_at_pericentre(int index_of_the_particle, bool value);
int get_check_for_RLOF_at_pericentre(int index_of_the_particle, bool* value);
int set_check_for_RLOF_at_pericentre_use_sepinsky_fit(int index_of_the_particle, bool value);
int get_check_for_RLOF_at_pericentre_use_sepinsky_fit(int index_of_the_particle, bool* value);

int set_root_finding_state(int index_of_the_particle, bool secular_breakdown_has_occurred, bool dynamical_instability_has_occurred, bool physical_collision_or_orbit_crossing_has_occurred, bool minimum_periapse_distance_has_occurred, bool RLOF_at_pericentre_has_occurred);
int get_root_finding_state(int index_of_the_particle, bool *secular_breakdown_has_occurred, bool *dynamical_instability_has_occurred, bool *physical_collision_or_orbit_crossing_has_occurred, bool* minimum_periapse_distance_has_occurred, bool *RLOF_at_pericentre_has_occurred);


/***********************
/* interface functions *
 ***********************/
int evolve_interface(double start_time, double time_step, double *output_time, double *hamiltonian, int *flag, int *error_code);
int determine_binary_parents_levels_and_masses_interface();
int apply_external_perturbation_assuming_integrated_orbits_interface();
int apply_user_specified_instantaneous_perturbation_interface();
int set_positions_and_velocities_interface();

/**********************************************
/* orbital element/vector conversion routines *
 **********************************************/
void compute_h_tot_vector(ParticlesMap* particlesMap, double h_tot_vector[3]);
int compute_orbital_vectors_from_orbital_elements(double child1_mass, double child2_mass, double semimajor_axis, double eccentricity, double inclination, double argument_of_pericenter,double longitude_of_ascending_node, double *e_vec_x, double *e_vec_y, double *e_vec_z, double *h_vec_x, double *h_vec_y, double *h_vec_z);
int compute_orbital_elements_from_orbital_vectors(double child1_mass, double child2_mass, double h_tot_vec[3], double e_vec_x, double e_vec_y, double e_vec_z, double h_vec_x, double h_vec_y, double h_vec_z, double *semimajor_axis, double *eccentricity, double *inclination, double *argument_of_pericenter,double *longitude_of_ascending_node);
int get_inclination_relative_to_parent(int index_of_the_particle, double *inclination_relative_to_parent);

void compute_eccentric_anomaly_from_mean_anomaly(double mean_anomaly, double eccentricity, double *cos_eccentric_anomaly, double *sin_eccentric_anomaly);
void compute_true_anomaly_from_eccentric_anomaly(double cos_eccentric_anomaly, double sin_eccentric_anomaly, double eccentricity, double *cos_true_anomaly, double *sin_true_anomaly);
double compute_true_anomaly_from_mean_anomaly(double mean_anomaly, double eccentricity);
double sample_random_true_anomaly(double eccentricity,int seed);

void from_orbital_vectors_to_cartesian(double child1_mass, double child2_mass, double e_vec[3], double h_vec[3], double true_anomaly, double r[3], double v[3]);
void from_cartesian_to_orbital_vectors(double child1_mass, double child2_mass, double r[3], double v[3], double e_vec[3], double h_vec[3]);

int get_de_dt(int index_of_the_particle, double *de_dt);

void get_position_and_velocity_vectors_from_particle(Particle *p, double r[3], double v[3]);
void set_position_and_velocity_vectors_in_particle(Particle *p,  double r[3], double v[3]);
void get_e_and_h_vectors_from_particle(Particle *p, double e_vec[3], double h_vec[3]);
void set_e_and_h_vectors_in_particle(Particle *p, double e_vec[3], double h_vec[3]);

/************************
/* interface parameters *
 ************************/
extern double relative_tolerance;
extern double absolute_tolerance_eccentricity_vectors;
extern bool include_quadrupole_order_terms;
extern bool include_octupole_order_binary_pair_terms;
extern bool include_octupole_order_binary_triplet_terms;
extern bool include_hexadecupole_order_binary_pair_terms;
extern bool include_dotriacontupole_order_binary_pair_terms;
extern int orbital_phases_random_seed;

int get_relative_tolerance(double *value);
int set_relative_tolerance(double value);

int get_absolute_tolerance_eccentricity_vectors(double *value);
int set_absolute_tolerance_eccentricity_vectors(double value);

int get_include_quadrupole_order_terms(bool *value);
int set_include_quadrupole_order_terms(bool value);

int get_include_octupole_order_binary_pair_terms(bool *value);
int set_include_octupole_order_binary_pair_terms(bool value);

int get_include_octupole_order_binary_triplet_terms(bool *value);
int set_include_octupole_order_binary_triplet_terms(bool value);

int get_include_hexadecupole_order_binary_pair_terms(bool *value);
int set_include_hexadecupole_order_binary_pair_terms(bool value);

int get_include_dotriacontupole_order_binary_pair_terms(bool *value);
int set_include_dotriacontupole_order_binary_pair_terms(bool value);

int get_orbital_phases_random_seed(int *value);
int set_orbital_phases_random_seed(int value);
