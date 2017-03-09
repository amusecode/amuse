#include "types.h"

void extract_ODE_variables(ParticlesMap *particlesMap, N_Vector &y, bool reset_ODE_quantities);
void write_ODE_variables_dots(ParticlesMap *particlesMap, N_Vector &y_dot);

int compute_y_dot(realtype time, N_Vector y, N_Vector y_dot, void *data_);
double compute_EOM_hierarchical_triple(ParticlesMap *particlesMap, int inner_binary_index, int outer_binary_index, int connecting_child_in_outer_binary, bool compute_hamiltonian_only);

double compute_EOM_binary_pairs(ParticlesMap *particlesMap, int inner_binary_index, int outer_binary_index, int connecting_child_in_outer_binary, bool compute_hamiltonian_only);
double compute_EOM_binary_triplets(ParticlesMap *particlesMap, int binary_A_index, int binary_B_index, int binary_C_index, int connecting_child_in_binary_B_to_binary_A, int connecting_child_in_binary_C_to_binary_B, bool compute_hamiltonian_only);

double compute_EOM_pairwise_1PN(ParticlesMap *particlesMap, int binary_index, bool compute_hamiltonian_only);
double compute_EOM_pairwise_25PN(ParticlesMap *particlesMap, int binary_index, bool compute_hamiltonian_only);

double compute_EOM_equilibrium_tide(ParticlesMap *particlesMap, int binary_index, int star_index, int companion_index, int include_tidal_friction_terms, int include_tidal_bulges_precession_terms, int include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession);
double f_tides1_function(double e_p2, double j_p10_inv, double j_p13_inv);
double f_tides2_function(double e_p2, double j_p10_inv, double j_p13_inv);
double f_tides3_function(double e_p2, double j_p10_inv, double j_p13_inv);
double f_tides4_function(double e_p2, double j_p10_inv, double j_p13_inv);
double f_tides5_function(double e_p2, double j_p10_inv, double j_p13_inv);

double compute_EOM_equilibrium_tide_BO(ParticlesMap *particlesMap, int binary_index, int star_index, int companion_index, int include_tidal_friction_terms, int include_tidal_bulges_precession_terms, int include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession);
double compute_EOM_equilibrium_tide_BO_full(ParticlesMap *particlesMap, int binary_index, int star_index, int companion_index, int include_tidal_friction_terms, int include_tidal_bulges_precession_terms, int include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession);
double f_tides1_function_BO(double e_p2, double j_p10_inv, double j_p13_inv);
double f_tides2_function_BO(double e_p2, double j_p10_inv, double j_p13_inv);
double f_tides3_function_BO(double e_p2, double j_p10_inv, double j_p13_inv);
double f_tides4_function_BO(double e_p2, double j_p10_inv, double j_p13_inv);
double f_tides5_function_BO(double e_p2, double j_p10_inv, double j_p13_inv);

int VWXYZ_tides_function
(
    int include_tidal_friction_terms, int include_tidal_bulges_precession_terms, int include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession, \
    double t_f_inv, double k_AM, \
    double m, double M, double mu, double n, double R_div_a_p5, \
    double e, double e_p2, double j_p4_inv, double j_p10_inv,double j_p13_inv, \
    double spin_vec_dot_h_vec_unit, double spin_vec_dot_e_vec_unit,double spin_vec_dot_q_vec_unit, \
    double* V, double* W, double* X, double* Y, double* Z
);


double compute_orbital_period(Particle *particle);
double compute_AM_time_scale(Particle *P_p);
void cross_section_function(Particle *p,double *cross_section);

int root_finding_functions(realtype t, N_Vector y, realtype *root_functions, void *data_);
double roche_radius_pericenter_eggleton(double rp, double q);
double roche_radius_pericenter_sepinsky(double rp, double q, double e, double f);
