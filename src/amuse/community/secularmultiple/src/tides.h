#include "types.h"

bool check_for_radiative_damping(int stellar_type, double mass, double convective_envelope_mass, double convective_envelope_radius);
bool check_for_convective_damping(int stellar_type);

double from_k_AM_div_T_to_t_V(double k_AM_div_T, double apsidal_motion_constant);

double compute_t_V(Particle *star, Particle *companion, double semimajor_axis);
double compute_t_V_hurley
(
    int stellar_type,
    double mass,
    double convective_envelope_mass,
    double companion_mass,
    double semimajor_axis,
    double radius,
    double convective_envelope_radius,
    double luminosity,
    double spin_angular_frequency,
    double gyration_radius,
    double apsidal_motion_constant
);

double compute_EOM_equilibrium_tide(ParticlesMap *particlesMap, int binary_index, int star_index, int companion_index, int include_tidal_friction_terms, int include_tidal_bulges_precession_terms, int include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession);
double f_tides1_function(double e_p2, double j_p10_inv, double j_p13_inv);
double f_tides2_function(double e_p2, double j_p10_inv, double j_p13_inv);
double f_tides3_function(double e_p2, double j_p10_inv, double j_p13_inv);
double f_tides4_function(double e_p2, double j_p10_inv, double j_p13_inv);
double f_tides5_function(double e_p2, double j_p10_inv, double j_p13_inv);

double compute_EOM_equilibrium_tide_BO_full(ParticlesMap *particlesMap, int binary_index, int star_index, int companion_index, int include_tidal_friction_terms, int include_tidal_bulges_precession_terms, int include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession, int tides_method);
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
