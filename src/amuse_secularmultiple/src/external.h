#include "types.h"

void apply_user_specified_instantaneous_perturbation(ParticlesMap *particlesMap);
void set_positions_and_velocities(ParticlesMap *particlesMap);
void update_masses_positions_and_velocities_of_all_bodies(ParticlesMap *particlesMap);
void update_masses_positions_and_velocities_of_all_binaries(ParticlesMap *particlesMap);
void update_orbital_vectors_in_binaries_from_positions_and_velocities(ParticlesMap *particlesMap);



void update_position_vectors_external_particles(ParticlesMap *particlesMap, External_ParticlesMap *external_particlesMap, double time);
void compute_position_vectors_external_particles(ParticlesMap *particlesMap, External_Particle *perturber, double time, double *r_per, double r_per_vec[3]);

double compute_EOM_binary_pairs_external_perturbation(ParticlesMap *particlesMap, External_ParticlesMap *external_particlesMap, int binary_index, int perturber_index, double time, bool compute_hamiltonian_only);
void compute_true_anomaly_from_mean_anomaly_hyperbolic_orbit(double mean_anomaly, double eccentricity,double *cos_true_anomaly,double *sin_true_anomaly);

int apply_external_perturbation_assuming_integrated_orbits(ParticlesMap *particlesMap, External_ParticlesMap *external_particlesMap);
void apply_external_perturbation_assuming_integrated_orbits_binary_pair(ParticlesMap *particlesMap, External_ParticlesMap *external_particlesMap, int binary_index, int perturber_index);

double retrieve_D_function(int function_index, double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef);
double retrieve_D_function_e_derivative(int function_index, double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef);
