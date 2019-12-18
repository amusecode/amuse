#include "types.h"

double compute_orbital_period(Particle *particle);

double compute_EOM_binary_pairs(ParticlesMap *particlesMap, int inner_binary_index, int outer_binary_index, int connecting_child_in_outer_binary, bool compute_hamiltonian_only);
double compute_EOM_binary_triplets(ParticlesMap *particlesMap, int binary_A_index, int binary_B_index, int binary_C_index, int connecting_child_in_binary_B_to_binary_A, int connecting_child_in_binary_C_to_binary_B, bool compute_hamiltonian_only);
