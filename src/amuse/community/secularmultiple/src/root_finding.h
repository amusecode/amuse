#include "types.h"

double compute_AM_time_scale(Particle *P_p);
void cross_section_function(Particle *p,double *cross_section);

int root_finding_functions(realtype t, N_Vector y, realtype *root_functions, void *data_);
double roche_radius_pericenter_eggleton(double rp, double q);
double roche_radius_pericenter_sepinsky(double rp, double q, double e, double f);

int read_root_finding_data(ParticlesMap *particlesMap, int *roots_found);
int check_for_initial_roots(ParticlesMap *particlesMap);
