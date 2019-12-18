#include "types.h"

void extract_ODE_variables(ParticlesMap *particlesMap, N_Vector &y, double delta_time, bool reset_ODE_quantities);
void write_ODE_variables_dots(ParticlesMap *particlesMap, N_Vector &y_dot);

int compute_y_dot(realtype time, N_Vector y, N_Vector y_dot, void *data_);

void set_initial_ODE_variables(ParticlesMap *particlesMap, N_Vector &y, N_Vector &y_abs_tol,double abs_tol_spin_vec, double abs_tol_e_vec, double abs_tol_h_vec);
void extract_final_ODE_variables(ParticlesMap *particlesMap, N_Vector &y_out);
