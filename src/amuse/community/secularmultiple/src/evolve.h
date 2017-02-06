#include <math.h>
#include <cstdlib>
#include <map>

#include "types.h"
#include "../interface.h"

#include "cvode/cvode.h"					/* prototypes for CVODE fcts., consts. */
#include "cvode/nvector_serial.h"			/* serial N_Vector types, fcts., macros */
#include "cvode/cvode_dense.h"				/* prototype for CVDense */
#include "cvode/sundials_dense.h"			/* definitions DlsMat DENSE_ELEM */
#include "cvode/sundials_types.h"			/* definition of type realtype */

#include "ODE_system.h"

void set_initial_ODE_variables(ParticlesMap *particlesMap, N_Vector &y, N_Vector &y_abs_tol,double abs_tol_spin_vec, double abs_tol_e_vec, double abs_tol_h_vec);
void extract_final_ODE_variables(ParticlesMap *particlesMap, N_Vector &y_out);

int read_root_finding_data(ParticlesMap *particlesMap, int *roots_found);
int check_for_initial_roots(ParticlesMap *particlesMap);

int evolve(ParticlesMap *particlesMap, double start_time, double time_step, double *output_time, double *hamiltonian, int *output_flag, int *error_code);
static int check_flag(void *flagvalue, char *funcname, int opt);
