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

#include "structure.h"
#include "ODE_system.h"
#include "root_finding.h"
#include "newtonian.h"
#include "postnewtonian.h"
#include "tides.h"
#include "external.h"


int evolve(ParticlesMap *particlesMap, External_ParticlesMap *external_particlesMap, double start_time, double time_step, double *output_time, double *hamiltonian, int *output_flag, int *error_code);
static int check_flag(void *flagvalue, char *funcname, int opt);
