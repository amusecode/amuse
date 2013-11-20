#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "gpu.h"
#include "aton_cpp.h"
#include "aton_fortran.h"

#define NCELL4 ((NCELLZ+NBOUND2)*(NCELLY+NBOUND2)*(NCELLX+NBOUND2))

// FIXME: Replace this with aton_cell_index.
static int CellIndex(int i, int j, int k) {
  return ((i + NBOUND) +
	  (j + NBOUND)*(NCELLX + NBOUND2) +
	  (k + NBOUND)*(NCELLX + NBOUND2)*(NCELLY + NBOUND2));
}


extern "C" void aton_validate_(
    const int *label, const double *c,
    const double *cpu_e, const double* cpu_d, const double* cpu_t,
    const double* cpu_x, const double *cpu_f) {
  for (int i = 0; i < NCELLX; i++) {
    for (int j = 0; j < NCELLY; j++) {
      for (int k = 0; k < NCELLZ; k++) {
	int b = CellIndex(i, j, k);

	bool bad = false;

	if (cpu_d[b] == 0.0) {
	  printf("(%d) zero cpu_d\n", *label);
	  bad = true;
	}

	if (isnan(cpu_e[b])) {
	  printf("(%d) nan cpu_e\n", *label);
	  bad = true;
	}
	if (isinf(cpu_e[b])) {
	  printf("(%d) inf cpu_e\n", *label);
	  bad = true;
	}
	if (cpu_e[b] <= 0.0) {
	  printf("(%d) cpu_e <= 0: %e\n", *label, cpu_e[b]);
	  bad = true;
	}

	if (isnan(cpu_t[b])) {
	  printf("(%d) nan cpu_t\n", *label);
	  bad = true;
	}
	if (isinf(cpu_t[b])) {
	  printf("(%d) inf cpu_t\n", *label);
	  bad = true;
	}
	if (cpu_t[b] <= 0.0) {
	  printf("(%d) cpu_t <= 0: %e\n", *label, cpu_t[b]);
	  bad = true;
	}

	if (isnan(cpu_x[b])) {
	  printf("(%d) nan cpu_x\n", *label);
	  bad = true;
	}
	if (isinf(cpu_x[b])) {
	  printf("(%d) inf cpu_x\n", *label);
	  bad = true;
	}
	if (cpu_x[b] < 0.0) {
	  printf("(%d) cpu_x < 0: %e\n", *label, cpu_x[b]);
	  bad = true;
	}
	if (cpu_x[b] > 1.0) {
	  printf("(%d) cpu_x > 1: %e\n", *label, cpu_x[b]);
	  bad = true;
	}

	if (isnan(cpu_f[b+0*NCELL4])) {
	  printf("(%d) nan cpu_f_x\n", *label);
	  bad = true;
	}
	if (isnan(cpu_f[b+1*NCELL4])) {
	  printf("(%d) nan cpu_f_y\n", *label);
	  bad = true;
	}
	if (isnan(cpu_f[b+2*NCELL4])) {
	  printf("(%d) nan cpu_f_z\n", *label);
	  bad = true;
	}
	if (isinf(cpu_f[b+0*NCELL4])) {
	  printf("(%d) inf cpu_f_x\n", *label);
	  bad = true;
	}
	if (isinf(cpu_f[b+1*NCELL4])) {
	  printf("(%d) inf cpu_f_y\n", *label);
	  bad = true;
	}
	if (isinf(cpu_f[b+2*NCELL4])) {
	  printf("(%d) inf cpu_f_z\n", *label);
	  bad = true;
	}

	double F = sqrtf(cpu_f[b+0*NCELL4]*cpu_f[b+0*NCELL4] +
			 cpu_f[b+1*NCELL4]*cpu_f[b+1*NCELL4] +
			 cpu_f[b+2*NCELL4]*cpu_f[b+2*NCELL4]);
	if (F > (*c)*cpu_e[b]) {
	  printf("(%d) F > cN. %e > %e * %e\n",
		 *label, F, *c, cpu_e[b]);
	  bad = true;
	}

	if (bad) {
	  printf("Error (%d): Assertion failed at %d,%d,%d, idx=%d\n",
		 *label, i, j, k, b);
	  return;
	}

      }
    }
  }
}

namespace aton {
  
bool validate(State state, double c_light) {
  int label = 0;
  aton_validate_(&label, &c_light, state.E, state.nH, state.T, state.xHII, state.F);
  return true;  // FIXME
}

}
