#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "gpu.h"

#define NCELLS3 (NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2)
#define NBUFF (std::max(NCELLX, std::max(NCELLY, NCELLZ)))

static void PackBoundary(int i, int j, int i1, int i2, int i3,
			 int fixed,
			 const double* cpu_e, const double* cpu_f,
			 double *boundary_values) {
  int bound_idx = i + NBUFF*j;

  int pos[3];
  pos[i1] = fixed;
  pos[i2] = i;
  pos[i3] = j;
  int idx = ((pos[0]+NBOUND) +
	     (pos[1]+NBOUND) * (NCELLX + 2*NBOUND) +
	     (pos[2]+NBOUND) * (NCELLX + 2*NBOUND) * (NCELLY + 2*NBOUND));

  boundary_values[bound_idx + 0*NBUFF*NBUFF] = cpu_e[idx];
  boundary_values[bound_idx + 1*NBUFF*NBUFF] = cpu_f[idx + 0*NCELLS3];
  boundary_values[bound_idx + 2*NBUFF*NBUFF] = cpu_f[idx + 1*NCELLS3];
  boundary_values[bound_idx + 3*NBUFF*NBUFF] = cpu_f[idx + 2*NCELLS3];
}

static void UnpackBoundary(int i, int j, int i1, int i2, int i3,
			   int fixed,
			   double* cpu_e, double* cpu_f,
			   const double *boundary_values) {
  int bound_idx = i + NBUFF*j;

  int pos[3];
  pos[i1] = fixed;
  pos[i2] = i;
  pos[i3] = j;
  int idx = ((pos[0]+NBOUND) +
	     (pos[1]+NBOUND) * (NCELLX + 2*NBOUND) +
	     (pos[2]+NBOUND) * (NCELLX + 2*NBOUND) * (NCELLY + 2*NBOUND));

  cpu_e[idx] = boundary_values[bound_idx + 0*NBUFF*NBUFF];
  cpu_f[idx + 0*NCELLS3] = boundary_values[bound_idx + 1*NBUFF*NBUFF];
  cpu_f[idx + 1*NCELLS3] = boundary_values[bound_idx + 2*NBUFF*NBUFF];
  cpu_f[idx + 2*NCELLS3] = boundary_values[bound_idx + 3*NBUFF*NBUFF];
}

extern "C" void cpu_pack_boundary_values_(
    const double *cpu_e, const double *cpu_f,
    double *boundary_values) {

  int stride = 4*NBUFF*NBUFF;

  for (int j=0; j<NCELLY; j++) {
    for (int k=0; k<NCELLZ; k++) {
      PackBoundary(j, k, 0, 1, 2, 0, cpu_e, cpu_f,
		   &boundary_values[0*stride]);
      PackBoundary(j, k, 0, 1, 2, NCELLX-1, cpu_e, cpu_f,
		   &boundary_values[1*stride]);
    }
  }

  for (int i=0; i<NCELLX; i++) {
    for (int k=0; k<NCELLZ; k++) {
      PackBoundary(i, k, 1, 0, 2, 0, cpu_e, cpu_f,
		   &boundary_values[2*stride]);
      PackBoundary(i, k, 1, 0, 2, NCELLY-1, cpu_e, cpu_f,
		   &boundary_values[3*stride]);
    }
  }

  for (int i=0; i<NCELLX; i++) {
    for (int j=0; j<NCELLY; j++) {
      PackBoundary(i, j, 2, 0, 1, 0, cpu_e, cpu_f,
		   &boundary_values[4*stride]);
      PackBoundary(i, j, 2, 0, 1, NCELLZ-1, cpu_e, cpu_f,
		   &boundary_values[5*stride]);
    }
  }
}

extern "C" void cpu_unpack_boundary_values_(
    const double *boundary_values,
    double *cpu_e, double *cpu_f) {

  int stride = 4*NBUFF*NBUFF;

  for (int j=0; j<NCELLY; j++) {
    for (int k=0; k<NCELLZ; k++) {
      UnpackBoundary(j, k, 0, 1, 2, -1, cpu_e, cpu_f,
		     &boundary_values[0*stride]);
      UnpackBoundary(j, k, 0, 1, 2, NCELLX, cpu_e, cpu_f,
		     &boundary_values[1*stride]);
    }
  }

  for (int i=0; i<NCELLX; i++) {
    for (int k=0; k<NCELLZ; k++) {
      UnpackBoundary(i, k, 1, 0, 2, -1, cpu_e, cpu_f,
		     &boundary_values[2*stride]);
      UnpackBoundary(i, k, 1, 0, 2, NCELLY, cpu_e, cpu_f,
		     &boundary_values[3*stride]);
    }
  }

  for (int i=0; i<NCELLX; i++) {
    for (int j=0; j<NCELLY; j++) {
      UnpackBoundary(i, j, 2, 0, 1, -1, cpu_e, cpu_f,
		     &boundary_values[4*stride]);
      UnpackBoundary(i, j, 2, 0, 1, NCELLZ, cpu_e, cpu_f,
		     &boundary_values[5*stride]);
    }
  }
}
