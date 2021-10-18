#include <math.h>
#include "ppmExtrap.h"

void ppmExtrap(const int nx, const double *dx, double *v, double *vL, 
	       double *vR, double *dv) {
  int i;

  /* C & W eqn. 1.7 - 1.8 */
  for (i=1; i<nx-2; i++) {

    dv[i] = dx[i] / (dx[i-1]+dx[i]+dx[i+1]) * 
      ( (2*dx[i-1]+dx[i]) * (v[i+i]-v[i]) / (dx[i+1]+dx[i]) +
	(dx[i]+2*dx[i+1]) * (v[i]-v[i-1]) / (dx[i-1]+dx[i]) );

    /* C & W eqn. 1.8 */
    if ((v[i+1]-v[i])*(v[i]-v[i-1]) > 0) {
      dv[i] = copysign(fmin(fabs(dv[i]), 2*fabs(v[i]-v[i-1])), dv[i]);
    } else {
      dv[i] = 0.0;
    }
  }

  /* C & W eqn. 1.6 */
  for (i=1; i<nx-2; i++) {
    vR[i] = v[i] + dx[i]*(v[i+1]-v[i])/(dx[i]+dx[i+1]) +
      ( (2*dx[i+1]*dx[i])/(dx[i]+dx[i+1]) * 
	( (dx[i-1]+dx[i])/(2*dx[i]+dx[i+1]) - 
	  (dx[i+2]+dx[i+1])/(2*dx[i+1]+dx[i]) ) * (v[i+1]-v[i])
	- dx[i]*(dx[i-1]+dx[i])*dv[i+1]/(2*dx[i]+dx[i+1])
	+ dx[i+1]*(dx[i+1]+dx[i+2])*dv[i]/(dx[i]+2*dx[i+1]) ) /
      (dx[i-1]+dx[i]+dx[i+1]+dx[i+2]);
    vL[i+1] = vR[i];
  }

  /* Assign L and R values in special cases (C & W 1.10) */
  for (i=2; i<nx-2; i++) {
    if ((vR[i]-v[i])*(v[i]-vL[i]) < 0) {
      /* Use flat interpolation if we're at a local min or max */
      vL[i] = vR[i] = v[i];
    } else if ((vR[i]-vL[i])*(v[i]-0.5*(vL[i]+vR[i])) >
	       (vR[i]-vL[i])*(vR[i]-vL[i])/6.0) {
      /* Enforce monotonicity on the left */
      vL[i] = 3*v[i]-2*vR[i];
    } else if (-(vR[i]-vL[i])*(vR[i]-vL[i])/6.0 >
	       (vR[i]-vL[i])*(v[i]-0.5*(vR[i]+vL[i]))) {
      /* Enforce monotonicity on the right */
      vR[i] = 3*v[i]-2*vL[i];
    }
  }

  /* Handle edge cases. */
  vL[1] = vR[0] = 0.5*(v[0]+v[1]);
  vL[nx-1] = vR[nx-2] = 0.5*(v[nx-2]+v[nx-1]);
  vL[0] = v[0];
  vR[nx-1] = v[nx-1];
}
