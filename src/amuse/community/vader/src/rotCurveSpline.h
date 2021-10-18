#ifndef _rotCurveSpline_h_
#define _rotCurveSpline_h_

/* Routine to use B splines to fit tabulated rotation curves */

void 
rotCurveSpline(const double *rTab, const double *vTab, 
	       const unsigned long nTab,
	       const unsigned long bspline_degree, 
	       const unsigned long bspline_breakpoints,
	       const double *r, const unsigned long nr, 
	       double *vphi, double *psi,
	       double *beta);
#endif
