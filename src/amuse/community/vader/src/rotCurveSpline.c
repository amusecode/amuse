/* This module contains the routine rotCurveSpline, which takes a
   tabulated rotaiton curve and calculates a B-spline fit to it and
   its derivatives. */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_multifit.h>
#include "rotCurveSpline.h"

struct dpsidr_params {
  gsl_bspline_workspace *wksp;
#ifdef GSL1
  gsl_bspline_deriv_workspace *dwksp;
#endif
  gsl_vector *B, *coef;
  gsl_matrix *dB, *cov;
};

double dpsidr(double r, void *params) {
  double vphi, err;
  struct dpsidr_params *p = (struct dpsidr_params *) params;

#ifdef GSL1
  gsl_bspline_deriv_eval(log(r), 2, p->dB, p->wksp, p->dwksp);
#else
  gsl_bspline_deriv_eval(log(r), 2, p->dB, p->wksp);
#endif
  gsl_matrix_get_col(p->B, p->dB, 0);
  gsl_multifit_linear_est(p->B, p->coef, p->cov, &vphi, &err);
  return(-vphi*vphi/r);
}

void 
rotCurveSpline(const double *rTab, const double *vTab,
	       const unsigned long nTab,
	       const unsigned long bspline_degree,
	       const unsigned long bspline_breakpoints,
	       const double *r, const unsigned long nr, 
	       double *vphi, double *psi,
	       double *beta) {

  gsl_vector *xInput, *vInput, *coef, *B, *bkpt_vec;
  gsl_matrix *dB, *cov, *pred;
  gsl_bspline_workspace *wksp;
#ifdef GSL1
  gsl_bspline_deriv_workspace *dwksp;
#endif
  gsl_multifit_linear_workspace *mw;
  gsl_integration_workspace *intwksp;
  gsl_function dpsidr_func;
  struct dpsidr_params params;
  double x0, x1, v0, bkpt_sum_target, sum, x, Bij, chisq;
  double dvphidlogr, err, psiInt;
  unsigned long bspline_ncoef, i, j, ptr;

  /* Safety check: make sure we have enough data points for the choice
     of degree and breakpoints */
  bspline_ncoef = bspline_breakpoints + bspline_degree - 2;
  if (bspline_ncoef > nTab) {
    fprintf(stderr, "Error: %lu data points for %lu coefficients; ndata must be >= ncoef = nbreakpoints + degree - 2\n", nTab, bspline_ncoef);
    exit(1);
  }

  /* Initialize the GSL vectors to hold the data table */
  xInput = gsl_vector_alloc(nTab);
  vInput = gsl_vector_alloc(nTab);
  for (i=0; i<nTab; i++) {
    gsl_vector_set(xInput, i, log(rTab[i]));
    gsl_vector_set(vInput, i, vTab[i]);
  }

  /* Allocate workspace for and associated stuff for B-spline solver */
  wksp = gsl_bspline_alloc(bspline_degree, bspline_breakpoints);
#ifdef GSL1
  dwksp = gsl_bspline_deriv_alloc(bspline_degree);
#endif
  B = gsl_vector_alloc(bspline_ncoef);
  dB = gsl_matrix_alloc(bspline_ncoef, 3);

  /* Compute breakpoint locations using the method of Gans & Gill,
     Applied Spectroscopy, 38, 297, 1984. The formula for
     uniformly-spaced data is that two successive breakpoints t_i and t_i+1
     should be separated by an amount such that
     sum_{x=t_i}^{x=t_i+1} f(x)^0.5 dx = S_T
     where
     S_T = (1/(m+1)) sum f(x_i)^0.5 dx.
     Here the first sum runs over all data points between t_i and
     t_i+1, and the second runs all data points; m is the number of
     knots, and f(x) is the function to be approximated. This formula
     assumes uniform spacing dx, but here it is modified for
     non-uniform data spacing such that dx is taken to be the mean
     size of the data interval, averaging over the left and right
     sides.
  */
  v0 = gsl_vector_get(vInput, 0);
  x0 = gsl_vector_get(xInput, 0);
  x1 = gsl_vector_get(xInput, 1);
  bkpt_sum_target = sqrt(v0) * (x1-x0);
  for (i=1; i<nTab-1; i++) {
    v0 = gsl_vector_get(vInput, i);
    x0 = gsl_vector_get(xInput, i-1);
    x1 = gsl_vector_get(xInput, i+1);
    bkpt_sum_target += sqrt(v0) * 0.5 * (x1-x0);
  }
  v0 = gsl_vector_get(vInput, nTab-1);
  x0 = gsl_vector_get(xInput, nTab-2);
  x1 = gsl_vector_get(xInput, nTab-1);
  bkpt_sum_target += sqrt(v0) * (x1-x0);
  bkpt_sum_target /= (bspline_breakpoints + 1);
  bkpt_vec = gsl_vector_alloc(bspline_breakpoints);
  gsl_vector_set(bkpt_vec, 0, gsl_vector_get(xInput, 0));
  gsl_vector_set(bkpt_vec, bspline_breakpoints-1, 
		 gsl_vector_get(xInput, nTab-1));
  ptr=0;
  for (i=1; i<bspline_breakpoints-2; i++) {
    sum = 0;
    while (sum < bkpt_sum_target) {
      if (ptr == 0) {
	v0 = gsl_vector_get(vInput, 0);
	x0 = gsl_vector_get(xInput, 0);
	x1 = gsl_vector_get(xInput, 1);
	sum += sqrt(v0) * (x1-x0);
      } else if (ptr == nTab-1) {
	v0 = gsl_vector_get(vInput, nTab-1);
	x0 = gsl_vector_get(xInput, nTab-2);
	x1 = gsl_vector_get(xInput, nTab-1);
	sum += sqrt(v0) * (x1-x0);
      } else {
	v0 = gsl_vector_get(vInput, ptr);
	x0 = gsl_vector_get(xInput, ptr-1);
	x1 = gsl_vector_get(xInput, ptr+1);
	sum += sqrt(v0) * 0.5 * (x1-x0);
      }
      ptr++;
      if (ptr == nTab) break;
    }
    if (ptr == nTab) break;
    gsl_vector_set(bkpt_vec, i, gsl_vector_get(xInput, ptr-1));
  }
  gsl_vector_set(bkpt_vec, i, gsl_vector_get(xInput, nTab-1));


  /* Construct knots vector */
  gsl_bspline_knots(bkpt_vec, wksp);

  /* Allocate workspace for the least squares fit to find the B spline
     coefficients; pred is the predictor variable matrix, cov is the
     covariance matrix, coef is the vector of best-fitting
     coefficients, and mw is the workspace used by the gsl fitter */
  pred = gsl_matrix_alloc(nTab, bspline_ncoef);
  mw = gsl_multifit_linear_alloc(nTab, bspline_ncoef);
  cov = gsl_matrix_alloc(bspline_ncoef, bspline_ncoef);
  coef = gsl_vector_alloc(bspline_ncoef);

  /* Load the predictor matrix with the B spline functions evaluated
     at the input positions */
  for (i=0; i<nTab; i++) {
    x = gsl_vector_get(xInput, i); /* Position of data point */
    gsl_bspline_eval(x, B, wksp); /* Evaluate B splines at this
				     position */
    /* Load predictor matrix */
    for (j=0; j<bspline_ncoef; j++) {
      Bij = gsl_vector_get(B, j);
      gsl_matrix_set(pred, i, j, Bij);
    }
  }
  
  /* Perform least squares fit to find coefficients */
  gsl_multifit_linear(pred, vInput, coef, cov, &chisq, mw);

  /* Compute B-spline fits at input positions */
  for (i=0; i<nr; i++) {
#ifdef GSL1
    gsl_bspline_deriv_eval(log(r[i]), 2, dB, wksp, dwksp);
#else
    gsl_bspline_deriv_eval(log(r[i]), 2, dB, wksp);
#endif
    gsl_matrix_get_col(B, dB, 0);
    gsl_multifit_linear_est(B, coef, cov, vphi+i, &err);
    gsl_matrix_get_col(B, dB, 1);
    gsl_multifit_linear_est(B, coef, cov, &dvphidlogr, &err);
    beta[i] = dvphidlogr / vphi[i];
  }

  /* Compute potential at cell centers and edges by numerically
     integrating from grid edge to center */
  intwksp = gsl_integration_workspace_alloc(4096);
  dpsidr_func.function = &dpsidr;
  dpsidr_func.params = &params;
  params.wksp = wksp;
#ifdef GSL1
  params.dwksp = dwksp;
#endif
  params.B = B;
  params.coef = coef;
  params.dB = dB;
  params.cov = cov;
  psi[nr-1] = 0.0;
  for (i=nr-1; i>0; i--) {
    gsl_integration_qag(&dpsidr_func, r[i], r[i-1],
			1e-6*vphi[i]*vphi[i]*
			(r[i+1]-r[i]),
			1e-6, 4096, 6, intwksp, &psiInt, &err);
    psi[i-1] = psi[i] - psiInt;
  }

  /* Free memory */
  gsl_matrix_free(pred);
  gsl_vector_free(bkpt_vec);
  gsl_vector_free(xInput);
  gsl_vector_free(vInput);
  gsl_multifit_linear_free(mw);
  gsl_vector_free(coef);
  gsl_vector_free(B);
  gsl_matrix_free(dB);
  gsl_matrix_free(cov);
  gsl_bspline_free(wksp);
#ifdef GSL1
  gsl_bspline_deriv_free(dwksp);
#endif
  gsl_integration_workspace_free(intwksp);
}
