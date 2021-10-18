#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include "userFunc.h"

#define G        6.67384e-8    /* 2010 CODATA value in CGS units */
#define QCRIT    1.0           /* Q value at which GI turns on */
#define QDAMP    10.0          /* Rate of damping of alpha above Qcrit */
#define ALPHAMIN 1.0e-6        /* Minimum allowed alpha value */

/**************************************************************************/
/* This defines userFunc routines for the Krumholz & Burkert (2010) GI    */
/* disk problem.                                                          */
/**************************************************************************/

void
userEOS(const double t, const double dt, const grid *grd, 
	const double *col, const double *pres, const double *eInt,
	void *params,
	double *gamma, double *delta) {
  fprintf(stderr, 
	  "Warning: userEOS function called but not implemented!\n");
  return;
}

void
userAlpha(const double t, const double dt, const grid *grd, 
	  const double *col, const double *pres, const double *eInt,
	  const double *gamma, const double *delta,
	  void *params,
	  double *alpha) {

  static int nr = 0;
  static double rmin = -1.0, rmax = -1.0;
  static double *x = NULL;
  static double *h1 = NULL;
  static double *h0 = NULL;
  static double *H = NULL;
  static double *Q = NULL;
  static gsl_vector *ud_g = NULL;
  static gsl_vector *ld_g = NULL;
  static gsl_vector *diag_g = NULL;
  static gsl_vector *rhs_g = NULL;
  static gsl_vector *tau_g = NULL;
  double eta = ((double *) params)[0];
  double chi = ((double *) params)[1];
  double tQ = ((double *) params)[2];
  double u, dlnx, sL, s, sR, dsdx, beta, dbetadx;
  double mdot;
  double dlnQdlnt;
  int i;

  /* Allocate memory if necessary */
  if (grd->nr != nr) {
    nr = grd->nr;
    rmin = -1.0;
    rmax = -1.0;
    /* If previously allocated, free */
    if (x != NULL) free(x);
    if (h1 != NULL) free(h1);
    if (h0 != NULL) free(h0);
    if (H != NULL) free(H);
    if (Q != NULL) free(Q);
    if (ud_g != NULL) gsl_vector_free(ud_g);
    if (ld_g != NULL) gsl_vector_free(ld_g);
    if (diag_g != NULL) gsl_vector_free(diag_g);
    if (rhs_g != NULL) gsl_vector_free(rhs_g);
    if (tau_g != NULL) gsl_vector_free(tau_g);
    /* Allocate new memory */
    if (!(x = calloc(grd->nr, sizeof(double)))) {
      fprintf(stderr, "Error: unable to allocate memory in userAlpha!\n");
      exit(1);
    }
    if (!(h1 = calloc(grd->nr, sizeof(double)))) {
      fprintf(stderr, "Error: unable to allocate memory in userAlpha!\n");
      exit(1);
    }
    if (!(h0 = calloc(grd->nr, sizeof(double)))) {
      fprintf(stderr, "Error: unable to allocate memory in userAlpha!\n");
      exit(1);
    }
    if (!(H = calloc(grd->nr, sizeof(double)))) {
      fprintf(stderr, "Error: unable to allocate memory in userAlpha!\n");
      exit(1);
    }
    if (!(Q = calloc(grd->nr, sizeof(double)))) {
      fprintf(stderr, "Error: unable to allocate memory in userAlpha!\n");
      exit(1);
    }
    if (!(ud_g = gsl_vector_alloc(grd->nr+1))) {
      fprintf(stderr, "Error: unable to allocate memory in userAlpha!\n");
      exit(1);
    }
    if (!(ld_g = gsl_vector_alloc(grd->nr+1))) {
      fprintf(stderr, "Error: unable to allocate memory in userAlpha!\n");
      exit(1);
    }
    if (!(diag_g = gsl_vector_alloc(grd->nr+2))) {
      fprintf(stderr, "Error: unable to allocate memory in userAlpha!\n");
      exit(1);
    }
    if (!(rhs_g = gsl_vector_alloc(grd->nr+2))) {
      fprintf(stderr, "Error: unable to allocate memory in userAlpha!\n");
      exit(1);
    }
    if (!(tau_g = gsl_vector_alloc(grd->nr+2))) {
      fprintf(stderr, "Error: unable to allocate memory in userAlpha!\n");
      exit(1);
    }
  }
  /* Initailize x if necessary */
  if ((rmin != grd->r_h[0]) || (rmax != grd->r_h[grd->nr+1])) {
    for (i=0; i<grd->nr; i++) x[i] = grd->r_g[i+1]/grd->r_h[grd->nr];
    rmin = grd->r_h[0];
    rmax = grd->r_h[grd->nr+1];
  }

  /* Compute coefficients of torque equation */
  dlnx = log(grd->r_g[1]/grd->r_g[0]);
  for (i=0; i<grd->nr; i++) {

    /* Define dimensionless velocity */
    u = grd->vphi_g[i+1] / grd->vphi_h[grd->nr];

    /* Note: compute derivatives of s using one-sided differences at
       edge cells, centered-limited differences elsewhere. */
    if (i==0) {
      s = sL = sqrt(pres[0]/col[0])/grd->vphi_h[grd->nr];
      sR = sqrt(pres[1]/col[1])/grd->vphi_h[grd->nr];
      dsdx = (sR - sL) / (x[i]*dlnx);
    } else if (i==grd->nr-1) {
      sL = sqrt(pres[grd->nr-2]/col[grd->nr-2])/grd->vphi_h[grd->nr];
      s = sR = sqrt(pres[grd->nr-1]/col[grd->nr-1])/grd->vphi_h[grd->nr];
      dsdx = (sR - sL) / (x[i]*dlnx);
    } else {
      sL = sqrt(pres[i+1]/col[i+1])/grd->vphi_h[grd->nr];
      s = sqrt(pres[i]/col[i])/grd->vphi_h[grd->nr];
      sR = sqrt(pres[i-1]/col[i-1])/grd->vphi_h[grd->nr];
      dsdx = ((sL-s)*(s-sR)) > 0 ? (sL-sR)/(2.0*x[i]*dlnx) : 0;
    }

    /* Compute derivative of beta using centered difference */
    beta = grd->beta_g[i+1];
    dbetadx = (grd->beta_h[i+1] - grd->beta_h[i]) / (x[i]*dlnx);

    /* Compute torque equation coefficients */
    h1[i] = -(5.0*(beta+1)*x[i]*dsdx + 2.0*s*(beta+SQR(beta)+x[i]*dbetadx)) /
      (2.0*(beta+1)*s*x[i]);
    h0[i] = (SQR(beta)-1) * SQR(u/(x[i]*s)) / 2.0;
    Q[i] = sqrt(2.0*(1+beta))*grd->vphi_g[i+1]/grd->r_g[i+1] * 
      s*grd->vphi_h[grd->nr]/(M_PI*G*col[i]);
    dlnQdlnt = Q[i] < QCRIT ? -(exp(QCRIT/Q[i]) - exp(1.0))*u/x[i] : 0.0;
    H[i] = G*grd->r_h[grd->nr]*u*(1+beta)*col[i] *
      (2*M_PI*u*eta - 3*dlnQdlnt/tQ) / (2.0*SQR(grd->vphi_h[grd->nr])*chi);
  }

  /* Load matrix elements for tridiagonal solve */
  for (i=0; i<grd->nr; i++) {
    gsl_vector_set(ld_g, i, 1.0/SQR(x[i]*dlnx) + 
		   1.0/(2.0*SQR(x[i])*dlnx) - h1[i]/(2.0*x[i]*dlnx));
    gsl_vector_set(diag_g, i+1, -2.0/SQR(x[i]*dlnx) + h0[i]);
    gsl_vector_set(ud_g, i+1, 1.0/SQR(x[i]*dlnx) - 
		   1.0/(2.0*SQR(x[i])*dlnx) + h1[i]/(2.0*x[i]*dlnx));
    gsl_vector_set(rhs_g, i+1, H[i]);
  }

  /* Load boundary conditions into matrix */
  /* Inner BC: tau = -x in ghost cell */
  gsl_vector_set(diag_g, 0, 1);
  gsl_vector_set(ud_g, 0, 0);
  gsl_vector_set(rhs_g, 0, -grd->r_g[0]/grd->r_h[grd->nr]);
  /* Outer BC: tau' = -beta-1 */
  gsl_vector_set(ld_g, grd->nr, -1.0/dlnx);
  gsl_vector_set(diag_g, grd->nr+1, 1.0/dlnx);
  gsl_vector_set(rhs_g, grd->nr+1, -grd->beta_g[grd->nr+1]-1);

  /* Solve matrix equation to get dimensionless torque */
  gsl_linalg_solve_tridiag(diag_g, ud_g, ld_g, rhs_g, tau_g);

  /* Use torque to compute alpha */
  mdot = chi*SQR(grd->vphi_h[grd->nr])*grd->vphi_h[grd->nr] / G;
  for (i=0; i<grd->nr; i++) {
    alpha[i] = -gsl_vector_get(tau_g, i+1)*mdot*grd->vphi_h[grd->nr] *
      grd->r_h[grd->nr] /
      (2.0*M_PI*SQR(grd->r_g[i+1])*pres[i]);
    if (Q[i] > QCRIT) alpha[i] *= exp(-QDAMP*(Q[i]-QCRIT));
    if (alpha[i] < 0.0) alpha[i] = ALPHAMIN;
  }
}

void
userMassSrc(const double t, const double dt, const grid *grd,
	    const double *col, const double *pres, const double *eInt,
	    const double *gamma, const double *delta,
	    void *params,
	    double *massSrc) {
  fprintf(stderr, 
	  "Warning: userMassSrc function called but not implemented!\n");
  return;
}

void
userIntEnSrc(const double t, const double dt, const grid *grd,
	     const double *col, const double *pres, const double *eInt,
	     const double *gamma, const double *delta,
	     void *params, 
	     double *intEnSrc) {
  /* Cooling rate = eta Sigma sigma^2 Omega = eta P vphi/r */
  int i;
  double eta = ((double *) params)[0];

  for (i=0; i<grd->nr; i++)
    intEnSrc[i] = -eta * pres[i] * grd->vphi_g[i+1] / grd->r_g[i+1];
}

void
userIBC(const double t, const double dt, const grid *grd,
	const double *col, const double *pres, const double *eInt,
	const double *gamma, const double *delta,
	const pres_bc_type ibc_pres, const enth_bc_type ibc_enth,
	void *params, 
	double *ibc_pres_val, double *ibc_enth_val) {

  /* Set torque. If Q <= QCRIT in innermost cell, so that disk is
     gravitationally unstable at the center, dimensionless torque =
     -x, dimensional torque = -(r_in/R) Mdot vphi/R. If Q > QCRIT in
     innermost cell, scale down torque in the same way as in
     userAlpha. */
  double chi = ((double *) params)[1];
  double x, mdot, s, Q;

  s = sqrt(pres[0]/col[0])/grd->vphi_h[grd->nr];
  Q = sqrt(2.0*(1+grd->beta_g[1]))*grd->vphi_g[1]/grd->r_g[1] * 
    s*grd->vphi_h[grd->nr]/(M_PI*G*col[0]);
  mdot = chi*SQR(grd->vphi_h[grd->nr])*grd->vphi_h[grd->nr] / G;
  x = grd->r_g[0] / grd->r_h[grd->nr];
  *ibc_pres_val = -x * mdot * grd->vphi_g[0] * grd->r_h[grd->nr];
  if (Q > QCRIT) *ibc_pres_val *= exp(-QDAMP*(Q-QCRIT));
  *ibc_enth_val = gamma[0]/(gamma[0]-1.0)*pres[0]/col[0];
}

void
userOBC(const double t, const double dt, const grid *grd,
	const double *col, const double *pres, const double *eInt,
	const double *gamma, const double *delta,
	const pres_bc_type obc_pres, const enth_bc_type obc_enth,
	void *params, 
	double *obc_pres_val, double *obc_enth_val) {
  fprintf(stderr, 
	  "Warning: userOBC function called but not implemented!\n");
  return;
}


void
userPreTimestep(const double t, const double dt,
		const grid *grd, double *col, double *pres,
		double *eInt, double *mBnd, double *eBnd,
		double *mSrc, double *eSrc,
		void *params, const unsigned long nUserOut,
		double *userOut) {
  fprintf(stderr,
	  "Warning: userPreTimestep function called but not implemented!\n");
  return;
}

void
userPostTimestep(const double t, const double dt,
		 const grid *grd, double *col, double *pres,
		 double *eInt, double *mBnd, double *eBnd,
		 double *mSrc, double *eSrc,
		 void *params, const unsigned long nUserOut,
		 double *userOut) {
  fprintf(stderr,
	  "Warning: userPostTimestep function called but not implemented!\n");
  return;
}

void
userCheckRead(
	      FILE *fp, grid *grd, const unsigned long nOut,
	      double *tOut, double *colOut,
	      double *presOut, double *eIntOut, double *mBndOut,
	      double *eBndOut, double *mSrcOut, double *eSrcOut,
	      const unsigned long nUserOut, double *userOut,
	      void *params
	      ) {
  fprintf(stderr,
	  "Warning: userCheckRead function called but not implemented!\n");
  return;
}

void
userCheckWrite(
	      FILE *fp,
	      const grid *grd, const unsigned long nOut,
	      const double *tOut, const double *colOut,
	      const double *presOut, const double *eIntOut,
	      const double *mBndOut, const double *eBndOut,
	      const double *mSrcOut, const double *eSrcOut,
	      const unsigned long nUserOut, const double *userOut,
	      const void *params
	      ) {
  fprintf(stderr,
	  "Warning: userCheckWrite function called but not implemented!\n");
  return;
}
