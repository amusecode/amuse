#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_poly.h>
#include "userFunc.h"

#define G        6.67384e-8    /* 2010 CODATA value in CGS units */
#define ALPHAMIN 1.0e-3        /* Minimum alpha allowed for numerical reasons */

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

  gsl_poly_complex_workspace *wksp;
  double coef[6], sol[10];
  double omega, kappa, T1, sigma, kcrit, J, Q, d, tfast, tgrowth;
  int i, j;
  const int m=2;
  gsl_complex nu2, nu;
  double alphacoef = ((double *) params)[4];

  /* Allocate workspace */
  wksp = gsl_poly_complex_workspace_alloc(6);

  /* Store the coefficiencts of the dispersion relation polynomial
     that don't change */
  coef[2] = -8.0;
  coef[3] = coef[4] = 0.0;

  /* Loop over cells */
  for (i=0; i<grd->nr; i++) {

    /* Handle negative column densities, which may arise if the time
       step overshoots. We just need to make sure the code doesn't
       barf here, as the negative column densities will be fixed by
       the iterative solver. */
    if (col[i] <= 0.0) {
      alpha[i] = ALPHAMIN;
      continue;
    }

    /* The various quantities needed for the stability analysis */
    sigma = sqrt(pres[i]/col[i]);
    omega = grd->vphi_g[i+1]/grd->r_g[i+1];
    kappa = sqrt(2.0*(1.0+grd->beta_g[i+1]))*omega;
    T1 = -(2.0*m*omega/(kappa*grd->r_g[i+1])) * 
      (2.0*m*omega/(kappa*grd->r_g[i+1])) * (grd->beta_g[i+1]-1);
    kcrit = kappa*kappa / (2.0*M_PI*G*col[i]);
    Q = kappa*sigma / (M_PI*G*col[i]);
    if (Q < 0)
      printf("Q < 0! kappa = %f, sigma = %f, col = %f\n", kappa, sigma, col[i]);
    J = sqrt(T1)/kcrit;
   
    /* Coefficients of the disperison relation polynomial */
    coef[0] = -Q*Q*Q*Q;
    coef[1] = 6.0*Q*Q;
    coef[5] = 16.0*J*J;

    /* Find roots, giving minima and maxima of D */
    gsl_poly_complex_solve(coef, 6, wksp, sol);

    /* For each root, compute the growth time of the instability */
    tfast = 1.0e30;
    for (j=0; j<5; j++) {

      /* Skip roots with negative real part, or imaginary part that is
	 greater than roundoff */
      if ((sol[2*j] < 0) || (fabs(sol[2*j+1]) > 1.0e-6)) continue;

      /* Compute D */
      d = (Q*Q/(4.0*sol[2*j]*sol[2*j]) - 1.0/sol[2*j]) *
	(Q*Q/(4.0*sol[2*j]*sol[2*j]) - 1.0/sol[2*j] - 
	 4.0*J*J*sol[2*j]*sol[2*j]);

      /* Ignore stable roots */
      if (d > 0) continue;

      /* Get frequency and growth time */
      GSL_SET_COMPLEX(&nu2, 
		      1.0 + 0.5*(Q*Q/(4.0*sol[2*j]*sol[2*j]) - 1.0/sol[2*j]),
		      0.5 * sqrt(-d));
      nu = gsl_complex_sqrt(nu2);
      tgrowth = 1.0 / (GSL_IMAG(nu)*kappa) / (2.0*M_PI/omega);

      /* Store minimum growth time */
      if (tgrowth < tfast) tfast = tgrowth;
    }

    /* Now check for gravitational instability */
    if (Q < 1) {
      GSL_SET_COMPLEX(&nu, 0.0, sqrt(1.0/(Q*Q)-1.0));
      tgrowth = 1.0 / (GSL_IMAG(nu)*kappa) / (2.0*M_PI/omega);
      if (tgrowth < tfast) tfast = tgrowth;
    }

    /* Compute alpha based on the fastest growing mode timescale */
    alpha[i] = alphacoef * exp(-tfast+1.0);
    if (alpha[i] > 1.0) alpha[i] = 1.0;
    if (alpha[i] < ALPHAMIN) alpha[i] = ALPHAMIN;
  }

  /* Free workspace */
  gsl_poly_complex_workspace_free(wksp);
}

void
userMassSrc(const double t, const double dt, const grid *grd,
	    const double *col, const double *pres, const double *eInt,
	    const double *gamma, const double *delta,
	    void *params,
	    double *massSrc) {
  /* Estimate Mdot by 
     dSigma/dt = - eta_ML Sigma_SFR 
               = - eta_ML eps_ff Sigma_g/t_ff, 
  */

  int i;
  double shapefac = ((double *) params)[2];
  double zetad = ((double *) params)[3];
  double etaML = ((double *) params)[5];
  double epsffmax = ((double *) params)[6];
  double rhostar, a, b, c, h;
  double rhog, tff, alphavir, epsff;

  for (i=0; i<grd->nr; i++) {
    /* Stellar density */
    rhostar = shapefac * SQR(grd->vphi_g[i+1]) * 
      (1.0+2.0*grd->beta_g[i+1]) / 
      (4.0*M_PI*G*SQR(grd->r_g[i+1]));
    /* Get scale height from OML model */
    if (rhostar > 0.0) {
      a = 2.0*M_PI*zetad*G*rhostar*col[i];
      b = M_PI/2.0*G*SQR(col[i]);
      c = -pres[i];
      h = (-b + sqrt(b*b-4.0*a*c))/(2.0*a);
    } else {
      h = pres[i] / (M_PI/2.0*G*SQR(col[i]));
    }
    /* Gas density and free-fall time */
    rhog = col[i]/(2*h);
    tff = sqrt(3*M_PI/(32*G*rhog));
    /* Virial ratio */
    alphavir = pres[i]/(M_PI/2.0*G*col[i]*col[i]*h);
    /* epsff */
    if (alphavir <= 1.0) epsff = epsffmax;
    else epsff = epsffmax*exp(-alphavir);
    /* Mdot */
    massSrc[i] = -(1.0+etaML) * epsff * col[i] / tff;
  }
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
  double sigmath = ((double *) params)[1];
  double shapefac = ((double *) params)[2];
  double zetad = ((double *) params)[3];
  double sigma2, sigmaNT, rhostar, a, b, c, h;

  for (i=0; i<grd->nr; i++) {
    /* Gas velocity dispersion */
    sigma2 = pres[i]/col[i];
    if (sigma2 > SQR(sigmath)) {
      sigmaNT = sqrt(sigma2 - SQR(sigmath));
      /* Stellar density */
      rhostar = shapefac * SQR(grd->vphi_g[i+1]) * 
	(1.0+2.0*grd->beta_g[i+1]) / 
	(4.0*M_PI*G*SQR(grd->r_g[i+1]));
      /* Get scale height from OML model */
      if (rhostar > 0.0) {
	a = 2.0*M_PI*zetad*G*rhostar*col[i];
	b = M_PI/2.0*G*SQR(col[i]);
	c = -pres[i];
	h = (-b + sqrt(b*b-4.0*a*c))/(2.0*a);
      } else {
	h = pres[i] / (M_PI/2.0*G*SQR(col[i]));
      }
      /* Cooling rate = eta Sigma sigmaNT^2 / (h/sigmaNT) */
      intEnSrc[i] = -eta * col[i] * SQR(sigmaNT) / (h/sigmaNT);
    } else {
      intEnSrc[i] = 0.0;
    }
  }
}

void
userIBC(const double t, const double dt, const grid *grd,
	const double *col, const double *pres, const double *eInt,
	const double *gamma, const double *delta,
	const pres_bc_type ibc_pres, const enth_bc_type ibc_enth,
	void *params, 
	double *ibc_pres_val, double *ibc_enth_val) {
  fprintf(stderr, 
	  "Warning: userIBC function called but not implemented!\n");
  return;
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
