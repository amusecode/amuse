#include "userFunc.h"
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

/**************************************************************************/
/* This defines userFunc routines for the Pringle (1981) ring test, with  */
/* a complex equation of state that includes both gas and radiation       */
/* presssure.                                                             */
/**************************************************************************/

void
userEOS(const double t, const double dt, const grid *grd, 
	const double *col, const double *pres, const double *eInt,
	void *params,
	double *gamma, double *delta) {
  int i;
  double gammaGas = 5.0/3.0;   /* Gas gamma value */

  /* Fill gamma and delta values */
  for (i=0; i<grd->nr; i++) {
    gamma[i] = ((16-3*gammaGas)*pres[i] + (16-15*gammaGas)*eInt[i]) /
      (9*pres[i] + (13-12*gammaGas)*eInt[i]);
    delta[i] = 4*(3*pres[i]-eInt[i])*((gammaGas-1)*eInt[i]-pres[i]) /
      (pres[i]*(3*(gammaGas-1)*eInt[i]+(3*gammaGas-7)*pres[i]));
  }
}

void
userAlpha(const double t, const double dt, const grid *grd, 
	  const double *col, const double *pres, const double *eInt,
	  const double *gamma, const double *delta,
	  void *params,
	  double *alpha) {
  /* alpha = nu col/pres vphi / r */

  int i;
  double nu = *((double *) params);

  for (i=0; i<grd->nr; i++)
    alpha[i] = nu * col[i]/pres[i] * grd->vphi_g[i+1] / grd->r_g[i+1];

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
  fprintf(stderr, 
	  "Warning: userIntEnSrc function called but not implemented!\n");
  return;
}

void
userIBC(const double t, const double dt, const grid *grd,
	const double *col, const double *pres, const double *eInt,
	const double *gamma, const double *delta,
	const pres_bc_type ibc_pres, const enth_bc_type ibc_enth,
	void *params, 
	double *ibc_pres_val, double *ibc_enth_val) {
  double nu = ((double *) params)[0];
  double r0 = ((double *) params)[1];
  double m0 = ((double *) params)[2];
  double colMin = ((double *) params)[3];
  double pOverCol = ((double *) params)[4];
  double x = grd->r_g[0]/r0;
  double tau = 12.0*nu*t/SQR(r0) + SMALL;
  double sigma0 = m0/(M_PI*SQR(r0));
  double sigma;

  sigma = sigma0/(pow(x, 0.25)*tau) * exp(-SQR(x-1.0)/tau) *
    gsl_sf_bessel_Inu_scaled(0.25, 2*x/tau);
  if (!(isfinite(sigma))) sigma=0.0;
  sigma = (sigma < colMin) ? colMin : sigma;
  *ibc_pres_val = -3.0*M_PI*grd->r_g[0]*nu*grd->vphi_g[0]*sigma;
  *ibc_enth_val = gamma[0]/(gamma[0]-1)*pOverCol;
}

void
userOBC(const double t, const double dt, const grid *grd,
	const double *col, const double *pres, const double *eInt,
	const double *gamma, const double *delta,
	const pres_bc_type obc_pres, const enth_bc_type obc_enth,
	void *params, 
	double *obc_pres_val, double *obc_enth_val) {
  double nu = ((double *) params)[0];
  double r0 = ((double *) params)[1];
  double m0 = ((double *) params)[2];
  double colMin = ((double *) params)[3];
  double pOverCol = ((double *) params)[4];
  double x = grd->r_g[grd->nr+1]/r0;
  double tau = 12.0*nu*t/SQR(r0) + SMALL;
  double sigma0 = m0/(M_PI*SQR(r0));
  double sigma;

  sigma = sigma0/(pow(x, 0.25)*tau) * exp(-SQR(x-1.0)/tau) *
    gsl_sf_bessel_Inu_scaled(0.25, 2*x/tau);
  if (!(isfinite(sigma))) sigma=0.0;
  sigma = (sigma < colMin) ? colMin : sigma;
  *obc_pres_val = -3.0*M_PI*grd->r_g[grd->nr+1]*nu*
    grd->vphi_g[grd->nr+1]*sigma;
  *obc_enth_val = gamma[0]/(gamma[0]-1)*pOverCol;
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
