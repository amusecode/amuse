#include <math.h>
#include "userFunc.h"

/**************************************************************************/
/* This defines userFunc routines for the Lynden-Bell & Pringle (1974,    */
/* MNRAS, 168, 603) self-similar disk evolution test                      */
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
  /* alpha = (nu0 vphi/R0) col/pres */

  int i;
  double nu0 = ((double *) params)[0];
  double R0 = ((double *) params)[1];

  for (i=0; i<grd->nr; i++)
    alpha[i] = nu0 * grd->vphi_g[i+1]/R0 * col[i]/pres[i];
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

  double nu0 = ((double *) params)[0];
  double R0 = ((double *) params)[1];
  double Mdot0 = ((double *) params)[2];
  double pOverCol = ((double *) params)[3];
  double ts = R0*R0/(3.0*nu0);
  double x = grd->r_g[0]/R0;
  double T = t/ts;

  *ibc_pres_val = -Mdot0*grd->vphi_g[0]*R0*x*exp(-x/T)/pow(T,1.5);
  *ibc_enth_val = gamma[0]/(gamma[0]-1)*pOverCol;
}

void
userOBC(const double t, const double dt, const grid *grd,
	const double *col, const double *pres, const double *eInt,
	const double *gamma, const double *delta,
	const pres_bc_type obc_pres, const enth_bc_type obc_enth,
	void *params, 
	double *obc_pres_val, double *obc_enth_val) {

  double nu0 = ((double *) params)[0];
  double R0 = ((double *) params)[1];
  double Mdot0 = ((double *) params)[2];
  double pOverCol = ((double *) params)[3];
  double ts = R0*R0/(3.0*nu0);
  double x = grd->r_g[grd->nr+1]/R0;
  double T = t/ts;

  *obc_pres_val = -Mdot0*grd->vphi_g[grd->nr+1]*R0*x*exp(-x/T)/pow(T,1.5);
  *obc_enth_val = gamma[grd->nr-1]/(gamma[grd->nr-1]-1)*pOverCol;
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

