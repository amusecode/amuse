#include "userFunc.h"

/**************************************************************************/
/* This defines userFunc routines that do nothing. Use only if you plan   */
/* to run the code with fixed values for boundary conditions, alpha, and  */
/* heating and cooling rates.                                             */
/**************************************************************************/

void
userAlpha(const double t, const double dt, const grid *grd, 
	  const double *col, const double *pres, const double *eInt,
	  const double *gamma, const double *delta,
	  void *params,
	  double *alpha) {
  fprintf(stderr, 
	  "Warning: userAlpha function called but not implemented!\n");
  return;
}

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
