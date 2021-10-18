#ifndef _userFunc_h_
#define _userFunc_h_

#include "vader_common.h"

/**************************************************************************/
/* Interface declarations for user-defined functions for alpha, source    */
/* terms, boundary conditions, pre- and post-timestep operations, and     */
/* checkpoint reading and writing operations.                             */
/**************************************************************************/

void
userAlpha(
	  /* Inputs */
	  const double t, const double dt, const grid *grd, 
	  const double *col, const double *pres,
	  const double *eInt, const double *gamma, 
	  const double *delta, void *params,
	  /* Outputs */
	  double *alpha
	  );

void
userEOS(
	/* Inputs */
	const double t, const double dt, const grid *grd, 
	const double *col, const double *pres, const double *eInt,
	void *params,
	/* Outputs */
	double *gamma, double *delta
	);

void
userMassSrc(
	    /* Inputs */
	    const double t, const double dt, const grid *grd,
	    const double *col, const double *pres,
	    const double *eInt, const double *gamma, 
	    const double *delta, void *params,
	    /* Outputs */
	    double *massSrc
	    );

void
userIntEnSrc(
	     /* Inputs */
	     const double t, const double dt, const grid *grd,
	     const double *col, const double *pres,
	     const double *eInt, const double *gamma, 
	     const double *delta, void *params, 
	     /* Outputs */
	     double *intEnSrc
	     );

void
userIBC(
	/* Inputs */
	const double t, const double dt, const grid *grd,
	const double *col, const double *pres, const double *eInt,
	const double *gamma, const double *delta,
	const pres_bc_type ibc_pres, const enth_bc_type ibc_enth,
	void *params, 
	/* Outputs */
	double *ibc_pres_val, double *ibc_enth_val
	);

void
userOBC(
	/* Inputs */
	const double t, const double dt, const grid *grd,
	const double *col, const double *pres, const double *eInt,
	const double *gamma, const double *delta,
	const pres_bc_type obc_pres, const enth_bc_type obc_enth,
	void *params, 
	/* Outputs */
	double *obc_pres_val, double *obc_enth_val
	);

void
userPreTimestep(
		/* Inputs / outputs */
		const double t, const double dt,
		const grid *grd, double *col, double *pres,
		double *eInt, double *mBnd, double *eBnd,
		double *mSrc, double *eSrc,
		void *params, const unsigned long nUserOut,
		double *userOut
		);

void
userPostTimestep(
		 /* Inputs / outputs */
		 const double t, const double dt,
		 const grid *grd, double *col, double *pres,
		 double *eInt, double *mBnd, double *eBnd,
		 double *mSrc, double *eSrc,
		 void *params, const unsigned long nUserOut,
		 double *userOut
		 );

void
userCheckRead(
	      FILE *fp, grid *grd, const unsigned long nOut,
	      double *tOut, double *colOut,
	      double *presOut, double *eIntOut, double *mBndOut,
	      double *eBndOut, double *mSrcOut, double *eSrcOut,
	      const unsigned long nUserOut, double *userOut,
	      void *params
	      );

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
	      );

#endif
