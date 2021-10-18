#ifndef _advance_h_
#define _advance_h_

#include "vader_common.h"


/**************************************************************************/
/* This routines advances the computation one time step using             */
/* backwards Euler time centering                                         */
/**************************************************************************/

double 
advanceBE(
	  /* Time step and grid */
	  const double t, const double dt, const grid *grd, 
	  /* Starting data */
	  double *col, double *pres, double *eInt,
	  /* Records of mass and energy transported off grid and added
	     by sources */
	  double *mBnd, double *eBnd, double *mSrc, double *eSrc,
	  /* Equation of state parameters */
	  const bool eos_func, const double gamma_val, 
	  const double delta_val,
	  /* Dimensionless viscosity parameters */
	  const bool alpha_func, const double alpha_val,
	  /* Inner boundary condition parameters */
	  const pres_bc_type pres_ibc, const enth_bc_type enth_ibc,
	  const bool ibc_func, const double ibc_pres_val, 
	  const double ibc_enth_val,
	  /* Outer boundary condition parameters */
	  const pres_bc_type pres_obc, const enth_bc_type enth_obc,
	  const bool obc_func, const double obc_pres_val, 
	  const double obc_enth_val,
	  /* Source function parameters */
	  const bool massSrc_func, const double massSrc_val,
	  const bool intEnSrc_func, const double intEnSrc_val,
	  /* Control and method parameters */
	  const double errTol, const double dtTol, 
	  const unsigned long maxIter, const unsigned long interpOrder,
	  const bool noUpdate, const bool verbose, 
	  const wksp *w, 
	  /* User-defined extra parameters */
	  void *params,
	  /* Diagnostic outputs */
	  unsigned long *itCount
#ifdef TESTING_MODE
	  , double *resid, unsigned long *rtype,
	  double *advanceTime, double *nextIterTime,
	  double *userTime
#endif
	  );

#endif
