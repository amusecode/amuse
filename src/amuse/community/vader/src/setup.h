#ifndef _setup_h_
#define _setup_h_

#include "driver.h"
#include "vader_common.h"

/**************************************************************************/
/* Routines to set up a calculation; called by driver                     */
/**************************************************************************/

/* Set up a calculation from a checkpoint */
setup_status
setup_checkpoint(
		 /* Restart checkpoint name */
		 const char *checkname,
		 /* Time parameters */
		 const double tStart, const double tEnd,
		 /* Equation of state parameters */
		 const bool eos_func, 
		 /* Source function parameters */
		 const bool massSrc_func,
		 const bool intEnSrc_func,
		 /* Control and method parameters */
		 const bool userReadCheckpoint,
		 const unsigned long verbosity,
		 /* Output control parameters  */
		 const unsigned long nSave,
		 const double *tSave,
		 const unsigned long nUserOut,
		 const bool *userOutCum,
		 /* Computational grid (should be NULL) */
		 grid **grd,
		 /* Current state (should be NULL) */
		 double **col, double **pres, double **eInt,
		 /* Diagnostic outputs */
		 unsigned long *nStep, unsigned long *nIter,
		 unsigned long *nFail,
		 /* Output values; double pointers should all be NULL */
		 double **tOut, double **colOut,
		 double **presOut, double **eIntOut,
		 double **mBndOut, double **eBndOut, 
		 double **mSrcOut, double **eSrcOut,
		 double **userOut,
		 /* User-defined parameters */
		 void *params,
		 /* Workspace, to be allocated */
		 wksp **w,
		 /* Current time and time step, to be read */
		 double *t, double *dt,
		 /* Pointers to where we are in the output block and
		    save times */
		 unsigned long *outPtr,
		 unsigned long *savePtr,
		 /* Total size of output block */
		 unsigned long *outSize,
		 /* Flag to indicate if we should write out on the
		    next time step */
		 bool *writeOut,
		 /* Checkpoint counter, to be read */
		 unsigned long *chkNum);

/* Set up a new calculation */
setup_status
setup_new(
	  /* Time parameters */
	  const double tStart, const double tEnd,
	  /* Equation of state parameters */
	  const bool eos_func, const double gamma_val, 
	  const double delta_val,
	  /* Dimensionless viscosity parameters */
	  const bool alpha_func, const double alpha_val,
	  /* Inner boundary condition parameters */
	  const pres_bc_type ibc_pres, const enth_bc_type ibc_enth,
	  const bool ibc_func, const double ibc_pres_val, 
	  const double ibc_enth_val,
	  /* Outer boundary condition parameters */
	  const pres_bc_type obc_pres, const enth_bc_type obc_enth,
	  const bool obc_func, const double obc_pres_val, 
	  const double obc_enth_val,
	  /* Source function parameters */
	  const bool massSrc_func, const double massSrc_val,
	  const bool intEnSrc_func, const double intEnSrc_val,
	  /* Control and method parameters */
	  const double dtStart,
	  const double dtMin, 
	  const double dtTol, 
	  const double errTol,
	  const double maxDtIncrease,
	  const unsigned long maxIter,
	  const unsigned long interpOrder,
	  const unsigned long maxStep,
	  const bool useBE,
	  const bool userWriteCheckpoint,
	  const unsigned long verbosity,
	  /* Output control parameters */
	  const unsigned long nSave,
	  const double *tSave,
	  const unsigned long nUserOut,
	  const bool *writeCheckpoint,
	  const char *checkname,
	  /* Computational grid and workspace */
	  grid *grd, wksp **w_ptr, 
	  /* Starting data */
	  double *col, double *pres, double *eInt,
	  /* User-defined extra parameters */
	  void *params,
	  /* Outputs; all arrays to be allocated */
	  double **tOut_ptr, double **colOut_ptr,
	  double **presOut_ptr, double **eIntOut_ptr,
	  double **mBndOut_ptr, double **eBndOut_ptr, 
	  double **mSrcOut_ptr, double **eSrcOut_ptr,
	  double **userOut_ptr,
	  /* Initial time step, to be computed */
	  double *dt,
	  /* Pointers to where we are in the output block and
	     save times */
	  unsigned long *outPtr,
	  unsigned long *savePtr,
	  /* Flag to indicate we should write on the next time step */
	  bool *writeOut,
	  /* Checkpoint counter */
	  unsigned long *chkNum
#ifdef TESTING_MODE
	  /* Parameters used in code tests */
	  , double *residSum, unsigned long *iterStep,
	  double *advanceTime, double *nextIterTime,
	  double *userTime
#endif
	  );

#endif
/* _setup_h_ */
