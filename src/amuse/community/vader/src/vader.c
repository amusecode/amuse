/**************************************************************************/
/* This is the main simulation routine. It initializes either from given  */
/* inputs or from a checkpoint, runs the simulation, and assembles the    */
/* final results.                                                         */
/**************************************************************************/

#include "init.h"
#include "vader.h"
#include "setup.h"
#include <string.h>
#include <gsl/gsl_errno.h>

double
vader(
      /* Restart checkpoint name (empty string for new start) */
      const char *restart_file,
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
      const bool preTimestep_func,
      const bool postTimestep_func,
      const unsigned long verbosity,
      /* Output control parameters */
      const unsigned long nSave, const double *tSave,
      const unsigned long nUserOut,
      const bool *userOutCum,
      const bool *writeCheckpoint,
      const char *checkname,
      const bool userReadCheckpoint,
      const bool userWriteCheckpoint,
      /* Computational grid (leave as NULL for a restart) */
      grid **grd_ptr,
      /* Starting data (leave as NULL for a restart) */
      double **col_init, double **pres_init, double **eInt_init,
      /* User-defined extra parameters */
      void *params,
      /* Diagnostic outputs */
      unsigned long *nStep, unsigned long *nIter,
      unsigned long *nFail,
      /* Output values; all the doubles should be passed in set to
	 NULL, and will be allocated appropriately */
      unsigned long *nOut,
      double **tOut_ptr, double **colOut_ptr,
      double **presOut_ptr, double **eIntOut_ptr,
      double **mBndOut_ptr, double **eBndOut_ptr, 
      double **mSrcOut_ptr, double **eSrcOut_ptr,
      double **userOut_ptr
#ifdef TESTING_MODE
      /* Parameters used in code tests */
      , double *residSum, unsigned long *iterStep,
      double *driverTime, double *advanceTime,
      double *nextIterTime, double *userTime
#endif
      ) {

  setup_status set_stat;
  gsl_error_handler_t *gslerr;
  double *col, *pres, *eInt;
  double *tOut, *colOut, *presOut, *eIntOut, *mBndOut, *eBndOut,
    *mSrcOut, *eSrcOut, *userOut;
  const double *tSave_tmp;
  const bool *writeCheckpoint_tmp, *userOutCum_tmp;
  grid *grd;
  unsigned long nOutMax, savePtr, chkNum;
  wksp *w;
  bool writeOut, restart;
  double t, dt, tSimEnd;

  /* Turn off GSL error handling, since we'll handle errors manually */
  gslerr = gsl_set_error_handler_off();

  /* Sanity check on the interpolation order */
  if ((interpOrder < 1) || (interpOrder > 3)) {
    fprintf(stderr, "vader: error: interpOrder must be 1, 2, or 3\n");
    return(-1.0);
  }

  /* Do the initial setup */
  if (!restart_file) restart = false;
  else if (strlen(restart_file) == 0) restart = false;
  else restart = true;
  if (!restart) {

    /* We're starting new, so the initial data and grid should already
       be present, so make our working pointers point to them */
    grd = *grd_ptr;
    col = *col_init;
    pres = *pres_init;
    if (eos_func) eInt = *eInt_init;
    else eInt = NULL;

    /* Initialize time and counters */
    t = tStart;
    *nStep = 1;
    *nIter = *nFail = 0;
    nOutMax = nSave;
    *nOut = savePtr = chkNum = 0;

    /* Call the setup routine; this stores the first outputs if
       necessary and computes the first time step */
    set_stat =
      setup_new(tStart, tEnd,
		eos_func, gamma_val, delta_val, alpha_func, alpha_val,
		ibc_pres, ibc_enth, ibc_func, ibc_pres_val, ibc_enth_val,
		obc_pres, obc_enth, obc_func, obc_pres_val, obc_enth_val,
		massSrc_func, massSrc_val, intEnSrc_func, intEnSrc_val,
		dtStart, dtMin, dtTol, errTol, maxDtIncrease,
		maxIter, interpOrder, maxStep, useBE,
		userWriteCheckpoint, verbosity,
		nSave, tSave, nUserOut, writeCheckpoint,
		checkname,
		grd, &w, col, pres, eInt, params,
		tOut_ptr, colOut_ptr, presOut_ptr, eIntOut_ptr, mBndOut_ptr,
		eBndOut_ptr, mSrcOut_ptr, eSrcOut_ptr, userOut_ptr,
		&dt, nOut, &savePtr, &writeOut, &chkNum
#ifdef TESTING_MODE
		, residSum, iterStep, advanceTime, nextIterTime,
		userTime
#endif		
		);
    if (set_stat != GOOD_START) return -LARGE;

  } else {

    /* Set up from the checkpoint */
    set_stat =
      setup_checkpoint(restart_file, tStart, tEnd,
		       eos_func, massSrc_func, intEnSrc_func,
		       userReadCheckpoint, verbosity, nSave, tSave,
		       nUserOut, userOutCum, &grd, &col, &pres, &eInt,
		       nStep, nIter, nFail,
		       tOut_ptr, colOut_ptr, presOut_ptr, eIntOut_ptr,
		       mBndOut_ptr, eBndOut_ptr, mSrcOut_ptr, eSrcOut_ptr,
		       userOut_ptr, params, &w, &t, &dt, nOut, &savePtr,
		       &nOutMax, &writeOut, &chkNum);
    if (set_stat != GOOD_START) return -LARGE;

    /* Point inputs to allocated memory */
    *grd_ptr = grd;
    *col_init = col;
    *pres_init = pres;
    if (eos_func) *eInt_init = eInt;
  }

  /* Set pointers to output blocks and output controls; note that some
     of the control flags might be NULL, and we have to make sure to
     preserve those values, because they are significant in driver */
  tOut = *tOut_ptr;
  colOut = *colOut_ptr;
  presOut = *presOut_ptr;
  if (eos_func) eIntOut = *eIntOut_ptr; else eIntOut = NULL;
  mBndOut = *mBndOut_ptr;
  eBndOut = *eBndOut_ptr;
  if (massSrc_func) mSrcOut = *mSrcOut_ptr; else mSrcOut = NULL;
  if (massSrc_func || intEnSrc_func) eSrcOut = *eSrcOut_ptr;
  else eSrcOut = NULL;
  if (nUserOut > 0) userOut = *userOut_ptr; else userOut = NULL;
  if (tSave) tSave_tmp = tSave + savePtr; else tSave_tmp = NULL;
  if (writeCheckpoint) writeCheckpoint_tmp = writeCheckpoint + savePtr;
  else writeCheckpoint_tmp = NULL;
  if (userOutCum) userOutCum_tmp = userOutCum + savePtr;
  else userOutCum_tmp = NULL;

  /* Now call the driver routine to run the simulation */
  tSimEnd = driver(t, tEnd, eos_func, gamma_val, delta_val,
		   alpha_func, alpha_val,
		   ibc_pres, ibc_enth, ibc_func,
		   ibc_pres_val, ibc_enth_val,
		   obc_pres, obc_enth, obc_func,
		   obc_pres_val, obc_enth_val,
		   massSrc_func, massSrc_val,
		   intEnSrc_func, intEnSrc_val,
		   dt, dtMin, dtTol, errTol, maxDtIncrease,
		   maxIter, interpOrder, maxStep,
		   useBE, preTimestep_func, postTimestep_func,
		   verbosity, nOutMax - *nOut,
		   tSave_tmp, nUserOut, userOutCum_tmp,
		   writeCheckpoint_tmp, checkname,
		   userWriteCheckpoint, writeOut, 
		   chkNum, grd, w, col, pres, eInt,
		   params, nStep, nIter, nFail,
		   nOut, tOut, colOut, presOut, eIntOut,
		   mBndOut, eBndOut, mSrcOut, eSrcOut,
		   userOut
#ifdef TESTING_MODE
		   , residSum, iterStep, driverTime, advanceTime,
		   nextIterTime, userTime
#endif
		   );

  /* If we didn't use up all the output slots because the calculation
     terminated early, free unneeded memory */
  if (*nOut < nOutMax) {
    outputResize(*nOut, eos_func, massSrc_func, intEnSrc_func,
		 nUserOut, grd, tOut_ptr, colOut_ptr, presOut_ptr,
		 eIntOut_ptr, mBndOut_ptr, eBndOut_ptr, mSrcOut_ptr,
		 eSrcOut_ptr, userOut_ptr);
  }

  /* Free the workspace */
  wkspFree(w);
  
  /* Return */
  return tSimEnd;
}
