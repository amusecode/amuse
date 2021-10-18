#include <math.h>
#include <gsl/gsl_errno.h>
#include "advanceCN.h"
#include "advanceBE.h"
#include "checkpoint.h"
#include "driver.h"
#include "setup.h"
#include "userFunc.h"

/**************************************************************************/
/* Driver routine                                                         */
/**************************************************************************/

double 
driver(
       /* Time parameters */
       const double tStart,
       const double tEnd,
       /* Equation of state parameters */
       const bool eos_func,
       const double gamma_val, 
       const double delta_val,
       /* Dimensionless viscosity parameters */
       const bool alpha_func,
       const double alpha_val,
       /* Inner boundary condition parameters */
       const pres_bc_type ibc_pres,
       const enth_bc_type ibc_enth,
       const bool ibc_func,
       const double ibc_pres_val, 
       const double ibc_enth_val,
       /* Outer boundary condition parameters */
       const pres_bc_type obc_pres,
       const enth_bc_type obc_enth,
       const bool obc_func,
       const double obc_pres_val, 
       const double obc_enth_val,
       /* Source function parameters */
       const bool massSrc_func,
       const double massSrc_val,
       const bool intEnSrc_func,
       const double intEnSrc_val,
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
       const unsigned long nSave,
       const double *tSave,
       const unsigned long nUserOut,
       const bool *userOutCum,
       const bool *writeCheckpoint,
       const char *checkname,
       const bool userWriteCheckpoint,
       const bool writeFirstStep,
       const unsigned long checknum,
       /* Computational grid and workspace */
       const grid *grd,
       const wksp *w, 
       /* Input data */
       double *col,
       double *pres,
       double *eInt,
       /* User-defined extra parameters */
       void *params,
       /* Diagnostic outputs */
       unsigned long *nStep,
       unsigned long *nIter,
       unsigned long *nFail,
       /* Storage for outputs */
       unsigned long *nOut,
       double *tOut,
       double *colOut,
       double *presOut,
       double *eIntOut,
       double *mBndOut,
       double *eBndOut,
       double *mSrcOut,
       double *eSrcOut,
       double *userOut
#ifdef TESTING_MODE
       /* Parameters used in code tests */
       , double *residSum, unsigned long *iterStep,
       double *driverTime, double *advanceTime,
       double *nextIterTime, double *userTime
#endif
       ) {
;
  status stat = RUNNING;
  double t = tStart;
  double dt = dtStart;
  bool writeOut = writeFirstStep;
  unsigned long savePtr = 0, chkPtr = checknum;
  double dtNew;
  unsigned long i, j;
  unsigned long nIterTmp, nreduce;
#ifdef TESTING_MODE
  unsigned long *residType;
  clock_t start_t, end_t;

  /* Start the clocks */
  start_t = clock();
  *advanceTime = *nextIterTime = *userTime = 0.0;

  /* In testing mode, allocate memory for diagnostic array */
  if (!(residType = calloc(maxIter, sizeof(unsigned long)))) {
    fprintf(stderr, "Error: unable to allocate memory for diagnostic arrays!\n");
    exit(1);
  }

  /* In testing mode, initialize sum of residuals and iteration
     counter arrays */
  for (i=0; i<maxIter; i++) residSum[i] = 0.0;
  for (i=0; i<maxStep; i++) iterStep[i] = 0;  
#endif

  /* Begin main loop */
  while (stat == RUNNING) {

    /* Print status */
    if (verbosity > 1) {
      printf("Step %ld, t = %e, dt = %e\n", *nStep, t, dt);
    } else if (verbosity > 0) {
      if (*nStep % 100 == 0) {
	printf("Step %ld, t = %e, dt = %e\n", *nStep, t, dt);
      }
    }

    /* Call user work function */
    if (preTimestep_func) {
		
      userPreTimestep(t, dt, grd, col, pres, eInt, mBndOut+2*(*nOut),
		      eBndOut+2*(*nOut), mSrcOut+grd->nr*(*nOut),
		      eSrcOut+grd->nr*(*nOut), params, nUserOut,
		      userOut+nUserOut*grd->nr*(*nOut));
    }

    /* Try to advance */
    dtNew = -1.0;
    nreduce = 0;
    while (dtNew < 0) {
      if (useBE == 0) { 

	dtNew = advanceCN(t, dt, grd, col, pres, eInt, 
			  mBndOut+2*(*nOut),
			  eBndOut+2*(*nOut),
			  mSrcOut+grd->nr*(*nOut),
			  eSrcOut+grd->nr*(*nOut),
			  eos_func, gamma_val, delta_val,
			  alpha_func, alpha_val, 
			  ibc_pres, ibc_enth,
			  ibc_func, ibc_pres_val, ibc_enth_val,
			  obc_pres, obc_enth,
			  obc_func, obc_pres_val, obc_enth_val,
			  massSrc_func, massSrc_val,
			  intEnSrc_func, intEnSrc_val,
			  errTol, dtTol, maxIter,
			  interpOrder, false, verbosity > 2, w, params,
			  &nIterTmp
#ifdef TESTING_MODE
			  , residSum, residType,
			  advanceTime, nextIterTime, userTime
#endif
			  );
      } else { 
	dtNew = advanceBE(t, dt, grd, col, pres, eInt, 
			  mBndOut+2*(*nOut),
			  eBndOut+2*(*nOut),
			  mSrcOut+grd->nr*(*nOut),
			  eSrcOut+grd->nr*(*nOut),
			  eos_func, gamma_val, delta_val,
			  alpha_func, alpha_val, 
			  ibc_pres, ibc_enth,
			  ibc_func, ibc_pres_val, ibc_enth_val,
			  obc_pres, obc_enth,
			  obc_func, obc_pres_val, obc_enth_val,
			  massSrc_func, massSrc_val,
			  intEnSrc_func, intEnSrc_val,
			  errTol, dtTol, maxIter,
			  interpOrder, false, verbosity > 2, w, params,
			  &nIterTmp
#ifdef TESTING_MODE
			  , residSum, residType,
			  advanceTime, nextIterTime, userTime
#endif
			  );
      }

      /* Upate iteration count. We do it here to ensure that we count
	 iterations even if they fail to converge. */
      *nIter += nIterTmp;
#ifdef TESTING_MODE
      iterStep[*nStep-1] += nIterTmp;
#endif
      if (dtNew < 0) {
	writeOut = false;
	dt = dt/2.0;
	nreduce++;
	(*nFail)++;
	if (verbosity > 1)
	  printf("   Iterative solver non-convergence! Reducing dt to %e.\n", dt);
	if (dt < dtMin*(tEnd - tStart)) {
	  stat = ZENO_ERROR;
	  dtNew = dt;
	  break;
	}
      }
    }
    for (i=0; i<nreduce; i++) dtNew = dtNew / 2;

    /* Update time and step counter */
    t += dt;
    (*nStep)++;

    /* Call user work function */
    if (postTimestep_func) {
      userPostTimestep(t, dt, grd, col, pres, eInt, mBndOut+2*(*nOut),
		       eBndOut+2*(*nOut), mSrcOut+grd->nr*(*nOut),
		       eSrcOut+grd->nr*(*nOut), params, nUserOut,
		       userOut+nUserOut*grd->nr*(*nOut));
    }

    /* If requested, scale back the next time step */
    if (dtNew/dt > maxDtIncrease) dtNew = maxDtIncrease * dt;

    /* Store state */
    if (writeOut) {
      if (verbosity > 0)
	printf("Storing output %lu at step %lu, time %e\n", (*nOut)+1,
	       (*nStep)-1, t);
      tOut[(*nOut)] = t;
      for (i=0; i<grd->nr; i++) {
	colOut[i + (*nOut)*grd->nr] = col[i];
	presOut[i + (*nOut)*grd->nr] = pres[i];
	if (eos_func) eIntOut[i + (*nOut)*grd->nr] = eInt[i];
      }
      /* Boundary and source values are cumulative, so just copy the
      final value during this output interval into the slot for the
      next interval */
      if (savePtr<nSave-1) {
	mBndOut[2*((*nOut)+1)] = mBndOut[2*(*nOut)];
	mBndOut[2*((*nOut)+1)+1] = mBndOut[2*(*nOut)+1];
	eBndOut[2*((*nOut)+1)] = eBndOut[2*(*nOut)];
	eBndOut[2*((*nOut)+1)+1] = eBndOut[2*(*nOut)+1];
	if (massSrc_func == 1) {
	  for (i=0; i<grd->nr; i++) {
	    mSrcOut[((*nOut)+1)*grd->nr+i] = mSrcOut[(*nOut)*grd->nr+i];
	    eSrcOut[((*nOut)+1)*grd->nr+i] = eSrcOut[(*nOut)*grd->nr+i];
	  }
	} else if (intEnSrc_func == 1) {
	  for (i=0; i<grd->nr; i++) {
	    eSrcOut[((*nOut)+1)*grd->nr+i] = eSrcOut[(*nOut)*grd->nr+i];
	  }
	}
	/* Handle any user-defined cumulative outputs in the same way;
	   note that user outputs are a 3D array, with the dimensions
	   being (output time, output variable number, cell number),
	   where cell number is the fastest varying index and output
	   time is the slowest varying index in memory */
	for (j=0; j<nUserOut; j++) {
	  if (userOutCum) {
	    if (userOutCum[j]) {
	      for (i=0; i<grd->nr; i++)
		userOut[((*nOut)+1)*nUserOut*grd->nr + j*grd->nr + i] =
		  userOut[(*nOut)*nUserOut*grd->nr + j*grd->nr + i];
	    }
	  }
	}
      }

      /* Write checkpoint if requested */
      if (0&&writeCheckpoint) {
	if (writeCheckpoint[savePtr]) {
	  saveCheckpoint(checkname, chkPtr,
			 eos_func, massSrc_func, intEnSrc_func,
			 nUserOut, t, dtNew, *nStep, *nIter, *nFail, grd,
			 *nOut+1, tOut, colOut, presOut, eIntOut,
			 mBndOut, eBndOut, mSrcOut, eSrcOut, userOut,
#if AA_M > 0
			 w->constraint,
#endif
			 params, userWriteCheckpoint, verbosity);
	  /* Update checkpoint counter */
	  chkPtr++;
	}
      }

      /* Update counters and flags */
      (*nOut)++;
      savePtr++;
      writeOut = false;
    }

    /* Check termination conditions */
    if ((*nStep > maxStep) && (maxStep > 0)) stat = TOO_MANY_STEPS;
    if (dtNew < dtMin*(tEnd-tStart)) stat = ZENO_ERROR;
    if (t >= tEnd) stat = NORMAL_EXIT;

    /* If necessary, reduce next time step to hit the next output
       time or the end time */
    if (savePtr < nSave) {
      if (t + dtNew > tSave[savePtr]) {
	dtNew = tSave[savePtr] - t;
	writeOut = true;
      }
    }
    if (t + dtNew >= tEnd) dtNew = (1.0+1.0e-10)*(tEnd - t);
    dt = dtNew;

  } /* Main loop end */
  
#ifdef TESTING_MODE
  /* In testing mode, free memory for diagnostic array */
  free(residType);
#endif

  /* On early exit, store final state in last output */
  if (stat != NORMAL_EXIT) {
    if (savePtr < nSave) {
      tOut[(*nOut)] = t;
      for (i=0; i<grd->nr; i++) {
	colOut[i + (*nOut)*grd->nr] = col[i];
	presOut[i + (*nOut)*grd->nr] = pres[i];
	if (eos_func) eIntOut[i + (*nOut)*grd->nr] = eInt[i];
      }
      (*nOut)++;
    }
  }

  /* Print final status and return */
  if (verbosity > 0) {
    *nStep = *nStep-1;
    if (stat == NORMAL_EXIT)
      printf("Finished computation in %ld steps, t = %e\n", *nStep, t);
    else if (stat == TOO_MANY_STEPS)
      printf("Reached maximum number of steps = %ld, t = %e\n", *nStep, t);
    else if (stat == ZENO_ERROR)
      printf("Time step too small after %ld steps, t = %e, dt = %e\n",
	     *nStep, t, dt);
  }
#ifdef TESTING_MODE
  end_t = clock();
  *driverTime = (double) (end_t - start_t) / CLOCKS_PER_SEC;
#endif
  return(t);
}


