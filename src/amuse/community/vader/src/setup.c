#include "setup.h"
#include "advanceCN.h"
#include "advanceBE.h"
#include "checkpoint.h"
#include "init.h"
#include <stdio.h>

/**************************************************************************/
/* Routine to set up a calculation from a checkpoint                      */
/**************************************************************************/

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
		 /* Total size of output block needed for calculation */
		 unsigned long *outSize,
		 /* Flag to indicate if we should write out on the
		    next time step */
		 bool *writeOut,
		 /* Checkpoint counter, to be read */
		 unsigned long *chkNum) {

  unsigned long nOut_chk, nUserOut_chk;
  bool eos_func_chk, massSrc_func_chk, intEnSrc_func_chk;
  iostatus readstat;
  unsigned long i, j;

  /* First step: read the checkpoint; bail if we fail */
  readstat = readCheckpoint(checkname, chkNum, &eos_func_chk,
			    &massSrc_func_chk, &intEnSrc_func_chk,
			    &nUserOut_chk, t, dt, nStep, nIter, nFail,
			    grd, w, &nOut_chk, tOut, colOut, presOut,
			    eIntOut, mBndOut, eBndOut,
			    mSrcOut, eSrcOut, userOut,
			    params, userReadCheckpoint, verbosity);
  if (readstat != GOOD_IO) return RESTART_ERROR;
  
  /* Sanity check: make sure that the flags we were passed are
     consistent with those in the checkpoint we just read. If not,
     deallocate and bail. */
  if ((nUserOut != nUserOut_chk) || (eos_func != eos_func_chk) ||
      (massSrc_func != massSrc_func_chk) ||
      (intEnSrc_func != intEnSrc_func_chk)) {
    fprintf(stderr, "vader: setup_checkpoint: error: checkpoint inconsistent with control flags -- values are (checkpoint, input): nUserOut = (%ld %ld), eos_func = (%d %d), massSrc_func = (%d %d), intEnSrcFunc = (%d %d)\n",
	    nUserOut_chk, nUserOut, eos_func_chk, eos_func, massSrc_func_chk,
	    massSrc_func, intEnSrc_func_chk, intEnSrc_func);
    outputFree(tOut, colOut, presOut, eIntOut, mBndOut, eBndOut,
	       mSrcOut, eSrcOut, userOut);
    wkspFree(*w);
    gridFree(*grd);
    *w = NULL;
    *grd = NULL;
    return RESTART_ERROR;
  }

  /* Sanity check: make sure the time in the checkpoint is between the
     start and end times */
  if ((*t < tStart) || (*t > tEnd)) {
    fprintf(stderr, "vader: setup_checkpoint: error: checkpoint time %e is not between the input start and end times %e - %e\n", *t, tStart, tEnd);
    outputFree(tOut, colOut, presOut, eIntOut, mBndOut, eBndOut,
	       mSrcOut, eSrcOut, userOut);
    wkspFree(*w);
    gridFree(*grd);
    *w = NULL;
    *grd = NULL;
    return RESTART_ERROR;
  }
  
  /* Figure out where we are in the list of output times that we were
     passed; these need not match up the times that are listed in the
     checkpoint, but we issue a warning the user requests a time that
     is before this checkpoint. */
  *outPtr = *savePtr = 0;
  while (*outPtr < nOut_chk) {
    if ((*tOut)[*outPtr] == tSave[*savePtr]) {
      /* Times match, so increment both pointers */
      (*savePtr)++;
      (*outPtr)++;
    } else if (tSave[*savePtr] < (*tOut)[*outPtr]) {
      /* Requested output at a time that is not present in the
	 checkpoint; issue a warning, increment the output pointer,
	 and continue */
      fprintf(stderr, "vader: setup_checkpoint: warning: requested output time %e not found in checkpoint file %s, will not be include in subsequent output\n",
	      tSave[*savePtr], checkname);
      (*savePtr)++;
    } else {
      /* Output time found in checkpoint was not requested; this is
	 not a problem, so just move the pointer in the checkpoint */
      (*outPtr)++;
    }
  }

  /* Compute the size of the output block we need, and allocate it */
  *outSize = nSave - *savePtr + *outPtr;
  if (!outputResize(*outSize, eos_func, massSrc_func, intEnSrc_func,
		    nUserOut, *grd, tOut, colOut, presOut, eIntOut,
		    mBndOut, eBndOut, mSrcOut, eSrcOut, userOut)) {
    wkspFree(*w);
    gridFree(*grd);
    *w = NULL;
    *grd = NULL;
    return MEMORY_ERROR;
  }  

  /* Allocate memory to hold the working column density, pressure, and
     internal energy arrays, and then fill them from the last output */
  *col = (double *) calloc((*grd)->nr, sizeof(double));
  *pres = (double *) calloc((*grd)->nr, sizeof(double));
  if (eos_func)
    *eInt = (double *) calloc((*grd)->nr, sizeof(double));
  else if (eInt)
    *eInt = NULL;
  if (!col || !pres || (eos_func && !eInt)) {
    fprintf(stderr, "vader: setup_checkpoint: error: unable to allocate memory\n");
    if (col) free(*col);
    if (pres) free(*pres);
    if (eInt) free(*eInt);
    return MEMORY_ERROR;
  }
  for (i=0; i<(*grd)->nr; i++) {
    (*col)[i] = (*colOut)[i+(*grd)->nr*(*outPtr-1)];
    (*pres)[i] = (*presOut)[i+(*grd)->nr*(*outPtr-1)];
    if (eos_func)
      (*eInt)[i] = (*eIntOut)[i+(*grd)->nr*(*outPtr-1)];
  }

  /* For output quantities that should be computed cumulatively,
     initialize the next output holders into which we'll be
     accumulating */
  if (*outPtr < *outSize-1) {
    (*mBndOut)[2*(*outPtr+1)] = (*mBndOut)[2*(*outPtr)];
    (*mBndOut)[2*(*outPtr+1)+1] = (*mBndOut)[2*(*outPtr)+1];
    (*eBndOut)[2*(*outPtr+1)] = (*eBndOut)[2*(*outPtr)];
    (*eBndOut)[2*(*outPtr+1)+1] = (*eBndOut)[2*(*outPtr)+1];
    if (massSrc_func) {
      for (i=0; i<(*grd)->nr; i++) {
	(*mSrcOut)[(*outPtr+1)*(*grd)->nr+i] =
	  (*mSrcOut)[(*outPtr)*(*grd)->nr+i];
	(*eSrcOut)[(*outPtr+1)*(*grd)->nr+i] =
	  (*eSrcOut)[(*outPtr)*(*grd)->nr+i];
      }
    } else if (intEnSrc_func) {
      for (i=0; i<(*grd)->nr; i++) {
	(*eSrcOut)[(*outPtr+1)*(*grd)->nr+i] =
	  (*eSrcOut)[(*outPtr)*(*grd)->nr+i];
      }
    }
    for (j=0; j<nUserOut; j++) {
      if (userOutCum) {
	if (userOutCum[j]) {
	  for (i=0; i<(*grd)->nr; i++) {
	    (*userOut)[(*outPtr+1)*(*grd)->nr+i] =
	      (*userOut)[(*outPtr)*(*grd)->nr+i];
	  }
	}
      }
    }
  }
	
  /* Scale back time step if necessary */
  *writeOut = 0;
  if (*outPtr < *outSize) {
    if (*t + *dt > tSave[*savePtr]) {
      *dt = tSave[*savePtr] - *t;
      *writeOut = true;
    }
  }
  if (*t + *dt > tEnd) *dt = tEnd - *t;

  /* Increment the checkpoint counter */
  (*chkNum)++;

  /* Return */
  return GOOD_START;
}


/**************************************************************************/
/* Routine to set up a new calculation                                    */
/**************************************************************************/
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
	  /* Pointer to where we are in the output block and save
	     times */
	  unsigned long *outPtr,
	  unsigned long *savePtr,
	  /* Flag to indicate if we should write out on the next time
	     step */
	  bool *writeOut,
	  /* Checkpoint counter */
	  unsigned long *chkNum
#ifdef TESTING_MODE
	  /* Parameters used in code tests */
	  , double *residSum, unsigned long *iterStep,
	  double *advanceTime, double *nextIterTime,
	  double *userTime
#endif
	  ) {

  unsigned long i, nreduce, nIterTmp;
  double dtNew;
  double *tOut, *colOut, *presOut, *eIntOut, *mBndOut, *eBndOut,
    *mSrcOut, *eSrcOut, *userOut;
  wksp *w;
#ifdef TESTING_MODE
  unsigned long *residType = NULL;
#endif

  /* Allocate memory for workspace and to store the outputs */
  *w_ptr = wkspAlloc(grd->nr);
  if (!(*w_ptr)) return MEMORY_ERROR;
  if (!outputAlloc(nSave, eos_func, massSrc_func, intEnSrc_func,
		   nUserOut, grd, tOut_ptr, colOut_ptr, presOut_ptr,
		   eIntOut_ptr, mBndOut_ptr, eBndOut_ptr, mSrcOut_ptr,
		   eSrcOut_ptr, userOut_ptr))
    return MEMORY_ERROR;

  /* Set some temporary pointers */
  w = *w_ptr;
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

  /* Initialize constraint vector */
#if AA_M > 0
  gsl_vector_set_zero(w->constraint);
#endif

  /* Determine initial time step */
  if (dtStart > 0) {
    dtNew = dtStart;
    *dt = dtStart;
  } else {
#if TESTING_MODE
    residType = (unsigned long *) calloc(maxIter, sizeof(unsigned long));
#endif
    *dt = 1.0e-4*grd->r_h[0]/grd->vphi_h[0];
    dtNew = -1.0;
    nreduce = 0;
    while (dtNew < 0) {
      if (useBE == 0) {
	dtNew = advanceCN(tStart, *dt, grd, col, pres, eInt, 
			  mBndOut, eBndOut, mSrcOut, eSrcOut,
			  eos_func, gamma_val, delta_val,
			  alpha_func, alpha_val, 
			  ibc_pres, ibc_enth,
			  ibc_func, ibc_pres_val, ibc_enth_val,
			  obc_pres, obc_enth,
			  obc_func, obc_pres_val, obc_enth_val,
			  massSrc_func, massSrc_val,
			  intEnSrc_func, intEnSrc_val,
			  errTol, dtTol, maxIter,
			  interpOrder, true, verbosity > 2, w, params,
			  &nIterTmp
#ifdef TESTING_MODE
			  , residSum, residType,
			  advanceTime, nextIterTime, userTime
#endif
			  );
      } else {
	dtNew = advanceBE(tStart, *dt, grd, col, pres, eInt, 
			  mBndOut, eBndOut, mSrcOut, eSrcOut,
			  eos_func, gamma_val, delta_val,
			  alpha_func, alpha_val, 
			  ibc_pres, ibc_enth,
			  ibc_func, ibc_pres_val, ibc_enth_val,
			  obc_pres, obc_enth,
			  obc_func, obc_pres_val, obc_enth_val,
			  massSrc_func, massSrc_val,
			  intEnSrc_func, intEnSrc_val,
			  errTol, dtTol, maxIter,
			  interpOrder, true, verbosity > 2, w, params,
			  &nIterTmp
#ifdef TESTING_MODE
			  , residSum, residType,
			  advanceTime, nextIterTime, userTime
#endif
			  );	
      }
      
      if (dtNew < 0) {
	*dt = *dt/2.0;
	nreduce++;
      }
      if (*dt < dtMin*(tEnd - tStart)) {
	fprintf(stderr, "vader: setup_new: initial time step too small, dt = %e\n",
		*dt);
#ifdef TESTING_MODE
	free(residType);
#endif
	return FIRST_DT_ERROR;
      }
    }
    *dt = dtNew;
    for (i=0; i<nreduce; i++) *dt = *dt/2;
  }
  
  /* If requested, store the initial state and write an initial checkpoint */
  if (nSave > 0) {
    if (tSave[0] == tStart) {
      if (verbosity > 0)
	printf("Storing output 1 at step 0, time %e\n", tStart);
      tOut[0] = tStart;
      for (i=0; i<grd->nr; i++) {
	colOut[i] = col[i];
	presOut[i] = pres[i];
	if (eos_func) eIntOut[i] = eInt[i];
      }
      (*outPtr)++;
      (*savePtr)++;
      if (writeCheckpoint != NULL) {
	if (writeCheckpoint[0]) {
	  saveCheckpoint(/* Name and number */
			 checkname, 0,
			 /* Control parameters */
			 eos_func,
			 massSrc_func,
			 intEnSrc_func,
			 /* Output control parameters */
			 nUserOut,
			 /* Time and timestep */
			 tStart, dtNew,
			 /* Diagnostic outputs */
			 0, 0, 0,
			 /* Computational grid */
			 grd,
			 /* Simulation outputs */
			 1, tOut, colOut, presOut,
			 eIntOut, mBndOut, eBndOut, mSrcOut, eSrcOut,
			 userOut,
#if AA_M > 0
			 w->constraint,
#endif
			 /* Control parameters */
			 params,
			 userWriteCheckpoint,
			 verbosity);
	  (*chkNum)++;
	}
      }
    }
  }

  /* Scale back time step if necessary */
  *writeOut = 0;
  if (*outPtr < nSave) {
    if (tStart + *dt > tSave[*savePtr]) {
      *dt = tSave[*savePtr] - tStart;
      *writeOut = 1;
    }
  }
  if (tStart + *dt > tEnd) *dt = tEnd - tStart;

  /* Initialize boundary and source values to zero */
  mBndOut[0] = mBndOut[1] = eBndOut[0] = eBndOut[1] = 0.0;
  if (massSrc_func == 1) {
    for (i=0; i<grd->nr; i++)
      mSrcOut[i] = eSrcOut[i] = 0.0;
  } else if (intEnSrc_func == 1) {
    for (i=0; i<grd->nr; i++)
      eSrcOut[i] = 0.0;
  }

  /* Free temporary array used in testing mode */
#ifdef TESTING_MODE
  free(residType);
#endif

  /* Return */
  return GOOD_START;
}
