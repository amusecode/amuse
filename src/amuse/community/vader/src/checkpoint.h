#ifndef _checkpoint_h_
#define _checkpoint_h_

#include "vader_common.h"

/**************************************************************************/
/* Routines to write and read checkpoints                                 */
/**************************************************************************/

iostatus readCheckpoint(
			/* Checkpoint file name and number */
			const char *checkname,
			unsigned long *checknum,
			/* Control parameters */
			bool *eos_func, bool *massSrc_func,
			bool *intEnSrc_func,
			/* Output control parameters */
			unsigned long *nUserOut, 
			/* Simulation state parameters */
			double *t, double *dt,
			/* Diagnostic outputs */
			unsigned long *nStep,
			unsigned long *nIter,
			unsigned long *nFail,
			/* The grid */
			grid **grd, wksp **w,
			/* Simulation outputs */
			unsigned long *nOut,
			double **tOut, double **colOut,
			double **presOut, double **eIntOut,
			double **mBndOut, double **eBndOut,
			double **mSrcOut, double **eSrcOut,
			double **userOut,
			/* Control options */
			void *params,
			const bool userReadCheckpoint,
			const unsigned long verbosity);

iostatus saveCheckpoint(
			/* Checkpoint file name */
			const char *checkname,
			const unsigned long checknum,
			/* Control parameters */
			const bool eos_func,
			const bool massSrc_func,
			const bool intEnSrc_func,
			/* Output control parameters */
			const unsigned long nUserOut, 
			/* Simulation state parameters */
			const double t,
			const double dt,
			/* Diagnostic outputs */
			const unsigned long nStep,
			const unsigned long nIter,
			const unsigned long nFail,
			/* The grid */
			const grid *grd,
			/* Simulation outputs */
			const unsigned long nOut,
			const double *tOut,
			const double *colOut,
			const double *presOut, 
			const double *eIntOut,
			const double *mBndOut,
			const double *eBndOut,
			const double *mSrcOut,
			const double *eSrcOut,
			const double *userOut,
#if AA_M > 0
			const gsl_vector *constraint,
#endif
			/* Conrol options */
			const void *params,
			const bool userWriteCheckpoint,		   
			const unsigned long verbosity);

#endif
/* _checkpoint_h_ */
