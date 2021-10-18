/**************************************************************************/
/* This is the main simulation routine. It initializes either from given  */
/* inputs or from a checkpoint, runs the simulation, and assembles the    */
/* final results.                                                         */
/**************************************************************************/

#ifndef _vader_h_
#define _vader_h_

#include "vader_common.h"
#include <stdbool.h>

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
      const unsigned long nSave,
      const double *tSave,
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
      double **tOut, double **colOut,
      double **presOut, double **eIntOut,
      double **mBndOut, double **eBndOut, 
      double **mSrcOut, double **eSrcOut,
      double **userOut
#ifdef TESTING_MODE
      /* Parameters used in code tests */
      , double *residSum, unsigned long *iterStep,
      double *driverTime, double *advanceTime,
      double *nextIterTime, double *userTime
#endif
      );

#endif 
/* end _vader_h_ */
