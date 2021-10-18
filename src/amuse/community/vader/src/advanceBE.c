#include <math.h>
#include <gsl/gsl_linalg.h>
#include "advanceBE.h"
#include "applyBC.h"
#include "userFunc.h"
#include "ppmExtrap.h"
#include "getNextIterate.h"

/**************************************************************************/
/* Time step advance routine with backwards Euler time differencing       */
/**************************************************************************/

double 
advanceBE(const double t, const double dt, const grid *grd, 
	  double *col, double *pres, double *eInt,
	  double *mBnd, double *eBnd, double *mSrc, double *eSrc,
	  const bool eos_func, const double gamma_val, 
	  const double delta_val,
	  const bool alpha_func, const double alpha_val,
	  const pres_bc_type ibc_pres, const enth_bc_type ibc_enth,
	  const bool ibc_func, const double ibc_pres_val, 
	  const double ibc_enth_val,
	  const pres_bc_type obc_pres, const enth_bc_type obc_enth,
	  const bool obc_func, const double obc_pres_val, 
	  const double obc_enth_val,
	  const bool massSrc_func, const double massSrc_val,
	  const bool intEnSrc_func, const double intEnSrc_val,
	  const double errTol, const double dtTol, 
	  const unsigned long maxIter, const unsigned long interpOrder, 
	  const bool noUpdate, const bool verbose,
	  const wksp *w, void *params, unsigned long *itCount
#ifdef TESTING_MODE
	  , double *resid, unsigned long *rtype,
	  double *advanceTime, double *nextIterTime,
	  double *userTime
#endif
	  ) {

  double pIn, qIn, pOut, qOut;
  double err, errCell;
  double dPres, dCol, dEInt, dtMin;
  double ibc_pres_val1, ibc_enth_val1, obc_pres_val1, obc_enth_val1;
  double h_up, h_up1;
  double mBndTemp[2], eBndTemp[2];
  unsigned long i, maxErrCell, converged, residType;
#if AA_M > 0
  long nHist;
  double residMax;
#endif
#ifdef TESTING_MODE
  clock_t start_t, end_t, user_start_t, user_end_t;
  clock_t next_iter_start_t, next_iter_end_t;

  start_t = clock();
#endif

  /****************************************************************/
  /* Step 1: initialize workspace arrays from starting conditions */
  /****************************************************************/

  for (i=0; i<grd->nr; i++) {
    w->pres_g[i+1] = w->presNew_g[i+1] = pres[i];
    w->colNew[i] = col[i];
    w->intEnSrc[i] = 0.0;
  }
  /* Get starting gamma and delta values */
  if (eos_func == 1) {
#if AA_M > 0
    nHist = 3*grd->nr;
#endif
#ifdef TESTING_MODE
    user_start_t = clock();
#endif
    userEOS(t, dt, grd, col, pres, eInt, params,
	    w->gammaLast, w->deltaLast);
#ifdef TESTING_MODE
    user_end_t = clock();
    *userTime += (double) (user_end_t - user_start_t) 
      / CLOCKS_PER_SEC;
#endif
    for (i=0; i<grd->nr; i++) {
      w->eIntTmp[i] = w->eIntNew[i] = eInt[i];
      w->gammaNew[i] = w->gammaLast[i];
      w->deltaNew[i] = w->deltaLast[i];
    }
  } else {
#if AA_M > 0
    nHist = 2*grd->nr;
#endif
    for (i=0; i<grd->nr; i++) {
      w->gammaLast[i] = w->gammaNew[i] = gamma_val;
      w->deltaLast[i] = w->deltaNew[i] = delta_val;
    }
  }
  /* Get starting viscosity */
  if (alpha_func == 1) {
#ifdef TESTING_MODE
    user_start_t = clock();
#endif
    userAlpha(t, dt, grd, col, pres, eInt, w->gammaLast, w->deltaLast, 
	      params, w->alpha_g+1);
#ifdef TESTING_MODE
    user_end_t = clock();
    *userTime += (double) (user_end_t - user_start_t) 
      / CLOCKS_PER_SEC;
#endif
  } else {
    for (i=1; i<=grd->nr; i++) w->alpha_g[i] = alpha_val;
  }
  /* Store starting mass and energy fluxes through boundary */
  mBndTemp[0] = mBnd[0];
  mBndTemp[1] = mBnd[1];
  eBndTemp[0] = eBnd[0];
  eBndTemp[1] = eBnd[1];
  /* Store starting source function values */
  if (massSrc_func == 1) {
    for (i=0; i<grd->nr; i++) {
      w->mSrc[i] = mSrc[i];
      w->eSrc[i] = eSrc[i];
    }
  } else if (intEnSrc_func == 1) {
    for (i=0; i<grd->nr; i++) w->eSrc[i] = eSrc[i];
  }

  /****************************************************************/
  /* Step 2: initialize RHS                                       */
  /****************************************************************/

  for (i=0; i<grd->nr; i++) {
    gsl_vector_set(w->rhs_g, i+1, pres[i]);
  }

  /****************************************************************/
  /* Step 3: get initial guess at mass fluxes; this is needed to  */
  /* upwind the enthalpy properly during the first loop iteration */
  /****************************************************************/

  /* Set alpha in ghost cells */
  w->alpha_g[0] = w->alpha_g[1];
  w->alpha_g[grd->nr+1] = w->alpha_g[grd->nr];

  /* Get pressure / enthalpy boundary values */
  if (ibc_func == 1) {
#ifdef TESTING_MODE
    user_start_t = clock();
#endif
    userIBC(t, dt, grd, col, pres, eInt, w->gammaLast, w->deltaLast, 
	    ibc_pres, ibc_enth, params, 
	    &ibc_pres_val1, &ibc_enth_val1);
#ifdef TESTING_MODE
    user_end_t = clock();
    *userTime += (double) (user_end_t - user_start_t) 
      / CLOCKS_PER_SEC;
#endif
  } else {
    ibc_pres_val1 = ibc_pres_val;
    ibc_enth_val1 = ibc_enth_val;
  }
  if (obc_func == 1) {
#ifdef TESTING_MODE
    user_start_t = clock();
#endif
    userOBC(t, dt, grd, col, pres, eInt, w->gammaLast, w->deltaLast, 
	    obc_pres, obc_enth, params, 
	    &obc_pres_val1, &obc_enth_val1);
#ifdef TESTING_MODE
    user_end_t = clock();
    *userTime += (double) (user_end_t - user_start_t) 
      / CLOCKS_PER_SEC;
#endif
  } else {
    obc_pres_val1 = obc_pres_val;
    obc_enth_val1 = obc_enth_val;
  }

  /* Set boundary conditions from user-specified parameters */
  applyBC(grd, w->alpha_g, 
	  ibc_pres, ibc_enth, ibc_pres_val1, ibc_enth_val1,
	  obc_pres, obc_enth, obc_pres_val1, obc_enth_val1,
	  w->pres_g, w->hint_g, &pIn, &qIn, &pOut, &qOut);

  /* Get old time mass flux */
  for (i=0; i<=grd->nr; i++) {
    w->fmNew_h[i] = -grd->g_h[i] *
      (SQR(grd->r_g[i+1])*(1-grd->beta_g[i+1]) * 
       w->pres_g[i+1]*w->alpha_g[i+1] - 
       SQR(grd->r_g[i])*(1-grd->beta_g[i]) * 
       w->pres_g[i]*w->alpha_g[i]);
  }

  /****************************************************************/
  /* Step 4: main iteration loop                                  */
  /****************************************************************/

  /* Status message */
  if (verbose != 0)
    printf("advanceBE: starting iteration, dt = %e\n", dt);

  /* Start loop */
  for (*itCount = 1, converged = 0; 1; (*itCount)++) {

    /**************************************************************/
    /* Step 4a: recompute alpha, enthalpy, BC's for current guess */
    /**************************************************************/

    /* Get new EOS parameters and cell-center enthalpies */
    if (eos_func == 1) {
#ifdef TESTING_MODE
      user_start_t = clock();
#endif
      userEOS(t+dt, dt, grd, w->colNew, w->presNew_g+1, 
	      w->eIntNew, params, w->gammaNew, w->deltaNew);
#ifdef TESTING_MODE
      user_end_t = clock();
      *userTime += (double) (user_end_t - user_start_t) 
	/ CLOCKS_PER_SEC;
#endif
      for (i=1; i<=grd->nr; i++) 
	w->hint_g[i] = (w->presNew_g[i] + w->eIntNew[i-1]) /
	  w->colNew[i-1];
    } else {
      for (i=1; i<=grd->nr; i++) {
	w->hint_g[i] = gamma_val/(gamma_val-1.0) * 
	  w->presNew_g[i]/w->colNew[i-1];
      }
    }

    /* If alpha is non-constant, get new values */
    if (alpha_func == 1) {
#ifdef TESTING_MODE
      user_start_t = clock();
#endif
      userAlpha(t+dt, dt, grd, w->colNew, w->presNew_g+1, 
		w->eIntNew, w->gammaNew, 
		w->deltaNew, params, w->alpha_g+1);
#ifdef TESTING_MODE
      user_end_t = clock();
      *userTime += (double) (user_end_t - user_start_t) 
	/ CLOCKS_PER_SEC;
#endif
      w->alpha_g[0] = w->alpha_g[1];
      w->alpha_g[grd->nr+1] = w->alpha_g[grd->nr];
    }

    /* Subtract off energy source term from previous iteration; we
       need to do this here because we're about to update gamma */
    for (i=0; i<grd->nr; i++) {
      gsl_vector_set(w->rhs_g, i+1, 
		     gsl_vector_get(w->rhs_g, i+1) - 
		     dt*(w->gammaNew[i]-1)*w->intEnSrc[i]);
    }

    /* Get new pressure / enthalpy boundary values */
    if (ibc_func == 1) {
#ifdef TESTING_MODE
      user_start_t = clock();
#endif
      userIBC(t+dt, dt, grd, w->colNew, w->presNew_g+1, 
	      w->eIntNew, w->gammaNew, w->deltaNew, 
	      ibc_pres, ibc_enth, params, 
	      &ibc_pres_val1, &ibc_enth_val1);
#ifdef TESTING_MODE
      user_end_t = clock();
      *userTime += (double) (user_end_t - user_start_t) 
	/ CLOCKS_PER_SEC;
#endif
    }
    if (obc_func == 1) {
#ifdef TESTING_MODE
      user_start_t = clock();
#endif
      userOBC(t+dt, dt, grd, w->colNew, w->presNew_g+1, 
	      w->eIntNew, w->gammaNew, w->deltaNew, 
	      obc_pres, obc_enth, params, 
	      &obc_pres_val1, &obc_enth_val1);
#ifdef TESTING_MODE
      user_end_t = clock();
      *userTime += (double) (user_end_t - user_start_t) 
	/ CLOCKS_PER_SEC;
#endif
    }

    /* Set boundary conditions from user-specified parameters, and
       also set ghost cell elements in RHS matrix */
    applyBC(grd, w->alpha_g, 
	    ibc_pres, ibc_enth, ibc_pres_val1, ibc_enth_val1,
	    obc_pres, obc_enth, obc_pres_val1, obc_enth_val1,
	    w->presNew_g, w->hint_g, &pIn, &qIn, &pOut, &qOut);
    gsl_vector_set(w->rhs_g, 0, pIn);
    gsl_vector_set(w->rhs_g, grd->nr+1, pOut);

    /* Get new enthalpies at cell centers and extrapolated to cell edges */
    if (interpOrder == 1) {
      for (i=0; i<=grd->nr+1; i++)
	w->hintL_g[i] = w->hintR_g[i] = w->hint_g[i];
    } else if (interpOrder == 2) {
      if (grd->linear == 1) {
	for (i=1; i<=grd->nr+1; i++)
	  w->hintR_g[i-1] = w->hintL_g[i] = 
	    ((grd->r_g[i]-grd->r_h[i-1]) * w->hint_g[i-1] +
	     (grd->r_h[i-1]-grd->r_g[i-1]) * w->hint_g[i]) /
	    (grd->r_g[i] - grd->r_g[i-1]);
      } else {
	for (i=1; i<=grd->nr+1; i++)
	  w->hintR_g[i-1] = w->hintL_g[i] = 
	    (log(grd->r_g[i]/grd->r_h[i-1]) * w->hint_g[i-1] +
	     log(grd->r_h[i-1]/grd->r_g[i-1]) * w->hint_g[i]) /
	    log(grd->r_g[i]/grd->r_g[i-1]);
      }
      /* Limit the slope */
      for (i=1; i<=grd->nr+1; i++) {
	if (w->hintR_g[i-1]/w->hint_g[i-1]-1.0 > SLOPELIMIT) {
	  w->hintR_g[i-1] = w->hint_g[i-1]*(1.0+SLOPELIMIT);
	} else if (w->hintR_g[i-1]/w->hint_g[i-1] < 1.0-SLOPELIMIT) {
	  w->hintR_g[i-1] = w->hint_g[i-1]*(1.0-SLOPELIMIT);
	}
	if (w->hintL_g[i]/w->hint_g[i]-1.0 > SLOPELIMIT) {
	  w->hintL_g[i] = w->hint_g[i]*(1.0+SLOPELIMIT);
	} else if (w->hintL_g[i]/w->hint_g[i] < 1.0-SLOPELIMIT) {
	  w->hintL_g[i] = w->hint_g[i]*(1.0-SLOPELIMIT);
	}
      }
    } else if (interpOrder == 3) {
      ppmExtrap(grd->nr+2, grd->dr_g, w->hint_g, w->hintL_g, w->hintR_g, 
		w->ppmwksp_g);
    }

    /**************************************************************/
    /* Step 4b: add energy source term to RHS                     */
    /**************************************************************/

    /* If source term can vary, compute its value for this 
       iteration */
    if (intEnSrc_func == 1) {
#ifdef TESTING_MODE
      user_start_t = clock();
#endif
      userIntEnSrc(t+dt, dt, grd, w->colNew, w->presNew_g+1, 
		   w->eIntNew, w->gammaNew, w->deltaNew, params, 
		   w->intEnSrc);
#ifdef TESTING_MODE
      user_end_t = clock();
      *userTime += (double) (user_end_t - user_start_t) 
	/ CLOCKS_PER_SEC;
#endif
    }
    for (i=0; i<grd->nr; i++) {
      gsl_vector_set(w->rhs_g, i+1, 
		     gsl_vector_get(w->rhs_g, i+1) +
		     dt*(w->gammaNew[i]-1)*w->intEnSrc[i]);
    }

    /**************************************************************/
    /* Step 4c: construct the tridiagonal matrix                  */
    /**************************************************************/

    /* Upper diagonal */
    gsl_vector_set(w->ud_g, 0, -qIn);
    for (i=1; i<=grd->nr; i++) {
      if (w->fmNew_h[i] > 0) h_up = w->hintR_g[i] + grd->psiEff_h[i];
      else h_up = w->hintL_g[i+1] + grd->psiEff_h[i];
      gsl_vector_set(w->ud_g, i,
		     dt*w->alpha_g[i+1]*(w->gammaNew[i-1]-1.0)
		     / grd->area[i-1] *
		     (grd->g_h[i]*SQR(grd->r_g[i+1]) *
		      (1-grd->beta_g[i+1]) *
		      (grd->psiEff_g[i]
		       + w->deltaNew[i-1] * 
		       w->presNew_g[i]/w->colNew[i-1] 
		       - h_up) +
		      M_PI*grd->r_h[i]*grd->vphi_h[i] *
		      (1-grd->beta_h[i])));
    }

    /* Lower diagonal */
    for (i=0; i<grd->nr; i++) {
      if (w->fmNew_h[i] > 0) h_up = w->hintR_g[i] + grd->psiEff_h[i];
      else h_up = w->hintL_g[i+1] + grd->psiEff_h[i];
      gsl_vector_set(w->ld_g, i,
		     dt*w->alpha_g[i]*(w->gammaNew[i]-1.0)
		     / grd->area[i] *
		     (grd->g_h[i]*SQR(grd->r_g[i]) *
		      (1-grd->beta_g[i]) *
		      (grd->psiEff_g[i+1] 
		       + w->deltaNew[i] * 
		       w->presNew_g[i+1]/w->colNew[i]
		       - h_up) -
		      M_PI*grd->r_h[i]*grd->vphi_h[i] *
		      (1-grd->beta_h[i])));
    }
    gsl_vector_set(w->ld_g, grd->nr, -qOut);

    /* Diagonal */
    for (i=1; i<=grd->nr; i++) {
      if (w->fmNew_h[i-1] > 0) h_up = w->hintR_g[i-1] + grd->psiEff_h[i-1];
      else h_up = w->hintL_g[i] + grd->psiEff_h[i-1];
      if (w->fmNew_h[i] > 0) h_up1 = w->hintR_g[i] + grd->psiEff_h[i];
      else h_up1 = w->hintL_g[i+1] + grd->psiEff_h[i];
      gsl_vector_set(w->diag_g, i,
		     1.0 + dt*w->alpha_g[i]*(w->gammaNew[i-1]-1.0) 
		     / grd->area[i-1] *
		     ((grd->g_h[i]*
		       (h_up1 - grd->psiEff_g[i]
			- w->deltaNew[i-1] * 
			w->presNew_g[i]/w->colNew[i-1]) +
		       grd->g_h[i-1]*
		       (h_up - grd->psiEff_g[i]
			- w->deltaNew[i-1] * 
			w->presNew_g[i]/w->colNew[i-1])) *
		      SQR(grd->r_g[i]) * (1-grd->beta_g[i]) +
		      M_PI*(grd->r_h[i]*grd->vphi_h[i]*
			    (1-grd->beta_h[i]) -
			    grd->r_h[i-1]*grd->vphi_h[i-1]*
			    (1-grd->beta_h[i-1]))));
    }

    /**************************************************************/
    /* Step 4d: solve the matrix equation to get the new pressure */
    /* guess                                                      */
    /**************************************************************/

    if (gsl_linalg_solve_tridiag(w->diag_g, w->ud_g, w->ld_g, 
				 w->rhs_g, w->presTmp_g)) {
      /* Non-zero return code means error has occured, probably due to
	 bad inputs data. Bail out, and let the driver try again */
      if (verbose > 0) 
	printf("   ...iteration %ld, detected NaN or Inf in tridiagonal inputs, exiting without convergence!\n", *itCount+1);
#ifdef TESTING_MODE
      end_t = clock();
      *advanceTime += (double) (end_t - start_t) / CLOCKS_PER_SEC;
#endif
      return(-1.0);
    }
    
    /**************************************************************/
    /* Step 4e: get the new column density guess                  */
    /**************************************************************/

    /* Get new mass fluxes */
    for (i=0; i<=grd->nr; i++)
      w->fmNew_h[i] = -grd->g_h[i] * 
	(SQR(grd->r_g[i+1])*
	 (1-grd->beta_g[i+1]) *
	 gsl_vector_get(w->presTmp_g, i+1)*w->alpha_g[i+1] - 
	 SQR(grd->r_g[i])*
	 (1-grd->beta_g[i]) * 
	 gsl_vector_get(w->presTmp_g, i)*w->alpha_g[i]);

    /* Get source terms at new time */
    if (massSrc_func == 1) {
#ifdef TESTING_MODE
      user_start_t = clock();
#endif
      userMassSrc(t+dt, dt, grd, w->colNew, w->presNew_g+1, 
		  w->eIntNew, w->gammaNew, 
		  w->deltaNew, params, w->massSrcNew);
#ifdef TESTING_MODE
      user_end_t = clock();
      *userTime += (double) (user_end_t - user_start_t) 
	/ CLOCKS_PER_SEC;
#endif
    } else {
      for (i=0; i<grd->nr; i++) w->massSrcNew[i] = massSrc_val;
    }

    /* Get new column density guess */
    for (i=0; i<grd->nr; i++) {
      w->colTmp[i] = col[i] - dt/grd->area[i] * 
	(w->fmNew_h[i+1] - w->fmNew_h[i]) +
	dt * w->massSrcNew[i];
    }

    /**************************************************************/
    /* Step 4f: if doing separate internal energy evolution, get  */
    /* new internal energy guess                                  */
    /**************************************************************/

    if (eos_func == 1) {
      for (i=0; i<grd->nr; i++) {
	w->eIntTmp[i] = eInt[i] + 1.0/(w->gammaNew[i]-1.0) *
	  (gsl_vector_get(w->presTmp_g, i+1) - pres[i]) +
	  w->deltaNew[i]*gsl_vector_get(w->presTmp_g, i+1) / 
	  w->colTmp[i] * (w->colTmp[i] - col[i]);
      }
    }

    /**************************************************************/
    /* Step 4g: check for convergence or loop termination; print  */
    /* status message                                             */
    /**************************************************************/

    /* Compute residuals on pressure */
    err = 0.0;
    maxErrCell = -1;
    residType = 0;
    for (i=0; i<grd->nr; i++) {
      if (!isfinite(gsl_vector_get(w->presTmp_g, i+1))) {
	/* If we've diverged to the point of getting inf's, we won't
	   converge, so bail out */
	if (verbose > 0) 
	  printf("   ...iteration %ld, detected NaN or Inf in tridiagonal solution, cell %ld, exiting without convergence!\n", *itCount, i);
#ifdef TESTING_MODE
	end_t = clock();
	*advanceTime += (double) (end_t - start_t) / CLOCKS_PER_SEC;
#endif
	return(-1.0);
      }
      errCell = (gsl_vector_get(w->presTmp_g, i+1) - 
		 w->presNew_g[i+1]) /
	gsl_vector_get(w->presTmp_g, i+1);
      if (fabs(errCell) > fabs(err)) {
	maxErrCell = i;
	err = errCell;
      }
    }

    /* Compute residuals on column density */
    for (i=0; i<grd->nr; i++) {
      errCell = (w->colTmp[i] - w->colNew[i]) / w->colTmp[i];
      if (fabs(errCell) > fabs(err)) {
	maxErrCell = i;
	err = errCell;
	residType = 1;
      }
    }

    /* Compute residuals on internal energy */
    for (i=0; i<grd->nr; i++) {
      errCell = (w->eIntTmp[i] - w->eIntNew[i]) / w->eIntTmp[i];
      if (fabs(errCell) > fabs(err)) {
	maxErrCell = i;
	err = errCell;
	residType = 2;
      }
    }

#ifdef TESTING_MODE
    /* In testing mode, store residuals and residual type */
    resid[*itCount-1] += err;
    rtype[*itCount-1] = residType;
#endif

    /* Print status, and bail out if we've exceed the iteration limit */
    if (fabs(err) < errTol) converged = 1;
    if (verbose > 0) {
      printf("   ...iteration %ld, max residual = %e in cell %ld",
	     *itCount, err, maxErrCell);
      if (residType==0) 
	printf(" pressure\n");
      else if (residType==1)
	printf(" density\n");
      else if (residType==2)
	printf(" internal energy\n");
      if (converged == 1) {
	printf("   ...converged!\n");
      } else if (*itCount == maxIter) {
	printf("   ...reached maximum iterations, exiting without convergence! Max residual = %e in cell %ld\n",
	       err, maxErrCell);
      }
    }

    /* If we have converged, shift temporaries into final arrays, then
       exit the loop */
    if (converged) {
      for (i=0; i<grd->nr; i++) {
	w->colNew[i] = w->colTmp[i];
	w->presNew_g[i+1] = gsl_vector_get(w->presTmp_g, i+1);
	if (eos_func == 1)
	  w->eIntNew[i] = w->eIntTmp[i];
      }
      w->presNew_g[0] = gsl_vector_get(w->presTmp_g, 0);
      w->presNew_g[grd->nr+1] = gsl_vector_get(w->presTmp_g, grd->nr+1);
      break;
    }

    /* If we've exceeded iteration limit, bail out */
    if (*itCount == maxIter) {
#ifdef TESTING_MODE
      end_t = clock();
      *advanceTime += (double) (end_t - start_t) / CLOCKS_PER_SEC;
#endif
      return(-1.0);
    }

    /**************************************************************/
    /* Step 4h: if we're here, we haven't converged yet, so use   */
    /* Anderson acceleration to generate new guesses at solution  */
    /**************************************************************/

#ifdef TESTING_MODE
    next_iter_start_t = clock();
#endif
    getNextIterate(grd, eos_func, w
#if AA_M > 0
		   , *itCount, nHist, &residMax
#endif
		   );
#ifdef TESTING_MODE
    next_iter_end_t = clock();
    *nextIterTime += (double) (next_iter_end_t - next_iter_start_t) /
      CLOCKS_PER_SEC;
#endif
    
  } /* End main iteration loop */


  /****************************************************************/
  /* Step 5: compute rate of change of pressure and column        */
  /* density; use this to estimate new time step                  */
  /****************************************************************/

  dtMin = LARGE;
  for (i=0; i<grd->nr; i++) {
    dCol = fabs(w->colNew[i] - col[i]) + SMALL;
    dPres = fabs(w->presNew_g[i+1] - pres[i]) + SMALL;
    dtMin = fmin( dtMin, dtTol*dt * 
		  fmin(col[i]/dCol, pres[i]/dPres) );
    if (eos_func == 1) {
      dEInt = fabs(w->eIntNew[i] - eInt[i]);
      dtMin = fmin( dtMin, dtTol*dt * eInt[i]/dEInt );
    }
  }
  if (verbose > 0)
    printf("advanceBE: estimated dtNew = %e\n", dtMin);

  /****************************************************************/
  /* Step 6: record the total mass and energy transported across  */
  /* the inner and outer boundaries during this time step.        */
  /****************************************************************/

  /* Mass fluxes */
  mBndTemp[0] += dt*w->fmNew_h[0];
  mBndTemp[1] += dt*w->fmNew_h[grd->nr];

  /* New time energy fluxes; need to get enthalpy to evaluate f_E */
  if (w->fmNew_h[0] > 0) h_up = w->hintR_g[0] + grd->psiEff_h[0];
  else h_up = w->hintL_g[1] + grd->psiEff_h[0];
  eBndTemp[0] += dt*h_up*w->fmNew_h[0];
  if (w->fmNew_h[grd->nr] > 0) 
    h_up = w->hintR_g[grd->nr] + grd->psiEff_h[grd->nr];
  else h_up = w->hintL_g[grd->nr+1] + grd->psiEff_h[grd->nr];
  eBndTemp[1] += dt*h_up*w->fmNew_h[grd->nr];

  /* New time torque fluxes */
  eBndTemp[0] += dt*M_PI*grd->r_h[0]*grd->vphi_h[0]*(1.0-grd->beta_h[0]) *
    (w->presNew_g[1]*w->alpha_g[1] + w->presNew_g[0]*w->alpha_g[0]);
  eBndTemp[1] += dt*M_PI*grd->r_h[grd->nr]*grd->vphi_h[grd->nr]
    *(1.0-grd->beta_h[grd->nr]) *
    (w->presNew_g[grd->nr+1]*w->alpha_g[grd->nr+1] 
     + w->presNew_g[grd->nr]*w->alpha_g[grd->nr]);

  /****************************************************************/
  /* Step 7: if using source functions, record total mass and     */
  /* energy added to every cell by sources                        */
  /****************************************************************/

  if (massSrc_func == 1) {
    for (i=0; i<grd->nr; i++) {
      w->mSrc[i] += dt * w->massSrcNew[i];
      w->eSrc[i] += dt * w->massSrcNew[i] *
	(grd->psiEff_g[i+1] + 
	 w->deltaNew[i]*w->presNew_g[i+1]/w->colNew[i]);
    }
  }
  if (intEnSrc_func == 1) {
    for (i=0; i<grd->nr; i++)
      w->eSrc[i] += dt * w->intEnSrc[i];
  }

  /****************************************************************/
  /* Step 8: store final values in arrays                         */
  /****************************************************************/

  if (noUpdate==0) {
    for (i=0; i<grd->nr; i++) {
      pres[i] = w->presNew_g[i+1];
      col[i] = w->colNew[i];
    }
    if (eos_func==1)
      for (i=0; i<grd->nr; i++) eInt[i] = w->eIntNew[i];
    mBnd[0] = mBndTemp[0];
    mBnd[1] = mBndTemp[1];
    eBnd[0] = eBndTemp[0];
    eBnd[1] = eBndTemp[1];
    if (massSrc_func == 1) {
      for (i=0; i<grd->nr; i++) {
	mSrc[i] = w->mSrc[i];
	eSrc[i] = w->eSrc[i];
      }
    } else if (intEnSrc_func == 1) {
      for (i=0; i<grd->nr; i++) eSrc[i] = w->eSrc[i];
    }
  }

  /****************************************************************/
  /* Step 8: return                                               */
  /****************************************************************/

#ifdef TESTING_MODE
  end_t = clock();
  *advanceTime += (double) (end_t - start_t) / CLOCKS_PER_SEC;
#endif
  return(dtMin);
}
