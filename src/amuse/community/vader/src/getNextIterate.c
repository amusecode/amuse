#include <math.h>
#include <gsl/gsl_linalg.h>
#include "getNextIterate.h"

void getNextIterate(const grid *grd, const bool eos_func, 
		    const wksp *w
#if AA_M > 0 
		    , const unsigned long itCount, const long nHist,
		    double *residMax
#endif
		    ) {

  long i;
#if AA_M > 0
  long j, mk, colptr;
  gsl_matrix *resid;
  gsl_vector *tau, *llsresid, *aa_wgts, *constraint;
  gsl_vector_view constraint_view;
#endif

#if AA_M == 0
/* No Anderson acceleration, so just move temporaries into next
   guess */
  for (i=0; i<grd->nr; i++) {
    w->colNew[i] = w->colTmp[i];
    w->presNew_g[i+1] = gsl_vector_get(w->presTmp_g, i+1);
    if (eos_func==1)
      w->eIntNew[i] = w->eIntTmp[i];
  }
#else

  /* Anderson acceleration code */

  /* Set up pointers */
  mk = itCount < AA_M+1 ? itCount : AA_M+1;
  colptr = (itCount - 1) % (AA_M+1);

  /* Store the current guesses in the history array */
  for (i=0; i<grd->nr; i++) {
    w->colHist[i + colptr*grd->nr] = w->colTmp[i];
    w->presHist[i + colptr*grd->nr] = gsl_vector_get(w->presTmp_g, i+1);
    if (eos_func)
      w->eIntHist[i + colptr*grd->nr] = w->eIntTmp[i];
  }

  /* Compute residuals and store them in the residual matrix */
  for (i=0; i<grd->nr; i++) {
    w->colResid[i + colptr*grd->nr] = 
      (w->colHist[i+colptr*grd->nr]-w->colNew[i]) / 
      w->colHist[i+colptr*grd->nr];
    if (fabs(w->colResid[i + colptr*grd->nr]) > *residMax) {
      *residMax = fabs(w->colResid[i + colptr*grd->nr]);
    }
    w->presResid[i + colptr*grd->nr] = 
      (w->presHist[i+colptr*grd->nr]-w->presNew_g[i+1]) / 
      w->presHist[i+colptr*grd->nr];
    if (fabs(w->presResid[i + colptr*grd->nr]) > *residMax) {
      *residMax = fabs(w->presResid[i + colptr*grd->nr]);
    }
    if (eos_func) {
      w->eIntResid[i + colptr*grd->nr] = 
	(w->eIntHist[i+colptr*grd->nr]-w->eIntNew[i]) / 
	w->eIntHist[i+colptr*grd->nr];
      if (fabs(w->eIntResid[i + colptr*grd->nr]) > *residMax)
	*residMax = fabs(w->eIntResid[i + colptr*grd->nr]);
    }
  }

  /* If we're in the first iteration, set move temporaries to next
     guess */
  if (itCount == 1) {
    for (i=0; i<grd->nr; i++) {
      w->colNew[i] = w->colHist[i+colptr*grd->nr];
      w->presNew_g[i+1] = w->presHist[i+colptr*grd->nr];
      if (eos_func)
	w->eIntNew[i] = w->eIntHist[i+colptr*grd->nr];
    }
  } else {
    /* On subsequent iterations, solve the linear least squares
       optimization problem to get a weighted sum of previous
       guesses for the new guess */

    /* Allocate workspace for this iteration */
    resid = gsl_matrix_alloc(nHist+1, mk);
    llsresid = gsl_vector_alloc(nHist+1);
    aa_wgts = gsl_vector_alloc(mk);
    if (nHist+1 < mk) {
      tau = gsl_vector_alloc(nHist+1);
    } else {
      tau = gsl_vector_alloc(mk);
    }

    /* Build matrix of residuals */
    for (i=0; i<grd->nr; i++) {
      for (j=0; j<mk; j++) {
	gsl_matrix_set(resid, i, j, w->colResid[i+j*grd->nr]);
	gsl_matrix_set(resid, i+grd->nr, j, w->presResid[i+j*grd->nr]);
	if (eos_func)
	  gsl_matrix_set(resid, i+2*grd->nr, j, w->eIntResid[i+j*grd->nr]);
      }
    }

    /* Set up constraints on coefficients */
    for (j=0; j<mk; j++) gsl_matrix_set(resid, nHist, j, 1.0e3*(*residMax));
    constraint_view = gsl_vector_subvector(w->constraint, 0, nHist+1);
    constraint = &(constraint_view.vector);
    gsl_vector_set(constraint, nHist, 1e3*(*residMax));

    /* Solve linear least squares system by QR decomposition */
    gsl_linalg_QR_decomp(resid, tau);
    gsl_linalg_QR_lssolve(resid, tau, constraint, aa_wgts, llsresid);

    /* Build the new guess from the computed weights */
    for (i=0; i<grd->nr; i++) {
      w->colNew[i] = 
	gsl_vector_get(aa_wgts, 0)*w->colHist[i];
      w->presNew_g[i+1] = 
	gsl_vector_get(aa_wgts, 0)*w->presHist[i];
      if (eos_func)
	w->eIntNew[i] = 
	  gsl_vector_get(aa_wgts, 0)*w->eIntHist[i];
    }
    for (j=1; j<mk; j++) {
      for (i=0; i<grd->nr; i++) {
	w->colNew[i] += gsl_vector_get(aa_wgts, j) *
	  w->colHist[i+j*grd->nr];
	w->presNew_g[i+1] += gsl_vector_get(aa_wgts, j) *
	  w->presHist[i+j*grd->nr];
	if (eos_func)
	  w->eIntNew[i] += gsl_vector_get(aa_wgts, j) *
	    w->eIntHist[i+j*grd->nr];
      }
    }

    /* Free memory for this iteration */
    gsl_vector_free(llsresid);
    gsl_vector_free(tau);
    gsl_vector_free(aa_wgts);
    gsl_matrix_free(resid);
  }

#endif
  /* End Anderson acceleration code */
}
