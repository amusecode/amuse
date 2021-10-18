#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "init.h"
#include "rotCurveSpline.h"

#define G 6.67384e-8   /* 2010 CODATA value */

grid *
gridAlloc(const unsigned long nr) {
  grid *grd;

  /* Allocate memory */
  if (!(grd = malloc(sizeof(grid)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  grd->nr = nr;
  if (!(grd->r_h = calloc(nr+1, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  if (!(grd->r_g = calloc(nr+2, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  if (!(grd->dr_g = calloc(nr+2, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  if (!(grd->area = calloc(nr, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  if (!(grd->vphi_g = calloc(nr+2, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  if (!(grd->vphi_h = calloc(nr+1, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  if (!(grd->beta_g = calloc(nr+2, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  if (!(grd->beta_h = calloc(nr+1, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  if (!(grd->psiEff_g = calloc(nr+2, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  if (!(grd->g_h = calloc(nr+1, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  if (!(grd->psiEff_h = calloc(nr+1, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }

  return(grd);
}


grid *
gridInit(const unsigned long nr, 
	 const double *r_g, const double *r_h, 
	 const double *vphi_g, const double *vphi_h, 
	 const double *beta_g, const double *beta_h,
	 const double *psiEff_g, const double *psiEff_h,
	 const double *g_h,
	 const bool linear) {

  int i;
  grid *grd;

  /* Allocate */
  grd = gridAlloc(nr);

  /* Initialize */
  grd->nr = nr;
  grd->linear = linear;
  for (i=0; i<nr; i++) {
    grd->r_g[i] = r_g[i];
    grd->r_h[i] = r_h[i];
    grd->vphi_g[i] = vphi_g[i];
    grd->vphi_h[i] = vphi_h[i];
    grd->beta_g[i] = beta_g[i];
    grd->beta_h[i] = beta_h[i];
    grd->psiEff_g[i] = psiEff_g[i];
    grd->psiEff_h[i] = psiEff_h[i];
    grd->g_h[i] = g_h[i];
  }
  grd->r_g[nr] = r_g[nr];
  grd->r_h[nr] = r_h[nr];
  grd->vphi_g[nr] = vphi_g[nr];
  grd->vphi_h[nr] = vphi_h[nr];
  grd->beta_g[nr] = beta_g[nr];
  grd->beta_h[nr] = beta_h[nr];
  grd->psiEff_g[nr] = psiEff_g[nr];
  grd->psiEff_h[nr] = psiEff_h[nr];
  grd->g_h[nr] = g_h[nr];
  grd->r_g[nr+1] = r_g[nr+1];
  grd->vphi_g[nr+1] = vphi_g[nr+1];
  grd->beta_g[nr+1] = beta_g[nr+1];
  grd->psiEff_g[nr+1] = psiEff_g[nr+1];
  for (i=0; i<nr; i++)
    grd->area[i] = M_PI*(SQR(grd->r_h[i+1]) - SQR(grd->r_h[i]));
  if (linear) {
    for (i=1; i<=grd->nr; i++) 
      grd->dr_g[i] = grd->r_h[i] - grd->r_h[i-1];
  } else {
    for (i=1; i<=grd->nr; i++) 
      grd->dr_g[i] = log(grd->r_h[i]/grd->r_h[i-1]);
  }
  grd->dr_g[0] = grd->dr_g[1];
  grd->dr_g[nr+1] = grd->dr_g[nr];

  return(grd);
}


grid *
gridInitTabulated(const unsigned long nr, const double *rTab, 
		  const double *vTab, const int nTab,
		  const double rmin, const double rmax,
		  const unsigned int bspline_degree,
		  const unsigned int bspline_breakpoints,
		  const bool linear) {

  int i;
  grid *grd;
  double dx;
  double *rtmp, *vphi, *psi, *beta;

  /* Allocate */
  grd = (grid *) gridAlloc(nr);

  /* Initialize */

  /* Grid structure */
  grd->nr = nr;
  grd->linear = linear;
  if (linear) {
    dx = (rmax-rmin) / nr;
    for (i=0; i<=nr; i++)
      grd->r_h[i] = rmin + dx*i;      
    for (i=0; i<nr+2; i++) {
      grd->r_g[i] = rmin + (i-0.5)*dx;
    }
  } else {
    dx = (log(rmax)-log(rmin)) / nr;
    for (i=0; i<=nr; i++)
      grd->r_h[i] = exp(log(rmin) + dx*i);      
    for (i=0; i<nr+2; i++) {
      grd->r_g[i] = exp(log(rmin) + (i-0.5)*dx);
    }
  }
  grd->dr_g[0] = grd->dr_g[1];
  grd->dr_g[nr+1] = grd->dr_g[nr];

  /* Safety check: make sure the table covers the range of the data */
  if ((grd->r_g[0] < rTab[0]) || (grd->r_g[nr+1] > rTab[nTab-1])) {
    gridFree(grd);
    fprintf(stderr, "Error: table too small to cover required range r = %e to %e\n", grd->r_g[0], grd->r_g[1]);
    return(NULL);
  }

  /* Area */
  for (i=0; i<nr; i++)
    grd->area[i] = M_PI*(SQR(grd->r_h[i+1]) - SQR(grd->r_h[i]));
  if (linear) {
    for (i=1; i<=grd->nr; i++) 
      grd->dr_g[i] = grd->r_h[i] - grd->r_h[i-1];
  } else {
    for (i=1; i<=grd->nr; i++) 
      grd->dr_g[i] = log(grd->r_h[i]/grd->r_h[i-1]);
  }
  grd->dr_g[0] = grd->dr_g[1];
  grd->dr_g[nr+1] = grd->dr_g[nr];

  /* Load temporary array of r values */
  if (!(rtmp = calloc(2*nr+3, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  for (i=0; i<=nr+1; i++) rtmp[2*i] = grd->r_g[i];
  for (i=0; i<=nr; i++) rtmp[2*i+1] = grd->r_h[i];

  /* Create temporary vphi, psi arrays */
  if (!(vphi = calloc(2*nr+3, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  if (!(psi = calloc(2*nr+3, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  if (!(beta = calloc(2*nr+3, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }

  /* Call interpolation routine */
  rotCurveSpline(rTab, vTab, nTab, bspline_degree, bspline_breakpoints,
		 rtmp, 2*nr+3, vphi, psi, beta);

  /* Unpack results into grid structure */
  for (i=0; i<=nr+1; i++) {
    grd->vphi_g[i] = vphi[2*i];
    grd->beta_g[i] = beta[2*i];
    grd->psiEff_g[i] = psi[2*i] + SQR(vphi[2*i])/2.0;
  }
  for (i=0; i<=nr; i++) {
    grd->vphi_h[i] = vphi[2*i+1];
    grd->beta_h[i] = beta[2*i+1];
    grd->psiEff_h[i] = psi[2*i+1] + SQR(vphi[2*i+1])/2.0;
  }

  /* Compute g_h */
  if (linear) {
    for (i=0; i<nr+1; i++)
      grd->g_h[i] = 2.0*M_PI / (grd->vphi_h[i]*(1+grd->beta_h[i])) *
	1.0/dx;
  } else {
    for (i=0; i<nr+1; i++)
      grd->g_h[i] = 2.0*M_PI / (grd->vphi_h[i]*(1+grd->beta_h[i])) *
	1.0/(dx*grd->r_h[i]);
  }

  /* Free temporary storage */
  free(rtmp);
  free(vphi);
  free(psi);
  free(beta);

  /* Return */
  return(grd);
}


grid *
gridInitTabulatedNonUniform(const unsigned long nr, const double *rTab, 
			    const double *vTab, const int nTab,
			    const double *r_g, const double *r_h,
			    const bool linear) {

  int i;
  grid *grd;
  double *rtmp, *vphi, *psi, *beta;

  /* Allocate */
  grd = (grid *) gridAlloc(nr);

  /* Initialize */

  /* Grid structure */
  grd->nr = nr;
  grd->linear = linear;
  for (i=0; i<nr+1; i++) {
    grd->r_h[i] = r_h[i];
    grd->r_g[i] = r_g[i];
  }
  grd->r_g[nr+1] = r_g[nr+1];

  /* Safety check: make sure the table covers the range of the data */
  if ((grd->r_g[0] < rTab[0]) || (grd->r_g[nr+1] > rTab[nTab-1])) {
    gridFree(grd);
    fprintf(stderr, "Error: table too small to cover required range r = %e to %e\n", grd->r_g[0], grd->r_g[1]);
    return(NULL);
  }

  /* Area */
  for (i=0; i<nr; i++)
    grd->area[i] = M_PI*(SQR(grd->r_h[i+1]) - SQR(grd->r_h[i]));
  if (linear) {
    for (i=1; i<=grd->nr; i++) 
      grd->dr_g[i] = grd->r_h[i] - grd->r_h[i-1];
  } else {
    for (i=1; i<=grd->nr; i++) 
      grd->dr_g[i] = log(grd->r_h[i]/grd->r_h[i-1]);
  }
  grd->dr_g[0] = grd->dr_g[1];
  grd->dr_g[nr+1] = grd->dr_g[nr];

  /* Load temporary array of r values */
  if (!(rtmp = calloc(2*nr+3, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  for (i=0; i<=nr+1; i++) rtmp[2*i] = grd->r_g[i];
  for (i=0; i<=nr; i++) rtmp[2*i+1] = grd->r_h[i];

  /* Create temporary vphi, psi arrays */
  if (!(vphi = calloc(2*nr+3, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  if (!(psi = calloc(2*nr+3, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }
  if (!(beta = calloc(2*nr+3, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new grid!\n");
    return(NULL);
  }

  /* Call interpolation routine */
  rotCurveSpline(rTab, vTab, nTab, 6, 15, rtmp, nr, vphi, psi, beta);

  /* Unpack results into grid structure */
  for (i=0; i<=nr+1; i++) {
    grd->vphi_g[i] = vphi[2*i];
    grd->beta_g[i] = beta[2*i];
    grd->psiEff_g[i] = psi[2*i] + SQR(vphi[2*i])/2.0;
  }
  for (i=0; i<=nr; i++) {
    grd->vphi_h[i] = vphi[2*i+1];
    grd->beta_h[i] = beta[2*i+1];
    grd->psiEff_h[i] = psi[2*i+1] + SQR(vphi[2*i]+1)/2.0;
  }

  /* Compute g_h */
  if (linear) {
    for (i=0; i<nr+1; i++)
      grd->g_h[i] = 2.0*M_PI / (grd->vphi_h[i]*(1+grd->beta_h[i])) *
	1.0/(grd->r_g[i+1]-grd->r_g[i]);
  } else {
    for (i=0; i<nr+1; i++)
      grd->g_h[i] = 2.0*M_PI / (grd->vphi_h[i]*(1+grd->beta_h[i])) *
	1.0/(log(grd->r_g[i+1]/grd->r_g[i])*grd->r_h[i]);
  }

  /* Free temporary storage */
  free(rtmp);
  free(vphi);
  free(psi);
  free(beta);

  /* Return */
  return(grd);
}


grid *
gridInitKeplerian(const unsigned long nr, const double rmin, 
		  const double rmax, const double m, 
		  const bool linear) {

  int i;
  grid *grd;
  double dx;

  /* Allocate */
  grd = (grid *) gridAlloc(nr);

  /* Initialize */

  /* Grid structure */
  grd->nr = nr;
  grd->linear = linear;
  if (linear) {
    dx = (rmax-rmin) / nr;
    for (i=0; i<=nr; i++)
      grd->r_h[i] = rmin + dx*i;      
    for (i=0; i<nr+2; i++) {
      grd->r_g[i] = rmin + (i-0.5)*dx;
    }
  } else {
    dx = (log(rmax)-log(rmin)) / nr;
    for (i=0; i<=nr; i++)
      grd->r_h[i] = exp(log(rmin) + dx*i);      
    for (i=0; i<nr+2; i++) {
      grd->r_g[i] = exp(log(rmin) + (i-0.5)*dx);
    }
  }

  /* Other stuff */
  for (i=0; i<nr; i++) {
    grd->area[i] = M_PI*(SQR(grd->r_h[i+1]) - SQR(grd->r_h[i]));
  }
  if (linear) {
    for (i=1; i<=grd->nr; i++) 
      grd->dr_g[i] = grd->r_h[i] - grd->r_h[i-1];
  } else {
    for (i=1; i<=grd->nr; i++) 
      grd->dr_g[i] = log(grd->r_h[i]/grd->r_h[i-1]);
  }
  grd->dr_g[0] = grd->dr_g[1];
  grd->dr_g[nr+1] = grd->dr_g[nr];
  for (i=0; i<nr+1; i++) {
    grd->vphi_h[i] = sqrt(G*m/grd->r_h[i]);
    grd->beta_h[i] = -0.5;
    grd->psiEff_h[i] = -G*m/(2.0*grd->r_h[i]);
  }
  for (i=0; i<nr+2; i++) {
    grd->vphi_g[i] = sqrt(G*m/grd->r_g[i]);
    grd->beta_g[i] = -0.5;
    grd->psiEff_g[i] = -G*m/(2.0*grd->r_g[i]);
  }
  if (linear) {
    for (i=0; i<nr+1; i++)
      grd->g_h[i] = 2.0*M_PI / (grd->vphi_h[i]*(1+grd->beta_h[i])) *
	1.0/dx;
  } else {
    for (i=0; i<nr+1; i++)
      grd->g_h[i] = 2.0*M_PI / (grd->vphi_h[i]*(1+grd->beta_h[i])) *
	1.0/(dx*grd->r_h[i]);
  }

  return(grd);
}


grid *
gridInitFlat(const unsigned long nr, const double rmin, 
	     const double rmax, const double vphi, 
	     const bool linear) {

  int i;
  grid *grd;
  double dx;

  /* Allocate */
  grd = (grid *) gridAlloc(nr);

  /* Initialize */

  /* Grid structure */
  grd->nr = nr;
  grd->linear = linear;
  if (linear) {
    dx = (rmax-rmin) / nr;
    for (i=0; i<=nr; i++)
      grd->r_h[i] = rmin + dx*i;      
    for (i=0; i<nr+2; i++) {
      grd->r_g[i] = rmin + (i-0.5)*dx;
    }
  } else {
    dx = log(rmax/rmin) / nr;
    for (i=0; i<=nr; i++)
      grd->r_h[i] = exp(log(rmin) + dx*i);      
    for (i=0; i<nr+2; i++) {
      grd->r_g[i] = exp(log(rmin) + (i-0.5)*dx);
    }
  }

  /* Other stuff */
  for (i=0; i<nr; i++) {
    grd->area[i] = M_PI*(SQR(grd->r_h[i+1]) - SQR(grd->r_h[i]));
  }
  for (i=0; i<nr+1; i++) {
    grd->vphi_h[i] = vphi;
    grd->beta_h[i] = 0.0;
    grd->psiEff_h[i] = (0.5+log(grd->r_h[i]/grd->r_h[nr]))*SQR(vphi);
  }
  if (linear) {
    for (i=1; i<=grd->nr; i++) 
      grd->dr_g[i] = grd->r_h[i] - grd->r_h[i-1];
  } else {
    for (i=1; i<=grd->nr; i++) 
      grd->dr_g[i] = log(grd->r_h[i]/grd->r_h[i-1]);
  }
  grd->dr_g[0] = grd->dr_g[1];
  grd->dr_g[nr+1] = grd->dr_g[nr];
  for (i=0; i<nr+2; i++) {
    grd->vphi_g[i] = vphi;
    grd->beta_g[i] = 0.0;
    grd->psiEff_g[i] = (0.5+log(grd->r_g[i]/grd->r_h[nr]))*SQR(vphi);
  }
  if (linear) {
    for (i=0; i<nr+1; i++)
      grd->g_h[i] = 2.0*M_PI / (grd->vphi_h[i]*(1+grd->beta_h[i])) *
	1.0/dx;
  } else {
    for (i=0; i<nr+1; i++)
      grd->g_h[i] = 2.0*M_PI / (grd->vphi_h[i]*(1+grd->beta_h[i])) *
	1.0/(dx*grd->r_h[i]);
  }

  return(grd);
}


bool outputAlloc(const unsigned long nOut, const bool eos_func,
		 const bool massSrc_func, const bool intEnSrc_func,
		 const unsigned long nUserOut, const grid *grd,
		 double **tOut, double **colOut, double **presOut,
		 double **eIntOut, double **mBndOut, double **eBndOut,
		 double **mSrcOut, double **eSrcOut, double **userOut) {
  *tOut = (double *) calloc(nOut, sizeof(double));
  *colOut = (double *) calloc(nOut*grd->nr, sizeof(double));
  *presOut = (double *) calloc(nOut*grd->nr, sizeof(double));
  if (eos_func)
    *eIntOut = (double *) calloc(nOut*grd->nr, sizeof(double));
  else if (eIntOut)
    *eIntOut = NULL;
  *mBndOut = (double *) calloc(nOut*2, sizeof(double));
  *eBndOut = (double *) calloc(nOut*2, sizeof(double));
  if (massSrc_func)
    *mSrcOut = (double *) calloc(nOut*grd->nr, sizeof(double));
  else if (mSrcOut)
    *mSrcOut = NULL;
  if (massSrc_func || intEnSrc_func)
    *eSrcOut = (double *) calloc(nOut*grd->nr, sizeof(double));
  else if (eSrcOut)
    *eSrcOut = NULL;
  if (nUserOut > 0)
    *userOut = (double *) calloc(nOut*grd->nr*nUserOut, sizeof(double));
  else if (userOut)
    *userOut = NULL;
  if (!*tOut || !*colOut || !*presOut || !*mBndOut || !*eBndOut ||
      (eos_func && !*eIntOut) || (massSrc_func && !*mSrcOut) ||
      ((massSrc_func || intEnSrc_func) && !*eSrcOut) ||
      ((nUserOut > 0) && !*userOut)) {
    fprintf(stderr, "vader: outputAlloc: cannot allocate memory to hold outputs!");
    outputFree(tOut, colOut, presOut, eIntOut, mBndOut, eBndOut,
	       mSrcOut, eSrcOut, userOut);
    return false;
  }
  return true;
}

void outputFree(double **tOut, double **colOut, double **presOut,
		double **eIntOut, double **mBndOut, double **eBndOut,
		double **mSrcOut, double **eSrcOut, double **userOut) {
  if (tOut) { free(*tOut); *tOut = NULL; }
  if (colOut) { free(*colOut); *colOut = NULL; }
  if (presOut) { free(*presOut); *presOut = NULL; }
  if (eIntOut) { free(*eIntOut); *eIntOut = NULL; }
  if (mBndOut) { free(*mBndOut); *mBndOut = NULL; }
  if (eBndOut) { free(*eBndOut); *eBndOut = NULL; }
  if (mSrcOut) { free(*mSrcOut); *mSrcOut = NULL; }
  if (eSrcOut) { free(*eSrcOut); *eSrcOut = NULL; }
  if (userOut) { free(*userOut); *userOut = NULL; }
}

bool outputResize(const unsigned long nOut, const bool eos_func,
		  const bool massSrc_func, const bool intEnSrc_func,
		  const unsigned long nUserOut, const grid *grd,
		  double **tOut, double **colOut, double **presOut,
		  double **eIntOut, double **mBndOut, double **eBndOut,
		  double **mSrcOut, double **eSrcOut, double **userOut) {
  *tOut = (double *) realloc(*tOut, nOut*sizeof(double));
  *colOut = (double *) realloc(*colOut, nOut*grd->nr*sizeof(double));
  *presOut = (double *) realloc(*presOut, nOut*grd->nr*sizeof(double));
  if (eos_func)
    *eIntOut = (double *) realloc(*eIntOut, nOut*grd->nr*sizeof(double));
  else if (eIntOut)
    *eIntOut = NULL;
  *mBndOut = (double *) realloc(*mBndOut, nOut*2*sizeof(double));
  *eBndOut = (double *) realloc(*eBndOut, nOut*2*sizeof(double));
  if (massSrc_func)
    *mSrcOut = (double *) realloc(*mSrcOut, nOut*grd->nr*sizeof(double));
  else if (mSrcOut)
    *mSrcOut = NULL;
  if (massSrc_func || intEnSrc_func)
    *eSrcOut = (double *) realloc(*eSrcOut, nOut*grd->nr*sizeof(double));
  else if (eSrcOut)
    *eSrcOut = NULL;
  if (nUserOut > 0)
    *userOut = (double *) realloc(*userOut,
				  nOut*grd->nr*nUserOut*sizeof(double));
  else if (userOut)
    *userOut = NULL;
  if (!*tOut || !*colOut || !*presOut || !*mBndOut || !*eBndOut ||
      (eos_func && !*eIntOut) || (massSrc_func && !*mSrcOut) ||
      ((massSrc_func || intEnSrc_func) && !*eSrcOut) ||
      ((nUserOut > 0) && !*userOut)) {
    fprintf(stderr, "vader: outputAlloc: cannot allocate memory to hold outputs!");
    outputFree(tOut, colOut, presOut, eIntOut, mBndOut, eBndOut,
	       mSrcOut, eSrcOut, userOut);
    return false;
  }
  return true;
}

wksp *
wkspAlloc(const unsigned long nr) {
  wksp *w;

  /* Allocate memory */
  if (!(w = malloc(sizeof(wksp)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->pres_g = calloc(nr+2, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->presNew_g = calloc(nr+2, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->colNew = calloc(nr, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->colTmp = calloc(nr, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->alpha_g = calloc(nr+2, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->hint_g = calloc(nr+2, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->hintL_g = calloc(nr+2, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->hintR_g = calloc(nr+2, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->ppmwksp_g = calloc(nr+2, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->fmLast_h = calloc(nr+1, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->fmNew_h = calloc(nr+1, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->ftLast_h = calloc(nr+1, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->feLast_h = calloc(nr+1, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->massSrcLast = calloc(nr, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->massSrcNew = calloc(nr, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->intEnSrc = calloc(nr, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->gammaLast = calloc(nr, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->deltaLast = calloc(nr, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->gammaNew = calloc(nr, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->deltaNew = calloc(nr, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->eIntTmp = calloc(nr, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->eIntNew = calloc(nr, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->mSrc = calloc(nr, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->eSrc = calloc(nr, sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->ud_g = gsl_vector_alloc(nr+1))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->ld_g = gsl_vector_alloc(nr+1))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->diag_g = gsl_vector_alloc(nr+2))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->rhs_g = gsl_vector_alloc(nr+2))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->presTmp_g = gsl_vector_alloc(nr+2))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }

  /* Since the ghost cell entries for the diagonal will always be
     unity, set them here */
  gsl_vector_set(w->diag_g, 0, 1.0);
  gsl_vector_set(w->diag_g, nr+1, 1.0);

  /* If using Anderson acceleration, allocate memory for it */
#if AA_M > 0
  if (!(w->colHist = calloc(nr*(AA_M+1), sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->presHist = calloc(nr*(AA_M+1), sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->eIntHist = calloc(nr*(AA_M+1), sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->colResid = calloc(nr*(AA_M+1), sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->presResid = calloc(nr*(AA_M+1), sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->eIntResid = calloc(nr*(AA_M+1), sizeof(double)))) {
    fprintf(stderr, "Error: unable to allocate memory for new workspace!\n");
    return(NULL);
  }
  if (!(w->constraint = gsl_vector_alloc(3*nr+1))) {
    fprintf(stderr, "Error: unable to allocate memory for workspace!\n");
    exit(1);
  }
#endif

  return(w);
}


void gridFree(grid *grd) {
  free(grd->r_h);
  free(grd->r_g);
  free(grd->dr_g);
  free(grd->area);
  free(grd->vphi_g);
  free(grd->vphi_h);
  free(grd->beta_g);
  free(grd->beta_h);
  free(grd->psiEff_g);
  free(grd->psiEff_h);
  free(grd->g_h);
  free(grd);
}

void wkspFree(wksp *w) {
  free(w->pres_g);
  free(w->presNew_g);
  free(w->colNew);
  free(w->colTmp);
  free(w->alpha_g);
  free(w->hint_g);
  free(w->hintL_g);
  free(w->hintR_g);
  free(w->ppmwksp_g);
  free(w->fmLast_h);
  free(w->fmNew_h);
  free(w->ftLast_h);
  free(w->feLast_h);
  free(w->massSrcLast);
  free(w->massSrcNew);
  free(w->intEnSrc);
  free(w->gammaLast);
  free(w->gammaNew);
  free(w->deltaLast);
  free(w->deltaNew);
  free(w->eIntTmp);
  free(w->eIntNew);
  free(w->mSrc);
  free(w->eSrc);
  gsl_vector_free(w->ud_g);
  gsl_vector_free(w->ld_g);
  gsl_vector_free(w->diag_g);
  gsl_vector_free(w->rhs_g);
  gsl_vector_free(w->presTmp_g);
#if AA_M > 0
  free(w->colHist);
  free(w->presHist);
  free(w->eIntHist);
  free(w->colResid);
  free(w->presResid);
  free(w->eIntResid);
  gsl_vector_free(w->constraint);
#endif
  free(w);
}
