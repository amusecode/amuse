#ifndef _init_h_
#define _init_h_

#include "vader_common.h"
#include <stdbool.h>

/**************************************************************************/
/* Routines to initialize and de-allocate memory in vader                 */
/*                                                                        */
/* gridAlloc : allocate a grid object                                     */
/* gridInit : allocate and initialize a grid object from specified        */
/*            values of all grid parameters                               */
/* gridInitTabulated : allocate and initialize a grid object using        */
/*             tabulated values of (r, v)                                 */
/* gridInitTabulatedNonUniform : allocate and initialize a grid object    */
/*            using tabulated values of (r, v) and a                      */
/*            non-uniformly-spaced mesh                                   */
/* gridInitKeplerian : allocate and initialize a grid object for          */
/*                     Keplerian rotation                                 */
/* gridInitFlat : allocate and initialize a grid object for a flat        */
/*                rotation curve                                          */
/* gridFree : deallocate a grid object                                    */
/* outputAlloc : allocate memory to hold vader outputs                    */
/* outputFree : de-allocate vader outputs                                 */
/* wkspAlloc : allocate a workspace object                                */
/* wkspFree : deallocate a workspace object                               */
/**************************************************************************/

grid *gridAlloc(const unsigned long nr);

grid *gridInit(const unsigned long nr,
	       const double *r_g, const double *r_h, 
	       const double *vphi_g, const double *vphi_h, 
	       const double *beta_g, const double *beta_h, 
	       const double *psiEff_g, const double *psiEff_h,
	       const double *g_h,
	       const bool linear);

grid *gridInitTabulated(const unsigned long nr, const double *rTab, 
			const double *vTab, const int nTab,
			const double rmin, const double rmax,
			const unsigned int bspline_degree,
			const unsigned int bspline_breakpoints,
			const bool linear);

grid *gridInitKeplerian(const unsigned long nr, const double rmin, 
			const double rmax, const double m, 
			const bool linear);

grid *gridInitFlat(const unsigned long nr, const double rmin, 
		   const double rmax, const double vphi, 
		   const bool linear);

bool outputAlloc(const unsigned long nOut, const bool eos_func,
		 const bool massSrc_func, const bool intEnSrc_func,
		 const unsigned long nUserOut, const grid *grd,
		 double **tOut, double **colOut, double **presOut,
		 double **eIntOut, double **mBndOut, double **eBndOut,
		 double **mSrcOut, double **eSrcOut, double **userOut);

void outputFree(double **tOut, double **colOut, double **presOut,
		double **eIntOut, double **mBndOut, double **eBndOut,
		double **mSrcOut, double **eSrcOut, double **userOut);

bool outputResize(const unsigned long nOut, const bool eos_func,
		  const bool massSrc_func, const bool intEnSrc_func,
		  const unsigned long nUserOut, const grid *grd,
		  double **tOut, double **colOut, double **presOut,
		  double **eIntOut, double **mBndOut, double **eBndOut,
		  double **mSrcOut, double **eSrcOut, double **userOut);

wksp *wkspAlloc(const unsigned long nr);

void gridFree(grid *grd);

void wkspFree(wksp *w);

#endif
