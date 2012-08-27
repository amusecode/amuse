#include "../copyright.h"
/*============================================================================*/
/*! \file lr_states_prim2.c
 *  \brief Second order (piecewise linear) spatial reconstruction in the
 *   primitive variables. 
 *
 * PURPOSE: Second order (piecewise linear) spatial reconstruction in the
 *   primitive variables. With the CTU integrator, a time-evolution
 *   (characteristic tracing) step is used to interpolate interface values
 *   to the half time level {n+1/2}.
 *
 *   Limiting is performed in the primitive (rather than characteristic)
 *   variables.  When used with the VL integrator, an eigenvalue decomposition
 *   is NOT needed.
 *
 * NOTATION: 
 * - W_{L,i-1/2} is reconstructed value on the left-side of interface at i-1/2
 * - W_{R,i-1/2} is reconstructed value on the right-side of interface at i-1/2
 *
 *   The L- and R-states at the left-interface in each cell are indexed i.
 * - W_{L,i-1/2} is denoted by Wl[i  ];   W_{R,i-1/2} is denoted by Wr[i  ]
 * - W_{L,i+1/2} is denoted by Wl[i+1];   W_{R,i+1/2} is denoted by Wr[i+1]
 *
 *   Internally, in this routine, Wlv and Wrv are the reconstructed values on
 *   the left-and right-side of cell center.  Thus (see Step 8),
 * -   W_{L,i-1/2} = Wrv(i-1);  W_{R,i-1/2} = Wlv(i)
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - lr_states()          - computes L/R states
 * - lr_states_init()     - initializes memory for static global arrays
 * - lr_states_destruct() - frees memory for static global arrays	      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef SECOND_ORDER_PRIM

static Real **pW=NULL;
#ifdef SPECIAL_RELATIVITY
static Real **vel=NULL;
#endif

/*----------------------------------------------------------------------------*/
/*! \fn void lr_states(const GridS *pG, const Prim1DS W[], const Real Bxc[],
 *               const Real dt, const Real dx, const int il, const int iu,
 *               Prim1DS Wl[], Prim1DS Wr[], const int dir)
 *  \brief Computes L/R states
 * Input Arguments:
 * - W = PRIMITIVE variables at cell centers along 1-D slice
 * - Bxc = B in direction of slice at cell center
 * - dtodx = dt/dx
 * - il,iu = lower and upper indices of zone centers in slice
 * W and Bxc must be initialized over [il-2:iu+2]
 *
 * Output Arguments:
 * - Wl,Wr = L/R-states of PRIMITIVE variables at interfaces over [il:iu+1]
 */

void lr_states(const GridS *pG, const Prim1DS W[], const Real Bxc[],
               const Real dt, const Real dx, const int il, const int iu,
               Prim1DS Wl[], Prim1DS Wr[], const int dir)
{
  int i,n,m;
  Real lim_slope1,lim_slope2,qa,qx;
  Real ev[NWAVE],rem[NWAVE][NWAVE],lem[NWAVE][NWAVE];
  Real dWc[NWAVE+NSCALARS],dWl[NWAVE+NSCALARS];
  Real dWr[NWAVE+NSCALARS],dWg[NWAVE+NSCALARS];
  Real Wlv[NWAVE+NSCALARS],Wrv[NWAVE+NSCALARS];
  Real dW[NWAVE+NSCALARS],dWm[NWAVE+NSCALARS];
  Real *pWl, *pWr;
  Real dtodx = dt/dx;

/* Set pointer to primitive variables */
  for (i=il-2; i<=iu+2; i++) pW[i] = (Real*)&(W[i]);

#ifdef SPECIAL_RELATIVITY /* calculate 4-velocity */
  for (i=il-2; i<=iu+2; i++) {
    vel[i][0] = 1.0 - (SQR(pW[i][1]) + SQR(pW[i][2]) + SQR(pW[i][3]));
    vel[i][0] = 1.0/sqrt(vel[i][0]);
    vel[i][1] = vel[i][0]*pW[i][1];
    vel[i][2] = vel[i][0]*pW[i][2];
    vel[i][3] = vel[i][0]*pW[i][3];
  }
#endif


#ifdef CTU_INTEGRATOR /* zero eigenmatrices if using CTU integrator */
  for (n=0; n<NWAVE; n++) {
    for (m=0; m<NWAVE; m++) {
      rem[n][m] = 0.0;
      lem[n][m] = 0.0;
    }
  }
#endif /* CTU_INTEGRATOR */

/*========================== START BIG LOOP OVER i =======================*/
  for (i=il-1; i<=iu+1; i++) {

/*--- Step 1. ------------------------------------------------------------------
 * Compute centered, L/R, and van Leer differences of primitive variables
 * Note we access contiguous array elements by indexing pointers for speed */

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      dWc[n] = pW[i+1][n] - pW[i-1][n];
      dWl[n] = pW[i][n]   - pW[i-1][n];
      dWr[n] = pW[i+1][n] - pW[i][n];
      if (dWl[n]*dWr[n] > 0.0) {
        dWg[n] = 2.0*dWl[n]*dWr[n]/(dWl[n]+dWr[n]);
      } else {
        dWg[n] = 0.0;
      }
    }

/*--- Step 2. ------------------------------------------------------------------
 * Apply monotonicity constraints to differences in primitive vars. */
    for (n=0; n<(NWAVE+NSCALARS); n++) {
      dWm[n] = 0.0;
      if (dWl[n]*dWr[n] > 0.0) {
        lim_slope1 = MIN(    fabs(dWl[n]),fabs(dWr[n]));
        lim_slope2 = MIN(0.5*fabs(dWc[n]),fabs(dWg[n]));
        dWm[n] = SIGN(dWc[n])*MIN(2.0*lim_slope1,lim_slope2);
      }
    }

/*--- Step 3. ------------------------------------------------------------------
 * Compute L/R values, ensure they lie between neighboring cell-centered vals */

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      Wlv[n] = pW[i][n] - 0.5*dWm[n];
      Wrv[n] = pW[i][n] + 0.5*dWm[n];
    }

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      Wlv[n] = MAX(MIN(pW[i][n],pW[i-1][n]),Wlv[n]);
      Wlv[n] = MIN(MAX(pW[i][n],pW[i-1][n]),Wlv[n]);
      Wrv[n] = MAX(MIN(pW[i][n],pW[i+1][n]),Wrv[n]);
      Wrv[n] = MIN(MAX(pW[i][n],pW[i+1][n]),Wrv[n]);
    }

/*--- Step 4. ------------------------------------------------------------------
 * Set L/R values */

    pWl = (Real *) &(Wl[i+1]);
    pWr = (Real *) &(Wr[i]);

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      pWl[n] = Wrv[n];
      pWr[n] = Wlv[n];
    }

#ifdef SPECIAL_RELATIVITY /* reconstruct on 4-velocity for robustness */
/*--- Step 1. ------------------------------------------------------------------
 * Compute centered, L/R, and van Leer differences of primitive variables
 * Note we access contiguous array elements by indexing pointers for speed */

    for (n=0; n==3; n++) {
      dWc[n] = vel[i+1][n] - vel[i-1][n];
      dWl[n] = vel[i][n]   - vel[i-1][n];
      dWr[n] = vel[i+1][n] - vel[i][n];
      if (dWl[n]*dWr[n] > 0.0) {
        dWg[n] = 2.0*dWl[n]*dWr[n]/(dWl[n]+dWr[n]);
      } else {
        dWg[n] = 0.0;
      }
    }

/*--- Step 2. ------------------------------------------------------------------
 * Apply monotonicity constraints to differences in primitive vars. */
    for (n=0; n==3; n++) {
      dWm[n] = 0.0;
      if (dWl[n]*dWr[n] > 0.0) {
        lim_slope1 = MIN(    fabs(dWl[n]),fabs(dWr[n]));
        lim_slope2 = MIN(0.5*fabs(dWc[n]),fabs(dWg[n]));
        dWm[n] = SIGN(dWc[n])*MIN(2.0*lim_slope1,lim_slope2);
      }
    }

/*--- Step 3. ------------------------------------------------------------------
 * Compute L/R values, ensure they lie between neighboring cell-centered vals */

    for (n=0; n==3; n++) {
      Wlv[n] = pW[i][n] - 0.5*dWm[n];
      Wrv[n] = pW[i][n] + 0.5*dWm[n];
    }

    for (n=0; n==3; n++) {
      Wlv[n] = MAX(MIN(vel[i][n],vel[i-1][n]),Wlv[n]);
      Wlv[n] = MIN(MAX(vel[i][n],vel[i-1][n]),Wlv[n]);
      Wrv[n] = MAX(MIN(vel[i][n],vel[i+1][n]),Wrv[n]);
      Wrv[n] = MIN(MAX(vel[i][n],vel[i+1][n]),Wrv[n]);
    }

/*--- Step 4. ------------------------------------------------------------------
 * Set L/R values */

    for (n=1; n==3; n++) {
      pWl[n] = Wrv[n]/Wrv[0];
      pWr[n] = Wlv[n]/Wlv[0];
    }

#endif

#ifdef CTU_INTEGRATOR /* only include steps below if CTU integrator */

/*--- Step 5. ------------------------------------------------------------------
 * Compute eigensystem in primitive variables.  */

#ifdef HYDRO
#ifdef ISOTHERMAL
    esys_prim_iso_hyd(W[i].d,W[i].Vx,             ev,rem,lem);
#else
    esys_prim_adb_hyd(W[i].d,W[i].Vx,Gamma*W[i].P,ev,rem,lem);
#endif /* ISOTHERMAL */
#endif /* HYDRO */

#ifdef MHD
#ifdef ISOTHERMAL
    esys_prim_iso_mhd(
      W[i].d,W[i].Vx,             Bxc[i],W[i].By,W[i].Bz,ev,rem,lem);
#else
    esys_prim_adb_mhd(
      W[i].d,W[i].Vx,Gamma*W[i].P,Bxc[i],W[i].By,W[i].Bz,ev,rem,lem);
#endif /* ISOTHERMAL */
#endif /* MHD */


/*--- Step 6. ------------------------------------------------------------------
 * Integrate linear interpolation function over domain of dependence defined by
 * max(min) eigenvalue
 */

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      dW[n] = Wrv[n] - Wlv[n];
    }

    qx = 0.5*MAX(ev[NWAVE-1],0.0)*dtodx;
    for (n=0; n<(NWAVE+NSCALARS); n++) {
      pWl[n] -= qx*dW[n];
    }

    qx = -0.5*MIN(ev[0],0.0)*dtodx;
    for (n=0; n<(NWAVE+NSCALARS); n++) {
      pWr[n] += qx*dW[n];
    }

/*--- Step 7. ------------------------------------------------------------------
 * Then subtract amount of each wave n that does not reach the interface
 * during timestep (CW eqn 3.5ff).  For HLL fluxes, must subtract waves that
 * move in both directions.
 */

    for (n=0; n<NWAVE; n++) {
      if (ev[n] > 0.) {
	qa  = 0.0;
	for (m=0; m<NWAVE; m++) {
	  qa += lem[n][m]*0.5*dtodx*(ev[NWAVE-1]-ev[n])*dW[m];
	}
	for (m=0; m<NWAVE; m++) pWl[m] += qa*rem[m][n];
/* For HLL fluxes, subtract wave moving away from interface as well. */
#if defined(HLLE_FLUX) || defined(HLLC_FLUX) || defined(HLLD_FLUX)
        qa = 0.0;
        for (m=0; m<NWAVE; m++) {
          qa += lem[n][m]*0.5*dtodx*(ev[n]-ev[0])*dW[m];
        }
        for (m=0; m<NWAVE; m++) pWr[m] -= qa*rem[m][n];
#endif /* HLL_FLUX */
      }
    }

    for (n=0; n<NWAVE; n++) {
      if (ev[n] < 0.) {
        qa = 0.0;
        for (m=0; m<NWAVE; m++) {
          qa += lem[n][m]*0.5*dtodx*(ev[0]-ev[n])*dW[m];
        }
        for (m=0; m<NWAVE; m++) pWr[m] += qa*rem[m][n];
/* For HLL fluxes, subtract wave moving away from interface as well. */
#if defined(HLLE_FLUX) || defined(HLLC_FLUX) || defined(HLLD_FLUX)
	qa  = 0.0;
	for (m=0; m<NWAVE; m++) {
	  qa += lem[n][m]*0.5*dtodx*(ev[n]-ev[NWAVE-1])*dW[m];
	}
	for (m=0; m<NWAVE; m++) pWl[m] -= qa*rem[m][n];
#endif /* HLL_FLUX */
      }
    }

/* Wave subtraction for advected quantities */
    for (n=NWAVE; n<(NWAVE+NSCALARS); n++) {
      if (W[i].Vx > 0.) {
        pWl[n] += 0.5*dtodx*(ev[NWAVE-1]-W[i].Vx)*dW[n];
      } else if (W[i].Vx < 0.) {
        pWr[n] += 0.5*dtodx*(ev[0]-W[i].Vx)*dW[n];
      }
    }

#endif /* CTU_INTEGRATOR */

  } /*===================== END BIG LOOP OVER i ===========================*/

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void lr_states_init(MeshS *pM)
 *  \brief Allocate enough memory for work arrays */

void lr_states_init(MeshS *pM)
{
  int nmax,size1=0,size2=0,size3=0,nl,nd,n4v=4;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2, Nx3 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if (pM->Domain[nl][nd].Grid->Nx[0] > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
        if (pM->Domain[nl][nd].Grid->Nx[1] > size2){
          size2 = pM->Domain[nl][nd].Grid->Nx[1];
        }
        if (pM->Domain[nl][nd].Grid->Nx[2] > size3){
          size3 = pM->Domain[nl][nd].Grid->Nx[2];
        }
      }
    }
  }

  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  size3 = size3 + 2*nghost;
  nmax = MAX((MAX(size1,size2)),size3);

  if ((pW = (Real**)malloc(nmax*sizeof(Real*))) == NULL) goto on_error;
#ifdef SPECIAL_RELATIVITY
  if ((vel = (Real**)calloc_2d_array(nmax, n4v, sizeof(Real))) == NULL)
    goto on_error;
#endif

  return;
  on_error:
    lr_states_destruct();
    ath_error("[lr_states_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void lr_states_destruct(void)
 *  \brief Free memory used by work arrays */

void lr_states_destruct(void)
{
  if (pW != NULL) free(pW);
#ifdef SPECIAL_RELATIVITY
  if (vel != NULL) free_2d_array(vel);
#endif
  return;
}

#endif /* SECOND_ORDER_PRIM */
