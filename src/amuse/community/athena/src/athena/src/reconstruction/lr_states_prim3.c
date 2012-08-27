#include "../copyright.h"
/*============================================================================*/
/*! \file lr_states_prim3.c
 *  \brief Third order (piecewise parabolic) spatial reconstruction in the
 *   primitive variables using the extremum-preserving limiters of
 *   Colella & Sekora.
 *
 * PURPOSE: Third order (piecewise parabolic) spatial reconstruction in the
 *   primitive variables using the extremum-preserving limiters of
 *   Colella & Sekora.  With the CTU integrator, a time-evolution
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
 *   the left-and right-side of cell center.  Thus (see Step 19),
 * -   W_{L,i-1/2} = Wrv(i-1);  W_{R,i-1/2} = Wlv(i)
 *
 * REFERENCE:
 * - P. Colella & P. Woodward, "The piecewise parabolic method (PPM) for
 *     gas-dynamical simulations", JCP, 54, 174 (1984).
 * - P. Colella & M. Sekora, "A limiter for PPM that preserves accuracy at
 *     smooth extrema", JCP, submitted (2007)
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

#ifdef THIRD_ORDER_PRIM

static Real **pW=NULL, **Whalf=NULL;

/*----------------------------------------------------------------------------*/
/*! \fn void lr_states(const GridS *pG, const Prim1DS W[], const Real Bxc[],
 *               const Real dt, const Real dx, const int il, const int iu,
 *               Prim1DS Wl[], Prim1DS Wr[], const int dir)
 *  \brief Computes L/R states
 *
 * Input Arguments:
 * - W = PRIMITIVE variables at cell centers along 1-D slice
 * - Bxc = B in direction of slice at cell center
 * - dtodx = dt/dx
 * - il,iu = lower and upper indices of zone centers in slice
 * W and Bxc must be initialized over [il-3:iu+3]
 *
 * Output Arguments:
 * - Wl,Wr = L/R-states of PRIMITIVE variables at interfaces over [il:iu+1]
 */

void lr_states(const GridS *pG, const Prim1DS W[], const Real Bxc[],
               const Real dt, const Real dx, const int il, const int iu,
               Prim1DS Wl[], Prim1DS Wr[], const int dir)
{
  int i,n,m;
  Real lim_slope,qa,qb,qc,qx;
  Real ev[NWAVE],rem[NWAVE][NWAVE],lem[NWAVE][NWAVE];
  Real d2Wc[NWAVE+NSCALARS],d2Wl[NWAVE+NSCALARS];
  Real d2Wr[NWAVE+NSCALARS],d2W [NWAVE+NSCALARS];
  Real d2Wlim[NWAVE+NSCALARS];
  Real Wlv[NWAVE+NSCALARS],Wrv[NWAVE+NSCALARS];
  Real dW[NWAVE+NSCALARS],W6[NWAVE+NSCALARS];
  Real *pWl, *pWr;
  Real dtodx = dt/dx;

/* Set pointer to primitive variables */
  for (i=il-3; i<=iu+3; i++) pW[i] = (Real*)&(W[i]);

#ifdef CTU_INTEGRATOR /* zero eigenmatrices if using CTU integrator */
  for (n=0; n<NWAVE; n++) {
    for (m=0; m<NWAVE; m++) {
      rem[n][m] = 0.0;
      lem[n][m] = 0.0;
    }
  }
#endif /* CTU_INTEGRATOR */

/*--- Step 1. ------------------------------------------------------------------
 * Compute interface states (CS eqns 12-15) over entire 1D pencil.  Using usual
 * Athena notation that index i for face-centered quantities denotes L-edge
 * (interface i-1/2), then Whalf[i] = W[i-1/2]. */

  for (i=il-1; i<=iu+2; i++) {
    for (n=0; n<(NWAVE+NSCALARS); n++) {
      Whalf[i][n]=(7.0*(pW[i-1][n]+pW[i][n]) - (pW[i-2][n]+pW[i+1][n]))/12.0;
    }
    for (n=0; n<(NWAVE+NSCALARS); n++) {
      d2Wc[n] = 3.0*(pW[i-1][n] - 2.0*Whalf[i][n] + pW[i][n]);
      d2Wl[n] = (pW[i-2][n] - 2.0*pW[i-1][n] + pW[i  ][n]);
      d2Wr[n] = (pW[i-1][n] - 2.0*pW[i  ][n] + pW[i+1][n]);
      d2Wlim[n] = 0.0;
      lim_slope = MIN(fabs(d2Wl[n]),fabs(d2Wr[n]));
      if (d2Wc[n] > 0.0 && d2Wl[n] > 0.0 && d2Wr[n] > 0.0) {
        d2Wlim[n] = SIGN(d2Wc[n])*MIN(1.25*lim_slope,fabs(d2Wc[n]));
      }
      if (d2Wc[n] < 0.0 && d2Wl[n] < 0.0 && d2Wr[n] < 0.0) {
        d2Wlim[n] = SIGN(d2Wc[n])*MIN(1.25*lim_slope,fabs(d2Wc[n]));
      }
    }
    for (n=0; n<(NWAVE+NSCALARS); n++) {
      Whalf[i][n] = 0.5*((pW[i-1][n]+pW[i][n]) - d2Wlim[n]/3.0);
    }
  }

/*====================== START BIG LOOP OVER i ============================*/
  for (i=il-1; i<=iu+1; i++) {

/*--- Step 2. ------------------------------------------------------------------
 * Compute L/R values
 * Wlv = W at left  side of cell-center = W[i-1/2] = a_{j,-} in CS
 * Wrv = W at right side of cell-center = W[i+1/2] = a_{j,+} in CS
 */

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      Wlv[n] = Whalf[i  ][n];
      Wrv[n] = Whalf[i+1][n];
    }

/*--- Step 3. ------------------------------------------------------------------
 * Construct parabolic interpolant (CS eqn 16-19) */

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      qa = (Wrv[n]-pW[i][n])*(pW[i][n]-Wlv[n]);
      qb = (pW[i-1][n]-pW[i][n])*(pW[i][n]-pW[i+1][n]);
      if (qa <= 0.0 && qb <= 0.0) {
        qc = 6.0*(pW[i][n] - 0.5*(Wlv[n]+Wrv[n]));
        d2W [n] = -2.0*qc;
        d2Wc[n] = (pW[i-1][n] - 2.0*pW[i  ][n] + pW[i+1][n]);
        d2Wl[n] = (pW[i-2][n] - 2.0*pW[i-1][n] + pW[i  ][n]);
        d2Wr[n] = (pW[i  ][n] - 2.0*pW[i+1][n] + pW[i+2][n]);
        d2Wlim[n] = 0.0;
        lim_slope = MIN(fabs(d2Wl[n]),fabs(d2Wr[n]));
        lim_slope = MIN(fabs(d2Wc[n]),lim_slope);
        if (d2Wc[n] > 0.0 && d2Wl[n] > 0.0 && d2Wr[n] > 0.0 && d2W[n] > 0.0) {
          d2Wlim[n] = SIGN(d2W[n])*MIN(1.25*lim_slope,fabs(d2W[n]));
        }
        if (d2Wc[n] < 0.0 && d2Wl[n] < 0.0 && d2Wr[n] < 0.0 && d2W[n] < 0.0) {
          d2Wlim[n] = SIGN(d2W[n])*MIN(1.25*lim_slope,fabs(d2W[n]));
        }
        if (d2W[n] == 0.0) {
          Wlv[n] = pW[i][n];
          Wrv[n] = pW[i][n];
        } else {
          Wlv[n] = pW[i][n] + (Wlv[n] - pW[i][n])*d2Wlim[n]/d2W[n];
          Wrv[n] = pW[i][n] + (Wrv[n] - pW[i][n])*d2Wlim[n]/d2W[n];
        }
      }
    }

/*--- Step 4. ------------------------------------------------------------------
 * Monotonize again (CS eqn 21-22) */

/*  First tests showed following does not work.  Not sure why, so leave out.
    for (n=0; n<(NWAVE+NSCALARS); n++) {
      s = SIGN(pW[i+1][n] - pW[i-1][n]);
      alphal = Wlv[n] - pW[i][n];
      alphar = Wrv[n] - pW[i][n];
      if (fabs(alphar)>=2.0*fabs(alphal)) {
        dalpha = pW[i+1][n] - pW[i][n]; 
        dI = -0.25*alphar*alphar/(alphar + alphal);
        if (s*dI >= s*dalpha) {
          Wlv[n] = pW[i][n]-2.0*(dalpha + s*sqrt(dalpha*dalpha-dalpha*alphar));
          Wrv[n] = pW[i][n]-2.0*(dalpha + s*sqrt(dalpha*dalpha-dalpha*alphal));
        }
      }
      if (fabs(alphal)>=2.0*fabs(alphar)){
        dalpha = pW[i-1][n] - pW[i][n]; 
        dI = -0.25*alphal*alphal/(alphar + alphal);
        if (s*dI >= s*dalpha) {
          Wlv[n] = pW[i][n]-2.0*(dalpha + s*sqrt(dalpha*dalpha-dalpha*alphar));
          Wrv[n] = pW[i][n]-2.0*(dalpha + s*sqrt(dalpha*dalpha-dalpha*alphal));
        }
      }
    }
*/

/* Monotonize again (CW eqn 1.10), ensure they lie between neighboring
 * cell-centered vals */

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      qa = (Wrv[n]-pW[i][n])*(pW[i][n]-Wlv[n]);
      qb = Wrv[n]-Wlv[n];
      qc = 6.0*(pW[i][n] - 0.5*(Wlv[n]+Wrv[n]));
      if (qa <= 0.0) {
        Wlv[n] = pW[i][n];
        Wrv[n] = pW[i][n];
      } else if ((qb*qc) > (qb*qb)) {
        Wlv[n] = 3.0*pW[i][n] - 2.0*Wrv[n];
      } else if ((qb*qc) < -(qb*qb)) {
        Wrv[n] = 3.0*pW[i][n] - 2.0*Wlv[n];
      }
    }

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      Wlv[n] = MAX(MIN(pW[i][n],pW[i-1][n]),Wlv[n]);
      Wlv[n] = MIN(MAX(pW[i][n],pW[i-1][n]),Wlv[n]);
      Wrv[n] = MAX(MIN(pW[i][n],pW[i+1][n]),Wrv[n]);
      Wrv[n] = MIN(MAX(pW[i][n],pW[i+1][n]),Wrv[n]);
    }

/*--- Step 5. ------------------------------------------------------------------
 * Set L/R values */

    pWl = (Real *) &(Wl[i+1]);
    pWr = (Real *) &(Wr[i]);

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      pWl[n] = Wrv[n];
      pWr[n] = Wlv[n];
    }

#ifdef CTU_INTEGRATOR /* only include steps below if CTU integrator */

/*--- Step 6. ------------------------------------------------------------------
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

/*--- Step 7. ------------------------------------------------------------------
 * Compute coefficients of interpolation parabolae (CW eqn 1.5) */

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      dW[n] = Wrv[n] - Wlv[n];
      W6[n] = 6.0*(pW[i][n] - 0.5*(Wlv[n] + Wrv[n]));
    }

/*--- Step 8. ------------------------------------------------------------------
 * Integrate linear interpolation function over domain of dependence defined by
 * max(min) eigenvalue (CW eqn 1.12)
 */

    qx = TWO_3RDS*MAX(ev[NWAVE-1],0.0)*dtodx;
    for (n=0; n<(NWAVE+NSCALARS); n++) {
      pWl[n] -= 0.75*qx*(dW[n] - (1.0 - qx)*W6[n]);
    }

    qx = -TWO_3RDS*MIN(ev[0],0.0)*dtodx;
    for (n=0; n<(NWAVE+NSCALARS); n++) {
      pWr[n] += 0.75*qx*(dW[n] + (1.0 - qx)*W6[n]);
    }

/*--- Step 9. ------------------------------------------------------------------
 * Then subtract amount of each wave m that does not reach the interface
 * during timestep (CW eqn 3.5ff).  For HLL fluxes, must subtract waves that
 * move in both directions, but only to 2nd order.
 */

    for (n=0; n<NWAVE; n++) {
      if (ev[n] > 0.) {
	qa  = 0.0;
        qb = 0.5*dtodx*(ev[NWAVE-1]-ev[n]);
        qc = 0.5*dtodx*dtodx*TWO_3RDS*(ev[NWAVE-1]*ev[NWAVE-1] - ev[n]*ev[n]);
	for (m=0; m<NWAVE; m++) {
	  qa += lem[n][m]*(qb*(dW[m]-W6[m]) + qc*W6[m]);
	}
	for (m=0; m<NWAVE; m++) pWl[m] += qa*rem[m][n];
/* For HLL fluxes, subtract wave moving away from interface as well. */
#if defined(HLLE_FLUX) || defined(HLLC_FLUX) || defined(HLLD_FLUX)
        qa  = 0.0;
        qb = 0.5*dtodx*(ev[0]-ev[n]);
        qc = 0.5*dtodx*dtodx*TWO_3RDS*(ev[0]*ev[0] - ev[n]*ev[n]);
        for (m=0; m<NWAVE; m++) {
          qa += lem[n][m]*(qb*(dW[m]+W6[m]) + qc*W6[m]);
        }

        for (m=0; m<NWAVE; m++) pWr[m] += qa*rem[m][n];
#endif /* HLL_FLUX */
      }
    }

    for (n=0; n<NWAVE; n++) {
      if (ev[n] < 0.) {
        qa = 0.0;
        qb = 0.5*dtodx*(ev[0]-ev[n]);
        qc = 0.5*dtodx*dtodx*TWO_3RDS*(ev[0]*ev[0] - ev[n]*ev[n]);
        for (m=0; m<NWAVE; m++) {
          qa += lem[n][m]*(qb*(dW[m]+W6[m]) + qc*W6[m]);
        }
        for (m=0; m<NWAVE; m++) pWr[m] += qa*rem[m][n];
/* For HLL fluxes, subtract wave moving away from interface as well. */
#if defined(HLLE_FLUX) || defined(HLLC_FLUX) || defined(HLLD_FLUX)
        qa = 0.0;
        qb = 0.5*dtodx*(ev[NWAVE-1]-ev[n]);
        qc = 0.5*dtodx*dtodx*TWO_3RDS*(ev[NWAVE-1]*ev[NWAVE-1] - ev[n]*ev[n]);
        for (m=0; m<NWAVE; m++) {
          qa += lem[n][m]*(qb*(dW[m]-W6[m]) + qc*W6[m]);
        }

        for (m=0; m<NWAVE; m++) pWl[m] += qa*rem[m][n];
#endif /* HLL_FLUX */
      }
    }

/* Wave subtraction for advected quantities */
    for (n=NWAVE; n<(NWAVE+NSCALARS); n++) {
      if (W[i].Vx > 0.) {
        qb = 0.5*dtodx*(ev[NWAVE-1]-W[i].Vx);
        qc = 0.5*dtodx*dtodx*TWO_3RDS*(SQR(ev[NWAVE-1]) - SQR(W[i].Vx));
        pWl[n] += qb*(dW[m]-W6[m]) + qc*W6[m];
      } else if (W[i].Vx < 0.) {
        qb = 0.5*dtodx*(ev[0]-W[i].Vx);
        qc = 0.5*dtodx*dtodx*TWO_3RDS*(ev[0]*ev[0] - W[i].Vx*W[i].Vx);
        pWr[n] += qb*(dW[m]+W6[m]) + qc*W6[m];
      }
    }

#endif /* CTU_INTEGRATOR */

  } /*======================= END BIG LOOP OVER i ========================*/

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void lr_states_init(MeshS *pM)
 *  \brief Allocate enough memory for work arrays */

void lr_states_init(MeshS *pM)
{
  int nmax,size1=0,size2=0,size3=0,nl,nd;

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

  if ((Whalf = (Real**)calloc_2d_array(nmax, NWAVE, sizeof(Real))) == NULL)
    goto on_error;

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
  if (Whalf != NULL) free_2d_array(Whalf);
  return;
}

#endif /* THIRD_ORDER_PRIM */
