#include "../copyright.h"
/*============================================================================*/
/*! \file lr_states_plm.c
 *  \brief Second order (piecewise linear) spatial reconstruction using
 *   characteristic interpolation in the primitive variables.  
 *
 * PURPOSE: Second order (piecewise linear) spatial reconstruction using
 *   characteristic interpolation in the primitive variables.  With the CTU
 *   integrator, a time-evolution (characteristic tracing) step is used to
 *   interpolate interface values to the half time level {n+1/2}. 
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
 * - lr_states_destruct() - frees memory for static global arrays
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef SECOND_ORDER_CHAR
#ifdef SPECIAL_RELATIVITY
#error : PLM reconstruction (order=2) cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */

static Real **pW=NULL;

/*----------------------------------------------------------------------------*/
/*! \fn void lr_states(const GridS *pG, const Prim1DS W[], const Real Bxc[], 
 *               const Real dt, const Real dx, const int il, const int iu, 
 *               Prim1DS Wl[], Prim1DS Wr[], const int dir)
 *  \brief Computes L/R states
 * Input Arguments:
 *   W = PRIMITIVE variables at cell centers along 1-D slice
 *   Bxc = B in direction of slice at cell center
 *   dtodx = dt/dx
 *   il,iu = lower and upper indices of zone centers in slice
 * W and Bxc must be initialized over [il-2:iu+2]
 *
 * Output Arguments:
 *   Wl,Wr = L/R-states of PRIMITIVE variables at interfaces over [il:iu+1]
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
  Real dac[NWAVE+NSCALARS],dal[NWAVE+NSCALARS];
  Real dar[NWAVE+NSCALARS],dag[NWAVE+NSCALARS],da[NWAVE+NSCALARS];
  Real Wlv[NWAVE+NSCALARS],Wrv[NWAVE+NSCALARS];
  Real dW[NWAVE+NSCALARS],dWm[NWAVE+NSCALARS];
  Real *pWl, *pWr;
  Real qx1,qx2,C;

  /* ADDITIONAL VARIABLES REQUIRED FOR CYLINDRICAL COORDINATES */
  Real zl,zr,zc,gamma_curv,opg,omg,beta,betai;
  int hllallwave_flag = 0;
  const Real dtodx = dt/dx;
#ifdef CYLINDRICAL
  const Real *r=pG->r, *ri=pG->ri;
#endif /* CYLINDRICAL */

/* Zero eigenmatrices, set pointer to primitive variables */
  for (n=0; n<NWAVE; n++) {
    for (m=0; m<NWAVE; m++) {
      rem[n][m] = 0.0;
      lem[n][m] = 0.0;
    }
  }
  for (i=il-2; i<=iu+2; i++) pW[i] = (Real*)&(W[i]);

/*========================== START BIG LOOP OVER i =======================*/
  for (i=il-1; i<=iu+1; i++) {

/*--- Step 1. ------------------------------------------------------------------
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

/*--- Step 2. ------------------------------------------------------------------
 * Compute centered, L/R, and van Leer differences of primitive variables
 * Note we access contiguous array elements by indexing pointers for speed */

#ifdef CYLINDRICAL
    if (dir==1) {
      /* compute cylindrical weighting factors */
     zc = 1.0/(1.0 - SQR(dx)/(12.0*r[i+1]*r[i-1]));
     zl = 1.0/(1.0 - SQR(dx)/(12.0*r[i  ]*r[i-1]));
     zr = 1.0/(1.0 - SQR(dx)/(12.0*r[i+1]*r[i  ]));
    }
#endif
    for (n=0; n<(NWAVE+NSCALARS); n++) {
      dWc[n] = pW[i+1][n] - pW[i-1][n];
      dWl[n] = pW[i][n]   - pW[i-1][n];
      dWr[n] = pW[i+1][n] - pW[i][n];
#ifdef CYLINDRICAL
      if (dir==1) {
        dWc[n] *= zc;
        dWl[n] *= zl;
        dWr[n] *= zr;
      }
#endif
      if (dWl[n]*dWr[n] > 0.0) {
        dWg[n] = 2.0*dWl[n]*dWr[n]/(dWl[n]+dWr[n]);
      } else {
        dWg[n] = 0.0;
      }
    }

/*--- Step 3. ------------------------------------------------------------------
 * Project differences in primitive variables along characteristics */

    for (n=0; n<NWAVE; n++) {
      dac[n] = lem[n][0]*dWc[0];
      dal[n] = lem[n][0]*dWl[0];
      dar[n] = lem[n][0]*dWr[0];
      dag[n] = lem[n][0]*dWg[0];
      for (m=1; m<NWAVE; m++) {
	dac[n] += lem[n][m]*dWc[m];
	dal[n] += lem[n][m]*dWl[m];
	dar[n] += lem[n][m]*dWr[m];
	dag[n] += lem[n][m]*dWg[m];
      }
    }

/* Advected variables are treated differently; for them the right and left
 * eigenmatrices are simply the identitiy matrix.
 */
#if (NSCALARS > 0)
    for (n=NWAVE; n<(NWAVE+NSCALARS); n++) {
      dac[n] = dWc[n];
      dal[n] = dWl[n];
      dar[n] = dWr[n];
      dag[n] = dWg[n];
    }
#endif

/*--- Step 4. ------------------------------------------------------------------
 * Apply monotonicity constraints to characteristic projections */

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      da[n] = 0.0;
      if (dal[n]*dar[n] > 0.0) {
        lim_slope1 = MIN(    fabs(dal[n]),fabs(dar[n]));
        lim_slope2 = MIN(0.5*fabs(dac[n]),fabs(dag[n]));
        da[n] = SIGN(dac[n])*MIN(2.0*lim_slope1,lim_slope2);
      }
    }

/*--- Step 5. ------------------------------------------------------------------
 * Project monotonic slopes in characteristic back to primitive variables  */

    for (n=0; n<NWAVE; n++) {
      dWm[n] = da[0]*rem[n][0];
      for (m=1; m<NWAVE; m++) {
        dWm[n] += da[m]*rem[n][m];
      }
    }

#if (NSCALARS > 0)
    for (n=NWAVE; n<(NWAVE+NSCALARS); n++) {
      dWm[n] = da[n];
    }
#endif

/*--- Step 6. ------------------------------------------------------------------
 * Limit velocity difference to sound speed (deleted).  Was added to make
 * turbulence runs more robust, but found added noise to Noh shocktube and 
 * other tests.  See r995 and earlier for this step */

/*--- Step 7. ------------------------------------------------------------------
 * Compute L/R values, ensure they lie between neighboring cell-centered vals */

    gamma_curv = 0.0;
#ifdef CYLINDRICAL
    if (dir==1) gamma_curv = dx/(6.0*r[i]);
#endif
    opg = 1.0 + gamma_curv;
    omg = 1.0 - gamma_curv;
    beta  = omg/opg;
    betai = opg/omg;

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      Wlv[n] = pW[i][n] - 0.5*dWm[n]*opg;
      Wrv[n] = pW[i][n] + 0.5*dWm[n]*omg;
    }

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      C = Wrv[n] + beta*Wlv[n];
      Wlv[n] = MAX(MIN(pW[i][n],pW[i-1][n]),Wlv[n]);
      Wlv[n] = MIN(MAX(pW[i][n],pW[i-1][n]),Wlv[n]);
      Wrv[n] = C - beta*Wlv[n];

      Wrv[n] = MAX(MIN(pW[i][n],pW[i+1][n]),Wrv[n]);
      Wrv[n] = MIN(MAX(pW[i][n],pW[i+1][n]),Wrv[n]);
      Wlv[n] = (C - Wrv[n])*betai;
    }

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      dW[n] = Wrv[n] - Wlv[n];
    }

/*--- Step 8. ------------------------------------------------------------------
 * Integrate linear interpolation function over domain of dependence defined by
 * max(min) eigenvalue
 */

    pWl = (Real *) &(Wl[i+1]);
    pWr = (Real *) &(Wr[i]);

#ifndef CTU_INTEGRATOR

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      pWl[n] = Wrv[n];
      pWr[n] = Wlv[n];
    }

#elif defined(HLLE_FLUX) || defined(HLLC_FLUX) || defined(HLLD_FLUX)

    for (n=0; n<NWAVE; n++) {
      pWl[n] = Wrv[n];
      pWr[n] = Wlv[n];
    }

#ifdef HLL_ALL_WAVE
    hllallwave_flag = 1;
#endif

    for (n=0; n<NWAVE; n++) {
      qa = 0.0;
      if (hllallwave_flag || ev[n] > 0.) {
        qx = 0.5*dtodx*ev[n];
#ifdef CYLINDRICAL
        if (dir==1) 
          qx *= 1.0 - dx*qx/(3.0*(ri[i+1]-dx*qx));
#endif
        for (m=0; m<NWAVE; m++) {
          qa += lem[n][m]*qx*dW[m];
        }
        for (m=0; m<NWAVE; m++) pWl[m] -= qa*rem[m][n];
      }

      qa = 0.0;
      if (hllallwave_flag || ev[n] < 0.) {
        qx = 0.5*dtodx*ev[n];
#ifdef CYLINDRICAL
        if (dir==1)
          qx *= 1.0 - dx*qx/(3.0*(ri[i]-dx*qx));
#endif
        for (m=0; m<NWAVE; m++) {
          qa += lem[n][m]*qx*dW[m];
        }
        for (m=0; m<NWAVE; m++) pWr[m] -= qa*rem[m][n];
      }
    }

#else  /* include steps 8-9 only if using CTU integrator (AND NOT HLL) */   
    qx = 0.5*MAX(ev[NWAVE-1],0.0)*dtodx;
#ifdef CYLINDRICAL
    if (dir==1) 
      qx *= 1.0 - dx*qx/(3.0*(ri[i+1]-dx*qx));
#endif
    for (n=0; n<(NWAVE+NSCALARS); n++) {
      pWl[n] = Wrv[n] - qx*dW[n];
    }

    qx = -0.5*MIN(ev[0],0.0)*dtodx;
#ifdef CYLINDRICAL
    if (dir==1) 
      qx *= 1.0 + dx*qx/(3.0*(ri[i]+dx*qx));
#endif
    for (n=0; n<(NWAVE+NSCALARS); n++) {
      pWr[n] = Wlv[n] + qx*dW[n];
    }


/*--- Step 9. ------------------------------------------------------------------
 * Then subtract amount of each wave n that does not reach the interface
 * during timestep (CW eqn 3.5ff).  For HLL fluxes, must subtract waves that
 * move in both directions.
 */

    for (n=0; n<NWAVE; n++) {
      if (ev[n] >= 0.0) {
        qa  = 0.0;
        qx1 = 0.5*dtodx*ev[NWAVE-1];
        qx2 = 0.5*dtodx*ev[n];
#ifdef CYLINDRICAL
        if (dir==1) {
          qx1 *= 1.0 - dx*qx1/(3.0*(ri[i+1]-dx*qx1));
          qx2 *= 1.0 - dx*qx2/(3.0*(ri[i+1]-dx*qx2));
        }
#endif
        qx = qx1 - qx2;
        for (m=0; m<NWAVE; m++) {
          qa += lem[n][m]*qx*dW[m];
        }
        for (m=0; m<NWAVE; m++) pWl[m] += qa*rem[m][n];
      }
    }

    for (n=0; n<NWAVE; n++) {
      if (ev[n] <= 0.0) {
        qa = 0.0;
        qx1 = -0.5*dtodx*ev[0];
        qx2 = -0.5*dtodx*ev[n];
#ifdef CYLINDRICAL
        if (dir==1) {
          qx1 *= 1.0 + dx*qx1/(3.0*(ri[i]+dx*qx1));
          qx2 *= 1.0 + dx*qx2/(3.0*(ri[i]+dx*qx2));
        }
#endif
        qx = -qx1 + qx2;
        for (m=0; m<NWAVE; m++) {
          qa += lem[n][m]*qx*dW[m];
        }
        for (m=0; m<NWAVE; m++) pWr[m] += qa*rem[m][n];
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
  return;
}

#endif /* SECOND_ORDER_CHAR */
