#include "../copyright.h"
/*============================================================================*/
/*! \file force.c
 *  \brief Computes 1D fluxes using a Riemann solver similar, but not
 *   identical, to Toro's FORCE (First-ORder-CEntred) flux.
 *
 * PURPOSE: Computes 1D fluxes using a Riemann solver similar, but not
 *   identical, to Toro's FORCE (First-ORder-CEntred) flux.  It uses the
 *   average of a Lax-Wendroff and HLLE flux for sub-sonic flow, and the
 *   appropriate upwind flux otherwise.
 *
 * HISTORY: -- TAG -- 3/4/2005
 *
 * REFERENCES:
 * - E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics",
 *   2nd ed., Springer-Verlag, Berlin, (1999), section 7.4.2.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - fluxes() - all Riemann solvers in Athena must have this function name and
 *              use the same argument list as defined in rsolvers/prototypes.h*/
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef FORCE_FLUX
/*----------------------------------------------------------------------------*/
/*! \fn void fluxes(const Cons1DS Ul, const Cons1DS Ur, 
 *            const Prim1DS Wl, const Prim1DS Wr,
 *            const Real Bxi, Cons1DS *pFlux)
 *  \brief Computes 1D fluxes
 *   Input Arguments:
 *  -  Bxi = B in direction of 1D slice at cell interface
 *  -  Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *   Output Arguments:
 *  -  pFlux = pointer to fluxes of CONSERVED variables at cell interface
 */

void fluxes(const Cons1DS Ul, const Cons1DS Ur, 
            const Prim1DS Wl, const Prim1DS Wr,
            const Real Bxi, Cons1DS *pFlux)
{
  Real sqrtdl,sqrtdr,isdlpdr,droe,v1roe,v2roe,v3roe,pbl=0.0,pbr=0.0;
  Real asq,vaxsq=0.0,qsq,cfsq,cfl,cfr,bp,bm,ct2=0.0,tmp;
#ifndef ISOTHERMAL
  Real hroe;
#endif
#ifdef MHD
  Real b2roe,b3roe,x,y;
#endif
  Real ev[NWAVE],al,ar;
  Real *pFl, *pFc, *pFr, *pUc, *pF;
  Prim1DS Wc;
  Cons1DS Fl, Fc, Fr, Uc;
  int n;

/* The first 5 steps are identical to those in hlle fluxes */
/*--- Step 1. ------------------------------------------------------------------
 * Convert left- and right- states in conserved to primitive variables.
 */
/*
  pbl = Cons1D_to_Prim1D(&Ul,&Wl,&Bxi);
  pbr = Cons1D_to_Prim1D(&Ur,&Wr,&Bxi);
*/

/*--- Step 2. ------------------------------------------------------------------
 * Compute Roe-averaged data from left- and right-states
 */

  sqrtdl = sqrt((double)Wl.d);
  sqrtdr = sqrt((double)Wr.d);
  isdlpdr = 1.0/(sqrtdl + sqrtdr);

  droe  = sqrtdl*sqrtdr;
  v1roe = (sqrtdl*Wl.Vx + sqrtdr*Wr.Vx)*isdlpdr;
  v2roe = (sqrtdl*Wl.Vy + sqrtdr*Wr.Vy)*isdlpdr;
  v3roe = (sqrtdl*Wl.Vz + sqrtdr*Wr.Vz)*isdlpdr;

/*
 * The Roe average of the magnetic field is defined differently.
 */

#ifdef MHD
  b2roe = (sqrtdr*Wl.By + sqrtdl*Wr.By)*isdlpdr;
  b3roe = (sqrtdr*Wl.Bz + sqrtdl*Wr.Bz)*isdlpdr;
  x = 0.5*(SQR(Wl.By - Wr.By) + SQR(Wl.Bz - Wr.Bz))/(SQR(sqrtdl + sqrtdr));
  y = 0.5*(Wl.d + Wr.d)/droe;
  pbl = 0.5*(SQR(Bxi) + SQR(Wl.By) + SQR(Wl.Bz));
  pbr = 0.5*(SQR(Bxi) + SQR(Wr.By) + SQR(Wr.Bz));
#endif

/*
 * Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
 * rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
 */

#ifndef ISOTHERMAL
  hroe  = ((Ul.E + Wl.P + pbl)/sqrtdl + (Ur.E + Wr.P + pbr)/sqrtdr)*isdlpdr;
#endif

/*--- Step 3. ------------------------------------------------------------------
 * Compute eigenvalues using Roe-averaged values
 */

#ifdef HYDRO
#ifdef ISOTHERMAL
  esys_roe_iso_hyd(v1roe, v2roe, v3roe, ev, NULL, NULL);
#else
  esys_roe_adb_hyd(v1roe, v2roe, v3roe, hroe, ev, NULL, NULL);
#endif /* ISOTHERMAL */
#endif /* HYDRO */

#ifdef MHD
#ifdef ISOTHERMAL
 esys_roe_iso_mhd(droe,v1roe,v2roe,v3roe,     Bxi,b2roe,b3roe,x,y,ev,NULL,NULL);
#else
 esys_roe_adb_mhd(droe,v1roe,v2roe,v3roe,hroe,Bxi,b2roe,b3roe,x,y,ev,NULL,NULL);
#endif /* ISOTHERMAL */
#endif /* MHD */

/*--- Step 4. ------------------------------------------------------------------
 * Compute the max and min wave speeds
 */

/* left state */
#ifdef ISOTHERMAL
  asq = Iso_csound2;
#else
  asq = Gamma*Wl.P/Wl.d;
#endif
#ifdef MHD
  vaxsq = Bxi*Bxi/Wl.d;
  ct2 = (Ul.By*Ul.By + Ul.Bz*Ul.Bz)/Wl.d;
#endif
  qsq = vaxsq + ct2 + asq;
  tmp = vaxsq + ct2 - asq;
  cfsq = 0.5*(qsq + sqrt((double)(tmp*tmp + 4.0*asq*ct2)));
  cfl = sqrt((double)cfsq);

/* right state */
#ifdef ISOTHERMAL
  asq = Iso_csound2;
#else
  asq = Gamma*Wr.P/Wr.d; 
#endif
#ifdef MHD
  vaxsq = Bxi*Bxi/Wr.d;
  ct2 = (Ur.By*Ur.By + Ur.Bz*Ur.Bz)/Wr.d;
#endif
  qsq = vaxsq + ct2 + asq;
  tmp = vaxsq + ct2 - asq;
  cfsq = 0.5*(qsq + sqrt((double)(tmp*tmp + 4.0*asq*ct2)));
  cfr = sqrt((double)cfsq);

/* take max/min of Roe eigenvalues and L/R state wave speeds */
  ar = MAX(ev[NWAVE-1],(Wr.Vx + cfr));
  al = MIN(ev[0]      ,(Wl.Vx - cfl));

  bp = MAX(ar, 0.0);
  bm = MIN(al, 0.0);

/*--- Step 5. ------------------------------------------------------------------
 * Compute L/R fluxes along the line bm, bp
 */

  Fl.d  = Ul.Mx - bm*Ul.d;
  Fr.d  = Ur.Mx - bp*Ur.d;

  Fl.Mx = Ul.Mx*(Wl.Vx - bm);
  Fr.Mx = Ur.Mx*(Wr.Vx - bp);

  Fl.My = Ul.My*(Wl.Vx - bm);
  Fr.My = Ur.My*(Wr.Vx - bp);

  Fl.Mz = Ul.Mz*(Wl.Vx - bm);
  Fr.Mz = Ur.Mz*(Wr.Vx - bp);

#ifdef ISOTHERMAL
  Fl.Mx += Wl.d*Iso_csound2;
  Fr.Mx += Wr.d*Iso_csound2;
#else
  Fl.Mx += Wl.P;
  Fr.Mx += Wr.P;

  Fl.E  = Ul.E*(Wl.Vx - bm) + Wl.P*Wl.Vx;
  Fr.E  = Ur.E*(Wr.Vx - bp) + Wr.P*Wr.Vx;
#endif /* ISOTHERMAL */

#ifdef MHD
  Fl.Mx -= 0.5*(Bxi*Bxi - SQR(Wl.By) - SQR(Wl.Bz));
  Fr.Mx -= 0.5*(Bxi*Bxi - SQR(Wr.By) - SQR(Wr.Bz));

  Fl.My -= Bxi*Wl.By;
  Fr.My -= Bxi*Wr.By;
    
  Fl.Mz -= Bxi*Wl.Bz;
  Fr.Mz -= Bxi*Wr.Bz;

#ifndef ISOTHERMAL
  Fl.E += (pbl*Wl.Vx - Bxi*(Bxi*Wl.Vx + Wl.By*Wl.Vy + Wl.Bz*Wl.Vz));
  Fr.E += (pbr*Wr.Vx - Bxi*(Bxi*Wr.Vx + Wr.By*Wr.Vy + Wr.Bz*Wr.Vz));
#endif /* ISOTHERMAL */

  Fl.By = Wl.By*(Wl.Vx - bm) - Bxi*Wl.Vy;
  Fr.By = Wr.By*(Wr.Vx - bp) - Bxi*Wr.Vy;

  Fl.Bz = Wl.Bz*(Wl.Vx - bm) - Bxi*Wl.Vz;
  Fr.Bz = Wr.Bz*(Wr.Vx - bp) - Bxi*Wr.Vz;
#endif /* MHD */

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) {
    Fl.s[n] = Fl.d*Wl.x[n];
    Fr.s[n] = Fr.d*Wr.x[n];
  }
#endif

/*--- Step 6. ------------------------------------------------------------------
 * For supersonic flow, return the upwind flux.
 */

  if(al >= 0.0){
    *pFlux = Fl;
    return;
  }

  if(ar <= 0.0){
    *pFlux = Fr;
    return;
  }

/*--- Step 7. ------------------------------------------------------------------
 * Compute the LW flux, start with the HLL mean state */

  pFl = (Real *)&(Fl);
  pFr = (Real *)&(Fr);
  pUc = (Real *)&(Uc);
  tmp = 1.0/(ar - al);
  for (n=0; n<(NWAVE+NSCALARS); n++){
    pUc[n] = (pFl[n] - pFr[n])*tmp;
  }

/* Convert the HLL mean state to primitive variables */
  Cons1D_to_Prim1D(&Uc,&Wc,&Bxi);

/* Compute the LW flux along the line dx/dt = 0 */

  Fc.d  = Uc.Mx;
  Fc.Mx = Uc.Mx*Wc.Vx;
  Fc.My = Uc.My*Wc.Vx;
  Fc.Mz = Uc.Mz*Wc.Vx;

#ifdef ISOTHERMAL
  Fc.Mx += Wc.d*Iso_csound2;
#else
  Fc.Mx += Wc.P;

  Fc.E  = Uc.E*Wc.Vx + Wc.P*Wc.Vx;
#endif /* ISOTHERMAL */

#ifdef MHD
  Fc.Mx -= 0.5*(Bxi*Bxi - SQR(Wc.By) - SQR(Wc.Bz));
  Fc.My -= Bxi*Wc.By;
  Fc.Mz -= Bxi*Wc.Bz;

#ifndef ISOTHERMAL
  Fc.E += (pbl*Wc.Vx - Bxi*(Bxi*Wc.Vx + Wc.By*Wc.Vy + Wc.Bz*Wc.Vz));
#endif /* ISOTHERMAL */

  Fc.By = Wc.By*Wc.Vx - Bxi*Wc.Vy;
  Fc.Bz = Wc.Bz*Wc.Vx - Bxi*Wc.Vz;
#endif /* MHD */

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) {
    Fc.s[n] = Fc.d*Wc.x[n];
  }
#endif


/*--- Step 8. ------------------------------------------------------------------
 * Compute the average of the Lax-Wendroff & HLLE flux
 */
  pFl = (Real *)&(Fl);
  pFc = (Real *)&(Fc);
  pFr = (Real *)&(Fr);
  pF  = (Real *)pFlux;
  tmp = 0.25*(bp + bm)/(bp - bm);
  for (n=0; n<(NWAVE+NSCALARS); n++){
    pF[n] = 0.5*pFc[n] + 0.25*(pFl[n] + pFr[n]) + (pFl[n] - pFr[n])*tmp;
  }

  return;
}
#endif /* FORCE_FLUX */
