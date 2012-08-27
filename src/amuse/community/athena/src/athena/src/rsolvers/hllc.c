#include "../copyright.h"
/*============================================================================*/
/*! \file hllc.c
 *  \brief Computes 1D fluxes using the HLLC Riemann solver.
 *
 * PURPOSE: Computes 1D fluxes using the HLLC Riemann solver, an extension of
 *   the HLLE fluxes to include the contact wave.  Currently only works for
 *   hydrodynamics.  For an extension to MHD, see hlld.c
 *
 * REFERENCES:
 * - E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics",
 *   2nd ed., Springer-Verlag, Berlin, (1999) chpt. 10.
 *
 * - P. Batten, N. Clarke, C. Lambert, and D. M. Causon,
 *   "On the Choice of Wavespeeds for the HLLC Riemann Solver", 
 *   SIAM J. Sci. & Stat. Comp. 18, 6, 1553-1570, (1997).
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

#ifdef HLLC_FLUX
#ifndef SPECIAL_RELATIVITY

#ifdef MHD
#error : The HLLC flux only works for hydro.
#endif /* MHD */

/*----------------------------------------------------------------------------*/
/*! \fn void fluxes(const Cons1DS Ul, const Cons1DS Ur,
 *            const Prim1DS Wl, const Prim1DS Wr,
 *            const Real Bxi, Cons1DS *pFlux)
 *  \brief Computes 1D fluxes
 *   Input Arguments:
 *   - Ul,Ur = L/R-states of CONSERVED variables at cell interface 
 *   Output Arguments:
 *   - pFlux = pointer to fluxes of CONSERVED variables at cell interface 
 */

void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr,
            const Real Bxi, Cons1DS *pFlux)
{
  Real sqrtdl,sqrtdr,isdlpdr,droe,v1roe,v2roe,v3roe;
#ifndef BAROTROPIC
  Real hroe;
#endif
  Real ev[NWAVE];
  Real *pFl, *pFr, *pF;
  Cons1DS Fl,Fr;
  int n;
  Real cfl,cfr,bp,bm,tmp;
  Real al,ar; /* Min and Max wave speeds */
  Real am,cp; /* Contact wave speed and pressure */
  Real tl,tr,dl,dr,sl,sm,sr;

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
 * Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
 * rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
 */

#ifndef ISOTHERMAL
  hroe = ((Ul.E + Wl.P)/sqrtdl + (Ur.E + Wr.P)/sqrtdr)*isdlpdr;
#endif

/*--- Step 3. ------------------------------------------------------------------
 * Compute eigenvalues using Roe-averaged values
 */

#ifdef ISOTHERMAL
  esys_roe_iso_hyd(v1roe, v2roe, v3roe, ev, NULL, NULL);
#else
  esys_roe_adb_hyd(v1roe, v2roe, v3roe, hroe, ev, NULL, NULL);
#endif /* ISOTHERMAL */

/*--- Step 4. ------------------------------------------------------------------
 * Compute the max and min wave speeds
 */

#ifdef ISOTHERMAL
  cfl = cfr = Iso_csound;
#else
  cfl = sqrt((double)(Gamma*Wl.P/Wl.d));
  cfr = sqrt((double)(Gamma*Wr.P/Wr.d));
#endif

  ar = MAX(ev[NWAVE-1],(Wr.Vx + cfr));
  al = MIN(ev[0]      ,(Wl.Vx - cfl));

  bp = ar > 0.0 ? ar : 0.0;
  bm = al < 0.0 ? al : 0.0;

/*--- Step 5. ------------------------------------------------------------------
 * Compute the contact wave speed and Pressure
 */

#ifdef ISOTHERMAL
  tl = Wl.d*Iso_csound2 + (Wl.Vx - al)*Ul.Mx;
  tr = Wr.d*Iso_csound2 + (Wr.Vx - ar)*Ur.Mx;
#else
  tl = Wl.P + (Wl.Vx - al)*Ul.Mx;
  tr = Wr.P + (Wr.Vx - ar)*Ur.Mx;
#endif

  dl =   Ul.Mx - Ul.d*al;
  dr = -(Ur.Mx - Ur.d*ar);

  tmp = 1.0/(dl + dr);
/* Determine the contact wave speed... */
  am = (tl - tr)*tmp;
/* ...and the pressure at the contact surface */
  cp = (dl*tr + dr*tl)*tmp;
  if(cp < 0.0) ath_perr(1,"[hllc flux]: Contact Pressure = %g\n",cp);
  cp = cp > 0.0 ? cp : 0.0;

/*--- Step 6. ------------------------------------------------------------------
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

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) {
    Fl.s[n] = Fl.d*Wl.r[n];
    Fr.s[n] = Fr.d*Wr.r[n];
  }
#endif

/*--- Step 7. ------------------------------------------------------------------
 * Compute flux weights or scales
 */

  if (am >= 0.0) {
    sl =  am/(am - bm);
    sr = 0.0;
    sm = -bm/(am - bm);
  }
  else {
    sl =  0.0;
    sr = -am/(bp - am);
    sm =  bp/(bp - am);
  }

/*--- Step 8. ------------------------------------------------------------------
 * Compute the HLLC flux at interface
 */
  pFl = (Real *)&(Fl);
  pFr = (Real *)&(Fr);
  pF  = (Real *)pFlux;
  for (n=0; n<(NWAVE+NSCALARS); n++) pF[n] = sl*pFl[n] + sr*pFr[n];

/* Add the weighted contribution of the flux along the contact */
  pFlux->Mx += sm*cp;
#ifndef ISOTHERMAL
  pFlux->E  += sm*cp*am;
#endif /* ISOTHERMAL */

#ifdef CYLINDRICAL
  if (al > 0.0) {
#ifndef ISOTHERMAL
    pFlux->Pflux = Wl.P;
#else /* ISOTHERMAL */
    pFlux->Pflux = Wl.d*Iso_csound2;
#endif /* ISOTHERMAL */
  }
  else if (ar < 0.0) {
#ifndef ISOTHERMAL
    pFlux->Pflux = Wr.P;
#else /* ISOTHERMAL */
    pFlux->Pflux = Wr.d*Iso_csound2;
#endif /* ISOTHERMAL */
  }
  else {
#ifndef ISOTHERMAL
    pFlux->Pflux = cp;
#else /* ISOTHERMAL */
    if (am >= 0.0) {
      pFlux->Pflux = Wl.d*(al-Wl.Vx)/(al-am);
    }
    else {
      pFlux->Pflux = Wr.d*(ar-Wr.Vx)/(ar-am);
    }
#endif /* ISOTHERMAL */
  }
#endif /* CYLINDRICAL */

  return;
}
#endif /* SPECIAL_RELATIVITY */
#endif /* HLLC_FLUX */
