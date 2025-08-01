#include "../copyright.h"
/*============================================================================*/
/*! \file two_shock.c
 *  \brief Computes 1D fluxes using simple two-shock Riemann solver.
 *
 * PURPOSE: Computes 1D fluxes using simple two-shock Riemann solver.
 *   Currently only isothermal hydrodynamics has been implemented.  
 *
 * REFERENCES:
 * - E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics",
 *   2nd ed., Springer-Verlag, Berlin, (1999).
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

#ifdef TWO_SHOCK_FLUX

#ifdef MHD
#error : The 2-shock flux for MHD has not been implemented.
#endif /* MHD */

#ifndef ISOTHERMAL
#error : The 2-shock flux for adiabatic EOS has not been implemented.
#endif /* ISOTHERMAL */

#if (NSCALARS > 0)
#error : Passive scalars have not been implemented in the 2-shock flux.
#endif /* NSCALARS */

/*----------------------------------------------------------------------------*/
/*! \fn void fluxes(const Cons1DS Ul, const Cons1DS Ur,
 *            const Prim1DS Wl, const Prim1DS Wr,
 *            const Real Bxi, Cons1DS *pFlux)
 *  \brief Computes 1D fluxes
 *   Input Arguments:
 *   - Bxi = B in direction of slice at cell interface
 *   - Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *   Output Arguments:
 *   - pFlux = pointer to fluxes of CONSERVED variables at cell interface
 */

void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr,
            const Real Bxi, Cons1DS *pFlux)
{
  Real zl, zc, zr, dc, Vxc, tmp;
  Real sl, sr;  /* Left and right going shock velocity */
  Real al, ar;  /* HLL a_l, a_r -> min and max signal velocity */

  if(!(Ul.d > 0.0)||!(Ur.d > 0.0))
    ath_error("[two-shock flux]: Non-positive densities: dl = %e  dr = %e\n",
	      Ul.d, Ur.d);

/*--- Step 1. ------------------------------------------------------------------
 * Convert left- and right- states in conserved to primitive variables.
 */

/*
  pbl = Cons1D_to_Prim1D(&Ul,&Wl,&Bxi);
  pbr = Cons1D_to_Prim1D(&Ur,&Wr,&Bxi);
*/

/*--- Step 2. ------------------------------------------------------------------
 * Compute the velocity and density of the contact
 */

  zl = sqrt((double)Wl.d);
  zr = sqrt((double)Wr.d);

  tmp = zl*zr*(Wl.Vx - Wr.Vx)/(2.0*Iso_csound*(zl + zr));
  zc = tmp + sqrt((double)(tmp*tmp + zl*zr));

  Vxc = (Wl.Vx*zl + Wr.Vx*zr)/(zl + zr) + Iso_csound*(zl - zr)/zc;

/* The L/R-going shock velocity */
  sl = Wl.Vx - Iso_csound*zc/zl;
  sr = Wr.Vx + Iso_csound*zc/zr;

/*--- Step 3. ------------------------------------------------------------------
 * For supersonic flow, return the upwind flux.
 */

  if(sr <= 0.0){
    pF->d  = Ur.Mx;
    pF->Mx = Ur.Mx*(Wr.Vx) + Wr.d*Iso_csound2;
    pF->My = Ur.My*(Wr.Vx);
    pF->Mz = Ur.Mz*(Wr.Vx);
    return;
  }

  if(sl >= 0.0){
    pF->d  = Ul.Mx;
    pF->Mx = Ul.Mx*(Wl.Vx) + Wl.d*Iso_csound2;
    pF->My = Ul.My*(Wl.Vx);
    pF->Mz = Ul.Mz*(Wl.Vx);
    return;
  }

/*--- Step 4. ------------------------------------------------------------------
 * Calculate the Interface Flux */

  dc = zc*zc;
  if(Vxc >= 0.0){
    pF->d  = dc*Vxc;
    pF->Mx = dc*Vxc*Vxc + dc*Iso_csound2;
    pF->My = dc*Vxc*Wl.Vy;
    pF->Mz = dc*Vxc*Wl.Vz;
  }
  else{
    pF->d  = dc*Vxc;
    pF->Mx = dc*Vxc*Vxc + dc*Iso_csound2;
    pF->My = dc*Vxc*Wr.Vy;
    pF->Mz = dc*Vxc*Wr.Vz;
  }

  return;
}
#endif /* TWO_SHOCK_FLUX */
