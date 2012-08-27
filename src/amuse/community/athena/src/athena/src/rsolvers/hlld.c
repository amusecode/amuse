#include "../copyright.h"
/*============================================================================*/
/*! \file hlld.c
 *  \brief Computes 1D fluxes using the HLLD Riemann solver.
 *
 * PURPOSE: Computes 1D fluxes using the HLLD Riemann solver, an extension of
 *   the HLLE solver to MHD.  Only works for MHD problems.  SEPARATE code
 *   blocks for adiabatic and isothermal equations of state.
 *
 * REFERENCES:
 * - T. Miyoshi & K. Kusano, "A multi-state HLL approximate Riemann solver
 *   for ideal MHD", JCP, 208, 315 (2005)
 * - A. Mignone, "A simple and accurate Riemann solver for isothermal MHD",
 *   JPC, 225, 1427 (2007)
 *
 * HISTORY: Adiabatic version written by Brian Biskeborn, May 8, 2006,
 *             COS sophmore project.
 *          Isothermal version written by Nicole Lemaster, May 1, 2008.
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

#ifdef HLLD_FLUX
#ifndef SPECIAL_RELATIVITY

#define SMALL_NUMBER 1e-8

#ifndef MHD
#error : The HLLD flux only works for mhd.
#endif /* MHD */

#ifndef ISOTHERMAL

/*----------------------------------------------------------------------------*/
/*! \fn void fluxes(const Cons1DS Ul, const Cons1DS Ur,
 *           const Prim1DS Wl, const Prim1DS Wr, const Real Bxi, Cons1DS *pFlux)
 *  \brief Compute 1D fluxes
 * Input Arguments:
 * - Bxi = B in direction of slice at cell interface
 * - Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *
 * Output Arguments:
 * - Flux = fluxes of CONSERVED variables at cell interface
 */

void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr, const Real Bxi, Cons1DS *pFlux)
{
  Cons1DS Ulst,Uldst,Urdst,Urst;       /* Conserved variable for all states */
  Prim1DS Wlst,Wrst;                   /* Primitive variables for all states */
  Cons1DS Fl,Fr;                       /* Fluxes for left & right states */
  Real spd[5];                        /* signal speeds, left to right */
/*  Real maxspd; */
  Real sdl,sdr,sdml,sdmr;             /* S_i-u_i, S_i-S_M (i=L or R) */
  Real pbl,pbr;                       /* Magnetic pressures */
  Real cfl,cfr,cfmax;                 /* Cf (left & right), max(cfl,cfr) */
  Real gpl,gpr,gpbl,gpbr;             /* gamma*P, gamma*P + B */
  Real sqrtdl,sqrtdr;                 /* sqrt of the L* & R* densities */
  Real invsumd;                       /* 1/(sqrtdl + sqrtdr) */
  Real ptl,ptr,ptst;                  /* total pressures */
  Real vbstl,vbstr;                   /* v_i* dot B_i* for i=L or R */
  Real Bxsig;                         /* sign(Bx) = 1 for Bx>0, -1 for Bx<0 */
  Real Bxsq = SQR(Bxi);               /* Bx^2 */
  Real tmp;                      /* Temp variable for repeated calculations */
#if (NSCALARS > 0)
  int n;
#endif



/*--- Step 1. ------------------------------------------------------------------
 * Convert left- and right- states in conserved to primitive variables.
 */

/*
  pbl = Cons1D_to_Prim1D(&Ul,&Wl,&Bxi);
  pbr = Cons1D_to_Prim1D(&Ur,&Wr,&Bxi);
*/

/*--- Step 2. ------------------------------------------------------------------
 * Compute left & right wave speeds according to Miyoshi & Kusano, eqn. (67)
 */

  pbl = 0.5*(SQR(Bxi) + SQR(Wl.By) + SQR(Wl.Bz));
  pbr = 0.5*(SQR(Bxi) + SQR(Wr.By) + SQR(Wr.Bz));
  gpl  = Gamma * Wl.P;
  gpr  = Gamma * Wr.P;
  gpbl = gpl + 2.0*pbl;
  gpbr = gpr + 2.0*pbr;

  cfl = sqrt((gpbl + sqrt(SQR(gpbl)-4*gpl*Bxsq))/(2.0*Wl.d));
  cfr = sqrt((gpbr + sqrt(SQR(gpbr)-4*gpr*Bxsq))/(2.0*Wr.d));
  cfmax = MAX(cfl,cfr);

  if(Wl.Vx <= Wr.Vx) {
    spd[0] = Wl.Vx - cfmax;
    spd[4] = Wr.Vx + cfmax;
  }
  else {
    spd[0] = Wr.Vx - cfmax;
    spd[4] = Wl.Vx + cfmax;
  }

/*  maxspd = MAX(fabs(spd[0]),fabs(spd[4])); */

/*--- Step 3. ------------------------------------------------------------------
 * Compute L/R fluxes
 */

  /* total pressure */
  ptl = Wl.P + pbl;
  ptr = Wr.P + pbr;

  Fl.d  = Ul.Mx;
  Fl.Mx = Ul.Mx*Wl.Vx + ptl - Bxsq;
  Fl.My = Ul.d*Wl.Vx*Wl.Vy - Bxi*Ul.By;
  Fl.Mz = Ul.d*Wl.Vx*Wl.Vz - Bxi*Ul.Bz;
  Fl.E  = Wl.Vx*(Ul.E + ptl - Bxsq) - Bxi*(Wl.Vy*Ul.By + Wl.Vz*Ul.Bz);
  Fl.By = Ul.By*Wl.Vx - Bxi*Wl.Vy;
  Fl.Bz = Ul.Bz*Wl.Vx - Bxi*Wl.Vz;

  Fr.d  = Ur.Mx;
  Fr.Mx = Ur.Mx*Wr.Vx + ptr - Bxsq;
  Fr.My = Ur.d*Wr.Vx*Wr.Vy - Bxi*Ur.By;
  Fr.Mz = Ur.d*Wr.Vx*Wr.Vz - Bxi*Ur.Bz;
  Fr.E  = Wr.Vx*(Ur.E + ptr - Bxsq) - Bxi*(Wr.Vy*Ur.By + Wr.Vz*Ur.Bz);
  Fr.By = Ur.By*Wr.Vx - Bxi*Wr.Vy;
  Fr.Bz = Ur.Bz*Wr.Vx - Bxi*Wr.Vz;

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) {
    Fl.s[n] = Fl.d*Wl.r[n];
    Fr.s[n] = Fr.d*Wr.r[n];
  }
#endif

/*--- Step 4. ------------------------------------------------------------------
 * Return upwind flux if flow is supersonic
 */

  if(spd[0] >= 0.0){
    *pFlux = Fl;
#if defined(CYLINDRICAL) && !defined(BAROTROPIC)
    pFlux->Pflux = ptl;
#endif 
    return;
  }

  if(spd[4] <= 0.0){
    *pFlux = Fr;
#if defined(CYLINDRICAL) && !defined(BAROTROPIC)
    pFlux->Pflux = ptr;
#endif 
    return;
  }

/*--- Step 5. ------------------------------------------------------------------
 * Compute middle and Alfven wave speeds
 */

  sdl = spd[0] - Wl.Vx;
  sdr = spd[4] - Wr.Vx;

  /* S_M: eqn (38) of Miyoshi & Kusano */
  spd[2] = (sdr*Wr.d*Wr.Vx - sdl*Wl.d*Wl.Vx - ptr + ptl) /
           (sdr*Wr.d-sdl*Wl.d);

  sdml   = spd[0] - spd[2];
  sdmr   = spd[4] - spd[2];
  /* eqn (43) of Miyoshi & Kusano */
  Ulst.d = Ul.d * sdl/sdml;
  Urst.d = Ur.d * sdr/sdmr;
  sqrtdl = sqrt(Ulst.d);
  sqrtdr = sqrt(Urst.d);

  /* eqn (51) of Miyoshi & Kusano */
  spd[1] = spd[2] - fabs(Bxi)/sqrtdl;
  spd[3] = spd[2] + fabs(Bxi)/sqrtdr;

/*--- Step 6. ------------------------------------------------------------------
 * Compute intermediate states
 */

  ptst = ptl + Ul.d*sdl*(sdl-sdml);
 
/* Ul* */
  /* eqn (39) of M&K */
  Ulst.Mx = Ulst.d * spd[2];
//   if((fabs(spd[2]/Wl.Vx-1.0)<SMALL_NUMBER) ||
//      (fabs(spd[2])/fabs(spd[0]) <= SMALL_NUMBER &&
//       fabs(Wl.Vx)/fabs(spd[0]) <= SMALL_NUMBER)) {
//     Ulst.My = Ulst.d * Wl.Vy;
//     Ulst.Mz = Ulst.d * Wl.Vz;
// 
//     Ulst.By = Ul.By;
//     Ulst.Bz = Ul.Bz;
//   }
  if (fabs(Ul.d*sdl*sdml/Bxsq-1.0) < SMALL_NUMBER) {
    /* Degenerate case */
    Ulst.My = Ulst.d * Wl.Vy;
    Ulst.Mz = Ulst.d * Wl.Vz;

    Ulst.By = Ul.By;
    Ulst.Bz = Ul.Bz;
  }
  else {
    /* eqns (44) and (46) of M&K */
    tmp = Bxi*(sdl-sdml)/(Ul.d*sdl*sdml-Bxsq);
    Ulst.My = Ulst.d * (Wl.Vy - Ul.By*tmp);
    Ulst.Mz = Ulst.d * (Wl.Vz - Ul.Bz*tmp);
//     if(Ul.By == 0.0 && Ul.Bz == 0.0) {
//       Ulst.By = 0.0;
//       Ulst.Bz = 0.0;
//     }
//     else {
//       /* eqns (45) and (47) of M&K */
//       tmp = (Ul.d*SQR(sdl)-Bxsq)/(Ul.d*sdl*sdml - Bxsq);
//       Ulst.By = Ul.By * tmp;
//       Ulst.Bz = Ul.Bz * tmp;
//     }

    /* eqns (45) and (47) of M&K */
    tmp = (Ul.d*SQR(sdl)-Bxsq)/(Ul.d*sdl*sdml - Bxsq);
    Ulst.By = Ul.By * tmp;
    Ulst.Bz = Ul.Bz * tmp;
  }
  vbstl = (Ulst.Mx*Bxi+Ulst.My*Ulst.By+Ulst.Mz*Ulst.Bz)/Ulst.d;
  /* eqn (48) of M&K */
  Ulst.E = (sdl*Ul.E - ptl*Wl.Vx + ptst*spd[2] +
            Bxi*(Wl.Vx*Bxi+Wl.Vy*Ul.By+Wl.Vz*Ul.Bz - vbstl))/sdml;
  Wlst = Cons1D_to_Prim1D(&Ulst,&Bxi);


/* Ur* */
  /* eqn (39) of M&K */
  Urst.Mx = Urst.d * spd[2];
//   if((fabs(spd[2]/Wr.Vx-1.0)<SMALL_NUMBER) ||
//      (fabs(spd[2])/fabs(spd[4]) <= SMALL_NUMBER &&
//       fabs(Wr.Vx)/fabs(spd[4]) <= SMALL_NUMBER)) {
//     Urst.My = Urst.d * Wr.Vy;
//     Urst.Mz = Urst.d * Wr.Vz;
// 
//     Urst.By = Ur.By;
//     Urst.Bz = Ur.Bz;
//   }
  if (fabs(Ur.d*sdr*sdmr/Bxsq-1.0) < SMALL_NUMBER) {
    /* Degenerate case */
    Urst.My = Urst.d * Wr.Vy;
    Urst.Mz = Urst.d * Wr.Vz;

    Urst.By = Ur.By;
    Urst.Bz = Ur.Bz;
  }
  else {
    /* eqns (44) and (46) of M&K */
    tmp = Bxi*(sdr-sdmr)/(Ur.d*sdr*sdmr-Bxsq);
    Urst.My = Urst.d * (Wr.Vy - Ur.By*tmp);
    Urst.Mz = Urst.d * (Wr.Vz - Ur.Bz*tmp);

//     if(Ur.By == 0.0 && Ur.Bz == 0.0) {
//       Urst.By = 0.0;
//       Urst.Bz = 0.0;
//     }
//     else {
//       /* eqns (45) and (47) of M&K */
//       tmp = (Ur.d*SQR(sdr)-Bxsq)/(Ur.d*sdr*sdmr - Bxsq);
//       Urst.By = Ur.By * tmp;
//       Urst.Bz = Ur.Bz * tmp;
//     }

    /* eqns (45) and (47) of M&K */
    tmp = (Ur.d*SQR(sdr)-Bxsq)/(Ur.d*sdr*sdmr - Bxsq);
    Urst.By = Ur.By * tmp;
    Urst.Bz = Ur.Bz * tmp;
  }
  vbstr = (Urst.Mx*Bxi+Urst.My*Urst.By+Urst.Mz*Urst.Bz)/Urst.d;
  /* eqn (48) of M&K */
  Urst.E = (sdr*Ur.E - ptr*Wr.Vx + ptst*spd[2] +
            Bxi*(Wr.Vx*Bxi+Wr.Vy*Ur.By+Wr.Vz*Ur.Bz - vbstr))/sdmr;
  Wrst = Cons1D_to_Prim1D(&Urst,&Bxi);


/* Ul** and Ur** - if Bx is zero, same as *-states */
//   if(Bxi == 0.0) {
  if(0.5*Bxsq/MIN(pbl,pbr) < SQR(SMALL_NUMBER)) {
    Uldst = Ulst;
    Urdst = Urst;
  }
  else {
    invsumd = 1.0/(sqrtdl + sqrtdr);
    if(Bxi > 0) Bxsig =  1;
    else        Bxsig = -1;

    Uldst.d = Ulst.d;
    Urdst.d = Urst.d;

    Uldst.Mx = Ulst.Mx;
    Urdst.Mx = Urst.Mx;

    /* eqn (59) of M&K */
    tmp = invsumd*(sqrtdl*Wlst.Vy + sqrtdr*Wrst.Vy + Bxsig*(Urst.By-Ulst.By));
    Uldst.My = Uldst.d * tmp;
    Urdst.My = Urdst.d * tmp;

    /* eqn (60) of M&K */
    tmp = invsumd*(sqrtdl*Wlst.Vz + sqrtdr*Wrst.Vz + Bxsig*(Urst.Bz-Ulst.Bz));
    Uldst.Mz = Uldst.d * tmp;
    Urdst.Mz = Urdst.d * tmp;

    /* eqn (61) of M&K */
    tmp = invsumd*(sqrtdl*Urst.By + sqrtdr*Ulst.By +
                   Bxsig*sqrtdl*sqrtdr*(Wrst.Vy-Wlst.Vy));
    Uldst.By = Urdst.By = tmp;

    /* eqn (62) of M&K */
    tmp = invsumd*(sqrtdl*Urst.Bz + sqrtdr*Ulst.Bz +
                   Bxsig*sqrtdl*sqrtdr*(Wrst.Vz-Wlst.Vz));
    Uldst.Bz = Urdst.Bz = tmp;

    /* eqn (63) of M&K */
    tmp = spd[2]*Bxi + (Uldst.My*Uldst.By + Uldst.Mz*Uldst.Bz)/Uldst.d;
    Uldst.E = Ulst.E - sqrtdl*Bxsig*(vbstl - tmp);
    Urdst.E = Urst.E + sqrtdr*Bxsig*(vbstr - tmp);
  }

/*--- Step 7. ------------------------------------------------------------------
 * Compute flux
 */

  if(spd[1] >= 0) {
/* return Fl* */
    pFlux->d  = Fl.d  + spd[0]*(Ulst.d  - Ul.d);
    pFlux->Mx = Fl.Mx + spd[0]*(Ulst.Mx - Ul.Mx);
    pFlux->My = Fl.My + spd[0]*(Ulst.My - Ul.My);
    pFlux->Mz = Fl.Mz + spd[0]*(Ulst.Mz - Ul.Mz);
    pFlux->E  = Fl.E  + spd[0]*(Ulst.E  - Ul.E);
    pFlux->By = Fl.By + spd[0]*(Ulst.By - Ul.By);
    pFlux->Bz = Fl.Bz + spd[0]*(Ulst.Bz - Ul.Bz);
  }
  else if(spd[2] >= 0) {
/* return Fl** */
    tmp = spd[1] - spd[0];
    pFlux->d  = Fl.d  - spd[0]*Ul.d  - tmp*Ulst.d  + spd[1]*Uldst.d;
    pFlux->Mx = Fl.Mx - spd[0]*Ul.Mx - tmp*Ulst.Mx + spd[1]*Uldst.Mx;
    pFlux->My = Fl.My - spd[0]*Ul.My - tmp*Ulst.My + spd[1]*Uldst.My;
    pFlux->Mz = Fl.Mz - spd[0]*Ul.Mz - tmp*Ulst.Mz + spd[1]*Uldst.Mz;
    pFlux->E  = Fl.E  - spd[0]*Ul.E  - tmp*Ulst.E  + spd[1]*Uldst.E;
    pFlux->By = Fl.By - spd[0]*Ul.By - tmp*Ulst.By + spd[1]*Uldst.By;
    pFlux->Bz = Fl.Bz - spd[0]*Ul.Bz - tmp*Ulst.Bz + spd[1]*Uldst.Bz;
  }
  else if(spd[3] > 0) {
/* return Fr** */
    tmp = spd[3] - spd[4];
    pFlux->d  = Fr.d  - spd[4]*Ur.d  - tmp*Urst.d  + spd[3]*Urdst.d;
    pFlux->Mx = Fr.Mx - spd[4]*Ur.Mx - tmp*Urst.Mx + spd[3]*Urdst.Mx;
    pFlux->My = Fr.My - spd[4]*Ur.My - tmp*Urst.My + spd[3]*Urdst.My;
    pFlux->Mz = Fr.Mz - spd[4]*Ur.Mz - tmp*Urst.Mz + spd[3]*Urdst.Mz;
    pFlux->E  = Fr.E  - spd[4]*Ur.E  - tmp*Urst.E  + spd[3]*Urdst.E;
    pFlux->By = Fr.By - spd[4]*Ur.By - tmp*Urst.By + spd[3]*Urdst.By;
    pFlux->Bz = Fr.Bz - spd[4]*Ur.Bz - tmp*Urst.Bz + spd[3]*Urdst.Bz;
  }
  else {
/* return Fr* */
    pFlux->d  = Fr.d  + spd[4]*(Urst.d  - Ur.d);
    pFlux->Mx = Fr.Mx + spd[4]*(Urst.Mx - Ur.Mx);
    pFlux->My = Fr.My + spd[4]*(Urst.My - Ur.My);
    pFlux->Mz = Fr.Mz + spd[4]*(Urst.Mz - Ur.Mz);
    pFlux->E  = Fr.E  + spd[4]*(Urst.E  - Ur.E);
    pFlux->By = Fr.By + spd[4]*(Urst.By - Ur.By);
    pFlux->Bz = Fr.Bz + spd[4]*(Urst.Bz - Ur.Bz);
  }

/* Fluxes of passively advected scalars, computed from density flux */
#if (NSCALARS > 0)
  if (pFlux->d >= 0.0) {
    for (n=0; n<NSCALARS; n++) pFlux->s[n] = pFlux->d*Wl.r[n];
  } else {
    for (n=0; n<NSCALARS; n++) pFlux->s[n] = pFlux->d*Wr.r[n];
  }
#endif

#if defined(CYLINDRICAL) && !defined(BAROTROPIC)
  pFlux->Pflux = ptst;
#endif
  return;
}

#else /* ISOTHERMAL */

/*----------------------------------------------------------------------------*/
/*! \fn void fluxes(const Cons1DS Ul, const Cons1DS Ur,
 *          const Prim1DS Wl, const Prim1DS Wr, const Real Bxi, Cons1DS *pFlux)
 *  \brief Compute 1D fluxes
 * Input Arguments:
 * - Bxi = B in direction of slice at cell interface
 * - Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *
 * Output Arguments:
 * - Flux = fluxes of CONSERVED variables at cell interface
 */

void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr, const Real Bxi, Cons1DS *pFlux)
{
  Cons1DS Ulst,Ucst,Urst;              /* Conserved variable for all states */
  Cons1DS Fl,Fr;                       /* Fluxes for left & right states */
  Real spd[5],idspd;                  /* signal speeds, left to right */
  Real pbl,pbr;                       /* Magnetic pressures */
  Real cfl,cfr;                       /* Cf (left & right) */
  Real gpl,gpr,gpbl,gpbr;             /* gamma*P, gamma*P + B */
  Real mxhll,ustar,dhll,sqrtdhll;     /* Mx, vel, and density in star states */
  Real fdhll,fmxhll;                  /* HLL fluxes (for star states) */
  Real ptl,ptr;                       /* total pressures */
  Real Bxsq = SQR(Bxi);               /* Bx^2 */
  Real tmp,mfact,bfact,X;             /* temporary variables */
  int n;

//   /* FROM ROE */
//   Real sqrtdl,sqrtdr,isdlpdr,idroe,v1roe,b2roe,b3roe,x,y;
//   Real bt_starsq,vaxsq,twid_csq,ct2,tsum,tdif,cf2_cs2,cfsq,cf;

/*--- Step 1. ------------------------------------------------------------------
 * Convert left- and right- states in conserved to primitive variables.
 */

/*
  pbl = Cons1D_to_Prim1D(&Ul,&Wl,&Bxi);
  pbr = Cons1D_to_Prim1D(&Ur,&Wr,&Bxi);
*/

/*--- Step 2. ------------------------------------------------------------------
 * Compute left & right wave speeds according to Mignone, eqns. (39), (8a),
 * and (9)
 */

  pbl = 0.5*(SQR(Bxi) + SQR(Wl.By) + SQR(Wl.Bz));
  pbr = 0.5*(SQR(Bxi) + SQR(Wr.By) + SQR(Wr.Bz));
  gpl  = Wl.d*Iso_csound2;
  gpr  = Wr.d*Iso_csound2;
  gpbl = gpl + 2.0*pbl;
  gpbr = gpr + 2.0*pbr;

  cfl = sqrt((gpbl + sqrt(SQR(gpbl)-4*gpl*Bxsq))/(2.0*Wl.d));
  cfr = sqrt((gpbr + sqrt(SQR(gpbr)-4*gpr*Bxsq))/(2.0*Wr.d));

  spd[0] = MIN(Wl.Vx-cfl,Wr.Vx-cfr);
  spd[4] = MAX(Wl.Vx+cfl,Wr.Vx+cfr);

//   /* COMPUTE ROE AVERAGES */
//   sqrtdl = sqrt((double)Wl.d);
//   sqrtdr = sqrt((double)Wr.d);
//   isdlpdr = 1.0/(sqrtdl + sqrtdr);
//   idroe  = 1.0/(sqrtdl*sqrtdr);
//   v1roe = (sqrtdl*Wl.Vx + sqrtdr*Wr.Vx)*isdlpdr;
//   b2roe = (sqrtdr*Wl.By + sqrtdl*Wr.By)*isdlpdr;
//   b3roe = (sqrtdr*Wl.Bz + sqrtdl*Wr.Bz)*isdlpdr;
//   x = 0.5*(SQR(Wl.By - Wr.By) + SQR(Wl.Bz - Wr.Bz))*SQR(isdlpdr);
//   y = 0.5*(Wl.d + Wr.d)*idroe;
// 
//   /* COMPUTE FAST MAGNETOSONIC SPEED */
//   bt_starsq = (SQR(b2roe) + SQR(b3roe))*y;
//   vaxsq = Bxsq*idroe;
//   twid_csq = Iso_csound2 + x;
//   ct2 = bt_starsq*idroe;
//   tsum = vaxsq + ct2 + twid_csq;
//   cfsq = 0.5*(tsum + sqrt((double)(SQR(tsum) - 4.0*twid_csq*vaxsq)));
//   cf = sqrt((double)cfsq);
// 
//   spd[0] = MIN(Wl.Vx-cfl,v1roe-cf);
//   spd[4] = MAX(v1roe+cf,Wr.Vx+cfr);

/*--- Step 3. ------------------------------------------------------------------
 * Compute L/R fluxes
 */

  /* total pressure */
  ptl = gpl + pbl;
  ptr = gpr + pbr;

  Fl.d  = Ul.Mx;
  Fl.Mx = Ul.Mx*Wl.Vx + ptl - Bxsq;
  Fl.My = Ul.d*Wl.Vx*Wl.Vy - Bxi*Ul.By;
  Fl.Mz = Ul.d*Wl.Vx*Wl.Vz - Bxi*Ul.Bz;
  Fl.By = Ul.By*Wl.Vx - Bxi*Wl.Vy;
  Fl.Bz = Ul.Bz*Wl.Vx - Bxi*Wl.Vz;

  Fr.d  = Ur.Mx;
  Fr.Mx = Ur.Mx*Wr.Vx + ptr - Bxsq;
  Fr.My = Ur.d*Wr.Vx*Wr.Vy - Bxi*Ur.By;
  Fr.Mz = Ur.d*Wr.Vx*Wr.Vz - Bxi*Ur.Bz;
  Fr.By = Ur.By*Wr.Vx - Bxi*Wr.Vy;
  Fr.Bz = Ur.Bz*Wr.Vx - Bxi*Wr.Vz;

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) {
    Fl.s[n] = Fl.d*Wl.r[n];
    Fr.s[n] = Fr.d*Wr.r[n];
  }
#endif

/*--- Step 4. ------------------------------------------------------------------
 * Return upwind flux if flow is supersonic
 */

  /* eqn. (38a) of Mignone */
  if(spd[0] >= 0.0){
    *pFlux = Fl;
    return;
  }

  /* eqn. (38e) of Mignone */
  if(spd[4] <= 0.0){
    *pFlux = Fr;
    return;
  }

/*--- Step 5. ------------------------------------------------------------------
 * Compute hll averages and Alfven wave speed
 */

  /* inverse of difference between right and left signal speeds */
  idspd = 1.0/(spd[4]-spd[0]);

  /* rho component of U^{hll} from Mignone eqn. (15);
   * uses F_L and F_R from eqn. (6) */
  dhll = (spd[4]*Ur.d-spd[0]*Ul.d-Fr.d+Fl.d)*idspd;
  sqrtdhll = sqrt(dhll);

  /* rho and mx components of F^{hll} from Mignone eqn. (17) */
  fdhll = (spd[4]*Fl.d-spd[0]*Fr.d+spd[4]*spd[0]*(Ur.d-Ul.d))*idspd;
  fmxhll = (spd[4]*Fl.Mx-spd[0]*Fr.Mx+spd[4]*spd[0]*(Ur.Mx-Ul.Mx))*idspd;

  /* ustar from paragraph between eqns. (23) and (24) */
  ustar = fdhll/dhll;

  /* mx component of U^{hll} from Mignone eqn. (15); paragraph referenced
   * above states that mxhll should NOT be used to compute ustar */
  mxhll = (spd[4]*Ur.Mx-spd[0]*Ul.Mx-Fr.Mx+Fl.Mx)*idspd;

  /* S*_L and S*_R from Mignone eqn. (29) */
  spd[1] = ustar - fabs(Bxi)/sqrtdhll;
  spd[3] = ustar + fabs(Bxi)/sqrtdhll;

/*--- Step 6. ------------------------------------------------------------------
 * Compute intermediate states
 */

/* Ul* */
  /* eqn. (20) of Mignone */
  Ulst.d = dhll;
  /* eqn. (24) of Mignone */
  Ulst.Mx = mxhll;

  tmp = (spd[0]-spd[1])*(spd[0]-spd[3]);
//   if (tmp == 0) {
  if ((fabs(spd[0]/spd[1]-1.0) < SMALL_NUMBER) 
        || (fabs(spd[0]/spd[3]-1.0) < SMALL_NUMBER)) {
    /* degenerate case described below eqn. (39) */
    Ulst.My = Ul.My;
    Ulst.Mz = Ul.Mz;
    Ulst.By = Ul.By;
    Ulst.Bz = Ul.Bz;
  } else {
    mfact = Bxi*(ustar-Wl.Vx)/tmp;
    bfact = (Ul.d*SQR(spd[0]-Wl.Vx)-Bxsq)/(dhll*tmp);

    /* eqn. (30) of Mignone */
    Ulst.My = dhll*Wl.Vy-Ul.By*mfact;
    /* eqn. (31) of Mignone */
    Ulst.Mz = dhll*Wl.Vz-Ul.Bz*mfact;
    /* eqn. (32) of Mignone */
    Ulst.By = Ul.By*bfact;
    /* eqn. (33) of Mignone */
    Ulst.Bz = Ul.Bz*bfact;
  }

/* Ur* */
  /* eqn. (20) of Mignone */
  Urst.d = dhll;
  /* eqn. (24) of Mignone */
  Urst.Mx = mxhll;

  tmp = (spd[4]-spd[1])*(spd[4]-spd[3]);
//   if (tmp == 0) {
  if ((fabs(spd[4]/spd[1]-1.0) < SMALL_NUMBER) 
        || (fabs(spd[4]/spd[3]-1.0) < SMALL_NUMBER)) {
    /* degenerate case described below eqn. (39) */
    Urst.My = Ur.My;
    Urst.Mz = Ur.Mz;
    Urst.By = Ur.By;
    Urst.Bz = Ur.Bz;
  } else {
    mfact = Bxi*(ustar-Wr.Vx)/tmp;
    bfact = (Ur.d*SQR(spd[4]-Wr.Vx)-Bxsq)/(dhll*tmp);

    /* eqn. (30) of Mignone */
    Urst.My = dhll*Wr.Vy-Ur.By*mfact;
    /* eqn. (31) of Mignone */
    Urst.Mz = dhll*Wr.Vz-Ur.Bz*mfact;
    /* eqn. (32) of Mignone */
    Urst.By = Ur.By*bfact;
    /* eqn. (33) of Mignone */
    Urst.Bz = Ur.Bz*bfact;
  }

/* Uc* */
  /* from below eqn. (37) of Mignone */
  X = sqrtdhll*SIGN(Bxi);
  /* eqn. (20) of Mignone */
  Ucst.d = dhll;
  /* eqn. (24) of Mignone */
  Ucst.Mx = mxhll;
  /* eqn. (34) of Mignone */
  Ucst.My = 0.5*(Ulst.My+Urst.My+X*(Urst.By-Ulst.By));
  /* eqn. (35) of Mignone */
  Ucst.Mz = 0.5*(Ulst.Mz+Urst.Mz+X*(Urst.Bz-Ulst.Bz));
  /* eqn. (36) of Mignone */
  Ucst.By = 0.5*(Ulst.By+Urst.By+(Urst.My-Ulst.My)/X);
  /* eqn. (37) of Mignone */
  Ucst.Bz = 0.5*(Ulst.Bz+Urst.Bz+(Urst.Mz-Ulst.Mz)/X);

/*--- Step 7. ------------------------------------------------------------------
 * Compute flux
 */

  if(spd[1] >= 0) {
/* return (Fl+Sl*(Ulst-Ul)), eqn. (38b) of Mignone */
    pFlux->d  = Fl.d  + spd[0]*(Ulst.d  - Ul.d);
    pFlux->Mx = Fl.Mx + spd[0]*(Ulst.Mx - Ul.Mx);
    pFlux->My = Fl.My + spd[0]*(Ulst.My - Ul.My);
    pFlux->Mz = Fl.Mz + spd[0]*(Ulst.Mz - Ul.Mz);
    pFlux->By = Fl.By + spd[0]*(Ulst.By - Ul.By);
    pFlux->Bz = Fl.Bz + spd[0]*(Ulst.Bz - Ul.Bz);
  }
  else if (spd[3] <= 0) {
/* return (Fr+Sr*(Urst-Ur)), eqn. (38d) of Mignone */
    pFlux->d  = Fr.d  + spd[4]*(Urst.d  - Ur.d);
    pFlux->Mx = Fr.Mx + spd[4]*(Urst.Mx - Ur.Mx);
    pFlux->My = Fr.My + spd[4]*(Urst.My - Ur.My);
    pFlux->Mz = Fr.Mz + spd[4]*(Urst.Mz - Ur.Mz);
    pFlux->By = Fr.By + spd[4]*(Urst.By - Ur.By);
    pFlux->Bz = Fr.Bz + spd[4]*(Urst.Bz - Ur.Bz);
  }
  else {
/* return Fcst, eqn. (38c) of Mignone, using eqn. (24) */
    pFlux->d = dhll*ustar;
    pFlux->Mx = fmxhll;
    pFlux->My = Ucst.My*ustar - Bxi*Ucst.By;
    pFlux->Mz = Ucst.Mz*ustar - Bxi*Ucst.Bz;
    pFlux->By = Ucst.By*ustar - Bxi*Ucst.My/Ucst.d;
    pFlux->Bz = Ucst.Bz*ustar - Bxi*Ucst.Mz/Ucst.d;
  }

/* Fluxes of passively advected scalars, computed from density flux */
#if (NSCALARS > 0)
  if (pFlux->d >= 0.0) {
    for (n=0; n<NSCALARS; n++) pFlux->s[n] = pFlux->d*Wl.r[n];
  } else {
    for (n=0; n<NSCALARS; n++) pFlux->s[n] = pFlux->d*Wr.r[n];
  }
#endif

  return;
}

#endif /* ISOTHERMAL */
#endif /* SPECIAL_RELATIVITY */
#endif /* HLLD_FLUX */
