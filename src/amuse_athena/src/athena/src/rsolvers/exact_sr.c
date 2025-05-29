#include "../copyright.h"
/*============================================================================*/
/*! \file exact_sr.c
 *  \brief Computes 1D fluxes using exact special relativistic Riemann solver.
 *
 * PURPOSE: Computes 1D fluxes using exact special relativistic Riemann solver.
 *   currently works only for hydrodynamics
 *
 * REFERENCES:
 * - Rezzolla, Zanotti, and Pons. "An Improved Exact Riemann Solver for
 *   Multidimensional Relativistic Flows." 2002. 
 *
 * HISTORY:
 * - April-2010:  Written by Nick Hand. 
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - fluxes() - all Riemann solvers in Athena must have this function name and
 *              use the same argument list as defined in rsolvers/prototypes.h
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef EXACT_FLUX
#ifdef SPECIAL_RELATIVITY

#ifdef MHD
#error : The exact flux for MHD has not been implemented.
#endif /* MHD */

#if (NSCALARS > 0)
#error : Passive scalars have not been implemented in the exact flux.
#endif /* NSCALARS */

#define EPS 3.0e-11 /* EPS is the relative precision. */
#define JMAX 40     /* Maximum number of allowed bisections */
enum WaveType {Two_S, RS, SR, Two_R, None}; 

static enum WaveType wave; /* the wave pattern of the Riemann problem */ 

void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr, const Real Bx, Cons1DS *pF);

void gauleg(double x1, double x2, Real x[], Real w[], int n);

Real rtbis_vel(Real (*func)(const Prim1DS Wl, const Prim1DS Wr, Real p), const Prim1DS Wl, 
	   const Prim1DS Wr, Real x1, Real x2, Real xacc);  

Real getDelVRel(const Prim1DS Wl, const Prim1DS Wr, Real p);

Real getVRel_2R(const Prim1DS Wl, const Prim1DS Wr, Real p); 

Real getVRel_RS(const Prim1DS Wl, const Prim1DS Wr, Real p); 

Real getVRel_2S(const Prim1DS Wl, const Prim1DS Wr, Real p);

Real getVlim_2S(const Prim1DS Wl, const Prim1DS Wr);

Real getVlim_RS(const Prim1DS Wl, const Prim1DS Wr); 

Real getVlim_2R(const Prim1DS Wl, const Prim1DS Wr);

Real integrateRaref(const Prim1DS W, Real a, Real b);

Real getP(const Prim1DS Wl, const Prim1DS Wr); 

void getShockVars(const Prim1DS Wa, Real Pb, char dir, Real *pJ, Real *pV_shock, Real *pd);
 
Real getVb_Shock(const Prim1DS Wa, Real Pb, char dir);

Real getVb_Raref(const Prim1DS Wa, Real Pb, char dir);

Real getXi(const Prim1DS Wa, Real P, Real vxc, char dir);

Real rtbis_xi(Real (*func)(const Prim1DS Wa, Real P, Real vx, char), const Prim1DS Wa, 
		 char dir, Real f_x1, Real f_x2, Real x1, Real x2, Real xacc);

void getVelT_Raref(const Prim1DS Wa, Real P, Real vxb, Real *pVy, Real *pVz);

void getVelT_Shock(const Prim1DS Wa, Real P, Real vxb, Real *pVy, Real *pVz);

void setFluxes(Real vx, Real vy, Real vz, Real P, Real d,Cons1DS *pF );

/*--------------------------------------------------------------------------*/
/*! \fn void fluxes(const Cons1DS Ul, const Cons1DS Ur,
 *            const Prim1DS Wl, const Prim1DS Wr, const Real Bx, Cons1DS *pF)
 *  \brief Computes 1D fluxes using exact special relativistic Riemann solver.
 *
 *   Input Arguments:
 *   - Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *   - Bx = B in direction of slice at cell interface
 *   Output Arguments:
 *   - pF = pointer to fluxes of CONSERVED variables at cell interface
 *
 *--------------------------------------------------------------------------*/
void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr, const Real Bx, Cons1DS *pF)
{
  Real pc, vxc, dcl, dcr; /* pressure, normal velocity, density in center region*/
  Real vl_shock, vr_shock;  /* left/right shock velocity */
  Real vr_hd, vr_tl, vl_hd, vl_tl;  /* raref head/tail velocities */
  Real eps, xacc; 
  Real vt, p, vx, vy, vz, C, d; 
  Real signY = 1.0; 
  Real signZ = 1.0;
  Real TOL = 1e-5;
 
  /* check if intial data states are the same */
  if ((fabs(Wl.P-Wr.P) <= TOL)  && (fabs(Wl.Vx-Wr.Vx) <= TOL)) {  
    if (vxc >= 0.0) {
      setFluxes(Wl.Vx, Wl.Vy, Wl.Vz, Wl.P, Wl.d, pF);
      return; 
    }
    else {
      setFluxes(Wl.Vx, Wr.Vy, Wr.Vz, Wl.P, Wr.d, pF);
      return;
    }
  }

  /* get the pressure in the intermediate state */
  pc = getP(Wl, Wr);

  /* determine density and velocity behind waves */
  if (wave == RS) {
    vxc = getVb_Shock(Wr, pc, 'R');
    getShockVars(Wr, pc, 'R', NULL, &vr_shock, &dcr); 
    dcl = Wl.d*pow((double)(pc/Wl.P), (double) (1.0/Gamma));
  }
  else if (wave == SR) {
    vxc = getVb_Shock(Wl, pc, 'L'); 
    getShockVars(Wl, pc, 'L', NULL, &vl_shock, &dcl); 
    dcr = Wr.d*pow((double)(pc/Wr.P), (double) (1.0/Gamma));
  }
  else if (wave == Two_R) {
    vxc = getVb_Raref(Wl, pc, 'L'); 
    dcl = Wl.d*pow((double)(pc/Wl.P), (double) (1.0/Gamma));
    dcr = Wr.d*pow((double)(pc/Wr.P), (double) (1.0/Gamma));
  }
  else if (wave == Two_S) {
    vxc = getVb_Shock(Wl, pc, 'L');
    getShockVars(Wl, pc, 'L', NULL, &vl_shock, &dcl); 
    getShockVars(Wr, pc, 'R', NULL, &vr_shock, &dcr);    
  }
  
/*-----------------------------------------------------------------
* Calculate the interface flux if the wave speeds are such that we aren't
* actually in the intermediate state */
  
  if (pc > Wl.P) {
    /* left shock wave */

    if (vl_shock >= 0.0) {
      /* to left of shock */   
      setFluxes(Wl.Vx, Wl.Vy, Wl.Vz, Wl.P, Wl.d, pF);
      return; 
    }
  }
  else {
    /* left rarefaction */
    
   /* calculate velocity at head and tail of raref */
    vl_hd = getXi(Wl, Wl.P, Wl.Vx, 'L');  
    vl_tl = getXi(Wl, pc, vxc, 'L'); 
    
    if (vl_hd >= 0.0) {
      /* To left of rarefaction */
      setFluxes(Wl.Vx, Wl.Vy, Wl.Vz, Wl.P, Wl.d, pF); 
      return;
    } 
    else if (vl_tl >= 0.0) {
      /* Inside rarefaction fan */
      
      /* machine precision */ 
      eps = 1.0f; 
      do eps /= 2.0f; 
      while ((float)(1.0 + (eps/2.0)) != 1.0);
      xacc = eps*0.5*(Wl.P + pc); 

      /* get vx, p, d in raref fan */
      p = rtbis_xi(getXi, Wl, 'L', vl_hd, vl_tl, Wl.P, pc, xacc); 
      vx = getVb_Raref(Wl, p, 'L'); 
      d = Wl.d*pow((double)(p/Wl.P), (double) (1.0/Gamma));
     
      getVelT_Raref(Wl, p, vx, &vy, &vz);     
      setFluxes(vx, vy, vz, p, d, pF); 
      return; 
    } 
  }

  if (pc > Wr.P) {
    /* right shock wave */
    
    /* right shock speed */
   
    if (vr_shock <= 0.0) {
      /* to right of shock */

      setFluxes(Wr.Vx, Wr.Vy, Wr.Vz, Wr.P, Wr.d, pF); 
      return; 
    }
  }
  else {
    /* right rarefaction */
   
    /* calculate velocity at head and tail of raref */
    vr_hd = getXi(Wr, Wr.P, Wr.Vx, 'R');  
    vr_tl = getXi(Wr, pc, vxc, 'R'); 
    
    if (vr_hd <= 0.0) {
      /* To right of rarefaction */
      
      setFluxes(Wr.Vx, Wr.Vy, Wr.Vz, Wr.P, Wr.d, pF); 
      return;
      
    }
    else if (vr_tl <= 0.0) {
      /* Inside rarefaction fan */
      
      /* machine precision */ 
      eps = 1.0f; 
      do eps /= 2.0f; 
      while ((float)(1.0 + (eps/2.0)) != 1.0); 
      xacc = eps*0.5*(Wr.P + pc); 

      /* get vx, p, d in raref fan */
      p = rtbis_xi(getXi, Wr, 'R', vr_hd, vr_tl, Wr.P, pc, xacc); 
      vx = getVb_Raref(Wr, p, 'R'); 
      d = Wr.d*pow((double)(p/Wr.P), (double) (1.0/Gamma));
     
      getVelT_Raref(Wr, p, vx, &vy, &vz);  
      setFluxes(vx, vy, vz, p, d, pF); 
      return;
    }
  }

/* We are in the intermediate state */
/*---------------------------------------------------------------------
 * Calculate the interface flux */

  if (vxc >= 0.0) {
    
    if (Wl.P > pc) {
      /* left rarefaction */
      getVelT_Raref(Wl, pc, vxc, &vy, &vz);
      setFluxes(vxc, vy, vz, pc, dcl, pF); 
      return;
    }
    else {
      /*left shock */
      getVelT_Shock(Wl, pc, vxc, &vy, &vz);
      setFluxes(vxc, vy, vz, pc, dcl, pF); 
      return;
    }
   
  }
  else {
    
    if (pc > Wr.P) {
      /* right shock */
      getVelT_Shock(Wr, pc, vxc, &vy, &vz);    
      setFluxes(vxc, vy, vz, pc, dcr, pF); 
      return;
    }
    else {
      /* right rarefaction */
      
      getVelT_Raref(Wr, pc, vxc, &vy, &vz);
      setFluxes(vxc, vy, vz, pc, dcr, pF); 
      return;
    }
  }
}
/*----------------------------------------------------------------------------*/
/*! \fn Real integrateRaref(const Prim1DS Wa, Real a, Real b)
 *  \brief Numerically integrate the integral for rarefactions (eq 3.22 in  RZP)
 * using Gaussian-Legendre quadrature */
Real integrateRaref(const Prim1DS Wa, Real a, Real b)
{
  
  int i;
  Real integral, sum, diff, f, xx; 
  Real * x; 
  Real * w; 
  int n = 10; /* number of integration points */
  Real k, h, vt, G, A;
  Real dd, ccs2, hh;


  /* allocate memory for absicssas and weights */
  x = (Real *)malloc((size_t) ((n+1)*sizeof(Real)));
  w = (Real *)malloc((size_t) ((n+1)*sizeof(Real)));
 
  /* get the abscissas and weights for integration */
  gauleg(-1.0, 1.0, x, w, n);

  sum = 0.5 * (b + a); 
  diff = 0.5 * (b - a); 
  integral = 0.0; 

  for (i = 1; i <= n; i++) {
    xx = diff*x[i] + sum; 

    h = 1.0 + Gamma*Wa.P/((Gamma-1.0)*Wa.d); 
    vt = sqrt((double)(Wa.Vy*Wa.Vy + Wa.Vz*Wa.Vz)); 
  
    G = 1.0/sqrt((double)(1.0 - Wa.Vx*Wa.Vx - Wa.Vy*Wa.Vy - Wa.Vz*Wa.Vz)); 
    A = h*G*vt;

    dd = Wa.d*pow((double)(xx/Wa.P), (double) (1.0/Gamma)); 
    ccs2 = Gamma*(Gamma - 1.)*xx / (Gamma*xx + (Gamma - 1.)*dd); 
    hh = 1.0 + xx*Gamma/(dd*(Gamma - 1.0)); 
  
    f = sqrt((double) (hh*hh + A*A*(1.0-ccs2)))/(dd*sqrt((double)ccs2)*(hh*hh + A*A));

    integral = integral + diff*w[i]*f; 
  }
  
  free(x); 
  free(w); 

  return integral; 
}

/*----------------------------------------------------------------------------*/
/*! \fn void getShockVars(const Prim1DS Wa, Real Pb, char dir, Real *pJ, 
 *                        Real *pV_s, Real *pD)
 *  \brief Determine the mass flux, shock velocity and density across a shock 
 *  wave */
void getShockVars(const Prim1DS Wa, Real Pb, char dir, Real *pJ, Real *pV_s, 
                  Real *pD)
{
  Real sign, ha, Ga;
  Real db, hb;  /* specific enthalpy and density behind wave */
  Real A, B, C,D; 
  Real d, J, v_s;
  Real TOL = 1e-5;

  if (dir == 'L') sign = -1.0; 
  if (dir == 'R') sign = 1.0;

  /* check if pressure is constant across wave */
  if (fabs(Wa.P - Pb) <= TOL) { 
    if (pJ != NULL)
      *pJ = 0.0; 
    if (pD != NULL)
      *pD = Wa.d;
    if (pV_s != NULL)
      *pV_s = Wa.Vx;
    
    return;
  }
  
  ha = 1.0 + Gamma*Wa.P/((Gamma-1.0)*Wa.d);
  Ga = 1.0/sqrt((double)(1.0 - Wa.Vx*Wa.Vx - Wa.Vy*Wa.Vy - Wa.Vz*Wa.Vz)); 

  /* calculate specific enthalpy */ 
  A = 1.0 + (Gamma - 1.0)*(Wa.P-Pb)/(Gamma*Pb); 
  B = 1.0 - A;
  C = ha*(Wa.P - Pb)/Wa.d - ha*ha; 

  /*check for unphysical specific enthalpies*/
  if (C > (B*B/(4.0*A)))
    ath_error("[exact flux]: Unphysical specific enthalpy in intermediate state");
    
  hb = (-B + sqrt((double)(B*B - 4.0*A*C)))/(2.0*A); 

  d = Gamma*Pb/((Gamma - 1.0)*(hb - 1.0));  
  if (pD != NULL)
    *pD = d; 
  
  J =  sign*sqrt((double)((Pb - Wa.P)/(ha/Wa.d - hb/d)));
  if (pJ != NULL)
    *pJ = J;
 
  A = Wa.d*Wa.d*Ga*Ga;
  if (pV_s != NULL) 
    *pV_s = (A*Wa.Vx +sign*fabs((float) J)*
	     sqrt((double)(J*J + A*(1.0 - Wa.Vx*Wa.Vx))))/(A + J*J);
  
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn Real getVlim_2R(const Prim1DS Wl, const Prim1DS Wr)
 *  \brief Determine vlim for two rarefactions (eq. 4.12 in RZP) */
Real getVlim_2R(const Prim1DS Wl, const Prim1DS Wr)
{
  
  Real vxlc, vxrc;  
  Real integral; 
  Real vlb, vrb; 

  /* compute vxl,c from RZP eq 4.13 */
  vlb = getVb_Raref(Wl, 0.0, 'L');
  vxlc = (Wl.Vx - vlb)/(1.0 - Wl.Vx*vlb);

  /* compute vxr,c from RZP eq 4.14 */
  vrb = getVb_Raref(Wr, 0.0, 'R');
  vxrc = (Wr.Vx - vrb)/(1.0 - Wr.Vx*vrb);

  return (vxlc - vxrc) / (1.0 - vxlc*vxrc); 
  
}

/*----------------------------------------------------------------------------*/
/*! \fn Real getVlim_RS(const Prim1DS Wl, const Prim1DS Wr)
 *  \brief Determine vlim for the case of rarefaction-shock (eq 4.10) */
Real getVlim_RS(const Prim1DS Wl, const Prim1DS Wr)
{
 
  Real a, b, integral; 
  Real vlb, vrb, vlc, vrc; 
  
  if (Wl.P < Wr.P) { 
    /* shock to left and rarefaction to right */
    /* limiting p in center is Wl.P */
    vlb = Wl.Vx;
    vrb = getVb_Raref(Wr, Wl.P, 'R');
  }
  else {
    /* shock to right and rarefaction to left */
    /* limiting p in center is Wr.P */
    vlb = getVb_Raref(Wl, Wr.P, 'L');
    vrb = Wr.Vx;
  }
     
  vlc = (Wl.Vx - vlb)/(1.0 - Wl.Vx*vlb);
  vrc = (Wr.Vx - vrb)/(1.0 - Wr.Vx*vrb);
    
  return (vlc - vrc)/(1.0 - vlc*vrc);
} 

/*--------------------------------------------------------------------------*/
/*! \fn Real getVlim_2S(const Prim1DS Wl, const Prim1DS Wr)
 *  \brief Determine vlim for the case of two shocks (eq 4.5) */
Real getVlim_2S(const Prim1DS Wl, const Prim1DS Wr)
{  
  if (Wl.P > Wr.P) {
    Real vlc, vrc, vrb; 

    vrb = getVb_Shock(Wr, Wl.P, 'R');  
    vrc = (Wr.Vx - vrb)/(1.0 - vrb*Wr.Vx); 
    vlc = 0.0; 
    return (vlc - vrc)/(1.0 - vlc*vrc);  
  }
  else {
    Real vlc, vrc, vlb; 

    vlb = getVb_Shock(Wl, Wr.P, 'L'); 
    vlc = (Wl.Vx - vlb)/(1.0 - vlb*Wl.Vx); 
    vrc = 0.0; 
    return (vlc - vrc)/(1.0 - vlc*vrc);
    
  }
}

/*----------------------------------------------------------------------------*/       
/*! \fn Real getDelVRel(const Prim1DS Wl, const Prim1DS Wr, Real p)
 *  \brief Determine vRel for the wave pattern indicated by wave */
Real getDelVRel(const Prim1DS Wl, const Prim1DS Wr, Real p)
{
  Real vRel; 
  Real vRel_0 = (Wl.Vx - Wr.Vx)/(1.0 - Wl.Vx*Wr.Vx); 
  
  if (wave == Two_R) {
    vRel = getVRel_2R(Wl, Wr, p); 
    return vRel - vRel_0; 
  }
  
  else if (wave == RS || wave == SR) {
    vRel = getVRel_RS(Wl, Wr, p);
    return vRel - vRel_0; 
  }
  else if (wave == Two_S) {
    vRel = getVRel_2S(Wl, Wr, p);
    return vRel - vRel_0; 
  }

  return -1.0; /* should never reach here*/
}
/*----------------------------------------------------------------------------*/
/*! \fn Real getP(const Prim1DS Wl, const Prim1DS Wr)
 *  \brief Determine the pressure in the intermediate state using interval 
 *  bisection */
Real getP(const Prim1DS Wl, const Prim1DS Wr)
{
  Real vRel, vRS, vSS, vRR; 
  Real pmin, pmax;  /* brackets for pressure */ 
  Real xacc; 
  Real pc;  /* pressure in intermediate state*/
  float tol, eps; 

  /* 1. Compute the relative x-velocity between Wl and Wr */
  vRel = (Wl.Vx - Wr.Vx) / (1.0 - Wl.Vx*Wr.Vx); 
  
  /* 2. Determine the wave pattern by comparing vRel with the limiting values of
     the relative velocity and bracket pressure in intermed state */
  vSS = getVlim_2S(Wl, Wr); 
  vRS = getVlim_RS(Wl, Wr); 

  if (vRel <= vRS) {
    wave = Two_R;
    pmin = 0.0; 
    if (Wl.P < Wr.P) pmax = Wl.P; 
    else pmax = Wr.P;
  }
  else if (vRel > vRS && vRel <= vSS) {
    if (Wl.P > Wr.P) {
      wave = RS; 
      pmin = Wr.P;
      pmax = Wl.P;
    }
    else { 
      wave = SR; 
      pmin = Wl.P; 
      pmax = Wr.P;
    } 
  }
  else {
    wave = Two_S; 
    if (Wl.P > Wr.P) pmin = Wl.P;
    else pmin = Wr.P;
    pmax = 1000 * 0.5*(Wl.P + Wr.P); 
  }

  /* 3. compute the machine precision */
  eps = 1.0f; 
  do eps /= 2.0f; 
  while ((float)(1.0 + (eps/2.0)) != 1.0); 
  xacc = eps*0.5*(pmin + pmax); 

  /* 4. use interval bisection to find p in the intermediate state */
  pc = rtbis_vel(getDelVRel, Wl, Wr, pmin, pmax, xacc); 
  return pc; 
  
}

/*---------------------------------------------------------------------------*/
/*! \fn Real getVRel_2R(const Prim1DS Wl, const Prim1DS Wr, Real p)
 *  \brief Determine vRel for the case of two rarefactions */
Real getVRel_2R(const Prim1DS Wl, const Prim1DS Wr, Real p)
{

  Real vlb, vrb; /* the normal velocity behind the left and right rarefactions */
  Real integral; 
  Real vlc, vrc; /* normal velocity of left/right waves in frame of contact */

  /* determine x-velocity behind left/right waves */
  vlb = getVb_Raref(Wl, p, 'L'); 
 
  vrb = getVb_Raref(Wr, p, 'R'); 

  /* determine velocities in frame of contact discontinuity */ 
  vlc = (Wl.Vx - vlb) / (1.0 - Wl.Vx*vlb); 
  vrc = (Wr.Vx - vrb) / (1.0 - Wr.Vx*vrb); 

  return (vlc - vrc)/(1.0 - vrc*vlc); 
}

/*---------------------------------------------------------------------------*/
/*! \fn Real getVRel_RS(const Prim1DS Wl, const Prim1DS Wr, Real p)
 *  \brief Determine vRel for the case of rarefaction-shock wave */
Real getVRel_RS(const Prim1DS Wl, const Prim1DS Wr, Real p)
{
  Real vlb, vrb; /* normal velocity behind the waves */
  Real vlc, vrc; /* normal velocity in frame of contact */
  Real v_shock, G_shock, J, ha, Ga; 
  Real num, denom, integral; 

  if (wave == SR) {
    /* shock wave moving to left and rarefaction to right */
    
    /* calculate vlb for shock*/
    vlb = getVb_Shock(Wl, p, 'L'); 

    /*calculate vrb for rarefaction */
    vrb = getVb_Raref(Wr, p, 'R');

    vlc = (Wl.Vx - vlb) / (1.0 - Wl.Vx*vlb); 
    vrc = (Wr.Vx - vrb) / (1.0 - Wr.Vx*vrb); 

    return (vlc - vrc)/(1.0 - vrc*vlc); 
  }
  else {
    /* shock wave moving to right and rarefaction to left */
    
    /* calculate vrb for shock  */
    vrb = getVb_Shock(Wr, p, 'R'); 

    /*calculate vlb for raref*/
    vlb = getVb_Raref(Wl, p, 'L'); 

    vlc = (Wl.Vx - vlb) / (1.0 - Wl.Vx*vlb); 
    vrc = (Wr.Vx - vrb) / (1.0 - Wr.Vx*vrb); 

    return (vlc - vrc)/(1.0 - vrc*vlc); 
  } 
}
/*---------------------------------------------------------------------------*/
/*! \fn Real getVRel_2S(const Prim1DS Wl, const Prim1DS Wr, Real p)
 *  \brief Determine vRel for the case of two shocks */
Real getVRel_2S(const Prim1DS Wl, const Prim1DS Wr, Real p)
{
  Real vlb, vrb; /* normal velocity behind the shock waves */
  Real vlc, vrc; /* normal velocity in frame of contact */
  char dir;
  Real v_shock, G_shock, J, ha, Ga; 
  Real num, denom; 

  /* calculate vlb */
  vlb = getVb_Shock(Wl, p, 'L');

  /* calculate vrb */ 
  vrb = getVb_Shock(Wr, p, 'R'); 

  /* determine velocities in frame of contact discontinuity */ 
  vlc = (Wl.Vx - vlb) / (1.0 - Wl.Vx*vlb); 
  vrc = (Wr.Vx - vrb) / (1.0 - Wr.Vx*vrb);

  return (vlc - vrc)/(1.0 - vlc*vrc); 
}
/*---------------------------------------------------------------------------*/
/*! \fn Real getVb_Shock(const Prim1DS Wa, Real Pb, char dir) 
 *  \brief Determine the normal velocity behind a shock wave */
Real getVb_Shock(const Prim1DS Wa, Real Pb, char dir)
{
  Real num, denom;
  Real G_shock, ha, Ga;
  Real J, v_shock; 
  
  getShockVars(Wa, Pb, dir, &J, &v_shock, NULL); 
  if (fabs(J) == 0.0) return Wa.Vx;

  G_shock = 1.0 / sqrt((double) (1.0 - v_shock*v_shock));
  
  /* specific enthalpy and gamma factor for left state */
  ha = 1.0 + Gamma*Wa.P/((Gamma-1.0)*Wa.d);
  Ga = 1.0/sqrt((double)(1.0 - Wa.Vx*Wa.Vx - Wa.Vy*Wa.Vy - Wa.Vz*Wa.Vz)); 

  num = ha*Ga*Wa.Vx + G_shock*(Pb - Wa.P)/J; 
  denom = ha*Ga + (Pb - Wa.P)*(G_shock*Wa.Vx/J + 1.0/(Wa.d*Ga));  

  return num/denom; 
}
/*----------------------------------------------------------------------------*/
/*! \fn Real getVb_Raref(const Prim1DS Wa, Real Pb, char dir)
 *  \brief Determine the normal velocity behind a rarefaction wave */
Real getVb_Raref(const Prim1DS Wa, Real Pb, char dir)
{
  Real integral, sign, B; 
  
  if (dir == 'L') sign = -1.0; 
  if (dir == 'R') sign = 1.0; 
  
  integral = integrateRaref(Wa, Wa.P, Pb);
  B = 0.5*log((double)((1. + Wa.Vx)/(1. - Wa.Vx))) + sign*integral; 
	      
  return tanh(B); 
}
/*----------------------------------------------------------------------------*/
/*! \fn Real getXi(const Prim1DS Wa, Real P, Real vxc, char dir) 
 *  \brief Determine the self-similarity variable xi for rarefactions (eq 3.15) 
*/
Real getXi(const Prim1DS Wa, Real P, Real vxc, char dir) 
{
  Real vt, h, G, A; /* vars associated with state ahead of
		      rarefaction wave */
  Real dc, hc, vtc, cs, v2, num, denom; /* vars inside wave */
  Real sign;  

  if (dir == 'L') sign = -1.0; 
  if (dir == 'R') sign = 1.0; 
  
  vt = sqrt((double)(Wa.Vy*Wa.Vy + Wa.Vz*Wa.Vz)); 
  h = 1.0 + Gamma*Wa.P/((Gamma-1.0)*Wa.d);  
  G = 1.0/sqrt((double)(1.0 - Wa.Vx*Wa.Vx - Wa.Vy*Wa.Vy - Wa.Vz*Wa.Vz)); 
  A = h*G*vt;

  dc = Wa.d*pow((double)(P/Wa.P), (double) (1.0/Gamma));
  hc =  1.0 + Gamma*P/((Gamma-1.0)*dc);  
  vtc = A*sqrt((double)((1.0 - vxc*vxc)/(hc*hc + A*A)));
  cs = sqrt((double)(Gamma*(Gamma-1.0)*P/((Gamma -1.0)*dc + Gamma*P)));

  v2 = vxc*vxc + vtc*vtc; 
  num = vxc*(1-cs*cs) + sign*cs*sqrt((double)((1-v2)*(1.-v2*cs*cs-vxc*vxc*(1.-cs*cs))));
  denom = 1.0 - v2*cs*cs; 
  
  return num/denom; 
}
/*---------------------------------------------------------------------------*/ 
/*! \fn void getVelT_Raref(const Prim1DS Wa, Real P, Real vxb, Real *pVy, 
 *			   Real *pVz) 
 *  \brief Determine the tangential velocities vy and vz behind a rarefaction 
 *  wave */
void getVelT_Raref(const Prim1DS Wa, Real P, Real vxb, Real *pVy, Real *pVz) 
{
  Real va_t, ha, Ga, A; 
  Real db, hb, vb_t, Gb, C;
  Real signY = 1.0; 
  Real signZ = 1.0;

  if (Wa.Vy == 0.0 && Wa.Vz == 0.0) { 
    *pVy = 0.0;  
    *pVz = 0.0;  
    return;
  }
  
  va_t = sqrt((double)(Wa.Vy*Wa.Vy + Wa.Vz*Wa.Vz)); 
  ha = 1.0 + Gamma*Wa.P/((Gamma-1.0)*Wa.d);  
  Ga = 1.0/sqrt((double)(1.0 - Wa.Vx*Wa.Vx - Wa.Vy*Wa.Vy - Wa.Vz*Wa.Vz)); 
  A = ha*Ga*va_t;

  db = Wa.d*pow((double)(P/Wa.P), (double) (1.0/Gamma));
  hb =  1.0 + Gamma*P/((Gamma-1.0)*db);  

  vb_t = A*sqrt((double)((1.0 - vxb*vxb)/(hb*hb + A*A)));
  Gb = 1.0/sqrt((double)(1.0 -vxb*vxb - vb_t*vb_t));
  
 
  if (Wa.Vy == 0.0) { 
    if (Wa.Vz < 0.0) signZ = -1.0; 
    *pVy = 0.0;
    *pVz = signZ*vb_t;
  }
  else if (Wa.Vz == 0.0) {
    if (Wa.Vy < 0.0) signY = -1.0; 
      *pVz = 0.0;
      *pVy = signY*vb_t;
  }
  else {
    C = (Wa.Vy*Wa.Vy) / (Wa.Vz*Wa.Vz);
    if (Wa.Vy < 0.0) signY = -1.0;
    if (Wa.Vz < 0.0) signZ = -1.0; 
    
    *pVz = signZ*sqrt((double)(vb_t*vb_t / (1.0 + C))); 
    *pVy = signY*sqrt((double)(vb_t*vb_t / (1.0 + 1.0/C))); 
  }
}
/*----------------------------------------------------------------------------*/
/*! \fn void getVelT_Shock(const Prim1DS Wa, Real P, Real vxb, Real *pVy, 
 *			   Real *pVz)
 *  \brief Determine the tangential velocitiies vy and vz behind a shock wave. 
 *
 *  Note that this function differs from Pons et al. (2000) equation 4.13, 
 *  which is only valid for one of vz or vy equal to zero. */
void getVelT_Shock(const Prim1DS Wa, Real P, Real vxb, Real *pVy, Real *pVz)
{
  Real va_t, ha, Ga, Ay, Az; 
  Real hb, db, vy, vz, Gb;
  Real Cy, Cz; 
  Real D; 

  if (Wa.Vy == 0.0 && Wa.Vz == 0.0) { 
    *pVy = 0.0;  
    *pVz = 0.0;  
    return;
  }

  ha = 1.0 + Gamma*Wa.P/((Gamma-1.0)*Wa.d);  
  Ga = 1.0/sqrt((double)(1.0 - Wa.Vx*Wa.Vx - Wa.Vy*Wa.Vy - Wa.Vz*Wa.Vz)); 
  
  Ay = ha*Ga*Wa.Vy;
  Az = ha*Ga*Wa.Vz; 
  db = Wa.d*pow((double)(P/Wa.P), (double) (1.0/Gamma));
  hb =  1.0 + Gamma*P/((Gamma-1.0)*db);  

  Cz = Az*Az/(hb*hb+Az*Az);
  Cy = Ay*Ay/(hb*hb+Ay*Ay);
  D = (1.0 - Cz*Cy);

  vz = sqrt((double)(Cz*(1-vxb*vxb)*(1-Cy)/D));
  vy = sqrt((double)(Cy*(1-vxb*vxb)*(1-Cz)/D));
  
  Gb = 1.0/sqrt((double)(1.0 - vxb*vxb - vy*vy - vz*vz));

  if (Wa.Vy >= 0.0) 
    *pVy = vy; 
  else 
    *pVy = -vy;
  
  if (Wa.Vz >= 0.0) 
    *pVz = vz; 
  else
    *pVz = -vz;
}
/*---------------------------------------------------------------------------*/ 
/*! \fn void setFluxes(Real vx, Real vy, Real vz, Real P, Real d, Cons1DS *pF)
 *  \brief Set the corresponding fluxes as the fields of pF */
void setFluxes(Real vx, Real vy, Real vz, Real P, Real d, Cons1DS *pF)
{
  Real G, h, D, Sx, Sy, Sz; 
  
  if ((vx*vx + vy*vy + vz*vz) >= 1.0)
    ath_error("[exact_flux]: Superluminal velocities vx = %f, vy = %f, vz = %f\n",
	      vx, vy, vz);

  G = 1.0/sqrt((double)(1.0 - vx*vx - vy*vy - vz*vz));
  h = 1.0 + Gamma*P/((Gamma-1.0)*d); 
  D = d*G; 
  Sx = d*h*G*G*vx;
  Sy = d*h*G*G*vy;
  Sz = d*h*G*G*vz; 

  pF->d  = D*vx;
  pF->Mx = Sx*vx + P;
  pF->My = Sy*vx;
  pF->Mz = Sz*vx;
  pF->E = Sx;
}
/*----------------------------------------------------------------------------*/
/*! \fn void gauleg(double x1, double x2, Real x[], Real w[], int n)
 *  \brief Given the lower and upper limits of integration x1 and x2, and 
 *   the abscissas and weights of the Gauss-Legendre n-point quadrature formula.
 *  
 *  Given the lower and upper limits of integration x1 and x2, and given n, 
 *   this routine returns arrays x[1..n] and w[1..n] of length n, containing 
 *   the abscissas and weights of the Gauss-Legendre n-point quadrature formula.
 *   See Numerical Recipes in C, Press et al. */
void gauleg(double x1, double x2, Real x[], Real w[], int n)
{
  int m,j,i;
  Real z1,z,xm,xl,pp,p3,p2,p1;
  m=(n+1)/2; /* The roots are symmetric, so we only find half of them. */
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  
  for (i=1; i<=m; i++) { /* Loop over the desired roots. */
    z=cos(3.141592654*(i-0.25)/(n+0.5));
    /* Starting with the above approximation to the ith root, we enter */
    /* the main loop of refinement by Newton's method.                 */
    do {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) { /* Recurrence to get Legendre polynomial. */
	p3=p2;
	p2=p1;
	p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      /* p1 is now the desired Legendre polynomial. We next compute */
      /* pp, its derivative, by a standard relation involving also  */
      /* p2, the polynomial of one lower order.                     */
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp; /* Newton's method. */
    } while (fabs(z-z1) > EPS);
    x[i]=xm-xl*z;      /* Scale the root to the desired interval, */
    x[n+1-i]=xm+xl*z;  /* and put in its symmetric counterpart.   */
    w[i]=2.0*xl/((1.0-z*z)*pp*pp); /* Compute the weight             */
    w[n+1-i]=w[i];                 /* and its symmetric counterpart. */
  }
}

/*-------------------------------------------------------------------------*/

/*! \fn Real rtbis_vel(Real (*func)(const Prim1DS Wl, const Prim1DS Wr, Real p),
 *	       const Prim1DS Wl, const Prim1DS Wr, Real x1, Real x2, Real xacc)
 *  \brief Implements bisection root finding algorithm to solve the velocity 
 *  eqn. 
 *
 * Assumes func(x1) and func(x2) have opposite signs without a check */
Real rtbis_vel(Real (*func)(const Prim1DS Wl, const Prim1DS Wr, Real p),
	       const Prim1DS Wl, const Prim1DS Wr, Real x1, Real x2, Real xacc)
{
  Real j; 
  Real dx, f, fmid, xmid, rtb; 
  Real vRel_0 = (Wl.Vx-Wr.Vx)/(1.0-Wl.Vx*Wr.Vx);

  if (wave == RS || wave == SR) {
    f = getVlim_RS(Wl, Wr) - vRel_0;
    fmid = getVlim_2S(Wl, Wr) - vRel_0; 
  }
 
  if (wave == Two_R) {
    f = getVlim_2R(Wl, Wr) - vRel_0; 
    fmid = getVlim_RS(Wl, Wr) - vRel_0; 
  }
  
  if (wave == Two_S) {
    f = getVlim_2S(Wl, Wr) - vRel_0; 
    fmid = (*func)(Wl, Wr, x2); 
  }
  if (f < 0.0) {
    dx = x2-x1;
    rtb =x1; 
  }
  else {
    dx = x1-x2; 
    rtb = x2; 
  }
  
  for (j = 1; j <= JMAX; j++) {
    dx = dx / 2.0; 
    xmid = rtb + dx;
    fmid = (*func)(Wl, Wr, xmid);  
    if (fmid <= 0.0) rtb = xmid;
    if (fabs(dx) < xacc || fmid == 0.0) {
      return rtb;
    }
  }   
}
/*----------------------------------------------------------------------------*/
/*! \fn Real rtbis_xi(Real (*func)(const Prim1DS Wa, Real P, Real vx, char), 
 *		      const Prim1DS Wa, char dir, Real f_x1, Real f_x2, Real x1,
 *		      Real x2, Real xacc)
 *  \brief Implements bisection root finding algorithm to find the rarefaction 
 *  head and tail velocity. 
 *  Assumes func(x1) and func(x2) have opposite signs without a check */ 
   
Real rtbis_xi(Real (*func)(const Prim1DS Wa, Real P, Real vx, char), 
	      const Prim1DS Wa, char dir, Real f_x1, Real f_x2, Real x1, 
	      Real x2, Real xacc) 
{
  Real j, vx; 
  Real dx, f, fmid, xmid, rtb; 

  f = f_x1; 
  fmid = f_x2; 

  if (f < 0.0) {
    dx = x2-x1;
    rtb =x1; 
  }
  else {
    dx = x1-x2; 
    rtb = x2; 
  }
  
  for (j = 1; j <= JMAX; j++) {
    dx = dx / 2.0; 
    xmid = rtb + dx;
    vx = getVb_Raref(Wa, xmid, dir); 
    fmid = (*func)(Wa, xmid, vx, dir);
 
    
    if (fmid <= 0.0) rtb = xmid;
    if (fabs(dx) < xacc || fmid == 0.0) {
      return rtb;
    }
  }
  
 
}
#endif /* Special Relativity */
#endif /* Exact flux */
