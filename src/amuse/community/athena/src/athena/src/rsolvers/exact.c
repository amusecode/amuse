#include "../copyright.h"
/*============================================================================*/
/*! \file exact.c
 *  \brief Computes 1D fluxes using exact nonlinear Riemann solver.
 *
 * PURPOSE: Computes 1D fluxes using exact nonlinear Riemann solver.
 *   Currently only isothermal hydrodynamics has been implemented.  
 *
 * REFERENCES:
 * - R.J. LeVeque, "Numerical Methods for Conservation Laws", 2nd ed.,
 *   Birkhauser Verlag, Basel, (1992).
 *
 * - E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics",
 *   2nd ed., Springer-Verlag, Berlin, (1999).
 *
 * HISTORY:
 * - dec-2006  Isothermal hydro version written by Nicole Lemaster
 * - mar-2010  Adiabatic hydro version written by Nick Hand
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - fluxes() - all Riemann solvers in Athena must have this function name and
 *              use the same argument list as defined in rsolvers/prototypes.h*/
/*============================================================================*/

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
#ifndef SPECIAL_RELATIVITY

#ifdef MHD
#error : The exact flux for MHD has not been implemented.
#endif /* MHD */

#if (NSCALARS > 0)
#error : Passive scalars have not been implemented in the exact flux.
#endif /* NSCALARS */

#ifdef ISOTHERMAL

static void srder(double dm, double vl, double vr, double dmin, double dmax, 
		  double *y, double *dydx);
static double rtsafe(void (*funcd)(double, double, double, double, double,
				   double *, double *), double x1, double x2, double xacc, 
		     double vl, double vr, double dmin, double dmax);

/*----------------------------------------------------------------------------*/
/*! \fn void fluxes(const Cons1DS Ul, const Cons1DS Ur,
 *            const Prim1DS Wl, const Prim1DS Wr,
 *            const Real Bxi, Cons1DS *pF)
 *  \brief Calculates fluxes of CONSERVED variables
 *
 *   Input Arguments:
 *   - Bxi = B in direction of slice at cell interface
 *   - Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *   Output Arguments:
 *   - pFlux = pointer to fluxes of CONSERVED variables at cell interface
 */


void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr,
            const Real Bxi, Cons1DS *pF)
{
  Real zl, zr, zm, dm, Vxm, Mxm, tmp, dmin, dmax;
  Real sl, sr;    /* Left and right going shock velocity */
  Real hdl, hdr;  /* Left and right going rarefaction head velocity */
  Real tll, tlr;  /* Left and right going rarefaction tail velocity */
  char soln;      /* two bits: 0=shock, 1=raref */

  if(!(Ul.d > 0.0)||!(Ur.d > 0.0))
    ath_error("[exact flux]: Non-positive densities: dl = %e  dr = %e\n",
	      Ul.d, Ur.d);

/*--- Step 1. ------------------------------------------------------------------
 * Compute the density and momentum of the intermediate state
 */

  zl = sqrt((double)Wl.d);
  zr = sqrt((double)Wr.d);

  /* --- 1-shock and 2-shock --- */
  soln = 0;

  /* Start by finding density if shocks on both left and right.
   * This will only be the case if dm > Wl.d and dm > Wr.d */
  tmp = zl*zr*(Wl.Vx - Wr.Vx)/(2.0*Iso_csound*(zl + zr));
  zm = tmp + sqrt((double)(tmp*tmp + zl*zr));
  dm = zm*zm;

  /* Get velocity from 1-shock formula */
  Vxm = Wl.Vx - Iso_csound*(dm-Wl.d)/(zm*zl);

  /* If left or right density is greater than intermediate density,
   * then at least one side has rarefaction instead of shock */
  dmin = MIN(Wl.d, Wr.d);
  dmax = MAX(Wl.d, Wr.d);
  if (dm < dmax) {
    /* --- 1-rarefaction and 2-rarefaction --- */
    soln = 3;

    /* Try rarefactions on both left and right, since it's a quicker
     * calculation than 1-shock+2-raref or 1-raref+2-shock */
    dm = zl*zr*exp((Wl.Vx-Wr.Vx)/(2.0*Iso_csound));

    /* Get velocity from 1-rarefaction formula */
    Vxm = Wl.Vx - Iso_csound*log(dm/Wl.d);

    /* If left or right density is smaller than intermediate density,
     * then we must instead have a combination of shock and rarefaction */
    if (dm > dmin) {
      /* --- EITHER 1-rarefaction and 2-shock ---
       * --- OR     1-shock and 2-rarefaction --- */

      /* Solve iteratively equation for shock and rarefaction
       * If Wl.d > Wr.d ==> 1-rarefaction and 2-shock
       * If Wr.d > Wl.d ==> 1-shock and 2-rarefaction */
      if (Wl.d > Wr.d) soln = 2; else soln = 1;

      dm = rtsafe(&srder,dmin,dmax, 2.0*DBL_EPSILON, Wl.Vx, Wr.Vx, dmin, dmax);

      /* Don't be foolish enough to take ln of zero */
      if ((dm > dmin) && (dm <= dmax)) {
        if (Wl.d > Wr.d) {
          /* Get velocity from 1-rarefaction formula */
          Vxm = Wl.Vx - Iso_csound*log(dm/Wl.d);
        } else {
          /* Get velocity from 2-rarefaction formula */
          Vxm = Wr.Vx + Iso_csound*log(dm/Wr.d);
        }
      } else {
        /* --- DEFAULT 1-rarefaction and 2-rarefaction --- */
        soln = 3;

        /* In the event that the intermediate density fails to fall between
         * the left and right densities (should only happen when left and
         * right densities differ only slightly and intermediate density
         * calculated in any step has significant truncation and/or roundoff
         * errors), default to rarefactions on both left and right */
        dm = zl*zr*exp((Wl.Vx-Wr.Vx)/(2.0*Iso_csound));

        /* Get velocity from 1-rarefaction formula */
        Vxm = Wl.Vx - Iso_csound*log(dm/Wl.d);
      }
    }
  }

  if (dm < 0.0)
    ath_error("[exact flux]: Solver finds negative density %5.4e\n", dm);

/*--- Step 2. ------------------------------------------------------------------
 * Calculate the Interface Flux if the wave speeds are such that we aren't
 * actually in the intermediate state
 */

  if (soln & 2) { /* left rarefaction */
    /* The L-going rarefaction head/tail velocity */
    hdl = Wl.Vx - Iso_csound;
    tll = Vxm - Iso_csound;

    if (hdl >= 0.0) {
      /* To left of rarefaction */
      pF->d  = Ul.Mx;
      pF->Mx = Ul.Mx*(Wl.Vx) + Wl.d*Iso_csound2;
      pF->My = Ul.My*(Wl.Vx);
      pF->Mz = Ul.Mz*(Wl.Vx);
      return;
    } else if (tll >= 0.0) {
      /* Inside rarefaction fan */
      dm = Ul.d*exp(hdl/Iso_csound);
      Mxm = Ul.d*Iso_csound*exp(hdl/Iso_csound);
      Vxm = (dm == 0.0 ? 0.0 : Mxm / dm);

      pF->d  = Mxm;
      pF->Mx = Mxm*Vxm + dm*Iso_csound2;
      pF->My = Mxm*Wl.Vy;
      pF->Mz = Mxm*Wl.Vz;
      return;
    }
  } else { /* left shock */
    /* The L-going shock velocity */
    sl = Wl.Vx - Iso_csound*sqrt(dm)/zl;

    if(sl >= 0.0) {
      /* To left of shock */
      pF->d  = Ul.Mx;
      pF->Mx = Ul.Mx*(Wl.Vx) + Wl.d*Iso_csound2;
      pF->My = Ul.My*(Wl.Vx);
      pF->Mz = Ul.Mz*(Wl.Vx);
      return;
    }
  }

  if (soln & 1) { /* right rarefaction */
    /* The R-going rarefaction head/tail velocity */
    hdr = Wr.Vx + Iso_csound;
    tlr = Vxm + Iso_csound;

    if (hdr <= 0.0) {
      /* To right of rarefaction */
      pF->d  = Ur.Mx;
      pF->Mx = Ur.Mx*(Wr.Vx) + Wr.d*Iso_csound2;
      pF->My = Ur.My*(Wr.Vx);
      pF->Mz = Ur.Mz*(Wr.Vx);
      return;
    } else if (tlr <= 0.0) {
      /* Inside rarefaction fan */
      tmp = dm;
      dm = tmp*exp(-tlr/Iso_csound);
      Mxm = -tmp*Iso_csound*exp(-tlr/Iso_csound);
      Vxm = (dm == 0.0 ? 0.0 : Mxm / dm);

      pF->d  = Mxm;
      pF->Mx = Mxm*Vxm + dm*Iso_csound2;
      pF->My = Mxm*Wr.Vy;
      pF->Mz = Mxm*Wr.Vz;
      return;
    }
  } else { /* right shock */
    /* The R-going shock velocity */
    sr = Wr.Vx + Iso_csound*sqrt(dm)/zr;

    if(sr <= 0.0) {
      /* To right of shock */
      pF->d  = Ur.Mx;
      pF->Mx = Ur.Mx*(Wr.Vx) + Wr.d*Iso_csound2;
      pF->My = Ur.My*(Wr.Vx);
      pF->Mz = Ur.Mz*(Wr.Vx);
      return;
    }
  }

/* If we make it this far, then we're in the intermediate state */

/*--- Step 3. ------------------------------------------------------------------
 * Calculate the Interface Flux */

  if(Vxm >= 0.0){
    pF->d  = dm*Vxm;
    pF->Mx = dm*Vxm*Vxm + dm*Iso_csound2;
    pF->My = dm*Vxm*Wl.Vy;
    pF->Mz = dm*Vxm*Wl.Vz;
  }
  else{
    pF->d  = dm*Vxm;
    pF->Mx = dm*Vxm*Vxm + dm*Iso_csound2;
    pF->My = dm*Vxm*Wr.Vy;
    pF->Mz = dm*Vxm*Wr.Vz;
  }

  return;
}

/*! \fn static void srder(double dm, double vl, double vr, double dmin, 
 *			  double dmax, double *y, double *dydx)
 *  \brief Equation to solve iteratively for shock and rarefaction as well
 * as its derivative.  Used by rtsafe() */

static void srder(double dm, double vl, double vr, double dmin, double dmax, 
			double *y, double *dydx)
{
  *y = (vr - vl) + Iso_csound*(log(dm/dmax) + (dm-dmin)/sqrt(dm*dmin));
  *dydx = Iso_csound/dm*(1.0 + 0.5*(dm+dmin)/sqrt(dm*dmin));

  return;
}

/*! \fn static double rtsafe(void (*funcd)(double, double, double, double, 
 *			    double, double *, double *), double x1, double x2, 
 *			    double xacc, double vl, double vr, double dmin, 
 *			    double dmax)
 *  \brief Numerical Recipes function rtsafe modified to use doubles and take
 * extra parameters to pass to the derivative function above */

#define MAXIT 100

static double rtsafe(void (*funcd)(double, double, double, double, double,
			double *, double *), double x1, double x2, double xacc, 
			double vl, double vr, double dmin, double dmax)
{
	int j;
	double df,dx,dxold,f,fh,fl;
	double temp,xh,xl,rts;

	(*funcd)(x1,vl,vr,dmin,dmax,&fl,&df);
	(*funcd)(x2,vl,vr,dmin,dmax,&fh,&df);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
		return 0.0;
	}
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	rts=0.5*(x1+x2);
	dxold=fabs(x2-x1);
	dx=dxold;
	(*funcd)(rts,vl,vr,dmin,dmax,&f,&df);
	for (j=1;j<=MAXIT;j++) {
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
				if (xl == rts) return rts;
		} else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts;
		(*funcd)(rts,vl,vr,dmin,dmax,&f,&df);
		if (f < 0.0)
			xl=rts;
		else
			xh=rts;
	}
        /* Cap the number of iterations but don't fail */
	return rts;
}

#undef MAXIT 

#elif defined ADIABATIC

/*------------------------------------------------------------------*/
/*! \fn static Real guessP(const Prim1DS Wl, const Prim1DS Wr)
 *  \brief  Provides a guess value for the pressure in the center region. 
 *
 * The choice is made according to the Two-Shock approximate
 * Riemman solver algorithm. Returns pressure guess value p_0, or 
 *   TOL if p_0 is less than zero */
 
static Real guessP(const Prim1DS Wl, const Prim1DS Wr)
{
 
  Real al = sqrt((double) Gamma*Wl.P/Wl.d);   /* left sound speed */
  Real ar = sqrt((double) Gamma*Wr.P/Wr.d);   /* right sound speed */
  Real ppv;  /* inital guess for pressure */
  Real gl, gr, p_0;
  Real TOL = 1.0e-6; 

  ppv = 0.5*(Wl.P + Wr.P)-0.125*(Wr.Vx-Wl.Vx)*(Wl.d+Wr.d)*(al+ar);

  if (ppv < 0.0)
    ppv = 0.0;

  gl = sqrt((double) (2.0/(Wl.d*(Gamma+1)))/
	    ((Gamma-1)*Wl.P/(Gamma+1) + ppv)); 
  gr = sqrt((double) (2.0/(Wr.d*(Gamma+1)))/
	    ((Gamma-1)*Wr.P/(Gamma+1) + ppv));
  p_0 = (gl*Wl.P + gr*Wr.P - (Wr.Vx-Wl.Vx))/(gr + gl); 

  if (p_0 < 0.0)
    p_0 = TOL; 

  return p_0; 
}  
/*------------------------------------------------------------------*/

/*! \fn static Real PFunc(const Prim1DS W, Real POld)
 *  \brief Evaluates the pressure function (see Toro) for the exact 
   riemann solver, given the pressure POld */
static Real PFunc(const Prim1DS W, Real POld)
{
  Real f; 
  Real a = sqrt((double) Gamma*W.P/W.d);
  Real Ak, Bk; 

  /* rarefaction wave */
  if (POld <= W.P)  
    f = 2*a/(Gamma-1) * (pow((double) POld/W.P, (double) 
			     (Gamma-1)/(2*Gamma)) - 1); 

  /* shock wave */
  else {
    Ak = 2/(W.d*(Gamma+1)); 
    Bk = W.P*(Gamma-1)/(Gamma+1); 

    f = (POld - W.P) * sqrt((double) Ak/(POld + Bk)); 
  }

  return f; 
}
   
/*------------------------------------------------------------------*/

/*! \fn static Real PFuncDeriv(const Prim1DS W, Real POld)
 *  \brief  evaluate the derivative of the pressure function (see Toro) 
 * at the pressure value POld */
static Real PFuncDeriv(const Prim1DS W, Real POld)
{
  Real fDer; 
  Real a = sqrt((double) Gamma*W.P/W.d);
  Real Ak, Bk; 

  /* rarefaction wave */
  if (POld <= W.P)  
    fDer = 1/(a * W.d) * pow((double) POld/W.P, (double) -
			     (Gamma+1)/(2*Gamma)); 

  /* shock wave */
  else {
    Ak = 2/(W.d*(Gamma+1)); 
    Bk = W.P*(Gamma-1)/(Gamma+1); 

    fDer = sqrt((double) Ak/(POld + Bk))
      *(1.0 - 0.5*(POld - W.P)/(Bk + POld)); 
  }

  return fDer; 
}
   
/*------------------------------------------------------------------*/ 
    
/*! \fn static Real getPC(const Prim1DS Wl, const Prim1DS Wr)
 *  \brief Uses Newton-Raphson iteration to find the alebraic root of the 
 * pressure function and determine the pressure solution in the 
 * center region, Pc. Fails if iteration diverges */
static Real getPC(const Prim1DS Wl, const Prim1DS Wr)
{
  Real POld; 
  Real VxDiff = Wr.Vx - Wl.Vx;
  int i = 0;   
  Real TOL = 1.0e-6;
  Real change, p, fr, fl, frder, flder; 
  Real MAX_ITER = 100; 

  POld = guessP(Wl, Wr);
  
  while (i < MAX_ITER) {

    fl = PFunc(Wl, POld); 
    fr = PFunc(Wr, POld); 
    flder = PFuncDeriv(Wl, POld); 
    frder = PFuncDeriv(Wr, POld); 
    
    p = POld - (fl + fr + VxDiff)/(flder + frder);

    change = 2.0*fabs((p-POld)/(p+POld)); 
    if (change <= TOL)
      return p;
    
    if (p < 0.0)
      p = TOL; 

    POld = p; 
    i++; 
  }
   
  /* exit failure due to divergence */
  ath_error("[exact_flux]: Divergence in Newton-Raphson \
iteration: p = %e\n", p); 
  
}
/*------------------------------------------------------------------*/
/*! \fn void fluxes(const Cons1DS Ul, const Cons1DS Ur, 
 *	    const Prim1DS Wl, const Prim1DS Wr, 
 *	    const Real Bxi, Cons1DS *pF)
 *  \brief Computes fluxes of CONSERVES variables
 *
 *   Input Arguments:
 *   - Bxi = B in direction of slice at cell interface
 *   - Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *   Output Arguments:
 *   - pF = pointer to fluxes of CONSERVED variables at cell interface
 */

void fluxes(const Cons1DS Ul, const Cons1DS Ur, 
	    const Prim1DS Wl, const Prim1DS Wr, 
	    const Real Bxi, Cons1DS *pF)
{
  Real pc, fr, fl; 
  Real dc, dcl, dcr, Vxc;
  Real sl, sr;    /* Left and right going shock velocity */
  Real hdl, hdr;  /* Left and right going rarefaction head velocity */
  Real tll, tlr;  /* Left and right going rarefaction tail velocity */
  Real al = sqrt((double) Gamma*Wl.P/Wl.d); /* left sound speed */
  Real ar = sqrt((double) Gamma*Wr.P/Wr.d); /* right sound speed */
  Real tmp1, tmp2;
  Real e, V, E; 
 
  if(!(Ul.d > 0.0)||!(Ur.d > 0.0))
    ath_error("[exact flux]: Non-positive densities: dl = %e  dr = %e\n", 
	      Ul.d, Ur.d);

  pc = getPC(Wl, Wr); 
  
  /*----------------------------------------------------------------*/
  /* calculate Vxc */
  fr = PFunc(Wr, pc); 
  fl = PFunc(Wl, pc); 
  Vxc = 0.5*(Wl.Vx + Wr.Vx) + 0.5*(fr - fl); 

  /*-----------------------------------------------------------------*/
  /* calucate density to left of contact (dcl) */
  if (pc > Wl.P) {

    /* left shock wave */
    Real tmp = (Gamma - 1.0) / (Gamma + 1.0);
    dcl = Wl.d*(pc/Wl.P + tmp)/(tmp*pc/Wl.P + 1); 
  }
  else {

    /* left rarefaction wave */
    dcl = Wl.d*pow((double) (pc/Wl.P), (double) (1/Gamma)); 
  }

  if (dcl < 0.0)
     ath_error("[exact flux]: Solver finds negative density %5.4e\n", dcl);
  
  /*-----------------------------------------------------------------*/
  /* calculate density to the right of contact (dcr) */
  if (pc > Wr.P) {

    /* right shock wave */
    Real tmp = (Gamma - 1)/(Gamma + 1);
    dcr = Wr.d*(pc/Wr.P + tmp)/(tmp*pc/Wr.P + 1); 
  }
  else {

    /* right rarefaction wave */
    dcr = Wr.d*pow((double) (pc/Wr.P), (double) (1/Gamma)); 
  }

  if (dcr < 0.0)
    ath_error("[exact flux]: Solver finds negative density %5.4e\n", dcr);
 /*-----------------------------------------------------------------
  * Calculate the Interface Flux if the wave speeds are such that we aren't
  * actually in the intermediate state */
  
  if (pc > Wl.P) {
    /* left shock wave */
    
    /* left shock speed */
    sl = Wl.Vx - al*sqrt((double)(pc*(Gamma+1)/(2*Gamma*Wl.P) + 
				  (Gamma-1)/(2*Gamma))); 
    if (sl >= 0.0) {
      /* to left of shock */
     
      e = Wl.P/(Wl.d*(Gamma-1));
      V = Wl.Vx*Wl.Vx + Wl.Vy*Wl.Vy + Wl.Vz*Wl.Vz;
      E = Wl.d*(0.5*V + e); 
      
      pF->E = Wl.Vx*(E + Wl.P);
      pF->d  = Ul.Mx;
      pF->Mx = Ul.Mx*(Wl.Vx) + Wl.P;
      pF->My = Ul.My*(Wl.Vx);
      pF->Mz = Ul.Mz*(Wl.Vx);

      return; 
    }
  }
  else {
    /* left rarefaction */
    
    Real alc = al*pow((double)(pc/Wl.P), (double)(Gamma-1)/(2*Gamma));
    
    hdl = Wl.Vx - al; 
    tll = Vxc - alc; 
    
    if (hdl >= 0.0) {
      /* To left of rarefaction */
	
      e = Wl.P/(Wl.d*(Gamma-1));
      V = Wl.Vx*Wl.Vx + Wl.Vy*Wl.Vy + Wl.Vz*Wl.Vz;
      E = Wl.d*(0.5*V + e); 
      
      pF->E =  Wl.Vx*(E + Wl.P);
      pF->d  = Ul.Mx;
      pF->Mx = Ul.Mx*(Wl.Vx) + Wl.P;
      pF->My = Ul.My*(Wl.Vx);
      pF->Mz = Ul.Mz*(Wl.Vx);
      return;
      
    } 
    else if (tll >= 0.0) {
      /* Inside rarefaction fan */
       
      
      tmp1 = 2/(Gamma + 1); 
      tmp2 = (Gamma - 1)/(al*(Gamma+1)); 
      
      dc = Wl.d*pow((double)(tmp1 + tmp2*Wl.Vx), (double)(2/(Gamma-1))); 
      Vxc = tmp1*(al + Wl.Vx*(Gamma-1)/2);
      pc = Wl.P*pow((double)(tmp1+tmp2*Wl.Vx),(double)(2*Gamma/(Gamma-1))); 
      
      e = pc/(dc*(Gamma-1));
      V = Vxc*Vxc + Wl.Vy*Wl.Vy + Wl.Vz*Wl.Vz;
      E = dc*(0.5*V + e); 	
      
      pF->E = Vxc*(E + pc); 
      pF->d  = dc*Vxc;
      pF->Mx = dc*Vxc*Vxc + pc;
      pF->My = dc*Vxc*Wl.Vy;
      pF->Mz = dc*Vxc*Wl.Vz;
      return;
    } 
  }

  if (pc > Wr.P) {
    /* right shock wave */
    
    /* right shock speed */
    sr = Wr.Vx + ar*sqrt((double)(pc*(Gamma+1)/(2*Gamma*Wr.P) + 
				  (Gamma-1)/(2*Gamma))); 
    if (sr <= 0.0) {
      /* to right of shock */

      e = Wr.P/(Wr.d*(Gamma-1));
      V = Wr.Vx*Wr.Vx + Wr.Vy*Wr.Vy + Wr.Vz*Wr.Vz;
      E = Wr.d*(0.5*V + e);     
      
      pF->E = Wr.Vx*(E + Wr.P);
      pF->d  = Ur.Mx;
      pF->Mx = Ur.Mx*(Wr.Vx) + Wr.P;
      pF->My = Ur.My*(Wr.Vx);
      pF->Mz = Ur.Mz*(Wr.Vx);
      return; 
    }
  }
  else {
    /* right rarefaction */
    
    Real arc = ar*pow((double)(pc/Wr.P), (double)(Gamma-1)/(2*Gamma)); 
      
    hdr = Wr.Vx + ar; 
    tlr = Vxc + arc; 
    
    if (hdr <= 0.0) {
      /* To right of rarefaction */

      e = Wr.P/(Wr.d*(Gamma-1));
      V = Wr.Vx*Wr.Vx + Wr.Vy*Wr.Vy + Wr.Vz*Wr.Vz;
      E = Wr.d*(0.5*V + e); 	
      
      pF->E = Wr.Vx*(E + Wr.P);
      pF->d  = Ur.Mx;
      pF->Mx = Ur.Mx*(Wr.Vx) + Wr.P;
      pF->My = Ur.My*(Wr.Vx);
      pF->Mz = Ur.Mz*(Wr.Vx);
      return;
      
    } else if (tlr <= 0.0) {
      /* Inside rarefaction fan */
      
      tmp1 = 2/(Gamma + 1); 
      tmp2 = (Gamma - 1)/(ar*(Gamma+1)); 
      
      dc = Wr.d*pow((double)(tmp1 - tmp2*Wr.Vx), (double)(2/(Gamma-1))); 
      Vxc = tmp1*(-ar + Wr.Vx*(Gamma-1)/2);
      pc = Wr.P*pow((double)(tmp1-tmp2*Wr.Vx), (double)(2*Gamma/(Gamma-1))); 
      
      e = pc/(dc*(Gamma-1));
      V = Vxc*Vxc + Wr.Vy*Wr.Vy + Wr.Vz*Wr.Vz;
      E = dc*(0.5*V + e); 
      
      pF->E = Vxc*(E + pc);
      pF->d  = dc*Vxc;
      pF->Mx = dc*Vxc*Vxc + pc;
      pF->My = dc*Vxc*Wr.Vy;
      pF->Mz = dc*Vxc*Wr.Vz;
      return;
    }
  }
    
/* We are in the intermediate state */

/*---------------------------------------------------------------------
 * Calculate the Interface Flux */
  if (Vxc >= 0.0) {

    e = pc/(dcl*(Gamma-1));
    V = Vxc*Vxc + Wl.Vy*Wl.Vy + Wl.Vz*Wl.Vz;
    E = dcl*(0.5*V + e); 
    
    pF->E = Vxc*(E + pc); 
    pF->d  = dcl*Vxc;
    pF->Mx = dcl*Vxc*Vxc + pc;
    pF->My = dcl*Vxc*Wl.Vy;
    pF->Mz = dcl*Vxc*Wl.Vz;
  }
  else {

    e = pc/(dcr*(Gamma-1));
    V = Vxc*Vxc + Wr.Vy*Wr.Vy + Wr.Vz*Wr.Vz;
    E = dcr*(0.5*V + e); 
    
    pF->E = Vxc*(E + pc); 
    pF->d  = dcr*Vxc;
    pF->Mx = dcr*Vxc*Vxc + pc;
    pF->My = dcr*Vxc*Wr.Vy;
    pF->Mz = dcr*Vxc*Wr.Vz;
  }
  
  return;
}
/*------------------------------------------------------------------*/

#endif /*ISOTHERMAL*/

#endif /* SPECIAL_RELATIVITY */
#endif /* EXACT_FLUX */
