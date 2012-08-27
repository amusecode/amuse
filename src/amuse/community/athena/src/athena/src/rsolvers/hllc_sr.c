#include "../copyright.h"
/*============================================================================*/
/*! \file hllc_sr.c
 *  \brief Computes 1D fluxes using the relativistic HLLC Riemann solver. 
 *
 * PURPOSE: Computes 1D fluxes using the relativistic HLLC Riemann solver, 
 *   an extension of the HLLE fluxes to include the contact wave.  Currently 
 *   only works for hydrodynamics.  For an extension to MHD, see hlld_sr.c
 *
 * REFERENCES:
 * - A. Mignone and G. Bodo, "An HLLC Riemann solver for relativistic flows",
 *   Mon. Not. R. Astron. Soc. 364, 126-136 (2005)
 *
 * - A. Mignone and G. Bodo, "An HLLC Solver for Relativistic Flows - II:
 *   Magnetohydrodynamics", arxiv:astro-ph/0601640v1 (2006)
 *
 * HISTORY: Written by Jonathan FUlton, February 2009
 *          Extended to MHD by Kris Beckwith, Spring 2010
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
#ifdef SPECIAL_RELATIVITY

#if (NSCALARS > 0)
#error : The SR HLLC flux does not work with passive scalars.
#endif

#ifdef MHD
void flux_LR(Cons1DS U, Prim1DS W, Cons1DS *flux, Real Bx, Real* p);
void getMaxSignalSpeeds_pluto(const Prim1DS Wl, const Prim1DS Wr,
                        const Real Bx, Real* low, Real* high);
void getMaxSignalSpeeds_echo(const Prim1DS Wl, const Prim1DS Wr,
			 const Real Bx, Real* low, Real* high);
void getVChar_pluto(const Prim1DS W, const Real Bx, Real* lml, Real* lmr);
void getVChar_echo (const Prim1DS W, const Real Bx, Real* lml, Real* lmr);
/* solves quartic equation defined by a and returns roots in root
 * returns the number of Real roots 
 * error specifies an accuracy
 * currently force four Real solutions b/c it's physical */
int QUARTIC (Real b, Real c, Real d, Real e, Real z[]);

/* solves cubic equation defined by a and stores roots in root
 * returns number of Real roots */
int CUBIC(Real b, Real c, Real d, Real z[]);
#endif

#ifdef HYDRO
/*----------------------------------------------------------------------------*/
/*! \fn void fluxes(const Cons1DS Ul, const Cons1DS Ur,
 *      const Prim1DS Wl, const Prim1DS Wr, const Real Bxi, Cons1DS *pFlux)
 *  \brief Computes 1D fluxes
 *   Input Arguments:
 *   - Ul,Ur = L/R-states of CONSERVED variables at cell interface 
 *   - Wl,Wr = L/R-states of PRIMITIVE variables at cell interface 
 *   Output Arguments:
 *   - pFlux = pointer to fluxes of CONSERVED variables at cell interface 
 */

void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr, const Real Bxi, Cons1DS *pFlux)
{
  Cons1DS Fl,Fr,Fhll,Uhll,Usl,Usr;
  Real rhl, rhr, csl, csr, cslsq, csrsq, vsql, vsqr, gammasql, gammasqr;
  Real ssl, ssr, radl, radr, lmdapl, lmdapr, lmdaml, lmdamr, lmdatlmda;
  Real lmdal,lmdar; /* Left and Right wave speeds */
  Real lmdas; /* Contact wave speed */
  Real ovlrmll;
  Real a,b,c,quad,rad;
  Real den,ps; /* PressUre in inner region */

/*--- Step 1. ------------------------------------------------------------------
 * Compute the max and min wave speeds used in Mignone 
 */

  rhl = Wl.d + Wl.P * Gamma / Gamma_1; /* Mignone Eq 3.5 */
  rhr = Wr.d + Wr.P * Gamma / Gamma_1;

  csl = sqrt(Gamma * Wl.P / rhl); /* Mignone Eq 4 */
  csr = sqrt(Gamma * Wr.P / rhr);

  cslsq = SQR(csl);
  csrsq = SQR(csr);

  vsql = SQR(Wl.Vx) + SQR(Wl.Vy) + SQR(Wl.Vz);
  vsqr = SQR(Wr.Vx) + SQR(Wr.Vy) + SQR(Wr.Vz);

  gammasql = 1.0 / (1.0 - vsql);
  gammasqr = 1.0 / (1.0 - vsqr);

  ssl = cslsq / ( gammasql * (1.0 - cslsq) ); /* Mignone Eq 22.5 */
  ssr = csrsq / ( gammasqr * (1.0 - csrsq) );

  radl = sqrt( ssl*(1.0-SQR(Wl.Vx)+ssl) ); /* Mignone Eq 23 (radical part) */
  radr = sqrt( ssr*(1.0-SQR(Wr.Vx)+ssr) );

  lmdapl = (Wl.Vx + radl) / (1.0 + ssl); /* Mignone Eq 23 */
  lmdapr = (Wr.Vx + radr) / (1.0 + ssr);
  lmdaml = (Wl.Vx - radl) / (1.0 + ssl);
  lmdamr = (Wr.Vx - radr) / (1.0 + ssr);

  lmdal = MIN(lmdaml, lmdamr); /* Mignone Eq 21 */
  lmdar = MAX(lmdapl, lmdapr);
  

/*--- Step 2. ------------------------------------------------------------------
 * Compute L/R fluxes according to Mignone 2
 */

  Fl.d  = Ul.d * Wl.Vx;
  Fl.Mx = Ul.Mx * Wl.Vx + Wl.P;
  Fl.My = Ul.My * Wl.Vx;
  Fl.Mz = Ul.Mz * Wl.Vx;
  Fl.E  = Ul.Mx;

  Fr.d  = Ur.d * Wr.Vx;
  Fr.Mx = Ur.Mx * Wr.Vx + Wr.P;
  Fr.My = Ur.My * Wr.Vx;
  Fr.Mz = Ur.Mz * Wr.Vx;
  Fr.E  = Ur.Mx;

/*--- Step 3. ------------------------------------------------------------------
 * Compute HLL flux using Mignone Eq 11 (necessary for computing lmdas (Eq 18)
 * Compute HLL conserved quantities using Mignone eq 9
 */

  ovlrmll = 1.0 / ( lmdar - lmdal );
  lmdatlmda = lmdal*lmdar;

  Fhll.d  = (lmdar*Fl.d  - lmdal*Fr.d  + lmdatlmda * (Ur.d  - Ul.d) ) * ovlrmll;
  Fhll.Mx = (lmdar*Fl.Mx - lmdal*Fr.Mx + lmdatlmda * (Ur.Mx - Ul.Mx)) * ovlrmll;
  Fhll.My = (lmdar*Fl.My - lmdal*Fr.My + lmdatlmda * (Ur.My - Ul.My)) * ovlrmll;
  Fhll.Mz = (lmdar*Fl.Mz - lmdal*Fr.Mz + lmdatlmda * (Ur.Mz - Ul.Mz)) * ovlrmll;
  Fhll.E  = (lmdar*Fl.E  - lmdal*Fr.E  + lmdatlmda * (Ur.E  - Ul.E )) * ovlrmll;

  Uhll.d  = (lmdar * Ur.d  - lmdal * Ul.d  + Fl.d  - Fr.d ) * ovlrmll;
  Uhll.Mx = (lmdar * Ur.Mx - lmdal * Ul.Mx + Fl.Mx - Fr.Mx) * ovlrmll;
  Uhll.My = (lmdar * Ur.My - lmdal * Ul.My + Fl.My - Fr.My) * ovlrmll;
  Uhll.Mz = (lmdar * Ur.Mz - lmdal * Ul.Mz + Fl.Mz - Fr.Mz) * ovlrmll;
  Uhll.E  = (lmdar * Ur.E  - lmdal * Ul.E  + Fl.E  - Fr.E ) * ovlrmll;

/*--- Step 4. ------------------------------------------------------------------
 * Compute contact wave speed using larger root from Mignone Eq 18
 * Physical root is the root with the minus sign
 */

  /* quadratic formUla calcUlation */

  a = Fhll.E;
  b = -(Uhll.E + Fhll.Mx);
  c = Uhll.Mx;


  quad = -0.5*(b + SIGN(b)*sqrt(b*b - 4.0*a*c));
  lmdas = c/quad;

/*--- Step 5. ------------------------------------------------------------------
 * Determine intercell flux according to Mignone 13
 */

  if( lmdal >= 0.0){ /* Fl */
    /* intercell flux is left flux */
    pFlux->d  = Fl.d;
    pFlux->Mx = Fl.Mx;
    pFlux->My = Fl.My;
    pFlux->Mz = Fl.Mz;
    pFlux->E  = Fl.E;

   return;
  }
  else if( lmdas >= 0.0){ /* Fls */

    /* Mignone 2006 Eq 48 */
    ps = -Fhll.E*lmdas + Fhll.Mx;

    /* now calcUlate Usl with Mignone Eq 16 */
    den = 1.0 / (lmdal - lmdas);

    Usl.d  = Ul.d * (lmdal - Wl.Vx) * den;
    Usl.Mx = (Ul.Mx * (lmdal - Wl.Vx) + ps - Wl.P) * den;
    Usl.My = Ul.My * (lmdal - Wl.Vx) * den;
    Usl.Mz = Ul.Mz * (lmdal - Wl.Vx) * den;
    Usl.E  = (Ul.E * (lmdal - Wl.Vx) + ps * lmdas - Wl.P * Wl.Vx) * den;

    /* now calcUlate Fsr using Mignone Eq 14 */

    pFlux->d  = lmdal*(Usl.d  - Ul.d ) + Fl.d;
    pFlux->Mx = lmdal*(Usl.Mx - Ul.Mx) + Fl.Mx;
    pFlux->My = lmdal*(Usl.My - Ul.My) + Fl.My;
    pFlux->Mz = lmdal*(Usl.Mz - Ul.Mz) + Fl.Mz;
    pFlux->E  = lmdal*(Usl.E  - Ul.E ) + Fl.E;

    return;
  }
  else if( lmdar >= 0.0){ /* Frs */

    /* Mignone 2006 Eq 48 */
    ps = -Fhll.E*lmdas + Fhll.Mx;

    /* now calcUlate Usr with Mignone Eq 16 */
    den = 1.0 / (lmdar - lmdas);

    Usr.d  = Ur.d * (lmdar - Wr.Vx) * den;
    Usr.Mx = (Ur.Mx * (lmdar - Wr.Vx) + ps - Wr.P) * den;
    Usr.My = Ur.My * (lmdar - Wr.Vx) * den;
    Usr.Mz = Ur.Mz * (lmdar - Wr.Vx) * den;
    Usr.E  = (Ur.E * (lmdar - Wr.Vx) + ps * lmdas - Wr.P * Wr.Vx) * den;

    /* now calcUlate Fsr using Mignone Eq 14 */

    pFlux->d  = lmdar*(Usr.d  - Ur.d ) + Fr.d;
    pFlux->Mx = lmdar*(Usr.Mx - Ur.Mx) + Fr.Mx;
    pFlux->My = lmdar*(Usr.My - Ur.My) + Fr.My;
    pFlux->Mz = lmdar*(Usr.Mz - Ur.Mz) + Fr.Mz;
    pFlux->E  = lmdar*(Usr.E  - Ur.E ) + Fr.E;

    return;
  }
  else{ /* Fr */
    /* intercell flux is right flux */
    pFlux->d  = Fr.d;
    pFlux->Mx = Fr.Mx;
    pFlux->My = Fr.My;
    pFlux->Mz = Fr.Mz;
    pFlux->E  = Fr.E;

    return;
  }
  
  /* need to deal with scalar fluxes */
}

#endif /* HYDRO */

#ifdef MHD
/*----------------------------------------------------------------------------*/
/*! \fn void fluxes(const Cons1DS Ul, const Cons1DS Ur,
 *          const Prim1DS Wl, const Prim1DS Wr, const Real Bxi, Cons1DS *pFlux)
 *  \brief Computes 1D fluxes
 *   Input Arguments:
 *  -  Ul,Ur = L/R-states of CONSERVED variables at cell interface 
 *  -  Wl,Wr = L/R-states of PRIMITIVE variables at cell interface 
 *   Output Arguments:
 *  -  pFlux = pointer to fluxes of CONSERVED variables at cell interface 
 */

void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr, const Real Bxi, Cons1DS *pFlux)
{
  Cons1DS Fl,Fr,Fhll,Uhll,Usl,Usr;
  Prim1DS Wsl,Wsr,Whll;
  Real Sl, Sr, Pl, Pr;
  Real Sla, Sra;
  Real dS_1, scrh;
  Real Bx, Bys, Bzs;
  Real BtFBt, Bt2, FBt2;
  Real a, b, c;
  Real ps, vxs, vys, vzs, gammas_2, vBs, V2l, V2r;
  Real vxl, vxr, alpha_l, alpha_r;
  int switch_to_hll,wave_speed_fail;

  wave_speed_fail = 0;
  switch_to_hll = 0;
	
/*--- Step 1. ------------------------------------------------------------------
 * Compute the max and min wave speeds used in Mignone 
 */
  getMaxSignalSpeeds_pluto(Wl,Wr,Bxi,&Sl,&Sr);
	
  if (Sl != Sl) {
    wave_speed_fail = 1;
    printf("[hllc_sr_mhd]: NaN in Sl %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr =  1.0;
  }
	
  if (Sr != Sr) {
    wave_speed_fail = 1;
    printf("[hllc_sr_mhd]: NaN in Sr %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }
	
  if (Sl < -1.0) {
    wave_speed_fail = 1;
    printf("[hllc_sr_mhd]: Superluminal Sl %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }
  if (Sr > 1.0) {
    wave_speed_fail = 1;
    printf("[hllc_sr_mhd]: Superluminal Sr %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }

/*--- Step 1a. -----------------------------------------------------------------
 * If PLUTO wavespeeds are bad, fall back to the estimate used in ECHO
 */
  if (wave_speed_fail){
    getMaxSignalSpeeds_echo (Wl,Wr,Bxi,&Sla,&Sra);
	
    if (Sla != Sla) {
      switch_to_hll = 1;
      printf("[hllc_sr_mhd]: NaN in Sl %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra =  1.0;
    }
	
    if (Sra != Sra) {
      switch_to_hll = 1;
      printf("[hllc_sr_mhd]: NaN in Sr %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra = 1.0;
    }
	
    if (Sla < -1.0) {
      switch_to_hll = 1;
      printf("[hllc_sr_mhd]: Superluminal Sl %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra = 1.0;
    }
    if (Sra > 1.0) {
      switch_to_hll = 1;
      printf("[hllc_sr_mhd]: Superluminal Sr %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra = 1.0;
    }

    Sl = Sla;
    Sr = Sra;

  }

  /* compute L/R fluxes */
  flux_LR(Ul,Wl,&Fl,Bxi,&Pl);
  flux_LR(Ur,Wr,&Fr,Bxi,&Pr);
	
/*--- Step 2. ------------------------------------------------------------------
 * Construct HLL fluxes & average state
 */
  dS_1 = 1.0/(Sr - Sl);

  Uhll.d  = (Sr*Ur.d  - Sl*Ul.d  + Fl.d  - Fr.d ) * dS_1;
  Uhll.Mx = (Sr*Ur.Mx - Sl*Ul.Mx + Fl.Mx - Fr.Mx) * dS_1;
  Uhll.My = (Sr*Ur.My - Sl*Ul.My + Fl.My - Fr.My) * dS_1;
  Uhll.Mz = (Sr*Ur.Mz - Sl*Ul.Mz + Fl.Mz - Fr.Mz) * dS_1;
  Uhll.E  = (Sr*Ur.E  - Sl*Ul.E  + Fl.E  - Fr.E ) * dS_1;
  Uhll.By = (Sr*Ur.By - Sl*Ul.By + Fl.By - Fr.By) * dS_1;
  Uhll.Bz = (Sr*Ur.Bz - Sl*Ul.Bz + Fl.Bz - Fr.Bz) * dS_1;
		
  Fhll.d  = (Sr*Fl.d  - Sl*Fr.d  + Sl*Sr*(Ur.d  - Ul.d )) * dS_1;
  Fhll.Mx = (Sr*Fl.Mx - Sl*Fr.Mx + Sl*Sr*(Ur.Mx - Ul.Mx)) * dS_1;
  Fhll.My = (Sr*Fl.My - Sl*Fr.My + Sl*Sr*(Ur.My - Ul.My)) * dS_1;
  Fhll.Mz = (Sr*Fl.Mz - Sl*Fr.Mz + Sl*Sr*(Ur.Mz - Ul.Mz)) * dS_1;
  Fhll.E  = (Sr*Fl.E  - Sl*Fr.E  + Sl*Sr*(Ur.E  - Ul.E )) * dS_1;
  Fhll.By = (Sr*Fl.By - Sl*Fr.By + Sl*Sr*(Ur.By - Ul.By)) * dS_1;
  Fhll.Bz = (Sr*Fl.Bz - Sl*Fr.Bz + Sl*Sr*(Ur.Bz - Ul.Bz)) * dS_1;


  if (switch_to_hll) {
    pFlux->d = Fhll.d;
    pFlux->Mx = Fhll.Mx;
    pFlux->My = Fhll.My;
    pFlux->Mz = Fhll.Mz;
    pFlux->E = Fhll.E;
    pFlux->By = Fhll.By;
    pFlux->Bz = Fhll.Bz;
		
    return;
  }

	
/*--- Step 3. ------------------------------------------------------------------
 * Compute fluxes based on wave speeds (Mignone et al. eqn 26)
 */
  if(Sl >= 0.0){
    pFlux->d  = Fl.d;
    pFlux->Mx = Fl.Mx;
    pFlux->My = Fl.My;
    pFlux->Mz = Fl.Mz;
    pFlux->E  = Fl.E;
    pFlux->By = Fl.By;
    pFlux->Bz = Fl.Bz;
		
    return;
   
  }
  else if(Sr <= 0.0){
    pFlux->d  = Fr.d;
    pFlux->Mx = Fr.Mx;
    pFlux->My = Fr.My;
    pFlux->Mz = Fr.Mz;
    pFlux->E  = Fr.E;
    pFlux->By = Fr.By;
    pFlux->Bz = Fr.Bz;
		
    return;
  }
  else {

    switch_to_hll = 0;
		
    /* Construct HLLC fluxes */
    vxl = Wl.Vx;
    vxr = Wr.Vx;
		
    Bx  = Bxi;
    Bys = Uhll.By;
    Bzs = Uhll.Bz;
		
    if (fabs(Bx) < 1.0e-12) {
      a  = Fhll.E;
      b  = - (Fhll.Mx + Uhll.E);
      c  = Uhll.Mx;
    } else {
      BtFBt = Uhll.By*Fhll.By + Uhll.Bz*Fhll.Bz;
      Bt2 = Uhll.By*Uhll.By + Uhll.Bz*Uhll.Bz;
      FBt2 = Fhll.By*Fhll.By + Fhll.Bz*Fhll.Bz;                
			
      a  = Fhll.E - BtFBt;
      b  = Bt2 + FBt2 - (Fhll.Mx + Uhll.E);
      c  = Uhll.Mx - BtFBt;
    }

    if (fabs(a) > 1.e-12){
      scrh = 1.0 + sqrt(1.0 - 4.0*a*c/(b*b));
      if (scrh != scrh) {
	switch_to_hll = 1;
      }
      vxs  = - 2.0*c/(b*scrh);
    } else {
      vxs = -c/b;
    }
    if ((vxs != vxs || vxs > 1.0) && (switch_to_hll == 0)) {
      switch_to_hll = 1;
    }
		
    if (switch_to_hll) {
      pFlux->d = Fhll.d;
      pFlux->Mx = Fhll.Mx;
      pFlux->My = Fhll.My;
      pFlux->Mz = Fhll.Mz;
      pFlux->E = Fhll.E;
      pFlux->By = Fhll.By;
      pFlux->Bz = Fhll.Bz;

      if (pFlux->d != pFlux->d) {
      }
		
      return;
    }
		
    if (fabs(Bx) < 1.0e-12) {
			
      /* -------------------------------
	 the value of vy and vz
	 is irrelevant in this case  
	 ------------------------------- */
			
      ps  = Fhll.Mx - Fhll.E*vxs;

      if (ps < 0) {
	switch_to_hll = 1;
      } else {			
	alpha_l = (Sl - vxl)/(Sl - vxs);
	alpha_r = (Sr - vxr)/(Sr - vxs);
			
	Usl.d = Ul.d*alpha_l;
	Usr.d = Ur.d*alpha_r;
			
	Usl.E = (Sl*Ul.E - Fl.E + ps*vxs)/(Sl - vxs);
	Usr.E = (Sr*Ur.E - Fr.E + ps*vxs)/(Sr - vxs);
			
	Usl.Mx = (Usl.E + ps)*vxs; 
	Usr.Mx = (Usr.E + ps)*vxs;
	Usl.My = Ul.My*alpha_l; 
	Usr.My = Ur.My*alpha_r; 
	Usl.Mz = Ul.Mz*alpha_l; 
	Usr.Mz = Ur.Mz*alpha_r;
			
	Usl.By = Ul.By*alpha_l;
	Usr.By = Ur.By*alpha_r;
	Usl.Bz = Ul.Bz*alpha_l;
	Usr.Bz = Ur.Bz*alpha_r;
      }
		
    } else {

      vys = (Bys*vxs - Fhll.By)/Bx;
      vzs = (Bzs*vxs - Fhll.Bz)/Bx;

      gammas_2 = vxs*vxs + vys*vys + vzs*vzs;
      gammas_2 = 1.0 - gammas_2;
      vBs = vxs*Bx + vys*Bys + vzs*Bzs;
			
      ps = (Bx*vBs - Fhll.E)*vxs + (Bx*Bx*gammas_2) + Fhll.Mx;
      if (ps < 0) {
	switch_to_hll = 1;
      } else {
			
	alpha_l = (Sl - vxl)/(Sl - vxs);
	alpha_r = (Sr - vxr)/(Sr - vxs);
			
	if (alpha_l != alpha_l) {
	  switch_to_hll = 1;
	}
			
	if (alpha_r != alpha_r) {
	  switch_to_hll = 1;

	}

	if (switch_to_hll == 0){
			
	  Usl.d = Ul.d*alpha_l;
	  Usr.d = Ur.d*alpha_r;
			
	  Usl.E = (Sl*Ul.E - Fl.E + ps*vxs - vBs*Bx)/(Sl - vxs);
	  Usr.E = (Sr*Ur.E - Fr.E + ps*vxs - vBs*Bx)/(Sr - vxs);
			
	  Usl.Mx = (Usl.E + ps)*vxs - vBs*Bx; 
	  Usr.Mx = (Usr.E + ps)*vxs - vBs*Bx;
	  Usl.My = (Sl*Ul.My - Fl.My - Bx*(Bys*gammas_2 + vBs*vys))/(Sl - vxs); 
	  Usr.My = (Sr*Ur.My - Fr.My - Bx*(Bys*gammas_2 + vBs*vys))/(Sr - vxs);  
	  Usl.Mz = (Sl*Ul.Mz - Fl.Mz - Bx*(Bzs*gammas_2 + vBs*vzs))/(Sl - vxs); 
	  Usr.Mz = (Sr*Ur.Mz - Fr.Mz - Bx*(Bzs*gammas_2 + vBs*vzs))/(Sr - vxs);
			
	  Usl.By = Usr.By = Bys;
	  Usl.Bz = Usr.Bz = Bzs;
	}
      }
    }

    if (switch_to_hll) {
      pFlux->d = Fhll.d;
      pFlux->Mx = Fhll.Mx;
      pFlux->My = Fhll.My;
      pFlux->Mz = Fhll.Mz;
      pFlux->E = Fhll.E;
      pFlux->By = Fhll.By;
      pFlux->Bz = Fhll.Bz;
		
      return;
    }
		
    /*  ----  Compute HLLC flux  ----  */
    if (vxs > 0.0) {
      pFlux->d  = Fl.d  + Sl*(Usl.d  - Ul.d );
      pFlux->Mx = Fl.Mx + Sl*(Usl.Mx - Ul.Mx);
      pFlux->My = Fl.My + Sl*(Usl.My - Ul.My);
      pFlux->Mz = Fl.Mz + Sl*(Usl.Mz - Ul.Mz);
      pFlux->E  = Fl.E  + Sl*(Usl.E  - Ul.E );
      pFlux->By = Fl.By + Sl*(Usl.By - Ul.By);
      pFlux->Bz = Fl.Bz + Sl*(Usl.Bz - Ul.Bz);
			
      if (pFlux->d != pFlux->d) {

	pFlux->d = Fhll.d;
	pFlux->Mx = Fhll.Mx;
	pFlux->My = Fhll.My;
	pFlux->Mz = Fhll.Mz;
	pFlux->E = Fhll.E;
	pFlux->By = Fhll.By;
	pFlux->Bz = Fhll.Bz;

      }
			
      return;
			
    } else {
	pFlux->d  = Fr.d  + Sr*(Usr.d  - Ur.d );
	pFlux->Mx = Fr.Mx + Sr*(Usr.Mx - Ur.Mx);
	pFlux->My = Fr.My + Sr*(Usr.My - Ur.My);
	pFlux->Mz = Fr.Mz + Sr*(Usr.Mz - Ur.Mz);
	pFlux->E  = Fr.E  + Sr*(Usr.E  - Ur.E );
	pFlux->By = Fr.By + Sr*(Usr.By - Ur.By);
	pFlux->Bz = Fr.Bz + Sr*(Usr.Bz - Ur.Bz);
			
      if (pFlux->d != pFlux->d) {

	pFlux->d = Fhll.d;
	pFlux->Mx = Fhll.Mx;
	pFlux->My = Fhll.My;
	pFlux->Mz = Fhll.Mz;
	pFlux->E = Fhll.E;
	pFlux->By = Fhll.By;
	pFlux->Bz = Fhll.Bz;

      }
			
      return;
    }

  }

}

/*! \fn void entropy_flux (const Cons1DS Ul, const Cons1DS Ur,
 *            const Prim1DS Wl, const Prim1DS Wr, const Real Bx, Real *pFlux)
 *  \brief Calculate entropy flux */
void entropy_flux (const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr, const Real Bx, Real *pFlux)
{
  Real Fl, Fr;
  Real USl, USr;
  Real WSl, WSr;
  Real Uhll, Fhll;
  Real Pl, Pr;
  Real Sl, Sr;
  Real Sla, Sra;
  Real dS_1;
  int wave_speed_fail;

  wave_speed_fail = 0;
	
/*--- Step 1. ------------------------------------------------------------------
 * Compute the max and min wave speeds used in Mignone 
 */
  getMaxSignalSpeeds_pluto(Wl,Wr,Bx,&Sl,&Sr);
	
  if (Sl != Sl) {
    wave_speed_fail = 1;
    printf("[hlle_sr_mhd]: NaN in Sl %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr =  1.0;
  }
	
  if (Sr != Sr) {
    wave_speed_fail = 1;
    printf("[hlle_sr_mhd]: NaN in Sr %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }
	
  if (Sl < -1.0) {
    wave_speed_fail = 1;
    printf("[hlle_sr_mhd]: Superluminal Sl %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }
  if (Sr > 1.0) {
    wave_speed_fail = 1;
    printf("[hlle_sr_mhd]: Superluminal Sr %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }

/*--- Step 1a. -----------------------------------------------------------------
 * If PLUTO wavespeeds are bad, fall back to the estimate used in ECHO
 */
  if (wave_speed_fail){
    getMaxSignalSpeeds_echo (Wl,Wr,Bx,&Sla,&Sra);
	
    if (Sla != Sla) {
      printf("[hlle_sr_mhd]: NaN in Sl %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra =  1.0;
    }
	
    if (Sra != Sra) {
      printf("[hlle_sr_mhd]: NaN in Sr %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra = 1.0;
    }
	
    if (Sla < -1.0) {
      printf("[hlle_sr_mhd]: Superluminal Sl %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra = 1.0;
    }
    if (Sra > 1.0) {
      printf("[hlle_sr_mhd]: Superluminal Sr %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra = 1.0;
    }

    Sl = Sla;
    Sr = Sra;

  }

  /* compute L/R fluxes */
  WSl = Wl.P*pow(Wl.d,1.0-Gamma);
  WSr = Wr.P*pow(Wr.d,1.0-Gamma);
  USl = WSl * Ul.d/Wl.d;
  USr = WSr * Ur.d/Wr.d;
  Fl = USl * Wl.Vx;
  Fr = USr * Wr.Vx;

  if(Sl >= 0.0){
    *pFlux = Fl;
    return;
  }
  else if(Sr <= 0.0){
    *pFlux = Fr;
    return;
  }
  else{
    /* Compute HLL average state */
    dS_1 = 1.0/(Sr - Sl);
    *pFlux = (Sr*Fl  - Sl*Fr  + Sl*Sr*(USr  - USl)) * dS_1;
    return;
  }
}

/*! \fn void flux_LR(Cons1DS U, Prim1DS W, Cons1DS *flux, Real Bx, Real* p)
 *  \brief Calculate LR flux. */
void flux_LR(Cons1DS U, Prim1DS W, Cons1DS *flux, Real Bx, Real* p)
{
  Real wtg2, pt, g, g2, g_1,g_2, h, gmmr, theta;
  Real bx, by, bz, vB, b2, Bmag2;
	
  /* calcUlate enthalpy */
	
  theta = W.P/W.d;
  gmmr = Gamma / Gamma_1;
	
  h = 1.0 + gmmr*theta;
	
  /* calcUlate gamma */
  g   = U.d/W.d;
  g2  = SQR(g);
  g_1 = 1.0/g;
  g_2 = 1.0/g2;
	
  pt = W.P;
  wtg2 = W.d*h*g2;
	
  vB = W.Vx*Bx + W.Vy*W.By + W.Vz*W.Bz;
  Bmag2 = SQR(Bx) + SQR(W.By) + SQR(W.Bz);
	
  bx = g*(  Bx*g_2 + vB*W.Vx);
  by = g*(W.By*g_2 + vB*W.Vy);
  bz = g*(W.Bz*g_2 + vB*W.Vz);
	
  b2 = Bmag2*g_2 + vB*vB;
	
  pt += 0.5*b2;
  wtg2 += b2*g2;
	
  flux->d  = U.d*W.Vx;
  flux->Mx = wtg2*W.Vx*W.Vx + pt;
  flux->My = wtg2*W.Vy*W.Vx;
  flux->Mz = wtg2*W.Vz*W.Vx;
  flux->E  = U.Mx;

  flux->Mx -= bx*bx;
  flux->My -= by*bx;
  flux->Mz -= bz*bx;
  flux->By = W.Vx*W.By - Bx*W.Vy;
  flux->Bz = W.Vx*W.Bz - Bx*W.Vz;
	
  *p = pt;
}

/*! \fn void getMaxSignalSpeeds_pluto(const Prim1DS Wl, const Prim1DS Wr,
 *			      const Real Bx, Real* low, Real* high)
 *  \brief
 */
void getMaxSignalSpeeds_pluto(const Prim1DS Wl, const Prim1DS Wr,
			      const Real Bx, Real* low, Real* high)
{
	
  Real lml,lmr;        /* smallest roots, Mignone Eq 55 */
  Real lpl,lpr;        /* largest roots, Mignone Eq 55 */
  Real al,ar;
	
  getVChar_pluto(Wl,Bx,&lml,&lpl);
  getVChar_pluto(Wr,Bx,&lmr,&lpr);
	
  *low =  MIN(lml, lmr);
  *high = MAX(lpl, lpr);
}

/*! \fn void getVChar_pluto(const Prim1DS W, const Real Bx, Real* lm, Real* lp)
 *  \brief
 */
void getVChar_pluto(const Prim1DS W, const Real Bx, Real* lm, Real* lp)
{
  Real rhoh,vsq,bsq;
  Real cssq,vasq,asq;
  Real Vx2,gamma,gamma2;
  Real Bx2,Bsq,vDotB,vDotBsq,b0,bx;
  Real w_1,a0,a1,a2,a3,a4,Q;
  Real scrh,scrh1,scrh2,eps2;
  Real a2_w,one_m_eps2,lambda[5];
  int iflag;
	
  rhoh = W.d + (Gamma/Gamma_1) * (W.P);
	
  Vx2 = SQR(W.Vx);
  vsq = Vx2 + SQR(W.Vy) + SQR(W.Vz);
  if (vsq > 1.0){
    /*printf("[getVChar]: |v|= %f > 1\n",vsq);*/	
    *lm = -1.0;
    *lp = 1.0;	
    return;		
  }
  gamma = 1.0 / sqrt(1 - vsq);    
  gamma2 = SQR(gamma);
	
  Bx2 = SQR(Bx);
  Bsq = Bx2 + SQR(W.By) + SQR(W.Bz);
  vDotB = W.Vx*Bx + W.Vy*W.By + W.Vz*W.Bz;
  vDotBsq = SQR(vDotB);
  b0 = gamma * vDotB;
  bx = Bx/gamma2 + W.Vx*vDotB;
  bsq = Bsq / gamma2 + SQR(vDotB);
	
  cssq = (Gamma * W.P) / (rhoh);
  vasq = bsq / (rhoh + bsq);
	
  if (bsq < 0.0) bsq = 0.0;
  if (cssq < 0.0) cssq = 0.0;
  if (cssq > 1.0) cssq = 1.0;
  if (vasq > 1.0) bsq = rhoh + bsq;
	
  if (vsq < 1.0e-12) {
    w_1  = 1.0/(rhoh + bsq);   
    eps2 = cssq + bsq*w_1*(1.0 - cssq);
    a0   = cssq*Bx*Bx*w_1;
    a1   = - a0 - eps2;
    scrh = a1*a1 - 4.0*a0;
    if (scrh < 0.0) scrh = 0.0;
		
    scrh = sqrt(0.5*(-a1 + sqrt(scrh)));
    *lp =  scrh;
    *lm = -scrh;
    return;
  }
	
  w_1 = 1.0/(rhoh + bsq);   
	
  if (Bx < 1.0e-14) {
		
    eps2  = cssq + bsq*w_1*(1.0 - cssq);
		
    scrh1 = (1.0 - eps2)*gamma2;
    scrh2 = cssq*vDotBsq*w_1 - eps2;
		
    a2  = scrh1 - scrh2;
    a1  = -2.0*W.Vx*scrh1;
    a0  = Vx2*scrh1 + scrh2;
		
    *lp = 0.5*(-a1 + sqrt(a1*a1 - 4.0*a2*a0))/a2;
    *lm = 0.5*(-a1 - sqrt(a1*a1 - 4.0*a2*a0))/a2;
		
    return;
  }
	
  scrh1 = bx;  /* -- this is bx/u0 -- */
  scrh2 = scrh1*scrh1;  
	
  a2_w       = cssq*w_1;
  eps2       = (cssq*rhoh + bsq)*w_1;
  one_m_eps2 = gamma2*rhoh*(1.0 - cssq)*w_1;
	
  /* ---------------------------------------
     Define coefficients for the quartic  
     --------------------------------------- */
	
  scrh = 2.0*(a2_w*vDotB*scrh1 - eps2*W.Vx);
  a4 = one_m_eps2 - a2_w*vDotBsq + eps2;
  a3 = - 4.0*W.Vx*one_m_eps2 + scrh;
  a2 =   6.0*Vx2*one_m_eps2 + a2_w*(vDotBsq - scrh2) + eps2*(Vx2 - 1.0);
  a1 = - 4.0*W.Vx*Vx2*one_m_eps2 - scrh;
  a0 = Vx2*Vx2*one_m_eps2 + a2_w*scrh2 - eps2*Vx2;
	
  if (a4 < 1.e-12){
    /*printPrim1D(W);*/
    printf("[MAX_CH_SPEED]: Can not divide by a4 in MAX_CH_SPEED\n");
		
    *lm = -1.0;
    *lp = 1.0;
		
    return;
  }
	
  scrh = 1.0/a4;
	
  a3 *= scrh;
  a2 *= scrh;
  a1 *= scrh;
  a0 *= scrh;
  iflag = QUARTIC(a3, a2, a1, a0, lambda);
	
  if (iflag){
    printf ("Can not find max speed:\n");
    /*SHOW(uprim,i);*/
    printf("QUARTIC: f(x) = %12.6e + x*(%12.6e + x*(%12.6e ",
	   a0*a4, a1*a4, a2*a4);
    printf("+ x*(%12.6e + x*%12.6e)))\n", a3*a4, a4);
    printf("[MAX_CH_SPEED]: Failed to find wave speeds");
		
    *lm = -1.0;
    *lp = 1.0;
		
    return;
  }
	
  *lp = MIN(1.0,MAX(lambda[3], lambda[2]));
  *lp = MIN(1.0,MAX(*lp, lambda[1]));
  *lp = MIN(1.0,MAX(*lp, lambda[0]));

  *lm = MAX(-1.0,MIN(lambda[3], lambda[2]));
  *lm = MAX(-1.0,MIN(*lm, lambda[1]));
  *lm = MAX(-1.0,MIN(*lm, lambda[0]));
	
  return;
	
}

/*! \fn void getMaxSignalSpeeds_echo (const Prim1DS Wl, const Prim1DS Wr,
 *			      const Real Bx, Real* low, Real* high)
 *  \brief
 */
void getMaxSignalSpeeds_echo (const Prim1DS Wl, const Prim1DS Wr,
			      const Real Bx, Real* low, Real* high)
{
	
  Real lml,lmr;        /* smallest roots, Mignone Eq 55 */
  Real lpl,lpr;        /* largest roots, Mignone Eq 55 */
  Real al,ar;
	
  getVChar_echo(Wl,Bx,&lml,&lpl);
  getVChar_echo(Wr,Bx,&lmr,&lpr);
	
  *low =  MIN(lml, lmr);
  *high = MAX(lpl, lpr);
}

/*! \fn void getVChar_echo(const Prim1DS W, const Real Bx, Real* lm, Real* lp)
 *  \brief
 */
void getVChar_echo(const Prim1DS W, const Real Bx, Real* lm, Real* lp)
{
  Real rhoh,vsq,bsq;
  Real cssq,vasq,asq;
  Real gamma,gamma2;
  Real Bsq,vDotB,b0,bx;
  Real tmp1,tmp2,tmp3,tmp4,tmp5;
  Real vm,vp;
	
  rhoh = W.d + (Gamma/Gamma_1) * (W.P);
	
  vsq = SQR(W.Vx) + SQR(W.Vy) + SQR(W.Vz);
  gamma = 1.0 / sqrt(1 - vsq);    
  gamma2 = SQR(gamma);
	
  Bsq = SQR(Bx) + SQR(W.By) + SQR(W.Bz);
  vDotB = W.Vx*Bx + W.Vy*W.By + W.Vz*W.Bz;
  b0 = gamma * vDotB;
  bx = Bx/gamma2 + W.Vx*vDotB;
  bsq = Bsq / gamma2 + SQR(vDotB);
	
  cssq = (Gamma * W.P) / (rhoh);
  vasq = bsq / (rhoh + bsq);
  asq = cssq + vasq - (cssq*vasq);
	
  if (cssq < 0.0) cssq = 0.0;
  if (vasq > 0.0) vasq = 0.0;
  if (asq < 0.0) asq = 0.0;
  if (cssq > 1.0) cssq = 1.0;
  if (vasq > 1.0) vasq = 1.0;
  if (asq > 1.0) asq = 1.0;
	
  tmp1 = (1.0 - asq);
  tmp2 = (1.0 - vsq);
  tmp3 = (1.0 - vsq*asq);
  tmp4 = SQR(W.Vx);
  tmp5 = 1.0 / tmp3;
	
  vm = tmp1*W.Vx - sqrt(asq*tmp2*(tmp3 - tmp1*tmp4));
  vp = tmp1*W.Vx + sqrt(asq*tmp2*(tmp3 - tmp1*tmp4));
  vm *=tmp5;
  vp *=tmp5;
	
  if (vp > vm) {
    *lm = vm;
    *lp = vp;
  } else {
    *lm = vp;
    *lp = vm;
  }
}

/* ******************************************** */
/*! \fn int QUARTIC (Real b, Real c, Real d, Real e, Real z[])
 *  \brief Solve a quartic equation */
/* 
 *
 * PURPOSE:
 *
 *   Solve a quartic equation in the form 
 *
 *    - z^4 + bz^3 + cz^2 + dz + e = 0
 *
 *   For its purpose, it is assumed that ALL 
 *   roots are real. This makes things faster.
 *
 *
 * ARGUMENTS
 *
 * - b, c,
 * - d, e  (IN)  = coefficient of the quartic
 *                 z^4 + bz^3 + cz^2 + dz + e = 0
 *
 * - z[]   (OUT) = a vector containing the 
 *                 (real) roots of the quartic
 *   
 *
 * REFERENCE:
 *
 * - http://www.1728.com/quartic2.htm 
 * 
 *
 */
/********************************************** */
int QUARTIC (Real b, Real c, Real d, 
             Real e, Real z[])
{
  int    n, ifail;
  Real b2, f, g, h;
  Real a2, a1, a0, u[4];
  Real p, q, r, s;
  static Real three_256 = 3.0/256.0;
  static Real one_64 = 1.0/64.0;
	
  b2 = b*b;
	
  f = c - b2*0.375;
  g = d + b2*b*0.125 - b*c*0.5;
  h = e - b2*b2*three_256 + 0.0625*b2*c - 0.25*b*d;
	
  a2 = 0.5*f;
  a1 = (f*f - 4.0*h)*0.0625;
  a0 = -g*g*one_64;
	
  ifail = CUBIC(a2, a1, a0, u);
	
  if (ifail)return(1);
	
  if (u[1] < 1.e-14){
		
    p = sqrt(u[2]);
    s = 0.25*b;
    z[0] = z[2] = - p - s;
    z[1] = z[3] = + p - s;
		
  }else{
		
    p = sqrt(u[1]);
    q = sqrt(u[2]);
		
    r = -0.125*g/(p*q);
    s =  0.25*b;
		
    z[0] = - p - q + r - s;
    z[1] =   p - q - r - s;
    z[2] = - p + q - r - s;
    z[3] =   p + q + r - s;
		
  }  
	
  /* ----------------------------------------------
     verify that cmax and cmin satisfy original 
     equation
     ---------------------------------------------- */  
	
  for (n = 0; n < 4; n++){
    s = e + z[n]*(d + z[n]*(c + z[n]*(b + z[n])));
    if (s != s) {
      printf ("Nan found in QUARTIC \n");
      return(1);
    }
    if (fabs(s) > 1.e-6) {
      printf ("Solution does not satisfy f(z) = 0; f(z) = %12.6e\n",s);
      return(1);
    }
  }
	
  return(0);
  /*  
      printf (" z: %f ; %f ; %f ; %f\n",z[0], z[1], z[2], z[3]);
  */
}
/* *************************************************** */
/*! \fn int CUBIC(Real b, Real c, Real d, Real z[])
 *  \brief Solve a cubic equation.
 *
 * PURPOSE:
 *
 *   Solve a cubic equation in the form 
 *
 *   -  z^3 + bz^2 + cz + d = 0
 *
 *   For its purpose, it is assumed that ALL 
 *   roots are real. This makes things faster.
 *
 *
 * ARGUMENTS
 *
 * - b, c, d (IN)  = coefficient of the cubic
 *                    z^3 + bz^2 + cz + d = 0
 *
 * - z[]   (OUT)   = a vector containing the 
 *                   (real) roots of the cubic.
 *                   Roots should be sorted
 *                   in increasing order.
 *   
 *
 * REFERENCE:
 *
 * - http://www.1728.com/cubic2.htm 
 *
 *
 *
 ***************************************************** */
int CUBIC(Real b, Real c, Real d, Real z[])
{
  Real b2, g2;
  Real f, g, h;
  Real i, i2, j, k, m, n, p;
  static Real one_3 = 1.0/3.0, one_27=1.0/27.0;
	
  b2 = b*b;
	
  /*  ----------------------------------------------
      the expression for f should be 
      f = c - b*b/3.0; however, to avoid negative
      round-off making h > 0.0 or g^2/4 - h < 0.0
      we let c --> c(1- 1.1e-16)
      ---------------------------------------------- */
	
  f  = c*(1.0 - 1.e-16) - b2*one_3;
  g  = b*(2.0*b2 - 9.0*c)*one_27 + d; 
  g2 = g*g;
  i2 = -f*f*f*one_27;
  h  = g2*0.25 - i2;
	
  /* --------------------------------------------
     Real roots are possible only when 
	 
     h <= 0 
     -------------------------------------------- */
	
  if (h > 1.e-12){
    printf ("Only one real root (%12.6e)!\n", h);
  }
  if (i2 < 0.0){
    /*
      printf ("i2 < 0.0 %12.6e\n",i2);
      return(1);
    */
    i2 = 0.0;
  }
	
  /* --------------------------------------
     i^2 must be >= g2*0.25
     -------------------------------------- */
	
  i = sqrt(i2);       /*  > 0   */
  j = pow(i, one_3);  /*  > 0   */
  k = -0.5*g/i;
	
  /*  this is to prevent unpleseant situation 
      where both g and i are close to zero       */
	
  k = (k < -1.0 ? -1.0:k);
  k = (k >  1.0 ?  1.0:k);
	
  k = acos(k)*one_3;       /*  pi/3 < k < 0 */
	
  m = cos(k);              /*   > 0   */
  n = sqrt(3.0)*sin(k);    /*   > 0   */
  p = -b*one_3;
	
  z[0] = -j*(m + n) + p;
  z[1] = -j*(m - n) + p;
  z[2] =  2.0*j*m + p;
	
  /* ------------------------------------------------------
     Since j, m, n > 0, it should follow that from
    
     z0 = -jm - jn + p
     z1 = -jm + jn + p
     z2 = 2jm + p
	 
     z2 is the greatest of the roots, while z0 is the 
     smallest one.
     ------------------------------------------------------ */
	
  return(0);
}

#endif /* MHD */

#endif /* HLLC_FLUX */
#endif /* SPECIAL_RELATIVITY */
