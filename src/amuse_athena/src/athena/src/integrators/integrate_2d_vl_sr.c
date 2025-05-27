#include "../copyright.h"
/*============================================================================*/
/*! \file integrate_2d_vl_sr.c
 *  \brief Integrate MHD equations using 2D version of the directionally
 *   unsplit MUSCL-Hancock (VL) integrator.
 *
 * PURPOSE: Integrate MHD equations using 2D version of the directionally
 *   unsplit MUSCL-Hancock (VL) integrator.  The variables updated are:
 *   -  U.[d,M1,M2,M3,E,B1c,B2c,B3c,s] -- where U is of type ConsS
 *   -  B1i, B2i  -- interface magnetic field
 *   Also adds gravitational source terms, self-gravity, and the H-correction
 *   of Sanders et al.
 *
 * REFERENCE: 
 * - J.M Stone & T.A. Gardiner, "A simple, unsplit Godunov method
 *   for multidimensional MHD", NewA 14, 139 (2009)
 *
 * - R. Sanders, E. Morano, & M.-C. Druguet, "Multidimensional dissipation for
 *   upwind schemes: stability and applications to gas dynamics", JCP, 145, 511
 *   (1998)
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - integrate_2d_vl()
 * - integrate_destruct_2d()
 * - integrate_init_2d() */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef SPECIAL_RELATIVITY

#if defined(VL_INTEGRATOR) && defined(CARTESIAN)

#ifdef MHD
#define USE_ENTROPY_FIX
#endif /* MHD */

/* The L/R states of primitive variables and fluxes at each cell face */
static Prim1DS **Wl_x1Face=NULL, **Wr_x1Face=NULL;
static Prim1DS **Wl_x2Face=NULL, **Wr_x2Face=NULL;
static Cons1DS **x1Flux=NULL, **x2Flux=NULL;
#ifdef FIRST_ORDER_FLUX_CORRECTION
static Cons1DS **x1FluxP=NULL, **x2FluxP=NULL;
static Real **emf3P=NULL;
#endif
#ifdef USE_ENTROPY_FIX
static Real **S=NULL,**Shalf=NULL;
static Real **x1FluxS=NULL,**x2FluxS=NULL;
#ifdef FIRST_ORDER_FLUX_CORRECTION
static Real **x1FluxSP=NULL,**x2FluxSP=NULL;
#endif
#endif

/* The interface magnetic fields and emfs */
#ifdef MHD
static Real **B1_x1Face=NULL, **B2_x2Face=NULL;
static Real **emf3=NULL, **emf3_cc=NULL;
#endif /* MHD */

/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W1d=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL, *Ul=NULL, *Ur=NULL;

/* primitive variables at t^{n} */
static PrimS **W=NULL;

/* conserved and primitive variables at t^{n+1/2} computed in predict step */
static ConsS **Uhalf=NULL;
static PrimS **Whalf=NULL;

/* variables needed for H-correction of Sanders et al (1998) */
extern Real etah;
#ifdef H_CORRECTION
static Real **eta1=NULL, **eta2=NULL;
#endif

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   integrate_emf3_corner() - the upwind CT method of Gardiner & Stone (2005) 
 *   FixCell() - apply first-order correction to one cell
 *============================================================================*/
#ifdef MHD
static void integrate_emf3_corner(GridS *pG);
#endif
#ifdef FIRST_ORDER_FLUX_CORRECTION
static void FixCell(GridS *pG, Int3Vect);
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void integrate_2d_vl(DomainS *pD)
 *  \brief van Leer unsplit integrator in 2D. 
 *
 *   The numbering of steps follows the numbering in the 3D version.
 *   NOT ALL STEPS ARE NEEDED IN 2D.
 */
void integrate_2d_vl(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  Real dtodx1=pG->dt/pG->dx1, dtodx2=pG->dt/pG->dx2;
  Real hdtodx1 = 0.5*dtodx1, hdtodx2 = 0.5*dtodx2;
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int ks = pG->ks;
  int cart_x1 = 1, cart_x2 = 2, cart_x3 = 3;
  Real pmhalf,pmnew;
  Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic,Bx;
#if (NSCALARS > 0)
  int n;
#endif
#ifdef SELF_GRAVITY
  Real gxl,gxr,gyl,gyr,flx_m1l,flx_m1r,flx_m2l,flx_m2r;
#endif
#ifdef H_CORRECTION
  Real cfr,cfl,lambdar,lambdal;
#endif
#ifdef STATIC_MESH_REFINEMENT
  int ncg,npg,dim;
  int ii,ics,ice,jj,jcs,jce,ips,ipe,jps,jpe;
#endif
#ifdef FIRST_ORDER_FLUX_CORRECTION
  ConsS Ucheck;
  PrimS Wcheck;
  int flag_cell=0,negd=0,negP=0,superl=0,NaNFlux=0;
  int entropy,final,fail;
  Real Vsq;
  Int3Vect BadCell;
#endif
  int il=is-(nghost-1), iu=ie+(nghost-1);
  int jl=js-(nghost-1), ju=je+(nghost-1);

/* Set etah=0 so first calls to flux functions do not use H-correction */
  etah = 0.0;

  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      Uhalf[j][i] = pG->U[ks][j][i];
      W[j][i] = Cons_to_Prim(&(pG->U[ks][j][i]));
#ifdef MHD
      B1_x1Face[j][i] = pG->B1i[ks][j][i]; 
      B2_x2Face[j][i] = pG->B2i[ks][j][i]; 
#endif /* MHD */
#ifdef USE_ENTROPY_FIX
      S[j][i] = W[j][i].P * pow(W[j][i].d,1.0-Gamma);
      S[j][i]*= pG->U[ks][j][i].d / W[j][i].d;
      Shalf[j][i] = S[j][i];
#endif
    }
  }

/*=== STEP 1: Compute first-order fluxes at t^{n} in x1-direction ============*/
/* No source terms are needed since there is no temporal evolution */

/*--- Step 1a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * W1d = (d, V1, V2, V3, P, B2c, B3c, s[n])
 */

  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      W1d[i].d  = W[j][i].d;
      W1d[i].Vx = W[j][i].V1;
      W1d[i].Vy = W[j][i].V2;
      W1d[i].Vz = W[j][i].V3;
      W1d[i].P  = W[j][i].P;
#ifdef MHD
      W1d[i].By = W[j][i].B2c;
      W1d[i].Bz = W[j][i].B3c;
      Bxc[i] = W[j][i].B1c;
      Bxi[i] = pG->B1i[ks][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) W1d[i].s[n] = W[j][i].s[n];
#endif
    }

/*--- Step 1b ------------------------------------------------------------------
 * Compute first-order L/R states */
    for (i=il; i<=ie+nghost; i++) {
      Wl[i] = W1d[i-1];
      Wr[i] = W1d[i  ];

/* Compute U from W in case Pfloor used in Cons1D_to_Prim1D */
      Ul[i] = Prim1D_to_Cons1D(&Wl[i], &Bxc[i-1]);
      Ur[i] = Prim1D_to_Cons1D(&Wr[i], &Bxc[i  ]);
    }

/*--- Step 1c ------------------------------------------------------------------
 * No source terms needed */

/*--- Step 1d ------------------------------------------------------------------
 * Compute flux in x1-direction */

    for (i=il; i<=ie+nghost; i++) {
      fluxes(Ul[i],Ur[i],Wl[i],Wr[i],Bxi[i],&x1Flux[j][i]);
#ifdef USE_ENTROPY_FIX
      entropy_flux(Ul[i],Ur[i],Wl[i],Wr[i],Bxi[i],&x1FluxS[j][i]);
#endif
    }
  }

/*=== STEP 2: Compute first-order fluxes at t^{n} in x2-direction ============*/
/* No source terms are needed since there is no temporal evolution */

/*--- Step 2a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;
 * W1d = (d, V2, V3, V1, P, B3c, B1c, s[n])
 */

  for (i=is-nghost; i<=ie+nghost; i++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      W1d[j].d  = W[j][i].d;
      W1d[j].Vx = W[j][i].V2;
      W1d[j].Vy = W[j][i].V3;
      W1d[j].Vz = W[j][i].V1;
      W1d[j].P  = W[j][i].P;
#ifdef MHD
      W1d[j].By = W[j][i].B3c;
      W1d[j].Bz = W[j][i].B1c;
      Bxc[j] = W[j][i].B2c;
      Bxi[j] = pG->B2i[ks][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) W1d[j].s[n] = W[j][i].s[n];
#endif
    }

/*--- Step 2b ------------------------------------------------------------------
 * Compute first-order L/R states */

    for (j=jl; j<=je+nghost; j++) {
      Wl[j] = W1d[j-1];
      Wr[j] = W1d[j  ];

/* Compute U from W in case Pfloor used in Cons1D_to_Prim1D */
      Ul[j] = Prim1D_to_Cons1D(&Wl[j], &Bxc[j-1]);
      Ur[j] = Prim1D_to_Cons1D(&Wr[j], &Bxc[j  ]);
    }

/*--- Step 2c ------------------------------------------------------------------
 * No source terms needed */

/*--- Step 2d ------------------------------------------------------------------
 * Compute flux in x2-direction */

    for (j=jl; j<=je+nghost; j++) {
      fluxes(Ul[j],Ur[j],Wl[j],Wr[j],Bxi[j],&x2Flux[j][i]);
#ifdef USE_ENTROPY_FIX
      entropy_flux(Ul[j],Ur[j],Wl[j],Wr[j],Bxi[j],&x2FluxS[j][i]);
#endif
    }
  }

/*=== STEP 3: Not needed in 2D ===*/

/*=== STEP 4:  Update face-centered B for 0.5*dt =============================*/

/*--- Step 4a ------------------------------------------------------------------
 * Calculate the cell centered value of emf1,2,3 at t^{n} and integrate
 * to corner.
 */

#ifdef MHD
  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      emf3_cc[j][i] = (W[j][i].B1c*W[j][i].V2 - W[j][i].B2c*W[j][i].V1);
    }
  }
  integrate_emf3_corner(pG);

/*--- Step 4b ------------------------------------------------------------------
 * Update the interface magnetic fields using CT for a half time step.
 */

  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      B1_x1Face[j][i] -= hdtodx2*(emf3[j+1][i  ] - emf3[j][i]);
      B2_x2Face[j][i] += hdtodx1*(emf3[j  ][i+1] - emf3[j][i]);
    }
    B1_x1Face[j][iu+1] -= hdtodx2*(emf3[j+1][iu+1]-emf3[j][iu+1]);
  }
  for (i=il; i<=iu; i++) {
    B2_x2Face[ju+1][i] += hdtodx1*(emf3[ju+1][i+1]-emf3[ju+1][i]);
  }

/*--- Step 4c ------------------------------------------------------------------
 * Compute cell-centered magnetic fields at half-timestep from average of
 * face-centered fields.
 */

  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      Uhalf[j][i].B1c = 0.5*(B1_x1Face[j][i] + B1_x1Face[j][i+1]);
      Uhalf[j][i].B2c = 0.5*(B2_x2Face[j][i] + B2_x2Face[j+1][i]);
    }
  }
#endif /* MHD */

/*=== STEP 5: Update cell-centered variables to half-timestep ================*/

/*--- Step 5a ------------------------------------------------------------------
 * Update cell-centered variables (including B3c) to half-timestep with x1Flux
 */

  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      Uhalf[j][i].d   -= hdtodx1*(x1Flux[j][i+1].d  - x1Flux[j][i].d );
      Uhalf[j][i].M1  -= hdtodx1*(x1Flux[j][i+1].Mx - x1Flux[j][i].Mx);
      Uhalf[j][i].M2  -= hdtodx1*(x1Flux[j][i+1].My - x1Flux[j][i].My);
      Uhalf[j][i].M3  -= hdtodx1*(x1Flux[j][i+1].Mz - x1Flux[j][i].Mz);
#ifndef BAROTROPIC
      Uhalf[j][i].E   -= hdtodx1*(x1Flux[j][i+1].E  - x1Flux[j][i].E );
#endif /* BAROTROPIC */
#ifdef MHD
      Uhalf[j][i].B3c -= hdtodx1*(x1Flux[j][i+1].Bz - x1Flux[j][i].Bz);
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        Uhalf[j][i].s[n] -= hdtodx1*(x1Flux[j][i+1].s[n] - x1Flux[j][i].s[n]);
#endif
#ifdef USE_ENTROPY_FIX
      Shalf[j][i]   -= hdtodx1*(x1FluxS[j][i+1]  - x1FluxS[j][i]);
#endif
    }
  }

/*--- Step 5b ------------------------------------------------------------------
 * Update cell-centered variables (including B3c) to half-timestep with x2Flux
 */

  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      Uhalf[j][i].d   -= hdtodx2*(x2Flux[j+1][i].d  - x2Flux[j][i].d );
      Uhalf[j][i].M1  -= hdtodx2*(x2Flux[j+1][i].Mz - x2Flux[j][i].Mz);
      Uhalf[j][i].M2  -= hdtodx2*(x2Flux[j+1][i].Mx - x2Flux[j][i].Mx);
      Uhalf[j][i].M3  -= hdtodx2*(x2Flux[j+1][i].My - x2Flux[j][i].My);
#ifndef BAROTROPIC
      Uhalf[j][i].E   -= hdtodx2*(x2Flux[j+1][i].E  - x2Flux[j][i].E );
#endif /* BAROTROPIC */
#ifdef MHD
      Uhalf[j][i].B3c -= hdtodx2*(x2Flux[j+1][i].By - x2Flux[j][i].By);
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        Uhalf[j][i].s[n] -= hdtodx2*(x2Flux[j+1][i].s[n] - x2Flux[j][i].s[n]);
#endif
#ifdef USE_ENTROPY_FIX
      Shalf[j][i]   -= hdtodx2*(x2FluxS[j+1][i]  - x2FluxS[j][i] );
#endif
    }
  }

#ifdef FIRST_ORDER_FLUX_CORRECTION
/*--- Step 5d ------------------------------------------------------------------
 * With first-order flux correction, check predict fluxes & emf3 for NaN's
 * and save predict fluxes & emf3 for first-order flux correction
 */

  NaNFlux = 0;
  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      if ((x1Flux[j][i].d  != x1Flux[j][i].d)  ||
#ifndef BAROTROPIC
          (x1Flux[j][i].E  != x1Flux[j][i].E)  ||
#endif
#ifdef MHD
          (x1Flux[j][i].By != x1Flux[j][i].By) ||
          (x1Flux[j][i].Bz != x1Flux[j][i].Bz) ||
#endif
          (x1Flux[j][i].Mx != x1Flux[j][i].Mx) ||
          (x1Flux[j][i].My != x1Flux[j][i].My) ||
          (x1Flux[j][i].Mz != x1Flux[j][i].Mz)) {
        NaNFlux++;
      }          
      if ((x2Flux[j][i].d  != x2Flux[j][i].d)  ||
#ifndef BAROTROPIC
          (x2Flux[j][i].E  != x2Flux[j][i].E)  ||
#endif
#ifdef MHD
          (x2Flux[j][i].By != x2Flux[j][i].By) ||
          (x2Flux[j][i].Bz != x2Flux[j][i].Bz) ||
#endif
          (x2Flux[j][i].Mx != x2Flux[j][i].Mx) ||
          (x2Flux[j][i].My != x2Flux[j][i].My) ||
          (x2Flux[j][i].Mz != x2Flux[j][i].Mz)) {
        NaNFlux++;
      }
      x1FluxP[j][i] = x1Flux[j][i];
      x2FluxP[j][i] = x2Flux[j][i];
#ifdef MHD
      emf3P[j][i] = emf3[j][i];
#endif
#ifdef USE_ENTROPY_FIX
      x1FluxSP[j][i] = x1FluxS[j][i];
      x2FluxSP[j][i] = x2FluxS[j][i];
#endif
    }
  }
  if (NaNFlux != 0) {
    ath_error("[Step5d] %i first-order fluxes are NaN!\n",NaNFlux);
  }
#endif

/*=== STEP 6: Add source terms to predict values at half-timestep ============*/

/*--- Step 6a ------------------------------------------------------------------
 * Add source terms from a static gravitational potential for 0.5*dt to predict
 * step.  To improve conservation of total energy, we average the energy
 * source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phic = (*StaticGravPot)( x1,             x2,x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

        Uhalf[j][i].M1 -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;
#ifndef BAROTROPIC
        Uhalf[j][i].E -= hdtodx1*(x1Flux[j][i  ].d*(phic - phil)
                           + x1Flux[j][i+1].d*(phir - phic));
#endif
        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        Uhalf[j][i].M2 -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;
#ifndef BAROTROPIC
        Uhalf[j][i].E -= hdtodx2*(x2Flux[j  ][i].d*(phic - phil)
                                + x2Flux[j+1][i].d*(phir - phic));
#endif
      }
    }
  }


/*=== STEP 7: Conserved->Primitive variable inversion at t^{n+1/2} ===========*/
        
/* Invert conserved variables at t^{n+1/2} to primitive variables. With FOFC, 
 * check if cell-centered d < 0, P< 0, or v^2 > 1. With Entropy fix, correct
 * by computing new primitive state using the entropy variable, otherwise
 * correct by switching back to values at beginning of step, rendering update
 * first order in time for that cell.
 */

#ifdef FIRST_ORDER_FLUX_CORRECTION        
  negd = 0;
  negP = 0;
  superl = 0;
  flag_cell = 0;
#endif
  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      Whalf[j][i] = check_Prim(&(Uhalf[j][i]));
#ifdef FIRST_ORDER_FLUX_CORRECTION   
      if (Whalf[j][i].d < 0.0) {
        flag_cell = 1;
        negd++;
      }
      if (Whalf[j][i].P < 0.0) {
        flag_cell = 1;
        negP++;
      }
      Vsq = SQR(Whalf[j][i].V1) + SQR(Whalf[j][i].V2) + SQR(Whalf[j][i].V3);
      if (Vsq > 1.0) {
        flag_cell = 1;
        superl++;
      }
      if (flag_cell != 0) {
#ifdef USE_ENTROPY_FIX
	Wcheck = entropy_fix (&(Uhalf[j][i]),&(Shalf[j][i]));
	Vsq = SQR(Wcheck.V1) + SQR(Wcheck.V2) + SQR(Wcheck.V3);
	if (Wcheck.d > 0.0 && Wcheck.P > 0.0 && Vsq < 1.0){
	  entropy++;
	  Whalf[j][i].d = Wcheck.d;
	  Whalf[j][i].P = Wcheck.P;
	  Whalf[j][i].V1 = Wcheck.V1;
	  Whalf[j][i].V2 = Wcheck.V2;
	  Whalf[j][i].V3 = Wcheck.V3;
	  flag_cell=0;
	} else {
#endif /* USE_ENTROPY_FIX */
	  Whalf[j][i].d = W[j][i].d;
	  Whalf[j][i].V1 = W[j][i].V1;
	  Whalf[j][i].V2 = W[j][i].V2;
	  Whalf[j][i].V3 = W[j][i].V3;
	  Whalf[j][i].P = W[j][i].P;
	  flag_cell=0;
#ifdef USE_ENTROPY_FIX
	}
#endif /* USE_ENTROPY_FIX */
      }
#endif
    }
  }

#ifdef FIRST_ORDER_FLUX_CORRECTION  
  if (negd > 0 || negP > 0 || superl > 0){
    printf("[Step7]: %i cells had d<0; %i cells had P<0\n",negd,negP);
    printf("[Step7]: %i cells had v>1 at t_half\n",superl);
#ifdef USE_ENTROPY_FIX
    printf("[Step7]: %i cells fixed using entropy at t_half\n",entropy);
#endif /* USE_ENTROPY_FIX */
  }
#endif 

/*=== STEP 8: Compute second-order L/R x1-interface states ===================*/

/*--- Step 8a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;
 * W1d = (d, V1, V2, V3, P, B2c, B3c, s[n])
 */

  for (j=js-1; j<=je+1; j++) {
    for (i=il; i<=iu; i++) {
      W1d[i].d  = Whalf[j][i].d;
      W1d[i].Vx = Whalf[j][i].V1;
      W1d[i].Vy = Whalf[j][i].V2;
      W1d[i].Vz = Whalf[j][i].V3;
#ifndef BAROTROPIC
      W1d[i].P  = Whalf[j][i].P;
#endif /* BAROTROPIC */
#ifdef MHD
      W1d[i].By = Whalf[j][i].B2c;
      W1d[i].Bz = Whalf[j][i].B3c;
      Bxc[i] = Whalf[j][i].B1c;
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) W1d[i].s[n] = Whalf[j][i].s[n];
#endif /* NSCALARS */
    }

/*--- Step 8b ------------------------------------------------------------------
 * Compute L/R states on x1-interfaces, store into arrays
 */

    lr_states(pG,W1d,Bxc,pG->dt,pG->dx1,is,ie,Wl,Wr,cart_x1);

#ifdef FIRST_ORDER_FLUX_CORRECTION
    for (i=il; i<=iu; i++) {
      Vsq = SQR(Wl[i].Vx) + SQR(Wl[i].Vy) + SQR(Wl[i].Vz);
      if (Vsq > 1.0){
        Wl[i] = W1d[i];
        Wr[i] = W1d[i];
      }
      Vsq = SQR(Wr[i].Vx) + SQR(Wr[i].Vy) + SQR(Wr[i].Vz);
      if (Vsq > 1.0){
        Wl[i] = W1d[i];
        Wr[i] = W1d[i];
      }
    }
#endif

    for (i=is; i<=ie+1; i++) {
      Wl_x1Face[j][i] = Wl[i];
      Wr_x1Face[j][i] = Wr[i];
    }
  }



/*=== STEP 9: Compute second-order L/R x2-interface states ===================*/

/*--- Step 9a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;
 * W1d = (d, V2, V3, V1, P, B3c, B1c, s[n])
 */

  for (i=is-1; i<=ie+1; i++) {
    for (j=jl; j<=ju; j++) {
      W1d[j].d  = Whalf[j][i].d;
      W1d[j].Vx = Whalf[j][i].V2;
      W1d[j].Vy = Whalf[j][i].V3;
      W1d[j].Vz = Whalf[j][i].V1;
#ifndef BAROTROPIC
      W1d[j].P  = Whalf[j][i].P;
#endif /* BAROTROPIC */
#ifdef MHD
      W1d[j].By = Whalf[j][i].B3c;
      W1d[j].Bz = Whalf[j][i].B1c;
      Bxc[j] = Whalf[j][i].B2c;
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) W1d[j].s[n] = Whalf[j][i].s[n];
#endif /* NSCALARS */
    }

/*--- Step 9b ------------------------------------------------------------------
 * Compute L/R states on x1-interfaces, store into arrays
 */

    lr_states(pG,W1d,Bxc,pG->dt,pG->dx2,js,je,Wl,Wr,cart_x2);

#ifdef FIRST_ORDER_FLUX_CORRECTION
    for (j=jl; j<=ju; j++) {
      Vsq = SQR(Wl[j].Vx) + SQR(Wl[j].Vy) + SQR(Wl[j].Vz);
      if (Vsq > 1.0){
        Wl[j] = W1d[j];
        Wr[j] = W1d[j];
      }
      Vsq = SQR(Wr[j].Vx) + SQR(Wr[j].Vy) + SQR(Wr[j].Vz);
      if (Vsq > 1.0){
        Wl[j] = W1d[j];
        Wr[j] = W1d[j];
      }
    }
#endif

    for (j=js; j<=je+1; j++) {
      Wl_x2Face[j][i] = Wl[j];
      Wr_x2Face[j][i] = Wr[j];
    }
  }

/*=== STEP 10: Not needed in 2D ===*/

/*=== STEP 11: Compute 2D x1-Flux, x2-Flux, ==================================*/

/*--- Step 11a -----------------------------------------------------------------
 * Compute maximum wavespeeds in multidimensions (eta in eq. 10 from Sanders et
 *  al. (1998)) for H-correction, if needed.
 */

#ifdef H_CORRECTION
  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=iu; i++) {
#ifdef MHD
      Bx = B1_x1Face[j][i];
#endif
      cfr = cfast(&(Ur_x1Face[j][i]),&Bx);
      cfl = cfast(&(Ul_x1Face[j][i]),&Bx);
      lambdar = Ur_x1Face[j][i].Mx/Ur_x1Face[j][i].d + cfr;
      lambdal = Ul_x1Face[j][i].Mx/Ul_x1Face[j][i].d - cfl;
      eta1[j][i] = 0.5*fabs(lambdar - lambdal);
    }
  }

  for (j=js-1; j<=ju; j++) {
    for (i=is-1; i<=ie+1; i++) {
#ifdef MHD
      Bx = B2_x2Face[j][i];
#endif
      cfr = cfast(&(Ur_x2Face[j][i]),&Bx);
      cfl = cfast(&(Ul_x2Face[j][i]),&Bx);
      lambdar = Ur_x2Face[j][i].Mx/Ur_x2Face[j][i].d + cfr;
      lambdal = Ul_x2Face[j][i].Mx/Ul_x2Face[j][i].d - cfl;
      eta2[j][i] = 0.5*fabs(lambdar - lambdal);
    }
  }
#endif /* H_CORRECTION */

/*--- Step 11b -----------------------------------------------------------------
 * Compute second-order fluxes in x1-direction
 */

  for (j=js-1; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
#ifdef H_CORRECTION
      etah = MAX(eta2[j][i-1],eta2[j][i]);
      etah = MAX(etah,eta2[j+1][i-1]);
      etah = MAX(etah,eta2[j+1][i  ]);
      etah = MAX(etah,eta1[j  ][i  ]);
#endif /* H_CORRECTION */
#ifdef MHD
      Bx = B1_x1Face[j][i];
#endif
      Ul[i] = Prim1D_to_Cons1D(&Wl_x1Face[j][i],&Bx);
      Ur[i] = Prim1D_to_Cons1D(&Wr_x1Face[j][i],&Bx);

      fluxes(Ul[i],Ur[i],Wl_x1Face[j][i],Wr_x1Face[j][i],Bx,&x1Flux[j][i]);

#ifdef USE_ENTROPY_FIX
      entropy_flux(Ul[i],          Ur[i],
		   Wl_x1Face[j][i],Wr_x1Face[j][i],
		   Bx,             &x1FluxS[j][i]);
#endif

#ifdef FIRST_ORDER_FLUX_CORRECTION
/* revert to predictor flux if this flux Nan'ed */
      if ((x1Flux[j][i].d  != x1Flux[j][i].d)  ||
#ifndef BAROTROPIC
          (x1Flux[j][i].E  != x1Flux[j][i].E)  ||
#endif
#ifdef MHD
          (x1Flux[j][i].By != x1Flux[j][i].By) ||
          (x1Flux[j][i].Bz != x1Flux[j][i].Bz) ||
#endif
          (x1Flux[j][i].Mx != x1Flux[j][i].Mx) ||
          (x1Flux[j][i].My != x1Flux[j][i].My) ||
          (x1Flux[j][i].Mz != x1Flux[j][i].Mz)) {
        x1Flux[j][i] = x1FluxP[j][i];
        NaNFlux++;
      }
#endif

    }
  }

/*--- Step 11c -----------------------------------------------------------------
 * Compute second-order fluxes in x2-direction
 */

  for (j=js; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
#ifdef H_CORRECTION
      etah = MAX(eta1[j-1][i],eta1[j][i]);
      etah = MAX(etah,eta1[j-1][i+1]);
      etah = MAX(etah,eta1[j  ][i+1]);
      etah = MAX(etah,eta2[j  ][i]);
#endif /* H_CORRECTION */
#ifdef MHD
      Bx = B2_x2Face[j][i];
#endif
      Ul[i] = Prim1D_to_Cons1D(&Wl_x2Face[j][i],&Bx);
      Ur[i] = Prim1D_to_Cons1D(&Wr_x2Face[j][i],&Bx);

      fluxes(Ul[i],Ur[i],Wl_x2Face[j][i],Wr_x2Face[j][i],Bx,&x2Flux[j][i]);

#ifdef USE_ENTROPY_FIX
      entropy_flux(Ul[i],          Ur[i],
		   Wl_x2Face[j][i],Wr_x2Face[j][i],
		   Bx,             &x2FluxS[j][i]);
#endif

#ifdef FIRST_ORDER_FLUX_CORRECTION
/* revert to predictor flux if this flux NaN'ed */
      if ((x2Flux[j][i].d  != x2Flux[j][i].d)  ||
#ifndef BAROTROPIC
          (x2Flux[j][i].E  != x2Flux[j][i].E)  ||
#endif
#ifdef MHD
          (x2Flux[j][i].By != x2Flux[j][i].By) ||
          (x2Flux[j][i].Bz != x2Flux[j][i].Bz) ||
#endif
          (x2Flux[j][i].Mx != x2Flux[j][i].Mx) ||
          (x2Flux[j][i].My != x2Flux[j][i].My) ||
          (x2Flux[j][i].Mz != x2Flux[j][i].Mz)) {
        x2Flux[j][i] = x2FluxP[j][i];
        NaNFlux++;
      }
#endif

    }
  }

#ifdef FIRST_ORDER_FLUX_CORRECTION
  if (NaNFlux != 0) {
    printf("[Step11] %i second-order fluxes replaced\n",NaNFlux);
    NaNFlux=0;
  }
#endif
  

/*=== STEP 12: Update face-centered B for a full timestep ====================*/
        
/*--- Step 12a -----------------------------------------------------------------
 * Calculate the cell centered value of emf1,2,3 at the half-time-step.
 */

#ifdef MHD
  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
      emf3_cc[j][i] = (Whalf[j][i].B1c*Whalf[j][i].V2 - 
		       Whalf[j][i].B2c*Whalf[j][i].V1);
    }
  }

/*--- Step 12b -----------------------------------------------------------------
 * Integrate emf*^{n+1/2} to the grid cell corners
 */

  integrate_emf3_corner(pG);

/*--- Step 12c -----------------------------------------------------------------
 * Update the interface magnetic fields using CT for a full time step.
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->B1i[ks][j][i] -= dtodx2*(emf3[j+1][i  ] - emf3[j][i]);
      pG->B2i[ks][j][i] += dtodx1*(emf3[j  ][i+1] - emf3[j][i]);
    }
    pG->B1i[ks][j][ie+1] -= dtodx2*(emf3[j+1][ie+1] - emf3[j][ie+1]);
  }
  for (i=is; i<=ie; i++) {
    pG->B2i[ks][je+1][i] += dtodx1*(emf3[je+1][i+1] - emf3[je+1][i]);
  }

/*--- Step 12d -----------------------------------------------------------------
 * Set cell centered magnetic fields to average of updated face centered fields.
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].B1c = 0.5*(pG->B1i[ks][j][i]+pG->B1i[ks][j][i+1]);
      pG->U[ks][j][i].B2c = 0.5*(pG->B2i[ks][j][i]+pG->B2i[ks][j+1][i]);
    }
  }
#endif /* MHD */

/*=== STEP 13: Add source terms for a full timestep using n+1/2 states =======*/
       
/*--- Step 13a -----------------------------------------------------------------
 * Add gravitational source terms due to a Static Potential
 * To improve conservation of total energy, we average the energy
 * source term computed at cell faces.
 *    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phic = (*StaticGravPot)( x1,             x2,x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

        pG->U[ks][j][i].M1 -= dtodx1*(phir-phil)*Uhalf[j][i].d;
#ifndef BAROTROPIC
        pG->U[ks][j][i].E -= dtodx1*(x1Flux[j][i  ].d*(phic - phil)
                                   + x1Flux[j][i+1].d*(phir - phic));
#endif
        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        pG->U[ks][j][i].M2 -= dtodx2*(phir-phil)*Uhalf[j][i].d;
#ifndef BAROTROPIC
        pG->U[ks][j][i].E -= dtodx2*(x2Flux[j  ][i].d*(phic - phil)
                                   + x2Flux[j+1][i].d*(phir - phic));
#endif
      }
    }
  }

/*=== STEP 14: Update cell-centered values for a full timestep ===============*/

/*--- Step 14a -----------------------------------------------------------------
 * Update cell-centered variables in pG (including B3c) using 2D x1-Fluxes
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].d   -= dtodx1*(x1Flux[j][i+1].d  - x1Flux[j][i].d );
      pG->U[ks][j][i].M1  -= dtodx1*(x1Flux[j][i+1].Mx - x1Flux[j][i].Mx);
      pG->U[ks][j][i].M2  -= dtodx1*(x1Flux[j][i+1].My - x1Flux[j][i].My);
      pG->U[ks][j][i].M3  -= dtodx1*(x1Flux[j][i+1].Mz - x1Flux[j][i].Mz);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E   -= dtodx1*(x1Flux[j][i+1].E  - x1Flux[j][i].E );
#endif /* BAROTROPIC */
#ifdef MHD
      pG->U[ks][j][i].B3c -= dtodx1*(x1Flux[j][i+1].Bz - x1Flux[j][i].Bz);
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        pG->U[ks][j][i].s[n] -= dtodx1*(x1Flux[j][i+1].s[n]
                                      - x1Flux[j][i  ].s[n]);
#endif
#ifdef USE_ENTROPY_FIX
      S[j][i] -= dtodx1*(x1FluxS[j][i+1]  - x1FluxS[j][i]);
#endif
    }
  }

/*--- Step 14b -----------------------------------------------------------------
 * Update cell-centered variables in pG (including B3c) using 2D x2-Fluxes
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].d   -= dtodx2*(x2Flux[j+1][i].d  - x2Flux[j][i].d );
      pG->U[ks][j][i].M1  -= dtodx2*(x2Flux[j+1][i].Mz - x2Flux[j][i].Mz);
      pG->U[ks][j][i].M2  -= dtodx2*(x2Flux[j+1][i].Mx - x2Flux[j][i].Mx);
      pG->U[ks][j][i].M3  -= dtodx2*(x2Flux[j+1][i].My - x2Flux[j][i].My);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E   -= dtodx2*(x2Flux[j+1][i].E  - x2Flux[j][i].E );
#endif /* BAROTROPIC */
#ifdef MHD
      pG->U[ks][j][i].B3c -= dtodx2*(x2Flux[j+1][i].By - x2Flux[j][i].By);
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        pG->U[ks][j][i].s[n] -= dtodx2*(x2Flux[j+1][i].s[n]
                                      - x2Flux[j  ][i].s[n]);
#endif
#ifdef USE_ENTROPY_FIX
      S[j][i] -= dtodx2*(x2FluxS[j+1][i]  - x2FluxS[j][i]);
#endif
    }
  }  

#ifdef FIRST_ORDER_FLUX_CORRECTION
/*=== STEP 15: First-order flux correction ===================================*/

/*--- Step 15a -----------------------------------------------------------------
 * If cell-centered d or P have gone negative, or if v^2 > 1 in SR, correct
 * by using 1st order predictor fluxes
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      Wcheck = check_Prim(&(pG->U[ks][j][i]));
      if (Wcheck.d < 0.0) {
        flag_cell = 1;
        BadCell.i = i;
        BadCell.j = j;
        BadCell.k = ks;
        negd++;
      }
      if (Wcheck.P < 0.0) {
        flag_cell = 1;
        BadCell.i = i;
        BadCell.j = j;
        BadCell.k = ks;
        negP++;
      }
      Vsq = SQR(Wcheck.V1) + SQR(Wcheck.V2) + SQR(Wcheck.V3);
      if (Vsq > 1.0) {
        flag_cell = 1;
        BadCell.i = i;
        BadCell.j = j;
        BadCell.k = ks;
        superl++;
      }
      if (flag_cell != 0) {
        FixCell(pG, BadCell);
        flag_cell=0;
      }
    }
  }

  if (negd > 0 || negP > 0 || superl > 0){
    printf("[Step15a]: %i cells had d<0; %i cells had P<0;\n",negd,negP);
    printf("[Step15a]: %i cells had v>1	at 1st order correction\n",superl);
  }

/*--- Step 15b -----------------------------------------------------------------
 * In SR the first-order flux correction can fail to fix an unphysical state.
 * We must fix these cells in order to avoid NaN's at the next timestep,
 * particuarly if v^2 > 1. We have 2 approaches; firstly, we use the entropy
 * equation (which we have not applied a 1st order flux correction to) to
 * calculate the pressure and the Lorentz factor of the gas. If this produces
 * and unphysical state, then we floor the pressure and iterate on v^2 until
 * v^2 < 1. Possibly could improved by averaging density and pressure from
 * adjacent cells and then calculating pressure.
 */

#ifdef MHD
  negd = 0;
  negP = 0;
  superl = 0;
  entropy = 0;
  final = 0;
  fail = 0;
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      flag_cell=0;
      Wcheck = check_Prim(&(pG->U[ks][j][i]));
      if (Wcheck.d < 0.0) {
        flag_cell = 1;
        negd++;
      }
      if (Wcheck.P < 0.0) {
        flag_cell = 1;
        negP++;
      }
      Vsq = SQR(Wcheck.V1) + SQR(Wcheck.V2) + SQR(Wcheck.V3);
      if (Vsq > 1.0) {
        flag_cell = 1;
        superl++;
      }
#ifdef USE_ENTROPY_FIX
      if (flag_cell != 0) {
	entropy++;
	Wcheck = entropy_fix (&(pG->U[ks][j][i]),&(S[j][i]));
	Ucheck = Prim_to_Cons(&Wcheck);
	Wcheck = check_Prim(&Ucheck);
	Vsq = SQR(Wcheck.V1) + SQR(Wcheck.V2) + SQR(Wcheck.V3);
	if (Wcheck.d > 0.0 && Wcheck.P > 0.0 && Vsq < 1.0){
          pG->U[ks][j][i].d = Ucheck.d;
          pG->U[ks][j][i].M1 = Ucheck.M1;
          pG->U[ks][j][i].M2 = Ucheck.M2;
          pG->U[ks][j][i].M3 = Ucheck.M3;
          pG->U[ks][j][i].E = Ucheck.E;
	  flag_cell = 0;
	}
      }
#endif
      if (flag_cell != 0) {
	final++;
	Wcheck = fix_vsq (&(pG->U[ks][j][i]));
	Ucheck = Prim_to_Cons(&Wcheck);
	pG->U[ks][j][i].d = Ucheck.d;
	pG->U[ks][j][i].M1 = Ucheck.M1;
	pG->U[ks][j][i].M2 = Ucheck.M2;
	pG->U[ks][j][i].M3 = Ucheck.M3;
	pG->U[ks][j][i].E = Ucheck.E;
	Wcheck = check_Prim(&(pG->U[ks][j][i]));
	Vsq = SQR(Wcheck.V1) + SQR(Wcheck.V2) + SQR(Wcheck.V3);
	if (Wcheck.d < 0.0 || Wcheck.P < 0.0 || Vsq > 1.0){
	  fail++;
	}  
      }
    }
  }

  if (negd > 0 || negP > 0 || superl > 0) {
    printf("[Step15b]: %i cells had d<0; %i cells had P<0;\n",negd,negP);
    printf("[Step15b]: %i cells had v>1	after 1st order correction\n",superl);
    printf("[Step15b]: %i cells required an entropy fix\n",entropy);
    printf("[Step15b]: %i cells required a  final fix\n",final);
    printf("[Step15b]: %i cells had an unphysical state\n",fail);
  }
#endif /* MHD */
#endif /* FIRST_ORDER_FLUX_CORRECTION */

#ifdef STATIC_MESH_REFINEMENT
/*=== STEP 16: With SMR, store fluxes at fine/coarse boundaries ==============*/

/*--- Step 16a -----------------------------------------------------------------
 * store fluxes at boundaries of child grids.
 */

  for (ncg=0; ncg<pG->NCGrid; ncg++) {

/* x1-boundaries of child Grids (interior to THIS Grid) */

    for (dim=0; dim<2; dim++){
      if (pG->CGrid[ncg].myFlx[dim] != NULL) {

        if (dim==0) i = pG->CGrid[ncg].ijks[0];
        if (dim==1) i = pG->CGrid[ncg].ijke[0] + 1;
        jcs = pG->CGrid[ncg].ijks[1];
        jce = pG->CGrid[ncg].ijke[1];

        for (j=jcs, jj=0; j<=jce; j++, jj++){
          pG->CGrid[ncg].myFlx[dim][ks][jj].d  = x1Flux[j][i].d;
          pG->CGrid[ncg].myFlx[dim][ks][jj].M1 = x1Flux[j][i].Mx;
          pG->CGrid[ncg].myFlx[dim][ks][jj].M2 = x1Flux[j][i].My;
          pG->CGrid[ncg].myFlx[dim][ks][jj].M3 = x1Flux[j][i].Mz;
#ifndef BAROTROPIC
          pG->CGrid[ncg].myFlx[dim][ks][jj].E  = x1Flux[j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
          pG->CGrid[ncg].myFlx[dim][ks][jj].B1c = 0.0;
          pG->CGrid[ncg].myFlx[dim][ks][jj].B2c = x1Flux[j][i].By;
          pG->CGrid[ncg].myFlx[dim][ks][jj].B3c = x1Flux[j][i].Bz;
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++)
            pG->CGrid[ncg].myFlx[dim][ks][jj].s[n]  = x1Flux[j][i].s[n];
#endif
        }
#ifdef MHD
        for (j=jcs, jj=0; j<=jce+1; j++, jj++){
          pG->CGrid[ncg].myEMF3[dim][ks][jj] = emf3[j][i];
        }
#endif /* MHD */
      }
    }

/* x2-boundaries of child Grids (interior to THIS Grid) */

    for (dim=2; dim<4; dim++){
      if (pG->CGrid[ncg].myFlx[dim] != NULL) {

        ics = pG->CGrid[ncg].ijks[0];
        ice = pG->CGrid[ncg].ijke[0];
        if (dim==2) j = pG->CGrid[ncg].ijks[1];
        if (dim==3) j = pG->CGrid[ncg].ijke[1] + 1;

        for (i=ics, ii=0; i<=ice; i++, ii++){
          pG->CGrid[ncg].myFlx[dim][ks][ii].d  = x2Flux[j][i].d;
          pG->CGrid[ncg].myFlx[dim][ks][ii].M1 = x2Flux[j][i].Mz;
          pG->CGrid[ncg].myFlx[dim][ks][ii].M2 = x2Flux[j][i].Mx;
          pG->CGrid[ncg].myFlx[dim][ks][ii].M3 = x2Flux[j][i].My;
#ifndef BAROTROPIC
          pG->CGrid[ncg].myFlx[dim][ks][ii].E  = x2Flux[j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
          pG->CGrid[ncg].myFlx[dim][ks][ii].B1c = x2Flux[j][i].Bz;
          pG->CGrid[ncg].myFlx[dim][ks][ii].B2c = 0.0;
          pG->CGrid[ncg].myFlx[dim][ks][ii].B3c = x2Flux[j][i].By;
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++)
            pG->CGrid[ncg].myFlx[dim][ks][ii].s[n]  = x2Flux[j][i].s[n];
#endif
        }
#ifdef MHD
        for (i=ics, ii=0; i<=ice+1; i++, ii++){
          pG->CGrid[ncg].myEMF3[dim][ks][ii] = emf3[j][i];
        }
#endif /* MHD */
      }
    }
  }

/*--- Step 16b -----------------------------------------------------------------
 * store fluxes at boundaries with parent grids.
 */

  for (npg=0; npg<pG->NPGrid; npg++) {

/* x1-boundaries of parent Grids (at boundaries of THIS Grid)  */

    for (dim=0; dim<2; dim++){
      if (pG->PGrid[npg].myFlx[dim] != NULL) {

        if (dim==0) i = pG->PGrid[npg].ijks[0];
        if (dim==1) i = pG->PGrid[npg].ijke[0] + 1;
        jps = pG->PGrid[npg].ijks[1];
        jpe = pG->PGrid[npg].ijke[1];

        for (j=jps, jj=0; j<=jpe; j++, jj++){
          pG->PGrid[npg].myFlx[dim][ks][jj].d  = x1Flux[j][i].d;
          pG->PGrid[npg].myFlx[dim][ks][jj].M1 = x1Flux[j][i].Mx;
          pG->PGrid[npg].myFlx[dim][ks][jj].M2 = x1Flux[j][i].My;
          pG->PGrid[npg].myFlx[dim][ks][jj].M3 = x1Flux[j][i].Mz;
#ifndef BAROTROPIC
          pG->PGrid[npg].myFlx[dim][ks][jj].E  = x1Flux[j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
          pG->PGrid[npg].myFlx[dim][ks][jj].B1c = 0.0;
          pG->PGrid[npg].myFlx[dim][ks][jj].B2c = x1Flux[j][i].By;
          pG->PGrid[npg].myFlx[dim][ks][jj].B3c = x1Flux[j][i].Bz;
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++)
            pG->PGrid[npg].myFlx[dim][ks][jj].s[n]  = x1Flux[j][i].s[n];
#endif
        }
#ifdef MHD
        for (j=jps, jj=0; j<=jpe+1; j++, jj++){
          pG->PGrid[npg].myEMF3[dim][ks][jj] = emf3[j][i];
        }
#endif /* MHD */
      }
    }

/* x2-boundaries of parent Grids (at boundaries of THIS Grid)  */

    for (dim=2; dim<4; dim++){
      if (pG->PGrid[npg].myFlx[dim] != NULL) {

        ips = pG->PGrid[npg].ijks[0];
        ipe = pG->PGrid[npg].ijke[0];
        if (dim==2) j = pG->PGrid[npg].ijks[1];
        if (dim==3) j = pG->PGrid[npg].ijke[1] + 1;

        for (i=ips, ii=0; i<=ipe; i++, ii++){
          pG->PGrid[npg].myFlx[dim][ks][ii].d  = x2Flux[j][i].d;
          pG->PGrid[npg].myFlx[dim][ks][ii].M1 = x2Flux[j][i].Mz;
          pG->PGrid[npg].myFlx[dim][ks][ii].M2 = x2Flux[j][i].Mx;
          pG->PGrid[npg].myFlx[dim][ks][ii].M3 = x2Flux[j][i].My;
#ifndef BAROTROPIC
          pG->PGrid[npg].myFlx[dim][ks][ii].E  = x2Flux[j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
          pG->PGrid[npg].myFlx[dim][ks][ii].B1c = x2Flux[j][i].Bz;
          pG->PGrid[npg].myFlx[dim][ks][ii].B2c = 0.0;
          pG->PGrid[npg].myFlx[dim][ks][ii].B3c = x2Flux[j][i].By;
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++)
            pG->PGrid[npg].myFlx[dim][ks][ii].s[n]  = x2Flux[j][i].s[n];
#endif
        }
#ifdef MHD
        for (i=ips, ii=0; i<=ipe+1; i++, ii++){
          pG->PGrid[npg].myEMF3[dim][ks][ii] = emf3[j][i];
        }
#endif /* MHD */
      }
    }
  }

#endif /* STATIC_MESH_REFINEMENT */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_init_2d(MeshS *pM)
 *  \brief Allocate temporary integration arrays */
void integrate_init_2d(MeshS *pM)
{
  int nmax,size1=0,size2=0,nl,nd;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if (pM->Domain[nl][nd].Grid->Nx[0] > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
        if (pM->Domain[nl][nd].Grid->Nx[1] > size2){
          size2 = pM->Domain[nl][nd].Grid->Nx[1];
        }
      }
    }
  }

  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  nmax = MAX(size1,size2);

#ifdef MHD
  if ((emf3 = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((emf3_cc = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* MHD */
#ifdef H_CORRECTION
  if ((eta1 = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((eta2 = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* H_CORRECTION */
  if ((Bxc = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
#ifdef MHD
  if ((B1_x1Face = (Real**)calloc_2d_array(size2,size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((B2_x2Face = (Real**)calloc_2d_array(size2,size1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* MHD */

  if ((U1d= (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ul = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ur = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W1d= (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wl = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wr = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((Wl_x1Face = (Prim1DS**)calloc_2d_array(size2,size1, sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wr_x1Face = (Prim1DS**)calloc_2d_array(size2,size1, sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wl_x2Face = (Prim1DS**)calloc_2d_array(size2,size1, sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wr_x2Face = (Prim1DS**)calloc_2d_array(size2,size1, sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((x1Flux    = (Cons1DS**)calloc_2d_array(size2,size1, sizeof(Cons1DS))) 
    == NULL) goto on_error;
  if ((x2Flux    = (Cons1DS**)calloc_2d_array(size2,size1, sizeof(Cons1DS))) 
    == NULL) goto on_error;
#ifdef FIRST_ORDER_FLUX_CORRECTION
  if ((x1FluxP = (Cons1DS**)calloc_2d_array(size2,size1, sizeof(Cons1DS))) 
    == NULL) goto on_error;
  if ((x2FluxP = (Cons1DS**)calloc_2d_array(size2,size1, sizeof(Cons1DS))) 
    == NULL) goto on_error;
#ifdef MHD
  if ((emf3P = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
#endif
#endif /* FIRST_ORDER_FLUX_CORRECTION */
#ifdef USE_ENTROPY_FIX
  if ((S = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((Shalf = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((x1FluxS = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((x2FluxS = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
#ifdef FIRST_ORDER_FLUX_CORRECTION
  if ((x1FluxSP = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((x2FluxSP = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
#endif
#endif

  if ((Uhalf = (ConsS**)calloc_2d_array(size2,size1,sizeof(ConsS)))==NULL)
    goto on_error;
  if ((Whalf = (PrimS**)calloc_2d_array(size2,size1,sizeof(PrimS)))==NULL)
    goto on_error;
  if ((W     = (PrimS**)calloc_2d_array(size2,size1,sizeof(PrimS)))==NULL)
    goto on_error;

  return;

  on_error:
    integrate_destruct();
    ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_destruct_2d(void)
 *  \brief Free temporary integration arrays */
void integrate_destruct_2d(void)
{
#ifdef MHD
  if (emf3    != NULL) free_2d_array(emf3);
  if (emf3_cc != NULL) free_2d_array(emf3_cc);
#endif /* MHD */
#ifdef H_CORRECTION
  if (eta1 != NULL) free_2d_array(eta1);
  if (eta2 != NULL) free_2d_array(eta2);
#endif /* H_CORRECTION */
  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);
#ifdef MHD
  if (B1_x1Face != NULL) free_2d_array(B1_x1Face);
  if (B2_x2Face != NULL) free_2d_array(B2_x2Face);
#endif /* MHD */

  if (U1d != NULL) free(U1d);
  if (Ul  != NULL) free(Ul);
  if (Ur  != NULL) free(Ur);
  if (W1d != NULL) free(W1d);
  if (Wl  != NULL) free(Wl);
  if (Wr  != NULL) free(Wr);

  if (Wl_x1Face != NULL) free_2d_array(Wl_x1Face);
  if (Wr_x1Face != NULL) free_2d_array(Wr_x1Face);
  if (Wl_x2Face != NULL) free_2d_array(Wl_x2Face);
  if (Wr_x2Face != NULL) free_2d_array(Wr_x2Face);
  if (x1Flux    != NULL) free_2d_array(x1Flux);
  if (x2Flux    != NULL) free_2d_array(x2Flux);
#ifdef FIRST_ORDER_FLUX_CORRECTION
  if (x1FluxP   != NULL) free_2d_array(x1FluxP);
  if (x2FluxP   != NULL) free_2d_array(x2FluxP);
#ifdef MHD
  if (emf3P     != NULL) free_2d_array(emf3P);
#endif
#endif
#ifdef USE_ENTROPY_FIX
  if (S         != NULL) free_2d_array(S);
  if (Shalf     != NULL) free_2d_array(Shalf);
  if (x1FluxS   != NULL) free_2d_array(x1FluxS);
  if (x2FluxS   != NULL) free_2d_array(x2FluxS);
#ifdef FIRST_ORDER_FLUX_CORRECTION
  if (x1FluxSP  != NULL) free_2d_array(x1FluxSP);
  if (x2FluxSP  != NULL) free_2d_array(x2FluxSP);
#endif
#endif

  if (Uhalf    != NULL) free_2d_array(Uhalf);
  if (Whalf    != NULL) free_2d_array(Whalf);
  if (W        != NULL) free_2d_array(W);

  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/*! \fn static void integrate_emf3_corner(GridS *pG)
 *  \brief Integrates face centered B-fluxes to compute corner EMFs. 
 *
 *   Note:
 * - x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 * - x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 * - x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 * - x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 */ 

#ifdef MHD
static void integrate_emf3_corner(GridS *pG)
{
  int i,il,iu,j,jl,ju;
  Real emf_l1, emf_r1, emf_l2, emf_r2;

  il = pG->is-(nghost-1);   iu = pG->ie+(nghost-1);
  jl = pG->js-(nghost-1);   ju = pG->je+(nghost-1);

/* NOTE: The x1-Flux of B2 is -E3.  The x2-Flux of B1 is +E3. */
  for (j=jl; j<=ju+1; j++) {
    for (i=il; i<=iu+1; i++) {
      if (x1Flux[j-1][i].d > 0.0) {
        emf_l2 = -x1Flux[j-1][i].By
          + (x2Flux[j][i-1].Bz - emf3_cc[j-1][i-1]);
      }
      else if (x1Flux[j-1][i].d < 0.0) {
        emf_l2 = -x1Flux[j-1][i].By
          + (x2Flux[j][i].Bz - emf3_cc[j-1][i]);

      } else {
        emf_l2 = -x1Flux[j-1][i].By
          + 0.5*(x2Flux[j][i-1].Bz - emf3_cc[j-1][i-1] +
                 x2Flux[j][i  ].Bz - emf3_cc[j-1][i  ] );
      }

      if (x1Flux[j][i].d > 0.0) {
        emf_r2 = -x1Flux[j][i].By
          + (x2Flux[j][i-1].Bz - emf3_cc[j][i-1]);
      }
      else if (x1Flux[j][i].d < 0.0) {
        emf_r2 = -x1Flux[j][i].By
          + (x2Flux[j][i].Bz - emf3_cc[j][i]);

      } else {
        emf_r2 = -x1Flux[j][i].By
          + 0.5*(x2Flux[j][i-1].Bz - emf3_cc[j][i-1] +
                 x2Flux[j][i  ].Bz - emf3_cc[j][i  ] );
      }

      if (x2Flux[j][i-1].d > 0.0) {
        emf_l1 = x2Flux[j][i-1].Bz
          + (-x1Flux[j-1][i].By - emf3_cc[j-1][i-1]);
      }
      else if (x2Flux[j][i-1].d < 0.0) {
        emf_l1 = x2Flux[j][i-1].Bz
          + (-x1Flux[j][i].By - emf3_cc[j][i-1]);
      } else {
        emf_l1 = x2Flux[j][i-1].Bz
          + 0.5*(-x1Flux[j-1][i].By - emf3_cc[j-1][i-1]
                 -x1Flux[j  ][i].By - emf3_cc[j  ][i-1] );
      }

      if (x2Flux[j][i].d > 0.0) {
        emf_r1 = x2Flux[j][i].Bz
          + (-x1Flux[j-1][i].By - emf3_cc[j-1][i]);
      }
      else if (x2Flux[j][i].d < 0.0) {
        emf_r1 = x2Flux[j][i].Bz
          + (-x1Flux[j][i].By - emf3_cc[j][i]);
      } else {
        emf_r1 = x2Flux[j][i].Bz
          + 0.5*(-x1Flux[j-1][i].By - emf3_cc[j-1][i]
                 -x1Flux[j  ][i].By - emf3_cc[j  ][i] );
      }

      emf3[j][i] = 0.25*(emf_l1 + emf_r1 + emf_l2 + emf_r2);
    }
  }

  return;
}
#endif /* MHD */

#ifdef FIRST_ORDER_FLUX_CORRECTION
/*----------------------------------------------------------------------------*/
/*! \fn static void FixCell(GridS *pG, Int3Vect ix)
 *  \brief Uses first order fluxes to fix negative d,P or superluminal v
 */ 
static void FixCell(GridS *pG, Int3Vect ix)
{
  int ks=pG->ks;
  Real dtodx1=pG->dt/pG->dx1, dtodx2=pG->dt/pG->dx2;
  Cons1DS x1FD_i, x1FD_ip1, x2FD_j, x2FD_jp1;
#ifdef MHD
  int i,j;
  Real emf3D_ji, emf3D_jip1, emf3D_jp1i, emf3D_jp1ip1;
#endif
#ifdef USE_ENTROPY_FIX
  Real x1SFD_i, x1SFD_ip1, x2SFD_j, x2SFD_jp1;
#endif


/* Compute difference of predictor and corrector fluxes at cell faces */

  x1FD_i.d = x1Flux[ix.j][ix.i].d - x1FluxP[ix.j][ix.i].d;
  x2FD_j.d = x2Flux[ix.j][ix.i].d - x2FluxP[ix.j][ix.i].d;

  x1FD_ip1.d = x1Flux[ix.j][ix.i+1].d - x1FluxP[ix.j][ix.i+1].d;
  x2FD_jp1.d = x2Flux[ix.j+1][ix.i].d - x2FluxP[ix.j+1][ix.i].d;

  x1FD_i.Mx = x1Flux[ix.j][ix.i].Mx - x1FluxP[ix.j][ix.i].Mx;
  x2FD_j.Mx = x2Flux[ix.j][ix.i].Mx - x2FluxP[ix.j][ix.i].Mx;

  x1FD_ip1.Mx = x1Flux[ix.j][ix.i+1].Mx - x1FluxP[ix.j][ix.i+1].Mx;
  x2FD_jp1.Mx = x2Flux[ix.j+1][ix.i].Mx - x2FluxP[ix.j+1][ix.i].Mx;

  x1FD_i.My = x1Flux[ix.j][ix.i].My - x1FluxP[ix.j][ix.i].My;
  x2FD_j.My = x2Flux[ix.j][ix.i].My - x2FluxP[ix.j][ix.i].My;

  x1FD_ip1.My = x1Flux[ix.j][ix.i+1].My - x1FluxP[ix.j][ix.i+1].My;
  x2FD_jp1.My = x2Flux[ix.j+1][ix.i].My - x2FluxP[ix.j+1][ix.i].My;

  x1FD_i.Mz = x1Flux[ix.j][ix.i].Mz - x1FluxP[ix.j][ix.i].Mz;
  x2FD_j.Mz = x2Flux[ix.j][ix.i].Mz - x2FluxP[ix.j][ix.i].Mz;

  x1FD_ip1.Mz = x1Flux[ix.j][ix.i+1].Mz - x1FluxP[ix.j][ix.i+1].Mz;
  x2FD_jp1.Mz = x2Flux[ix.j+1][ix.i].Mz - x2FluxP[ix.j+1][ix.i].Mz;

#ifndef BAROTROPIC
  x1FD_i.E = x1Flux[ix.j][ix.i].E - x1FluxP[ix.j][ix.i].E;
  x2FD_j.E = x2Flux[ix.j][ix.i].E - x2FluxP[ix.j][ix.i].E;

  x1FD_ip1.E = x1Flux[ix.j][ix.i+1].E - x1FluxP[ix.j][ix.i+1].E;
  x2FD_jp1.E = x2Flux[ix.j+1][ix.i].E - x2FluxP[ix.j+1][ix.i].E;
#endif /* BAROTROPIC */

#ifdef MHD
  x1FD_i.Bz = x1Flux[ix.j][ix.i].Bz - x1FluxP[ix.j][ix.i].Bz;
  x2FD_j.By = x2Flux[ix.j][ix.i].By - x2FluxP[ix.j][ix.i].By;

  x1FD_ip1.Bz = x1Flux[ix.j][ix.i+1].Bz - x1FluxP[ix.j][ix.i+1].Bz;
  x2FD_jp1.By = x2Flux[ix.j+1][ix.i].By - x2FluxP[ix.j+1][ix.i].By;

  emf3D_ji     = emf3[ix.j  ][ix.i  ] - emf3P[ix.j  ][ix.i  ];
  emf3D_jip1   = emf3[ix.j  ][ix.i+1] - emf3P[ix.j  ][ix.i+1];
  emf3D_jp1i   = emf3[ix.j+1][ix.i  ] - emf3P[ix.j+1][ix.i  ];
  emf3D_jp1ip1 = emf3[ix.j+1][ix.i+1] - emf3P[ix.j+1][ix.i+1];
#endif /* MHD */

  /*#ifdef USE_ENTROPY_FIX
  x1SFD_i = x1FluxS[ix.j][ix.i] - x1FluxSP[ix.j][ix.i];
  x2SFD_j = x2FluxS[ix.j][ix.i] - x2FluxSP[ix.j][ix.i];

  x1SFD_ip1 = x1FluxS[ix.j][ix.i+1] - x1FluxSP[ix.j][ix.i+1];
  x2SFD_jp1 = x2FluxS[ix.j+1][ix.i] - x2FluxSP[ix.j+1][ix.i];
  #endif*/

/* Use flux differences to correct bad cell-centered quantities */

  pG->U[ks][ix.j][ix.i].d  += dtodx1*(x1FD_ip1.d  - x1FD_i.d );
  pG->U[ks][ix.j][ix.i].M1 += dtodx1*(x1FD_ip1.Mx - x1FD_i.Mx);
  pG->U[ks][ix.j][ix.i].M2 += dtodx1*(x1FD_ip1.My - x1FD_i.My);
  pG->U[ks][ix.j][ix.i].M3 += dtodx1*(x1FD_ip1.Mz - x1FD_i.Mz);
#ifndef BAROTROPIC
  pG->U[ks][ix.j][ix.i].E  += dtodx1*(x1FD_ip1.E  - x1FD_i.E );
#endif /* BAROTROPIC */
#ifdef MHD
  pG->U[ks][ix.j][ix.i].B3c += dtodx1*(x1FD_ip1.Bz - x1FD_i.Bz);
#endif /* MHD */
  /*#ifdef USE_ENTROPY_FIX
  S[ix.j][ix.i]  += dtodx1*(x1SFD_ip1  - x1SFD_i );
  #endif*/

  pG->U[ks][ix.j][ix.i].d  += dtodx2*(x2FD_jp1.d  - x2FD_j.d );
  pG->U[ks][ix.j][ix.i].M1 += dtodx2*(x2FD_jp1.Mz - x2FD_j.Mz);
  pG->U[ks][ix.j][ix.i].M2 += dtodx2*(x2FD_jp1.Mx - x2FD_j.Mx);
  pG->U[ks][ix.j][ix.i].M3 += dtodx2*(x2FD_jp1.My - x2FD_j.My);
#ifndef BAROTROPIC
  pG->U[ks][ix.j][ix.i].E  += dtodx2*(x2FD_jp1.E  - x2FD_j.E );
#endif /* BAROTROPIC */
#ifdef MHD
  pG->U[ks][ix.j][ix.i].B3c += dtodx2*(x2FD_jp1.By - x2FD_j.By);
#endif /* MHD */
  /*#ifdef USE_ENTROPY_FIX
  S[ix.j][ix.i]  += dtodx2*(x2SFD_jp1  - x2SFD_j );
  #endif*/

/* Use flux differences to correct cell-centered values at i-1 and i+1 */

  if (ix.i > pG->is) {
    pG->U[ks][ix.j][ix.i-1].d  += dtodx1*(x1FD_i.d );
    pG->U[ks][ix.j][ix.i-1].M1 += dtodx1*(x1FD_i.Mx);
    pG->U[ks][ix.j][ix.i-1].M2 += dtodx1*(x1FD_i.My);
    pG->U[ks][ix.j][ix.i-1].M3 += dtodx1*(x1FD_i.Mz);
#ifndef BAROTROPIC
    pG->U[ks][ix.j][ix.i-1].E  += dtodx1*(x1FD_i.E );
#endif /* BAROTROPIC */
#ifdef MHD
    pG->U[ks][ix.j][ix.i-1].B3c += dtodx1*(x1FD_i.Bz);
#endif /* MHD */
    /*#ifdef USE_ENTROPY_FIX
    S[ix.j][ix.i-1]  += dtodx1*(x1SFD_i );
    #endif*/
  }

  if (ix.i < pG->ie) {
    pG->U[ks][ix.j][ix.i+1].d  -= dtodx1*(x1FD_ip1.d );
    pG->U[ks][ix.j][ix.i+1].M1 -= dtodx1*(x1FD_ip1.Mx);
    pG->U[ks][ix.j][ix.i+1].M2 -= dtodx1*(x1FD_ip1.My);
    pG->U[ks][ix.j][ix.i+1].M3 -= dtodx1*(x1FD_ip1.Mz);
#ifndef BAROTROPIC
    pG->U[ks][ix.j][ix.i+1].E  -= dtodx1*(x1FD_ip1.E );
#endif /* BAROTROPIC */
#ifdef MHD
    pG->U[ks][ix.j][ix.i+1].B3c -= dtodx1*(x1FD_ip1.Bz);
#endif /* MHD */
    /*#ifdef USE_ENTROPY_FIX
    S[ix.j][ix.i+1]  -= dtodx1*(x1SFD_ip1 );
    #endif*/
  }

/* Use flux differences to correct cell-centered values at j-1 and j+1 */

  if (ix.j > pG->js) {
    pG->U[ks][ix.j-1][ix.i].d  += dtodx2*(x2FD_j.d );
    pG->U[ks][ix.j-1][ix.i].M1 += dtodx2*(x2FD_j.Mz);
    pG->U[ks][ix.j-1][ix.i].M2 += dtodx2*(x2FD_j.Mx);
    pG->U[ks][ix.j-1][ix.i].M3 += dtodx2*(x2FD_j.My);
#ifndef BAROTROPIC
    pG->U[ks][ix.j-1][ix.i].E  += dtodx2*(x2FD_j.E );
#endif /* BAROTROPIC */
#ifdef MHD
    pG->U[ks][ix.j-1][ix.i].B3c += dtodx2*(x2FD_j.By);
#endif /* MHD */
    /*#ifdef USE_ENTROPY_FIX
   S[ix.j-1][ix.i]  += dtodx2*(x2SFD_j );
   #endif*/
  }

  if (ix.j < pG->je) {
    pG->U[ks][ix.j+1][ix.i].d  -= dtodx2*(x2FD_jp1.d );
    pG->U[ks][ix.j+1][ix.i].M1 -= dtodx2*(x2FD_jp1.Mz);
    pG->U[ks][ix.j+1][ix.i].M2 -= dtodx2*(x2FD_jp1.Mx);
    pG->U[ks][ix.j+1][ix.i].M3 -= dtodx2*(x2FD_jp1.My);
#ifndef BAROTROPIC
    pG->U[ks][ix.j+1][ix.i].E  -= dtodx2*(x2FD_jp1.E );
#endif /* BAROTROPIC */
#ifdef MHD
    pG->U[ks][ix.j+1][ix.i].B3c -= dtodx2*(x2FD_jp1.By);
#endif /* MHD */
    /*#ifdef USE_ENTROPY_FIX
    S[ix.j+1][ix.i]  -= dtodx2*(x2SFD_jp1 );
    #endif*/
  }

/* Use emf3 difference to correct face-centered fields on edges of bad cell */
#ifdef MHD
  pG->B1i[ks][ix.j][ix.i  ] += dtodx2*(emf3D_jp1i   - emf3D_ji);
  pG->B1i[ks][ix.j][ix.i+1] += dtodx2*(emf3D_jp1ip1 - emf3D_jip1);
  pG->B2i[ks][ix.j  ][ix.i] -= dtodx1*(emf3D_jip1   - emf3D_ji);
  pG->B2i[ks][ix.j+1][ix.i] -= dtodx1*(emf3D_jp1ip1 - emf3D_jp1i);

/* Use emf3 difference to correct face-centered fields around bad cell */
  if (ix.j > pG->js) {
    pG->B1i[ks][ix.j-1][ix.i  ] -= dtodx2*(emf3D_ji);
    pG->B1i[ks][ix.j-1][ix.i+1] -= dtodx2*(emf3D_jip1);
  }
  if (ix.j < pG->je) {
    pG->B1i[ks][ix.j+1][ix.i  ] += dtodx2*(emf3D_jp1i);
    pG->B1i[ks][ix.j+1][ix.i+1] += dtodx2*(emf3D_jp1ip1);
  }

  if (ix.i > pG->is) {
    pG->B2i[ks][ix.j  ][ix.i-1] += dtodx1*(emf3D_ji);
    pG->B2i[ks][ix.j+1][ix.i-1] += dtodx1*(emf3D_jp1i);
  }
  if (ix.i < pG->ie) {
    pG->B2i[ks][ix.j  ][ix.i+1] -= dtodx1*(emf3D_jip1);
    pG->B2i[ks][ix.j+1][ix.i+1] -= dtodx1*(emf3D_jp1ip1);
  }

/* compute new cell-centered fields */
  for (j=(ix.j-1); j<=(ix.j+1); j++) {
  for (i=(ix.i-1); i<=(ix.i+1); i++) {
    pG->U[ks][j][i].B1c = 0.5*(pG->B1i[ks][j][i] + pG->B1i[ks][j][i+1]);
    pG->U[ks][j][i].B2c = 0.5*(pG->B2i[ks][j][i] + pG->B2i[ks][j+1][i]);
  }}
#endif /* MHD */

#ifdef STATIC_MESH_REFINEMENT
/* With SMR, replace higher-order fluxes with predict fluxes in case they are
 * used at fine/coarse grid boundaries */ 

  x1Flux[ix.j][ix.i] = x1FluxP[ix.j][ix.i];
  x2Flux[ix.j][ix.i] = x2FluxP[ix.j][ix.i];
#endif /* STATIC_MESH_REFINEMENT */

}
#endif /* FIRST_ORDER_FLUX_CORRECTION */

#endif /* VL_INTEGRATOR */

#endif /* SPECIAL_RELATIVITY */
