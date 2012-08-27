#include "../copyright.h"
/*============================================================================*/
/*! \file integrate_3d_vl_sr.c
 *  \brief Integrate MHD equations using 3D version of the directionally
 *   unsplit MUSCL-Hancock (VL) integrator. 
 *
 * PURPOSE: Integrate MHD equations using 3D version of the directionally
 *   unsplit MUSCL-Hancock (VL) integrator.  The variables updated are:
 *    - U.[d,M1,M2,M3,E,B1c,B2c,B3c,s] -- where U is of type ConsS
 *    - B1i, B2i, B3i  -- interface magnetic field
 *   Also adds gravitational source terms, self-gravity, and the H-correction
 *   of Sanders et al.
 *   - For adb hydro, requires (9*Cons1DS + 3*Real + 1*ConsS) = 53 3D arrays
 *   - For adb mhd, requires   (9*Cons1DS + 9*Real + 1*ConsS) = 80 3D arrays
 *
 * REFERENCE: 
 * - J.M Stone & T.A. Gardiner, "A simple, unsplit Godunov method
 *   for multidimensional MHD", NewA 14, 139 (2009)
 *
 * - R. Sanders, E. Morano, & M.-C. Druguet, "Multidimensinal dissipation for
 *   upwind schemes: stability and applications to gas dynamics", JCP, 145, 511
 *   (1998)
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - integrate_3d_vl()
 * - integrate_destruct_3d()
 * - integrate_init_3d() */
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
static Prim1DS ***Wl_x1Face=NULL, ***Wr_x1Face=NULL;
static Prim1DS ***Wl_x2Face=NULL, ***Wr_x2Face=NULL;
static Prim1DS ***Wl_x3Face=NULL, ***Wr_x3Face=NULL;
static Cons1DS ***x1Flux=NULL, ***x2Flux=NULL, ***x3Flux=NULL;
#ifdef FIRST_ORDER_FLUX_CORRECTION
static Cons1DS ***x1FluxP=NULL, ***x2FluxP=NULL, ***x3FluxP=NULL;
static Real ***emf1P=NULL, ***emf2P=NULL, ***emf3P=NULL;
#endif
#ifdef USE_ENTROPY_FIX
static Real ***S=NULL,***Shalf=NULL;
static Real ***x1FluxS=NULL,***x2FluxS=NULL,***x3FluxS=NULL;
#ifdef FIRST_ORDER_FLUX_CORRECTION
static Real ***x1FluxSP=NULL,***x2FluxSP=NULL,***x3FluxSP=NULL;
#endif
#endif


/* The interface magnetic fields and emfs */
#ifdef MHD
static Real ***B1_x1Face=NULL, ***B2_x2Face=NULL, ***B3_x3Face=NULL;
static Real ***emf1=NULL, ***emf2=NULL, ***emf3=NULL;
static Real ***emf1_cc=NULL, ***emf2_cc=NULL, ***emf3_cc=NULL;
#endif /* MHD */

/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W1d=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL, *Ul=NULL, *Ur=NULL;

/* primitive variables at t^{n} computed in predict step */
static PrimS ***W=NULL;

/* conserved & primitive variables at t^{n+1/2} computed in predict step */
static ConsS ***Uhalf=NULL;
static PrimS ***Whalf=NULL;

/* variables needed for H-correction of Sanders et al (1998) */
extern Real etah;
#ifdef H_CORRECTION
static Real ***eta1=NULL, ***eta2=NULL, ***eta3=NULL;
#endif


/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   integrate_emf1_corner() - upwind CT method of GS (2005) for emf1
 *   integrate_emf2_corner() - upwind CT method of GS (2005) for emf2 
 *   integrate_emf3_corner() - upwind CT method of GS (2005) for emf3
 *   FixCell() - apply first-order correction to one cell
 *============================================================================*/
#ifdef MHD
static void integrate_emf1_corner(const GridS *pG);
static void integrate_emf2_corner(const GridS *pG);
static void integrate_emf3_corner(const GridS *pG);
#endif
#ifdef FIRST_ORDER_FLUX_CORRECTION
static void FixCell(GridS *pG, Int3Vect);
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void integrate_3d_vl(DomainS *pD) 
 *  \brief 3D van Leer unsplit integrator for MHD. 
 */
void integrate_3d_vl(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  Real dtodx1=pG->dt/pG->dx1, dtodx2=pG->dt/pG->dx2, dtodx3=pG->dt/pG->dx3;
  Real q1 = 0.5*dtodx1, q2 = 0.5*dtodx2, q3 = 0.5*dtodx3;
  Real dt = pG->dt, hdt = 0.5*pG->dt;
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int cart_x1 = 1, cart_x2 = 2, cart_x3 = 3;
  Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic,Bx;
#if (NSCALARS > 0)
  int n;
#endif
#ifdef SELF_GRAVITY
  Real gxl,gxr,gyl,gyr,gzl,gzr,flx_m1l,flx_m1r,flx_m2l,flx_m2r,flx_m3l,flx_m3r;
#endif
#ifdef H_CORRECTION
  Real cfr,cfl,lambdar,lambdal;
#endif
#ifdef STATIC_MESH_REFINEMENT
  int ncg,npg,dim;
  int ii,ics,ice,jj,jcs,jce,kk,kcs,kce,ips,ipe,jps,jpe,kps,kpe;
#endif
#ifdef FIRST_ORDER_FLUX_CORRECTION
  int flag_cell=0,negd=0,negP=0,superl=0,NaNFlux=0;
  int entropy,final,fail;
  Real Vsq;
  ConsS Ucheck;
  PrimS Wcheck;
  Int3Vect BadCell;
#endif
  int il=is-(nghost-1), iu=ie+(nghost-1);
  int jl=js-(nghost-1), ju=je+(nghost-1);
  int kl=ks-(nghost-1), ku=ke+(nghost-1);

/* Set etah=0 so first calls to flux functions do not use H-correction */
  etah = 0.0;

  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        Uhalf[k][j][i] = pG->U[k][j][i];
        W[k][j][i] = Cons_to_Prim(&(pG->U[k][j][i]));
#ifdef USE_ENTROPY_FIX
	S[k][j][i] = W[k][j][i].P * pow(W[k][j][i].d,1.0-Gamma);
	S[k][j][i]*= pG->U[k][j][i].d / W[k][j][i].d;
	Shalf[k][j][i] = S[k][j][i];
#endif
#ifdef MHD
        B1_x1Face[k][j][i] = pG->B1i[k][j][i];
        B2_x2Face[k][j][i] = pG->B2i[k][j][i];
        B3_x3Face[k][j][i] = pG->B3i[k][j][i];
#endif /* MHD */
      }
    }
  }

/*=== STEP 1: Compute first-order fluxes at t^{n} in x1-direction ============*/
/* No source terms are needed since there is no temporal evolution */

/*--- Step 1a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;
 * W1d = (d, V1, V2, V3, P, B2c, B3c, s[n])
 */

  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
	W1d[i].d  = W[k][j][i].d;
	W1d[i].Vx = W[k][j][i].V1;
	W1d[i].Vy = W[k][j][i].V2;
	W1d[i].Vz = W[k][j][i].V3;
#ifndef BAROTROPIC
	W1d[i].P  = W[k][j][i].P;
#endif /* BAROTROPIC */
#ifdef MHD
	W1d[i].By = W[k][j][i].B2c;
	W1d[i].Bz = W[k][j][i].B3c;
        Bxc[i] = W[k][j][i].B1c;
        Bxi[i] = pG->B1i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) W1d[i].s[n] = W[k][j][i].s[n];
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
 * No source terms needed.
 */

/*--- Step 1d ------------------------------------------------------------------
 * Compute flux in x1-direction */

      for (i=il; i<=ie+nghost; i++) {
        fluxes(Ul[i],Ur[i],Wl[i],Wr[i],Bxi[i],&x1Flux[k][j][i]);
#ifdef USE_ENTROPY_FIX
	entropy_flux(Ul[i],Ur[i],Wl[i],Wr[i],Bxi[i],&x1FluxS[k][j][i]);
#endif
      }
    }
  }

/*=== STEP 2: Compute first-order fluxes at t^{n} in x2-direction ============*/
/* No source terms are needed since there is no temporal evolution */

/*--- Step 2a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;
 * W1d = (d, V2, V3, V1, P, B3c, B1c, s[n])
 */

  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      for (j=js-nghost; j<=je+nghost; j++) {
	W1d[j].d  = W[k][j][i].d;
	W1d[j].Vx = W[k][j][i].V2;
	W1d[j].Vy = W[k][j][i].V3;
	W1d[j].Vz = W[k][j][i].V1;
#ifndef BAROTROPIC
	W1d[j].P  = W[k][j][i].P;
#endif /* BAROTROPIC */
#ifdef MHD
	W1d[j].By = W[k][j][i].B3c;
	W1d[j].Bz = W[k][j][i].B1c;
        Bxc[j] = W[k][j][i].B2c;
        Bxi[j] = pG->B2i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) W1d[j].s[n] = W[k][j][i].s[n];
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
 * No source terms needed
 */

/*--- Step 2d ------------------------------------------------------------------
 * Compute flux in x2-direction */

      for (j=jl; j<=je+nghost; j++) {
        fluxes(Ul[j],Ur[j],Wl[j],Wr[j],Bxi[j],&x2Flux[k][j][i]);
#ifdef USE_ENTROPY_FIX
	entropy_flux(Ul[j],Ur[j],Wl[j],Wr[j],Bxi[j],&x2FluxS[k][j][i]);
#endif
      }
    }
  }


/*=== STEP 3: Compute first-order fluxes at t^{n} in x3-direction ============*/
/* No source terms are needed since there is no temporal evolution */

/*--- Step 3a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;
 * W1d = (d, V3, V1, V2, P, B1c, B2c, s[n])
 */

  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      for (k=ks-nghost; k<=ke+nghost; k++) {
	W1d[k].d  = W[k][j][i].d;
	W1d[k].Vx = W[k][j][i].V3;
	W1d[k].Vy = W[k][j][i].V1;
	W1d[k].Vz = W[k][j][i].V2;
#ifndef BAROTROPIC
	W1d[k].P  = W[k][j][i].P;
#endif /* BAROTROPIC */
#ifdef MHD
	W1d[k].By = W[k][j][i].B1c;
	W1d[k].Bz = W[k][j][i].B2c;
        Bxc[k] = W[k][j][i].B3c;
        Bxi[k] = pG->B3i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) W1d[k].s[n] = W[k][j][i].s[n];
#endif
      }

/*--- Step 3b ------------------------------------------------------------------
 * Compute first-order L/R states */      
      for (k=kl; k<=ke+nghost; k++) { 
        Wl[k] = W1d[k-1];
        Wr[k] = W1d[k  ]; 

/* Compute U from W in case Pfloor used in Cons1D_to_Prim1D */
        Ul[k] = Prim1D_to_Cons1D(&Wl[k], &Bxc[k-1]);
        Ur[k] = Prim1D_to_Cons1D(&Wr[k], &Bxc[k  ]);
      }

/*--- Step 3c ------------------------------------------------------------------
 * No source terms needed. 
 */

/*--- Step 3d ------------------------------------------------------------------
 * Compute flux in x1-direction */

      for (k=kl; k<=ke+nghost; k++) {
        fluxes(Ul[k],Ur[k],Wl[k],Wr[k],Bxi[k],&x3Flux[k][j][i]);
#ifdef USE_ENTROPY_FIX
	entropy_flux(Ul[k],Ur[k],Wl[k],Wr[k],Bxi[k],&x3FluxS[k][j][i]);
#endif
      }
    }
  }

/*=== STEP 4:  Update face-centered B for 0.5*dt =============================*/

/*--- Step 4a ------------------------------------------------------------------
 * Calculate the cell centered value of emf1,2,3 at t^{n} and integrate
 * to corner.
 */

#ifdef MHD
  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        emf1_cc[k][j][i] = (W[k][j][i].B2c*W[k][j][i].V3 - 
			    W[k][j][i].B3c*W[k][j][i].V2);
        emf2_cc[k][j][i] = (W[k][j][i].B3c*W[k][j][i].V1 -
			    W[k][j][i].B1c*W[k][j][i].V3);
        emf3_cc[k][j][i] = (W[k][j][i].B1c*W[k][j][i].V2 -
			    W[k][j][i].B2c*W[k][j][i].V1);
      }
    }
  }
  integrate_emf1_corner(pG);
  integrate_emf2_corner(pG);
  integrate_emf3_corner(pG);

/*--- Step 4b ------------------------------------------------------------------
 * Update the interface magnetic fields using CT for a half time step.
 */

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        B1_x1Face[k][j][i] += q3*(emf2[k+1][j  ][i  ] - emf2[k][j][i]) -
                              q2*(emf3[k  ][j+1][i  ] - emf3[k][j][i]);
        B2_x2Face[k][j][i] += q1*(emf3[k  ][j  ][i+1] - emf3[k][j][i]) -
                              q3*(emf1[k+1][j  ][i  ] - emf1[k][j][i]);
        B3_x3Face[k][j][i] += q2*(emf1[k  ][j+1][i  ] - emf1[k][j][i]) -
                              q1*(emf2[k  ][j  ][i+1] - emf2[k][j][i]);
      }
      B1_x1Face[k][j][iu+1] += q3*(emf2[k+1][j  ][iu+1]-emf2[k][j][iu+1]) -
                               q2*(emf3[k  ][j+1][iu+1]-emf3[k][j][iu+1]);
    }
    for (i=il; i<=iu; i++) {
      B2_x2Face[k][ju+1][i] += q1*(emf3[k  ][ju+1][i+1]-emf3[k][ju+1][i]) -
                               q3*(emf1[k+1][ju+1][i  ]-emf1[k][ju+1][i]);
    }
  }
  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      B3_x3Face[ku+1][j][i] += q2*(emf1[ku+1][j+1][i  ]-emf1[ku+1][j][i]) -
                               q1*(emf2[ku+1][j  ][i+1]-emf2[ku+1][j][i]);
    }
  }

/*--- Step 4c ------------------------------------------------------------------
 * Compute cell-centered magnetic fields at half-timestep from average of
 * face-centered fields.
 */

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        Uhalf[k][j][i].B1c = 0.5*(B1_x1Face[k][j][i] + B1_x1Face[k][j][i+1]);
        Uhalf[k][j][i].B2c = 0.5*(B2_x2Face[k][j][i] + B2_x2Face[k][j+1][i]);
        Uhalf[k][j][i].B3c = 0.5*(B3_x3Face[k][j][i] + B3_x3Face[k+1][j][i]);
      }
    }
  }
#endif /* MHD */

/*=== STEP 5: Update cell-centered variables to half-timestep ================*/

/*--- Step 5a ------------------------------------------------------------------
 * Update cell-centered variables to half-timestep using x1-fluxes
 */

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        Uhalf[k][j][i].d   -= q1*(x1Flux[k][j][i+1].d  - x1Flux[k][j][i].d );
        Uhalf[k][j][i].M1  -= q1*(x1Flux[k][j][i+1].Mx - x1Flux[k][j][i].Mx);
        Uhalf[k][j][i].M2  -= q1*(x1Flux[k][j][i+1].My - x1Flux[k][j][i].My);
        Uhalf[k][j][i].M3  -= q1*(x1Flux[k][j][i+1].Mz - x1Flux[k][j][i].Mz);
#ifndef BAROTROPIC
        Uhalf[k][j][i].E   -= q1*(x1Flux[k][j][i+1].E  - x1Flux[k][j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          Uhalf[k][j][i].s[n] -= 
            q1*(x1Flux[k][j][i+1].s[n] - x1Flux[k][j][i  ].s[n]);
#endif
#ifdef USE_ENTROPY_FIX
      Shalf[k][j][i]   -= q1*(x1FluxS[k][j][i+1]  - x1FluxS[k][j][i]);
#endif
      }
    }
  }

/*--- Step 5b ------------------------------------------------------------------
 * Update cell-centered variables to half-timestep using x2-fluxes
 */

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        Uhalf[k][j][i].d   -= q2*(x2Flux[k][j+1][i].d  - x2Flux[k][j][i].d );
        Uhalf[k][j][i].M1  -= q2*(x2Flux[k][j+1][i].Mz - x2Flux[k][j][i].Mz);
        Uhalf[k][j][i].M2  -= q2*(x2Flux[k][j+1][i].Mx - x2Flux[k][j][i].Mx);
        Uhalf[k][j][i].M3  -= q2*(x2Flux[k][j+1][i].My - x2Flux[k][j][i].My);
#ifndef BAROTROPIC
        Uhalf[k][j][i].E   -= q2*(x2Flux[k][j+1][i].E  - x2Flux[k][j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          Uhalf[k][j][i].s[n] -= 
            q2*(x2Flux[k][j+1][i].s[n] - x2Flux[k][j  ][i].s[n]);
#endif
#ifdef USE_ENTROPY_FIX
      Shalf[k][j][i]   -= q2*(x2FluxS[k][j+1][i]  - x2FluxS[k][j][i] );
#endif
      }
    }
  }

/*--- Step 5c ------------------------------------------------------------------
 * Update cell-centered variables to half-timestep using x3-fluxes
 */

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        Uhalf[k][j][i].d   -= q3*(x3Flux[k+1][j][i].d  - x3Flux[k][j][i].d );
        Uhalf[k][j][i].M1  -= q3*(x3Flux[k+1][j][i].My - x3Flux[k][j][i].My);
        Uhalf[k][j][i].M2  -= q3*(x3Flux[k+1][j][i].Mz - x3Flux[k][j][i].Mz);
        Uhalf[k][j][i].M3  -= q3*(x3Flux[k+1][j][i].Mx - x3Flux[k][j][i].Mx);
#ifndef BAROTROPIC
        Uhalf[k][j][i].E   -= q3*(x3Flux[k+1][j][i].E  - x3Flux[k][j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          Uhalf[k][j][i].s[n] -=
            q3*(x3Flux[k+1][j][i].s[n] - x3Flux[k  ][j][i].s[n]);
#endif
#ifdef USE_ENTROPY_FIX
      Shalf[k][j][i]   -= q3*(x3FluxS[k+1][j][i]  - x3FluxS[k][j][i] );
#endif
      }
    }
  }

#ifdef FIRST_ORDER_FLUX_CORRECTION
/*--- Step 5d ------------------------------------------------------------------
 * With first-order flux correction, save predict fluxes and emf3
 */

  NaNFlux = 0;
  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        x1FluxP[k][j][i] = x1Flux[k][j][i];
        x2FluxP[k][j][i] = x2Flux[k][j][i];
        x3FluxP[k][j][i] = x3Flux[k][j][i];
#ifdef MHD
        emf1P[k][j][i] = emf1[k][j][i];
        emf2P[k][j][i] = emf2[k][j][i];
        emf3P[k][j][i] = emf3[k][j][i];
#endif
#ifdef USE_ENTROPY_FIX
      x1FluxSP[k][j][i] = x1FluxS[k][j][i];
      x2FluxSP[k][j][i] = x2FluxS[k][j][i];
      x3FluxSP[k][j][i] = x3FluxS[k][j][i];
#endif
        if ((x1Flux[k][j][i].d  != x1Flux[k][j][i].d)  ||
#ifndef BAROTROPIC
            (x1Flux[k][j][i].E  != x1Flux[k][j][i].E)  ||
#endif
#ifdef MHD
            (x1Flux[k][j][i].By != x1Flux[k][j][i].By) ||
            (x1Flux[k][j][i].Bz != x1Flux[k][j][i].Bz) ||
#endif
            (x1Flux[k][j][i].Mx != x1Flux[k][j][i].Mx) ||
            (x1Flux[k][j][i].My != x1Flux[k][j][i].My) ||
            (x1Flux[k][j][i].Mz != x1Flux[k][j][i].Mz)) {
          x1Flux[k][j][i] = x1FluxP[k][j][i];
          NaNFlux++;
        }
        if ((x2Flux[k][j][i].d  != x2Flux[k][j][i].d)  ||
#ifndef BAROTROPIC
            (x2Flux[k][j][i].E  != x2Flux[k][j][i].E)  ||
#endif
#ifdef MHD
            (x2Flux[k][j][i].By != x2Flux[k][j][i].By) ||
            (x2Flux[k][j][i].Bz != x2Flux[k][j][i].Bz) ||
#endif
            (x2Flux[k][j][i].Mx != x2Flux[k][j][i].Mx) ||
            (x2Flux[k][j][i].My != x2Flux[k][j][i].My) ||
            (x2Flux[k][j][i].Mz != x2Flux[k][j][i].Mz)) {
          x2Flux[k][j][i] = x2FluxP[k][j][i];
          NaNFlux++;
        }
        if ((x3Flux[k][j][i].d  != x3Flux[k][j][i].d)  ||
#ifndef BAROTROPIC
            (x3Flux[k][j][i].E  != x3Flux[k][j][i].E)  ||
#endif
#ifdef MHD
            (x3Flux[k][j][i].By != x3Flux[k][j][i].By) ||
            (x3Flux[k][j][i].Bz != x3Flux[k][j][i].Bz) ||
#endif
            (x3Flux[k][j][i].Mx != x3Flux[k][j][i].Mx) ||
            (x3Flux[k][j][i].My != x3Flux[k][j][i].My) ||
            (x3Flux[k][j][i].Mz != x3Flux[k][j][i].Mz)) {
          x3Flux[k][j][i] = x3FluxP[k][j][i];
          NaNFlux++;
        }

      }
    }
  }
  if (NaNFlux != 0) {
    ath_error("[Step5d] %i first-order fluxes are NaN!\n",NaNFlux);
  }
  NaNFlux = 0;
#endif

/*=== STEP 6: Add source terms to predict values at half-timestep ============*/

/*--- Step 6a ------------------------------------------------------------------
 * Add source terms from a static gravitational potential for 0.5*dt to predict
 * step.  To improve conservation of total energy, we average the energy
 * source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          cc_pos(pG,i,j,k,&x1,&x2,&x3);
          phic = (*StaticGravPot)(x1,x2,x3);
          phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
          phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

          Uhalf[k][j][i].M1 -= q1*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
          Uhalf[k][j][i].E -= q1*(x1Flux[k][j][i  ].d*(phic - phil)
                                + x1Flux[k][j][i+1].d*(phir - phic));
#endif
          phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
          phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

          Uhalf[k][j][i].M2 -= q2*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
          Uhalf[k][j][i].E -= q2*(x2Flux[k][j  ][i].d*(phic - phil)
                                + x2Flux[k][j+1][i].d*(phir - phic));
#endif
          phir = (*StaticGravPot)(x1,x2,(x3+0.5*pG->dx3));
          phil = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));

          Uhalf[k][j][i].M3 -= q3*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
          Uhalf[k][j][i].E -= q3*(x3Flux[k  ][j][i].d*(phic - phil)
                                + x3Flux[k+1][j][i].d*(phir - phic));
#endif
        }
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
  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      for (j=js-nghost; j<=je+nghost; j++) {
	Whalf[k][j][i] = check_Prim(&(Uhalf[k][j][i]));
#ifdef FIRST_ORDER_FLUX_CORRECTION
	if (Whalf[k][j][i].d < 0.0) {
	  flag_cell = 1;
	  BadCell.i = i;
	  BadCell.j = j;
	  BadCell.k = ks;
	  negd++;
	}
	if (Whalf[k][j][i].P < 0.0) {
	  flag_cell = 1;
	  BadCell.i = i;
	  BadCell.j = j;
	  BadCell.k = ks;
	  negP++;
	}
	Vsq = SQR(Whalf[k][j][i].V1) +
	      SQR(Whalf[k][j][i].V2) + 
	      SQR(Whalf[k][j][i].V3);
	if (Vsq > 1.0) {
	  flag_cell = 1;
	  BadCell.i = i;
	  BadCell.j = j;
	  BadCell.k = ks;
	  superl++;
	}
	if (flag_cell != 0) {
#ifdef USE_ENTROPY_FIX
	Wcheck = entropy_fix (&(Uhalf[k][j][i]),&(Shalf[k][j][i]));
	Vsq = SQR(Wcheck.V1) + SQR(Wcheck.V2) + SQR(Wcheck.V3);
	if (Wcheck.d > 0.0 && Wcheck.P > 0.0 && Vsq < 1.0){
	  entropy++;
	  Whalf[k][j][i].d = Wcheck.d;
	  Whalf[k][j][i].P = Wcheck.P;
	  Whalf[k][j][i].V1 = Wcheck.V1;
	  Whalf[k][j][i].V2 = Wcheck.V2;
	  Whalf[k][j][i].V3 = Wcheck.V3;
	  flag_cell=0;
	} else {
#endif /* USE_ENTROPY_FIX */
	  Whalf[k][j][i].d = W[k][j][i].d;
	  Whalf[k][j][i].V1 = W[k][j][i].V1;
	  Whalf[k][j][i].V2 = W[k][j][i].V2;
	  Whalf[k][j][i].V3 = W[k][j][i].V3;
	  Whalf[k][j][i].P = W[k][j][i].P;
	  flag_cell=0;
#ifdef USE_ENTROPY_FIX
	}
#endif /* USE_ENTROPY_FIX */
	}
#endif
      }
    }
  }

#ifdef FIRST_ORDER_FLUX_CORRECTION
  if (negd > 0 || negP > 0 || superl > 0){
    printf("[Step7]: %i cells had d<0; %i cells had P<0;\n",negd,negP); 
    printf("[Step7]: %i cells had v>1 at t_half\n",superl);
#ifdef USE_ENTROPY_FIX
    printf("[Step7]: %i cells fixed using entropy at t_half\n",entropy);
#endif /* USE_ENTROPY_FIX */
  }
#endif 

/*=== STEP 8: Compute second-order L/R x1-interface states ===================*/

/*--- Step 8a ------------------------------------------------------------------
 * Load 1D vector of primitve variables;
 * W1d = (d, V1, V2, V3, P, B2c, B3c, s[n])
 */

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=il; i<=iu; i++) {
        W1d[i].d  = Whalf[k][j][i].d;
        W1d[i].Vx = Whalf[k][j][i].V1;
        W1d[i].Vy = Whalf[k][j][i].V2;
        W1d[i].Vz = Whalf[k][j][i].V3;
#ifndef BAROTROPIC
        W1d[i].P  = Whalf[k][j][i].P;
#endif /* BAROTROPIC */
#ifdef MHD
        W1d[i].By = Whalf[k][j][i].B2c;
        W1d[i].Bz = Whalf[k][j][i].B3c;
        Bxc[i] = Whalf[k][j][i].B1c;
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) W1d[i].s[n] = Whalf[k][j][i].s[n];
#endif
      }

/*--- Step 8b ------------------------------------------------------------------
 * Compute L and R states at x1-interfaces, store in 3D array
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
        Wl_x1Face[k][j][i] = Wl[i];
        Wr_x1Face[k][j][i] = Wr[i];
      }
    }
  }

/*=== STEP 9: Compute second-order L/R x2-interface states ===================*/

/*--- Step 9a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * W1d = (d, V2, V3, V1, P, B3c, B1c, s[n])
 */

  for (k=ks-1; k<=ke+1; k++) {
    for (i=is-1; i<=ie+1; i++) {
      for (j=jl; j<=ju; j++) {
        W1d[j].d  = Whalf[k][j][i].d;
        W1d[j].Vx = Whalf[k][j][i].V2;
        W1d[j].Vy = Whalf[k][j][i].V3;
        W1d[j].Vz = Whalf[k][j][i].V1;
#ifndef BAROTROPIC
        W1d[j].P  = Whalf[k][j][i].P;
#endif /* BAROTROPIC */
#ifdef MHD
        W1d[j].By = Whalf[k][j][i].B3c;
        W1d[j].Bz = Whalf[k][j][i].B1c;
        Bxc[j] = Whalf[k][j][i].B2c;
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) W1d[j].s[n] = Whalf[k][j][i].s[n];
#endif
      }

/*--- Step 9b ------------------------------------------------------------------
 * Compute L and R states at x2-interfaces, store in 3D array
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
        Wl_x2Face[k][j][i] = Wl[j];
        Wr_x2Face[k][j][i] = Wr[j];
      }
    }
  }

/*=== STEP 10: Compute second-order L/R x3-interface states ==================*/

/*--- Step 9a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * W1d = (d, V3, V1, V2, P, B1c, B2c, s[n])
 */

  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
      for (k=kl; k<=ku; k++) {
        W1d[k].d  = Whalf[k][j][i].d;
        W1d[k].Vx = Whalf[k][j][i].V3;
        W1d[k].Vy = Whalf[k][j][i].V1;
        W1d[k].Vz = Whalf[k][j][i].V2;
#ifndef BAROTROPIC
        W1d[k].P  = Whalf[k][j][i].P;
#endif /* BAROTROPIC */
#ifdef MHD
        W1d[k].By = Whalf[k][j][i].B1c;
        W1d[k].Bz = Whalf[k][j][i].B2c;
        Bxc[k] = Whalf[k][j][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) W1d[k].s[n] = Whalf[k][j][i].s[n];
#endif
      }

/*--- Step 9b ------------------------------------------------------------------
 * Compute L and R states at x3-interfaces, store in 3D array
 */
      lr_states(pG,W1d,Bxc,pG->dt,pG->dx3,ks,ke,Wl,Wr,cart_x3);

#ifdef FIRST_ORDER_FLUX_CORRECTION
    for (k=kl; k<=ku; k++) {
      Vsq = SQR(Wl[k].Vx) + SQR(Wl[k].Vy) + SQR(Wl[k].Vz);
      if (Vsq > 1.0){
        Wl[k] = W1d[k];
        Wr[k] = W1d[k];
      }
      Vsq = SQR(Wr[k].Vx) + SQR(Wr[k].Vy) + SQR(Wr[k].Vz);
      if (Vsq > 1.0){
        Wl[k] = W1d[k];
        Wr[k] = W1d[k];
      }
    }
#endif

      for (k=ks; k<=ke+1; k++) {
        Wl_x3Face[k][j][i] = Wl[k];
        Wr_x3Face[k][j][i] = Wr[k];
      }
    }
  }

/*=== STEP 11: Compute 3D x1-Flux, x2-Flux, x3-Flux ==========================*/

/*--- Step 11a -----------------------------------------------------------------
 * Compute maximum wavespeeds in multidimensions (eta in eq. 10 from Sanders et
 *  al. (1998)) for H-correction, if needed.
 */

#ifdef H_CORRECTION
  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=iu; i++) {
#ifdef MHD
        Bx = B1_x1Face[k][j][i];
#endif
        cfr = cfast(&(Ur_x1Face[k][j][i]),&Bx);
        cfl = cfast(&(Ul_x1Face[k][j][i]),&Bx);
        lambdar = Ur_x1Face[k][j][i].Mx/Ur_x1Face[k][j][i].d + cfr;
        lambdal = Ul_x1Face[k][j][i].Mx/Ul_x1Face[k][j][i].d - cfl;
        eta1[k][j][i] = 0.5*fabs(lambdar - lambdal);
      }
    }
  }

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=ju; j++) {
      for (i=is-1; i<=ie+1; i++) {
#ifdef MHD
        Bx = B2_x2Face[k][j][i];
#endif
        cfr = cfast(&(Ur_x2Face[k][j][i]),&Bx);
        cfl = cfast(&(Ul_x2Face[k][j][i]),&Bx);
        lambdar = Ur_x2Face[k][j][i].Mx/Ur_x2Face[k][j][i].d + cfr;
        lambdal = Ul_x2Face[k][j][i].Mx/Ul_x2Face[k][j][i].d - cfl;
        eta2[k][j][i] = 0.5*fabs(lambdar - lambdal);
      }
    }
  }

  for (k=ks-1; k<=ku; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
#ifdef MHD
        Bx = B3_x3Face[k][j][i];
#endif
        cfr = cfast(&(Ur_x3Face[k][j][i]),&Bx);
        cfl = cfast(&(Ul_x3Face[k][j][i]),&Bx);
        lambdar = Ur_x3Face[k][j][i].Mx/Ur_x3Face[k][j][i].d + cfr;
        lambdal = Ul_x3Face[k][j][i].Mx/Ul_x3Face[k][j][i].d - cfl;
        eta3[k][j][i] = 0.5*fabs(lambdar - lambdal);
      }
    }
  }
#endif /* H_CORRECTION */

/*--- Step 11b -----------------------------------------------------------------
 * Compute second-order fluxes in x1-direction
 */

#ifdef FIRST_ORDER_FLUX_CORRECTION
  NaNFlux = 0.0;
#endif
  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
#ifdef H_CORRECTION
        etah = MAX(eta2[k][j][i-1],eta2[k][j][i]);
        etah = MAX(etah,eta2[k][j+1][i-1]);
        etah = MAX(etah,eta2[k][j+1][i  ]);

        etah = MAX(etah,eta3[k  ][j][i-1]);
        etah = MAX(etah,eta3[k  ][j][i  ]);
        etah = MAX(etah,eta3[k+1][j][i-1]);
        etah = MAX(etah,eta3[k+1][j][i  ]);

        etah = MAX(etah,eta1[k  ][j][i  ]);
#endif /* H_CORRECTION */
#ifdef MHD
        Bx = B1_x1Face[k][j][i];
#endif
        Ul[i] = Prim1D_to_Cons1D(&Wl_x1Face[k][j][i],&Bx);
        Ur[i] = Prim1D_to_Cons1D(&Wr_x1Face[k][j][i],&Bx);

        fluxes(Ul[i],Ur[i],Wl_x1Face[k][j][i],Wr_x1Face[k][j][i],Bx,
               &x1Flux[k][j][i]);
#ifdef USE_ENTROPY_FIX
	entropy_flux(Ul[i],             Ur[i],
		     Wl_x1Face[k][j][i],Wr_x1Face[k][j][i],
		     Bx,                &x1FluxS[k][j][i]);
#endif


#ifdef FIRST_ORDER_FLUX_CORRECTION
/* revert to predictor flux if this flux Nan'ed */
        if ((x1Flux[k][j][i].d  != x1Flux[k][j][i].d)  ||
#ifndef BAROTROPIC
            (x1Flux[k][j][i].E  != x1Flux[k][j][i].E)  ||
#endif
#ifdef MHD
            (x1Flux[k][j][i].By != x1Flux[k][j][i].By) ||
            (x1Flux[k][j][i].Bz != x1Flux[k][j][i].Bz) ||
#endif
            (x1Flux[k][j][i].Mx != x1Flux[k][j][i].Mx) ||
            (x1Flux[k][j][i].My != x1Flux[k][j][i].My) ||
            (x1Flux[k][j][i].Mz != x1Flux[k][j][i].Mz)) {
          x1Flux[k][j][i] = x1FluxP[k][j][i];
#ifdef USE_ENTROPY_FIX
	  x1FluxS[k][j][i] = x1FluxSP[k][j][i];
#endif
          NaNFlux++;
        }
#endif

      }
    }
  }

#ifdef FIRST_ORDER_FLUX_CORRECTION
  if (NaNFlux != 0) {
    printf("[Step11b] %i second-order fluxes replaced\n",NaNFlux);
    NaNFlux=0;
  }
#endif

/*--- Step 11c -----------------------------------------------------------------
 * Compute second-order fluxes in x2-direction
 */

#ifdef FIRST_ORDER_FLUX_CORRECTION
  NaNFlux = 0.0;
#endif
  for (k=ks-1; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
#ifdef H_CORRECTION
        etah = MAX(eta1[k][j-1][i],eta1[k][j][i]);
        etah = MAX(etah,eta1[k][j-1][i+1]);
        etah = MAX(etah,eta1[k][j  ][i+1]);

        etah = MAX(etah,eta3[k  ][j-1][i]);
        etah = MAX(etah,eta3[k  ][j  ][i]);
        etah = MAX(etah,eta3[k+1][j-1][i]);
        etah = MAX(etah,eta3[k+1][j  ][i]);

        etah = MAX(etah,eta2[k  ][j  ][i]);
#endif /* H_CORRECTION */
#ifdef MHD
        Bx = B2_x2Face[k][j][i];
#endif
        Ul[i] = Prim1D_to_Cons1D(&Wl_x2Face[k][j][i],&Bx);
        Ur[i] = Prim1D_to_Cons1D(&Wr_x2Face[k][j][i],&Bx);

        fluxes(Ul[i],Ur[i],Wl_x2Face[k][j][i],Wr_x2Face[k][j][i],Bx,
               &x2Flux[k][j][i]);

#ifdef USE_ENTROPY_FIX
	entropy_flux(Ul[i],          Ur[i],
		     Wl_x2Face[k][j][i],Wr_x2Face[k][j][i],
		     Bx,                &x2FluxS[k][j][i]);
#endif


#ifdef FIRST_ORDER_FLUX_CORRECTION
/* revert to predictor flux if this flux NaN'ed */
        if ((x2Flux[k][j][i].d  != x2Flux[k][j][i].d)  ||
#ifndef BAROTROPIC
            (x2Flux[k][j][i].E  != x2Flux[k][j][i].E)  ||
#endif
#ifdef MHD
            (x2Flux[k][j][i].By != x2Flux[k][j][i].By) ||
            (x2Flux[k][j][i].Bz != x2Flux[k][j][i].Bz) ||
#endif
            (x2Flux[k][j][i].Mx != x2Flux[k][j][i].Mx) ||
            (x2Flux[k][j][i].My != x2Flux[k][j][i].My) ||
            (x2Flux[k][j][i].Mz != x2Flux[k][j][i].Mz)) {
          x2Flux[k][j][i] = x2FluxP[k][j][i];
#ifdef USE_ENTROPY_FIX
	  x2FluxS[k][j][i] = x2FluxSP[k][j][i];
#endif
          NaNFlux++;
        }
#endif

      }
    }
  }

#ifdef FIRST_ORDER_FLUX_CORRECTION
  if (NaNFlux != 0) {
    printf("[Step11c] %i second-order fluxes replaced\n",NaNFlux);
    NaNFlux=0;
  }
#endif

/*--- Step 11d -----------------------------------------------------------------
 * Compute second-order fluxes in x3-direction
 */

#ifdef FIRST_ORDER_FLUX_CORRECTION
  NaNFlux = 0.0;
#endif
  for (k=ks; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
#ifdef H_CORRECTION
        etah = MAX(eta1[k-1][j][i],eta1[k][j][i]);
        etah = MAX(etah,eta1[k-1][j][i+1]);
        etah = MAX(etah,eta1[k][j  ][i+1]);

        etah = MAX(etah,eta2[k-1][j  ][i]);
        etah = MAX(etah,eta2[k  ][j  ][i]);
        etah = MAX(etah,eta2[k-1][j+1][i]);
        etah = MAX(etah,eta2[k  ][j+1][i]);

        etah = MAX(etah,eta3[k  ][j  ][i]);
#endif /* H_CORRECTION */
#ifdef MHD
        Bx = B3_x3Face[k][j][i];
#endif
        Ul[i] = Prim1D_to_Cons1D(&Wl_x3Face[k][j][i],&Bx);
        Ur[i] = Prim1D_to_Cons1D(&Wr_x3Face[k][j][i],&Bx);

        fluxes(Ul[i],Ur[i],Wl_x3Face[k][j][i],Wr_x3Face[k][j][i],Bx,
               &x3Flux[k][j][i]);

#ifdef USE_ENTROPY_FIX
	entropy_flux(Ul[i],          Ur[i],
		     Wl_x3Face[k][j][i],Wr_x3Face[k][j][i],
		     Bx,                &x3FluxS[k][j][i]);
#endif

#ifdef FIRST_ORDER_FLUX_CORRECTION
/* revert to predictor flux if this flux NaN'ed */
        if ((x3Flux[k][j][i].d  != x3Flux[k][j][i].d)  ||
#ifndef BAROTROPIC
            (x3Flux[k][j][i].E  != x3Flux[k][j][i].E)  ||
#endif
#ifdef MHD
            (x3Flux[k][j][i].By != x3Flux[k][j][i].By) ||
            (x3Flux[k][j][i].Bz != x3Flux[k][j][i].Bz) ||
#endif
            (x3Flux[k][j][i].Mx != x3Flux[k][j][i].Mx) ||
            (x3Flux[k][j][i].My != x3Flux[k][j][i].My) ||
            (x3Flux[k][j][i].Mz != x3Flux[k][j][i].Mz)) {
          x3Flux[k][j][i] = x3FluxP[k][j][i];
#ifdef USE_ENTROPY_FIX
	  x3FluxS[k][j][i] = x3FluxSP[k][j][i];
#endif
          NaNFlux++;
        }
#endif

      }
    }
  }

#ifdef FIRST_ORDER_FLUX_CORRECTION
  if (NaNFlux != 0) {
    printf("[Step11d] %i second-order fluxes replaced\n",NaNFlux);
    NaNFlux=0;
  }
#endif

/*=== STEP 12: Update face-centered B for a full timestep ====================*/
        
/*--- Step 12a -----------------------------------------------------------------
 * Calculate the cell centered value of emf1,2,3 at the half-time-step.
 */

#ifdef MHD
  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
        emf1_cc[k][j][i] = (Whalf[k][j][i].B2c*Whalf[k][j][i].V3 - 
			    Whalf[k][j][i].B3c*Whalf[k][j][i].V2);
        emf2_cc[k][j][i] = (Whalf[k][j][i].B3c*Whalf[k][j][i].V1 -
			    Whalf[k][j][i].B1c*Whalf[k][j][i].V3);
        emf3_cc[k][j][i] = (Whalf[k][j][i].B1c*Whalf[k][j][i].V2 -
			    Whalf[k][j][i].B2c*Whalf[k][j][i].V1);
      }
    }
  }
#endif

/*--- Step 12b -----------------------------------------------------------------
 * Integrate emf*^{n+1/2} to the grid cell corners
 */

#ifdef MHD
  integrate_emf1_corner(pG);
  integrate_emf2_corner(pG);
  integrate_emf3_corner(pG);

/*--- Step 12c -----------------------------------------------------------------
 * Update the interface magnetic fields using CT for a full time step.
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->B1i[k][j][i] += dtodx3*(emf2[k+1][j  ][i  ] - emf2[k][j][i]) -
                            dtodx2*(emf3[k  ][j+1][i  ] - emf3[k][j][i]);
        pG->B2i[k][j][i] += dtodx1*(emf3[k  ][j  ][i+1] - emf3[k][j][i]) -
                            dtodx3*(emf1[k+1][j  ][i  ] - emf1[k][j][i]);
        pG->B3i[k][j][i] += dtodx2*(emf1[k  ][j+1][i  ] - emf1[k][j][i]) -
                            dtodx1*(emf2[k  ][j  ][i+1] - emf2[k][j][i]);
      }
      pG->B1i[k][j][ie+1] +=
        dtodx3*(emf2[k+1][j  ][ie+1] - emf2[k][j][ie+1]) -
        dtodx2*(emf3[k  ][j+1][ie+1] - emf3[k][j][ie+1]);
    }
    for (i=is; i<=ie; i++) {
      pG->B2i[k][je+1][i] +=
        dtodx1*(emf3[k  ][je+1][i+1] - emf3[k][je+1][i]) -
        dtodx3*(emf1[k+1][je+1][i  ] - emf1[k][je+1][i]);
    }
  }
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->B3i[ke+1][j][i] += 
        dtodx2*(emf1[ke+1][j+1][i  ] - emf1[ke+1][j][i]) -
        dtodx1*(emf2[ke+1][j  ][i+1] - emf2[ke+1][j][i]);
    }
  }

/*--- Step 12d -----------------------------------------------------------------
 * Set cell centered magnetic fields to average of updated face centered fields.
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].B1c = 0.5*(pG->B1i[k][j][i]+pG->B1i[k][j][i+1]);
        pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i]+pG->B2i[k][j+1][i]);
        pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i]+pG->B3i[k+1][j][i]);
      }
    }
  }
#endif

/*=== STEP 13: Add source terms for a full timestep using n+1/2 states =======*/
       
/*--- Step 13a -----------------------------------------------------------------
 * Add gravitational source terms due to a Static Potential
 * To improve conservation of total energy, we average the energy
 * source term computed at cell faces.
 *    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          cc_pos(pG,i,j,k,&x1,&x2,&x3);
          phic = (*StaticGravPot)(x1,x2,x3);
          phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
          phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

          pG->U[k][j][i].M1 -= dtodx1*(phir-phil)*Uhalf[k][j][i].d;
#ifndef BAROTROPIC
          pG->U[k][j][i].E -= dtodx1*(x1Flux[k][j][i  ].d*(phic - phil)
                                    + x1Flux[k][j][i+1].d*(phir - phic));
#endif
          phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
          phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

          pG->U[k][j][i].M2 -= dtodx2*(phir-phil)*Uhalf[k][j][i].d;
#ifndef BAROTROPIC
          pG->U[k][j][i].E -= dtodx2*(x2Flux[k][j  ][i].d*(phic - phil)
                                    + x2Flux[k][j+1][i].d*(phir - phic));
#endif
          phir = (*StaticGravPot)(x1,x2,(x3+0.5*pG->dx3));
          phil = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));

          pG->U[k][j][i].M3 -= dtodx3*(phir-phil)*Uhalf[k][j][i].d;
#ifndef BAROTROPIC
          pG->U[k][j][i].E -= dtodx3*(x3Flux[k  ][j][i].d*(phic - phil)
                                    + x3Flux[k+1][j][i].d*(phir - phic));
#endif
        }
      }
    }
  }

/*=== STEP 14: Update cell-centered values for a full timestep ===============*/

/*--- Step 14a -----------------------------------------------------------------
 * Update cell-centered variables in pG using 3D x1-Fluxes
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].d -=dtodx1*(x1Flux[k][j][i+1].d -x1Flux[k][j][i].d );
        pG->U[k][j][i].M1-=dtodx1*(x1Flux[k][j][i+1].Mx-x1Flux[k][j][i].Mx);
        pG->U[k][j][i].M2-=dtodx1*(x1Flux[k][j][i+1].My-x1Flux[k][j][i].My);
        pG->U[k][j][i].M3-=dtodx1*(x1Flux[k][j][i+1].Mz-x1Flux[k][j][i].Mz);
#ifndef BAROTROPIC
        pG->U[k][j][i].E -=dtodx1*(x1Flux[k][j][i+1].E -x1Flux[k][j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          pG->U[k][j][i].s[n] -= dtodx1*(x1Flux[k][j][i+1].s[n]
                                       - x1Flux[k][j][i  ].s[n]);
#endif
#ifdef USE_ENTROPY_FIX
	S[k][j][i] -= dtodx1*(x1FluxS[k][j][i+1]  - x1FluxS[k][j][i]);
#endif
      }
    }
  }

/*--- Step 14b -----------------------------------------------------------------
 * Update cell-centered variables in pG using 3D x2-Fluxes
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].d -=dtodx2*(x2Flux[k][j+1][i].d -x2Flux[k][j][i].d );
        pG->U[k][j][i].M1-=dtodx2*(x2Flux[k][j+1][i].Mz-x2Flux[k][j][i].Mz);
        pG->U[k][j][i].M2-=dtodx2*(x2Flux[k][j+1][i].Mx-x2Flux[k][j][i].Mx);
        pG->U[k][j][i].M3-=dtodx2*(x2Flux[k][j+1][i].My-x2Flux[k][j][i].My);
#ifndef BAROTROPIC
        pG->U[k][j][i].E -=dtodx2*(x2Flux[k][j+1][i].E -x2Flux[k][j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          pG->U[k][j][i].s[n] -= dtodx2*(x2Flux[k][j+1][i].s[n]
                                       - x2Flux[k][j  ][i].s[n]);
#endif
#ifdef USE_ENTROPY_FIX
      S[k][j][i] -= dtodx2*(x2FluxS[k][j+1][i]  - x2FluxS[k][j][i]);
#endif
      }
    }
  }

/*--- Step 14c -----------------------------------------------------------------
 * Update cell-centered variables in pG using 3D x3-Fluxes
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].d -=dtodx3*(x3Flux[k+1][j][i].d -x3Flux[k][j][i].d );
        pG->U[k][j][i].M1-=dtodx3*(x3Flux[k+1][j][i].My-x3Flux[k][j][i].My);
        pG->U[k][j][i].M2-=dtodx3*(x3Flux[k+1][j][i].Mz-x3Flux[k][j][i].Mz);
        pG->U[k][j][i].M3-=dtodx3*(x3Flux[k+1][j][i].Mx-x3Flux[k][j][i].Mx);
#ifndef BAROTROPIC
        pG->U[k][j][i].E -=dtodx3*(x3Flux[k+1][j][i].E -x3Flux[k][j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          pG->U[k][j][i].s[n] -= dtodx3*(x3Flux[k+1][j][i].s[n]
                                       - x3Flux[k  ][j][i].s[n]);
#endif
#ifdef USE_ENTROPY_FIX
      S[k][j][i] -= dtodx3*(x3FluxS[k+1][j][i]  - x3FluxS[k][j][i]);
#endif
      }
    }
  }

#ifdef FIRST_ORDER_FLUX_CORRECTION
/*=== STEP 15: First-order flux correction ===================================*/
/*--- Step 15a -----------------------------------------------------------------
 * If cell-centered d or P have gone negative, or if v^2 > 1 in SR, correct
 * by using 1st order predictor fluxes */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        Wcheck = check_Prim(&(pG->U[k][j][i]));
        if (Wcheck.d < 0.0) {
          flag_cell = 1;
          BadCell.i = i;
          BadCell.j = j;
          BadCell.k = k;
          negd++;
        }
        if (Wcheck.P < 0.0) {
          flag_cell = 1;
          BadCell.i = i;
          BadCell.j = j;
          BadCell.k = k;
          negP++;
        }
        Vsq = SQR(Wcheck.V1) + SQR(Wcheck.V2) + SQR(Wcheck.V3);
        if (Vsq > 1.0) {
          flag_cell = 1;
          BadCell.i = i;
          BadCell.j = j;
          BadCell.k = k;
          superl++; 
        }
        if (flag_cell != 0) {
          FixCell(pG, BadCell);
          flag_cell=0;
        }
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
 * adjacent cells and then calculating pressure. */

#ifdef MHD
  negd = 0;
  negP = 0;
  superl = 0;
  entropy = 0;
  final = 0;
  fail = 0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	flag_cell = 0;
        Wcheck = check_Prim(&(pG->U[k][j][i]));
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
	  Wcheck = entropy_fix (&(pG->U[k][j][i]),&(S[k][j][i]));
	  Ucheck = Prim_to_Cons(&Wcheck);
	  Wcheck = check_Prim(&Ucheck);
	  Vsq = SQR(Wcheck.V1) + SQR(Wcheck.V2) + SQR(Wcheck.V3);
	  if (Wcheck.d > 0.0 && Wcheck.P > 0.0 && Vsq < 1.0){
	    pG->U[k][j][i].d = Ucheck.d;
	    pG->U[k][j][i].M1 = Ucheck.M1;
	    pG->U[k][j][i].M2 = Ucheck.M2;
	    pG->U[k][j][i].M3 = Ucheck.M3;
	    pG->U[k][j][i].E = Ucheck.E;
	    flag_cell = 0;
	  }
	}
#endif /* USE_ENTROPY_FIX */
	if (flag_cell != 0) {
	  final++;
	  Wcheck = fix_vsq (&(pG->U[k][j][i]));
	  Ucheck = Prim_to_Cons(&Wcheck);
	  pG->U[k][j][i].d = Ucheck.d;
	  pG->U[k][j][i].M1 = Ucheck.M1;
	  pG->U[k][j][i].M2 = Ucheck.M2;
	  pG->U[k][j][i].M3 = Ucheck.M3;
	  pG->U[k][j][i].E = Ucheck.E;
	  Wcheck = check_Prim(&(pG->U[k][j][i]));
	  Vsq = SQR(Wcheck.V1) + SQR(Wcheck.V2) + SQR(Wcheck.V3);
	  if (Wcheck.d < 0.0 || Wcheck.P < 0.0 || Vsq > 1.0){
	    fail++;
	  }
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
        kcs = pG->CGrid[ncg].ijks[2];
        kce = pG->CGrid[ncg].ijke[2];

        for (k=kcs, kk=0; k<=kce; k++, kk++){
          for (j=jcs, jj=0; j<=jce; j++, jj++){
            pG->CGrid[ncg].myFlx[dim][kk][jj].d  = x1Flux[k][j][i].d;
            pG->CGrid[ncg].myFlx[dim][kk][jj].M1 = x1Flux[k][j][i].Mx;
            pG->CGrid[ncg].myFlx[dim][kk][jj].M2 = x1Flux[k][j][i].My;
            pG->CGrid[ncg].myFlx[dim][kk][jj].M3 = x1Flux[k][j][i].Mz;
#ifndef BAROTROPIC
            pG->CGrid[ncg].myFlx[dim][kk][jj].E  = x1Flux[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
            pG->CGrid[ncg].myFlx[dim][kk][jj].B1c = 0.0;
            pG->CGrid[ncg].myFlx[dim][kk][jj].B2c = x1Flux[k][j][i].By;
            pG->CGrid[ncg].myFlx[dim][kk][jj].B3c = x1Flux[k][j][i].Bz;
#endif /* MHD */
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++)
              pG->CGrid[ncg].myFlx[dim][kk][jj].s[n]  = x1Flux[k][j][i].s[n];
#endif
          }
        }
#ifdef MHD
        for (k=kcs, kk=0; k<=kce+1; k++, kk++){
          for (j=jcs, jj=0; j<=jce; j++, jj++){
            pG->CGrid[ncg].myEMF2[dim][kk][jj] = emf2[k][j][i];
          }
        }
        for (k=kcs, kk=0; k<=kce; k++, kk++){
          for (j=jcs, jj=0; j<=jce+1; j++, jj++){
            pG->CGrid[ncg].myEMF3[dim][kk][jj] = emf3[k][j][i];
          }
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
        kcs = pG->CGrid[ncg].ijks[2];
        kce = pG->CGrid[ncg].ijke[2];

        for (k=kcs, kk=0; k<=kce; k++, kk++){
          for (i=ics, ii=0; i<=ice; i++, ii++){
            pG->CGrid[ncg].myFlx[dim][kk][ii].d  = x2Flux[k][j][i].d;
            pG->CGrid[ncg].myFlx[dim][kk][ii].M1 = x2Flux[k][j][i].Mz;
            pG->CGrid[ncg].myFlx[dim][kk][ii].M2 = x2Flux[k][j][i].Mx;
            pG->CGrid[ncg].myFlx[dim][kk][ii].M3 = x2Flux[k][j][i].My;
#ifndef BAROTROPIC
            pG->CGrid[ncg].myFlx[dim][kk][ii].E  = x2Flux[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
            pG->CGrid[ncg].myFlx[dim][kk][ii].B1c = x2Flux[k][j][i].Bz;
            pG->CGrid[ncg].myFlx[dim][kk][ii].B2c = 0.0;
            pG->CGrid[ncg].myFlx[dim][kk][ii].B3c = x2Flux[k][j][i].By;
#endif /* MHD */
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++)
              pG->CGrid[ncg].myFlx[dim][kk][ii].s[n]  = x2Flux[k][j][i].s[n];
#endif
          }
        }
#ifdef MHD
        for (k=kcs, kk=0; k<=kce+1; k++, kk++){
          for (i=ics, ii=0; i<=ice; i++, ii++){
            pG->CGrid[ncg].myEMF1[dim][kk][ii] = emf1[k][j][i];
          }
        }
        for (k=kcs, kk=0; k<=kce; k++, kk++){
          for (i=ics, ii=0; i<=ice+1; i++, ii++){
            pG->CGrid[ncg].myEMF3[dim][kk][ii] = emf3[k][j][i];
          }
        }
#endif /* MHD */
      }
    }

/* x3-boundaries of child Grids (interior to THIS Grid) */

    for (dim=4; dim<6; dim++){
      if (pG->CGrid[ncg].myFlx[dim] != NULL) {

        ics = pG->CGrid[ncg].ijks[0];
        ice = pG->CGrid[ncg].ijke[0];
        jcs = pG->CGrid[ncg].ijks[1];
        jce = pG->CGrid[ncg].ijke[1];
        if (dim==4) k = pG->CGrid[ncg].ijks[2];
        if (dim==5) k = pG->CGrid[ncg].ijke[2] + 1;

        for (j=jcs, jj=0; j<=jce; j++, jj++){
          for (i=ics, ii=0; i<=ice; i++, ii++){
            pG->CGrid[ncg].myFlx[dim][jj][ii].d  = x3Flux[k][j][i].d;
            pG->CGrid[ncg].myFlx[dim][jj][ii].M1 = x3Flux[k][j][i].My;
            pG->CGrid[ncg].myFlx[dim][jj][ii].M2 = x3Flux[k][j][i].Mz;
            pG->CGrid[ncg].myFlx[dim][jj][ii].M3 = x3Flux[k][j][i].Mx;
#ifndef BAROTROPIC
            pG->CGrid[ncg].myFlx[dim][jj][ii].E  = x3Flux[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
            pG->CGrid[ncg].myFlx[dim][jj][ii].B1c = x3Flux[k][j][i].By;
            pG->CGrid[ncg].myFlx[dim][jj][ii].B2c = x3Flux[k][j][i].Bz;
            pG->CGrid[ncg].myFlx[dim][jj][ii].B3c = 0.0;
#endif /* MHD */
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++)
              pG->CGrid[ncg].myFlx[dim][jj][ii].s[n]  = x3Flux[k][j][i].s[n];
#endif
          }
        }
#ifdef MHD
        for (j=jcs, jj=0; j<=jce+1; j++, jj++){
          for (i=ics, ii=0; i<=ice; i++, ii++){
            pG->CGrid[ncg].myEMF1[dim][jj][ii] = emf1[k][j][i];
          }
        }
        for (j=jcs, jj=0; j<=jce; j++, jj++){
          for (i=ics, ii=0; i<=ice+1; i++, ii++){
            pG->CGrid[ncg].myEMF2[dim][jj][ii] = emf2[k][j][i];
          }
        }
#endif /* MHD */
      }
    }
  } /* end loop over child Grids */

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
        kps = pG->PGrid[npg].ijks[2];
        kpe = pG->PGrid[npg].ijke[2];

        for (k=kps, kk=0; k<=kpe; k++, kk++){
          for (j=jps, jj=0; j<=jpe; j++, jj++){
            pG->PGrid[npg].myFlx[dim][kk][jj].d  = x1Flux[k][j][i].d;
            pG->PGrid[npg].myFlx[dim][kk][jj].M1 = x1Flux[k][j][i].Mx;
            pG->PGrid[npg].myFlx[dim][kk][jj].M2 = x1Flux[k][j][i].My;
            pG->PGrid[npg].myFlx[dim][kk][jj].M3 = x1Flux[k][j][i].Mz;
#ifndef BAROTROPIC
            pG->PGrid[npg].myFlx[dim][kk][jj].E  = x1Flux[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
            pG->PGrid[npg].myFlx[dim][kk][jj].B1c = 0.0;
            pG->PGrid[npg].myFlx[dim][kk][jj].B2c = x1Flux[k][j][i].By;
            pG->PGrid[npg].myFlx[dim][kk][jj].B3c = x1Flux[k][j][i].Bz;
#endif /* MHD */
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++)
              pG->PGrid[npg].myFlx[dim][kk][jj].s[n]  = x1Flux[k][j][i].s[n];
#endif
          }
        }
#ifdef MHD
        for (k=kps, kk=0; k<=kpe+1; k++, kk++){
          for (j=jps, jj=0; j<=jpe; j++, jj++){
            pG->PGrid[npg].myEMF2[dim][kk][jj] = emf2[k][j][i];
          }
        }
        for (k=kps, kk=0; k<=kpe; k++, kk++){
          for (j=jps, jj=0; j<=jpe+1; j++, jj++){
            pG->PGrid[npg].myEMF3[dim][kk][jj] = emf3[k][j][i];
          }
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
        kps = pG->PGrid[npg].ijks[2];
        kpe = pG->PGrid[npg].ijke[2];

        for (k=kps, kk=0; k<=kpe; k++, kk++){
          for (i=ips, ii=0; i<=ipe; i++, ii++){
            pG->PGrid[npg].myFlx[dim][kk][ii].d  = x2Flux[k][j][i].d;
            pG->PGrid[npg].myFlx[dim][kk][ii].M1 = x2Flux[k][j][i].Mz;
            pG->PGrid[npg].myFlx[dim][kk][ii].M2 = x2Flux[k][j][i].Mx;
            pG->PGrid[npg].myFlx[dim][kk][ii].M3 = x2Flux[k][j][i].My;
#ifndef BAROTROPIC
            pG->PGrid[npg].myFlx[dim][kk][ii].E  = x2Flux[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
            pG->PGrid[npg].myFlx[dim][kk][ii].B1c = x2Flux[k][j][i].Bz;
            pG->PGrid[npg].myFlx[dim][kk][ii].B2c = 0.0;
            pG->PGrid[npg].myFlx[dim][kk][ii].B3c = x2Flux[k][j][i].By;
#endif /* MHD */
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++)
              pG->PGrid[npg].myFlx[dim][kk][ii].s[n]  = x2Flux[k][j][i].s[n];
#endif
          }
        }
#ifdef MHD
        for (k=kps, kk=0; k<=kpe+1; k++, kk++){
          for (i=ips, ii=0; i<=ipe; i++, ii++){
            pG->PGrid[npg].myEMF1[dim][kk][ii] = emf1[k][j][i];
          }
        }
        for (k=kps, kk=0; k<=kpe; k++, kk++){
          for (i=ips, ii=0; i<=ipe+1; i++, ii++){
            pG->PGrid[npg].myEMF3[dim][kk][ii] = emf3[k][j][i];
          }
        }
#endif /* MHD */
      }
    }

/* x3-boundaries of parent Grids (at boundaries of THIS Grid)  */
 
    for (dim=4; dim<6; dim++){
      if (pG->PGrid[npg].myFlx[dim] != NULL) {

        ips = pG->PGrid[npg].ijks[0];
        ipe = pG->PGrid[npg].ijke[0];
        jps = pG->PGrid[npg].ijks[1];
        jpe = pG->PGrid[npg].ijke[1];
        if (dim==4) k = pG->PGrid[npg].ijks[2];
        if (dim==5) k = pG->PGrid[npg].ijke[2] + 1;

        for (j=jps, jj=0; j<=jpe; j++, jj++){
          for (i=ips, ii=0; i<=ipe; i++, ii++){
            pG->PGrid[npg].myFlx[dim][jj][ii].d  = x3Flux[k][j][i].d;
            pG->PGrid[npg].myFlx[dim][jj][ii].M1 = x3Flux[k][j][i].My;
            pG->PGrid[npg].myFlx[dim][jj][ii].M2 = x3Flux[k][j][i].Mz;
            pG->PGrid[npg].myFlx[dim][jj][ii].M3 = x3Flux[k][j][i].Mx;
#ifndef BAROTROPIC
            pG->PGrid[npg].myFlx[dim][jj][ii].E  = x3Flux[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
            pG->PGrid[npg].myFlx[dim][jj][ii].B1c = x3Flux[k][j][i].By;
            pG->PGrid[npg].myFlx[dim][jj][ii].B2c = x3Flux[k][j][i].Bz;
            pG->PGrid[npg].myFlx[dim][jj][ii].B3c = 0.0;
#endif /* MHD */
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++)
              pG->PGrid[npg].myFlx[dim][jj][ii].s[n]  = x3Flux[k][j][i].s[n];
#endif
          }
        }
#ifdef MHD
        for (j=jps, jj=0; j<=jpe+1; j++, jj++){
          for (i=ips, ii=0; i<=ipe; i++, ii++){
            pG->PGrid[npg].myEMF1[dim][jj][ii] = emf1[k][j][i];
          }
        }
        for (j=jps, jj=0; j<=jpe; j++, jj++){
          for (i=ips, ii=0; i<=ipe+1; i++, ii++){
            pG->PGrid[npg].myEMF2[dim][jj][ii] = emf2[k][j][i];
          }
        }
#endif /* MHD */
      }
    }
  }

#endif /* STATIC_MESH_REFINEMENT */

  return;
}


/*----------------------------------------------------------------------------*/
/*! \fn void integrate_init_3d(MeshS *pM)
 *  \brief Allocate temporary integration arrays */
void integrate_init_3d(MeshS *pM)
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

#ifdef MHD
  if ((emf1 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf2 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf3 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;

  if ((emf1_cc=(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf2_cc=(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf3_cc=(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
#endif /* MHD */
#ifdef H_CORRECTION
  if ((eta1 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((eta2 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((eta3 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
#endif /* H_CORRECTION */
  if ((Wl_x1Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wr_x1Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wl_x2Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wr_x2Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wl_x3Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wr_x3Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;

  if ((Bxc = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;

#ifdef MHD
  if ((B1_x1Face = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))
    == NULL) goto on_error;
  if ((B2_x2Face = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))
    == NULL) goto on_error;
  if ((B3_x3Face = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))
    == NULL) goto on_error;
#endif /* MHD */

  if ((U1d = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ul  = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ur  = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W1d = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wl  = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wr  = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((x1Flux = (Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x2Flux = (Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x3Flux = (Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
#ifdef FIRST_ORDER_FLUX_CORRECTION
  if ((x1FluxP =(Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x2FluxP =(Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x3FluxP =(Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
#ifdef MHD
  if ((emf1P = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
  if ((emf2P = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
  if ((emf3P = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
#endif
#endif /* FIRST_ORDER_FLUX_CORRECTION */
#ifdef USE_ENTROPY_FIX
  if ((S = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
      == NULL) goto on_error;
  if ((Shalf = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
      == NULL) goto on_error;
  if ((x1FluxS = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
      == NULL) goto on_error;
  if ((x2FluxS = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
      == NULL) goto on_error;
  if ((x3FluxS = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
      == NULL) goto on_error;
#ifdef FIRST_ORDER_FLUX_CORRECTION
  if ((x1FluxSP = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
      == NULL) goto on_error;
  if ((x2FluxSP = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
      == NULL) goto on_error;
  if ((x3FluxSP = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
      == NULL) goto on_error;
#endif
#endif

  if ((Uhalf = (ConsS***)calloc_3d_array(size3,size2,size1,sizeof(ConsS)))
    == NULL) goto on_error;
  if ((Whalf = (PrimS***)calloc_3d_array(size3,size2,size1,sizeof(PrimS)))
    == NULL) goto on_error;
  if ((W     = (PrimS***)calloc_3d_array(size3,size2,size1,sizeof(PrimS)))
    == NULL) goto on_error;

  return;

  on_error:
    integrate_destruct();
    ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_destruct_3d(void)
 *  \brief Free temporary integration arrays */
void integrate_destruct_3d(void)
{
#ifdef MHD
  if (emf1 != NULL) free_3d_array(emf1);
  if (emf2 != NULL) free_3d_array(emf2);
  if (emf3 != NULL) free_3d_array(emf3);
  if (emf1_cc != NULL) free_3d_array(emf1_cc);
  if (emf2_cc != NULL) free_3d_array(emf2_cc);
  if (emf3_cc != NULL) free_3d_array(emf3_cc);
#endif /* MHD */
#ifdef H_CORRECTION
  if (eta1 != NULL) free_3d_array(eta1);
  if (eta2 != NULL) free_3d_array(eta2);
  if (eta3 != NULL) free_3d_array(eta3);
#endif /* H_CORRECTION */
  if (Wl_x1Face != NULL) free_3d_array(Wl_x1Face);
  if (Wr_x1Face != NULL) free_3d_array(Wr_x1Face);
  if (Wl_x2Face != NULL) free_3d_array(Wl_x2Face);
  if (Wr_x2Face != NULL) free_3d_array(Wr_x2Face);
  if (Wl_x3Face != NULL) free_3d_array(Wl_x3Face);
  if (Wr_x3Face != NULL) free_3d_array(Wr_x3Face);

  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);
#ifdef MHD
  if (B1_x1Face != NULL) free_3d_array(B1_x1Face);
  if (B2_x2Face != NULL) free_3d_array(B2_x2Face);
  if (B3_x3Face != NULL) free_3d_array(B3_x3Face);
#endif /* MHD */

  if (U1d != NULL) free(U1d);
  if (Ul  != NULL) free(Ul);
  if (Ur  != NULL) free(Ur);
  if (W1d != NULL) free(W1d);
  if (Wl  != NULL) free(Wl);
  if (Wr  != NULL) free(Wr);

  if (x1Flux  != NULL) free_3d_array(x1Flux);
  if (x2Flux  != NULL) free_3d_array(x2Flux);
  if (x3Flux  != NULL) free_3d_array(x3Flux);
#ifdef FIRST_ORDER_FLUX_CORRECTION
  if (x1FluxP != NULL) free_3d_array(x1FluxP);
  if (x2FluxP != NULL) free_3d_array(x2FluxP);
  if (x3FluxP != NULL) free_3d_array(x3FluxP);
#ifdef MHD
  if (emf1P   != NULL) free_3d_array(emf1P);
  if (emf2P   != NULL) free_3d_array(emf2P);
  if (emf3P   != NULL) free_3d_array(emf3P);
#endif
#endif /* FIRST_ORDER_FLUX_CORRECTION */
#ifdef USE_ENTROPY_FIX
  if (S       != NULL) free_3d_array(S);
  if (Shalf   != NULL) free_3d_array(Shalf);
  if (x1FluxS != NULL) free_3d_array(x1FluxS);
  if (x2FluxS != NULL) free_3d_array(x2FluxS);
  if (x3FluxS != NULL) free_3d_array(x3FluxS);
#ifdef FIRST_ORDER_FLUX_CORRECTION
  if (x1FluxSP!= NULL) free_3d_array(x1FluxSP);
  if (x2FluxSP!= NULL) free_3d_array(x2FluxSP);
  if (x3FluxSP!= NULL) free_3d_array(x3FluxSP);
#endif
#endif

  if (Uhalf  != NULL) free_3d_array(Uhalf);
  if (Whalf  != NULL) free_3d_array(Whalf);
  if (W      != NULL) free_3d_array(W);

  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/

#ifdef MHD
/*! \fn static void integrate_emf1_corner(const GridS *pG)
 *  \brief Integrates face centered B-fluxes to compute corner EMFs.  
 *
 *  Note:
 * - x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 * - x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 * - x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 * - x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 * - x3Flux.By = VxBy - BxVy = v3*b1-b3*v1 = -EMFY
 * - x3Flux.Bz = VxBz - BxVz = v3*b2-b3*v2 = EMFX
 *- first_order_correction() - Added by N. Lemaster to run supersonic turbulence
 */ 


static void integrate_emf1_corner(const GridS *pG)
{
  int i,il,iu,j,jl,ju,k,kl,ku;
  Real de1_l2, de1_r2, de1_l3, de1_r3;

  il = pG->is-(nghost-1);   iu = pG->ie+(nghost-1);
  jl = pG->js-(nghost-1);   ju = pG->je+(nghost-1);
  kl = pG->ks-(nghost-1);   ku = pG->ke+(nghost-1);

  for (k=kl; k<=ku+1; k++) {
    for (j=jl; j<=ju+1; j++) {
      for (i=il; i<=iu; i++) {
/* NOTE: The x2-Flux of By is -E1. */
/*       The x3-Flux of Bz is +E1. */
	if (x2Flux[k-1][j][i].d > 0.0)
	  de1_l3 = x3Flux[k][j-1][i].Bz - emf1_cc[k-1][j-1][i];
	else if (x2Flux[k-1][j][i].d < 0.0)
	  de1_l3 = x3Flux[k][j][i].Bz - emf1_cc[k-1][j][i];
	else {
	  de1_l3 = 0.5*(x3Flux[k][j-1][i].Bz - emf1_cc[k-1][j-1][i] +
			x3Flux[k][j  ][i].Bz - emf1_cc[k-1][j  ][i] );
	}

	if (x2Flux[k][j][i].d > 0.0)
	  de1_r3 = x3Flux[k][j-1][i].Bz - emf1_cc[k][j-1][i];
	else if (x2Flux[k][j][i].d < 0.0)
	  de1_r3 = x3Flux[k][j][i].Bz - emf1_cc[k][j][i];
	else {
	  de1_r3 = 0.5*(x3Flux[k][j-1][i].Bz - emf1_cc[k][j-1][i] +
			x3Flux[k][j  ][i].Bz - emf1_cc[k][j  ][i] );
	}

	if (x3Flux[k][j-1][i].d > 0.0)
	  de1_l2 = -x2Flux[k-1][j][i].By - emf1_cc[k-1][j-1][i];
	else if (x3Flux[k][j-1][i].d < 0.0)
	  de1_l2 = -x2Flux[k][j][i].By - emf1_cc[k][j-1][i];
	else {
	  de1_l2 = 0.5*(-x2Flux[k-1][j][i].By - emf1_cc[k-1][j-1][i]
			-x2Flux[k  ][j][i].By - emf1_cc[k  ][j-1][i] );
	}

	if (x3Flux[k][j][i].d > 0.0)
	  de1_r2 = -x2Flux[k-1][j][i].By - emf1_cc[k-1][j][i];
	else if (x3Flux[k][j][i].d < 0.0)
	  de1_r2 = -x2Flux[k][j][i].By - emf1_cc[k][j][i];
	else {
	  de1_r2 = 0.5*(-x2Flux[k-1][j][i].By - emf1_cc[k-1][j][i]
			-x2Flux[k  ][j][i].By - emf1_cc[k  ][j][i] );
	}

        emf1[k][j][i] = 0.25*(  x3Flux[k][j][i].Bz + x3Flux[k][j-1][i].Bz
                              - x2Flux[k][j][i].By - x2Flux[k-1][j][i].By 
			      + de1_l2 + de1_r2 + de1_l3 + de1_r3);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void integrate_emf2_corner(const GridS *pG)
 *  \brief Integrates face centered B-fluxes to compute corner EMFs.  
 *
 *  Note:
 * - x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 * - x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 * - x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 * - x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 * - x3Flux.By = VxBy - BxVy = v3*b1-b3*v1 = -EMFY
 * - x3Flux.Bz = VxBz - BxVz = v3*b2-b3*v2 = EMFX
 *- first_order_correction() - Added by N. Lemaster to run supersonic turbulence
 */ 

static void integrate_emf2_corner(const GridS *pG)
{
  int i,il,iu,j,jl,ju,k,kl,ku;
  Real de2_l1, de2_r1, de2_l3, de2_r3;

  il = pG->is-(nghost-1);   iu = pG->ie+(nghost-1);
  jl = pG->js-(nghost-1);   ju = pG->je+(nghost-1);
  kl = pG->ks-(nghost-1);   ku = pG->ke+(nghost-1);

  for (k=kl; k<=ku+1; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu+1; i++) {
/* NOTE: The x1-Flux of Bz is +E2. */
/*       The x3-Flux of By is -E2. */
	if (x1Flux[k-1][j][i].d > 0.0)
	  de2_l3 = -x3Flux[k][j][i-1].By - emf2_cc[k-1][j][i-1];
	else if (x1Flux[k-1][j][i].d < 0.0)
	  de2_l3 = -x3Flux[k][j][i].By - emf2_cc[k-1][j][i];
	else {
	  de2_l3 = 0.5*(-x3Flux[k][j][i-1].By - emf2_cc[k-1][j][i-1] 
			-x3Flux[k][j][i  ].By - emf2_cc[k-1][j][i  ] );
	}

	if (x1Flux[k][j][i].d > 0.0)
	  de2_r3 = -x3Flux[k][j][i-1].By - emf2_cc[k][j][i-1];
	else if (x1Flux[k][j][i].d < 0.0)
	  de2_r3 = -x3Flux[k][j][i].By - emf2_cc[k][j][i];
	else {
	  de2_r3 = 0.5*(-x3Flux[k][j][i-1].By - emf2_cc[k][j][i-1] 
			-x3Flux[k][j][i  ].By - emf2_cc[k][j][i  ] );
	}

	if (x3Flux[k][j][i-1].d > 0.0)
	  de2_l1 = x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i-1];
	else if (x3Flux[k][j][i-1].d < 0.0)
	  de2_l1 = x1Flux[k][j][i].Bz - emf2_cc[k][j][i-1];
	else {
	  de2_l1 = 0.5*(x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i-1] +
			x1Flux[k  ][j][i].Bz - emf2_cc[k  ][j][i-1] );
	}

	if (x3Flux[k][j][i].d > 0.0)
	  de2_r1 = x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i];
	else if (x3Flux[k][j][i].d < 0.0)
	  de2_r1 = x1Flux[k][j][i].Bz - emf2_cc[k][j][i];
	else {
	  de2_r1 = 0.5*(x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i] +
			x1Flux[k  ][j][i].Bz - emf2_cc[k  ][j][i] );
	}

	emf2[k][j][i] = 0.25*(  x1Flux[k][j][i].Bz + x1Flux[k-1][j][i  ].Bz
                              - x3Flux[k][j][i].By - x3Flux[k  ][j][i-1].By
			      + de2_l1 + de2_r1 + de2_l3 + de2_r3);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/

/*! \fn static void integrate_emf3_corner(const GridS *pG)
 *  \brief Integrates face centered B-fluxes to compute corner EMFs.  
 *
 *  Note:
 * - x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 * - x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 * - x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 * - x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 * - x3Flux.By = VxBy - BxVy = v3*b1-b3*v1 = -EMFY
 * - x3Flux.Bz = VxBz - BxVz = v3*b2-b3*v2 = EMFX
 *- first_order_correction() - Added by N. Lemaster to run supersonic turbulence
 */ 
static void integrate_emf3_corner(const GridS *pG)
{
  int i,il,iu,j,jl,ju,k,kl,ku;
  Real de3_l1, de3_r1, de3_l2, de3_r2;

  il = pG->is-(nghost-1);   iu = pG->ie+(nghost-1);
  jl = pG->js-(nghost-1);   ju = pG->je+(nghost-1);
  kl = pG->ks-(nghost-1);   ku = pG->ke+(nghost-1);

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju+1; j++) {
      for (i=il; i<=iu+1; i++) {
/* NOTE: The x1-Flux of By is -E3. */
/*       The x2-Flux of Bx is +E3. */
	if (x1Flux[k][j-1][i].d > 0.0)
	  de3_l2 = x2Flux[k][j][i-1].Bz - emf3_cc[k][j-1][i-1];
	else if (x1Flux[k][j-1][i].d < 0.0)
	  de3_l2 = x2Flux[k][j][i].Bz - emf3_cc[k][j-1][i];
	else {
	  de3_l2 = 0.5*(x2Flux[k][j][i-1].Bz - emf3_cc[k][j-1][i-1] + 
			x2Flux[k][j][i  ].Bz - emf3_cc[k][j-1][i  ] );
	}

	if (x1Flux[k][j][i].d > 0.0)
	  de3_r2 = x2Flux[k][j][i-1].Bz - emf3_cc[k][j][i-1];
	else if (x1Flux[k][j][i].d < 0.0)
	  de3_r2 = x2Flux[k][j][i].Bz - emf3_cc[k][j][i];
	else {
	  de3_r2 = 0.5*(x2Flux[k][j][i-1].Bz - emf3_cc[k][j][i-1] + 
			x2Flux[k][j][i  ].Bz - emf3_cc[k][j][i  ] );
	}

	if (x2Flux[k][j][i-1].d > 0.0)
	  de3_l1 = -x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i-1];
	else if (x2Flux[k][j][i-1].d < 0.0)
	  de3_l1 = -x1Flux[k][j][i].By - emf3_cc[k][j][i-1];
	else {
	  de3_l1 = 0.5*(-x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i-1]
			-x1Flux[k][j  ][i].By - emf3_cc[k][j  ][i-1] );
	}

	if (x2Flux[k][j][i].d > 0.0)
	  de3_r1 = -x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i];
	else if (x2Flux[k][j][i].d < 0.0)
	  de3_r1 = -x1Flux[k][j][i].By - emf3_cc[k][j][i];
	else {
	  de3_r1 = 0.5*(-x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i]
			-x1Flux[k][j  ][i].By - emf3_cc[k][j  ][i] );
	}

	emf3[k][j][i] = 0.25*(  x2Flux[k][j  ][i-1].Bz + x2Flux[k][j][i].Bz
			      - x1Flux[k][j-1][i  ].By - x1Flux[k][j][i].By
			      + de3_l1 + de3_r1 + de3_l2 + de3_r2);
      }
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
  Real dtodx1=pG->dt/pG->dx1, dtodx2=pG->dt/pG->dx2, dtodx3=pG->dt/pG->dx3;
  Cons1DS x1FD_i, x1FD_ip1, x2FD_j, x2FD_jp1, x3FD_k, x3FD_kp1;
#ifdef MHD
  int i,j,k;
  Real emf1D_kj, emf1D_kp1j, emf1D_kjp1, emf1D_kp1jp1;
  Real emf2D_ki, emf2D_kip1, emf2D_kp1i, emf2D_kp1ip1;
  Real emf3D_ji, emf3D_jip1, emf3D_jp1i, emf3D_jp1ip1;
#endif

/* Compute difference of predictor and corrector fluxes at cell faces */

  x1FD_i.d = x1Flux[ix.k][ix.j][ix.i].d - x1FluxP[ix.k][ix.j][ix.i].d;
  x2FD_j.d = x2Flux[ix.k][ix.j][ix.i].d - x2FluxP[ix.k][ix.j][ix.i].d;
  x3FD_k.d = x3Flux[ix.k][ix.j][ix.i].d - x3FluxP[ix.k][ix.j][ix.i].d;

  x1FD_ip1.d = x1Flux[ix.k][ix.j][ix.i+1].d - x1FluxP[ix.k][ix.j][ix.i+1].d;
  x2FD_jp1.d = x2Flux[ix.k][ix.j+1][ix.i].d - x2FluxP[ix.k][ix.j+1][ix.i].d;
  x3FD_kp1.d = x2Flux[ix.k+1][ix.j][ix.i].d - x3FluxP[ix.k+1][ix.j][ix.i].d;

  x1FD_i.Mx = x1Flux[ix.k][ix.j][ix.i].Mx - x1FluxP[ix.k][ix.j][ix.i].Mx;
  x2FD_j.Mx = x2Flux[ix.k][ix.j][ix.i].Mx - x2FluxP[ix.k][ix.j][ix.i].Mx;
  x3FD_k.Mx = x3Flux[ix.k][ix.j][ix.i].Mx - x3FluxP[ix.k][ix.j][ix.i].Mx;

  x1FD_ip1.Mx = x1Flux[ix.k][ix.j][ix.i+1].Mx - x1FluxP[ix.k][ix.j][ix.i+1].Mx;
  x2FD_jp1.Mx = x2Flux[ix.k][ix.j+1][ix.i].Mx - x2FluxP[ix.k][ix.j+1][ix.i].Mx;
  x3FD_kp1.Mx = x2Flux[ix.k+1][ix.j][ix.i].Mx - x3FluxP[ix.k+1][ix.j][ix.i].Mx;

  x1FD_i.My = x1Flux[ix.k][ix.j][ix.i].My - x1FluxP[ix.k][ix.j][ix.i].My;
  x2FD_j.My = x2Flux[ix.k][ix.j][ix.i].My - x2FluxP[ix.k][ix.j][ix.i].My;
  x3FD_k.My = x3Flux[ix.k][ix.j][ix.i].My - x3FluxP[ix.k][ix.j][ix.i].My;

  x1FD_ip1.My = x1Flux[ix.k][ix.j][ix.i+1].My - x1FluxP[ix.k][ix.j][ix.i+1].My;
  x2FD_jp1.My = x2Flux[ix.k][ix.j+1][ix.i].My - x2FluxP[ix.k][ix.j+1][ix.i].My;
  x3FD_kp1.My = x2Flux[ix.k+1][ix.j][ix.i].My - x3FluxP[ix.k+1][ix.j][ix.i].My;

  x1FD_i.Mz = x1Flux[ix.k][ix.j][ix.i].Mz - x1FluxP[ix.k][ix.j][ix.i].Mz;
  x2FD_j.Mz = x2Flux[ix.k][ix.j][ix.i].Mz - x2FluxP[ix.k][ix.j][ix.i].Mz;
  x3FD_k.Mz = x3Flux[ix.k][ix.j][ix.i].Mz - x3FluxP[ix.k][ix.j][ix.i].Mz;

  x1FD_ip1.Mz = x1Flux[ix.k][ix.j][ix.i+1].Mz - x1FluxP[ix.k][ix.j][ix.i+1].Mz;
  x2FD_jp1.Mz = x2Flux[ix.k][ix.j+1][ix.i].Mz - x2FluxP[ix.k][ix.j+1][ix.i].Mz;
  x3FD_kp1.Mz = x2Flux[ix.k+1][ix.j][ix.i].Mz - x3FluxP[ix.k+1][ix.j][ix.i].Mz;

#ifndef BAROTROPIC
  x1FD_i.E = x1Flux[ix.k][ix.j][ix.i].E - x1FluxP[ix.k][ix.j][ix.i].E;
  x2FD_j.E = x2Flux[ix.k][ix.j][ix.i].E - x2FluxP[ix.k][ix.j][ix.i].E;
  x3FD_k.E = x3Flux[ix.k][ix.j][ix.i].E - x3FluxP[ix.k][ix.j][ix.i].E;

  x1FD_ip1.E = x1Flux[ix.k][ix.j][ix.i+1].E - x1FluxP[ix.k][ix.j][ix.i+1].E;
  x2FD_jp1.E = x2Flux[ix.k][ix.j+1][ix.i].E - x2FluxP[ix.k][ix.j+1][ix.i].E;
  x3FD_kp1.E = x2Flux[ix.k+1][ix.j][ix.i].E - x3FluxP[ix.k+1][ix.j][ix.i].E;
#endif /* BAROTROPIC */

#ifdef MHD
  emf1D_kj     = emf1[ix.k  ][ix.j  ][ix.i] - emf1P[ix.k  ][ix.j  ][ix.i];
  emf1D_kjp1   = emf1[ix.k  ][ix.j+1][ix.i] - emf1P[ix.k  ][ix.j+1][ix.i];
  emf1D_kp1j   = emf1[ix.k+1][ix.j  ][ix.i] - emf1P[ix.k+1][ix.j  ][ix.i];
  emf1D_kp1jp1 = emf1[ix.k+1][ix.j+1][ix.i] - emf1P[ix.k+1][ix.j+1][ix.i];

  emf2D_ki     = emf2[ix.k  ][ix.j][ix.i  ] - emf2P[ix.k  ][ix.j][ix.i  ];
  emf2D_kip1   = emf2[ix.k  ][ix.j][ix.i+1] - emf2P[ix.k  ][ix.j][ix.i+1];
  emf2D_kp1i   = emf2[ix.k+1][ix.j][ix.i  ] - emf2P[ix.k+1][ix.j][ix.i  ];
  emf2D_kp1ip1 = emf2[ix.k+1][ix.j][ix.i+1] - emf2P[ix.k+1][ix.j][ix.i+1];

  emf3D_ji     = emf3[ix.k][ix.j  ][ix.i  ] - emf3P[ix.k][ix.j  ][ix.i  ];
  emf3D_jip1   = emf3[ix.k][ix.j  ][ix.i+1] - emf3P[ix.k][ix.j  ][ix.i+1];
  emf3D_jp1i   = emf3[ix.k][ix.j+1][ix.i  ] - emf3P[ix.k][ix.j+1][ix.i  ];
  emf3D_jp1ip1 = emf3[ix.k][ix.j+1][ix.i+1] - emf3P[ix.k][ix.j+1][ix.i+1];
#endif /* MHD */

/* Use flux differences to correct bad cell-centered quantities */
  
  pG->U[ix.k][ix.j][ix.i].d  += dtodx1*(x1FD_ip1.d  - x1FD_i.d );
  pG->U[ix.k][ix.j][ix.i].M1 += dtodx1*(x1FD_ip1.Mx - x1FD_i.Mx);
  pG->U[ix.k][ix.j][ix.i].M2 += dtodx1*(x1FD_ip1.My - x1FD_i.My);
  pG->U[ix.k][ix.j][ix.i].M3 += dtodx1*(x1FD_ip1.Mz - x1FD_i.Mz);
#ifndef BAROTROPIC
  pG->U[ix.k][ix.j][ix.i].E  += dtodx1*(x1FD_ip1.E  - x1FD_i.E );
#endif /* BAROTROPIC */
  
  pG->U[ix.k][ix.j][ix.i].d  += dtodx2*(x2FD_jp1.d  - x2FD_j.d );
  pG->U[ix.k][ix.j][ix.i].M1 += dtodx2*(x2FD_jp1.Mz - x2FD_j.Mz);
  pG->U[ix.k][ix.j][ix.i].M2 += dtodx2*(x2FD_jp1.Mx - x2FD_j.Mx);
  pG->U[ix.k][ix.j][ix.i].M3 += dtodx2*(x2FD_jp1.My - x2FD_j.My);
#ifndef BAROTROPIC
  pG->U[ix.k][ix.j][ix.i].E  += dtodx2*(x2FD_jp1.E  - x2FD_j.E );
#endif /* BAROTROPIC */

  pG->U[ix.k][ix.j][ix.i].d  += dtodx3*(x3FD_kp1.d  - x3FD_k.d );
  pG->U[ix.k][ix.j][ix.i].M1 += dtodx3*(x3FD_kp1.Mz - x3FD_k.Mz);
  pG->U[ix.k][ix.j][ix.i].M2 += dtodx3*(x3FD_kp1.Mx - x3FD_k.Mx);
  pG->U[ix.k][ix.j][ix.i].M3 += dtodx3*(x3FD_kp1.My - x3FD_k.My);
#ifndef BAROTROPIC
  pG->U[ix.k][ix.j][ix.i].E  += dtodx3*(x3FD_kp1.E  - x3FD_k.E );
#endif /* BAROTROPIC */

/* Use flux differences to correct cell-centered values at i-1 and i+1 */
  
  if (ix.i > pG->is) {
    pG->U[ix.k][ix.j][ix.i-1].d  += dtodx1*(x1FD_i.d );
    pG->U[ix.k][ix.j][ix.i-1].M1 += dtodx1*(x1FD_i.Mx);
    pG->U[ix.k][ix.j][ix.i-1].M2 += dtodx1*(x1FD_i.My);
    pG->U[ix.k][ix.j][ix.i-1].M3 += dtodx1*(x1FD_i.Mz);
#ifndef BAROTROPIC
    pG->U[ix.k][ix.j][ix.i-1].E  += dtodx1*(x1FD_i.E );
#endif /* BAROTROPIC */
  }
  
  if (ix.i < pG->ie) {
    pG->U[ix.k][ix.j][ix.i+1].d  -= dtodx1*(x1FD_ip1.d );
    pG->U[ix.k][ix.j][ix.i+1].M1 -= dtodx1*(x1FD_ip1.Mx);
    pG->U[ix.k][ix.j][ix.i+1].M2 -= dtodx1*(x1FD_ip1.My);
    pG->U[ix.k][ix.j][ix.i+1].M3 -= dtodx1*(x1FD_ip1.Mz);
#ifndef BAROTROPIC
    pG->U[ix.k][ix.j][ix.i+1].E  -= dtodx1*(x1FD_ip1.E );
#endif /* BAROTROPIC */
  }

/* Use flux differences to correct cell-centered values at j-1 and j+1 */
  
  if (ix.j > pG->js) {
    pG->U[ix.k][ix.j-1][ix.i].d  += dtodx2*(x2FD_j.d );
    pG->U[ix.k][ix.j-1][ix.i].M1 += dtodx2*(x2FD_j.Mz);
    pG->U[ix.k][ix.j-1][ix.i].M2 += dtodx2*(x2FD_j.Mx);
    pG->U[ix.k][ix.j-1][ix.i].M3 += dtodx2*(x2FD_j.My);
#ifndef BAROTROPIC
    pG->U[ix.k][ix.j-1][ix.i].E  += dtodx2*(x2FD_j.E );
#endif /* BAROTROPIC */
  }
  
  if (ix.j < pG->je) {
    pG->U[ix.k][ix.j+1][ix.i].d  -= dtodx2*(x2FD_jp1.d );
    pG->U[ix.k][ix.j+1][ix.i].M1 -= dtodx2*(x2FD_jp1.Mz);
    pG->U[ix.k][ix.j+1][ix.i].M2 -= dtodx2*(x2FD_jp1.Mx);
    pG->U[ix.k][ix.j+1][ix.i].M3 -= dtodx2*(x2FD_jp1.My);
#ifndef BAROTROPIC
    pG->U[ix.k][ix.j+1][ix.i].E  -= dtodx2*(x2FD_jp1.E );
#endif /* BAROTROPIC */
  }

/* Use flux differences to correct cell-centered values at k-1 and k+1 */
 
  if (ix.k > pG->ks) {
    pG->U[ix.k-1][ix.j][ix.i].d  += dtodx3*(x3FD_k.d );
    pG->U[ix.k-1][ix.j][ix.i].M1 += dtodx3*(x3FD_k.Mz);
    pG->U[ix.k-1][ix.j][ix.i].M2 += dtodx3*(x3FD_k.Mx);
    pG->U[ix.k-1][ix.j][ix.i].M3 += dtodx3*(x3FD_k.My);
#ifndef BAROTROPIC
    pG->U[ix.k-1][ix.j][ix.i].E  += dtodx3*(x3FD_k.E );
#endif /* BAROTROPIC */
  }
 
  if (ix.k < pG->ke) {
    pG->U[ix.k+1][ix.j][ix.i].d  -= dtodx3*(x3FD_kp1.d );
    pG->U[ix.k+1][ix.j][ix.i].M1 -= dtodx3*(x3FD_kp1.Mz);
    pG->U[ix.k+1][ix.j][ix.i].M2 -= dtodx3*(x3FD_kp1.Mx);
    pG->U[ix.k+1][ix.j][ix.i].M3 -= dtodx3*(x3FD_kp1.My);
#ifndef BAROTROPIC
    pG->U[ix.k+1][ix.j][ix.i].E  -= dtodx3*(x3FD_kp1.E );
#endif /* BAROTROPIC */
  }

/* Use emf differences to correct face-centered fields on edges of bad cell */
#ifdef MHD
  pG->B1i[ix.k][ix.j][ix.i  ] -= dtodx3*(emf2D_kp1i   - emf2D_ki) -
                                 dtodx2*(emf3D_jp1i   - emf3D_ji);
  pG->B1i[ix.k][ix.j][ix.i+1] -= dtodx3*(emf2D_kp1ip1 - emf2D_kip1) -
                                 dtodx2*(emf3D_jp1ip1 - emf3D_jip1);
  pG->B2i[ix.k][ix.j  ][ix.i] -= dtodx1*(emf3D_jip1   - emf3D_ji) -
                                 dtodx3*(emf1D_kp1j   - emf1D_kj);
  pG->B2i[ix.k][ix.j+1][ix.i] -= dtodx1*(emf3D_jp1ip1 - emf3D_jp1i) -
                                 dtodx3*(emf1D_kp1jp1 - emf1D_kjp1);
  pG->B3i[ix.k  ][ix.j][ix.i] -= dtodx2*(emf1D_kjp1   - emf1D_kj) -
                                 dtodx1*(emf2D_kip1   - emf2D_ki);
  pG->B3i[ix.k+1][ix.j][ix.i] -= dtodx2*(emf1D_kp1jp1 - emf1D_kp1j) -
                                 dtodx1*(emf2D_kp1ip1 - emf2D_kp1i);
  
/* Use emf differences to correct face-centered fields around bad cell */
  if (ix.i > pG->is) {
    pG->B2i[ix.k  ][ix.j  ][ix.i-1] += dtodx1*(emf3D_ji);
    pG->B2i[ix.k  ][ix.j+1][ix.i-1] += dtodx1*(emf3D_jp1i);
    pG->B3i[ix.k  ][ix.j  ][ix.i-1] -= dtodx1*(emf2D_ki);
    pG->B3i[ix.k+1][ix.j  ][ix.i-1] -= dtodx1*(emf2D_kp1i);
  }
  if (ix.i < pG->ie) {
    pG->B2i[ix.k  ][ix.j  ][ix.i+1] -= dtodx1*(emf3D_jip1);
    pG->B2i[ix.k  ][ix.j+1][ix.i+1] -= dtodx1*(emf3D_jp1ip1);
    pG->B3i[ix.k  ][ix.j  ][ix.i+1] += dtodx1*(emf2D_kip1);
    pG->B3i[ix.k+1][ix.j  ][ix.i+1] += dtodx1*(emf2D_kp1ip1);
  }
  
  if (ix.j > pG->js) {
    pG->B3i[ix.k  ][ix.j-1][ix.i  ] += dtodx2*(emf1D_kj);
    pG->B3i[ix.k+1][ix.j-1][ix.i  ] += dtodx2*(emf1D_kp1j);
    pG->B1i[ix.k  ][ix.j-1][ix.i  ] -= dtodx2*(emf3D_ji);
    pG->B1i[ix.k  ][ix.j-1][ix.i+1] -= dtodx2*(emf3D_jip1);
  }
  if (ix.j < pG->je) {
    pG->B3i[ix.k  ][ix.j+1][ix.i  ] -= dtodx2*(emf1D_kjp1);
    pG->B3i[ix.k+1][ix.j+1][ix.i  ] -= dtodx2*(emf1D_kp1jp1);
    pG->B1i[ix.k  ][ix.j+1][ix.i  ] += dtodx2*(emf3D_jp1i);
    pG->B1i[ix.k  ][ix.j+1][ix.i+1] += dtodx2*(emf3D_jp1ip1);
  }
  
  if (ix.k > pG->ks) {
    pG->B1i[ix.k-1][ix.j  ][ix.i  ] += dtodx3*(emf2D_ki);
    pG->B1i[ix.k-1][ix.j  ][ix.i+1] += dtodx3*(emf2D_kip1);
    pG->B2i[ix.k-1][ix.j  ][ix.i  ] -= dtodx3*(emf1D_kj);
    pG->B2i[ix.k-1][ix.j+1][ix.i  ] -= dtodx3*(emf1D_kjp1);
  }
  if (ix.k < pG->ke) {
    pG->B1i[ix.k+1][ix.j  ][ix.i  ] -= dtodx3*(emf2D_kp1i);
    pG->B1i[ix.k+1][ix.j  ][ix.i+1] -= dtodx3*(emf2D_kp1ip1);
    pG->B2i[ix.k+1][ix.j  ][ix.i  ] += dtodx3*(emf1D_kp1j);
    pG->B2i[ix.k+1][ix.j+1][ix.i  ] += dtodx3*(emf1D_kp1jp1);
  }
  
/* Compute new cell-centered fields */
  for (k=(ix.k-1); k<=(ix.k+1); k++) {
  for (j=(ix.j-1); j<=(ix.j+1); j++) {
  for (i=(ix.i-1); i<=(ix.i+1); i++) {
    pG->U[k][j][i].B1c = 0.5*(pG->B1i[k][j][i] + pG->B1i[k][j][i+1]);
    pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j+1][i]);
    pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i] + pG->B3i[k+1][j][i]);
  }}}
#endif /* MHD */

/* With SMR, replace higher-order fluxes with predict fluxes in case they are
 * used at fine/coarse grid boundaries */
#ifdef STATIC_MESH_REFINEMENT
  x1Flux[ix.k][ix.j][ix.i] = x1FluxP[ix.k][ix.j][ix.i];
  x2Flux[ix.k][ix.j][ix.i] = x2FluxP[ix.k][ix.j][ix.i];
  x3Flux[ix.k][ix.j][ix.i] = x3FluxP[ix.k][ix.j][ix.i];
#endif /* STATIC_MESH_REFINEMENT */

}
#endif /* FIRST_ORDER_FLUX_CORRECTION */

#endif /* VL_INTEGRATOR */

#endif /* SPECIAL_RELATIVITY */
