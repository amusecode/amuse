#include "../copyright.h"
/*============================================================================*/
/*! \file integrate_1d_vl.c
 *  \brief Integrate MHD equations using 1D version of MUSCL-Hancock (VL)
 *   integrator.
 *
 * PURPOSE: Integrate MHD equations using 1D version of MUSCL-Hancock (VL)
 *   integrator.  Updates U.[d,M1,M2,M3,E,B2c,B3c,s] in Grid structure.
 *   Adds gravitational source terms, self-gravity.
 *        
 * CONTAINS PUBLIC FUNCTIONS: 
 * - integrate_1d_vl()
 * - integrate_init_1d()
 * - integrate_destruct_1d() */
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

#ifndef SPECIAL_RELATIVITY

#if defined(VL_INTEGRATOR) && defined(CARTESIAN)

/* The L/R states of primitive variables and fluxes at each cell face */
static Prim1DS *Wl_x1Face=NULL, *Wr_x1Face=NULL;
static Cons1DS *x1Flux=NULL;

/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL, *Ul=NULL, *Ur=NULL;

/* conserved variables at t^{n+1/2} computed in predict step */
static ConsS *Uhalf=NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void integrate_1d_vl(DomainS *pD)
 *  \brief 1D version of van Leer unsplit integrator for MHD.  
 *
 *   The numbering of steps follows the numbering in the 3D version.
 *   NOT ALL STEPS ARE NEEDED IN 1D.
 */

void integrate_1d_vl(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  Real dtodx1=pG->dt/pG->dx1, hdtodx1=0.5*pG->dt/pG->dx1;
  int i, is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks;
  Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic;
#if (NSCALARS > 0)
  int n;
#endif
#ifdef SELF_GRAVITY
  Real gxl,gxr,flx_m1l,flx_m1r;
#endif
#ifdef STATIC_MESH_REFINEMENT
  int ncg,npg,dim;
#endif

  int il=is-(nghost-1), iu=ie+(nghost-1);

  for (i=is-nghost; i<=ie+nghost; i++) {
    Uhalf[i] = pG->U[ks][js][i];
  }

/*=== STEP 1: Compute first-order fluxes at t^{n} in x1-direction ============*/
/* No source terms are needed since there is no temporal evolution */

/*--- Step 1a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;  
 * U1d = (d, M1, M2, M3, E, B2c, B3c, s[n])
 */

  for (i=is-nghost; i<=ie+nghost; i++) {
    U1d[i].d  = pG->U[ks][js][i].d;
    U1d[i].Mx = pG->U[ks][js][i].M1;
    U1d[i].My = pG->U[ks][js][i].M2;
    U1d[i].Mz = pG->U[ks][js][i].M3;
#ifndef BAROTROPIC
    U1d[i].E  = pG->U[ks][js][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
    U1d[i].By = pG->U[ks][js][i].B2c;
    U1d[i].Bz = pG->U[ks][js][i].B3c;
    Bxc[i] = pG->U[ks][js][i].B1c;
    Bxi[i] = pG->B1i[ks][js][i];
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) U1d[i].s[n] = pG->U[ks][js][i].s[n];
#endif
  }

/*--- Step 1b ------------------------------------------------------------------
 * Compute first-order L/R states */

  for (i=is-nghost; i<=ie+nghost; i++) {
    W[i] = Cons1D_to_Prim1D(&U1d[i],&Bxc[i]);
  }

  for (i=il; i<=ie+nghost; i++) {
    Wl[i] = W[i-1];
    Wr[i] = W[i  ];

    Ul[i] = U1d[i-1];
    Ur[i] = U1d[i  ];
  }

/*--- Step 1c ------------------------------------------------------------------
 * No source terms needed */

/*--- Step 1d ------------------------------------------------------------------
 * Compute flux in x1-direction */

  for (i=il; i<=ie+nghost; i++) {
    fluxes(Ul[i],Ur[i],Wl[i],Wr[i],Bxi[i],&x1Flux[i]);
  }

/*=== STEPS 2-4: Not needed in 1D ===*/

/*=== STEP 5: Update cell-centered variables to half-timestep ================*/

/*--- Step 5a ------------------------------------------------------------------
 * Update cell-centered variables (including B2c and B3c) to half-timestep
 */

  for (i=il; i<=iu; i++) {
    Uhalf[i].d   -= hdtodx1*(x1Flux[i+1].d  - x1Flux[i].d );
    Uhalf[i].M1  -= hdtodx1*(x1Flux[i+1].Mx - x1Flux[i].Mx);
    Uhalf[i].M2  -= hdtodx1*(x1Flux[i+1].My - x1Flux[i].My);
    Uhalf[i].M3  -= hdtodx1*(x1Flux[i+1].Mz - x1Flux[i].Mz);
#ifndef BAROTROPIC
    Uhalf[i].E   -= hdtodx1*(x1Flux[i+1].E  - x1Flux[i].E );
#endif /* BAROTROPIC */
#ifdef MHD
    Uhalf[i].B2c -= hdtodx1*(x1Flux[i+1].By - x1Flux[i].By);
    Uhalf[i].B3c -= hdtodx1*(x1Flux[i+1].Bz - x1Flux[i].Bz);
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++)
      Uhalf[i].s[n] -= hdtodx1*(x1Flux[i+1].s[n] - x1Flux[i].s[n]);
#endif
  }

/*=== STEP 6: Add source terms to predict values at half-timestep ============*/

/*--- Step 6a ------------------------------------------------------------------
 * Add source terms from a static gravitational potential for 0.5*dt to predict
 * step.  To improve conservation of total energy, we average the energy
 * source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (i=il; i<=iu; i++) {
      cc_pos(pG,i,js,ks,&x1,&x2,&x3);
      phic = (*StaticGravPot)((x1            ),x2,x3);
      phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
      phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

      Uhalf[i].M1 -= hdtodx1*pG->U[ks][js][i].d*(phir-phil);
#ifndef BAROTROPIC
      Uhalf[i].E -= hdtodx1*(x1Flux[i  ].d*(phic - phil) +
                             x1Flux[i+1].d*(phir - phic));
#endif
    }
  }

/*--- Step 6b ------------------------------------------------------------------
 * Add source terms for self gravity for 0.5*dt to predict step.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

#ifdef SELF_GRAVITY
  for (i=il; i<=iu; i++) {
    phic = pG->Phi[ks][js][i];
    phir = 0.5*(pG->Phi[ks][js][i] + pG->Phi[ks][js][i+1]);
    phil = 0.5*(pG->Phi[ks][js][i] + pG->Phi[ks][js][i-1]);

    Uhalf[i].M1 -= hdtodx1*pG->U[ks][js][i].d*(phir-phil);
#ifndef BAROTROPIC
    Uhalf[i].E -= hdtodx1*(x1Flux[i  ].d*(phic - phil) +
                           x1Flux[i+1].d*(phir - phic));
#endif
  }
#endif /* SELF_GRAVITY */

/*=== STEP 7: Compute second-order L/R x1-interface states ===================*/

/*--- Step 7a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U = (d, M1, M2, M3, E, B2c, B3c, s[n])
 */

  for (i=il; i<=iu; i++) {
    U1d[i].d  = Uhalf[i].d;
    U1d[i].Mx = Uhalf[i].M1;
    U1d[i].My = Uhalf[i].M2;
    U1d[i].Mz = Uhalf[i].M3;
#ifndef BAROTROPIC
    U1d[i].E  = Uhalf[i].E;
#endif /* BAROTROPIC */
#ifdef MHD
    U1d[i].By = Uhalf[i].B2c;
    U1d[i].Bz = Uhalf[i].B3c;
    Bxc[i]  = Uhalf[i].B1c;
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) U1d[i].s[n] = Uhalf[i].s[n];
#endif /* NSCALARS */
  }

/*--- Step 7b ------------------------------------------------------------------
 * Compute L/R states on x1-interfaces, store into arrays
 */

  for (i=il; i<=iu; i++) {
    W[i] = Cons1D_to_Prim1D(&U1d[i],&Bxc[i]);
  }

  lr_states(pG,W,Bxc,pG->dt,pG->dx1,is,ie,Wl,Wr,1);

  for (i=is; i<=ie+1; i++) {
    Wl_x1Face[i] = Wl[i];
    Wr_x1Face[i] = Wr[i];
  }

/*=== STEPS 8-9: Not needed in 1D ===*/

/*=== STEP 10: Compute x1-Flux ===============================================*/

/*--- Step 10b -----------------------------------------------------------------
 * Compute second-order fluxes in x1-direction
 */

  for (i=is; i<=ie+1; i++) {
    Ul[i] = Prim1D_to_Cons1D(&Wl_x1Face[i],&Bxi[i]);
    Ur[i] = Prim1D_to_Cons1D(&Wr_x1Face[i],&Bxi[i]);

    fluxes(Ul[i],Ur[i],Wl_x1Face[i],Wr_x1Face[i],Bxi[i],&x1Flux[i]);
  }

/*=== STEP 11: Not needed in 1D ===*/
        
/*=== STEP 12: Add source terms for a full timestep using n+1/2 states =======*/
       
/*--- Step 12a -----------------------------------------------------------------
 * Add gravitational source terms due to a Static Potential
 * To improve conservation of total energy, we average the energy
 * source term computed at cell faces.
 *    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (i=is; i<=ie; i++) {
      cc_pos(pG,i,js,ks,&x1,&x2,&x3);
      phic = (*StaticGravPot)((x1            ),x2,x3);
      phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
      phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

      pG->U[ks][js][i].M1 -= dtodx1*Uhalf[i].d*(phir-phil);
#ifndef BAROTROPIC
      pG->U[ks][js][i].E -= dtodx1*(x1Flux[i  ].d*(phic - phil) +
                                    x1Flux[i+1].d*(phir - phic));
#endif
    }
  }

/*--- Step 12b -----------------------------------------------------------------
 * Add gravitational source terms for self-gravity.
 * A flux correction using Phi^{n+1} in the main loop is required to make
 * the source terms 2nd order: see selfg_flux_correction().
 */

#ifdef SELF_GRAVITY
/* Add fluxes and source terms due to (d/dx1) terms  */

  for (i=is; i<=ie; i++){
    phic = pG->Phi[ks][js][i];
    phil = 0.5*(pG->Phi[ks][js][i-1] + pG->Phi[ks][js][i  ]);
    phir = 0.5*(pG->Phi[ks][js][i  ] + pG->Phi[ks][js][i+1]);

/* gx, gy and gz centered at L and R x1-faces */
    gxl = (pG->Phi[ks][js][i-1] - pG->Phi[ks][js][i  ])/(pG->dx1);
    gxr = (pG->Phi[ks][js][i  ] - pG->Phi[ks][js][i+1])/(pG->dx1);

/* momentum fluxes in x1.  2nd term is needed only if Jean's swindle used */
    flx_m1l = 0.5*(gxl*gxl)/four_pi_G + grav_mean_rho*phil;
    flx_m1r = 0.5*(gxr*gxr)/four_pi_G + grav_mean_rho*phir;

/* Update momenta and energy with d/dx1 terms  */
    pG->U[ks][js][i].M1 -= dtodx1*(flx_m1r - flx_m1l);
#ifndef BAROTROPIC
    pG->U[ks][js][i].E -= dtodx1*(x1Flux[i  ].d*(phic - phil) +
                                  x1Flux[i+1].d*(phir - phic));
#endif /* BAROTROPIC */
  }

/* Save mass fluxes in Grid structure for source term correction in main loop */

  for (i=is; i<=ie+1; i++) {
    pG->x1MassFlux[ks][js][i] = x1Flux[i].d;
  }
#endif /* SELF_GRAVITY */

/*=== STEP 13: Update cell-centered values for a full timestep ===============*/

/*--- Step 13a -----------------------------------------------------------------
 * Update cell-centered variables in pG (including B2c and B3c) using x1-Fluxes
 */

  for (i=is; i<=ie; i++) {
    pG->U[ks][js][i].d  -= dtodx1*(x1Flux[i+1].d  - x1Flux[i].d );
    pG->U[ks][js][i].M1 -= dtodx1*(x1Flux[i+1].Mx - x1Flux[i].Mx);
    pG->U[ks][js][i].M2 -= dtodx1*(x1Flux[i+1].My - x1Flux[i].My);
    pG->U[ks][js][i].M3 -= dtodx1*(x1Flux[i+1].Mz - x1Flux[i].Mz);
#ifndef BAROTROPIC
    pG->U[ks][js][i].E  -= dtodx1*(x1Flux[i+1].E  - x1Flux[i].E );
#endif /* BAROTROPIC */
#ifdef MHD
    pG->U[ks][js][i].B2c -= dtodx1*(x1Flux[i+1].By - x1Flux[i].By);
    pG->U[ks][js][i].B3c -= dtodx1*(x1Flux[i+1].Bz - x1Flux[i].Bz);
/* For consistency, set B2i and B3i to cell-centered values.  */
    pG->B2i[ks][js][i] = pG->U[ks][js][i].B2c;
    pG->B3i[ks][js][i] = pG->U[ks][js][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++)
      pG->U[ks][js][i].s[n] -= dtodx1*(x1Flux[i+1].s[n] - x1Flux[i].s[n]);
#endif
  }

#ifdef STATIC_MESH_REFINEMENT
/*--- Step 13d -----------------------------------------------------------------
 * With SMR, store fluxes at boundaries of child and parent grids.
 */

/* x1-boundaries of child Grids (interior to THIS Grid) */
  for (ncg=0; ncg<pG->NCGrid; ncg++) {
    for (dim=0; dim<2; dim++){
      if (pG->CGrid[ncg].myFlx[dim] != NULL) {

        if (dim==0) i = pG->CGrid[ncg].ijks[0];
        if (dim==1) i = pG->CGrid[ncg].ijke[0] + 1;

        pG->CGrid[ncg].myFlx[dim][ks][js].d  = x1Flux[i].d;
        pG->CGrid[ncg].myFlx[dim][ks][js].M1 = x1Flux[i].Mx;
        pG->CGrid[ncg].myFlx[dim][ks][js].M2 = x1Flux[i].My;
        pG->CGrid[ncg].myFlx[dim][ks][js].M3 = x1Flux[i].Mz;
#ifndef BAROTROPIC
        pG->CGrid[ncg].myFlx[dim][ks][js].E  = x1Flux[i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        pG->CGrid[ncg].myFlx[dim][ks][js].B1c = 0.0;
        pG->CGrid[ncg].myFlx[dim][ks][js].B2c = x1Flux[i].By;
        pG->CGrid[ncg].myFlx[dim][ks][js].B3c = x1Flux[i].Bz;
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          pG->CGrid[ncg].myFlx[dim][ks][js].s[n]  = x1Flux[i].s[n];
#endif
      }
    }
  }

/* x1-boundaries of parent Grids (at boundaries of THIS Grid)  */
  for (npg=0; npg<pG->NPGrid; npg++) {
    for (dim=0; dim<2; dim++){
      if (pG->PGrid[npg].myFlx[dim] != NULL) {

        if (dim==0) i = pG->PGrid[npg].ijks[0];
        if (dim==1) i = pG->PGrid[npg].ijke[0] + 1;

        pG->PGrid[npg].myFlx[dim][ks][js].d  = x1Flux[i].d;
        pG->PGrid[npg].myFlx[dim][ks][js].M1 = x1Flux[i].Mx;
        pG->PGrid[npg].myFlx[dim][ks][js].M2 = x1Flux[i].My;
        pG->PGrid[npg].myFlx[dim][ks][js].M3 = x1Flux[i].Mz;
#ifndef BAROTROPIC
        pG->PGrid[npg].myFlx[dim][ks][js].E  = x1Flux[i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        pG->PGrid[npg].myFlx[dim][ks][js].B1c = 0.0;
        pG->PGrid[npg].myFlx[dim][ks][js].B2c = x1Flux[i].By;
        pG->PGrid[npg].myFlx[dim][ks][js].B3c = x1Flux[i].Bz;
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          pG->PGrid[npg].myFlx[dim][ks][js].s[n]  = x1Flux[i].s[n];
#endif
      }
    }
  }

#endif /* STATIC_MESH_REFINEMENT */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_init_1d(MeshS *pM)
 *  \brief Allocate temporary integration arrays */
void integrate_init_1d(MeshS *pM)
{
  int size1=0,nl,nd;

/* Cycle over all Grids on this processor to find maximum Nx1 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if (pM->Domain[nl][nd].Grid->Nx[0] > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
      }
    }
  }

  size1 = size1 + 2*nghost;

  if ((Wl_x1Face=(Prim1DS*)malloc(size1*sizeof(Prim1DS))) ==NULL) goto on_error;
  if ((Wr_x1Face=(Prim1DS*)malloc(size1*sizeof(Prim1DS))) ==NULL) goto on_error;
  if ((x1Flux   =(Cons1DS*)malloc(size1*sizeof(Cons1DS))) ==NULL) goto on_error;

  if ((Bxc = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;

  if ((U1d= (Cons1DS*)malloc(size1*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ul = (Cons1DS*)malloc(size1*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ur = (Cons1DS*)malloc(size1*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W  = (Prim1DS*)malloc(size1*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wl = (Prim1DS*)malloc(size1*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wr = (Prim1DS*)malloc(size1*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((Uhalf = (ConsS*)malloc(size1*sizeof(ConsS)))==NULL) goto on_error;

  return;

  on_error:
    integrate_destruct();
    ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_destruct_1d(void)
 *  \brief Free temporary integration arrays */
void integrate_destruct_1d(void)
{
  if (Wl_x1Face != NULL) free(Wl_x1Face);
  if (Wr_x1Face != NULL) free(Wr_x1Face);
  if (x1Flux    != NULL) free(x1Flux);

  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);

  if (U1d      != NULL) free(U1d);
  if (Ul       != NULL) free(Ul);
  if (Ur       != NULL) free(Ur);
  if (W        != NULL) free(W);
  if (Wl       != NULL) free(Wl);
  if (Wr       != NULL) free(Wr);

  if (Uhalf    != NULL) free(Uhalf);

  return;
}
#endif /* VL_INTEGRATOR */

#endif /* SPECIAL_RELATIVITY */
