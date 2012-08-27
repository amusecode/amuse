#include "../copyright.h"
/*============================================================================*/
/*! \file integrate_1d_vl_sr.c
 *  \brief Integrate SRMHD equations using 1D version of MUSCL-Hancock (VL)
 *   integrator.
 *
 * PURPOSE: Integrate SRMHD equations using 1D version of MUSCL-Hancock (VL)
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

#ifdef SPECIAL_RELATIVITY

#if defined(VL_INTEGRATOR) && defined(CARTESIAN)

/* The L/R states of primitive variables and fluxes at each cell face */
static Prim1DS *Wl_x1Face=NULL, *Wr_x1Face=NULL;
static Cons1DS *x1Flux=NULL;

/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W1d=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL, *Ul=NULL, *Ur=NULL;
#ifdef FIRST_ORDER_FLUX_CORRECTION
static Cons1DS *x1FluxP=NULL;
#endif

/* conserved & primitive variables at t^{n} */
static PrimS *W=NULL;

/* conserved variables at t^{n+1/2} computed in predict step */
static ConsS *Uhalf=NULL;
static PrimS *Whalf=NULL;

#ifdef FIRST_ORDER_FLUX_CORRECTION
static void FixCell(GridS *pG, Int3Vect);
#endif


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
  ConsS U;
  Real dtodx1=pG->dt/pG->dx1, hdtodx1=0.5*pG->dt/pG->dx1;
  int i, is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks;
  int cart_x1 = 1, cart_x2 = 2, cart_x3 = 3;
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
#ifdef FIRST_ORDER_FLUX_CORRECTION
  int flag_cell=0,negd=0,negP=0,superl=0,NaNFlux=0;
  int fail=0,final=0;
  Real Vsq,Bx;
  PrimS Wcheck;
  Int3Vect BadCell;
#endif


  int il=is-(nghost-1), iu=ie+(nghost-1);

  for (i=is-nghost; i<=ie+nghost; i++) {
    Uhalf[i] = pG->U[ks][js][i];
    W[i] = Cons_to_Prim(&(pG->U[ks][js][i]));
  }

/*=== STEP 1: Compute first-order fluxes at t^{n} in x1-direction ============*/
/* No source terms are needed since there is no temporal evolution */

/*--- Step 1a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;  
 * W1d = (d, V1, V2, V3, P, B2c, B3c, s[n])
 */

  for (i=is-nghost; i<=ie+nghost; i++) {
    W1d[i].d  = W[i].d;
    W1d[i].Vx = W[i].V1;
    W1d[i].Vy = W[i].V2;
    W1d[i].Vz = W[i].V3;
    W1d[i].P  = W[i].P;
#ifdef MHD
    W1d[i].By = W[i].B2c;
    W1d[i].Bz = W[i].B3c;
    Bxc[i] = W[i].B1c;
    Bxi[i] = pG->B1i[ks][js][i];
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) W1d[i].s[n] = W[i].s[n];
#endif
  }

/*--- Step 1b ------------------------------------------------------------------
 * Compute first-order L/R states */

/* Ensure that W & U are consistent */
  for (i=is-nghost; i<=ie+nghost; i++) {
    U1d[i] = Prim1D_to_Cons1D(&W1d[i],&Bxc[i]);
  }

  for (i=il; i<=ie+nghost; i++) {
    Wl[i] = W1d[i-1];
    Wr[i] = W1d[i  ];

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
    Uhalf[i].E   -= hdtodx1*(x1Flux[i+1].E  - x1Flux[i].E );
#ifdef MHD
    Uhalf[i].B2c -= hdtodx1*(x1Flux[i+1].By - x1Flux[i].By);
    Uhalf[i].B3c -= hdtodx1*(x1Flux[i+1].Bz - x1Flux[i].Bz);
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++)
      Uhalf[i].s[n] -= hdtodx1*(x1Flux[i+1].s[n] - x1Flux[i].s[n]);
#endif
#ifdef FIRST_ORDER_FLUX_CORRECTION
          x1FluxP[i] = x1Flux[i];
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
      Uhalf[i].E -= hdtodx1*(x1Flux[i  ].d*(phic - phil) +
                             x1Flux[i+1].d*(phir - phic));
    }
  }

/*=== STEP 7: Conserved->Primitive variable inversion at t^{n+1/2} ===========*/
        
/* Invert conserved variables at t^{n+1/2} to primitive variables. With FOFC, 
 * if cell-centered d < 0, P< 0, or v^2 > 1, correct by switching back to 
 * values at beginning of step, rendering update first order in time for that
 * cell.
 */

#ifdef FIRST_ORDER_FLUX_CORRECTION        
  negd = 0;
  negP = 0;
  superl = 0;
  flag_cell = 0;
#endif
  for (i=il; i<=iu; i++) {
    Whalf[i] = Cons_to_Prim(&Uhalf[i]);
#ifdef FIRST_ORDER_FLUX_CORRECTION   
    if (Whalf[i].d < 0.0) {
      flag_cell = 1;
      negd++;
    }
    if (Whalf[i].P < 0.0) {
      flag_cell = 1;
      negP++;
    }
    Vsq = SQR(Whalf[i].V1) + SQR(Whalf[i].V2) + SQR(Whalf[i].V3);
    if (Vsq > 1.0) {
      flag_cell = 1;
      superl++;
    }
    if (flag_cell != 0) {
      Whalf[i].d = W[i].d;
      Whalf[i].V1 = W[i].V1;
      Whalf[i].V2 = W[i].V2;
      Whalf[i].V3 = W[i].V3;
      Whalf[i].P = W[i].P;
      flag_cell=0;
    }
#endif
  }

#ifdef FIRST_ORDER_FLUX_CORRECTION  
  if (negd > 0 || negP > 0 || superl > 0)
    printf("[Step7]: %i cells had d<0; %i cells had P<0; %i cells had v>1\n"
                                 ,negd,negP,superl);
#endif


/*=== STEP 8: Compute second-order L/R x1-interface states ===================*/

/*--- Step 8a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;
 * W = (d, V1, V2, V3, P, B2c, B3c, s[n])
 */

  for (i=il; i<=iu; i++) {
    W1d[i].d  = Whalf[i].d;
    W1d[i].Vx = Whalf[i].V1;
    W1d[i].Vy = Whalf[i].V2;
    W1d[i].Vz = Whalf[i].V3;
    W1d[i].P  = Whalf[i].P;
#ifdef MHD
    W1d[i].By = Whalf[i].B2c;
    W1d[i].Bz = Whalf[i].B3c;
    Bxc[i]  = Whalf[i].B1c;
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) W1d[i].s[n] = Whalf[i].s[n];
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
    Wl_x1Face[i] = Wl[i];
    Wr_x1Face[i] = Wr[i];
  }

/*=== STEPS 9-10: Not needed in 1D ===*/

/*=== STEP 11: Compute x1-Flux ===============================================*/

/*--- Step 11b -----------------------------------------------------------------
 * Compute second-order fluxes in x1-direction
 */

  for (i=is; i<=ie+1; i++) {
    Ul[i] = Prim1D_to_Cons1D(&Wl_x1Face[i],&Bxi[i]);
    Ur[i] = Prim1D_to_Cons1D(&Wr_x1Face[i],&Bxi[i]);

    fluxes(Ul[i],Ur[i],Wl_x1Face[i],Wr_x1Face[i],Bxi[i],&x1Flux[i]);
  }

/*=== STEP 12: Not needed in 1D ===*/
        
/*=== STEP 13: Add source terms for a full timestep using n+1/2 states =======*/
       
/*--- Step 13a -----------------------------------------------------------------
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

/*=== STEP 14: Update cell-centered values for a full timestep ===============*/

/*--- Step 14a -----------------------------------------------------------------
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

#ifdef FIRST_ORDER_FLUX_CORRECTION
/*=== STEP 15: First-order flux correction ===================================*/

/*--- Step 15a -----------------------------------------------------------------
 * If cell-centered d or P have gone negative, or if v^2 > 1, correct
 * by using 1st order predictor fluxes */
        
  for (i=is; i<=ie; i++) {
      Wcheck = check_Prim(&(pG->U[ks][js][i]));
      if (Wcheck.d < 0.0) {
        flag_cell = 1;
        BadCell.i = i;
        BadCell.j = js;
        BadCell.k = ks;
        negd++;
      }
      if (Wcheck.P < 0.0) {
        flag_cell = 1;
        BadCell.i = i;
        BadCell.j = js;
        BadCell.k = ks;
        negP++;
      }
      Vsq = SQR(Wcheck.V1) + SQR(Wcheck.V2) + SQR(Wcheck.V3);
      if (Vsq > 1.0) {
        flag_cell = 1;
        BadCell.i = i;
        BadCell.j = js;
        BadCell.k = ks;
        superl++;
      }
      if (flag_cell != 0) {
        FixCell(pG, BadCell);
        flag_cell=0;
      }

  }

  if (negd > 0 || negP > 0 || superl > 0){
    printf("[Step15a]: %i cells had d<0; %i cells had P<0;\n",negd,negP);
    printf("[Step15a]: %i cells had v>1 at 1st correction\n",superl);
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
  fail = 0;
  negd = 0;
  negP = 0;
  final = 0;
  superl = 0;
  for (i=is; i<=ie; i++) {
    flag_cell=0;
    Wcheck = check_Prim(&(pG->U[ks][js][i]));
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
    if (flag_cell != 0) {
      final++;
      Wcheck = fix_vsq (&(pG->U[ks][js][i]));
      U = Prim_to_Cons(&Wcheck);
      pG->U[ks][js][i].d = U.d;
      pG->U[ks][js][i].M1 = U.M1;
      pG->U[ks][js][i].M2 = U.M2;
      pG->U[ks][js][i].M3 = U.M3;
      pG->U[ks][js][i].E = U.E;
      Wcheck = check_Prim(&(pG->U[ks][js][i]));
      Vsq = SQR(Wcheck.V1) + SQR(Wcheck.V2) + SQR(Wcheck.V3);
      if (Wcheck.d < 0.0 || Wcheck.P < 0.0 || Vsq > 1.0){
	fail++;
      }
    }
  }

  if (negd > 0 || negP > 0 || superl > 0) {
    printf("[Step15b]: %i cells had d<0; %i cells had P<0;\n",negd,negP);
    printf("[Step15b]: %i cells had v>1	after 1st order correction\n",superl);
    printf("[Step15b]: %i cells required a  final fix\n",final);
    printf("[Step15b]: %i cells had an unphysical state\n",fail);
  }
#endif /* MHD */
#endif /* FIRST_ORDER_FLUX_CORRECTION */



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
#ifdef FIRST_ORDER_FLUX_CORRECTION
  if ((x1FluxP  =(Cons1DS*)malloc(size1*sizeof(Cons1DS))) ==NULL) goto on_error;
#endif /* FIRST_ORDER_FLUX_CORRECTION */

  if ((Bxc = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;

  if ((U1d= (Cons1DS*)malloc(size1*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ul = (Cons1DS*)malloc(size1*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ur = (Cons1DS*)malloc(size1*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W1d= (Prim1DS*)malloc(size1*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wl = (Prim1DS*)malloc(size1*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wr = (Prim1DS*)malloc(size1*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((Uhalf = (ConsS*)malloc(size1*sizeof(ConsS)))==NULL) goto on_error;
  if ((Whalf = (PrimS*)malloc(size1*sizeof(PrimS)))==NULL) goto on_error;
  if ((W     = (PrimS*)malloc(size1*sizeof(PrimS)))==NULL) goto on_error;

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
#ifdef FIRST_ORDER_FLUX_CORRECTION
  if (x1FluxP   != NULL) free(x1FluxP);
#endif


  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);

  if (U1d      != NULL) free(U1d);
  if (Ul       != NULL) free(Ul);
  if (Ur       != NULL) free(Ur);
  if (W1d      != NULL) free(W1d);
  if (Wl       != NULL) free(Wl);
  if (Wr       != NULL) free(Wr);

  if (Uhalf    != NULL) free(Uhalf);
  if (Whalf    != NULL) free(Whalf);
  if (W        != NULL) free(W);

  return;
}

#ifdef FIRST_ORDER_FLUX_CORRECTION
/*----------------------------------------------------------------------------*/
/*! \fn static void FixCell(GridS *pG, Int3Vect indx)
 *  \brief Uses first order fluxes to fix negative d,P or superluminal v
 */ 
static void FixCell(GridS *pG, Int3Vect indx)
{
  int ks=pG->ks,js=pG->js;
  Real dtodx1=pG->dt/pG->dx1;
  Cons1DS x1FD_i, x1FD_ip1, x2FD_j, x2FD_jp1;
  
  /* Compute difference of predictor and corrector fluxes at cell faces */
  
  x1FD_i.d    = x1Flux[indx.i  ].d - x1FluxP[indx.i  ].d;
  x1FD_ip1.d  = x1Flux[indx.i+1].d - x1FluxP[indx.i+1].d;
        
  x1FD_i.Mx    = x1Flux[indx.i  ].Mx - x1FluxP[indx.i  ].Mx;
  x1FD_ip1.Mx  = x1Flux[indx.i+1].Mx - x1FluxP[indx.i+1].Mx;
        
  x1FD_i.My    = x1Flux[indx.i  ].My - x1FluxP[indx.i  ].My;
  x1FD_ip1.My  = x1Flux[indx.i+1].My - x1FluxP[indx.i+1].My;
        
  x1FD_i.Mz    = x1Flux[indx.i  ].Mz - x1FluxP[indx.i  ].Mz;
  x1FD_ip1.Mz  = x1Flux[indx.i+1].Mz - x1FluxP[indx.i+1].Mz;
        
#ifdef MHD
  x1FD_i.By    = x1Flux[indx.i  ].By - x1FluxP[indx.i  ].By;
  x1FD_ip1.By  = x1Flux[indx.i+1].By - x1FluxP[indx.i+1].By;
        
  x1FD_i.Bz    = x1Flux[indx.i  ].Bz - x1FluxP[indx.i  ].Bz;
  x1FD_ip1.Bz  = x1Flux[indx.i+1].Bz - x1FluxP[indx.i+1].Bz;
#endif
        
#ifndef BAROTROPIC
  x1FD_i.E    = x1Flux[indx.i  ].E - x1FluxP[indx.i  ].E;
  x1FD_ip1.E  = x1Flux[indx.i+1].E - x1FluxP[indx.i+1].E;
#endif /* BAROTROPIC */
        
  /* Use flux differences to correct bad cell */
  pG->U[ks][js][indx.i].d  += dtodx1*(x1FD_ip1.d  - x1FD_i.d );
  pG->U[ks][js][indx.i].M1 += dtodx1*(x1FD_ip1.Mx - x1FD_i.Mx);
  pG->U[ks][js][indx.i].M2 += dtodx1*(x1FD_ip1.My - x1FD_i.My);
  pG->U[ks][js][indx.i].M3 += dtodx1*(x1FD_ip1.Mz - x1FD_i.Mz);
#ifdef MHD
  pG->U[ks][js][indx.i].B2c += dtodx1*(x1FD_ip1.By - x1FD_i.By);
  pG->U[ks][js][indx.i].B3c += dtodx1*(x1FD_ip1.Bz - x1FD_i.Bz);
  /* For consistency, set B2i and B3i to cell-centered values.  */
  pG->B2i[ks][js][indx.i] = pG->U[ks][js][indx.i].B2c;
  pG->B3i[ks][js][indx.i] = pG->U[ks][js][indx.i].B3c;
#endif
#ifndef BAROTROPIC
  pG->U[ks][js][indx.i].E  += dtodx1*(x1FD_ip1.E  - x1FD_i.E );
#endif /* BAROTROPIC */
        
        
  /* Use flux differences to correct bad cell neighbors at i-1 and i+1 */      
  if (indx.i > pG->is) {
    pG->U[ks][js][indx.i-1].d  += dtodx1*(x1FD_i.d );
    pG->U[ks][js][indx.i-1].M1 += dtodx1*(x1FD_i.Mx);
    pG->U[ks][js][indx.i-1].M2 += dtodx1*(x1FD_i.My);
    pG->U[ks][js][indx.i-1].M3 += dtodx1*(x1FD_i.Mz);
#ifdef MHD
    pG->U[ks][js][indx.i-1].B2c += dtodx1*(x1FD_i.By);
    pG->U[ks][js][indx.i-1].B3c += dtodx1*(x1FD_i.Bz);
    /* For consistency, set B2i and B3i to cell-centered values.  */
    pG->B2i[ks][js][indx.i-1] = pG->U[ks][js][indx.i-1].B2c;
    pG->B3i[ks][js][indx.i-1] = pG->U[ks][js][indx.i-1].B3c;
#endif
#ifndef BAROTROPIC
    pG->U[ks][js][indx.i-1].E  += dtodx1*(x1FD_i.E );
#endif /* BAROTROPIC */
  }
        
  if (indx.i < pG->ie) {
    pG->U[ks][js][indx.i+1].d  -= dtodx1*(x1FD_ip1.d );
    pG->U[ks][js][indx.i+1].M1 -= dtodx1*(x1FD_ip1.Mx);
    pG->U[ks][js][indx.i+1].M2 -= dtodx1*(x1FD_ip1.My);
    pG->U[ks][js][indx.i+1].M3 -= dtodx1*(x1FD_ip1.Mz);
#ifdef MHD
    pG->U[ks][js][indx.i+1].B2c -= dtodx1*(x1FD_ip1.By);
    pG->U[ks][js][indx.i+1].B3c -= dtodx1*(x1FD_ip1.Bz);
    /* For consistency, set B2i and B3i to cell-centered values.  */
    pG->B2i[ks][js][indx.i+1] = pG->U[ks][js][indx.i+1].B2c;
    pG->B3i[ks][js][indx.i+1] = pG->U[ks][js][indx.i+1].B3c;
#endif MHD
#ifndef BAROTROPIC
    pG->U[ks][js][indx.i+1].E  -= dtodx1*(x1FD_ip1.E );
#endif /* BAROTROPIC */
  }
        
}
#endif /* FIRST_ORDER_FLUX_CORRECTION */

#endif /* VL_INTEGRATOR */

#endif /* SPECIAL_RELATIVITY */
