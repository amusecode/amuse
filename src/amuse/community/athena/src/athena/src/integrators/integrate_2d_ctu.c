#include "../copyright.h"
/*============================================================================*/
/*! \file integrate_2d_ctu.c
 *  \brief Integrate MHD equations in 2D using the directionally unsplit CTU
 *   method of Colella (1990). 
 *
 * PURPOSE: Integrate MHD equations in 2D using the directionally unsplit CTU
 *   method of Colella (1990).  The variables updated are:
 *   -  U.[d,M1,M2,M3,E,B1c,B2c,B3c,s] -- where U is of type ConsS
 *   -  B1i, B2i -- interface magnetic field
 *   Also adds gravitational source terms, self-gravity, optically-thin cooling,
 *   and H-correction of Sanders et al.
 *
 * REFERENCES:
 * - P. Colella, "Multidimensional upwind methods for hyperbolic conservation
 *   laws", JCP, 87, 171 (1990)
 *
 * - T. Gardiner & J.M. Stone, "An unsplit Godunov method for ideal MHD via
 *   constrained transport", JCP, 205, 509 (2005)
 *
 * - R. Sanders, E. Morano, & M.-C. Druguet, "Multidimensinal dissipation for
 *   upwind schemes: stability and applications to gas dynamics", JCP, 145, 511
 *   (1998)
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - integrate_2d_ctu()
 * - integrate_init_2d()
 * - integrate_destruct_2d() */
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
#ifdef PARTICLES
#include "../particles/particle.h"
#endif

#if defined(CTU_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The CTU integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */

/* The L/R states of conserved variables and fluxes at each cell face */
static Cons1DS **Ul_x1Face=NULL, **Ur_x1Face=NULL;
static Cons1DS **Ul_x2Face=NULL, **Ur_x2Face=NULL;
static Cons1DS **x1Flux=NULL, **x2Flux=NULL;

/* The interface magnetic fields and emfs */
#ifdef MHD
static Real **B1_x1Face=NULL, **B2_x2Face=NULL;
static Real **emf3=NULL, **emf3_cc=NULL;
#endif /* MHD */

/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL, *Ul=NULL, *Ur=NULL;

/* density and Pressure at t^{n+1/2} needed by MHD, cooling, and gravity */
static Real **dhalf = NULL,**phalf = NULL;

/* variables needed for H-correction of Sanders et al (1998) */
extern Real etah;
#ifdef H_CORRECTION
static Real **eta1=NULL, **eta2=NULL;
#endif

/* VARIABLES NEEDED FOR CYLINDRICAL COORDINATES */
#ifdef CYLINDRICAL
static Real **geom_src=NULL;
#endif

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   integrate_emf3_corner() - the upwind CT method in Gardiner & Stone (2005) 
 *============================================================================*/

#ifdef MHD
static void integrate_emf3_corner(GridS *pG);
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void integrate_2d_ctu(DomainS *pD)
 *  \brief CTU integrator in 2D.
 *
 *   The numbering of steps follows the numbering in the 3D version.
 *   NOT ALL STEPS ARE NEEDED IN 2D.
 */

void integrate_2d_ctu(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  Real dtodx1 = pG->dt/pG->dx1, dtodx2 = pG->dt/pG->dx2;
  Real hdtodx1 = 0.5*dtodx1, hdtodx2 = 0.5*dtodx2;
  Real hdt = 0.5*pG->dt, dx2 = pG->dx2;
  int i,il,iu,is=pG->is, ie=pG->ie;
  int j,jl,ju,js=pG->js, je=pG->je;
  int ks=pG->ks;
  Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic,M1h,M2h,M3h,Bx=0.0;
#ifndef BAROTROPIC
  Real coolfl,coolfr,coolf,Eh=0.0;
#endif
#ifdef MHD
  Real MHD_src,dbx,dby,B1,B2,B3,V3;
  Real B1ch, B2ch, B3ch;
#endif
#ifdef H_CORRECTION
  Real cfr,cfl,lambdar,lambdal;
#endif
#if (NSCALARS > 0)
  int n;
#endif
#ifdef SELF_GRAVITY
  Real gxl,gxr,gyl,gyr,flux_m1l,flux_m1r,flux_m2l,flux_m2r;
#endif
#ifdef FEEDBACK
  Real dt1 = 1.0/pG->dt;
#endif
#ifdef SHEARING_BOX
/* in XY shearing sheet 2=phi; in XZ shearing sheet 2=Z and 3=phi */
  Real Vphi, Mphi;
  Real M1n, dM2n, dM3n;   /* M1, dM2/3=(My+d*q*Omega_0*x) at time n */
  Real M1e, dM2e, dM3e;   /* M1, dM2/3 evolved by dt/2  */
  Real flx1_dM2, frx1_dM2, flx2_dM2, frx2_dM2;
  Real flx1_dM3, frx1_dM3, flx2_dM3, frx2_dM3;
  Real fact, qom, om_dt = Omega_0*pG->dt;
#endif /* SHEARING_BOX */
#ifdef STATIC_MESH_REFINEMENT
  int ncg,npg,dim;
  int ii,ics,ice,jj,jcs,jce,ips,ipe,jps,jpe;
#endif

  /* VARIABLES NEEDED FOR CYLINDRICAL COORDINATES */
#ifdef CYLINDRICAL
#ifndef ISOTHERMAL
  Real Pavgh;
#endif
  Real g,gl,gr,rinv;
  Real geom_src_d,geom_src_Vx,geom_src_Vy,geom_src_P,geom_src_By,geom_src_Bz;
  const Real *r=pG->r, *ri=pG->ri;
#endif /* CYLINDRICAL */
  Real lsf=1.0, rsf=1.0;

/* With particles, one more ghost cell must be updated in predict step */
#ifdef PARTICLES
  Real d1;
  il = is - 3;
  iu = ie + 3;
  jl = js - 3;
  ju = je + 3;
#else
  il = is - 2;
  iu = ie + 2;
  jl = js - 2;
  ju = je + 2;
#endif

/* Set etah=0 so first calls to flux functions do not use H-correction */
  etah = 0.0;

/* Compute predictor feedback from particle drag */
#ifdef FEEDBACK
  feedback_predictor(pG);
#endif

/*=== STEP 1: Compute L/R x1-interface states and 1D x1-Fluxes ===============*/

/*--- Step 1a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M1, M2, M3, E, B2c, B3c, s[n])
 */

  for (j=jl; j<=ju; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      U1d[i].d  = pG->U[ks][j][i].d;
      U1d[i].Mx = pG->U[ks][j][i].M1;
      U1d[i].My = pG->U[ks][j][i].M2;
      U1d[i].Mz = pG->U[ks][j][i].M3;
#ifndef BAROTROPIC
      U1d[i].E  = pG->U[ks][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
      U1d[i].By = pG->U[ks][j][i].B2c;
      U1d[i].Bz = pG->U[ks][j][i].B3c;
      Bxc[i] = pG->U[ks][j][i].B1c;
      Bxi[i] = pG->B1i[ks][j][i];
      B1_x1Face[j][i] = pG->B1i[ks][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) U1d[i].s[n] = pG->U[ks][j][i].s[n];
#endif
    }

/*--- Step 1b ------------------------------------------------------------------
 * Compute L and R states at X1-interfaces, add "MHD source terms" for 0.5*dt
 */

    for (i=is-nghost; i<=ie+nghost; i++) {
      W[i] = Cons1D_to_Prim1D(&U1d[i],&Bxc[i]);

      /* CALCULATE THE CELL-CENTERED GEOMETRIC SOURCE VECTOR NOW USING U^{n}
      * THIS WILL BE USED AT THE END OF THE INTEGRATION STEP AS A SOURCE TERM
      * FOR THE CELL-CENTERED CONSERVED VARIABLES (STEPS 6D,8B) */
#ifdef CYLINDRICAL 
      geom_src[j][i]  = W[i].d*SQR(W[i].Vy);
#ifdef MHD
      geom_src[j][i] += 0.5*(SQR(Bxc[i]) - SQR(W[i].By) + SQR(W[i].Bz));
#endif
#ifdef ISOTHERMAL
      geom_src[j][i] += Iso_csound2*W[i].d;
#else
      geom_src[j][i] += W[i].P;
#endif
      geom_src[j][i] /= x1vc(pG,i);
#endif /* CYLINDRICAL */
    }

    lr_states(pG,W,Bxc,pG->dt,pG->dx1,il+1,iu-1,Wl,Wr,1);

#ifdef MHD
    for (i=il+1; i<=iu; i++) {
#ifdef CYLINDRICAL
      rsf = ri[i]/r[i-1];  lsf = ri[i-1]/r[i-1];
#endif
      MHD_src = (pG->U[ks][j][i-1].M2/pG->U[ks][j][i-1].d)*
               (rsf*pG->B1i[ks][j][i] - lsf*pG->B1i[ks][j][i-1])/pG->dx1;
      Wl[i].By += hdt*MHD_src;

#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
      MHD_src = (pG->U[ks][j][i].M2/pG->U[ks][j][i].d)*
               (rsf*pG->B1i[ks][j][i+1] - lsf*pG->B1i[ks][j][i])/pG->dx1;
      Wr[i].By += hdt*MHD_src;
    }
#endif /* MHD */

/*--- Step 1c ------------------------------------------------------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */

    if (StaticGravPot != NULL){
      for (i=il+1; i<=iu; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
#ifdef CYLINDRICAL
        gl = (*x1GravAcc)(x1vc(pG,i-1),x2,x3);
        gr = (*x1GravAcc)(x1vc(pG,i),x2,x3);

        /* APPLY GRAV. SOURCE TERMS TO VELOCITY USING ACCELERATION FOR (dt/2) */
        Wl[i].Vx -= hdt*gl;
        Wr[i].Vx -= hdt*gr;
#else
        phicr = (*StaticGravPot)( x1             ,x2,x3);
        phicl = (*StaticGravPot)((x1-    pG->dx1),x2,x3);
        phifc = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

        Wl[i].Vx -= dtodx1*(phifc - phicl);
        Wr[i].Vx -= dtodx1*(phicr - phifc);
#endif
      }
    }

/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms for self-gravity for 0.5*dt to L/R states
 */

#ifdef SELF_GRAVITY
    for (i=il+1; i<=iu; i++) {
      Wl[i].Vx -= hdtodx1*(pG->Phi[ks][j][i] - pG->Phi[ks][j][i-1]);
      Wr[i].Vx -= hdtodx1*(pG->Phi[ks][j][i] - pG->Phi[ks][j][i-1]);
    }
#endif

/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms from optically-thin cooling for 0.5*dt to L/R states
 */

#ifndef BAROTROPIC
    if (CoolingFunc != NULL){
      for (i=il+1; i<=iu; i++) {
        coolfl = (*CoolingFunc)(Wl[i].d,Wl[i].P,(0.5*pG->dt));
        coolfr = (*CoolingFunc)(Wr[i].d,Wr[i].P,(0.5*pG->dt));

        Wl[i].P -= 0.5*pG->dt*Gamma_1*coolfl;
        Wr[i].P -= 0.5*pG->dt*Gamma_1*coolfr;
      }
    }
#endif /* BAROTROPIC */

/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms for shearing box (Coriolis forces) for 0.5*dt to L/R states
 * starting with tidal gravity terms added through the ShearingBoxPot
 *    Vx source term = (dt/2)*( 2 Omega_0 Vy)
 *    Vy source term = (dt/2)*(-2 Omega_0 Vx)
 *    Vy source term = (dt/2)*((q-2) Omega_0 Vx) (with FARGO)
 * (x1,x2,x3) in code = (X,Z,Y) in 2D shearing sheet
 */

#ifdef SHEARING_BOX
    if (ShearingBoxPot != NULL){
      for (i=il+1; i<=iu; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phicr = (*ShearingBoxPot)( x1             ,x2,x3);
        phicl = (*ShearingBoxPot)((x1-    pG->dx1),x2,x3);
        phifc = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);

        Wl[i].Vx -= dtodx1*(phifc - phicl);
        Wr[i].Vx -= dtodx1*(phicr - phifc);
      }
    }

    if (ShBoxCoord == xz){
      for (i=il+1; i<=iu; i++) {
        Wl[i].Vx += pG->dt*Omega_0*W[i-1].Vz;
        Wr[i].Vx += pG->dt*Omega_0*W[i].Vz;
#ifdef FARGO
        Wl[i].Vz += hdt*(qshear-2.)*Omega_0*W[i-1].Vx;
        Wr[i].Vz += hdt*(qshear-2.)*Omega_0*W[i].Vx;
#else
        Wl[i].Vz -= pG->dt*Omega_0*W[i-1].Vx;
        Wr[i].Vz -= pG->dt*Omega_0*W[i].Vx;
#endif
      }
    }

    if (ShBoxCoord == xy) {
      for (i=il+1; i<=iu; i++) {
        Wl[i].Vx += pG->dt*Omega_0*W[i-1].Vy;
        Wr[i].Vx += pG->dt*Omega_0*W[i].Vy;
#ifdef FARGO
        Wl[i].Vy += hdt*(qshear-2.)*Omega_0*W[i-1].Vx;
        Wr[i].Vy += hdt*(qshear-2.)*Omega_0*W[i].Vx;
#else
        Wl[i].Vy -= pG->dt*Omega_0*W[i-1].Vx;
        Wr[i].Vy -= pG->dt*Omega_0*W[i].Vx;
#endif
      }
    }
#endif /* SHEARING_BOX */

/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms for particle feedback for 0.5*dt to L/R states
 */

#ifdef FEEDBACK
    for (i=il+1; i<=iu; i++) {
      d1 = 1.0/W[i-1].d;
      Wl[i].Vx -= pG->Coup[ks][j][i-1].fb1*d1;
      Wl[i].Vy -= pG->Coup[ks][j][i-1].fb2*d1;
      Wl[i].Vz -= pG->Coup[ks][j][i-1].fb3*d1;

      d1 = 1.0/W[i].d;
      Wr[i].Vx -= pG->Coup[ks][j][i].fb1*d1;
      Wr[i].Vy -= pG->Coup[ks][j][i].fb2*d1;
      Wr[i].Vz -= pG->Coup[ks][j][i].fb3*d1;

#ifndef BAROTROPIC
      Wl[i].P += pG->Coup[ks][j][i-1].Eloss*Gamma_1;
      Wr[i].P += pG->Coup[ks][j][i].Eloss*Gamma_1;
#endif
    }
#endif /* FEEDBACK */

/*--- Step 1c (cont) -----------------------------------------------------------
 * ADD THE GEOMETRIC SOURCE-TERMS NOW USING CELL-CENTERED PRIMITIVE
 * VARIABLES AT TIME t^n
 */
#ifdef CYLINDRICAL
      for (i=il+1; i<=iu; i++) {

        // LEFT STATE GEOMETRIC SOURCE TERM (USES W[i-1])
        rinv = 1.0/x1vc(pG,i-1);
        geom_src_d  = -W[i-1].d*W[i-1].Vx*rinv;
        geom_src_Vx =  SQR(W[i-1].Vy);
        geom_src_Vy = -W[i-1].Vx*W[i-1].Vy;
#ifdef MHD
        geom_src_Vx -= SQR(W[i-1].By)/W[i-1].d;
        geom_src_Vy += Bxc[i-1]*W[i-1].By/W[i-1].d;
        geom_src_By =  -W[i-1].Vy*Bxc[i-1]*rinv;
        geom_src_Bz =  -W[i-1].Vx*W[i-1].Bz*rinv;
#endif /* MHD */
        geom_src_Vx *= rinv;
        geom_src_Vy *= rinv;
#ifndef ISOTHERMAL
        geom_src_P  = -Gamma*W[i-1].P*W[i-1].Vx*rinv;
#endif /* ISOTHERMAL */

        // ADD SOURCE TERM TO LEFT STATE
        Wl[i].d  += hdt*geom_src_d;
        Wl[i].Vx += hdt*geom_src_Vx;
        Wl[i].Vy += hdt*geom_src_Vy;
#ifdef MHD
        Wl[i].By += hdt*geom_src_By;
        Wl[i].Bz += hdt*geom_src_Bz;
#endif /* MHD */
#ifndef ISOTHERMAL
        Wl[i].P  += hdt*geom_src_P;
#endif /* ISOTHERMAL */

        // RIGHT STATE GEOMETRIC SOURCE TERM (USES W[i])
        rinv = 1.0/x1vc(pG,i);
        geom_src_d  = -W[i].d*W[i].Vx*rinv;
        geom_src_Vx =  SQR(W[i].Vy);
        geom_src_Vy = -W[i].Vx*W[i].Vy;
#ifdef MHD
        geom_src_Vx -= SQR(W[i].By)/W[i].d;
        geom_src_Vy += Bxc[i]*W[i].By/W[i].d;
        geom_src_By =  -W[i].Vy*Bxc[i]*rinv;
        geom_src_Bz =  -W[i].Vx*W[i].Bz*rinv;
#endif /* MHD */
        geom_src_Vx *= rinv;
        geom_src_Vy *= rinv;
#ifndef ISOTHERMAL
        geom_src_P  = -Gamma*W[i].P*W[i].Vx*rinv;
#endif /* ISOTHERMAL */

        // ADD SOURCE TERM TO RIGHT STATE
        Wr[i].d  += hdt*geom_src_d;
        Wr[i].Vx += hdt*geom_src_Vx;
        Wr[i].Vy += hdt*geom_src_Vy;
#ifdef MHD
        Wr[i].By += hdt*geom_src_By;
        Wr[i].Bz += hdt*geom_src_Bz;
#endif /* MHD */
#ifndef ISOTHERMAL
        Wr[i].P  += hdt*geom_src_P;
#endif /* ISOTHERMAL */
      }
#endif /* CYLINDRICAL */

/*--- Step 1d ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction, storing into 2D array
 */

    for (i=il+1; i<=iu; i++) {
      Ul_x1Face[j][i] = Prim1D_to_Cons1D(&Wl[i],&Bxi[i]);
      Ur_x1Face[j][i] = Prim1D_to_Cons1D(&Wr[i],&Bxi[i]);

#ifdef MHD
      Bx = B1_x1Face[j][i];
#endif
      fluxes(Ul_x1Face[j][i],Ur_x1Face[j][i],Wl[i],Wr[i],Bx,&x1Flux[j][i]);
    }
  }

/*=== STEP 2: Compute L/R x2-interface states and 1D x2-Fluxes ===============*/

/*--- Step 2a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M2, M3, M1, E, B3c, B1c, s[n])
 */

  for (i=il; i<=iu; i++) {
#ifdef CYLINDRICAL
    dx2 = r[i]*pG->dx2;
    dtodx2 = pG->dt/dx2;
    hdtodx2 = 0.5*dtodx2;
#endif
    for (j=js-nghost; j<=je+nghost; j++) {
      U1d[j].d  = pG->U[ks][j][i].d;
      U1d[j].Mx = pG->U[ks][j][i].M2;
      U1d[j].My = pG->U[ks][j][i].M3;
      U1d[j].Mz = pG->U[ks][j][i].M1;
#ifndef BAROTROPIC
      U1d[j].E  = pG->U[ks][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
      U1d[j].By = pG->U[ks][j][i].B3c;
      U1d[j].Bz = pG->U[ks][j][i].B1c;
      Bxc[j] = pG->U[ks][j][i].B2c;
      Bxi[j] = pG->B2i[ks][j][i];
      B2_x2Face[j][i] = pG->B2i[ks][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) U1d[j].s[n] = pG->U[ks][j][i].s[n];
#endif
    }

/*--- Step 2b ------------------------------------------------------------------
 * Compute L and R states at X2-interfaces, add "MHD source terms" for 0.5*dt
 */

    for (j=js-nghost; j<=je+nghost; j++) {
      W[j] = Cons1D_to_Prim1D(&U1d[j],&Bxc[j]);
    }

    lr_states(pG,W,Bxc,pG->dt,dx2,jl+1,ju-1,Wl,Wr,2);

#ifdef MHD
    for (j=jl+1; j<=ju; j++) {
      MHD_src = (pG->U[ks][j-1][i].M1/pG->U[ks][j-1][i].d)*
        (pG->B2i[ks][j][i] - pG->B2i[ks][j-1][i])/dx2;
      Wl[j].Bz += hdt*MHD_src;

      MHD_src = (pG->U[ks][j][i].M1/pG->U[ks][j][i].d)*
        (pG->B2i[ks][j+1][i] - pG->B2i[ks][j][i])/dx2;
      Wr[j].Bz += hdt*MHD_src;
    }
#endif

/*--- Step 2c ------------------------------------------------------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */
  
    if (StaticGravPot != NULL){
      for (j=jl+1; j<=ju; j++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phicr = (*StaticGravPot)(x1, x2             ,x3);
        phicl = (*StaticGravPot)(x1,(x2-    pG->dx2),x3);
        phifc = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        Wl[j].Vx -= dtodx2*(phifc - phicl);
        Wr[j].Vx -= dtodx2*(phicr - phifc);
      }
    }

/*--- Step 2c (cont) -----------------------------------------------------------
 * Add source terms for self-gravity for 0.5*dt to L/R states
 */

#ifdef SELF_GRAVITY
    for (j=jl+1; j<=ju; j++) {
      Wl[j].Vx -= hdtodx2*(pG->Phi[ks][j][i] - pG->Phi[ks][j-1][i]);
      Wr[j].Vx -= hdtodx2*(pG->Phi[ks][j][i] - pG->Phi[ks][j-1][i]);
    }
#endif

/*--- Step 2c (cont) -----------------------------------------------------------
 * Add source terms from optically-thin cooling for 0.5*dt to L/R states
 */

#ifndef BAROTROPIC
    if (CoolingFunc != NULL){
      for (j=jl+1; j<=ju; j++) {
        coolfl = (*CoolingFunc)(Wl[j].d,Wl[j].P,(0.5*pG->dt));
        coolfr = (*CoolingFunc)(Wr[j].d,Wr[j].P,(0.5*pG->dt));

        Wl[j].P -= 0.5*pG->dt*Gamma_1*coolfl;
        Wr[j].P -= 0.5*pG->dt*Gamma_1*coolfr;
      }
    }
#endif /* BAROTROPIC */

/*--- Step 2c (cont) -----------------------------------------------------------
 * Add source terms for particle feedback for 0.5*dt to L/R states
 */

#ifdef FEEDBACK
   for (j=jl+1; j<=ju; j++) {
      d1 = 1.0/W[j-1].d;
      Wl[j].Vx -= pG->Coup[ks][j-1][i].fb2*d1;
      Wl[j].Vy -= pG->Coup[ks][j-1][i].fb3*d1;
      Wl[j].Vz -= pG->Coup[ks][j-1][i].fb1*d1;

      d1 = 1.0/W[j].d;
      Wr[j].Vx -= pG->Coup[ks][j][i].fb2*d1;
      Wr[j].Vy -= pG->Coup[ks][j][i].fb3*d1;
      Wr[j].Vz -= pG->Coup[ks][j][i].fb1*d1;

#ifndef BAROTROPIC
      Wl[i].P += pG->Coup[ks][j-1][i].Eloss*Gamma_1;
      Wr[i].P += pG->Coup[ks][j][i].Eloss*Gamma_1;
#endif
    }
#endif /* FEEDBACK */

/*--- Step 2d ------------------------------------------------------------------
 * Compute 1D fluxes in x2-direction, storing into 2D array
 */

    for (j=jl+1; j<=ju; j++) {
      Ul_x2Face[j][i] = Prim1D_to_Cons1D(&Wl[j],&Bxi[j]);
      Ur_x2Face[j][i] = Prim1D_to_Cons1D(&Wr[j],&Bxi[j]);
#ifdef MHD
      Bx = B2_x2Face[j][i];
#endif
      fluxes(Ul_x2Face[j][i],Ur_x2Face[j][i],Wl[j],Wr[j],Bx,&x2Flux[j][i]);
    }
  }

/*=== STEP 3: Not needed in 2D ===*/

/*=== STEP 4:  Update face-centered B for 0.5*dt =============================*/

/*--- Step 4a ------------------------------------------------------------------
 * Calculate the cell centered value of emf_3 at t^{n} and integrate
 * to corner.
 */

#ifdef MHD
  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      emf3_cc[j][i] =
	(pG->U[ks][j][i].B1c*pG->U[ks][j][i].M2 -
	 pG->U[ks][j][i].B2c*pG->U[ks][j][i].M1 )/pG->U[ks][j][i].d;
    }
  }
  integrate_emf3_corner(pG);

/*--- Step 4b ------------------------------------------------------------------
 * Update the interface magnetic fields using CT for a half time step.
 */

  for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu-1; i++) {
#ifdef CYLINDRICAL
      hdtodx2 = hdt/(ri[i]*pG->dx2);
#endif
      B1_x1Face[j][i] -= hdtodx2*(emf3[j+1][i  ] - emf3[j][i]);
      B2_x2Face[j][i] += hdtodx1*(emf3[j  ][i+1] - emf3[j][i]);
    }
#ifdef CYLINDRICAL
    hdtodx2 = hdt/(ri[iu]*pG->dx2);
#endif
    B1_x1Face[j][iu] -= hdtodx2*(emf3[j+1][iu] - emf3[j][iu]);
  }
  for (i=il+1; i<=iu-1; i++) {
    B2_x2Face[ju][i] += hdtodx1*(emf3[ju][i+1] - emf3[ju][i]);
  }
#endif /* MHD */

/*=== STEP 5: Correct x1-interface states with transverse flux gradients =====*/

/*--- Step 5a ------------------------------------------------------------------
 * Correct x1-interface states using x2-fluxes computed in Step 2d.
 * Since the fluxes come from an x2-sweep, (x,y,z) on RHS -> (z,x,y) on LHS
 */

  for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu; i++) {
#ifdef CYLINDRICAL
      hdtodx2 = hdt/(r[i-1]*pG->dx2);
#endif
      Ul_x1Face[j][i].d  -= hdtodx2*(x2Flux[j+1][i-1].d  - x2Flux[j][i-1].d );
      Ul_x1Face[j][i].Mx -= hdtodx2*(x2Flux[j+1][i-1].Mz - x2Flux[j][i-1].Mz);
      Ul_x1Face[j][i].My -= hdtodx2*(x2Flux[j+1][i-1].Mx - x2Flux[j][i-1].Mx);
      Ul_x1Face[j][i].Mz -= hdtodx2*(x2Flux[j+1][i-1].My - x2Flux[j][i-1].My);
#ifndef BAROTROPIC
      Ul_x1Face[j][i].E  -= hdtodx2*(x2Flux[j+1][i-1].E  - x2Flux[j][i-1].E );
#endif /* BAROTROPIC */
#ifdef MHD
      Ul_x1Face[j][i].Bz -= hdtodx2*(x2Flux[j+1][i-1].By - x2Flux[j][i-1].By);
#endif
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) {
      Ul_x1Face[j][i].s[n]-=hdtodx2*(x2Flux[j+1][i-1].s[n]-x2Flux[j][i-1].s[n]);
      }
#endif

#ifdef CYLINDRICAL
      hdtodx2 = hdt/(r[i]*pG->dx2);
#endif
      Ur_x1Face[j][i].d  -= hdtodx2*(x2Flux[j+1][i  ].d  - x2Flux[j][i  ].d );
      Ur_x1Face[j][i].Mx -= hdtodx2*(x2Flux[j+1][i  ].Mz - x2Flux[j][i  ].Mz);
      Ur_x1Face[j][i].My -= hdtodx2*(x2Flux[j+1][i  ].Mx - x2Flux[j][i  ].Mx);
      Ur_x1Face[j][i].Mz -= hdtodx2*(x2Flux[j+1][i  ].My - x2Flux[j][i  ].My);
#ifndef BAROTROPIC
      Ur_x1Face[j][i].E  -= hdtodx2*(x2Flux[j+1][i  ].E  - x2Flux[j][i  ].E );
#endif /* BAROTROPIC */
#ifdef MHD
      Ur_x1Face[j][i].Bz -= hdtodx2*(x2Flux[j+1][i  ].By - x2Flux[j][i  ].By);
#endif
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) {
      Ur_x1Face[j][i].s[n]-=hdtodx2*(x2Flux[j+1][i  ].s[n]-x2Flux[j][i  ].s[n]);
      }
#endif
    }
  }

/*--- Step 5b: Not needed in 2D ---*/
/*--- Step 5c ------------------------------------------------------------------
 * Add the "MHD source terms" from the x2-flux gradient to the conservative
 * variables on the x1Face..
 */

#ifdef MHD
  for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu; i++) {
#ifdef CYLINDRICAL
      dbx = (ri[i]*pG->B1i[ks][j][i] - ri[i-1]*pG->B1i[ks][j][i-1])/r[i-1];
#else
      dbx = pG->B1i[ks][j][i] - pG->B1i[ks][j][i-1];
#endif
      B1 = pG->U[ks][j][i-1].B1c;
      B2 = pG->U[ks][j][i-1].B2c;
      B3 = pG->U[ks][j][i-1].B3c;
      V3 = pG->U[ks][j][i-1].M3/pG->U[ks][j][i-1].d;

      Ul_x1Face[j][i].Mx += hdtodx1*B1*dbx;
      Ul_x1Face[j][i].My += hdtodx1*B2*dbx;
      Ul_x1Face[j][i].Mz += hdtodx1*B3*dbx;
      Ul_x1Face[j][i].Bz += hdtodx1*V3*dbx;
#ifndef BAROTROPIC
      Ul_x1Face[j][i].E  += hdtodx1*B3*V3*dbx;
#endif /* BAROTROPIC */

#ifdef CYLINDRICAL
      dbx = (ri[i+1]*pG->B1i[ks][j][i+1] - ri[i]*pG->B1i[ks][j][i])/r[i];
#else
      dbx = pG->B1i[ks][j][i+1] - pG->B1i[ks][j][i];
#endif
      B1 = pG->U[ks][j][i].B1c;
      B2 = pG->U[ks][j][i].B2c;
      B3 = pG->U[ks][j][i].B3c;
      V3 = pG->U[ks][j][i].M3/pG->U[ks][j][i].d;

      Ur_x1Face[j][i].Mx += hdtodx1*B1*dbx;
      Ur_x1Face[j][i].My += hdtodx1*B2*dbx;
      Ur_x1Face[j][i].Mz += hdtodx1*B3*dbx;
      Ur_x1Face[j][i].Bz += hdtodx1*V3*dbx;
#ifndef BAROTROPIC
      Ur_x1Face[j][i].E  += hdtodx1*B3*V3*dbx;
#endif /* BAROTROPIC */
    }
  }
#endif /* MHD */

/*--- Step 5d ------------------------------------------------------------------
 * Add source terms for a static gravitational potential arising from x2-Flux
 * gradients.  To improve conservation of total energy, average
 * the energy source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (j=jl+1; j<=ju-1; j++) {
      for (i=il+1; i<=iu; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phic = (*StaticGravPot)(x1, x2             ,x3);
        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

#ifdef CYLINDRICAL
        hdtodx2 = hdt/(r[i]*pG->dx2);
#endif
        Ur_x1Face[j][i].My -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;
#ifndef BAROTROPIC
        Ur_x1Face[j][i].E -= hdtodx2*(x2Flux[j  ][i  ].d*(phic - phil) +
                                      x2Flux[j+1][i  ].d*(phir - phic));
#endif

        phic = (*StaticGravPot)((x1-pG->dx1), x2             ,x3);
        phir = (*StaticGravPot)((x1-pG->dx1),(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)((x1-pG->dx1),(x2-0.5*pG->dx2),x3);

#ifdef CYLINDRICAL
        hdtodx2 = hdt/(r[i-1]*pG->dx2);
#endif
        Ul_x1Face[j][i].My -= hdtodx2*(phir-phil)*pG->U[ks][j][i-1].d;
#ifndef BAROTROPIC
        Ul_x1Face[j][i].E -= hdtodx2*(x2Flux[j  ][i-1].d*(phic - phil) +
                                      x2Flux[j+1][i-1].d*(phir - phic));
#endif
      }
    }
  }

/*--- Step 5d (cont) -----------------------------------------------------------
 * Add source terms for self gravity arising from x2-Flux gradient
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

#ifdef SELF_GRAVITY
  for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu; i++) {
      phic = pG->Phi[ks][j][i];
      phir = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j+1][i]);
      phil = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j-1][i]);

      Ur_x1Face[j][i].My -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;
#ifndef BAROTROPIC
      Ur_x1Face[j][i].E -= hdtodx2*(x2Flux[j  ][i  ].d*(phic - phil) +
                                    x2Flux[j+1][i  ].d*(phir - phic));
#endif

      phic = pG->Phi[ks][j][i-1];
      phir = 0.5*(pG->Phi[ks][j][i-1] + pG->Phi[ks][j+1][i-1]);
      phil = 0.5*(pG->Phi[ks][j][i-1] + pG->Phi[ks][j-1][i-1]);

      Ul_x1Face[j][i].My -= hdtodx2*(phir-phil)*pG->U[ks][j][i-1].d;
#ifndef BAROTROPIC
      Ul_x1Face[j][i].E -= hdtodx2*(x2Flux[j  ][i-1].d*(phic - phil) +
                                    x2Flux[j+1][i-1].d*(phir - phic));
#endif
    }
  }
#endif /* SELF_GRAVITY */

/*=== STEP 6: Correct x2-interface states with transverse flux gradients =====*/

/*--- Step 6a ------------------------------------------------------------------
 * Correct x2-interface states using x1-fluxes computed in Step 1d.
 * Since the fluxes come from an x1-sweep, (x,y,z) on RHS -> (y,z,x) on LHS
 */

  for (j=jl+1; j<=ju; j++) {
    for (i=il+1; i<=iu-1; i++) {
#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
      Ul_x2Face[j][i].d  -= hdtodx1*(rsf*x1Flux[j-1][i+1].d  - lsf*x1Flux[j-1][i].d );
      Ul_x2Face[j][i].Mx -= hdtodx1*(SQR(rsf)*x1Flux[j-1][i+1].My - SQR(lsf)*x1Flux[j-1][i].My);
      Ul_x2Face[j][i].My -= hdtodx1*(rsf*x1Flux[j-1][i+1].Mz - lsf*x1Flux[j-1][i].Mz);
      Ul_x2Face[j][i].Mz -= hdtodx1*(rsf*x1Flux[j-1][i+1].Mx - lsf*x1Flux[j-1][i].Mx);
#ifndef BAROTROPIC
      Ul_x2Face[j][i].E  -= hdtodx1*(rsf*x1Flux[j-1][i+1].E  - lsf*x1Flux[j-1][i].E );
#endif /* BAROTROPIC */
#ifdef MHD
      Ul_x2Face[j][i].By -= hdtodx1*(rsf*x1Flux[j-1][i+1].Bz - lsf*x1Flux[j-1][i].Bz);
#endif

      Ur_x2Face[j][i].d  -= hdtodx1*(rsf*x1Flux[j  ][i+1].d  - lsf*x1Flux[j  ][i].d );
      Ur_x2Face[j][i].Mx -= hdtodx1*(SQR(rsf)*x1Flux[j  ][i+1].My - SQR(lsf)*x1Flux[j  ][i].My);
      Ur_x2Face[j][i].My -= hdtodx1*(rsf*x1Flux[j  ][i+1].Mz - lsf*x1Flux[j  ][i].Mz);
      Ur_x2Face[j][i].Mz -= hdtodx1*(rsf*x1Flux[j  ][i+1].Mx - lsf*x1Flux[j  ][i].Mx);
#ifndef BAROTROPIC
      Ur_x2Face[j][i].E  -= hdtodx1*(rsf*x1Flux[j  ][i+1].E  - lsf*x1Flux[j  ][i].E );
#endif /* BAROTROPIC */
#ifdef MHD
      Ur_x2Face[j][i].By -= hdtodx1*(rsf*x1Flux[j  ][i+1].Bz - lsf*x1Flux[j  ][i].Bz);
#endif
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) {
      Ul_x2Face[j][i].s[n]-=hdtodx1*(rsf*x1Flux[j-1][i+1].s[n]-lsf*x1Flux[j-1][i].s[n]);
      Ur_x2Face[j][i].s[n]-=hdtodx1*(rsf*x1Flux[j  ][i+1].s[n]-lsf*x1Flux[j  ][i].s[n]);
      }
#endif
    }
  }

/*--- Step 6b: Not needed in 2D ---*/
/*--- Step 6c ------------------------------------------------------------------
 * Add the "MHD source terms" from the x1-flux-gradients to the
 * conservative variables on the x2Face.
 */

#ifdef MHD
  for (j=jl+1; j<=ju; j++) {
    for (i=il+1; i<=iu-1; i++) {
#ifdef CYLINDRICAL
      hdtodx2 = hdt/(r[i]*pG->dx2);
#endif
      dby = pG->B2i[ks][j][i] - pG->B2i[ks][j-1][i];
      B1 = pG->U[ks][j-1][i].B1c;
      B2 = pG->U[ks][j-1][i].B2c;
      B3 = pG->U[ks][j-1][i].B3c;
      V3 = pG->U[ks][j-1][i].M3/pG->U[ks][j-1][i].d;

      Ul_x2Face[j][i].Mz += hdtodx2*B1*dby;
      Ul_x2Face[j][i].Mx += hdtodx2*B2*dby;
      Ul_x2Face[j][i].My += hdtodx2*B3*dby;
      Ul_x2Face[j][i].By += hdtodx2*V3*dby;
#ifndef BAROTROPIC
      Ul_x2Face[j][i].E  += hdtodx2*B3*V3*dby;
#endif /* BAROTROPIC */

      dby = pG->B2i[ks][j+1][i] - pG->B2i[ks][j][i];
      B1 = pG->U[ks][j][i].B1c;
      B2 = pG->U[ks][j][i].B2c;
      B3 = pG->U[ks][j][i].B3c;
      V3 = pG->U[ks][j][i].M3/pG->U[ks][j][i].d;

      Ur_x2Face[j][i].Mz += hdtodx2*B1*dby;
      Ur_x2Face[j][i].Mx += hdtodx2*B2*dby;
      Ur_x2Face[j][i].My += hdtodx2*B3*dby;
      Ur_x2Face[j][i].By += hdtodx2*V3*dby;
#ifndef BAROTROPIC
      Ur_x2Face[j][i].E  += hdtodx2*B3*V3*dby;
#endif /* BAROTROPIC */
    }
  }
#endif /* MHD */

/*--- Step 6d ------------------------------------------------------------------
 * Add source terms for a static gravitational potential arising from x1-Flux
 * gradients. To improve conservation of total energy, average the energy
 * source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (j=jl+1; j<=ju; j++) {
      for (i=il+1; i<=iu-1; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phic = (*StaticGravPot)((x1            ),x2,x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

#ifdef CYLINDRICAL
        g = (*x1GravAcc)(x1vc(pG,i),x2,x3);
        rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
        Ur_x2Face[j][i].Mz -= hdt*pG->U[ks][j][i].d*g;
#else
        Ur_x2Face[j][i].Mz -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;
#endif
#ifndef BAROTROPIC
        Ur_x2Face[j][i].E -= hdtodx1*(lsf*x1Flux[j  ][i  ].d*(phic - phil) +
                                      rsf*x1Flux[j  ][i+1].d*(phir - phic));
#endif

        phic = (*StaticGravPot)((x1            ),(x2-pG->dx2),x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),(x2-pG->dx2),x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),(x2-pG->dx2),x3);

#ifdef CYLINDRICAL
        g = (*x1GravAcc)(x1vc(pG,i),(x2-pG->dx2),x3);
        Ul_x2Face[j][i].Mz -= hdt*pG->U[ks][j-1][i].d*g;
#else
        Ul_x2Face[j][i].Mz -= hdtodx1*(phir-phil)*pG->U[ks][j-1][i].d;
#endif
#ifndef BAROTROPIC
        Ul_x2Face[j][i].E -= hdtodx1*(lsf*x1Flux[j-1][i  ].d*(phic - phil) +
                                      rsf*x1Flux[j-1][i+1].d*(phir - phic));
#endif
      }
    }
  }

/*--- Step 6d (cont) -----------------------------------------------------------
 * Add source terms for self gravity arising from x1-Flux gradients
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

#ifdef SELF_GRAVITY
  for (j=jl+1; j<=ju; j++) {
    for (i=il+1; i<=iu-1; i++) {
      phic = pG->Phi[ks][j][i];
      phir = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j][i+1]);
      phil = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j][i-1]);

      Ur_x2Face[j][i].Mz -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;
#ifndef BAROTROPIC
      Ur_x2Face[j][i].E -= hdtodx1*(x1Flux[j  ][i  ].d*(phic - phil) +
                                    x1Flux[j  ][i+1].d*(phir - phic));
#endif

      phic = pG->Phi[ks][j-1][i];
      phir = 0.5*(pG->Phi[ks][j-1][i] + pG->Phi[ks][j-1][i+1]);
      phil = 0.5*(pG->Phi[ks][j-1][i] + pG->Phi[ks][j-1][i-1]);

      Ul_x2Face[j][i].Mz -= hdtodx1*(phir-phil)*pG->U[ks][j-1][i].d;
#ifndef BAROTROPIC
      Ul_x2Face[j][i].E -= hdtodx1*(x1Flux[j-1][i  ].d*(phic - phil) +
                                    x1Flux[j-1][i+1].d*(phir - phic));
#endif
    }
  }

#endif /* SELF_GRAVITY */

/*--- Step 6d (cont) -----------------------------------------------------------
 * Add source terms for shearing box (Coriolis forces) for 0.5*dt arising from
 * x1-Flux gradient.  The tidal gravity terms are added via ShearingBoxPot
 *    Vx source term is (dt/2)( 2 Omega_0 Vy)
 *    Vy source term is (dt/2)(-2 Omega_0 Vx)
 *    Vy source term is (dt/2)((q-2) Omega_0 Vx) (with FARGO)
 * (x1,x2,x3) in code = (X,Z,Y) in shearing sheet
 */

#ifdef SHEARING_BOX
  if (ShearingBoxPot != NULL){
    for (j=jl+1; j<=ju; j++) {
      for (i=il+1; i<=iu-1; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phic = (*ShearingBoxPot)((x1            ),x2,x3);
        phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);

        Ur_x2Face[j][i].Mz -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;
#ifndef BAROTROPIC
        Ur_x2Face[j][i].E -= hdtodx1*(x1Flux[j  ][i  ].d*(phic - phil) +
                                      x1Flux[j  ][i+1].d*(phir - phic));
#endif

        phic = (*ShearingBoxPot)((x1            ),(x2-pG->dx2),x3);
        phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),(x2-pG->dx2),x3);
        phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),(x2-pG->dx2),x3);

        Ul_x2Face[j][i].Mz -= hdtodx1*(phir-phil)*pG->U[ks][j-1][i].d;
#ifndef BAROTROPIC
        Ul_x2Face[j][i].E -= hdtodx1*(x1Flux[j-1][i  ].d*(phic - phil) +
                                      x1Flux[j-1][i+1].d*(phir - phic));
#endif
      }
    }
  }

  if (ShBoxCoord == xz){
    for (j=jl+1; j<=ju; j++) {
      for (i=il+1; i<=iu-1; i++) {
        Ur_x2Face[j][i].Mz += pG->dt*Omega_0*pG->U[ks][j][i].M3;
        Ul_x2Face[j][i].Mz += pG->dt*Omega_0*pG->U[ks][j-1][i].M3;
#ifdef FARGO
        Ur_x2Face[j][i].My += hdt*(qshear-2.)*Omega_0*pG->U[ks][j][i].M1;
        Ul_x2Face[j][i].My += hdt*(qshear-2.)*Omega_0*pG->U[ks][j-1][i].M1;
#else
        Ur_x2Face[j][i].My -= pG->dt*Omega_0*pG->U[ks][j][i].M1;
        Ul_x2Face[j][i].My -= pG->dt*Omega_0*pG->U[ks][j-1][i].M1;
#endif
      }
    }
  }

  if (ShBoxCoord == xy){
    for (j=jl+1; j<=ju; j++) {
      for (i=il+1; i<=iu-1; i++) {
        Ur_x2Face[j][i].Mz += pG->dt*Omega_0*pG->U[ks][j][i].M2;
        Ul_x2Face[j][i].Mz += pG->dt*Omega_0*pG->U[ks][j-1][i].M2;
#ifdef FARGO
        Ur_x2Face[j][i].Mx += hdt*(qshear-2.)*Omega_0*pG->U[ks][j][i].M1;
        Ul_x2Face[j][i].Mx += hdt*(qshear-2.)*Omega_0*pG->U[ks][j-1][i].M1;
#else
        Ur_x2Face[j][i].Mx -= pG->dt*Omega_0*pG->U[ks][j][i].M1;
        Ul_x2Face[j][i].Mx -= pG->dt*Omega_0*pG->U[ks][j-1][i].M1;
#endif
      }
    }
  }
#endif /* SHEARING_BOX */

/*--- Step 6d (cont) -----------------------------------------------------------
 * ADD THE GEOMETRIC SOURCE TERMS IN THE x1-DIRECTION FOR dt/2
 */
#ifdef CYLINDRICAL
  for (j=js-1; j<=ju; j++) {
    for (i=is-1; i<=ie+1; i++) {
      Ur_x2Face[j][i].Mz += hdt*geom_src[j  ][i];
      Ul_x2Face[j][i].Mz += hdt*geom_src[j-1][i];
    }
  }
#endif

/*=== STEP 7: Not needed in 2D ===*/

/*=== STEP 8: Compute cell-centered values at n+1/2 ==========================*/

/*--- Step 8a ------------------------------------------------------------------
 * Calculate d^{n+1/2} (needed with static potential, cooling, or MHD)
 */

#ifndef MHD
#ifndef PARTICLES
  if ((StaticGravPot != NULL) || (CoolingFunc != NULL)) 
#endif
#endif
  {
    for (j=jl+1; j<=ju-1; j++) {
      for (i=il+1; i<=iu-1; i++) {
#ifdef CYLINDRICAL
        rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
        hdtodx2 = hdt/(r[i]*pG->dx2);
#endif
        dhalf[j][i] = pG->U[ks][j][i].d
          - hdtodx1*(rsf*x1Flux[j  ][i+1].d - lsf*x1Flux[j][i].d)
          - hdtodx2*(    x2Flux[j+1][i  ].d -     x2Flux[j][i].d);
#ifdef PARTICLES
        pG->Coup[ks][j][i].grid_d = dhalf[j][i];
#endif
      }
    }
  }

/*--- Step 8b ------------------------------------------------------------------
 * Calculate P^{n+1/2} (needed with cooling), and cell centered emf3_cc^{n+1/2}
 */

#ifndef MHD
#ifndef PARTICLES
  if (CoolingFunc != NULL) 
#endif /* PARTICLES */
#endif /* MHD */
  {
  for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu-1; i++) {
#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
      hdtodx2 = hdt/(r[i]*pG->dx2);
#endif
      M1h = pG->U[ks][j][i].M1
        - hdtodx1*(rsf*x1Flux[j][i+1].Mx - lsf*x1Flux[j][i].Mx)
        - hdtodx2*(    x2Flux[j+1][i].Mz -     x2Flux[j][i].Mz);

      M2h = pG->U[ks][j][i].M2
        - hdtodx1*(SQR(rsf)*x1Flux[j][i+1].My - SQR(lsf)*x1Flux[j][i].My)
        - hdtodx2*(         x2Flux[j+1][i].Mx -          x2Flux[j][i].Mx);

      M3h = pG->U[ks][j][i].M3
        - hdtodx1*(rsf*x1Flux[j][i+1].Mz - lsf*x1Flux[j][i].Mz)
        - hdtodx2*(    x2Flux[j+1][i].My -     x2Flux[j][i].My);

#ifndef BAROTROPIC
      Eh = pG->U[ks][j][i].E
        - hdtodx1*(rsf*x1Flux[j][i+1].E - lsf*x1Flux[j][i].E)
        - hdtodx2*(    x2Flux[j+1][i].E -     x2Flux[j][i].E);
#endif

/* Add source terms for fixed gravitational potential */
      if (StaticGravPot != NULL){
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
#ifdef CYLINDRICAL
        g = (*x1GravAcc)(x1vc(pG,i),x2,x3);
        M1h -= hdt*pG->U[ks][j][i].d*g;
#else
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);
        M1h -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;
#endif
        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);
        M2h -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;
      }

/* Add source terms due to self-gravity  */
#ifdef SELF_GRAVITY
      phir = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j][i+1]);
      phil = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j][i-1]);
      M1h -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;

      phir = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j+1][i]);
      phil = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j-1][i]);
      M2h -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;
#endif /* SELF_GRAVITY */

/* Add the tidal gravity and Coriolis terms for shearing box. */
#ifdef SHEARING_BOX
      if (ShearingBoxPot != NULL){
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);
        M1h -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;

        phir = (*ShearingBoxPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*ShearingBoxPot)(x1,(x2-0.5*pG->dx2),x3);
        M2h -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;
      }

      if (ShBoxCoord == xy) M1h += pG->dt*Omega_0*pG->U[ks][j][i].M2;
      if (ShBoxCoord == xz) M1h += pG->dt*Omega_0*pG->U[ks][j][i].M3;
#ifdef FARGO
      if (ShBoxCoord == xy) M2h += hdt*(qshear-2.)*Omega_0*pG->U[ks][j][i].M1;
      if (ShBoxCoord == xz) M3h += hdt*(qshear-2.)*Omega_0*pG->U[ks][j][i].M1;
#else
      if (ShBoxCoord == xy) M2h -= pG->dt*Omega_0*pG->U[ks][j][i].M1;
      if (ShBoxCoord == xz) M3h -= pG->dt*Omega_0*pG->U[ks][j][i].M1;
#endif
#endif /* SHEARING_BOX */

/* Add the particle feedback terms */
#ifdef FEEDBACK
      M1h -= pG->Coup[ks][j][i].fb1;
      M2h -= pG->Coup[ks][j][i].fb2;
      M3h -= pG->Coup[ks][j][i].fb3;
#endif /* FEEDBACK */


/* Add the geometric source term */
#ifdef CYLINDRICAL
      M1h += hdt*geom_src[j][i];
#endif

#ifndef BAROTROPIC
      phalf[j][i] = Eh - 0.5*(M1h*M1h + M2h*M2h + M3h*M3h)/dhalf[j][i];
#endif

#ifdef MHD
      B1ch = 0.5*(lsf*B1_x1Face[j][i] + rsf*B1_x1Face[j][i+1]);
      B2ch = 0.5*(    B2_x2Face[j][i] +     B2_x2Face[j+1][i]);
      B3ch = pG->U[ks][j][i].B3c 
        - hdtodx1*(rsf*x1Flux[j][i+1].Bz - lsf*x1Flux[j][i].Bz)
        - hdtodx2*(    x2Flux[j+1][i].By -     x2Flux[j][i].By);
      emf3_cc[j][i] = (B1ch*M2h - B2ch*M1h)/dhalf[j][i];

#ifndef BAROTROPIC
      phalf[j][i] -= 0.5*(B1ch*B1ch + B2ch*B2ch + B3ch*B3ch);
#endif
#endif /* MHD */

#ifndef BAROTROPIC
      phalf[j][i] *= Gamma_1;
#endif

#ifdef PARTICLES
      d1 = 1.0/dhalf[j][i];
      pG->Coup[ks][j][i].grid_v1 = M1h*d1;
      pG->Coup[ks][j][i].grid_v2 = M2h*d1;
      pG->Coup[ks][j][i].grid_v3 = M3h*d1;
#ifndef BAROTROPIC
      pG->Coup[ks][j][i].grid_cs = sqrt(Gamma*phalf[j][i]*d1);
#endif  /* BAROTROPIC */
#endif /* PARTICLES */
    }
  }
  }

/*=== STEP 8.5: Integrate the particles, compute the feedback ================*/

#ifdef PARTICLES
  Integrate_Particles(pG,pD);
#ifdef FEEDBACK
  exchange_feedback(pG, pD);
#endif
#endif

/*=== STEP 9: Compute 2D x1-Flux, x2-Flux ====================================*/

/*--- Step 9a ------------------------------------------------------------------
 * Compute maximum wavespeeds in multidimensions (eta in eq. 10 from Sanders et
 *  al. (1998)) for H-correction
 */

#ifdef H_CORRECTION
  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=ie+2; i++) {
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

  for (j=js-1; j<=je+2; j++) {
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

/*--- Step 9b ------------------------------------------------------------------
 * Compute 2D x1-fluxes from corrected L/R states.
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
      Wl[i] = Cons1D_to_Prim1D(&Ul_x1Face[j][i],&Bx);
      Wr[i] = Cons1D_to_Prim1D(&Ur_x1Face[j][i],&Bx);

      fluxes(Ul_x1Face[j][i],Ur_x1Face[j][i],Wl[i],Wr[i],Bx,&x1Flux[j][i]);
    }
  }

/*--- Step 9c ------------------------------------------------------------------
 * Compute 2D x2-fluxes from corrected L/R states.
 */

  for (j=js; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
#ifdef H_CORRECTION
      etah = MAX(eta1[j-1][i],eta1[j][i]);
      etah = MAX(etah,eta1[j-1][i+1]);
      etah = MAX(etah,eta1[j  ][i+1]);
      etah = MAX(etah,eta2[j  ][i  ]);
#endif /* H_CORRECTION */
#ifdef MHD
      Bx = B2_x2Face[j][i];
#endif
      Wl[i] = Cons1D_to_Prim1D(&Ul_x2Face[j][i],&Bx);
      Wr[i] = Cons1D_to_Prim1D(&Ur_x2Face[j][i],&Bx);

      fluxes(Ul_x2Face[j][i],Ur_x2Face[j][i],Wl[i],Wr[i],Bx,&x2Flux[j][i]);
    }
  }

/*=== STEP 10: Update face-centered B for a full timestep ====================*/

/*--- Step 10a -----------------------------------------------------------------
 * Integrate emf*^{n+1/2} to the grid cell corners
 */

#ifdef MHD
  integrate_emf3_corner(pG);

/*--- Step 10b -----------------------------------------------------------------
 * Update the interface magnetic fields using CT for a full time step.
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
      dtodx2 = pG->dt/(ri[i]*pG->dx2);
#endif
      pG->B1i[ks][j][i] -= dtodx2*(emf3[j+1][i  ] - emf3[j][i]);
      pG->B2i[ks][j][i] += dtodx1*(emf3[j  ][i+1] - emf3[j][i]);
    }
#ifdef CYLINDRICAL
    dtodx2 = pG->dt/(ri[ie+1]*pG->dx2);
#endif
    pG->B1i[ks][j][ie+1] -= dtodx2*(emf3[j+1][ie+1] - emf3[j][ie+1]);
  }
  for (i=is; i<=ie; i++) {
    pG->B2i[ks][je+1][i] += dtodx1*(emf3[je+1][i+1] - emf3[je+1][i]);
  }
#endif /* MHD */

/*=== STEP 11: Add source terms for a full timestep using n+1/2 states =======*/

/*--- Step 11a -----------------------------------------------------------------
 * ADD GEOMETRIC SOURCE TERMS.
 */

#ifdef CYLINDRICAL
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
      hdtodx2 = hdt/(r[i]*pG->dx2);

      /* CALCULATE d AT TIME n+1/2 */
      dhalf[j][i] = pG->U[ks][j][i].d
        - hdtodx1*(rsf*x1Flux[j  ][i+1].d - lsf*x1Flux[j][i].d)
        - hdtodx2*(    x2Flux[j+1][i  ].d -     x2Flux[j][i].d);

      /* CALCULATE M2 AT TIME n+1/2 */
      M2h = pG->U[ks][j][i].M2
        - hdtodx1*(SQR(rsf)*x1Flux[j][i+1].My - SQR(lsf)*x1Flux[j][i].My)
        - hdtodx2*(         x2Flux[j+1][i].Mx -          x2Flux[j][i].Mx);

      /* ADD SOURCE TERM FOR FIXED GRAVITATIONAL POTENTIAL FOR 0.5*dt */
      if (StaticGravPot != NULL){
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);
        M2h -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;
      }

      /* COMPUTE GEOMETRIC SOURCE TERM AT TIME n+1/2 */
      geom_src[j][i] = SQR(M2h)/dhalf[j][i];
#ifdef MHD
      B2ch = 0.5*(B2_x2Face[j][i] + B2_x2Face[j+1][i]);
      geom_src[j][i] -= SQR(B2ch);
#endif
#ifdef ISOTHERMAL
      geom_src[j][i] += Iso_csound2*dhalf[j][i];
#ifdef MHD
      B1ch = 0.5*(rsf*B1_x1Face[j][i+1] + lsf*B1_x1Face[j][i]);
      B3ch = pG->U[ks][j][i].B3c
        - hdtodx1*(rsf*x1Flux[j  ][i+1].Bz - lsf*x1Flux[j][i].Bz)
        - hdtodx2*(    x2Flux[j+1][i  ].By -     x2Flux[j][i].By);
      geom_src[j][i] += 0.5*(SQR(B1ch)+SQR(B2ch)+SQR(B3ch));
#endif
#else /* ISOTHERMAL */
      Pavgh = 0.5*(lsf*x1Flux[j][i].Pflux + rsf*x1Flux[j][i+1].Pflux);
      geom_src[j][i] += Pavgh;
#endif
      geom_src[j][i] /= x1vc(pG,i);

      /* ADD TIME-CENTERED GEOMETRIC SOURCE TERM FOR FULL dt */
      pG->U[ks][j][i].M1 += pG->dt*geom_src[j][i];
    }
  }
#endif /* CYLINDRICAL */

/*--- Step 11a -----------------------------------------------------------------
 * Add gravitational (or shearing box) source terms.
 * Remember with the 2D shearing box (1,2,3) = (x,z,y)
 *   A Crank-Nicholson update is used for shearing box terms.
 *   The energy source terms computed at cell faces are averaged to improve
 * conservation of total energy.
 *    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

#ifdef SHEARING_BOX
  fact = om_dt/(2. + (2.-qshear)*om_dt*om_dt);
  qom = qshear*Omega_0;
  for(j=js; j<=je; j++) {
    for(i=is; i<=ie; i++) {
      cc_pos(pG,i,j,ks,&x1,&x2,&x3);

/* Store the current state */
      M1n  = pG->U[ks][j][i].M1;
#ifdef FARGO
      if (ShBoxCoord==xy) dM2n = pG->U[ks][j][i].M2;
      if (ShBoxCoord==xz) dM3n = pG->U[ks][j][i].M3;
#else
      if (ShBoxCoord==xy) dM2n = pG->U[ks][j][i].M2 + qom*x1*pG->U[ks][j][i].d;
      if (ShBoxCoord==xz) dM3n = pG->U[ks][j][i].M3 + qom*x1*pG->U[ks][j][i].d;
#endif

/* Calculate the flux for the y-momentum fluctuation (M3 in 2D) */
      if (ShBoxCoord==xy){
        frx1_dM2 = x1Flux[j][i+1].My;
        flx1_dM2 = x1Flux[j][i  ].My;
        frx2_dM2 = x2Flux[j+1][i].Mx;
        flx2_dM2 = x2Flux[j  ][i].Mx;
      }
      if (ShBoxCoord==xz){
        frx1_dM3 = x1Flux[j][i+1].Mz;
        flx1_dM3 = x1Flux[j][i  ].Mz;
        frx2_dM3 = x2Flux[j+1][i].My;
        flx2_dM3 = x2Flux[j  ][i].My;
      }
#ifndef FARGO
      if (ShBoxCoord==xy){
        frx1_dM2 += qom*(x1+0.5*pG->dx1)*x1Flux[j][i+1].d;
        flx1_dM2 += qom*(x1-0.5*pG->dx1)*x1Flux[j][i  ].d;
        frx2_dM2 += qom*(x1            )*x2Flux[j+1][i].d;
        flx2_dM2 += qom*(x1            )*x2Flux[j  ][i].d;
      }
      if (ShBoxCoord==xz){
        frx1_dM3 += qom*(x1+0.5*pG->dx1)*x1Flux[j][i+1].d;
        flx1_dM3 += qom*(x1-0.5*pG->dx1)*x1Flux[j][i  ].d;
        frx2_dM3 += qom*(x1            )*x2Flux[j+1][i].d;
        flx2_dM3 += qom*(x1            )*x2Flux[j  ][i].d;
      }
#endif

/* evolve M1n and dM3n by dt/2 using flux gradients */
      M1e = M1n - hdtodx1*(x1Flux[j][i+1].Mx - x1Flux[j][i].Mx)
                - hdtodx2*(x2Flux[j+1][i].Mz - x2Flux[j][i].Mz);

      if (ShBoxCoord==xy){
        dM2e = dM2n - hdtodx1*(frx1_dM2 - flx1_dM2)
                    - hdtodx2*(frx2_dM2 - flx2_dM2);
      }
      if (ShBoxCoord==xz){
        dM3e = dM3n - hdtodx1*(frx1_dM3 - flx1_dM3)
                    - hdtodx2*(frx2_dM3 - flx2_dM3);
      }

#ifdef FEEDBACK
      M1e -= 0.5*pG->Coup[ks][j][i].fb1;
      if (ShBoxCoord==xy) dM2e -= 0.5*pG->Coup[ks][j][i].fb2;
      if (ShBoxCoord==xz) dM3e -= 0.5*pG->Coup[ks][j][i].fb3;
#endif

/* Update the 1- and 2-momentum (or 1- and 3-momentum in XZ 2D shearing box)
 * for the Coriolis and tidal potential source terms using a Crank-Nicholson
 * discretization for the momentum fluctuation equation. */

      if (ShBoxCoord==xy){
        pG->U[ks][j][i].M1 += (4.0*dM2e + 2.0*(qshear-2.)*om_dt*M1e)*fact;
        pG->U[ks][j][i].M2 += 2.0*(qshear-2.)*(M1e + om_dt*dM2e)*fact;
#ifndef FARGO
        pG->U[ks][j][i].M2 -=0.5*qshear*om_dt*(x1Flux[j][i].d+x1Flux[j][i+1].d);
#endif
      }
      if (ShBoxCoord==xz){
        pG->U[ks][j][i].M1 += (4.0*dM3e + 2.0*(qshear-2.)*om_dt*M1e)*fact;
        pG->U[ks][j][i].M3 += 2.0*(qshear-2.)*(M1e + om_dt*dM3e)*fact;
#ifndef FARGO
        pG->U[ks][j][i].M3 -=0.5*qshear*om_dt*(x1Flux[j][i].d+x1Flux[j][i+1].d);
#endif
      }

/* Update the energy for a fixed potential.
 * This update is identical to non-SHEARING_BOX below  */

      phic = (*ShearingBoxPot)((x1            ),x2,x3);
      phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
      phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E -= dtodx1*(x1Flux[j][i  ].d*(phic - phil) +
                                   x1Flux[j][i+1].d*(phir - phic));
#endif

      phir = (*ShearingBoxPot)(x1,(x2+0.5*pG->dx2),x3);
      phil = (*ShearingBoxPot)(x1,(x2-0.5*pG->dx2),x3);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E -= dtodx2*(x2Flux[j  ][i].d*(phic - phil) +
                                   x2Flux[j+1][i].d*(phir - phic));
#endif
    }
  }

#endif /* SHEARING_BOX */

  if (StaticGravPot != NULL){
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phic = (*StaticGravPot)((x1            ),x2,x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

#ifdef CYLINDRICAL
        g = (*x1GravAcc)(x1vc(pG,i),x2,x3);
        rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
        dtodx2 = pG->dt/(r[i]*pG->dx2);
        pG->U[ks][j][i].M1 -= pG->dt*dhalf[j][i]*g;
#else
        pG->U[ks][j][i].M1 -= dtodx1*dhalf[j][i]*(phir-phil);
#endif

#ifndef BAROTROPIC
        pG->U[ks][j][i].E -= dtodx1*(lsf*x1Flux[j][i  ].d*(phic - phil) +
                                     rsf*x1Flux[j][i+1].d*(phir - phic));
#endif
        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        pG->U[ks][j][i].M2 -= dtodx2*dhalf[j][i]*(phir-phil);

#ifndef BAROTROPIC
        pG->U[ks][j][i].E -= dtodx2*(x2Flux[j  ][i].d*(phic - phil) +
                                     x2Flux[j+1][i].d*(phir - phic));
#endif
      }
    }
  }


/*--- Step 11b -----------------------------------------------------------------
 * Add source terms for self-gravity.
 * A flux correction using Phi^{n+1} in the main loop is required to make
 * the source terms 2nd order: see selfg_flux_correction().
 */

#ifdef SELF_GRAVITY
/* Add fluxes and source terms due to (d/dx1) terms  */

  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      phic = pG->Phi[ks][j][i];
      phil = 0.5*(pG->Phi[ks][j][i-1] + pG->Phi[ks][j][i  ]);
      phir = 0.5*(pG->Phi[ks][j][i  ] + pG->Phi[ks][j][i+1]);

/* gx and gy centered at L and R x1-faces */
      gxl = (pG->Phi[ks][j][i-1] - pG->Phi[ks][j][i  ])/(pG->dx1);
      gxr = (pG->Phi[ks][j][i  ] - pG->Phi[ks][j][i+1])/(pG->dx1);

      gyl = 0.25*((pG->Phi[ks][j-1][i-1] - pG->Phi[ks][j+1][i-1]) +
                  (pG->Phi[ks][j-1][i  ] - pG->Phi[ks][j+1][i  ]) )/(pG->dx2);
      gyr = 0.25*((pG->Phi[ks][j-1][i  ] - pG->Phi[ks][j+1][i  ]) +
                  (pG->Phi[ks][j-1][i+1] - pG->Phi[ks][j+1][i+1]) )/(pG->dx2);

/* momentum fluxes in x1.  2nd term is needed only if Jean's swindle used */
      flux_m1l = 0.5*(gxl*gxl-gyl*gyl)/four_pi_G + grav_mean_rho*phil;
      flux_m1r = 0.5*(gxr*gxr-gyr*gyr)/four_pi_G + grav_mean_rho*phir;

      flux_m2l = gxl*gyl/four_pi_G;
      flux_m2r = gxr*gyr/four_pi_G;

/* Update momenta and energy with d/dx1 terms  */
      pG->U[ks][j][i].M1 -= dtodx1*(flux_m1r - flux_m1l);
      pG->U[ks][j][i].M2 -= dtodx1*(flux_m2r - flux_m2l);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E -= dtodx1*(x1Flux[j][i  ].d*(phic - phil) +
                                   x1Flux[j][i+1].d*(phir - phic));
#endif
    }
  }

/* Add fluxes and source terms due to (d/dx2) terms  */

  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      phic = pG->Phi[ks][j][i];
      phil = 0.5*(pG->Phi[ks][j-1][i] + pG->Phi[ks][j  ][i]);
      phir = 0.5*(pG->Phi[ks][j  ][i] + pG->Phi[ks][j+1][i]);

/* gx and gy centered at L and R x2-faces */
      gxl = 0.25*((pG->Phi[ks][j-1][i-1] - pG->Phi[ks][j-1][i+1]) +
                  (pG->Phi[ks][j  ][i-1] - pG->Phi[ks][j  ][i+1]) )/(pG->dx1);
      gxr = 0.25*((pG->Phi[ks][j  ][i-1] - pG->Phi[ks][j  ][i+1]) +
                  (pG->Phi[ks][j+1][i-1] - pG->Phi[ks][j+1][i+1]) )/(pG->dx1);

      gyl = (pG->Phi[ks][j-1][i] - pG->Phi[ks][j  ][i])/(pG->dx2);
      gyr = (pG->Phi[ks][j  ][i] - pG->Phi[ks][j+1][i])/(pG->dx2);

/* momentum fluxes in x2.  2nd term is needed only if Jean's swindle used */
      flux_m1l = gyl*gxl/four_pi_G;
      flux_m1r = gyr*gxr/four_pi_G;

      flux_m2l = 0.5*(gyl*gyl-gxl*gxl)/four_pi_G + grav_mean_rho*phil;
      flux_m2r = 0.5*(gyr*gyr-gxr*gxr)/four_pi_G + grav_mean_rho*phir;

/* Update momenta and energy with d/dx2 terms  */
      pG->U[ks][j][i].M1 -= dtodx2*(flux_m1r - flux_m1l);
      pG->U[ks][j][i].M2 -= dtodx2*(flux_m2r - flux_m2l);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E -= dtodx2*(x2Flux[j  ][i].d*(phic - phil) +
                                   x2Flux[j+1][i].d*(phir - phic));
#endif
    }
  }

/* Save mass fluxes in Grid structure for source term correction in main loop */
  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      pG->x1MassFlux[ks][j][i] = x1Flux[j][i].d;
      pG->x2MassFlux[ks][j][i] = x2Flux[j][i].d;
    }
  }
#endif /* SELF_GRAVITY */

/*--- Step 11c -----------------------------------------------------------------
 * Add source terms for optically thin cooling
 */

#ifndef BAROTROPIC
  if (CoolingFunc != NULL){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        coolf = (*CoolingFunc)(dhalf[j][i],phalf[j][i],pG->dt);
        pG->U[ks][j][i].E -= pG->dt*coolf;
      }
    }
  }
#endif /* BAROTROPIC */

/*--- Step 11d -----------------------------------------------------------------
 * Add source terms for particle feedback
 */

#ifdef FEEDBACK
  for (j=js; j<=je; j++)
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].M1 -= pG->Coup[ks][j][i].fb1;
      pG->U[ks][j][i].M2 -= pG->Coup[ks][j][i].fb2;
      pG->U[ks][j][i].M3 -= pG->Coup[ks][j][i].fb3;
#ifndef BAROTROPIC
      pG->U[ks][j][i].E += pG->Coup[ks][j][i].Eloss;
      pG->Coup[ks][j][i].Eloss *= dt1; /* for history output purpose */
      pG->Eloss[ks][j][i] *= dt1; /* for history output purpose */
#endif
    }
#endif

/*=== STEP 12: Update cell-centered values for a full timestep ===============*/

/*--- Step 12a -----------------------------------------------------------------
 * Update cell-centered variables in pG using 2D x1-fluxes
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
      pG->U[ks][j][i].d  -= dtodx1*(rsf*x1Flux[j][i+1].d  - lsf*x1Flux[j][i].d );
      pG->U[ks][j][i].M1 -= dtodx1*(rsf*x1Flux[j][i+1].Mx - lsf*x1Flux[j][i].Mx);
      pG->U[ks][j][i].M2 -= dtodx1*(SQR(rsf)*x1Flux[j][i+1].My - SQR(lsf)*x1Flux[j][i].My);
      pG->U[ks][j][i].M3 -= dtodx1*(rsf*x1Flux[j][i+1].Mz - lsf*x1Flux[j][i].Mz);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E  -= dtodx1*(rsf*x1Flux[j][i+1].E  - lsf*x1Flux[j][i].E );
#endif /* BAROTROPIC */
#ifdef MHD
      pG->U[ks][j][i].B2c -= dtodx1*(x1Flux[j][i+1].By - x1Flux[j][i].By);
      pG->U[ks][j][i].B3c -= dtodx1*(rsf*x1Flux[j][i+1].Bz - lsf*x1Flux[j][i].Bz);
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        pG->U[ks][j][i].s[n] -= dtodx1*(rsf*x1Flux[j][i+1].s[n] 
                                      - lsf*x1Flux[j][i  ].s[n]);
#endif
    }
  }

/*--- Step 12b -----------------------------------------------------------------
 * Update cell-centered variables in pG using 2D x2-fluxes
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
      dtodx2 = pG->dt/(r[i]*pG->dx2);
#endif
      pG->U[ks][j][i].d  -= dtodx2*(x2Flux[j+1][i].d  - x2Flux[j][i].d );
      pG->U[ks][j][i].M1 -= dtodx2*(x2Flux[j+1][i].Mz - x2Flux[j][i].Mz);
      pG->U[ks][j][i].M2 -= dtodx2*(x2Flux[j+1][i].Mx - x2Flux[j][i].Mx);
      pG->U[ks][j][i].M3 -= dtodx2*(x2Flux[j+1][i].My - x2Flux[j][i].My);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E  -= dtodx2*(x2Flux[j+1][i].E  - x2Flux[j][i].E );
#endif /* BAROTROPIC */
#ifdef MHD
      pG->U[ks][j][i].B3c -= dtodx2*(x2Flux[j+1][i].By - x2Flux[j][i].By);
      pG->U[ks][j][i].B1c -= dtodx2*(x2Flux[j+1][i].Bz - x2Flux[j][i].Bz);
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        pG->U[ks][j][i].s[n] -= dtodx2*(x2Flux[j+1][i].s[n] 
                                         - x2Flux[j  ][i].s[n]);
#endif
    }
  }

/*--- Step 12c: Not needed in 2D ---*/
/*--- Step 12d -----------------------------------------------------------------
 * Set cell centered magnetic fields to average of updated face centered fields.
 */

#ifdef MHD
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
      pG->U[ks][j][i].B1c =0.5*(lsf*pG->B1i[ks][j][i] + rsf*pG->B1i[ks][j][i+1]);
      pG->U[ks][j][i].B2c =0.5*(    pG->B2i[ks][j][i] +     pG->B2i[ks][j+1][i]);
/* Set the 3-interface magnetic field equal to the cell center field. */
      pG->B3i[ks][j][i] = pG->U[ks][j][i].B3c;
    }
  }
#endif /* MHD */

#ifdef STATIC_MESH_REFINEMENT
/*--- Step 12e -----------------------------------------------------------------
 * With SMR, store fluxes at boundaries of child and parent grids.
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
  if ((B1_x1Face = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((B2_x2Face = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* MHD */

  if ((U1d= (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ul = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ur = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W  = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wl = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wr = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((Ul_x1Face=(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;
  if ((Ur_x1Face=(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;
  if ((Ul_x2Face=(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;
  if ((Ur_x2Face=(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;

  if ((x1Flux   =(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;
  if ((x2Flux   =(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;

#ifndef CYLINDRICAL
#ifndef MHD
#ifndef PARTICLES
  if((StaticGravPot != NULL) || (CoolingFunc != NULL))
#endif
#endif
#endif
  {
  if ((dhalf = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((phalf = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  }

  /* DATA STRUCTURES FOR CYLINDRICAL COORDINATES */
#ifdef CYLINDRICAL
  if ((geom_src = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL) 
    goto on_error;
#endif

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

  if (U1d      != NULL) free(U1d);
  if (Ul       != NULL) free(Ul);
  if (Ur       != NULL) free(Ur);
  if (W        != NULL) free(W);
  if (Wl       != NULL) free(Wl);
  if (Wr       != NULL) free(Wr);

  if (Ul_x1Face != NULL) free_2d_array(Ul_x1Face);
  if (Ur_x1Face != NULL) free_2d_array(Ur_x1Face);
  if (Ul_x2Face != NULL) free_2d_array(Ul_x2Face);
  if (Ur_x2Face != NULL) free_2d_array(Ur_x2Face);
  if (x1Flux    != NULL) free_2d_array(x1Flux);
  if (x2Flux    != NULL) free_2d_array(x2Flux);
  if (dhalf     != NULL) free_2d_array(dhalf);
  if (phalf     != NULL) free_2d_array(phalf);

  /* DATA STRUCTURES FOR CYLINDRICAL COORDINATES */
#ifdef CYLINDRICAL
  if (geom_src  != NULL) free_2d_array(geom_src);
#endif

  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/*! \fn static void integrate_emf3_corner(GridS *pG)
 *  \brief Corner EMF terms in the upwind CT method */
#ifdef MHD
static void integrate_emf3_corner(GridS *pG)
{
  int i,is,ie,j,js,je;
  Real emf_l1, emf_r1, emf_l2, emf_r2;
  Real rsf=1.0,lsf=1.0;

  is = pG->is;   ie = pG->ie;
  js = pG->js;   je = pG->je;

/* NOTE: The x1-Flux of B2 is -E3.  The x2-Flux of B1 is +E3. */
  for (j=js-1; j<=je+2; j++) {
    for (i=is-1; i<=ie+2; i++) {
#ifdef CYLINDRICAL
      rsf = pG->ri[i]/pG->r[i];  lsf = pG->ri[i]/pG->r[i-1];
#endif
      if (x1Flux[j-1][i].d > 0.0) {
	emf_l2 = -x1Flux[j-1][i].By
	  + (x2Flux[j][i-1].Bz - emf3_cc[j-1][i-1])*lsf;
      }
      else if (x1Flux[j-1][i].d < 0.0) {
	emf_l2 = -x1Flux[j-1][i].By
	  + (x2Flux[j][i].Bz - emf3_cc[j-1][i])*rsf;

      } else {
	emf_l2 = -x1Flux[j-1][i].By
	  + 0.5*((x2Flux[j][i-1].Bz - emf3_cc[j-1][i-1])*lsf + 
		 (x2Flux[j][i  ].Bz - emf3_cc[j-1][i  ])*rsf );
      }

      if (x1Flux[j][i].d > 0.0) {
	emf_r2 = -x1Flux[j][i].By 
	  + (x2Flux[j][i-1].Bz - emf3_cc[j][i-1])*lsf;
      }
      else if (x1Flux[j][i].d < 0.0) {
	emf_r2 = -x1Flux[j][i].By 
	  + (x2Flux[j][i].Bz - emf3_cc[j][i])*rsf;

      } else {
	emf_r2 = -x1Flux[j][i].By 
	  + 0.5*((x2Flux[j][i-1].Bz - emf3_cc[j][i-1])*lsf + 
		 (x2Flux[j][i  ].Bz - emf3_cc[j][i  ])*rsf );
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

#endif /* CTU_INTEGRATOR */
