#include "../copyright.h"
/*============================================================================*/
/*! \file conduction.c
 *  \brief Adds explicit thermal conduction term to the energy equation,
 *      dE/dt = Div(Q)
 *
 *   where 
 *    - Q = kappa_iso Grad(T) + kappa_aniso([b Dot Grad(T)]b) = heat flux
 *    -    T = (P/d)*(mbar/k_B) = temperature
 *    -    b = magnetic field unit vector
 *
 *   Here 
 *    - kappa_iso   is the coeffcient for   isotropic conduction
 *    -   kappa_aniso is the coeffcient for anisotropic conduction
 *
 *   The heat flux Q is calculated by calls to HeatFlux_* functions.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - conduction() - updates energy equation with thermal conduction
 * - conduction_init() - allocates memory needed
 * - conduction_destruct() - frees memory used */
/*============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef THERMAL_CONDUCTION

#ifdef BAROTROPIC
#error : Thermal conduction requires an adiabatic EOS
#endif

/* Arrays for the temperature and heat fluxes */
static Real ***Temp=NULL;
static Real3Vect ***Q=NULL;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   HeatFlux_iso   - computes   isotropic heat flux
 *   HeatFlux_aniso - computes anisotropic heat flux
 *============================================================================*/

void HeatFlux_iso(DomainS *pD);
void HeatFlux_aniso(DomainS *pD);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void conduction(DomainS *pD)
 *  \brief Explicit thermal conduction
 */
void conduction(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, jl, ju, js = pG->js, je = pG->je;
  int k, kl, ku, ks = pG->ks, ke = pG->ke;
  Real dtodx1=pG->dt/pG->dx1, dtodx2=0.0, dtodx3=0.0;

  if (pG->Nx[1] > 1){
    jl = js - 1;
    ju = je + 1;
    dtodx2 = pG->dt/pG->dx2;
  } else {
    jl = js;
    ju = je;
  }
  if (pG->Nx[2] > 1){
    kl = ks - 1;
    ku = ke + 1;
    dtodx3 = pG->dt/pG->dx3;
  } else {
    kl = ks;
    ku = ke;
  }

/* Zero heat heat flux array; compute temperature at cell centers.  Temperature
 * is dimensionless (in units of [mbar/k_B]).  For cgs units, the value of
 * kappa_iso or kappa_aniso must contain a factor (be multiplied by) mbar/k_B.
 */

  for (k=kl; k<=ku; k++) {
  for (j=jl; j<=ju; j++) {
  for (i=is-1; i<=ie+1; i++) {

    Q[k][j][i].x = 0.0;
    Q[k][j][i].y = 0.0;
    Q[k][j][i].z = 0.0;

    Temp[k][j][i] = pG->U[k][j][i].E - (0.5/pG->U[k][j][i].d)*
      (SQR(pG->U[k][j][i].M1) +SQR(pG->U[k][j][i].M2) +SQR(pG->U[k][j][i].M3));
#ifdef MHD
    Temp[k][j][i] -= (0.5)*(SQR(pG->U[k][j][i].B1c) +
      SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
#endif
    Temp[k][j][i] *= (Gamma_1/pG->U[k][j][i].d);

  }}}

/* Compute isotropic and anisotropic heat fluxes.  Heat fluxes and temperature
 * are global variables in this file. */

  if (kappa_iso > 0.0)   HeatFlux_iso(pD);
  if (kappa_aniso > 0.0) HeatFlux_aniso(pD);

/* Update energy using x1-fluxes */

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[k][j][i].E += dtodx1*(Q[k][j][i+1].x - Q[k][j][i].x);
    }
  }}

/* Update energy using x2-fluxes */

  if (pG->Nx[1] > 1){
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E += dtodx2*(Q[k][j+1][i].y - Q[k][j][i].y);
      }
    }}
  }

/* Update energy using x3-fluxes */

  if (pG->Nx[2] > 1){
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E += dtodx3*(Q[k+1][j][i].z - Q[k][j][i].z);
      }
    }}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void HeatFlux_iso(DomainS *pD)
 *  \brief Calculate heat fluxes with isotropic conduction
 */

void HeatFlux_iso(DomainS *pD)
{ 
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;

/* Add heat fluxes in 1-direction */ 

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie+1; i++) {
      Q[k][j][i].x += kappa_iso*(Temp[k][j][i] - Temp[k][j][i-1])/pG->dx1;
    }
  }}

/* Add heat fluxes in 2-direction */ 

  if (pG->Nx[1] > 1) {
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie; i++) {
        Q[k][j][i].y += kappa_iso*(Temp[k][j][i] - Temp[k][j-1][i])/pG->dx2;
      }
    }}
  }

/* Add heat fluxes in 3-direction */

  if (pG->Nx[2] > 1) {
    for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        Q[k][j][i].z += kappa_iso*(Temp[k][j][i] - Temp[k-1][j][i])/pG->dx3;
      }
    }}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void HeatFlux_aniso(DomainS *pD)
 *  \brief Calculate heat fluxes with anisotropic conduction
 */

void HeatFlux_aniso(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real Bx,By,Bz,B02,dTc,dTl,dTr,lim_slope,dTdx,dTdy,dTdz,bDotGradT;

#ifdef MHD
  if (pD->Nx[1] == 1) return;  /* problem must be at least 2D */

/* Compute heat fluxes in 1-direction  --------------------------------------*/

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie+1; i++) {

/* Monotonized temperature difference dT/dy */
      dTr = 0.5*((Temp[k][j+1][i-1] + Temp[k][j+1][i]) -
                 (Temp[k][j  ][i-1] + Temp[k][j  ][i]));
      dTl = 0.5*((Temp[k][j  ][i-1] + Temp[k][j  ][i]) -
                 (Temp[k][j-1][i-1] + Temp[k][j-1][i]));
      dTc = dTr + dTl;

      dTdy = 0.0;
      if (dTl*dTr > 0.0) {
        lim_slope = MIN(fabs(dTl),fabs(dTr));
        dTdy = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx2;
      }

/* Monotonized temperature difference dT/dz, 3D problem ONLY */
      if (pD->Nx[2] > 1) {
        dTr = 0.5*((Temp[k+1][j][i-1] + Temp[k+1][j][i]) -
                   (Temp[k  ][j][i-1] + Temp[k  ][j][i]));
        dTl = 0.5*((Temp[k  ][j][i-1] + Temp[k  ][j][i]) -
                   (Temp[k-1][j][i-1] + Temp[k-1][j][i]));
        dTc = dTr + dTl; 

        dTdz = 0.0;
        if (dTl*dTr > 0.0) {
          lim_slope = MIN(fabs(dTl),fabs(dTr));
          dTdz = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx3;
        }
      }

/* Add flux at x1-interface, 2D PROBLEM */

      if (pD->Nx[2] == 1) {
        By = 0.5*(pG->U[k][j][i-1].B2c + pG->U[k][j][i].B2c);
        B02 = SQR(pG->B1i[k][j][i]) + SQR(By);
        B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */
        bDotGradT = pG->B1i[k][j][i]*(Temp[k][j][i]-Temp[k][j][i-1])/pG->dx1
           + By*dTdy;
        Q[k][j][i].x += kappa_aniso*(pG->B1i[k][j][i]*bDotGradT)/B02;

/* Add flux at x1-interface, 3D PROBLEM */

      } else {
        By = 0.5*(pG->U[k][j][i-1].B2c + pG->U[k][j][i].B2c);
        Bz = 0.5*(pG->U[k][j][i-1].B3c + pG->U[k][j][i].B3c);
        B02 = SQR(pG->B1i[k][j][i]) + SQR(By) + SQR(Bz);
        B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */
        bDotGradT = pG->B1i[k][j][i]*(Temp[k][j][i]-Temp[k][j][i-1])/pG->dx1
           + By*dTdy + Bz*dTdz;
        Q[k][j][i].x += kappa_aniso*(pG->B1i[k][j][i]*bDotGradT)/B02;
      }
    }
  }}

/* Compute heat fluxes in 2-direction  --------------------------------------*/

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie; i++) {

/* Monotonized temperature difference dT/dx */
      dTr = 0.5*((Temp[k][j-1][i+1] + Temp[k][j][i+1]) -
                 (Temp[k][j-1][i  ] + Temp[k][j][i  ]));
      dTl = 0.5*((Temp[k][j-1][i  ] + Temp[k][j][i  ]) -
                 (Temp[k][j-1][i-1] + Temp[k][j][i-1]));
      dTc = dTr + dTl;

      dTdx = 0.0;
      if (dTl*dTr > 0.0) {
        lim_slope = MIN(fabs(dTl),fabs(dTr));
        dTdx = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx1;
      }

/* Monotonized temperature difference dT/dz, 3D problem ONLY */
      if (pD->Nx[2] > 1) {
        dTr = 0.5*((Temp[k+1][j-1][i] + Temp[k+1][j][i]) -
                   (Temp[k  ][j-1][i] + Temp[k  ][j][i]));
        dTl = 0.5*((Temp[k  ][j-1][i] + Temp[k  ][j][i]) -
                   (Temp[k-1][j-1][i] + Temp[k-1][j][i]));
        dTc = dTr + dTl;

        dTdz = 0.0;
        if (dTl*dTr > 0.0) {
          lim_slope = MIN(fabs(dTl),fabs(dTr));
          dTdz = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx3;
        }
      }

/* Add flux at x2-interface, 2D PROBLEM */

      if (pD->Nx[2] == 1) {
        Bx = 0.5*(pG->U[k][j-1][i].B1c + pG->U[k][j][i].B1c);
        B02 = SQR(Bx) + SQR(pG->B2i[k][j][i]);
        B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */

        bDotGradT = pG->B2i[k][j][i]*(Temp[k][j][i]-Temp[k][j-1][i])/pG->dx2
           + Bx*dTdx;
        Q[k][j][i].y += kappa_aniso*(pG->B2i[k][j][i]*bDotGradT)/B02;

/* Add flux at x2-interface, 3D PROBLEM */

      } else {
        Bx = 0.5*(pG->U[k][j-1][i].B1c + pG->U[k][j][i].B1c);
        Bz = 0.5*(pG->U[k][j-1][i].B3c + pG->U[k][j][i].B3c);
        B02 = SQR(Bx) + SQR(pG->B2i[k][j][i]) + SQR(Bz);
        B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */
        bDotGradT = pG->B2i[k][j][i]*(Temp[k][j][i]-Temp[k][j-1][i])/pG->dx2
           + Bx*dTdx + Bz*dTdz;
        Q[k][j][i].y += kappa_aniso*(pG->B2i[k][j][i]*bDotGradT)/B02;
      }
    }
  }}

/* Compute heat fluxes in 3-direction, 3D problem ONLY  ---------------------*/

  if (pD->Nx[2] > 1) {
    for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {

/* Monotonized temperature difference dT/dx */
        dTr = 0.5*((Temp[k-1][j][i+1] + Temp[k][j][i+1]) -
                   (Temp[k-1][j][i  ] + Temp[k][j][i  ]));
        dTl = 0.5*((Temp[k-1][j][i  ] + Temp[k][j][i  ]) -
                   (Temp[k-1][j][i-1] + Temp[k][j][i-1]));
        dTc = dTr + dTl;

        dTdx = 0.0; 
        if (dTl*dTr > 0.0) {
          lim_slope = MIN(fabs(dTl),fabs(dTr));
          dTdx = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx1;
        }

/* Monotonized temperature difference dT/dy */
        dTr = 0.5*((Temp[k-1][j+1][i] + Temp[k][j+1][i]) -
                   (Temp[k-1][j  ][i] + Temp[k][j  ][i]));
        dTl = 0.5*((Temp[k-1][j  ][i] + Temp[k][j  ][i]) -
                   (Temp[k-1][j-1][i] + Temp[k][j-1][i]));
        dTc = dTr + dTl;

        dTdy = 0.0; 
        if (dTl*dTr > 0.0) {
          lim_slope = MIN(fabs(dTl),fabs(dTr));
          dTdy = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx2;
        }

/* Add flux at x3-interface, 3D PROBLEM */

        Bx = 0.5*(pG->U[k-1][j][i].B1c + pG->U[k][j][i].B1c);
        By = 0.5*(pG->U[k-1][j][i].B2c + pG->U[k][j][i].B2c);
        B02 = SQR(Bx) + SQR(By) + SQR(pG->B3i[k][j][i]);
        B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */
        bDotGradT = pG->B3i[k][j][i]*(Temp[k][j][i]-Temp[k-1][j][i])/pG->dx3
           + Bx*dTdx + By*dTdy;
        Q[k][j][i].z += kappa_aniso*(pG->B3i[k][j][i]*bDotGradT)/B02;
      }
    }}
  }
#endif /* MHD */

  return;
}


/*----------------------------------------------------------------------------*/
/*! \fn void conduction_init(MeshS *pM) 
 *  \brief Allocate temporary arrays
 */

void conduction_init(MeshS *pM)
{
  int nl,nd,size1=0,size2=0,size3=0,Nx1,Nx2,Nx3;

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

  Nx1 = size1 + 2*nghost;

  if (pM->Nx[1] > 1){
    Nx2 = size2 + 2*nghost;
  } else {
    Nx2 = size2;
  }

  if (pM->Nx[2] > 1){
    Nx3 = size3 + 2*nghost;
  } else {
    Nx3 = size3;
  }
  if ((Temp = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((Q = (Real3Vect***)calloc_3d_array(Nx3,Nx2,Nx1,sizeof(Real3Vect)))==NULL)
    goto on_error;
  return;

  on_error:
  conduction_destruct();
  ath_error("[conduct_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void conduction_destruct(void)
 *  \brief Free temporary arrays
 */

void conduction_destruct(void)
{
  if (Temp != NULL) free_3d_array(Temp);
  if (Q != NULL) free_3d_array(Q);
  return;
}
#endif /* THERMAL_CONDUCTION */
