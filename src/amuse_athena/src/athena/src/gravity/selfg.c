#include "../copyright.h"
/*============================================================================*/
/*! \file selfg.c
 *  \brief Contains functions to control solution of Poisson's equation for
 *   self-gravity.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - selfg_fc()   - 2nd order flux-corrections for self-gravity terms
 * - selfg_init() - sets pointer to appropriate self-gravity function
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef SELF_GRAVITY
/*----------------------------------------------------------------------------*/
/*! \fn void selfg_fc(DomainS *pD)
 *  \brief Adds flux-correction to make the integration algorithms for the
 *   source terms in the momentum and energy equations second-order.  
 *
 *   This
 *   requires subtracting 1/2 the source terms computed with the old potential,
 *   and adding 1/2 the source terms computed with the new potential.
 *
 *   The source terms for momentum are computed using the divergence of the
 *   gravitational stress tensor to conserve momentum exactly.
 *   - dM/dt = -Div(G);   G = (gg - 0.5g^2)/4\piG;   g=-Grad(Phi);
 *
 *   The source terms for the energy are added using the mass fluxes at cell
 *   faces, to improve conservation.
 *   - S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

void selfg_fc(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int dim=0;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;
  Real dtodx3 = pG->dt/pG->dx3;
  Real phic,phil,phir,phil_old,phir_old,dphic,dphil,dphir;
  Real gxl,gxr,gyl,gyr,gzl,gzr;
  Real flx_m1r,flx_m1l,flx_m2r,flx_m2l,flx_m3r,flx_m3l;
  
/* Calculate the dimensions  */
  if(pG->Nx[0] > 1) dim++;
  if(pG->Nx[1] > 1) dim++;
  if(pG->Nx[2] > 1) dim++;


/* The divergence of the gravitational stress tensor depends on the dimensions
 * of the problem.
 */

  switch(dim){
/*------------------------- 1D problem ---------------------------------------*/
  case 1:

/* Step 1 for 1D.  Add fluxes and source terms due to (d/dx1) terms  */

    for (i=is; i<=ie; i++){
      phic = pG->Phi[ks][js][i];
      phil = 0.5*(pG->Phi[ks][js][i-1] + pG->Phi[ks][js][i  ]);
      phir = 0.5*(pG->Phi[ks][js][i  ] + pG->Phi[ks][js][i+1]);
      phil_old = 0.5*(pG->Phi_old[ks][js][i-1] + pG->Phi_old[ks][js][i  ]);
      phir_old = 0.5*(pG->Phi_old[ks][js][i  ] + pG->Phi_old[ks][js][i+1]);

      dphic = phic - pG->Phi_old[ks][js][i];
      dphil = phil - phil_old;
      dphir = phir - phir_old;

/* momentum fluxes in x1.  gx centered at L and R x1-faces */
      gxl = (pG->Phi[ks][js][i-1] - pG->Phi[ks][js][i  ])/(pG->dx1);
      gxr = (pG->Phi[ks][js][i  ] - pG->Phi[ks][js][i+1])/(pG->dx1);

      flx_m1l = 0.5*(gxl*gxl)/four_pi_G + grav_mean_rho*phil;
      flx_m1r = 0.5*(gxr*gxr)/four_pi_G + grav_mean_rho*phir;

/*  subtract off momentum fluxes from old potential */
      gxl = (pG->Phi_old[ks][js][i-1] - pG->Phi_old[ks][js][i  ])/(pG->dx1);
      gxr = (pG->Phi_old[ks][js][i  ] - pG->Phi_old[ks][js][i+1])/(pG->dx1);

      flx_m1l -= 0.5*(gxl*gxl)/four_pi_G + grav_mean_rho*phil_old;
      flx_m1r -= 0.5*(gxr*gxr)/four_pi_G + grav_mean_rho*phir_old;

/* Update momenta and energy with d/dx1 terms  */
      pG->U[ks][js][i].M1 -= 0.5*dtodx1*(flx_m1r-flx_m1l);
#ifndef ISOTHERMAL
      pG->U[ks][js][i].E -=
         0.5*dtodx1*(pG->x1MassFlux[ks][js][i  ]*(dphic - dphil) +
                     pG->x1MassFlux[ks][js][i+1]*(dphir - dphic));
#endif
    }
    break;

/*------------------------- 2D problem ---------------------------------------*/
  case 2:

/* Step 1 for 2D.  Add fluxes and source terms due to (d/dx1) terms  */

    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        phic = pG->Phi[ks][j][i];
        phil = 0.5*(pG->Phi[ks][j][i-1] + pG->Phi[ks][j][i  ]);
        phir = 0.5*(pG->Phi[ks][j][i  ] + pG->Phi[ks][j][i+1]);
        phil_old = 0.5*(pG->Phi_old[ks][j][i-1] + pG->Phi_old[ks][j][i  ]);
        phir_old = 0.5*(pG->Phi_old[ks][j][i  ] + pG->Phi_old[ks][j][i+1]);

        dphic = phic - pG->Phi_old[ks][j][i];
        dphil = phil - phil_old;
        dphir = phir - phir_old;

/*  momentum fluxes in x1.  gx and gy centered at L and R x1-faces */
        gxl = (pG->Phi[ks][j][i-1] - pG->Phi[ks][j][i  ])/(pG->dx1);
        gxr = (pG->Phi[ks][j][i  ] - pG->Phi[ks][j][i+1])/(pG->dx1);

        gyl = (pG->Phi[ks][j-1][i-1] - pG->Phi[ks][j+1][i-1]) +
              (pG->Phi[ks][j-1][i  ] - pG->Phi[ks][j+1][i  ]);
        gyl *= (0.25/pG->dx2);

        gyr = (pG->Phi[ks][j-1][i  ] - pG->Phi[ks][j+1][i  ]) +
              (pG->Phi[ks][j-1][i+1] - pG->Phi[ks][j+1][i+1]);
        gyr *= (0.25/pG->dx2);

        flx_m1l = 0.5*(gxl*gxl-gyl*gyl)/four_pi_G + grav_mean_rho*phil;
        flx_m1r = 0.5*(gxr*gxr-gyr*gyr)/four_pi_G + grav_mean_rho*phir;

        flx_m2l = gxl*gyl/four_pi_G;
        flx_m2r = gxr*gyr/four_pi_G;

/*  subtract off momentum fluxes from old potential */
	gxl = (pG->Phi_old[ks][j][i-1] - pG->Phi_old[ks][j][i  ])/(pG->dx1);
	gxr = (pG->Phi_old[ks][j][i  ] - pG->Phi_old[ks][j][i+1])/(pG->dx1);

	gyl = (pG->Phi_old[ks][j-1][i-1] - pG->Phi_old[ks][j+1][i-1]) +
              (pG->Phi_old[ks][j-1][i  ] - pG->Phi_old[ks][j+1][i  ]);
        gyl *= (0.25/pG->dx2);

	gyr = (pG->Phi_old[ks][j-1][i  ] - pG->Phi_old[ks][j+1][i  ]) +
              (pG->Phi_old[ks][j-1][i+1] - pG->Phi_old[ks][j+1][i+1]);
        gyr *= (0.25/pG->dx2);

	flx_m1l -= 0.5*(gxl*gxl-gyl*gyl)/four_pi_G + grav_mean_rho*phil_old;
	flx_m1r -= 0.5*(gxr*gxr-gyr*gyr)/four_pi_G + grav_mean_rho*phir_old;

	flx_m2l -= gxl*gyl/four_pi_G;
	flx_m2r -= gxr*gyr/four_pi_G;

/* Update momenta and energy with d/dx1 terms  */
        pG->U[ks][j][i].M1 -= 0.5*dtodx1*(flx_m1r - flx_m1l);
        pG->U[ks][j][i].M2 -= 0.5*dtodx1*(flx_m2r - flx_m2l);
#ifndef ISOTHERMAL
        pG->U[ks][j][i].E -=
           0.5*dtodx1*(pG->x1MassFlux[ks][j][i  ]*(dphic - dphil) +
                       pG->x1MassFlux[ks][j][i+1]*(dphir - dphic));
#endif
      }
    }

/* Step 2 for 2D.  Add fluxes and source terms due to (d/dx2) terms  */

    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        phic = pG->Phi[ks][j][i];
        phil = 0.5*(pG->Phi[ks][j-1][i] + pG->Phi[ks][j  ][i]);
        phir = 0.5*(pG->Phi[ks][j  ][i] + pG->Phi[ks][j+1][i]);
        phil_old = 0.5*(pG->Phi_old[ks][j-1][i] + pG->Phi_old[ks][j  ][i]);
        phir_old = 0.5*(pG->Phi_old[ks][j  ][i] + pG->Phi_old[ks][j+1][i]);

        dphic = phic - pG->Phi[ks][j][i];
        dphil = phil - phil_old;
        dphir = phir - phir_old;

/*  momentum fluxes in x2.  gx and gy centered at L and R x2-faces */
        gxl = (pG->Phi[ks][j-1][i-1] - pG->Phi[ks][j-1][i+1]) +
              (pG->Phi[ks][j  ][i-1] - pG->Phi[ks][j  ][i+1]);
        gxl *= (0.25/pG->dx1);

        gxr = (pG->Phi[ks][j  ][i-1] - pG->Phi[ks][j  ][i+1]) +
              (pG->Phi[ks][j+1][i-1] - pG->Phi[ks][j+1][i+1]);
        gxr *= (0.25/pG->dx1);

        gyl = (pG->Phi[ks][j-1][i] - pG->Phi[ks][j  ][i])/(pG->dx2);
        gyr = (pG->Phi[ks][j  ][i] - pG->Phi[ks][j+1][i])/(pG->dx2);

        flx_m1l = gyl*gxl/four_pi_G;
        flx_m1r = gyr*gxr/four_pi_G;

        flx_m2l = 0.5*(gyl*gyl-gxl*gxl)/four_pi_G + grav_mean_rho*phil;
        flx_m2r = 0.5*(gyr*gyr-gxr*gxr)/four_pi_G + grav_mean_rho*phir;

/*  subtract off momentum fluxes from old potential */
        gxl = (pG->Phi_old[ks][j-1][i-1] - pG->Phi_old[ks][j-1][i+1]) +
              (pG->Phi_old[ks][j  ][i-1] - pG->Phi_old[ks][j  ][i+1]);
        gxl *= (0.25/pG->dx1);

        gxr = (pG->Phi_old[ks][j  ][i-1] - pG->Phi_old[ks][j  ][i+1]) +
              (pG->Phi_old[ks][j+1][i-1] - pG->Phi_old[ks][j+1][i+1]);
        gxr *= (0.25/pG->dx1);

        gyl = (pG->Phi_old[ks][j-1][i] - pG->Phi_old[ks][j  ][i])/(pG->dx2);
        gyr = (pG->Phi_old[ks][j  ][i] - pG->Phi_old[ks][j+1][i])/(pG->dx2);

        flx_m1l -= gyl*gxl/four_pi_G;
        flx_m1r -= gyr*gxr/four_pi_G;

        flx_m2l -= 0.5*(gyl*gyl-gxl*gxl)/four_pi_G + grav_mean_rho*phil_old;
        flx_m2r -= 0.5*(gyr*gyr-gxr*gxr)/four_pi_G + grav_mean_rho*phir_old;

/* Update momenta and energy with d/dx2 terms  */
        pG->U[ks][j][i].M1 -= 0.5*dtodx2*(flx_m1r - flx_m1l);
        pG->U[ks][j][i].M2 -= 0.5*dtodx2*(flx_m2r - flx_m2l);
#ifndef ISOTHERMAL
        pG->U[ks][j][i].E -=
           0.5*dtodx2*(pG->x2MassFlux[ks][j  ][i]*(dphic - dphil) +
                       pG->x2MassFlux[ks][j+1][i]*(dphir - dphic));
#endif
      }
    }

    break;

/*------------------------- 3D problem ---------------------------------------*/
  case 3:

/* Step 1 for 3D.  Add fluxes and source terms due to (d/dx1) terms  */

    for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        phic = pG->Phi[k][j][i];
        phil = 0.5*(pG->Phi[k][j][i-1] + pG->Phi[k][j][i  ]);
        phir = 0.5*(pG->Phi[k][j][i  ] + pG->Phi[k][j][i+1]);
        phil_old = 0.5*(pG->Phi_old[k][j][i-1] + pG->Phi_old[k][j][i  ]);
        phir_old = 0.5*(pG->Phi_old[k][j][i  ] + pG->Phi_old[k][j][i+1]);

        dphic = phic - pG->Phi_old[k][j][i];
        dphil = phil - phil_old;
        dphir = phir - phir_old;

/*  momentum fluxes in x1. gx, gy and gz centered at L and R x1-faces */
        gxl = (pG->Phi[k][j][i-1] - pG->Phi[k][j][i  ])/(pG->dx1);
        gxr = (pG->Phi[k][j][i  ] - pG->Phi[k][j][i+1])/(pG->dx1);

        gyl = (pG->Phi[k][j-1][i-1] - pG->Phi[k][j+1][i-1]) +
              (pG->Phi[k][j-1][i  ] - pG->Phi[k][j+1][i  ]);
        gyl *= (0.25/pG->dx2);

        gyr = (pG->Phi[k][j-1][i  ] - pG->Phi[k][j+1][i  ]) +
              (pG->Phi[k][j-1][i+1] - pG->Phi[k][j+1][i+1]);
        gyr *= (0.25/pG->dx2);

        gzl = (pG->Phi[k-1][j][i-1] - pG->Phi[k+1][j][i-1]) +
              (pG->Phi[k-1][j][i  ] - pG->Phi[k+1][j][i  ]);
        gzl *= (0.25/pG->dx3);

        gzr = (pG->Phi[k-1][j][i  ] - pG->Phi[k+1][j][i  ]) +
              (pG->Phi[k-1][j][i+1] - pG->Phi[k+1][j][i+1]);
        gzr *= (0.25/pG->dx3);

        flx_m1l = 0.5*(gxl*gxl-gyl*gyl-gzl*gzl)/four_pi_G + grav_mean_rho*phil;
        flx_m1r = 0.5*(gxr*gxr-gyr*gyr-gzr*gzr)/four_pi_G + grav_mean_rho*phir;

        flx_m2l = gxl*gyl/four_pi_G;
        flx_m2r = gxr*gyr/four_pi_G;

        flx_m3l = gxl*gzl/four_pi_G;
        flx_m3r = gxr*gzr/four_pi_G;

/*  subtract off momentum fluxes from old potential */
        gxl = (pG->Phi_old[k][j][i-1] - pG->Phi_old[k][j][i  ])/(pG->dx1);
        gxr = (pG->Phi_old[k][j][i  ] - pG->Phi_old[k][j][i+1])/(pG->dx1);

        gyl = (pG->Phi_old[k][j-1][i-1] - pG->Phi_old[k][j+1][i-1]) +
              (pG->Phi_old[k][j-1][i  ] - pG->Phi_old[k][j+1][i  ]);
        gyl *= (0.25/pG->dx2);

        gyr = (pG->Phi_old[k][j-1][i  ] - pG->Phi_old[k][j+1][i  ]) +
              (pG->Phi_old[k][j-1][i+1] - pG->Phi_old[k][j+1][i+1]);
        gyr *= (0.25/pG->dx2);

        gzl = (pG->Phi_old[k-1][j][i-1] - pG->Phi_old[k+1][j][i-1]) +
              (pG->Phi_old[k-1][j][i  ] - pG->Phi_old[k+1][j][i  ]);
        gzl *= (0.25/pG->dx3);

        gzr = (pG->Phi_old[k-1][j][i  ] - pG->Phi_old[k+1][j][i  ]) +
              (pG->Phi_old[k-1][j][i+1] - pG->Phi_old[k+1][j][i+1]);
        gzr *= (0.25/pG->dx3);

        flx_m1l -= 0.5*(gxl*gxl-gyl*gyl-gzl*gzl)/four_pi_G 
                   + grav_mean_rho*phil_old;
        flx_m1r -= 0.5*(gxr*gxr-gyr*gyr-gzr*gzr)/four_pi_G
                   + grav_mean_rho*phir_old;

        flx_m2l -= gxl*gyl/four_pi_G;
        flx_m2r -= gxr*gyr/four_pi_G;

        flx_m3l -= gxl*gzl/four_pi_G;
        flx_m3r -= gxr*gzr/four_pi_G;
/* Update momenta and energy with d/dx1 terms  */
        pG->U[k][j][i].M1 -= 0.5*dtodx1*(flx_m1r - flx_m1l);
        pG->U[k][j][i].M2 -= 0.5*dtodx1*(flx_m2r - flx_m2l);
        pG->U[k][j][i].M3 -= 0.5*dtodx1*(flx_m3r - flx_m3l);
#ifdef ADIABATIC
        pG->U[k][j][i].E -= 0.5*dtodx1*
          (pG->x1MassFlux[k][j][i  ]*(dphic - dphil) +
           pG->x1MassFlux[k][j][i+1]*(dphir - dphic));
#endif /* ADIABATIC */
      }
    }}

/* Step 2 for 3D.  Add fluxes and source terms due to (d/dx2) terms  */

    for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        phic = pG->Phi[k][j][i];
        phil = 0.5*(pG->Phi[k][j-1][i] + pG->Phi[k][j  ][i]);
        phir = 0.5*(pG->Phi[k][j  ][i] + pG->Phi[k][j+1][i]);
        phil_old = 0.5*(pG->Phi_old[k][j-1][i] + pG->Phi_old[k][j  ][i]);
        phir_old = 0.5*(pG->Phi_old[k][j  ][i] + pG->Phi_old[k][j+1][i]);

        dphic = phic - pG->Phi_old[k][j][i];
        dphil = phil - phil_old;
        dphir = phir - phir_old;

/* gx, gy and gz centered at L and R x2-faces */
        gxl = (pG->Phi[k][j-1][i-1] - pG->Phi[k][j-1][i+1]) +
              (pG->Phi[k][j  ][i-1] - pG->Phi[k][j  ][i+1]);
        gxl *= (0.25/pG->dx1);

        gxr = (pG->Phi[k][j  ][i-1] - pG->Phi[k][j  ][i+1]) +
              (pG->Phi[k][j+1][i-1] - pG->Phi[k][j+1][i+1]);
        gxr *= (0.25/pG->dx1);

        gyl = (pG->Phi[k][j-1][i] - pG->Phi[k][j  ][i])/(pG->dx2);
        gyr = (pG->Phi[k][j  ][i] - pG->Phi[k][j+1][i])/(pG->dx2);

        gzl = (pG->Phi[k-1][j-1][i] - pG->Phi[k+1][j-1][i]) +
              (pG->Phi[k-1][j  ][i] - pG->Phi[k+1][j  ][i]);
        gzl *= (0.25/pG->dx3);

        gzr = (pG->Phi[k-1][j  ][i] - pG->Phi[k+1][j  ][i]) +
              (pG->Phi[k-1][j+1][i] - pG->Phi[k+1][j+1][i]);
        gzr *= (0.25/pG->dx3);

        flx_m1l = gyl*gxl/four_pi_G;
        flx_m1r = gyr*gxr/four_pi_G;

        flx_m2l = 0.5*(gyl*gyl-gxl*gxl-gzl*gzl)/four_pi_G + grav_mean_rho*phil;
        flx_m2r = 0.5*(gyr*gyr-gxr*gxr-gzr*gzr)/four_pi_G + grav_mean_rho*phir;

        flx_m3l = gyl*gzl/four_pi_G;
        flx_m3r = gyr*gzr/four_pi_G;

/*  subtract off momentum fluxes from old potential */
        gxl = (pG->Phi_old[k][j-1][i-1] - pG->Phi_old[k][j-1][i+1]) +
              (pG->Phi_old[k][j  ][i-1] - pG->Phi_old[k][j  ][i+1]);
        gxl *= (0.25/pG->dx1);

        gxr = (pG->Phi_old[k][j  ][i-1] - pG->Phi_old[k][j  ][i+1]) +
              (pG->Phi_old[k][j+1][i-1] - pG->Phi_old[k][j+1][i+1]);
        gxr *= (0.25/pG->dx1);

        gyl = (pG->Phi_old[k][j-1][i] - pG->Phi_old[k][j  ][i])/(pG->dx2);
        gyr = (pG->Phi_old[k][j  ][i] - pG->Phi_old[k][j+1][i])/(pG->dx2);

        gzl = (pG->Phi_old[k-1][j-1][i] - pG->Phi_old[k+1][j-1][i]) +
              (pG->Phi_old[k-1][j  ][i] - pG->Phi_old[k+1][j  ][i]);
        gzl *= (0.25/pG->dx3);

        gzr = (pG->Phi_old[k-1][j  ][i] - pG->Phi_old[k+1][j  ][i]) +
              (pG->Phi_old[k-1][j+1][i] - pG->Phi_old[k+1][j+1][i]);
        gzr *= (0.25/pG->dx3);

        flx_m1l -= gyl*gxl/four_pi_G;
        flx_m1r -= gyr*gxr/four_pi_G;

        flx_m2l -= 0.5*(gyl*gyl-gxl*gxl-gzl*gzl)/four_pi_G
                 + grav_mean_rho*phil_old;
        flx_m2r -= 0.5*(gyr*gyr-gxr*gxr-gzr*gzr)/four_pi_G 
                 + grav_mean_rho*phir_old;

        flx_m3l -= gyl*gzl/four_pi_G;
        flx_m3r -= gyr*gzr/four_pi_G;

/* Update momenta and energy with d/dx2 terms  */
        pG->U[k][j][i].M1 -= 0.5*dtodx2*(flx_m1r - flx_m1l);
        pG->U[k][j][i].M2 -= 0.5*dtodx2*(flx_m2r - flx_m2l);
        pG->U[k][j][i].M3 -= 0.5*dtodx2*(flx_m3r - flx_m3l);
#ifdef ADIABATIC
        pG->U[k][j][i].E -= 0.5*dtodx2*
          (pG->x2MassFlux[k][j  ][i]*(dphic - dphil) +
           pG->x2MassFlux[k][j+1][i]*(dphir - dphic));
#endif /* ADIABATIC */
      }
    }}

/* Step 3 for 3D.  Add fluxes and source terms due to (d/dx3) terms  */

    for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        phic = pG->Phi[k][j][i];
        phil = 0.5*(pG->Phi[k-1][j][i] + pG->Phi[k  ][j][i]);
        phir = 0.5*(pG->Phi[k  ][j][i] + pG->Phi[k+1][j][i]);
        phil_old = 0.5*(pG->Phi_old[k-1][j][i] + pG->Phi_old[k  ][j][i]);
        phir_old = 0.5*(pG->Phi_old[k  ][j][i] + pG->Phi_old[k+1][j][i]);

        dphic = phic - pG->Phi_old[k][j][i];
        dphil = phil - phil_old;
        dphir = phir - phir_old;

/*  momentum fluxes in x3. gx, gy and gz centered at L and R x3-faces */
        gxl = (pG->Phi[k-1][j][i-1] - pG->Phi[k-1][j][i+1]) +
              (pG->Phi[k  ][j][i-1] - pG->Phi[k  ][j][i+1]);
        gxl *= (0.25/pG->dx1);

        gxr = (pG->Phi[k  ][j][i-1] - pG->Phi[k  ][j][i+1]) +
              (pG->Phi[k+1][j][i-1] - pG->Phi[k+1][j][i+1]);
        gxr *= (0.25/pG->dx1);

        gyl = (pG->Phi[k-1][j-1][i] - pG->Phi[k-1][j+1][i]) +
              (pG->Phi[k  ][j-1][i] - pG->Phi[k  ][j+1][i]);
        gyl *= (0.25/pG->dx2);

        gyr = (pG->Phi[k  ][j-1][i] - pG->Phi[k  ][j+1][i]) +
              (pG->Phi[k+1][j-1][i] - pG->Phi[k+1][j+1][i]);
        gyr *= (0.25/pG->dx2);

        gzl = (pG->Phi[k-1][j][i] - pG->Phi[k  ][j][i])/(pG->dx3);
        gzr = (pG->Phi[k  ][j][i] - pG->Phi[k+1][j][i])/(pG->dx3);

        flx_m1l = gzl*gxl/four_pi_G;
        flx_m1r = gzr*gxr/four_pi_G;

        flx_m2l = gzl*gyl/four_pi_G;
        flx_m2r = gzr*gyr/four_pi_G;

        flx_m3l = 0.5*(gzl*gzl-gxl*gxl-gyl*gyl)/four_pi_G + grav_mean_rho*phil;
        flx_m3r = 0.5*(gzr*gzr-gxr*gxr-gyr*gyr)/four_pi_G + grav_mean_rho*phir;

/*  subtract off momentum fluxes from old potential */
        gxl = (pG->Phi_old[k-1][j][i-1] - pG->Phi_old[k-1][j][i+1]) +
              (pG->Phi_old[k  ][j][i-1] - pG->Phi_old[k  ][j][i+1]);
        gxl *= (0.25/pG->dx1);

        gxr = (pG->Phi_old[k  ][j][i-1] - pG->Phi_old[k  ][j][i+1]) +
              (pG->Phi_old[k+1][j][i-1] - pG->Phi_old[k+1][j][i+1]);
        gxr *= (0.25/pG->dx1);

        gyl = (pG->Phi_old[k-1][j-1][i] - pG->Phi_old[k-1][j+1][i]) +
              (pG->Phi_old[k  ][j-1][i] - pG->Phi_old[k  ][j+1][i]);
        gyl *= (0.25/pG->dx2);

        gyr = (pG->Phi_old[k  ][j-1][i] - pG->Phi_old[k  ][j+1][i]) +
              (pG->Phi_old[k+1][j-1][i] - pG->Phi_old[k+1][j+1][i]);
        gyr *= (0.25/pG->dx2);

        gzl = (pG->Phi_old[k-1][j][i] - pG->Phi_old[k  ][j][i])/(pG->dx3);
        gzr = (pG->Phi_old[k  ][j][i] - pG->Phi_old[k+1][j][i])/(pG->dx3);

        flx_m1l -= gzl*gxl/four_pi_G;
        flx_m1r -= gzr*gxr/four_pi_G;

        flx_m2l -= gzl*gyl/four_pi_G;
        flx_m2r -= gzr*gyr/four_pi_G;

        flx_m3l -= 0.5*(gzl*gzl-gxl*gxl-gyl*gyl)/four_pi_G 
                 + grav_mean_rho*phil_old;
        flx_m3r -= 0.5*(gzr*gzr-gxr*gxr-gyr*gyr)/four_pi_G
                 + grav_mean_rho*phir_old;

/* Update momenta and energy with d/dx3 terms  */
        pG->U[k][j][i].M1 -= 0.5*dtodx3*(flx_m1r - flx_m1l);
        pG->U[k][j][i].M2 -= 0.5*dtodx3*(flx_m2r - flx_m2l);
        pG->U[k][j][i].M3 -= 0.5*dtodx3*(flx_m3r - flx_m3l);
#ifdef ADIABATIC
        pG->U[k][j][i].E -= 0.5*dtodx3*
          (pG->x3MassFlux[k  ][j][i]*(dphic - dphil) +
           pG->x3MassFlux[k+1][j][i]*(dphir - dphic));
#endif /* ADIABATIC */
      }
    }}

    break;

  } /* end of switch statement */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn VDFun_t selfg_init(MeshS *pM)
 *  \brief Initialize pointer to appropriate self-gravity f'n, 
 *   allocates memory for dPhi array used for flux correction.
 */

VDFun_t selfg_init(MeshS *pM)
{
  int dim = 0;

/* Calculate the dimensions  */
  if(pM->Nx[0] > 1) dim++;
  if(pM->Nx[1] > 1) dim++;
  if(pM->Nx[2] > 1) dim++;

/* test that user set values for constants */
  if (grav_mean_rho < 0.0)
    ath_error("[selfg_init] grav_mean_rho must be set >0 in prob generator\n");
  if (four_pi_G < 0.0)
    ath_error("[selfg_init] four_pi_G must be set >0 in prob generator\n");

/* Return function pointer based on dimensions and algorithm */

  switch(dim){
#ifdef SELF_GRAVITY_USING_MULTIGRID
  case 1:
    return selfg_multig_1d;
  case 2:
    return selfg_multig_2d;
  case 3:
    selfg_multig_3d_init(pM);
    return selfg_multig_3d;
#endif

/* for gravity using FFTs, also initialize plans and data for FFTW */
#ifdef SELF_GRAVITY_USING_FFT
  case 1:
    return selfg_fft_1d;
  case 2:
    selfg_fft_2d_init(pM);
    return selfg_fft_2d;
  case 3:
    selfg_fft_3d_init(pM);
    return selfg_fft_3d;
#endif
/* for gravity using FFTs with open BC, also initialize plans and data for FFTW */
#ifdef SELF_GRAVITY_USING_FFT_OBC
  case 1:
    ath_error("[selfg_init] FFT with open BC not defined for 1D \n");
  case 2:
    ath_error("[selfg_init] FFT with open BC not defined for 2D \n");
  case 3:
    selfg_fft_obc_3d_init(pM);
    return selfg_fft_obc_3d;
#endif



  }

  return NULL;
}
#endif /* SELF_GRAVITY */
