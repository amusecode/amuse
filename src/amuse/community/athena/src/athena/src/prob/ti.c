#include "copyright.h"
/*============================================================================*/
/*! \file ti.c
 *  \brief Problem generator for thermal instability.
 *
 *  PURPOSE: Problem generator for thermal instability.  Setup to use the 
 *   Koyama & Inutsuka cooling function, so must use cgs units for calculations
 *
 *  TI test w/ Isotropic Conduction - hydro
 *   The coefficient of conductivity kappa should be defined in cgs unit
 *
 *  TI test w/ Anisotropic Conduction - mhd
 *   The coefficient of Spitzer Coulombic conductivity should be given.
 *
 *  Various test cases are possible
 *- (iprob=1) : TI test w/ isotropic conduction in 1D - sinusoidal fluctuations
 *- (iprob=2) : TI test w/ isotropic conduction in 2D
 *                       with sinusoidal fluctuations along the diagonal        
 *- (iprob=3) : TI test w/ conduction in 2D with random pressure purtubation
 *              in case of anisotropic conduction : mhd - Diagonal field line
 *- (iprob=4) : TI test w/ conduction in 2D with random pressure purtubation
 *              in case of anisotropic conduction : mhd - Tangled field line  */
/*============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/* These constants are in cgs */
static const Real mbar = (1.27)*(1.6733e-24);
static const Real kb = 1.380658e-16;
static const Real pi = 3.14159265358979323846;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()
 *============================================================================*/

static double ran2(long int *idum);
static Real logd(const Grid *pG, const int i, const int j, const int k);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pGrid, Domain *pDomain)
{
  int i=0,j=0,k=0,n=0,m=0,l=0;
  int is,ie,js,je,ks,ke,iprob,
      nkx=16,nky=8,nkz=8,nx1,nx2,nx3;
  Real P_k, n0, T0, kappa,chi, x1L, x2L, x3L, x12L, num,
			 drho, dp,
       krho, krho_i, krho_ij, angle, phi, phi2, phi3, ksi, tmp,
			 rho0, rho1, n_anal, r, r_ij, r_proj, 
       beta, b0, Btot=0.0,Bsum=0.0,
       ***az,***ax,***ay,amp,kx[16],ky[8],kz[8],xa,ya,za;
  long int iseed = -1, rseed;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

/* Read problem parameters */
  P_k   = par_getd("problem","P_k");
  n0    = par_getd("problem","n0");
  x1L		= par_getd("grid","x1max");
	x2L   = par_getd("grid","x2max");
  x3L   = par_getd("grid","x3max");
	iprob = par_geti("problem","iprob");

#ifdef MHD
	beta  = par_getd("problem","beta");
//  b0    = par_getd("problem","b0");
	b0 = sqrt(2.0*P_k*kb/beta);
#endif
#ifdef ISOTROPIC_CONDUCTION
  kappa = par_getd("problem","kappa");
  kappa_T = (mbar/kb)*kappa;
#endif
#ifdef ANISOTROPIC_CONDUCTION
  chi = par_getd("problem","chi");
  chi_C = (mbar/kb)*chi;
#endif


/* Set up Initial Temperature and density */
	T0 = P_k / n0;
  rho0 = n0 * mbar;

/* For (iprob=1) -- TI w/ isotropic conduction - growth rate test in 1D
 *   2 with sinusoidal fluctuations in density, pressure and velocity  
 *    drho (initial fluctuation amplitude in density) should be given. 
 *    krho = 2pi/x1L * n */
if(iprob==1) {
  drho  = par_getd("problem","drho");
  n_anal = par_getd("problem","n_anal"); // analytic growth rate - Field 1965
  num   = par_getd("problem","num"); // wavenumber n

  krho_i = 2.0 * pi * num / (float)(ie-is+1);
  krho = 2.0 * pi * num / x1L;
	rho1 = rho0 * drho;
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pGrid->U[k][j][i].d = rho0 + rho1*cos(krho_i*(float)(i-is));
      pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d * n_anal * rho1 * sin(krho_i*(float)(i-is)) / rho0 / krho * (-1.0);
      pGrid->U[k][j][i].M2 = 0.0;
      pGrid->U[k][j][i].M3 = 0.0;
      pGrid->U[k][j][i].E = (n0*kb*T0)/Gamma_1 - n_anal*n_anal/ (krho * krho * Gamma_1) * rho1*cos(krho_i*(float)(i-is))
             + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
             + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
       }
     }
   } 
 }


/* For (iprob=2) -- TI w/ isotropic conduction - growth rate test in 2D
 *    with sinusoidal fluctuations in density, pressure and velocity
 *    drho (initial fluctuation amplitude in density) should be given.
 *    krho = 2pi/x12L * n  - x12L (diagonal length of the 2D box)      */
if(iprob==2) {
	angle = atan(x2L/x1L);
	x12L = pow(x1L*x1L + x2L*x2L, 0.5);	
	r_ij = pow(pow((float)(ie-is),2)+pow((float)(je-js),2),0.5);
	krho = 2.0 * pi * num / x12L;
	krho_ij = 2.0 * pi * num / r_ij;
	rho1 = rho0 * drho;
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      r = pow(pow((float)(i-is),2)+pow((float)(j-js),2),0.5);
      if ((i-is) == 0) {
      phi = 0;
      }
      else{
      phi = atan((float)(j-js)/(float)(i-is));
      }
      r_proj = r * sin(phi + angle) / sin ( 2 * angle);
      pGrid->U[k][j][i].d = rho0 + rho1 * cos(krho_ij * r_proj);
      pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d * n_anal * rho1 * sin(krho_ij * r_proj) / rho0 / krho * cos(angle) * (-1.0);
      pGrid->U[k][j][i].M2 = pGrid->U[k][j][i].d * n_anal * rho1 * sin(krho_ij * r_proj) / rho0 / krho * sin(angle) * (-1.0);
      pGrid->U[k][j][i].M3 = 0.0;
      pGrid->U[k][j][i].E = (n0*kb*T0)/Gamma_1 
             - n_anal*n_anal/ (krho * krho * Gamma_1) * rho1 * cos(krho_ij * r_proj)
             + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
             + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
       }
     }
   }
 }

/* For (iprob=3) -- TI w/ (an)isotropic conduction - growth rate test in 2D
 *    with RANDOM fluctuations in pressure
 *    dp (initial fluctuation amplitude in pressure) should be given.
 *    iprob=3 - if MHD - diagonal streight field line / 
                if HYDRO - isotropic conduction test with random perturbation in pressure  */
 if(iprob==3) {
  dp  = par_getd("problem","dp");
  angle = atan(x2L/x1L);
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pGrid->U[k][j][i].d = rho0; 
      pGrid->U[k][j][i].M1 = 0.0; 
      pGrid->U[k][j][i].M2 = 0.0;  
      pGrid->U[k][j][i].M3 = 0.0;
      pGrid->U[k][j][i].E = (n0*kb*T0)/Gamma_1 + (n0*kb*T0)/Gamma_1 * (ran2(&rseed) - 0.5) * dp
             + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
             + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#ifdef MHD
      pGrid->B1i[k][j][i] = b0 * cos(angle);
      if (i == ie) pGrid->B1i[k][j][i+1] = b0 * cos(angle);
      pGrid->U[k][j][i].B1c = b0 * cos(angle);

      pGrid->B2i[k][j][i] = b0 * sin(angle);
      if (j==je) pGrid->B2i[k][j+1][i] = b0 * sin(angle);
      pGrid->U[k][j][i].B2c = b0 * sin(angle);

      pGrid->B3i[k][j][i] = 0.0;
      if (k==ke && pGrid->Nx3 > 1) pGrid->B3i[k+1][j][i] = 0.0;
      pGrid->U[k][j][i].B3c = 0.0;
      pGrid->U[k][j][i].E += 0.5*(SQR(pGrid->U[k][j][i].B1c)
                           + SQR(pGrid->U[k][j][i].B2c) + SQR(pGrid->U[k][j][i].B3c));
#endif /* MHD */
      }
    }
   }
}

/* For (iprob=4) -- TI w/ (an)isotropic conduction - growth rate test in 2D
 *    with RANDOM fluctuations in pressure
 *    dp (initial fluctuation amplitude in pressure) should be given.
 *    iprob=4 - tangled field line */
if(iprob==4){

  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  if ((nx1 == 1) || (nx2 == 1)) {
    ath_error("[TI-tangledB]: This problem can only be run with Nx1>1\n");
  }
  if ((az = (Real***)calloc_3d_array(nx3,nx2, nx1, sizeof(Real))) == NULL) {
    ath_error("[TI-tangledB]: Error allocating memory for vector pot\n");
  }

/* Vector Potential Generation*/
  for (n=0; n<=nkx-1; n++){
    kx[n] = (1.+(float)(n))*(2.0*pi)/2.0;
  }
	kx[nkx-1] = kx[0] * (-1.0);
  for (m=0; m<=nky-1; m++){
    ky[m] = (1.+(float)(m))*(2.0*pi);
  }
  ky[nky-1] = ky[0] * (-1.0);

  for (k=ks; k=ks; k++) {
  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      az[k][j][i] = 0.0;
     }
   }
  }
  for (n=0; n<=nkx-1; n++){
    for (m=0; m<=nky-1; m++){
      phi = 2.0 * pi * ran2(&rseed);
//			printf("phi=%f",phi);
      amp = (1.0/sqrt( kx[n]*kx[n] + ky[m]*ky[m]));
     for (j=js; j<=je+1; j++) {
       ya=(float)(j-js)*(1.0/(float)(je-js+1));
/*			 if(j == js | j==je+1){
				printf("ya=%f\n",ya);
       }*/
      for (i=is; i<=ie+1; i++) {
       xa=(float)(i-is)*(2.0/(float)(ie-is+1));
        az[ks][j][i] += amp*sin(kx[n]*xa + ky[m]*ya + phi);
       }
     }
    }
  }
/*Vector Potential for 2D generation end*/

#ifdef MHD
  dp  = par_getd("problem","dp");
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pGrid->U[k][j][i].d = rho0;
      pGrid->U[k][j][i].M1 = 0.0;
      pGrid->U[k][j][i].M2 = 0.0;
      pGrid->U[k][j][i].M3 = 0.0;
      pGrid->U[k][j][i].E = (n0*kb*T0)/Gamma_1 + (n0*kb*T0)/Gamma_1 * (ran2(&rseed) - 0.5) * dp
             + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
             + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
      pGrid->B1i[k][j][i] = (az[ks][j+1][i] - az[ks][j][i]);
      pGrid->B2i[k][j][i] =-(az[ks][j][i+1] - az[ks][j][i]);

      pGrid->B3i[k][j][i] = 0.0;
      if (k==ke && pGrid->Nx3 > 1) pGrid->B3i[k+1][j][i] = 0.0;

      Bsum += (pGrid->B1i[k][j][i]*pGrid->B1i[k][j][i] + 
             pGrid->B2i[k][j][i]*pGrid->B2i[k][j][i] + pGrid->B3i[k][j][i]*pGrid->B3i[k][j][i]);
       }
     }
   }
	Btot = sqrt(Bsum/((float)(je-js+1) * (float)(ie-is+1)));
  printf("Btot=%f\n",Btot);
	Bsum = 0.0;

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pGrid->B1i[k][j][i] = pGrid->B1i[k][j][i] * b0 / Btot;
      if (i == ie) pGrid->B1i[k][j][i+1] = pGrid->B1i[k][j][is];
      pGrid->B2i[k][j][i] = pGrid->B2i[k][j][i] * b0 / Btot;
      if (j == je) pGrid->B2i[k][j+1][i] = pGrid->B2i[k][js][i];
      Bsum += (pGrid->B1i[k][j][i]*pGrid->B1i[k][j][i] +
             pGrid->B2i[k][j][i]*pGrid->B2i[k][j][i] + pGrid->B3i[k][j][i]*pGrid->B3i[k][j][i]);
       }
     }
   }

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i]+pGrid->B1i[k][j][i+1]);
      pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i]+pGrid->B2i[k][j+1][i]);
      pGrid->U[k][j][i].B3c = 0.0;
      pGrid->U[k][j][i].E += 0.5*(SQR(pGrid->U[k][j][i].B1c)
                           + SQR(pGrid->U[k][j][i].B2c) + SQR(pGrid->U[k][j][i].B3c));
       }
     }
   }
  Btot = Bsum/((float)(je-js+1) * (float)(ie-is+1));
  printf("B^2tot,b0^2=%e  %e \n",Btot,b0*b0);
#endif //MHD
 }

/* For (iprob=5) -- TI w/ (an)isotropic conduction - growth rate test in 3D
 *    with RANDOM fluctuations in pressure
 *    dp (initial fluctuation amplitude in pressure) should be given.
 *    iprob=5 - if MHD - diagonal streight field line /
                if HYDRO - isotropic conduction test with random perturbation in pressure  */
 if(iprob==5) {
  dp  = par_getd("problem","dp");
  angle = atan(x2L/x1L);
  x12L = pow(x1L*x1L + x2L*x2L, 0.5);
	ksi = atan(x3L/x12L);
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pGrid->U[k][j][i].d = rho0;
      pGrid->U[k][j][i].M1 = 0.0;
      pGrid->U[k][j][i].M2 = 0.0;
      pGrid->U[k][j][i].M3 = 0.0;
      pGrid->U[k][j][i].E = (n0*kb*T0)/Gamma_1 + (n0*kb*T0)/Gamma_1 * (ran2(&rseed) - 0.5) * dp
             + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
             + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#ifdef MHD
      pGrid->B1i[k][j][i] = b0 * cos(angle) * cos(ksi);
      if (i == ie) pGrid->B1i[k][j][i+1] = b0 * cos(angle) * cos(ksi);
      pGrid->U[k][j][i].B1c = b0 * cos(angle) * cos(ksi);

      pGrid->B2i[k][j][i] = b0 * sin(angle) * cos(ksi);
      if (j==je) pGrid->B2i[k][j+1][i] = b0 * sin(angle) * cos(ksi);
      pGrid->U[k][j][i].B2c = b0 * sin(angle) * cos(ksi);

      pGrid->B3i[k][j][i] = b0 * sin(ksi);
      if (k==ke && pGrid->Nx3 > 1) pGrid->B3i[k+1][j][i] = b0 * sin(ksi);
      pGrid->U[k][j][i].B3c = b0 * sin(ksi);
      pGrid->U[k][j][i].E += 0.5*(SQR(pGrid->U[k][j][i].B1c)
                           + SQR(pGrid->U[k][j][i].B2c) + SQR(pGrid->U[k][j][i].B3c));
#endif /* MHD */
      }
    }
   }
}


/* For (iprob=6) -- TI w/ (an)isotropic conduction - growth rate test in *3D*
 *    with RANDOM fluctuations in pressure
 *    dp (initial fluctuation amplitude in pressure) should be given.
 *    iprob=6 - in mhd, tangled field line */
if(iprob==6){
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;
  if ((nx1 == 1) || (nx2 == 1) || (nx3 == 1)) {
    ath_error("[TI-tangledB]: This problem can only be run with Nx1>1\n");
  }
  if ((ay = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL) {
    ath_error("[TI-tangledB]: Error allocating memory for vector pot\n");
  }
  if ((az = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL) {
    ath_error("[TI-tangledB]: Error allocating memory for vector pot\n");
  }
  if ((ax = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL) {
    ath_error("[TI-tangledB]: Error allocating memory for vector pot\n");
  }

/* Vector Potential Generation*/
  for (n=0; n<=nkx-1; n++){
    kx[n] = (1.+(float)(n))*(2.0*pi)/2.0;
  }
  kx[nkx-1] = kx[0] * (-1.0);
  for (m=0; m<=nky-1; m++){
    ky[m] = (1.+(float)(m))*(2.0*pi);
  }
  ky[nky-1] = ky[0] * (-1.0);
  for (l=0; l<=nkz-1; l++){
    kz[l] = (1.+(float)(l))*(2.0*pi);
  }
  kz[nkz-1] = kz[0] * (-1.0);

  for (k=ks; k<=ke+1; k++) {
   for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      ax[k][j][i] = 0.0;
      ay[k][j][i] = 0.0;
      az[k][j][i] = 0.0;
    }
   }
  }

  for (n=0; n<=nkx-1; n++){
   for (m=0; m<=nky-1; m++){
    for (l=0; l<=nkz-1; l++){
      phi = 2.0 * pi * ran2(&rseed);
      phi2 = 2.0 * pi * ran2(&rseed);
      phi3 = 2.0 * pi * ran2(&rseed);
      amp = (1.0/sqrt(kx[n]*kx[n] + ky[m]*ky[m] + kz[l]*kz[l]));
     for (k=ks; k<=ke+1; k++) {
       za=(float)(k-ks)*(1.0/(float)(ke-ks+1));
     for (j=js; j<=je+1; j++) {
       ya=(float)(j-js)*(1.0/(float)(je-js+1));
      for (i=is; i<=ie+1; i++) {
       xa=(float)(i-is)*(2.0/(float)(ie-is+1));
       tmp = kx[n]*xa + ky[m]*ya + kz[l]*za;
        az[k][j][i] += amp*sin(tmp + phi);
        ay[k][j][i] += amp*sin(tmp + phi2);
        ax[k][j][i] += amp*sin(tmp + phi3);
       }
      }
     }
    }
   }
  }

  printf("za=%f\n",za);
  Bsum = 0.0;
#ifdef MHD
  dp  = par_getd("problem","dp");
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pGrid->U[k][j][i].d = rho0;
      pGrid->U[k][j][i].M1 = 0.0;
      pGrid->U[k][j][i].M2 = 0.0;
      pGrid->U[k][j][i].M3 = 0.0;
      pGrid->U[k][j][i].E = (n0*kb*T0)/Gamma_1 + (n0*kb*T0)/Gamma_1 * (ran2(&rseed) - 0.5) * dp
             + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
             + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
      pGrid->B1i[k][j][i] = (az[k][j+1][i] - az[k][j][i]) -
                            (ay[k+1][j][i] - ay[k][j][i]);
      pGrid->B2i[k][j][i] = (ax[k+1][j][i] - ax[k][j][i]) -
                            (az[k][j][i+1] - az[k][j][i]);
      pGrid->B3i[k][j][i] = (ay[k][j][i+1] - ay[k][j][i]) -
                            (ax[k][j+1][i] - ax[k][j][i]);

      Bsum += sqrt(pGrid->B1i[k][j][i]*pGrid->B1i[k][j][i] +
             pGrid->B2i[k][j][i]*pGrid->B2i[k][j][i] + pGrid->B3i[k][j][i]*pGrid->B3i[k][j][i]);
       }
     }
   }
  Btot = Bsum/((float)(je-js+1) * (float)(ie-is+1) * (float)(ke-ks+1));
  printf("Btot=%f\n",Btot);
  Bsum = 0.0;
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pGrid->B1i[k][j][i] = pGrid->B1i[k][j][i] * b0 / Btot;
      if (i == ie) pGrid->B1i[k][j][i+1] = pGrid->B1i[k][j][is];
      pGrid->B2i[k][j][i] = pGrid->B2i[k][j][i] * b0 / Btot;
      if (j == je) pGrid->B2i[k][j+1][i] = pGrid->B2i[k][js][i];
      pGrid->B3i[k][j][i] = pGrid->B3i[k][j][i] * b0 / Btot;
      if (k == ke) pGrid->B3i[k+1][j][i] = pGrid->B3i[ks][j][i];

      Bsum += sqrt(pGrid->B1i[k][j][i]*pGrid->B1i[k][j][i] +
             pGrid->B2i[k][j][i]*pGrid->B2i[k][j][i] + pGrid->B3i[k][j][i]*pGrid->B3i[k][j][i]);
       }
     }
   }

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i]+pGrid->B1i[k][j][i+1]);
      pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i]+pGrid->B2i[k][j+1][i]);
      pGrid->U[k][j][i].B3c = 0.5*(pGrid->B3i[k][j][i]+pGrid->B3i[k+1][j][i]);
      pGrid->U[k][j][i].E += 0.5*(SQR(pGrid->U[k][j][i].B1c)
                           + SQR(pGrid->U[k][j][i].B2c) + SQR(pGrid->U[k][j][i].B3c));
       }
     }
   }
  Btot = Bsum/((float)(je-js+1) * (float)(ie-is+1) * (float)(ke-ks+1));
  printf("Btot/b0=%e  %e\n",Btot,b0);
#endif //MHD
 }
//printf("Btot=%f\n",xa);
  CoolingFunc = KoyInut;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(Grid *pG, Domain *pD, FILE *fp)
{
  return;
}

void problem_read_restart(Grid *pG, Domain *pD, FILE *fp)
{
  return;
}

/*! \fn static Real logd(const Grid *pG, const int i, const int j, const int k)
 *  \brief Log10 of density */
static Real logd(const Grid *pG, const int i, const int j, const int k)
{
  return log10(pG->U[k][j][i].d);
}

Gasfun_t get_usr_expr(const char *expr)
{
  if(strcmp(expr,"logd")==0) return logd;
  return NULL;
}

VGFunout_t get_usr_out_fun(const char *name){
  return NULL;
}


void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
  return;
}

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/*------------------------------------------------------------------------------
 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/*! \fn double ran2(long int *idum){
 *  \brief Extracted from the Numerical Recipes in C (version 2) code.  Modified
 *   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003
 *
 * Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1. 

 */
double ran2(long int *idum){
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX
