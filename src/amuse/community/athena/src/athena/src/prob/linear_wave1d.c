#include "copyright.h"
/*============================================================================*/
/*! \file linear_wave1d.c
 *  \brief Problem generator for plane-parallel, grid-aligned linear wave tests.
 *
 * PURPOSE: Problem generator for plane-parallel, grid-aligned linear wave
 *   tests.  If wavevector is in x2 (x3) direction, then the grid must be
 *   2D (3D). Only works for wavevector parallel to grid.  Tests in 2D and 3D
 *   with wavevector at an arbitrary angle are initialized with different
 *   functions (linear_wave2d/3d).
 *
 *   Can be used for either standing (problem/vflow=1.0) or travelling
 *   (problem/vflow=0.0) waves.
 *
 *   Configure --with-gravity=fft to check Jeans stability of plane waves
 *   propagating parallel to grid.
 *
 *   Configure --enable-resistivity and/or --with-viscosity=ns to check damping
 *   of linear waves by resistivity/ambipolar diffusion and/or viscosity
 *
 * USERWORK_AFTER_LOOP function computes L1 error norm in solution by comparing
 *   to initial conditions.  Problem must be evolved for an integer number of
 *   wave periods for this to work (only works for ideal MHD).		      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/* Initial solution on root level, shared with Userwork_after_loop  */
static ConsS ***RootSoln=NULL;
static int wave_flag;


/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  ConsS ***Soln;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,n,m,nx1,nx2,nx3,wave_dir;
  Real amp,vflow;
  Real d0,p0,u0,v0,w0,h0;
  Real x1,x2,x3,r,ev[NWAVE],rem[NWAVE][NWAVE],lem[NWAVE][NWAVE];
#ifdef MHD
  Real bx0,by0,bz0,xfact,yfact;
#endif /* MHD */
  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;

/* Read initial conditions  */

  wave_flag = par_geti("problem","wave_flag");
  amp = par_getd("problem","amp");
  vflow = par_getd("problem","vflow");
  wave_dir = par_getd("problem","wave_dir");

/* allocate memory for solution on this level.  If this is root level
 * also allocate memory for RootSoln */

  if ((Soln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS)))==NULL)
    ath_error("[problem]: Error allocating memory for Soln\n");
  if (pDomain->Level == 0){
    if ((RootSoln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS)))
      == NULL) ath_error("[problem]: Error alloc memory for RootSoln\n");
  }

/* Get eigenmatrix, where the quantities u0 and bx0 are parallel to the
 * wavevector, and v0,w0,by0,bz0 are perpendicular.  Background state chosen
 * carefully so wave speeds are well-separated: Va=1, Vf=2, Vs=0.5 */

  d0 = 1.0;
#ifndef ISOTHERMAL
  p0 = 1.0/Gamma;
  u0 = vflow*sqrt(Gamma*p0/d0);
#else
  u0 = vflow*Iso_csound;
#endif
  v0 = 0.0;
  w0 = 0.0;
#ifdef MHD
  bx0 = 1.0;
  by0 = sqrt(2.0);
  bz0 = 0.5;
  xfact = 0.0;
  yfact = 1.0;
#endif

  for (n=0; n<NWAVE; n++) {
    for (m=0; m<NWAVE; m++) {
      rem[n][m] = 0.0;
      lem[n][m] = 0.0;
    }
  }

#ifdef HYDRO
#ifdef ISOTHERMAL
  esys_roe_iso_hyd(u0,v0,w0,   ev,rem,lem);
#else
  h0 = ((p0/Gamma_1 + 0.5*d0*(u0*u0+v0*v0+w0*w0)) + p0)/d0;
  esys_roe_adb_hyd(u0,v0,w0,h0,ev,rem,lem);
  ath_pout(0,"Ux - Cs = %e, %e\n",ev[0],rem[0][wave_flag]);
  ath_pout(0,"Ux      = %e, %e\n",ev[1],rem[1][wave_flag]);
  ath_pout(0,"Ux + Cs = %e, %e\n",ev[4],rem[4][wave_flag]);
#endif /* ISOTHERMAL */
#endif /* HYDRO */

#ifdef MHD
#if defined(ISOTHERMAL)
  esys_roe_iso_mhd(d0,u0,v0,w0,bx0,by0,bz0,xfact,yfact,ev,rem,lem);
  ath_pout(0,"Ux - Cf = %e, %e\n",ev[0],rem[0][wave_flag]);
  ath_pout(0,"Ux - Ca = %e, %e\n",ev[1],rem[1][wave_flag]);
  ath_pout(0,"Ux - Cs = %e, %e\n",ev[2],rem[2][wave_flag]);
  ath_pout(0,"Ux + Cs = %e, %e\n",ev[3],rem[3][wave_flag]);
  ath_pout(0,"Ux + Ca = %e, %e\n",ev[4],rem[4][wave_flag]);
  ath_pout(0,"Ux + Cf = %e, %e\n",ev[5],rem[5][wave_flag]);
#else
  h0 = ((p0/Gamma_1+0.5*(bx0*bx0+by0*by0+bz0*bz0)+0.5*d0*(u0*u0+v0*v0+w0*w0))
               + (p0+0.5*(bx0*bx0+by0*by0+bz0*bz0)))/d0;
  esys_roe_adb_mhd(d0,u0,v0,w0,h0,bx0,by0,bz0,xfact,yfact,ev,rem,lem);
  ath_pout(0,"Ux - Cf = %e, %e\n",ev[0],rem[0][wave_flag]);
  ath_pout(0,"Ux - Ca = %e, %e\n",ev[1],rem[1][wave_flag]);
  ath_pout(0,"Ux - Cs = %e, %e\n",ev[2],rem[2][wave_flag]);
  ath_pout(0,"Ux      = %e, %e\n",ev[3],rem[3][wave_flag]);
  ath_pout(0,"Ux + Cs = %e, %e\n",ev[4],rem[4][wave_flag]);
  ath_pout(0,"Ux + Ca = %e, %e\n",ev[5],rem[5][wave_flag]);
  ath_pout(0,"Ux + Cf = %e, %e\n",ev[6],rem[6][wave_flag]);
#endif /* ISOTHERMAL */
#endif /* MHD */

/* Now initialize 1D solution vector  */

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie; i++) {

/* Set background state */

    cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
    Soln[k][j][i].d = d0;
#ifndef ISOTHERMAL
    Soln[k][j][i].E = p0/Gamma_1 + 0.5*d0*u0*u0;
#ifdef MHD
    Soln[k][j][i].E += 0.5*(bx0*bx0 + by0*by0 + bz0*bz0);
#endif /* MHD */
#endif /* ISOTHERMAL */

/* Select appropriate solution based on direction of wavevector */

    switch(wave_dir){

/* wave in x1-direction */

    case 1:
      Soln[k][j][i].d += amp*sin(2.0*PI*x1)*rem[0][wave_flag];
#ifndef ISOTHERMAL
      Soln[k][j][i].E += amp*sin(2.0*PI*x1)*rem[4][wave_flag];
#endif /* ISOTHERMAL */
      Soln[k][j][i].M1 = d0*vflow;
      Soln[k][j][i].M2 = 0.0;
      Soln[k][j][i].M3 = 0.0;
      Soln[k][j][i].M1 += amp*sin(2.0*PI*x1)*rem[1][wave_flag];
      Soln[k][j][i].M2 += amp*sin(2.0*PI*x1)*rem[2][wave_flag];
      Soln[k][j][i].M3 += amp*sin(2.0*PI*x1)*rem[3][wave_flag];
#ifdef MHD
      Soln[k][j][i].B1c = bx0;
      Soln[k][j][i].B2c = by0 + amp*sin(2.0*PI*x1)*rem[NWAVE-2][wave_flag];
      Soln[k][j][i].B3c = bz0 + amp*sin(2.0*PI*x1)*rem[NWAVE-1][wave_flag];
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        Soln[k][j][i].s[n] = amp*(1.0 + sin(2.0*PI*x1));
#endif
      break;

/* Wave in x2-direction */

    case 2:
      Soln[k][j][i].d += amp*sin(2.0*PI*x2)*rem[0][wave_flag];
#ifndef ISOTHERMAL
      Soln[k][j][i].E += amp*sin(2.0*PI*x2)*rem[4][wave_flag];
#endif /* ISOTHERMAL */
      Soln[k][j][i].M1 = 0.0;
      Soln[k][j][i].M2 = d0*vflow;
      Soln[k][j][i].M3 = 0.0;
      Soln[k][j][i].M1 += amp*sin(2.0*PI*x2)*rem[3][wave_flag];
      Soln[k][j][i].M2 += amp*sin(2.0*PI*x2)*rem[1][wave_flag];
      Soln[k][j][i].M3 += amp*sin(2.0*PI*x2)*rem[2][wave_flag];
#ifdef MHD
      Soln[k][j][i].B1c = bz0 + amp*sin(2.0*PI*x2)*rem[NWAVE-1][wave_flag];
      Soln[k][j][i].B2c = bx0;
      Soln[k][j][i].B3c = by0 + amp*sin(2.0*PI*x2)*rem[NWAVE-2][wave_flag];
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        Soln[k][j][i].s[n] = amp*(1.0 + sin(2.0*PI*x2));
#endif
      break;

/* Wave in x3-direction */

    case 3:
      Soln[k][j][i].d += amp*sin(2.0*PI*x3)*rem[0][wave_flag];
#ifndef ISOTHERMAL
      Soln[k][j][i].E += amp*sin(2.0*PI*x3)*rem[4][wave_flag];
#endif /* ISOTHERMAL */
      Soln[k][j][i].M1 = 0.0;
      Soln[k][j][i].M2 = 0.0;
      Soln[k][j][i].M3 = d0*vflow;
      Soln[k][j][i].M1 += amp*sin(2.0*PI*x3)*rem[2][wave_flag];
      Soln[k][j][i].M2 += amp*sin(2.0*PI*x3)*rem[3][wave_flag];
      Soln[k][j][i].M3 += amp*sin(2.0*PI*x3)*rem[1][wave_flag];
#ifdef MHD
      Soln[k][j][i].B1c = by0 + amp*sin(2.0*PI*x3)*rem[NWAVE-2][wave_flag];
      Soln[k][j][i].B2c = bz0 + amp*sin(2.0*PI*x3)*rem[NWAVE-1][wave_flag];
      Soln[k][j][i].B3c = bx0;
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        Soln[k][j][i].s[n] = amp*(1.0 + sin(2.0*PI*x3));
#endif
      break;
    default:
      ath_error("[linear_wave1d]: wave direction %d not allowed\n",wave_dir);
    }
  }}}

/* Now set initial conditions to wave solution */ 

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie; i++) {
    pGrid->U[k][j][i].d = Soln[k][j][i].d;
#ifndef ISOTHERMAL
    pGrid->U[k][j][i].E = Soln[k][j][i].E;
#endif /* ISOTHERMAL */
    pGrid->U[k][j][i].M1 = Soln[k][j][i].M1;
    pGrid->U[k][j][i].M2 = Soln[k][j][i].M2;
    pGrid->U[k][j][i].M3 = Soln[k][j][i].M3;
#ifdef MHD
    pGrid->U[k][j][i].B1c = Soln[k][j][i].B1c;
    pGrid->U[k][j][i].B2c = Soln[k][j][i].B2c;
    pGrid->U[k][j][i].B3c = Soln[k][j][i].B3c;
    pGrid->B1i[k][j][i] = Soln[k][j][i].B1c;
    pGrid->B2i[k][j][i] = Soln[k][j][i].B2c;
    pGrid->B3i[k][j][i] = Soln[k][j][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        pGrid->U[k][j][i].s[n] = Soln[k][j][i].s[n]; 
#endif
  }}}
#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      pGrid->B1i[k][j][ie+1] = pGrid->B1i[k][j][ie];
    }
  }
  if (pGrid->Nx[1] > 1) {
    for (k=ks; k<=ke; k++) {
      for (i=is; i<=ie; i++) {
        pGrid->B2i[k][je+1][i] = pGrid->B2i[k][je][i];
      }
    }
  }
  if (pGrid->Nx[2] > 1) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->B3i[ke+1][j][i] = pGrid->B3i[ke][j][i];
      }
    }
  }
#endif /* MHD */

/* For self-gravitating problems, read 4\piG and compute mean density */

#ifdef SELF_GRAVITY
  four_pi_G = par_getd("problem","four_pi_G");
  grav_mean_rho = d0;
#endif /* SELF_GRAVITY */

/* With viscosity and/or resistivity, read eta_R and nu_V */

#ifdef RESISTIVITY
  eta_Ohm = par_getd("problem","eta_O");
  Q_AD    = par_getd_def("problem","Q_AD",0.0);
  d_ind   = par_getd_def("problem","d_ind",0.0);
#endif
#ifdef NAVIER_STOKES
  nu_V = par_getd("problem","nu");
#endif

/* save solution on root grid */

  if (pDomain->Level == 0) {
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      RootSoln[k][j][i].d  = Soln[k][j][i].d ;
      RootSoln[k][j][i].M1 = Soln[k][j][i].M1;
      RootSoln[k][j][i].M2 = Soln[k][j][i].M2;
      RootSoln[k][j][i].M3 = Soln[k][j][i].M3;
#ifndef ISOTHERMAL
      RootSoln[k][j][i].E  = Soln[k][j][i].E ;
#endif /* ISOTHERMAL */
#ifdef MHD
      RootSoln[k][j][i].B1c = Soln[k][j][i].B1c;
      RootSoln[k][j][i].B2c = Soln[k][j][i].B2c;
      RootSoln[k][j][i].B3c = Soln[k][j][i].B3c;
#endif
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) 
        RootSoln[k][j][i].s[n] = Soln[k][j][i].s[n];
#endif
    }}}
  }

  free_3d_array(Soln);

  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
#ifdef SELF_GRAVITY
  Real d0 = 1.0;
  four_pi_G = par_getd("problem","four_pi_G");
  grav_mean_rho = d0;
#endif /* SELF_GRAVITY */
#ifdef RESISTIVITY
  eta_Ohm = par_getd_def("problem","eta_O",0.0);
  Q_Hall  = par_getd_def("problem","Q_H",0.0);
  Q_AD    = par_getd_def("problem","Q_AD",0.0);
#endif
#ifdef NAVIER_STOKES
  nu_V = par_getd("problem","nu");
#endif
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

#ifdef RESISTIVITY
/*! \fn void get_eta_user(GridS *pG, int i, int j, int k,
 *                           Real *eta_O, Real *eta_H, Real *eta_A)
 *  \brief Get user defined resistivity. */
void get_eta_user(GridS *pG, int i, int j, int k,
                             Real *eta_O, Real *eta_H, Real *eta_A)
{
  *eta_O = 0.0;
  *eta_H = 0.0;
  *eta_A = 0.0;

  return;
}
#endif

void Userwork_in_loop(MeshS *pM)
{
}

/*------------------------------------------------------------------------------
 * Userwork_after_loop: computes L1-error in linear waves,
 * ASSUMING WAVE HAS PROPAGATED AN INTEGER NUMBER OF PERIODS
 * Must set parameters (tlim, wavespeeds, box size) in input file appropriately
 * so that this is true
 * With nested grids, computes L1 error only on root grid (level=0).
 */

void Userwork_after_loop(MeshS *pM)
{
  GridS *pGrid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke;
#if (NSCALARS > 0)
   int n;
#endif

  Real rms_error=0.0;
  ConsS error,total_error;
  FILE *fp;
  char *fname;
  int Nx1, Nx2, Nx3, count;

  total_error.d = 0.0;
  total_error.M1 = 0.0;
  total_error.M2 = 0.0;
  total_error.M3 = 0.0;
#ifdef MHD
  total_error.B1c = 0.0;
  total_error.B2c = 0.0;
  total_error.B3c = 0.0;
#endif /* MHD */
#ifndef ISOTHERMAL
  total_error.E = 0.0;
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) total_error.s[n] = 0.0;
#endif

/* Compute error only on root Grid, which is in Domain[0][0] */

  pGrid=pM->Domain[0][0].Grid;
  if (pGrid == NULL) return;

/* compute L1 error in each variable, and rms total error */

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    error.d = 0.0;
    error.M1 = 0.0;
    error.M2 = 0.0;
    error.M3 = 0.0;
#ifdef MHD
    error.B1c = 0.0;
    error.B2c = 0.0;
    error.B3c = 0.0;
#endif /* MHD */
#ifndef ISOTHERMAL
    error.E = 0.0;
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) error.s[n] = 0.0;
#endif

    for (i=is; i<=ie; i++) {
      error.d   += fabs(pGrid->U[k][j][i].d   - RootSoln[k][j][i].d );
      error.M1  += fabs(pGrid->U[k][j][i].M1  - RootSoln[k][j][i].M1);
      error.M2  += fabs(pGrid->U[k][j][i].M2  - RootSoln[k][j][i].M2);
      error.M3  += fabs(pGrid->U[k][j][i].M3  - RootSoln[k][j][i].M3); 
#ifdef MHD
      error.B1c += fabs(pGrid->U[k][j][i].B1c - RootSoln[k][j][i].B1c);
      error.B2c += fabs(pGrid->U[k][j][i].B2c - RootSoln[k][j][i].B2c);
      error.B3c += fabs(pGrid->U[k][j][i].B3c - RootSoln[k][j][i].B3c);
#endif /* MHD */
#ifndef ISOTHERMAL
      error.E   += fabs(pGrid->U[k][j][i].E   - RootSoln[k][j][i].E );
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) 
        error.s[n] += fabs(pGrid->U[k][j][i].s[n] - RootSoln[k][j][i].s[n]);;
#endif
    }

    total_error.d += error.d;
    total_error.M1 += error.M1;
    total_error.M2 += error.M2;
    total_error.M3 += error.M3;
#ifdef MHD
    total_error.B1c += error.B1c;
    total_error.B2c += error.B2c;
    total_error.B3c += error.B3c;
#endif /* MHD */
#ifndef ISOTHERMAL
    total_error.E += error.E;
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) total_error.s[n] += error.s[n];
#endif
  }}

  Nx1 = ie - is + 1;
  Nx2 = je - js + 1;
  Nx3 = ke - ks + 1;

  count = Nx1*Nx2*Nx3;

/* Compute RMS error over all variables, and print out */

  rms_error = SQR(total_error.d) + SQR(total_error.M1) + SQR(total_error.M2)
                + SQR(total_error.M3);
#ifdef MHD
  rms_error += SQR(total_error.B1c) + SQR(total_error.B2c) 
               + SQR(total_error.B3c);
#endif /* MHD */
#ifndef ISOTHERMAL
  rms_error += SQR(total_error.E);
#endif /* ISOTHERMAL */
  rms_error = sqrt(rms_error)/(double)count;

/* Print error to file "LinWave-errors.#.dat", where #=wave_flag  */

  fname = ath_fname(NULL,"LinWave-errors",NULL,NULL,1,wave_flag,NULL,"dat");
/* The file exists -- reopen the file in append mode */
  if((fp=fopen(fname,"r")) != NULL){
    if((fp = freopen(fname,"a",fp)) == NULL){
      ath_error("[Userwork_after_loop]: Unable to reopen file.\n");
      return;
    }
  }
/* The file does not exist -- open the file in write mode */
  else{
    if((fp = fopen(fname,"w")) == NULL){
      ath_error("[Userwork_after_loop]: Unable to open file.\n");
      return;
    }
/* Now write out some header information */
    fprintf(fp,"# Nx1  Nx2  Nx3  RMS-Error  d  M1  M2  M3");
#ifndef ISOTHERMAL
    fprintf(fp,"  E");
#endif /* ISOTHERMAL */
#ifdef MHD
    fprintf(fp,"  B1c  B2c  B3c");
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) {
      fprintf(fp,"  S[ %d ]",n);
    }
#endif

    fprintf(fp,"\n#\n");
  }
 
  fprintf(fp,"%d  %d  %d  %e",Nx1,Nx2,Nx3,rms_error);

  fprintf(fp,"  %e  %e  %e  %e",
	  (total_error.d/(double) count),
	  (total_error.M1/(double)count),
	  (total_error.M2/(double)count),
	  (total_error.M3/(double)count));
#ifndef ISOTHERMAL
  fprintf(fp,"  %e",(total_error.E/(double)count));
#endif /* ISOTHERMAL */
#ifdef MHD
  fprintf(fp,"  %e  %e  %e",
	  (total_error.B1c/(double)count),
	  (total_error.B2c/(double)count),
	  (total_error.B3c/(double)count));
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) {
      fprintf(fp,"  %e",total_error.s[n]/(double)count);
    }
#endif
  fprintf(fp,"\n");

  fclose(fp);

  return;
}
