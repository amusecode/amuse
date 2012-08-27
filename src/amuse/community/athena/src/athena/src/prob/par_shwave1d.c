#include "copyright.h"
/*============================================================================*/
/*! \file par_shwave1d.c
 *  \brief Problem generator radial epicyclic wave test with Lagrangian
 *   particles.
 *
 * PURPOSE: Problem generator radial epicyclic wave test with Lagrangian
 *   particles.
 *
 *   Only works for 2D, and the wavevector is in x1 direction. Particles can be
 *   initialized either with uniform density or with density profile same as the
 *   gas. SHEARING_BOX must be turned on, FARGO is optional.
 *   ipert: 1 for linear wave; 2 for non-linear wave from Fromang & Papaloizou.
 *
 * Reference: Fromang & Papaloizou, A&A, 468, 1-18 (2007).
 */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "particles/particle.h"

#ifndef SHEARING_BOX
#error : The shear wave problem requires shearing-box to be enabled.
#endif /* SHEARING_BOX */


/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * GetPosition()    - get particle status (grid/crossing)
 * ShearingBoxPot() - shearing box tidal potential
 *============================================================================*/
#ifdef PARTICLES
int GetPosition(Grain *gr);
#endif
static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3);

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(Grid *pGrid, Domain *pDomain)
{
  FILE *fp;
  Real xFP[160],dFP[160],vxFP[160],vyFP[160];
  int i=0,j=0,k=0,ipert;
  int is,ie,js,je,ks,ke,n,m,nwave,samp;
  Real x1,x2,x3,x1max,x1min,x2max,x2min;
  Real kx,omg,omg2,amp,v1x,v1y;
#ifdef PARTICLES
  long p;
  int Npar,ip,jp;
  Real x1p,x2p,x3p,x1l,x1u,x2l,x2u,par_amp,factor2;
#endif

  if ((par_geti("grid","Nx2") == 1) || (par_geti("grid","Nx3") > 1)) {
    ath_error("[par_shwave1d]: par_shwave1d must work in 2D grid.\n");
  }

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

/* Read initial conditions  */
  amp = par_getd("problem","amp");
  Omega_0 = par_getd_def("problem","omega",1.0);
  qshear  = par_getd_def("problem","qshear",1.5);
  ipert = par_geti_def("problem","ipert", 1);
  nwave = par_geti("problem","nwave");
  samp = par_geti("problem","sample");
  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");
  if (ipert != 1) samp = 2;

#ifdef PARTICLES
/* Read initial conditions for the particles */
  /* basic parameters */
  if (par_geti("particle","partypes") != 1)
    ath_error("[par_shwave1d]: This test only allows ONE particle species!\n");

  Npar = (int)(sqrt(par_geti("particle","parnumcell")));
  pGrid->nparticle = Npar*Npar*pGrid->Nx1*pGrid->Nx2;
  pGrid->grproperty[0].num = pGrid->nparticle;
  if (pGrid->nparticle+2 > pGrid->arrsize)
    particle_realloc(pGrid, pGrid->nparticle+2);

  /* particle stopping time */
  tstop0[0] = par_getd_def("problem","tstop",0.0); /* in code unit */
  if (par_geti("particle","tsmode") != 3)
    ath_error("[par_shwave1d]: This test only allows fixed stopping time!\n");

  par_amp = amp;
  factor2 = 0.5;
#endif /* PARTICLES */

/* Initialize gas and particles */
  if (ipert == 1)
  {
    if (nwave <= 0)
      ath_error("[par_shwave1d]: nwave must be positive!\n");

  /* calculate dispersion relation and find eigen vectors */
    kx = 2.0*(PI)*nwave/(x1max-x1min);
    omg2 = SQR(Iso_csound*kx)+2.0*(2.0-qshear)*SQR(Omega_0);
    omg = sqrt(omg2);
    v1x = omg*amp/kx;
    v1y = (2.0-qshear)*Omega_0/omg*v1x;
#ifdef PARTICLES
    /* particle perturbation amplitude */
//    par_amp = amp*kx*pGrid->dx1/sin(kx*pGrid->dx1);
//    factor2 = 0.5*tan(kx*pGrid->dx1)/(kx*pGrid->dx1);
#endif

  /* Now set initial conditions to wave solution */ 
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
      pGrid->U[k][j][i].d = 1.0+amp*cos(kx*x1);
      pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d*v1x*cos(kx*x1);
      pGrid->U[k][j][i].M3 = pGrid->U[k][j][i].d*v1y*sin(kx*x1);
#ifndef FARGO
      pGrid->U[k][j][i].M3 -= pGrid->U[k][j][i].d*qshear*Omega_0*x1;
#endif
      pGrid->U[k][j][i].M2 = 0.0;
#if (NSCALARS > 0)
      if (samp == 1)
        for (n=0; n<NSCALARS; n++)
          pGrid->U[k][j][i].s[n] = pGrid->U[k][j][i].d;
      else
        for (n=0; n<NSCALARS; n++)
          pGrid->U[k][j][i].s[n] = 1.0;
#endif
    }}}
  }
  else /* ipert != 1 (read FP solution */
  {
    if (pGrid->Nx1 == 160) {
      if((fp = fopen("Data-160-FPwave.dat","r")) == NULL)
         ath_error("Error opening Data-160-FPwave.dat\n");
      for (i=0; i<160; i++) {
        fscanf(fp,"%lf %lf %lf %lf",&xFP[i],&dFP[i],&vxFP[i],&vyFP[i]);
      }
    }

    if (pGrid->Nx1 == 40) {
      if((fp = fopen("Data-40-FPwave.dat","r")) == NULL)
         ath_error("Error opening Data-40-FPwave.dat\n");
      for (i=0; i<40; i++) {
        fscanf(fp,"%lf %lf %lf %lf",&xFP[i],&dFP[i],&vxFP[i],&vyFP[i]);
      }
    }

    if (x1min != -4.7965)
      ath_error("[par_shwave1d]: ipert=2 requires x1min=-4.7965\n");
    if (x1max != 4.7965)
      ath_error("[par_shwave1d]: ipert=2 requires x1max=4.7965\n");

  /* Now set initial conditions to wave solution */

    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
      pGrid->U[k][j][i].d = dFP[i+pGrid->idisp];
      pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d*vxFP[i+pGrid->idisp];
      pGrid->U[k][j][i].M3 = pGrid->U[k][j][i].d*vyFP[i+pGrid->idisp];
#ifdef FARGO
      pGrid->U[k][j][i].M3 += pGrid->U[k][j][i].d*qshear*Omega_0*x1;
#endif
      pGrid->U[k][j][i].M2 = 0.0;
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        pGrid->U[k][j][i].s[n] = 1.0;
#endif
    }}}
  }

#ifdef PARTICLES
  /* Now set initial conditions for the particles */
  p = 0;
  x3p = pGrid->x3_0 + (pGrid->ks+pGrid->kdisp)*pGrid->dx3;

  for (j=pGrid->js; j<=pGrid->je; j++)
  {
    x2l = pGrid->x2_0 + (j+pGrid->jdisp)*pGrid->dx2;
    x2u = pGrid->x2_0 + ((j+pGrid->jdisp)+1.0)*pGrid->dx2;

    for (i=pGrid->is; i<=pGrid->ie; i++)
    {
      x1l = pGrid->x1_0 + (i + pGrid->idisp)*pGrid->dx1;
      x1u = pGrid->x1_0 + ((i + pGrid->idisp) + 1.0)*pGrid->dx1;

      for (ip=0;ip<Npar;ip++)
      {
        x1p = x1l+(x1u-x1l)/Npar*(ip+0.5);

        for (jp=0;jp<Npar;jp++)
        {
          x2p = x2l+(x2u-x2l)/Npar*(jp+0.5);

          pGrid->particle[p].property = 0;

          pGrid->particle[p].x1 = x1p;
          pGrid->particle[p].x2 = x2p;
          if (samp == 1)
            pGrid->particle[p].x1 += -par_amp*sin(kx*x1p)/kx
                                     +factor2*SQR(amp)*sin(2.0*kx*x1p)/kx;
          pGrid->particle[p].x3 = x3p;

          pGrid->particle[p].v1 = v1x*cos(kx*x1p);
          pGrid->particle[p].v3 = v1y*sin(kx*x1p);
#ifndef FARGO
          pGrid->particle[p].v3 -= qshear*Omega_0*x1;
#endif
          pGrid->particle[p].v2 = 0.0;

          pGrid->particle[p].pos = GetPosition(&pGrid->particle[p]);

          pGrid->particle[p].my_id = p;
#ifdef MPI_PARALLEL
          pGrid->particle[p].init_id = pGrid->my_id;
#endif
          p += 1;
        }
      }
    }
  }
#endif /* PARTICLES */

/* enroll gravitational potential function, shearing sheet BC functions */

  StaticGravPot = ShearingBoxPot;

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

void problem_write_restart(Grid *pG, Domain *pD, FILE *fp)
{
  return;
}

void problem_read_restart(Grid *pG, Domain *pD, FILE *fp)
{
  Omega_0 = par_getd_def("problem","omega",1.0);
  qshear  = par_getd_def("problem","qshear",1.5);

  StaticGravPot = ShearingBoxPot;
  return;
}

#if (NSCALARS > 0)
/*! \fn static Real ScalarDen(const Grid *pG,const int i,
 *			      const int j,const int k)
 *  \brief Scalar density */
static Real ScalarDen(const Grid *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].s[0];
}
#endif

#ifdef PARTICLES
/*! \fn static Real dratio(const Grid *pG, const int i, const int j,const int k)
 *  \brief Density ratio */
static Real dratio(const Grid *pG, const int i, const int j, const int k)
{
#if (NSCALARS > 0)
  return pG->Coup[k][j][i].grid_d/pG->U[k][j][i].s[0];
#else
  return pG->Coup[k][j][i].grid_d/pG->U[k][j][i].d;
#endif
}
#endif

/*! \fn static Real expr_dV3(const Grid *pG,const int i,const int j,const int k)
 *  \brief 3-component of velocity */
static Real expr_dV3(const Grid *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  return pG->U[k][j][i].M3/pG->U[k][j][i].d;
#else
  return (pG->U[k][j][i].M3/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
}

#ifdef PARTICLES
/*! \fn static Real expr_dV3par(const Grid *pG, const int i, const int j, 
 *				const int k)
 *  \brief 3-component of particle velocity */
static Real expr_dV3par(const Grid *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  return expr_V3par(pG, i, j, k);
#else
  return expr_V3par(pG, i, j, k) + qshear*Omega_0*x1;
#endif
}
#endif /* PARTICLES */

Gasfun_t get_usr_expr(const char *expr)
{
#if (NSCALARS > 0)
  if(strcmp(expr,"scalar")==0) return ScalarDen;
#endif
#ifdef PARTICLES
  if(strcmp(expr,"dratio")==0) return dratio;
  if(strcmp(expr,"dVypar")==0) return expr_dV3par;
#endif
  if(strcmp(expr,"dVy")==0) return expr_dV3;
  return NULL;
}

VGFunout_t get_usr_out_fun(const char *name){
  return NULL;
}

#ifdef PARTICLES
PropFun_t get_usr_par_prop(const char *name)
{
  return NULL;
}

void gasvshift(const Real x1, const Real x2, const Real x3,
                                    Real *u1, Real *u2, Real *u3)
{
  return;
}

void Userforce_particle(Vector *ft, const Real x1, const Real x2, const Real x3,
                                    const Real v1, const Real v2, const Real v3)
{
  return;
}
#endif

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
  return;
}

/*---------------------------------------------------------------------------
 * Userwork_after_loop: computes L1-error in linear waves,
 * ASSUMING WAVE HAS PROPAGATED AN INTEGER NUMBER OF PERIODS
 * Must set parameters in input file appropriately so that this is true
 */

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,n,m,nwave,samp;
  Real x1,x2,x3,x1max,x1min,x2max,x2min;
  Real time,kx,omg,omg2,v1x,v1y,amp,SolGasd,SolLagd,SolGasM1,SolGasM3;
  char *fname;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  /* Read initial conditions  */
  amp = par_getd("problem","amp");
  Omega_0 = par_getd_def("problem","omega",1.0);
  qshear  = par_getd_def("problem","qshear",1.5);
  nwave = par_geti("problem","nwave");
  samp = par_geti("problem","sample");
  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");

/* calculate dispersion relation and find eigen vectors */
  kx = 2.0*(PI)*nwave/(x1max-x1min);
  omg2 = SQR(Iso_csound*kx)+2.0*(2.0-qshear)*SQR(Omega_0);
  omg = sqrt(omg2);
  v1x = omg*amp/kx;
  v1y = (2.0-qshear)*Omega_0/omg*v1x;

  time = pGrid->time;

#ifdef PARTICLES
  /* Bin particles to grid */
  particle_to_grid(pGrid, pDomain, property_all);
#endif

  /* Print error to file "Par_LinWave-errors.#.dat", where #=wavedir  */
  fname = ath_fname(NULL,"Par_Shwave1d-errors",0,0,NULL,"dat");
  /* Open output file in write mode */
  FILE *fid = fopen(fname,"w");

  /* Calculate the error and output */
  fprintf(fid, "%f	%d\n", time, ie-is+1);
  for (i=is; i<=ie; i++) {
    /* calculate analytic solution */
    cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);
    SolGasd = 1.0+amp*cos(kx*x1-omg*time);
    SolGasM1 = SolGasd*v1x*cos(kx*x1);
    SolGasM3 = SolGasd*v1y*sin(kx*x1);
    if (samp == 1)
      SolLagd = SolGasd;
    else
      SolLagd = SolGasd-amp*cos(kx*x1);
    /* output */
    fprintf(fid,"%f	",x1);
    fprintf(fid,"%e	%e	%e	",pGrid->U[ks][js][i].d-1.0,
                                SolGasd-1.0,pGrid->U[ks][js][i].d-SolGasd);
#ifdef PARTICLES
    fprintf(fid,"%e	%e	%e	", pG->Coup[ks][js][i].grid_d-1.0,
                        SolLagd-1.0, pG->Coup[ks][js][i].grid_d-SolLagd);
#endif
#ifdef FARGO
//    fprintf(fid,"%e     %e      %e      ",pGrid->U[ks][js][i].d, pGrid->U[ks][js][i].M1, pGrid->U[ks][js][i].M3);
#else
//    fprintf(fid,"%e     %e      %e      ",pGrid->U[ks][js][i].d, pGrid->U[ks][js][i].M1, pGrid->U[ks][js][i].M3+pGrid->U[ks][js][i].d*qshear*Omega_0*x1);
#endif
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++)
      fprintf(fid,"%e	%e	",pGrid->U[ks][js][i].s[n]-1.0,
                                  pGrid->U[ks][js][i].s[n]-SolLagd);
#endif
    fprintf(fid,"\n");
  }

  fclose(fid);

  return;
}


/*=========================== PRIVATE FUNCTIONS ==============================*/
/*--------------------------------------------------------------------------- */

#ifdef PARTICLES
/*! \fn int GetPosition(Grain *gr)
 *  \brief Get particle status (grid/crossing) */
int GetPosition(Grain *gr)
{
  if ((gr->x1>=x1upar) || (gr->x1<x1lpar) || (gr->x2>=x2upar) || (gr->x2<x2lpar) || (gr->x3>=x3upar) || (gr->x3<x3lpar))
    return 10; /* crossing particle */
  else
    return 1;  /* grid particle */
}
#endif

/*------------------------------------------------------------------------------
 * ShearingBoxPot: 
 */

/*! \fn static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3)
 *  \brief shearing box tidal potential */
static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifndef FARGO
  phi -= qshear*SQR(Omega_0*x1);
#endif
  return phi;
}
