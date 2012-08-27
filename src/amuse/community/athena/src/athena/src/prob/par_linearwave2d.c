#include "copyright.h"
/*============================================================================*/
/*! \file par_linearwave2d.c
 *  \brief Problem generator for oblique linear waves with Lagrangian particles.
 *
 * PURPOSE: Problem generator for oblique linear waves with Lagrangian
 *   particles. Only works for 2D grid. The angle the wave propagates to the
 *   grid is automatically computed to be tan^{-1} (X/Y), so that periodic
 *   boundary conditions can be used. Particles can be initialized either with
 *   uniform density or with the same density profile as gas.
 *
 *   Note angle=0 or 90 [grid aligned waves] is not allowed in this routine.
 *   Use linearwave1d for grid aligned wave on 1D/2D/3D grids.
 *
 *   Can be used for either standing (problem/vflow=1.0) or travelling
 *   (problem/vflow=0.0) waves.
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

#ifndef PARTICLES
#error : The par_linearwave problem requires particles to be enabled.
#endif /* PARTICLES */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * GetPosition()    - get particle status (grid/crossing)
 *============================================================================*/
int GetPosition(Grain *gr);

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(Grid *pGrid, Domain *pDomain)
{
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,n,m,nwave1,nwave2,samp;
  Real x1,x2,x3,x1max,x1min,x2max,x2min,amp,vflow,kx1,kx2,ktot,kr,angle;
#ifdef PARTICLES
  long p;
  int Npar,ip,jp;
  Real x1p,x2p,x3p,x1l,x1u,x2l,x2u;
  Real kd1,kd2,par_amp,factor2;
#endif

  if ((par_geti("grid","Nx2") == 1) || (par_geti("grid","Nx3") > 1)) {
    ath_error("[par_linearwave2d]: par_linearwave1d must work in 2D grid.\n");
  }

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

/* Read initial conditions  */
  amp = par_getd("problem","amp");
  vflow = par_getd("problem","vflow");
  nwave1 = par_geti("problem","nwave1");
  nwave2 = par_geti("problem","nwave2");
  samp = par_geti("problem","sample");
  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");

  if ((nwave1 <= 0) || (nwave2 <= 0))
    ath_error("[par_linearwave2d]: nwave must be positive!\n");

  kx1 = 2.0*(PI)*nwave1/(x1max-x1min);
  kx2 = 2.0*(PI)*nwave2/(x2max-x2min);
  angle = atan(((x1max-x1min)/nwave1)/((x2max-x2min)/nwave2));
  ktot = sqrt(SQR(kx1)+SQR(kx2));

/* Now set initial conditions to wave solution */ 

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie; i++) {
    cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
    kr = kx1*x1+kx2*x2;
    pGrid->U[k][j][i].d = 1.0+amp*sin(kr);
    pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d*
                           (vflow + amp*Iso_csound*sin(kr))*cos(angle);
    pGrid->U[k][j][i].M2 = pGrid->U[k][j][i].d*
                           (vflow + amp*Iso_csound*sin(kr))*sin(angle);
    pGrid->U[k][j][i].M3 = 0.0;
#if (NSCALARS > 0)
    if (samp == 1)
      for (n=0; n<NSCALARS; n++)
        pGrid->U[k][j][i].s[n] = pGrid->U[k][j][i].d;
    else
      for (n=0; n<NSCALARS; n++)
        pGrid->U[k][j][i].s[n] = 1.0;
#endif
  }}}

/* Read initial conditions for the particles */
#ifdef PARTICLES

  /* basic parameters */
  if (par_geti("particle","partypes") != 1)
    ath_error("[par_linwave1d]: This test only allows ONE particle species!\n");

  Npar = (int)(sqrt(par_geti("particle","parnumcell")));
  pGrid->nparticle = Npar*Npar*pGrid->Nx1*pGrid->Nx2;
  pGrid->grproperty[0].num = pGrid->nparticle;
  if (pGrid->nparticle+2 > pGrid->arrsize)
    particle_realloc(pGrid, pGrid->nparticle+2);

  /* particle stopping time */
  tstop0[0] = par_getd_def("problem","tstop",0.0); /* in code unit */
  if (par_geti("particle","tsmode") != 3)
    ath_error("[par_linwave1d]: This test only allows fixed stopping time!\n");

  /* particle perturbation amplitude */
  kd1 = kx1*pGrid->dx1;
  kd2 = kx2*pGrid->dx2;

  par_amp = 4.0*amp*SQR(ktot)/(SQR(kx1)*sin(kd1)/(kd1)*(3.0+cos(kd2))
                               +SQR(kx2)*sin(kd2)/(kd2)*(3.0+cos(kd1)));

  factor2 = (SQR(SQR(kx1)*sin(kd1)/kd1)*(1.0+SQR(cos(kd2)))
             +SQR(kx1*kx2)*sin(2*kd1)*sin(2*kd2)/(kd1*kd2)
             +SQR(SQR(kx2)*sin(kd2)/kd2)*(1.0+SQR(cos(kd1))))
            /(SQR(kx1)*sin(2.0*kd1)/(2.0*kd1)*(3.0+cos(2.0*kd2))
             +SQR(kx2)*sin(2.0*kd2)/(2.0*kd2)*(3.0+cos(2.0*kd1)))/SQR(ktot);
//par_amp=amp;
//factor2 = 0.5;

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

            kr = kx1*x1p + kx2*x2p;

            pGrid->particle[p].property = 0;

            pGrid->particle[p].x1 = x1p;
            pGrid->particle[p].x2 = x2p;
            if (samp == 1)
            {
              pGrid->particle[p].x1 += (par_amp*cos(kr)
                            - factor2*SQR(par_amp)*sin(2.0*kr))*cos(angle)/ktot;
              pGrid->particle[p].x2 += (par_amp*cos(kr)
                            - factor2*SQR(par_amp)*sin(2.0*kr))*sin(angle)/ktot;
            }
            pGrid->particle[p].x3 = x3p;

            pGrid->particle[p].v1 = (vflow+amp*Iso_csound*sin(kr))*cos(angle);
            pGrid->particle[p].v2 = (vflow+amp*Iso_csound*sin(kr))*sin(angle);
            pGrid->particle[p].v3 = 0.0;

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
  return;
}

#if (NSCALARS > 0)
/*! |fn static Real ScalarDen(const Grid *pG, const int i, const int j, 
 *			      const int k)
  * \brief Scalar ratio */
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

Gasfun_t get_usr_expr(const char *expr)
{
#if (NSCALARS > 0)
  if(strcmp(expr,"scalar")==0) return ScalarDen;
#endif
#ifdef PARTICLES
  if(strcmp(expr,"dratio")==0) return dratio;
#endif
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
  int is,ie,js,je,ks,ke,n,nwave1,nwave2,samp,Npar;
  Real x1,x2,x3,x1max,x1min,x2max,x2min,amp,vflow,kx1,kx2,ktot,kr,angle;
  Real time,omega,SolGasd,SolLagd;
  char *fname;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

/* Read initial conditions  */
  amp = par_getd("problem","amp");
  vflow = par_getd("problem","vflow");
  nwave1 = par_geti("problem","nwave1");
  nwave2 = par_geti("problem","nwave2");
  samp = par_geti("problem","sample");
  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");
  Npar = (int)(sqrt(par_geti("particle","parnumcell")));

  /* calculate dispersion relation */
  kx1 = 2.0*(PI)*nwave1/(x1max-x1min);
  kx2 = 2.0*(PI)*nwave2/(x2max-x2min);
  angle = atan(((x1max-x1min)/nwave1)/((x2max-x2min)/nwave2));
  ktot = sqrt(SQR(kx1)+SQR(kx2));

  time = pGrid->time;
  omega = ktot*Iso_csound;

  /* Bin particles to grid */
  particle_to_grid(pGrid, pDomain, property_all);

  /* Print error to file "Par_LinWave-errors.#.dat", where #=wavedir  */
  fname = ath_fname(NULL,"Par_LinWave2d-errors",0,0,NULL,"dat");
  /* Open output file in write mode */
  FILE *fid = fopen(fname,"w");

  fprintf(fid, "%f	%d\n", time, ie-is+1);
  for (i=is; i<=ie; i++) {
    /* calculate analytic solution */
    cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);
    kr = kx1*x1+kx2*x2;
    SolGasd = 1.0+amp*sin(kr-ktot*vflow*time-omega*time);
    if (samp == 1)
      SolLagd = SolGasd;
    else
      SolLagd = SolGasd-amp*sin(kr-ktot*vflow*time);
    /* output */
    fprintf(fid,"%f	",x1);
    fprintf(fid,"%e	%e	%e	",pGrid->U[ks][js][i].d-1.0,
                                SolGasd-1.0,pGrid->U[ks][js][i].d-SolGasd);
    fprintf(fid,"%e	%e	%e	",
                 pG->Coup[ks][js][i].grid_d/SQR(Npar)-1.0,
                 SolLagd-1.0,pG->Coup[ks][js][i].grid_d/SQR(Npar)-SolLagd);
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++)
      fprintf(fid,"%e	%e	",pGrid->U[ks][js][i].s[n]-1.0,
                                  pGrid->U[ks][js][i].s[n]-SolLagd);
#endif
    fprintf(fid,"\n");
  }

  return;
}


/*=========================== PRIVATE FUNCTIONS ==============================*/
/*--------------------------------------------------------------------------- */
/*! \fn int GetPosition(Grain *gr)
 *  \brief get particle status (grid/crossing)*/
int GetPosition(Grain *gr)
{
  if ((gr->x1>=x1upar) || (gr->x1<x1lpar) || (gr->x2>=x2upar) || (gr->x2<x2lpar)
                       || (gr->x3>=x3upar) || (gr->x3<x3lpar))
    return 10; /* crossing particle */
  else
    return 1;  /* grid particle */
}
