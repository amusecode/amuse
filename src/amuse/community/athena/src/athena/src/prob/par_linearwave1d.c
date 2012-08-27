#include "copyright.h"
/*============================================================================*/
/*! \file par_linearwave1d.c
 *  \brief Problem generator for plane-parallel, grid-aligned linear wave tests
 *   with Lagrangian particles.
 *
 * PURPOSE: Problem generator for plane-parallel, grid-aligned linear wave tests
 *   with Lagrangian particles. Works only for 2D grid and with wavevector
 *   parallel to grid. Particles can be initialized either with uniform density
 *   or with the same density profile as gas.
 *
 *  Can be used for either standing (problem/vflow=1.0) or travelling
 *  (problem/vflow=0.0) waves.
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
  int is,ie,js,je,ks,ke,n,wavedir,nwave,samp;
  Real x1,x2,x3,x1max,x1min,x2max,x2min,amp,vflow,kw;
#ifdef PARTICLES
  long p;
  int Npar,ip,jp;
  Real x1p,x2p,x3p,x1l,x1u,x2l,x2u;
  Real par_amp, factor2;
#endif

  if ((par_geti("grid","Nx2") == 1) || (par_geti("grid","Nx3") > 1)) {
    ath_error("[par_linearwave1d]: par_linearwave1d must work in 2D grid.\n");
  }

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

/* Read initial conditions  */
  amp = par_getd("problem","amp");
  wavedir = par_geti("problem","wavedir");
  vflow = par_getd("problem","vflow");
  nwave = par_geti("problem","nwave");
  samp = par_geti("problem","sample");
  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");

  if (wavedir == 1)
    kw = 2.0*(PI)*nwave/(x1max-x1min);
  else if (wavedir == 2)
    kw = 2.0*(PI)*nwave/(x2max-x2min);

/* Now set initial conditions to wave solution */ 

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie; i++) {
    cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

    switch(wavedir){
    case 1:
      pGrid->U[k][j][i].d = 1.0+amp*sin(kw*x1);
      pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d*
                             (vflow+amp*Iso_csound*sin(kw*x1));
      pGrid->U[k][j][i].M2 = 0.0;
      break;

    case 2:
      pGrid->U[k][j][i].d = 1.0+amp*sin(kw*x2);
      pGrid->U[k][j][i].M1 = 0.0;
      pGrid->U[k][j][i].M2 = pGrid->U[k][j][i].d*
                             (vflow+amp*Iso_csound*sin(kw*x2));
      break;

    default:
      ath_error("[par_linearwave1d]: wavedir must be either 1 or 2!\n");
    }

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
  switch(wavedir){
  case 1:
    par_amp = amp*kw*pGrid->dx1/sin(kw*pGrid->dx1);
    factor2 = 0.5*tan(kw*pGrid->dx1)/(kw*pGrid->dx1);
    break;
  case 2:
    par_amp = amp*kw*pGrid->dx2/sin(kw*pGrid->dx2);
    factor2 = 0.5*tan(kw*pGrid->dx2)/(kw*pGrid->dx2);
    break;
  default:
   ath_error("[par_linearwave1d]: wavedir must be either 1 or 2!\n");
  }

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

            pGrid->particle[p].property = 0;

            switch(wavedir){
            case 1:
              pGrid->particle[p].x1 = x1p;
              if (samp == 1) {
                pGrid->particle[p].x1 += par_amp*cos(kw*x1p)/kw
                                      - factor2*SQR(par_amp)*sin(2.0*kw*x1p)/kw;
              }
              pGrid->particle[p].x2 = x2p;
              pGrid->particle[p].v1 = vflow+amp*Iso_csound*sin(kw*x1p);
              pGrid->particle[p].v2 = 0.0;
              break;

            case 2:
              pGrid->particle[p].x1 = x1p;
              pGrid->particle[p].x2 = x2p;
              if (samp == 1) {
                pGrid->particle[p].x2 += par_amp*cos(kw*x2p)/kw
                                      - factor2*SQR(par_amp)*sin(2.0*kw*x2p)/kw;
              }
              pGrid->particle[p].v1 = 0.0;
              pGrid->particle[p].v2 = vflow+amp*Iso_csound*sin(kw*x2p);
              break;

            default:
              ath_error("[par_linearwave1d]: wavedir must be either 1 or 2!\n");
            }

            pGrid->particle[p].x3 = x3p;
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
/*! \fn static Real ScalarDen(const Grid *pG, const int i, const int j, 
 *			      const int k)
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
  int is,ie,js,je,ks,ke,n,wavedir,nwave,samp;
  Real x1,x2,x3,x1max,x1min,x2max,x2min,amp,vflow,kw;
  Real time,omega,SolGasd,SolLagd;
  char *fname;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  /* Read initial conditions  */
  amp = par_getd("problem","amp");
  wavedir = par_geti("problem","wavedir");
  vflow = par_getd("problem","vflow");
  nwave = par_geti("problem","nwave");
  samp = par_geti("problem","sample");
  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");

  /* calculate dispersion relation */
  if (wavedir == 1)
    kw = 2.0*(PI)*nwave/(x1max-x1min);
  else if (wavedir == 2)
    kw = 2.0*(PI)*nwave/(x2max-x2min);

  time = pGrid->time;
  omega = kw*Iso_csound;

  /* Bin particles to grid */
  particle_to_grid(pGrid, pDomain, property_all);

  /* Print error to file "Par_LinWave-errors.#.dat", where #=wavedir  */
  fname = ath_fname(NULL,"Par_LinWave1d-errors",1,wavedir,NULL,"dat");

  /* Open output file in write mode */
  FILE *fid = fopen(fname,"w");

  /* Calculate the error and output */
  switch(wavedir){
  case 1:
    fprintf(fid, "%f	%d\n", time, ie-is+1);
    for (i=is; i<=ie; i++) {
      /* calculate analytic solution */
      cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);
      SolGasd = 1.0+amp*sin(kw*(x1-vflow*time)-omega*time);
      if (samp == 1)
        SolLagd = SolGasd;
      else
        SolLagd = SolGasd-amp*sin(kw*(x1-vflow*time));

      /* output */
      fprintf(fid,"%f	",x1);
      fprintf(fid,"%e	%e	%e	",pGrid->U[ks][js][i].d-1.0,
                                SolGasd-1.0,pGrid->U[ks][js][i].d-SolGasd);
      fprintf(fid,"%e	%e	%e	",pG->Coup[ks][js][i].grid_d-1.0,
                                SolLagd-1.0,pG->Coup[ks][js][i].grid_d-SolLagd);
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        fprintf(fid,"%e	%e	",pGrid->U[ks][js][i].s[n]-1.0,
                                  pGrid->U[ks][js][i].s[n]-SolLagd);
#endif
      fprintf(fid,"\n");
    }
    break;

  case 2:
    fprintf(fid, "%f	%d\n", time, je-js+1);
    for (j=js; j<=je; j++) {
      /* calculate analytic solution */
      cc_pos(pGrid,is,j,ks,&x1,&x2,&x3);
      SolGasd = 1.0+amp*sin(kw*(x2-vflow*time)-omega*time);
      if (samp == 1)
        SolLagd = SolGasd;
      else
        SolLagd = SolGasd-amp*sin(kw*(x2-vflow*time));

      /* output */
      fprintf(fid,"%f   ",x2);
      fprintf(fid,"%e   %e      %e      ",pGrid->U[ks][j][is].d-1.0,
                                SolGasd-1.0,pGrid->U[ks][j][is].d-SolGasd);
      fprintf(fid,"%e   %e      %e      ",pG->Coup[ks][j][is].grid_d-1.0,
                                SolLagd-1.0,pG->Coup[ks][j][is].grid_d-SolLagd);
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        fprintf(fid,"%e %e      ",pGrid->U[ks][j][is].s[n]-1.0,
                        pGrid->U[ks][j][is].s[n]-SolLagd); 
#endif
      fprintf(fid,"\n");
    }
    break; 

  }

  fclose(fid);

  return;
}


/*=========================== PRIVATE FUNCTIONS ==============================*/
/*--------------------------------------------------------------------------- */

/*! \fn int GetPosition(Grain *gr)
 *  \brief get particle status (grid/crossing) */
int GetPosition(Grain *gr)
{
  if ((gr->x1>=x1upar) || (gr->x1<x1lpar) || (gr->x2>=x2upar) || (gr->x2<x2lpar) || (gr->x3>=x3upar) || (gr->x3<x3lpar))
    return 10; /* crossing particle */
  else
    return 1;  /* grid particle */
}
