#include "copyright.h"
/*============================================================================*/
/*! \file par_collision.c
 *  \brief Problem generator for particle feedback test in 2D. 
 *
 * PURPOSE: Problem generator for particle feedback test in 2D. The particles
 *   and gas are initialized to have the same total momentum in the opposite
 *   direction. Both have uniform density. One test particle is picked up for
 *   precision test.
 *
 * - Configure --with-particle=feedback --with-eos=isothermal
 *
 * USERWORK_IN_LOOP function is used to output particle positions.
 */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "particles/particle.h"

#ifndef PARTICLES
#error : The feedback problem requires particles to be enabled.
#endif /* PARTICLES */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ParticleTroj()   - analytical particle trajectory
 * ParticleVel()    - analytical particle velocity
 *============================================================================*/
static Vector ParticleTroj(Real t);
static Vector ParticleVel(Real t);

/*------------------------ filewide global variables -------------------------*/
Real mratio,wx,wy,wz;
Real x1min,x1max,x2min,x2max;
Real x0[3]; /* test particle initial position */
int Npar;
long idlab; /* test particle id */
int cpuid;  /* test particle cpu id */
char name[50];

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(Grid *pGrid, Domain *pDomain)
{
  int i,j,ip,jp;
  long p,lab;
  Real x1,x2,x3,x1l,x1u,x2l,x2u,x1p,x2p,x3p,cosphi,sinphi;
  Real vdif_p,vdif_v,w,u,ux,uy,uz;

  if (par_geti("grid","Nx2") == 1) {
    ath_error("[par_collision]: this test only works for 2D problem.\n");
  }
  if (par_geti("grid","Nx3") > 1) {
    ath_error("[par_collision]: this test only works for 2D problem.\n");
  }

  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");

/* Read the velocity difference */
  vdif_p = par_getd("problem","vdif_p");
  vdif_v = par_getd("problem","vdif_v");
  cosphi = par_getd("problem","oblique");
  sinphi = sqrt(1.0-SQR(cosphi));

/* Read particle parameters */

  /* basic parameters */
  if (par_geti("particle","partypes") != 1)
    ath_error("[par_collision]: This test only allows ONE particle species!\n");

  Npar = (int)(sqrt(par_geti("particle","parnumcell")));
  pGrid->nparticle = Npar*Npar*pGrid->Nx1*pGrid->Nx2;
  pGrid->grproperty[0].num = pGrid->nparticle;
  if (pGrid->nparticle+2 > pGrid->arrsize)
    particle_realloc(pGrid, pGrid->nparticle+2);

  /* particle stopping time */
  tstop0[0] = par_getd("problem","tstop"); /* in code unit */
  if (par_geti("particle","tsmode") != 3)
    ath_error("[par_collision]: This test only allow fixed stopping time!\n");

  /* assign particle effective mass */
#ifdef FEEDBACK
  mratio = par_getd_def("problem","mratio",0.0); /* total mass fraction */
  if (mratio < 0.0)
    ath_error("[par_collision]: mratio must be positive!\n");
  pGrid->grproperty[0].m = mratio/SQR(Npar);
#else
  mratio = 0.0;
#endif

  w = vdif_p/(1.0+mratio); /* dust velocity */
  u = -mratio*w;           /* gas velocity */

  wx = w * cosphi;
  wy = w * sinphi;
  ux = u * cosphi;
  uy = u * sinphi;
  wz = vdif_v/(1.0+mratio);
  uz = -mratio*wz;

/* Now set initial conditions for the gas */

  for (j=pGrid->js; j<=pGrid->je; j++) {
  for (i=pGrid->is; i<=pGrid->ie; i++) {
    cc_pos(pGrid,i,j,pGrid->ks,&x1,&x2,&x3);
    pGrid->U[pGrid->ks][j][i].d = 1.0;

    pGrid->U[pGrid->ks][j][i].M1 = ux;
    pGrid->U[pGrid->ks][j][i].M2 = uy;
    pGrid->U[pGrid->ks][j][i].M3 = uz;
  }}

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
            pGrid->particle[p].x3 = x3p;

            pGrid->particle[p].v1 = wx;
            pGrid->particle[p].v2 = wy;
            pGrid->particle[p].v3 = wz;

            pGrid->particle[p].pos = 1; /* grid particle */
            pGrid->particle[p].my_id = p;
#ifdef MPI_PARALLEL
            pGrid->particle[p].init_id = pGrid->my_id;
#endif
            p += 1;
          }
        }
    }
  }

  /* label a test particle */
  if (pGrid->my_id == 0) {
    lab = 0;
    x0[0] = pGrid->particle[lab].x1;
    x0[1] = pGrid->particle[lab].x2;
    x0[2] = pGrid->particle[lab].x3;
    idlab = pGrid->particle[lab].my_id;
    cpuid = 0;
  }

  if (pGrid->my_id == 0) {
  /* flush output file */
    sprintf(name, "%s_partroj.dat", pGrid->outfilename);
    FILE *fid = fopen(name,"w");
    fclose(fid);
#ifdef MPI_PARALLEL
    sprintf(name, "../%s_partroj.dat", pGrid->outfilename);
#else
    sprintf(name, "%s_partroj.dat", pGrid->outfilename);
#endif
  }

#ifdef MPI_PARALLEL
  MPI_Bcast(name,50,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(&x0,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&idlab,1,MPI_LONG,0,MPI_COMM_WORLD);
  MPI_Bcast(&cpuid,1,MPI_INT,0,MPI_COMM_WORLD);
#endif

  return;
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
  fwrite(&wx, sizeof(Real),1,fp);
  fwrite(&wy, sizeof(Real),1,fp);
  fwrite(&wz, sizeof(Real),1,fp);
  fwrite(x0, sizeof(Real),3,fp);
  fwrite(&idlab, sizeof(long),1,fp);
  fwrite(&cpuid, sizeof(int),1,fp);
  return;
}

void problem_read_restart(Grid *pG, Domain *pD, FILE *fp)
{
  fread(&wx, sizeof(Real),1,fp);
  fread(&wy, sizeof(Real),1,fp);
  fread(&wz, sizeof(Real),1,fp);
  fread(x0, sizeof(Real),3,fp);
  fread(&idlab, sizeof(long),1,fp);
  fread(&cpuid, sizeof(int),1,fp);

  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");
  Npar = (int)(sqrt(par_geti("particle","parnumcell")));
#ifdef FEEDBACK
  mratio = par_getd_def("problem","mratio",0.0); /* total mass fraction */
#else
  mratio = 0.0;
#endif

  sprintf(name, "%s_partroj.dat", pG->outfilename);
  return;
}

Gasfun_t get_usr_expr(const char *expr)
{
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

/* calculate particle analytical positions and compare */
void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
  long p, lab;
  Grain *gr;
  Real Mg1,Mg2,Mg3,mp,t,s,ds,v,dv,Mtot;
  Vector pos0, vel0;
  FILE *fid;

  p = 0;  lab = -1;
  while (p<pGrid->nparticle) {
#ifdef MPI_PARALLEL
    if ((pGrid->particle[p].my_id == idlab) && 
        (pGrid->particle[p].init_id == cpuid))
#else
    if (pGrid->particle[p].my_id == idlab)
#endif
    {
      lab = p;
      break;
    }
    p++;
  }

  if (lab >= 0)
  {/* if particle is in the grid */

    /* analytic position */
    t = pGrid->time + pGrid->dt;
    pos0 = ParticleTroj(t);
    vel0 = ParticleVel(t);
    Mg1 = pGrid->U[pGrid->ks][pGrid->js][pGrid->is].M1;
    Mg2 = pGrid->U[pGrid->ks][pGrid->js][pGrid->is].M2;
    Mg3 = pGrid->U[pGrid->ks][pGrid->js][pGrid->is].M3;
    gr = &(pGrid->particle[lab]);
#ifdef FEEDBACK
    mp = pGrid->grproperty[0].m*SQR(Npar);
#else
    mp = SQR(Npar);
#endif

    /* output quantities */
    s  = sqrt(SQR(gr->x1-x0[0])+SQR(gr->x2-x0[1])+SQR(gr->x3-x0[2]));
    ds = sqrt(SQR(gr->x1-pos0.x1)+SQR(gr->x2-pos0.x2)+SQR(gr->x3-pos0.x3));
    v  = sqrt(SQR(gr->v1)+SQR(gr->v2)+SQR(gr->v3));
    dv = sqrt(SQR(gr->v1-vel0.x1)+SQR(gr->v2-vel0.x2)+SQR(gr->v3-vel0.x3));
    Mtot = sqrt(SQR(mp*gr->v1+Mg1)+SQR(mp*gr->v2+Mg2)+SQR(mp*gr->v3+Mg3));

    /* output */
    fid = fopen(name,"a+");
    fprintf(fid,"%e	%e	%e	%e	%e	%e\n", 
                  t,     s,     ds,     v,       dv,   Mtot);
    fclose(fid);
  }

  return;
}

/*---------------------------------------------------------------------------
 * Userwork_after_loop: computes L1-error in linear waves,
 * ASSUMING WAVE HAS PROPAGATED AN INTEGER NUMBER OF PERIODS
 * Must set parameters in input file appropriately so that this is true
 */

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
  return;
}


/*=========================== PRIVATE FUNCTIONS ==============================*/
/*--------------------------------------------------------------------------- */
/*! \fn static Vector ParticleTroj(Real t)
 *  \brief Compute particle trajectory */
static Vector ParticleTroj(Real t)
{
  Vector pos;
  Real L1,L2,L3;
  Real factor=tstop0[0]/(1.0+mratio)*(1.0-exp(-(1.0+mratio)*t/tstop0[0]));

  pos.x1 = x0[0] + wx*factor;
  pos.x2 = x0[1] + wy*factor;
  pos.x3 = x0[2];

  L1 = x1max-x1min;	L2 = x2max-x2min;

  pos.x1 = pos.x1-floor((pos.x1-x1min)/L1)*L1;
  pos.x2 = pos.x2-floor((pos.x2-x2min)/L2)*L2;

  return pos;
}

/*! \fn static Vector ParticleVel(Real t)
 *  \brief Compute particle velocity */
static Vector ParticleVel(Real t)
{
  Vector vel;
  Real factor = exp(-(1.0+mratio)*t/tstop0[0]);

  vel.x1 = wx*factor;
  vel.x2 = wy*factor;
  vel.x3 = wz*factor;

  return vel;
}
