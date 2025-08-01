#include "copyright.h"
/*============================================================================*/
/*! \file par_friction.c
 *  \brief Problem generator for particle code test, works for 2D and 3D.
 *
 * PURPOSE: Problem generator for particle code test, works for 2D and 3D. The
 *   fluid is set to be at rest. One test particle with initial velocity v0 is
 *   then stopped by the gas. This problem is used to test particle integrator
 *   performance in weak coupling regime.
 *
 * - Configure --with-particle=passive --with-eos=isothermal
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
#error : The friction problem requires particles to be enabled.
#endif /* PARTICLES */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ParticleTroj()    - analytical particle trajectory
 * ParticleVel()     - analytical particle velocity
 * ParticleLocator() - locate the particles (for mpi)
 *============================================================================*/

static Vector ParticleTroj(Real t);
static Vector ParticleVel(Real t);
static int ParticleLocator(Real x1, Real x2, Real x3);

/*------------------------ filewide global variables -------------------------*/
Real x1c,x2c,x3c,v01,v02,v03;
Real x1min,x1max,x2min,x2max,x3min,x3max;
char name[50];

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(Grid *pGrid, Domain *pDomain)
{
  int i,j,k;
  long p,in;


  if (par_geti("grid","Nx1") == 1 || par_geti("grid","Nx2") == 1) {
    ath_error("[par_fric]: this test only works with Nx1 & Nx2 > 1\n");
  }

/* Initialize boxsize */
  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");
  x3min = par_getd("grid","x3min");
  x3max = par_getd("grid","x3max");
  x1c = 0.5*(x1min+x1max);
  x2c = 0.5*(x2min+x2max);
  x3c = 0.5*(x3min+x3max);

/* Read initial conditions for the gas */
  v01 = par_getd("problem","v1");
  v02 = par_getd("problem","v2");
  v03 = par_getd("problem","v3");

/* particle type */
  if (par_geti("particle","partypes") != 1)
    ath_error("[par_fric]: number of particle types must be 1!\n");

/* particle stopping time */
  tstop0 = par_getd("problem","tstop"); /* in code unit */
  if (par_geti("particle","tsmode") != 3)
    ath_error("[par_fric]: This test works only for fixed stopping time!\n");

/* initial particle position */
  in = ParticleLocator(x1c, x2c, x3c);

  pGrid->nparticle         = in;
  pGrid->grproperty[0].num = in;

  if (pGrid->nparticle+2 > pGrid->arrsize)
    particle_realloc(pGrid, pGrid->nparticle+2);

/* Now set initial conditions for the gas */

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
  for (j=pGrid->js; j<=pGrid->je; j++) {
  for (i=pGrid->is; i<=pGrid->ie; i++) {
    pGrid->U[k][j][i].d = 1.0;
    pGrid->U[k][j][i].M1 = 0.0;
    pGrid->U[k][j][i].M2 = 0.0;
    pGrid->U[k][j][i].M3 = 0.0;
  }}}

/* Now set initial conditions for the particles */
  for (p=0; p<in; p++)
  {
    pGrid->particle[p].property = 0;
    pGrid->particle[p].x1 = x1c;
    pGrid->particle[p].x2 = x2c;
    pGrid->particle[p].x3 = x3c;
    pGrid->particle[p].v1 = v01;
    pGrid->particle[p].v2 = v02;
    pGrid->particle[p].v3 = v03;
    pGrid->particle[p].pos = 1; /* grid particle */
    pGrid->particle[p].my_id = p;
#ifdef MPI_PARALLEL
    pGrid->particle[p].init_id = pGrid->my_id;
#endif
  }

  if (pGrid->my_id == 0) {
  /* flush output file */
    sprintf(name, "%s_Err.dat", pGrid->outfilename);
    FILE *fid = fopen(name,"w");
    fclose(fid);
#ifdef MPI_PARALLEL
    sprintf(name, "../%s_Err.dat", pGrid->outfilename);
#else
    sprintf(name, "%s_Err.dat", pGrid->outfilename);
#endif
  }

#ifdef MPI_PARALLEL
  MPI_Bcast(name,50,MPI_CHAR,0,MPI_COMM_WORLD);
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
  fwrite(name, sizeof(char),50,fp);
  return;
}

void problem_read_restart(Grid *pG, Domain *pD, FILE *fp)
{
/* Initialize boxsize */
  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");
  x3min = par_getd("grid","x3min");
  x3max = par_getd("grid","x3max");
  x1c = 0.5*(x1min+x1max);
  x2c = 0.5*(x2min+x2max);
  x3c = 0.5*(x3min+x3max);

/* Read initial conditions for the gas */
  v01 = par_getd("problem","v1");
  v02 = par_getd("problem","v2");
  v03 = par_getd("problem","v3");

  fread(name, sizeof(char),50,fp);

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

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
  long p;
  Real t,ds,dv;
  Vector pos0,vel0;
  Grain *gr;
  FILE *fid;

  t = pGrid->time + pGrid->dt;

  pos0 = ParticleTroj(t);
  vel0 = ParticleVel(t);

  for (p=0; p<pGrid->nparticle; p++)
  {
    if (pGrid->particle[p].pos == 1) /* grid particle */
    {
      gr = &(pGrid->particle[p]);
      ds = sqrt(SQR(pos0.x1-gr->x1)+SQR(pos0.x2-gr->x2)+SQR(pos0.x3-gr->x3));
      dv = sqrt(SQR(vel0.x1-gr->v1)+SQR(vel0.x2-gr->v2)+SQR(vel0.x3-gr->v3));

      fid = fopen(name,"a+");
      fprintf(fid,"%e	%e	%e\n", t, ds, dv);
      fclose(fid);
    }
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

  pos.x1 = x1c+v01*tstop0[0]*(1.0-exp(-t/tstop0[0]));
  pos.x2 = x2c+v02*tstop0[0]*(1.0-exp(-t/tstop0[0]));
  pos.x3 = x3c+v03*tstop0[0]*(1.0-exp(-t/tstop0[0]));

  L1 = x1max-x1min;	L2 = x2max-x2min;	L3 = x3max-x3min;

  pos.x1 = pos.x1-floor((pos.x1-x1min)/L1)*L1;
  pos.x2 = pos.x2-floor((pos.x2-x2min)/L2)*L2;
  pos.x3 = pos.x3-floor((pos.x3-x3min)/L3)*L3;

  return pos;
}

/*! \fn static Vector ParticleVel(Real t)
 *  \brief Compute particle velocity */
static Vector ParticleVel(Real t)
{
  Vector vel;

  vel.x1 = v01*exp(-t/tstop0[0]);
  vel.x2 = v02*exp(-t/tstop0[0]);
  vel.x3 = v03*exp(-t/tstop0[0]);

  return vel;
}

/*! \fn static int ParticleLocator(Real x1, Real x2, Real x3)
 *  \brief Judge if the particle is in this cpu */
static int ParticleLocator(Real x1, Real x2, Real x3)
{
  if ((x1<x1upar) && (x1>=x1lpar) && (x2<x2upar) && (x2>=x2lpar) &&(x3<x3upar) && (x3>=x3lpar))
    return 1;	/* yes */
  else
    return 0;	/* no */
}
