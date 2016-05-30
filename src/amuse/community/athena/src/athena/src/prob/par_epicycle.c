#include "copyright.h"
/*============================================================================*/
/*! \file par_epicycle.c
 *  \brief Problem generator for particle epicycle trajectory presicion test.
 *
 * PURPOSE: Problem generator for particle epicycle trajectory presicion test.
 *   This code works for both 2D and 3D. No gas is involved, but gas has to be
 *   initialized anyway. The particle stopping time should be set to be
 *   sufficiently large so that only shear terms affect particle motion.
 *
 *  Should be configured using --enable-shearing-box and --with-eos=isothermal.
 *  Optional choices are --enable-fargo and --enable-mpi
 */
/*============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "particles/particle.h"

#ifndef SHEARING_BOX
#error : The epicycle problem requires shearing-box to be enabled.
#endif /* SHEARING_BOX */

#ifndef PARTICLES
#error : The epicycle problem requires particles to be enabled.
#endif /* PARTICLES */

#ifndef ISOTHERMAL
#error : The epicycle problem requires isothermal equation of state.
#endif /* ISOTHERMAL */


/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ShearingBoxPot()   - shearing box tidal gravitational potential
 * ParticlePosition() - analytical particle trajectory
 * ParticleVelocity() - analytical particle velocity
 * ParticleLocator()  - locate the particles (for mpi)
 *============================================================================*/

static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3);
static Vector ParticlePosition(const Real t);
static Vector ParticleVelocity(const Vector pos, const Real t);
static int ParticleLocator(const Vector pos);

/*------------------------ filewide global variables -------------------------*/
char name[50];
/* trajectory variables */
Real amp, Lx, Ly, x1min, x1max, x2min, x2max, omg;


/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(Grid *pGrid, Domain *pDomain)
{
  int in,i,j,k;
  Real x1,x2,x3;
  long p;
  Vector parpos, parvel;

  if (par_geti("grid","Nx2") == 1) {
    ath_error("[par_epicycle]: par_epicycle must work in 2D or 3D.\n");
  }

/* Initialize boxsize */
  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  Lx = x1max - x1min;

  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");
  Ly = x2max - x2min;	/* for 3D problem */
  if (par_geti("grid","Nx3") == 1) {
    Ly = 0.0;
  }

/* Read initial conditions */
  Omega_0 = par_getd("problem","omega");
  qshear  = par_getd_def("problem","qshear",1.5);
  amp = par_getd("problem","amp");
  omg = sqrt(2.0*(2.0-qshear))*Omega_0;

/* particle type */
  if (par_geti("particle","partypes") != 1)
    ath_error("[par_epicycle]: This test only allows ONE particle species!\n");

/* particle stopping time */
  tstop0[0] = par_getd_def("problem","tstop",1.0e20); /* in code unit */
  if (par_geti("particle","tsmode") != 3)
    ath_error("[par_epicycle]: This test only allows fixed stopping time!\n");

/* particle position */
  parpos = ParticlePosition(0.0);
  parvel = ParticleVelocity(parpos, 0.0);
  in = ParticleLocator(parpos);

  pGrid->nparticle         = in;
  pGrid->grproperty[0].num = in;

  if (pGrid->nparticle+2 > pGrid->arrsize)
    particle_realloc(pGrid, pGrid->nparticle+2);

/* Now set initial conditions for the gas */
  for (k=pGrid->ks; k<=pGrid->ke; k++) {
  for (j=pGrid->js; j<=pGrid->je; j++) {
  for (i=pGrid->is; i<=pGrid->ie; i++) {
    cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
    pGrid->U[k][j][i].d = 1.0;
    pGrid->U[k][j][i].M1 = 0.0;
    pGrid->U[k][j][i].M2 = 0.0;
    pGrid->U[k][j][i].M3 = 0.0;
#ifndef FARGO
    if (Ly>0.0) /* 3D */
      pGrid->U[k][j][i].M2 -= qshear*Omega_0*x1;
    else /* 2D */
      pGrid->U[k][j][i].M3 -= qshear*Omega_0*x1;
#endif
  }}}

/* Now set initial conditions for the particles */
  for (p=0; p<in; p++)
  {
    pGrid->particle[p].property = 0;
    pGrid->particle[p].x1 = parpos.x1;
    pGrid->particle[p].x2 = parpos.x2;
    pGrid->particle[p].x3 = parpos.x3;
    pGrid->particle[p].v1 = parvel.x1;
    pGrid->particle[p].v2 = parvel.x2;
    pGrid->particle[p].v3 = parvel.x3;
    pGrid->particle[p].pos = 1; /* grid particle */
    pGrid->particle[p].my_id = p;
#ifdef MPI_PARALLEL
    pGrid->particle[p].init_id = pGrid->my_id;
#endif
  }

/* enroll gravitational potential function, shearing sheet BC functions */
  StaticGravPot = ShearingBoxPot;

  if (pGrid->my_id == 0) {
  /* flush output file */
    sprintf(name, "%s_Traj.dat", pGrid->outfilename);
    FILE *fid = fopen(name,"w");
    fclose(fid);
#ifdef MPI_PARALLEL
    sprintf(name, "../%s_Traj.dat", pGrid->outfilename);
#else
    sprintf(name, "%s_Traj.dat", pGrid->outfilename);
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
  Omega_0 = par_getd("problem","omega");
  qshear  = par_getd_def("problem","qshear",1.5);

  StaticGravPot = ShearingBoxPot;

  Ly = x2max - x2min;	/* for 3D problem */
  if (par_geti("grid","Nx3") == 1) {
    Ly = 0.0;
  }

  amp = par_getd("problem","amp");
  omg = sqrt(2.0*(2.0-qshear))*Omega_0;

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
  Real t, ds, E, vshift;
  Vector pos0;
  Grain *gr;
  FILE *fid;

  t = pGrid->time+pGrid->dt;
  pos0 = ParticlePosition(t);
  for (p=0; p<pGrid->nparticle; p++)
  {
    gr = &(pGrid->particle[p]);
    if (gr->pos == 1) /* grid particle */
    {
      /* position error */
      ds = sqrt(SQR(gr->x1-pos0.x1)+SQR(gr->x2-pos0.x2)+SQR(gr->x3-pos0.x3));

      /* total energy */
      E = 0.5*SQR(gr->v1) - qshear*SQR(Omega_0*gr->x1);
#ifdef FARGO
      vshift = -qshear*Omega_0*gr->x1;
#else
      vshift = 0.0;
#endif
      if (Ly>0.0)
        E += 0.5*SQR(gr->v2+vshift);
      else
        E += 0.5*SQR(gr->v3+vshift);

      /* output */
      fid = fopen(name,"a+");
      fprintf(fid,"%e	%e	%e	%e	%e	%e	%e	%e	%e\n", t, ds, E, gr->x1, gr->x2, gr->x3, pos0.x1, pos0.x2, pos0.x3);
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
/*! \fn static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3)
 *  \brief shearing box tidal gravitational potential*/
static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifndef FARGO
  phi -= qshear*Omega_0*Omega_0*x1*x1;
#endif
  return phi;
}

/*! \fn static Vector ParticlePosition(const Real t)
 *  \brief Calculate the particle position */
static Vector ParticlePosition(const Real t)
{
  Real x,y;
  Vector pos;
  x = amp*cos(omg*t);
  y = -2.0*amp*Omega_0/omg*sin(omg*t);

  x = x-floor((x-x1min)/Lx)*Lx;
  if (Ly > 0.0)
    y = y-floor((y-x2min)/Ly)*Ly;
  else
    y = 0.0;

  pos.x1 = x;	pos.x2 = y;	pos.x3 = 0.0;
  return pos;
}

/*! \fn static Vector ParticleVelocity(const Vector pos, const Real t)
 *  \brief Calculate the particle velocity */
static Vector ParticleVelocity(const Vector pos, const Real t)
{
  Real vx,vy;
  Vector vel;
  vx = -amp*omg*sin(omg*t);
  vy = -2.0*amp*Omega_0*cos(omg*t);
#ifdef FARGO
  vy = vy + qshear*Omega_0*pos.x1;
#endif
  if (Ly>0.0) {
    vel.x1 = vx;	vel.x2 = vy;	vel.x3 = 0.0;
  } else {
    vel.x1 = vx;	vel.x3 = vy;	vel.x2 = 0.0;
  }
  return vel;
}

/*! \fn static int ParticleLocator(const Vector pos)
 *  \brief Judge if the particle is in this cpu */
static int ParticleLocator(const Vector pos)
{
  if ((pos.x1<x1upar) && (pos.x1>=x1lpar) && (pos.x2<x2upar)
      && (pos.x2>=x2lpar) &&(pos.x3<x3upar) && (pos.x3>=x3lpar))
    return 1;	/* yes */
  else
    return 0;	/* no */
}
