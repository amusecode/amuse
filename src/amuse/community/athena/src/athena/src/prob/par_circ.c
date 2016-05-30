#include "copyright.h"
/*============================================================================*/
/*! \file par_circ.c
 *  \brief Problem generator for particle code test, works for 2D or 3D.
 *
 * PURPOSE: Problem generator for particle code test, works for 2D or 3D. The 
 *   gas is set tobe steady state, with a circular motion around the center of
 *   the domain. This is not done by playing with gas dynamics, but by setting
 *   the gas at rest, and properly treat gasvshift function. A particle with
 *   zero stopping is initiated and serve as Lagrangian particle to trace fluid
 *   element trajectories. This test is used to test particle integrator
 *   performance in strong coupling regime.
 *
 *   Configure --with-particle=passive --with-eos=isothermal
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
#error : The circular motion problem requires particles to be enabled.
#endif /* PARTICLES */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ParticleTroj()    - analytical particle trajectory
 * ParticleVel()     - analytical particle velocity
 * ParticleLocator() - locate the particles (for mpi)
 * ran2()            - random number generator
 *============================================================================*/

static Vector ParticleTroj(Real t);
static Vector ParticleVel(Vector pos);
static int ParticleLocator(Real x1, Real x2, Real x3);
double ran2(long int *idum);
extern Vector Get_Term(Grid *pG, int type, Real x1, Real x2, Real x3,
                                 Vector cell1,      Real *tstop);
/*------------------------ filewide global variables -------------------------*/
Real x1c,x2c,x3c;
Real x01,x02,x03;
Real omg, omgx1, omgx2, omgx3;	/* angular velocity and orientation */
Real r01,r02,r03,r0dn,r0cn1,r0cn2,r0cn3,n1,n2,n3;
char name[50];

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(Grid *pGrid, Domain *pDomain)
{
  int i,j,k,in;
  long p, seed;
  Real x1min,x1max,x2min,x2max,x3min,x3max;
  Real ranv1, ranv2, ranv3, ranvnorm, tstop;
  Real theta, phi, rad, vran, thetaran, phiran;
  Vector cell1, vterm;

  if (par_geti("grid","Nx1") == 1 || par_geti("grid","Nx2") == 1) {
    ath_error("[par_circ]: this test only works with Nx1 & Nx2 > 1\n");
  }

  /* cell1 is a shortcut expressions as well as dimension indicator */
  if (pGrid->Nx1 > 1)  cell1.x1 = 1.0/pGrid->dx1;  else cell1.x1 = 0.0;
  if (pGrid->Nx2 > 1)  cell1.x2 = 1.0/pGrid->dx2;  else cell1.x2 = 0.0;
  if (pGrid->Nx3 > 1)  cell1.x3 = 1.0/pGrid->dx3;  else cell1.x3 = 0.0;

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

/* Read initial conditions for gas */
  omg = par_getd("problem","omega");
  theta = par_getd("problem","theta");
  phi = par_getd("problem","phi");
  rad = par_getd("problem","rad");
  vran = par_getd("problem","vran");
  seed = par_geti_def("problem","seed",0);

  omgx1 = omg*sin(theta)*cos(phi);
  omgx2 = omg*sin(theta)*sin(phi);
  omgx3 = omg*cos(theta);
  n1 = omgx1/omg;
  n2 = omgx2/omg;
  n3 = omgx3/omg;

/* particle type */
  if (par_geti("particle","partypes") != 1)
    ath_error("[par_circ]: number of particle types must be 1!\n");

/* particle stopping time */
  tstop0[0] = par_getd_def("problem","tstop",0.0); /* in code unit */
  if (par_geti("particle","tsmode") != 3)
    ath_error("[par_circ]: This test works only for fixed stopping time!\n");

/* particle position relative to domain center */
  r01 = rad*cos(theta)*cos(phi);
  r02 = rad*cos(theta)*sin(phi);
  r03 = -rad*sin(theta);

/* initial particle position */
  x01 = x1c+r01;
  x02 = x2c+r02;
  x03 = x3c+r03;

  r0dn = r01*n1+r02*n2+r03*n3;/* r0 dot n */
  r0cn1 = n2*r03-n3*r02;      /* r0 cross n */
  r0cn2 = n3*r01-n1*r03;
  r0cn3 = n1*r02-n2*r01;

  in = ParticleLocator(x01, x02, x03);

  pGrid->nparticle         = in;
  pGrid->grproperty[0].num = in;

  if (pGrid->nparticle+2 > pGrid->arrsize)
    particle_realloc(pGrid, pGrid->nparticle+2);

/* initial particle random velocity */
  ranv1 = ran2(&seed)-0.5;
  ranv2 = ran2(&seed)-0.5;
  ranv3 = ran2(&seed)-0.5;
  ranvnorm = vran/sqrt(SQR(ranv1)+SQR(ranv2)+SQR(ranv3)+1.0e-8);
  ranv1 = ranvnorm*ranv1;
  ranv2 = ranvnorm*ranv2;
  ranv3 = ranvnorm*ranv3;

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
    pGrid->particle[p].x1 = x01;
    pGrid->particle[p].x2 = x02;
    pGrid->particle[p].x3 = x03;

    vterm = Get_Term(pGrid,0,x01,x02,x03,cell1,&tstop);

//    pGrid->particle[p].v1 = omgx2*r03-omgx3*r02+ranv1;
//    pGrid->particle[p].v2 = omgx3*r01-omgx1*r03+ranv2;
//    pGrid->particle[p].v3 = omgx1*r02-omgx2*r01+ranv3;
    pGrid->particle[p].v1 = ranv1+vterm.x1;
    pGrid->particle[p].v2 = ranv2+vterm.x2;
    pGrid->particle[p].v3 = ranv3+vterm.x3;

    pGrid->particle[p].pos = 1; /* grid particle */
    pGrid->particle[p].my_id = p;
#ifdef MPI_PARALLEL
    pGrid->particle[p].init_id = pGrid->my_id;
#endif
  }

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
 * get_usr_par_prop()      - returns a user defined particle selection function
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
  Real x1min,x1max,x2min,x2max,x3min,x3max;
  Real theta, phi, rad;

  x1min = par_getd("grid","x1min");  x1max = par_getd("grid","x1max");
  x2min = par_getd("grid","x2min");  x2max = par_getd("grid","x2max");
  x3min = par_getd("grid","x3min");  x3max = par_getd("grid","x3max");
  x1c = 0.5*(x1min+x1max);  x2c = 0.5*(x2min+x2max);  x3c = 0.5*(x3min+x3max);

  omg = par_getd("problem","omega");
  theta = par_getd("problem","theta");
  phi = par_getd("problem","phi");
  rad = par_getd("problem","rad");

  omgx1 = omg*sin(theta)*cos(phi);  r01 = rad*cos(theta)*cos(phi);
  omgx2 = omg*sin(theta)*sin(phi);  r02 = rad*cos(theta)*sin(phi);
  omgx3 = omg*cos(theta);           r03 = -rad*sin(theta);

  n1 = omgx1/omg;  n2 = omgx2/omg;  n3 = omgx3/omg;
  x01 = x1c+r01;  x02 = x2c+r02;  x03 = x3c+r03;

  r0dn = r01*n1+r02*n2+r03*n3;/* r0 dot n */
  r0cn1 = n2*r03-n3*r02;      /* r0 cross n */
  r0cn2 = n3*r01-n1*r03;
  r0cn3 = n1*r02-n2*r01;

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

/*! \fn void gasvshift(const Real x1, const Real x2, const Real x3,
 *                                           Real *u1, Real *u2, Real *u3)
 *  \brief Gas velocity shift */
void gasvshift(const Real x1, const Real x2, const Real x3,
                                             Real *u1, Real *u2, Real *u3)
{
  Real dx1, dx2, dx3;
  dx1 = x1-x1c;          dx2 = x2-x2c;          dx3 = x3-x3c;
  *u1 = omgx2*dx3-omgx3*dx2;
  *u2 = omgx3*dx1-omgx1*dx3;
  *u3 = omgx1*dx2-omgx2*dx1;

  return;
}

/*! \fn void Userforce_particle(Vector *ft, const Real x1, const Real x2, 
 *		      const Real x3,
 *  \brief User-supplied particle force */
void Userforce_particle(Vector *ft, const Real x1, const Real x2, const Real x3,
                                    const Real v1, const Real v2, const Real v3)
{
  return;
}
#endif

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
  long p;
  Real t,ds,dv,r1,r0,dr;
  Vector pos0, vel0, pos;
  Grain *gr;
  FILE *fid;

  t = pGrid->time+pGrid->dt;
  pos0 = ParticleTroj(t);

  for (p=0; p<pGrid->nparticle; p++)
  {
    gr = &(pGrid->particle[p]);
    if (gr->pos == 1) /* grid particle */
    {
      pos.x1 = gr->x1;
      pos.x2 = gr->x2;
      pos.x3 = gr->x3;
      vel0 = ParticleVel(pos);

      ds = sqrt(SQR(pos0.x1-gr->x1)+SQR(pos0.x2-gr->x2)+SQR(pos0.x3-gr->x3));
      dv = sqrt(SQR(vel0.x1-gr->v1)+SQR(vel0.x2-gr->v2)+SQR(vel0.x3-gr->v3));
      r0 = sqrt(SQR(pos0.x1-x1c)+SQR(pos0.x2-x2c)+SQR(pos0.x3-x3c));
      r1 = sqrt(SQR(gr->x1-x1c)+SQR(gr->x2-x2c)+SQR(gr->x3-x3c));
      dr = r1-r0;

      fid = fopen(name,"a+");
//      fprintf(fid,"%e	%e	%e	%e	%e	%e	%e	%e	%e\n", t, ds, dv, gr->x1, gr->x2, gr->x3, pos0.x1, pos0.x2, pos0.x3);
      fprintf(fid,"%e   %e      %e\n", t, ds, dr);
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
  Real r1,r2,r3,omgt,sinot, cosot;

  omgt = omg*t;
  sinot = sin(omgt);	cosot = cos(omgt);
  r1 = r0dn*n1 + sinot*r0cn1 + cosot*(r01-r0dn*n1);
  r2 = r0dn*n2 + sinot*r0cn2 + cosot*(r02-r0dn*n2);
  r3 = r0dn*n3 + sinot*r0cn3 + cosot*(r03-r0dn*n3);

  pos.x1 = x1c+r1;     pos.x2 = x2c+r2;        pos.x3 = x3c+r3;
  return pos;
}

/*! \fn static Vector ParticleVel(Vector pos)
 *  \brief Compute particle velocity */
static Vector ParticleVel(Vector pos)
{
  Vector vel;
  Real r1,r2,r3;

  r1 = pos.x1-x1c;	r2 = pos.x2-x2c;	r3 = pos.x3-x3c;
  vel.x1 = omgx2*r3-omgx3*r2;
  vel.x2 = omgx3*r1-omgx1*r3;
  vel.x3 = omgx1*r2-omgx2*r1;

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

/*----------------------------------------------------------------------------*/

#define DBL_EPSILON 1.0e-16
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

/*! \fn double ran2(long int *idum)
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

double ran2(long int *idum)
{
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
#undef DBL_EPSILON
