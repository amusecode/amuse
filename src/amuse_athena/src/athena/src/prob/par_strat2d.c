#include "copyright.h"
/*============================================================================*/
/*! \file par_strat2d.c
 *  \brief Problem generator for non-linear streaming instability in
 *   stratified disks.
 *
 * PURPOSE: Problem generator for non-linear streaming instability in
 *   stratified disks. This code works in 2D ONLY. Isothermal eos is assumed,
 *   and the value etavk/iso_sound is fixed. MPI domain decomposition in x is
 *   allowed, but not in z.
 *
 * Perturbation modes:
 *  - ipert = 0: multi-nsh equilibtium
 *  - ipert = 1: white noise within the entire grid
 *  - ipert = 2: non-nsh velocity
 *
 *  Must be configured using --enable-shearing-box and --with-eos=isothermal.
 *  FARGO is recommended.
 *
 * Reference:
 *   Johansen & Youdin, 2007, ApJ, 662, 627
 *   Bai & Stone, 2009, in preparation					      */
/*============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "particles/particle.h"

#ifndef SHEARING_BOX
#error : The streaming2d problem requires shearing-box to be enabled.
#endif /* SHEARING_BOX */

#ifndef PARTICLES
#error : The streaming2d problem requires particles to be enabled.
#endif /* PARTICLES */

#ifndef ISOTHERMAL
#error : The streaming2d problem requires isothermal equation of state.
#endif /* ISOTHERMAL */

/*------------------------ filewide global variables -------------------------*/
/* flow properties */
Real vsc1,vsc2;
int ipert;
/* domain size variables */
Real x1min,x1max,x2min,x2max,Lx,Lz,Lg;
long Npar;
/* output variables */
long ntrack;   /* number of tracer particles */
long nlis;     /* number of output particles for list output */
int mytype;    /* specific particle type for output */
Real dpar_thresh; /* threshold particle density */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()            - random number generator
 * Normal()          - normal distribution generator
 * Erf()             - error function
 * MultiNSH()        - multiple component NSH equilibrium solver
 * ShearingBoxPot()  - shearing box tidal gravitational potential
 * hst_rho_Vx_dVy()  - total Reynolds stress for history dump
 * property_???()    - particle property selection function
 *============================================================================*/
double ran2(long int *idum);
double Normal(long int *idum);
Real Erf(Real z);

void MultiNSH(int n, Real *tstop, Real *mratio, Real etavk,
                     Real *uxNSH, Real *uyNSH,  Real *wxNSH, Real *wyNSH);
static Real hst_rho_Vx_dVy(const Grid *pG,const int i,const int j,const int k);
static void close_ix2(Grid *pGrid);
static void close_ox2(Grid *pGrid);
static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3);

static int property_limit(const Grain *gr, const GrainAux *grsub);
static int property_trace(const Grain *gr, const GrainAux *grsub);
static int property_type(const Grain *gr, const GrainAux *grsub);
static int property_dense(const Grain *gr, const GrainAux *grsub);

extern Real expr_dpar(const Grid *pG, const int i, const int j, const int k);
extern Real expr_V1par(const Grid *pG, const int i, const int j, const int k);
extern Real expr_V2par(const Grid *pG, const int i, const int j, const int k);
extern Real expr_V3par(const Grid *pG, const int i, const int j, const int k);
extern Real expr_V2(const Grid *pG, const int i, const int j, const int k);

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(Grid *pGrid, Domain *pDomain)
{
  int i,j,js,pt,tsmode,vertical_bc;
  long p,q;
  Real ScaleHg,tsmin,tsmax,tscrit,amin,amax,Hparmin,Hparmax;
  Real *ep,*ScaleHpar,epsum,mratio,pwind,rhoaconv,etavk;
  Real *epsilon,*uxNSH,*uyNSH,**wxNSH,**wyNSH;
  Real rhog,h,x1,x2,x3,t,x1p,x2p,x3p,zmin,zmax,dx2_1,b;
  long int iseed = -1; /* Initialize on the first call to ran2 */

  if (par_geti("grid","Nx2") == 1) {
    ath_error("[par_strat2d]: par_strat2d only works for 2D problem.\n");
  }
  if (par_geti("grid","Nx3") > 1) {
    ath_error("[par_strat2d]: par_strat2d only works for 2D problem.\n");
  }
#ifdef MPI_PARALLEL
  if (par_geti("parallel","NGrid_x2") > 2) {
    ath_error("[par_strat2d]: Domain decomposition in x2 (z) must be no more\
 than 2!\n");
  }
#endif

/* Initialize boxsize */
  x1min = pGrid->x1_0 + (pGrid->is + pGrid->idisp)*pGrid->dx1;
  x1max = pGrid->x1_0 + (pGrid->ie + pGrid->idisp + 1.0)*pGrid->dx1;
  Lx = x1max - x1min;

  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");
  Lz = x2max - x2min;

  Lg = nghost*pGrid->dx2; /* size of the ghost zone */

  js = pGrid->js;

/* Read initial conditions */
  Omega_0 = par_getd("problem","omega");
  qshear = par_getd_def("problem","qshear",1.5);
  ipert = par_geti_def("problem","ipert",1);
  vsc1 = par_getd_def("problem","vsc1",0.05); /* in unit of iso_sound (N.B.!) */
  vsc2 = par_getd_def("problem","vsc2",0.0);

  vsc1 = vsc1 * Iso_csound;
  vsc2 = vsc2 * Iso_csound;

  ScaleHg = Iso_csound/Omega_0;

  /* particle number */
  Npar  = (long)(par_geti("particle","parnumgrid"));

  pGrid->nparticle = Npar*pGrid->partypes;
  for (i=0; i<pGrid->partypes; i++)
    pGrid->grproperty[i].num = Npar;

  if (pGrid->nparticle+2 > pGrid->arrsize)
    particle_realloc(pGrid, pGrid->nparticle+2);

  ep = (Real*)calloc_1d_array(pGrid->partypes, sizeof(Real));
  ScaleHpar = (Real*)calloc_1d_array(pGrid->partypes, sizeof(Real));

  epsilon = (Real*)calloc_1d_array(pGrid->partypes, sizeof(Real));
  wxNSH   = (Real**)calloc_2d_array(pGrid->Nx2+1, pGrid->partypes,sizeof(Real));
  wyNSH   = (Real**)calloc_2d_array(pGrid->Nx2+1, pGrid->partypes,sizeof(Real));
  uxNSH   = (Real*)calloc_1d_array(pGrid->Nx2+1, sizeof(Real));
  uyNSH   = (Real*)calloc_1d_array(pGrid->Nx2+1, sizeof(Real));

  /* particle stopping time */
  tsmode = par_geti("particle","tsmode");
  if (tsmode == 3) {/* fixed stopping time */
    tsmin = par_getd("problem","tsmin"); /* in code unit */
    tsmax = par_getd("problem","tsmax");
    tscrit= par_getd("problem","tscrit");

    for (i=0; i<pGrid->partypes; i++) {
      tstop0[i] = tsmin*exp(i*log(tsmax/tsmin)/MAX(pGrid->partypes-1,1.0));
      pGrid->grproperty[i].rad = tstop0[i];
      /* use fully implicit integrator for well coupled particles */
      if (tstop0[i] < tscrit) pGrid->grproperty[i].integrator = 3;
    }
  }
  else { 
    amin = par_getd("problem","amin");
    amax = par_getd("problem","amax");

    for (i=0; i<pGrid->partypes; i++)
      pGrid->grproperty[i].rad = amin*
                           exp(i*log(amax/amin)/MAX(pGrid->partypes-1,1.0));

    if (tsmode >= 2) {/* Epstein/General regime */
      /* conversion factor for rhoa */
      rhoaconv = par_getd_def("problem","rhoaconv",1.0);
      for (i=0; i<pGrid->partypes; i++)
        grrhoa[i]=pGrid->grproperty[i].rad*rhoaconv;
    }

    if (tsmode == 3)  /* General drag formula */
      alamcoeff = par_getd("problem","alamcoeff");
  }

  /* particle scale height */
  Hparmax = par_getd("problem","hparmax"); /* in unit of gas scale height */
  Hparmin = par_getd("problem","hparmin");
  for (i=0; i<pGrid->partypes; i++) 
    ScaleHpar[i] = Hparmax*
                   exp(-i*log(Hparmax/Hparmin)/MAX(pGrid->partypes-1,1.0));

#ifdef FEEDBACK
  mratio = par_getd_def("problem","mratio",0.0); /* total mass fraction */
  pwind = par_getd_def("problem","pwind",0.0);   /* power law index */
  if (mratio < 0.0)
    ath_error("[par_strat2d]: mratio must be positive!\n");

  epsum = 0.0;
  for (i=0; i<pGrid->partypes; i++)
  {
    ep[i] = pow(pGrid->grproperty[i].rad,pwind);
    epsum += ep[i];
  }

  for (i=0; i<pGrid->partypes; i++)
  {
    ep[i] = mratio*ep[i]/epsum;
    pGrid->grproperty[i].m = sqrt(2.0*PI)*ScaleHg/Lz*ep[i]*
                                                  pGrid->Nx1*pGrid->Nx2/Npar;
  }
#else
  mratio = 0.0;
  for (i=0; i<pGrid->partypes; i++)
    ep[i] = 0.0;
#endif

  /* NSH equilibrium */
  for (j=pGrid->js; j<=pGrid->je+1; j++) {

    h = pGrid->x2_0 + (j+pGrid->jdisp)*pGrid->dx2;
    q = j - js;
    etavk = fabs(vsc1+vsc2*SQR(h));

    for (i=0; i<pGrid->partypes; i++) {
      epsilon[i] = ep[i]/ScaleHpar[i]/erf(Lz/(sqrt(8.0)*ScaleHpar[i]*ScaleHg))*
                   exp(-0.5*SQR(h/ScaleHg)*(SQR(1.0/ScaleHpar[i])-1.0));

      if (tsmode != 3)
        tstop0[i] = get_ts(pGrid,i,exp(-0.5*SQR(h/ScaleHg)),Iso_csound,etavk);
    }

    MultiNSH(pGrid->partypes, tstop0, epsilon, etavk,
                              &uxNSH[q], &uyNSH[q], wxNSH[q], wyNSH[q]);
  }

/* Now set initial conditions for the gas */
  for (j=pGrid->js; j<=pGrid->je; j++) {
  for (i=pGrid->is; i<=pGrid->ie; i++) {
    cc_pos(pGrid,i,j,pGrid->ks,&x1,&x2,&x3);

    rhog = exp(-0.5*SQR(x2/ScaleHg));
    pGrid->U[pGrid->ks][j][i].d = rhog;

    if (ipert != 1) {/* NSH velocity */
      pGrid->U[pGrid->ks][j][i].M1 = 0.5*rhog*(uxNSH[j-js]+uxNSH[j-js+1]);
      pGrid->U[pGrid->ks][j][i].M3 = 0.5*rhog*(uyNSH[j-js]+uyNSH[j-js+1]);
    } else {
      pGrid->U[pGrid->ks][j][i].M1 = 0.0;
      pGrid->U[pGrid->ks][j][i].M3 = 0.0;
    }

    pGrid->U[pGrid->ks][j][i].M2 = 0.0;
#ifndef FARGO
    pGrid->U[pGrid->ks][j][i].M3 -= qshear*rhog*Omega*x1;
#endif

  }}

/* Now set initial conditions for the particles */
  p = 0;
  dx2_1 = 1.0/pGrid->dx2;
  x3p = pGrid->x3_0 + (pGrid->ks+pGrid->kdisp)*pGrid->dx3;
  zmin = pGrid->x2_0 + (pGrid->js+pGrid->jdisp)*pGrid->dx2;
  zmax = pGrid->x2_0 + (pGrid->je+pGrid->jdisp+1.0)*pGrid->dx2;

  for (q=0; q<Npar; q++) {

    for (pt=0; pt<pGrid->partypes; pt++) {

      x1p = x1min + Lx*ran2(&iseed);
      x2p = ScaleHpar[pt]*ScaleHg*Normal(&iseed);
      while ((x2p >= zmax) || (x2p < zmin))
        x2p = ScaleHpar[pt]*ScaleHg*Normal(&iseed);

      pGrid->particle[p].property = pt;
      pGrid->particle[p].x1 = x1p;
      pGrid->particle[p].x2 = x2p;
      pGrid->particle[p].x3 = x3p;

      if (ipert != 1) {/* NSH velocity */

        cellj(pGrid, x2p, dx2_1, &j, &b);
        j = j-pGrid->js;  b = b - pGrid->js;

        pGrid->particle[p].v1 = (j+1-b)*wxNSH[j][pt]+(b-j)*wxNSH[j+1][pt];
        pGrid->particle[p].v3 = (j+1-b)*wyNSH[j][pt]+(b-j)*wyNSH[j+1][pt];

      } else {

        pGrid->particle[p].v1 = 0.0;
        pGrid->particle[p].v3 = vsc1+vsc2*SQR(x2p);

      }

      pGrid->particle[p].v2 = 0.0;
#ifndef FARGO
      pGrid->particle[p].v3 -= qshear*Omega*x1p;
#endif

      pGrid->particle[p].pos = 1; /* grid particle */
      pGrid->particle[p].my_id = p;
#ifdef MPI_PARALLEL
      pGrid->particle[p].init_id = pGrid->my_id;
#endif
      p++;
  }}

/* enroll gravitational potential function, shearing sheet BC functions */
  StaticGravPot = ShearingBoxPot;

  dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");

  /* set the # of the particles in list output
   * (by default, output 1 particle per cell)
   */
  nlis = par_geti_def("problem","nlis",pGrid->Nx1*pGrid->Nx2);

  /* set the number of particles to keep track of */
  ntrack = par_geti_def("problem","ntrack",2000);

  /* set the threshold particle density */
  dpar_thresh = par_geti_def("problem","dpar_thresh",10.0);

  /* set vertical boundary conditions (by default, periodic) */
  vertical_bc = par_geti_def("problem","vertical_bc",0);

  if (vertical_bc == 1) /* closed BC */
  {
    set_bvals_mhd_fun(left_x2,  close_ix2);
    set_bvals_mhd_fun(right_x2, close_ox2);
  }

  /* finalize */
  free(ep);  free(ScaleHpar);
  free(epsilon);
  free_2d_array(wxNSH);  free_2d_array(wyNSH);
  free(uxNSH);           free(uyNSH);

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
  return;
}

void problem_read_restart(Grid *pG, Domain *pD, FILE *fp)
{

  StaticGravPot = ShearingBoxPot;

  Omega_0 = par_getd("problem","omega");
  qshear = par_getd_def("problem","qshear",1.5);
  ipert = par_geti_def("problem","ipert",1);

  x1min = pG->x1_0 + (pG->is + pG->idisp)*pG->dx1;
  x1max = pG->x1_0 + (pG->ie + pG->idisp + 1.0)*pG->dx1;
  Lx = x1max - x1min;

  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");
  Lz = x2max - x2min;

  Lg = nghost*pG->dx2; /* size of the ghost zone */

  vsc1 = par_getd_def("problem","vsc1",0.05); /* in unit of iso_sound (N.B.!) */
  vsc2 = par_getd_def("problem","vsc2",0.0);

  vsc1 = vsc1 * Iso_csound;
  vsc2 = vsc2 * Iso_csound;

  Npar  = (int)(sqrt(par_geti("particle","parnumgrid")));
  nlis = par_geti_def("problem","nlis",pG->Nx1*pG->Nx2);
  ntrack = par_geti_def("problem","ntrack",2000);

  dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");

  return;
}

static Real hst_rho_Vx_dVy(const Grid *pG,
                           const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  return pG->U[k][j][i].M1*pG->U[k][j][i].M3/pG->U[k][j][i].d;
#else
  return pG->U[k][j][i].M1*(pG->U[k][j][i].M3/pG->U[k][j][i].d
                             + qshear*Omega_0*x1);
#endif
}

static void close_ix2(Grid *pGrid)
{
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  for (k=ks; k<=ke; k++)
    for (j=1; j<=nghost; j++)
      for (i=il; i<=iu; i++) {
        pGrid->U[k][js-j][i] = pGrid->U[k][js][i];
        pGrid->U[k][js-j][i].M2 = 0.0;
      }

  return;
}

static void close_ox2(Grid *pGrid)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  for (k=ks; k<=ke; k++)
    for (j=1; j<=nghost; j++)
      for (i=il; i<=iu; i++) {
        pGrid->U[k][je+j][i] = pGrid->U[k][je][i];
        pGrid->U[k][je+j][i].M2 = 0.0;
      }

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
  if (strcmp(name,"limit")==0)    return property_limit;
  if (strcmp(name,"trace")==0)    return property_trace;
  if (strcmp(name,"type")==0)    return property_type;

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
  Real z,fac;

  ft->x1 -= 2.0*(vsc1 + vsc2*SQR(x2))*Omega_0;

  if(x2 > x2max)
    z = x2-Lz;
  else if (x2 < x2min)
    z = x2+Lz;
  else
    z = x2;

  fac = Lg/(0.5*Lz+Lg-fabs(z));
  ft->x2 -= SQR(Omega_0)*z*(1.0-SQR(fac)*fac); /* 3rd order sharp */
//  *ft.x2 -= SQR(Omega_0)*z*(1.0-SQR(fac));  /* 2nd order sharp */

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
  return;
}
 
/*=========================== PRIVATE FUNCTIONS ==============================*/
/*--------------------------------------------------------------------------- */
/*! \fn static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3)
 *  \brief Shearing box tidal gravitational potential */
static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3)
{
  Real z,z0,phi=0.0;

#ifndef FARGO
  phi -= qshear*SQR(Omega_0*x1);
#endif
  /* Ensure vertical periodicity in ghost zones */
  if(x2 > x2max)
    z=fabs(x2-Lz);
  else if (x2 < x2min)
    z=fabs(x2+Lz);
  else
    z=fabs(x2);

  phi += 0.5*SQR(Omega_0*z);

  /* smooth the potential at domain edges */
  z0 = 0.5*Lz+Lg;
  phi -= SQR(Omega_0*Lg)*Lg*(2.0*z-z0)/(2.0*SQR(z0-z)); /* 3rd order sharp */

//  phi -= SQR(Omega_0*Lg)*(z/(z0-z)+log(z0/(z0-z))); /* 2nd order sharp */

  return phi;
}

/*! \fn static int property_limit(const Grain *gr, const GrainAux *grsub)
 *  \brief User defined particle selection function (1: true; 0: false) */
static int property_limit(const Grain *gr, const GrainAux *grsub)
{
  if ((gr->pos == 1) && (gr->my_id<nlis))
    return 1;
  else
    return 0;
}

/*! \fn static int property_trace(const Grain *gr, const GrainAux *grsub)
 *  \brief User defined particle selection function (1: true; 0: false) */
static int property_trace(const Grain *gr, const GrainAux *grsub)
{
  if ((gr->pos == 1) && (gr->my_id<ntrack))
    return 1;
  else
    return 0;
}

/*! \fn static int property_type(const Grain *gr, const GrainAux *grsub)
 *  \brief User defined particle selection function (1: true; 0: false) */
static int property_type(const Grain *gr, const GrainAux *grsub)
{
  if ((gr->pos == 1) && (gr->property == mytype))
    return 1;
  else
    return 0;
}

/*! \fn static int property_dense(const Grain *gr, const GrainAux *grsub)
 *  \brief User defined particle selection function (1: true; 0: false) */
static int property_dense(const Grain *gr, const GrainAux *grsub)
{
  if ((gr->pos == 1) && (grsub->dpar > dpar_thresh))
    return 1;
  else
    return 0;
}

/*----------------------------------------------------------------------------*/
/*! \fn void MultiNSH(int n, Real *tstop, Real *mratio, Real etavk,
 *                   Real *uxNSH, Real *uyNSH, Real *wxNSH, Real *wyNSH)
 *  \brief Multi-species NSH equilibrium
 *
 * Input: # of particle types (n), dust stopping time and mass ratio array, and
 *        drift speed etavk.
 * Output: gas NSH equlibrium velocity (u), and dust NSH equilibrium velocity
 *         array (w).
 */
void MultiNSH(int n, Real *tstop, Real *mratio, Real etavk,
                     Real *uxNSH, Real *uyNSH, Real *wxNSH, Real *wyNSH)
{
  int i,j;
  Real *Lambda1,**Lam1GamP1, **A, **B, **Tmp;

  Lambda1 = (Real*)calloc_1d_array(n, sizeof(Real));     /* Lambda^{-1} */
  Lam1GamP1=(Real**)calloc_2d_array(n, n, sizeof(Real)); /* Lambda1*(1+Gamma) */
  A       = (Real**)calloc_2d_array(n, n, sizeof(Real));
  B       = (Real**)calloc_2d_array(n, n, sizeof(Real));
  Tmp     = (Real**)calloc_2d_array(n, n, sizeof(Real));

  /* definitions */
  for (i=0; i<n; i++){
    for (j=0; j<n; j++)
      Lam1GamP1[i][j] = mratio[j];
    Lam1GamP1[i][i] += 1.0;
    Lambda1[i] = 1.0/(tstop[i]+1.0e-16);
    for (j=0; j<n; j++)
      Lam1GamP1[i][j] *= Lambda1[i];
  }

  /* Calculate A and B */
  MatrixMult(Lam1GamP1, Lam1GamP1, n,n,n, Tmp);
  for (i=0; i<n; i++) Tmp[i][i] += 1.0;
  InverseMatrix(Tmp, n, B);
  for (i=0; i<n; i++)
  for (j=0; j<n; j++)
    B[i][j] *= Lambda1[j];
  MatrixMult(Lam1GamP1, B, n,n,n, A);

  /* Obtain NSH velocities */
  *uxNSH = 0.0;  *uyNSH = 0.0;
  for (i=0; i<n; i++){
    wxNSH[i] = 0.0;
    wyNSH[i] = 0.0;
    for (j=0; j<n; j++){
      wxNSH[i] -= B[i][j];
      wyNSH[i] -= A[i][j];
    }
    wxNSH[i] *= 2.0*etavk;
    wyNSH[i] *= etavk;
    *uxNSH -= mratio[i]*wxNSH[i];
    *uyNSH -= mratio[i]*wyNSH[i];
    wyNSH[i] += etavk;
  }

  free(Lambda1);
  free_2d_array(A);         free_2d_array(B);
  free_2d_array(Lam1GamP1); free_2d_array(Tmp);

  return;
}

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

/*--------------------------------------------------------------------------- */
/*! \fn double Normal(long int *idum)
 *  \brief Normal distribution random number generator */
double Normal(long int *idum)
{
  double Y,X1,X2;

  X1 = ran2(idum);
  X2 = ran2(idum);

  Y = sqrt(-2.0*log(X1+TINY_NUMBER))*cos(2*PI*X2);

  return Y;
}

/*--------------------------------------------------------------------------- */
/*! \fn Real Erf(Real z)
 *  \brief Error function  */
Real Erf(Real z)
{
  /* arrays of the error function y=erf(x) */
  static double x[101]={
        0.000000e+000,  3.783387e-003,  7.709914e-003,  1.178500e-002,  1.601425e-002,  2.040352e-002,
        2.495885e-002,  2.968653e-002,  3.459307e-002,  3.968525e-002,  4.497008e-002,  5.045486e-002,
        5.614715e-002,  6.205480e-002,  6.818596e-002,  7.454909e-002,  8.115295e-002,  8.800667e-002,
        9.511969e-002,  1.025018e-001,  1.101632e-001,  1.181145e-001,  1.263667e-001,  1.349310e-001,
        1.438193e-001,  1.530440e-001,  1.626176e-001,  1.725534e-001,  1.828652e-001,  1.935671e-001,
        2.046738e-001,  2.162008e-001,  2.281639e-001,  2.405796e-001,  2.534651e-001,  2.668380e-001,
        2.807169e-001,  2.951209e-001,  3.100699e-001,  3.255844e-001,  3.416859e-001,  3.583966e-001,
        3.757395e-001,  3.937386e-001,  4.124186e-001,  4.318054e-001,  4.519256e-001,  4.728071e-001,
        4.944786e-001,  5.169701e-001,  5.403124e-001,  5.645379e-001,  5.896800e-001,  6.157732e-001,
        6.428537e-001,  6.709587e-001,  7.001271e-001,  7.303990e-001,  7.618162e-001,  7.944220e-001,
        8.282614e-001,  8.633812e-001,  8.998296e-001,  9.376570e-001,  9.769156e-001,  1.017659e+000,
        1.059945e+000,  1.103830e+000,  1.149376e+000,  1.196644e+000,  1.245701e+000,  1.296614e+000,
        1.349454e+000,  1.404292e+000,  1.461205e+000,  1.520272e+000,  1.581573e+000,  1.645193e+000,
        1.711221e+000,  1.779746e+000,  1.850864e+000,  1.924673e+000,  2.001274e+000,  2.080774e+000,
        2.163281e+000,  2.248909e+000,  2.337778e+000,  2.430008e+000,  2.525728e+000,  2.625070e+000,
        2.728170e+000,  2.835170e+000,  2.946219e+000,  3.061469e+000,  3.181080e+000,  3.305216e+000,
        3.434048e+000,  3.567755e+000,  3.706521e+000,  3.850536e+000,  4.000000e+000
  };
  static double y[101]={
        0.00000000e+000,        4.26907434e-003,        8.69953340e-003,        1.32973284e-002,        1.80686067e-002,        2.30197153e-002,
        2.81572033e-002,        3.34878242e-002,        3.90185379e-002,        4.47565113e-002,        5.07091186e-002,        5.68839404e-002,
        6.32887618e-002,        6.99315688e-002,        7.68205444e-002,        8.39640613e-002,        9.13706742e-002,        9.90491090e-002,
        1.07008250e-001,        1.15257124e-001,        1.23804883e-001,        1.32660778e-001,        1.41834139e-001,        1.51334337e-001,
        1.61170754e-001,        1.71352743e-001,        1.81889576e-001,        1.92790394e-001,        2.04064148e-001,        2.15719527e-001,
        2.27764884e-001,        2.40208149e-001,        2.53056730e-001,        2.66317410e-001,        2.79996226e-001,        2.94098338e-001,
        3.08627885e-001,        3.23587825e-001,        3.38979770e-001,        3.54803790e-001,        3.71058224e-001,        3.87739454e-001,
        4.04841688e-001,        4.22356710e-001,        4.40273635e-001,        4.58578645e-001,        4.77254725e-001,        4.96281391e-001,
        5.15634428e-001,        5.35285634e-001,        5.55202571e-001,        5.75348359e-001,        5.95681482e-001,        6.16155658e-001,
        6.36719759e-001,        6.57317799e-001,        6.77889021e-001,        6.98368078e-001,        7.18685336e-001,        7.38767318e-001,
        7.58537287e-001,        7.77916009e-001,        7.96822665e-001,        8.15175962e-001,        8.32895397e-001,        8.49902691e-001,
        8.66123358e-001,        8.81488386e-001,        8.95935967e-001,        9.09413237e-001,        9.21877939e-001,        9.33299942e-001,
        9.43662512e-001,        9.52963249e-001,        9.61214608e-001,        9.68443923e-001,        9.74692870e-001,        9.80016358e-001,
        9.84480847e-001,        9.88162149e-001,        9.91142807e-001,        9.93509200e-001,        9.95348535e-001,        9.96745927e-001,
        9.97781755e-001,        9.98529475e-001,        9.99054014e-001,        9.99410828e-001,        9.99645625e-001,        9.99794704e-001,
        9.99885782e-001,        9.99939162e-001,        9.99969080e-001,        9.99985060e-001,        9.99993164e-001,        9.99997050e-001,
        9.99998805e-001,        9.99999548e-001,        9.99999841e-001,        9.99999948e-001,        9.99999985e-001
  };
  double z1, g;
  int il, iu, im;
  z1 = fabs(z);
  /* look for the location of z1 in the x array */
  il=0; iu=100;
  while(iu-il>1) {
    im = (iu+il)/2;
    if (x[im] > z1) iu = im;
    else il = im;
  }
  /* linear interpolation */
  g = (y[iu]*(z1-x[il])+y[il]*(x[iu]-z1))/(x[iu]-x[il]);

  g = MIN(g,1.0);
  if (z<0.0) g = -g;
  return g;
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
