#include "copyright.h"
/*============================================================================*/
/*! \file turb.c
 *  \brief Problem generator for driven and decaying turbulence.
 *
 * PURPOSE: Problem generator for driven and decaying turbulence. Only works
 *   in 3D with periodic BC.  Arbitrary power spectrum specified using ispect:
 *  -  ispect=1: power law - original form
 *  -  ispect=2: form from Gammie&Ostriker
 *  Driving specified using idrive
 *  -  idrive=0: driven turbulence (de=dedt*dt before each time step,
 *                                     unless IMPULSIVE_DRIVING enabled)
 *  -  idrive=1: decaying turbulence (de=dedt before first time step ONLY;
 *                                     this mode uses less memory!)
 *
 *  HISTORY:
 *  - Original ZEUS version (gmc.f) written by J. Stone, 24 Jan 1996
 *  - First Athena version written by J. Stone, 9 June 2004
 *  - Major rewrite to add MPI and use FFTW by Nicole Lemaster, 28 Sept 2006
 *
 *  - Last updated May 11, 2007
 *
 *  REFERENCE: "Dissipation in Compressible MHD Turbulence", by J. Stone,
 *    E. Ostriker, & C. Gammie, ApJ 508, L99 (1998)			      */
/*============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"
#include "globals.h"

#ifndef ISOTHERMAL
#error Problem generator only works for isothermal turbulence
#endif /* ISOTHERMAL */

#ifdef MPI_PARALLEL
#include "mpi.h"
#ifdef DOUBLE_PREC
#define MPI_RL MPI_DOUBLE
#else /* DOUBLE_PREC */
#define MPI_RL MPI_FLOAT
#endif /* DOUBLE_PREC */
#endif /* MPI_PARALLEL */

/* Uncomment the following define to drive the flow in an impulsive manner
   as was done originally.  Restarts for this mode not yet implemented! */
/* #define IMPULSIVE_DRIVING */

/* KEEP SEMI-COLONS OUT OF THESE PRE-PROCESSOR DIRECTIVES! */
/* FFT indexing Nfast=k, Nmid=j, Nslow=i (opposite to Athena)
 * For OFST, i,j,k,nx2,nx3 reference the local grid */
#define OFST(i, j, k) ((k) + nx3*((j) + nx2*(i)))
/* KWVM: magnitude of wavenumber k in units of dkx */
#define KWVM(i, j, k) (sqrt(SQR(KCOMP(i,gis,gnx1))+ \
                            SQR(KCOMP(j,gjs,gnx2))+SQR(KCOMP(k,gks,gnx3))))

/* FFTW - Variables, Plan, etc. */
/* These are made static global variables so that they need not be
   allocated AND destroyed with each call to pspect! */
static struct ath_3d_fft_plan *plan;
/* Between calls to generate(), these have unshifted, unnormalized
 * velocity perturbations. */
static ath_fft_data *fv1=NULL, *fv2=NULL, *fv3=NULL;

/* Normalized, shifted velocity perturbations */
static Real ***dv1=NULL, ***dv2=NULL, ***dv3=NULL;
/* Cutoff wavenumbers, G&O spect peak, power law spect exponent, 2 pi/L */
static Real klow,khigh,kpeak,expo,dkx;
/* Energy injection rate, last planned driving time, driving interval.
 * If not using impulsive driving, then the time quantities above are for 
 * computing a new spectrum, not driving */
static Real dedt,tdrive,dtdrive;
/* Driving properties */
static int ispect,idrive;
/* Number of cells in local grid, number of cells in global grid */
static int nx1,nx2,nx3,gnx1,gnx2,gnx3;
/* Starting and ending indices for global grid */
static int gis,gie,gjs,gje,gks,gke;
/* Seed for random number generator */
long int rseed;
#ifdef MHD
/* beta = isothermal pressure / magnetic pressure
 * B0 = sqrt(2.0*Iso_csound2*rhobar/beta) is init magnetic field strength */
static Real beta,B0;
#endif /* MHD */
/* Initial density (will be average density throughout simulation) */
static const Real rhobar = 1.0;

/* Functions appear in this file in the same order that they appear in the
 * prototypes below */

/* Function prototypes for generating velocity perturbations */
static void pspect(ath_fft_data *ampl);
static void project();
static inline void transform();
static inline void generate();
static void perturb(Grid *pGrid, Real dt);

/* Function prototypes for initializing and interfacing with Athena */
static void initialize(Grid *pGrid, Domain *pD);
/* void problem(Grid *pGrid, Domain *pD); */
/* void Userwork_in_loop(Grid *pGrid, Domain *pD); */
/* void Userwork_after_loop(Grid *pGrid, Domain *pD); */
/* void problem_write_restart(Grid *pG, Domain *pD, FILE *fp); */
/* void problem_read_restart(Grid *pG, Domain *pD, FILE *fp); */
/* Gasfun_t get_usr_expr(const char *expr); */

/* Function prototypes for analysis and outputs */
static Real hst_dEk(const Grid *pG, const int i, const int j, const int k);
static Real hst_dEb(const Grid *pG, const int i, const int j, const int k);

/* Function prototypes for Numerical Recipes functions */
static double ran2(long int *idum);

/* ========================================================================== */

/*! \fn static void pspect(ath_fft_data *ampl)
 *  \brief computes component of velocity with specific power
 *  spectrum in Fourier space determined by ispect
 *
 *  Velocity power spectrum returned in ampl
 *  - klow   = multiple of 2 pi/L for cut-off at low  wavenumbers
 *  - khigh  = multiple of 2 pi/L for cut-off at high wavenumbers
 *  - expo   = exponent of power law
 *  - ispect = integer flag which specifies spectrum
 *
 *  Note that the fourier amplitudes are stored in an array with no
 *  ghost zones
 */
static void pspect(ath_fft_data *ampl)
{
  int i,j,k;
  double q1,q2,q3;

  /* set random amplitudes with gaussian deviation */
  for (k=0; k<nx3; k++) {
    for (j=0; j<nx2; j++) {
      for (i=0; i<nx1; i++) {
        q1 = ran2(&rseed);
        q2 = ran2(&rseed);
        q3 = sqrt(-2.0*log(q1+1.0e-20))*cos(2.0*PI*q2);
        q1 = ran2(&rseed);
        ampl[OFST(i,j,k)][0] = q3*cos(2.0*PI*q1);
        ampl[OFST(i,j,k)][1] = q3*sin(2.0*PI*q1);
      }
    }
  }

  /* set power spectrum
   *   ispect=1: power law - original form
   *   ispect=2: form from Gammie&Ostriker
   */
  for (k=0; k<nx3; k++) {
    for (j=0; j<nx2; j++) {
      for (i=0; i<nx1; i++) {
        /* compute k/dkx */
        q3 = KWVM(i,j,k);
        if ((q3 > klow) && (q3 < khigh)) {
          q3 *= dkx; /* multiply by 2 pi/L */
          if (ispect == 1) {
            /* decreasing power law */
            ampl[OFST(i,j,k)][0] /= pow(q3,(expo+2.0)/2.0);
            ampl[OFST(i,j,k)][1] /= pow(q3,(expo+2.0)/2.0);
          } else if (ispect == 2) {
            /* G&O form */
            ampl[OFST(i,j,k)][0] *= pow(q3,3.0)*exp(-4.0*q3/kpeak);
            ampl[OFST(i,j,k)][1] *= pow(q3,3.0)*exp(-4.0*q3/kpeak);
          }
        } else {
          /* introduce cut-offs at klow and khigh */
          ampl[OFST(i,j,k)][0] = 0.0;
          ampl[OFST(i,j,k)][1] = 0.0;
        }
      }
    }
  }
  ampl[0][0] = 0.0;
  ampl[0][1] = 0.0;

  return;
}

/* ========================================================================== */

/*! \fn static void project()
 *  \brief Makes velocity perturbations divergence free
 */
static void project()
{
  int i,j,k,m,ind;
  double kap[3], kapn[3], mag;
  ath_fft_data dot;
  
  /* Project off non-solenoidal component of velocity */
  for (k=0; k<nx3; k++) {
    kap[2] = sin(2.0*PI*(gks+k)/gnx3);
    for (j=0; j<nx2; j++) {
      kap[1] = sin(2.0*PI*(gjs+j)/gnx2);
      for (i=0; i<nx1; i++) {
        if (((gis+i)+(gjs+j)+(gks+k)) != 0) {
          kap[0] = sin(2.0*PI*(gis+i)/gnx1);
          ind = OFST(i,j,k);

          /* make kapn a unit vector */
          mag = sqrt(SQR(kap[0]) + SQR(kap[1]) + SQR(kap[2]));
          for (m=0; m<3; m++) kapn[m] = kap[m] / mag;

          /* find fv_0 dot kapn */
          dot[0] = fv1[ind][0]*kapn[0]+fv2[ind][0]*kapn[1]+fv3[ind][0]*kapn[2];
          dot[1] = fv1[ind][1]*kapn[0]+fv2[ind][1]*kapn[1]+fv3[ind][1]*kapn[2];

          /* fv = fv_0 - (fv_0 dot kapn) * kapn */
          fv1[ind][0] -= dot[0]*kapn[0];
          fv2[ind][0] -= dot[0]*kapn[1];
          fv3[ind][0] -= dot[0]*kapn[2];

          fv1[ind][1] -= dot[1]*kapn[0];
          fv2[ind][1] -= dot[1]*kapn[1];
          fv3[ind][1] -= dot[1]*kapn[2];
        }
      }
    }
  }

  return;
}

/* ========================================================================== */

/*! \fn static inline void transform()
 *  \brief Generate velocities from fourier transform
 */
static inline void transform()
{
  /* Transform velocities from k space to physical space */
  ath_3d_fft(plan, fv1);
  ath_3d_fft(plan, fv2);
  ath_3d_fft(plan, fv3);

  /* Should technically renormalize (divide by gnx1*gnx2*gnx3) here, but
   * since we're going to renormalize to get the desired energy injection
   * rate anyway, there's no point */
 
  return;
}

/* ========================================================================== */

/*! \fn static inline void generate()
 *  \brief Generate the velocity perturbations
 */
static inline void generate()
{
  /* Generate new perturbations following appropriate power spectrum */
  pspect(fv1);
  pspect(fv2);
  pspect(fv3);

  /* Require div V = 0 */
  project();

  /* Transform perturbations to real space, but don't normalize until
   * just before we apply them in perturb() */
  transform();

  return;
}

/* ========================================================================== */

/*! \fn static void perturb(Grid *pGrid, Real dt)
 *  \brief  Shifts velocities so no net momentum change, normalizes to keep
 *  dedt fixed, and then sets velocities
 */
static void perturb(Grid *pGrid, Real dt)
{
  int i, is=pGrid->is, ie = pGrid->ie;
  int j, js=pGrid->js, je = pGrid->je;
  int k, ks=pGrid->ks, ke = pGrid->ke;
  int ind, mpierr;
  Real dvol, aa, b, c, s, de, qa, v1, v2, v3;
  Real t0, t0ij, t0i, t1, t1ij, t1i;
  Real t2, t2ij, t2i, t3, t3ij, t3i;
  Real m[4], gm[4];

  /* Set the velocities in real space */
  dvol = 1.0/((Real)(gnx1*gnx2*gnx3));
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        ind = OFST(i-is,j-js,k-ks);
        dv1[k][j][i] = fv1[ind][0]*dvol;
        dv2[k][j][i] = fv2[ind][0]*dvol;
        dv3[k][j][i] = fv3[ind][0]*dvol;
      }
    }
  }

  /* Calculate net momentum pertubation components t1, t2, t3 */
  t0 = 0.0;  t1 = 0.0;  t2 = 0.0;  t3 = 0.0;
  for (k=ks; k<=ke; k++) {
    t0ij = 0.0;  t1ij = 0.0;  t2ij = 0.0;  t3ij = 0.0;
    for (j=js; j<=je; j++) {
      t0i = 0.0;  t1i = 0.0;  t2i = 0.0;  t3i = 0.0;
      for (i=is; i<=ie; i++) {
        t0i += pGrid->U[k][j][i].d;

	/* The net momentum perturbation */
        t1i += pGrid->U[k][j][i].d * dv1[k][j][i];
        t2i += pGrid->U[k][j][i].d * dv2[k][j][i];
        t3i += pGrid->U[k][j][i].d * dv3[k][j][i];
      }
      t0ij += t0i;  t1ij += t1i;  t2ij += t2i;  t3ij += t3i;
    }
    t0 += t0ij;  t1 += t1ij;  t2 += t2ij;  t3 += t3ij;
  }

#ifdef MPI_PARALLEL
  /* Sum the perturbations over all processors */
  m[0] = t0;  m[1] = t1;  m[2] = t2;  m[3] = t3;
  mpierr = MPI_Allreduce(m, gm, 4, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) ath_error("[normalize]: MPI_Allreduce error = %d\n", mpierr);
  t0 = gm[0];  t1 = gm[1];  t2 = gm[2];  t3 = gm[3];
#endif /* MPI_PARALLEL */

  /* Subtract the mean velocity perturbation so that the net momentum
   * perturbation is zero. */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        dv1[k][j][i] -= t1/t0;
        dv2[k][j][i] -= t2/t0;
        dv3[k][j][i] -= t3/t0;
      }
    }
  }

  /* Calculate unscaled energy of perturbations */
  t1 = 0.0;  t2 = 0.0;
  for (k=ks; k<=ke; k++) {
    t1ij = 0.0;  t2ij = 0.0;
    for (j=js; j<=je; j++) {
      t1i = 0.0;  t2i = 0.0;
      for (i=is; i<=ie; i++) {
        /* Calculate velocity pertubation at cell center from
         * perturbations at cell faces */
	v1 = dv1[k][j][i];
	v2 = dv2[k][j][i];
	v3 = dv3[k][j][i];

        t1i += (pGrid->U[k][j][i].d)*(SQR(v1) + SQR(v2) + SQR(v3));
	t2i +=  (pGrid->U[k][j][i].M1)*v1 + (pGrid->U[k][j][i].M2)*v2 +
                     (pGrid->U[k][j][i].M3)*v3;
      }
      t1ij += t1i;  t2ij += t2i;
    }
    t1 += t1ij;  t2 += t2ij;
  }

#ifdef MPI_PARALLEL
  /* Sum the perturbations over all processors */
  m[0] = t1;  m[1] = t2;
  mpierr = MPI_Allreduce(m, gm, 2, MPI_RL, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) ath_error("[normalize]: MPI_Allreduce error = %d\n", mpierr);
  t1 = gm[0];  t2 = gm[1];
#endif /* MPI_PARALLEL */

  /* Rescale to give the correct energy injection rate */
  dvol = pGrid->dx1*pGrid->dx2*pGrid->dx3;
  if (idrive == 0) {
    /* driven turbulence */
    de = dedt*dt;
  } else {
    /* decaying turbulence (all in one shot) */
    de = dedt;
  }
  aa = 0.5*t1;
  aa = MAX(aa,1.0e-20);
  b = t2;
  c = -de/dvol;
  if(b >= 0.0)
    s = (-2.0*c)/(b + sqrt(b*b - 4.0*aa*c));
  else
    s = (-b + sqrt(b*b - 4.0*aa*c))/(2.0*aa);

  if (isnan(s)) ath_error("[perturb]: s is NaN!\n");

  /* Apply momentum pertubations */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        qa = s*pGrid->U[k][j][i].d;
        pGrid->U[k][j][i].M1 += qa*dv1[k][j][i];
        pGrid->U[k][j][i].M2 += qa*dv2[k][j][i];
        pGrid->U[k][j][i].M3 += qa*dv3[k][j][i];
      }
    }
  }

  return;
}

/* ========================================================================== */
/*! \fn static void initialize(Grid *pGrid, Domain *pD)
 *  \brief  Allocate memory and initialize FFT plans */
static void initialize(Grid *pGrid, Domain *pD)
{
  int i, is=pGrid->is, ie = pGrid->ie;
  int j, js=pGrid->js, je = pGrid->je;
  int k, ks=pGrid->ks, ke = pGrid->ke;
  int nbuf, mpierr, nx1gh, nx2gh, nx3gh;
  float kwv, kpara, kperp;
  char donedrive = 0;

/* -----------------------------------------------------------
 * Variables within this block are stored globally, and used
 * within preprocessor macros.  Don't create variables with
 * these names within your function if you are going to use
 * OFST(), KCOMP(), or KWVM() within the function! */

  /* Get local grid size */
  nx1 = (ie-is+1);
  nx2 = (je-js+1);
  nx3 = (ke-ks+1);

  /* Get global grid size */
  gnx1 = pD->ide - pD->ids + 1;
  gnx2 = pD->jde - pD->jds + 1;
  gnx3 = pD->kde - pD->kds + 1;

  /* Get extents of local FFT grid in global coordinates */
  gis=is+pGrid->idisp;  gie=ie+pGrid->idisp;
  gjs=js+pGrid->jdisp;  gje=je+pGrid->jdisp;
  gks=ks+pGrid->kdisp;  gke=ke+pGrid->kdisp;
/* ----------------------------------------------------------- */

  /* Get size of arrays with ghost cells */
  nx1gh = nx1 + 2*nghost;
  nx2gh = nx2 + 2*nghost;
  nx3gh = nx3 + 2*nghost;

  /* Get input parameters */

  /* interval for generating new driving spectrum; also interval for
   * driving when IMPULSIVE_DRIVING is used */
  dtdrive = par_getd("problem","dtdrive");
#ifdef MHD
  /* magnetic field strength */
  beta = par_getd("problem","beta");
  /* beta = isothermal pressure/magnetic pressure */
  B0 = sqrt(2.0*Iso_csound2*rhobar/beta);
#endif /* MHD */
  /* energy injection rate */
  dedt = par_getd("problem","dedt");

  /* parameters for spectrum */
  ispect = par_geti("problem","ispect");
  if (ispect == 1) {
    expo = par_getd("problem","expo");
  } else if (ispect == 2) {
    kpeak = par_getd("problem","kpeak")*2.0*PI;
  } else {
    ath_error("Invalid value for ispect\n");
  }
  /* Cutoff wavenumbers of spectrum */
  klow = par_getd("problem","klow"); /* in integer units */
  khigh = par_getd("problem","khigh"); /* in integer units */
  dkx = 2.0*PI/(pGrid->dx1*gnx1); /* convert k from integer */

  /* Driven or decaying */
  idrive = par_geti("problem","idrive");
  if ((idrive < 0) || (idrive > 1)) ath_error("Invalid value for idrive\n");
  /* If restarting with decaying turbulence, no driving necessary. */
  if ((idrive == 1) && (pGrid->nstep > 0)) {
    donedrive = 1;
  }

  if (donedrive == 0) {
    /* Allocate memory for components of velocity perturbation */
    if ((dv1=(Real***)calloc_3d_array(nx3gh,nx2gh,nx1gh,sizeof(Real)))==NULL) {
      ath_error("[problem]: Error allocating memory for vel pert\n");
    }
    if ((dv2=(Real***)calloc_3d_array(nx3gh,nx2gh,nx1gh,sizeof(Real)))==NULL) {
      ath_error("[problem]: Error allocating memory for vel pert\n");
    }
    if ((dv3=(Real***)calloc_3d_array(nx3gh,nx2gh,nx1gh,sizeof(Real)))==NULL) {
      ath_error("[problem]: Error allocating memory for vel pert\n");
    }
  }

  /* Initialize the FFT plan */
  plan = ath_3d_fft_quick_plan(pGrid, pD, NULL, ATH_FFT_BACKWARD);

  /* Allocate memory for FFTs */
  if (donedrive == 0) {
    fv1 = ath_3d_fft_malloc(plan);
    fv2 = ath_3d_fft_malloc(plan);
    fv3 = ath_3d_fft_malloc(plan);
  }

  /* Enroll outputs */
  dump_history_enroll(hst_dEk,"<dE_K>");
  dump_history_enroll(hst_dEb,"<dE_B>");

  return;
}

/* ========================================================================== */

/*
 *  Function problem
 *
 *  Set up initial conditions, allocate memory, and initialize FFT plans
 */

void problem(Grid *pGrid, Domain *pD)
{
  int i, is=pGrid->is, ie = pGrid->ie;
  int j, js=pGrid->js, je = pGrid->je;
  int k, ks=pGrid->ks, ke = pGrid->ke;

  rseed = (pGrid->my_id+1);
  initialize(pGrid, pD);
  tdrive = 0.0;

  /* Initialize uniform density and momenta */
  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][j][i].d = rhobar;
        pGrid->U[k][j][i].M1 = 0.0;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;
      }
    }
  }

#ifdef MHD
  /* Initialize uniform magnetic field */
  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][j][i].B1c  = B0;
        pGrid->U[k][j][i].B2c  = 0.0;
        pGrid->U[k][j][i].B3c  = 0.0;
        pGrid->B1i[k][j][i] = B0;
        pGrid->B2i[k][j][i] = 0.0;
        pGrid->B3i[k][j][i] = 0.0;
      }
    }
  }
#endif /* MHD */

  /* Set the initial perturbations.  Note that we're putting in too much
   * energy this time.  This is okay since we're only interested in the
   * saturated state. */
  generate();
  perturb(pGrid, dtdrive);

  /* If decaying turbulence, no longer need the driving memory */
  if (idrive == 1) {
    ath_pout(0,"De-allocating driving memory.\n");

    /* Free Athena-style arrays */
    free_3d_array(dv1);
    free_3d_array(dv2);
    free_3d_array(dv3);

    /* Free FFTW-style arrays */
    ath_3d_fft_free(fv1);
    ath_3d_fft_free(fv2);
    ath_3d_fft_free(fv3);
  }

  return;
}

/* ========================================================================== */

/*
 *  Function Userwork_in_loop
 *
 *  Drive velocity field for turbulence in GMC problems
 */

void Userwork_in_loop(Grid *pGrid, Domain *pD)
{
  int i, is=pGrid->is, ie = pGrid->ie;
  int j, js=pGrid->js, je = pGrid->je;
  int k, ks=pGrid->ks, ke = pGrid->ke;
  Real newtime;

  if (isnan(pGrid->dt)) ath_error("Time step is NaN!");

  if (idrive == 0) {  /* driven turbulence */
    /* Integration has already been done, but time not yet updated */
    newtime = pGrid->time + pGrid->dt;

#ifndef IMPULSIVE_DRIVING
    /* Drive on every time step */
    perturb(pGrid, pGrid->dt);
#endif /* IMPULSIVE_DRIVING */

    if (newtime >= (tdrive+dtdrive)) {
      /* If we start with large time steps so that tdrive would get way
       * behind newtime, this makes sure we don't keep generating after
       * dropping down to smaller time steps */
      while ((tdrive+dtdrive) <= newtime) tdrive += dtdrive;

#ifdef IMPULSIVE_DRIVING
      /* Only drive at intervals of dtdrive */
      perturb(pGrid, dtdrive);
#endif /* IMPULSIVE_DRIVING */

      /* Compute new spectrum after dtdrive.  Putting this after perturb()
       * means we won't be applying perturbations from a new power spectrum
       * just before writing outputs.  At the very beginning, we'll go a
       * little longer before regenerating, but the energy injection rate
       * was off on the very first timestep anyway.  When studying driven
       * turbulence, all we care about is the saturated state. */
      generate();
    }
  }

  return;
}

/* ========================================================================== */

void Userwork_after_loop(Grid *pGrid, Domain *pD)
{
  /* Don't free memory here if doing any analysis because final
   * output hasn't been written yet!! */
  return;
}

void problem_write_restart(Grid *pG, Domain *pD, FILE *fp)
{  return;  }

void problem_read_restart(Grid *pG, Domain *pD, FILE *fp)
{  
  /* Allocate memory and initialize everything */
  rseed  = (pG->my_id+1);
  initialize(pG, pD);
  tdrive = pG->time;

  /* Generate a new power spectrum */
  if (idrive == 0) generate();

  return;
}

Gasfun_t get_usr_expr(const char *expr)
{  return NULL;  }

VGFunout_t get_usr_out_fun(const char *name){
  return NULL;
}

/* ========================================================================== */

/*
 *  Function hst_*
 *
 *  Dumps to history file
 */

/*! \fn static Real hst_dEk(const Grid *pG, const int i,const int j,const int k)
 *  \brief Dump kinetic energy in perturbations */
static Real hst_dEk(const Grid *pG, const int i, const int j, const int k)
{ /* The kinetic energy in perturbations is 0.5*d*V^2 */
  return 0.5*(pG->U[k][j][i].M1*pG->U[k][j][i].M1 +
	      pG->U[k][j][i].M2*pG->U[k][j][i].M2 +
	      pG->U[k][j][i].M3*pG->U[k][j][i].M3)/pG->U[k][j][i].d;
}

/*! \fn static Real hst_dEb(const Grid *pG, const int i,const int j,const int k)
 *  \brief Dump magnetic energy in perturbations */
static Real hst_dEb(const Grid *pG, const int i, const int j, const int k)
{ /* The magnetic energy in perturbations is 0.5*B^2 - 0.5*B0^2 */
#ifdef MHD
  return 0.5*((pG->U[k][j][i].B1c*pG->U[k][j][i].B1c +
	       pG->U[k][j][i].B2c*pG->U[k][j][i].B2c +
	       pG->U[k][j][i].B3c*pG->U[k][j][i].B3c)-B0*B0);
#else /* MHD */
  return 0.0;
#endif /* MHD */
}

/* ========================================================================== */

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
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)
#define NTAB 32

/*! \fn double ran2(long int *idum){
 *  \brief The routine ran2() is extracted from the Numerical Recipes in C 
 *
 * The routine ran2() is extracted from the Numerical Recipes in C
 * (version 2) code.  I've modified it to use doubles instead of
 * floats. -- T. A. Gardiner -- Aug. 12, 2003 
 *
 * Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1. */

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

#undef OFST
#undef KCOMP
#undef KWVM
