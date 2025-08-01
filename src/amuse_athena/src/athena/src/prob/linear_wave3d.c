#include "copyright.h"
/*============================================================================*/
/*! \file linear_wave3d.c
 *  \brief Linear wave problem generator for 3D problems.
 *
 * PURPOSE: Linear wave problem generator for 3D problems. The generalization
 *   of linear_wave.c to handle 1D, 2D, and 3D problems seemed more complex
 *   than neccessary, so this special routine has been written to handle only
 *   3D problems.  This routine automatically sets the wavevector along the
 *   domain diagonal (cannot run grid-aligned waves).
 *
 *   Can be used for either standing (problem/vflow=1.0) or travelling
 *   (problem/vflow=0.0) waves.
 *
 *   Configure --with-gravity=fft to check Jeans stability of plane waves
 *   propagating parallel to grid.
 *
 *   Configure --enable-resistivity and/or --with-viscosity=ns to check
 *   damping of linear waves by resistivity/ambipolar diffusion and/or viscosity
 *
 * This code is most easily understood in terms of a one dimensional
 * problem in the coordinate system (x,y,z).  Two coordinate rotations are
 * applied to obtain a new wave vector in a 3D space in the (x1,x2,x3)
 * coordinate system.
 *
 *   First rotate about the y axis:
 *   - x' = x*cos(ang_2) - z*sin(ang_2)
 *   - y' = y
 *   - z' = x*sin(ang_2) + z*cos(ang_2)
 *
 *   Next rotate about the z' axis:
 *   - x = x'*cos(ang_3) - y'*sin(ang_3)
 *   - y = x'*sin(ang_3) + y'*cos(ang_3)
 *   - z = z'
 *
 *   Expanding this out we get:
 *   - x1 = x*cos(ang_2)*cos(ang_3) - y*sin(ang_3) - z*sin(ang_2)*cos(ang_3)
 *   - x2 = x*cos(ang_2)*sin(ang_3) - y*cos(ang_3) - z*sin(ang_2)*sin(ang_3)
 *   - x3 = x*sin(ang_2)                           + z*cos(ang_2)
 *
 *   This inverts to:
 *   - x =  x1*cos(ang_2)*cos(ang_3) + x2*cos(ang_2)*sin(ang_3) + x3*sin(ang_2)
 *   - y = -x1*sin(ang_3)            + x2*cos(ang_3)
 *   - z = -x1*sin(ang_2)*cos(ang_3) - x2*sin(ang_2)*sin(ang_3) + x3*cos(ang_2)
 *
 *   The magnetic field is given by:
 *   - B_x = b_par
 *   - B_y = b_perp*cos(k*x)
 *   - B_z = b_perp*sin(k*x)
 *   where k = 2.0*PI/lambda
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - A1() - 1-component of vector potential for initial conditions
 * - A2() - 2-component of vector potential for initial conditions
 * - A3() - 3-component of vector potential for initial conditions
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

/* Initial solution, shared with Userwork_after_loop to compute L1 error */
static ConsS ***RootSoln=NULL;

/* Parameters which define initial solution -- made global so that they can be
 * shared with functions A1,2,3 which compute vector potentials */
#ifdef MHD
static Real bx0, by0, bz0, dby, dbz;
#endif
static int wave_flag;
static Real ang_2, ang_3; /* Rotation angles about the y and z' axis */
static Real sin_a2, cos_a2, sin_a3, cos_a3;
static Real lambda, k_par; /* Wavelength, 2*PI/wavelength */

/*=============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * A1() - 1-component of vector potential for initial conditions
 * A2() - 2-component of vector potential for initial conditions
 * A3() - 3-component of vector potential for initial conditions
 *============================================================================*/

#ifdef MHD
static Real A1(const Real x1, const Real x2, const Real x3);
static Real A2(const Real x1, const Real x2, const Real x3);
static Real A3(const Real x1, const Real x2, const Real x3);
#endif /* MHD */

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  ConsS ***Soln;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,n,m,nx1,nx2,nx3;
  Real amp,vflow;
  Real d0,p0,u0,v0,w0,h0;
  Real x,x1,x2,x3,sn;
  Real Mx, My, Mz;
  Real ev[NWAVE],rem[NWAVE][NWAVE],lem[NWAVE][NWAVE];
  Real x1size, x2size, x3size;
#ifdef MHD
  Real xfact,yfact;
  Real ***a1,***a2,***a3;
#endif /* MHD */

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;

  if (nx3 <= 1)
    ath_error("[linear_wave3d]: Only true 3D problems can be run\n");

/* allocate memory for solution and vector potential */

#ifdef MHD
  if ((a1 = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL)
    ath_error("[problem]: Error allocating memory for \"a1\"\n");

  if ((a2 = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL)
    ath_error("[problem]: Error allocating memory for \"a2\"\n");

  if ((a3 = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL)
    ath_error("[problem]: Error allocating memory for \"a3\"\n");
#endif /* MHD */

  if ((Soln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS)))==NULL)
    ath_error("[problem]: Error allocating memory for \"Soln\"\n");

  if (pDomain->Level == 0){
    if ((RootSoln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS)))
      == NULL) ath_error("[problem]: Error alloc memory for RootSoln\n");
  }


/* Read initial conditions */

  wave_flag = par_geti("problem","wave_flag");
  amp = par_getd("problem","amp");
  vflow = par_getd("problem","vflow");

/* Imposing periodicity and one wavelength along each grid direction */

  x1size = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  x2size = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  x3size = pDomain->RootMaxX[2] - pDomain->RootMinX[2];

  ang_3 = atan(x1size/x2size);
  sin_a3 = sin(ang_3);
  cos_a3 = cos(ang_3);

  ang_2 = atan(0.5*(x1size*cos_a3 + x2size*sin_a3)/x3size);
  sin_a2 = sin(ang_2);
  cos_a2 = cos(ang_2);

  x1 = x1size*cos_a2*cos_a3;
  x2 = x2size*cos_a2*sin_a3;
  x3 = x3size*sin_a2;

/* For lambda choose the smaller of the 3 */
  lambda = x1;
  lambda = MIN(lambda,x2);
  lambda = MIN(lambda,x3);

/* Initialize k_parallel */
  k_par = 2.0*PI/lambda;

/* Get eigenmatrix, where the quantities u0 and bx0 are parallel to the
 * wavevector, and v0,w0,by0,bz0 are perpendicular. */

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


/* Now initialize wave at appropriate angle to grid */
/* Fields are initialized using vector potentials */

#ifdef MHD
/* Initialize the wave amplitudes for By and Bz */
  dby = amp*rem[NWAVE-2][wave_flag];
  dbz = amp*rem[NWAVE-1][wave_flag];

/* Initialize the Vector potential */
  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

	a1[k][j][i] = A1(x1, x2-0.5*pGrid->dx2, x3-0.5*pGrid->dx3);

	a2[k][j][i] = A2(x1-0.5*pGrid->dx1, x2, x3-0.5*pGrid->dx3);

	a3[k][j][i] = A3(x1-0.5*pGrid->dx1, x2-0.5*pGrid->dx2, x3);
      }
    }
  }

/* Initialize the interface fields in the Grid */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
	pGrid->B1i[k][j][i] = (a3[k  ][j+1][i] - a3[k][j][i])/pGrid->dx2 -
	                      (a2[k+1][j  ][i] - a2[k][j][i])/pGrid->dx3;
      }
    }
  }

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie; i++) {
	pGrid->B2i[k][j][i] = (a1[k+1][j][i  ] - a1[k][j][i])/pGrid->dx3 -
	                      (a3[k  ][j][i+1] - a3[k][j][i])/pGrid->dx1;
      }
    }
  }

  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	pGrid->B3i[k][j][i] = (a2[k][j  ][i+1] - a2[k][j][i])/pGrid->dx1 -
	                      (a1[k][j+1][i  ] - a1[k][j][i])/pGrid->dx2;
      }
    }
  }

/* compute cell-centered B-fields */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	Soln[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i]+pGrid->B1i[k  ][j  ][i+1]);
	Soln[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i]+pGrid->B2i[k  ][j+1][i  ]);
	Soln[k][j][i].B3c = 0.5*(pGrid->B3i[k][j][i]+pGrid->B3i[k+1][j  ][i  ]);
      }
    }
  }
#endif /* MHD */

/* The initial solution is stored into global 3D arrays so the error norm can
 * be computed at the end of the run.  After calculating these arrays, the
 * conserved variables are set equal to them */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

	x = cos_a2*(x1*cos_a3 + x2*sin_a3) + x3*sin_a2;

	sn = sin(k_par*x);

	Soln[k][j][i].d = d0 + amp*sn*rem[0][wave_flag];

	Mx = d0*vflow + amp*sn*rem[1][wave_flag];
	My = amp*sn*rem[2][wave_flag];
	Mz = amp*sn*rem[3][wave_flag];

	Soln[k][j][i].M1 = Mx*cos_a2*cos_a3 - My*sin_a3 - Mz*sin_a2*cos_a3;
	Soln[k][j][i].M2 = Mx*cos_a2*sin_a3 + My*cos_a3 - Mz*sin_a2*sin_a3;
	Soln[k][j][i].M3 = Mx*sin_a2                    + Mz*cos_a2;

#ifndef ISOTHERMAL
	Soln[k][j][i].E = p0/Gamma_1 + 0.5*d0*u0*u0 +
#ifdef MHD
	  0.5*(bx0*bx0 + by0*by0 + bz0*bz0) +
#endif /* MHD */
	  amp*sn*rem[4][wave_flag];
#endif /* ISOTHERMAL */
      }
    }
  }

/* Now set initial conditions to 3d wave solution */ 

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
#endif /* MHD */
      }
    }
  }

#ifdef MHD
  free_3d_array((void***)a1);
  free_3d_array((void***)a2);
  free_3d_array((void***)a3);
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

/*---------------------------------------------------------------------------
 * Userwork_after_loop: computes L1-error in linear waves,
 * ASSUMING WAVE HAS PROPAGATED AN INTEGER NUMBER OF PERIODS
 * Must set parameters in input file appropriately so that this is true
 */

void Userwork_after_loop(MeshS *pM)
{
  GridS *pGrid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke;
  Real rms_error=0.0;
  ConsS error,total_error;
  FILE *fp;
  char *fname;
  int Nx1, Nx2, Nx3, count;
#if defined MPI_PARALLEL
  double err[8], tot_err[8];
  int ierr,myID;
#endif

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
  }}

#ifdef MPI_PARALLEL
  Nx1 = pM->Domain[0][0].Nx[0];
  Nx2 = pM->Domain[0][0].Nx[1];
  Nx3 = pM->Domain[0][0].Nx[2];
#else
  Nx1 = ie - is + 1;
  Nx2 = je - js + 1;
  Nx3 = ke - ks + 1;
#endif
  count = Nx1*Nx2*Nx3;

#ifdef MPI_PARALLEL 
/* Now we have to use an All_Reduce to get the total error over all the MPI
 * grids.  Begin by copying the error into the err[] array */

  err[0] = total_error.d;
  err[1] = total_error.M1;
  err[2] = total_error.M2;
  err[3] = total_error.M3;
#ifdef MHD
  err[4] = total_error.B1c;
  err[5] = total_error.B2c;
  err[6] = total_error.B3c;
#endif /* MHD */
#ifndef ISOTHERMAL
  err[7] = total_error.E;
#endif /* ISOTHERMAL */

  /* Sum up the Computed Error */
  ierr = MPI_Reduce(err,tot_err,8,MPI_DOUBLE,MPI_SUM,0,
    pM->Domain[0][0].Comm_Domain);

/* If I'm the parent, copy the sum back to the total_error variable */

  ierr = MPI_Comm_rank(pM->Domain[0][0].Comm_Domain, &myID);
  if(myID == 0){ /* I'm the parent */
    total_error.d   = tot_err[0];
    total_error.M1  = tot_err[1];
    total_error.M2  = tot_err[2];
    total_error.M3  = tot_err[3];
#ifdef MHD
    total_error.B1c = tot_err[4];
    total_error.B2c = tot_err[5];
    total_error.B3c = tot_err[6];
#endif /* MHD */
#ifndef ISOTHERMAL
    total_error.E   = tot_err[7];
#endif /* ISOTHERMAL */
  }
  else return; /* The child grids do not do any of the following code */

#endif /* MPI_PARALLEL */

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

#ifdef MPI_PARALLEL
  fname = ath_fname("../","LinWave-errors",NULL,NULL,1,wave_flag,NULL,"dat");
#else
  fname = ath_fname(NULL,"LinWave-errors",NULL,NULL,1,wave_flag,NULL,"dat");
#endif

/* The file exists -- reopen the file in append mode */
  if((fp=fopen(fname,"r")) != NULL){
    if((fp = freopen(fname,"a",fp)) == NULL){
      ath_perr(-1,"[Userwork_after_loop]: Unable to reopen file.\n");
      free(fname);
      return;
    }
  }
/* The file does not exist -- open the file in write mode */
  else{
    if((fp = fopen(fname,"w")) == NULL){
      ath_perr(-1,"[Userwork_after_loop]: Unable to open file.\n");
      free(fname);
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
    fprintf(fp,"\n#\n");
  }

  fprintf(fp,"%d  %d  %d  %e",Nx1,Nx2,Nx3,rms_error);

  fprintf(fp,"  %e  %e  %e  %e",
	  (total_error.d/(double)count),
	  (total_error.M1/(double)count),
	  (total_error.M2/(double)count),
	  (total_error.M3/(double)count) );

#ifndef ISOTHERMAL
  fprintf(fp,"  %e",(total_error.E/(double)count) );
#endif /* ISOTHERMAL */

#ifdef MHD
  fprintf(fp,"  %e  %e  %e",
	  (total_error.B1c/(double)count),
	  (total_error.B2c/(double)count),
	  (total_error.B3c/(double)count));
#endif /* MHD */

  fprintf(fp,"\n");

  fclose(fp);
  free(fname);

  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/*! \fn static Real A1(const Real x1, const Real x2, const Real x3)
 *  \brief A1: 1-component of vector potential, using a gauge such that Ax = 0,
 * and Ay, Az are functions of x and y alone.
 */

#ifdef MHD
static Real A1(const Real x1, const Real x2, const Real x3)
{
  Real x, y;
  Real Ay, Az;

  x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  y = -x1*sin_a3        + x2*cos_a3;

  Ay =  bz0*x - (dbz/k_par)*cos(k_par*x);
  Az = -by0*x + (dby/k_par)*cos(k_par*x) + bx0*y;

  return -Ay*sin_a3 - Az*sin_a2*cos_a3;
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real A2(const Real x1, const Real x2, const Real x3)
 *  \brief A2: 2-component of vector potential
 */

static Real A2(const Real x1, const Real x2, const Real x3)
{
  Real x, y;
  Real Ay, Az;

  x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  y = -x1*sin_a3        + x2*cos_a3;

  Ay =  bz0*x - (dbz/k_par)*cos(k_par*x);
  Az = -by0*x + (dby/k_par)*cos(k_par*x) + bx0*y;

  return Ay*cos_a3 - Az*sin_a2*sin_a3;
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real A3(const Real x1, const Real x2, const Real x3)
 *  \brief A3: 3-component of vector potential
 */

static Real A3(const Real x1, const Real x2, const Real x3)
{
  Real x, y;
  Real Az;

  x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  y = -x1*sin_a3        + x2*cos_a3;

  Az = -by0*x + (dby/k_par)*cos(k_par*x) + bx0*y;

  return Az*cos_a2;
}
#endif /* MHD */
