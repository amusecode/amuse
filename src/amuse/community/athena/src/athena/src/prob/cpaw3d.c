#include "copyright.h"
/*============================================================================*/
/*! \file cpaw3d.c
 *  \brief Problem generator for circularly polarized Alfven wave (CPAW) in 3D
 *   test.  
 *
 * PURPOSE: Problem generator for circularly polarized Alfven wave (CPAW) in 3D
 *   test.  The angles the wave propagates to the grid is automatically computed
 *   to be alpha12 = tan^{-1} (Y/X) and alpha23 = tan^{-1} (Z/Y). 
 *
 * The wave is defined with reference to a coordinate system (x,y,z) with
 *  transformation rules to the code coordinate system (x1,x2,x3)
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
 * Note these transformations are the same used in linear_wave3d.c
 *
 * Can be used for either standing (problem/v_par=1.0) or travelling
 * (problem/v_par=0.0) waves.
 *
 * USERWORK_AFTER_LOOP function computes L1 error norm in solution by comparing
 *   to initial conditions.  Problem must be evolved for an integer number of
 *   wave periods for this to work.
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - A1() - 1-component of vector potential for initial conditions
 * - A2() - 2-component of vector potential for initial conditions
 * - A3() - 3-component of vector potential for initial conditions
 *
 * REFERENCE: G. Toth,  "The div(B)=0 constraint in shock capturing MHD codes",
 *   JCP, 161, 605 (2000)						      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifndef MHD
#error : The cpaw3d test only works for MHD.
#endif

/* Initial solution, shared with Userwork_after_loop to compute L1 error */
static ConsS ***RootSoln=NULL;

/* Parameters which define initial solution -- made global so that they can be
 * shared with functions A1,2,3 which compute vector potentials */
static Real b_par, b_perp;
static Real ang_2, ang_3; /* Rotation angles about the y and z' axis */
static Real fac, sin_a2, cos_a2, sin_a3, cos_a3;
static Real lambda, k_par; /* Wavelength, 2*PI/wavelength */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * A1() - 1-component of vector potential for initial conditions
 * A2() - 2-component of vector potential for initial conditions
 * A3() - 3-component of vector potential for initial conditions
 *============================================================================*/

static Real A1(const Real x1, const Real x2, const Real x3);
static Real A2(const Real x1, const Real x2, const Real x3);
static Real A3(const Real x1, const Real x2, const Real x3);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  ConsS ***Soln;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,nx1,nx2,nx3,dir;
  Real dx1 = pGrid->dx1;
  Real dx2 = pGrid->dx2;
  Real dx3 = pGrid->dx3;
  Real hdx1 = 0.5*pGrid->dx1;
  Real hdx2 = 0.5*pGrid->dx2;
  Real hdx3 = 0.5*pGrid->dx3;
  Real x1size, x2size, x3size;
  Real x,x1,x2,x3,cs,sn;
  Real v_par, v_perp, den, pres;
#ifdef RESISTIVITY
  Real v_A, kva, omega_h, omega_l, omega_r;
#endif

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;

  if(pGrid->Nx[2] <= 1)
    ath_error("[problem]: cp_alfven3d assumes a 3D grid\n");

  if ((Soln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS)))==NULL)
    ath_error("[problem]: Error allocating memory for Soln\n");

  if (pDomain->Level == 0){
    if ((RootSoln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS)))
      == NULL) ath_error("[problem]: Error allocating memory for RootSoln\n");
  }

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
  b_par = par_getd("problem","b_par");
  den = 1.0;

  ath_pout(0,"va_parallel = %g\n",b_par/sqrt(den));

  b_perp = par_getd("problem","b_perp");
  v_perp = b_perp/sqrt((double)den);

  dir    = par_geti_def("problem","dir",1); /* right(1)/left(2) polarization */
  if (dir == 1) /* right polarization */
    fac = 1.0;
  else          /* left polarization */
    fac = -1.0;

#ifdef RESISTIVITY
  Q_Hall = par_getd("problem","Q_H");
  d_ind  = 0.0;
  v_A    = b_par/sqrt((double)den);
  if (Q_Hall > 0.0)
  {
    kva     = k_par*v_A;
    omega_h = 1.0/Q_Hall;

    omega_r = 0.5*SQR(kva)/omega_h*(sqrt(1.0+SQR(2.0*omega_h/kva)) + 1.0);
    omega_l = 0.5*SQR(kva)/omega_h*(sqrt(1.0+SQR(2.0*omega_h/kva)) - 1.0);

    if (dir == 1) /* right polarization (whistler wave) */
      v_perp = v_perp * kva / omega_r;
    else          /* left polarization */
      v_perp = v_perp * kva / omega_l;
  }
#endif

/* The gas pressure and parallel velocity are free parameters. */

  pres = par_getd("problem","pres");
  v_par = par_getd("problem","v_par");

/* Use the vector potential to initialize the interface magnetic fields */

  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	x1 -= 0.5*pGrid->dx1;
	x2 -= 0.5*pGrid->dx2;
	x3 -= 0.5*pGrid->dx3;

	pGrid->B1i[k][j][i] = (A3(x1,x2+dx2 ,x3+hdx3) - A3(x1,x2,x3+hdx3))/dx2 -
	                      (A2(x1,x2+hdx2,x3+dx3 ) - A2(x1,x2+hdx2,x3))/dx3;

	pGrid->B2i[k][j][i] = (A1(x1+hdx1,x2,x3+dx3 ) - A1(x1+hdx1,x2,x3))/dx3 -
	                      (A3(x1+dx1 ,x2,x3+hdx3) - A3(x1,x2,x3+hdx3))/dx1;

	pGrid->B3i[k][j][i] = (A2(x1+dx1,x2+hdx2,x3) - A2(x1,x2+hdx2,x3))/dx1 -
	                      (A1(x1+hdx1,x2+dx2,x3) - A1(x1+hdx1,x2,x3))/dx2;
      }
    }
  }

/* Now initialize the cell centered quantities */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	x = cos_a2*(x1*cos_a3 + x2*sin_a3) + x3*sin_a2;
	sn = sin(k_par*x);
	cs = fac*cos(k_par*x);

	pGrid->U[k][j][i].d  = den;
	pGrid->U[k][j][i].M1 = den*(v_par*cos_a2*cos_a3 +
				    v_perp*sn*sin_a3 +
				    v_perp*cs*sin_a2*cos_a3);

	pGrid->U[k][j][i].M2 = den*(v_par*cos_a2*sin_a3 -
				    v_perp*sn*cos_a3 +
				    v_perp*cs*sin_a2*sin_a3);

	pGrid->U[k][j][i].M3 = den*(v_par*sin_a2 -
				    v_perp*cs*cos_a2);

	pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i  ] +
			             pGrid->B1i[k][j][i+1]);
	pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j  ][i] +
			             pGrid->B2i[k][j+1][i]);
	pGrid->U[k][j][i].B3c = 0.5*(pGrid->B3i[k  ][j][i] +
			             pGrid->B3i[k+1][j][i]);

#ifndef ISOTHERMAL
	pGrid->U[k][j][i].E = pres/Gamma_1 
	  + 0.5*(SQR(pGrid->U[k][j][i].B1c) +
		 SQR(pGrid->U[k][j][i].B2c) +
		 SQR(pGrid->U[k][j][i].B3c) )
	  + 0.5*(SQR(pGrid->U[k][j][i].M1) +
		 SQR(pGrid->U[k][j][i].M2) +
		 SQR(pGrid->U[k][j][i].M3) )/den;
#endif

/* And store the initial state */
	Soln[k][j][i] = pGrid->U[k][j][i];
      }
    }
  }

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
      RootSoln[k][j][i].B1c = Soln[k][j][i].B1c;
      RootSoln[k][j][i].B2c = Soln[k][j][i].B2c;
      RootSoln[k][j][i].B3c = Soln[k][j][i].B3c;
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


ConsFun_t get_usr_expr(const char *expr){
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
 * Userwork_after_loop: computes L1-error in CPAW,
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
  total_error.B1c = 0.0;
  total_error.B2c = 0.0;
  total_error.B3c = 0.0;
#ifndef ISOTHERMAL
  total_error.E = 0.0;
#endif /* ISOTHERMAL */

/* Compute error only on root Grid, which is in Domain[0][0] */

  pGrid=pM->Domain[0][0].Grid;
  if (pGrid == NULL) return;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

/* compute L1 error in each variable, and rms total error */

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    error.d = 0.0;
    error.M1 = 0.0;
    error.M2 = 0.0;
    error.M3 = 0.0;
    error.B1c = 0.0;
    error.B2c = 0.0;
    error.B3c = 0.0;
#ifndef ISOTHERMAL
    error.E = 0.0;
#endif /* ISOTHERMAL */

    for (i=is; i<=ie; i++) {
      error.d   += fabs(pGrid->U[k][j][i].d   - RootSoln[k][j][i].d);
      error.M1  += fabs(pGrid->U[k][j][i].M1  - RootSoln[k][j][i].M1);
      error.M2  += fabs(pGrid->U[k][j][i].M2  - RootSoln[k][j][i].M2);
      error.M3  += fabs(pGrid->U[k][j][i].M3  - RootSoln[k][j][i].M3); 
      error.B1c += fabs(pGrid->U[k][j][i].B1c - RootSoln[k][j][i].B1c);
      error.B2c += fabs(pGrid->U[k][j][i].B2c - RootSoln[k][j][i].B2c);
      error.B3c += fabs(pGrid->U[k][j][i].B3c - RootSoln[k][j][i].B3c);
#ifndef ISOTHERMAL
      error.E   += fabs(pGrid->U[k][j][i].E   -  RootSoln[k][j][i].E);
#endif /* ISOTHERMAL */
    }

    total_error.d += error.d;
    total_error.M1 += error.M1;
    total_error.M2 += error.M2;
    total_error.M3 += error.M3;
    total_error.B1c += error.B1c;
    total_error.B2c += error.B2c;
    total_error.B3c += error.B3c;
#ifndef ISOTHERMAL
    total_error.E += error.E;
#endif /* ISOTHERMAL */
  }}

#if defined MPI_PARALLEL
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
  err[4] = total_error.B1c;
  err[5] = total_error.B2c;
  err[6] = total_error.B3c;
#ifndef ISOTHERMAL
  err[7] = total_error.E;
#endif /* ISOTHERMAL */

  ierr = MPI_Reduce(err,tot_err,8,MPI_DOUBLE,MPI_SUM,0,
    pM->Domain[0][0].Comm_Domain);

/* If I'm the parent, copy the sum back to the total_error variable */

  ierr = MPI_Comm_rank(pM->Domain[0][0].Comm_Domain, &myID);
  if(myID == 0){ /* I'm the parent */
    total_error.d   = tot_err[0];
    total_error.M1  = tot_err[1];
    total_error.M2  = tot_err[2];
    total_error.M3  = tot_err[3];
    total_error.B1c = tot_err[4];
    total_error.B2c = tot_err[5];
    total_error.B3c = tot_err[6];
#ifndef ISOTHERMAL
    total_error.E   = tot_err[7];
#endif /* ISOTHERMAL */
  }
  else return; /* The child grids do not do any of the following code */

#endif /* MPI_PARALLEL */

/* Compute RMS error over all variables, and print out */

  rms_error = SQR(total_error.d) + SQR(total_error.M1) + SQR(total_error.M2)
                + SQR(total_error.M3);
  rms_error += SQR(total_error.B1c) + SQR(total_error.B2c) 
               + SQR(total_error.B3c);
#ifndef ISOTHERMAL
  rms_error += SQR(total_error.E);
#endif /* ISOTHERMAL */
  rms_error = sqrt(rms_error)/(double)count;

/* Print error to file "cpaw3d-errors.dat" */
  
#ifdef MPI_PARALLEL
  fname = "../cpaw3d-errors.dat";
#else
  fname = "cpaw3d-errors.dat";
#endif
/* The file exists -- reopen the file in append mode */
  if((fp=fopen(fname,"r")) != NULL){
    if((fp = freopen(fname,"a",fp)) == NULL){
      ath_error("[Userwork_after_loop]: Unable to reopen file.\n");
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
    fprintf(fp,"  B1c  B2c  B3c");
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
  fprintf(fp,"  %e  %e  %e",
	  (total_error.B1c/(double)count),
	  (total_error.B2c/(double)count),
	  (total_error.B3c/(double)count));
  fprintf(fp,"\n");

  fclose(fp);

  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
  
/*----------------------------------------------------------------------------*/
/*! \fn static Real A1(const Real x1, const Real x2, const Real x3)
 *  \brief A1: 1-component of vector potential, using a gauge such that Ax = 0,
 * and Ay, Az are functions of x and y alone.
 */

static Real A1(const Real x1, const Real x2, const Real x3)
{
  Real x, y;
  Real Ay, Az;

  x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  y = -x1*sin_a3        + x2*cos_a3;

  Ay = fac*(b_perp/k_par)*sin(k_par*x);
  Az = (b_perp/k_par)*cos(k_par*x) + b_par*y;

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

  Ay = fac*(b_perp/k_par)*sin(k_par*x);
  Az = (b_perp/k_par)*cos(k_par*x) + b_par*y;

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

  Az = (b_perp/k_par)*cos(k_par*x) + b_par*y;

  return Az*cos_a2;
}
