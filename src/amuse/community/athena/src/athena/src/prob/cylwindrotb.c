#include "copyright.h"
/*============================================================================*/
/*! \file cylwindrotb.c
 *  \brief The cylindrical analogue of the Bondi accretion (Parker wind) 
 *  problem with rotation and magnetic field.  Axisymmetric.
 */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * grav_pot() - gravitational potential
 * grav_acc() - gravitational acceleration
 * myfunc()   - used to compute transonic solution
 *============================================================================*/

static Real theta,omega,eta,E,vz;
const Real rho_A = 1.0;
const Real R_A   = 1.0;
const Real GM    = 1.0;
static Real grav_pot(const Real x1, const Real x2, const Real x3);
static Real grav_acc(const Real x1, const Real x2, const Real x3);
Real myfunc(const Real x, const Real y);
static ConsS ***RootSoln=NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */
void problem(DomainS *pDomain)
{
  GridS *pG = pDomain->Grid;
  int i,j,k,n,converged;
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  int nx1, nx2, nx3;
  Real x1, x2, x3;
  Real a,b,c,d,xmin,xmax,ymin,ymax;
  Real x,y,xslow,yslow,xfast,yfast;
  Real R0,R1,R2,rho,Mdot,K,Omega,Pgas,beta,vR,BR,vphi,Bphi;
  ConsS *Wind=NULL;
  Real *pU=NULL,*pUl=NULL,*pUr=NULL;
  Real lsf,rsf;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifndef CYLINDRICAL
  ath_error("[cylwindrotb]: This problem only works in cylindrical!\n");
#endif
#ifndef MHD
  ath_error("[cylwindrotb]: This problem only works in MHD!\n");
#endif

  if (nx1==1) {
    ath_error("[cylwindrotb]: Only R can be used in 1D!\n");
  }
  else if (nx2==1 && nx3>1) {
    ath_error("[cylwindrotb]: Only (R,phi) can be used in 2D!\n");
  }

  /* ALLOCATE MEMORY FOR WIND SOLUTION */
  if ((Wind = (ConsS*)calloc_1d_array(nx1+1,sizeof(ConsS))) == NULL)
    ath_error("[cylwindrotb]: Error allocating memory\n");

  /* ALLOCATE MEMORY FOR GRID SOLUTION */
  if ((RootSoln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS))) == NULL)
    ath_error("[cylwindrotb]: Error allocating memory\n");

  theta = par_getd("problem","theta");
  omega = par_getd("problem","omega");
  vz    = par_getd("problem","vz");

  /* THIS NUMERICAL SOLUTION WAS OBTAINED FROM MATLAB.
   * IDEALLY, WE REPLACE THIS WITH A NONLINEAR SOLVER... */
  xslow = 0.5243264128;
  yslow = 2.4985859152;
  xfast = 1.6383327831;
  yfast = 0.5373957134;
  E     = 7.8744739104;
  eta   = 2.3608500383;

  xmin = par_getd("domain1","x1min")/R_A;
  xmax = par_getd("domain1","x1max")/R_A;
  ymin = 0.45/rho_A;
  ymax = 2.6/rho_A;

  printf("theta = %f,\t omega = %f,\t eta = %f,\t E = %f\n", theta,omega,eta,E);
  printf("xslow = %f,\t yslow = %f,\t xfast = %f,\t yfast = %f\n", xslow,yslow,xfast,yfast);
  printf("xmin = %f,\t ymin = %f,\t xmax = %f,\t ymax = %f\n", xmin,ymin,xmax,ymax);


  /* CALCULATE THE 1D WIND SOLUTION AT CELL-INTERFACES */
  for (i=il; i<=iu+1; i++) {
    memset(&(Wind[i]),0.0,sizeof(ConsS));
    cc_pos(pG,i,js,ks,&x1,&x2,&x3);

    /* WANT THE SOLUTION AT R-INTERFACES */
    R0 = x1 - 0.5*pG->dx1;
    x = R0/R_A;

    /* LOOK FOR A SIGN CHANGE INTERVAL */
    if (x < xslow) {
      sign_change(myfunc,yslow,10.0*ymax,x,&a,&b);
      sign_change(myfunc,b,10.0*ymax,x,&a,&b);
    } else if (x < 1.0) {
      sign_change(myfunc,1.0+TINY_NUMBER,yslow,x,&a,&b);
    } else if (x < xfast) {
      sign_change(myfunc,yfast,1.0-TINY_NUMBER,x,&a,&b);
      if (!sign_change(myfunc,b,1.0-TINY_NUMBER,x,&a,&b)) {
        a = yfast;
        b = 1.0-TINY_NUMBER;
      }
    } else {
      sign_change(myfunc,0.5*ymin,yfast,x,&a,&b);
    }

    /* USE BISECTION TO FIND THE ROOT */
    converged = bisection(myfunc,a,b,x,&y);
    if(!converged) {
      ath_error("[cylwindrotb]:  Bisection did not converge!\n");
    }

    /* CONSTRUCT THE SOLUTION */
    rho = rho_A*y;
    Mdot = sqrt(R_A*SQR(rho_A)*GM*eta);
    Omega = sqrt((GM*omega)/pow(R_A,3));
    K = (GM*theta)/(Gamma*pow(rho_A,Gamma_1)*R_A);
    Pgas = K*pow(rho,Gamma);
    vR = Mdot/(R0*rho);
    beta = sqrt(1.0/rho_A);
    BR = beta*rho*vR;
    vphi = R0*Omega*(1.0/SQR(x)-y)/(1.0-y);
    Bphi = beta*rho*(vphi-R0*Omega);

    Wind[i].d   = rho;
    Wind[i].M1  = rho*vR;
    Wind[i].M2  = rho*vphi;
    Wind[i].M3  = rho*vz;
    Wind[i].B1c = BR;
    Wind[i].B2c = Bphi;
    Wind[i].B3c = 0.0;
    Wind[i].E   = Pgas/Gamma_1
      + 0.5*(SQR(Wind[i].B1c) + SQR(Wind[i].B2c) + SQR(Wind[i].B3c))
      + 0.5*(SQR(Wind[i].M1 ) + SQR(Wind[i].M2 ) + SQR(Wind[i].M3 ))/Wind[i].d;
  }

  /* AVERAGE THE WIND SOLUTION ACROSS THE ZONE FOR CC VARIABLES */
  for (i=il; i<=iu; i++) {
    memset(&(pG->U[ks][js][i]),0.0,sizeof(ConsS));
    cc_pos(pG,i,js,ks,&x1,&x2,&x3);
    lsf = (x1 - 0.5*pG->dx1)/x1;
    rsf = (x1 + 0.5*pG->dx1)/x1;

    pU  = (Real*)&(pG->U[ks][js][i]);
    pUl = (Real*)&(Wind[i]);
    pUr = (Real*)&(Wind[i+1]);
    for (n=0; n<NWAVE; n++) {
      pU[n] = 0.5*(lsf*pUl[n] + rsf*pUr[n]);
    }

    pG->B1i[ks][js][i]   = Wind[i].B1c;
    pG->B2i[ks][js][i]   = Wind[i].B2c;
    pG->B3i[ks][js][i]   = Wind[i].B3c;
  }

  /* COPY 1D SOLUTION ACROSS THE GRID AND SAVE */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pG->U[k][j][i] = pG->U[ks][js][i];
        RootSoln[k][j][i]  = pG->U[ks][js][i];
      }
    }
  }

  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  bvals_mhd_fun(pDomain,left_x1,do_nothing_bc);
  bvals_mhd_fun(pDomain,right_x1,do_nothing_bc);

  free_1d_array((void *)Wind);

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

void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
  compute_l1_error("CylWindRotB", pM, RootSoln, 1);
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*! \fn static Real grav_pot(const Real x1, const Real x2, const Real x3) 
 *  \brief Gravitational potential */
static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  return -GM/x1;
}

/*! \fn static Real grav_acc(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational acceleration */
static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  return GM/SQR(x1);
}

/*----------------------------------------------------------------------------*/
/*! \fn Real myfunc(const Real x, const Real y) 
 *  \brief This function is used to calculate y (ie. rho) as a function of x 
 *  (ie. R), gamma, eta, theta, omega, and E using the bisection method.
 *
 */

Real myfunc(const Real x, const Real y) 
{
  return eta/(2.0*pow(x,2)*pow(y,2)) + (theta/Gamma_1)*pow(y,Gamma_1) 
    - 1.0/x + 0.5*omega*(pow(x-1.0/x,2)/pow(y-1,2) - pow(x,2)) - E;
}
