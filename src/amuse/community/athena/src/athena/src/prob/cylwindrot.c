#include "copyright.h"
/*============================================================================*/
/*! \file cylwindrot.c
 *  \brief The cylindrical analogue of the Bondi accretion (Parker wind) 
 *  problem with rotation.  Hydrodynamic, axisymmetric. 
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

static Real ang_mom,c_infty,lambda_s,vz0;
static int iprob;
static Real grav_pot(const Real x1, const Real x2, const Real x3);
static Real grav_acc(const Real x1, const Real x2, const Real x3);
Real myfunc(Real x, Real v);
static ConsS ***RootSoln=NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */
void problem(DomainS *pDomain)
{
  GridS *pG = pDomain->Grid;
  int i,j,k;
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  int nx1,nx2,nx3;
  Real x1,x2,x3;
  Real xs,vs,v,pgas0,pgas,alpha,beta,a,b,converged;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifdef MHD
  ath_error("[cylwindrot]: This problem only works in hydro!\n");
#endif

#ifndef CYLINDRICAL
  ath_error("[cylwindrot]: This problem only works in cylindrical!\n");
#endif

  if (nx1==1) {
    ath_error("[cylwindrot]: This problem can only be run in 2D or 3D!\n");
  }
  else if (nx2==1 && nx3>1) {
    ath_error("[cylwindrot]: Only (R,phi) can be used in 2D!\n");
  }

  /* ALLOCATE MEMORY FOR SOLUTION */
  if ((RootSoln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS))) == NULL)
    ath_error("[cylwindrot]: Error allocating memory for solution\n");

  ang_mom = par_getd("problem","ang_mom");
  c_infty = par_getd("problem","c_infty");
  vz0     = par_getd("problem","vz0");
  iprob   = par_geti("problem","iprob");
  printf("gamma = %f,\t ang_mom = %f,\t c_infty = %f\n", Gamma, ang_mom, c_infty);

  beta = 2.0*Gamma_1/(Gamma+1.0);
  xs = (3.0-Gamma+sqrt(SQR(Gamma-3.0)-16.0*SQR(ang_mom)))/4.0;
  lambda_s = 1.0/Gamma_1*pow(xs,beta)+pow(xs,beta-1.0)-0.5*SQR(ang_mom)*pow(xs,beta-2.0);
  lambda_s = pow(lambda_s/(0.5+1.0/Gamma_1),1.0/beta);
  vs = c_infty*pow(lambda_s/xs,0.5*beta);
  printf("xs = %13.10f,\t lambda_s = %13.10f,\t vs = %13.10f\n", xs, lambda_s, vs);

  // COMPUTE 1D WIND/ACCRETION SOLUTION
  for (i=il; i<=iu; i++) {
    cc_pos(pG,i,j,k,&x1,&x2,&x3);
    memset(&(pG->U[ks][js][i]),0.0,sizeof(ConsS));
    vs = pow(lambda_s/x1,0.5*beta);

    switch(iprob) {
      case 1: /* WIND */
              if (x1 < xs) {
                a = TINY_NUMBER;  b = vs;
              }
              else {
                a = vs;           b = HUGE_NUMBER;
              }
              break;
      case 2: /* ACCRETION */
              if (x1 < xs) {
                a = vs;           b = HUGE_NUMBER;
              }
              else {
                a = TINY_NUMBER;  b = vs;
              }
              break;
      default:  ath_error("[cylwindrot]:  Not an accepted problem number!\n");
    }

    converged = bisection(myfunc,a,b,x1,&v);
    if (!converged) ath_error("[cylwindrot]:  Bisection did not converge!\n");

    pG->U[ks][js][i].d   = lambda_s/(x1*v);
    pG->U[ks][js][i].M1  = lambda_s/x1;
    if (iprob==2)
      pG->U[ks][js][i].M1  *= -1.0;
    pG->U[ks][js][i].M2  = pG->U[ks][js][i].d*ang_mom/x1;
    pG->U[ks][js][i].M3  = pG->U[ks][js][i].d*vz0;

        /* INITIALIZE TOTAL ENERGY */
#ifndef ISOTHERMAL
        pgas0 = 1.0/Gamma;
        pgas = pgas0*pow(pG->U[ks][js][i].d,Gamma);
        pG->U[ks][js][i].E = pgas/Gamma_1
          + 0.5*(SQR(pG->U[ks][js][i].M1) + SQR(pG->U[ks][js][i].M2) + SQR(pG->U[ks][js][i].M3))/pG->U[ks][js][i].d;
#endif /* ISOTHERMAL */
  }

  /* COPY 1D SOLUTION AND SAVE */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pG->U[k][j][i] = pG->U[ks][js][i];
        RootSoln[k][j][i] = pG->U[ks][js][i];
      }
    }
  }

  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  bvals_mhd_fun(pDomain,left_x1,do_nothing_bc);
  bvals_mhd_fun(pDomain,right_x1,do_nothing_bc);

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
  compute_l1_error("CylWindRot", pM, RootSoln, 1);
}


/*=========================== PRIVATE FUNCTIONS ==============================*/

/*! \fn static Real grav_pot(const Real x1, const Real x2, const Real x3) 
 *  \brief Gravitational potential */
static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  return -SQR(c_infty)/x1;
}

/*! \fn static Real grav_acc(const Real x1, const Real x2, const Real x3)
  * \brief Gravitational acceleration */
static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  return SQR(c_infty/x1);
}

/*----------------------------------------------------------------------------*/
/*! \fn Real myfunc(Real x, Real v) 
 * \brief This function is used to calculate v as a function of x, gamma, and 
 *  lambda  using the bisection method.   */
/*----------------------------------------------------------------------------*/
Real myfunc(Real x, Real v) 
{
  return Gamma_1*(1/x + 1/Gamma_1 - 0.5*(SQR(v/c_infty)+SQR(ang_mom/x)))*pow(v*x/c_infty,Gamma_1) - pow(lambda_s,Gamma_1);
}
