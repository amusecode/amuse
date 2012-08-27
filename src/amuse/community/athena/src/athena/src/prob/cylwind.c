#include "copyright.h"
/*============================================================================*/
/*! \file cylwind.c
 * \brief The cylindrical analogue of the Bondi accretion (Parker wind) problem.
 *
 * Hydrodynamic, rotationless, axisymmetric.
 *
 * - REFERENCE: F. Shu, "The Physics of Astrophysics, Vol. II:  Gas Dynamics",
 *   1992.  ISBN 0935702652
 * - REFERENCE: L. Spitzer, "Physical Processes in the Interstellar Medium",
 *   1998.  ISBN 0471293350
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

static Real b0,lambda_s;
static int iprob;
Real grav_pot(const Real x1, const Real x2, const Real x3);
Real grav_acc(const Real x1, const Real x2, const Real x3);
Real myfunc(const Real x, const Real v);
static ConsS ***RootSoln=NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */
void problem(DomainS *pDomain)
{
  GridS *pG = pDomain->Grid;
  int i,j,k,converged;
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  int nx1,nx2,nx3; 
  Real x1,x2,x3,a,b;
  Real xs,vs,v,pgas0,pgas,beta;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifndef CYLINDRICAL
  ath_error("[cylwind]: This problem only works in cylindrical!\n");
#endif

  if (nx1==1) {
    ath_error("[cylwind]: Only R can be used in 1D!\n");
  }
  else if (nx2==1 && nx3>1) {
    ath_error("[cylwind]: Only (R,phi) can be used in 2D!\n");
  }

  if ((Soln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS))) == NULL)
    ath_error("[cylwind]: Error allocating memory for solution\n");

  b0    = par_getd("problem","b0");
  iprob = par_geti("problem","iprob");

  beta = 2.0*Gamma_1/(Gamma+1.0);
  xs = (3.0-Gamma)/2.0;
  lambda_s = pow(xs,(beta-1)/beta);
  vs = sqrt(2.0/(3.0-Gamma));
  printf("xs = %f, \tlambda_s = %f, \tvs = %f\n", xs, lambda_s, vs);

  for (i=il; i<=iu; i++) {
    cc_pos(pG,i,js,ks,&x1,&x2,&x3);
    memset(&(pG->U[ks][js][i]),0.0,sizeof(ConsS));

    switch(iprob) {
      case 1: /* PARKER WIND */
              if (x1 < xs)
                a = TINY_NUMBER;  b = vs;
              else
                a = vs;           b = HUGE_NUMBER;
              break;
      case 2: /* BONDI ACCRETION */
              if (x1 < xs)
                a = vs;           b = HUGE_NUMBER;
              else
                a = TINY_NUMBER;  b = vs;
              break;
      default:  ath_error("[cylwind]:  Not an accepted problem number!\n");
    }

    converged = bisection(myfunc,a,b,x1,&v);
    if (!converged) ath_error("[cylwind]:  Bisection did not converge!\n");

    pG->U[ks][js][i].d   = lambda_s/(x1*v);
    pG->U[ks][js][i].M1  = lambda_s/x1;
    if (iprob==2)
      pG->U[ks][js][i].M1  *= -1.0;

#ifdef MHD
    pG->U[ks][js][i].B1c = b0/x1;
    pG->B1i[ks][js][i]   = b0/(x1-0.5*pG->dx1);
#endif /* MHD */

    /* INITIALIZE TOTAL ENERGY */
#ifndef ISOTHERMAL
    pgas0 = 1.0/Gamma;
    pgas = pgas0*pow(pG->U[ks][js][i].d,Gamma);
    pG->U[ks][js][i].E = pgas/Gamma_1 + 0.5*SQR(pG->U[ks][js][i].M1)/pG->U[ks][js][i].d;
#endif /* ISOTHERMAL */
  }

  /* COPY 1D SOLUTION AND SAVE */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pG->U[k][j][i] = pG->U[ks][js][i];
        RootSoln[k][j][i] = pG->U[k][j][i];
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
  compute_l1_error("CylWind", pM, RootSoln, 1);
}


/*=========================== PRIVATE FUNCTIONS ==============================*/

/*! \fn Real grav_pot(const Real x1, const Real x2, const Real x3) 
 *  \brief Gravitational potential  */
Real grav_pot(const Real x1, const Real x2, const Real x3) {
  return -1.0/x1;
}

/*! \fn Real grav_acc(const Real x1, const Real x2, const Real x3) 
 *  \brief Gravitational acceleration */
Real grav_acc(const Real x1, const Real x2, const Real x3) {
  return 1.0/SQR(x1);
}

/*----------------------------------------------------------------------------*/
/*! \fn Real myfunc(const Real x, const Real v)
 *  \brief This funciton is used to calculate velocity v as a function of 
 *  position x
 *  using lambda_c, the critical value of the dimensionless mass wind/accretion
 *  rate.  Standard bisection is used to find the root(s). */
/*----------------------------------------------------------------------------*/
Real myfunc(const Real x, const Real v)
{
  return Gamma_1*(1/x + 1/Gamma_1 - 0.5*SQR(v))*pow(v*x,Gamma_1)
    - pow(lambda_s,Gamma_1);
}
