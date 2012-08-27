#include "copyright.h"
/*============================================================================*/
/*! \file cylunif.c
 *  \brief A test of conservation using uniform initial conditions.  
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
#include "cyl.h"

static Real br, bphi, omega, rho, pgas;
static int iprob;


/*! \fn static Real grav_pot(const Real x1, const Real x2, const Real x3) 
 *  \brief Gravitational potential */
static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  switch(iprob) {
    case 0:
            return 0.0;
            break;
    case 1:
            return 0.5*SQR(x1*omega);
            break;
    case 2:
            return -1.0/(rho*x1) + 0.5*SQR(x1*omega);
            break;
    case 3:
            return x1/rho + 0.5*SQR(x1*omega);
            break;
    case 4:
            return 0.0;
            break;
    default:
            ath_error("[cylunif]:  Invalid problem number!\n");
  }
}

/*! \fn static Real grav_acc(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational acceleration */
static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  switch(iprob) {
    case 0:
            return 0.0;
            break;
    case 1:
            return x1*SQR(omega);
            break;
    case 2:
            return 1.0/(rho*SQR(x1)) + x1*SQR(omega);
            break;
    case 3:
            return 1.0/rho + x1*SQR(omega);
            break;
    case 4:
            return 0.0;
            break;
    default:
            ath_error("[cylunif]:  Invalid problem number!\n");
  }
}

static Gas ***Soln=NULL;


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pG, Domain *pDomain)
{
  int i,j,k;
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  int nx1, nx2, nx3;
  Real x1, x2, x3, y1, y2, y3, r, noise, r1, r2;
  Real x1min, x1max, x2min, x2max, pb;
  Real ptotal;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifndef CYLINDRICAL
  ath_error("[cylunif]: This problem only works in cylindrical!\n");
#endif

  /* ALLOCATE MEMORY FOR SOLUTION */
  if ((Soln = (Gas***)calloc_3d_array(nx3,nx2,nx1,sizeof(Gas))) == NULL)
    ath_error("[cylunif]: Error allocating memory for solution\n");


  omega = par_getd("problem", "omega");
  br    = par_getd("problem", "br");
  bphi  = par_getd("problem", "bphi");
  rho   = par_getd("problem", "rho");
  pgas  = par_getd("problem", "pgas");
  iprob = par_geti("problem","iprob");


/* INITIALIZE DENSITY, MOMENTUM, AND B-FIELDS */


  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        vc_pos(pG,i,j,k,&y1,&y2,&y3);
        memset(&(pG->U[k][j][i]),0.0,sizeof(Gas));


        switch(iprob) {
          case 0:  // COMPLETELY STATIC
                  pG->U[k][j][i].d  = rho;
                  pG->U[k][j][i].E  = pgas/Gamma_1;
                  break;
          case 1:  // ROTATION, GRAVITY BALANCE (UNIF. PRESS.)
                  pG->U[k][j][i].d  = rho;
                  pG->U[k][j][i].M2 = pG->U[k][j][i].d*y1*omega;
                  pG->U[k][j][i].E  = pgas/Gamma_1;
                  break;
          case 2:  // ROTATION, PRESSURE, GRAVITY BALANCE (1/R PRESS.)
                  pG->U[k][j][i].d  = rho;
                  pG->U[k][j][i].M2 = pG->U[k][j][i].d*y1*omega;
                  pG->U[k][j][i].E  = (pgas/x1)/Gamma_1 + 0.5*SQR(pG->U[k][j][i].M2)/pG->U[k][j][i].d;
                  break;
          case 3:  // ROTATION, PRESSURE, GRAVITY BALANCE (R PRESS.)
                  pG->U[k][j][i].d  = rho;
                  pG->U[k][j][i].M2 = pG->U[k][j][i].d*y1*omega;
//                   pG->U[k][j][i].E  = (1.0+x1max-y1)/Gamma_1 + 0.5*SQR(pG->U[k][j][i].M2)/pG->U[k][j][i].d;
                  pG->U[k][j][i].E  = (pgas+5.0-y1)/Gamma_1 + 0.5*SQR(pG->U[k][j][i].M2)/pG->U[k][j][i].d;
#ifdef MHD
                  pG->U[k][j][i].B1c = br/x1;
                  pG->B1i[k][j][i]   = br/(x1-0.5*pG->dx1);
                  pG->U[k][j][i].E  += 0.5*SQR(pG->U[k][j][i].B1c);
#endif
                  break;
          case 4:  // TESTING...
                  pG->U[k][j][i].d  = rho*(SQR(x1)+0.25*SQR(pG->dx1));
                  pG->U[k][j][i].E  = pgas/Gamma_1;
                  break;
          default:
                  ath_error("[cylunif]:  Invalid problem number!\n");
        }

        // SAVE SOLUTION
        Soln[k][j][i] = pG->U[k][j][i];
      }
    }
  }



  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  set_bvals_mhd_fun(left_x1,do_nothing_bc);
  set_bvals_mhd_fun(right_x1,do_nothing_bc);

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
  return;
}

void problem_read_restart(Grid *pG, Domain *pD, FILE *fp)
{
  return;
}

Gasfun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VGFunout_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
}

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
  compute_l1_error("CylUnif", pGrid, pDomain, Soln, 0);
}
