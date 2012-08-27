#define SEED 661979

#include "copyright.h"
/*============================================================================*/
/*! \file cylrayleigh.c
 *  \brief A test of the Rayleigh instability using omega(R) = omega_0/R^q.
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

static Real rho0,omega0,q;
static Real grav_pot(const Real x1, const Real x2, const Real x3);
static Real grav_acc(const Real x1, const Real x2, const Real x3);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */
void problem(DomainS *pDomain)
{
  GridS *pG = pDomain->Grid;
  int i,j,k;
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  int nx1,nx2,nx3,myid=0;
  Real x1,x2,x3,R1,R2;
  Real r,noise,omega,bphi0,pgas0,noise_level;
  Real Eint,Emag,Ekin;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifndef CYLINDRICAL
  ath_error("[cylrayleigh]: This problem only works in cylindrical!\n");
#endif

  if (nx1==1) {
    ath_error("[cylrayleigh]: This problem can only be run in 2D or 3D!\n");
  }
  else if (nx2==1 && nx3>1) {
    ath_error("[cylrayleigh]: Only (R,phi) can be used in 2D!\n");
  }

  /* SEED THE RANDOM NUMBER GENERATOR */
  srand(SEED + pG->my_id);

  omega0      = par_getd("problem", "omega0");
  bphi0       = par_getd("problem", "bphi0");
  rho0        = par_getd("problem", "rho0");
  pgas0       = par_getd("problem", "pgas0");
  q           = par_getd("problem", "q");
  noise_level = par_getd("problem", "noise_level");


  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        memset(&(pG->U[k][j][i]),0.0,sizeof(Gas));
        R1 = x1 - 0.5*pG->dx1;
        R2 = x1 + 0.5*pG->dx1;

        // RANDOM NUMBER BETWEEN 0 AND 1
        r = ((double) rand()/((double)RAND_MAX + 1.0));
        // RANDOM NUMBER BETWEEN +/- noise_level
        noise = noise_level*(2.0*r-1.0);

        pG->U[k][j][i].d  = rho0;
        pG->U[k][j][i].M2 = avg1d(M2,pG,i,j,k);
        // NOW PERTURB v_phi
        if ((i>=is) && (i<=ie)) {
          pG->U[k][j][i].M2 *= (1.0 + noise);
        }

#ifdef MHD
        pG->U[k][j][i].B2c = bphi0/x1;
        pG->B2i[k][j][i]   = bphi0/x1;
#endif

#ifndef ISOTHERMAL
        Eint = pgas0/Gamma_1;
        Emag = 0.0;
#ifdef MHD
        Emag = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
#endif
        Ekin = 0.5*(SQR(pG->U[k][j][i].M1 ) + SQR(pG->U[k][j][i].M2 ) + SQR(pG->U[k][j][i].M3 ))/pG->U[k][j][i].d;
        pG->U[k][j][i].E = Eint + Emag + Ekin;
#endif
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
  rho0   = par_getd("problem","rho0");
  omega0 = par_getd("problem","omega0");
  q      = par_getd("problem","q");

  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  bvals_mhd_fun(pDomain,left_x1,do_nothing_bc);
  bvals_mhd_fun(pDomain,right_x1,do_nothing_bc);
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
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*! \fn static Real grav_pot(const Real x1, const Real x2, const Real x3) 
 *  \brief Gravitational potential */
static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  if (q == 1.0) {
    return SQR(omega0)*log(x1);
  }
  else {
    Real omega = omega0/pow(x1,q);
    return 0.5*SQR(x1*omega)/(1.0-q);
  }
}

/*! \fn static Real grav_acc(const Real x1, const Real x2, const Real x3) {
 *  \brief Gravitational acceleration */
static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  Real omega = omega0/pow(x1,q);
  return x1*SQR(omega);
}

/*! \fn Real M2(const Real x1, const Real x2, const Real x3) 
 *  \brief 2-component of momentum */
Real M2(const Real x1, const Real x2, const Real x3) {
  return rho0*omega0*pow(x1,1.0-q);
}
