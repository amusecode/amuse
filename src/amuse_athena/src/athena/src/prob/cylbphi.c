#include "copyright.h"
/*============================================================================*/
/*! \file cylbphi.c
 *  \brief A simple magnetostatic test of pressure balance using a B-field with 
 *  uniform phi-component. */
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
 * M2()       - phi-momentum
 * A3()       - magnetic vector potential
 *============================================================================*/

static Real omega0,rho0;
static int iprob;
static Real grav_pot(const Real x1, const Real x2, const Real x3);
static Real grav_acc(const Real x1, const Real x2, const Real x3);
Real M2(const Real x1, const Real x2, const Real x3);
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
  Real Eint,Emag,Ekin,vphi0,vz0,pgas0;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifndef MHD
  ath_error("[cylbphi]: This problem only works in MHD!\n");
#endif
#ifndef CYLINDRICAL
  ath_error("[cylbphi]: This problem only works in cylindrical!\n");
#endif

  if (nx1==1) {
    ath_error("[cylbphi]: This problem can only be run in 2D or 3D!\n");
  }
  else if (nx2==1 && nx3>1) {
    ath_error("[cylbphi]: Only (R,phi) can be used in 2D!\n");
  }

  omega0 = par_getd("problem", "omega0");
  vz0    = par_getd("problem", "vz0");
  bphi0  = par_getd("problem", "bphi0");
  rho0   = par_getd("problem", "rho0");
  pgas0  = par_getd("problem", "pgas0");
  iprob  = par_geti("problem", "iprob");


  /* ALLOCATE MEMORY FOR SOLUTION */
  if ((RootSoln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS))) == NULL)
    ath_error("[cylbphi]: Error allocating memory for solution!\n");


  /* SET DENSITY, MOMENTUM, AND MAGNETIC FIELDS
     iprob = 1, CONSTANT B-PHI, CONSTANT PRESSURE
     iprob = 2, B-PHI GOES AS 1/R, CONSTANT PRESSURE
  */

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        memset(&(pG->U[k][j][i]),0.0,sizeof(Gas));

        pG->U[k][j][i].d  = rho0;
        pG->U[k][j][i].M1 = 0.0;
        pG->U[k][j][i].M2 = avg1d(M2,pG,i,j,k);
        pG->U[k][j][i].M3 = pG->U[k][j][i].d*vz0;
        switch (iprob) {
          case 1:   // CONSTANT B_phi
                    pG->B2i[k][j][i]   = bphi0;
                    pG->U[k][j][i].B2c = bphi0;
                    break;
          case 2:   // B_phi GOES AS 1/R
                    pG->B2i[k][j][i]   = bphi0/x1;
                    pG->U[k][j][i].B2c = bphi0/x1;
                    break;
          default:  printf("[cylbphi]:  Not an accepted problem number\n");
        }

        /* INITIALIZE TOTAL ENERGY */
        Eint = pgas0/Gamma_1;
        Emag = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c)
                  + SQR(pG->U[k][j][i].B3c));
        Ekin = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2)
                  + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
        pG->U[k][j][i].E = Eint + Emag + Ekin;

        /* SAVE SOLUTION */
        RootSoln[k][j][i] = pG->U[ks][j][i];
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
//   printf("Max divB = %1.10e\n", compute_div_b(pM->Domain[0][0].Grid));
}

void Userwork_after_loop(MeshS *pM)
{
  compute_l1_error("CylBPhi", pM, RootSoln, 1);
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*! \fn static Real grav_pot(const Real x1, const Real x2, const Real x3) 
 *  \brief  Gravitational potential. */
static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  switch (iprob) {
    case 1:   return 0.5*SQR(x1*omega0) - (SQR(bphi0)/rho0)*log(x1);
              break;
    case 2:   return 0.5*SQR(x1*omega0);
              break;
    default:  return 0.0;
  }
}

/*! \fn static Real grav_acc(const Real x1, const Real x2, const Real x3) 
 *  \brief  Gravitational acceleration */
static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  switch (iprob) {
    case 1:   return x1*SQR(omega0) - SQR(bphi0)/(rho0*x1);
              break;
    case 2:   return x1*SQR(omega0);
              break;
    default:  return 0.0;
  }
}

/*! \fn Real M2(const Real x1, const Real x2, const Real x3) 
 *  \brief 2-component of momentum  */
Real M2(const Real x1, const Real x2, const Real x3) {
  return rho0*omega0*x1;
}
