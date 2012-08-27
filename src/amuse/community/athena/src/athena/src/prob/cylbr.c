#include "copyright.h"
/*============================================================================*/
/*! \file cylbr.c
 *  \brief A dynamic test of force balance using a B_R-only, time-dependent,
 * non-axisymmetric magnetic field.
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
 * d()        - density
 * M2()       - phi-momentum
 * B1()       - magnetic field
 * B1i()      - R-interface magnetic field
 * Pgas()     - gas pressure
 * cylbr_ix1  - inner-R boundary condition
 * cylbr_ox1  - outer-R boundary condition
 *============================================================================*/

static Real br0,omega0,vz0,rho0,pgas0,a,x2min,x2max,phi0,x1save,x3save;
Real d(const Real x1, const Real x2, const Real x3);
Real M2(const Real x1, const Real x2, const Real x3);
Real B1(const Real x1, const Real x2, const Real x3);
Real B1i(const Real x2);
Real Pgas(const Real x1, const Real x2, const Real x3);
static Real grav_pot(const Real x1, const Real x2, const Real x3);
static Real grav_acc(const Real x1, const Real x2, const Real x3);
void cylbr_ix1(GridS *pG);
void cylbr_ox1(GridS *pG);
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
  Real x1,x2,x3,x1a,x2a,x2b;
  Real Eint,Ekin,Emag;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifndef CYLINDRICAL
  ath_error("[cylbr]: This problem only works in cylindrical!\n");
#endif

  if (nx1==1) {
    ath_error("[cylbr]: This problem can only be run in 2D or 3D!\n");
  }
  else if (nx2==1 && nx3>1) {
    ath_error("[cylbr]: Only (R,phi) can be used in 2D!\n");
  }

  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");
  a = 2.0*PI/(x2max-x2min);
  phi0 = a*(omega0*pG->time + x2min);

  omega0 = par_getd("problem", "omega0");
  vz0    = par_getd("problem", "vz0");
  br0    = par_getd("problem", "br0");
  rho0   = par_getd("problem", "rho0");
  pgas0  = par_getd("problem", "pgas0");

  /* ALLOCATE MEMORY FOR SOLUTION */
  if ((RootSoln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS))) == NULL)
    ath_error("[cylbr]: Error allocating memory for solution\n");

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        x1a = x1 - 0.5*pG->dx1;
        x2a = x2 - 0.5*pG->dx2;
        x2b = x2 + 0.5*pG->dx2;
        x1save = x1a;
        x3save = x3;

        pG->U[k][j][i].d   = avg2d(d,pG,i,j,k);
        pG->U[k][j][i].M2  = avg2d(M2,pG,i,j,k);
        pG->U[k][j][i].M3  = pG->U[k][j][i].d*vz0;
        pG->U[k][j][i].B1c = avg2d(B1,pG,i,j,k);
        pG->B1i[k][j][i]   = qsimp(B1i,x2a,x2b)/pG->dx2;

        Eint = avg2d(Pgas,pG,i,j,k)/Gamma_1;
        Emag = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
        Ekin = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d; 
        pG->U[k][j][i].E = Eint + Emag + Ekin;

        /* SAVE SOLUTION */
        RootSoln[k][j][i] = pG->U[k][j][i];
      }
    }
  }

  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  bvals_mhd_fun(pDomain,left_x1,cylbr_ix1);
  bvals_mhd_fun(pDomain,right_x1,cylbr_ox1);

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
  printf("Max divB = %1.10e\n", compute_div_b(pM->Domain[0][0].Grid));
}

void Userwork_after_loop(MeshS *pM)
{
  compute_l1_error("CylBR", pM, RootSoln, 1);
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*! \fn static Real grav_pot(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential */
static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  return 0.5*SQR(x1*omega0) - 0.5*SQR(br0)/SQR(x1);
}

/*! \fn static Real grav_acc(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational acceleration. */
static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  return x1*SQR(omega0) + SQR(br0)/pow(x1,3);
}

/*! \fn Real d(const Real x1, const Real x2, const Real x3) 
 *  \brief Density */
Real d(const Real x1, const Real x2, const Real x3) {
  return rho0 + SQR(sin(a*x2-phi0));
}

/*! \fn Real M2(const Real x1, const Real x2, const Real x3) 
 *  \brief 2-component of momentum */
Real M2(const Real x1, const Real x2, const Real x3) {
  return (rho0 + SQR(sin(a*x2-phi0)))*omega0*x1;
}

/*! \fn Real B1(const Real x1, const Real x2, const Real x3)
 *  \brief 1-component of B field */
Real B1(const Real x1, const Real x2, const Real x3) {
  return br0*(cos(a*x2-phi0))/x1;
}

/*! \fn Real B1i(const Real x2) 
 *  \brief 1-component of B field at interface */
Real B1i(const Real x2) {
  return B1(x1save,x2,x3save);
}

/*! \fn Real Pgas(const Real x1, const Real x2, const Real x3) 
 *  \brief gas pressure */
Real Pgas(const Real x1, const Real x2, const Real x3) {
  return pgas0 + 0.5*SQR(br0)*(rho0 + SQR(sin(a*x2-phi0)))/SQR(x1);
}

/*----------------------------------------------------------------------------*/
/*! \fn void cylbr_ix1(GridS *pG)
 *  \brief  Inner-R boundary conditions.  d, M2, B1, B1i, and P are all
 *   functions of R, phi, and t.
 */

void cylbr_ix1(GridS *pG)
{
  int is = pG->is;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  Real x1,x2,x3,x1a,x2a,x2b;
  Real Eint,Emag,Ekin;

  phi0 = a*(omega0*pG->time + x2min);
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        cc_pos(pG,is-i,j,k,&x1,&x2,&x3);
        x1a = x1 - 0.5*pG->dx1;
        x2a = x2 - 0.5*pG->dx2;
        x2b = x2 + 0.5*pG->dx2;
        x1save = x1a;
        x3save = x3;

        pG->U[k][j][is-i].d   = avg2d(d,pG,is-i,j,k);
        pG->U[k][j][is-i].M2  = avg2d(M2,pG,is-i,j,k);
        pG->U[k][j][is-i].M3  = pG->U[k][j][is-i].d*vz0;
        pG->U[k][j][is-i].B1c = avg2d(B1,pG,is-i,j,k);
        pG->B1i[k][j][is-i]   = qsimp(B1i,x2a,x2b)/pG->dx2;

        Eint = avg2d(Pgas,pG,is-i,j,k)/Gamma_1;
        Emag = 0.5*(SQR(pG->U[k][j][is-i].B1c) + SQR(pG->U[k][j][is-i].B2c) + SQR(pG->U[k][j][is-i].B3c));
        Ekin = 0.5*(SQR(pG->U[k][j][is-i].M1) + SQR(pG->U[k][j][is-i].M2) + SQR(pG->U[k][j][is-i].M3))/pG->U[k][j][is-i].d; 
        pG->U[k][j][is-i].E = Eint + Emag + Ekin;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void cylbr_ox1(GridS *pG)
 *  \brief B_R = B_0/R boundary conditions, Outer x1 boundary
 */
void cylbr_ox1(GridS *pG)
{
  int ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  Real x1,x2,x3,x1a,x2a,x2b;
  Real Eint,Emag,Ekin;

  phi0 = a*(omega0*pG->time + x2min);
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        cc_pos(pG,ie+i,j,k,&x1,&x2,&x3);
        x1a = x1 - 0.5*pG->dx1;
        x2a = x2 - 0.5*pG->dx2;  
        x2b = x2 + 0.5*pG->dx2;
        x1save = x1a;  
        x3save = x3;

        pG->U[k][j][ie+i].d   = avg2d(d,pG,ie+i,j,k);
        pG->U[k][j][ie+i].M2  = avg2d(M2,pG,ie+i,j,k);
        pG->U[k][j][ie+i].M3  = pG->U[k][j][ie+i].d*vz0;
        pG->U[k][j][ie+i].B1c = avg2d(B1,pG,ie+i,j,k);
        pG->B1i[k][j][ie+i]   = qsimp(B1i,x2a,x2b)/pG->dx2;

        Eint = avg2d(Pgas,pG,ie+i,j,k)/Gamma_1;
        Emag = 0.5*(SQR(pG->U[k][j][ie+i].B1c) + SQR(pG->U[k][j][ie+i].B2c) + SQR(pG->U[k][j][ie+i].B3c));
        Ekin = 0.5*(SQR(pG->U[k][j][ie+i].M1) + SQR(pG->U[k][j][ie+i].M2) + SQR(pG->U[k][j][ie+i].M3))/pG->U[k][j][ie+i].d; 
        pG->U[k][j][ie+i].E = Eint + Emag + Ekin;
      }
    }
  }

  return;
}
