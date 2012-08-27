#include "copyright.h"
/*============================================================================*/
/*! \file planet-disk.c
 *  \brief Problem generator for planet embedded in a disk, using the
 *   shearing sheet approximation.
 *
 * Code must be configured using --enable-shearing-box */
/*============================================================================*/
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * PlanetPot()   - static gravitational potential of planet
 * UnstratifiedDisk() - tidal potential in shearing box
 * expr_dV2()    - computes delta(Vy)
 * hst_*         - new history variables
 *============================================================================*/


void constant_iib(GridS *pGrid);
void constant_oib(GridS *pGrid);

static Real Mp,Rsoft,Xplanet,Yplanet,Zplanet;
static Real ramp_time,insert_time;
static Real PlanetPot(const Real x1, const Real x2, const Real x3);
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3);
static Real expr_dV2(const GridS *pG, const int i, const int j, const int k);
static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j,const int k);
static Real hst_rho_dVy2(const GridS *pG,const int i, const int j, const int k);
#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k);
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,BCFlag;
  Real x1,x2,x3;
  Real den = 1.0, pres = 1.0e-6;
  static int frst=1;  /* flag so new history variables enrolled only once */

#ifdef SHEARING_BOX
/* specify xy (r-phi) plane */
  ShBoxCoord = xy;
#endif

/* Read problem parameters.  Note Omega_0 set to 10^{-3} by default */
#ifdef SHEARING_BOX
  Omega_0 = par_getd_def("problem","omega",1.0e-3);
  qshear  = par_getd_def("problem","qshear",1.5);
#endif
  Mp      = par_getd_def("problem","Mplanet",0.0);
  Xplanet = par_getd_def("problem","Xplanet",0.0);
  Yplanet = par_getd_def("problem","Yplanet",0.0);
  Zplanet = par_getd_def("problem","Zplanet",0.0);
  Rsoft   = par_getd_def("problem","Rsoft",0.1);
  ramp_time = 0.0;
  insert_time = par_getd_def("problem","insert_time",0.0);

/* Compute field strength based on beta.  */
#ifdef ISOTHERMAL
  pres = Iso_csound2;
#endif

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

/* Initialize d, M, and P.  With FARGO do not initialize the background shear */

      pGrid->U[k][j][i].d  = den;
      pGrid->U[k][j][i].M1 = 0.0;
      pGrid->U[k][j][i].M2 = 0.0;
#ifdef SHEARING_BOX
#ifndef FARGO
      pGrid->U[k][j][i].M2 -= den*(qshear*Omega_0*x1);
#endif
#endif
      pGrid->U[k][j][i].M3 = 0.0;
#ifdef ADIABATIC
      pGrid->U[k][j][i].E = pres/Gamma_1
        + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2) 
             + SQR(pGrid->U[k][j][i].M3))/den;
#endif

    }
  }}

/* enroll gravitational potential of planet & shearing-box potential fns */

  StaticGravPot = PlanetPot;
  ShearingBoxPot = UnstratifiedDisk;

/* enroll new history variables, only once  */

  if (frst == 1) {
    dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
    dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
#ifdef ADIABATIC
    dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif
    frst = 0;
  }

/* With viscosity and/or resistivity, read eta_Ohm and nu_V */
#ifdef NAVIER_STOKES
  nu_V = par_getd("problem","nu");
#endif

/* Enroll outflow BCs if perdiodic BCs NOT selected.  This assumes the root
 * level grid is specified by the <domain1> block in the input file */

  BCFlag = par_geti_def("domain1","bc_ix1",0);
  if (BCFlag != 4) {
    if (pDomain->Disp[0] == 0) bvals_mhd_fun(pDomain, left_x1,  constant_iib);
  }
  BCFlag = par_geti_def("domain1","bc_ox1",0);
  if (BCFlag != 4) {
    if (pDomain->MaxX[0] == pDomain->RootMaxX[0])
      bvals_mhd_fun(pDomain, right_x1, constant_oib);
  }

  return;
}

/*==============================================================================
 * PUBLIC PROBLEM USER FUNCTIONS:
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

/* 'problem_read_restart' must enroll gravity on restarts */

void problem_read_restart(MeshS *pM, FILE *fp)
{
  int nl,nd,BCFlag_ix1,BCFlag_ox1;
/* Read Omega, and with viscosity and/or resistivity, read eta_Ohm and nu_V */

#ifdef SHEARING_BOX
  Omega_0 = par_getd_def("problem","omega",1.0e-3);
  qshear  = par_getd_def("problem","qshear",1.5);
#endif
  Mp      = par_getd_def("problem","Mplanet",0.0);
  Xplanet = par_getd_def("problem","Xplanet",0.0);
  Yplanet = par_getd_def("problem","Yplanet",0.0);
  Zplanet = par_getd_def("problem","Zplanet",0.0);
  Rsoft   = par_getd_def("problem","Rsoft",0.1);
  ramp_time = 0.0;
  insert_time = par_getd_def("problem","insert_time",0.0);
#ifdef NAVIER_STOKES
  nu_V = par_getd("problem","nu");
#endif

/* enroll gravitational potential of planet & shearing-box potential fns */

  StaticGravPot = PlanetPot;
  ShearingBoxPot = UnstratifiedDisk;

/* enroll new history variables */

  dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
  dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
#ifdef ADIABATIC
  dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif

  BCFlag_ix1 = par_geti_def("domain1","bc_ix1",0);
  BCFlag_ox1 = par_geti_def("domain1","bc_ox1",0);
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Disp[0] == 0 && BCFlag_ix1 != 4) 
        bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x1,  constant_iib);
      if (pM->Domain[nl][nd].MaxX[0] == pM->Domain[nl][nd].RootMaxX[0] 
          && BCFlag_ox1 != 4)
        bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x1, constant_oib);
    }
  }

  return;
}

/* Get_user_expression computes dVy */
ConsFun_t get_usr_expr(const char *expr)
{
  if(strcmp(expr,"dVy")==0) return expr_dV2;
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  ramp_time = pM->time;
}

void Userwork_after_loop(MeshS *pM)
{
}

/*------------------------------------------------------------------------------
 * PlanetPot:
 */
/*! \fn static Real PlanetPot(const Real x1, const Real x2, const Real x3)
 *  \brief static gravitational potential of planet */
static Real PlanetPot(const Real x1, const Real x2, const Real x3)
{
  Real rad,phi=0.0;
  rad = sqrt(SQR(x1-Xplanet) + SQR(x2-Yplanet) + SQR(x3-Zplanet));
  phi = -1.0*MIN(1.0,(ramp_time/(insert_time+0.0001)))*Mp/(rad+Rsoft);
  return phi;
}

/*------------------------------------------------------------------------------
 *! \fn static Real UnstratifiedDisk(const Real x1, const Real x2,const Real x3)
 *  \brief tidal potential in shearing box
 */

static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifdef SHEARING_BOX
#ifndef FARGO
  phi -= qshear*Omega_0*Omega_0*x1*x1;
#endif
#endif
  return phi;
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real expr_dV2(const GridS *pG, const int i, const int j, 
 *			     const int k)
 *  \brief Computes delta(Vy) 
 */

static Real expr_dV2(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef SHEARING_BOX
#ifdef FARGO
  return (pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  return (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
#endif
}

/*---------------------------------------------------------------------------*/
/* Hydro history variables:
 * hst_rho_Vx_dVy: Reynolds stress, added as history variable.
 * hst_rho_dVy2: KE in y-velocity fluctuations
 * hst_E_total: total energy (including tidal potential).
 */

/*! \fn static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j, 
 *				   const int k)
 *  \brief Reynolds stress, added as history variable.*/
static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef SHEARING_BOX
#ifdef FARGO
  return pG->U[k][j][i].M1*(pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  return pG->U[k][j][i].M1*
    (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
#endif
}

/*! \fn static Real hst_rho_dVy2(const GridS *pG, const int i, const int j,  
 *				const int k)
 *  \brief KE in y-velocity fluctuations */
static Real hst_rho_dVy2(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,dVy;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef SHEARING_BOX
#ifdef FARGO
  dVy = (pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  dVy = (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
#endif
  return pG->U[k][j][i].d*dVy*dVy;
}

#ifdef ADIABATIC
/*! \fn static Real hst_E_total(const GridS *pG, const int i, const int j, 
 *			        const int k)
 *  \brief total energy (including tidal potential). */
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,phi;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  phi = UnstratifiedDisk(x1, x2, x3);

  return pG->U[k][j][i].E + pG->U[k][j][i].d*phi;
}
#endif /* ADIABATIC */

/*---------------------------------------------------------------------------*/
/*! \fn void constant_iib(GridS *pGrid)
 *  \brief Sets boundary condition on left X boundary (iib) to constant
 * state (initial values).
 */

void constant_iib(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
  Real x1,x2,x3;
  Real den = 1.0, pres = 1.0e-6;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
      cc_pos(pGrid,(is-i),j,k,&x1,&x2,&x3);

/* Initialize d, M, and P.  With FARGO do not initialize the background shear */

      pGrid->U[k][j][is-i].d  = den;
      pGrid->U[k][j][is-i].M1 = 0.0;
      pGrid->U[k][j][is-i].M2 = 0.0;
#ifdef SHEARING_BOX
#ifndef FARGO
      pGrid->U[k][j][is-i].M2 -= den*(qshear*Omega_0*x1);
#endif
#endif
      pGrid->U[k][j][is-i].M3 = 0.0;
#ifdef ADIABATIC
      pGrid->U[k][j][is-i].E = pres/Gamma_1
        + 0.5*(SQR(pGrid->U[k][j][is-i].M1) + SQR(pGrid->U[k][j][is-i].M2)
             + SQR(pGrid->U[k][j][is-i].M3))/den;
#endif
      }
    }
  } 
}

/*---------------------------------------------------------------------------*/
/*! \fn void constant_oib(GridS *pGrid)
 *  \brief  Sets boundary condition on right X boundary (oib) to constant
 * state (initial values).
 */

void constant_oib(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
  Real x1,x2,x3;
  Real den = 1.0, pres = 1.0e-6;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
      cc_pos(pGrid,(ie+i),j,k,&x1,&x2,&x3);

/* Initialize d, M, and P.  With FARGO do not initialize the background shear */

      pGrid->U[k][j][ie+i].d  = den;
      pGrid->U[k][j][ie+i].M1 = 0.0;
      pGrid->U[k][j][ie+i].M2 = 0.0;
#ifdef SHEARING_BOX
#ifndef FARGO
      pGrid->U[k][j][ie+i].M2 -= den*(qshear*Omega_0*x1);
#endif
#endif
      pGrid->U[k][j][ie+i].M3 = 0.0;
#ifdef ADIABATIC
      pGrid->U[k][j][ie+i].E = pres/Gamma_1
        + 0.5*(SQR(pGrid->U[k][j][ie+i].M1) + SQR(pGrid->U[k][j][ie+i].M2)
             + SQR(pGrid->U[k][j][ie+i].M3))/den;
#endif
      }
    }
  }
}

