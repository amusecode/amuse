#include "../copyright.h"
/*==============================================================================
 * FILE: get_eta.c
 *
 * PURPOSE: Functions to calculate the diffusion coefficients for resistivity.
 *   There are two standard prescriptions:
 *
 *   A. Single ion prescription (eq. 21 of Balbus & Terquem 2001, ApJ):
 *      E = eta_Ohm*J + (JxB)/(4pi*n_e*e/c) - [(JxB)xB]/(4pi*gamma*rho_i*rho)
 *
 *   B. Generalized prescription (eq. 25-31 of Wardle 2007, ApSS):
 *      E = eta_Ohm*J + eta_Hall*(JxB)/B + eta_AD*J_perp
 *   where eta_Ohm, eta_Hall and eta_AD are determined by sigma_O, sigma_H and
 *   sigma_P.
 *
 *   We provide two options for the user to set the diffusion coefficients:
 *
 *      CASE 1: eta_Ohm is constant and n_e and/or rho_i are proportional to
 *   some power (d_ind) of the gas density rho. The user needs to provide
 *   eta_Ohm and two proportional coefficients Q_Hall and Q_AD, which are
 *   defined in global.h.
 *
 *      CASE 2: The user provides his/her own function to evaluate the etas,
 *   e.g., spatially variable diffussivity, or adopting some chemical model,
 *   If the generalized prescription is used, the user can call the function
 *   convert_diffusion to convert the conductivitis to resistivities. The
 *   (blank) function get_eta_user() is given in the problem generator.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *  get_eta()          - main function call to get magnetic diffusivities
 *  eta_single_const() - constant diffusivities for single ion prescription
 *  eta_general()      - user defined diffusivities
 *  convert_diffusion() - convert conductivities to diffusion coefficients
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef RESISTIVITY

#ifdef HYDRO
#error : resistivity only works for MHD.
#endif /* HYDRO */

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* Get magnetic diffusivities
 */
void get_eta(GridS *pG)
{
  int i, il, iu, is = pG->is, ie = pG->ie;
  int j, jl, ju, js = pG->js, je = pG->je;
  int k, kl, ku, ks = pG->ks, ke = pG->ke;

  il = is - 2;
  iu = ie + 2;
  if (pG->Nx[1] > 1){
    jl = js - 2;
    ju = je + 2;
  } else {
    jl = js;
    ju = je;
  }
  if (pG->Nx[2] > 1){
    kl = ks - 2;
    ku = ke + 2;
  } else {
    kl = ks;
    ku = ke;
  }

  for (k=kl; k<=ku; k++) {
  for (j=jl; j<=ju; j++) {
  for (i=il; i<=iu; i++) {

     get_myeta(pG, i,j,k, &(pG->eta_Ohm[k][j][i]), 
                          &(pG->eta_Hall[k][j][i]), &(pG->eta_AD[k][j][i]));

  }}}

  return;
}

/*----------------------------------------------------------------------------*/
/* Single-ion prescription with constant diffusivities
 * NEED global parameters: eta_Ohm, Q_Hall and Q_AD.
 * ONLY eta_Ohm has the dimension of diffusivity.
 * This correspond to CASE 1.
 */

void eta_single_const(GridS *pG, int i, int j, int k,
                               Real *eta_O, Real *eta_H, Real *eta_A)
{
  Real Bsq, Bmag;

  *eta_O = eta_Ohm;

  if ((Q_Hall > 0.0) || (Q_AD > 0.0))
  {
    Bsq  = SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c)
                                   + SQR(pG->U[k][j][i].B3c);
    Bmag = sqrt(Bsq);

    *eta_H = Q_Hall * Bmag / pow(pG->U[k][j][i].d, d_ind);
    *eta_A = Q_AD * Bsq / pow(pG->U[k][j][i].d, 1.0+d_ind);
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* User defined resistivities.
 * This correspond to CASE 2.
 */

void eta_general(GridS *pG, int i, int j, int k,
                            Real *eta_O, Real *eta_H, Real *eta_A)
{
  get_eta_user(pG, i,j,k, eta_O, eta_H, eta_A);

  return;
}

/*----------------------------------------------------------------------------*/
/* Convert conductivities to diffusivities
 * NOTE: c^2/4pi factor is NOT included, so (eta x sigma) is dimensionless.
 * The conductivity coefficients should have the right dimension.
 */

void convert_diffusion(Real sigma_O, Real sigma_H, Real sigma_P,
                       Real *eta_O,  Real *eta_H,  Real *eta_A )
{
  Real sigma_perp2;

  sigma_perp2 = SQR(sigma_H) + SQR(sigma_P);

  *eta_O = 1.0/sigma_O;
  *eta_H = sigma_H / sigma_perp2;
  *eta_A = sigma_P / sigma_perp2 - *eta_O;

  return;
}

#endif /* RESISTIVITY */
