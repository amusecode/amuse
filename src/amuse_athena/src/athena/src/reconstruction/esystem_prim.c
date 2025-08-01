#include "../copyright.h"
/*============================================================================*/
/*! \file esystem_prim.c
 *  \brief Functions to evaluate the eigenvalues, and left- and
 * right-eigenvectors for the linearized system in the 
 * PRIMITIVE variables.
 *
 * PURPOSE: Functions to evaluate the eigenvalues, and left- and
 * right-eigenvectors for the linearized system in the 
 * PRIMITIVE variables, i.e. W,t = AW,x, where W=(d,vx,vy,vz,[P],[By,Bz]).
 * The eigenvalues are returned through the argument list as a vector of length
 * NWAVE.  The eigenvectors are returned as matrices of size (NWAVE)x(NWAVE),
 * with right-eigenvectors stored as COLUMNS (R_i[*] = right_eigenmatrix[*][i]),
 * and left-eigenvectors stored as ROWS (L_i[*] = left_eigenmatrix[i][*]).
 *
 * To improve performance components of the eigenvectors which are zero
 * are not set here (eigenmatrices must be initialized to zero in calling
 * routine).   However, for completeness statements which set these values
 * are included, but are commented out.
 *
 * REFERENCES:
 * - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new
 *   code for astrophysical MHD", ApJS, (2008), Appendix A
 *   Equation numbers refer to this paper.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - esys_prim_iso_hyd() - isothermal hydrodynamics
 * - esys_prim_adb_hyd() - adiabatic hydrodynamics
 * - esys_prim_iso_mhd() - isothermal MHD
 * - esys_prim_adb_mhd() - adiabatic MHD
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

/*----------------------------------------------------------------------------*/
/*! \fn void esys_prim_iso_hyd(const Real d, const Real v1, 
 *  Real eigenvalues[],
 *  Real right_eigenmatrix[][4], Real left_eigenmatrix[][4])
 *  \brief ISOTHERMAL HYDRO
 *   Input: d, v1 = primitive variables
 *   Output: eigenvalues[4], right_eigenmatrix[4,4], left_eigenmatrix[4,4];
 */

#if defined(BAROTROPIC) && defined(HYDRO)
void esys_prim_iso_hyd(const Real d, const Real v1, 
  Real eigenvalues[],
  Real right_eigenmatrix[][4], Real left_eigenmatrix[][4])
{

/* Compute eigenvalues (eq. A6) */

  eigenvalues[0] = v1 - Iso_csound;
  eigenvalues[1] = v1;
  eigenvalues[2] = v1;
  eigenvalues[3] = v1 + Iso_csound;

/* Right-eigenvectors, stored as COLUMNS (eq. A3) */

  right_eigenmatrix[0][0] = 1.0;
  right_eigenmatrix[1][0] = -Iso_csound/d;
/*right_eigenmatrix[2][0] = 0.0; */
/*right_eigenmatrix[3][0] = 0.0; */

/*right_eigenmatrix[0][1] = 0.0; */
/*right_eigenmatrix[1][1] = 0.0; */
  right_eigenmatrix[2][1] = 1.0;
/*right_eigenmatrix[3][1] = 0.0; */

/*right_eigenmatrix[0][2] = 0.0; */
/*right_eigenmatrix[1][2] = 0.0; */
/*right_eigenmatrix[2][2] = 0.0; */
  right_eigenmatrix[3][2] = 1.0;

  right_eigenmatrix[0][3] = 1.0;
  right_eigenmatrix[1][3] = Iso_csound/d;
/*right_eigenmatrix[2][3] = 0.0; */
/*right_eigenmatrix[3][3] = 0.0; */

/* Left-eigenvectors, stored as ROWS (eq. A7) */

  left_eigenmatrix[0][0] = 0.5;
  left_eigenmatrix[0][1] = -0.5*d/Iso_csound;
/*left_eigenmatrix[0][2] = 0.0; */
/*left_eigenmatrix[0][3] = 0.0; */

/*left_eigenmatrix[1][0] = 0.0; */
/*left_eigenmatrix[1][1] = 0.0; */
  left_eigenmatrix[1][2] = 1.0;
/*left_eigenmatrix[1][3] = 0.0; */

/*left_eigenmatrix[2][0] = 0.0; */
/*left_eigenmatrix[2][1] = 0.0; */
/*left_eigenmatrix[2][2] = 0.0; */
  left_eigenmatrix[2][3] = 1.0;

  left_eigenmatrix[3][0] = 0.5;
  left_eigenmatrix[3][1] = 0.5*d/Iso_csound;
/*left_eigenmatrix[3][2] = 0.0; */
/*left_eigenmatrix[3][3] = 0.0; */

}
#endif

/*----------------------------------------------------------------------------*/
/*! \fn void esys_prim_adb_hyd(const Real d, const Real v1, const Real rho_a2,
 *  Real eigenvalues[],
 *  Real right_eigenmatrix[][5], Real left_eigenmatrix[][5])
 *  \brief ADIABATIC HYDRO
 *   Input: d, v1, p = primitive variables
 *   Output: eigenvalues[5], right_eigenmatrix[5,5], left_eigenmatrix[5,5];
 */

#if !defined(BAROTROPIC) && defined(HYDRO)
void esys_prim_adb_hyd(const Real d, const Real v1, const Real rho_a2,
  Real eigenvalues[],
  Real right_eigenmatrix[][5], Real left_eigenmatrix[][5])
{
  Real asq,a;
  asq = rho_a2/d;
  a = sqrt(asq);

/* Compute eigenvalues (eq. A2) */

  eigenvalues[0] = v1 - a;
  eigenvalues[1] = v1;
  eigenvalues[2] = v1;
  eigenvalues[3] = v1;
  eigenvalues[4] = v1 + a;

/* Right-eigenvectors, stored as COLUMNS (eq. A3) */

  right_eigenmatrix[0][0] = 1.0;
  right_eigenmatrix[1][0] = -a/d;
/*right_eigenmatrix[2][0] = 0.0; */
/*right_eigenmatrix[3][0] = 0.0; */
  right_eigenmatrix[4][0] = asq;

  right_eigenmatrix[0][1] = 1.0;
/*right_eigenmatrix[1][1] = 0.0; */
/*right_eigenmatrix[2][1] = 0.0; */
/*right_eigenmatrix[3][1] = 0.0; */
/*right_eigenmatrix[4][1] = 0.0; */

/*right_eigenmatrix[0][2] = 0.0; */
/*right_eigenmatrix[1][2] = 0.0; */
  right_eigenmatrix[2][2] = 1.0;
/*right_eigenmatrix[3][2] = 0.0; */
/*right_eigenmatrix[4][2] = 0.0; */

/*right_eigenmatrix[0][3] = 0.0; */
/*right_eigenmatrix[1][3] = 0.0; */
/*right_eigenmatrix[2][3] = 0.0; */
  right_eigenmatrix[3][3] = 1.0;
/*right_eigenmatrix[4][3] = 0.0; */

  right_eigenmatrix[0][4] = 1.0;
  right_eigenmatrix[1][4] = -right_eigenmatrix[1][0];
/*right_eigenmatrix[2][4] = 0.0; */
/*right_eigenmatrix[3][4] = 0.0; */
  right_eigenmatrix[4][4] = asq;

/* Left-eigenvectors, stored as ROWS (eq. A4) */

/*left_eigenmatrix[0][0] = 0.0; */
  left_eigenmatrix[0][1] = -0.5*d/a;
/*left_eigenmatrix[0][2] = 0.0; */
/*left_eigenmatrix[0][3] = 0.0; */
  left_eigenmatrix[0][4] = 0.5/asq;

  left_eigenmatrix[1][0] = 1.0;
/*left_eigenmatrix[1][1] = 0.0; */
/*left_eigenmatrix[1][2] = 0.0; */
/*left_eigenmatrix[1][3] = 0.0; */
  left_eigenmatrix[1][4] = -1.0/asq;

/*left_eigenmatrix[2][0] = 0.0; */
/*left_eigenmatrix[2][1] = 0.0; */
  left_eigenmatrix[2][2] = 1.0;
/*left_eigenmatrix[2][3] = 0.0; */
/*left_eigenmatrix[2][4] = 0.0; */

/*left_eigenmatrix[3][0] = 0.0; */
/*left_eigenmatrix[3][1] = 0.0; */
/*left_eigenmatrix[3][2] = 0.0; */
  left_eigenmatrix[3][3] = 1.0;
/*left_eigenmatrix[3][4] = 0.0; */

/*left_eigenmatrix[4][0] = 0.0; */
  left_eigenmatrix[4][1] = -left_eigenmatrix[0][1];
/*left_eigenmatrix[4][2] = 0.0; */
/*left_eigenmatrix[4][3] = 0.0; */
  left_eigenmatrix[4][4] = left_eigenmatrix[0][4];
}
#endif

/*----------------------------------------------------------------------------*/
/*! \fn void esys_prim_iso_mhd(const Real d, const Real v1, const Real b1, 
 *  const Real b2, const Real b3, Real eigenvalues[],
 *  Real right_eigenmatrix[][6], Real left_eigenmatrix[][6])
 *  \brief ISOTHERMAL MHD
 *   Input: d, v1, b1, b2, b3 = density, velocities, and B field
 *   Output: eigenvalues[6], right_eigenmatrix[6,6], left_eigenmatrix[6,6];
 */

#if defined(BAROTROPIC) && defined(MHD)
void esys_prim_iso_mhd(const Real d, const Real v1, const Real b1,
  const Real b2, const Real b3, Real eigenvalues[],
  Real right_eigenmatrix[][6], Real left_eigenmatrix[][6])
{
  Real btsq,vaxsq,cfsq,cf,cssq,cs,bt,bet2,bet3,alpha_f,alpha_s;
  Real sqrtd,s,qf,qs,af,as,vax,norm,af_prime,as_prime;
  Real ct2,tsum,tdif,cf2_cs2;
  Real di = 1.0/d;
  btsq  = b2*b2 + b3*b3;
  vaxsq = b1*b1*di;

/* Compute fast- and slow-magnetosonic speeds (eq. A10) */

  ct2 = btsq*di;
  tsum = vaxsq + ct2 + (Iso_csound2);
  tdif = vaxsq + ct2 - (Iso_csound2);
  cf2_cs2 = sqrt((double)(tdif*tdif + 4.0*(Iso_csound2)*ct2));

  cfsq = 0.5*(tsum + cf2_cs2);
  cf = sqrt((double)cfsq);

  cssq = (Iso_csound2)*vaxsq/cfsq;
  cs = sqrt((double)cssq);

/* Compute beta(s) (eq A17) */

  bt  = sqrt(btsq);
  if (bt == 0.0) {
    bet2 = 1.0;
    bet3 = 0.0;
  } else {
    bet2 = b2/bt;
    bet3 = b3/bt;
  }

/* Compute alpha(s) (eq A16) */

  if ((cfsq-cssq) == 0.0) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  } else if ( (Iso_csound2 - cssq) <= 0.0) {
    alpha_f = 0.0;
    alpha_s = 1.0;
  } else if ( (cfsq - Iso_csound2) <= 0.0) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  } else {
    alpha_f = sqrt((Iso_csound2 - cssq)/(cfsq - cssq));
    alpha_s = sqrt((cfsq - Iso_csound2)/(cfsq - cssq));
  }

/* Compute Q(s) and A(s) (eq. A14-15), etc. */

  sqrtd = sqrt(d);
  s = SIGN(b1);
  qf = cf*alpha_f*s;
  qs = cs*alpha_s*s;
  af = Iso_csound*alpha_f*sqrtd;
  as = Iso_csound*alpha_s*sqrtd;

/* Compute eigenvalues (eq. A21) */

  vax = sqrt(vaxsq);
  eigenvalues[0] = v1 - cf;
  eigenvalues[1] = v1 - vax;
  eigenvalues[2] = v1 - cs;
  eigenvalues[3] = v1 + cs;
  eigenvalues[4] = v1 + vax;
  eigenvalues[5] = v1 + cf;

/* Right-eigenvectors, stored as COLUMNS (eq A12) */

  right_eigenmatrix[0][0] = d*alpha_f;
  right_eigenmatrix[1][0] = -cf*alpha_f;
  right_eigenmatrix[2][0] = qs*bet2;
  right_eigenmatrix[3][0] = qs*bet3;
  right_eigenmatrix[4][0] = as*bet2;
  right_eigenmatrix[5][0] = as*bet3;

/*right_eigenmatrix[0][1] = 0.0; */
/*right_eigenmatrix[1][1] = 0.0; */
  right_eigenmatrix[2][1] = -bet3;
  right_eigenmatrix[3][1] = bet2;
  right_eigenmatrix[4][1] = -bet3*s*sqrtd;
  right_eigenmatrix[5][1] = bet2*s*sqrtd;

  right_eigenmatrix[0][2] = d*alpha_s;
  right_eigenmatrix[1][2] = -cs*alpha_s;
  right_eigenmatrix[2][2] = -qf*bet2;
  right_eigenmatrix[3][2] = -qf*bet3;
  right_eigenmatrix[4][2] = -af*bet2;
  right_eigenmatrix[5][2] = -af*bet3;

  right_eigenmatrix[0][3] = right_eigenmatrix[0][2];
  right_eigenmatrix[1][3] = -right_eigenmatrix[1][2];
  right_eigenmatrix[2][3] = -right_eigenmatrix[2][2];
  right_eigenmatrix[3][3] = -right_eigenmatrix[3][2];
  right_eigenmatrix[4][3] = right_eigenmatrix[4][2];
  right_eigenmatrix[5][3] = right_eigenmatrix[5][2];

/*right_eigenmatrix[0][4] = 0.0; */
/*right_eigenmatrix[1][4] = 0.0; */
  right_eigenmatrix[2][4] = bet3;
  right_eigenmatrix[3][4] = -bet2;
  right_eigenmatrix[4][4] = right_eigenmatrix[4][1];
  right_eigenmatrix[5][4] = right_eigenmatrix[5][1];

  right_eigenmatrix[0][5] = right_eigenmatrix[0][0];
  right_eigenmatrix[1][5] = -right_eigenmatrix[1][0];
  right_eigenmatrix[2][5] = -right_eigenmatrix[2][0];
  right_eigenmatrix[3][5] = -right_eigenmatrix[3][0];
  right_eigenmatrix[4][5] = right_eigenmatrix[4][0];
  right_eigenmatrix[5][5] = right_eigenmatrix[5][0];

/* Left-eigenvectors, stored as ROWS (eq A22) */

  norm = 0.5/Iso_csound2;
  qf = norm*qf;
  qs = norm*qs;
  af_prime = norm*af*di;
  as_prime = norm*as*di;

  left_eigenmatrix[0][0] = norm*alpha_f*Iso_csound2*di;
  left_eigenmatrix[0][1] = -norm*cf*alpha_f;
  left_eigenmatrix[0][2] = qs*bet2;
  left_eigenmatrix[0][3] = qs*bet3;
  left_eigenmatrix[0][4] = as_prime*bet2;
  left_eigenmatrix[0][5] = as_prime*bet3;

/*left_eigenmatrix[1][0] = 0.0; */
/*left_eigenmatrix[1][1] = 0.0; */
  left_eigenmatrix[1][2] = -0.5*bet3;
  left_eigenmatrix[1][3] = 0.5*bet2;
  left_eigenmatrix[1][4] = -0.5*bet3*s/sqrtd;
  left_eigenmatrix[1][5] = 0.5*bet2*s/sqrtd;

  left_eigenmatrix[2][0] = norm*alpha_s*Iso_csound2*di;
  left_eigenmatrix[2][1] = -norm*cs*alpha_s;
  left_eigenmatrix[2][2] = -qf*bet2;
  left_eigenmatrix[2][3] = -qf*bet3;
  left_eigenmatrix[2][4] = -af_prime*bet2;
  left_eigenmatrix[2][5] = -af_prime*bet3;

  left_eigenmatrix[3][0] = left_eigenmatrix[2][0];
  left_eigenmatrix[3][1] = -left_eigenmatrix[2][1];
  left_eigenmatrix[3][2] = -left_eigenmatrix[2][2];
  left_eigenmatrix[3][3] = -left_eigenmatrix[2][3];
  left_eigenmatrix[3][4] = left_eigenmatrix[2][4];
  left_eigenmatrix[3][5] = left_eigenmatrix[2][5];

/*left_eigenmatrix[4][0] = 0.0; */
/*left_eigenmatrix[4][1] = 0.0; */
  left_eigenmatrix[4][2] = -left_eigenmatrix[1][2];
  left_eigenmatrix[4][3] = -left_eigenmatrix[1][3];
  left_eigenmatrix[4][4] = left_eigenmatrix[1][4];
  left_eigenmatrix[4][5] = left_eigenmatrix[1][5];

  left_eigenmatrix[5][0] = left_eigenmatrix[0][0];
  left_eigenmatrix[5][1] = -left_eigenmatrix[0][1];
  left_eigenmatrix[5][2] = -left_eigenmatrix[0][2];
  left_eigenmatrix[5][3] = -left_eigenmatrix[0][3];
  left_eigenmatrix[5][4] = left_eigenmatrix[0][4];
  left_eigenmatrix[5][5] = left_eigenmatrix[0][5];
}
#endif

/*----------------------------------------------------------------------------*/
/*! \fn void esys_prim_adb_mhd(const Real d, const Real v1, const Real rho_a2,
 *  const Real b1, const Real b2, const Real b3, Real eigenvalues[],
 *  Real right_eigenmatrix[][7], Real left_eigenmatrix[][7])
 *  \brief ADIABATIC MHD
 *   Input: d, v1, p, b1, b2, b3 = density, velocities, pressure, and B field
 *   Output: eigenvalues[7], right_eigenmatrix[7,7], left_eigenmatrix[7,7];
 */

#if !defined(BAROTROPIC) && defined(MHD)
void esys_prim_adb_mhd(const Real d, const Real v1, const Real rho_a2,
  const Real b1, const Real b2, const Real b3, Real eigenvalues[],
  Real right_eigenmatrix[][7], Real left_eigenmatrix[][7])
{
  Real di,cfsq,cf,cssq,cs,bt,bet2,bet3,alpha_f,alpha_s;
  Real sqrtd,s,a,qf,qs,af,as,vax,na,af_prime,as_prime;
  Real tsum,tdif,cf2_cs2,ct2;
  Real btsq,vaxsq,asq;
  di = 1.0/d;
  btsq  = b2*b2 + b3*b3;
  vaxsq = b1*b1*di;
  asq   = rho_a2*di;

/* Compute fast- and slow-magnetosonic speeds (eq. A10) */

  ct2 = btsq*di;
  tsum = vaxsq + ct2 + asq;
  tdif = vaxsq + ct2 - asq;
  cf2_cs2 = sqrt((double)(tdif*tdif + 4.0*asq*ct2));

  cfsq = 0.5*(tsum + cf2_cs2);
  cf = sqrt((double)cfsq);

  cssq = asq*vaxsq/cfsq;
  cs = sqrt((double)cssq);

/* Compute beta(s) (eq A17) */

  bt  = sqrt(btsq);
  if (bt == 0.0) {
    bet2 = 1.0;
    bet3 = 0.0;
  } else {
    bet2 = b2/bt;
    bet3 = b3/bt;
  }

/* Compute alpha(s) (eq A16) */

  if (cf2_cs2 == 0.0) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  } else if ( (asq - cssq) <= 0.0) {
    alpha_f = 0.0;
    alpha_s = 1.0;
  } else if ( (cfsq - asq) <= 0.0) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  } else {
    alpha_f = sqrt((asq - cssq)/cf2_cs2);
    alpha_s = sqrt((cfsq - asq)/cf2_cs2);
  }

/* Compute Q(s) and A(s) (eq. A14-15), etc. */

  sqrtd = sqrt(d);
  s = SIGN(b1);
  a = sqrt(asq);
  qf = cf*alpha_f*s;
  qs = cs*alpha_s*s;
  af = a*alpha_f*sqrtd;
  as = a*alpha_s*sqrtd;

/* Compute eigenvalues (eq. A9) */

  vax = sqrt(vaxsq);
  eigenvalues[0] = v1 - cf;
  eigenvalues[1] = v1 - vax;
  eigenvalues[2] = v1 - cs;
  eigenvalues[3] = v1;
  eigenvalues[4] = v1 + cs;
  eigenvalues[5] = v1 + vax;
  eigenvalues[6] = v1 + cf;

/* Right-eigenvectors, stored as COLUMNS (eq. A12) */
/* Note statements are grouped in ROWS for optimization, even though rem[*][n]
 * is the nth right eigenvector */


  right_eigenmatrix[0][0] = d*alpha_f;
/*right_eigenmatrix[0][1] = 0.0; */
  right_eigenmatrix[0][2] = d*alpha_s;
  right_eigenmatrix[0][3] = 1.0;
  right_eigenmatrix[0][4] = right_eigenmatrix[0][2];
/*right_eigenmatrix[0][5] = 0.0; */
  right_eigenmatrix[0][6] = right_eigenmatrix[0][0];

  right_eigenmatrix[1][0] = -cf*alpha_f;
/*right_eigenmatrix[1][1] = 0.0; */
  right_eigenmatrix[1][2] = -cs*alpha_s;
/*right_eigenmatrix[1][3] = 0.0; */
  right_eigenmatrix[1][4] = -right_eigenmatrix[1][2];
/*right_eigenmatrix[1][5] = 0.0; */
  right_eigenmatrix[1][6] = -right_eigenmatrix[1][0];

  right_eigenmatrix[2][0] = qs*bet2;
  right_eigenmatrix[2][1] = -bet3;
  right_eigenmatrix[2][2] = -qf*bet2;
/*right_eigenmatrix[2][3] = 0.0; */
  right_eigenmatrix[2][4] = -right_eigenmatrix[2][2];
  right_eigenmatrix[2][5] = bet3;
  right_eigenmatrix[2][6] = -right_eigenmatrix[2][0];

  right_eigenmatrix[3][0] = qs*bet3;
  right_eigenmatrix[3][1] = bet2;
  right_eigenmatrix[3][2] = -qf*bet3;
/*right_eigenmatrix[3][3] = 0.0; */
  right_eigenmatrix[3][4] = -right_eigenmatrix[3][2];
  right_eigenmatrix[3][5] = -bet2;
  right_eigenmatrix[3][6] = -right_eigenmatrix[3][0];

  right_eigenmatrix[4][0] = d*asq*alpha_f;
/*right_eigenmatrix[4][1] = 0.0; */
  right_eigenmatrix[4][2] = d*asq*alpha_s;
/*right_eigenmatrix[4][3] = 0.0; */
  right_eigenmatrix[4][4] = right_eigenmatrix[4][2];
/*right_eigenmatrix[4][5] = 0.0; */
  right_eigenmatrix[4][6] = right_eigenmatrix[4][0];

  right_eigenmatrix[5][0] = as*bet2;
  right_eigenmatrix[5][1] = -bet3*s*sqrtd;
  right_eigenmatrix[5][2] = -af*bet2;
/*right_eigenmatrix[5][3] = 0.0; */
  right_eigenmatrix[5][4] = right_eigenmatrix[5][2];
  right_eigenmatrix[5][5] = right_eigenmatrix[5][1];
  right_eigenmatrix[5][6] = right_eigenmatrix[5][0];

  right_eigenmatrix[6][0] = as*bet3;
  right_eigenmatrix[6][1] = bet2*s*sqrtd;
  right_eigenmatrix[6][2] = -af*bet3;
/*right_eigenmatrix[6][3] = 0.0; */
  right_eigenmatrix[6][4] = right_eigenmatrix[6][2];
  right_eigenmatrix[6][5] = right_eigenmatrix[6][1];
  right_eigenmatrix[6][6] = right_eigenmatrix[6][0];

/* Left-eigenvectors, stored as ROWS (eq. A18) */

  na = 0.5/asq;
  qf = na*qf;
  qs = na*qs;
  af_prime = na*af*di;
  as_prime = na*as*di;

/*left_eigenmatrix[0][0] = 0.0; */
  left_eigenmatrix[0][1] = -na*cf*alpha_f;
  left_eigenmatrix[0][2] = qs*bet2;
  left_eigenmatrix[0][3] = qs*bet3;
  left_eigenmatrix[0][4] = na*alpha_f*di;
  left_eigenmatrix[0][5] = as_prime*bet2;
  left_eigenmatrix[0][6] = as_prime*bet3;

/*left_eigenmatrix[1][0] = 0.0; */
/*left_eigenmatrix[1][1] = 0.0; */
  left_eigenmatrix[1][2] = -0.5*bet3;
  left_eigenmatrix[1][3] = 0.5*bet2;
/*left_eigenmatrix[1][4] = 0.0; */
  left_eigenmatrix[1][5] = -0.5*bet3*s/sqrtd;
  left_eigenmatrix[1][6] = 0.5*bet2*s/sqrtd;

/*left_eigenmatrix[2][0] = 0.0; */
  left_eigenmatrix[2][1] = -na*cs*alpha_s;
  left_eigenmatrix[2][2] = -qf*bet2;
  left_eigenmatrix[2][3] = -qf*bet3;
  left_eigenmatrix[2][4] = na*alpha_s*di;
  left_eigenmatrix[2][5] = -af_prime*bet2;
  left_eigenmatrix[2][6] = -af_prime*bet3;

  left_eigenmatrix[3][0] = 1.0;
/*left_eigenmatrix[3][1] = 0.0; */
/*left_eigenmatrix[3][2] = 0.0; */
/*left_eigenmatrix[3][3] = 0.0; */
  left_eigenmatrix[3][4] = -1.0/asq;
/*left_eigenmatrix[3][5] = 0.0; */
/*left_eigenmatrix[3][6] = 0.0; */

/*left_eigenmatrix[4][0] = 0.0; */
  left_eigenmatrix[4][1] = -left_eigenmatrix[2][1];
  left_eigenmatrix[4][2] = -left_eigenmatrix[2][2];
  left_eigenmatrix[4][3] = -left_eigenmatrix[2][3];
  left_eigenmatrix[4][4] = left_eigenmatrix[2][4];
  left_eigenmatrix[4][5] = left_eigenmatrix[2][5];
  left_eigenmatrix[4][6] = left_eigenmatrix[2][6];

/*left_eigenmatrix[5][0] = 0.0; */
/*left_eigenmatrix[5][1] = 0.0; */
  left_eigenmatrix[5][2] = -left_eigenmatrix[1][2];
  left_eigenmatrix[5][3] = -left_eigenmatrix[1][3];
/*left_eigenmatrix[5][4] = 0.0; */
  left_eigenmatrix[5][5] = left_eigenmatrix[1][5];
  left_eigenmatrix[5][6] = left_eigenmatrix[1][6];

/*left_eigenmatrix[6][0] = 0.0; */
  left_eigenmatrix[6][1] = -left_eigenmatrix[0][1];
  left_eigenmatrix[6][2] = -left_eigenmatrix[0][2];
  left_eigenmatrix[6][3] = -left_eigenmatrix[0][3];
  left_eigenmatrix[6][4] = left_eigenmatrix[0][4];
  left_eigenmatrix[6][5] = left_eigenmatrix[0][5];
  left_eigenmatrix[6][6] = left_eigenmatrix[0][6];
}
#endif
