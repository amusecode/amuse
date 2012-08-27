#include "../copyright.h"
/*============================================================================*/
/*! \file esystem_roe.c
 *  \brief Functions to evaluate the eigenvalues, and left- and
 * right-eigenvectors of "Roe's matrix A" for the linearized system in the 
 * CONSERVED variables.
 *
 * PURPOSE: Functions to evaluate the eigenvalues, and left- and
 * right-eigenvectors of "Roe's matrix A" for the linearized system in the 
 * CONSERVED variables, i.e. U,t = AU,x, where U=(d,d*vx,d*vy,d*vz,[E],[By,Bz]).
 * The eigenvalues are returned through the argument list as a vector of length
 * NWAVE.  The eigenvectors are returned as matrices of size (NWAVE)x(NWAVE),
 * with right-eigenvectors stored as COLUMNS (so R_i = right_eigenmatrix[*][i]),
 * and left-eigenvectors stored as ROWS (so L_i = left_eigenmatrix[i][*]).
 *
 * To improve performance components of the eigenvectors which are zero
 * are not set here (eigenmatrices must be initialized to zero in calling
 * routine).  However, for completeness, statements which set these values
 * are included but are commented out.
 *
 * If the input argument for the left- or right-eigenmatrix is set to the NULL
 * pointer, only the eigenvalues are returned
 *
 * The "Roe-averaging" of the L/R states must be performed in the calling funct
 *
 * REFERENCES:
 * - P. Cargo & G. Gallice, "Roe matrices for ideal MHD and systematic
 *   construction of Roe matrices for systems of conservation laws",
 *   JCP, 136, 446 (1997)
 *
 * - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new
 *   code for astrophysical MHD", ApJS, (2008), Appendix B
 *   Equation numbers refer to this paper.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - esys_roe_iso_hyd()
 * - esys_roe_adb_hyd()
 * - esys_roe_iso_mhd()
 * - esys_roe_adb_mhd()							      */
/*============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

/*----------------------------------------------------------------------------*/
/*! \fn  void esys_roe_iso_hyd(const Real v1, const Real v2, const Real v3,
 *  Real eigenvalues[],
 *  Real right_eigenmatrix[][4], Real left_eigenmatrix[][4])
 *  \brief ISOTHERMAL HYDRO
 *
 *  - Input: v1,v2,v3 = Roe averaged components of velocity
 *  - Output: eigenvalues[4],right_eigenmatrix[4,4],left_eigenmatrix[4,4]
 */

#if defined(ISOTHERMAL) && defined(HYDRO)
void esys_roe_iso_hyd(const Real v1, const Real v2, const Real v3,
  Real eigenvalues[],
  Real right_eigenmatrix[][4], Real left_eigenmatrix[][4])
{

/* Compute eigenvalues (eq. B6) */

  eigenvalues[0] = v1 - Iso_csound;
  eigenvalues[1] = v1;
  eigenvalues[2] = v1;
  eigenvalues[3] = v1 + Iso_csound;
  if (right_eigenmatrix == NULL || left_eigenmatrix == NULL) return;

/* Right-eigenvectors, stored as COLUMNS (eq. B3) */

  right_eigenmatrix[0][0] = 1.0;
  right_eigenmatrix[1][0] = v1 - Iso_csound;
  right_eigenmatrix[2][0] = v2;
  right_eigenmatrix[3][0] = v3;

/*right_eigenmatrix[0][1] = 0.0; */
/*right_eigenmatrix[1][1] = 0.0; */
  right_eigenmatrix[2][1] = 1.0;
/*right_eigenmatrix[3][1] = 0.0; */

/*right_eigenmatrix[0][2] = 0.0; */
/*right_eigenmatrix[1][2] = 0.0; */
/*right_eigenmatrix[2][2] = 0.0; */
  right_eigenmatrix[3][2] = 1.0;

  right_eigenmatrix[0][3] = 1.0;
  right_eigenmatrix[1][3] = v1 + Iso_csound;
  right_eigenmatrix[2][3] = v2;
  right_eigenmatrix[3][3] = v3;

/* Left-eigenvectors, stored as ROWS (eq. B7) */

  left_eigenmatrix[0][0] = 0.5*(1.0 + v1/Iso_csound);
  left_eigenmatrix[0][1] = -0.5/Iso_csound;
/*left_eigenmatrix[0][2] = 0.0; */
/*left_eigenmatrix[0][3] = 0.0; */

  left_eigenmatrix[1][0] = -v2;
/*left_eigenmatrix[1][1] = 0.0; */
  left_eigenmatrix[1][2] = 1.0;
/*left_eigenmatrix[1][3] = 0.0; */

  left_eigenmatrix[2][0] = -v3;
/*left_eigenmatrix[2][1] = 0.0; */
/*left_eigenmatrix[2][2] = 0.0; */
  left_eigenmatrix[2][3] = 1.0;

  left_eigenmatrix[3][0] = 0.5*(1.0 - v1/Iso_csound);
  left_eigenmatrix[3][1] = 0.5/Iso_csound;
/*left_eigenmatrix[3][2] = 0.0; */
/*left_eigenmatrix[3][3] = 0.0; */
}
#endif

/*----------------------------------------------------------------------------*/
/*! \fn void esys_roe_adb_hyd(const Real v1, const Real v2, const Real v3, 
 *			      const Real h, Real right_eigenmatrix[][5], 
 *			      Real left_eigenmatrix[][5])
 *  \brief ADIABATIC HYDRO
 *
 * - Input: v1,v2,v3,h = Roe averaged velocities and enthalpy
 * - Output: eigenvalues[5], right_eigenmatrix[5,5], left_eigenmatrix[5,5];
 */

#if defined(ADIABATIC) && defined(HYDRO)
void esys_roe_adb_hyd(const Real v1, const Real v2, const Real v3, const Real h,
  Real eigenvalues[],
  Real right_eigenmatrix[][5], Real left_eigenmatrix[][5])
{
  Real vsq,asq,a,na,qa;
  vsq = v1*v1 + v2*v2 + v3*v3;
  asq = Gamma_1*MAX((h-0.5*vsq), TINY_NUMBER);
  a = sqrt(asq);

/* Compute eigenvalues (eq. B2) */

  eigenvalues[0] = v1 - a;
  eigenvalues[1] = v1;
  eigenvalues[2] = v1;
  eigenvalues[3] = v1;
  eigenvalues[4] = v1 + a;
  if (right_eigenmatrix == NULL || left_eigenmatrix == NULL) return;

/* Right-eigenvectors, stored as COLUMNS (eq. B3) */

  right_eigenmatrix[0][0] = 1.0;
  right_eigenmatrix[1][0] = v1 - a;
  right_eigenmatrix[2][0] = v2;
  right_eigenmatrix[3][0] = v3;
  right_eigenmatrix[4][0] = h - v1*a;

/*right_eigenmatrix[0][1] = 0.0; */
/*right_eigenmatrix[1][1] = 0.0; */
  right_eigenmatrix[2][1] = 1.0;
/*right_eigenmatrix[3][1] = 0.0; */
  right_eigenmatrix[4][1] = v2;

/*right_eigenmatrix[0][2] = 0.0; */
/*right_eigenmatrix[1][2] = 0.0; */
/*right_eigenmatrix[2][2] = 0.0; */
  right_eigenmatrix[3][2] = 1.0;
  right_eigenmatrix[4][2] = v3;

  right_eigenmatrix[0][3] = 1.0;
  right_eigenmatrix[1][3] = v1;
  right_eigenmatrix[2][3] = v2;
  right_eigenmatrix[3][3] = v3;
  right_eigenmatrix[4][3] = 0.5*vsq;

  right_eigenmatrix[0][4] = 1.0;
  right_eigenmatrix[1][4] = v1 + a;
  right_eigenmatrix[2][4] = v2;
  right_eigenmatrix[3][4] = v3;
  right_eigenmatrix[4][4] = h + v1*a;

/* Left-eigenvectors, stored as ROWS (eq. B4) */

  na = 0.5/asq;
  left_eigenmatrix[0][0] = na*(0.5*Gamma_1*vsq + v1*a);
  left_eigenmatrix[0][1] = -na*(Gamma_1*v1 + a);
  left_eigenmatrix[0][2] = -na*Gamma_1*v2;
  left_eigenmatrix[0][3] = -na*Gamma_1*v3;
  left_eigenmatrix[0][4] = na*Gamma_1;

  left_eigenmatrix[1][0] = -v2;
/*left_eigenmatrix[1][1] = 0.0; */
  left_eigenmatrix[1][2] = 1.0;
/*left_eigenmatrix[1][3] = 0.0; */
/*left_eigenmatrix[1][4] = 0.0; */

  left_eigenmatrix[2][0] = -v3;
/*left_eigenmatrix[2][1] = 0.0; */
/*left_eigenmatrix[2][2] = 0.0; */
  left_eigenmatrix[2][3] = 1.0;
/*left_eigenmatrix[2][4] = 0.0; */

  qa = Gamma_1/asq;
  left_eigenmatrix[3][0] = 1.0 - na*Gamma_1*vsq;
  left_eigenmatrix[3][1] = qa*v1;
  left_eigenmatrix[3][2] = qa*v2;
  left_eigenmatrix[3][3] = qa*v3;
  left_eigenmatrix[3][4] = -qa;

  left_eigenmatrix[4][0] = na*(0.5*Gamma_1*vsq - v1*a);
  left_eigenmatrix[4][1] = -na*(Gamma_1*v1 - a);
  left_eigenmatrix[4][2] = left_eigenmatrix[0][2];
  left_eigenmatrix[4][3] = left_eigenmatrix[0][3];
  left_eigenmatrix[4][4] = left_eigenmatrix[0][4];
}
#endif

/*----------------------------------------------------------------------------*/
/*! \fn void esys_roe_iso_mhd(const Real d, const Real v1, const Real v2, 
 *			const Real v3, const Real b1, const Real b2, 
 *			const Real b3, const Real x, const Real y, 
 *                      Real eigenvalues[],
 *			Real right_eigenmatrix[][6], Real left_eigenmatrix[][6])
 *  \brief ISOTHERMAL MHD
 *
 * - Input: d,v1,v2,v3,b1,b2,b3 = Roe averaged density, velocities, and B field
 *          x,y = numerical factors (eqs. )
 * - Output: eigenvalues[6], right_eigenmatrix[6,6], left_eigenmatrix[6,6];
 */

#if defined(ISOTHERMAL) && defined(MHD)
void esys_roe_iso_mhd(const Real d, const Real v1, const Real v2, const Real v3,
  const Real b1, const Real b2, const Real b3, const Real x, const Real y, 
  Real eigenvalues[],
  Real right_eigenmatrix[][6], Real left_eigenmatrix[][6])
{
  Real btsq,bt_starsq,vaxsq,twid_csq,cfsq,cf,cssq,cs;
  Real bt,bt_star,bet2,bet3,bet2_star,bet3_star,bet_starsq,alpha_f,alpha_s;
  Real sqrtd,s,twid_c,qf,qs,af_prime,as_prime,vax;
  Real norm,cff,css,af,as,afpb,aspb,q2_star,q3_star,vqstr;
  Real ct2,tsum,tdif,cf2_cs2;
  Real di = 1.0/d;
  btsq = b2*b2 + b3*b3;
  bt_starsq = btsq*y;
  vaxsq = b1*b1*di;
  twid_csq = Iso_csound2 + x;

/* Compute fast- and slow-magnetosonic speeds (eq. B39) */

  ct2 = bt_starsq*di;
  tsum = vaxsq + ct2 + twid_csq;
  tdif = vaxsq + ct2 - twid_csq;
  cf2_cs2 = sqrt((double)(tdif*tdif + 4.0*twid_csq*ct2));

  cfsq = 0.5*(tsum + cf2_cs2);
  cf = sqrt((double)cfsq);

  cssq = twid_csq*vaxsq/cfsq;
  cs = sqrt((double)cssq);

/* Compute beta's (eqs. A17, B28, B40) */

  bt = sqrt(btsq);
  bt_star = sqrt(bt_starsq);
  if (bt == 0.0) {
    bet2 = 1.0;
    bet3 = 0.0;
  } 
  else {
    bet2 = b2/bt;
    bet3 = b3/bt;
  }
  bet2_star = bet2/sqrt(y);
  bet3_star = bet3/sqrt(y);
  bet_starsq = bet2_star*bet2_star + bet3_star*bet3_star;

/* Compute alpha's (eq. A16) */

  if ((cfsq-cssq) == 0.0) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  } else if ((twid_csq - cssq) <= 0.0) {
    alpha_f = 0.0;
    alpha_s = 1.0;
  } else if ((cfsq - twid_csq) <= 0.0) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  } else {
    alpha_f = sqrt((twid_csq - cssq)/(cfsq - cssq));
    alpha_s = sqrt((cfsq - twid_csq)/(cfsq - cssq));
  }

/* Compute Q's (eq. A14-15), etc. */

  sqrtd = sqrt(d);
  s = SIGN(b1);
  twid_c = sqrt(twid_csq);
  qf = cf*alpha_f*s;
  qs = cs*alpha_s*s;
  af_prime = twid_c*alpha_f/sqrtd;
  as_prime = twid_c*alpha_s/sqrtd;

/* Compute eigenvalues (eq. B38) */

  vax  = sqrt(vaxsq);
  eigenvalues[0] = v1 - cf;
  eigenvalues[1] = v1 - vax;
  eigenvalues[2] = v1 - cs;
  eigenvalues[3] = v1 + cs;
  eigenvalues[4] = v1 + vax;
  eigenvalues[5] = v1 + cf;
  if (right_eigenmatrix == NULL || left_eigenmatrix == NULL) return;

/* Right-eigenvectors, stored as COLUMNS (eq. B21) */

  right_eigenmatrix[0][0] = alpha_f;
  right_eigenmatrix[1][0] = alpha_f*(v1 - cf);
  right_eigenmatrix[2][0] = alpha_f*v2 + qs*bet2_star;
  right_eigenmatrix[3][0] = alpha_f*v3 + qs*bet3_star;
  right_eigenmatrix[4][0] = as_prime*bet2_star;
  right_eigenmatrix[5][0] = as_prime*bet3_star;

/*right_eigenmatrix[0][1] = 0.0; */
/*right_eigenmatrix[1][1] = 0.0; */
  right_eigenmatrix[2][1] = -bet3;
  right_eigenmatrix[3][1] = bet2;
  right_eigenmatrix[4][1] = -bet3*s/sqrtd;
  right_eigenmatrix[5][1] = bet2*s/sqrtd;

  right_eigenmatrix[0][2] = alpha_s;
  right_eigenmatrix[1][2] = alpha_s*(v1 - cs);
  right_eigenmatrix[2][2] = alpha_s*v2 - qf*bet2_star;
  right_eigenmatrix[3][2] = alpha_s*v3 - qf*bet3_star;
  right_eigenmatrix[4][2] = -af_prime*bet2_star;
  right_eigenmatrix[5][2] = -af_prime*bet3_star;

  right_eigenmatrix[0][3] = alpha_s;
  right_eigenmatrix[1][3] = alpha_s*(v1 + cs);
  right_eigenmatrix[2][3] = alpha_s*v2 + qf*bet2_star;
  right_eigenmatrix[3][3] = alpha_s*v3 + qf*bet3_star;
  right_eigenmatrix[4][3] = right_eigenmatrix[4][2];
  right_eigenmatrix[5][3] = right_eigenmatrix[5][2];

/*right_eigenmatrix[0][4] = 0.0; */
/*right_eigenmatrix[1][4] = 0.0; */
  right_eigenmatrix[2][4] = bet3;
  right_eigenmatrix[3][4] = -bet2;
  right_eigenmatrix[4][4] = right_eigenmatrix[4][1];
  right_eigenmatrix[5][4] = right_eigenmatrix[5][1];

  right_eigenmatrix[0][5] = alpha_f;
  right_eigenmatrix[1][5] = alpha_f*(v1 + cf);
  right_eigenmatrix[2][5] = alpha_f*v2 - qs*bet2_star;
  right_eigenmatrix[3][5] = alpha_f*v3 - qs*bet3_star;
  right_eigenmatrix[4][5] = right_eigenmatrix[4][0];
  right_eigenmatrix[5][5] = right_eigenmatrix[5][0];

/* Left-eigenvectors, stored as ROWS (eq. B41) */
/* Normalize by 1/2a^{2}: quantities denoted by \hat{f} */

  norm = 0.5/twid_csq;
  cff = norm*alpha_f*cf;
  css = norm*alpha_s*cs;
  qf *= norm;
  qs *= norm;
  af = norm*af_prime*d;
  as = norm*as_prime*d;
  afpb = norm*af_prime*bt_star;
  aspb = norm*as_prime*bt_star;

  q2_star = bet2_star/bet_starsq;
  q3_star = bet3_star/bet_starsq;
  vqstr = (v2*q2_star + v3*q3_star);

  left_eigenmatrix[0][0] = cff*(cf+v1) - qs*vqstr - aspb;
  left_eigenmatrix[0][1] = -cff;
  left_eigenmatrix[0][2] = qs*q2_star;
  left_eigenmatrix[0][3] = qs*q3_star;
  left_eigenmatrix[0][4] = as*q2_star;
  left_eigenmatrix[0][5] = as*q3_star;

  left_eigenmatrix[1][0] = 0.5*(v2*bet3 - v3*bet2);
/*left_eigenmatrix[1][1] = 0.0; */
  left_eigenmatrix[1][2] = -0.5*bet3;
  left_eigenmatrix[1][3] = 0.5*bet2;
  left_eigenmatrix[1][4] = -0.5*sqrtd*bet3*s;
  left_eigenmatrix[1][5] = 0.5*sqrtd*bet2*s;

  left_eigenmatrix[2][0] = css*(cs+v1) + qf*vqstr + afpb;
  left_eigenmatrix[2][1] = -css;
  left_eigenmatrix[2][2] = -qf*q2_star;
  left_eigenmatrix[2][3] = -qf*q3_star;
  left_eigenmatrix[2][4] = -af*q2_star;
  left_eigenmatrix[2][5] = -af*q3_star;

  left_eigenmatrix[3][0] = css*(cs-v1) - qf*vqstr + afpb;
  left_eigenmatrix[3][1] = css;
  left_eigenmatrix[3][2] = -left_eigenmatrix[2][2];
  left_eigenmatrix[3][3] = -left_eigenmatrix[2][3];
  left_eigenmatrix[3][4] = left_eigenmatrix[2][4];
  left_eigenmatrix[3][5] = left_eigenmatrix[2][5];

  left_eigenmatrix[4][0] = -left_eigenmatrix[1][0];
/*left_eigenmatrix[4][1] = 0.0; */
  left_eigenmatrix[4][2] = -left_eigenmatrix[1][2];
  left_eigenmatrix[4][3] = -left_eigenmatrix[1][3];
  left_eigenmatrix[4][4] = left_eigenmatrix[1][4];
  left_eigenmatrix[4][5] = left_eigenmatrix[1][5];

  left_eigenmatrix[5][0] = cff*(cf-v1) + qs*vqstr - aspb;
  left_eigenmatrix[5][1] = cff;
  left_eigenmatrix[5][2] = -left_eigenmatrix[0][2];
  left_eigenmatrix[5][3] = -left_eigenmatrix[0][3];
  left_eigenmatrix[5][4] = left_eigenmatrix[0][4];
  left_eigenmatrix[5][5] = left_eigenmatrix[0][5];
}
#endif

/*----------------------------------------------------------------------------*/
/*! \fn void esys_roe_adb_mhd(const Real d, const Real v1, const Real v2, 
 *			const Real v3, const Real h, const Real b1, 
 *			const Real b2, const Real b3, Real eigenvalues[], 
 *			Real right_eigenmatrix[][7], Real left_eigenmatrix[][7])
 *  \brief ADIABATIC MHD
 *
 * - Input: d,v1,v2,v3,h,b1,b2,b3=Roe averaged density, velocities, enthalpy, B
 *          x,y = numerical factors (see eqn XX)
 * - Output: eigenvalues[7], right_eigenmatrix[7,7], left_eigenmatrix[7,7];
 */

#if defined(ADIABATIC) && defined(MHD)
void esys_roe_adb_mhd(const Real d, const Real v1, const Real v2, const Real v3,
  const Real h, const Real b1, const Real b2, const Real b3, 
  const Real x, const Real y,
  Real eigenvalues[],
  Real right_eigenmatrix[][7], Real left_eigenmatrix[][7])
{
  Real di,vsq,btsq,bt_starsq,vaxsq,hp,twid_asq,cfsq,cf,cssq,cs;
  Real bt,bt_star,bet2,bet3,bet2_star,bet3_star,bet_starsq,vbet,alpha_f,alpha_s;
  Real isqrtd,sqrtd,s,twid_a,qf,qs,af_prime,as_prime,afpbb,aspbb,vax;
  Real norm,cff,css,af,as,afpb,aspb,q2_star,q3_star,vqstr;
  Real ct2,tsum,tdif,cf2_cs2;
  Real qa,qb,qc,qd;
  di = 1.0/d;
  vsq = v1*v1 + v2*v2 + v3*v3;
  btsq = b2*b2 + b3*b3;
  bt_starsq = (Gamma_1 - Gamma_2*y)*btsq;
  vaxsq = b1*b1*di;
  hp = h - (vaxsq + btsq*di);
  twid_asq = MAX((Gamma_1*(hp-0.5*vsq)-Gamma_2*x), TINY_NUMBER);

/* Compute fast- and slow-magnetosonic speeds (eq. B18) */

  ct2 = bt_starsq*di;
  tsum = vaxsq + ct2 + twid_asq;
  tdif = vaxsq + ct2 - twid_asq;
  cf2_cs2 = sqrt((double)(tdif*tdif + 4.0*twid_asq*ct2));

  cfsq = 0.5*(tsum + cf2_cs2);
  cf = sqrt((double)cfsq);

  cssq = twid_asq*vaxsq/cfsq;
  cs = sqrt((double)cssq);

/* Compute beta(s) (eqs. A17, B20, B28) */

  bt = sqrt(btsq);
  bt_star = sqrt(bt_starsq);
  if (bt == 0.0) {
    bet2 = 1.0;
    bet3 = 0.0;
  } else {
    bet2 = b2/bt;
    bet3 = b3/bt;
  }
  bet2_star = bet2/sqrt(Gamma_1 - Gamma_2*y);
  bet3_star = bet3/sqrt(Gamma_1 - Gamma_2*y);
  bet_starsq = bet2_star*bet2_star + bet3_star*bet3_star;
  vbet = v2*bet2_star + v3*bet3_star;

/* Compute alpha(s) (eq. A16) */

  if ((cfsq-cssq) == 0.0) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  } else if ( (twid_asq - cssq) <= 0.0) {
    alpha_f = 0.0;
    alpha_s = 1.0;
  } else if ( (cfsq - twid_asq) <= 0.0) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  } else {
    alpha_f = sqrt((twid_asq - cssq)/(cfsq - cssq));
    alpha_s = sqrt((cfsq - twid_asq)/(cfsq - cssq));
  }

/* Compute Q(s) and A(s) (eq. A14-15), etc. */

  sqrtd = sqrt(d);
  isqrtd = 1.0/sqrtd;
  s = SIGN(b1);
  twid_a = sqrt(twid_asq);
  qf = cf*alpha_f*s;
  qs = cs*alpha_s*s;
  af_prime = twid_a*alpha_f*isqrtd;
  as_prime = twid_a*alpha_s*isqrtd;
  afpbb = af_prime*bt_star*bet_starsq;
  aspbb = as_prime*bt_star*bet_starsq;

/* Compute eigenvalues (eq. B17) */

  vax = sqrt(vaxsq);
  eigenvalues[0] = v1 - cf;
  eigenvalues[1] = v1 - vax;
  eigenvalues[2] = v1 - cs;
  eigenvalues[3] = v1;
  eigenvalues[4] = v1 + cs;
  eigenvalues[5] = v1 + vax;
  eigenvalues[6] = v1 + cf;
  if (right_eigenmatrix == NULL || left_eigenmatrix == NULL) return;

/* Right-eigenvectors, stored as COLUMNS (eq. B21) */
/* Note statements are grouped in ROWS for optimization, even though rem[*][n]
 * is the nth right eigenvector */

  right_eigenmatrix[0][0] = alpha_f;
/*right_eigenmatrix[0][1] = 0.0; */
  right_eigenmatrix[0][2] = alpha_s;
  right_eigenmatrix[0][3] = 1.0;
  right_eigenmatrix[0][4] = alpha_s;
/*right_eigenmatrix[0][5] = 0.0; */
  right_eigenmatrix[0][6] = alpha_f;

  right_eigenmatrix[1][0] = alpha_f*eigenvalues[0];
/*right_eigenmatrix[1][1] = 0.0; */
  right_eigenmatrix[1][2] = alpha_s*eigenvalues[2];
  right_eigenmatrix[1][3] = v1;
  right_eigenmatrix[1][4] = alpha_s*eigenvalues[4];
/*right_eigenmatrix[1][5] = 0.0; */
  right_eigenmatrix[1][6] = alpha_f*eigenvalues[6];

  qa = alpha_f*v2;
  qb = alpha_s*v2;
  qc = qs*bet2_star;
  qd = qf*bet2_star;
  right_eigenmatrix[2][0] = qa + qc;
  right_eigenmatrix[2][1] = -bet3;
  right_eigenmatrix[2][2] = qb - qd;
  right_eigenmatrix[2][3] = v2;
  right_eigenmatrix[2][4] = qb + qd;
  right_eigenmatrix[2][5] = bet3;
  right_eigenmatrix[2][6] = qa - qc;

  qa = alpha_f*v3;
  qb = alpha_s*v3;
  qc = qs*bet3_star;
  qd = qf*bet3_star;
  right_eigenmatrix[3][0] = qa + qc;
  right_eigenmatrix[3][1] = bet2;
  right_eigenmatrix[3][2] = qb - qd;
  right_eigenmatrix[3][3] = v3;
  right_eigenmatrix[3][4] = qb + qd;
  right_eigenmatrix[3][5] = -bet2;
  right_eigenmatrix[3][6] = qa - qc;

  right_eigenmatrix[4][0] = alpha_f*(hp - v1*cf) + qs*vbet + aspbb;
  right_eigenmatrix[4][1] = -(v2*bet3 - v3*bet2);
  right_eigenmatrix[4][2] = alpha_s*(hp - v1*cs) - qf*vbet - afpbb;
  right_eigenmatrix[4][3] = 0.5*vsq + Gamma_2*x/Gamma_1;
  right_eigenmatrix[4][4] = alpha_s*(hp + v1*cs) + qf*vbet - afpbb;
  right_eigenmatrix[4][5] = -right_eigenmatrix[4][1];
  right_eigenmatrix[4][6] = alpha_f*(hp + v1*cf) - qs*vbet + aspbb;

  right_eigenmatrix[5][0] = as_prime*bet2_star;
  right_eigenmatrix[5][1] = -bet3*s*isqrtd;
  right_eigenmatrix[5][2] = -af_prime*bet2_star;
/*right_eigenmatrix[5][3] = 0.0; */
  right_eigenmatrix[5][4] = right_eigenmatrix[5][2];
  right_eigenmatrix[5][5] = right_eigenmatrix[5][1];
  right_eigenmatrix[5][6] = right_eigenmatrix[5][0];

  right_eigenmatrix[6][0] = as_prime*bet3_star;
  right_eigenmatrix[6][1] = bet2*s*isqrtd;
  right_eigenmatrix[6][2] = -af_prime*bet3_star;
/*right_eigenmatrix[6][3] = 0.0; */
  right_eigenmatrix[6][4] = right_eigenmatrix[6][2];
  right_eigenmatrix[6][5] = right_eigenmatrix[6][1];
  right_eigenmatrix[6][6] = right_eigenmatrix[6][0];

/* Left-eigenvectors, stored as ROWS (eq. B29) */

/* Normalize by 1/2a^{2}: quantities denoted by \hat{f} */
  norm = 0.5/twid_asq;
  cff = norm*alpha_f*cf;
  css = norm*alpha_s*cs;
  qf *= norm;
  qs *= norm;
  af = norm*af_prime*d;
  as = norm*as_prime*d;
  afpb = norm*af_prime*bt_star;
  aspb = norm*as_prime*bt_star;

/* Normalize by (gamma-1)/2a^{2}: quantities denoted by \bar{f} */
  norm *= Gamma_1;
  alpha_f *= norm;
  alpha_s *= norm;
  q2_star = bet2_star/bet_starsq;
  q3_star = bet3_star/bet_starsq;
  vqstr = (v2*q2_star + v3*q3_star);
  norm *= 2.0;

  left_eigenmatrix[0][0] = alpha_f*(vsq-hp) + cff*(cf+v1) - qs*vqstr - aspb;
  left_eigenmatrix[0][1] = -alpha_f*v1 - cff;
  left_eigenmatrix[0][2] = -alpha_f*v2 + qs*q2_star;
  left_eigenmatrix[0][3] = -alpha_f*v3 + qs*q3_star;
  left_eigenmatrix[0][4] = alpha_f;
  left_eigenmatrix[0][5] = as*q2_star - alpha_f*b2;
  left_eigenmatrix[0][6] = as*q3_star - alpha_f*b3;

  left_eigenmatrix[1][0] = 0.5*(v2*bet3 - v3*bet2);
/*left_eigenmatrix[1][1] = 0.0; */
  left_eigenmatrix[1][2] = -0.5*bet3;
  left_eigenmatrix[1][3] = 0.5*bet2;
/*left_eigenmatrix[1][4] = 0.0; */
  left_eigenmatrix[1][5] = -0.5*sqrtd*bet3*s;
  left_eigenmatrix[1][6] = 0.5*sqrtd*bet2*s;

  left_eigenmatrix[2][0] = alpha_s*(vsq-hp) + css*(cs+v1) + qf*vqstr + afpb;
  left_eigenmatrix[2][1] = -alpha_s*v1 - css;
  left_eigenmatrix[2][2] = -alpha_s*v2 - qf*q2_star;
  left_eigenmatrix[2][3] = -alpha_s*v3 - qf*q3_star;
  left_eigenmatrix[2][4] = alpha_s;
  left_eigenmatrix[2][5] = -af*q2_star - alpha_s*b2;
  left_eigenmatrix[2][6] = -af*q3_star - alpha_s*b3;

  left_eigenmatrix[3][0] = 1.0 - norm*(0.5*vsq - Gamma_2*x/Gamma_1); 
  left_eigenmatrix[3][1] = norm*v1;
  left_eigenmatrix[3][2] = norm*v2;
  left_eigenmatrix[3][3] = norm*v3;
  left_eigenmatrix[3][4] = -norm;
  left_eigenmatrix[3][5] = norm*b2;
  left_eigenmatrix[3][6] = norm*b3;

  left_eigenmatrix[4][0] = alpha_s*(vsq-hp) + css*(cs-v1) - qf*vqstr + afpb;
  left_eigenmatrix[4][1] = -alpha_s*v1 + css;
  left_eigenmatrix[4][2] = -alpha_s*v2 + qf*q2_star;
  left_eigenmatrix[4][3] = -alpha_s*v3 + qf*q3_star;
  left_eigenmatrix[4][4] = alpha_s;
  left_eigenmatrix[4][5] = left_eigenmatrix[2][5];
  left_eigenmatrix[4][6] = left_eigenmatrix[2][6];

  left_eigenmatrix[5][0] = -left_eigenmatrix[1][0];
/*left_eigenmatrix[5][1] = 0.0; */
  left_eigenmatrix[5][2] = -left_eigenmatrix[1][2];
  left_eigenmatrix[5][3] = -left_eigenmatrix[1][3];
/*left_eigenmatrix[5][4] = 0.0; */
  left_eigenmatrix[5][5] = left_eigenmatrix[1][5];
  left_eigenmatrix[5][6] = left_eigenmatrix[1][6];

  left_eigenmatrix[6][0] = alpha_f*(vsq-hp) + cff*(cf-v1) + qs*vqstr - aspb;
  left_eigenmatrix[6][1] = -alpha_f*v1 + cff;
  left_eigenmatrix[6][2] = -alpha_f*v2 - qs*q2_star;
  left_eigenmatrix[6][3] = -alpha_f*v3 - qs*q3_star;
  left_eigenmatrix[6][4] = alpha_f;
  left_eigenmatrix[6][5] = left_eigenmatrix[0][5];
  left_eigenmatrix[6][6] = left_eigenmatrix[0][6];
}
#endif
