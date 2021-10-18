#include <math.h>
#include "applyBC.h"

/**************************************************************************/
/* Boundary conditions filler routine                                     */
/**************************************************************************/

void applyBC(const grid *grd, const double *alpha_g,
	     const pres_bc_type ibc_pres, const enth_bc_type ibc_enth,
	     const double ibc_pres_val, const double ibc_enth_val,
	     const pres_bc_type obc_pres, const enth_bc_type obc_enth,
	     const double obc_pres_val, const double obc_enth_val,
	     double *pres_g, double *hint_g, double *pIn, double *qIn, 
	     double *pOut, double *qOut) {

  /* Inner pressure BC */
  if (ibc_pres == FIXED_MASS_FLUX) {
    *pIn = ibc_pres_val / (grd->g_h[0]*alpha_g[0]*
			   (1-grd->beta_g[0])*SQR(grd->r_g[0]));
    *qIn = SQR(grd->r_g[1]/grd->r_g[0])*
      (1.0-grd->beta_g[1])/(1.0-grd->beta_g[0]);
  } else if (ibc_pres == FIXED_TORQUE_FLUX) {
    *qIn = -1.0;
    *pIn = ibc_pres_val / 
      (M_PI*alpha_g[0]*grd->r_h[0]*grd->vphi_h[0]*(1-grd->beta_h[0]));
  } else if (ibc_pres == FIXED_TORQUE) {
    *qIn = 0.0;
    *pIn = -ibc_pres_val / 
      (2.0*M_PI*SQR(grd->r_g[0])*alpha_g[0]*(1.0-grd->beta_g[0]));
  } else {
    fprintf(stderr, "Error: unknown inner pressure boundary condition\n");
    exit(1);
  }
  pres_g[0] = (*qIn) * pres_g[1] + (*pIn);

  /* Outer pressure BC */
  if (obc_pres == FIXED_MASS_FLUX) {
    *qOut = SQR(grd->r_g[grd->nr] / grd->r_g[grd->nr+1]) *
      (1.0-grd->beta_g[grd->nr]) / (1.0-grd->beta_g[grd->nr+1]);
    *pOut = -obc_pres_val
      / (grd->g_h[grd->nr]*alpha_g[grd->nr+1]*
	 (1-grd->beta_g[grd->nr+1])*SQR(grd->r_g[grd->nr+1]));
  } else if (obc_pres == FIXED_TORQUE_FLUX) {
    *qOut = -1.0;
    *pOut = obc_pres_val /
      (M_PI*alpha_g[grd->nr+1]*grd->r_h[grd->nr]*grd->vphi_h[grd->nr]
       *(1-grd->beta_h[grd->nr]));
  } else if (obc_pres == FIXED_TORQUE) {
    *qOut = 0.0;
    *pOut = -obc_pres_val /
      (2.0*M_PI*SQR(grd->r_g[grd->nr+1])*alpha_g[grd->nr+1]
       *(1.0-grd->beta_g[grd->nr+1]));
  } else {
    fprintf(stderr, "Error: unknown outer pressure boundary condition\n");
    exit(1);
  }
  pres_g[grd->nr+1] = (*qOut) * pres_g[grd->nr] + (*pOut);

  /* Inner enthalpy BC */
  if (ibc_enth == FIXED_ENTHALPY_VALUE) {
    hint_g[0] = ibc_enth_val;
  } else if (ibc_enth == FIXED_ENTHALPY_GRADIENT) {
    hint_g[0] = hint_g[1] - ibc_enth_val *
      (grd->r_g[1]-grd->r_g[0]);
  } else {
    fprintf(stderr, "Error: unknown inner enthalpy boundary condition\n");
    exit(1);
  }

  /* Outer enthalpy BC */
  if (obc_enth == FIXED_ENTHALPY_VALUE) {
    hint_g[grd->nr+1] = obc_enth_val;
  } else if (obc_enth == FIXED_ENTHALPY_GRADIENT) {
    hint_g[grd->nr+1] = hint_g[grd->nr] + obc_enth_val *
      (grd->r_g[grd->nr+1]-grd->r_g[grd->nr]);
  } else {
    fprintf(stderr, "Error: unknown outer enthalpy boundary condition\n");
    exit(1);
  }
}
