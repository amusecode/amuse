#ifndef _applyBC_h_
#define _applyBC_h_

#include "vader_common.h"

/**************************************************************************/
/* Boundary conditions filler routine                                     */
/**************************************************************************/

void 
applyBC(
	/* Grid structure and viscosity values on grid */
	const grid *grd, const double *alpha_g,
	/* Inner boundary condition specifications */
	const pres_bc_type ibc_pres, const enth_bc_type ibc_enth,
	const double ibc_pres_val, const double ibc_enth_val,
	/* Outer boundary condition specifications */
	const pres_bc_type obc_pres, const enth_bc_type obc_enth,
	const double obc_pres_val, const double obc_enth_val,
	/* Outputs */
	double *pres_g, double *hint_g, double *pIn, double *qIn, 
	double *pOut, double *qOut
	);
#endif
