#include "../copyright.h"
/*============================================================================*/
/*! \file integrate.c
 *  \brief Contains public functions to set integrator.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - integrate_init()        - set pointer to integrate function based on dim
 * - integrate_destruct()    - call destruct integrate function based on dim */
/*============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "prototypes.h"
#include "../prototypes.h"

/* dimension of calculation (determined at runtime) */
static int dim=0;

/*----------------------------------------------------------------------------*/
/*! \fn VDFun_t integrate_init(MeshS *pM)
 *  \brief Initialize integrator; VGDFun_t defined in athena.h   */
VDFun_t integrate_init(MeshS *pM)
{
  int i;
  Real cfl;
/* Calculate the dimensions (using root Domain)  */
  dim = 0;
  for (i=0; i<3; i++) if(pM->Nx[i] > 1) dim++;

/* set function pointer to appropriate integrator based on dimensions */
  switch(dim){

  case 1:
    if(pM->Nx[0] <= 1) break;
    integrate_init_1d(pM);
#if defined(CTU_INTEGRATOR)
    return integrate_1d_ctu;
#elif defined(VL_INTEGRATOR)
    cfl = par_getd("time","cour_no");
    if (cfl > 0.5)
      ath_error("<time>cour_no=%e, must be <= 0.5 with 1D VL integrator\n",cfl);
    return integrate_1d_vl;
#else
    ath_err("[integrate_init]: Invalid integrator defined for 1D problem");
#endif

  case 2:
    if(pM->Nx[2] > 1) break;
    integrate_init_2d(pM);
#if defined(CTU_INTEGRATOR)
    return integrate_2d_ctu;
#elif defined(VL_INTEGRATOR)
    cfl = par_getd("time","cour_no");
    if (cfl > 0.5)
      ath_error("<time>cour_no=%e, must be <= 0.5 with 2D VL integrator\n",cfl);
    return integrate_2d_vl;
#else
    ath_err("[integrate_init]: Invalid integrator defined for 2D problem");
#endif

  case 3:
    integrate_init_3d(pM);
#if defined(CTU_INTEGRATOR)
    cfl = par_getd("time","cour_no");
    if (cfl > 0.5)
      ath_error("<time>cour_no=%e, must be <= 0.5 with 3D CTU integrator\n",cfl);
    return integrate_3d_ctu;
#elif defined(VL_INTEGRATOR)
    cfl = par_getd("time","cour_no");
    if (cfl > 0.5)
      ath_error("<time>cour_no=%e, must be <= 0.5 with 3D VL integrator\n",cfl);
    return integrate_3d_vl;
#else
    ath_err("[integrate_init]: Invalid integrator defined for 3D problem");
#endif

  }

  if (dim == 1)
    ath_error("[integrate_init]: 1D problem must have Nx1 > 1: Nx1=%d, Nx2=%d, Nx3=%d\n",pM->Nx[0],pM->Nx[1],pM->Nx[2]);
  if (dim == 2)
     ath_error("[integrate_init]: 2D problem must have Nx1 and Nx2 > 1: Nx1=%d, Nx2=%d, Nx3=%d\n",pM->Nx[0],pM->Nx[1],pM->Nx[2]);

/* This is never executed, but generates a warning on some compilers. */
  return NULL;
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_destruct()
 *  \brief Free memory */
void integrate_destruct()
{
  switch(dim){
  case 1:
    integrate_destruct_1d();
    return;
  case 2:
    integrate_destruct_2d();
    return;
  case 3:
    integrate_destruct_3d();
    return;
  }

  ath_error("[integrate_destruct]: Grid dimension = %d\n",dim);
}
