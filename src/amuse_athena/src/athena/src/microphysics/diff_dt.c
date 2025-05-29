#include "../copyright.h"
/*============================================================================*/
/*! \file diff_dt.c
 *  \brief Computes diffusion timestep using CFL condition, for all diffusive
 *   processes currently implemented in code. 
 *
 *  These include:
 *     - Ohmic dissipation, Hall effect, ambipolar diffusion
 *     - Navier-Stokes and Braginskii viscosity
 *     - isotropic and anisotropic thermal conduction
 *   With MPI parallel jobs, finds minimum dt across all processors.
 *   Function returns minimum diffusion dt.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - diff_dt()  - computes dt */
/*============================================================================*/

#include <stdio.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

/*----------------------------------------------------------------------------*/
/*! \fn Real diff_dt(MeshS *pM)
 *  \brief Computes diffusion timestep */
Real diff_dt(MeshS *pM)
{
  int irefine, ir;
  Real dtmin_diffusion=(HUGE_NUMBER);
  Real dxmin,qa,fac;
#ifdef RESISTIVITY
  int i,j,k,nl,nd;
  int il,iu,jl,ju,kl,ku;
  Real qb;
  GridS *pG;
#endif
#ifdef MPI_PARALLEL
  double my_dt, dt;
  int ierr;
#endif

/* Calculate minimum dx.  Always given by Grid on highest level of refinement */

  irefine = 1;
  for (ir=1; ir<(pM->NLevels); ir++) irefine *= 2;

  dxmin = pM->dx[0]/(Real)(irefine);
  if (pM->Nx[1] > 1) dxmin = MIN( dxmin, (pM->dx[1]/(Real)(irefine)) );
  if (pM->Nx[2] > 1) dxmin = MIN( dxmin, (pM->dx[2]/(Real)(irefine)) );

  qa = CourNo*(dxmin*dxmin)/2.0;

  fac = 1.0;
  if (pM->Nx[1] > 1) fac = 2.0;
  if (pM->Nx[2] > 1) fac = 3.0;

  qa = qa / fac;

#ifdef THERMAL_CONDUCTION
  dtmin_diffusion = MIN(dtmin_diffusion,(qa/(kappa_iso + kappa_aniso)));
#endif
#ifdef VISCOSITY
  dtmin_diffusion = MIN(dtmin_diffusion,(qa/(nu_iso + nu_aniso)));
#endif

#ifdef RESISTIVITY
  qb = 0.25*qa;
  for (nl=pM->NLevels-1; nl>=0; nl--){
    qb *= 4.0;
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){

        pG = pM->Domain[nl][nd].Grid;

        il = pG->is - 2;    iu = pG->ie + 2;
        if (pG->Nx[1] > 1){
          jl = pG->js - 2;  ju = pG->je + 2;
        } else {
          jl = pG->js;      ju = pG->je;
        }
        if (pG->Nx[2] > 1){
          kl = pG->ks - 2;  ku = pG->ke + 2;
        } else {
          kl = pG->ks;      ku = pG->ke;
        }

        for (k=kl; k<=ku; k++) {
        for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {

          if (eta_Ohm > 0.0){
            dtmin_diffusion = MIN(dtmin_diffusion,(qb/pG->eta_Ohm[k][j][i]));
          }
          if (Q_Hall > 0.0)
            dtmin_diffusion = MIN(dtmin_diffusion,
                                         (0.5*fac*qb/pG->eta_Hall[k][j][i]));

          if (Q_AD > 0.0)
            dtmin_diffusion = MIN(dtmin_diffusion, (qb/pG->eta_AD[k][j][i]));
        }}}
      }
    }
  }
#endif

/* Find minimum timestep over all processors */
#ifdef MPI_PARALLEL
  my_dt = dtmin_diffusion;
  ierr = MPI_Allreduce(&my_dt, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  dtmin_diffusion = dt;
#endif /* MPI_PARALLEL */

  return dtmin_diffusion;
}
