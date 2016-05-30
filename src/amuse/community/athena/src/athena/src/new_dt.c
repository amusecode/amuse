#include "copyright.h"
/*============================================================================*/
/*! \file new_dt.c
 *  \brief Computes timestep using CFL condition.
 *
 * PURPOSE: Computes timestep using CFL condition on cell-centered velocities
 *   and sound speed, and Alfven speed from face-centered B, across all Grids
 *   being updated on this processor.  With MPI parallel jobs, also finds
 *   minimum dt across all processors.
 *
 * For special relativity, the time step limit is just (1/dx), since the fastest
 * wave speed is never larger than c=1.
 *
 * A CFL condition is also applied using particle velocities if PARTICLES is
 * defined.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - new_dt() - computes dt						      */
/*============================================================================*/

#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/*! \fn void new_dt(MeshS *pM)
 *  \brief Computes timestep using CFL condition. */ 

void new_dt(MeshS *pM)
{
  GridS *pGrid;
#ifndef SPECIAL_RELATIVITY
  int i,j,k;
  Real di,v1,v2,v3,qsq,asq,cf1sq,cf2sq,cf3sq;
#ifdef ADIABATIC
  Real p;
#endif
#ifdef MHD
  Real b1,b2,b3,bsq,tsum,tdif;
#endif /* MHD */
#ifdef PARTICLES
  long q;
#endif /* PARTICLES */
#endif /* SPECIAL RELATIVITY */
#ifdef MPI_PARALLEL
  double dt, my_dt;
  int ierr;
#endif
  int nl,nd;
  Real tlim,max_v1=0.0,max_v2=0.0,max_v3=0.0,max_dti = 0.0;
  Real x1,x2,x3;

/* Loop over all Domains with a Grid on this processor -----------------------*/

  for (nl=0; nl<(pM->NLevels); nl++){
  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){

  if (pM->Domain[nl][nd].Grid != NULL) {
    pGrid=(pM->Domain[nl][nd].Grid);

/* Maximum velocity is always c with special relativity */
#ifdef SPECIAL_RELATIVITY
    max_v1 = max_v2 = max_v3 = 1.0;
#else

    for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {
        di = 1.0/(pGrid->U[k][j][i].d);
        v1 = pGrid->U[k][j][i].M1*di;
        v2 = pGrid->U[k][j][i].M2*di;
        v3 = pGrid->U[k][j][i].M3*di;
        qsq = v1*v1 + v2*v2 + v3*v3;

#ifdef MHD

/* Use maximum of face-centered fields (always larger than cell-centered B) */
        b1 = pGrid->U[k][j][i].B1c 
          + fabs((double)(pGrid->B1i[k][j][i] - pGrid->U[k][j][i].B1c));
        b2 = pGrid->U[k][j][i].B2c 
          + fabs((double)(pGrid->B2i[k][j][i] - pGrid->U[k][j][i].B2c));
        b3 = pGrid->U[k][j][i].B3c 
          + fabs((double)(pGrid->B3i[k][j][i] - pGrid->U[k][j][i].B3c));
        bsq = b1*b1 + b2*b2 + b3*b3;
/* compute sound speed squared */
#ifdef ADIABATIC
        p = MAX(Gamma_1*(pGrid->U[k][j][i].E - 0.5*pGrid->U[k][j][i].d*qsq
                - 0.5*bsq), TINY_NUMBER);
        asq = Gamma*p*di;
#elif defined ISOTHERMAL
        asq = Iso_csound2;
#endif /* EOS */
/* compute fast magnetosonic speed squared in each direction */
        tsum = bsq*di + asq;
        tdif = bsq*di - asq;
        cf1sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b2*b2+b3*b3)*di));
        cf2sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b1*b1+b3*b3)*di));
        cf3sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b1*b1+b2*b2)*di));

#else /* MHD */

/* compute sound speed squared */
#ifdef ADIABATIC
        p = MAX(Gamma_1*(pGrid->U[k][j][i].E - 0.5*pGrid->U[k][j][i].d*qsq),
                TINY_NUMBER);
        asq = Gamma*p*di;
#elif defined ISOTHERMAL
        asq = Iso_csound2;
#endif /* EOS */
/* compute fast magnetosonic speed squared in each direction */
        cf1sq = asq;
        cf2sq = asq;
        cf3sq = asq;

#endif /* MHD */

/* compute maximum cfl velocity (corresponding to minimum dt) */
        if (pGrid->Nx[0] > 1)
          max_v1 = MAX(max_v1,fabs(v1)+sqrt((double)cf1sq));
        if (pGrid->Nx[1] > 1)
#ifdef CYLINDRICAL
          cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
          max_v2 = MAX(max_v2,(fabs(v2)+sqrt((double)cf2sq))/x1);
#else
          max_v2 = MAX(max_v2,fabs(v2)+sqrt((double)cf2sq));
#endif
        if (pGrid->Nx[2] > 1)
          max_v3 = MAX(max_v3,fabs(v3)+sqrt((double)cf3sq));
  
      }
    }}

#endif /* SPECIAL_RELATIVITY */

/* compute maximum velocity with particles */
#ifdef PARTICLES
    for (q=0; q<pGrid->nparticle; q++) {
      if (pGrid->Nx[0] > 1)
        max_v1 = MAX(max_v1, pGrid->particle[q].v1);
      if (pGrid->Nx[1] > 1)
        max_v2 = MAX(max_v2, pGrid->particle[q].v2);
      if (pGrid->Nx[2] > 1)
        max_v3 = MAX(max_v3, pGrid->particle[q].v3);
    }
#endif /* PARTICLES */

/* compute maximum inverse of dt (corresponding to minimum dt) */
    if (pGrid->Nx[0] > 1)
      max_dti = MAX(max_dti, max_v1/pGrid->dx1);
    if (pGrid->Nx[1] > 1)
      max_dti = MAX(max_dti, max_v2/pGrid->dx2);
    if (pGrid->Nx[2] > 1)
      max_dti = MAX(max_dti, max_v3/pGrid->dx3);

  }}} /*--- End loop over Domains --------------------------------------------*/

/* new timestep.  Limit increase to 2x old value */

  if (pM->nstep == 0) {
    pM->dt = CourNo/max_dti;
  } else {
    pM->dt = MIN(2.0*pM->dt, CourNo/max_dti);
  }

/* Find minimum timestep over all processors */

#ifdef MPI_PARALLEL
  my_dt = pM->dt;
  ierr = MPI_Allreduce(&my_dt, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  pM->dt = dt;
#endif /* MPI_PARALLEL */

/* modify timestep so loop finishes at t=tlim exactly */

  tlim = par_getd("time","tlim");
  if ((pM->time < tlim) && ((tlim - pM->time) < pM->dt))
    pM->dt = tlim - pM->time;

/* Spread timestep across all Grid structures in all Domains */

  for (nl=0; nl<=(pM->NLevels)-1; nl++){
    for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        pM->Domain[nl][nd].Grid->dt = pM->dt;
      }
    }
  }

  return;
}
