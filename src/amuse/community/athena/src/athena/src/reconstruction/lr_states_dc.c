#include "../copyright.h"
/*============================================================================*/
/*! \file lr_states_dc.c
 *  \brief First order (donor cell, piecewise constant) spatial reconstruction.
 *
 * PURPOSE: First order (donor cell, piecewise constant) spatial reconstruction.
 * - The L/R-states at the left-interface in each cell are indexed i.
 * - W_{L,i-1/2} is denoted by Wl[i  ];   W_{R,i-1/2} is denoted by Wr[i  ]
 * - W_{L,i+1/2} is denoted by Wl[i+1];   W_{R,i+1/2} is denoted by Wr[i+1]
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - lr_states()          - computes L/R states
 * - lr_states_init()     - NoOp function in this case
 * - lr_states_destruct() - NoOp function in this case			      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef FIRST_ORDER

/*----------------------------------------------------------------------------*/
/*! \fn void lr_states(const GridS *pG, const Prim1DS W[], const Real Bxc[],
 *               const Real dt, const Real dx, const int il, const int iu,
 *               Prim1DS Wl[], Prim1DS Wr[], const enum DIRECTION dir)
 * \brief Computes L/R states 
 *
 * Input Arguments:
 * - W = PRIMITIVE variables at cell centers along 1-D slice
 * - Bxc = B in direction of slice at cell centers
 * - dtodx = dt/dx
 * - il,iu = lower and upper indices of zone centers in slice
 * W must be initialized over [il-1:iu+1]
 *
 * Output Arguments:
 * - Wl,Wr = L/R-states of PRIMITIVE variables at interfaces over [il:iu+1]
 */

void lr_states(const GridS *pG, const Prim1DS W[], const Real Bxc[],
               const Real dt, const Real dx, const int il, const int iu,
               Prim1DS Wl[], Prim1DS Wr[], const enum DIRECTION dir)
{
  int i;

  for (i=il; i<=iu+1; i++) {
    Wl[i] = W[i-1];
    Wr[i] = W[i  ];
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void lr_states_init(MeshS *pM)
 *  \brief NoOp for first order, but included for compatibility
 *   with integrator (needed for 2nd and 3rd order). */

void lr_states_init(MeshS *pM)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void lr_states_destruct(void)
 *  \brief NoOp for first order, but included for compatibility
 *   with integrator (needed for 2nd and 3rd order). */

void lr_states_destruct(void)
{
  return;
}

#endif /* FIRST_ORDER */
