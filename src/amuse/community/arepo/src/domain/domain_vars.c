/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/domain_vars.c
 * \date        05/2018
 * \brief       Variables and memory allocation functions for domain
 *              decomposition.
 * \details     contains functions:
 *                void domain_allocate_lists(void)
 *                void domain_free_lists(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 05.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../mesh/voronoi/voronoi.h"
#include "domain.h"

struct domain_peano_hilbert_data *mp;

struct local_topnode_data *topNodes, *branchNodes; /*!< points to the root node of the top-level tree */

double totgravcost, totpartcount, gravcost, totsphcost, sphcost;

struct domain_cost_data *DomainLeaveNode;

double fac_work, fac_load, fac_worksph;
double normsum_work, normsum_load, normsum_worksph;

int Nbranch;

/*! toGo[partner] gives the number of particles on the current task that have to go to task 'partner'
 */
int *toGo, *toGoSph;
int *toGet, *toGetSph;
int *list_NumPart;
int *list_NumGas;
int *list_load;
int *list_loadsph;
double *list_work;
double *list_worksph;

/*! \brief Allocates lists needed for communication in domain decomposition.
 *
 *  These lists are holding information about other tasks (number of particles,
 *  load, work, etc.).
 *
 *  \return void
 */
void domain_allocate_lists(void)
{
  Key             = (peanokey *)mymalloc_movable(&Key, "domain_key", (sizeof(peanokey) * All.MaxPart));
  toGo            = (int *)mymalloc_movable(&toGo, "toGo", (sizeof(int) * NTask));
  toGoSph         = (int *)mymalloc_movable(&toGoSph, "toGoSph", (sizeof(int) * NTask));
  toGet           = (int *)mymalloc_movable(&toGet, "toGet", (sizeof(int) * NTask));
  toGetSph        = (int *)mymalloc_movable(&toGetSph, "toGetSph", (sizeof(int) * NTask));
  list_NumPart    = (int *)mymalloc_movable(&list_NumPart, "list_NumPart", (sizeof(int) * NTask));
  list_NumGas     = (int *)mymalloc_movable(&list_NumGas, "list_NumGas", (sizeof(int) * NTask));
  list_load       = (int *)mymalloc_movable(&list_load, "list_load", (sizeof(int) * NTask));
  list_loadsph    = (int *)mymalloc_movable(&list_loadsph, "list_loadsph", (sizeof(int) * NTask));
  list_work       = (double *)mymalloc_movable(&list_work, "list_work", (sizeof(double) * NTask));
  list_worksph    = (double *)mymalloc_movable(&list_worksph, "list_worksph", (sizeof(double) * NTask));
  DomainLeaveNode = (struct domain_cost_data *)mymalloc_movable(&DomainLeaveNode, "DomainLeaveNode",
                                                                (MaxTopNodes * sizeof(struct domain_cost_data)));
}

/*! \brief Frees lists needed for communication in domain decomposition.
 *
 *  This routine is the counterpart of domain_allocate_lists(void).
 *  Frees memory of all arrays allocated there, except Key, which is freed
 *  elsewhere (in void domain_Decomposition(void); see domain.c).
 *
 *  \return void
 */
void domain_free_lists(void)
{
  myfree(DomainLeaveNode);
  myfree(list_worksph);
  myfree(list_work);
  myfree(list_loadsph);
  myfree(list_load);
  myfree(list_NumGas);
  myfree(list_NumPart);
  myfree(toGetSph);
  myfree(toGet);
  myfree(toGoSph);
  myfree(toGo);
}
