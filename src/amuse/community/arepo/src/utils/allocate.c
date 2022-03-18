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
 * \file        src/utils/allocate.c
 * \date        05/2018
 * \brief       Functions to allocate and reallocate global arrays.
 * \details     contains functions
 *                void allocate_memory(void)
 *                void reallocate_memory_maxpart(void)
 *                void reallocate_memory_maxpartsph(void)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 03.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/*! \brief Allocates memory for global arrays.
 *
 *  This routine allocates memory for
 *  - particle storage, both the collisionless and the cells (SPH particles),
 *  - the ordered binary tree of the timeline,
 *  - communication buffers.
 *
 *  \return void
 */
void allocate_memory(void)
{
  int NTaskTimesThreads;

  NTaskTimesThreads = MaxThreads * NTask;

  Exportflag      = (int *)mymalloc("Exportflag", NTaskTimesThreads * sizeof(int));
  Exportindex     = (int *)mymalloc("Exportindex", NTaskTimesThreads * sizeof(int));
  Exportnodecount = (int *)mymalloc("Exportnodecount", NTaskTimesThreads * sizeof(int));

  Send = (struct send_recv_counts *)mymalloc("Send", sizeof(struct send_recv_counts) * NTask);
  Recv = (struct send_recv_counts *)mymalloc("Recv", sizeof(struct send_recv_counts) * NTask);

  TasksThatSend = (int *)mymalloc("TasksThatSend", sizeof(int) * NTask);
  TasksThatRecv = (int *)mymalloc("TasksThatRecv", sizeof(int) * NTask);

  Send_count  = (int *)mymalloc("Send_count", sizeof(int) * NTaskTimesThreads);
  Send_offset = (int *)mymalloc("Send_offset", sizeof(int) * NTaskTimesThreads);
  Recv_count  = (int *)mymalloc("Recv_count", sizeof(int) * NTask);
  Recv_offset = (int *)mymalloc("Recv_offset", sizeof(int) * NTask);

  Send_count_nodes  = (int *)mymalloc("Send_count_nodes", sizeof(int) * NTask);
  Send_offset_nodes = (int *)mymalloc("Send_offset_nodes", sizeof(int) * NTask);
  Recv_count_nodes  = (int *)mymalloc("Recv_count_nodes", sizeof(int) * NTask);
  Recv_offset_nodes = (int *)mymalloc("Recv_offset_nodes", sizeof(int) * NTask);

  Mesh_Send_count  = (int *)mymalloc("Mesh_Send_count", sizeof(int) * NTask);
  Mesh_Send_offset = (int *)mymalloc("Mesh_Send_offset", sizeof(int) * NTask);
  Mesh_Recv_count  = (int *)mymalloc("Mesh_Recv_count", sizeof(int) * NTask);
  Mesh_Recv_offset = (int *)mymalloc("Mesh_Recv_offset", sizeof(int) * NTask);

  Force_Send_count  = (int *)mymalloc("Force_Send_count", sizeof(int) * NTask);
  Force_Send_offset = (int *)mymalloc("Force_Send_offset", sizeof(int) * NTask);
  Force_Recv_count  = (int *)mymalloc("Force_Recv_count", sizeof(int) * NTask);
  Force_Recv_offset = (int *)mymalloc("Force_Recv_offset", sizeof(int) * NTask);

  mpi_printf("ALLOCATE: initial allocation for MaxPart = %d\n", All.MaxPart);
  P = (struct particle_data *)mymalloc_movable(&P, "P", All.MaxPart * sizeof(struct particle_data));

  mpi_printf("ALLOCATE: initial allocation for MaxPartSph = %d\n", All.MaxPartSph);
  SphP = (struct sph_particle_data *)mymalloc_movable(&SphP, "SphP", All.MaxPartSph * sizeof(struct sph_particle_data));

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
  PartSpecialListGlobal = (struct special_particle_data *)mymalloc_movable(&PartSpecialListGlobal, "PartSpecialListGlobal",
                                                                           All.MaxPartSpecial * sizeof(struct special_particle_data));
#endif /* #ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE */

  timebins_allocate(&TimeBinsHydro);
  timebins_allocate(&TimeBinsGravity);

  /* set to zero */
  memset(P, 0, All.MaxPart * sizeof(struct particle_data));
  memset(SphP, 0, All.MaxPartSph * sizeof(struct sph_particle_data));
}

/*! \brief Reallocates memory for particle data.
 *
 *  Reallocates memory for P and TimeBinsGravity arrays.
 *
 *  \return void
 */
void reallocate_memory_maxpart(void)
{
  mpi_printf("ALLOCATE: Changing to MaxPart = %d\n", All.MaxPart);

  P = (struct particle_data *)myrealloc_movable(P, All.MaxPart * sizeof(struct particle_data));
  timebins_reallocate(&TimeBinsGravity);
}

/*! \brief Reallocate memory for cell data.
 *
 *  Reallocates memory for cells in SphP and TimeBinsHydro arrays.
 *
 *  \return void
 */
void reallocate_memory_maxpartsph(void)
{
  mpi_printf("ALLOCATE: Changing to MaxPartSph = %d\n", All.MaxPartSph);

  SphP = (struct sph_particle_data *)myrealloc_movable(SphP, All.MaxPartSph * sizeof(struct sph_particle_data));
  timebins_reallocate(&TimeBinsHydro);
}
