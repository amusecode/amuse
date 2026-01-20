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
 * \file        src/add_backgroundgrid/add_bggrid.c
 * \date        05/2018
 * \brief       Re-gridding of ICs to ensure that the entire computational
 *              domain contains gas cells.
 * \details     Can be used to convert SPH ICs to Arepo ICs.
 *              contains functions:
 *                int add_backgroundgrid(void)
 *                void modify_boxsize(double new_val)
 *                void prepare_domain_backgroundgrid(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 11.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"
#include "add_bggrid.h"

#ifdef ADDBACKGROUNDGRID

static void modify_boxsize(double new_val);

MyIDType IDNew;

/*! \brief Re-gridding of ICs onto oct-tree nodes.
 *
 *  If this is active, no simulation is performed.
 *
 *  \return void
 */
int add_backgroundgrid(void)
{
  int i, no, numnodes;
  long long ngas_count_all_old;
  double vol, voltot, mgas, mtot;
  int flag_all, flag = 0;

  mpi_printf("\n\nADD BACKGROUND GRID: Adding background grid to IC file\n\n");

  for(i = 0, mgas = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      mgas += P[i].Mass;

  MPI_Allreduce(&mgas, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("ADD BACKGROUND GRID: Total gas mass before remap=%g\n", mtot);

  ngas_count_all_old = All.TotNumGas;

  ngb_treefree();

  domain_free();

  domain_Decomposition(); /* do new domain decomposition, will also make a new chained-list of synchronized particles */

  numnodes = construct_forcetree(1, 1, 0, 0); /* build tree only with gas cells */

  for(i = Tree_MaxPart, vol = 0; i < numnodes + Tree_MaxPart; i++)
    {
      if(Nodes[i].u.d.sibling == Nodes[i].u.d.nextnode) /* node is a leave */
        {
          vol += Nodes[i].len * Nodes[i].len * Nodes[i].len;
        }
    }

  for(i = 0; i < NumGas; i++)
    {
      no = Father[i];
      vol += Nodes[no].len * Nodes[no].len * Nodes[no].len / 8;
    }

  MPI_Allreduce(&vol, &voltot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("\nADD BACKGROUND GRID: voltot=%g  %g\n", voltot, pow(DomainLen, 3));

  int count_leaves = 0, count_leaves_all;

  for(i = Tree_MaxPart, vol = 0; i < numnodes + Tree_MaxPart; i++)
    {
      if(Nodes[i].u.d.sibling == Nodes[i].u.d.nextnode) /* node is a leave */
        {
          if(Nodes[i].center[0] > 0 && Nodes[i].center[0] < All.BoxSize)
            if(Nodes[i].center[1] > 0 && Nodes[i].center[1] < All.BoxSize)
              if(Nodes[i].center[2] > 0 && Nodes[i].center[2] < All.BoxSize)
                count_leaves++;
        }
    }

  MPI_Allreduce(&count_leaves, &count_leaves_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("ADD BACKGROUND GRID: count_leaves_all=%d\n\n", count_leaves_all);

  if((NumGas + count_leaves >= All.MaxPartSph) || (NumPart + count_leaves >= All.MaxPart))
    flag = 1;

  MPI_Allreduce(&flag, &flag_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  /*Increase storage for newly added gas particles */
  if(flag_all)
    domain_resize_storage(count_leaves, count_leaves, 0);

  /* determine maximum ID */
  MyIDType maxid, newid, *tmp;
  int *list;

  for(i = 0, maxid = 0; i < NumPart; i++)
    if(P[i].ID > maxid)
      maxid = P[i].ID;

  tmp = mymalloc("tmp", NTask * sizeof(MyIDType));

  MPI_Allgather(&maxid, sizeof(MyIDType), MPI_BYTE, tmp, sizeof(MyIDType), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < NTask; i++)
    if(tmp[i] > maxid)
      maxid = tmp[i];

  myfree(tmp);
  // maxid is now the total maximum ID number of all particles

  list = mymalloc("list", NTask * sizeof(int));

  MPI_Allgather(&count_leaves, 1, MPI_INT, list, 1, MPI_INT, MPI_COMM_WORLD);

  newid = maxid + 1;

  for(i = 0; i < ThisTask; i++)
    newid += list[i];

  myfree(list);

  // newid is now the maxid+total of count_leaves over all previous tasks

  IDNew = maxid + 1; /* old gas particles will have IDs below this */

  // move all particle and sph particle data down the arrays by
  // count_leaves.

  memmove(P + count_leaves, P, sizeof(struct particle_data) * NumPart);
  memmove(SphP + count_leaves, SphP, sizeof(struct sph_particle_data) * NumGas);

  NumPart += count_leaves;
  NumGas += count_leaves;

  // this is the same loop as determined count_leaves above, so
  // it will be applied count_leaves times again.
  count_leaves = 0;
  for(i = Tree_MaxPart, vol = 0; i < numnodes + Tree_MaxPart; i++)
    {
      if(Nodes[i].u.d.sibling == Nodes[i].u.d.nextnode) /* node is a leave */
        {
          if(Nodes[i].center[0] > 0 && Nodes[i].center[0] < All.BoxSize)
            if(Nodes[i].center[1] > 0 && Nodes[i].center[1] < All.BoxSize)
              if(Nodes[i].center[2] > 0 && Nodes[i].center[2] < All.BoxSize)
                {
                  P[count_leaves].Pos[0] = Nodes[i].center[0];
                  P[count_leaves].Pos[1] = Nodes[i].center[1];
                  P[count_leaves].Pos[2] = Nodes[i].center[2];
                  P[count_leaves].Vel[0] = 0;
                  P[count_leaves].Vel[1] = 0;
                  P[count_leaves].Vel[2] = 0;

                  P[count_leaves].Mass         = 0;
                  P[count_leaves].TimeBinHydro = 0;
                  P[count_leaves].TimeBinGrav  = 0;

                  P[count_leaves].Ti_Current = All.Ti_Current;

#ifdef MHD
                  SphP[count_leaves].B[0] = 0;
                  SphP[count_leaves].B[1] = 0;
                  SphP[count_leaves].B[2] = 0;
                  SphP[count_leaves].DivB = 0;
#endif /* #ifdef MHD */

                  P[count_leaves].Type          = 0;
                  P[count_leaves].SofteningType = All.SofteningTypeOfPartType[0];

                  // this puts the new ID at the right spot
                  P[count_leaves].ID = newid++;

                  SphP[count_leaves].Volume      = Nodes[i].len * Nodes[i].len * Nodes[i].len;
                  SphP[count_leaves].Utherm      = 0;
                  SphP[count_leaves].Energy      = 0;
                  SphP[count_leaves].Momentum[0] = 0;
                  SphP[count_leaves].Momentum[1] = 0;
                  SphP[count_leaves].Momentum[2] = 0;

                  count_leaves++;
                }
        }
    }

  /* Delete the force tree */
  myfree(Father);
  myfree(Nextnode);
  myfree(Tree_Points);
  force_treefree();

  calculate_weights();
  distribute_particles();

  int count_elim = 0, count_elim_all;

  for(i = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      {
        if(P[i].ID <= maxid)
          {
            // remove particle i by swapping in the last sph particle
            // and then swap the last particle to that spot
            P[i]          = P[NumGas - 1];
            P[NumGas - 1] = P[NumPart - 1];

            SphP[i] = SphP[NumGas - 1];

            NumPart--;
            NumGas--;
            i--;

            count_elim++;
          }
        else
          {
            if(P[i].Mass > 0)
              {
                SphP[i].Utherm = SphP[i].Energy / P[i].Mass;
                P[i].Vel[0]    = SphP[i].Momentum[0] / P[i].Mass;
                P[i].Vel[1]    = SphP[i].Momentum[1] / P[i].Mass;
                P[i].Vel[2]    = SphP[i].Momentum[2] / P[i].Mass;
              }
          }
      }

  MPI_Allreduce(&count_elim, &count_elim_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  sumup_large_ints(1, &NumPart, &All.TotNumPart);
  sumup_large_ints(1, &NumGas, &All.TotNumGas);

  mpi_printf("\nADD BACKGROUND GRID: count_elim_all=%d  IDNew=%d\n", count_elim_all, IDNew);
  mpi_printf("ADD BACKGROUND GRID: added particles=%d  (task 0: NumGas=%d)\n", count_leaves_all - count_elim_all, NumGas);
  mpi_printf("ADD BACKGROUND GRID: new particle number=%d\n", All.TotNumPart);
  mpi_printf("ADD BACKGROUND GRID: new gas particle number=%d\n\n", All.TotNumGas);

  for(i = 0, mgas = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      mgas += P[i].Mass;

  MPI_Allreduce(&mgas, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("ADD BACKGROUND GRID: Total gas mass after remap=%g\n", mtot);

  savepositions(0, 0);

  mpi_printf("\nADD BACKGROUND GRID: GridSize = %d\n", All.GridSize);
  mpi_printf(
      "ADD BACKGROUND GRID: Suggested value for MeanVolume = %g\nADD BACKGROUND GRID: Suggested value for ReferenceGasPartMass = %g\n",
      pow(All.BoxSize / All.GridSize, 3), mtot / ngas_count_all_old);
  mpi_printf("ADD BACKGROUND GRID: Suggested value for BoxSize = %g\n", All.BoxSize);
  mpi_printf("ADD BACKGROUND GRID: Done!\n\n");

  return 0;
}

/*! \brief Changes the box size to a new value.
 *
 *  LONG_X, LONG_Y and LONG_Z are still active as specified in Config file.
 *
 *  \param[in] new_val New box size.
 *
 *  \return void
 */
void modify_boxsize(double new_val)
{
  All.BoxSize = new_val;

  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
#endif /* #ifdef LONG_X */
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
#endif /* #ifdef LONG_Y */
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
#endif /* #ifdef LONG_Z */
}

/*! \brief Prepares computational box; makes sure simulation volume is large
 *         enough.
 *
 *  \return void
 */
void prepare_domain_backgroundgrid(void)
{
  int i, j, shift_half_box = 0, min_topleave_num = 0, set_grid_size_flag = 0;
  unsigned int size, bit_num;
  double len, xmin[3], xmax[3], xmin_glob[3], xmax_glob[3];
  double len_gas, xmin_gas[3], xmax_gas[3], xmin_gas_glob[3], xmax_gas_glob[3];
  double min_box_size, max_box_size;

  mpi_printf("\n\nADD BACKGROUND GRID: preparing domain for first domain decomposition\n");

  /* Checking GridSize limits */
  if(All.GridSize < 0)
    terminate("GridSize = %d is less than zero. This is not allowed.", All.GridSize);

  if(All.GridSize > ADDBACKGROUNDGRIDMAX)
    terminate("GridSize = %d is exceeding the max grid size = %d", All.GridSize, ADDBACKGROUNDGRIDMAX);

  if(All.GridSize > 0)
    set_grid_size_flag = 1;

  /* Now checking it is a power of two. If not assign the closest value (is this required?) */
  bit_num = 0;
  size    = ADDBACKGROUNDGRIDMAX;
  while(((size & 1) == 0) && size > 1)
    {
      size >>= 1;
      bit_num++;
    }

  for(j = 1; j < bit_num; j++)
    {
      size = All.GridSize;
      size >>= (bit_num - j);
      if((size & 1) == 1)
        break;
    }

  mpi_printf("ADD BACKGROUND GRID: original value of GridSize =  %d\n", All.GridSize);

  All.GridSize = (size << (bit_num - j - 1));

  if(All.GridSize < 1)
    All.GridSize = 1;

  mpi_printf("ADD BACKGROUND GRID: closest power of two corresponding to GridSize = %d is taken as initial guess\n", 2 * All.GridSize);

  /* determine local extension */
  for(j = 0; j < 3; j++)
    {
      xmin[j]     = MAX_REAL_NUMBER;
      xmax[j]     = -MAX_REAL_NUMBER;
      xmin_gas[j] = MAX_REAL_NUMBER;
      xmax_gas[j] = -MAX_REAL_NUMBER;
    }

  for(i = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
        {
          if(xmin[j] > P[i].Pos[j])
            xmin[j] = P[i].Pos[j];

          if(xmax[j] < P[i].Pos[j])
            xmax[j] = P[i].Pos[j];
        }
    }

  for(i = 0; i < NumGas; i++)
    {
      for(j = 0; j < 3; j++)
        {
          if(xmin_gas[j] > P[i].Pos[j])
            xmin_gas[j] = P[i].Pos[j];

          if(xmax_gas[j] < P[i].Pos[j])
            xmax_gas[j] = P[i].Pos[j];
        }
    }

  MPI_Allreduce(xmin, xmin_glob, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(xmax, xmax_glob, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(xmin_gas, xmin_gas_glob, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(xmax_gas, xmax_gas_glob, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  mpi_printf("ADD BACKGROUND GRID: Min and max coordinates.\n");
  mpi_printf("ADD BACKGROUND GRID: xmin|ymin|zmin=% g|% g|% g.\n", xmin_glob[0], xmin_glob[1], xmin_glob[2]);
  mpi_printf("ADD BACKGROUND GRID: xmax|ymax|zmax=% g|% g|% g.\n", xmax_glob[0], xmax_glob[1], xmax_glob[2]);
  mpi_printf("ADD BACKGROUND GRID: xmin_gas|ymin_gas|zmin_gas=% g|% g|% g.\n", xmin_gas_glob[0], xmin_gas_glob[1], xmin_gas_glob[2]);
  mpi_printf("ADD BACKGROUND GRID: xmax_gas|ymax_gas|zmax_gas=% g|% g|% g.\n", xmax_gas_glob[0], xmax_gas_glob[1], xmax_gas_glob[2]);

  len     = 0;
  len_gas = 0;
  for(j = 0; j < 3; j++)
    {
      if(xmax_glob[j] - xmin_glob[j] > len)
        len = xmax_glob[j] - xmin_glob[j];

      if(xmax_gas_glob[j] - xmin_gas_glob[j] > len_gas)
        len_gas = xmax_gas_glob[j] - xmin_gas_glob[j];

      if(xmin_glob[j] < 0)
        shift_half_box = 1;
    }

  max_box_size = FACTOR_MAX_BOX_SIZE * len_gas;
  min_box_size = FACTOR_MIN_BOX_SIZE * len_gas;

  if(All.BoxSize < min_box_size)
    {
      mpi_printf("ADD BACKGROUND GRID: Need to increase the BoxSize. Old value = %g, new value = %g\n", All.BoxSize, min_box_size);
      modify_boxsize(min_box_size);
    }
  if(All.BoxSize > max_box_size)
    {
      mpi_printf("ADD BACKGROUND GRID: Need to decrease the BoxSize. Old value = %g, new value = %g\n", All.BoxSize, max_box_size);
      modify_boxsize(max_box_size);
    }

  mpi_printf("ADD BACKGROUND GRID: Domain extent %g, BoxSize = %g, ratio = %g\n", len, All.BoxSize, len / All.BoxSize);
  mpi_printf("ADD BACKGROUND GRID: Gas extent %g, BoxSize = %g, ratio = %g\n", len_gas, All.BoxSize, len_gas / All.BoxSize);

  /* the terminate condition must be checked properly */
  if(!set_grid_size_flag)
    {
      while(min_topleave_num < NTask && (All.BoxSize / len_gas) > All.GridSize && All.GridSize < ADDBACKGROUNDGRIDMAX)
        {
          All.GridSize <<= 1;
          min_topleave_num = (int)pow(len_gas * All.GridSize / All.BoxSize, 3.0);
          mpi_printf("ADD BACKGROUND GRID: GridSize=%3d, min_topleave_num=%6d, NTask=%6d, BoxSize/GridSize=%g, len_gas/GridSize=%g\n",
                     All.GridSize, min_topleave_num, NTask, All.BoxSize / All.GridSize, len_gas / All.BoxSize);
        }
    }
  else
    {
      All.GridSize <<= 1;
      min_topleave_num = (int)pow(len_gas * All.GridSize / All.BoxSize, 3.0);
      mpi_printf("ADD BACKGROUND GRID: GridSize=%3d, min_topleave_num=%6d, NTask=%6d, BoxSize/GridSize=%g, len_gas/GridSize=%g\n",
                 All.GridSize, min_topleave_num, NTask, All.BoxSize / All.GridSize, len_gas / All.BoxSize);
    }

  if(min_topleave_num < NTask)
    {
      char buf[500];
      sprintf(buf,
              "min_topleave_num=%d < NTask=%d, MaxGridSize=%d. Try either to run with less task or to set the BoxSize to a smaller "
              "value\n",
              min_topleave_num, NTask, ADDBACKGROUNDGRIDMAX);
      terminate(buf);
    }

  if(len_gas / All.BoxSize > All.GridSize)
    {
      char buf[500];
      sprintf(buf, "len_gas/BoxSize=%g > GridSize=%d, MaxGridSize=%d. GridSize should be increased if possible\n",
              len_gas / All.BoxSize, All.GridSize, ADDBACKGROUNDGRIDMAX);
      terminate(buf);
    }

  if(shift_half_box)
    {
      mpi_printf("ADD BACKGROUND GRID: Need to shift particles by half box size\n\n");
      for(i = 0; i < NumPart; i++)
        {
          P[i].Pos[0] += 0.5 * All.BoxSize;
          P[i].Pos[1] += 0.5 * All.BoxSize;
          P[i].Pos[2] += 0.5 * All.BoxSize;
        }
    }
}

#endif /* #ifdef ADDBACKGROUNDGRID */
