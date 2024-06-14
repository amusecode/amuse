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
 * \file        src/mesh/voronoi/voronoi_refinement.c
 * \date        05/2018
 * \brief       Contains routines for refinement.
 * \details     contains functions:
 *                static void refine_add_ngb(int i, int j)
 *                int do_refinements(void)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 23.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include "../../main/allvars.h"
#include "../../main/proto.h"

#ifdef REFINEMENT_SPLIT_CELLS

static int *ref_SphP_dp_index;
static int *first_ngb, *last_ngb, first_free_ngb;

/*! \brief Linked list for neighbor data.
 *
 */
static struct ngb_data
{
  int index;
  int next_ngb;
} * ngb;

/*! \brief Add element to linked neighbor list.
 *
 *  \param[in] i Index of existing cell.
 *  \param[in] j Index of new cell.
 *
 *  \return void
 */
static void refine_add_ngb(int i, int j)
{
  if(i >= 0 && j >= 0)
    {
      if(i >= Mesh.Ndp || j >= Mesh.Ndp)
        {
          terminate("i>= Ndp || j>= Ndp");
        }

      if(first_ngb[i] >= 0)
        {
          ngb[last_ngb[i]].next_ngb = first_free_ngb;
          last_ngb[i]               = first_free_ngb;
        }
      else
        {
          first_ngb[i] = last_ngb[i] = first_free_ngb;
        }

      ngb[first_free_ngb].index    = j;
      ngb[first_free_ngb].next_ngb = -1;
      first_free_ngb++;
    }
}

/*! \brief Loops through active cells and refine cells if needed.
 *
 *  Splits the cell in random direction; moves mesh-generating point by
 *  0.025 cell radius and inserts a second mesh-generating point opposite to
 *  split the cell into two.
 *
 *  \return Number of cells that were refined.
 */
int do_refinements(void)
{
  char buf[1000];
  int idx, i, j, k, count, countall;
  double rad, fac;
  MyIDType newid = 0;

  TIMER_START(CPU_REFINE);

  ref_SphP_dp_index = mymalloc_movable(&ref_SphP_dp_index, "ref_SphP_dp_index", NumGas * sizeof(int));

  int NActiveParticles = TimeBinsHydro.NActiveParticles; /* save this since refinement is going to change it */
  for(idx = 0, count = 0; idx < NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(should_this_cell_be_split(i))
        {
          ref_SphP_dp_index[i] = -1;
          count++;
        }
    }

  MPI_Allreduce(&count, &countall, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  mpi_printf("REFINE: want to refine %d cells\n", countall);

  if(countall)
    {
      domain_resize_storage(count, count, 2);

      if(NumPart + count >= All.MaxPart)
        {
          sprintf(buf, "On Task=%d with NumPart=%d we try to produce %d cells. Sorry, no space left...(All.MaxPart=%d)\n", ThisTask,
                  NumPart, count, All.MaxPart);
          terminate(buf);
        }

      if(NumGas + count >= All.MaxPartSph)
        {
          sprintf(buf, "On Task=%d with NumGas=%d we try to produce %d cells. Sorry, no space left...(All.MaxPartSph=%d)\n", ThisTask,
                  NumGas, count, All.MaxPartSph);
          terminate(buf);
        }

      if(All.MaxID == 0) /* MaxID not calculated yet */
        calculate_maxid();

      int *list = mymalloc("list", NTask * sizeof(int));

      MPI_Allgather(&count, 1, MPI_INT, list, 1, MPI_INT, MPI_COMM_WORLD);

      newid = All.MaxID + 1;

      for(i = 0; i < ThisTask; i++)
        newid += list[i];

      All.MaxID += countall;

      myfree(list);

      Ngb_MarkerValue++;
      int nchanged  = 0;
      int *nodelist = (int *)mymalloc("nodelist", NTopleaves * sizeof(int));

      /*  create explicit list of neighbors */

      first_ngb = mymalloc("first_ngb", Mesh.Ndp * sizeof(int));
      ngb       = mymalloc("ngbs", 2 * Mesh.Nvf * sizeof(struct ngb_data));
      last_ngb  = mymalloc("last_ngb", Mesh.Ndp * sizeof(int));

      for(i = 0; i < Mesh.Ndp; i++)
        {
          first_ngb[i] = last_ngb[i] = -1;

          if(Mesh.DP[i].task == ThisTask)
            {
              int li = Mesh.DP[i].index;
              if(li >= 0 && li < NumGas)
                if(ref_SphP_dp_index[li] < 0)
                  ref_SphP_dp_index[li] = i; /* only guaranteed to be set for active cells */
            }
        }

      for(i = 0, first_free_ngb = 0; i < Mesh.Nvf; i++)
        {
          refine_add_ngb(Mesh.VF[i].p1, Mesh.VF[i].p2);
          refine_add_ngb(Mesh.VF[i].p2, Mesh.VF[i].p1);
        }

      myfree(last_ngb);

      int NActiveParticles = TimeBinsHydro.NActiveParticles;
      for(idx = 0, count = 0; idx < NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(should_this_cell_be_split(i))
            {
              int addToGravList = TimeBinSynchronized[P[i].TimeBinGrav];
              if(NumPart > NumGas)
                {
                  move_collisionless_particle(NumPart + count, NumGas + count);
                  if(TimeBinSynchronized[P[NumPart + count].TimeBinGrav] && P[i].Mass > 0)
                    addToGravList = 0;

                  /* there is already an entry in the list of active particles for
                     gravity that points to the index that we will use for our new cell */
                }

              /* now split the gas cell */

              j = NumGas + count;

              P[j]    = P[i];
              SphP[j] = SphP[i];

              P[j].ID = newid++;

              rad = get_cell_radius(i);

              double dir[3];
#ifdef TWODIMS
              double phi = 2 * M_PI * get_random_number();

              dir[0] = cos(phi);
              dir[1] = sin(phi);
              dir[2] = 0;
#else  /* #ifdef TWODIMS */
              double theta = acos(2 * get_random_number() - 1);
              double phi   = 2 * M_PI * get_random_number();

              dir[0] = sin(theta) * cos(phi);
              dir[1] = sin(theta) * sin(phi);
              dir[2] = cos(theta);
#endif /* #ifdef TWODIMS */
              fac = 0.025 * rad;

              P[j].Pos[0] = P[i].Pos[0] + fac * dir[0];
              P[j].Pos[1] = P[i].Pos[1] + fac * dir[1];
              P[j].Pos[2] = P[i].Pos[2] + fac * dir[2];

              SphP[j].SepVector[0] = SphP[i].SepVector[0] = dir[0];
              SphP[j].SepVector[1] = SphP[i].SepVector[1] = dir[1];
              SphP[j].SepVector[2] = SphP[i].SepVector[2] = dir[2];

              /**** create the voronoi cell of i as an auxiliary mesh */

              int jj = ref_SphP_dp_index[i]; /* this is the delaunay point of this cell */
              if(jj < 0)
                terminate("jj < 0");

              initialize_and_create_first_tetra(&DeRefMesh);

              DeRefMesh.DTC = mymalloc_movable(&DeRefMesh.DTC, "DeRefDTC", DeRefMesh.MaxNdt * sizeof(tetra_center));
              DeRefMesh.DTF = mymalloc_movable(&DeRefMesh.DTF, "DeRefDTF", DeRefMesh.MaxNdt * sizeof(char));
              for(k = 0; k < DeRefMesh.Ndt; k++)
                DeRefMesh.DTF[k] = 0;

              int tlast = 0;

              k = first_ngb[jj];
              while(k >= 0)
                {
                  int q = ngb[k].index;

                  if(DeRefMesh.Ndp + 2 >= DeRefMesh.MaxNdp)
                    {
                      DeRefMesh.Indi.AllocFacNdp *= ALLOC_INCREASE_FACTOR;
                      DeRefMesh.MaxNdp = DeRefMesh.Indi.AllocFacNdp;
#ifdef VERBOSE
                      printf("Task=%d: increase memory allocation, MaxNdp=%d Indi.AllocFacNdp=%g\n", ThisTask, DeRefMesh.MaxNdp,
                             DeRefMesh.Indi.AllocFacNdp);
#endif /* #ifdef VERBOSE */
                      DeRefMesh.DP -= 5;
                      DeRefMesh.DP = myrealloc_movable(DeRefMesh.DP, (DeRefMesh.MaxNdp + 5) * sizeof(point));
                      DeRefMesh.DP += 5;
                    }

                  DeRefMesh.DP[DeRefMesh.Ndp] = Mesh.DP[q];

                  double r =
                      sqrt(pow(DeRefMesh.DP[DeRefMesh.Ndp].x - P[i].Pos[0], 2) + pow(DeRefMesh.DP[DeRefMesh.Ndp].y - P[i].Pos[1], 2) +
                           pow(DeRefMesh.DP[DeRefMesh.Ndp].z - P[i].Pos[2], 2));

                  if(r < 2 * fac)
                    terminate("We are trying to split a heavily distorted cell... We better stop. Check your refinement criterion.");

#ifndef OPTIMIZE_MEMORY_USAGE
                  set_integers_for_point(&DeRefMesh, DeRefMesh.Ndp);
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE */
                  tlast = insert_point(&DeRefMesh, DeRefMesh.Ndp, tlast);

                  DeRefMesh.Ndp++;
                  k = ngb[k].next_ngb;
                }

              /* now add also the point jj itself (the one that is to be split */

              DeRefMesh.DP[DeRefMesh.Ndp] = Mesh.DP[jj];
#ifndef OPTIMIZE_MEMORY_USAGE
              set_integers_for_point(&DeRefMesh, DeRefMesh.Ndp);
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE */
              tlast = insert_point(&DeRefMesh, DeRefMesh.Ndp, tlast);
              DeRefMesh.Ndp++;

              /* and finally, add the newly generated point */

              DeRefMesh.DP[DeRefMesh.Ndp].x  = P[j].Pos[0];
              DeRefMesh.DP[DeRefMesh.Ndp].y  = P[j].Pos[1];
              DeRefMesh.DP[DeRefMesh.Ndp].z  = P[j].Pos[2];
              DeRefMesh.DP[DeRefMesh.Ndp].ID = P[j].ID;
#ifndef OPTIMIZE_MEMORY_USAGE
              set_integers_for_point(&DeRefMesh, DeRefMesh.Ndp);
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE */
              tlast = insert_point(&DeRefMesh, DeRefMesh.Ndp, tlast);
              DeRefMesh.Ndp++;

              /* compute circumcircles */
              compute_circumcircles(&DeRefMesh);

              double *Volume = mymalloc("Volume", DeRefMesh.Ndp * sizeof(double));

              derefine_refine_compute_volumes(Volume);

              double voli = Volume[DeRefMesh.Ndp - 2];
              double volj = Volume[DeRefMesh.Ndp - 1];

              myfree(Volume);

              myfree(DeRefMesh.DTF);
              myfree(DeRefMesh.DTC);
              DeRefMesh.DTC = NULL;

              myfree(DeRefMesh.DT);
              myfree(DeRefMesh.DP - 5);
              myfree(DeRefMesh.VF);

              /* now split the conserved variables according to the volume ratio of the split */

              double faci = voli / (voli + volj);
              double facj = volj / (voli + volj);

              P[i].Mass *= faci;
              P[j].Mass *= facj;
              SphP[i].OldMass *= faci;
              SphP[j].OldMass *= facj;

              SphP[i].Energy *= faci;
              SphP[j].Energy *= facj;

#ifdef MHD
              for(k = 0; k < 3; k++)
                {
                  SphP[i].B[k] = SphP[i].BConserved[k] / (voli + volj);
                  SphP[j].B[k] =
                      SphP[i].B[k] + SphP[i].Grad.dB[k][0] * (P[j].Pos[0] - P[i].Pos[0]) +
                      SphP[i].Grad.dB[k][1] * (P[j].Pos[1] - P[i].Pos[1]) +
                      SphP[i].Grad.dB[k][2] * (P[j].Pos[2] - P[i].Pos[2]); /* extrapolate B to the position of the new cell */

                  /* update conserved variables */
                  SphP[i].BConserved[k] = SphP[i].B[k] * voli;
                  SphP[j].BConserved[k] = SphP[j].B[k] * volj;
                }
#endif /* #ifdef MHD */

              for(k = 0; k < 3; k++)
                {
                  SphP[i].Momentum[k] *= faci;
                  SphP[j].Momentum[k] *= facj;
                }

#ifdef USE_SFR
              SphP[i].Sfr *= faci;
              SphP[j].Sfr *= facj;
#endif /* #ifdef USE_SFR */

#ifdef MAXSCALARS
              for(int s = 0; s < N_Scalar;
                  s++) /* Note, the changes in MATERIALS, HIGHRESGASMASS, etc., are treated as part of the Scalars */
                {
                  *(MyFloat *)(((char *)(&SphP[i])) + scalar_elements[s].offset_mass) *= faci;
                  *(MyFloat *)(((char *)(&SphP[j])) + scalar_elements[s].offset_mass) *= facj;
                }
#endif /* #ifdef MAXSCALARS */

#ifdef REFINEMENT_HIGH_RES_GAS
              /* the change in the SphP[].HighResMass is treated as part of the Scalars loop above */
              SphP[i].AllowRefinement += 2; /* increment the refinement "generation" of both cells */
              SphP[j].AllowRefinement += 2;
#endif /* #ifdef REFINEMENT_HIGH_RES_GAS */

              /* add the new particle into the neighbour tree */
              int no          = Ngb_Nextnode[i];
              Ngb_Nextnode[i] = j;
              Ngb_Nextnode[j] = no;
              Ngb_Father[j]   = Ngb_Father[i];

              ngb_update_rangebounds(j, &nchanged, nodelist);

              /* now add the new particle into the link-lists for the time integration */

              timebin_add_particle(&TimeBinsHydro, j, i, P[i].TimeBinHydro, 1);
              timebin_add_particle(&TimeBinsGravity, j, i, P[i].TimeBinGrav, addToGravList);

              SphP[j].first_connection = -1;
              SphP[j].last_connection  = -1;

              count++;
            }
        }

      NumPart += count;
      NumGas += count;
      All.TotNumPart += countall;
      All.TotNumGas += countall;

      myfree(ngb);
      myfree(first_ngb);

      ngb_finish_rangebounds_update(nchanged, nodelist);

      myfree(nodelist);
    }

  myfree(ref_SphP_dp_index);

  TIMER_STOP(CPU_REFINE);

  return countall;
}

#endif /* REFINEMENT_SPLIT_CELLS */
