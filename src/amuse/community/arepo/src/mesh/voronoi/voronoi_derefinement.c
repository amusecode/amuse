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
 * \file        src/mesh/voronoi/voronoi_derefinement.c
 * \date        05/2018
 * \brief       Contains routines for de-refinement.
 * \details     contains functions:
 *                static void derefine_add_ngb(int edge, int i, int j, double
 *                  area, int t, int nr)
 *                int do_derefinements(void)
 *                static void derefine_apply_probe_list(void)
 *                static void derefine_apply_flux_list(void)
 *                static int derefine_flux_list_data_compare(const void *a,
 *                  const void *b)
 *                static int derefine_probe_list_data_compare_task(const
 *                  void *a, const void *b)
 *                static int derefine_compare_seq_DP_ID(const void *a,
 *                  const void *b)
 *                static void derefine_exchange_flag(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 22.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../main/allvars.h"
#include "../../main/proto.h"

#include "voronoi.h"

#if defined(REFINEMENT_MERGE_CELLS) && !defined(ONEDIMS)
#define DEREF_SA_FAC 1.0e-4

int do_derefinements(void);
static void derefine_add_ngb(int edge, int i, int j, double area, int tt, int nr);
static int derefine_compare_seq_DP_ID(const void *a, const void *b);
static int derefine_flux_list_data_compare(const void *a, const void *b);
static void derefine_apply_flux_list(void);
static void derefine_exchange_flag(void);
static void derefine_apply_probe_list(void);
static int derefine_probe_list_data_compare_task(const void *a, const void *b);

/*! \brief Data for derefinement: flag for de-refinement and index of cell.
 */
static struct derefine_particle_data
{
  int Flag;
  int dp_index;
} * deref_SphP;

/*! \brief Data structure for communicating de-refinement flags.
 */
static struct flagexch
{
  int Flag;
  MyIDType ID;
} * FlagExch;

/*! \brief Data structure to flag Delaunay data.
 */
static struct flag_delaunay_data
{
  int Flag;
} * flag_DP;

/*! \brief Structure defining auxiliary Delaunay data (for sorting).
 */
static struct seq_delaunay_data
{
  MyFloat rnd;
  int rank, index;
  MyIDType ID;
} * seq_DP;

/*! \brief Structure defining probe list element.
 */
static struct probe_list_data
{
  int task, index;
  int sendpart;
  int flag;
} * ProbeList;

/*! \brief Structure defining flux list element.
 */
static struct flux_list_data
{
  int task, index;
  double dM, dP[3];
#ifdef MHD
  double dB[3];
#endif /* #ifdef MHD */

#ifndef ISOTHERM_EQS
  double dEnergy;
#endif /* #ifndef ISOTHERM_EQS */

#ifdef MAXSCALARS
  double dConservedScalars[MAXSCALARS];
#endif /* #ifdef MAXSCALARS */
} * FluxList;

static int Nflux, MaxNflux;

static int *first_ngb, *last_ngb, first_free_ngb;

/*! \brief Structure defining neighbour data.
 */
static struct ngb_data
{
#ifdef OPTIMIZE_MEMORY_USAGE
  MyFloat area;
#else  /* #ifdef OPTIMIZE_MEMORY_USAGE */
  double area;
#endif /* #ifdef OPTIMIZE_MEMORY_USAGE #else */
  int index;
  int edge;
  int next_ngb;
  int t, nr; /* delaunay tetra and edge number that generated this face */
} * ngb;

static int n_tri, max_n_tri;
static triangle *trilist;

#ifdef REFINEMENT_SPLIT_CELLS
extern char *FlagDoNotRefine;
#endif /* #ifdef REFINEMENT_SPLIT_CELLS */

/*! \brief Adds cell in list ngb.
 *
 *  \param[in] edge Element 'edge' in ngb.
 *  \param[in] i Index in first_ngb and last_ngb lists.
 *  \param[in] j Element 'index' in ngb.
 *  \param[in] area Element 'area' in ngb.
 *  \param[in] t Element 't' in ngb.
 *  \param[in] nr Element 'nr' in ngb.
 *
 *  \return void
 */
static void derefine_add_ngb(int edge, int i, int j, double area, int t, int nr)
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

      ngb[first_free_ngb].area     = area;
      ngb[first_free_ngb].edge     = edge;
      ngb[first_free_ngb].t        = t;
      ngb[first_free_ngb].nr       = nr;
      ngb[first_free_ngb].index    = j;
      ngb[first_free_ngb].next_ngb = -1;
      first_free_ngb++;
    }
}

/*! \brief Loop over all active cells and derefine the ones that need to be
 *         derefined.
 *
 *  \return Number of derefined cells.
 */
int do_derefinements(void)
{
  int idx, i, j, k, count, countall;

  TIMER_START(CPU_DEREFINE);

  deref_SphP = mymalloc_movable(&deref_SphP, "deref_SphP", NumGas * sizeof(struct derefine_particle_data));

  FlagExch = mymalloc_movable(&FlagExch, "FlagExch", Mesh_nimport * sizeof(struct flagexch));

  /* first, check whether we have cells to derefine */
  for(idx = 0, count = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
#ifdef REFINEMENT_SPLIT_CELLS
      FlagDoNotRefine[i] = 0;
#endif /* #ifdef REFINEMENT_SPLIT_CELLS */

      if(i >= NumGas)
        terminate("index of gas cell greater than NumGas");

      deref_SphP[i].Flag     = 0;
      deref_SphP[i].dp_index = -1;

      if(derefine_should_this_cell_be_merged(i, deref_SphP[i].Flag))
        {
          deref_SphP[i].Flag = 1;
          count++;
        }
    }

  MPI_Allreduce(&count, &countall, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("DEREFINE: Number of cells that want to be de-refined: %d\n", countall);

  if(countall)
    {
      int max_assumed_ntri = 0;

      /* tell the ghost cells whether they want to be refined or not */
      derefine_exchange_flag();

      /* let's create an explicit list of the neighbors of each cell */

      first_ngb = mymalloc("first_ngb", Mesh.Ndp * sizeof(int));
      ngb       = mymalloc("ngb", 2 * Mesh.Nvf * sizeof(struct ngb_data));

      last_ngb = mymalloc("last_ngb", Mesh.Ndp * sizeof(int));

      for(i = 0; i < Mesh.Ndp; i++)
        first_ngb[i] = last_ngb[i] = -1;

      for(i = 0, first_free_ngb = 0; i < Mesh.Nvf; i++)
        {
          derefine_add_ngb(i, Mesh.VF[i].p1, Mesh.VF[i].p2, Mesh.VF[i].area, Mesh.VF[i].t, Mesh.VF[i].nr);
          derefine_add_ngb(i, Mesh.VF[i].p2, Mesh.VF[i].p1, Mesh.VF[i].area, Mesh.VF[i].t, Mesh.VF[i].nr);
        }

      myfree(last_ngb);

      /* we now make a list of the delaunay points that we can sort in a globally unique way */
      flag_DP = mymalloc_movable(&flag_DP, "flag_DP", Mesh.Ndp * sizeof(struct flag_delaunay_data));
      seq_DP  = mymalloc("seq_DP", Mesh.Ndp * sizeof(struct seq_delaunay_data));

      for(i = 0; i < Mesh.Ndp; i++)
        {
          seq_DP[i].rank  = i;
          seq_DP[i].index = Mesh.DP[i].index;

          if(Mesh.DP[i].task == ThisTask)
            {
              int li = Mesh.DP[i].index;
              if(li < 0)
                {
                  flag_DP[i].Flag = 0;
                  seq_DP[i].ID    = 0;
                  seq_DP[i].rnd   = 0;
                }
              else
                {
                  if(li < NumGas)
                    if(deref_SphP[li].dp_index < 0)
                      deref_SphP[li].dp_index = i; /* only guaranteed to be set for active cells */

                  if(li >= NumGas)
                    li -= NumGas;

                  flag_DP[i].Flag = deref_SphP[li].Flag;
                  seq_DP[i].ID    = P[li].ID;
                  seq_DP[i].rnd   = get_random_number();
                }
            }
          else
            {
              flag_DP[i].Flag = FlagExch[Mesh.DP[i].index].Flag;
              seq_DP[i].ID    = FlagExch[Mesh.DP[i].index].ID;
              seq_DP[i].rnd   = get_random_number();
            }
        }

      /* sort according to ID */
      mysort(seq_DP, Mesh.Ndp, sizeof(struct seq_delaunay_data), derefine_compare_seq_DP_ID);

      /* now let's go through in sorted order. For each cell that is supposed to be refined, check whether any of the
       * neighbors is already refined. If yes, don't allow it to be refined.
       * Also, if there is a neighbour with the same ID, don't refine it, because this must be a mirrored particle
       */

      for(i = 0; i < Mesh.Ndp; i++)
        {
          if(seq_DP[i].ID != 0)
            {
              j = seq_DP[i].rank;

              if(flag_DP[j].Flag == 1) /* this cell is still eligible for derefinement */
                {
                  /* go through its neighbours and check whether one of them is already up for derefinement */

                  int n = 0;
                  k     = first_ngb[j];
                  while(k >= 0)
                    {
                      /* we only need to consider neighboring cells if they are active */
                      int q = ngb[k].index;

                      if(q >= 0)
                        {
                          int timebin;

                          if(Mesh.DP[q].task == ThisTask)
                            {
                              if(Mesh.DP[q].index < NumGas)
                                timebin = P[Mesh.DP[q].index].TimeBinHydro;
                              else
                                timebin = P[Mesh.DP[q].index - NumGas].TimeBinHydro;
                            }
                          else
                            {
#ifndef OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT
                              timebin = PrimExch[Mesh.DP[q].index].TimeBinHydro;
#else  /* #ifndef OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT */
                              timebin = RefExch[Mesh.DP[q].index].TimeBinHydro;
#endif /* #ifndef OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT #else */
                            }

                          if(TimeBinSynchronized[timebin])
                            {
                              if(flag_DP[q].Flag == 2 || flag_DP[q].Flag == 3)
                                n++;

                              if(Mesh.DP[q].ID == seq_DP[i].ID) /* same ID, so we have a mirrored particle */
                                n++;
                            }
                        }

                      k = ngb[k].next_ngb;
                    }

                  if(n == 0) /* ok, none have been found. This means this cell is allowed to be refined */
                    flag_DP[j].Flag = 2;
                  else
                    flag_DP[j].Flag = 3;
                }
            }
        }

      myfree(seq_DP);

      /* copy of the refinement flags to the cell structure */
      for(i = 0; i < Mesh.Ndp; i++)
        if(Mesh.DP[i].task == ThisTask && Mesh.DP[i].index >= 0 && Mesh.DP[i].index < NumGas)
          deref_SphP[Mesh.DP[i].index].Flag = flag_DP[i].Flag;

      myfree(flag_DP);

      /* now let's count again how many cells we would like to derefine */

      for(idx = 0, count = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(deref_SphP[i].Flag == 2)
            count++;
        }

      int in[2], out[2];
      in[0] = count;

      /* now we carry out an auxiliary check to make sure that we really
         avoid de-refining two neighboring cells.  If such a pair is
         found, both cells will not be derefined. */

      MaxNflux  = Mesh.Indi.AllocFacNflux;
      Nflux     = 0;
      ProbeList = mymalloc_movable(&ProbeList, "ProbeList", MaxNflux * sizeof(struct probe_list_data));

      count = 0;

      for(idx = 0, count = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(deref_SphP[i].Flag == 2)
            {
              j = deref_SphP[i].dp_index; /* this is the delaunay point of this cell */
              if(j < 0)
                terminate("j < 0");

              k = first_ngb[j];

              int flag = 0;

              while(k >= 0)
                {
                  if(ngb[k].area > DEREF_SA_FAC * SphP[i].SurfaceArea)
                    {
                      int q = ngb[k].index;

                      if(Mesh.DP[q].task == ThisTask)
                        {
                          int p = Mesh.DP[q].index;

                          if(p < 0)
                            terminate("p < 0");

                          if(p >= NumGas) /* this is a local ghost point */
                            p -= NumGas;

                          if(TimeBinSynchronized[P[p].TimeBinHydro])
                            if(deref_SphP[p].Flag == 2)
                              flag++;
                        }
                      else
                        {
                          /* here we have a foreign ghost point */
                          if(Nflux >= MaxNflux)
                            {
                              Mesh.Indi.AllocFacNflux *= ALLOC_INCREASE_FACTOR;
                              MaxNflux = Mesh.Indi.AllocFacNflux;
#ifdef VERBOSE
                              printf("Task=%d: increase memory allocation, MaxNflux=%d Indi.AllocFacNflux=%g\n", ThisTask, MaxNflux,
                                     Mesh.Indi.AllocFacNflux);
#endif /* #ifdef VERBOSE */
                              ProbeList = myrealloc_movable(ProbeList, MaxNflux * sizeof(struct probe_list_data));

                              if(Nflux >= MaxNflux)
                                terminate("Nflux >= MaxNflux");
                            }

                          ProbeList[Nflux].task     = Mesh.DP[q].task;
                          ProbeList[Nflux].index    = Mesh.DP[q].originalindex;
                          ProbeList[Nflux].sendpart = i;
                          ProbeList[Nflux].flag     = 0;

                          Nflux++;
                        }
                    }
                  k = ngb[k].next_ngb;
                }

              if(flag)
                {
                  /* ups. It looks like a neigboring point is also about to be dissolved. We hence do not
                     dissolve the current point
                   */
                  deref_SphP[i].Flag = 0;
                  count++;
                }
            }
        }

      /* now let's probe on other tasks */

      derefine_apply_probe_list();

      for(i = 0; i < Nflux; i++)
        {
          if(ProbeList[i].flag)
            if(deref_SphP[ProbeList[i].sendpart].Flag == 2)
              {
                deref_SphP[ProbeList[i].sendpart].Flag = 0;
                count++;
              }
        }

      myfree(ProbeList);

      in[1] = count;
      MPI_Reduce(in, out, 2, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      mpi_printf("DEREFINE: Number of cells that we could de-refine: %d, number of cells we exclude from this set:  %d\n", out[0],
                 out[1]);

      /* we now distribute the conserved quantities of the cell among the neighbours */

      MaxNflux = Mesh.Indi.AllocFacNflux;
      Nflux    = 0;
      FluxList = mymalloc_movable(&FluxList, "FluxList", MaxNflux * sizeof(struct flux_list_data));

      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(deref_SphP[i].Flag == 2)
            {
              j = deref_SphP[i].dp_index; /* this is the delaunay point of this cell */
              if(j < 0)
                terminate("j < 0");

              max_n_tri = 300000;
              n_tri     = 0;

              trilist = mymalloc("trilist", max_n_tri * sizeof(triangle));

              /* get a list of all the triangles that make up the Voronoi cell of j */
              k = first_ngb[j];
              while(k >= 0)
                {
                  n_tri = derefine_refine_get_triangles(&Mesh, ngb[k].t, ngb[k].nr, &Mesh.DP[j], trilist, n_tri, max_n_tri);

                  k = ngb[k].next_ngb;
                }

              /* assign the first point as owner to all tetras */
              k     = first_ngb[j];
              int q = ngb[k].index;
              int t;
              for(t = 0; t < n_tri; t++)
                trilist[t].owner = q;

              double vol = 0;
              for(k = 0; k < n_tri; k++)
                vol += get_tri_volume(k, trilist);

              /* now consider all the other points and split the triangles if needed */
              k = first_ngb[j];
              k = ngb[k].next_ngb;
              while(k >= 0)
                {
                  int q = ngb[k].index;
                  n_tri = derefine_add_point_and_split_tri(q, trilist, n_tri, max_n_tri, vol);
                  k     = ngb[k].next_ngb;
                }

              if(n_tri > max_assumed_ntri)
                max_assumed_ntri = n_tri;

              double *volume = mymalloc("volume", Mesh.Ndp * sizeof(double));

              /* clear the volume entries of the neighbors */
              k = first_ngb[j];
              while(k >= 0)
                {
                  int q     = ngb[k].index;
                  volume[q] = 0;
                  k         = ngb[k].next_ngb;
                }

              /* now assign the volume of the triangles to the neighbors */
              for(k = 0; k < n_tri; k++)
                {
                  if(trilist[k].owner < 0 || trilist[k].owner >= Mesh.Ndp)
                    terminate("strange owner");

                  volume[trilist[k].owner] += get_tri_volume(k, trilist);
                }

              /* first, let's establish the surface area sum for this cell */
              double voltot = 0;
              k             = first_ngb[j];
              while(k >= 0)
                {
                  if(ngb[k].area > DEREF_SA_FAC * SphP[i].SurfaceArea)
                    {
                      int q = ngb[k].index;
                      voltot += volume[q];
                    }
                  k = ngb[k].next_ngb;
                }

              /* now, distribute conserved quantities proportional to the gained volume */
              double facsum = 0;
              k             = first_ngb[j];
              while(k >= 0)
                {
                  if(ngb[k].area > DEREF_SA_FAC * SphP[i].SurfaceArea)
                    {
                      int q = ngb[k].index;

                      double fac = volume[q] / voltot;

                      if(fac < 0)
                        {
                          warn("strange: fac=%g\n", fac);
                          fac = 0;
                          // terminate("strange");
                        }
                      facsum += fac;

                      if(Mesh.DP[q].task == ThisTask)
                        {
                          int p = Mesh.DP[q].index;

                          if(p < 0)
                            terminate("p < 0");

                          if(p >= NumGas) /* this is a local ghost point */
                            p -= NumGas;
                          P[p].Mass += fac * P[i].Mass;
                          SphP[p].Momentum[0] += fac * SphP[i].Momentum[0];
                          SphP[p].Momentum[1] += fac * SphP[i].Momentum[1];
                          SphP[p].Momentum[2] += fac * SphP[i].Momentum[2];

#ifdef MHD
                          SphP[p].BConserved[0] += fac * SphP[i].BConserved[0];
                          SphP[p].BConserved[1] += fac * SphP[i].BConserved[1];
                          SphP[p].BConserved[2] += fac * SphP[i].BConserved[2];
#endif /* #ifdef MHD */

#ifndef ISOTHERM_EQS
                          SphP[p].Energy += fac * SphP[i].Energy;
#endif /* #ifndef ISOTHERM_EQS */

#ifdef MAXSCALARS
                          for(int s = 0; s < N_Scalar; s++)
                            *(MyFloat *)(((char *)(&SphP[p])) + scalar_elements[s].offset_mass) +=
                                fac * (*(MyFloat *)(((char *)(&SphP[i])) + scalar_elements[s].offset_mass));
#endif /* #ifdef MAXSCALARS */

#ifdef REFINEMENT_SPLIT_CELLS
                          FlagDoNotRefine[p] = 1;
#endif /* #ifdef REFINEMENT_SPLIT_CELLS */
                        }
                      else
                        {
                          /* here we have a foreign ghost point */
                          if(Mesh.DP[q].originalindex < 0)
                            {
                              char buf[1000];
                              sprintf(buf, "---> task=%d  q=%d j=%d Ndp=%d\n", ThisTask, q, j, Mesh.Ndp);
                              terminate(buf);
                            }

                          if(Nflux >= MaxNflux)
                            {
                              Mesh.Indi.AllocFacNflux *= ALLOC_INCREASE_FACTOR;
                              MaxNflux = Mesh.Indi.AllocFacNflux;
#ifdef VERBOSE
                              printf("Task=%d: increase memory allocation, MaxNflux=%d Indi.AllocFacNflux=%g\n", ThisTask, MaxNflux,
                                     Mesh.Indi.AllocFacNflux);
#endif /* #ifdef VERBOSE */
                              FluxList = myrealloc_movable(FluxList, MaxNflux * sizeof(struct flux_list_data));

                              if(Nflux >= MaxNflux)
                                terminate("Nflux >= MaxNflux");
                            }

                          FluxList[Nflux].task  = Mesh.DP[q].task;
                          FluxList[Nflux].index = Mesh.DP[q].originalindex;
                          FluxList[Nflux].dM    = fac * P[i].Mass;
                          FluxList[Nflux].dP[0] = fac * SphP[i].Momentum[0];
                          FluxList[Nflux].dP[1] = fac * SphP[i].Momentum[1];
                          FluxList[Nflux].dP[2] = fac * SphP[i].Momentum[2];
#ifdef MHD
                          FluxList[Nflux].dB[0] = fac * SphP[i].BConserved[0];
                          FluxList[Nflux].dB[1] = fac * SphP[i].BConserved[1];
                          FluxList[Nflux].dB[2] = fac * SphP[i].BConserved[2];
#endif /* #ifdef MHD */

#ifndef ISOTHERM_EQS
                          FluxList[Nflux].dEnergy = fac * SphP[i].Energy;
#endif /* #ifndef ISOTHERM_EQS */

#ifdef MAXSCALARS
                          for(int s = 0; s < N_Scalar; s++)
                            FluxList[Nflux].dConservedScalars[s] =
                                fac * (*(MyFloat *)(((char *)(&SphP[i])) + scalar_elements[s].offset_mass));
#endif /* #ifdef MAXSCALARS */
                          Nflux++;
                        }
                    }

                  k = ngb[k].next_ngb;
                }

              if(fabs(facsum - 1) > 1.0e-3)
                {
                  char buf[1000];
                  sprintf(buf, "facsum=%g\n", facsum);
                  terminate(buf);
                }

              myfree(volume);
              myfree(trilist);

              /* we set the dissolved cell to zero mass and zero ID. It will be eliminated from the list
               * of cells in the next domain decomposition
               */
              P[i].Mass   = 0;
              P[i].ID     = 0;
              P[i].Vel[0] = 0;
              P[i].Vel[1] = 0;
              P[i].Vel[2] = 0;

              SphP[i].VelVertex[0] = 0;
              SphP[i].VelVertex[1] = 0;
              SphP[i].VelVertex[2] = 0;

              timebin_remove_particle(&TimeBinsHydro, idx, P[i].TimeBinHydro);

              voronoi_remove_connection(i);
            }
        }

      /* now let's apply the flux-list */
      derefine_apply_flux_list();
      myfree(FluxList);

      myfree(ngb);
      myfree(first_ngb);

#ifdef VERBOSE
      MPI_Reduce(&max_assumed_ntri, &n_tri, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
      if(ThisTask == 0)
        printf("DEREFINE: maximum assumed n_tri = %d\n", n_tri);
#endif /* #ifdef VERBOSE */
    }

  myfree(FlagExch);
  myfree(deref_SphP);

  /* remove removed cells from list of active gravity cells */
  timebin_cleanup_list_of_active_particles(&TimeBinsGravity);

  TIMER_STOP(CPU_DEREFINE);

  return countall;
}

/*! \brief Communicates probe list data if needed.
 *
 *  \return void
 */
static void derefine_apply_probe_list(void)
{
  int i, j, p, nimport, ngrp, recvTask;

  /* now exchange the probe-list and apply it where needed */

  mysort(ProbeList, Nflux, sizeof(struct probe_list_data), derefine_probe_list_data_compare_task);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < Nflux; i++)
    Send_count[ProbeList[i].task]++;

  if(Send_count[ThisTask] > 0)
    terminate("Send_count[ThisTask]");

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct probe_list_data *ProbeListGet = (struct probe_list_data *)mymalloc("ProbeListGet", nimport * sizeof(struct probe_list_data));

  /* exchange particle data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&ProbeList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct probe_list_data), MPI_BYTE,
                           recvTask, TAG_DENS_A, &ProbeListGet[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct probe_list_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  /* apply the probes */

  for(i = 0; i < nimport; i++)
    {
      p = ProbeListGet[i].index;

      if(TimeBinSynchronized[P[p].TimeBinHydro])
        if(deref_SphP[p].Flag == 2)
          ProbeListGet[i].flag = 1;
    }

  /* send results back */

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&ProbeListGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct probe_list_data), MPI_BYTE,
                           recvTask, TAG_DENS_A, &ProbeList[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct probe_list_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(ProbeListGet);
}

/*! \brief Communicate flux list data if needed.
 *
 *  \return void
 */
static void derefine_apply_flux_list(void)
{
  int i, j, p, nimport, ngrp, recvTask;

  /* now exchange the flux-list and apply it when needed */

  mysort(FluxList, Nflux, sizeof(struct flux_list_data), derefine_flux_list_data_compare);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < Nflux; i++)
    Send_count[FluxList[i].task]++;

  if(Send_count[ThisTask] > 0)
    terminate("Send_count[ThisTask]");

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct flux_list_data *FluxListGet = (struct flux_list_data *)mymalloc("FluxListGet", nimport * sizeof(struct flux_list_data));

  /* exchange particle data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&FluxList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct flux_list_data), MPI_BYTE, recvTask,
                           TAG_DENS_A, &FluxListGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct flux_list_data),
                           MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  /* apply the fluxes */

  for(i = 0; i < nimport; i++)
    {
      p = FluxListGet[i].index;

      if(P[p].ID == 0)
        {
          char buf[1000];
#ifndef LONGIDS
          printf("On task=%d flux to ID=%d, but this is already deleted (index p=%d)\n", ThisTask, P[p].ID, p);
#else  /* #ifndef LONGIDS */
          printf("On task=%d flux to ID=%llu, but this is already deleted (index p=%d)\n", ThisTask, P[p].ID, p);
#endif /* #ifndef LONGIDS #else */
          terminate(buf);
        }

      P[p].Mass += FluxListGet[i].dM;
      SphP[p].Momentum[0] += FluxListGet[i].dP[0];
      SphP[p].Momentum[1] += FluxListGet[i].dP[1];
      SphP[p].Momentum[2] += FluxListGet[i].dP[2];
#ifdef MHD
      SphP[p].BConserved[0] += FluxListGet[i].dB[0];
      SphP[p].BConserved[1] += FluxListGet[i].dB[1];
      SphP[p].BConserved[2] += FluxListGet[i].dB[2];
#endif /* #ifdef MHD */

#ifdef MAXSCALARS
      int k;
      for(k = 0; k < N_Scalar; k++)
        *(MyFloat *)(((char *)(&SphP[p])) + scalar_elements[k].offset_mass) += FluxListGet[i].dConservedScalars[k];
#endif /* #ifdef MAXSCALARS */

#ifndef ISOTHERM_EQS
      SphP[p].Energy += FluxListGet[i].dEnergy;
#endif /* #ifndef ISOTHERM_EQS */

#ifdef REFINEMENT_SPLIT_CELLS
      FlagDoNotRefine[p] = 1;
#endif /* #ifdef REFINEMENT_SPLIT_CELLS */
    }

  myfree(FluxListGet);
}

/*! \brief Compares flux list data task of two elements.
 *
 *  \param[in] a Pointer to first flux list data object.
 *  \param[in] b Pointer to second flux list data object.
 *
 *  \return (-1,0,1); -1 if a->task < b->task.
 */
static int derefine_flux_list_data_compare(const void *a, const void *b)
{
  if(((struct flux_list_data *)a)->task < (((struct flux_list_data *)b)->task))
    return -1;

  if(((struct flux_list_data *)a)->task > (((struct flux_list_data *)b)->task))
    return +1;

  return 0;
}

/*! \brief Compares probe list data task of two elements.
 *
 *  \param[in] a Pointer to first probe list data object.
 *  \param[in] b Pointer to second probe list data object.
 *
 *  \return (-1,0,1); -1 if a->task < b->task.
 */
static int derefine_probe_list_data_compare_task(const void *a, const void *b)
{
  if(((struct probe_list_data *)a)->task < (((struct probe_list_data *)b)->task))
    return -1;

  if(((struct probe_list_data *)a)->task > (((struct probe_list_data *)b)->task))
    return +1;

  return 0;
}

/*! \brief Compares seq delaunay data task of two elements.
 *
 *  Comparison criteria (most important first)
 *    rnd
 *    ID
 *    index
 *    rank
 *
 *  \param[in] a Pointer to first seq delaunay data object.
 *  \param[in] b Pointer to second seq delaunay data object.
 *
 *  \return (-1,0,1); -1 if a < b.
 */
static int derefine_compare_seq_DP_ID(const void *a, const void *b)
{
  if(((struct seq_delaunay_data *)a)->rnd < (((struct seq_delaunay_data *)b)->rnd))
    return -1;

  if(((struct seq_delaunay_data *)a)->rnd > (((struct seq_delaunay_data *)b)->rnd))
    return +1;

  if(((struct seq_delaunay_data *)a)->ID < (((struct seq_delaunay_data *)b)->ID))
    return -1;

  if(((struct seq_delaunay_data *)a)->ID > (((struct seq_delaunay_data *)b)->ID))
    return +1;

  if(((struct seq_delaunay_data *)a)->index < (((struct seq_delaunay_data *)b)->index))
    return -1;

  if(((struct seq_delaunay_data *)a)->index > (((struct seq_delaunay_data *)b)->index))
    return +1;

  if(((struct seq_delaunay_data *)a)->rank < (((struct seq_delaunay_data *)b)->rank))
    return -1;

  if(((struct seq_delaunay_data *)a)->rank > (((struct seq_delaunay_data *)b)->rank))
    return +1;

  return 0;
}

/*! \brief Sets exchange flag in de-refinement algorithm.
 *
 *  Loops through gas cells in mesh, sets set export flag and communicates this
 *  information to the appropriate tasks.
 *
 *  \return void
 */
static void derefine_exchange_flag(void)
{
  int listp;
  int i, j, p, task, off;
  int ngrp, recvTask, place;

  struct exchange_data
  {
    int Flag;
    MyIDType ID;
  } * tmpExch, *tmpRecv;

  tmpExch = (struct exchange_data *)mymalloc("tmpExch", Mesh_nexport * sizeof(struct exchange_data));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off   = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

              tmpExch[off].Flag = 0;
              tmpExch[off].ID   = P[place].ID;

              if(P[place].Type == 0)
                if(TimeBinSynchronized[P[place].TimeBinHydro])
                  if(!(P[place].Mass == 0 && P[place].ID == 0))
                    tmpExch[off].Flag = deref_SphP[place].Flag;
            }
          listp = ListExports[listp].nextexport;
        }
    }

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              tmpRecv = (struct exchange_data *)mymalloc("tmpRecv", Mesh_Recv_count[recvTask] * sizeof(struct exchange_data));

              /* get the values */
              MPI_Sendrecv(&tmpExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] * sizeof(struct exchange_data), MPI_BYTE,
                           recvTask, TAG_DENS_A, tmpRecv, Mesh_Recv_count[recvTask] * sizeof(struct exchange_data), MPI_BYTE, recvTask,
                           TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              for(i = 0; i < Mesh_Recv_count[recvTask]; i++)
                {
                  if(Mesh_Recv_offset[recvTask] + i >= Mesh_nimport)
                    terminate("number of imported mesh points grater than Mesh_nimport");
                  FlagExch[Mesh_Recv_offset[recvTask] + i].Flag = tmpRecv[i].Flag;
                  FlagExch[Mesh_Recv_offset[recvTask] + i].ID   = tmpRecv[i].ID;
                }

              myfree(tmpRecv);
            }
        }
    }

  myfree(tmpExch);
}

#endif /* #if defined(REFINEMENT_MERGE_CELLS) && !defined(ONEDIMS) */
