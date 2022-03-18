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
 * \file        src/subfind/subfind_collective.c
 * \date        05/2018
 * \brief       Subfind algorithm running collectively on all tasks.
 * \details     contains functions:
 *                void subfind_process_group_collectively(int nsubgroups_cat)
 *                void subfind_fof_calc_am_collective(int snapnr, int
 *                  ngroups_cat)
 *                void subfind_col_find_coll_candidates(int totgrouplen)
 *                void subfind_unbind_independent_ones(int count_cand)
 *                int subfind_col_unbind(struct unbind_data *d, int num, int
 *                  *num_non_gas)
 *                void subfind_poll_for_requests(void)
 *                long long subfind_distlinklist_setrank_and_get_next(
 *                  long long index, long long *rank)
 *                void subfind_distlinklist_set_next(long long index,
 *                  long long next)
 *                void subfind_distlinklist_add_particle(long long index)
 *                void subfind_distlinklist_mark_particle(long long index,
 *                  int target, int submark)
 *                void subfind_distlinklist_add_bound_particles(
 *                  long long index, int nsub)
 *                long long subfind_distlinklist_get_next(long long index)
 *                long long subfind_distlinklist_get_rank(long long index)
 *                long long subfind_distlinklist_get_head(long long index)
 *                void subfind_distlinklist_get_two_heads(long long ngb_index1,
 *                  long long ngb_index2, long long *head, long long
 *                  *head_attach)
 *                void subfind_distlinklist_set_headandnext(long long index,
 *                  long long head, long long next)
 *                int subfind_distlinklist_get_tail_set_tail_increaselen(
 *                  long long index, long long *tail, long long newtail)
 *                void subfind_distlinklist_set_tailandlen(long long index,
 *                  long long tail, int len)
 *                void subfind_distlinklist_get_tailandlen(long long index,
 *                  long long *tail, int *len)
 *                void subfind_distlinklist_set_all(long long index,
 *                  long long head, long long tail, int len, long long next)
 *                int subfind_compare_densities(const void *a, const void *b)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 15.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef SUBFIND
#include "../fof/fof.h"
#include "subfind.h"

#define TAG_POLLING_DONE 201
#define TAG_SET_ALL 202
#define TAG_GET_NGB_INDICES 204
#define TAG_GET_TAILANDLEN 205
#define TAG_GET_TAILANDLEN_DATA 206
#define TAG_SET_TAILANDLEN 207
#define TAG_SET_HEADANDNEXT 209
#define TAG_SETHEADGETNEXT_DATA 210
#define TAG_SET_NEXT 211
#define TAG_SETHEADGETNEXT 213
#define TAG_GET_NEXT 215
#define TAG_GET_NEXT_DATA 216
#define TAG_GET_HEAD 217
#define TAG_GET_HEAD_DATA 218
#define TAG_ADD_PARTICLE 219
#define TAG_ADDBOUND 220
#define TAG_NID 222
#define TAG_NID_DATA 223
#define TAG_SETRANK 224
#define TAG_SETRANK_OUT 226
#define TAG_GET_RANK 227
#define TAG_GET_RANK_DATA 228
#define TAG_MARK_PARTICLE 229
#define TAG_SET_NEWTAIL 230
#define TAG_GET_OLDTAIL 231
#define TAG_GET_TWOHEADS 232
#define TAG_GET_TWOHEADS_DATA 233

#define MASK ((((long long)1) << 32) - 1)
#define HIGHBIT (1 << 30)

static long long *Head, *Next, *Tail;
static int *Len;
static int LocalLen;
static int count_cand, max_coll_candidates;

static struct unbind_data *ud;

/*! \brief Data structure for sorting density data.
 */
static struct sort_density_data
{
  MyFloat density;
  int ngbcount;
  long long index; /* this will store the task in the upper word */
  long long ngb_index1, ngb_index2;
} * sd;

/*! \brief Processes a group collectively on all MPI tasks.
 *
 *  \param[in] nsubgroups_cat (unused)
 *
 *  \return void
 */
void subfind_process_group_collectively(int nsubgroups_cat)
{
  int totgrouplen1, totgrouplen2;

  /* make a sanity check: We should have exactly 1 group, stored on the root of the processor subset */
  if(SubThisTask == 0)
    {
      if(Ngroups != 1)
        terminate("Ngroups=%d != 1  SubNTask=%d SubThisTask=%d", Ngroups, SubNTask, SubThisTask);
    }
  else
    {
      if(Ngroups != 0)
        terminate("Ngroups=%d != 0  SubNTask=%d SubThisTask=%d", Ngroups, SubNTask, SubThisTask);
    }

  if(SubThisTask == 0)
    {
      printf("SUBFIND-COLLECTIVE, root-task=%d: Collectively doing halo %d of length  %d  on  %d  processors.\n", ThisTask,
             Group[0].GrNr, Group[0].Len, SubNTask);

      GrNr         = Group[0].GrNr;
      totgrouplen2 = Group[0].Len;
      for(int j = 0; j < 3; j++)
        GrCM[j] = Group[0].CM[j];
    }

  /* tell everybody in the set the group number, the center of mass, and the grouplen */
  MPI_Bcast(&GrNr, 1, MPI_INT, 0, SubComm);
  MPI_Bcast(&GrCM[0], 3 * sizeof(MyDouble), MPI_BYTE, 0, SubComm);
  MPI_Bcast(&totgrouplen2, 1, MPI_INT, 0, SubComm);

  NumPartGroup = 0;
  for(int i = 0; i < NumPart; i++)
    if(PS[i].GrNr == GrNr)
      NumPartGroup++;

  MPI_Allreduce(&NumPartGroup, &totgrouplen1, 1, MPI_INT, MPI_SUM, SubComm);

  /* sanity check that we actually have all the right particles on the processor subset */
  if(totgrouplen1 != totgrouplen2)
    terminate("totgrouplen1=%d != totgrouplen2=%d", totgrouplen1, totgrouplen2); /* inconsistency */

  /* do a domain decomposition just for this halo */
  subfind_coll_domain_decomposition();

  /* copy over the domain dimensions to serial tree code, as this may be used in the collective unbinding */
  subfind_loctree_copyExtent();

  /* now let us sort according to GrNr and Density. This step will temporarily break the association with SphP[] and other arrays! */
  submp = (struct submp_data *)mymalloc("submp", sizeof(struct submp_data) * NumPart);
  for(int i = 0; i < NumPart; i++)
    {
      PS[i].SubNr         = TotNgroups + 1; /* set a default that is larger than reasonable group number */
      PS[i].OldIndex      = i;
      submp[i].index      = i;
      submp[i].GrNr       = PS[i].GrNr;
      submp[i].DM_Density = PS[i].Density;
    }
  qsort(submp, NumPart, sizeof(struct submp_data), subfind_compare_submp_GrNr_DM_Density);
  subfind_reorder_according_to_submp();
  myfree(submp);

  /* note: now we have the particles of the group at the beginning, but SPH particles are not aligned.
     They can however be accessed via SphP[PS[i].OldIndex] */

  /* re-determine the number of local group particles, which has changed due to domain decomposition */
  NumPartGroup = 0;
  for(int i = 0; i < NumPart; i++)
    if(PS[i].GrNr == GrNr)
      NumPartGroup++;

  /* allocate some storage for the halo */
  subfind_coll_treeallocate(NumPart, All.MaxPart);

  /* construct a tree for the halo */
  subfind_coll_treebuild(NumPartGroup, NULL);

#ifdef SUBFIND_EXTENDED_PROPERTIES
  // calculate binding energy of full fof group
  {
    struct unbind_data *ud = (struct unbind_data *)mymalloc_movable(&ud, "ud", NumPartGroup * sizeof(struct unbind_data));

    NumPartGroup = 0;
    for(int i = 0; i < NumPart; i++)
      if(PS[i].GrNr == GrNr)
        ud[NumPartGroup++].index = i;

    subfind_potential_compute(NumPartGroup, ud, 0, 0);

    double binding_energy_local = 0, binding_energy_global;

    for(int i = 0; i < NumPartGroup; i++)
      binding_energy_local += 0.5 * P[ud[i].index].Mass * PS[ud[i].index].Potential;

    MPI_Allreduce(&binding_energy_local, &binding_energy_global, 1, MPI_DOUBLE, MPI_SUM, SubComm);
    Group[0].Epot = binding_energy_global;

    myfree(ud);
    ud = NULL;
  }
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

  long long p;
  int len;
  int ncand, parent, totcand, nremaining;
  int max_loc_length, max_length;
  int count, countall, *countlist, *offset;
  int i, j, k, nr, grindex = 0, nsubs, subnr;
  int count_leaves, tot_count_leaves, master;
  struct coll_cand_dat *tmp_coll_candidates = 0;
  double t0, t1, tt0, tt1;

  /* determine the radius that encloses a certain number of link particles */
  subfind_find_linkngb();

  sd = (struct sort_density_data *)mymalloc_movable(&sd, "sd", NumPartGroup * sizeof(struct sort_density_data));

  /* determine the indices of the nearest two denser neighbours within the link region */
  NgbLoc = (struct nearest_ngb_data *)mymalloc("NgbLoc", NumPartGroup * sizeof(struct nearest_ngb_data));
  R2Loc  = (struct nearest_r2_data *)mymalloc("R2Loc", NumPartGroup * sizeof(struct nearest_r2_data));

  subfind_find_nearesttwo();

  for(i = 0; i < NumPartGroup; i++)
    {
      sd[i].density    = PS[i].Density;
      sd[i].ngbcount   = NgbLoc[i].count;
      sd[i].index      = (((long long)SubThisTask) << 32) + i;
      sd[i].ngb_index1 = NgbLoc[i].index[0];
      sd[i].ngb_index2 = NgbLoc[i].index[1];
    }
  myfree(R2Loc);
  myfree(NgbLoc);

  if(SubThisTask == 0)
    {
      printf("SUBFIND-COLLECTIVE, root-task=%d: before parallel sort of 'sd'.\n", ThisTask);
      fflush(stdout);
    }

  /* sort the densities */
  parallel_sort_comm(sd, NumPartGroup, sizeof(struct sort_density_data), subfind_compare_densities, SubComm);

  if(SubThisTask == 0)
    {
      printf("SUBFIND-COLLECTIVE, root-task=%d: parallel sort of 'sd' done.\n", ThisTask);
      fflush(stdout);
    }

  /* allocate and initialize distributed link list */
  Head = (long long *)mymalloc_movable(&Head, "Head", NumPartGroup * sizeof(long long));
  Next = (long long *)mymalloc_movable(&Next, "Next", NumPartGroup * sizeof(long long));
  Tail = (long long *)mymalloc_movable(&Tail, "Tail", NumPartGroup * sizeof(long long));
  Len  = (int *)mymalloc_movable(&Len, "Len", NumPartGroup * sizeof(int));

  for(i = 0; i < NumPartGroup; i++)
    {
      Head[i] = Next[i] = Tail[i] = -1;
      Len[i]                      = 0;
    }

  /* allocate a list to store subhalo coll_candidates */
  max_coll_candidates = imax((NumPartGroup / 50), 200);
  coll_candidates     = (struct coll_cand_dat *)mymalloc_movable(&coll_candidates, "coll_candidates",
                                                             max_coll_candidates * sizeof(struct coll_cand_dat));
  count_cand          = 0;

  subfind_col_find_coll_candidates(totgrouplen1);

  /* establish total number of coll_candidates */
  MPI_Allreduce(&count_cand, &totcand, 1, MPI_INT, MPI_SUM, SubComm);
  if(SubThisTask == 0)
    {
      printf("SUBFIND-COLLECTIVE, root-task=%d: total number of subhalo coll_candidates=%d\n", ThisTask, totcand);
      fflush(stdout);
    }

  nremaining = totcand;

  for(i = 0; i < NumPartGroup; i++)
    Tail[i] = -1;

  for(i = 0; i < count_cand; i++)
    coll_candidates[i].parent = 0;

  do
    {
      /* Let's see which coll_candidates can be unbound independent from each other.
         We identify them with those coll_candidates that have no embedded other candidate */
      t0 = second();
      if(SubThisTask == 0)
        tmp_coll_candidates = (struct coll_cand_dat *)mymalloc("tmp_coll_candidates", totcand * sizeof(struct coll_cand_dat));

      count = count_cand;
      count *= sizeof(struct coll_cand_dat);

      countlist = (int *)mymalloc("countlist", SubNTask * sizeof(int));
      offset    = (int *)mymalloc("offset", SubNTask * sizeof(int));

      MPI_Allgather(&count, 1, MPI_INT, countlist, 1, MPI_INT, SubComm);

      for(i = 1, offset[0] = 0; i < SubNTask; i++)
        offset[i] = offset[i - 1] + countlist[i - 1];

      MPI_Gatherv(coll_candidates, countlist[SubThisTask], MPI_BYTE, tmp_coll_candidates, countlist, offset, MPI_BYTE, 0, SubComm);

      if(SubThisTask == 0)
        {
          for(k = 0; k < totcand; k++)
            {
              tmp_coll_candidates[k].nsub  = k;
              tmp_coll_candidates[k].subnr = k;
            }

          qsort(tmp_coll_candidates, totcand, sizeof(struct coll_cand_dat), subfind_compare_coll_candidates_rank);
          for(k = 0; k < totcand; k++)
            {
              if(tmp_coll_candidates[k].parent >= 0)
                {
                  tmp_coll_candidates[k].parent = 0;

                  for(j = k + 1; j < totcand; j++)
                    {
                      if(tmp_coll_candidates[j].rank > tmp_coll_candidates[k].rank + tmp_coll_candidates[k].len)
                        break;

                      if(tmp_coll_candidates[j].parent < 0) /* ignore these */
                        continue;

                      if(tmp_coll_candidates[k].rank + tmp_coll_candidates[k].len >=
                         tmp_coll_candidates[j].rank + tmp_coll_candidates[j].len)
                        {
                          tmp_coll_candidates[k].parent++; /* we here count the number of subhalos that are enclosed */
                        }
                      else
                        {
                          terminate("k=%d|%d has rank=%d and len=%d.  j=%d has rank=%d and len=%d\n", k, totcand,
                                    (int)tmp_coll_candidates[k].rank, (int)tmp_coll_candidates[k].len, j,
                                    (int)tmp_coll_candidates[j].rank, (int)tmp_coll_candidates[j].len);
                        }
                    }
                }
            }

          qsort(tmp_coll_candidates, totcand, sizeof(struct coll_cand_dat), subfind_compare_coll_candidates_subnr);
        }

      MPI_Scatterv(tmp_coll_candidates, countlist, offset, MPI_BYTE, coll_candidates, countlist[SubThisTask], MPI_BYTE, 0, SubComm);

      myfree(offset);
      myfree(countlist);

      if(SubThisTask == 0)
        myfree(tmp_coll_candidates);

      for(i = 0, count_leaves = 0, max_loc_length = 0; i < count_cand; i++)
        if(coll_candidates[i].parent == 0)
          {
            if(coll_candidates[i].len > max_loc_length)
              max_loc_length = coll_candidates[i].len;

            if(coll_candidates[i].len > 0.20 * All.TotNumPart / NTask) /* seems large, let's rather do it collectively */
              {
                coll_candidates[i].parent++; /* this will ensure that it is not considered in this round */
              }
            else
              {
                count_leaves++;
              }
          }

      MPI_Allreduce(&count_leaves, &tot_count_leaves, 1, MPI_INT, MPI_SUM, SubComm);
      MPI_Allreduce(&max_loc_length, &max_length, 1, MPI_INT, MPI_MAX, SubComm);

      t1 = second();
      if(SubThisTask == 0)
        printf(
            "SUBFIND-COLLECTIVE, root-task=%d: number of subhalo coll_candidates that can be done independently=%d. (Largest size %d, "
            "finding took %g sec)\n",
            ThisTask, tot_count_leaves, max_length, timediff(t0, t1));

      if(tot_count_leaves <= 0) /* if there are none left, we break and do the reset collectively */
        {
          if(SubThisTask == 0)
            printf("SUBFIND-COLLECTIVE, root-task=%d: too few, I do the rest of %d collectively\n", ThisTask, nremaining);
          break;
        }

      nremaining -= tot_count_leaves;

      for(i = 0; i < NumPart; i++)
        {
          PS[i].origintask = PS[i].TargetTask = SubThisTask;
          PS[i].originindex                   = i;
          PS[i].submark                       = HIGHBIT;
          if(i < NumPartGroup)
            if(Tail[i] >= 0) /* this means this particle is already bound to a substructure */
              PS[i].origintask |= HIGHBIT;
        }

      /* we now mark the particles that are in subhalo coll_candidates that can be processed independently in parallel */
      nsubs = 0;
      t0    = second();
      for(master = 0; master < SubNTask; master++)
        {
          ncand = count_cand;

          MPI_Bcast(&ncand, sizeof(ncand), MPI_BYTE, master, SubComm);

          for(k = 0; k < ncand; k++)
            {
              if(SubThisTask == master)
                {
                  len    = coll_candidates[k].len;
                  parent = coll_candidates[k].parent; /* this is here actually the daughter count */
                }

              MPI_Bcast(&len, sizeof(len), MPI_BYTE, master, SubComm);
              MPI_Bcast(&parent, sizeof(parent), MPI_BYTE, master, SubComm);
              MPI_Barrier(SubComm);

              if(parent == 0)
                {
                  if(SubThisTask != master)
                    subfind_poll_for_requests();
                  else
                    {
                      for(i = 0, p = coll_candidates[k].head; i < coll_candidates[k].len; i++)
                        {
                          subfind_distlinklist_mark_particle(p, master, nsubs);

                          if(p < 0)
                            terminate("Bummer i=%d \n", i);

                          p = subfind_distlinklist_get_next(p);
                        }

                      /* now tell the others to stop polling */
                      for(i = 0; i < SubNTask; i++)
                        if(i != SubThisTask)
                          MPI_Send(&i, 1, MPI_INT, i, TAG_POLLING_DONE, SubComm);
                    }

                  MPI_Barrier(SubComm);
                }

              nsubs++;
            }
        }
      t1 = second();
      if(SubThisTask == 0)
        {
          printf("SUBFIND-COLLECTIVE, root-task=%d: particles are marked (took %g)\n", ThisTask, timediff(t0, t1));
          fflush(stdout);
        }

      for(i = 0; i < NumPart; i++)
        PS[i].TargetIndex = PS[i].submark; /* this will make sure that the particles are grouped by submark on the target task */

      t0 = second();
      subfind_distribute_particles(SubComm); /* assemble the particles on individual processors */
      t1 = second();
      if(SubThisTask == 0)
        {
          printf("SUBFIND-COLLECTIVE, root-task=%d: distribution of independent ones took %g sec\n", ThisTask, timediff(t0, t1));
          fflush(stdout);
        }

      MPI_Barrier(SubComm);
      t0 = second();

      subfind_unbind_independent_ones(count_cand);

      MPI_Barrier(SubComm);
      t1 = second();

      if(SubThisTask == 0)
        {
          printf("SUBFIND-COLLECTIVE, root-task=%d: unbinding of independent ones took %g sec\n", ThisTask, timediff(t0, t1));
          fflush(stdout);
        }

      for(i = 0; i < NumPart; i++)
        {
          PS[i].origintask &= (HIGHBIT - 1); /* clear high bit if set */
          PS[i].TargetTask  = PS[i].origintask;
          PS[i].TargetIndex = PS[i].originindex;
        }

      t0 = second();
      subfind_distribute_particles(SubComm); /* bring them back to their original processor */

      t1 = second();
      if(SubThisTask == 0)
        {
          printf("SUBFIND-COLLECTIVE, root-task=%d: bringing the independent ones back took %g sec\n", ThisTask, timediff(t0, t1));
          fflush(stdout);
        }

      /* now mark the bound particles */
      for(i = 0; i < NumPartGroup; i++)
        if(PS[i].submark >= 0 && PS[i].submark < nsubs)
          Tail[i] = PS[i].submark; /* we use this to flag bound parts of substructures */

      for(i = 0; i < count_cand; i++)
        if(coll_candidates[i].parent == 0)
          coll_candidates[i].parent = -1;
    }
  while(tot_count_leaves > 0);

  /**** now we do the collective unbinding of the subhalo coll_candidates that contain other subhalo coll_candidates ****/
  ud = (struct unbind_data *)mymalloc_movable(&ud, "ud", NumPartGroup * sizeof(struct unbind_data));

  t0 = second();
  for(master = 0, nr = 0; master < SubNTask; master++)
    {
      ncand = count_cand;

      MPI_Bcast(&ncand, sizeof(ncand), MPI_BYTE, master, SubComm);

      for(k = 0; k < ncand; k++)
        {
          if(SubThisTask == master)
            {
              len    = coll_candidates[k].len;
              nsubs  = coll_candidates[k].nsub;
              parent = coll_candidates[k].parent; /* this is here actually the daughter count */
            }

          MPI_Bcast(&parent, sizeof(parent), MPI_BYTE, master, SubComm);
          MPI_Barrier(SubComm);

          if(parent >= 0)
            {
              MPI_Bcast(&len, sizeof(len), MPI_BYTE, master, SubComm);
              MPI_Bcast(&nsubs, sizeof(nsubs), MPI_BYTE, master, SubComm);

              if(SubThisTask == 0)
                {
                  printf("SUBFIND-COLLECTIVE, root-task=%d: collective unbinding of nr=%d (%d) of length=%d\n", ThisTask, nr,
                         nremaining, (int)len);
                  fflush(stdout);
                }

              nr++;

              LocalLen = 0;

              tt0 = second();

              if(SubThisTask != master)
                subfind_poll_for_requests();
              else
                {
                  for(i = 0, p = coll_candidates[k].head; i < coll_candidates[k].len; i++)
                    {
                      subfind_distlinklist_add_particle(p);
                      if(p < 0)
                        terminate("Bummer i=%d \n", i);

                      p = subfind_distlinklist_get_next(p);
                    }

                  /* now tell the others to stop polling */
                  for(i = 0; i < SubNTask; i++)
                    if(i != SubThisTask)
                      MPI_Send(&i, 1, MPI_INT, i, TAG_POLLING_DONE, SubComm);
                }

              int LocalNonGasLen;

              LocalLen = subfind_col_unbind(ud, LocalLen, &LocalNonGasLen);

              tt1 = second();
              if(SubThisTask == 0)
                {
                  printf("SUBFIND-COLLECTIVE, root-task=%d: took %g sec\n", ThisTask, timediff(tt0, tt1));
                  fflush(stdout);
                }

              MPI_Allreduce(&LocalLen, &len, 1, MPI_INT, MPI_SUM, SubComm);

              if(len >= All.DesLinkNgb)
                {
                  /* ok, we found a substructure */

                  for(i = 0; i < LocalLen; i++)
                    Tail[ud[i].index] = nsubs; /* we use this to flag the substructures */

                  if(SubThisTask == master)
                    {
                      coll_candidates[k].bound_length = len;
                    }
                }
              else
                {
                  if(SubThisTask == master)
                    {
                      coll_candidates[k].bound_length = 0;
                    }
                }
            }
        }
    }
  t1 = second();

  if(SubThisTask == 0)
    {
      printf("SUBFIND-COLLECTIVE, root-task=%d: the collective unbinding of remaining halos took %g sec\n", ThisTask,
             timediff(t0, t1));
      fflush(stdout);
    }

  for(k = 0, count = 0; k < count_cand; k++)
    if(coll_candidates[k].bound_length >= All.DesLinkNgb)
      {
        if(coll_candidates[k].len < All.DesLinkNgb)
          terminate("coll_candidates[k=%d].len=%d bound=%d\n", k, coll_candidates[k].len, coll_candidates[k].bound_length);

        count++;
      }

  MPI_Allreduce(&count, &countall, 1, MPI_INT, MPI_SUM, SubComm);

  if(SubThisTask == 0)
    {
      printf("SUBFIND-COLLECTIVE, root-task=%d: found %d bound substructures in FoF group of length %d\n", ThisTask, countall,
             totgrouplen1);
      fflush(stdout);
    }

  /* now determine the parent subhalo for each candidate */
  t0 = second();
  parallel_sort_comm(coll_candidates, count_cand, sizeof(struct coll_cand_dat), subfind_compare_coll_candidates_boundlength, SubComm);

  if(SubThisTask == 0)
    tmp_coll_candidates = (struct coll_cand_dat *)mymalloc("tmp_coll_candidates", totcand * sizeof(struct coll_cand_dat));

  count = count_cand;
  count *= sizeof(struct coll_cand_dat);

  countlist = (int *)mymalloc("countlist", SubNTask * sizeof(int));
  offset    = (int *)mymalloc("offset", SubNTask * sizeof(int));

  MPI_Allgather(&count, 1, MPI_INT, countlist, 1, MPI_INT, SubComm);

  for(i = 1, offset[0] = 0; i < SubNTask; i++)
    offset[i] = offset[i - 1] + countlist[i - 1];

  MPI_Gatherv(coll_candidates, countlist[SubThisTask], MPI_BYTE, tmp_coll_candidates, countlist, offset, MPI_BYTE, 0, SubComm);

  if(SubThisTask == 0)
    {
      for(k = 0; k < totcand; k++)
        {
          tmp_coll_candidates[k].subnr  = k;
          tmp_coll_candidates[k].parent = 0;
        }

      qsort(tmp_coll_candidates, totcand, sizeof(struct coll_cand_dat), subfind_compare_coll_candidates_rank);

      for(k = 0; k < totcand; k++)
        {
          for(j = k + 1; j < totcand; j++)
            {
              if(tmp_coll_candidates[j].rank > tmp_coll_candidates[k].rank + tmp_coll_candidates[k].len)
                break;

              if(tmp_coll_candidates[k].rank + tmp_coll_candidates[k].len >= tmp_coll_candidates[j].rank + tmp_coll_candidates[j].len)
                {
                  if(tmp_coll_candidates[k].bound_length >= All.DesLinkNgb)
                    tmp_coll_candidates[j].parent = tmp_coll_candidates[k].subnr;
                }
              else
                {
                  terminate("k=%d|%d has rank=%d and len=%d.  j=%d has rank=%d and len=%d bound=%d\n", k, countall,
                            (int)tmp_coll_candidates[k].rank, (int)tmp_coll_candidates[k].len,
                            (int)tmp_coll_candidates[k].bound_length, (int)tmp_coll_candidates[j].rank,
                            (int)tmp_coll_candidates[j].len, (int)tmp_coll_candidates[j].bound_length);
                }
            }
        }

      qsort(tmp_coll_candidates, totcand, sizeof(struct coll_cand_dat), subfind_compare_coll_candidates_subnr);
    }

  MPI_Scatterv(tmp_coll_candidates, countlist, offset, MPI_BYTE, coll_candidates, countlist[SubThisTask], MPI_BYTE, 0, SubComm);

  myfree(offset);
  myfree(countlist);

  if(SubThisTask == 0)
    myfree(tmp_coll_candidates);

  t1 = second();
  if(SubThisTask == 0)
    {
      printf("SUBFIND-COLLECTIVE, root-task=%d: determination of parent subhalo took %g sec (presently allocated %g MB)\n", ThisTask,
             timediff(t0, t1), AllocatedBytes / (1024.0 * 1024.0));
      fflush(stdout);
    }

  /* Now let's save some properties of the substructures */
  if(SubThisTask == 0)
    Group[0].Nsubs = countall;

  t0 = second();
  for(master = 0, subnr = 0; master < SubNTask; master++)
    {
      ncand = count_cand;
      MPI_Bcast(&ncand, sizeof(ncand), MPI_BYTE, master, SubComm);

      for(k = 0; k < ncand; k++)
        {
          if(SubThisTask == master)
            {
              len    = coll_candidates[k].bound_length;
              nsubs  = coll_candidates[k].nsub;
              parent = coll_candidates[k].parent;
            }

          MPI_Bcast(&len, sizeof(len), MPI_BYTE, master, SubComm);
          MPI_Barrier(SubComm);

          if(len > 0)
            {
              MPI_Bcast(&nsubs, sizeof(nsubs), MPI_BYTE, master, SubComm);
              MPI_Bcast(&parent, sizeof(parent), MPI_BYTE, master, SubComm);

              LocalLen = 0;

              if(SubThisTask != master)
                subfind_poll_for_requests();
              else
                {
                  for(i = 0, p = coll_candidates[k].head; i < coll_candidates[k].len; i++)
                    {
                      subfind_distlinklist_add_bound_particles(p, nsubs);
                      p = subfind_distlinklist_get_next(p);
                    }

                  /* now tell the others to stop polling */
                  for(i = 0; i < SubNTask; i++)
                    if(i != SubThisTask)
                      MPI_Send(&i, 1, MPI_INT, i, TAG_POLLING_DONE, SubComm);
                }

              MPI_Barrier(SubComm);

              if(SubThisTask == 0)
                {
                  if(Nsubgroups >= MaxNsubgroups)
                    terminate("Nsubgroups=%d >= MaxNsubgroups=%d", Nsubgroups, MaxNsubgroups);
                }

              tt0 = second();
              subfind_determine_sub_halo_properties(ud, LocalLen, &SubGroup[Nsubgroups], GrNr, subnr, 1, nsubgroups_cat);
              tt1 = second();

              /* we have filled into ud the binding energy and the particle ID return */

              if(SubThisTask == 0)
                {
                  if(Nsubgroups >= MaxNsubgroups)
                    terminate("Nsubgroups >= MaxNsubgroups");

                  if(subnr == 0)
                    {
                      for(j = 0; j < 3; j++)
                        Group[grindex].Pos[j] = SubGroup[Nsubgroups].Pos[j];
                    }

                  SubGroup[Nsubgroups].GrNr      = GrNr;
                  SubGroup[Nsubgroups].SubNr     = subnr;
                  SubGroup[Nsubgroups].SubParent = parent;

                  Nsubgroups++;
                }

              /* Let's now assign the subgroup number */
              for(i = 0; i < LocalLen; i++)
                PS[ud[i].index].SubNr = subnr;

              subnr++;
            }
        }
    }

  t1 = second();
  if(SubThisTask == 0)
    {
      printf("SUBFIND-COLLECTIVE, root-task=%d: determining substructure properties took %g sec (presently allocated %g MB)\n",
             ThisTask, timediff(t0, t1), AllocatedBytes / (1024.0 * 1024.0));
      fflush(stdout);
    }

  myfree(ud);
  ud = NULL;
  myfree(coll_candidates);
  myfree(Len);
  myfree(Tail);
  myfree(Next);
  myfree(Head);
  myfree(sd);

  subfind_coll_treefree();
  subfind_coll_domain_free();

  /* undo local rearrangement that made group consecutive. After that, the association of SphP[] will be correct again */
  submp = (struct submp_data *)mymalloc("submp", sizeof(struct submp_data) * NumPart);
  for(int i = 0; i < NumPart; i++)
    {
      submp[i].index    = i;
      submp[i].OldIndex = PS[i].OldIndex;
    }
  qsort(submp, NumPart, sizeof(struct submp_data), subfind_compare_submp_OldIndex);
  subfind_reorder_according_to_submp();
  myfree(submp);
}

#ifdef SUBFIND_EXTENDED_PROPERTIES
/*! \brief Calculates angualar momentum collectively on all MPI tasks.
 *
 *  \param[in] snapnr (unused)
 *  \param[in] ngroups_cat (unused)
 *
 *  \return void
 */
void subfind_fof_calc_am_collective(int snapnr, int ngroups_cat)
{
  int len, totgrouplen1, totgrouplen2;
  long long index;

  int grindex = 0, i, k, ptype;
  double Pos_pbc[3], Vel_tot[3], gr_pos[3], gr_vel[3];
  double gr_Jtot[3], gr_Jdm[3], gr_Jgas[3], gr_Jstars[3], jpart[3];
  double gr_CMFrac, gr_CMFracType[NTYPES];
  int gr_len_dm;
  double gr_mass, gr_mass_gas, gr_mass_stars;  // gr_mass_dm,
  double gr_Ekin, gr_Ethr;

  /* make a sanity check: We should have exactly 1 group, stored on the root of the processor subset */
  if(SubThisTask == 0)
    {
      if(Ngroups != 1)
        terminate("Ngroups=%d != 1  SubNTask=%d SubThisTask=%d", Ngroups, SubNTask, SubThisTask);
    }
  else
    {
      if(Ngroups != 0)
        terminate("Ngroups=%d != 0  SubNTask=%d SubThisTask=%d", Ngroups, SubNTask, SubThisTask);
    }

  if(SubThisTask == 0)
    {
      printf("SUBFIND-COLLECTIVE, root-task=%d: Collectively doing AM of halo %d of length %d on %d processors.\n", ThisTask,
             Group[0].GrNr, Group[0].Len, SubNTask);

      totgrouplen2 = Group[0].Len;
    }

  /* tell everybody in the set the group number and the grouplen */
  MPI_Bcast(&GrNr, 1, MPI_INT, 0, SubComm);
  MPI_Bcast(&totgrouplen2, 1, MPI_INT, 0, SubComm);

  for(i = 0, NumPartGroup = 0; i < NumPart; i++)
    if(PS[i].GrNr == GrNr)
      NumPartGroup++;

  MPI_Allreduce(&NumPartGroup, &totgrouplen1, 1, MPI_INT, MPI_SUM, SubComm);

  /* sanity check that we actually have all the right particles on the processor subset */
  if(totgrouplen1 != totgrouplen2)
    terminate("totgrouplen1 != totgrouplen2"); /* inconsistency */

  /* do a domain decomposition just for this halo */
  subfind_coll_domain_decomposition();

  /* copy over the domain dimensions to serial tree code, as this may be used in the collective unbinding */
  subfind_loctree_copyExtent();

  /* now let us sort according to GrNr and Density. This step will temporarily break the association with SphP[] and other arrays! */
  submp = (struct submp_data *)mymalloc("submp", sizeof(struct submp_data) * NumPart);
  for(i = 0; i < NumPart; i++)
    {
      PS[i].OldIndex      = i;
      submp[i].index      = i;
      submp[i].GrNr       = PS[i].GrNr;
      submp[i].DM_Density = PS[i].Density;
    }
  qsort(submp, NumPart, sizeof(struct submp_data), subfind_compare_submp_GrNr_DM_Density);
  subfind_reorder_according_to_submp();
  myfree(submp);

  /* note: now we have the particles of the group at the beginning, but SPH particles are not aligned.
     They can however be accessed via SphP[PS[i].OldIndex] */

  /* re-determine the number of local group particles, which has changed due to domain decomposition */
  for(i = 0, NumPartGroup = 0; i < NumPart; i++)
    if(PS[i].GrNr == GrNr)
      NumPartGroup++;

  ud  = (struct unbind_data *)mymalloc("ud", NumPartGroup * sizeof(struct unbind_data));
  len = NumPartGroup;

  // pick my particles
  for(i = 0; i < len; i++)
    ud[i].index = i;

  // initialize
  gr_CMFrac = 0;
  gr_Ekin   = 0;
  gr_Ethr   = 0;
  for(k = 0; k < 3; k++)
    {
      gr_Jtot[k]   = 0;
      gr_Jdm[k]    = 0;
      gr_Jgas[k]   = 0;
      gr_Jstars[k] = 0;
    }
  for(k = 0; k < NTYPES; k++)
    {
      gr_CMFracType[k] = 0;
    }

  if(SubThisTask == 0)
    {
      for(k = 0; k < 3; k++)
        {
          gr_pos[k] = Group[grindex].Pos[k];
          gr_vel[k] = Group[grindex].Vel[k];
        }
    }

  // send group properties stored only on root task to all participating tasks
  MPI_Bcast(gr_pos, 3, MPI_DOUBLE, 0, SubComm);
  MPI_Bcast(gr_vel, 3, MPI_DOUBLE, 0, SubComm);

  for(k = 0; k < len; k++)
    {
      index = ud[k].index;
      ptype = P[index].Type;

      for(i = 0; i < 3; i++)
        Pos_pbc[i] = P[index].Pos[i] - gr_pos[i];

      for(i = 0; i < 3; i++)
        Pos_pbc[i] = fof_periodic(Pos_pbc[i]);

      for(i = 0; i < 3; i++)
        Pos_pbc[i] = Pos_pbc[i] * All.cf_atime; /* convert to physical length */

      for(i = 0; i < 3; i++)
        Vel_tot[i] = P[index].Vel[i] / All.cf_atime - gr_vel[i] / All.cf_atime + All.cf_Hrate * Pos_pbc[i];

      gr_Ekin += (P[index].Mass / 2) * (Vel_tot[0] * Vel_tot[0] + Vel_tot[1] * Vel_tot[1] + Vel_tot[2] * Vel_tot[2]);
      if(P[index].Type == 0)
        gr_Ethr += P[index].Mass * SphP[PS[index].OldIndex].Utherm;

      gr_Jtot[0] += P[index].Mass * (Pos_pbc[1] * Vel_tot[2] - Pos_pbc[2] * Vel_tot[1]);
      gr_Jtot[1] += P[index].Mass * (Pos_pbc[2] * Vel_tot[0] - Pos_pbc[0] * Vel_tot[2]);
      gr_Jtot[2] += P[index].Mass * (Pos_pbc[0] * Vel_tot[1] - Pos_pbc[1] * Vel_tot[0]);

      if(ptype == 1)  // dm illustris
        {
          gr_Jdm[0] += P[index].Mass * (Pos_pbc[1] * Vel_tot[2] - Pos_pbc[2] * Vel_tot[1]);
          gr_Jdm[1] += P[index].Mass * (Pos_pbc[2] * Vel_tot[0] - Pos_pbc[0] * Vel_tot[2]);
          gr_Jdm[2] += P[index].Mass * (Pos_pbc[0] * Vel_tot[1] - Pos_pbc[1] * Vel_tot[0]);
        }
      if(ptype == 0)  // gas (incl. winds)
        {
          gr_Jgas[0] += P[index].Mass * (Pos_pbc[1] * Vel_tot[2] - Pos_pbc[2] * Vel_tot[1]);
          gr_Jgas[1] += P[index].Mass * (Pos_pbc[2] * Vel_tot[0] - Pos_pbc[0] * Vel_tot[2]);
          gr_Jgas[2] += P[index].Mass * (Pos_pbc[0] * Vel_tot[1] - Pos_pbc[1] * Vel_tot[0]);
        }
      if(ptype == 4)  // stars
        {
          gr_Jstars[0] += P[index].Mass * (Pos_pbc[1] * Vel_tot[2] - Pos_pbc[2] * Vel_tot[1]);
          gr_Jstars[1] += P[index].Mass * (Pos_pbc[2] * Vel_tot[0] - Pos_pbc[0] * Vel_tot[2]);
          gr_Jstars[2] += P[index].Mass * (Pos_pbc[0] * Vel_tot[1] - Pos_pbc[1] * Vel_tot[0]);
        }
    }

  MPI_Allreduce(MPI_IN_PLACE, gr_Jtot, 3, MPI_DOUBLE, MPI_SUM, SubComm);
  MPI_Allreduce(MPI_IN_PLACE, gr_Jdm, 3, MPI_DOUBLE, MPI_SUM, SubComm);
  MPI_Allreduce(MPI_IN_PLACE, gr_Jgas, 3, MPI_DOUBLE, MPI_SUM, SubComm);
  MPI_Allreduce(MPI_IN_PLACE, gr_Jstars, 3, MPI_DOUBLE, MPI_SUM, SubComm);
  MPI_Allreduce(MPI_IN_PLACE, &gr_Ekin, 1, MPI_DOUBLE, MPI_SUM, SubComm);
  MPI_Allreduce(MPI_IN_PLACE, &gr_Ethr, 1, MPI_DOUBLE, MPI_SUM, SubComm);

  // save the properties
  if(SubThisTask == 0)
    {
      Group[grindex].Ekin = gr_Ekin;
      Group[grindex].Ethr = gr_Ethr;
      for(i = 0; i < 3; i++)
        {
          Group[grindex].J[i]      = gr_Jtot[i];
          Group[grindex].JDM[i]    = gr_Jdm[i];
          Group[grindex].JGas[i]   = gr_Jgas[i];
          Group[grindex].JStars[i] = gr_Jstars[i];
        }
    }

  // calculate counter-rotating fractions
  gr_len_dm = 0;
  gr_mass = gr_mass_gas = gr_mass_stars = 0;

  for(k = 0; k < len; k++)
    {
      index = ud[k].index;
      ptype = P[index].Type;

      for(i = 0; i < 3; i++)
        Pos_pbc[i] = P[index].Pos[i] - gr_pos[i];

      for(i = 0; i < 3; i++)
        Pos_pbc[i] = fof_periodic(Pos_pbc[i]);

      for(i = 0; i < 3; i++)
        Pos_pbc[i] = Pos_pbc[i] * All.cf_atime;  // units: phys kpc/h

      for(i = 0; i < 3; i++)
        Vel_tot[i] = P[index].Vel[i] / All.cf_atime - gr_vel[i] / All.cf_atime + All.cf_Hrate * Pos_pbc[i];

      jpart[0] = P[index].Mass * (Pos_pbc[1] * Vel_tot[2] - Pos_pbc[2] * Vel_tot[1]);
      jpart[1] = P[index].Mass * (Pos_pbc[2] * Vel_tot[0] - Pos_pbc[0] * Vel_tot[2]);
      jpart[2] = P[index].Mass * (Pos_pbc[0] * Vel_tot[1] - Pos_pbc[1] * Vel_tot[0]);

      gr_mass += P[index].Mass;
      if((gr_Jtot[0] * jpart[0] + gr_Jtot[1] * jpart[1] + gr_Jtot[2] * jpart[2]) < 0.)
        gr_CMFrac += P[index].Mass;  // / gr_mass;

      if(ptype == 1)  // dm illustris
        {
          gr_len_dm++;
          if((gr_Jdm[0] * jpart[0] + gr_Jdm[1] * jpart[1] + gr_Jdm[2] * jpart[2]) < 0.)
            gr_CMFracType[1]++;  //= P[index].Mass / gr_mass_dm;
        }
      if(ptype == 0)  // gas (incl. winds)
        {
          gr_mass_gas += P[index].Mass;
          if((gr_Jgas[0] * jpart[0] + gr_Jgas[1] * jpart[1] + gr_Jgas[2] * jpart[2]) < 0.)
            gr_CMFracType[0] += P[index].Mass;  // / gr_mass_gas;
        }
      if(ptype == 4)  // stars
        {
          gr_mass_stars += P[index].Mass;
          if((gr_Jstars[0] * jpart[0] + gr_Jstars[1] * jpart[1] + gr_Jstars[2] * jpart[2]) < 0.)
            gr_CMFracType[4] += P[index].Mass;  // / gr_mass_stars;
        }
    }

  MPI_Allreduce(MPI_IN_PLACE, &gr_mass, 1, MPI_DOUBLE, MPI_SUM, SubComm);
  MPI_Allreduce(MPI_IN_PLACE, &gr_len_dm, 1, MPI_INT, MPI_SUM, SubComm);
  MPI_Allreduce(MPI_IN_PLACE, &gr_mass_gas, 1, MPI_DOUBLE, MPI_SUM, SubComm);
  MPI_Allreduce(MPI_IN_PLACE, &gr_mass_stars, 1, MPI_DOUBLE, MPI_SUM, SubComm);
  MPI_Allreduce(MPI_IN_PLACE, &gr_CMFrac, 1, MPI_DOUBLE, MPI_SUM, SubComm);
  MPI_Allreduce(MPI_IN_PLACE, gr_CMFracType, NTYPES, MPI_DOUBLE, MPI_SUM, SubComm);

  // save the properties
  if(SubThisTask == 0)
    {
      gr_CMFrac /= gr_mass;
      gr_CMFracType[1] /= gr_len_dm;
      gr_CMFracType[0] /= gr_mass_gas;
      gr_CMFracType[4] /= gr_mass_stars;

      Group[grindex].CMFrac = gr_CMFrac;
      for(i = 0; i < NTYPES; i++)
        Group[grindex].CMFracType[i] = gr_CMFracType[i];
    }

  myfree(ud);

  if(SubThisTask == 0)
    printf("SUBFIND-COLLECTIVE: root-task = %d AM done.\n", ThisTask);

  subfind_coll_domain_free();

  /* undo local rearrangement that made group consecutive. After that, the association of SphP[] will be correct again */
  submp = (struct submp_data *)mymalloc("submp", sizeof(struct submp_data) * NumPart);
  for(i = 0; i < NumPart; i++)
    {
      submp[i].index    = i;
      submp[i].OldIndex = PS[i].OldIndex;
    }
  qsort(submp, NumPart, sizeof(struct submp_data), subfind_compare_submp_OldIndex);
  subfind_reorder_according_to_submp();
  myfree(submp);
}
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

/*! \brief Finds candidates for subfind collective.
 *
 *  \param[in] totgrouplen Length of group.
 *
 *  \return void
 */
void subfind_col_find_coll_candidates(int totgrouplen)
{
  int ngbcount, retcode, len_attach;
  int i, k, len, master;
  long long prev, tail, tail_attach, tmp, next, index;
  long long p, ss, head, head_attach, ngb_index1, ngb_index2, rank;
  double t0, t1, tt0, tt1;

  if(SubThisTask == 0)
    {
      printf("SUBFIND-COLLECTIVE, root-task=%d: building distributed linked list. (presently allocated %g MB)\n", ThisTask,
             AllocatedBytes / (1024.0 * 1024.0));
      fflush(stdout);
    }

  /* now find the subhalo coll_candidates by building up link lists from high density to low density */
  t0 = second();
  for(master = 0; master < SubNTask; master++)
    {
      tt0 = second();
      if(SubThisTask != master)
        subfind_poll_for_requests();
      else
        {
          for(k = 0; k < NumPartGroup; k++)
            {
              ngbcount   = sd[k].ngbcount;
              ngb_index1 = sd[k].ngb_index1;
              ngb_index2 = sd[k].ngb_index2;

              switch(ngbcount) /* treat the different possible cases */
                {
                  case 0: /* this appears to be a lonely maximum -> new group */
                    subfind_distlinklist_set_all(sd[k].index, sd[k].index, sd[k].index, 1, -1);
                    break;

                  case 1: /* the particle is attached to exactly one group */
                    head = subfind_distlinklist_get_head(ngb_index1);

                    if(head == -1)
                      terminate("We have a problem!  head=%d/%d for k=%d on task=%d\n", (int)(head >> 32), (int)head, k, SubThisTask);

                    retcode = subfind_distlinklist_get_tail_set_tail_increaselen(head, &tail, sd[k].index);

                    if(!(retcode & 1))
                      subfind_distlinklist_set_headandnext(sd[k].index, head, -1);
                    if(!(retcode & 2))
                      subfind_distlinklist_set_next(tail, sd[k].index);
                    break;

                  case 2: /* the particle merges two groups together */
                    if((ngb_index1 >> 32) == (ngb_index2 >> 32))
                      {
                        subfind_distlinklist_get_two_heads(ngb_index1, ngb_index2, &head, &head_attach);
                      }
                    else
                      {
                        head        = subfind_distlinklist_get_head(ngb_index1);
                        head_attach = subfind_distlinklist_get_head(ngb_index2);
                      }

                    if(head == -1 || head_attach == -1)
                      terminate("We have a problem!  head=%d/%d head_attach=%d/%d for k=%d on task=%d\n", (int)(head >> 32), (int)head,
                                (int)(head_attach >> 32), (int)head_attach, k, SubThisTask);

                    if(head != head_attach)
                      {
                        subfind_distlinklist_get_tailandlen(head, &tail, &len);
                        subfind_distlinklist_get_tailandlen(head_attach, &tail_attach, &len_attach);

                        if(len_attach > len ||
                           (len_attach == len &&
                            head_attach < head)) /* other group is longer, swap them. for equal length, take the larger head value */
                          {
                            tmp         = head;
                            head        = head_attach;
                            head_attach = tmp;
                            tmp         = tail;
                            tail        = tail_attach;
                            tail_attach = tmp;
                            tmp         = len;
                            len         = len_attach;
                            len_attach  = tmp;
                          }

                        /* only in case the attached group is long enough we bother to register it
                           as a subhalo candidate */

                        if(len_attach >= All.DesLinkNgb)
                          {
                            if(count_cand < max_coll_candidates)
                              {
                                coll_candidates[count_cand].len  = len_attach;
                                coll_candidates[count_cand].head = head_attach;
                                count_cand++;
                              }
                            else
                              terminate("Task %d: count=%d, max=%d, npartgroup=%d\n", SubThisTask, count_cand, max_coll_candidates,
                                        NumPartGroup);
                          }

                        /* now join the two groups */
                        subfind_distlinklist_set_tailandlen(head, tail_attach, len + len_attach);
                        subfind_distlinklist_set_next(tail, head_attach);

                        ss = head_attach;
                        do
                          {
                            ss = subfind_distlinklist_set_head_get_next(ss, head);
                          }
                        while(ss >= 0);
                      }

                    /* finally, attach the particle to 'head' */
                    retcode = subfind_distlinklist_get_tail_set_tail_increaselen(head, &tail, sd[k].index);

                    if(!(retcode & 1))
                      subfind_distlinklist_set_headandnext(sd[k].index, head, -1);
                    if(!(retcode & 2))
                      subfind_distlinklist_set_next(tail, sd[k].index);
                    break;
                }
            }

          fflush(stdout);

          /* now tell the others to stop polling */
          for(k = 0; k < SubNTask; k++)
            if(k != SubThisTask)
              MPI_Send(&k, 1, MPI_INT, k, TAG_POLLING_DONE, SubComm);
        }

      MPI_Barrier(SubComm);
      tt1 = second();
      if(SubThisTask == 0)
        {
          printf("SUBFIND-COLLECTIVE, root-task=%d: ma=%d/%d took %g sec\n", ThisTask, master, SubNTask, timediff(tt0, tt1));
          fflush(stdout);
        }
    }
  t1 = second();
  if(SubThisTask == 0)
    printf("SUBFIND-COLLECTIVE, root-task=%d: identification of primary coll_candidates took %g sec\n", ThisTask, timediff(t0, t1));

  /* add the full thing as a subhalo candidate */
  t0 = second();
  for(master = 0, head = -1, prev = -1; master < SubNTask; master++)
    {
      if(SubThisTask != master)
        subfind_poll_for_requests();
      else
        {
          for(i = 0; i < NumPartGroup; i++)
            {
              index = (((long long)SubThisTask) << 32) + i;

              if(Head[i] == index)
                {
                  subfind_distlinklist_get_tailandlen(Head[i], &tail, &len);
                  next = subfind_distlinklist_get_next(tail);
                  if(next == -1)
                    {
                      if(prev < 0)
                        head = index;

                      if(prev >= 0)
                        subfind_distlinklist_set_next(prev, index);

                      prev = tail;
                    }
                }
            }

          /* now tell the others to stop polling */
          for(k = 0; k < SubNTask; k++)
            if(k != SubThisTask)
              MPI_Send(&k, 1, MPI_INT, k, TAG_POLLING_DONE, SubComm);
        }

      MPI_Barrier(SubComm);
      MPI_Bcast(&head, sizeof(head), MPI_BYTE, master, SubComm);
      MPI_Bcast(&prev, sizeof(prev), MPI_BYTE, master, SubComm);
    }

  if(SubThisTask == SubNTask - 1)
    {
      if(count_cand < max_coll_candidates)
        {
          coll_candidates[count_cand].len  = totgrouplen;
          coll_candidates[count_cand].head = head;
          count_cand++;
        }
      else
        terminate("count_cand=%d >= max_coll_candidates=%d", count_cand, max_coll_candidates);
    }
  t1 = second();
  if(SubThisTask == 0)
    printf("SUBFIND-COLLECTIVE, root-task=%d: adding background as candidate took %g sec\n", ThisTask, timediff(t0, t1));

  /* go through the whole chain once to establish a rank order. For the rank we use Len[] */
  t0 = second();

  master = (head >> 32);

  if(SubThisTask != master)
    subfind_poll_for_requests();
  else
    {
      p    = head;
      rank = 0;

      while(p >= 0)
        {
          p = subfind_distlinklist_setrank_and_get_next(p, &rank);
        }

      /* now tell the others to stop polling */
      for(i = 0; i < SubNTask; i++)
        if(i != master)
          MPI_Send(&i, 1, MPI_INT, i, TAG_POLLING_DONE, SubComm);
    }

  MPI_Barrier(SubComm);
  MPI_Bcast(&rank, sizeof(rank), MPI_BYTE, master, SubComm); /* just for testing */

  /* for each candidate, we now pull out the rank of its head */
  for(master = 0; master < SubNTask; master++)
    {
      if(SubThisTask != master)
        subfind_poll_for_requests();
      else
        {
          for(k = 0; k < count_cand; k++)
            coll_candidates[k].rank = subfind_distlinklist_get_rank(coll_candidates[k].head);

          /* now tell the others to stop polling */
          for(i = 0; i < SubNTask; i++)
            if(i != SubThisTask)
              MPI_Send(&i, 1, MPI_INT, i, TAG_POLLING_DONE, SubComm);
        }
    }
  MPI_Barrier(SubComm);

  t1 = second();
  if(SubThisTask == 0)
    printf("SUBFIND-COLLECTIVE, root-task=%d: establishing of rank order took %g sec  (p=%d, grouplen=%d) presently allocated %g MB\n",
           ThisTask, timediff(t0, t1), (int)rank, totgrouplen, AllocatedBytes / (1024.0 * 1024.0));

  if(((int)rank) != totgrouplen)
    terminate("mismatch\n");
}

/*! \brief Unbinding for independent subgroups.
 *
 *  \param[in] cont_cand Number of subgroup candidates.
 *
 *  \return void
 */
void subfind_unbind_independent_ones(int count_cand)
{
  int i, j, k, len, nsubs, len_non_gas;

  ud = (struct unbind_data *)mymalloc("ud", NumPart * sizeof(struct unbind_data));

  subfind_loctree_treeallocate(All.TreeAllocFactor * NumPart, NumPart);

  qsort(coll_candidates, count_cand, sizeof(struct coll_cand_dat), subfind_compare_coll_candidates_nsubs);

  for(k = 0, i = 0; k < count_cand; k++)
    if(coll_candidates[k].parent == 0)
      {
        while(PS[i].submark < coll_candidates[k].nsub)
          {
            i++;
            if(i >= NumPart)
              terminate("i >= NumPart");
          }

        if(PS[i].submark >= 0 && PS[i].submark < HIGHBIT)
          {
            len   = 0;
            nsubs = PS[i].submark;

            if(nsubs != coll_candidates[k].nsub)
              {
                terminate("TASK=%d i=%d k=%d nsubs=%d coll_candidates[k].nsub=%d\n", SubThisTask, i, k, nsubs,
                          coll_candidates[k].nsub);
              }

            while(i < NumPart)
              {
                if(PS[i].submark == nsubs)
                  {
                    PS[i].submark = HIGHBIT;
                    if((PS[i].origintask & HIGHBIT) == 0)
                      {
                        ud[len].index = i;
                        len++;
                      }
                    i++;
                  }
                else
                  break;
              }

            /* call the serial unbind function */
            len = subfind_unbind(ud, len, &len_non_gas);

            if(len >= All.DesLinkNgb)
              {
                /* ok, we found a substructure */
                coll_candidates[k].bound_length = len;

                for(j = 0; j < len; j++)
                  PS[ud[j].index].submark = nsubs; /* we use this to flag the substructures */
              }
            else
              coll_candidates[k].bound_length = 0;
          }
      }

  subfind_loctree_treefree();

  myfree(ud);
}

/*! \brief Unbinding for subfind collective.
 *
 *  \param[in] d Unbind data.
 *  \param[in] num Number of particles in subgroup.
 *  \param[out] num_non_gas Number of particles which are not gas cells.
 *
 *  \return
 */
int subfind_col_unbind(struct unbind_data *d, int num, int *num_non_gas)
{
  int iter = 0;
  int i, j, p, part_index, minindex, task;
  int unbound, totunbound, numleft, mincpu;
  int *npart, *offset, *nbu_count, count_bound_unbound, phaseflag;
  double s[3], dx[3], ddxx, v[3], dv[3], sloc[3], vloc[3], pos[3];
  double vel_to_phys, atime;
  MyFloat minpot, *potlist;
  double boxsize, xtmp;
  double mass, massloc;
  double *bnd_energy, energy_limit, energy_limit_local, weakly_bound_limit_local, weakly_bound_limit = 0;

  if(SubThisTask == 0)
    {
      printf("SUBFIND-COLLECTIVE, root-task=%d: beginning of subfind_col_unbind()\n", ThisTask);
      fflush(stdout);
    }

  boxsize = All.BoxSize;

  vel_to_phys = 1.0 / All.cf_atime;
  atime       = All.cf_atime;

  phaseflag = 0; /* this means we will recompute the potential for all particles */

  do
    {
      subfind_coll_treebuild(num, d);

      /* let's compute the potential energy */

      subfind_potential_compute(num, d, phaseflag, weakly_bound_limit);

      if(phaseflag == 0)
        {
          potlist = (MyFloat *)mymalloc("potlist", SubNTask * sizeof(MyFloat));

          for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
            {
              if(gsl_isnan(PS[d[i].index].Potential))
                terminate("pot is nan");

              if(PS[d[i].index].Potential < minpot || minindex == -1)
                {
                  minpot   = PS[d[i].index].Potential;
                  minindex = d[i].index;
                }
            }

          MPI_Allgather(&minpot, sizeof(MyFloat), MPI_BYTE, potlist, sizeof(MyFloat), MPI_BYTE, SubComm);

          for(i = 0, mincpu = -1, minpot = 1.0e30; i < SubNTask; i++)
            if(potlist[i] < minpot)
              {
                mincpu = i;
                minpot = potlist[i];
              }

          if(mincpu < 0)
            terminate("mincpu < 0");

          myfree(potlist);

          if(SubThisTask == mincpu)
            {
#ifdef CELL_CENTER_GRAVITY
              if(P[minindex].Type == 0)
                {
                  for(j = 0; j < 3; j++)
                    pos[j] = PS[minindex].Center[j];
                }
              else
#endif /* #ifdef CELL_CENTER_GRAVITY */
                {
                  for(j = 0; j < 3; j++)
                    pos[j] = P[minindex].Pos[j];
                }
            }

          MPI_Bcast(&pos[0], 3, MPI_DOUBLE, mincpu, SubComm);
          /* pos[] now holds the position of minimum potential */
          /* we take that as the center */
        }

      /* let's get bulk velocity and the center-of-mass */

      for(j = 0; j < 3; j++)
        sloc[j] = vloc[j] = 0;

      for(i = 0, massloc = 0; i < num; i++)
        {
          part_index = d[i].index;

          for(j = 0; j < 3; j++)
            {
#ifdef CELL_CENTER_GRAVITY
              if(P[part_index].Type == 0)
                ddxx = GRAVITY_NEAREST_X(PS[part_index].Center[j] - pos[j]);
              else
#endif /* #ifdef CELL_CENTER_GRAVITY */
                ddxx = GRAVITY_NEAREST_X(P[part_index].Pos[j] - pos[j]);

              sloc[j] += P[part_index].Mass * ddxx;
              vloc[j] += P[part_index].Mass * P[part_index].Vel[j];
            }
          massloc += P[part_index].Mass;
        }

      MPI_Allreduce(sloc, s, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(vloc, v, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&massloc, &mass, 1, MPI_DOUBLE, MPI_SUM, SubComm);

      for(j = 0; j < 3; j++)
        {
          s[j] /= mass; /* center of mass */
          v[j] /= mass;

          s[j] += pos[j];

          while(s[j] < 0)
            s[j] += boxsize;
          while(s[j] >= boxsize)
            s[j] -= boxsize;
        }

      bnd_energy = (double *)mymalloc("bnd_energy", num * sizeof(double));

      for(i = 0; i < num; i++)
        {
          part_index = d[i].index;

          for(j = 0; j < 3; j++)
            {
              dv[j] = vel_to_phys * (P[part_index].Vel[j] - v[j]);

#ifdef CELL_CENTER_GRAVITY
              if(P[part_index].Type == 0)
                dx[j] = atime * GRAVITY_NEAREST_X(PS[part_index].Center[j] - s[j]);
              else
#endif /* #ifdef CELL_CENTER_GRAVITY */
                dx[j] = atime * GRAVITY_NEAREST_X(P[part_index].Pos[j] - s[j]);

              dv[j] += All.cf_Hrate * dx[j];
            }

          PS[part_index].BindingEnergy = PS[part_index].Potential + 0.5 * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
          PS[part_index].BindingEnergy += All.G / All.cf_atime * P[part_index].Mass /
                                          (All.ForceSoftening[P[part_index].SofteningType] / 2.8); /* add self-energy */

          if(P[part_index].Type == 0)
            PS[part_index].BindingEnergy += PS[part_index].Utherm;

          bnd_energy[i] = PS[part_index].BindingEnergy;
        }

      parallel_sort_comm(bnd_energy, num, sizeof(double), subfind_compare_binding_energy, SubComm);

      npart     = (int *)mymalloc("npart", SubNTask * sizeof(int));
      nbu_count = (int *)mymalloc("nbu_count", SubNTask * sizeof(int));
      offset    = (int *)mymalloc("offset", SubNTask * sizeof(int));

      MPI_Allgather(&num, 1, MPI_INT, npart, 1, MPI_INT, SubComm);
      MPI_Allreduce(&num, &numleft, 1, MPI_INT, MPI_SUM, SubComm);

      for(i = 1, offset[0] = 0; i < SubNTask; i++)
        offset[i] = offset[i - 1] + npart[i - 1];

      j = (int)(0.25 * numleft); /* index of limiting energy value */

      task = 0;
      while(j >= npart[task])
        {
          j -= npart[task];
          task++;
        }

      if(SubThisTask == task)
        energy_limit_local = bnd_energy[j];
      else
        energy_limit_local = 1.0e30;

      MPI_Allreduce(&energy_limit_local, &energy_limit, 1, MPI_DOUBLE, MPI_MIN, SubComm);

      for(i = 0, count_bound_unbound = 0; i < num; i++)
        {
          if(bnd_energy[i] > 0)
            count_bound_unbound++;
          else
            count_bound_unbound--;
        }

      MPI_Allgather(&count_bound_unbound, 1, MPI_INT, nbu_count, 1, MPI_INT, SubComm);

      for(i = 0, count_bound_unbound = 0; i < SubThisTask; i++)
        count_bound_unbound += nbu_count[i];

      for(i = 0; i < num - 1; i++)
        {
          if(bnd_energy[i] > 0)
            count_bound_unbound++;
          else
            count_bound_unbound--;
          if(count_bound_unbound <= 0)
            break;
        }

      if(num > 0 && count_bound_unbound <= 0)
        weakly_bound_limit_local = bnd_energy[i];
      else
        weakly_bound_limit_local = -1.0e30;

      MPI_Allreduce(&weakly_bound_limit_local, &weakly_bound_limit, 1, MPI_DOUBLE, MPI_MAX, SubComm);

      for(i = 0, unbound = 0; i < num; i++)
        {
          p = d[i].index;

          if(PS[p].BindingEnergy > 0 && PS[p].BindingEnergy > energy_limit)
            {
              unbound++;

              d[i] = d[num - 1];
              num--;
              i--;
            }
          else if(P[p].Type != 0)
            (*num_non_gas)++;
        }

      myfree(offset);
      myfree(nbu_count);
      myfree(npart);
      myfree(bnd_energy);

      MPI_Allreduce(&unbound, &totunbound, 1, MPI_INT, MPI_SUM, SubComm);
      MPI_Allreduce(&num, &numleft, 1, MPI_INT, MPI_SUM, SubComm);

      if(phaseflag == 0)
        {
          if(totunbound > 0)
            phaseflag = 1;
        }
      else
        {
          if(totunbound == 0)
            {
              phaseflag  = 0; /* this will make us repeat everything once more for all particles */
              totunbound = 1;
            }
        }

      iter++;
    }
  while(totunbound > 0 && numleft >= All.DesLinkNgb);

  return num;
}

/*! \brief Gets new request from other task.
 *
 *  \return void
 */
void subfind_poll_for_requests(void)
{
  int index, nsub, source, tag, ibuf[3], target, submark, task;
  long long head, next, rank, buf[5];
  long long oldtail, newtail;
  int task_newtail, i_newtail, task_oldtail, i_oldtail;
  char msg[200];
  MPI_Status status;

  do
    {
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, SubComm, &status);

      source = status.MPI_SOURCE;
      tag    = status.MPI_TAG;

      /* MPI_Get_count(&status, MPI_BYTE, &count); */
      switch(tag)
        {
          case TAG_GET_TWOHEADS:
            MPI_Recv(ibuf, 2, MPI_INT, source, TAG_GET_TWOHEADS, SubComm, MPI_STATUS_IGNORE);
            buf[0] = Head[ibuf[0]];
            buf[1] = Head[ibuf[1]];
            MPI_Send(buf, 2 * sizeof(long long), MPI_BYTE, source, TAG_GET_TWOHEADS_DATA, SubComm);
            break;
          case TAG_SET_NEWTAIL:
            MPI_Recv(buf, 2 * sizeof(long long), MPI_BYTE, source, TAG_SET_NEWTAIL, SubComm, MPI_STATUS_IGNORE);
            index       = buf[0];
            newtail     = buf[1];
            oldtail     = Tail[index]; /* return old tail */
            Tail[index] = newtail;
            Len[index]++;

            task_newtail = (newtail >> 32);
            if(task_newtail == SubThisTask)
              {
                i_newtail       = (newtail & MASK);
                Head[i_newtail] = (((long long)SubThisTask) << 32) + index;
                Next[i_newtail] = -1;
              }
            task_oldtail = (oldtail >> 32);
            if(task_oldtail == SubThisTask)
              {
                i_oldtail       = (oldtail & MASK);
                Next[i_oldtail] = newtail;
              }

            buf[0] = oldtail;
            MPI_Send(buf, 1 * sizeof(long long), MPI_BYTE, source, TAG_GET_OLDTAIL, SubComm);
            break;
          case TAG_SET_ALL:
            MPI_Recv(buf, 5 * sizeof(long long), MPI_BYTE, source, TAG_SET_ALL, SubComm, MPI_STATUS_IGNORE);
            index       = buf[0];
            Head[index] = buf[1];
            Tail[index] = buf[2];
            Len[index]  = buf[3];
            Next[index] = buf[4];
            break;
          case TAG_GET_TAILANDLEN:
            MPI_Recv(&index, 1, MPI_INT, source, tag, SubComm, &status);
            buf[0] = Tail[index];
            buf[1] = Len[index];
            MPI_Send(buf, 2 * sizeof(long long), MPI_BYTE, source, TAG_GET_TAILANDLEN_DATA, SubComm);
            break;
          case TAG_SET_TAILANDLEN:
            MPI_Recv(buf, 3 * sizeof(long long), MPI_BYTE, source, TAG_SET_TAILANDLEN, SubComm, MPI_STATUS_IGNORE);
            index       = buf[0];
            Tail[index] = buf[1];
            Len[index]  = buf[2];
            break;
          case TAG_SET_HEADANDNEXT:
            MPI_Recv(buf, 3 * sizeof(long long), MPI_BYTE, source, TAG_SET_HEADANDNEXT, SubComm, MPI_STATUS_IGNORE);
            index       = buf[0];
            Head[index] = buf[1];
            Next[index] = buf[2];
            break;
          case TAG_SET_NEXT:
            MPI_Recv(buf, 2 * sizeof(long long), MPI_BYTE, source, TAG_SET_NEXT, SubComm, MPI_STATUS_IGNORE);
            index       = buf[0];
            Next[index] = buf[1];
            break;
          case TAG_SETHEADGETNEXT:
            MPI_Recv(buf, 2 * sizeof(long long), MPI_BYTE, source, TAG_SETHEADGETNEXT, SubComm, MPI_STATUS_IGNORE);
            index = buf[0];
            head  = buf[1];
            do
              {
                Head[index] = head;
                next        = Next[index];
                task        = (next >> 32);
                index       = (next & MASK);
              }
            while(next >= 0 && task == SubThisTask);
            MPI_Send(&next, 1 * sizeof(long long), MPI_BYTE, source, TAG_SETHEADGETNEXT_DATA, SubComm);
            break;
          case TAG_GET_NEXT:
            MPI_Recv(&index, 1, MPI_INT, source, tag, SubComm, &status);
            MPI_Send(&Next[index], 1 * sizeof(long long), MPI_BYTE, source, TAG_GET_NEXT_DATA, SubComm);
            break;
          case TAG_GET_HEAD:
            MPI_Recv(&index, 1, MPI_INT, source, tag, SubComm, &status);
            MPI_Send(&Head[index], 1 * sizeof(long long), MPI_BYTE, source, TAG_GET_HEAD_DATA, SubComm);
            break;
          case TAG_ADD_PARTICLE:
            MPI_Recv(&index, 1, MPI_INT, source, tag, SubComm, &status);
            if(Tail[index] < 0) /* consider only particles not already in substructures */
              {
                ud[LocalLen].index = index;
                if(index >= NumPartGroup)
                  {
                    sprintf(msg, "What: index=%d NumPartGroup=%d\n", index, NumPartGroup);
                    terminate(msg);
                  }
                LocalLen++;
              }
            break;
          case TAG_MARK_PARTICLE:
            MPI_Recv(ibuf, 3, MPI_INT, source, TAG_MARK_PARTICLE, SubComm, MPI_STATUS_IGNORE);
            index   = ibuf[0];
            target  = ibuf[1];
            submark = ibuf[2];

            if(PS[index].submark != HIGHBIT)
              terminate("TasK=%d i=%d P[i].submark=%d?\n", SubThisTask, index, PS[index].submark);

            PS[index].TargetTask = target;
            PS[index].submark    = submark;
            break;
          case TAG_ADDBOUND:
            MPI_Recv(ibuf, 2, MPI_INT, source, TAG_ADDBOUND, SubComm, &status);
            index = ibuf[0];
            nsub  = ibuf[1];
            if(Tail[index] == nsub) /* consider only particles in this substructure */
              {
                ud[LocalLen].index = index;
                LocalLen++;
              }
            break;
          case TAG_SETRANK:
            MPI_Recv(buf, 2 * sizeof(long long), MPI_BYTE, source, TAG_SETRANK, SubComm, MPI_STATUS_IGNORE);
            index = buf[0];
            rank  = buf[1];
            do
              {
                Len[index] = rank++;
                next       = Next[index];
                if(next < 0)
                  break;
                index = (next & MASK);
              }
            while((next >> 32) == SubThisTask);
            buf[0] = next;
            buf[1] = rank;
            MPI_Send(buf, 2 * sizeof(long long), MPI_BYTE, source, TAG_SETRANK_OUT, SubComm);
            break;
          case TAG_GET_RANK:
            MPI_Recv(&index, 1, MPI_INT, source, tag, SubComm, &status);
            rank = Len[index];
            MPI_Send(&rank, 1 * sizeof(long long), MPI_BYTE, source, TAG_GET_RANK_DATA, SubComm);
            break;

          case TAG_POLLING_DONE:
            MPI_Recv(&index, 1, MPI_INT, source, tag, SubComm, &status);
            break;

          default:
            terminate("tag not present in the switch");
            break;
        }
    }
  while(tag != TAG_POLLING_DONE);
}

/*! \brief Sets rank in global linked list and gets next entry.
 *
 *  \param[in] index Index in global linked list.
 *  \param[in, out] rank Rank to be set in linked list.
 *
 *  \return Next entry
 */
long long subfind_distlinklist_setrank_and_get_next(long long index, long long *rank)
{
  int task, i;
  long long next;
  long long buf[2];

  task = (index >> 32);
  i    = (index & MASK);

  if(SubThisTask == task)
    {
      Len[i] = *rank;
      *rank  = *rank + 1;
      next   = Next[i];
    }
  else
    {
      buf[0] = i;
      buf[1] = *rank;

      MPI_Send(buf, 2 * sizeof(long long), MPI_BYTE, task, TAG_SETRANK, SubComm);
      MPI_Recv(buf, 2 * sizeof(long long), MPI_BYTE, task, TAG_SETRANK_OUT, SubComm, MPI_STATUS_IGNORE);
      next  = buf[0];
      *rank = buf[1];
    }
  return next;
}

/*! \brief Sets head in global linked list and gets next
 *
 *  \param[in] index Index in global linked list.
 *  \param[in] head Head value to be set.
 *
 *  \return Next value.
 */
long long subfind_distlinklist_set_head_get_next(long long index, long long head)
{
  int task, i;
  long long buf[2];
  long long next;

  task = (index >> 32);
  i    = (index & MASK);

  if(SubThisTask == task)
    {
      Head[i] = head;
      next    = Next[i];
    }
  else
    {
      buf[0] = i;
      buf[1] = head;
      MPI_Send(buf, 2 * sizeof(long long), MPI_BYTE, task, TAG_SETHEADGETNEXT, SubComm);
      MPI_Recv(&next, 1 * sizeof(long long), MPI_BYTE, task, TAG_SETHEADGETNEXT_DATA, SubComm, MPI_STATUS_IGNORE);
    }

  return next;
}

/*! \brief Sets next value in global linked list.
 *
 *  \param[in] index Index in global linked list.
 *  \param[in] next Next value to be set.
 *
 *  \return void
 */
void subfind_distlinklist_set_next(long long index, long long next)
{
  int task, i;
  long long buf[2];

  task = (index >> 32);
  i    = (index & MASK);

  if(SubThisTask == task)
    {
      Next[i] = next;
    }
  else
    {
      buf[0] = i;
      buf[1] = next;
      MPI_Send(buf, 2 * sizeof(long long), MPI_BYTE, task, TAG_SET_NEXT, SubComm);
    }
}

/*! \brief Adds particle to 'ud' list if not already in substructure.
 *
 *  \param[in] index Index in global linked list.
 *
 *  \return void
 */
void subfind_distlinklist_add_particle(long long index)
{
  int task, i;
  char msg[200];

  task = (index >> 32);
  i    = (index & MASK);

  if(SubThisTask == task)
    {
      if(Tail[i] < 0) /* consider only particles not already in substructures */
        {
          ud[LocalLen].index = i;
          if(i >= NumPartGroup)
            {
              sprintf(msg, "What: index=%d NumPartGroup=%d\n", i, NumPartGroup);
              terminate(msg);
            }

          LocalLen++;
        }
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_ADD_PARTICLE, SubComm);
    }
}

/*! \brief Sets target task and submark field in 'PS' structure.
 *
 *  \param[in] index Index in global linked list
 *  \param[in] target Value for TargetTask field.
 *  \param[in] submark Value for submark field.
 *
 *  \return void
 */
void subfind_distlinklist_mark_particle(long long index, int target, int submark)
{
  int task, i, ibuf[3];

  task = (index >> 32);
  i    = (index & MASK);

  if(SubThisTask == task)
    {
      if(PS[i].submark != HIGHBIT)
        terminate("Tas=%d i=%d P[i].submark=%d?\n", SubThisTask, i, PS[i].submark);

      PS[i].TargetTask = target;
      PS[i].submark    = submark;
    }
  else
    {
      ibuf[0] = i;
      ibuf[1] = target;
      ibuf[2] = submark;
      MPI_Send(ibuf, 3, MPI_INT, task, TAG_MARK_PARTICLE, SubComm);
    }
}

/*! \brief Add bound particle to 'ud' array.
 *
 *  \param[in] index Index in global linked list.
 *  \param[in] nsub Number of subgroups (i.e. if Tail index the same, not yet
 *             in a substructrue).
 *
 *  \return void
 */
void subfind_distlinklist_add_bound_particles(long long index, int nsub)
{
  int task, i, ibuf[2];

  task = (index >> 32);
  i    = (index & MASK);

  if(SubThisTask == task)
    {
      if(Tail[i] == nsub) /* consider only particles not already in substructures */
        {
          ud[LocalLen].index = i;
          LocalLen++;
        }
    }
  else
    {
      ibuf[0] = i;
      ibuf[1] = nsub;
      MPI_Send(ibuf, 2, MPI_INT, task, TAG_ADDBOUND, SubComm);
    }
}

/*! \brief Get Next value from global linked list.
 *
 *  \param[in] index Index in global linked list.
 *
 *  \return
 */
long long subfind_distlinklist_get_next(long long index)
{
  int task, i;
  long long next;

  task = (index >> 32);
  i    = (index & MASK);

  if(SubThisTask == task)
    {
      next = Next[i];
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_GET_NEXT, SubComm);
      MPI_Recv(&next, 1 * sizeof(long long), MPI_BYTE, task, TAG_GET_NEXT_DATA, SubComm, MPI_STATUS_IGNORE);
    }

  return next;
}

/*! \brief Get rank value from global linked list.
 *
 *  \param[in] index Index in global linked list.
 *
 *  \return Rank value.
 */
long long subfind_distlinklist_get_rank(long long index)
{
  int task, i;
  long long rank;

  task = (index >> 32);
  i    = (index & MASK);

  if(SubThisTask == task)
    {
      rank = Len[i];
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_GET_RANK, SubComm);
      MPI_Recv(&rank, 1 * sizeof(long long), MPI_BYTE, task, TAG_GET_RANK_DATA, SubComm, MPI_STATUS_IGNORE);
    }

  return rank;
}

/*! \brief Get the head value of global linked list.
 *
 *  \param[in] index Index in the global linked list.
 *
 *  \return Head value.
 */
long long subfind_distlinklist_get_head(long long index)
{
  int task, i;
  long long head;

  task = (index >> 32);
  i    = (index & MASK);

  if(SubThisTask == task)
    {
      head = Head[i];
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_GET_HEAD, SubComm);
      MPI_Recv(&head, 1 * sizeof(long long), MPI_BYTE, task, TAG_GET_HEAD_DATA, SubComm, MPI_STATUS_IGNORE);
    }

  return head;
}

/*! \brief Gets the head value of two entries in linked list.
 *
 *  \param[in] ngb_index1 Index of first subgroup.
 *  \param[in] ngb_index2 Index of second subgroup.
 *  \param[out] head Head value of first subgroup.
 *  \param[out] head_attach head value of second subgroup.
 *
 *  \return void
 */
void subfind_distlinklist_get_two_heads(long long ngb_index1, long long ngb_index2, long long *head, long long *head_attach)
{
  int task, i1, i2, ibuf[2];
  long long buf[2];

  task = (ngb_index1 >> 32);
  i1   = (ngb_index1 & MASK);
  i2   = (ngb_index2 & MASK);

  if(SubThisTask == task)
    {
      *head        = Head[i1];
      *head_attach = Head[i2];
    }
  else
    {
      ibuf[0] = i1;
      ibuf[1] = i2;
      MPI_Send(ibuf, 2, MPI_INT, task, TAG_GET_TWOHEADS, SubComm);
      MPI_Recv(buf, 2 * sizeof(long long), MPI_BYTE, task, TAG_GET_TWOHEADS_DATA, SubComm, MPI_STATUS_IGNORE);
      *head        = buf[0];
      *head_attach = buf[1];
    }
}

/*! \brief Sets Head and Next entries in global linked list.
 *
 *  \param[in] index Index in global linked list.
 *  \param[in] head Value for Head.
 *  \param[in] next Value for Next.
 *
 *  \return void
 */
void subfind_distlinklist_set_headandnext(long long index, long long head, long long next)
{
  int task, i;
  long long buf[3];

  task = (index >> 32);
  i    = (index & MASK);

  if(SubThisTask == task)
    {
      Head[i] = head;
      Next[i] = next;
    }
  else
    {
      buf[0] = i;
      buf[1] = head;
      buf[2] = next;
      MPI_Send(buf, 3 * sizeof(long long), MPI_BYTE, task, TAG_SET_HEADANDNEXT, SubComm);
    }
}

/*! \brief Returns old tail, sets a new tail, increases length of linked list.
 *
 *  \param[in] index Index of the subgroup.
 *  \param[out] tail Old value for tail.
 *  \param[in] newtail New value for tail.
 *
 *  \return return code
 */
int subfind_distlinklist_get_tail_set_tail_increaselen(long long index, long long *tail, long long newtail)
{
  int task, i, task_newtail, i_newtail, task_oldtail, i_oldtail, retcode;
  long long oldtail;
  long long buf[2];

  task = (index >> 32);
  i    = (index & MASK);

  retcode = 0;

  if(SubThisTask == task)
    {
      oldtail = Tail[i];
      Tail[i] = newtail;
      Len[i]++;
      *tail = oldtail;

      task_newtail = (newtail >> 32);
      if(task_newtail == SubThisTask)
        {
          i_newtail       = (newtail & MASK);
          Head[i_newtail] = index;
          Next[i_newtail] = -1;
          retcode |= 1;
        }
      task_oldtail = (oldtail >> 32);
      if(task_oldtail == SubThisTask)
        {
          i_oldtail       = (oldtail & MASK);
          Next[i_oldtail] = newtail;
          retcode |= 2;
        }
    }
  else
    {
      buf[0] = i;
      buf[1] = newtail;
      MPI_Send(buf, 2 * sizeof(long long), MPI_BYTE, task, TAG_SET_NEWTAIL, SubComm);
      MPI_Recv(&oldtail, 1 * sizeof(long long), MPI_BYTE, task, TAG_GET_OLDTAIL, SubComm, MPI_STATUS_IGNORE);
      *tail = oldtail;

      if((newtail >> 32) == task)
        retcode |= 1;
      if((oldtail >> 32) == task)
        retcode |= 2;
    }

  return retcode;
}

/*! \brief Set tail and len in global linked list.
 *
 *  \param[in] index Index in global linked list.
 *  \param[in] tail Value to be set in 'Tail'.
 *  \param[in] len Value to be set in 'Len'.
 *
 *  \return void
 */
void subfind_distlinklist_set_tailandlen(long long index, long long tail, int len)
{
  int task, i;
  long long buf[3];

  task = (index >> 32);
  i    = (index & MASK);

  if(SubThisTask == task)
    {
      Tail[i] = tail;
      Len[i]  = len;
    }
  else
    {
      buf[0] = i;
      buf[1] = tail;
      buf[2] = len;
      MPI_Send(buf, 3 * sizeof(long long), MPI_BYTE, task, TAG_SET_TAILANDLEN, SubComm);
    }
}

/*! \brief Get tail and len in global linked list.
 *
 *  \param[in] index Index in global linked list.
 *  \param[out] tail 'Tail' value.
 *  \param[out] len 'Len' value.
 *
 *  \return void
 */
void subfind_distlinklist_get_tailandlen(long long index, long long *tail, int *len)
{
  int task, i;
  long long buf[2];

  task = (index >> 32);
  i    = (index & MASK);

  if(SubThisTask == task)
    {
      *tail = Tail[i];
      *len  = Len[i];
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_GET_TAILANDLEN, SubComm);
      MPI_Recv(buf, 2 * sizeof(long long), MPI_BYTE, task, TAG_GET_TAILANDLEN_DATA, SubComm, MPI_STATUS_IGNORE);
      *tail = buf[0];
      *len  = buf[1];
    }
}

/*! \brief Sets head, tail, len and next in global linked list
 *
 *  \param[in] index Index in global linked list.
 *  \param[in] head Value for 'Head'.
 *  \param[in] tail Value for 'Tail'.
 *  \param[in] len Value for 'Len'.
 *  \param[in] next Value for 'Next'.
 *
 *  \return void
 */
void subfind_distlinklist_set_all(long long index, long long head, long long tail, int len, long long next)
{
  int task, i;
  long long buf[5];

  task = (index >> 32);
  i    = (index & MASK);

  if(SubThisTask == task)
    {
      Head[i] = head;
      Tail[i] = tail;
      Len[i]  = len;
      Next[i] = next;
    }
  else
    {
      buf[0] = i;
      buf[1] = head;
      buf[2] = tail;
      buf[3] = len;
      buf[4] = next;
      MPI_Send(buf, 5 * sizeof(long long), MPI_BYTE, task, TAG_SET_ALL, SubComm);
    }
}

/*! \brief Comparison function of sort_density_data objects.
 *
 *  Compares element density.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1); -1 if a > b
 */
int subfind_compare_densities(const void *a, const void *b) /* largest density first */
{
  if(((struct sort_density_data *)a)->density > (((struct sort_density_data *)b)->density))
    return -1;

  if(((struct sort_density_data *)a)->density < (((struct sort_density_data *)b)->density))
    return +1;

  return 0;
}

#endif
