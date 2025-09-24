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
 * \file        src/subfind/subfind_coll_treewalk.c
 * \date        05/2018
 * \brief       Algorithm for collective tree walk; computes gravitational
 *              binding energy.
 * \details     contains functions:
 *                static void particle2in(data_in * in, int i, int firstnode)
 *                static void out2particle(data_out * out, int i, int mode)
 *                static void kernel_local(void)
 *                static void kernel_imported(void)
 *                void subfind_potential_compute(int num, struct unbind_data
 *                  *darg, int phasearg, double weakly_bound_limit_arg)
 *                static int subfind_force_treeevaluate_potential(int target,
 *                  int mode, int threadid)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 15.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef SUBFIND
#include "../fof/fof.h"
#include "subfind.h"

static int subfind_force_treeevaluate_potential(int target, int mode, int threadid);

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  unsigned char SofteningType;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/*! \brief Routine that fills the relevant particle/cell data into the input
 *         structure defined above. Needed by generic_comm_helpers2.
 *
 *  \param[out] in Data structure to fill.
 *  \param[in] i Index of particle in P and SphP arrays.
 *  \param[in] firstnode First note of communication.
 *
 *  \return void
 */
static void particle2in(data_in *in, int i, int firstnode)
{
#ifdef CELL_CENTER_GRAVITY
  if(P[i].Type == 0)
    {
      for(int k = 0; k < 3; k++)
        in->Pos[k] = PS[i].Center[k];
    }
  else
#endif /* #ifdef CELL_CENTER_GRAVITY */
    {
      for(int k = 0; k < 3; k++)
        in->Pos[k] = P[i].Pos[k];
    }

  in->SofteningType = P[i].SofteningType;

  in->Firstnode = firstnode;
}

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
  MyFloat Potential;
} data_out;

static data_out *DataResult, *DataOut;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *  particle and cell data arrays (P, SphP,...)
 *  \param[in] i Index of particle in P and SphP arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *  communicated from other tasks and has to be added locally?
 *
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      PS[i].Potential = out->Potential;
    }
  else /* combine */
    {
      PS[i].Potential += out->Potential;
    }
}

#define USE_SUBCOMM_COMMUNICATOR
#include "../utils/generic_comm_helpers2.h"

static int Num;
static struct unbind_data *d;
static int phase;
static double weakly_bound_limit;

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int i, idx;

  {
    int j, threadid = get_thread_num();

    for(j = 0; j < SubNTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

        idx = NextParticle++;

        if(idx >= Num)
          break;

        i = d[idx].index;

        if(phase == 1)
          if(PS[i].BindingEnergy <= weakly_bound_limit)
            continue;

        subfind_force_treeevaluate_potential(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

/*! \brief Routine that defines what to do with imported particles.
 *
 *  Calls the *_evaluate function in MODE_IMPORTED_PARTICLES.
 *
 *  \return void
 */
static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
  {
    int threadid = get_thread_num();

    while(1)
      {
        i = cnt++;

        if(i >= Nimport)
          break;

        subfind_force_treeevaluate_potential(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*! \brief Computes potential energy.
 *
 *  \param[in] num Number of elements.
 *  \param[in] darg Unbind data.
 *  \param[in] phasearg Which phase are we in? 1:ignore weakly bound particles.
 *  \param[in] weakly_bound_limit_arg Minimum binding energy between two
 *             particles that is accounted for.
 *
 *  \return void
 */
void subfind_potential_compute(int num, struct unbind_data *darg, int phasearg, double weakly_bound_limit_arg)
{
  generic_set_MaxNexport();

  Num                = num;
  d                  = darg;
  phase              = phasearg;
  weakly_bound_limit = weakly_bound_limit_arg;

  generic_comm_pattern(Num, kernel_local, kernel_imported);

  double atime;

  if(All.ComovingIntegrationOn)
    atime = All.Time;
  else
    atime = 1;

  for(int i = 0; i < num; i++)
    {
      if(phase == 1)
        if(PS[d[i].index].BindingEnergy <= weakly_bound_limit)
          continue;

      PS[d[i].index].Potential *= All.G / atime;
    }
}

/*! \brief Evaluate function of potential calculation.
 *
 *  \param[in] target Index of particle/cell/imported data.
 *  \param[in] mode Flag if it operates on local or imported data.
 *  \param[in] threadid ID of thread.
 *
 *  \return 0
 */
static int subfind_force_treeevaluate_potential(int target, int mode, int threadid)
{
  struct NODE *nop = 0;
  int no, numnodes, *firstnode, k;
  double r2, dx, dy, dz, mass, r, u, h_i, h_j, hmax, h_inv, wp;
  double pos_x, pos_y, pos_z;
#ifdef MULTIPLE_NODE_SOFTENING
  struct ExtNODE *extnop = 0;
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
#if !defined(GRAVITY_NOT_PERIODIC)
  double xtmp, ytmp, ztmp;
#endif

  data_in local, *in;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos_x = in->Pos[0];
  pos_y = in->Pos[1];
  pos_z = in->Pos[2];
  h_i   = All.ForceSoftening[in->SofteningType];

  double pot = 0;

  for(k = 0; k < numnodes; k++)
    {
      if(mode == MODE_LOCAL_PARTICLES)
        no = SubTree_MaxPart; /* root node */
      else
        {
          no = firstnode[k];
          no = SubNodes[no].u.d.nextnode; /* open it */
        }

      while(no >= 0)
        {
#ifdef MULTIPLE_NODE_SOFTENING
          int indi_flag1 = -1, indi_flag2 = 0;
#endif                             /* #ifdef MULTIPLE_NODE_SOFTENING */
          if(no < SubTree_MaxPart) /* single particle */
            {
              dx = GRAVITY_NEAREST_X(SubTree_Pos_list[3 * no + 0] - pos_x);
              dy = GRAVITY_NEAREST_Y(SubTree_Pos_list[3 * no + 1] - pos_y);
              dz = GRAVITY_NEAREST_Z(SubTree_Pos_list[3 * no + 2] - pos_z);
              r2 = dx * dx + dy * dy + dz * dz;

              mass = P[no].Mass;

              h_j = All.ForceSoftening[P[no].SofteningType];

              if(h_j > h_i)
                hmax = h_j;
              else
                hmax = h_i;

              no = SubNextnode[no];
            }
          else if(no < SubTree_MaxPart + SubTree_MaxNodes) /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no < SubTree_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the
                                                           branch */
                    break;
                }

              nop  = &SubNodes[no];
              mass = nop->u.d.mass;

              dx = GRAVITY_NEAREST_X(nop->u.d.s[0] - pos_x);
              dy = GRAVITY_NEAREST_Y(nop->u.d.s[1] - pos_y);
              dz = GRAVITY_NEAREST_Z(nop->u.d.s[2] - pos_z);

              r2 = dx * dx + dy * dy + dz * dz;

              /* check Barnes-Hut opening criterion */
              if(nop->len * nop->len > r2 * All.ErrTolThetaSubfind * All.ErrTolThetaSubfind)
                {
                  /* open cell */
                  if(mass)
                    {
                      no = nop->u.d.nextnode;
                      continue;
                    }
                }

              h_j = All.ForceSoftening[nop->u.d.maxsofttype];

              if(h_j > h_i)
                {
#ifdef MULTIPLE_NODE_SOFTENING
#ifdef ADAPTIVE_HYDRO_SOFTENING
                  if(nop->u.d.maxhydrosofttype != nop->u.d.minhydrosofttype)
                    if(SubExtNodes[no].mass_per_type[0] > 0)
                      if(r2 < All.ForceSoftening[nop->u.d.maxhydrosofttype] * All.ForceSoftening[nop->u.d.maxhydrosofttype])
                        {
                          /* open cell */
                          no = nop->u.d.nextnode;
                          continue;
                        }
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
                  indi_flag1 = 0;
                  indi_flag2 = NSOFTTYPES;
#else  /* #ifdef MULTIPLE_NODE_SOFTENING */
                  if(r2 < h_j * h_j)
                    {
                      /* open cell */
                      no = nop->u.d.nextnode;
                      continue;
                    }
#endif /* #ifdef MULTIPLE_NODE_SOFTENING #else */
                  hmax = h_j;
                }
              else
                hmax = h_i;

                /* node can be used */
#ifdef MULTIPLE_NODE_SOFTENING
              extnop = &SubExtNodes[no];
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
              no = nop->u.d.sibling;
            }
          else if(no >= SubTree_ImportedNodeOffset) /* point from imported nodelist */
            {
              terminate("this is not expected here");
            }
          else
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES");

              subfind_treefind_collective_export_node_threads(no, target, threadid);

              no = SubNextnode[no - SubTree_MaxNodes];
              continue;
            }

          /* now evaluate the potential contribution */
          r = sqrt(r2);

#ifdef MULTIPLE_NODE_SOFTENING
          int type;
          for(type = indi_flag1; type < indi_flag2; type++)
            {
              if(type >= 0)
                {
                  mass = extnop->mass_per_type[type];

#ifdef ADAPTIVE_HYDRO_SOFTENING
                  if(type == 0)
                    h_j = All.ForceSoftening[nop->u.d.maxhydrosofttype];
                  else
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
                    h_j = All.ForceSoftening[type];

                  if(h_j > h_i)
                    hmax = h_j;
                  else
                    hmax = h_i;
                }

              if(mass)
                {
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

                  if(r >= hmax)
                    pot += FLT(-mass / r);
                  else
                    {
                      h_inv = 1.0 / hmax;

                      u = r * h_inv;
                      if(u < 0.5)
                        wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
                      else
                        wp = -3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u)));

                      pot += FLT(mass * h_inv * wp);
                    }
#ifdef MULTIPLE_NODE_SOFTENING
                }
            }
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
        }
    }

  out.Potential = pot;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif /* #ifdef SUBFIND */
