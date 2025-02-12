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
 * \file        src/time_integration/timestep_treebased.c
 * \date        05/2018
 * \brief       Algorithm to compute non-local time-step criterion.
 * \details     This is necessary for local time-stepping if material that
 *              would require a short time-step is arriving in cells that
 *              would formally be integrated at a large time-step.
 *              contains functions:
 *                static void particle2in(data_in * in, int i, int firstnode)
 *                static void out2particle(data_out * out, int i, int mode)
 *                static void kernel_local(void)
 *                static void kernel_imported(void)
 *                void tree_based_timesteps(void)
 *                int timestep_evaluate(int target, int mode, int threadid)
 *                void tree_based_timesteps_setsoundspeeds(void)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 11.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef TREE_BASED_TIMESTEPS
/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat Csnd;
  MyFloat cellrad;
  MyFloat CurrentMaxTiStep;

  int Firstnode; /* this is needed as part of the communication alogorithm */
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
  int k;

  for(k = 0; k < 3; k++)
    {
      in->Pos[k] = P[i].Pos[k];
      in->Vel[k] = P[i].Vel[k];
    }

  in->Csnd             = SphP[i].Csnd;
  in->cellrad          = get_cell_radius(i);
  in->CurrentMaxTiStep = SphP[i].CurrentMaxTiStep;

  in->Firstnode = firstnode;
}

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
  MyFloat CurrentMaxTiStep;
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
      SphP[i].CurrentMaxTiStep = out->CurrentMaxTiStep;
    }
  else /* combine */
    {
      if(SphP[i].CurrentMaxTiStep > out->CurrentMaxTiStep)
        SphP[i].CurrentMaxTiStep = out->CurrentMaxTiStep;
    }
}

#include "../utils/generic_comm_helpers2.h"

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int idx;

  /* do local particles */
  {
    int j, threadid = get_thread_num();

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

        idx = NextParticle++;

        if(idx >= TimeBinsHydro.NActiveParticles)
          break;

        int i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;

        if(P[i].Mass == 0 && P[i].ID == 0)
          continue;

        timestep_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        timestep_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*! \brief Main function to call tree-based timesteps.
 *
 *  This function is called in find_timesteps_without_gravity() (timestep.c).
 *
 *  \return void
 */
void tree_based_timesteps(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  tree_based_timesteps_setsoundspeeds();

  generic_set_MaxNexport();

  double t0 = second();

  generic_comm_pattern(TimeBinsHydro.NActiveParticles, kernel_local, kernel_imported);

  double t1 = second();

  mpi_printf("TIMESTEPS: timestep-treewalk: sec=%g\n", timediff(t0, t1));

  CPU_Step[CPU_TREE_TIMESTEPS] += measure_time();
}

/*! \brief The 'core' of the tree-based timestep computation.
 *
 *  A target particle is specified which may either be local, or reside in the
 *  communication buffer.
 *
 *  \param[in] target Index of particle/cell.
 *  \param[in] mode Flag if it operates on local or imported data.
 *  \param[in] threadid ID of thread.
 *
 *  \return cost, i.e. number of nodes that had to be opened.
 */
int timestep_evaluate(int target, int mode, int threadid)
{
  int k, cost = 0, numnodes, *firstnode;
  MyDouble *pos;
  MyFloat *vel;
  double dxp, dxm, dyp, dym, dzp, dzm, pos_m[3], pos_p[3];
  int no, p;
  struct NgbNODE *current;
  double dx, dy, dz, dist, csnd, cellrad, xtmp, ytmp, ztmp;

  data_out out;
  data_in local, *target_data;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos     = target_data->Pos;
  vel     = target_data->Vel;
  csnd    = target_data->Csnd;
  cellrad = target_data->cellrad;

  out.CurrentMaxTiStep = target_data->CurrentMaxTiStep;

  pos_m[0] = pos[0] - boxSize_X;
  pos_p[0] = pos[0] + boxSize_X;
  pos_m[1] = pos[1] - boxSize_Y;
  pos_p[1] = pos[1] + boxSize_Y;
  pos_m[2] = pos[2] - boxSize_Z;
  pos_p[2] = pos[2] + boxSize_Z;

  double atimeinv;
  if(All.ComovingIntegrationOn)
    atimeinv = 1 / All.Time;
  else
    atimeinv = 1.0;

  /* Now start the actual tree-walk computation for this particle */

  for(k = 0; k < numnodes; k++)
    {
      if(mode == MODE_LOCAL_PARTICLES)
        {
          no = Ngb_MaxPart; /* root node */
        }
      else
        {
          no = firstnode[k];
          no = Ngb_Nodes[no].u.d.nextnode; /* open it */
        }

      while(no >= 0)
        {
          cost++;

          if(no < Ngb_MaxPart) /* single particle */
            {
              p  = no;
              no = Ngb_Nextnode[no];

              if(P[p].Type > 0)
                continue;

              if(P[p].Mass == 0 && P[p].ID == 0) /* skip eliminated cells */
                continue;

              if(P[p].Ti_Current != All.Ti_Current)
                {
                  drift_particle(p, All.Ti_Current);
                }

              dx = NEAREST_X(P[p].Pos[0] - pos[0]);
              dy = NEAREST_Y(P[p].Pos[1] - pos[1]);
              dz = NEAREST_Z(P[p].Pos[2] - pos[2]);

              dist = sqrt(dx * dx + dy * dy + dz * dz);

              if(dist > 0)
                {
                  double vsig = csnd + SphP[p].Csnd -
                                ((P[p].Vel[0] - vel[0]) * dx + (P[p].Vel[1] - vel[1]) * dy + (P[p].Vel[2] - vel[2]) * dz) / dist;

                  if(vsig > 0)
                    {
                      dist += cellrad; /* take one cell radius as minimum distance in order to protect against unreasonably small steps
                                          if two mesh-generating points are extremely close */
                      if(out.CurrentMaxTiStep > dist / vsig)
                        out.CurrentMaxTiStep = dist / vsig;
                    }
                }
            }
          else if(no < Ngb_MaxPart + Ngb_MaxNodes) /* internal  */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no <
                     Ngb_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
                    break;
                }

              current = &Ngb_Nodes[no];

              if(current->Ti_Current != All.Ti_Current)
                {
                  drift_node(current, All.Ti_Current);
                }

              if(pos[0] > current->u.d.range_max[0] && pos_m[0] < current->u.d.range_min[0])
                {
                  dxp = pos[0] - current->u.d.range_max[0];
                  dxm = pos_m[0] - current->u.d.range_min[0]; /* negative */
                }
              else if(pos_p[0] > current->u.d.range_max[0] && pos[0] < current->u.d.range_min[0])
                {
                  dxp = pos_p[0] - current->u.d.range_max[0];
                  dxm = pos[0] - current->u.d.range_min[0]; /* negative */
                }
              else
                dxp = dxm = 0;

              if(pos[1] > current->u.d.range_max[1] && pos_m[1] < current->u.d.range_min[1])
                {
                  dyp = pos[1] - current->u.d.range_max[1];
                  dym = pos_m[1] - current->u.d.range_min[1]; /* negative */
                }
              else if(pos_p[1] > current->u.d.range_max[1] && pos[1] < current->u.d.range_min[1])
                {
                  dyp = pos_p[1] - current->u.d.range_max[1];
                  dym = pos[1] - current->u.d.range_min[1]; /* negative */
                }
              else
                dyp = dym = 0;

              if(pos[2] > current->u.d.range_max[2] && pos_m[2] < current->u.d.range_min[2])
                {
                  dzp = pos[2] - current->u.d.range_max[2];
                  dzm = pos_m[2] - current->u.d.range_min[2]; /* negative */
                }
              else if(pos_p[2] > current->u.d.range_max[2] && pos[2] < current->u.d.range_min[2])
                {
                  dzp = pos_p[2] - current->u.d.range_max[2];
                  dzm = pos[2] - current->u.d.range_min[2]; /* negative */
                }
              else
                dzp = dzm = 0;

              double vsig = csnd + ExtNgb_Nodes[no].MaxCsnd;

              int flag = 0;

              if(dxp + cellrad < out.CurrentMaxTiStep * (vsig + (ExtNgb_Nodes[no].vmax[0] * atimeinv - vel[0])))
                flag++;
              else if(-dxm + cellrad < out.CurrentMaxTiStep * (vsig - (ExtNgb_Nodes[no].vmin[0] * atimeinv - vel[0])))
                flag++;

              if(dyp + cellrad < out.CurrentMaxTiStep * (vsig + (ExtNgb_Nodes[no].vmax[1] * atimeinv - vel[1])))
                flag++;
              else if(-dym + cellrad < out.CurrentMaxTiStep * (vsig - (ExtNgb_Nodes[no].vmin[1] * atimeinv - vel[1])))
                flag++;

              if(dzp + cellrad < out.CurrentMaxTiStep * (vsig + (ExtNgb_Nodes[no].vmax[2] * atimeinv - vel[2])))
                flag++;
              else if(-dzm + cellrad < out.CurrentMaxTiStep * (vsig - (ExtNgb_Nodes[no].vmin[2] * atimeinv - vel[2])))
                flag++;

              if(flag >= 3)
                {
                  /* need to open */
                  no = current->u.d.nextnode;
                  continue;
                }

              /* in this case the node can be discarded */
              no = current->u.d.sibling;
              continue;
            }
          else /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == 1");

              if(target >= 0) /* if no target is given, export will not occur */
                ngb_treefind_export_node_threads(no, target, threadid, 0);

              no = Ngb_Nextnode[no - Ngb_MaxNodes];
              continue;
            }
        }
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return cost;
}

/*! \brief Sets local sound speed and time-step limits from local conditions.
 *
 *  This is a sort of initialization of the tree-based time-steps algorithm.
 *
 *  \return void
 */
void tree_based_timesteps_setsoundspeeds(void)
{
  int idx, i;
  double rad, csnd;
  double hubble_a, atime;

  if(All.ComovingIntegrationOn)
    {
      hubble_a = hubble_function(All.Time);
      atime    = All.Time;
    }
  else
    hubble_a = atime = 1;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      csnd = get_sound_speed(i);

      if(csnd <= 1.0e-30)
        csnd = 1.0e-30;

      SphP[i].Csnd = csnd;

      rad = get_cell_radius(i);

#ifdef VORONOI_STATIC_MESH
      csnd += sqrt(P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]) / All.cf_atime;
#else  /* #ifdef VORONOI_STATIC_MESH */
      csnd += sqrt((P[i].Vel[0] - SphP[i].VelVertex[0]) * (P[i].Vel[0] - SphP[i].VelVertex[0]) +
                   (P[i].Vel[1] - SphP[i].VelVertex[1]) * (P[i].Vel[1] - SphP[i].VelVertex[1]) +
                   (P[i].Vel[2] - SphP[i].VelVertex[2]) * (P[i].Vel[2] - SphP[i].VelVertex[2])) /
              All.cf_atime;
#endif /* #ifdef VORONOI_STATIC_MESH #else */

      SphP[i].CurrentMaxTiStep = rad / csnd;

      /* note: for cosmological integration, CurrentMaxTiStep stores  1/a times the maximum allowed physical timestep */

      if(SphP[i].CurrentMaxTiStep >= All.MaxSizeTimestep / (atime * hubble_a) / All.CourantFac)
        SphP[i].CurrentMaxTiStep = All.MaxSizeTimestep / (atime * hubble_a) / All.CourantFac;
    }
}

#endif /* #ifdef TREE_BASED_TIMESTEPS */
