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
 * \file        src/ngbtree/ngbtree_walk.c
 * \date        05/2018
 * \brief       Routines to walk the ngb tree.
 * \details     contains functions:
 *                int ngb_treefind_variable_threads(MyDouble searchcenter[3],
 *                  MyFloat hsml, int target, int mode, int thread_id, int
 *                  numnodes, int *firstnode)
 *                int ngb_treefind_export_node_threads(int no, int target, int
 *                  thread_id, int image_flag)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 16.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/*! \brief Finds all cells around seearchcenter in region with radius hsml.
 *
 *  This function returns the number of neighbors with distance <= hsml, and
 *  returns the particle indices in the global buffer Ngblist.
 *  The tree traversal starts at startnode.
 *  Keep in mind that this is usually called within an *_evaluate function
 *  within the generic communication pattern. This means that first, the local
 *  (bound to this task) search is performed and the local neighbors written
 *  to the array, then communication happens and afterwards, the function is
 *  called again in imported mode, finding particles on other tasks.
 *
 *  \param[in] searchcenter Center of the neighbor search.
 *  \param[in] hsml Radius of the search.
 *  \param[in] target Index of the particle around which the search is
 *             performed; needed for parallel search. If < 0, only local search
 *             is performed.
 *  \param[in] mode Mode for local or imported particle search.
 *  \param[in] thread_id ID of thread (always 0 in our case).
 *  \param[in] numnodes Number of nodes on this task (1 for mode local;
 *             for mode imported: given by generic_get_numnodes(...) ).
 *  \param[in] firstnode Node to start with (in case of mode imported).
 *
 *  \return The number of neighbors found.
 */
int ngb_treefind_variable_threads(MyDouble searchcenter[3], MyFloat hsml, int target, int mode, int thread_id, int numnodes,
                                  int *firstnode)
{
  MyDouble search_min[3], search_max[3], search_max_Lsub[3], search_min_Ladd[3];

  for(int i = 0; i < 3; i++)
    {
      search_min[i] = searchcenter[i] - 1.001 * hsml;
      search_max[i] = searchcenter[i] + 1.001 * hsml;
    }

  search_max_Lsub[0] = search_max[0] - boxSize_X;
  search_max_Lsub[1] = search_max[1] - boxSize_Y;
  search_max_Lsub[2] = search_max[2] - boxSize_Z;

  search_min_Ladd[0] = search_min[0] + boxSize_X;
  search_min_Ladd[1] = search_min[1] + boxSize_Y;
  search_min_Ladd[2] = search_min[2] + boxSize_Z;

  int numngb = 0;
  double xtmp, ytmp, ztmp;
  double hsml2 = hsml * hsml;

  for(int k = 0; k < numnodes; k++)
    {
      int no;

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
          if(no < Ngb_MaxPart) /* single particle */
            {
              int p = no;
              no    = Ngb_Nextnode[no];

              if(P[p].Type > 0)
                continue;

              if(P[p].Ti_Current != All.Ti_Current)
                {
                  drift_particle(p, All.Ti_Current);
                }

              double dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
              if(dx > hsml)
                continue;
              double dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
              if(dy > hsml)
                continue;
              double dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
              if(dz > hsml)
                continue;

              double r2 = dx * dx + dy * dy + dz * dz;
              if(r2 > hsml2)
                continue;

              Thread[thread_id].R2list[numngb]    = r2;
              Thread[thread_id].Ngblist[numngb++] = p;
            }
          else if(no < Ngb_MaxPart + Ngb_MaxNodes) /* internal node */
            {
              struct NgbNODE *current = &Ngb_Nodes[no];

              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no <
                     Ngb_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
                    break;
                }

              no = current->u.d.sibling; /* in case the node can be discarded */

              if(current->Ti_Current != All.Ti_Current)
                {
                  drift_node(current, All.Ti_Current);
                }

              if(search_min[0] > current->u.d.range_max[0] && search_max_Lsub[0] < current->u.d.range_min[0])
                continue;
              if(search_min_Ladd[0] > current->u.d.range_max[0] && search_max[0] < current->u.d.range_min[0])
                continue;

              if(search_min[1] > current->u.d.range_max[1] && search_max_Lsub[1] < current->u.d.range_min[1])
                continue;
              if(search_min_Ladd[1] > current->u.d.range_max[1] && search_max[1] < current->u.d.range_min[1])
                continue;

              if(search_min[2] > current->u.d.range_max[2] && search_max_Lsub[2] < current->u.d.range_min[2])
                continue;
              if(search_min_Ladd[2] > current->u.d.range_max[2] && search_max[2] < current->u.d.range_min[2])
                continue;

              no = current->u.d.nextnode; /* ok, we need to open the node */
            }
          else /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES should not occur here");

              if(target >= 0) /* if no target is given, export will not occur */
                if(ngb_treefind_export_node_threads(no, target, thread_id, 0))
                  return -1;

              no = Ngb_Nextnode[no - Ngb_MaxNodes];
              continue;
            }
        }
    }
  return numngb;
}

/*! \brief Prepares export of ngb-tree node.
 *
 *  \param[in] no Pseudoparticle node to be exported.
 *  \param[in] target (Local) index to identify what it refers to.
 *  \param[in] thread_id ID of thread (0 in our case).
 *  \param[in] image_flag Bit flag used in EXTENDED_GHOST_SEARCH.
 *
 *  \return 0
 */
int ngb_treefind_export_node_threads(int no, int target, int thread_id, int image_flag)
{
  /* The task indicated by the pseudoparticle node */
  int task = DomainTask[no - (Ngb_MaxPart + Ngb_MaxNodes)];

  if(Thread[thread_id].Exportflag[task] != target)
    {
      Thread[thread_id].Exportflag[task]     = target;
      int nexp                               = Thread[thread_id].Nexport++;
      Thread[thread_id].PartList[nexp].Task  = task;
      Thread[thread_id].PartList[nexp].Index = target;
      Thread[thread_id].ExportSpace -= Thread[thread_id].ItemSize;
    }

  int nexp                      = Thread[thread_id].NexportNodes++;
  nexp                          = -1 - nexp;
  struct datanodelist *nodelist = (struct datanodelist *)(((char *)Thread[thread_id].PartList) + Thread[thread_id].InitialSpace);
  nodelist[nexp].Task           = task;
  nodelist[nexp].Index          = target;
  nodelist[nexp].Node           = Ngb_DomainNodeIndex[no - (Ngb_MaxPart + Ngb_MaxNodes)];
#ifdef EXTENDED_GHOST_SEARCH
  nodelist[nexp].BitFlags = image_flag;
#endif /* #ifdef EXTENDED_GHOST_SEARCH */
  Thread[thread_id].ExportSpace -= sizeof(struct datanodelist) + sizeof(int);
  return 0;
}
