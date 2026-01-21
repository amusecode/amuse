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
 * \file        src/subfind/subfind_coll_tree.c
 * \date        05/2018
 * \brief       Functions for tree-construction for subfind collective.
 * \details     contains functions:
 *                int subfind_coll_treebuild(int npart, struct unbind_data *mp)
 *                int subfind_coll_treebuild_construct(int npart, struct
 *                  unbind_data *mp)
 *                int subfind_coll_treebuild_insert_single_point(int i,
 *                  unsigned long long *intpos, int th, unsigned char levels)
 *                int subfind_coll_create_empty_nodes(int no, int topnode,
 *                  int bits, int x, int y, int z, unsigned long long xc,
 *                  unsigned long long yc, unsigned long long zc,
 *                  unsigned long long ilen)
 *                void subfind_coll_insert_pseudo_particles(void)
 *                void subfind_coll_update_node_recursive(int no, int sib,
 *                  int father, int *last)
 *                void subfind_coll_exchange_topleafdata(void)
 *                void subfind_coll_treeupdate_toplevel(int no, int topnode,
 *                  int bits, int x, int y, int z)
 *                void subfind_coll_treeallocate(int maxpart, int maxindex)
 *                void subfind_coll_treefree(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 04.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef SUBFIND
#include "../gravity/forcetree.h"
#include "subfind.h"

/*! \brief Main function to build subfind collective tree.
 *
 *  \param[in] npart Number of particles.
 *  \param[in] mp Unbind data.
 *
 *  \return Number of nodes in tree.
 */
int subfind_coll_treebuild(int npart, struct unbind_data *mp)
{
  int flag;

  do
    {
      int flag_single = subfind_coll_treebuild_construct(npart, mp);

      MPI_Allreduce(&flag_single, &flag, 1, MPI_INT, MPI_MIN, SubComm);

      if(flag < 0)
        {
          subfind_coll_treefree();

          SubTreeAllocFactor *= 1.15;

          printf("SUBFIND-COLLECTIVE, root-task=%d: Increasing TreeAllocFactor, new value=%g\n", ThisTask, SubTreeAllocFactor);
          fflush(stdout);

          subfind_coll_treeallocate(NumPart, All.MaxPart);
        }
    }
  while(flag < 0);

  /* insert the pseudo particles that represent the mass distribution of other domains */
  subfind_coll_insert_pseudo_particles();

  /* now compute the multipole moments recursively */
  int last = -1;

  subfind_coll_update_node_recursive(SubTree_MaxPart, -1, -1, &last);

  if(last >= SubTree_MaxPart)
    {
      if(last >= SubTree_MaxPart + SubTree_MaxNodes) /* a pseudo-particle or imported particle */
        SubNextnode[last - SubTree_MaxNodes] = -1;
      else
        SubNodes[last].u.d.nextnode = -1;
    }
  else
    SubNextnode[last] = -1;

  subfind_coll_exchange_topleafdata();

  SubTree_NextFreeNode = SubTree_MaxPart + 1;

  subfind_coll_treeupdate_toplevel(SubTree_MaxPart, 0, 1, 0, 0, 0);

  return SubTree_NumNodes;
}

/*! \brief Constructs the collective subfind oct-tree.
 *
 *  The index convention for accessing tree nodes is the following:
 *  node index
 *  [0...SubTree_MaxPart-1]   references single particles, the indices
 *  [SubTree_MaxPart...SubTree_MaxPart+SubTree_MaxNodes-1] references tree
 *  nodes.
 *  [SubTree_MaxPart+SubTree_MaxNodes...
 *  SubTree_MaxPart+SubTree_MaxNodes+NTopleaves-1] references "pseudo
 *  particles", i.e. mark branches on foreign CPUs
 *  [SubTree_MaxPart+SubTree_MaxNodes+NTopleaves...
 *  SubTree_MaxPart+SubTree_MaxNodes+NTopleaves+0-1] references imported points
 *
 *  `Nodes_base' points to the first tree node, while `Nodes' is shifted such
 *  that SubNodes[SubTree_MaxPart] gives the root tree node.
 *
 *  \param[in] npart Number of particles.
 *  \param[in] mp Unbind data.
 *
 *  \return Number of nodes.
 */
int subfind_coll_treebuild_construct(int npart, struct unbind_data *mp)
{
  int i, j, k, no, flag_full = 0;
  unsigned long long *intposp;
  MyDouble *posp;
  unsigned long long ibaselen = ((unsigned long long)1) << 52;

  /* create an empty root node  */
  SubTree_NextFreeNode = SubTree_MaxPart;                 /* index of first free node */
  struct NODE *nfreep  = &SubNodes[SubTree_NextFreeNode]; /* select first node        */

  for(j = 0; j < 8; j++)
    nfreep->u.suns[j] = -1;

  nfreep->len = SubDomainLen;
  for(j = 0; j < 3; j++)
    nfreep->center[j] = SubDomainCenter[j];

  SubTree_NumNodes = 1;
  SubTree_NextFreeNode++;

  /* create a set of empty nodes corresponding to the top-level domain
   * grid. We need to generate these nodes first to make sure that we have a
   * complete top-level tree which allows the easy insertion of the
   * pseudo-particles at the right place
   */
  if(subfind_coll_create_empty_nodes(SubTree_MaxPart, 0, 1, 0, 0, 0, 0, 0, 0, ibaselen) < 0)
    return -1;

  SubTree_FirstNonTopLevelNode = SubTree_NextFreeNode;

  /* if a high-resolution region in a global tree is used, we need to generate
   * an additional set empty nodes to make sure that we have a complete
   * top-level tree for the high-resolution inset
   */

  SubTree_IntPos_list =
      (unsigned long long *)mymalloc_movable(&SubTree_IntPos_list, "SubTree_IntPos_list", 3 * NumPart * sizeof(unsigned long long));

  SubTree_ImportedNodeOffset = SubTree_MaxPart + SubTree_MaxNodes + SubNTopleaves;

  /* now we determine for each point the insertion top-level node, and the task on which this lies */
  for(i = 0; i < npart; i++)
    {
      for(j = 0; j < 3; j++)
        {
          if(mp)
            k = mp[i].index;
          else
            k = i;

#ifdef CELL_CENTER_GRAVITY
          if(P[k].Type == 0)
            posp = &PS[k].Center[j];
          else
#endif /* #ifdef CELL_CENTER_GRAVITY */
            posp = &P[k].Pos[j];

          if(*posp < SubDomainCorner[j] || *posp >= SubDomainCorner[j] + SubDomainLen)
            {
              terminate("out of box i=%d j=%d coord=%g SubDomainCorner=(%g|%g|%g) SubDomainLen=%g", i, j, *posp, SubDomainCorner[0],
                        SubDomainCorner[1], SubDomainCorner[2], SubDomainLen);
            }

          SubTree_Pos_list[3 * k + j] = *posp;
        }
    }

  for(i = 0; i < npart; i++)
    {
      if(mp)
        k = mp[i].index;
      else
        k = i;

      posp = &SubTree_Pos_list[3 * k];

      unsigned long long xxb  = force_double_to_int(((*posp++ - SubDomainCorner[0]) * SubDomainInverseLen) + 1.0);
      unsigned long long yyb  = force_double_to_int(((*posp++ - SubDomainCorner[1]) * SubDomainInverseLen) + 1.0);
      unsigned long long zzb  = force_double_to_int(((*posp++ - SubDomainCorner[2]) * SubDomainInverseLen) + 1.0);
      unsigned long long mask = ((unsigned long long)1) << (52 - 1);
      unsigned char shiftx    = (52 - 1);
      unsigned char shifty    = (52 - 2);
      unsigned char shiftz    = (52 - 3);
      unsigned char levels    = 0;

      intposp = &SubTree_IntPos_list[3 * k];

      *intposp++ = xxb;
      *intposp++ = yyb;
      *intposp++ = zzb;

      no = 0;
      while(SubTopNodes[no].Daughter >= 0)
        {
          unsigned char subnode = (((unsigned char)((xxb & mask) >> (shiftx--))) | ((unsigned char)((yyb & mask) >> (shifty--))) |
                                   ((unsigned char)((zzb & mask) >> (shiftz--))));

          mask >>= 1;
          levels++;

          no = SubTopNodes[no].Daughter + SubTopNodes[no].MortonToPeanoSubnode[subnode];
        }

      no = SubTopNodes[no].Leaf;

      if(no >= SubTree_ImportedNodeOffset)
        terminate("i=%d: no=%d SubTree_ImportedNodeOffset=%d", i, no, SubTree_ImportedNodeOffset);

      if(subfind_coll_treebuild_insert_single_point(k, &SubTree_IntPos_list[3 * k], SubDomainNodeIndex[no], levels) < 0)
        {
          flag_full = 1;
          break;
        }
    }

  myfree_movable(SubTree_IntPos_list);

  if(flag_full)
    return -1;

  return SubTree_NumNodes;
}

/*! \brief Inserts single point in tree.
 *
 *  \param[in] i Index of particle.
 *  \param[in] intpos Integer position.
 *  \param[in] th Index in SubNodes.
 *  \param[in] levels Level corresponding to subnode.
 *
 *  \return void
 */
int subfind_coll_treebuild_insert_single_point(int i, unsigned long long *intpos, int th, unsigned char levels)
{
  int j, parent = -1;
  unsigned char subnode       = 0;
  unsigned long long xxb      = intpos[0];
  unsigned long long yyb      = intpos[1];
  unsigned long long zzb      = intpos[2];
  unsigned long long mask     = ((unsigned long long)1) << ((52 - 1) - levels);
  unsigned char shiftx        = (52 - 1) - levels;
  unsigned char shifty        = (52 - 2) - levels;
  unsigned char shiftz        = (52 - 3) - levels;
  signed long long centermask = (0xFFF0000000000000llu);
  unsigned long long *intppos;
  centermask >>= levels;

  while(1)
    {
      if(th >= SubTree_MaxPart && th < SubTree_ImportedNodeOffset) /* we are dealing with an internal node */
        {
          subnode = (((unsigned char)((xxb & mask) >> (shiftx--))) | ((unsigned char)((yyb & mask) >> (shifty--))) |
                     ((unsigned char)((zzb & mask) >> (shiftz--))));

          centermask >>= 1;
          mask >>= 1;
          levels++;

          if(levels > MAX_TREE_LEVEL)
            {
              /* seems like we're dealing with particles at identical (or extremely close)
               * locations. Shift subnode index to allow tree construction. Note: Multipole moments
               * of tree are still correct, but one should MAX_TREE_LEVEL large enough to have
               *      DomainLen/2^MAX_TREE_LEEL  < gravitational softening length
               */
              for(j = 0; j < 8; j++)
                {
                  if(SubNodes[th].u.suns[subnode] < 0)
                    break;

                  subnode++;
                  if(subnode >= 8)
                    subnode = 7;
                }
            }

          int nn = SubNodes[th].u.suns[subnode];

          if(nn >= 0) /* ok, something is in the daughter slot already, need to continue */
            {
              parent = th;
              th     = nn;
            }
          else
            {
              /* here we have found an empty slot where we can attach
               * the new particle as a leaf.
               */
              SubNodes[th].u.suns[subnode] = i;
              break; /* done for this particle */
            }
        }
      else
        {
          /* We try to insert into a leaf with a single particle.  Need
           * to generate a new internal node at this point.
           */
          SubNodes[parent].u.suns[subnode] = SubTree_NextFreeNode;
          struct NODE *nfreep              = &SubNodes[SubTree_NextFreeNode];

          /* the other is: */
          double len = ((double)(mask << 1)) * SubDomainBigFac;
          double cx  = ((double)((xxb & centermask) | mask)) * SubDomainBigFac + SubDomainCorner[0];
          double cy  = ((double)((yyb & centermask) | mask)) * SubDomainBigFac + SubDomainCorner[1];
          double cz  = ((double)((zzb & centermask) | mask)) * SubDomainBigFac + SubDomainCorner[2];

          nfreep->len       = len;
          nfreep->center[0] = cx;
          nfreep->center[1] = cy;
          nfreep->center[2] = cz;

          for(j = 0; j < 8; j++)
            nfreep->u.suns[j] = -1;

          if(th >= SubTree_ImportedNodeOffset)
            {
              terminate("unexpected here: th=%d SubTree_ImportedNodeOffset=%d", th, SubTree_ImportedNodeOffset);
            }
          else
            intppos = &SubTree_IntPos_list[3 * th];

          subnode = (((unsigned char)((intppos[0] & mask) >> shiftx)) | ((unsigned char)((intppos[1] & mask) >> shifty)) |
                     ((unsigned char)((intppos[2] & mask) >> shiftz)));

          nfreep->u.suns[subnode] = th;

          th = SubTree_NextFreeNode; /* resume trying to insert the new particle the newly created internal node */
          SubTree_NumNodes++;
          SubTree_NextFreeNode++;

          if(SubTree_NumNodes >= SubTree_MaxNodes)
            {
              if(SubTreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
                {
                  char buf[500];
                  sprintf(buf,
                          "task %d: looks like a serious problem for particle %d, stopping with particle dump.  SubTree_NumNodes=%d "
                          "SubTree_MaxNodes=%d  0=%d NumPart=%d\n",
                          SubThisTask, i, SubTree_NumNodes, SubTree_MaxNodes, 0, NumPart);
                  dump_particles();
                  terminate(buf);
                }

              return -1;
            }
        }
    }

  return 0;
}

/*! \brief Recursively creates a set of empty tree nodes which corresponds to
 *         the top-level tree for the domain grid. This is done to ensure that
 *         this top-level tree is always "complete" so that we can easily
 *         associate the pseudo-particles of other CPUs with tree-nodes at a
 *         given level in the tree, even when the particle population is so
 *         sparse that some of these nodes are actually empty.
 *
 *  \param[in] no Index of node.
 *  \param[in] topnode Index of topnode.
 *  \param[in] bits Number of bits used for Peano key.
 *  \param[in] x Integer x position.
 *  \param[in] y Integer y position.
 *  \param[in] z Integer z position.
 *  \param[in] xc X position of corner.
 *  \param[in] yc Y position of corner.
 *  \param[in] zc Z position of corner.
 *  \param[in] ilen Sidelength.
 *
 *  \return 0: success; -1 Number of nodes exceeds maximum number of nodes.
 */
int subfind_coll_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z, unsigned long long xc, unsigned long long yc,
                                    unsigned long long zc, unsigned long long ilen)
{
  int i, j, k, n, sub, count;
  unsigned long long xxc, yyc, zzc, ilenhalf;

  ilen >>= 1;

  if(SubTopNodes[topnode].Daughter >= 0)
    {
      for(i = 0; i < 2; i++)
        for(j = 0; j < 2; j++)
          for(k = 0; k < 2; k++)
            {
              if(SubTree_NumNodes >= SubTree_MaxNodes)
                {
                  if(SubTreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
                    {
                      char buf[500];
                      sprintf(buf, "task %d: looks like a serious problem (NTopnodes=%d), stopping with particle dump.\n", SubThisTask,
                              NTopnodes);
                      dump_particles();
                      terminate(buf);
                    }
                  return -1;
                }

              sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);

              count = i + 2 * j + 4 * k;

              SubNodes[no].u.suns[count] = SubTree_NextFreeNode;

              xxc      = xc + i * ilen;
              yyc      = yc + j * ilen;
              zzc      = zc + k * ilen;
              ilenhalf = ilen >> 1;

              double len = ((double)ilen) * SubDomainBigFac;
              double cx  = ((double)(xxc + ilenhalf)) * SubDomainBigFac + SubDomainCorner[0];
              double cy  = ((double)(yyc + ilenhalf)) * SubDomainBigFac + SubDomainCorner[1];
              double cz  = ((double)(zzc + ilenhalf)) * SubDomainBigFac + SubDomainCorner[2];

              SubNodes[SubTree_NextFreeNode].len       = len;
              SubNodes[SubTree_NextFreeNode].center[0] = cx;
              SubNodes[SubTree_NextFreeNode].center[1] = cy;
              SubNodes[SubTree_NextFreeNode].center[2] = cz;

              for(n = 0; n < 8; n++)
                SubNodes[SubTree_NextFreeNode].u.suns[n] = -1;

              if(SubTopNodes[SubTopNodes[topnode].Daughter + sub].Daughter == -1)
                SubDomainNodeIndex[SubTopNodes[SubTopNodes[topnode].Daughter + sub].Leaf] = SubTree_NextFreeNode;

              SubTree_NextFreeNode++;
              SubTree_NumNodes++;

              if(subfind_coll_create_empty_nodes(SubTree_NextFreeNode - 1, SubTopNodes[topnode].Daughter + sub, bits + 1, 2 * x + i,
                                                 2 * y + j, 2 * z + k, xxc, yyc, zzc, ilen) < 0)
                return -1;
            }
    }

  return 0;
}

/*! \brief Inserts pseudo-particles which will represent the mass
 *         distribution of the other CPUs. Initially, the mass of the
 *         pseudo-particles is set to zero, and their coordinate is set to the
 *        center of the domain-cell they correspond to. These quantities will
 *        be updated later on.
 *
 *  \return void
 */
void subfind_coll_insert_pseudo_particles(void)
{
  int i, index;

  for(i = 0; i < SubNTopleaves; i++)
    {
      index = SubDomainNodeIndex[i];

      if(SubDomainTask[i] != SubThisTask)
        SubNodes[index].u.suns[0] = SubTree_MaxPart + SubTree_MaxNodes + i;
    }
}

/*! \brief Determines the multipole moments for a given internal node
 *         and all its subnodes using a recursive computation.  The result is
 *         stored in the SubNodes structure in the sequence of this tree-walk.
 *
 *  \param[in] no Index of node.
 *  \param[in] sib Index of sibling.
 *  \param[in] father Index of parent node.
 *  \param[in, out] last Node index of last call.
 *
 *  \return void
 */
void subfind_coll_update_node_recursive(int no, int sib, int father, int *last)
{
  int j, jj, p, pp, nextsib, suns[8];
  double s[3], mass;
  unsigned char maxsofttype;
#ifdef MULTIPLE_NODE_SOFTENING
  double mass_per_type[NSOFTTYPES];
#ifdef ADAPTIVE_HYDRO_SOFTENING
  unsigned char maxhydrosofttype;
  unsigned char minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

  if(no >= SubTree_MaxPart && no < SubTree_MaxPart + SubTree_MaxNodes) /* internal node */
    {
      for(j = 0; j < 8; j++)
        suns[j] = SubNodes[no].u.suns[j]; /* this "backup" is necessary because the nextnode entry will
                                             overwrite one element (union!) */
      if(*last >= 0)
        {
          if(*last >= SubTree_MaxPart)
            {
              if(*last >= SubTree_MaxPart + SubTree_MaxNodes)
                SubNextnode[*last - SubTree_MaxNodes] = no; /* a pseudo-particle or imported point */
              else
                SubNodes[*last].u.d.nextnode = no;
            }
          else
            SubNextnode[*last] = no;
        }

      *last = no;

      mass        = 0;
      s[0]        = 0;
      s[1]        = 0;
      s[2]        = 0;
      maxsofttype = NSOFTTYPES + NSOFTTYPES_HYDRO;

#ifdef MULTIPLE_NODE_SOFTENING
      for(j = 0; j < NSOFTTYPES; j++)
        mass_per_type[j] = 0;

#ifdef ADAPTIVE_HYDRO_SOFTENING
      maxhydrosofttype = NSOFTTYPES;
      minhydrosofttype = NSOFTTYPES + NSOFTTYPES_HYDRO - 1;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

      for(j = 0; j < 8; j++)
        {
          if((p = suns[j]) >= 0)
            {
              /* check if we have a sibling on the same level */
              for(jj = j + 1; jj < 8; jj++)
                if((pp = suns[jj]) >= 0)
                  break;

              if(jj < 8) /* yes, we do */
                nextsib = pp;
              else
                nextsib = sib;

              subfind_coll_update_node_recursive(p, nextsib, no, last);

              if(p < SubTree_MaxPart) /* a particle */
                {
                  MyDouble *pos = &SubTree_Pos_list[3 * p];

                  mass += P[p].Mass;
                  s[0] += P[p].Mass * pos[0];
                  s[1] += P[p].Mass * pos[1];
                  s[2] += P[p].Mass * pos[2];

                  if(All.ForceSoftening[maxsofttype] < All.ForceSoftening[P[p].SofteningType])
                    maxsofttype = P[p].SofteningType;

#ifdef MULTIPLE_NODE_SOFTENING
#ifdef ADAPTIVE_HYDRO_SOFTENING
                  mass_per_type[P[p].Type == 0 ? 0 : P[p].SofteningType] += P[p].Mass;

                  if(P[p].Type == 0)
                    {
                      if(maxhydrosofttype < P[p].SofteningType)
                        maxhydrosofttype = P[p].SofteningType;
                      if(minhydrosofttype > P[p].SofteningType)
                        minhydrosofttype = P[p].SofteningType;
                    }
#else  /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
                  mass_per_type[P[p].SofteningType] += P[p].Mass;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING #else */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
                }
              else if(p < SubTree_MaxPart + SubTree_MaxNodes) /* an internal node  */
                {
                  mass += SubNodes[p].u.d.mass;
                  s[0] += SubNodes[p].u.d.mass * SubNodes[p].u.d.s[0];
                  s[1] += SubNodes[p].u.d.mass * SubNodes[p].u.d.s[1];
                  s[2] += SubNodes[p].u.d.mass * SubNodes[p].u.d.s[2];

                  if(All.ForceSoftening[maxsofttype] < All.ForceSoftening[SubNodes[p].u.d.maxsofttype])
                    maxsofttype = SubNodes[p].u.d.maxsofttype;

#ifdef MULTIPLE_NODE_SOFTENING
                  int k;
                  for(k = 0; k < NSOFTTYPES; k++)
                    mass_per_type[k] += SubExtNodes[p].mass_per_type[k];

#ifdef ADAPTIVE_HYDRO_SOFTENING
                  if(maxhydrosofttype < SubNodes[p].u.d.maxhydrosofttype)
                    maxhydrosofttype = SubNodes[p].u.d.maxhydrosofttype;
                  if(minhydrosofttype > SubNodes[p].u.d.minhydrosofttype)
                    minhydrosofttype = SubNodes[p].u.d.minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
                }
              else if(p < SubTree_MaxPart + SubTree_MaxNodes + SubNTopleaves) /* a pseudo particle */
                {
                  /* nothing to be done here because the mass of the
                   *  pseudo-particle is still zero. This will be changed
                   * later.
                   */
                }
              else
                {
                  /* an imported point */
                  terminate("should not occur here");
                }
            }
        }

      if(mass)
        {
          s[0] /= mass;
          s[1] /= mass;
          s[2] /= mass;
        }
      else
        {
          s[0] = SubNodes[no].center[0];
          s[1] = SubNodes[no].center[1];
          s[2] = SubNodes[no].center[2];
        }

      SubNodes[no].u.d.mass        = mass;
      SubNodes[no].u.d.s[0]        = s[0];
      SubNodes[no].u.d.s[1]        = s[1];
      SubNodes[no].u.d.s[2]        = s[2];
      SubNodes[no].u.d.maxsofttype = maxsofttype;
      SubNodes[no].u.d.sibling     = sib;
      SubNodes[no].u.d.father      = father;

#ifdef MULTIPLE_NODE_SOFTENING
      int k;
      for(k = 0; k < NSOFTTYPES; k++)
        SubExtNodes[no].mass_per_type[k] = mass_per_type[k];

#ifdef ADAPTIVE_HYDRO_SOFTENING
      SubNodes[no].u.d.maxhydrosofttype = maxhydrosofttype;
      SubNodes[no].u.d.minhydrosofttype = minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
    }
  else /* single particle or pseudo particle */
    {
      if(*last >= 0)
        {
          if(*last >= SubTree_MaxPart)
            {
              if(*last >= SubTree_MaxPart + SubTree_MaxNodes)
                SubNextnode[*last - SubTree_MaxNodes] = no; /* a pseudo-particle or an imported point */
              else
                SubNodes[*last].u.d.nextnode = no;
            }
          else
            SubNextnode[*last] = no;
        }

      *last = no;
    }
}

/*! \brief This function communicates the values of the multipole moments of
 *         the top-level tree-nodes of the domain grid.  This data can then be
 *         used to update the pseudo-particles on each CPU accordingly.
 *
 *  \return void
 */
void subfind_coll_exchange_topleafdata(void)
{
  int n, no, idx, task;
  int *recvcounts, *recvoffset, *bytecounts, *byteoffset;
  struct DomainNODE
  {
    MyFloat s[3];
    MyFloat mass;
#ifdef MULTIPLE_NODE_SOFTENING
    MyDouble mass_per_type[NSOFTTYPES];
#ifdef ADAPTIVE_HYDRO_SOFTENING
    unsigned char maxhydrosofttype;
    unsigned char minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
    unsigned char maxsofttype;
  } * DomainMoment, *loc_DomainMoment;

  DomainMoment = (struct DomainNODE *)mymalloc("DomainMoment", SubNTopleaves * sizeof(struct DomainNODE));

  /* share the pseudo-particle data accross CPUs */
  recvcounts = (int *)mymalloc("recvcounts", sizeof(int) * SubNTask);
  recvoffset = (int *)mymalloc("recvoffset", sizeof(int) * SubNTask);
  bytecounts = (int *)mymalloc("bytecounts", sizeof(int) * SubNTask);
  byteoffset = (int *)mymalloc("byteoffset", sizeof(int) * SubNTask);

  for(task = 0; task < SubNTask; task++)
    recvcounts[task] = 0;

  for(n = 0; n < SubNTopleaves; n++)
    {
      if(SubDomainTask[n] < 0 || SubDomainTask[n] >= SubNTask)
        terminate("n=%d|%d: SubDomainTask[n]=%d", n, SubNTopleaves, SubDomainTask[n]);

      recvcounts[SubDomainTask[n]]++;
    }

  for(task = 0; task < SubNTask; task++)
    bytecounts[task] = recvcounts[task] * sizeof(struct DomainNODE);

  for(task = 1, recvoffset[0] = 0, byteoffset[0] = 0; task < SubNTask; task++)
    {
      recvoffset[task] = recvoffset[task - 1] + recvcounts[task - 1];
      byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];
    }

  loc_DomainMoment = (struct DomainNODE *)mymalloc("loc_DomainMoment", recvcounts[SubThisTask] * sizeof(struct DomainNODE));

  for(n = 0, idx = 0; n < SubNTopleaves; n++)
    {
      if(SubDomainTask[n] == SubThisTask)
        {
          no = SubDomainNodeIndex[n];

          /* read out the multipole moments from the local base cells */
          loc_DomainMoment[idx].s[0]        = SubNodes[no].u.d.s[0];
          loc_DomainMoment[idx].s[1]        = SubNodes[no].u.d.s[1];
          loc_DomainMoment[idx].s[2]        = SubNodes[no].u.d.s[2];
          loc_DomainMoment[idx].mass        = SubNodes[no].u.d.mass;
          loc_DomainMoment[idx].maxsofttype = SubNodes[no].u.d.maxsofttype;
#ifdef MULTIPLE_NODE_SOFTENING
          int k;
          for(k = 0; k < NSOFTTYPES; k++)
            loc_DomainMoment[idx].mass_per_type[k] = SubExtNodes[no].mass_per_type[k];

#ifdef ADAPTIVE_HYDRO_SOFTENING
          loc_DomainMoment[idx].maxhydrosofttype = SubNodes[no].u.d.maxhydrosofttype;
          loc_DomainMoment[idx].minhydrosofttype = SubNodes[no].u.d.minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
          idx++;
        }
    }

  MPI_Allgatherv(loc_DomainMoment, bytecounts[SubThisTask], MPI_BYTE, DomainMoment, bytecounts, byteoffset, MPI_BYTE, SubComm);

  for(task = 0; task < SubNTask; task++)
    recvcounts[task] = 0;

  for(n = 0; n < SubNTopleaves; n++)
    {
      task = SubDomainTask[n];
      if(task != SubThisTask)
        {
          no  = SubDomainNodeIndex[n];
          idx = recvoffset[task] + recvcounts[task]++;

          SubNodes[no].u.d.s[0]        = DomainMoment[idx].s[0];
          SubNodes[no].u.d.s[1]        = DomainMoment[idx].s[1];
          SubNodes[no].u.d.s[2]        = DomainMoment[idx].s[2];
          SubNodes[no].u.d.mass        = DomainMoment[idx].mass;
          SubNodes[no].u.d.maxsofttype = DomainMoment[idx].maxsofttype;
#ifdef MULTIPLE_NODE_SOFTENING
          int k;
          for(k = 0; k < NSOFTTYPES; k++)
            SubExtNodes[no].mass_per_type[k] = DomainMoment[idx].mass_per_type[k];
#ifdef ADAPTIVE_HYDRO_SOFTENING
          SubNodes[no].u.d.maxhydrosofttype = DomainMoment[idx].maxhydrosofttype;
          SubNodes[no].u.d.minhydrosofttype = DomainMoment[idx].minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
        }
    }

  myfree(loc_DomainMoment);
  myfree(byteoffset);
  myfree(bytecounts);
  myfree(recvoffset);
  myfree(recvcounts);
  myfree(DomainMoment);
}

/*! \brief This function updates the top-level tree after the multipole
 *         moments of the pseudo-particles have been updated.
 *
 *  \param[in] no Index of node.
 *  \param[in] topnode Index of topnode.
 *  \param[in] bits Number of bits used.
 *  \param[in] x Integer x position.
 *  \param[in] y Integer y position.
 *  \param[in] z Integer z position.
 *
 *  \return void
 */
void subfind_coll_treeupdate_toplevel(int no, int topnode, int bits, int x, int y, int z)
{
  int i, j, k, sub;
  int p;
  double s[3], mass;
  unsigned char maxsofttype;
#ifdef MULTIPLE_NODE_SOFTENING
  double mass_per_type[NSOFTTYPES];
#ifdef ADAPTIVE_HYDRO_SOFTENING
  unsigned char maxhydrosofttype;
  unsigned char minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

  if(SubTopNodes[topnode].Daughter >= 0)
    {
      for(i = 0; i < 2; i++)
        for(j = 0; j < 2; j++)
          for(k = 0; k < 2; k++)
            {
              sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);

              SubTree_NextFreeNode++;
              subfind_coll_treeupdate_toplevel(SubTree_NextFreeNode - 1, SubTopNodes[topnode].Daughter + sub, bits + 1, 2 * x + i,
                                               2 * y + j, 2 * z + k);
            }

      mass        = 0;
      s[0]        = 0;
      s[1]        = 0;
      s[2]        = 0;
      maxsofttype = NSOFTTYPES + NSOFTTYPES_HYDRO;
#ifdef MULTIPLE_NODE_SOFTENING
      for(j = 0; j < NSOFTTYPES; j++)
        mass_per_type[j] = 0;

#ifdef ADAPTIVE_HYDRO_SOFTENING
      maxhydrosofttype = NSOFTTYPES;
      minhydrosofttype = NSOFTTYPES + NSOFTTYPES_HYDRO - 1;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

      p = SubNodes[no].u.d.nextnode;

      for(j = 0; j < 8; j++) /* since we are dealing with top-level nodes, we know that there are 8 consecutive daughter nodes */
        {
          if(p >= SubTree_MaxPart && p < SubTree_MaxPart + SubTree_MaxNodes) /* internal node */
            {
              mass += SubNodes[p].u.d.mass;
              s[0] += SubNodes[p].u.d.mass * SubNodes[p].u.d.s[0];
              s[1] += SubNodes[p].u.d.mass * SubNodes[p].u.d.s[1];
              s[2] += SubNodes[p].u.d.mass * SubNodes[p].u.d.s[2];
              if(All.ForceSoftening[maxsofttype] < All.ForceSoftening[SubNodes[p].u.d.maxsofttype])
                maxsofttype = SubNodes[p].u.d.maxsofttype;
#ifdef MULTIPLE_NODE_SOFTENING
              int k;
              for(k = 0; k < NSOFTTYPES; k++)
                mass_per_type[k] += SubExtNodes[p].mass_per_type[k];

#ifdef ADAPTIVE_HYDRO_SOFTENING
              if(maxhydrosofttype < SubNodes[p].u.d.maxhydrosofttype)
                maxhydrosofttype = SubNodes[p].u.d.maxhydrosofttype;
              if(minhydrosofttype > SubNodes[p].u.d.minhydrosofttype)
                minhydrosofttype = SubNodes[p].u.d.minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
            }
          else
            terminate("may not happen");

          p = SubNodes[p].u.d.sibling;
        }

      if(mass)
        {
          s[0] /= mass;
          s[1] /= mass;
          s[2] /= mass;
        }
      else
        {
          s[0] = SubNodes[no].center[0];
          s[1] = SubNodes[no].center[1];
          s[2] = SubNodes[no].center[2];
        }

      SubNodes[no].u.d.s[0]        = s[0];
      SubNodes[no].u.d.s[1]        = s[1];
      SubNodes[no].u.d.s[2]        = s[2];
      SubNodes[no].u.d.mass        = mass;
      SubNodes[no].u.d.maxsofttype = maxsofttype;
#ifdef MULTIPLE_NODE_SOFTENING
      int k;
      for(k = 0; k < NSOFTTYPES; k++)
        SubExtNodes[no].mass_per_type[k] = mass_per_type[k];
#ifdef ADAPTIVE_HYDRO_SOFTENING
      SubNodes[no].u.d.maxhydrosofttype = maxhydrosofttype;
      SubNodes[no].u.d.minhydrosofttype = minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
    }
}

/*! \brief Allocates tree arrays.
 *
 *  This function allocates the memory used for storage of the tree nodes.
 *  Usually, the number of required nodes is of order 0.7*maxpart, but if this
 *  is insufficient, the code will try to allocated more space.
 *
 *  \param[in] maxpart Maximum number of nodes.
 *  \param[in] maxindex Maximum number of particles.
 *
 *  \return void
 */
void subfind_coll_treeallocate(int maxpart, int maxindex)
{
  if(SubNodes)
    terminate("already allocated");

  SubTree_MaxPart  = maxindex;
  SubTree_MaxNodes = (int)(SubTreeAllocFactor * maxpart) + SubNTopnodes;

  SubDomainNodeIndex = (int *)mymalloc_movable(&SubDomainNodeIndex, "SubDomainNodeIndex", SubNTopleaves * sizeof(int));

  SubTree_Pos_list = (MyDouble *)mymalloc_movable(&SubTree_Pos_list, "SubTree_Pos_list", 3 * maxpart * sizeof(MyDouble));

  SubNodes = (struct NODE *)mymalloc_movable(&SubNodes, "SubNodes", (SubTree_MaxNodes + 1) * sizeof(struct NODE));
  SubNodes -= SubTree_MaxPart;

#ifdef MULTIPLE_NODE_SOFTENING
  SubExtNodes = (struct ExtNODE *)mymalloc_movable(&SubExtNodes, "SubExtNodes", (SubTree_MaxNodes + 1) * sizeof(struct ExtNODE));
  SubExtNodes -= SubTree_MaxPart;
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

  SubNextnode = (int *)mymalloc_movable(&SubNextnode, "SubNextnode", (SubTree_MaxPart + SubNTopleaves) * sizeof(int));
}

/*! \brief Free tree arrays.
 *
 *  This function frees the memory allocated for the tree, i.e. it frees
 *  the space allocated by the function subfind_coll_treeallocate().
 *
 *  \return void
 */
void subfind_coll_treefree(void)
{
  if(SubNodes)
    {
      myfree(SubNextnode);

#ifdef MULTIPLE_NODE_SOFTENING
      myfree(SubExtNodes + SubTree_MaxPart);
      SubExtNodes = NULL;
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

      myfree(SubNodes + SubTree_MaxPart);
      myfree(SubTree_Pos_list);
      myfree(SubDomainNodeIndex);

      SubNodes           = NULL;
      SubDomainNodeIndex = NULL;
      SubNextnode        = NULL;
      SubTree_Pos_list   = NULL;
    }
  else
    terminate("trying to free the tree even though it's not allocated");
}

#endif /* #ifdef SUBFIND */
