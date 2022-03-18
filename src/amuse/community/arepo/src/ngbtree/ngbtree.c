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
 * \file        src/ngbtree/ngbtree.c
 * \date        05/2018
 * \brief       Construct neighbor tree.
 * \details     This file contains the neighbor tree construction. This is a
 *              tree structure that includes all gas cells, but no other
 *              particle types.
 *              contains functions:
 *                int ngb_treebuild(int npart)
 *                static inline unsigned long long ngb_double_to_int(double d)
 *                int ngb_treebuild_construct(int npart)
 *                int ngb_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z)
 *                void ngb_update_node_recursive(int no, int sib, int father, int *last, int mode)
 *                void ngb_record_topnode_siblings(int no, int sib)
 *                void ngb_exchange_topleafdata(void)
 *                void drift_node(struct NgbNODE *current, integertime time1)
 *                void ngb_update_velocities(void)
 *                void ngb_update_vbounds(int i, int *nchanged, int *nodelist)
 *                void ngb_finish_vounds_update(int nchanged, int *nodelist)
 *                void ngb_update_rangebounds(int i, int *nchanged, int *nodelist)
 *                void ngb_finish_rangebounds_update(int nchanged, int *nodelist)
 *                void ngb_treemodifylength(int delta_NgbMaxPart)
 *                void ngb_treeallocate(void)
 *                void ngb_treefree(void)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 21.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"
#include "../gravity/forcetree.h"

static void ngb_record_topnode_siblings(int no, int sib);
static int ngb_treebuild_construct(int npart);
static void ngb_update_node_recursive(int no, int sib, int father, int *last, int mode);
static void ngb_exchange_topleafdata(void);
static int ngb_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z);
static void ngb_update_vbounds(int i, int *nchanged, int *nodelist);
static void ngb_finish_vounds_update(int nchanged, int *nodelist);

static int *Ngb_Node_Tmp_Sibling;

/*! \brief This function is a driver routine for constructing the neighbor
 *         oct-tree, which is done by calling a small number of other
 *         functions.
 *
 *  Does not build a tree if All.TotNumGas == 0.
 *
 *  \param[in] npart Number of particles in tree.
 *
 *  \return Number of nodes in the tree.
 */
int ngb_treebuild(int npart)
{
  if(All.TotNumGas == 0)
    return 0;

  TIMER_START(CPU_NGBTREEBUILD);

  mpi_printf("NGBTREE: Ngb-tree construction.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  double t0 = second();

  int flag;
  do
    {
      int flag_single = ngb_treebuild_construct(npart);

      MPI_Allreduce(&flag_single, &flag, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      if(flag == -1)
        {
          myfree(Ngb_Node_Tmp_Sibling + Ngb_MaxPart);
          ngb_treefree();

          All.NgbTreeAllocFactor *= 1.15;
          mpi_printf("Increasing NgbTreeAllocFactor, new value=%g\n", All.NgbTreeAllocFactor);

          ngb_treeallocate();
        }
    }
  while(flag == -1);

  int ntopleaves = DomainNLocalTopleave[ThisTask];
  int *list      = DomainListOfLocalTopleaves + DomainFirstLocTopleave[ThisTask];

  for(int i = 0; i < ntopleaves; i++)
    {
      int last = -1;
      int no   = Ngb_DomainNodeIndex[list[i]];

      if(no < Ngb_MaxPart || no >= Ngb_MaxPart + Ngb_MaxNodes)
        terminate("i=%d no=%d  task=%d \n", i, no, DomainTask[list[i]]);

      ngb_update_node_recursive(no, Ngb_Node_Tmp_Sibling[no], no, &last, 0);

      /* if there was no particle in the node, we need to initialize nextnode of the node */
      if(no == last)
        Ngb_Nodes[no].u.d.nextnode = -1;

      Ngb_Nodes[no].u.d.sibling = last; /* we temporarily store this here and will later restore this sibling pointer,
                                           which is anyway equal to Ngb_Node_Tmp_Sibling[index] */
    }

  ngb_exchange_topleafdata();

  /* now put in "pseudo" particles as nextnode in non-local topleaves */
  for(int i = 0; i < NTopleaves; i++)
    {
      if(DomainTask[i] != ThisTask)
        {
          int index                     = Ngb_DomainNodeIndex[i];
          Ngb_Nodes[index].u.d.nextnode = Ngb_MaxPart + Ngb_MaxNodes + i;
        }
    }

  /* now update the top-level tree nodes */
  int last = -1;
  ngb_update_node_recursive(Ngb_MaxPart, -1, -1, &last, 1);

  if(last >= Ngb_MaxPart)
    {
      if(last >= Ngb_MaxPart + Ngb_MaxNodes) /* a pseudo-particle */
        Ngb_Nextnode[last - Ngb_MaxNodes] = -1;
      else
        Ngb_Nodes[last].u.d.nextnode = -1;
    }
  else
    Ngb_Nextnode[last] = -1;

  TIMER_STOPSTART(CPU_NGBTREEBUILD, CPU_LOGS);

  double numnodes = Ngb_NumNodes, tot_numnodes;
  MPI_Reduce(&numnodes, &tot_numnodes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  double t1 = second();
  mpi_printf("NGBTREE: Ngb-tree construction done. took %g sec  <numnodes>=%g  NTopnodes=%d NTopleaves=%d\n", timediff(t0, t1),
             tot_numnodes / NTask, NTopnodes, NTopleaves);

  myfree(Ngb_Node_Tmp_Sibling + Ngb_MaxPart);

  Ngb_MarkerValue = 0;
  memset(Ngb_Marker, 0, (Ngb_MaxPart + Ngb_NumNodes) * sizeof(int));

  TIMER_STOP(CPU_LOGS);

  return Ngb_NumNodes;
}

/*! \brief Converts double precision coordinate to unsigned long long int.
 *
 *  \param[in] d Double precision coordinate that is to be converted.
 *
 *  \return Unsigned long long int represenation of d.
 */
static inline unsigned long long ngb_double_to_int(double d)
{
  union
  {
    double d;
    unsigned long long ull;
  } u;
  u.d = d;
  return (u.ull & 0xFFFFFFFFFFFFFllu);
}

/*! \brief Constructs the neighbor oct-tree.
 *
 *  The index convention for accessing tree nodes is the following:
 *
 *  0...NumPart-1                             reference single particles.
 *  Ngb_MaxPart.... Ngb_MaxPart+Numnodes-1    references tree nodes.
 *  Ngb_MaxPart + All.MaxNgb_Nodes....                reference "pseudo
 *     particles", i.e. the marker that indicates a top-node lying on
 *     another CPU.
 *
 *  `Ngb_Nodes_base' points to the first tree node,
 *  `Ngb_Nodes' is shifted such that Ngb_Nodes[Ngb_MaxPart] gives the first
 *     tree node.
 *
 *  \param[in] npart Number of particles involved.
 *
 *  \return status: 0 (default) -1: too many nodes.
 */
int ngb_treebuild_construct(int npart)
{
  /* create an empty root node  */
  Ngb_NextFreeNode = Ngb_MaxPart; /* index of first free node */

  for(int i = 0; i < 8; i++)
    Ngb_Nodes[Ngb_NextFreeNode].u.suns[i] = -1;

  Ngb_NumNodes = 1;
  Ngb_NextFreeNode++;

  /* create a set of empty nodes corresponding to the top-level domain
   * grid. We need to generate these nodes first to make sure that we have a
   * complete top-level tree which allows the easy insertion of the
   * pseudo-particles at the right place
   */
  if(ngb_create_empty_nodes(Ngb_MaxPart, 0, 1, 0, 0, 0) < 0)
    return -1;

  Ngb_FirstNonTopLevelNode = Ngb_NextFreeNode;

  Ngb_Node_Tmp_Sibling = (int *)mymalloc("Ngb_Node_Tmp_Sibling", (Ngb_MaxNodes + 1) * sizeof(int));
  Ngb_Node_Tmp_Sibling -= Ngb_MaxPart;

  ngb_record_topnode_siblings(Ngb_MaxPart, -1);

  unsigned long long *ngbTree_IntPos_list =
      (unsigned long long *)mymalloc("ngbTree_IntPos_list", 3 * npart * sizeof(unsigned long long));

  /* now we insert all particles */
  {
    int out_of_space = 0;

    int threadid = get_thread_num();
    int start, end, size;

    int first_empty_slot = Ngb_NextFreeNode + threadid * TAKE_NSLOTS_IN_ONE_GO;
    int count_empty_slot = TAKE_NSLOTS_IN_ONE_GO;

    if(threadid == 0)
      Ngb_NextFreeNode += NUM_THREADS * TAKE_NSLOTS_IN_ONE_GO;

    size  = (npart - 1) / NUM_THREADS + 1;
    start = threadid * size;
    end   = (threadid + 1) * size - 1;
    if(end >= npart)
      end = npart - 1;

    for(int i = start; i <= end && out_of_space == 0; i++)
      {
        unsigned long long xxb  = ngb_double_to_int(((P[i].Pos[0] - DomainCorner[0]) * DomainInverseLen) + 1.0);
        unsigned long long yyb  = ngb_double_to_int(((P[i].Pos[1] - DomainCorner[1]) * DomainInverseLen) + 1.0);
        unsigned long long zzb  = ngb_double_to_int(((P[i].Pos[2] - DomainCorner[2]) * DomainInverseLen) + 1.0);
        unsigned long long mask = ((unsigned long long)1) << (52 - 1);
        unsigned char shiftx    = (52 - 1);
        unsigned char shifty    = (52 - 2);
        unsigned char shiftz    = (52 - 3);
        unsigned char levels    = 0;

        ngbTree_IntPos_list[3 * i + 0] = xxb;
        ngbTree_IntPos_list[3 * i + 1] = yyb;
        ngbTree_IntPos_list[3 * i + 2] = zzb;

        int no = 0;
        while(TopNodes[no].Daughter >= 0) /* walk down top tree to find correct leaf */
          {
            unsigned char subnode = (((unsigned char)((xxb & mask) >> (shiftx--))) | ((unsigned char)((yyb & mask) >> (shifty--))) |
                                     ((unsigned char)((zzb & mask) >> (shiftz--))));

            mask >>= 1;
            levels++;

            no = TopNodes[no].Daughter + TopNodes[no].MortonToPeanoSubnode[subnode];
          }

        no = TopNodes[no].Leaf;

        if(DomainTask[no] != ThisTask)
          terminate("STOP!  ID=%lld of type=%d is inserted into task=%d, but should be on task=%d no=%d\n", (long long)P[i].ID,
                    P[i].Type, ThisTask, DomainTask[no], no);

        int th = Ngb_DomainNodeIndex[no];

        signed long long centermask = (0xFFF0000000000000llu) >> levels;

        int parent            = -1; /* note: will not be used below before it is changed */
        unsigned char subnode = 0;

        while(1)
          {
            if(th >= Ngb_MaxPart) /* we are dealing with an internal node */
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
                    for(int j = 0; j < 8; j++)
                      {
                        if(Ngb_Nodes[th].u.suns[subnode] < 0)
                          break;

                        subnode++;
                        if(subnode >= 8)
                          subnode = 7;
                      }
                  }

                int nn = Ngb_Nodes[th].u.suns[subnode];

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
                    Ngb_Nodes[th].u.suns[subnode] = i;
                    break; /* done for this particle */
                  }
              }
            else
              {
                /* We try to insert into a leaf with a single particle.  Need
                 * to generate a new internal node at this point.
                 * Then resume trying to insert the new particle at
                 * the newly created internal node
                 */
                int thold = th;

                if(count_empty_slot)
                  {
                    th = first_empty_slot + (TAKE_NSLOTS_IN_ONE_GO - count_empty_slot);
                    count_empty_slot--;
                  }
                else
                  {
                    {
                      th = Ngb_NextFreeNode;
                      Ngb_NextFreeNode += TAKE_NSLOTS_IN_ONE_GO;
                    }

                    first_empty_slot = th;
                    count_empty_slot = (TAKE_NSLOTS_IN_ONE_GO - 1);

                    if(first_empty_slot + TAKE_NSLOTS_IN_ONE_GO - Ngb_MaxPart >= Ngb_MaxNodes)
                      {
                        out_of_space = 1;
                        break;
                      }
                  }

                Ngb_Nodes[parent].u.suns[subnode] = th;
                struct NgbNODE *nfreep            = &Ngb_Nodes[th];

                for(int j = 0; j < 8; j++)
                  nfreep->u.suns[j] = -1;

                unsigned long long *intppos = &ngbTree_IntPos_list[3 * thold];

                subnode = (((unsigned char)((intppos[0] & mask) >> shiftx)) | ((unsigned char)((intppos[1] & mask) >> shifty)) |
                           ((unsigned char)((intppos[2] & mask) >> shiftz)));

                nfreep->u.suns[subnode] = thold;
              }
          }
      }
  }

  myfree(ngbTree_IntPos_list);

  if((Ngb_NumNodes = Ngb_NextFreeNode - Ngb_MaxPart) >= Ngb_MaxNodes)
    {
      if(All.NgbTreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
        {
          dump_particles();
          terminate("task %d: out of space for neighbor tree, stopping with particle dump.\n", ThisTask);
        }
      else
        return -1;
    }

  return 0;
}

/*! \brief Create empty ngb-tree node.
 *
 *  This function recursively creates a set of empty tree nodes which
 *  corresponds to the top-level tree for the domain grid. This is done to
 *  ensure that this top-level tree is always "complete" so that we can easily
 *  associate the pseudo-particles of other CPUs with tree-nodes at a given
 *  level in the tree, even when the particle population is so sparse that
 *  some of these nodes are actually empty.
 *
 *  \param[in] no Index of node in Ngb_Nodes array.
 *  \param[in] topnode Index in TopNodes.
 *  \param[in] bits Number of bits used.
 *  \param[in] x Integer coordinate X.
 *  \param[in] y Integer coordinate Y.
 *  \param[in] z Integer coordinate Z.
 *
 *  \return Status: 0 success; -1 error.
 */
int ngb_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z)
{
  if(TopNodes[topnode].Daughter >= 0)
    {
      for(int i = 0; i < 2; i++)
        for(int j = 0; j < 2; j++)
          for(int k = 0; k < 2; k++)
            {
              if(Ngb_NumNodes >= Ngb_MaxNodes)
                {
                  if(All.NgbTreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
                    {
                      dump_particles();
                      terminate("task %d: looks like a serious problem (NTopnodes=%d), stopping with particle dump.\n", ThisTask,
                                NTopnodes);
                    }
                  return -1;
                }

              int sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);

              int count = i + 2 * j + 4 * k;

              Ngb_Nodes[no].u.suns[count] = Ngb_NextFreeNode;

              for(int n = 0; n < 8; n++)
                Ngb_Nodes[Ngb_NextFreeNode].u.suns[n] = -1;

              if(TopNodes[TopNodes[topnode].Daughter + sub].Daughter == -1)
                Ngb_DomainNodeIndex[TopNodes[TopNodes[topnode].Daughter + sub].Leaf] = Ngb_NextFreeNode;

              Ngb_NextFreeNode++;
              Ngb_NumNodes++;

              if(ngb_create_empty_nodes(Ngb_NextFreeNode - 1, TopNodes[topnode].Daughter + sub, bits + 1, 2 * x + i, 2 * y + j,
                                        2 * z + k) < 0)
                return -1;
            }
    }

  return 0;
}

/*! \brief Determine node ranges.
 *
 *  This routine determines the node ranges a given internal node
 *  and all its subnodes using a recursive computation.  The result is
 *  stored in the Ngb_Nodes[] structure in the sequence of this tree-walk.
 *
 *
 *  \param[in] no Index of node.
 *  \param[in] sib Sibling node of no.
 *  \param[in] father Parent node of no.
 *  \param[in, out] last Pointer to last node for which this function was
 *                  called.
 *  \param[in] mode 0: process a leave branch; 1: process top-level nodes.
 *
 *  \return void
 */
void ngb_update_node_recursive(int no, int sib, int father, int *last, int mode)
{
  int j, jj, k, p, pp, nextsib, suns[8];
  MyNgbTreeFloat range_min[3];
  MyNgbTreeFloat range_max[3];
  MyNgbTreeFloat vertex_vmin[3];
  MyNgbTreeFloat vertex_vmax[3];
#ifdef TREE_BASED_TIMESTEPS
  MyNgbTreeFloat vmin[3], vmax[3], maxcsnd;
#endif /* #ifdef TREE_BASED_TIMESTEPS */

  if(no >= Ngb_MaxPart && no < Ngb_MaxPart + Ngb_MaxNodes) /* internal node */
    {
      if(*last >= 0)
        {
          if(*last >= Ngb_MaxPart)
            {
              if(*last == no)
                terminate("as");

              if(*last >= Ngb_MaxPart + Ngb_MaxNodes) /* a pseudo-particle */
                Ngb_Nextnode[*last - Ngb_MaxNodes] = no;
              else
                Ngb_Nodes[*last].u.d.nextnode = no;
            }
          else
            Ngb_Nextnode[*last] = no;
        }

      *last = no;

      int not_interal_top_level = 0;

      if(mode == 1)
        {
          if(!(no >= Ngb_MaxPart && no < Ngb_FirstNonTopLevelNode))
            terminate("can't be");

          if(Ngb_Node_Tmp_Sibling[no] != -2)
            not_interal_top_level = 1;
        }

      if(not_interal_top_level)
        {
          p = Ngb_Nodes[no].u.d.nextnode;

          if(p >= Ngb_MaxPart + Ngb_MaxNodes &&
             p < Ngb_MaxPart + Ngb_MaxNodes + NTopleaves) /* a pseudo-particle, i.e. we are dealing with a non-local top-leave */
            ngb_update_node_recursive(p, sib, no, last, mode);
          else
            {
              /* this is local toplevel node */
              *last = Ngb_Nodes[no].u.d.sibling;
            }

          if(Ngb_Node_Tmp_Sibling[no] != sib)
            terminate("Ngb_Node_Tmp_Sibling[no] != sib");

          /* restore the sibling pointer for local toplevel nodes (we had temporarily stored the last element in this branch */
          Ngb_Nodes[no].u.d.sibling = sib;
          Ngb_Nodes[no].father      = father;
        }
      else
        {
          for(j = 0; j < 8; j++)
            suns[j] = Ngb_Nodes[no].u.suns[j]; /* this "backup" is necessary because the nextnode entry will
                                                  overwrite one element (union!) */

#ifdef TREE_BASED_TIMESTEPS
          maxcsnd = 0;
#endif /* #ifdef TREE_BASED_TIMESTEPS */
          for(k = 0; k < 3; k++)
            {
              range_min[k] = MAX_NGBRANGE_NUMBER;
              range_max[k] = -MAX_NGBRANGE_NUMBER;

              vertex_vmin[k] = MAX_NGBRANGE_NUMBER;
              vertex_vmax[k] = -MAX_NGBRANGE_NUMBER;

#ifdef TREE_BASED_TIMESTEPS
              vmin[k] = MAX_NGBRANGE_NUMBER;
              vmax[k] = -MAX_NGBRANGE_NUMBER;
#endif /* #ifdef TREE_BASED_TIMESTEPS */
            }

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

                  ngb_update_node_recursive(p, nextsib, no, last, mode);

                  if(p >= Ngb_MaxPart) /* an internal node or pseudo particle */
                    {
                      if(p >= Ngb_MaxPart + Ngb_MaxNodes) /* a pseudo particle */
                        {
                          /* nothing to be done here because the mass of the
                           * pseudo-particle is still zero. This will be changed
                           * later.
                           */
                        }
                      else
                        {
#ifdef TREE_BASED_TIMESTEPS
                          if(maxcsnd < ExtNgb_Nodes[p].MaxCsnd)
                            maxcsnd = ExtNgb_Nodes[p].MaxCsnd;
#endif /* #ifdef TREE_BASED_TIMESTEPS */
                          for(k = 0; k < 3; k++)
                            {
                              if(range_min[k] > Ngb_Nodes[p].u.d.range_min[k])
                                range_min[k] = Ngb_Nodes[p].u.d.range_min[k];

                              if(range_max[k] < Ngb_Nodes[p].u.d.range_max[k])
                                range_max[k] = Ngb_Nodes[p].u.d.range_max[k];

                              if(vertex_vmin[k] > Ngb_Nodes[p].vertex_vmin[k])
                                vertex_vmin[k] = Ngb_Nodes[p].vertex_vmin[k];

                              if(vertex_vmax[k] < Ngb_Nodes[p].vertex_vmax[k])
                                vertex_vmax[k] = Ngb_Nodes[p].vertex_vmax[k];

#ifdef TREE_BASED_TIMESTEPS
                              if(vmin[k] > ExtNgb_Nodes[p].vmin[k])
                                vmin[k] = ExtNgb_Nodes[p].vmin[k];

                              if(vmax[k] < ExtNgb_Nodes[p].vmax[k])
                                vmax[k] = ExtNgb_Nodes[p].vmax[k];
#endif /* #ifdef TREE_BASED_TIMESTEPS */
                            }
                        }
                    }
                  else /* a particle */
                    {
#ifdef TREE_BASED_TIMESTEPS
                      if(maxcsnd < SphP[p].Csnd)
                        maxcsnd = SphP[p].Csnd;
#endif /* #ifdef TREE_BASED_TIMESTEPS */
                      for(k = 0; k < 3; k++)
                        {
                          if(range_min[k] > P[p].Pos[k])
                            range_min[k] = P[p].Pos[k];

                          if(range_max[k] < P[p].Pos[k])
                            range_max[k] = P[p].Pos[k];

                          if(P[p].Type == 0)
                            {
                              if(vertex_vmin[k] > SphP[p].VelVertex[k])
                                vertex_vmin[k] = SphP[p].VelVertex[k];

                              if(vertex_vmax[k] < SphP[p].VelVertex[k])
                                vertex_vmax[k] = SphP[p].VelVertex[k];
                            }

#ifdef TREE_BASED_TIMESTEPS
                          if(vmin[k] > P[p].Vel[k])
                            vmin[k] = P[p].Vel[k];

                          if(vmax[k] < P[p].Vel[k])
                            vmax[k] = P[p].Vel[k];
#endif /* #ifdef TREE_BASED_TIMESTEPS */
                        }
                    }
                }
            }

#ifdef TREE_BASED_TIMESTEPS
          ExtNgb_Nodes[no].MaxCsnd = maxcsnd;
#endif /* #ifdef TREE_BASED_TIMESTEPS */

          for(k = 0; k < 3; k++)
            {
              Ngb_Nodes[no].u.d.range_min[k] = range_min[k];
              Ngb_Nodes[no].u.d.range_max[k] = range_max[k];
              Ngb_Nodes[no].vertex_vmin[k]   = vertex_vmin[k];
              Ngb_Nodes[no].vertex_vmax[k]   = vertex_vmax[k];
#ifdef TREE_BASED_TIMESTEPS
              ExtNgb_Nodes[no].vmin[k] = vmin[k];
              ExtNgb_Nodes[no].vmax[k] = vmax[k];
#endif /* #ifdef TREE_BASED_TIMESTEPS */
            }

          Ngb_Nodes[no].u.d.sibling = sib;
          Ngb_Nodes[no].father      = father;

          Ngb_Nodes[no].Ti_Current = All.Ti_Current;
        }
    }
  else /* single particle or pseudo particle */
    {
      if(*last >= 0)
        {
          if(*last >= Ngb_MaxPart)
            {
              if(*last >= Ngb_MaxPart + Ngb_MaxNodes) /* a pseudo-particle */
                Ngb_Nextnode[*last - Ngb_MaxNodes] = no;
              else
                Ngb_Nodes[*last].u.d.nextnode = no;
            }
          else
            {
              Ngb_Nextnode[*last] = no;
            }
        }
      if(no < Ngb_MaxPart) /* only set it for single particles... */
        {
          if(father < Ngb_MaxPart)
            terminate("no=%d father=%d\n", no, father);

          Ngb_Father[no] = father;
        }

      *last = no;
    }
}

/*! \brief Sets sibling information in u.suns for node no.
 *
 *  \param[in] no Index of node.
 *  \param[in] sib Index of sibling.
 *
 *  \return void
 */
void ngb_record_topnode_siblings(int no, int sib)
{
  /* note: when this routine is called, only toplevel tree nodes are present */

  if(Ngb_Nodes[no].u.suns[0] >= 0)
    {
      /* marker value to designate internal nodes in the top-level tree */
      Ngb_Node_Tmp_Sibling[no] = -2;

      if(Ngb_Nodes[no].u.suns[0] >= 0)
        for(int j = 0; j < 8; j++)
          {
            int p = Ngb_Nodes[no].u.suns[j];
            int nextsib;

            if(j < 7)
              nextsib = Ngb_Nodes[no].u.suns[j + 1];
            else
              nextsib = sib;

            ngb_record_topnode_siblings(p, nextsib);
          }
    }
  else
    Ngb_Node_Tmp_Sibling[no] = sib; /* a top-level leave node */
}

/*! \brief Communicates top leaf data.
 *
 *  \return void
 */
void ngb_exchange_topleafdata(void)
{
  struct DomainNODE
  {
    MyNgbTreeFloat range_min[3];
    MyNgbTreeFloat range_max[3];
    MyNgbTreeFloat vertex_vmin[3];
    MyNgbTreeFloat vertex_vmax[3];
#ifdef TREE_BASED_TIMESTEPS
    MyNgbTreeFloat MaxCsnd, vmin[3], vmax[3];
#endif /* #ifdef TREE_BASED_TIMESTEPS */
  };

  struct DomainNODE *DomainMoment = (struct DomainNODE *)mymalloc("DomainMoment", NTopleaves * sizeof(struct DomainNODE));

  /* share the pseudo-particle data accross CPUs */
  int *recvcounts = (int *)mymalloc("recvcounts", sizeof(int) * NTask);
  int *recvoffset = (int *)mymalloc("recvoffset", sizeof(int) * NTask);
  int *bytecounts = (int *)mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *)mymalloc("byteoffset", sizeof(int) * NTask);

  for(int task = 0; task < NTask; task++)
    recvcounts[task] = 0;

  for(int n = 0; n < NTopleaves; n++)
    recvcounts[DomainTask[n]]++;

  for(int task = 0; task < NTask; task++)
    bytecounts[task] = recvcounts[task] * sizeof(struct DomainNODE);

  recvoffset[0] = 0, byteoffset[0] = 0;
  for(int task = 1; task < NTask; task++)
    {
      recvoffset[task] = recvoffset[task - 1] + recvcounts[task - 1];
      byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];
    }

  struct DomainNODE *loc_DomainMoment =
      (struct DomainNODE *)mymalloc("loc_DomainMoment", recvcounts[ThisTask] * sizeof(struct DomainNODE));

  int idx = 0;
  for(int n = 0; n < NTopleaves; n++)
    {
      if(DomainTask[n] == ThisTask)
        {
          int no = Ngb_DomainNodeIndex[n];

          /* read out the multipole moments from the local base cells */
#ifdef TREE_BASED_TIMESTEPS
          loc_DomainMoment[idx].MaxCsnd = ExtNgb_Nodes[no].MaxCsnd;
#endif /* #ifdef TREE_BASED_TIMESTEPS */
          for(int k = 0; k < 3; k++)
            {
              loc_DomainMoment[idx].range_min[k]   = Ngb_Nodes[no].u.d.range_min[k];
              loc_DomainMoment[idx].range_max[k]   = Ngb_Nodes[no].u.d.range_max[k];
              loc_DomainMoment[idx].vertex_vmin[k] = Ngb_Nodes[no].vertex_vmin[k];
              loc_DomainMoment[idx].vertex_vmax[k] = Ngb_Nodes[no].vertex_vmax[k];
#ifdef TREE_BASED_TIMESTEPS
              loc_DomainMoment[idx].vmin[k] = ExtNgb_Nodes[no].vmin[k];
              loc_DomainMoment[idx].vmax[k] = ExtNgb_Nodes[no].vmax[k];
#endif /* #ifdef TREE_BASED_TIMESTEPS */
            }
          idx++;
        }
    }

  MPI_Allgatherv(loc_DomainMoment, bytecounts[ThisTask], MPI_BYTE, DomainMoment, bytecounts, byteoffset, MPI_BYTE, MPI_COMM_WORLD);

  for(int task = 0; task < NTask; task++)
    recvcounts[task] = 0;

  for(int n = 0; n < NTopleaves; n++)
    {
      int task = DomainTask[n];
      if(task != ThisTask)
        {
          int no  = Ngb_DomainNodeIndex[n];
          int idx = recvoffset[task] + recvcounts[task]++;

#ifdef TREE_BASED_TIMESTEPS
          ExtNgb_Nodes[no].MaxCsnd = DomainMoment[idx].MaxCsnd;
#endif /* #ifdef TREE_BASED_TIMESTEPS */
          for(int k = 0; k < 3; k++)
            {
              Ngb_Nodes[no].u.d.range_min[k] = DomainMoment[idx].range_min[k];
              Ngb_Nodes[no].u.d.range_max[k] = DomainMoment[idx].range_max[k];
              Ngb_Nodes[no].vertex_vmin[k]   = DomainMoment[idx].vertex_vmin[k];
              Ngb_Nodes[no].vertex_vmax[k]   = DomainMoment[idx].vertex_vmax[k];
#ifdef TREE_BASED_TIMESTEPS
              ExtNgb_Nodes[no].vmin[k] = DomainMoment[idx].vmin[k];
              ExtNgb_Nodes[no].vmax[k] = DomainMoment[idx].vmax[k];
#endif /* #ifdef TREE_BASED_TIMESTEPS */
            }
          Ngb_Nodes[no].Ti_Current = All.Ti_Current;
        }
    }

  myfree(loc_DomainMoment);
  myfree(byteoffset);
  myfree(bytecounts);
  myfree(recvoffset);
  myfree(recvcounts);
  myfree(DomainMoment);
}

/*! \brief Drifts a node to time time1.
 *
 *  \param[in] current Current node.
 *  \param[in] time1 Time to be drifted to.
 *
 *  \return void
 */
void drift_node(struct NgbNODE *current, integertime time1)
{
  double dt_drift;

  if(All.ComovingIntegrationOn)
    dt_drift = get_drift_factor(current->Ti_Current, time1);
  else
    dt_drift = (time1 - current->Ti_Current) * All.Timebase_interval;

  for(int j = 0; j < 3; j++)
    {
      current->u.d.range_min[j] += current->vertex_vmin[j] * dt_drift;
      current->u.d.range_max[j] += current->vertex_vmax[j] * dt_drift;
    }

  current->Ti_Current = time1;
}

/*! \brief Updates velocity informataion in ngb node data.
 *
 *  \return void
 */
void ngb_update_velocities(void)
{
  TIMER_START(CPU_NGBTREEUPDATEVEL);

  Ngb_MarkerValue++;

  int nchanged  = 0;
  int *nodelist = (int *)mymalloc("nodelist", NTopleaves * sizeof(int));

  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int target = TimeBinsHydro.ActiveParticleList[idx];
      if(target >= 0)
        if(P[target].Type == 0)
          ngb_update_vbounds(target, &nchanged, nodelist);
    }

  for(int timebin = All.HighestSynchronizedTimeBin; timebin >= 0; timebin--)
    {
      for(int target = TimeBinsGravity.FirstInTimeBin[timebin]; target >= 0; target = TimeBinsGravity.NextInTimeBin[target])
        if(target >= 0)
          if(P[target].Type == 0)
            ngb_update_vbounds(target, &nchanged, nodelist);
    }

  ngb_finish_vounds_update(nchanged, nodelist);

  myfree(nodelist);

  TIMER_STOP(CPU_NGBTREEUPDATEVEL);
}

/*! \brief Updates vmin and vmax in ngb nodes.
 *
 *  Inverse tree walk.
 *
 *  \param[in] i Index of particle.
 *  \param[in, out] nchanged Number of changed top level nodes.
 *  \param[out] nodelist Top level nodes that were changed.
 *
 *  \return void
 */
void ngb_update_vbounds(int i, int *nchanged, int *nodelist)
{
  int no = Ngb_Father[i];

  while(no >= 0)
    {
      if(Ngb_Nodes[no].Ti_Current != All.Ti_Current)
        drift_node(&Ngb_Nodes[no], All.Ti_Current);

      int flag_changed = 0;

      for(int j = 0; j < 3; j++)
        {
          if(Ngb_Nodes[no].vertex_vmin[j] > SphP[i].VelVertex[j])
            {
              Ngb_Nodes[no].vertex_vmin[j] = SphP[i].VelVertex[j];
              flag_changed                 = 1;
            }

          if(Ngb_Nodes[no].vertex_vmax[j] < SphP[i].VelVertex[j])
            {
              Ngb_Nodes[no].vertex_vmax[j] = SphP[i].VelVertex[j];
              flag_changed                 = 1;
            }

#ifdef TREE_BASED_TIMESTEPS
          if(ExtNgb_Nodes[no].vmin[j] > P[i].Vel[j])
            {
              ExtNgb_Nodes[no].vmin[j] = P[i].Vel[j];
              flag_changed             = 1;
            }

          if(ExtNgb_Nodes[no].vmax[j] < P[i].Vel[j])
            {
              ExtNgb_Nodes[no].vmax[j] = P[i].Vel[j];
              flag_changed             = 1;
            }
#endif /* #ifdef TREE_BASED_TIMESTEPS */
        }

      if(flag_changed == 0)
        break;

      if(no < Ngb_FirstNonTopLevelNode) /* top-level tree-node reached */
        {
          if(Ngb_Marker[no] != Ngb_MarkerValue)
            {
              Ngb_Marker[no]      = Ngb_MarkerValue;
              nodelist[*nchanged] = no;
              *nchanged           = *nchanged + 1;
            }
          break;
        }

      no = Ngb_Nodes[no].father;
    }
}

/*! \brief Finalizes velocity bounds update.
 *
 *  Exchanges changed information in top level nodes to all tasks.
 *
 *  \param[in] nchanged Number of changed top level nodes.
 *  \param[in] list of changed top level nodes
 *
 *  \return void
 */
void ngb_finish_vounds_update(int nchanged, int *nodelist)
{
  struct DomainNODE
  {
    int node;
    MyNgbTreeFloat vertex_vmin[3];
    MyNgbTreeFloat vertex_vmax[3];
#ifdef TREE_BASED_TIMESTEPS
    MyNgbTreeFloat vmin[3];
    MyNgbTreeFloat vmax[3];
#endif /* #ifdef TREE_BASED_TIMESTEPS */
  };

  /* share the pseudo-particle data accross CPUs */
  int *recvcounts = (int *)mymalloc("recvcounts", sizeof(int) * NTask);
  int *bytecounts = (int *)mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *)mymalloc("byteoffset", sizeof(int) * NTask);

  MPI_Allgather(&nchanged, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);

  for(int task = 0; task < NTask; task++)
    bytecounts[task] = recvcounts[task] * sizeof(struct DomainNODE);

  byteoffset[0] = 0;
  for(int task = 1; task < NTask; task++)
    byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];

  struct DomainNODE *loc_DomainMoment =
      (struct DomainNODE *)mymalloc("loc_DomainMoment", recvcounts[ThisTask] * sizeof(struct DomainNODE));

  for(int i = 0; i < nchanged; i++)
    {
      int no                   = nodelist[i];
      loc_DomainMoment[i].node = no;

      for(int j = 0; j < 3; j++)
        {
          loc_DomainMoment[i].vertex_vmin[j] = Ngb_Nodes[no].vertex_vmin[j];
          loc_DomainMoment[i].vertex_vmax[j] = Ngb_Nodes[no].vertex_vmax[j];
#ifdef TREE_BASED_TIMESTEPS
          loc_DomainMoment[i].vmin[j] = ExtNgb_Nodes[no].vmin[j];
          loc_DomainMoment[i].vmax[j] = ExtNgb_Nodes[no].vmax[j];
#endif /* #ifdef TREE_BASED_TIMESTEPS */
        }
    }

  int tot_nchanged = 0;
  for(int task = 0; task < NTask; task++)
    tot_nchanged += recvcounts[task];

  struct DomainNODE *tot_DomainMoment = (struct DomainNODE *)mymalloc("tot_DomainMoment", tot_nchanged * sizeof(struct DomainNODE));

  MPI_Allgatherv(loc_DomainMoment, bytecounts[ThisTask], MPI_BYTE, tot_DomainMoment, bytecounts, byteoffset, MPI_BYTE, MPI_COMM_WORLD);

  for(int i = 0; i < tot_nchanged; i++)
    {
      int no = tot_DomainMoment[i].node;

      if(Ngb_Nodes[no].Ti_Current != All.Ti_Current)
        drift_node(&Ngb_Nodes[no], All.Ti_Current);

      for(int j = 0; j < 3; j++)
        {
          Ngb_Nodes[no].vertex_vmin[j] = tot_DomainMoment[i].vertex_vmin[j];
          Ngb_Nodes[no].vertex_vmax[j] = tot_DomainMoment[i].vertex_vmax[j];
#ifdef TREE_BASED_TIMESTEPS
          ExtNgb_Nodes[no].vmin[j] = tot_DomainMoment[i].vmin[j];
          ExtNgb_Nodes[no].vmax[j] = tot_DomainMoment[i].vmax[j];
#endif /* #ifdef TREE_BASED_TIMESTEPS */
        }

      no = Ngb_Nodes[no].father;

      while(no >= 0)
        {
          if(Ngb_Nodes[no].Ti_Current != All.Ti_Current)
            drift_node(&Ngb_Nodes[no], All.Ti_Current);

          int flag_changed = 0;

          for(int j = 0; j < 3; j++)
            {
              if(Ngb_Nodes[no].vertex_vmin[j] > tot_DomainMoment[i].vertex_vmin[j])
                {
                  Ngb_Nodes[no].vertex_vmin[j] = tot_DomainMoment[i].vertex_vmin[j];
                  flag_changed                 = 1;
                }

              if(Ngb_Nodes[no].vertex_vmax[j] < tot_DomainMoment[i].vertex_vmax[j])
                {
                  Ngb_Nodes[no].vertex_vmax[j] = tot_DomainMoment[i].vertex_vmax[j];
                  flag_changed                 = 1;
                }
#ifdef TREE_BASED_TIMESTEPS
              if(ExtNgb_Nodes[no].vmin[j] > tot_DomainMoment[i].vmin[j])
                {
                  ExtNgb_Nodes[no].vmin[j] = tot_DomainMoment[i].vmin[j];
                  flag_changed             = 1;
                }

              if(ExtNgb_Nodes[no].vmax[j] < tot_DomainMoment[i].vmax[j])
                {
                  ExtNgb_Nodes[no].vmax[j] = tot_DomainMoment[i].vmax[j];
                  flag_changed             = 1;
                }
#endif /* #ifdef TREE_BASED_TIMESTEPS */
            }

          if(flag_changed == 0)
            break;

          no = Ngb_Nodes[no].father;
        }
    }

  myfree(tot_DomainMoment);
  myfree(loc_DomainMoment);
  myfree(byteoffset);
  myfree(bytecounts);
  myfree(recvcounts);
}

/*! \brief Updates min and max position in ngb nodes.
 *
 *  Inverse tree walk.
 *
 *  \param[in] i Index of particle.
 *  \param[in, out] nchanged Number of changed top level nodes.
 *  \param[out] nodelist Top level nodes that were changed.
 *
 *  \return void
 */
void ngb_update_rangebounds(int i, int *nchanged, int *nodelist)
{
  int no = Ngb_Father[i];

  while(no >= 0)
    {
      if(Ngb_Nodes[no].Ti_Current != All.Ti_Current)
        drift_node(&Ngb_Nodes[no], All.Ti_Current);

      int flag_changed = 0;

      for(int j = 0; j < 3; j++)
        {
          if(Ngb_Nodes[no].u.d.range_min[j] > P[i].Pos[j])
            {
              Ngb_Nodes[no].u.d.range_min[j] = P[i].Pos[j];
              flag_changed                   = 1;
            }

          if(Ngb_Nodes[no].u.d.range_max[j] < P[i].Pos[j])
            {
              Ngb_Nodes[no].u.d.range_max[j] = P[i].Pos[j];
              flag_changed                   = 1;
            }
        }

      if(flag_changed == 0)
        break;

      if(no < Ngb_FirstNonTopLevelNode) /* top-level tree-node reached */
        {
          if(Ngb_Marker[no] != Ngb_MarkerValue)
            {
              Ngb_Marker[no]      = Ngb_MarkerValue;
              nodelist[*nchanged] = no;
              *nchanged           = *nchanged + 1;
            }
          break;
        }

      no = Ngb_Nodes[no].father;
    }
}

/*! \brief Finalizes position bounds update.
 *
 *  Exchanges changed information in top level nodes to all tasks.
 *
 *  \param[in] nchanged Number of changed top level nodes.
 *  \param[in] nodelist List of changed top level nodes.
 *
 *  \return void
 */
void ngb_finish_rangebounds_update(int nchanged, int *nodelist)
{
  struct DomainNODE
  {
    int node;
    MyNgbTreeFloat range_min[3];
    MyNgbTreeFloat range_max[3];
  };

  /* share the pseudo-particle data accross CPUs */
  int *recvcounts = (int *)mymalloc("recvcounts", sizeof(int) * NTask);
  int *bytecounts = (int *)mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *)mymalloc("byteoffset", sizeof(int) * NTask);

  MPI_Allgather(&nchanged, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);

  for(int task = 0; task < NTask; task++)
    bytecounts[task] = recvcounts[task] * sizeof(struct DomainNODE);

  byteoffset[0] = 0;
  for(int task = 1; task < NTask; task++)
    byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];

  struct DomainNODE *loc_DomainMoment =
      (struct DomainNODE *)mymalloc("loc_DomainMoment", recvcounts[ThisTask] * sizeof(struct DomainNODE));

  for(int i = 0; i < nchanged; i++)
    {
      int no                   = nodelist[i];
      loc_DomainMoment[i].node = no;

      for(int j = 0; j < 3; j++)
        {
          loc_DomainMoment[i].range_min[j] = Ngb_Nodes[no].u.d.range_min[j];
          loc_DomainMoment[i].range_max[j] = Ngb_Nodes[no].u.d.range_max[j];
        }
    }

  int tot_nchanged = 0;
  for(int task = 0; task < NTask; task++)
    tot_nchanged += recvcounts[task];

  struct DomainNODE *tot_DomainMoment = (struct DomainNODE *)mymalloc("tot_DomainMoment", tot_nchanged * sizeof(struct DomainNODE));

  MPI_Allgatherv(loc_DomainMoment, bytecounts[ThisTask], MPI_BYTE, tot_DomainMoment, bytecounts, byteoffset, MPI_BYTE, MPI_COMM_WORLD);

  for(int i = 0; i < tot_nchanged; i++)
    {
      int no = tot_DomainMoment[i].node;

      if(Ngb_Nodes[no].Ti_Current != All.Ti_Current)
        drift_node(&Ngb_Nodes[no], All.Ti_Current);

      for(int j = 0; j < 3; j++)
        {
          Ngb_Nodes[no].u.d.range_min[j] = tot_DomainMoment[i].range_min[j];
          Ngb_Nodes[no].u.d.range_max[j] = tot_DomainMoment[i].range_max[j];
        }

      no = Ngb_Nodes[no].father;

      while(no >= 0)
        {
          if(Ngb_Nodes[no].Ti_Current != All.Ti_Current)
            drift_node(&Ngb_Nodes[no], All.Ti_Current);

          int flag_changed = 0;

          for(int j = 0; j < 3; j++)
            {
              if(Ngb_Nodes[no].u.d.range_min[j] > tot_DomainMoment[i].range_min[j])
                {
                  Ngb_Nodes[no].u.d.range_min[j] = tot_DomainMoment[i].range_min[j];
                  flag_changed                   = 1;
                }

              if(Ngb_Nodes[no].u.d.range_max[j] < tot_DomainMoment[i].range_max[j])
                {
                  Ngb_Nodes[no].u.d.range_max[j] = tot_DomainMoment[i].range_max[j];
                  flag_changed                   = 1;
                }
            }

          if(flag_changed == 0)
            break;

          no = Ngb_Nodes[no].father;
        }
    }

  myfree(tot_DomainMoment);
  myfree(loc_DomainMoment);
  myfree(byteoffset);
  myfree(bytecounts);
  myfree(recvcounts);
}

/*! \brief Adjust ngb-tree structures due to a change in number of gas cells.
 *
 *  \param[in] delta_NgbMaxPart Difference in number of cells.
 *
 *  \return void
 */
void ngb_treemodifylength(int delta_NgbMaxPart)
{
  mpi_printf("ALLOCATE: Need to adjust NgbTree because Ngb_MaxPart needs to grow by %d\n", delta_NgbMaxPart);

  for(int i = 0; i < Ngb_MaxPart + NTopleaves; i++) /* check for particles and pseudo particles */
    if(Ngb_Nextnode[i] >= Ngb_MaxPart)              /* internal node or pseudo particle */
      Ngb_Nextnode[i] += delta_NgbMaxPart;

  for(int i = 0; i < Ngb_MaxPart; i++)
    if(Ngb_Father[i] >= Ngb_MaxPart) /* internal node or pseudo particle */
      Ngb_Father[i] += delta_NgbMaxPart;

  for(int i = 0; i < Ngb_MaxNodes; i++)
    {
      if(Ngb_Nodes[i + Ngb_MaxPart].u.d.nextnode >= Ngb_MaxPart) /* internal node or pseudo particle */
        Ngb_Nodes[i + Ngb_MaxPart].u.d.nextnode += delta_NgbMaxPart;

      if(Ngb_Nodes[i + Ngb_MaxPart].u.d.sibling >= Ngb_MaxPart) /* internal node or pseudo particle */
        Ngb_Nodes[i + Ngb_MaxPart].u.d.sibling += delta_NgbMaxPart;

      if(Ngb_Nodes[i + Ngb_MaxPart].father >= Ngb_MaxPart)
        Ngb_Nodes[i + Ngb_MaxPart].father += delta_NgbMaxPart;
    }

  for(int i = 0; i < NTopleaves; i++)
    Ngb_DomainNodeIndex[i] += delta_NgbMaxPart;

  Ngb_Nextnode = (int *)myrealloc_movable(Ngb_Nextnode, (Ngb_MaxPart + delta_NgbMaxPart + NTopleaves) * sizeof(int));

  memmove(&Ngb_Nextnode[Ngb_MaxPart + delta_NgbMaxPart], &Ngb_Nextnode[Ngb_MaxPart], NTopleaves * sizeof(int));

  Ngb_MaxPart += delta_NgbMaxPart;

  Ngb_FirstNonTopLevelNode += delta_NgbMaxPart;

  Ngb_Nodes -= delta_NgbMaxPart;

#ifdef TREE_BASED_TIMESTEPS
  ExtNgb_Nodes -= delta_NgbMaxPart;
#endif /* #ifdef TREE_BASED_TIMESTEPS */

  Ngb_Father = (int *)myrealloc_movable(Ngb_Father, Ngb_MaxPart * sizeof(int));

  Ngb_Marker = (int *)myrealloc_movable(Ngb_Marker, (Ngb_MaxNodes + Ngb_MaxPart) * sizeof(int));
  memmove(Ngb_Marker + Ngb_MaxPart, Ngb_Marker + Ngb_MaxPart - delta_NgbMaxPart, Ngb_MaxNodes * sizeof(int));
  memset(Ngb_Marker + Ngb_MaxPart - delta_NgbMaxPart, -1, delta_NgbMaxPart * sizeof(int));
}

/*! \brief Allocates arrays for neighbor tree.
 *
 *  \return void
 */
void ngb_treeallocate(void)
{
  if(Ngb_MaxPart == 0)
    {
      Ngb_MaxPart  = All.MaxPartSph;
      Ngb_MaxNodes = (int)(All.NgbTreeAllocFactor * (All.MaxPartSph + BASENUMBER)) + NTopnodes;
    }

  if(All.TotNumGas == 0)
    return;

  if(Ngb_Nodes)
    terminate("already allocated");

  Ngb_DomainNodeIndex = (int *)mymalloc_movable(&Ngb_DomainNodeIndex, "Ngb_DomainNodeIndex", NTopleaves * sizeof(int));

  Ngb_Nodes = (struct NgbNODE *)mymalloc_movable(&Ngb_Nodes, "Ngb_Nodes", (Ngb_MaxNodes + 1) * sizeof(struct NgbNODE));
  Ngb_Nodes -= Ngb_MaxPart;

#ifdef TREE_BASED_TIMESTEPS
  ExtNgb_Nodes = (struct ExtNgbNODE *)mymalloc_movable(&ExtNgb_Nodes, "ExtNgb_Nodes", (Ngb_MaxNodes + 1) * sizeof(struct ExtNgbNODE));
  ExtNgb_Nodes -= Ngb_MaxPart;
#endif /* #ifdef TREE_BASED_TIMESTEPS */
  Ngb_Nextnode = (int *)mymalloc_movable(&Ngb_Nextnode, "Ngb_Nextnode", (Ngb_MaxPart + NTopleaves) * sizeof(int));
  Ngb_Father   = (int *)mymalloc_movable(&Ngb_Father, "Ngb_Father", Ngb_MaxPart * sizeof(int));

  Ngb_Marker = (int *)mymalloc_movable(&Ngb_Marker, "Ngb_Marker", (Ngb_MaxNodes + Ngb_MaxPart) * sizeof(int));
}

/*! \brief This function frees the memory allocated for the neighbor tree.
 *
 *  \return void
 */
void ngb_treefree(void)
{
  if(All.TotNumGas == 0)
    return;

  if(Ngb_Nodes)
    {
      myfree_movable(Ngb_Marker);
      myfree_movable(Ngb_Father);
      myfree_movable(Ngb_Nextnode);
#ifdef TREE_BASED_TIMESTEPS
      myfree_movable(ExtNgb_Nodes + Ngb_MaxPart);
      ExtNgb_Nodes = NULL;
#endif /* #ifdef TREE_BASED_TIMESTEPS */
      myfree_movable(Ngb_Nodes + Ngb_MaxPart);
      myfree_movable(Ngb_DomainNodeIndex);

      Ngb_Marker          = NULL;
      Ngb_Father          = NULL;
      Ngb_Nodes           = NULL;
      Ngb_DomainNodeIndex = NULL;
      Ngb_Nextnode        = NULL;
      Ngb_MaxPart         = 0;
      Ngb_MaxNodes        = 0;
    }
  else
    terminate("trying to free the tree even though it's not allocated");
}
