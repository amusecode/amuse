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
 * \file        src/subfind/subfind_so_potegy.c
 * \date        05/2018
 * \brief       Calculates the the potential energy.
 * \details     contains functions:
 *                static void subfind_so_potegy_loctree_findExtent(int npart,
 *                  int start)
 *                static int subfind_so_potegy_loctree_treebuild(int npart,
 *                  int start)
 *                static void subfind_so_potegy_loctree_update_node_recursive(
 *                  int no, int sib, int father)
 *                double subfind_so_potegy_loctree_treeevaluate_potential(int
 *                  target)
 *                static size_t subfind_so_potegy_loctree_treeallocate(int
 *                  maxnodes, int maxpart)
 *                static void subfind_so_potegy_loctree_treefree(void)
 *                static int subfind_compare_Paux_LocGrIndex(const void *a,
 *                  const void *b)
 *                double subfind_so_potegy(double *egypot)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 14.05.2018 Prepared file for public release -- Rainer Weinberger
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

#if defined(SUBFIND) && defined(SUBFIND_EXTENDED_PROPERTIES)

#include "../fof/fof.h"
#include "subfind.h"

static double RootLen, RootFac, RootBigFac, RootInverseLen, RootCenter[3], RootCorner[3];
static int LocMaxPart;
static int MaxNodes, last;
static int *LocNextNode;
static unsigned long long *LocTree_IntPos_list;
static struct paux_data *LocPaux;

static void subfind_so_potegy_loctree_update_node_recursive(int no, int sib, int father);

/*! \brief Node structure for local tree.
 */
static struct LocNODE
{
  union
  {
    int suns[8]; /*!< temporary pointers to daughter nodes */
    struct
    {
      MyDouble s[3]; /*!< center of mass of node */
      MyDouble mass; /*!< mass of node */
      unsigned char maxsofttype;
#if defined(MULTIPLE_NODE_SOFTENING) && defined(ADAPTIVE_HYDRO_SOFTENING)
      unsigned char maxhydrosofttype;
      unsigned char minhydrosofttype;
#endif              /* #if defined(MULTIPLE_NODE_SOFTENING) && defined(ADAPTIVE_HYDRO_SOFTENING) */
      int sibling;  /*!< this gives the next node in the walk in case the current node can be used */
      int nextnode; /*!< this gives the next node in case the current node needs to be opened */
    } d;
  } u;

  MyDouble center[3]; /*!< geometrical center of node */
  MyFloat len;        /*!< sidelength of treenode */

#ifdef MULTIPLE_NODE_SOFTENING
  MyDouble mass_per_type[NSOFTTYPES];
#endif
} * LocNodes_base, /*!< points to the actual memory allocted for the nodes */
    *LocNodes;     /*!< this is a pointer used to access the nodes which is shifted such that Nodes[LocMaxPart]
                      gives the first allocated node */

/*! \brief Finds spatial extent of local particles.
 *
 *  Sets global 'Root*' variables that determine root node properties.
 *
 *  \param[in] npart Number of particles.
 *  \param[in] start Start index.
 *
 *  \return void
 */
static void subfind_so_potegy_loctree_findExtent(int npart, int start)
{
  double len, xmin[3], xmax[3];

  /* determine extension */
  for(int i = 0; i < 3; i++)
    {
      xmin[i] = MAX_REAL_NUMBER;
      xmax[i] = -MAX_REAL_NUMBER;
    }

  for(int k = 0; k < npart; k++)
    {
      int i = start + k;

      for(int j = 0; j < 3; j++)
        {
          if(xmin[j] > LocPaux[i].Pos[j])
            xmin[j] = LocPaux[i].Pos[j];

          if(xmax[j] < LocPaux[i].Pos[j])
            xmax[j] = LocPaux[i].Pos[j];
        }
    }

  len = 0;
  for(int j = 0; j < 3; j++)
    if(xmax[j] - xmin[j] > len)
      len = xmax[j] - xmin[j];

  len *= 1.001;

  RootLen        = len;
  RootInverseLen = 1.0 / RootLen;
  RootFac        = 1.0 / len * (((peanokey)1) << (BITS_PER_DIMENSION));
  RootBigFac     = (RootLen / (((long long)1) << 52));

  for(int j = 0; j < 3; j++)
    {
      RootCenter[j] = 0.5 * (xmin[j] + xmax[j]);
      RootCorner[j] = 0.5 * (xmin[j] + xmax[j]) - 0.5 * len;
    }
}

/*! \brief Builds local tree.
 *
 *  \param[in] npart Number of particles.
 *  \param[in] start Start index.
 *
 *  \return Number of nodes in tree.
 */
static int subfind_so_potegy_loctree_treebuild(int npart, int start)
{
  int subnode = 0, parent = -1, numnodes;
  int nfree, th, nn;
  struct LocNODE *nfreep;

  /* select first node */
  nfree  = LocMaxPart;
  nfreep = &LocNodes[nfree];

  /* create an empty  root node  */
  nfreep->len = (MyFloat)RootLen;
  for(int i = 0; i < 3; i++)
    nfreep->center[i] = (MyFloat)RootCenter[i];

  for(int i = 0; i < 8; i++)
    nfreep->u.suns[i] = -1;

  numnodes = 1;
  nfreep++;
  nfree++;

  /* insert all particles */

  LocTree_IntPos_list =
      (unsigned long long *)mymalloc_movable(&LocTree_IntPos_list, "LocTree_IntPos_list", 3 * LocMaxPart * sizeof(unsigned long long));

  for(int k = 0; k < npart; k++)
    {
      int i = start + k;

      MyDouble *posp;

      posp = &LocPaux[i].Pos[0];

      unsigned long long xxb      = force_double_to_int(((posp[0] - RootCorner[0]) * RootInverseLen) + 1.0);
      unsigned long long yyb      = force_double_to_int(((posp[1] - RootCorner[1]) * RootInverseLen) + 1.0);
      unsigned long long zzb      = force_double_to_int(((posp[2] - RootCorner[2]) * RootInverseLen) + 1.0);
      unsigned long long mask     = ((unsigned long long)1) << (52 - 1);
      unsigned char shiftx        = (52 - 1);
      unsigned char shifty        = (52 - 2);
      unsigned char shiftz        = (52 - 3);
      signed long long centermask = (0xFFF0000000000000llu);
      unsigned char levels        = 0;

      unsigned long long *intposp = &LocTree_IntPos_list[3 * i];

      *intposp++ = xxb;
      *intposp++ = yyb;
      *intposp++ = zzb;

      th = LocMaxPart;

      while(1)
        {
          if(th >= LocMaxPart) /* we are dealing with an internal node */
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
                   *      DomainLen/2^MAX_TREE_LEVEL  < gravitational softening length
                   */
                  for(int j = 0; j < 8; j++)
                    {
                      if(LocNodes[th].u.suns[subnode] < 0)
                        break;

                      subnode++;
                      if(subnode >= 8)
                        subnode = 7;
                    }
                }

              nn = LocNodes[th].u.suns[subnode];

              if(nn >= 0) /* ok, something is in the daughter slot already, need to continue */
                {
                  parent = th; /* note: subnode can still be used in the next step of the walk */
                  th     = nn;
                }
              else
                {
                  /* here we have found an empty slot where we can
                   * attach the new particle as a leaf
                   */
                  LocNodes[th].u.suns[subnode] = i;
                  break; /* done for this particle */
                }
            }
          else
            {
              /* we try to insert into a leaf with a single particle
               * need to generate a new internal node at this point
               */
              LocNodes[parent].u.suns[subnode] = nfree;

              /* the other is: */
              double len = ((double)(mask << 1)) * RootBigFac;
              double cx  = ((double)((xxb & centermask) | mask)) * RootBigFac + RootCorner[0];
              double cy  = ((double)((yyb & centermask) | mask)) * RootBigFac + RootCorner[1];
              double cz  = ((double)((zzb & centermask) | mask)) * RootBigFac + RootCorner[2];

              nfreep->len       = len;
              nfreep->center[0] = cx;
              nfreep->center[1] = cy;
              nfreep->center[2] = cz;

              nfreep->u.suns[0] = -1;
              nfreep->u.suns[1] = -1;
              nfreep->u.suns[2] = -1;
              nfreep->u.suns[3] = -1;
              nfreep->u.suns[4] = -1;
              nfreep->u.suns[5] = -1;
              nfreep->u.suns[6] = -1;
              nfreep->u.suns[7] = -1;

              unsigned long long *intppos = &LocTree_IntPos_list[3 * th];

              subnode = (((unsigned char)((intppos[0] & mask) >> shiftx)) | ((unsigned char)((intppos[1] & mask) >> shifty)) |
                         ((unsigned char)((intppos[2] & mask) >> shiftz)));

              nfreep->u.suns[subnode] = th;

              th = nfree; /* resume trying to insert the new particle at
                             the newly created internal node */

              numnodes++;
              nfree++;
              nfreep++;

              if(numnodes >= MaxNodes)
                {
                  MaxNodes *= 1.2;

                  LocNodes_base = (struct LocNODE *)myrealloc_movable(LocNodes_base, (MaxNodes + 1) * sizeof(struct LocNODE));
                  LocNodes      = LocNodes_base - LocMaxPart;
                  nfreep        = &LocNodes[nfree];

                  if(numnodes > MaxNodes)
                    {
                      char buf[1000];

                      sprintf(buf, "maximum number %d of tree-nodes reached., for particle %d  %g %g %g", MaxNodes, i,
                              LocPaux[i].Pos[0], LocPaux[i].Pos[1], LocPaux[i].Pos[2]);
                      terminate(buf);
                    }
                }
            }
        }
    }

  myfree(LocTree_IntPos_list);

  /* now compute the multipole moments recursively */
  last = -1;
  subfind_so_potegy_loctree_update_node_recursive(LocMaxPart, -1, -1);

  if(last >= LocMaxPart)
    LocNodes[last].u.d.nextnode = -1;
  else
    LocNextNode[last] = -1;

  return numnodes;
}

/*! \brief Walk the tree and update node data recursively.
 *
 *  This routine computes the multipole moments for a given internal node and
 *  all its subnodes using a recursive computation. Note that this switches
 *  the information stored in LocNodes[no].u from suns to d!
 *
 *
 *  \param[in] no Node index.
 *  \param[in] sib Sibling index.
 *  \param[in] father Parent index.
 *
 *  \return void
 */
static void subfind_so_potegy_loctree_update_node_recursive(int no, int sib, int father)
{
  int j, jj, p, pp = 0, nextsib, suns[8];
  unsigned char maxsofttype;
#ifdef MULTIPLE_NODE_SOFTENING
  double mass_per_type[NSOFTTYPES];
#ifdef ADAPTIVE_HYDRO_SOFTENING
  unsigned char maxhydrosofttype;
  unsigned char minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

  double mass;
  double s[3];

  if(no >= LocMaxPart)
    {
      for(j = 0; j < 8; j++)
        suns[j] = LocNodes[no].u.suns[j]; /* this "backup" is necessary because the nextnode entry will
                                             overwrite one element (union!) */
      if(last >= 0)
        {
          if(last >= LocMaxPart)
            LocNodes[last].u.d.nextnode = no;
          else
            LocNextNode[last] = no;
        }

      last = no;

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

              subfind_so_potegy_loctree_update_node_recursive(p, nextsib, no);

              if(p >= LocMaxPart) /* an internal node  */
                {
                  mass += LocNodes[p].u.d.mass; /* we assume a fixed particle mass */
                  s[0] += LocNodes[p].u.d.mass * LocNodes[p].u.d.s[0];
                  s[1] += LocNodes[p].u.d.mass * LocNodes[p].u.d.s[1];
                  s[2] += LocNodes[p].u.d.mass * LocNodes[p].u.d.s[2];

                  if(All.ForceSoftening[maxsofttype] < All.ForceSoftening[LocNodes[p].u.d.maxsofttype])
                    maxsofttype = LocNodes[p].u.d.maxsofttype;

#ifdef MULTIPLE_NODE_SOFTENING
                  int k;
                  for(k = 0; k < NSOFTTYPES; k++)
                    mass_per_type[k] += LocNodes[p].mass_per_type[k];

#ifdef ADAPTIVE_HYDRO_SOFTENING
                  if(maxhydrosofttype < LocNodes[p].u.d.maxhydrosofttype)
                    maxhydrosofttype = LocNodes[p].u.d.maxhydrosofttype;
                  if(minhydrosofttype > LocNodes[p].u.d.minhydrosofttype)
                    minhydrosofttype = LocNodes[p].u.d.minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
                }
              else /* a particle */
                {
                  mass += LocPaux[p].Mass;

                  s[0] += LocPaux[p].Mass * LocPaux[p].Pos[0];
                  s[1] += LocPaux[p].Mass * LocPaux[p].Pos[1];
                  s[2] += LocPaux[p].Mass * LocPaux[p].Pos[2];

                  if(All.ForceSoftening[maxsofttype] < All.ForceSoftening[LocPaux[p].SofteningType])
                    maxsofttype = LocPaux[p].SofteningType;
#ifdef MULTIPLE_NODE_SOFTENING
#ifdef ADAPTIVE_HYDRO_SOFTENING
                  mass_per_type[LocPaux[p].Type == 0 ? 0 : LocPaux[p].SofteningType] += LocPaux[p].Mass;

                  if(LocPaux[p].Type == 0)
                    {
                      if(maxhydrosofttype < LocPaux[p].SofteningType)
                        maxhydrosofttype = LocPaux[p].SofteningType;
                      if(minhydrosofttype > LocPaux[p].SofteningType)
                        minhydrosofttype = LocPaux[p].SofteningType;
                    }
#else  /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
                  mass_per_type[LocPaux[p].SofteningType] += LocPaux[p].Mass;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING #else */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
                }
            }
        }

      if(mass > 0)
        {
          s[0] /= mass;
          s[1] /= mass;
          s[2] /= mass;
        }
      else
        {
          s[0] = LocNodes[no].center[0];
          s[1] = LocNodes[no].center[1];
          s[2] = LocNodes[no].center[2];
        }

      LocNodes[no].u.d.s[0]        = (MyFloat)s[0];
      LocNodes[no].u.d.s[1]        = (MyFloat)s[1];
      LocNodes[no].u.d.s[2]        = (MyFloat)s[2];
      LocNodes[no].u.d.mass        = (MyFloat)mass;
      LocNodes[no].u.d.maxsofttype = maxsofttype;
#ifdef MULTIPLE_NODE_SOFTENING
      int k;
      for(k = 0; k < NSOFTTYPES; k++)
        LocNodes[no].mass_per_type[k] = mass_per_type[k];

#ifdef ADAPTIVE_HYDRO_SOFTENING
      LocNodes[no].u.d.maxhydrosofttype = maxhydrosofttype;
      LocNodes[no].u.d.minhydrosofttype = minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

      LocNodes[no].u.d.sibling = sib;
    }
  else /* single particle or pseudo particle */
    {
      if(last >= 0)
        {
          if(last >= LocMaxPart)
            LocNodes[last].u.d.nextnode = no;
          else
            LocNextNode[last] = no;
        }

      last = no;
    }
}

/*! \brief Calculates the gravitational potential energy of single particle.
 *
 *  \pararm[in] target Target particle index (in LocPaux).
 *
 *  \return Gravitational potential.
 */
double subfind_so_potegy_loctree_treeevaluate_potential(int target)
{
  struct LocNODE *nop = 0;
  int no;
  double r2, dx, dy, dz, mass, r, u, h_i, h_j, hmax, h_inv, wp;
  double pot, pos_x, pos_y, pos_z, xtmp, ytmp, ztmp;

  pos_x = LocPaux[target].Pos[0];
  pos_y = LocPaux[target].Pos[1];
  pos_z = LocPaux[target].Pos[2];

  h_i = All.ForceSoftening[LocPaux[target].SofteningType];

  pot = 0;

  no = LocMaxPart;

  while(no >= 0)
    {
#ifdef MULTIPLE_NODE_SOFTENING
      int indi_flag1 = -1, indi_flag2 = 0;
#endif                    /* #ifdef MULTIPLE_NODE_SOFTENING */
      if(no < LocMaxPart) /* single particle */
        {
          dx = GRAVITY_NEAREST_X(LocPaux[no].Pos[0] - pos_x);
          dy = GRAVITY_NEAREST_Y(LocPaux[no].Pos[1] - pos_y);
          dz = GRAVITY_NEAREST_Z(LocPaux[no].Pos[2] - pos_z);

          r2 = dx * dx + dy * dy + dz * dz;

          mass = LocPaux[no].Mass;

          h_j = All.ForceSoftening[LocPaux[no].SofteningType];

          if(h_j > h_i)
            hmax = h_j;
          else
            hmax = h_i;

          no = LocNextNode[no];
        }
      else
        {
          nop  = &LocNodes[no];
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
                if(LocNodes[no].mass_per_type[0] > 0)
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

          no = nop->u.d.sibling; /* node can be used */
        }

      r = sqrt(r2);
#ifdef MULTIPLE_NODE_SOFTENING
      int type;
      for(type = indi_flag1; type < indi_flag2; type++)
        {
          if(type >= 0)
            {
              mass = nop->mass_per_type[type];

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
                pot -= mass / r;
              else
                {
                  h_inv = 1.0 / hmax;

                  u = r * h_inv;

                  if(u < 0.5)
                    wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
                  else
                    wp = -3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u)));

                  pot += mass * h_inv * wp;
#ifdef MULTIPLE_NODE_SOFTENING
                }
            }
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
        }
    }

  return pot;
}

/*! \brief Allocates memory used for storage of the tree and auxiliary arrays
 *         for tree-walk and link-lists.
 *
 *  \param[in] maxnodes Maximum number of nodes.
 *  \param[in] maxpart Maximum number of particles.
 *
 *  \return Number of allocated bytes.
 */
static size_t subfind_so_potegy_loctree_treeallocate(int maxnodes, int maxpart)
{
  size_t bytes, allbytes = 0;

  if(LocNextNode)
    terminate("loctree already allocated");

  MaxNodes   = maxnodes;
  LocMaxPart = maxpart;

  LocNextNode = (int *)mymalloc("LocNextNode", bytes = maxpart * sizeof(int));
  allbytes += bytes;

  R2list = (r2type *)mymalloc("R2list", bytes = maxpart * sizeof(r2type));
  allbytes += bytes;

  LocNodes_base = (struct LocNODE *)mymalloc_movable(&LocNodes_base, "LocNodes_base", bytes = (MaxNodes + 1) * sizeof(struct LocNODE));
  LocNodes      = LocNodes_base - LocMaxPart;
  allbytes += bytes;

  return allbytes;
}

/*! \brief Frees the allocated memory.
 *
 *  \return void
 */
static void subfind_so_potegy_loctree_treefree(void)
{
  myfree(LocNodes_base);
  myfree(R2list);
  myfree(LocNextNode);

  LocNextNode   = NULL;
  R2list        = NULL;
  LocNodes_base = NULL;
}

/*! \brief Comparison function for paux_data objects.
 *
 *  Compares field LocGrIndex.
 *
 *  \param[in] a First object to be compared.
 *  \param[in] b Second object to be compared.
 *
 *  \return (-1,0,1); -1 if a < b.
 */
static int subfind_compare_Paux_LocGrIndex(const void *a, const void *b)
{
  if(((struct paux_data *)a)->LocGrIndex < ((struct paux_data *)b)->LocGrIndex)
    return -1;

  if(((struct paux_data *)a)->LocGrIndex > ((struct paux_data *)b)->LocGrIndex)
    return +1;

  return 0;
}

/*! \brief Calculates potential energy of spherical overdensity groups.
 *
 *  \param[out] egypot Array with potential energies in each group.
 *
 *  \return Time this routine took.
 */
double subfind_so_potegy(double *egypot)
{
  double t0 = second();
  mpi_printf("SUBFIND: Starting SO potential energy computation\n");

  size_t *count_send  = (size_t *)mymalloc_movable(&count_send, "count_send", NTask * sizeof(size_t));
  size_t *offset_send = (size_t *)mymalloc_movable(&offset_send, "offset_send", NTask * sizeof(size_t));
  size_t *count_recv  = (size_t *)mymalloc_movable(&count_recv, "count_recv", NTask * sizeof(size_t));
  size_t *offset_recv = (size_t *)mymalloc_movable(&offset_recv, "offset_recv", NTask * sizeof(size_t));

  for(int i = 0; i < NTask; i++)
    count_send[i] = 0;

  for(int i = 0; i < NumPaux; i++)
    count_send[Paux[i].TaskOfGr]++;

  MPI_Alltoall(count_send, sizeof(size_t), MPI_BYTE, count_recv, sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);

  offset_send[0] = offset_recv[0] = 0;

  for(int i = 1; i < NTask; i++)
    {
      offset_send[i] = offset_send[i - 1] + count_send[i - 1];
      offset_recv[i] = offset_recv[i - 1] + count_recv[i - 1];
    }

  struct paux_data *PauxTmp = (struct paux_data *)mymalloc_movable(&PauxTmp, "PauxTmp", NumPaux * sizeof(struct paux_data));

  for(int i = 0; i < NTask; i++)
    count_send[i] = 0;

  for(int i = 0; i < NumPaux; i++)
    {
      int task     = Paux[i].TaskOfGr;
      int loc      = offset_send[task] + count_send[task]++;
      PauxTmp[loc] = Paux[i];
    }

  int NumPauxRecv = 0;

  for(int i = 0; i < NTask; i++)
    NumPauxRecv += count_recv[i];

  LocPaux = (struct paux_data *)mymalloc_movable(&LocPaux, "LocPaux", NumPauxRecv * sizeof(struct paux_data));

  myMPI_Alltoallv(PauxTmp, count_send, offset_send, LocPaux, count_recv, offset_recv, sizeof(struct paux_data), 1, MPI_COMM_WORLD);

  myfree_movable(PauxTmp);

  qsort(LocPaux, NumPauxRecv, sizeof(struct paux_data), subfind_compare_Paux_LocGrIndex);

  int *group_len = (int *)mymalloc("group_len", Ngroups * sizeof(int));
  int *group_off = (int *)mymalloc("group_off", Ngroups * sizeof(int));

  for(int i = 0; i < Ngroups; i++)
    group_len[i] = 0;

  for(int i = 0; i < NumPauxRecv; i++)
    {
      int j = LocPaux[i].LocGrIndex;
      if(j < 0 || j >= Ngroups)
        terminate("j=%d Ngroups=%d", j, Ngroups);

      group_len[j]++;
    }

  group_off[0] = 0;

  for(int i = 1; i < Ngroups; i++)
    group_off[i] = group_off[i - 1] + group_len[i - 1];

  int MaxAllocPart = NumPart;
  // extend in case a single group holds more particles than NumPart
  for(int i = 0; i < Ngroups; i++)
    if(group_len[i] > MaxAllocPart)
      MaxAllocPart = group_len[i];

  subfind_so_potegy_loctree_treeallocate((int)(All.TreeAllocFactor * MaxAllocPart) + NTopnodes, MaxAllocPart);

  /* now do the actual potential calculation */
  for(int i = 0; i < Ngroups; i++)
    {
      subfind_so_potegy_loctree_findExtent(group_len[i], group_off[i]);
      subfind_so_potegy_loctree_treebuild(group_len[i], group_off[i]);

      egypot[i] = 0;

      for(int j = 0; j < group_len[i]; j++)
        {
          int target = group_off[i] + j;

          double pot = subfind_so_potegy_loctree_treeevaluate_potential(target);

          /* remove self-potential */
          pot += LocPaux[target].Mass / (All.ForceSoftening[LocPaux[target].SofteningType] / 2.8);

          pot *= All.G / All.cf_atime;

          egypot[i] += 0.5 * pot * LocPaux[target].Mass;
        }
    }

  subfind_so_potegy_loctree_treefree();

  myfree(group_off);
  myfree(group_len);

  myfree(LocPaux);

  myfree(offset_recv);
  myfree(count_recv);
  myfree(offset_send);
  myfree(count_send);

  double t1 = second();
  mpi_printf("SUBFIND: SO potential energy computation took %g sec\n", timediff(t0, t1));

  return timediff(t0, t1);
}

#endif /* #if defined(SUBFIND) && defined(SUBFIND_EXTENDED_PROPERTIES) */
