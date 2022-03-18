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
 * \file        src/subfind/subfind_loctree.c
 * \date        05/2018
 * \brief       Algorithms for local tree in subfind.
 * \details     contains functions:
 *                void subfind_loctree_findExtent(int npart, struct unbind_data *mp)
 *                void subfind_loctree_copyExtent(void)
 *                int subfind_loctree_treebuild(int npart, struct unbind_data **udp)
 *                void subfind_loctree_update_node_recursive(int no, int sib, int father)
 *                double subfind_loctree_treeevaluate_potential(int target)
 *                int subfind_locngb_compare_key(const void *a, const void *b)
 *                double subfind_locngb_treefind(MyDouble xyz[3], int desngb, double hguess)
 *                int subfind_locngb_treefind_variable(MyDouble searchcenter[3], double hguess)
 *                size_t subfind_loctree_treeallocate(int maxnodes, int maxpart)
 *                void subfind_loctree_treefree(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 14.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"
#include "../gravity/forcetree.h"
#include "subfind.h"

#ifdef SUBFIND
static double RootLen, RootFac, RootBigFac, RootInverseLen, RootCenter[3], RootCorner[3];
static int LocMaxPart;
static int MaxNodes, last;
static int *LocNextNode;

static unsigned long long *LocTree_IntPos_list;

/*! \brief Node structure for subfind tree.
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
#endif             /* #ifdef MULTIPLE_NODE_SOFTENING */
} * LocNodes_base, /*!< points to the actual memory allocted for the nodes */
    *LocNodes;     /*!< this is a pointer used to access the nodes which is shifted such that Nodes[LocMaxPart]
                      gives the first allocated node */

/*! \brief Calculates min/max coordinate of particles in unbind data.
 *
 *  \param[in] npart Number of local particles (in unbind_data).
 *  \param[in] mp Pointer to unbind data.
 *
 *  \return void
 */
void subfind_loctree_findExtent(int npart, struct unbind_data *mp)
{
  int i, j, k;
  double len, xmin[3], xmax[3];

  /* determine extension */
  for(i = 0; i < 3; i++)
    {
      xmin[i] = MAX_REAL_NUMBER;
      xmax[i] = -MAX_REAL_NUMBER;
    }

  for(k = 0; k < npart; k++)
    {
      if(mp)
        i = mp[k].index;
      else
        terminate("what?");

#ifdef CELL_CENTER_GRAVITY
      if(P[i].Type == 0)
        {
          for(j = 0; j < 3; j++)
            {
              if(xmin[j] > PS[i].Center[j])
                xmin[j] = PS[i].Center[j];

              if(xmax[j] < PS[i].Center[j])
                xmax[j] = PS[i].Center[j];
            }
        }
      else
#endif /* #ifdef CELL_CENTER_GRAVITY */
        {
          for(j = 0; j < 3; j++)
            {
              if(xmin[j] > P[i].Pos[j])
                xmin[j] = P[i].Pos[j];

              if(xmax[j] < P[i].Pos[j])
                xmax[j] = P[i].Pos[j];
            }
        }
    }

  len = 0;
  for(j = 0; j < 3; j++)
    if(xmax[j] - xmin[j] > len)
      len = xmax[j] - xmin[j];

  len *= 1.001;

  RootLen        = len;
  RootInverseLen = 1.0 / RootLen;
  RootFac        = 1.0 / len * (((peanokey)1) << (BITS_PER_DIMENSION));
  RootBigFac     = (RootLen / (((long long)1) << 52));

  for(j = 0; j < 3; j++)
    {
      RootCenter[j] = 0.5 * (xmin[j] + xmax[j]);
      RootCorner[j] = 0.5 * (xmin[j] + xmax[j]) - 0.5 * len;
    }
}

/*! \brief Copy extent information from SubDomain to Root.
 *
 *  This is called from the collective subfind code.
 *
 *  \return void
 */
void subfind_loctree_copyExtent(void)
{
  int j;
  for(j = 0; j < 3; j++)
    {
      RootCenter[j] = SubDomainCenter[j];
      RootCorner[j] = SubDomainCorner[j];
    }
  RootLen        = SubDomainLen;
  RootInverseLen = SubDomainInverseLen;
  RootFac        = SubDomainFac;
  RootBigFac     = SubDomainBigFac;
}

/*! \brief Construct the subfind tree.
 *
 *  \param[in] npart Number of particles involved.
 *  \param[in] udp Unbind data.
 *
 *  \return Number of nodes.
 */
int subfind_loctree_treebuild(int npart, struct unbind_data **udp)
{
  int i, j, k, subnode = 0, parent = -1, numnodes;
  int nfree, th, nn;
  struct LocNODE *nfreep;
  struct unbind_data *mp;

  /* select first node */
  nfree  = LocMaxPart;
  nfreep = &LocNodes[nfree];

  mp = *udp;

  /* create an empty  root node  */
  nfreep->len = (MyFloat)RootLen;
  for(i = 0; i < 3; i++)
    nfreep->center[i] = (MyFloat)RootCenter[i];

  for(i = 0; i < 8; i++)
    nfreep->u.suns[i] = -1;

  numnodes = 1;
  nfreep++;
  nfree++;

  /* insert all particles */

  LocTree_IntPos_list =
      (unsigned long long *)mymalloc_movable(&LocTree_IntPos_list, "LocTree_IntPos_list", 3 * NumPart * sizeof(unsigned long long));

  for(k = 0; k < npart; k++)
    {
      if(mp)
        i = mp[k].index;
      else
        terminate("what?");

      MyDouble *posp;

#ifdef CELL_CENTER_GRAVITY
      if(P[i].Type == 0)
        posp = &PS[i].Center[0];
      else
#endif /* #ifdef CELL_CENTER_GRAVITY */
        posp = &P[i].Pos[0];

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
                   *      DomainLen/2^MAX_TREE_LEEL  < gravitational softening length
                   */
                  for(j = 0; j < 8; j++)
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
                  mp            = *udp;

                  if(numnodes > MaxNodes)
                    {
                      char buf[1000];

                      sprintf(buf, "maximum number %d of tree-nodes reached., for particle %d  %g %g %g", MaxNodes, i, P[i].Pos[0],
                              P[i].Pos[1], P[i].Pos[2]);
                      terminate(buf);
                    }
                }
            }
        }
    }

  myfree(LocTree_IntPos_list);

  /* now compute the multipole moments recursively */
  last = -1;
  subfind_loctree_update_node_recursive(LocMaxPart, -1, -1);

  if(last >= LocMaxPart)
    LocNodes[last].u.d.nextnode = -1;
  else
    LocNextNode[last] = -1;

  return numnodes;
}

/*! \brief Compute multipole moments.
 *
 *  This routine computes the multipole moments for a given internal node and
 *  all its subnodes using a recursive computation.
 *
 *  \param[in] no Node that we are in.
 *  \param[in] sib Sibling of the node.
 *  \param[in] father Parent node.
 *
 *  \return void
 */
void subfind_loctree_update_node_recursive(int no, int sib, int father)
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

              subfind_loctree_update_node_recursive(p, nextsib, no);

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
                  mass += P[p].Mass;
#ifdef CELL_CENTER_GRAVITY
                  if(P[p].Type == 0)
                    {
                      s[0] += P[p].Mass * PS[p].Center[0];
                      s[1] += P[p].Mass * PS[p].Center[1];
                      s[2] += P[p].Mass * PS[p].Center[2];
                    }
                  else
#endif /* #ifdef CELL_CENTER_GRAVITY */
                    {
                      s[0] += P[p].Mass * P[p].Pos[0];
                      s[1] += P[p].Mass * P[p].Pos[1];
                      s[2] += P[p].Mass * P[p].Pos[2];
                    }

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

/*! \brief Evaluates the potential by walking the subfind local tree.
 *
 *  \param[in] target Index of the particle.
 *
 *  \return Gravitational potiential.
 */
double subfind_loctree_treeevaluate_potential(int target)
{
  struct LocNODE *nop = 0;
  int no;
  double r2, dx, dy, dz, mass, r, u, h_i, h_j, hmax, h_inv, wp;
  double pot, pos_x, pos_y, pos_z, xtmp, ytmp, ztmp;

#ifdef CELL_CENTER_GRAVITY
  if(P[target].Type == 0)
    {
      pos_x = PS[target].Center[0];
      pos_y = PS[target].Center[1];
      pos_z = PS[target].Center[2];
    }
  else
#endif /* #ifdef CELL_CENTER_GRAVITY */
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
    }

  h_i = All.ForceSoftening[P[target].SofteningType];

  pot = 0;

  no = LocMaxPart;

  while(no >= 0)
    {
#ifdef MULTIPLE_NODE_SOFTENING
      int indi_flag1 = -1, indi_flag2 = 0;
#endif                    /* #ifdef MULTIPLE_NODE_SOFTENING */
      if(no < LocMaxPart) /* single particle */
        {
#ifdef CELL_CENTER_GRAVITY
          if(P[no].Type == 0)
            {
              dx = GRAVITY_NEAREST_X(PS[no].Center[0] - pos_x);
              dy = GRAVITY_NEAREST_Y(PS[no].Center[1] - pos_y);
              dz = GRAVITY_NEAREST_Z(PS[no].Center[2] - pos_z);
            }
          else
#endif /* #ifdef CELL_CENTER_GRAVITY */
            {
              dx = GRAVITY_NEAREST_X(P[no].Pos[0] - pos_x);
              dy = GRAVITY_NEAREST_Y(P[no].Pos[1] - pos_y);
              dz = GRAVITY_NEAREST_Z(P[no].Pos[2] - pos_z);
            }

          r2 = dx * dx + dy * dy + dz * dz;

          mass = P[no].Mass;

          h_j = All.ForceSoftening[P[no].SofteningType];

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

/*! \brief Comparison function for r2type objects.
 *
 *  Compares element r2.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a->r2 < b->r2.
 */
int subfind_locngb_compare_key(const void *a, const void *b)
{
  if(((r2type *)a)->r2 < (((r2type *)b)->r2))
    return -1;

  if(((r2type *)a)->r2 > (((r2type *)b)->r2))
    return +1;

  return 0;
}

/*! \brief Iterates on smoothing length of  neighbor search to get a desired
 *         number of neighbors.
 *
 *  \param[in] xyz Search center of neighbor search.
 *  \param[in] desngb Desired number of neighbors.
 *  \param[in] hguess Initial guess of smoothing length.
 *
 *  \return Distance of the outermost particle to seearch center.
 */
double subfind_locngb_treefind(MyDouble xyz[3], int desngb, double hguess)
{
  int numngb;
  double h2max;

  if(hguess == 0)
    terminate("hguess needed");

  while(1)
    {
      numngb = subfind_locngb_treefind_variable(xyz, hguess);

      if(numngb < desngb)
        {
          hguess *= 1.26;
          continue;
        }

      if(numngb >= desngb)
        {
          qsort(R2list, numngb, sizeof(r2type), subfind_locngb_compare_key);
          h2max = R2list[desngb - 1].r2;
          break;
        }

      hguess *= 1.26;
    }

  return sqrt(h2max);
}

/*! \brief (Local) tree-search in subfind tree.
 *
 *  Adds these cells to R2list.
 *
 *  \param[in] searchcenter Center around which particles are searched.
 *  \param[in] hguess Distance up to which particles are searched.
 *
 *  \return Number of neighbors found.
 */
int subfind_locngb_treefind_variable(MyDouble searchcenter[3], double hguess)
{
  int numngb, no, p;
  double dx, dy, dz, r2, h2;
  struct LocNODE *thisnode;
  double xtmp, ytmp, ztmp;

  h2 = hguess * hguess;

  numngb = 0;
  no     = LocMaxPart;

  while(no >= 0)
    {
      if(no < LocMaxPart) /* single particle */
        {
          p  = no;
          no = LocNextNode[no];
#ifdef CELL_CENTER_GRAVITY
          if(P[p].Type == 0)
            {
              dx = GRAVITY_NEAREST_X(PS[p].Center[0] - searchcenter[0]);
              dy = GRAVITY_NEAREST_Y(PS[p].Center[1] - searchcenter[1]);
              dz = GRAVITY_NEAREST_Z(PS[p].Center[2] - searchcenter[2]);
            }
          else
#endif /* #ifdef CELL_CENTER_GRAVITY */
            {
              dx = GRAVITY_NEAREST_X(P[p].Pos[0] - searchcenter[0]);
              dy = GRAVITY_NEAREST_Y(P[p].Pos[1] - searchcenter[1]);
              dz = GRAVITY_NEAREST_Z(P[p].Pos[2] - searchcenter[2]);
            }

          if(dx < -hguess)
            continue;
          if(dx > hguess)
            continue;

          if(dy < -hguess)
            continue;
          if(dy > hguess)
            continue;

          if(dz < -hguess)
            continue;
          if(dz > hguess)
            continue;

          r2 = dx * dx + dy * dy + dz * dz;

          if(r2 <= h2)
            {
              R2list[numngb].r2    = r2;
              R2list[numngb].index = p;
              numngb++;
            }
        }
      else
        {
          thisnode = &LocNodes[no];

          no = LocNodes[no].u.d.sibling; /* in case the node can be discarded */

          if((GRAVITY_NEAREST_X(thisnode->center[0] - searchcenter[0]) + 0.5 * thisnode->len) < -hguess)
            continue;
          if((GRAVITY_NEAREST_X(thisnode->center[0] - searchcenter[0]) - 0.5 * thisnode->len) > hguess)
            continue;
          if((GRAVITY_NEAREST_Y(thisnode->center[1] - searchcenter[1]) + 0.5 * thisnode->len) < -hguess)
            continue;
          if((GRAVITY_NEAREST_Y(thisnode->center[1] - searchcenter[1]) - 0.5 * thisnode->len) > hguess)
            continue;
          if((GRAVITY_NEAREST_Z(thisnode->center[2] - searchcenter[2]) + 0.5 * thisnode->len) < -hguess)
            continue;
          if((GRAVITY_NEAREST_Z(thisnode->center[2] - searchcenter[2]) - 0.5 * thisnode->len) > hguess)
            continue;

          no = thisnode->u.d.nextnode; /* ok, we need to open the node */
        }
    }

  return numngb;
}

/*! \brief Allocates memory used for storage of the tree
 *         and auxiliary arrays for tree-walk and link-lists.
 *
 *  \param[in] maxnodes Maximum number of nodes.
 *  \param[in] maxpart Maximum number of particles.
 *
 *  \return Size of allocated memory in bytes.
 */
size_t subfind_loctree_treeallocate(int maxnodes, int maxpart)
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

/*! \brief Frees the memory allocated for subfind_loctree.
 *
 *  \return void
 */
void subfind_loctree_treefree(void)
{
  myfree(LocNodes_base);
  myfree(R2list);
  myfree(LocNextNode);

  LocNextNode   = NULL;
  R2list        = NULL;
  LocNodes_base = NULL;
}

#endif /* #ifdef SUBFIND */
