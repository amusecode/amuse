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
 * \file        src/gravity/forcetree_walk.c
 * \date        05/2018
 * \brief       Gravitational tree walk code.
 * \details     This file contains the various gravitational tree walks.
 *              contains functions:
 *                void force_short_range_init(void)
 *                int force_treeevaluate(gravdata_in * in, gravdata_out * out,
 *                  int target, int mode, int thread_id, int numnodes, int
 *                  *firstnode, int measure_cost_flag)
 *                int tree_treefind_export_node_threads(int no, int i, int
 *                  thread_id)
 *                void force_evaluate_direct(int target, int result_idx,
 *                  int nimport)
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

/*! \brief Variable for short-range lookup table.
 *
 *  Contains the factor needed for the short range
 *  contribution of the tree to the gravity force.
 */
static float shortrange_table[NTAB + 1];

/*! \brief Variable for short-range lookup table.
 *
 *  Contains the factor needed for the short range
 *  contribution of the tree to the potential energy.
 */
static float shortrange_table_potential[NTAB + 1];

/*! \brief Initializes the short range table.
 *
 *  The short range table contains the complementary error function
 *  needed for the computation of the short range part of the gravity
 *  force/potential in case of the TreePM algorithm.
 *
 *  \return void
 */
void force_short_range_init(void)
{
  for(int i = 0; i <= NTAB; i++)
    {
      double u = ((RCUT / 2.0) / NTAB) * i;

      shortrange_table_potential[i] = -erfc(u); /* -r * g(r) */

      if(u > 0)
        shortrange_table[i] = (erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u) - 1.0) / (u * u); /* -g'(r) - 1/r^2 */
      else
        shortrange_table[i] = 0;
    }
}

/*! \brief This routine calculates the (short range) force contribution
 *   for a given particle in case the Tree(PM) algorithm is used.
 *
 *  In the TreePM algorithm, the tree is walked only locally around the
 *  target coordinate.  Tree nodes that fall outside a box of half
 *  side-length Rcut= RCUT*ASMTH*MeshSize can be discarded. The short-range
 *  potential is modified by a complementary error function, multiplied
 *  with the Newtonian form. The resulting short-range suppression compared
 *  to the Newtonian force is tabulated, because looking up from this table
 *  is faster than recomputing the corresponding factor, despite the
 *  memory-access penalty (which reduces cache performance) incurred by the
 *  table.
 *
 *  Depending on the value of TypeOfOpeningCriterion, either the geometrical BH
 *  cell-opening criterion, or the `relative' opening criterion is used.
 *
 *  \param[in] in Gravdata communicated into function.
 *  \param[in, out] out Gravdata communicated from function.
 *  \param[in] target Index of the particle to be processed.
 *  \param[in] mode 0: process local particle (phase 1), 1: process imported
 *             particle (phase 2).
 *  \param[in] thread_id Id of this thread.
 *  \param[in, out] firstnode First node involved in this algorithm.
 *  \param[in] measure_cost_flag Whether the cost of the tree walk should be
 *             measured.
 *
 *  \return Number of interactions processed for particle i.
 */
int force_treeevaluate(gravdata_in *in, gravdata_out *out, int target, int mode, int thread_id, int numnodes, int *firstnode,
                       int measure_cost_flag)
{
  struct NODE *nop = NULL;
#ifdef MULTIPLE_NODE_SOFTENING
  struct ExtNODE *extnop = 0;
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
#if !defined(GRAVITY_NOT_PERIODIC)
  double xtmp, ytmp, ztmp;
#endif /* #if !defined(GRAVITY_NOT_PERIODIC) */

  double acc_x = 0;
  double acc_y = 0;
  double acc_z = 0;
#ifdef EVALPOTENTIAL
  double pot = 0.0;
#endif /* #ifdef EVALPOTENTIAL */

  int ninteractions = 0;

  double pos_x = in->Pos[0];
  double pos_y = in->Pos[1];
  double pos_z = in->Pos[2];
  double aold  = All.ErrTolForceAcc * in->OldAcc;
  double h_i   = All.ForceSoftening[in->SofteningType];

#ifdef PMGRID
  double rcut  = All.Rcut[0];
  double asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
  if(pmforce_is_particle_high_res(in->Type, in->Pos))
    {
      rcut  = All.Rcut[1];
      asmth = All.Asmth[1];
    }
#endif /* #ifdef PLACEHIGHRESREGION */

  double rcut2     = rcut * rcut;
  double asmthinv  = 0.5 / asmth;
  double asmthinv2 = asmthinv * asmthinv;
  double asmthfac  = asmthinv * (NTAB / (RCUT / 2.0));
#endif /* #ifdef PMGRID */

  for(int k = 0; k < numnodes; k++)
    {
      int no;

      if(mode == 0)
        no = Tree_MaxPart; /* root node */
      else
        {
          no = firstnode[k];
          no = Nodes[no].u.d.nextnode; /* open it */
        }

      while(no >= 0)
        {
          double dx, dy, dz, r2, mass, hmax;

#ifdef MULTIPLE_NODE_SOFTENING
          int indi_flag1 = -1, indi_flag2 = 0;
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

          if(no < Tree_MaxPart) /* single particle */
            {
              dx = GRAVITY_NEAREST_X(Tree_Pos_list[3 * no + 0] - pos_x);
              dy = GRAVITY_NEAREST_Y(Tree_Pos_list[3 * no + 1] - pos_y);
              dz = GRAVITY_NEAREST_Z(Tree_Pos_list[3 * no + 2] - pos_z);
              r2 = dx * dx + dy * dy + dz * dz;

              mass = P[no].Mass;

              if(measure_cost_flag)
                Thread[thread_id].P_CostCount[no]++;

              double h_j = All.ForceSoftening[P[no].SofteningType];

              hmax = (h_j > h_i) ? h_j : h_i;

              no = Nextnode[no];
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes) /* we have an  internal node */
            {
              if(mode == 1)
                {
                  if(no <
                     Tree_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
                    {
                      no = -1;
                      continue;
                    }
                }

              nop = &Nodes[no];

              mass = nop->u.d.mass;
              dx   = GRAVITY_NEAREST_X(nop->u.d.s[0] - pos_x);
              dy   = GRAVITY_NEAREST_Y(nop->u.d.s[1] - pos_y);
              dz   = GRAVITY_NEAREST_Z(nop->u.d.s[2] - pos_z);

              r2 = dx * dx + dy * dy + dz * dz;

#if defined(PMGRID)
              if(r2 > rcut2)
                {
                  /* check whether we can stop walking along this branch */
                  double eff_dist = rcut + 0.5 * nop->len;

                  double dist = GRAVITY_NEAREST_X(nop->center[0] - pos_x);
                  if(dist < -eff_dist || dist > eff_dist)
                    {
                      no = nop->u.d.sibling;
                      continue;
                    }

                  dist = GRAVITY_NEAREST_Y(nop->center[1] - pos_y);
                  if(dist < -eff_dist || dist > eff_dist)
                    {
                      no = nop->u.d.sibling;
                      continue;
                    }

                  dist = GRAVITY_NEAREST_Z(nop->center[2] - pos_z);
                  if(dist < -eff_dist || dist > eff_dist)
                    {
                      no = nop->u.d.sibling;
                      continue;
                    }
                }
#endif /* #if defined(PMGRID) */

              if(All.ErrTolTheta) /* check Barnes-Hut opening criterion */
                {
                  if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
                    {
                      /* open cell */
                      no = nop->u.d.nextnode;
                      continue;
                    }
                }
              else /* check relative opening criterion */
                {
                  double len2 = nop->len * nop->len;

                  if(len2 > r2 * (1.2 * 1.2)) /* add a worst case protection */
                    {
                      /* open cell */
                      no = nop->u.d.nextnode;
                      continue;
                    }

                    // note that aold is strictly speaking |acceleration| / G
#ifdef ACTIVATE_MINIMUM_OPENING_ANGLE
                  if(mass * len2 > r2 * r2 * aold && len2 > r2 * (0.4 * 0.4))
#else  /* #ifdef ACTIVATE_MINIMUM_OPENING_ANGLE */
                  if(mass * len2 > r2 * r2 * aold)
#endif /* #ifdef ACTIVATE_MINIMUM_OPENING_ANGLE #else */
                    {
                      /* open cell */
                      no = nop->u.d.nextnode;
                      continue;
                    }

                  /* check in addition whether we lie inside or very close to the cell */
                  if(fabs(GRAVITY_NEAREST_X(nop->center[0] - pos_x)) < 0.60 * nop->len)
                    {
                      if(fabs(GRAVITY_NEAREST_Y(nop->center[1] - pos_y)) < 0.60 * nop->len)
                        {
                          if(fabs(GRAVITY_NEAREST_Z(nop->center[2] - pos_z)) < 0.60 * nop->len)
                            {
                              no = nop->u.d.nextnode;
                              continue;
                            }
                        }
                    }
                }

              double h_j = All.ForceSoftening[nop->u.d.maxsofttype];

              if(h_j > h_i)
                {
#ifdef MULTIPLE_NODE_SOFTENING
#ifdef ADAPTIVE_HYDRO_SOFTENING
                  if(nop->u.d.maxhydrosofttype != nop->u.d.minhydrosofttype)
                    if(ExtNodes[no].mass_per_type[0] > 0)
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

                /* ok, node can be used */
#ifdef MULTIPLE_NODE_SOFTENING
              extnop = &ExtNodes[no];
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
              if(measure_cost_flag && mass)
                Thread[thread_id].Node_CostCount[no]++;

              no = nop->u.d.sibling;
            }
          else if(no >= Tree_ImportedNodeOffset) /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;

              dx = GRAVITY_NEAREST_X(Tree_Points[n].Pos[0] - pos_x);
              dy = GRAVITY_NEAREST_Y(Tree_Points[n].Pos[1] - pos_y);
              dz = GRAVITY_NEAREST_Z(Tree_Points[n].Pos[2] - pos_z);

              r2 = dx * dx + dy * dy + dz * dz;

              mass = Tree_Points[n].Mass;

              if(measure_cost_flag)
                Thread[thread_id].TreePoints_CostCount[n]++;

              double h_j = All.ForceSoftening[Tree_Points[n].SofteningType];

              hmax = (h_j > h_i) ? h_j : h_i;

              no = Nextnode[no - Tree_MaxNodes];
            }
          else /* pseudo particle */
            {
              if(mode == 0)
                {
                  tree_treefind_export_node_threads(no, target, thread_id);
                }

              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }

          /* now evaluate the multipole moment */
          if(mass)
            {
              double r = sqrt(r2);

#ifdef PMGRID
              double tabentry = asmthfac * r;
              int tabindex    = (int)tabentry;

              if(tabindex < NTAB)
                {
                  double tabweight    = tabentry - tabindex;
                  double factor_force = (1.0 - tabweight) * shortrange_table[tabindex] + tabweight * shortrange_table[tabindex + 1];
#ifdef EVALPOTENTIAL
                  double factor_pot =
                      (1.0 - tabweight) * shortrange_table_potential[tabindex] + tabweight * shortrange_table_potential[tabindex + 1];
#endif /* #ifdef EVALPOTENTIAL */
#endif /* #ifdef PMGRID */

#ifdef MULTIPLE_NODE_SOFTENING
                  for(int type = indi_flag1; type < indi_flag2; type++)
                    {
                      if(type >= 0)
                        {
                          mass = extnop->mass_per_type[type];
                          double h_j;
#ifdef ADAPTIVE_HYDRO_SOFTENING
                          if(type == 0)
                            h_j = All.ForceSoftening[nop->u.d.maxhydrosofttype];
                          else
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
                            h_j = All.ForceSoftening[type];

                          hmax = (h_j > h_i) ? h_j : h_i;
                        }

                      if(mass)
                        {
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
                          double fac;
#ifdef EVALPOTENTIAL
                          double wp;
#endif /* #ifdef EVALPOTENTIAL */

                          if(r >= hmax)
                            {
                              double rinv  = 1.0 / r;
                              double rinv3 = rinv * rinv * rinv;
#ifdef PMGRID
                              fac = rinv3 + rinv * factor_force * asmthinv2; /* fac  = -g'(r)/r */
#ifdef EVALPOTENTIAL
                              wp = rinv * factor_pot; /* wp   = -g(r)    */
#endif                                                /* #ifdef EVALPOTENTIAL */
#else                                                 /* #ifdef PMGRID */
                  fac = rinv3;
#ifdef EVALPOTENTIAL
                  wp  = -rinv;
#endif /* #ifdef EVALPOTENTIAL */
#endif /* #ifdef PMGRID #else */
                            }
                          else
                            {
                              double h_inv  = 1.0 / hmax;
                              double h3_inv = h_inv * h_inv * h_inv;
                              double u      = r * h_inv;

                              if(u < 0.5)
                                {
                                  double u2 = u * u;
                                  fac       = h3_inv * (SOFTFAC1 + u2 * (SOFTFAC2 * u + SOFTFAC3));
#ifdef EVALPOTENTIAL
                                  wp = h_inv * (SOFTFAC4 + u2 * (SOFTFAC5 + u2 * (SOFTFAC6 * u + SOFTFAC7)));
#endif /* #ifdef EVALPOTENTIAL */
                                }
                              else
                                {
                                  double u2 = u * u;
                                  double u3 = u2 * u;
                                  fac       = h3_inv * (SOFTFAC8 + SOFTFAC9 * u + SOFTFAC10 * u2 + SOFTFAC11 * u3 + SOFTFAC12 / u3);
#ifdef EVALPOTENTIAL
                                  wp = h_inv * (SOFTFAC13 + SOFTFAC14 / u +
                                                u2 * (SOFTFAC1 + u * (SOFTFAC15 + u * (SOFTFAC16 + SOFTFAC17 * u))));
#endif /* #ifdef EVALPOTENTIAL */
                                }

#ifdef PMGRID
                              if(r > 0)
                                {
                                  double rinv = 1.0 / r;
                                  fac += rinv * factor_force * asmthinv2; /* fac  = -g'(r)/r */
#ifdef EVALPOTENTIAL
                                  wp += rinv * (factor_pot + 1.0); /* wp   = -g(r)    */
#endif                                                             /* #ifdef EVALPOTENTIAL */
                                }
#endif /* #ifdef PMGRID */
                            }

#ifdef EVALPOTENTIAL
                          pot += mass * wp;
#endif /* #ifdef EVALPOTENTIAL */
                          fac *= mass;

                          acc_x += dx * fac;
                          acc_y += dy * fac;
                          acc_z += dz * fac;

#if !defined(PMGRID) && defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) && !defined(ONEDIMS_SPHERICAL)
                          double fcorr[3];
                          ewald_corr(dx, dy, dz, fcorr);
                          acc_x += mass * fcorr[0];
                          acc_y += mass * fcorr[1];
                          acc_z += mass * fcorr[2];
#ifdef EVALPOTENTIAL
                          pot += mass * ewald_pot_corr(dx, dy, dz);
#endif /* #ifdef EVALPOTENTIAL */
#endif /* #if !defined(PMGRID) && defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) && !defined(ONEDIMS_SPHERICAL) */

#ifdef MULTIPLE_NODE_SOFTENING
                        }
                    }
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
                  ninteractions++;
#ifdef PMGRID
                }
#endif /* #ifdef PMGRID */
            }
        }
    }

  out->Acc[0] = acc_x;
  out->Acc[1] = acc_y;
  out->Acc[2] = acc_z;
#ifdef EVALPOTENTIAL
  out->Potential = pot;
#endif /* #ifdef EVALPOTENTIAL */
#ifdef NO_GRAVITY_TYPE
  if(in->Type == NO_GRAVITY_TYPE)
    {
      out->Acc[0] = 0.0;
      out->Acc[1] = 0.0;
      out->Acc[2] = 0.0;
#ifdef EVALPOTENTIAL
      out->Potential = 0.0;
#endif /* #ifdef EVALPOTENTIAL */
    }
#endif /* #ifdef NO_GRAVITY_TYPE */
#ifdef OUTPUTGRAVINTERACTIONS
  out->GravInteractions = ninteractions;
#endif /* #ifdef OUTPUTGRAVINTERACTIONS */

  return ninteractions;
}

/*! \brief Prepares node to be exported.
 *
 *  \param[in] no Index of node.
 *  \param[in] i Index of particle.
 *  \param[in] thread_id ID of thread.
 *
 *  \return 0
 */
int tree_treefind_export_node_threads(int no, int i, int thread_id)
{
  /* The task indicated by the pseudoparticle node */
  int task = DomainNewTask[no - (Tree_MaxPart + Tree_MaxNodes)];

  if(Thread[thread_id].Exportflag[task] != i)
    {
      Thread[thread_id].Exportflag[task]     = i;
      int nexp                               = Thread[thread_id].Nexport++;
      Thread[thread_id].PartList[nexp].Task  = task;
      Thread[thread_id].PartList[nexp].Index = i;
      Thread[thread_id].ExportSpace -= Thread[thread_id].ItemSize;
    }

  int nexp                      = Thread[thread_id].NexportNodes++;
  nexp                          = -1 - nexp;
  struct datanodelist *nodelist = (struct datanodelist *)(((char *)Thread[thread_id].PartList) + Thread[thread_id].InitialSpace);
  nodelist[nexp].Task           = task;
  nodelist[nexp].Index          = i;
  nodelist[nexp].Node           = DomainNodeIndex[no - (Tree_MaxPart + Tree_MaxNodes)];
  Thread[thread_id].ExportSpace -= sizeof(struct datanodelist) + sizeof(int);
  return 0;
}

#ifdef ALLOW_DIRECT_SUMMATION
/*! \brief Kernel of direct summation force calculation.
 *
 *  \param[in] target Index of particle in import array.
 *  \param[in] result_idx Index in result array.
 *  \param[in] nimport number of imported particles.
 *
 *  \return void
 */
void force_evaluate_direct(int target, int result_idx, int nimport)
{
#if !defined(GRAVITY_NOT_PERIODIC)
  double xtmp, ytmp, ztmp;
#endif /* #if !defined(GRAVITY_NOT_PERIODIC) */

  double acc_x = 0;
  double acc_y = 0;
  double acc_z = 0;
#ifdef EVALPOTENTIAL
  double pot = 0.0;
#endif /* #ifdef EVALPOTENTIAL */

  double pos_x = DirectDataAll[target].Pos[0];
  double pos_y = DirectDataAll[target].Pos[1];
  double pos_z = DirectDataAll[target].Pos[2];
  double h_i   = All.ForceSoftening[DirectDataAll[target].SofteningType];

#ifdef PMGRID
  double asmth = All.Asmth[0];
#if defined(PLACEHIGHRESREGION)
  int ptype_i = DirectDataAll[target].Type;
  if(pmforce_is_particle_high_res(ptype_i, DirectDataAll[target].Pos))
    asmth = All.Asmth[1];
#endif /* #if defined(PLACEHIGHRESREGION) */
  double asmthinv  = 0.5 / asmth;
  double asmthinv2 = asmthinv * asmthinv;
  double asmthfac  = asmthinv * (NTAB / (RCUT / 2.0));
#endif /* #ifdef PMGRID */

  for(int j = 0; j < nimport; j++)
    {
      double h_j = All.ForceSoftening[DirectDataAll[j].SofteningType];

      double hmax = (h_j > h_i) ? h_j : h_i;

      double dx = GRAVITY_NEAREST_X(DirectDataAll[j].Pos[0] - pos_x);
      double dy = GRAVITY_NEAREST_Y(DirectDataAll[j].Pos[1] - pos_y);
      double dz = GRAVITY_NEAREST_Z(DirectDataAll[j].Pos[2] - pos_z);

      double r2 = dx * dx + dy * dy + dz * dz;

      double mass = DirectDataAll[j].Mass;

      /* now evaluate the force component */

      double r = sqrt(r2);

#ifdef PMGRID
      double tabentry = asmthfac * r;
      int tabindex    = (int)tabentry;

      if(tabindex < NTAB)
        {
          double tabweight    = tabentry - tabindex;
          double factor_force = (1.0 - tabweight) * shortrange_table[tabindex] + tabweight * shortrange_table[tabindex + 1];
#ifdef EVALPOTENTIAL
          double factor_pot =
              (1.0 - tabweight) * shortrange_table_potential[tabindex] + tabweight * shortrange_table_potential[tabindex + 1];
#endif /* #ifdef EVALPOTENTIAL */
#endif /* #ifdef PMGRID */

          double fac;
#ifdef EVALPOTENTIAL
          double wp;
#endif /* #ifdef EVALPOTENTIAL */

          if(r >= hmax)
            {
              double rinv  = 1.0 / r;
              double rinv3 = rinv * rinv * rinv;
#ifdef PMGRID
              fac = rinv3 + rinv * factor_force * asmthinv2; /* fac  = -g'(r)/r */
#ifdef EVALPOTENTIAL
              wp = rinv * factor_pot; /* wp   = -g(r)    */
#endif                                /* #ifdef EVALPOTENTIAL */
#else                                 /* #ifdef PMGRID */
          fac = rinv3;
#ifdef EVALPOTENTIAL
          wp  = -rinv;
#endif /* #ifdef EVALPOTENTIAL */
#endif /* #ifdef PMGRID #else */
            }
          else
            {
              double h_inv  = 1.0 / hmax;
              double h3_inv = h_inv * h_inv * h_inv;
              double u      = r * h_inv;

              if(u < 0.5)
                {
                  double u2 = u * u;
                  fac       = h3_inv * (SOFTFAC1 + u2 * (SOFTFAC2 * u + SOFTFAC3));
#ifdef EVALPOTENTIAL
                  wp = h_inv * (SOFTFAC4 + u2 * (SOFTFAC5 + u2 * (SOFTFAC6 * u + SOFTFAC7)));
#endif /* #ifdef EVALPOTENTIAL */
                }
              else
                {
                  double u2 = u * u;
                  double u3 = u2 * u;
                  fac       = h3_inv * (SOFTFAC8 + SOFTFAC9 * u + SOFTFAC10 * u2 + SOFTFAC11 * u3 + SOFTFAC12 / u3);
#ifdef EVALPOTENTIAL
                  wp = h_inv * (SOFTFAC13 + SOFTFAC14 / u + u2 * (SOFTFAC1 + u * (SOFTFAC15 + u * (SOFTFAC16 + SOFTFAC17 * u))));
#endif /* #ifdef EVALPOTENTIAL */
                }
#ifdef PMGRID
              if(r > 0)
                {
                  double rinv = 1.0 / r;
                  fac += rinv * factor_force * asmthinv2; /* fac  = -g'(r)/r */
#ifdef EVALPOTENTIAL
                  wp += rinv * (factor_pot + 1.0); /* wp   = -g(r)    */
#endif                                             /* #ifdef EVALPOTENTIAL */
                }
#endif /* #ifdef PMGRID */
            }

#ifdef EVALPOTENTIAL
          pot += mass * wp;
#endif /* #ifdef EVALPOTENTIAL */
          fac *= mass;

          acc_x += dx * fac;
          acc_y += dy * fac;
          acc_z += dz * fac;

#if !defined(PMGRID) && defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) && !defined(ONEDIMS_SPHERICAL)
          {
            double fcorr[3];
            ewald_corr(dx, dy, dz, fcorr);
            acc_x += mass * fcorr[0];
            acc_y += mass * fcorr[1];
            acc_z += mass * fcorr[2];
#if defined(EVALPOTENTIAL)
            pot += mass * ewald_pot_corr(dx, dy, dz);
#endif /* #if defined(EVALPOTENTIAL) */
          }
#endif /* #if !defined(PMGRID) && defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) && !defined(ONEDIMS_SPHERICAL) */

#ifdef PMGRID
        }
#endif /* #ifdef PMGRID */
    }

  DirectAccOut[result_idx].Acc[0] = acc_x;
  DirectAccOut[result_idx].Acc[1] = acc_y;
  DirectAccOut[result_idx].Acc[2] = acc_z;
#ifdef EVALPOTENTIAL
  DirectAccOut[result_idx].Potential = pot;
#endif /* #ifdef EVALPOTENTIAL */
}
#endif /* #ifdef ALLOW_DIRECT_SUMMATION */
