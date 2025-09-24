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
 * \file        src/domain.h
 * \date        05/2018
 * \brief       Header for domain decomposition.
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 28.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifndef ALLVARS_H
#include "../main/allvars.h"
#endif /* #ifndef ALLVARS_H */

#ifndef DOMAIN_H
#define DOMAIN_H

#define MASK_ACTIVE_FLAG_IN_TYPE 127
#define SET_ACTIVE_FLAG_IN_TYPE 128

enum domain_displace_mode
{
  DISPLACE_POSITION_FORWARD,
  DISPLACE_POSITION_BACKWARD
};

extern struct local_topnode_data
{
  peanokey Size;     /*!< number of Peano-Hilbert mesh-cells represented by top-level node */
  peanokey StartKey; /*!< first Peano-Hilbert key in top-level node */
  long long Count;   /*!< counts the number of particles in this top-level node */
  double Cost;
  double SphCost;
  int Daughter; /*!< index of first daughter cell (out of 8) of top-level node */
  int Leaf;     /*!< if the node is a leaf, this gives its number when all leaves are traversed in Peano-Hilbert order */
  int Parent;
  int PIndex; /*!< first particle in node */

} * topNodes, *branchNodes; /*!< points to the root node of the top-level tree */

struct domain_count_data
{
  int task;
  int count;
  int origintask;
};

extern struct domain_peano_hilbert_data
{
  peanokey key;
  int index;
} * mp;

extern struct trans_data
{
  MyIDType ID;
  int new_task;
  int new_index;
  int wrapped;
} * trans_table;

extern int N_trans;

extern int Nbranch;

extern double fac_work, fac_load, fac_worksph;
extern double normsum_work, normsum_load, normsum_worksph;

extern double totgravcost, totpartcount, gravcost, totsphcost, sphcost;

extern struct domain_cost_data
{
  int no;
  float Work;    /*!< total "work" due to the particles stored by a leave node */
  float WorkSph; /*!< total "work" due to the particles stored by a leave node */
  int Count;     /*!< a table that gives the total number of particles held by each processor */
  int CountSph;  /*!< a table that gives the total number of SPH particles held by each processor */
} * DomainLeaveNode;

/* toGo[partner] gives the number of particles on the current task that have to go to task 'partner'
 */
extern int *toGo, *toGoSph;
extern int *toGet, *toGetSph;
extern int *list_NumPart;
extern int *list_NumGas;
extern int *list_load;
extern int *list_loadsph;
extern double *list_work;
extern double *list_worksph;

/* functions for domain decomposition */
peano1D domain_double_to_int(double d);
double domain_grav_tot_costfactor(int i);
double domain_hydro_tot_costfactor(int i);
void domain_init_sum_cost(void);
void domain_printf(char *buf);
void domain_report_balance(void);
int domain_sort_load(const void *a, const void *b);
int domain_compare_count(const void *a, const void *b);
int domain_sort_task(const void *a, const void *b);
int domain_compare_count(const void *a, const void *b);
void domain_rearrange_particle_sequence(void);
void domain_combine_topleaves_to_domains(int ncpu, int ndomain);
void domain_combine_multipledomains(void);
void domain_allocate(void);
void domain_Decomposition(void);
int domain_compare_key(const void *a, const void *b);
int domain_countToGo(void);
int domain_determineTopTree(void);
void domain_exchange(void);
void domain_findExtent(void);
void domain_free(void);
void domain_sumCost(void);
void domain_walktoptree(int no);
void domain_optimize_domain_to_task_mapping(void);
int domain_compare_count(const void *a, const void *b);
void domain_allocate_lists(void);
void domain_free_lists(void);
int domain_unpack_tree_branch(int no, int parent);
void domain_do_local_refine(int n, int *list);
void domain_preserve_relevant_topnode_data(void);
void domain_find_total_cost(void);
void domain_voronoi_dynamic_update_execute(void);
void domain_prepare_voronoi_dynamic_update(void);
void domain_voronoi_dynamic_flag_particles(void);
void domain_mark_in_trans_table(int i, int task);
void domain_exchange_and_update_DC(void);
int domain_compare_connection_ID(const void *a, const void *b);
int domain_compare_local_trans_data_ID(const void *a, const void *b);
int domain_compare_recv_trans_data_ID(const void *a, const void *b);
int domain_compare_recv_trans_data_oldtask(const void *a, const void *b);
void mysort_domain(void *b, size_t n, size_t s);
void domain_displacePosition(MyDouble *pos, enum domain_displace_mode mode);

#endif /* #ifndef DOMAIN_H */
