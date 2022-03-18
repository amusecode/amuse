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
 * \file        src/gravity/forcetree.h
 * \date        05/2018
 * \brief       Functions and data structurer for forcetree.
 * \details
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 28.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifndef FORCETREE_H
#define FORCETREE_H

#ifndef INLINE_FUNC
#define INLINE_FUNC
#endif /* #ifndef INLINE_FUNC */

typedef struct
{
  MyDouble Pos[3];
  float OldAcc;
  unsigned char Type;
  unsigned char SofteningType;

  int Firstnode;
} gravdata_in;

typedef struct
{
  MyFloat Acc[3];
#ifdef EVALPOTENTIAL
  MyFloat Potential;
#endif /* #ifdef EVALPOTENTIAL */
#ifdef OUTPUTGRAVINTERACTIONS
  int GravInteractions;
#endif /* #ifdef OUTPUTGRAVINTERACTIONS */

} gravdata_out;

#ifdef LONG_X
#define STRETCHX (LONG_X)
#else /* #ifdef LONG_X */
#define STRETCHX 1
#endif /* #ifdef LONG_X #else */

#ifdef LONG_Y
#define STRETCHY (LONG_Y)
#else /* #ifdef LONG_Y */
#define STRETCHY 1
#endif /* #ifdef LONG_Y #else */

#ifdef LONG_Z
#define STRETCHZ (LONG_Z)
#else /* #ifdef LONG_Z */
#define STRETCHZ 1
#endif /* #ifdef LONG_Z #else */

#define DBX 1
#define DBY 1
#define DBZ 1
#define DBX_EXTRA 0
#define DBY_EXTRA 0
#define DBZ_EXTRA 0

/*! length of lock-up table for short-range force kernel in TreePM algorithm */
#define NTAB 127

#if defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC)

#define EN 64

#define ENX (DBX * STRETCHX * EN)
#define ENY (DBY * STRETCHY * EN)
#define ENZ (DBZ * STRETCHZ * EN)

extern MyFloat Ewd_fcorrx[ENX + 1][ENY + 1][ENZ + 1];
extern MyFloat Ewd_fcorry[ENX + 1][ENY + 1][ENZ + 1];
extern MyFloat Ewd_fcorrz[ENX + 1][ENY + 1][ENZ + 1];
extern MyFloat Ewd_potcorr[ENX + 1][ENY + 1][ENZ + 1];
extern double Ewd_fac_intp;

extern int NTreeInsert;

#endif /* #if defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) */

#define MAX_TREE_LEVEL 30
#define MAX_TREE_ALLOC_FACTOR 30.0

#define TAKE_NSLOTS_IN_ONE_GO 32

#define MAX_IMPACT_BEFORE_OPTIMIZATION 1.03

#define BITFLAG_TOPLEVEL 0
#define BITFLAG_DEPENDS_ON_LOCAL_MASS 1
#define BITFLAG_DEPENDS_ON_EXTERN_MASS 2
#define BITFLAG_INTERNAL_TOPLEVEL 6
#define BITFLAG_MULTIPLEPARTICLES 7
#define BITFLAG_CONTAINS_GAS 10

#define BITFLAG_MASK ((1 << BITFLAG_CONTAINS_GAS) + (1 << BITFLAG_MULTIPLEPARTICLES))

static inline unsigned long long force_double_to_int(double d)
{
  union
  {
    double d;
    unsigned long long ull;
  } u;
  u.d = d;
  return (u.ull & 0xFFFFFFFFFFFFFllu);
}

static inline double force_int_to_double(unsigned long long x)
{
  union
  {
    double d;
    unsigned long long ull;
  } u;
  u.d = 1.0;
  u.ull |= x;
  return u.d;
}

int tree_treefind_export_node_threads(int no, int target, int thread_id);
int construct_forcetree(int mode, int optimized_domain_mapping, int insert_only_primary, int timebin);
int force_treebuild(int npart, int optimized_domain_mapping, int insert_only_primary, int timebin);
int force_treebuild_construct(int npart, int optimized_domain_mapping, int insert_only_primary, int timebin);
int force_treebuild_insert_single_point(int i, unsigned long long *intpos, int th, unsigned char level);
int force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z);
void force_insert_pseudo_particles(void);
void force_update_node_recursive(int no, int sib, int father, int *last);
void force_exchange_topleafdata(void);
void force_treeupdate_toplevel(int no, int topnode, int bits, int x, int y, int z);
void force_treeallocate(int maxpart, int maxindex);
void force_treefree(void);
void dump_particles(void);
int force_add_empty_nodes(void);
void force_short_range_init(void);
int force_treeevaluate(gravdata_in *in, gravdata_out *out, int target, int mode, int thread_id, int numnodes, int *firstnode,
                       int measure_cost_flag);
void force_assign_cost_values(void);
void force_optimize_domain_mapping(void);
double force_get_current_balance(double *impact);
void force_get_global_cost_for_leavenodes(int nexport);
void forcetest_ewald_init(void);

#endif /* #ifndef FORCETREE_H */
