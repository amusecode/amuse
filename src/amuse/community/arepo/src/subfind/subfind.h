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
 * \file        src/subfind/subfind.h
 * \date        05/2018
 * \brief       Header for subfind algorithm.
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 27.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifndef SUBFIND_H
#define SUBFIND_H

#include "../domain/domain.h"
#include "../main/allvars.h"

#define FIND_SMOOTHING_LENGTHS 0
#define FIND_TOTAL_DENSITIES 1
#define SUBFIND_SO_POT_CALCULATION_PARTICLE_NUMBER 10000
#define SUBFIND_GAL_RADIUS_FAC 2.0 /* for subfind metal calculation */

#if defined(SUBFIND) && defined(SUBFIND_EXTENDED_PROPERTIES)
extern int *NodeGrNr;
#endif /* #if defined(SUBFIND) && defined(SUBFIND_EXTENDED_PROPERTIES) */

extern int GrNr;
extern int NumPartGroup;

extern struct topnode_data *SubTopNodes;
extern struct local_topnode_data *Sub_LocTopNodes;

extern int *SubDomainTask;
extern int *SubDomainNodeIndex;
extern int *SubNextnode;
extern int SubNTopleaves;
extern int SubNTopnodes;

extern int SubTree_MaxPart;
extern int SubTree_NumNodes;
extern int SubTree_MaxNodes;
extern int SubTree_FirstNonTopLevelNode;
extern int SubTree_NumPartImported;
extern int SubTree_NumPartExported;
extern int SubTree_ImportedNodeOffset;
extern int SubTree_NextFreeNode;
extern MyDouble *SubTree_Pos_list;
extern struct NODE *SubNodes;
extern struct ExtNODE *SubExtNodes;

extern double SubTreeAllocFactor;

extern int *SubTree_ResultIndexList;
extern int *SubTree_Task_list;
extern unsigned long long *SubTree_IntPos_list;

extern double SubDomainCorner[3], SubDomainCenter[3], SubDomainLen, SubDomainFac;
extern double SubDomainInverseLen, SubDomainBigFac;

extern MyDouble GrCM[3];

extern int Ncollective;
extern int NprocsCollective;
extern int MaxNsubgroups;
extern int MaxNgbs;
extern int MaxSerialGroupLen;
extern r2type *R2list;

extern int CommSplitColor;
extern MPI_Comm SubComm;

extern int SubNTask, SubThisTask;
extern int SubTagOffset;

extern struct proc_assign_data
{
  int GrNr;
  int Len;
  int FirstTask;
  int NTask;
} * ProcAssign;

extern struct subgroup_properties
{
  int Len;
  int LenType[NTYPES];
  int GrNr;
  int SubNr;
  int SubParent;
  MyIDType SubMostBoundID;
  MyFloat Mass;
  MyFloat MassType[NTYPES];
  MyFloat SubVelDisp;
  MyFloat SubVmax;
  MyFloat SubVmaxRad;
  MyFloat SubHalfMassRad;
  MyFloat SubHalfMassRadType[NTYPES];
  MyFloat SubMassInRad;
  MyFloat SubMassInRadType[NTYPES];
  MyFloat SubMassInHalfRad;
  MyFloat SubMassInHalfRadType[NTYPES];
  MyFloat SubMassInMaxRad;
  MyFloat SubMassInMaxRadType[NTYPES];
  MyFloat Pos[3];
  MyFloat CM[3];
  MyFloat Vel[3];
  MyFloat Spin[3];

#ifdef MHD
  MyFloat Bfld_Halo, Bfld_Disk;
#endif /* #ifdef MHD */

#ifdef SUBFIND_EXTENDED_PROPERTIES
  MyFloat Ekin, Epot, Ethr;
  MyFloat J[3], Jdm[3], Jgas[3], Jstars[3], CMFrac, CMFracType[NTYPES];
  MyFloat J_inRad[3], Jdm_inRad[3], Jgas_inRad[3], Jstars_inRad[3], CMFrac_inRad, CMFracType_inRad[NTYPES];
  MyFloat J_inHalfRad[3], Jdm_inHalfRad[3], Jgas_inHalfRad[3], Jstars_inHalfRad[3], CMFrac_inHalfRad, CMFracType_inHalfRad[NTYPES];
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

#ifdef USE_SFR
  MyFloat Sfr, SfrInRad, SfrInHalfRad, SfrInMaxRad, GasMassSfr;
#endif /* #ifdef USE_SFR */
} * SubGroup;

extern struct nearest_r2_data
{
  double dist[2];
} * R2Loc;

extern struct nearest_ngb_data
{
  long long index[2];
  int count;
} * NgbLoc;

extern int NumPaux;

extern struct paux_data
{
  int TaskOfGr;
  int LocGrIndex;
  unsigned char Type;
  unsigned char SofteningType;
  MyDouble Pos[3];
  MyDouble Mass;
} * Paux;

extern struct submp_data
{
  int index;
  int GrNr;
  int OldIndex;
  MyFloat DM_Density;
} * submp;

extern struct cand_dat
{
  int head;
  int len;
  int nsub;
  int rank, subnr, parent;
  int bound_length;
} * candidates;

extern struct coll_cand_dat
{
  long long head;
  long long rank;
  int len;
  int nsub;
  int subnr, parent;
  int bound_length;
} * coll_candidates;

typedef struct
{
  double rho;
#ifdef SUBFIND_CALC_MORE
  double vx, vy, vz;
  double v2;
#endif
} SubDMData;

void subfind_determine_sub_halo_properties(struct unbind_data *d, int num, struct subgroup_properties *subgroup, int grnr, int subnr,
                                           int parallel_flag, int nsubgroups_cat);
int subfind_ngb_treefind_density(MyDouble searchcenter[3], double hsml, int target, int *startnode, int mode, int *exportflag,
                                 int *exportnodecount, int *exportindex, SubDMData *sub_dm_data);
int subfind_treefind_collective_export_node_threads(int no, int i, int thread_id);
void subfind_domain_do_local_refine(int n, int *list);
void assign_group_numbers_based_on_catalogue(int ngroups_cat, int nsubgroups_cat);
int subfind_compare_rlist_mhd(const void *a, const void *b);

#endif /* #ifndef SUBFIND_H */
