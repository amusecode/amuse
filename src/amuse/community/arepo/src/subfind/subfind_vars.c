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
 * \file        src/subfind/subfind_vars.c
 * \date        05/2018
 * \brief       Variables for the subfind algorithm.
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 14.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include "../main/allvars.h"

#ifdef SUBFIND

#include "../domain/domain.h"
#include "../fof/fof.h"
#include "subfind.h"

double SubDomainCorner[3], SubDomainCenter[3], SubDomainLen, SubDomainFac;
double SubDomainInverseLen, SubDomainBigFac;

MyDouble GrCM[3];

int GrNr;
int NumPartGroup;

MPI_Comm SubComm;
int CommSplitColor;
int SubNTask, SubThisTask;
int SubTagOffset;

struct topnode_data *SubTopNodes;
struct local_topnode_data *Sub_LocTopNodes;

double SubTreeAllocFactor;

#if defined(SUBFIND) && defined(SUBFIND_EXTENDED_PROPERTIES)
int *NodeGrNr;
#endif

int *SubDomainTask;
int *SubDomainNodeIndex;
int *SubNextnode;
int SubNTopleaves;
int SubNTopnodes;

int SubTree_MaxPart;
int SubTree_NumNodes;
int SubTree_MaxNodes;
int SubTree_FirstNonTopLevelNode;
int SubTree_NumPartImported;
int SubTree_NumPartExported;
int SubTree_ImportedNodeOffset;
int SubTree_NextFreeNode;
struct NODE *SubNodes;
struct ExtNODE *SubExtNodes;
int *SubTree_ResultIndexList;
int *SubTree_Task_list;
unsigned long long *SubTree_IntPos_list;
MyDouble *SubTree_Pos_list;

int Ncollective;
int NprocsCollective;
int MaxNsubgroups = 0;
int MaxNgbs;
int MaxSerialGroupLen;

r2type *R2list;

int NumPaux;

struct paux_data *Paux;
struct proc_assign_data *ProcAssign;
struct subgroup_properties *SubGroup;
struct nearest_r2_data *R2Loc;
struct nearest_ngb_data *NgbLoc;
struct submp_data *submp;
struct cand_dat *candidates;
struct coll_cand_dat *coll_candidates;

#endif /* #ifdef SUBFIND */
