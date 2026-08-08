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
 * \file        src/fof/fof.h
 * \date        05/2018
 * \brief       Header for Friend-of-Friends halo finder.
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 27.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifndef FOF_H
#define FOF_H

#include "../main/allvars.h"

extern int Ngroups, NgroupsExt, MaxNgroups, TotNgroups, Nsubgroups, TotNsubgroups;
extern int Nids;
extern long long TotNids;

extern int fof_OldMaxPart;
extern int fof_OldMaxPartSph;

extern double LinkL;
extern unsigned char *flag_node_inside_linkinglength;

#define BITFLAG_INSIDE_LINKINGLENGTH 1

#ifndef FOF_SECONDARY_LINK_TARGET_TYPES
#define FOF_SECONDARY_LINK_TARGET_TYPES FOF_PRIMARY_LINK_TYPES
#endif

extern struct group_properties
{
  int Len;
  MyIDType MinID;
  MyIDType MinIDTask;
  int GrNr;
  int LenType[NTYPES];
  MyFloat MassType[NTYPES];
  MyFloat Mass;
  MyDouble CM[3];
  MyFloat Vel[3];
  MyDouble Pos[3];

  MyDouble FirstPos[3];
#ifdef USE_SFR
  MyFloat Sfr;
#endif /* #ifdef USE_SFR */

#ifdef SUBFIND
  int TargetTask; /* primary CPU responsible for this group */
  int Nsubs;
  int FirstSub;
  MyFloat M_Mean200, R_Mean200;
  MyFloat M_Crit200, R_Crit200;
  MyFloat M_Crit500, R_Crit500;
  MyFloat M_TopHat200, R_TopHat200;
#ifdef SUBFIND_EXTENDED_PROPERTIES
  MyFloat J_Mean200[3], JDM_Mean200[3], JGas_Mean200[3], JStars_Mean200[3], MassType_Mean200[NTYPES], CMFrac_Mean200,
      CMFracType_Mean200[NTYPES];
  MyFloat J_Crit200[3], JDM_Crit200[3], JGas_Crit200[3], JStars_Crit200[3], MassType_Crit200[NTYPES], CMFrac_Crit200,
      CMFracType_Crit200[NTYPES];
  MyFloat J_Crit500[3], JDM_Crit500[3], JGas_Crit500[3], JStars_Crit500[3], MassType_Crit500[NTYPES], CMFrac_Crit500,
      CMFracType_Crit500[NTYPES];
  MyFloat J_TopHat200[3], JDM_TopHat200[3], JGas_TopHat200[3], JStars_TopHat200[3], MassType_TopHat200[NTYPES], CMFrac_TopHat200,
      CMFracType_TopHat200[NTYPES];
  int LenType_Mean200[NTYPES], LenType_Crit200[NTYPES], LenType_Crit500[NTYPES], LenType_TopHat200[NTYPES];
  MyFloat J[3], JDM[3], JGas[3], JStars[3], CMFrac, CMFracType[NTYPES];
  MyFloat Ekin, Epot, Ethr;
  MyFloat Ekin_Crit200, Epot_Crit200, Ethr_Crit200;
  MyFloat Ekin_Crit500, Epot_Crit500, Ethr_Crit500;
  MyFloat Ekin_Mean200, Epot_Mean200, Ethr_Mean200;
  MyFloat Ekin_TopHat200, Epot_TopHat200, Ethr_TopHat200;
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
#endif /* #ifdef SUBFIND */

} * Group;

struct data_aux_sort
{
  int OriginTask, OriginIndex;
  int TargetTask, TargetIndex;
  int GrNr;
  int Type;
  MyIDType ID;
#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT)
  MyIDType FileOrder;
#endif /* #if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) */
#ifdef SUBFIND
  int SubNr;
  MyFloat DM_BindingEnergy;
#endif /* #ifdef SUBFIND */
};

extern struct fof_particle_list
{
  MyIDType MinID;
  int MinIDTask;
  int Pindex;
} * FOF_PList;

extern struct fof_group_list
{
  MyIDType MinID;
  int MinIDTask;
  int LocCount;
  int ExtCount;
  int GrNr;
} * FOF_GList;

extern struct id_list
{
  MyIDType ID;
  int GrNr;
  int Type;
#ifdef SUBFIND
  int SubNr;
  MyFloat BindingEgy;
#endif /* #ifdef SUBFIND */
} * ID_list;

extern struct bit_flags
{
  unsigned char Nonlocal : 2, MinIDChanged : 2, Marked : 2, Changed : 2;
} * Flags;

struct fof_local_sort_data
{
  int targetindex;
  int index;
};

extern struct fof_subfind_header
{
  int Ngroups;
  int Nsubgroups;
  int Nids;
  int TotNgroups;
  int TotNsubgroups;
  long long TotNids;
  int num_files;
  double time;
  double redshift;
  double HubbleParam;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  int flag_doubleprecision;
} catalogue_header;

enum fof_subfind_iofields
{
  IO_FOF_LEN,
  IO_FOF_MTOT,
  IO_FOF_POS,
  IO_FOF_CM,
  IO_FOF_VEL,
  IO_FOF_LENTYPE,
  IO_FOF_MASSTYPE,
  IO_FOF_SFR,

  IO_FOF_M_MEAN200,
  IO_FOF_R_MEAN200,
  IO_FOF_M_CRIT200,
  IO_FOF_R_CRIT200,
  IO_FOF_M_TOPHAT200,
  IO_FOF_R_TOPHAT200,
  IO_FOF_M_CRIT500,
  IO_FOF_R_CRIT500,

#ifdef SUBFIND_EXTENDED_PROPERTIES
  IO_FOF_J_MEAN200,
  IO_FOF_JDM_MEAN200,
  IO_FOF_JGAS_MEAN200,
  IO_FOF_JSTARS_MEAN200,
  IO_FOF_MASSTYPE_MEAN200,
  IO_FOF_LENTYPE_MEAN200,
  IO_FOF_CMFRAC_MEAN200,
  IO_FOF_CMFRACTYPE_MEAN200,
  IO_FOF_J_CRIT200,
  IO_FOF_JDM_CRIT200,
  IO_FOF_JGAS_CRIT200,
  IO_FOF_JSTARS_CRIT200,
  IO_FOF_MASSTYPE_CRIT200,
  IO_FOF_LENTYPE_CRIT200,
  IO_FOF_CMFRAC_CRIT200,
  IO_FOF_CMFRACTYPE_CRIT200,
  IO_FOF_J_TOPHAT200,
  IO_FOF_JDM_TOPHAT200,
  IO_FOF_JGAS_TOPHAT200,
  IO_FOF_JSTARS_TOPHAT200,
  IO_FOF_MASSTYPE_TOPHAT200,
  IO_FOF_LENTYPE_TOPHAT200,
  IO_FOF_CMFRAC_TOPHAT200,
  IO_FOF_CMFRACTYPE_TOPHAT200,
  IO_FOF_J_CRIT500,
  IO_FOF_JDM_CRIT500,
  IO_FOF_JGAS_CRIT500,
  IO_FOF_JSTARS_CRIT500,
  IO_FOF_MASSTYPE_CRIT500,
  IO_FOF_LENTYPE_CRIT500,
  IO_FOF_CMFRAC_CRIT500,
  IO_FOF_CMFRACTYPE_CRIT500,
  IO_FOF_J,
  IO_FOF_JDM,
  IO_FOF_JGAS,
  IO_FOF_JSTARS,
  IO_FOF_CMFRAC,
  IO_FOF_CMFRACTYPE,
  IO_FOF_EKIN,
  IO_FOF_ETHR,
  IO_FOF_EPOT,
  IO_FOF_EPOT_CRIT200,
  IO_FOF_EKIN_CRIT200,
  IO_FOF_ETHR_CRIT200,
  IO_FOF_EPOT_MEAN200,
  IO_FOF_EKIN_MEAN200,
  IO_FOF_ETHR_MEAN200,
  IO_FOF_EPOT_TOPHAT200,
  IO_FOF_EKIN_TOPHAT200,
  IO_FOF_ETHR_TOPHAT200,
  IO_FOF_EPOT_CRIT500,
  IO_FOF_EKIN_CRIT500,
  IO_FOF_ETHR_CRIT500,
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

  IO_FOF_NSUBS,
  IO_FOF_FIRSTSUB,
  IO_FOF_FUZZOFFTYPE,

  IO_SUB_LEN,
  IO_SUB_MTOT,
  IO_SUB_POS,
  IO_SUB_VEL,
  IO_SUB_LENTYPE,
  IO_SUB_MASSTYPE,
  IO_SUB_CM,
  IO_SUB_SPIN,
  IO_SUB_BFLD_HALO,
  IO_SUB_BFLD_DISK,

#ifdef SUBFIND_EXTENDED_PROPERTIES
  IO_SUB_EKIN,
  IO_SUB_ETHR,
  IO_SUB_EPOT,
  IO_SUB_J,
  IO_SUB_JDM,
  IO_SUB_JGAS,
  IO_SUB_JSTARS,
  IO_SUB_JINHALFRAD,
  IO_SUB_JDMINHALFRAD,
  IO_SUB_JGASINHALFRAD,
  IO_SUB_JSTARSINHALFRAD,
  IO_SUB_JINRAD,
  IO_SUB_JDMINRAD,
  IO_SUB_JGASINRAD,
  IO_SUB_JSTARSINRAD,
  IO_SUB_CMFRAC,
  IO_SUB_CMFRACTYPE,
  IO_SUB_CMFRACINHALFRAD,
  IO_SUB_CMFRACTYPEINHALFRAD,
  IO_SUB_CMFRACINRAD,
  IO_SUB_CMFRACTYPEINRAD,
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

  IO_SUB_VELDISP,
  IO_SUB_VMAX,
  IO_SUB_VMAXRAD,
  IO_SUB_HALFMASSRAD,
  IO_SUB_HALFMASSRADTYPE,
  IO_SUB_MASSINRAD,
  IO_SUB_MASSINHALFRAD,
  IO_SUB_MASSINMAXRAD,
  IO_SUB_MASSINRADTYPE,
  IO_SUB_MASSINHALFRADTYPE,
  IO_SUB_MASSINMAXRADTYPE,
  IO_SUB_IDMOSTBOUND,
  IO_SUB_GRNR,
  IO_SUB_PARENT,
  IO_SUB_SFR,
  IO_SUB_SFRINRAD,
  IO_SUB_SFRINHALFRAD,
  IO_SUB_SFRINMAXRAD,
  IO_FOFSUB_IDS,
  IO_FOF_LASTENTRY
};

int fof_subfind_blockpresent(enum fof_subfind_iofields blocknr);
int fof_subfind_get_datatype(enum fof_subfind_iofields blocknr);
int fof_subfind_get_bytes_per_blockelement(enum fof_subfind_iofields blocknr);
int fof_subfind_get_particles_in_block(enum fof_subfind_iofields blocknr);
void fof_subfind_get_dataset_name(enum fof_subfind_iofields blocknr, char *label);
void fof_subfind_get_Tab_IO_Label(enum fof_subfind_iofields blocknr, char *label);
int fof_subfind_get_dataset_group(enum fof_subfind_iofields blocknr);
void fof_subfind_fill_write_buffer(enum fof_subfind_iofields blocknr, int *startindex, int pc);
int fof_subfind_get_values_per_blockelement(enum fof_subfind_iofields blocknr);

#endif /* #ifndef FOF_H */
