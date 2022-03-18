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
 * \file        src/fof/fof_vars.c
 * \date        05/2018
 * \brief       Iinstances for the global variables used by FOF, which are
 *              declared in fof.h
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 24.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_math.h>
#include <inttypes.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"
#include "../subfind/subfind.h"
#include "fof.h"

#ifdef FOF

int Ngroups, NgroupsExt, MaxNgroups, TotNgroups, Nsubgroups, TotNsubgroups;
int Nids;
long long TotNids;

double LinkL = 0;

int fof_OldMaxPart;
int fof_OldMaxPartSph;

unsigned char *flag_node_inside_linkinglength;

struct group_properties *Group;

struct fofdata_in *FoFDataIn, *FoFDataGet;

struct fofdata_out *FoFDataResult, *FoFDataOut;

struct fof_particle_list *FOF_PList;

struct fof_group_list *FOF_GList;

struct id_list *ID_list;

struct bit_flags *Flags;

struct fof_subfind_header catalogue_header;

#endif /* #ifdef FOF */
