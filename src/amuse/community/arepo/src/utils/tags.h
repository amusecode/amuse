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
 * \file        src/utils/tags.h
 * \date        05/2018
 * \brief       Tag defines.
 * \details     Choice of numbers for historic reasons.
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 28.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#define TAG_N 10 /*!< Various tags used for labelling MPI messages */
#define TAG_HEADER 11
#define TAG_PDATA 12
#define TAG_SPHDATA 13
#define TAG_KEY 14
#define TAG_GRAV_B 19
#define TAG_HYDRO_A 22
#define TAG_HYDRO_B 23
#define TAG_NFORTHISTASK 24
#define TAG_NONPERIOD_A 29
#define TAG_NONPERIOD_B 30
#define TAG_NONPERIOD_C 31
#define TAG_DENS_A 35
#define TAG_DENS_B 36
#define TAG_LOCALN 37
#define TAG_FOF_A 45
#define TAG_PDATA_SPH 70
#define TAG_KEY_SPH 71
#define TAG_BARRIER 85
#define TAG_NODE_DATA 105
