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
 * \file        src/add_backgroundgrid/add_bggrid.h
 * \date        05/2018
 * \brief       Re-gridding of ICs to ensure that the entire computational
 *              domain contains gas cells.
 * \details     Can be used to convert SPH ICs to Arepo ICs.
 *              Interface functions:
 *                int add_backgroundgrid(void);
 *                void prepare_domain_backgroundgrid(void);
 *              Functions of this module called in:
 *                init() (init.c)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 11.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifndef ADD_BGGRID_H
#define ADD_BGGRID_H

#include "../main/allvars.h"

#ifdef ADDBACKGROUNDGRID

#define ADDBACKGROUNDGRIDMAX 256
#define FACTOR_MAX_BOX_SIZE 15.0
#define FACTOR_MIN_BOX_SIZE 2.0

extern MyIDType IDNew;

int add_backgroundgrid(void);
void prepare_domain_backgroundgrid(void);
void calculate_weights();
void distribute_particles();

#endif /* #ifdef ADDBACKGROUNDGRID */

#endif /* ADD_BGGRID_H */
