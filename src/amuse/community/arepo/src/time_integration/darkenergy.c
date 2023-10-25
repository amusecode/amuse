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
 * \file        src/time_integration/darkenergy.c
 * \date        05/2018
 * \brief       Contains the hubble function for a LCDM cosmology.
 * \details     Using Dark Energy instead of a cosmological constant can be
 *              archived by replacing Lambda by Lambda * a^(-3*(1+w)) in the
 *              Hubble function. w = -1 gives back a  standard cosmological
 *              constant! Also w = -1/3 gives Lambda / a^2 which then cancel
 *              within the Hubble function and is then equal to the dynamics
 *              of a universe with Lambda = 0 !
 *
 *              For a time varying w once has to replace Lambda * a^(-3*(1+w))
 *              by Lambda * exp(Integral(a,1,3*(1+w)/a))
 *
 *              Dark Energy does not alter the powerspectrum of initial
 *              conditions. To get the same cluster for various values or
 *              functions of w, once has do assign a new redshift to the
 *              initial conditions to match the linear growth factors, so
 *              g(z=0)/g(z_ini) == g_w(z=0)/g_w(z_ini^new). Also the initial
 *              velocities field has to be scaled by
 *(Hubble_w(z_ini^new)*Omega_w(z_ini^new)^0.6)/(Hubble(z_ini)*Omega(z_ini)^0.6)
 *              where _w means the according functions including the terms for
 *              Dark Energy.
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 04.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/*! \brief Hubble function.
 *
 *  Returns the Hubble function at a given scalefactor for a LCDM cosmology.
 *
 *  \param[in] a Scalefactor.
 *
 *  \return Hubble parameter in internal units.
 */
double INLINE_FUNC hubble_function(double a)
{
  double hubble_a;

  hubble_a = All.Omega0 / (a * a * a) + (1 - All.Omega0 - All.OmegaLambda) / (a * a) + All.OmegaLambda;
  hubble_a = All.Hubble * sqrt(hubble_a);

  return (hubble_a);
}
