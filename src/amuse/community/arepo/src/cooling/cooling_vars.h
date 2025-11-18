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
 * \file        src/cooling/cooling_vars.h
 * \date        05/2018
 * \brief       Variables for cooling.
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 27.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#define NCOOLTAB 2000
#define SMALLNUM 1.0e-60
#define COOLLIM 0.1
#define HEATLIM 20.0
#define eV_to_K 11606.0
#define eV_to_erg 1.60184e-12
#define MAX_TABLESIZE 250 /* Max # of lines in TREECOOL */

/* data for gas state */
typedef struct
{
  double ne, necgs, nHcgs;
  double bH0, bHep, bff, aHp, aHep, aHepp, ad, geH0, geHe0, geHep;
  double gJH0ne, gJHe0ne, gJHepne;
  double nH0, nHp, nHep, nHe0, nHepp;
  double XH, yhelium;
  double mhboltz;
  double ethmin; /* minimum internal energy for neutral gas */
  double mu;
} GasState;

/* tabulated rates */
typedef struct
{
  double BetaH0, BetaHep, Betaff;
  double AlphaHp, AlphaHep, Alphad, AlphaHepp;
  double GammaeH0, GammaeHe0, GammaeHep;
} RateTable;

/* photo-ionization/heating rate table */
typedef struct
{
  float variable;       /* logz for UVB */
  float gH0, gHe, gHep; /* photo-ionization rates */
  float eH0, eHe, eHep; /* photo-heating rates */
} PhotoTable;

/* current interpolated photo-ionization/heating rates */
typedef struct
{
  char J_UV;
  double gJH0, gJHep, gJHe0, epsH0, epsHep, epsHe0;
} PhotoCurrent;

/* cooling data */
typedef struct
{
  double u_old_input, rho_input, dt_input, ne_guess_input;
} DoCoolData;
