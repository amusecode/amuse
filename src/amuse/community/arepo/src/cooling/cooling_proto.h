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
 * \file        src/cooling/cooling_proto.h
 * \date        05/2018
 * \brief       Header for cooling functions.
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 27.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifndef INLINE_FUNC
#define INLINE_FUNC
#endif /* #ifndef INLINE_FUNC */

void SetOutputGasState(int i, double *ne_guess, double *nH0, double *coolrate);

double convert_u_to_temp(double u, double rho, double *ne_guess);
double CoolingRate(double logT, double rho, double *nelec);
double CoolingRateFromU(double u, double rho, double *ne_guess);
double DoCooling(double u_old, double rho, double dt, double *ne_guess);
double GetCoolingTime(double u_old, double rho, double *ne_guess);

void find_abundances_and_rates(double logT, double rho, double *ne_guess);
void InitCool(void);
void IonizeParamsUVB(void);
void IonizeParams(void);
void ReadIonizeParams(char *fname, int which);
void SetZeroIonization(void);
