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
 * \file        src/cooling/cooling.c
 * \date        05/2018
 * \brief       Module for gas radiative cooling
 * \details     contains functions:
 *                double DoCooling(double u_old, double rho, double dt, double
 *                  *ne_guess)
 *                double GetCoolingTime(double u_old, double rho, double
 *                  *ne_guess)
 *                double convert_u_to_temp(double u, double rho, double
 *                  *ne_guess)
 *                void find_abundances_and_rates(double logT, double rho,
 *                  double *ne_guess)
 *                double CoolingRateFromU(double u, double rho, double
 *                  *ne_guess)
 *                void SetOutputGasState(int i, double *ne_guess, double *nH0,
 *                  double *coolrate)
 *                double CoolingRate(double logT, double rho, double *nelec)
 *                void MakeRateTable(void)
 *                void ReadIonizeParams(char *fname, int which)
 *                void IonizeParamsUVB(void)
 *                void SetZeroIonization(void)
 *                void IonizeParams(void)
 *                void InitCool(void)
 *                void cooling_only(void)
 *                void cool_cell(int i)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 24.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef COOLING

static double Tmin = 0.0;     /*!< min temperature in log10 */
static double Tmax = 9.0;     /*!< max temperature in log10 */
static double deltaT;         /*!< log10 of temperature spacing in the interpolation tables */
static GasState gs;           /*!< gas state */
static RateTable *RateT;      /*!< tabulated rates */
static PhotoTable *PhotoTUVB; /*!< photo-ionization/heating rate table for UV background */
static PhotoCurrent pc;       /*!< current interpolated photo rates */
static int NheattabUVB;       /*!< length of UVB photo table */
static DoCoolData DoCool;     /*!< cooling data */

/*! \brief Computes the new internal energy per unit mass.
 *
 *  The function solves for the new internal energy per unit mass of the gas
 *  by integrating the equation for the internal energy with an implicit
 *  Euler scheme. The root of resulting non linear equation,
 *  which gives tnew internal energy, is found with the bisection method.
 *  Arguments are passed in code units.
 *
 *  \param[in] u_old the initial (before cooling is applied) internal energy
 *             per unit mass of the gas cell.
 *  \param[in] rho   the proper density of the gas cell.
 *  \param[in] dt    the duration of the time step.
 *  \param[in] ne_guess electron number density relative to hydrogen number
 *             density (for molecular weight computation).
 *
 *  \return The new internal energy per unit mass of the gas cell.
 */
double DoCooling(double u_old, double rho, double dt, double *ne_guess)
{
  double u, du;
  double u_lower, u_upper;
  double ratefact;
  double LambdaNet;

  int iter = 0;

  DoCool.u_old_input    = u_old;
  DoCool.rho_input      = rho;
  DoCool.dt_input       = dt;
  DoCool.ne_guess_input = *ne_guess;

  if(!gsl_finite(u_old))
    terminate("invalid input: u_old=%g\n", u_old);

  if(u_old < 0 || rho < 0)
    terminate("invalid input: task=%d u_old=%g  rho=%g  dt=%g  All.MinEgySpec=%g\n", ThisTask, u_old, rho, dt, All.MinEgySpec);

  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam; /* convert to physical cgs units */
  u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
  dt *= All.UnitTime_in_s / All.HubbleParam;

  gs.nHcgs = gs.XH * rho / PROTONMASS; /* hydrogen number dens in cgs units */
  ratefact = gs.nHcgs * gs.nHcgs / rho;

  u       = u_old;
  u_lower = u;
  u_upper = u;

  LambdaNet = CoolingRateFromU(u, rho, ne_guess);

  /* bracketing */
  if(u - u_old - ratefact * LambdaNet * dt < 0) /* heating */
    {
      u_upper *= sqrt(1.1);
      u_lower /= sqrt(1.1);
      while(u_upper - u_old - ratefact * CoolingRateFromU(u_upper, rho, ne_guess) * dt < 0)
        {
          u_upper *= 1.1;
          u_lower *= 1.1;
        }
    }

  if(u - u_old - ratefact * LambdaNet * dt > 0)
    {
      u_lower /= sqrt(1.1);
      u_upper *= sqrt(1.1);
      while(u_lower - u_old - ratefact * CoolingRateFromU(u_lower, rho, ne_guess) * dt > 0)
        {
          u_upper /= 1.1;
          u_lower /= 1.1;
        }
    }

  do
    {
      u = 0.5 * (u_lower + u_upper);

      LambdaNet = CoolingRateFromU(u, rho, ne_guess);

      if(u - u_old - ratefact * LambdaNet * dt > 0)
        {
          u_upper = u;
        }
      else
        {
          u_lower = u;
        }

      du = u_upper - u_lower;

      iter++;

      if(iter >= (MAXITER - 10))
        printf("u= %g\n", u);
    }
  while(fabs(du / u) > 1.0e-6 && iter < MAXITER);

  if(iter >= MAXITER)
    terminate(
        "failed to converge in DoCooling(): DoCool.u_old_input=%g\nDoCool.rho_input= %g\nDoCool.dt_input= %g\nDoCool.ne_guess_input= "
        "%g\n",
        DoCool.u_old_input, DoCool.rho_input, DoCool.dt_input, DoCool.ne_guess_input);

  u *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs; /* to internal units */

  return u;
}

/*! \brief Returns the cooling time.
 *
 *  If we actually have heating, a cooling time of 0 is returned.
 *
 *  \param[in] u_old The initial (before cooling is applied) internal energy
 *             per unit mass of the gas cell.
 *  \param[in] rho The proper density of the gas cell.
 *  \param[in] ne_guess Electron number density relative to hydrogen number
 *             density (for molecular weight computation).
 *
 *  \return Cooling time; 0 if heating.
 */
double GetCoolingTime(double u_old, double rho, double *ne_guess)
{
  double u;
  double ratefact;
  double LambdaNet, coolingtime;

  DoCool.u_old_input    = u_old;
  DoCool.rho_input      = rho;
  DoCool.ne_guess_input = *ne_guess;

  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam; /* convert to physical cgs units */
  u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

  gs.nHcgs = gs.XH * rho / PROTONMASS; /* hydrogen number dens in cgs units */
  ratefact = gs.nHcgs * gs.nHcgs / rho;

  u = u_old;

  LambdaNet = CoolingRateFromU(u, rho, ne_guess);

  if(LambdaNet >= 0) /* ups, we have actually heating due to UV background */
    return 0;

  coolingtime = u_old / (-ratefact * LambdaNet);

  coolingtime *= All.HubbleParam / All.UnitTime_in_s;

  return coolingtime;
}

/*! \brief Compute gas temperature from internal energy per unit mass.
 *
 *   This function determines the electron fraction, and hence the mean
 *   molecular weight. With it arrives at a self-consistent temperature.
 *   Element abundances and the rates for the emission are also computed.
 *
 *  \param[in] u   internal energy per unit mass.
 *  \param[in] rho gas density.
 *  \param[in, out] ne_guess electron number density relative to hydrogen
 *                  number density
 *
 *  \return The gas temperature.
 */
double convert_u_to_temp(double u, double rho, double *ne_guess)
{
  double temp, temp_old, temp_new, max = 0, ne_old;
  double mu;
  int iter = 0;

  double u_input, rho_input, ne_input;

  u_input   = u;
  rho_input = rho;
  ne_input  = *ne_guess;

  mu   = (1 + 4 * gs.yhelium) / (1 + gs.yhelium + *ne_guess);
  temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

  do
    {
      ne_old = *ne_guess;

      find_abundances_and_rates(log10(temp), rho, ne_guess);
      temp_old = temp;

      mu = (1 + 4 * gs.yhelium) / (1 + gs.yhelium + *ne_guess);

      temp_new = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

      max = dmax(max, temp_new / (1 + gs.yhelium + *ne_guess) * fabs((*ne_guess - ne_old) / (temp_new - temp_old + 1.0)));

      temp = temp_old + (temp_new - temp_old) / (1 + max);
      iter++;

      if(iter > (MAXITER - 10))
        printf("-> temp= %g ne=%g\n", temp, *ne_guess);
    }
  while(fabs(temp - temp_old) > 1.0e-3 * temp && iter < MAXITER);

  if(iter >= MAXITER)
    {
      printf("failed to converge in convert_u_to_temp()\n");
      printf("u_input= %g\nrho_input=%g\n ne_input=%g\n", u_input, rho_input, ne_input);
      printf("DoCool.u_old_input=%g\nDoCool.rho_input= %g\nDoCool.dt_input= %g\nDoCool.ne_guess_input= %g\n", DoCool.u_old_input,
             DoCool.rho_input, DoCool.dt_input, DoCool.ne_guess_input);
      terminate("convergence failure");
    }

  gs.mu = mu;

  return temp;
}

/*! \brief Computes the actual abundance ratios.
 *
 *  The chemical composition of the gas is primordial (no metals are present).
 *
 *  \param[in] logT log10 of gas temperature.
 *  \param[in] rho Gas density.
 *  \param[in, out] ne_guess Electron number density relative to hydrogen
 *                  number density.
 *
 *  \return void
 */
void find_abundances_and_rates(double logT, double rho, double *ne_guess)
{
  double neold, nenew;
  int j, niter;
  double flow, fhi, t;

  double logT_input, rho_input, ne_input;

  logT_input = logT;
  rho_input  = rho;
  ne_input   = *ne_guess;

  if(!gsl_finite(logT))
    terminate("logT=%g\n", logT);

  if(logT <= Tmin) /* everything neutral */
    {
      gs.nH0    = 1.0;
      gs.nHe0   = gs.yhelium;
      gs.nHp    = 0;
      gs.nHep   = 0;
      gs.nHepp  = 0;
      gs.ne     = 0;
      *ne_guess = 0;
      return;
    }

  if(logT >= Tmax) /* everything is ionized */
    {
      gs.nH0    = 0;
      gs.nHe0   = 0;
      gs.nHp    = 1.0;
      gs.nHep   = 0;
      gs.nHepp  = gs.yhelium;
      gs.ne     = gs.nHp + 2.0 * gs.nHepp;
      *ne_guess = gs.ne; /* note: in units of the hydrogen number density */
      return;
    }

  t    = (logT - Tmin) / deltaT;
  j    = (int)t;
  fhi  = t - j;
  flow = 1 - fhi;

  if(*ne_guess == 0)
    *ne_guess = 1.0;

  gs.nHcgs = gs.XH * rho / PROTONMASS; /* hydrogen number dens in cgs units */

  gs.ne    = *ne_guess;
  neold    = gs.ne;
  niter    = 0;
  gs.necgs = gs.ne * gs.nHcgs;

  /* evaluate number densities iteratively (cf KWH eqns 33-38) in units of nH */
  do
    {
      niter++;

      gs.aHp   = flow * RateT[j].AlphaHp + fhi * RateT[j + 1].AlphaHp;
      gs.aHep  = flow * RateT[j].AlphaHep + fhi * RateT[j + 1].AlphaHep;
      gs.aHepp = flow * RateT[j].AlphaHepp + fhi * RateT[j + 1].AlphaHepp;
      gs.ad    = flow * RateT[j].Alphad + fhi * RateT[j + 1].Alphad;
      gs.geH0  = flow * RateT[j].GammaeH0 + fhi * RateT[j + 1].GammaeH0;
      gs.geHe0 = flow * RateT[j].GammaeHe0 + fhi * RateT[j + 1].GammaeHe0;
      gs.geHep = flow * RateT[j].GammaeHep + fhi * RateT[j + 1].GammaeHep;

      if(gs.necgs <= 1.e-25 || pc.J_UV == 0)
        {
          gs.gJH0ne = gs.gJHe0ne = gs.gJHepne = 0;
        }
      else
        {
          gs.gJH0ne  = pc.gJH0 / gs.necgs;
          gs.gJHe0ne = pc.gJHe0 / gs.necgs;
          gs.gJHepne = pc.gJHep / gs.necgs;
        }

      gs.nH0 = gs.aHp / (gs.aHp + gs.geH0 + gs.gJH0ne); /* eqn (33) */
      gs.nHp = 1.0 - gs.nH0;                            /* eqn (34) */

      if((gs.gJHe0ne + gs.geHe0) <= SMALLNUM) /* no ionization at all */
        {
          gs.nHep  = 0.0;
          gs.nHepp = 0.0;
          gs.nHe0  = gs.yhelium;
        }
      else
        {
          gs.nHep =
              gs.yhelium / (1.0 + (gs.aHep + gs.ad) / (gs.geHe0 + gs.gJHe0ne) + (gs.geHep + gs.gJHepne) / gs.aHepp); /* eqn (35) */
          gs.nHe0  = gs.nHep * (gs.aHep + gs.ad) / (gs.geHe0 + gs.gJHe0ne);                                          /* eqn (36) */
          gs.nHepp = gs.nHep * (gs.geHep + gs.gJHepne) / gs.aHepp;                                                   /* eqn (37) */
        }

      neold = gs.ne;

      gs.ne    = gs.nHp + gs.nHep + 2 * gs.nHepp; /* eqn (38) */
      gs.necgs = gs.ne * gs.nHcgs;

      if(pc.J_UV == 0)
        break;

      nenew    = 0.5 * (gs.ne + neold);
      gs.ne    = nenew;
      gs.necgs = gs.ne * gs.nHcgs;

      if(fabs(gs.ne - neold) < 1.0e-4)
        break;

      if(niter > (MAXITER - 10))
        printf("ne= %g  niter=%d\n", gs.ne, niter);
    }
  while(niter < MAXITER);

  if(niter >= MAXITER)
    {
      printf("gs.aHp = %le\n", gs.aHp);
      char buff[1000];
      sprintf(buff, "%s/cooling_task%d.dat", All.OutputDir, ThisTask);
      FILE *fp = fopen(buff, "w");
      fwrite(&All.Time, sizeof(double), 1, fp);
      fwrite(&logT_input, sizeof(double), 1, fp);
      fwrite(&rho_input, sizeof(double), 1, fp);
      fwrite(&ne_input, sizeof(double), 1, fp);
      fclose(fp);
      terminate(
          "no convergence reached in find_abundances_and_rates(): logT_input= %g  rho_input= %g  ne_input= %g "
          "DoCool.u_old_input=%g\nDoCool.rho_input= %g\nDoCool.dt_input= %g\nDoCool.ne_guess_input= %g\n",
          logT_input, rho_input, ne_input, DoCool.u_old_input, DoCool.rho_input, DoCool.dt_input, DoCool.ne_guess_input);
    }
  gs.bH0  = flow * RateT[j].BetaH0 + fhi * RateT[j + 1].BetaH0;
  gs.bHep = flow * RateT[j].BetaHep + fhi * RateT[j + 1].BetaHep;
  gs.bff  = flow * RateT[j].Betaff + fhi * RateT[j + 1].Betaff;

  *ne_guess = gs.ne;
}

/*! \brief Get cooling rate from gas internal energy.
 *
 *  This function first computes the self-consistent temperature
 *  and abundance ratios, and then it calculates
 *  (heating rate-cooling rate)/n_h^2 in cgs units.
 *
 *  \param[in] u Gas internal energy per unit mass.
 *  \param[in] rho Gas density.
 *  \param[in, out] ne_guess Electron number density relative to hydrogen
 *                  number density.
 *
 *  \return Cooling rate.
 */
double CoolingRateFromU(double u, double rho, double *ne_guess)
{
  double temp;

  temp = convert_u_to_temp(u, rho, ne_guess);

  return CoolingRate(log10(temp), rho, ne_guess);
}

/*! \brief  This function computes the self-consistent temperature and
 *          abundance ratios.
 *
 *  Used only in io_fields.c for calculating output fields.
 *
 *  \param[in] i index into SphP for gas cell to consider.
 *  \param[in, out] ne_guess pointer to electron number density relative to
 *                  hydrogen number density (modified).
 *  \param[out] nH0 Pointer to the neutral hydrogen fraction (set to value in
 *              the GasState struct).
 *  \param[out] coolrate Pointer to cooling rate (set to value from
 *              CoolingRateFromU).
 *
 *  \return void
 */
void SetOutputGasState(int i, double *ne_guess, double *nH0, double *coolrate)
{
  double sfr = 0;
  double rho = SphP[i].Density * All.cf_a3inv;
  double u   = dmax(All.MinEgySpec, SphP[i].Utherm);

  /* update GasState as appropriate given compile-time options and cell properties */
#if defined(USE_SFR)
  sfr = get_starformation_rate(i);
#endif /* #if defined(USE_SFR) */

  /* update DoCool */
  DoCool.u_old_input    = u;
  DoCool.rho_input      = rho;
  DoCool.ne_guess_input = *ne_guess;

  /* convert to physical cgs units */
  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  u *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

  /* calculate cooling rate (and so ne_guess and all of gs including nH0, nHeII) */
  *coolrate = CoolingRateFromU(u, rho, ne_guess);

  *nH0 = gs.nH0;
}

/*! \brief  Calculate (heating rate-cooling rate)/n_h^2 in cgs units.
 *
 *  \param[in] logT log10 of gas temperature.
 *  \param[in] rho Gas density.
 *  \param[in, out] nelec Electron number density relative to hydrogen number
 *                  density.
 *
 *  \return (heating rate-cooling rate)/n_h^2.
 */
double CoolingRate(double logT, double rho, double *nelec)
{
  double Lambda, Heat;
  double LambdaExc, LambdaIon, LambdaRec, LambdaFF, LambdaCmptn = 0.0;
  double LambdaExcH0, LambdaExcHep, LambdaIonH0, LambdaIonHe0, LambdaIonHep;
  double LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd;
  double redshift;
  double T;
  double LambdaPrim = 0.0, LambdaMet = 0.0, LambdaDust = 0.0, LambdaMol = 0.0;

  if(logT <= Tmin)
    logT = Tmin + 0.5 * deltaT; /* floor at Tmin */

  gs.nHcgs = gs.XH * rho / PROTONMASS; /* hydrogen number dens in cgs units */

  if(logT < Tmax)
    {
      find_abundances_and_rates(logT, rho, nelec);

      /* Compute cooling and heating rate (cf KWH Table 1) in units of nH**2 */
      T = pow(10.0, logT);

      LambdaExcH0   = gs.bH0 * gs.ne * gs.nH0;
      LambdaExcHep  = gs.bHep * gs.ne * gs.nHep;
      LambdaExc     = LambdaExcH0 + LambdaExcHep; /* excitation */
      LambdaIonH0   = 2.18e-11 * gs.geH0 * gs.ne * gs.nH0;
      LambdaIonHe0  = 3.94e-11 * gs.geHe0 * gs.ne * gs.nHe0;
      LambdaIonHep  = 8.72e-11 * gs.geHep * gs.ne * gs.nHep;
      LambdaIon     = LambdaIonH0 + LambdaIonHe0 + LambdaIonHep; /* ionization */
      LambdaRecHp   = 1.036e-16 * T * gs.ne * (gs.aHp * gs.nHp);
      LambdaRecHep  = 1.036e-16 * T * gs.ne * (gs.aHep * gs.nHep);
      LambdaRecHepp = 1.036e-16 * T * gs.ne * (gs.aHepp * gs.nHepp);
      LambdaRecHepd = 6.526e-11 * gs.ad * gs.ne * gs.nHep;
      LambdaRec     = LambdaRecHp + LambdaRecHep + LambdaRecHepp + LambdaRecHepd;
      LambdaFF      = gs.bff * (gs.nHp + gs.nHep + 4 * gs.nHepp) * gs.ne;
      LambdaPrim    = LambdaExc + LambdaIon + LambdaRec + LambdaFF;

      if(All.ComovingIntegrationOn)
        {
          redshift    = 1 / All.Time - 1;
          LambdaCmptn = 5.65e-36 * gs.ne * (T - 2.73 * (1. + redshift)) * pow(1. + redshift, 4.) / gs.nHcgs;
        }
      else
        LambdaCmptn = 0;

      Lambda = LambdaPrim + LambdaMet + LambdaDust + LambdaCmptn + LambdaMol;

      Heat = 0;
      if(pc.J_UV != 0)
        Heat += (gs.nH0 * pc.epsH0 + gs.nHe0 * pc.epsHe0 + gs.nHep * pc.epsHep) / gs.nHcgs;
    }
  else /* here we're outside of tabulated rates, T>Tmax K */
    {
      /* at high T (fully ionized); only free-free and Compton cooling are present. Assumes no heating. */
      Heat = 0;

      LambdaExcH0 = LambdaExcHep = LambdaIonH0 = LambdaIonHe0 = LambdaIonHep = LambdaRecHp = LambdaRecHep = LambdaRecHepp =
          LambdaRecHepd                                                                                   = 0;

      /* very hot: H and He both fully ionized */
      gs.nHp   = 1.0;
      gs.nHep  = 0;
      gs.nHepp = gs.yhelium;
      gs.ne    = gs.nHp + 2.0 * gs.nHepp;
      *nelec   = gs.ne; /* note: in units of the hydrogen number density */

      T        = pow(10.0, logT);
      LambdaFF = 1.42e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - logT) * (5.5 - logT) / 3)) * (gs.nHp + 4 * gs.nHepp) * gs.ne;

      if(All.ComovingIntegrationOn)
        {
          redshift = 1 / All.Time - 1;
          /* add inverse Compton cooling off the microwave background */
          LambdaCmptn = 5.65e-36 * gs.ne * (T - 2.73 * (1. + redshift)) * pow(1. + redshift, 4.) / gs.nHcgs;
        }
      else
        LambdaCmptn = 0;

      Lambda = LambdaFF + LambdaCmptn;
    }

  return (Heat - Lambda);
}

/*! \brief Make cooling rates interpolation table.
 *
 *  Set up interpolation tables in T for cooling rates given in
 *  KWH, ApJS, 105, 19.
 *
 *  \return void
 */
void MakeRateTable(void)
{
  int i;
  double T;
  double Tfact;

  gs.yhelium = (1 - gs.XH) / (4 * gs.XH);
  gs.mhboltz = PROTONMASS / BOLTZMANN;
  if(All.MinGasTemp > 0.0)
    Tmin = log10(0.1 * All.MinGasTemp);
  else
    Tmin = 1.0;
  deltaT    = (Tmax - Tmin) / NCOOLTAB;
  gs.ethmin = pow(10.0, Tmin) * (1. + gs.yhelium) / ((1. + 4. * gs.yhelium) * gs.mhboltz * GAMMA_MINUS1);
  /* minimum internal energy for neutral gas */

  for(i = 0; i <= NCOOLTAB; i++)
    {
      RateT[i].BetaH0 = RateT[i].BetaHep = RateT[i].Betaff = RateT[i].AlphaHp = RateT[i].AlphaHep = RateT[i].AlphaHepp =
          RateT[i].Alphad = RateT[i].GammaeH0 = RateT[i].GammaeHe0 = RateT[i].GammaeHep = 0;

      T     = pow(10.0, Tmin + deltaT * i);
      Tfact = 1.0 / (1 + sqrt(T / 1.0e5));

      /* collisional excitation */
      /* Cen 1992 */
      if(118348 / T < 70)
        RateT[i].BetaH0 = 7.5e-19 * exp(-118348 / T) * Tfact;
      if(473638 / T < 70)
        RateT[i].BetaHep = 5.54e-17 * pow(T, -0.397) * exp(-473638 / T) * Tfact;

      /* free-free */
      RateT[i].Betaff = 1.43e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - log10(T)) * (5.5 - log10(T)) / 3));

      /* recombination */
      /* Cen 1992 */
      /* Hydrogen II */
      RateT[i].AlphaHp = 8.4e-11 * pow(T / 1000, -0.2) / (1. + pow(T / 1.0e6, 0.7)) / sqrt(T);
      /* Helium II */
      RateT[i].AlphaHep = 1.5e-10 * pow(T, -0.6353);
      /* Helium III */
      RateT[i].AlphaHepp = 4. * RateT[i].AlphaHp;

      /* Cen 1992 */
      /* dielectric recombination */
      if(470000 / T < 70)
        RateT[i].Alphad = 1.9e-3 * pow(T, -1.5) * exp(-470000 / T) * (1. + 0.3 * exp(-94000 / T));

      /* collisional ionization */
      /* Cen 1992 */
      /* Hydrogen */
      if(157809.1 / T < 70)
        RateT[i].GammaeH0 = 5.85e-11 * sqrt(T) * exp(-157809.1 / T) * Tfact;
      /* Helium */
      if(285335.4 / T < 70)
        RateT[i].GammaeHe0 = 2.38e-11 * sqrt(T) * exp(-285335.4 / T) * Tfact;
      /* Hellium II */
      if(631515.0 / T < 70)
        RateT[i].GammaeHep = 5.68e-12 * sqrt(T) * exp(-631515.0 / T) * Tfact;
    }
}

/*! \brief Read table input for ionizing parameters.
 *
 *  \param[in] fname Name of file that contains the tabulated parameters.
 *  \param[in] which Flag used to identify the type of the ionizing background
 *                   (0 = UV background, 1 = AGN background, 2=RADCOOL).
 *
 *  \return void
 */
void ReadIonizeParams(char *fname, int which)
{
  int iter, i;
  FILE *fdcool;
  float dummy;

  if(which == 0)
    {
      NheattabUVB = 0;

      for(iter = 0, i = 0; iter < 2; iter++)
        {
          if(!(fdcool = fopen(fname, "r")))
            terminate("COOLING: cannot read ionization table in file `%s'\n", fname);
          if(iter == 0)
            while(fscanf(fdcool, "%g %g %g %g %g %g %g", &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy) != EOF)
              NheattabUVB++;
          if(iter == 1)
            while(fscanf(fdcool, "%g %g %g %g %g %g %g", &PhotoTUVB[i].variable, &PhotoTUVB[i].gH0, &PhotoTUVB[i].gHe,
                         &PhotoTUVB[i].gHep, &PhotoTUVB[i].eH0, &PhotoTUVB[i].eHe, &PhotoTUVB[i].eHep) != EOF)
              i++;
          fclose(fdcool);

          if(iter == 0)
            {
              PhotoTUVB = (PhotoTable *)mymalloc("PhotoT", NheattabUVB * sizeof(PhotoTable));
              mpi_printf("COOLING: read ionization table with %d entries in file `%s'.\n", NheattabUVB, fname);
            }
        }
      /* ignore zeros at end of treecool file */
      for(i = 0; i < NheattabUVB; ++i)
        if(PhotoTUVB[i].gH0 == 0.0)
          break;

      NheattabUVB = i;
      mpi_printf("COOLING: using %d ionization table entries from file `%s'.\n", NheattabUVB, fname);
    }
}

/*! \brief Set the ionization parameters for the UV background.
 *
 *  \return void
 */
void IonizeParamsUVB(void)
{
  int i, ilow;
  double logz, dzlow, dzhi;
  double redshift;

  if(All.ComovingIntegrationOn)
    redshift = 1 / All.Time - 1;
  else
    {
      redshift = 0.0;
    }

  logz = log10(redshift + 1.0);
  ilow = 0;
  for(i = 0; i < NheattabUVB; i++)
    {
      if(PhotoTUVB[i].variable < logz)
        ilow = i;
      else
        break;
    }

  dzlow = logz - PhotoTUVB[ilow].variable;
  dzhi  = PhotoTUVB[ilow + 1].variable - logz;

  if(NheattabUVB == 0 || logz > PhotoTUVB[NheattabUVB - 1].variable || PhotoTUVB[ilow].gH0 == 0 || PhotoTUVB[ilow + 1].gH0 == 0)
    {
      SetZeroIonization();
      return;
    }
  else
    pc.J_UV = 1;

  pc.gJH0   = pow(10., (dzhi * log10(PhotoTUVB[ilow].gH0) + dzlow * log10(PhotoTUVB[ilow + 1].gH0)) / (dzlow + dzhi));
  pc.gJHe0  = pow(10., (dzhi * log10(PhotoTUVB[ilow].gHe) + dzlow * log10(PhotoTUVB[ilow + 1].gHe)) / (dzlow + dzhi));
  pc.gJHep  = pow(10., (dzhi * log10(PhotoTUVB[ilow].gHep) + dzlow * log10(PhotoTUVB[ilow + 1].gHep)) / (dzlow + dzhi));
  pc.epsH0  = pow(10., (dzhi * log10(PhotoTUVB[ilow].eH0) + dzlow * log10(PhotoTUVB[ilow + 1].eH0)) / (dzlow + dzhi));
  pc.epsHe0 = pow(10., (dzhi * log10(PhotoTUVB[ilow].eHe) + dzlow * log10(PhotoTUVB[ilow + 1].eHe)) / (dzlow + dzhi));
  pc.epsHep = pow(10., (dzhi * log10(PhotoTUVB[ilow].eHep) + dzlow * log10(PhotoTUVB[ilow + 1].eHep)) / (dzlow + dzhi));

  return;
}

/*! \brief Reset the ionization parameters.
 *
 *  \return void
 */
void SetZeroIonization(void) { memset(&pc, 0, sizeof(PhotoCurrent)); }

/*! \brief Wrapper function to set the ionizing background.
 *
 *  \return void
 */
void IonizeParams(void) { IonizeParamsUVB(); }

/*! \brief Initialize the cooling module.
 *
 *  This function initializes the cooling module. In particular,
 *  it allocates the memory for the cooling rate and ionization tables
 *  and initializes them.
 *
 *  \return void
 */
void InitCool(void)
{
  /* set default hydrogen mass fraction */
  gs.XH = HYDROGEN_MASSFRAC;

  /* zero photo-ionization/heating rates */
  SetZeroIonization();

  /* allocate and construct rate table */
  RateT = (RateTable *)mymalloc("RateT", (NCOOLTAB + 1) * sizeof(RateTable));
  ;
  MakeRateTable();

  /* read photo tables */
  ReadIonizeParams(All.TreecoolFile, 0);

  mpi_printf("GFM_COOLING: time, time begin = %le\t%le\n", All.Time, All.TimeBegin);
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();

  IonizeParams();
}

/*! \brief Apply the isochoric cooling to all the active gas cells.
 *
 *  \return void
 */
void cooling_only(void) /* normal cooling routine when star formation is disabled */
{
  int idx, i;

  CPU_Step[CPU_MISC] += measure_time();

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i >= 0)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue; /* skip cells that have been swallowed or eliminated */

          cool_cell(i);
        }
    }
  CPU_Step[CPU_COOLINGSFR] += measure_time();
}

/*! \brief Apply the isochoric cooling to a given gas cell.
 *
 *  This function applies the normal isochoric cooling to a single gas cell.
 *  Once the cooling has been applied according to one of the cooling models
 *  implemented, the internal energy per unit mass, the total energy and the
 *  pressure of the cell are updated.
 *
 *  \param[in] i Index of the gas cell to which cooling is applied.
 *
 *  \return void
 */
void cool_cell(int i)
{
  double dt, dtime, ne = 1;
  double unew, dens, dtcool;

  dens = SphP[i].Density;

  dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

  dtime = All.cf_atime * dt / All.cf_time_hubble_a;

  dtcool = dtime;

  ne         = SphP[i].Ne; /* electron abundance (gives ionization state and mean molecular weight) */
  unew       = DoCooling(dmax(All.MinEgySpec, SphP[i].Utherm), dens * All.cf_a3inv, dtcool, &ne);
  SphP[i].Ne = ne;

  if(unew < 0)
    terminate("invalid temperature: Thistask=%d i=%d unew=%g\n", ThisTask, i, unew);

  double du = unew - SphP[i].Utherm;

  if(unew < All.MinEgySpec)
    du = All.MinEgySpec - SphP[i].Utherm;

  SphP[i].Utherm += du;
  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

#ifdef OUTPUT_COOLHEAT
  if(dtime > 0)
    SphP[i].CoolHeat = du * P[i].Mass / dtime;
#endif /* #ifdef OUTPUT_COOLHEAT */

  set_pressure_of_cell(i);
}

#endif /* #ifdef COOLING */
