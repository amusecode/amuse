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
 * \file        src/star_formation/sfr_eEOS.c
 * \date        05/2018
 * \brief       Star formation rate routines for the effective multi-phase
 *              model.
 * \details     contains functions:
 *                void cooling_and_starformation(void)
 *                double get_starformation_rate(int i)
 *                void init_clouds(void)
 *                void integrate_sfr(void)
 *                void set_units_sfr(void)
 *                double calc_egyeff(int i, double gasdens, double *ne,
 *                  double *x, double *tsfr, double *factorEVP)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../gravity/forcetree.h"

#ifdef USE_SFR

/*! \brief Main driver for star formation and gas cooling.
 *
 *  This function loops over all the active gas cells. If a given cell
 *  meets the criteria for star formation to be active the multi-phase
 *  model is activated, the properties of the cell are updated according to
 *  the latter and the star formation rate computed. In the other case, the
 *  standard isochoric cooling is applied to the gas cell by calling the
 *  function cool_cell() and the star formation rate is set to 0.
 *
 *  \return void
 */
void cooling_and_starformation(void)
{
  TIMER_START(CPU_COOLINGSFR);

  int idx, i, bin, flag;
  double dt, dtime, ne = 1;
  double unew, du;
  double cloudmass;
  double factorEVP, dens;
  double tsfr;
  double egyeff, x;

  double eos_dens_threshold = All.PhysDensThresh;

  /* note: assuming FULL ionization */
  double u_to_temp_fac =
      (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  /* clear the SFR stored in the active timebins */
  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinSynchronized[bin])
      TimeBinSfr[bin] = 0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Mass == 0 && P[i].ID == 0)
        continue; /* skip cells that have been swallowed or eliminated */

      dens = SphP[i].Density;

      dt    = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      dtime = All.cf_atime * dt / All.cf_time_hubble_a;

      /* apply the temperature floor */

      unew = dmax(All.MinEgySpec, SphP[i].Utherm);

      if(unew < 0)
        terminate("Invalid Temperature: Task=%d i=%d unew=%g\n", ThisTask, i, unew);

      du = unew - SphP[i].Utherm;
      SphP[i].Utherm += du;
      SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

      egyeff = 0.;
      /* calculate the effective equation of state for gas above the density threshold */
      if(dens * All.cf_a3inv >= eos_dens_threshold)
        {
          ne     = SphP[i].Ne;
          egyeff = calc_egyeff(i, dens * All.cf_a3inv, &ne, &x, &tsfr, &factorEVP);
        }

      /* do cooling, except for gas above the EOS density threshold that is colder than the eEOS */
      if(dens * All.cf_a3inv < eos_dens_threshold || (dens * All.cf_a3inv >= eos_dens_threshold && SphP[i].Utherm > egyeff))
        {
          cool_cell(i);
        }

      /* check whether conditions for star formation are fulfilled.
       * f=1  normal cooling
       * f=0  star formation
       */

      flag = 1; /* default is normal cooling */

      /* enable star formation if gas is above SF density threshold */
      if(dens * All.cf_a3inv >= eos_dens_threshold)
        if(SphP[i].Utherm <= egyeff || u_to_temp_fac * SphP[i].Utherm <= All.TemperatureThresh)
          flag = 0;

      if(All.ComovingIntegrationOn)
        if(dens < All.OverDensThresh)
          flag = 1;

      if(P[i].Mass == 0) /* tracer particles don't form stars */
        flag = 1;

      if(flag == 1)
        SphP[i].Sfr = 0;

      /* active star formation */
      if(flag == 0)
        {
          SphP[i].Ne = (HYDROGEN_MASSFRAC + 1) / 2 / HYDROGEN_MASSFRAC; /* note: assuming FULL ionization */

          cloudmass = x * P[i].Mass;

          if(tsfr < dtime)
            tsfr = dtime;

          if(dt > 0)
            {
              if(P[i].TimeBinHydro) /* upon start-up, we need to protect against dt==0 */
                {
                  unew = SphP[i].Utherm;

                  // put (cold) star forming cells on the effective equation of state
                  if(SphP[i].Utherm < egyeff)
                    {
                      unew = egyeff;
                    }

                  du = unew - SphP[i].Utherm;
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
            }

          SphP[i].Sfr = (1 - All.FactorSN) * cloudmass / tsfr * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

          TimeBinSfr[P[i].TimeBinHydro] += SphP[i].Sfr;
        }
    } /* end of main loop over active particles */

  TIMER_STOP(CPU_COOLINGSFR);
}

/*! \brief Return the star formation rate associated with the gas cell i.
 *
 *  \param[in] i the index of the gas cell.
 *
 *  \return star formation rate in solar masses / yr.
 */
double get_starformation_rate(int i)
{
  if(RestartFlag == 3)
    return SphP[i].Sfr;

  double rateOfSF;
  int flag;
  double tsfr;
  double factorEVP, egyeff, ne, x, cloudmass;
  /* note: assuming FULL ionization */
  double u_to_temp_fac =
      (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  double eos_dens_threshold = All.PhysDensThresh;

  flag   = 1; /* default is normal cooling */
  egyeff = 0.0;

  if(SphP[i].Density * All.cf_a3inv >= eos_dens_threshold)
    {
      ne     = SphP[i].Ne;
      egyeff = calc_egyeff(i, SphP[i].Density * All.cf_a3inv, &ne, &x, &tsfr, &factorEVP);
    }

  if(SphP[i].Density * All.cf_a3inv >= All.PhysDensThresh)
    if(SphP[i].Utherm <= 1.01 * egyeff || u_to_temp_fac * SphP[i].Utherm <= All.TemperatureThresh)
      flag = 0;

  if(All.ComovingIntegrationOn)
    if(SphP[i].Density < All.OverDensThresh)
      flag = 1;

  if(flag == 1)
    return 0;

  cloudmass = x * P[i].Mass;

  rateOfSF = (1 - All.FactorSN) * cloudmass / tsfr;

  /* convert to solar masses per yr */
  rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  return rateOfSF;
}

/*! \brief Initialize the parameters of effective multi-phase model.
 *
 *   In particular this function computes the value of PhysDensThresh, that is
 *   the physical density threshold above which star formation is active, if
 *   its value was set to 0 in the parameter file.
 *
 *   \return void
 */
void init_clouds(void)
{
  double A0, dens, tcool, ne, coolrate, egyhot, x, u4, meanweight;
  double tsfr, peff, fac, neff, egyeff, factorEVP, sigma, thresholdStarburst;

  if(All.PhysDensThresh == 0)
    {
      A0 = All.FactorEVP;

      egyhot = All.EgySpecSN / A0;

      meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)); /* note: assuming FULL ionization */
      u4         = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
      u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

      /* choose a high reference density to avoid that we pick up a compton cooling contribution */
      if(All.ComovingIntegrationOn)
        dens = 1.0e10 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
      else
        dens = 1.0e10 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

      if(All.ComovingIntegrationOn)
        {
          All.Time = 1.0; /* to be guaranteed to get z=0 rate */
          set_cosmo_factors_for_current_time();
          IonizeParams();
        }

      ne = 1.0;
      SetZeroIonization();

      tcool = GetCoolingTime(egyhot, dens, &ne);

      coolrate = egyhot / tcool / dens;

      x = (egyhot - u4) / (egyhot - All.EgySpecCold);

      All.PhysDensThresh =
          x / pow(1 - x, 2) * (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold) / (All.MaxSfrTimescale * coolrate);

      mpi_printf(
          "USE_SFR: A0=%g   PhysDensThresh=%g (int units) %g h^2 cm^-3   expected fraction of cold gas at threshold=%g   tcool=%g   "
          "dens=%g   egyhot=%g\n",
          A0, All.PhysDensThresh, All.PhysDensThresh / (PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs), x, tcool, dens,
          egyhot);

      dens = All.PhysDensThresh;

      do
        {
          ne     = 0.5;
          egyeff = calc_egyeff(-1, dens, &ne, &x, &tsfr, &factorEVP);
          peff   = GAMMA_MINUS1 * dens * egyeff;

          fac = 1 / (log(dens * 1.025) - log(dens));
          dens *= 1.025;

          neff = -log(peff) * fac;

          ne     = 0.5;
          egyeff = calc_egyeff(-1, dens, &ne, &x, &tsfr, &factorEVP);
          peff   = GAMMA_MINUS1 * dens * egyeff;

          neff += log(peff) * fac;
        }
      while(neff > 4.0 / 3);

      thresholdStarburst = dens;

      mpi_printf("USE_SFR: run-away sets in for dens=%g   dynamic range for quiescent star formation=%g\n", thresholdStarburst,
                 thresholdStarburst / All.PhysDensThresh);

      integrate_sfr();

      if(ThisTask == 0)
        {
          sigma = 10.0 / All.Hubble * 1.0e-10 / pow(1.0e-3, 2);

          printf("USE_SFR: isotherm sheet central density=%g   z0=%g\n", M_PI * All.G * sigma * sigma / (2 * GAMMA_MINUS1) / u4,
                 GAMMA_MINUS1 * u4 / (2 * M_PI * All.G * sigma));
          myflush(stdout);
        }

      mpi_printf("USE_SFR: SNII energy=%g [internal units] = %g [erg/M_sun] = %g [1e51 erg/Msun]\n", All.FactorSN * All.EgySpecSN,
                 All.FactorSN * All.EgySpecSN / (1 - All.FactorSN) / (All.UnitMass_in_g / All.UnitEnergy_in_cgs) * SOLAR_MASS,
                 All.FactorSN * All.EgySpecSN / (1 - All.FactorSN) / (All.UnitMass_in_g / All.UnitEnergy_in_cgs) * SOLAR_MASS / 1e51);

      if(All.ComovingIntegrationOn)
        {
          All.Time = All.TimeBegin;
          set_cosmo_factors_for_current_time();
          IonizeParams();
        }
    }
}

/*! \brief Compute the effective equation of state for the gas and
 *         the integrated SFR per unit area.
 *
 *  This function computes the effective equation of state for the gas and
 *  the integrated SFR per unit area. It saves the results into two files:
 *  eos.txt for the equation of state and sfrrate.txt for the integrated SFR.
 *  In the latter case, the SFR is determined by integrating along the vertical
 *  direction the gas density of an infinite self-gravitating isothermal sheet.
 *  The integrated gas density is saved as well, so effectively sfrrate.txt
 *  contains the Kennicutt-Schmidt law of the star formation model.
 *
 *  \return void
 */
void integrate_sfr(void)
{
  double rho0, rho, rho2, q, dz, gam, sigma = 0, sigma_u4, sigmasfr = 0, ne, P1;
  double x = 0, P, P2, x2, tsfr2, factorEVP2, drho, dq;
  double meanweight, u4, tsfr, factorEVP, egyeff, egyeff2;
  FILE *fd;

  double eos_dens_threshold = All.PhysDensThresh;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)); /* note: assuming FULL ionization */
  u4         = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  if(All.ComovingIntegrationOn)
    {
      All.Time = 1.0; /* to be guaranteed to get z=0 rate */
      set_cosmo_factors_for_current_time();
      IonizeParams();
    }

  if(WriteMiscFiles && (ThisTask == 0))
    fd = fopen("eos.txt", "w");
  else
    fd = 0;

  for(rho = eos_dens_threshold; rho <= 1000 * eos_dens_threshold; rho *= 1.1)
    {
      ne     = 1.0;
      egyeff = calc_egyeff(-1, rho, &ne, &x, &tsfr, &factorEVP);

      P = GAMMA_MINUS1 * rho * egyeff;

      if(WriteMiscFiles && (ThisTask == 0))
        {
          fprintf(fd, "%g %g %g\n", rho, P, x);
        }
    }

  if(WriteMiscFiles && (ThisTask == 0))
    fclose(fd);

  if(WriteMiscFiles && (ThisTask == 0))
    fd = fopen("sfrrate.txt", "w");
  else
    fd = 0;

  for(rho0 = eos_dens_threshold; rho0 <= 10000 * eos_dens_threshold; rho0 *= 1.02)
    {
      rho = rho0;
      q   = 0;
      dz  = 0.001;

      sigma = sigmasfr = sigma_u4 = 0;

      while(rho > 0.0001 * rho0)
        {
          if(rho > All.PhysDensThresh)
            {
              ne     = 1.0;
              egyeff = calc_egyeff(-1, rho, &ne, &x, &tsfr, &factorEVP);

              P = P1 = GAMMA_MINUS1 * rho * egyeff;

              rho2 = 1.1 * rho;

              egyeff2 = calc_egyeff(-1, rho2, &ne, &x2, &tsfr2, &factorEVP2);

              P2 = GAMMA_MINUS1 * rho2 * egyeff2;

              gam = log(P2 / P1) / log(rho2 / rho);
            }
          else
            {
              tsfr = 0;

              P   = GAMMA_MINUS1 * rho * u4;
              gam = 1.0;

              sigma_u4 += rho * dz;
            }

          drho = q;
          dq   = -(gam - 2) / rho * q * q - 4 * M_PI * All.G / (gam * P) * rho * rho * rho;

          sigma += rho * dz;
          if(tsfr > 0)
            {
              sigmasfr += (1 - All.FactorSN) * rho * x / tsfr * dz;
            }

          rho += drho * dz;
          q += dq * dz;
        }

      sigma *= 2; /* to include the other side */
      sigmasfr *= 2;
      sigma_u4 *= 2;

      sigma *= All.HubbleParam * (All.UnitMass_in_g / SOLAR_MASS) * PARSEC * PARSEC / (All.UnitLength_in_cm * All.UnitLength_in_cm);
      sigmasfr *= All.HubbleParam * All.HubbleParam * (All.UnitMass_in_g / SOLAR_MASS) * (SEC_PER_YEAR / All.UnitTime_in_s) * 1.0e6 *
                  PARSEC * PARSEC / (All.UnitLength_in_cm * All.UnitLength_in_cm);
      sigma_u4 *= All.HubbleParam * (All.UnitMass_in_g / SOLAR_MASS) * PARSEC * PARSEC / (All.UnitLength_in_cm * All.UnitLength_in_cm);

      if(WriteMiscFiles && (ThisTask == 0))
        {
          fprintf(fd, "%g %g %g %g\n", rho0, sigma, sigmasfr, sigma_u4);
        }
    }

  if(All.ComovingIntegrationOn)
    {
      All.Time = All.TimeBegin;
      set_cosmo_factors_for_current_time();
      IonizeParams();
    }

  if(WriteMiscFiles && (ThisTask == 0))
    fclose(fd);
}

/*! \brief Set the appropriate units for the parameters of the multi-phase
 *         model.
 *
 *  \return void
 */
void set_units_sfr(void)
{
  double meanweight;

  All.OverDensThresh = All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs;

  meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC); /* note: assuming NEUTRAL GAS */

  All.EgySpecCold = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempClouds;
  All.EgySpecCold *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)); /* note: assuming FULL ionization */

  All.EgySpecSN = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempSupernova;
  All.EgySpecSN *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
}

/*! \brief Calculate the effective energy of the multi-phase model.
 *
 *  \param[in] i (unused)
 *  \param[in] gasdens gas density.
 *  \param[in, out] ne Fractional electron density.
 *  \param[out] x Fraction cold gas within model.
 *  \param[out] tsfr Star formation timescale.
 *  \param[out] factorEVP Supernova evaporation factor for given density.
 */
double calc_egyeff(int i, double gasdens, double *ne, double *x, double *tsfr, double *factorEVP)
{
  double egyhot, egyeff, tcool, y;
  double rho = gasdens;

  rho = dmax(rho, All.PhysDensThresh);

  *tsfr = sqrt(All.PhysDensThresh / rho) * All.MaxSfrTimescale;

  *factorEVP = pow(rho / All.PhysDensThresh, -0.8) * All.FactorEVP;

  egyhot = All.EgySpecSN / (1 + *factorEVP) + All.EgySpecCold;

  tcool = GetCoolingTime(egyhot, rho, ne);

  y = *tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);

  *x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

  egyeff = egyhot * (1 - *x) + All.EgySpecCold * (*x);

  return egyeff;
}

#endif /* #ifdef USE_SFR */
