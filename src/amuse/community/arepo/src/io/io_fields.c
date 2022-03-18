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
 * \file        src/io/io_fields.c
 * \date        05/2018
 * \brief       User defined functions for output; needed for all
 *              quantities that are not stored in a global array
 * \details     contains functions:
 *                static void io_func_task(int particle, int components,
 *                  void *out_buffer, int mode)
 *                static void io_func_timebin_hydro(int particle, int
 *                  components, void *out_buffer, int mode)
 *                static void io_func_timestep(int particle, int components,
 *                  void *out_buffer, int mode)
 *                static void io_func_softenings(int particle, int components,
 *                  void *out_buffer, int mode)
 *                void io_func_pos(int particle, int components, void *buffer,
 *                  int mode)
 *                static void io_func_vel(int particle, int components, void
 *                  *buffer, int mode)
 *                static void io_func_coolrate(int particle, int components,
 *                  void *buffer, int mode)
 *                static void io_func_ne(int particle, int components, void
 *                  *buffer, int mode)
 *                static void io_func_nh(int particle, int components, void
 *                  *buffer, int mode)
 *                static void io_func_curlvel(int particle, int components,
 *                  void *out_buffer, int mode)
 *                static void io_func_vorticity(int particle, int components,
 *                  void *out_buffer, int mode)
 *                static void io_func_cell_spin(int particle, int components,
 *                  void *out_buffer, int mode)
 *                static void io_func_bfield(int particle, int components,
 *                  void *out_buffer, int mode)
 *                void init_io_fields()
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 07.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <errno.h>
#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef OUTPUT_TASK
/*! \brief Output of the task the particles are at.
 *
 *  \param[in] particle (unused)
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode (unused)
 *
 *  \return void
 */
static void io_func_task(int particle, int components, void *out_buffer, int mode) { ((int *)out_buffer)[0] = ThisTask; }
#endif /* #ifdef OUTPUT_TASK */

#ifdef OUTPUT_TIMEBIN_HYDRO
/*! \brief Output function of the timebin corresponding to the hydrodynamic
 *         timestep.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode (unused)
 *
 *  \return void
 */
static void io_func_timebin_hydro(int particle, int components, void *out_buffer, int mode)
{
  ((int *)out_buffer)[0] = P[particle].TimeBinHydro;
}
#endif /* #ifdef OUTPUT_TIMEBIN_HYDRO */

#ifdef OUTPUTTIMESTEP
/*! \brief Output function of the hydrodynamic timestep.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode (unused)
 *
 *  \return void
 */
static void io_func_timestep(int particle, int components, void *out_buffer, int mode)
{
  ((MyOutputFloat *)out_buffer)[0] =
      (P[particle].TimeBinHydro ? (((integertime)1) << P[particle].TimeBinHydro) : 0) * All.Timebase_interval;
}
#endif /* #ifdef OUTPUTTIMESTEP */

#ifdef OUTPUT_SOFTENINGS
/*! \brief Output function of the force softening.
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode Mode (output)
 *
 *  \return void
 */
static void io_func_softenings(int particle, int components, void *out_buffer, int mode)
{
  ((MyOutputFloat *)out_buffer)[0] = All.ForceSoftening[P[particle].SofteningType];
}
#endif /* #ifdef OUTPUT_SOFTENINGS */

/*! \brief IO function of the particle positions.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode Mode (0: output, 1: input).
 *
 *  \return void
 */
void io_func_pos(int particle, int components, void *buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      if(DumpFlag != 3)  // TODO: clean up this code duplication
        {
#ifdef OUTPUT_COORDINATES_IN_DOUBLEPRECISION
          double *pp = buffer;
#else  /* #ifdef OUTPUT_COORDINATES_IN_DOUBLEPRECISION */
          MyOutputFloat *pp = buffer;
#endif /* #ifdef OUTPUT_COORDINATES_IN_DOUBLEPRECISION #else */

          for(k = 0; k < 3; k++)
            {
              pp[k] = P[particle].Pos[k] - All.GlobalDisplacementVector[k];

#if defined(GRAVITY_NOT_PERIODIC)
              if(P[particle].Type != 0)
                continue;
#endif /* #if defined(GRAVITY_NOT_PERIODIC) */
              double boxSize = All.BoxSize;
#ifdef LONG_X
              if(k == 0)
                boxSize = All.BoxSize * LONG_X;
#endif /* #ifdef LONG_X */
#ifdef LONG_Y
              if(k == 1)
                boxSize = All.BoxSize * LONG_Y;
#endif /* #ifdef LONG_Y */
#ifdef LONG_Z
              if(k == 2)
                boxSize = All.BoxSize * LONG_Z;
#endif /* #ifdef LONG_Z */
              while(pp[k] < 0)
                pp[k] += boxSize;
              while(pp[k] >= boxSize)
                pp[k] -= boxSize;
            }
        }
      else
        {
          MyOutputFloat *pp = buffer;

          for(k = 0; k < 3; k++)
            {
              pp[k] = P[particle].Pos[k] - All.GlobalDisplacementVector[k];

#if defined(GRAVITY_NOT_PERIODIC)
              if(P[particle].Type != 0)
                continue;
#endif /* #if defined(GRAVITY_NOT_PERIODIC) */
              double boxSize = All.BoxSize;
#ifdef LONG_X
              if(k == 0)
                boxSize = All.BoxSize * LONG_X;
#endif /* #ifdef LONG_X */
#ifdef LONG_Y
              if(k == 1)
                boxSize = All.BoxSize * LONG_Y;
#endif /* #ifdef LONG_Y */
#ifdef LONG_Z
              if(k == 2)
                boxSize = All.BoxSize * LONG_Z;
#endif /* #ifdef LONG_Z */
              while(pp[k] < 0)
                pp[k] += boxSize;
              while(pp[k] >= boxSize)
                pp[k] -= boxSize;
            }
        }
    }
  else
    {
#ifdef READ_COORDINATES_IN_DOUBLE
      double *in_buffer = buffer;
#else  /* #ifdef READ_COORDINATES_IN_DOUBLE */
      MyInputFloat *in_buffer = buffer;
#endif /* #ifdef READ_COORDINATES_IN_DOUBLE #else */

      for(k = 0; k < components; k++)
        {
          P[particle].Pos[k] = in_buffer[k] + All.GlobalDisplacementVector[k];
        }
    }
}

/*! \brief IO function for velocities.
 *
 *  Note the different factors of scalefactor in the output than in the code!
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components Number of entries in array.
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode Mode 0: output, 1: input.
 *
 *  \return void
 */
static void io_func_vel(int particle, int components, void *buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      for(k = 0; k < components; k++)
        {
          ((MyOutputFloat *)buffer)[k] = P[particle].Vel[k];
          ((MyOutputFloat *)buffer)[k] *= sqrt(All.cf_a3inv); /* we are dealing with p = a^2 * xdot */
        }
    }
  else
    {
      for(k = 0; k < components; k++)
        {
          P[particle].Vel[k] = ((MyInputFloat *)buffer)[k];
        }
    }
}

#ifdef OUTPUTACCELERATION
/*! \brief IO function for gravitational accelerations.
 *
 *  Note different a factors in output than in code.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components Number of entries in array.
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode Mode 0: output, 1: input.
 *
 *  \return void
 */
static void io_func_accel(int particle, int components, void *out_buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      if(RestartFlag != 6)
        for(k = 0; k < 3; k++)
          ((MyOutputFloat *)out_buffer)[k] = All.cf_a2inv * P[particle].GravAccel[k];
      else
        for(k = 0; k < 3; k++)
          ((MyOutputFloat *)out_buffer)[k] = P[particle].GravAccel[k];
#ifdef PMGRID
      if(RestartFlag != 6)
        for(k = 0; k < 3; k++)
          ((MyOutputFloat *)out_buffer)[k] += All.cf_a2inv * P[particle].GravPM[k];
      else
        for(k = 0; k < 3; k++)
          ((MyOutputFloat *)out_buffer)[k] += P[particle].GravPM[k];
#endif /* #ifdef PMGRID */
    }
  else
    {
      for(k = 0; k < 3; k++)
        P[particle].GravAccel[k] = ((MyOutputFloat *)out_buffer)[k];
    }
}
#endif /* #ifdef OUTPUTACCELERATION */

/* -- user defined functions: additional physics -- */
#ifdef OUTPUTCOOLRATE
/*! \brief Output function of cooling rate.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode (unused)
 *
 *  \return void
 */
static void io_func_coolrate(int particle, int components, void *buffer, int mode)
{
  double tcool, ne, nh0, coolrate;

  ne = SphP[particle].Ne;
  SetOutputGasState(particle, &ne, &nh0, &coolrate);

  /* get cooling time */
  tcool = GetCoolingTime(SphP[particle].Utherm, SphP[particle].Density * All.cf_a3inv, &ne);

  /* convert cooling time with current thermal energy to du/dt */
  if(tcool != 0)
    ((MyOutputFloat *)buffer)[0] = SphP[particle].Utherm / tcool;
  else
    ((MyOutputFloat *)buffer)[0] = 0;
}
#endif /* #ifdef OUTPUTCOOLRATE */

/* -- user defined functions: gas properties -- */
#if defined(COOLING)
/*! \brief IO function of the electron number density.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File IO buffer.
 *  \param[in] mode Mode 0: output, 1: input.
 *
 *  \return void
 */
static void io_func_ne(int particle, int components, void *buffer, int mode)
{
  if(mode == 0)
    {
      // normal code path: calculate Ne accounting for GFM options and USE_SFR
      double ne = SphP[particle].Ne;

#if defined(USE_SFR)
      // reproduces previous behavior that Ne is updated prior to output only for Sfr>0 cells
      // if this is unwanted (or redundant) this if() condition should be removed
      double nh0, coolrate;
      if(get_starformation_rate(particle) > 0)
        SetOutputGasState(particle, &ne, &nh0, &coolrate);
#endif /* #if defined(USE_SFR) */

      ((MyOutputFloat *)buffer)[0] = ne;
    }
  else
    {
      SphP[particle].Ne = ((MyInputFloat *)buffer)[0];
    }
}
#endif /* #if defined(COOLING) */

#if defined(COOLING)
/*! \brief Output function for neutral hydrogen fraction.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode (unused)
 *
 *  \return void
 */
static void io_func_nh(int particle, int components, void *buffer, int mode)
{
  double ne, nh0, coolrate;

  ne = SphP[particle].Ne;
  SetOutputGasState(particle, &ne, &nh0, &coolrate);

  ((MyOutputFloat *)buffer)[0] = nh0;
}
#endif /* #if defined(COOLING) */

#ifdef USE_SFR
/*! \brief IO function for star formation rate.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File output buffer.
 *  \param[in] mode Mode 0: output, 1: input.
 *
 *  \return void
 */
static void io_func_sfr(int particle, int components, void *buffer, int mode)
{
  if(mode == 0)
    {
      ((MyOutputFloat *)buffer)[0] = get_starformation_rate(particle);
    }
  else
    {
      SphP[particle].Sfr = ((MyOutputFloat *)buffer)[0];
    }
}
#endif

/* -- user defined functions: other -- */
#if defined(OUTPUT_CURLVEL)
/*! \brief Output function for curl of velocity field.
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File IO buffer.
 *  \param[in] mode Mode 0: output.
 *
 *  \return void
 */
static void io_func_curlvel(int particle, int components, void *out_buffer, int mode)
{
  if(mode == 0)
    {
      ((MyOutputFloat *)out_buffer)[0] = SphP[particle].CurlVel;
    }
}
#endif /* #if defined(OUTPUT_CURLVEL) */

#ifdef OUTPUT_VORTICITY
/*! \brief Output function of vorticity (calculated from velocity spatial
 *         derivatives).
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File IO buffer.
 *  \param[in] mode Mode 0: output
 *
 *  \return void
 */
static void io_func_vorticity(int particle, int components, void *out_buffer, int mode)
{
  if(mode == 0)
    {
      ((MyOutputFloat *)out_buffer)[0] = SphP[particle].Grad.dvel[2][1] - SphP[particle].Grad.dvel[1][2];
      ((MyOutputFloat *)out_buffer)[1] = SphP[particle].Grad.dvel[0][2] - SphP[particle].Grad.dvel[2][0];
      ((MyOutputFloat *)out_buffer)[2] = SphP[particle].Grad.dvel[1][0] - SphP[particle].Grad.dvel[0][1];
    }
}
#endif /* #ifdef OUTPUT_VORTICITY */

#ifdef MHD
/*! \brief IO function for magnetic field.
 *
 *  Note that the output is in Gauss unit system (in code units) while the
 *  internal B-field is in Heaviside-Lorentz system (FACTOR of sqrt(4 PI)!).
 *
 *  \param[in] particle Index of particle/cell.
 *  \param[in] components (unused)
 *  \param[out] out_buffer File IO buffer.
 *  \param[in] mode Mode 0: output, 1: input.
 *
 *  \return void
 */
static void io_func_bfield(int particle, int components, void *out_buffer, int mode)
{
  int k;

  if(mode == 0)
    {
      /* writing: convert from Heavyside-Lorentz to Gauss */
      for(k = 0; k < 3; k++)
        ((MyOutputFloat *)out_buffer)[k] = SphP[particle].B[k] * sqrt(4. * M_PI);
    }
  else
    {
      /* reading: convert from Gauss to Heavyside-Lorentz */
      for(k = 0; k < 3; k++)
        SphP[particle].B[k] = ((MyInputFloat *)out_buffer)[k] / sqrt(4. * M_PI);
    }
}
#endif /* #ifdef MHD */

/*! \brief Function for field registering.
 *
 *  For init_field arguments read the description of init_field.
 *  Don't forget to add the new IO_FLAG to allvars.h.
 *
 *  \return void
 */
void init_io_fields()
{
  /* ALL TYPES */

#ifdef OUTPUT_COORDINATES_IN_DOUBLEPRECISION
  enum types_in_file pos_out = FILE_DOUBLE;
#else  /* #ifdef  OUTPUT_COORDINATES_IN_DOUBLEPRECISION */
  enum types_in_file pos_out = FILE_MY_IO_FLOAT;
#endif /* #ifdef  OUTPUT_COORDINATES_IN_DOUBLEPRECISION #else */
#ifdef READ_COORDINATES_IN_DOUBLE
  enum types_in_file pos_in = FILE_DOUBLE;
#else  /* #ifdef  READ_COORDINATES_IN_DOUBLE */
  enum types_in_file pos_in = FILE_MY_IO_FLOAT;
#endif /* #ifdef  READ_COORDINATES_IN_DOUBLE #else */
  init_field(IO_POS, "POS ", "Coordinates", MEM_MY_DOUBLE, pos_out, pos_in, 3, A_NONE, 0, io_func_pos, ALL_TYPES);
  init_units(IO_POS, 1., -1., 1., 0., 0., All.UnitLength_in_cm);

  init_field(IO_POS_MINI, "POS ", "Coordinates", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_NONE, 3, A_NONE, 0, io_func_pos, ALL_TYPES);
  init_units(IO_POS_MINI, 1., -1., 1., 0., 0., All.UnitLength_in_cm);
  init_snapshot_type(IO_POS_MINI, SN_MINI_ONLY); /* second IO tag output to mini-snaps always in single precision */

  init_field(IO_VEL, "VEL ", "Velocities", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_NONE, 0, io_func_vel,
             ALL_TYPES);                                                 /* particle velocities */
  init_units(IO_VEL, 0.5, 0., 0., 0., 1., All.UnitVelocity_in_cm_per_s); /* sqrt(a)*km/s */
  init_snapshot_type(IO_VEL, SN_MINI);

  init_field(IO_ID, "ID  ", "ParticleIDs", MEM_MY_ID_TYPE, FILE_MY_ID_TYPE, FILE_MY_ID_TYPE, 1, A_P, &P[0].ID, 0, ALL_TYPES);
  init_units(IO_ID, 0, 0, 0, 0, 0, 0);
  init_snapshot_type(IO_ID, SN_MINI);

  init_field(IO_MASS, "MASS", "Masses", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_P, &P[0].Mass, 0,
             SET_IN_GET_PARTICLES_IN_BLOCK); /* particle mass */
  init_units(IO_MASS, 0., -1., 0., 1., 0., All.UnitMass_in_g);
  init_snapshot_type(IO_MASS, SN_MINI);

#ifdef OUTPUTPOTENTIAL
  init_field(IO_POT, "POT ", "Potential", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_P, &P[0].Potential, 0,
             ALL_TYPES); /* gravitational potential */
  init_units(IO_POT, -1.0, 0.0, 0.0, 0.0, 2.0, All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s); /* (km/s)^2/a */

  init_field(IO_POT_MINI, "POT ", "Potential", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_P, &P[0].Potential, 0,
             STARS_ONLY | BHS_ONLY);
  init_units(IO_POT_MINI, -1.0, 0.0, 0.0, 0.0, 2.0, All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
  init_snapshot_type(IO_POT_MINI, SN_MINI_ONLY); /* second IO tag output to mini-snaps for stars/BHs only */
#endif                                           /* #ifdef OUTPUTPOTENTIAL */

  /* GAS CELLS */

  init_field(IO_U, "U   ", "InternalEnergy", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].Utherm, 0,
             GAS_ONLY); /* internal energy */
  init_units(IO_U, 0., 0., 0., 0., 2., All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
  init_snapshot_type(IO_U, SN_MINI);

  init_field(IO_RHO, "RHO ", "Density", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].Density, 0,
             GAS_ONLY); /* particle density */
  init_units(IO_RHO, -3., 2., -3., 1., 0., All.UnitDensity_in_cgs);
  init_snapshot_type(IO_RHO, SN_MINI);

#ifdef OUTPUT_PRESSURE
  init_field(IO_PRESSURE, "PRES", "Pressure", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].Pressure, 0, GAS_ONLY);
  init_units(IO_PRESSURE, -3.0, 2.0, -3.0, 1.0, 2.0,
             All.UnitDensity_in_cgs * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
#endif /* #ifdef OUTPUT_PRESSURE */

#ifdef OUTPUT_CSND
  init_field(IO_CSND, "CSND", "SoundSpeed", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].Csnd, 0, GAS_ONLY);
  init_units(IO_CSND, 0., 0., 0., 0., 1., All.UnitVelocity_in_cm_per_s);
#endif /* #ifdef OUTPUT_CSND */

#if defined(COOLING)
  init_field(IO_NE, "NE  ", "ElectronAbundance", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_NONE, 0, io_func_ne,
             GAS_ONLY);                /* electron abundance */
  init_units(IO_NE, 0, 0, 0, 0, 0, 0); /* dimensionless fraction */
  init_snapshot_type(IO_NE, SN_MINI);

  init_field(IO_NH, "NH  ", "NeutralHydrogenAbundance", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_nh,
             GAS_ONLY);                /* neutral hydrogen fraction */
  init_units(IO_NH, 0, 0, 0, 0, 0, 0); /* dimensionless fraction */
#endif                                 /* #if defined(COOLING) */

#ifdef USE_SFR
  init_field(IO_SFR, "SFR ", "StarFormationRate", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_NONE, 0, io_func_sfr,
             GAS_ONLY);                                                    /* star formation rate */
  init_units(IO_SFR, 0.0, 0.0, -1.0, 1.0, 1.0, SOLAR_MASS / SEC_PER_YEAR); /* Msun/yr */
  init_snapshot_type(IO_SFR, SN_MINI);
#endif /* #ifdef USE_SFR */

#ifdef OUTPUT_DIVVEL
  init_field(IO_DIVVEL, "DIVV", "VelocityDivergence", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].DivVel, 0,
             GAS_ONLY);
  init_units(IO_DIVVEL, 0.0, 1.0, -1.0, 0.0, 1.0, All.UnitVelocity_in_cm_per_s / All.UnitLength_in_cm);
#endif /* #ifdef OUTPUT_DIVVEL */

#if defined(OUTPUT_CURLVEL)
  init_field(IO_CURLVEL, "ROTV", "VelocityCurl", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_NONE, 0, io_func_curlvel,
             GAS_ONLY);
  init_units(IO_CURLVEL, 0.0, 1.0, -1.0, 0.0, 1.0, All.UnitVelocity_in_cm_per_s / All.UnitLength_in_cm);
#endif /* #if defined(OUTPUT_CURLVEL) */

#ifdef OUTPUT_COOLHEAT
  init_field(IO_COOLHEAT, "COHE", "CoolingHeatingEnergy", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].CoolHeat, 0,
             GAS_ONLY);
  init_units(IO_COOLHEAT, 0.0, 0.0, -1.0, 1.0, 3.0, All.UnitEnergy_in_cgs / All.UnitTime_in_s);
#endif /* #ifdef OUTPUT_COOLHEAT */

#ifdef OUTPUT_SURFACE_AREA
  init_field(IO_SAREA, "AREA", "SurfaceArea", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].SurfaceArea, 0,
             GAS_ONLY);
  init_units(IO_SAREA, 2.0, -2.0, 2.0, 0.0, 0.0, All.UnitLength_in_cm * All.UnitLength_in_cm);

  init_field(IO_NFACES, "NFAC", "NumFacesCell", MEM_INT, FILE_INT, FILE_INT, 1, A_SPHP, &SphP[0].CountFaces, 0, GAS_ONLY);
  init_units(IO_NFACES, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
#endif /* #ifdef OUTPUT_SURFACE_AREA */

#ifdef OUTPUTCOOLRATE
  init_field(IO_COOLRATE, "COOR", "CoolingRate", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_coolrate, GAS_ONLY);
  init_units(IO_COOLRATE, 0.0, 0.0, -1.0, 1.0, 3.0, 1.0);
#endif /* #ifdef OUTPUTCOOLRATE */

#ifdef OUTPUT_VORTICITY
  init_field(IO_VORT, "VORT", "Vorticity", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_NONE, 0, io_func_vorticity, GAS_ONLY);
  init_units(IO_VORT, 0.0, 1.0, -1.0, 0.0, 1.0, All.UnitVelocity_in_cm_per_s / All.UnitLength_in_cm);
#endif /* #ifdef OUTPUT_VORTICITY */

  /* GAS CELLS GRADIENTS */

#ifdef OUTPUT_PRESSURE_GRADIENT
  init_field(IO_GRADP, "GRAP", "PressureGradient", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 3, A_SPHP, &SphP[0].Grad.dpress[0], 0,
             GAS_ONLY);
  init_units(IO_GRADP, -4.0, 3.0, -4.0, 1.0, 2.0,
             All.UnitDensity_in_cgs * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s / All.UnitLength_in_cm);
#endif /* #ifdef OUTPUT_PRESSURE_GRADIENT */

#ifdef OUTPUT_DENSITY_GRADIENT
  init_field(IO_GRADR, "GRAR", "DensityGradient", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 3, A_SPHP, &SphP[0].Grad.drho[0], 0,
             GAS_ONLY);
  init_units(IO_GRADR, -4., 3., -4., 1., 0., All.UnitDensity_in_cgs / All.UnitLength_in_cm);
#endif /* #ifdef OUTPUT_DENSITY_GRADIENT */

#ifdef OUTPUT_VELOCITY_GRADIENT
  init_field(IO_GRADV, "GRAV", "VelocityGradient", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 9, A_SPHP, &SphP[0].Grad.dvel[0][0], 0,
             GAS_ONLY);
  init_units(IO_GRADV, 0., 1., -1., 0., 1., All.UnitVelocity_in_cm_per_s / All.UnitLength_in_cm); /* sqrt(a)*km/s */
#endif                                                                                            /* #ifdef OUTPUT_VELOCITY_GRADIENT */

#ifdef OUTPUT_BFIELD_GRADIENT
  init_field(IO_GRADB, "GRAB", "BfieldGradient", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 9, A_SPHP, &SphP[0].Grad.dB[0][0], 0,
             GAS_ONLY);
  init_units(IO_GRADB, -3., 2., -2.5, 0.5, 1., pow(All.UnitPressure_in_cgs, 0.5) / All.UnitLength_in_cm);
#endif /* #ifdef OUTPUT_BFIELD_GRADIENT */

  /* GAS CELLS (MESH PROPERTIES) */

#ifdef OUTPUT_VOLUME
  init_field(IO_VOL, "VOL ", "Volume", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].Volume, 0, GAS_ONLY);
  init_units(IO_VOL, 3., -3., 3., 0., 0., All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm);
#endif /* #ifdef OUTPUT_VOLUME */

#ifdef OUTPUT_VERTEX_VELOCITY
  init_field(IO_VERTEXVEL, "VEVE", "VertexVelocity", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_SPHP,
             &SphP[0].VelVertex[0], 0, GAS_ONLY);
  init_units(IO_VERTEXVEL, 1., 0., 0., 0., 1., All.UnitVelocity_in_cm_per_s);
#endif /* #ifdef OUTPUT_VERTEX_VELOCITY */

#ifdef OUTPUT_MESH_FACE_ANGLE
  init_field(IO_FACEANGLE, "FACA", "MaxFaceAngle", MEM_MY_SINGLE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].MaxFaceAngle, 0,
             GAS_ONLY);
  init_units(IO_FACEANGLE, 0., 0., 0., 0., 0., 0.0);
#endif /* #ifdef OUTPUT_MESH_FACE_ANGLE */

#ifdef OUTPUT_CENTER_OF_MASS
  init_field(IO_CM, "CMCE", "CenterOfMass", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_SPHP, &SphP[0].Center[0], 0,
             GAS_ONLY);
  init_units(IO_CM, 1., -1., 1., 0., 0., All.UnitLength_in_cm);
#endif /* #ifdef OUTPUT_CENTER_OF_MASS */

  /* DIAGNOSTIC */

#ifdef OUTPUT_TASK
  init_field(IO_TASK, "TASK", "task", MEM_INT, FILE_INT, FILE_NONE, 1, A_NONE, 0, io_func_task, GAS_ONLY);
  init_units(IO_TASK, 0., 0., 0., 0., 0., 0.0);
#endif /* #ifdef OUTPUT_TASK */

#ifdef OUTPUT_TIMEBIN_HYDRO
  init_field(IO_TIMEBIN_HYDRO, "TBH", "TimebinHydro", MEM_NONE, FILE_INT, FILE_NONE, 1, A_NONE, 0, io_func_timebin_hydro, GAS_ONLY);
  init_units(IO_TIMEBIN_HYDRO, 0., 0., 0., 0., 0., 0.0);
#endif /* #ifdef OUTPUT_TIMEBIN_HYDRO */

#ifdef OUTPUTTIMESTEP
  init_field(IO_TSTP, "TSTP", "TimeStep", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_NONE, 0, io_func_timestep, ALL_TYPES);
  init_units(IO_TSTP, 0., -1., 1., 0., -1., All.UnitTime_in_s);
#endif /* #ifdef OUTPUTTIMESTEP */

#ifdef OUTPUTACCELERATION
  init_field(IO_ACCEL, "ACCE", "Acceleration", MEM_NONE, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 3, A_NONE, 0, io_func_accel, ALL_TYPES);
  init_units(IO_ACCEL, -1., 1., -1., 0., 2., All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s / All.UnitLength_in_cm);
#endif /* #ifdef OUTPUTACCELERATION */

#ifdef OUTPUT_SOFTENINGS
  init_field(IO_SOFTENING, "SOFT", "Softenings", MEM_NONE, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_NONE, 0, io_func_softenings, ALL_TYPES);
  init_units(IO_SOFTENING, 1., -1., 1., 0., 0., All.UnitLength_in_cm);
#endif /* #ifdef OUTPUT_SOFTENINGS */

#ifdef OUTPUTGRAVINTERACTIONS
  init_field(IO_GRAVITERACTIONS, "GINT", "GravityInteractions", MEM_INT, FILE_INT, FILE_NONE, 1, A_SPHP, &SphP[0].GravInteractions, 0,
             ALL_TYPES);
  init_units(IO_GRAVITERACTIONS, 0., 0., 0., 0., 0., 0.0);
#endif /* #ifdef OUTPUTGRAVINTERACTIONS */

  /* MHD */

#ifdef MHD
  enum types_in_file mhd_read = FILE_MY_IO_FLOAT;
#if defined(MHD_SEEDFIELD)
  if(RestartFlag == 0)
    mhd_read = FILE_NONE; /* magnetic field not expected in ICs */
#endif                    /* #if defined(MHD_SEEDFIELD) */

  init_field(IO_BFLD, "BFLD", "MagneticField", MEM_NONE, FILE_MY_IO_FLOAT, mhd_read, 3, A_NONE, 0, io_func_bfield,
             GAS_ONLY); /* magnetic field  */
  init_units(IO_BFLD, -2., 1., -1.5, 0.5, 1., pow(All.UnitPressure_in_cgs, 0.5));

  init_field(IO_DIVB, "DIVB", "MagneticFieldDivergence", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, 1, A_SPHP, &SphP[0].DivB, 0,
             GAS_ONLY); /* divergence of magnetic field  */
  init_units(IO_DIVB, -3., 2., -2.5, 0.5, 1., pow(All.UnitPressure_in_cgs, 0.5) / All.UnitLength_in_cm);
#endif /* #ifdef MHD */

  /* Scalars */

#ifdef PASSIVE_SCALARS
  init_field(IO_PASS, "PASS", "PassiveScalars", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_MY_IO_FLOAT, PASSIVE_SCALARS, A_SPHP,
             &SphP[0].PScalars[0], 0, GAS_ONLY);
  init_units(IO_PASS, 0., 0., 0., 0., 0., 0.0);
#endif /* #ifdef PASSIVE_SCALARS */

  /* OTHER */

#ifdef SAVE_HSML_IN_SNAPSHOT
  init_field(IO_SUBFINDDENSITY, "SFDE", "SubfindDensity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_PS, &PS[0].SubfindDensity, 0,
             ALL_TYPES);
  init_units(IO_SUBFINDDENSITY, -3., 2., -3., 1., 0., All.UnitDensity_in_cgs);
  init_snapshot_type(IO_SUBFINDDENSITY, SN_NO_SUBBOX);

  init_field(IO_SUBFINDDMDENSITY, "SFDD", "SubfindDMDensity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_PS,
             &PS[0].SubfindDMDensity, 0, ALL_TYPES);
  init_units(IO_SUBFINDDMDENSITY, -3., 2., -3., 1., 0., All.UnitDensity_in_cgs);
  init_snapshot_type(IO_SUBFINDDMDENSITY, SN_NO_SUBBOX);

  init_field(IO_SUBFINDHSML, "SFHS", "SubfindHsml", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_PS, &PS[0].SubfindHsml, 0,
             ALL_TYPES);
  init_units(IO_SUBFINDHSML, 1., -1., 1., 0., 0., All.UnitLength_in_cm);
  init_snapshot_type(IO_SUBFINDHSML, SN_NO_SUBBOX);

  init_field(IO_SUBFINDVELDISP, "SFVD", "SubfindVelDisp", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_PS, &PS[0].SubfindVelDisp, 0,
             ALL_TYPES);
  init_units(IO_SUBFINDVELDISP, 0.0, 0.0, 0.0, 0.0, 1.0, All.UnitVelocity_in_cm_per_s);
  init_snapshot_type(IO_SUBFINDVELDISP, SN_NO_SUBBOX);
#endif /* #ifdef SAVE_HSML_IN_SNAPSHOT */

#if defined(REFINEMENT_HIGH_RES_GAS)
  init_field(IO_HIGHRESMASS, "HRGM", "HighResGasMass", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, FILE_NONE, 1, A_SPHP, &SphP[0].HighResMass, 0,
             GAS_ONLY);
  init_units(IO_HIGHRESMASS, 0, -1, 0, 1, 0, All.UnitMass_in_g);

  init_field(IO_ALLOWREFINEMENT, "REF ", "AllowRefinement", MEM_INT, FILE_INT, FILE_INT, 1, A_SPHP, &SphP[0].AllowRefinement, 0,
             GAS_ONLY);
  init_units(IO_ALLOWREFINEMENT, 0, 0, 0, 0, 0, 0);
#endif /* #if defined(REFINEMENT_HIGH_RES_GAS) */
}
