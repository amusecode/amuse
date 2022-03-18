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
 * \file        src/init/begrun.c
 * \date        05/2018
 * \brief       Initial set-up of a simulation run
 * \details     This file contains various functions to initialize a simulation
 *              run. In particular, the parameter file is read in and parsed
 *              and global variables are initialized to their proper values.
 *              contains functions:
 *                void hello(void)
 *                void begrun0(void)
 *                void begrun1(void)
 *                void begrun2(void)
 *                void set_units(void)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 03.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"
#include "../mesh/voronoi/voronoi.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
herr_t my_hdf5_error_handler(void *unused);
#endif

static void delete_end_file(void);

/*! \brief Prints a welcome message.
 *
 *  \return void
 */
void hello(void)
{
  mpi_printf(
      "\n   __    ____  ____  ____  _____\n  /__\\  (  _ \\( ___)(  _ \\(  _  )\n /(__)\\  )   / )__)  )___/ "
      ")(_)(\n(__)(__)(_)\\_)(____)(__)  (_____)\n\n");
}

/*! \brief Prints used compile options.
 *
 *  \return void
 */
void begrun0(void)
{
  mpi_printf(
      "\nThis is Arepo, version %s.\n\nRunning with %d MPI tasks.\n\nApparently we're using %d compute nodes (we have a minimum of %d "
      "MPI tasks per node, and a maximum of %d)\n\nCode was compiled with settings:\n\n",
      AREPO_VERSION, NTask, NumNodes, MinTasksPerNode, MaxTasksPerNode);

  if(ThisTask == 0)
    {
      output_compile_time_options();
    }
}

/*! \brief Initial setup of the simulation.
 *
 *  First, the parameter file is read by read_parameter_file(),
 *  then routines for setting units, etc are called. This function only does
 *  the setup necessary to load the IC file. After the IC file has been loaded
 *  and prepared by init(), setup continues with begrun2(). This splitting is
 *  done so that we can return cleanly from operations that don't actually
 *  start the simulation (converting snapshots, making projected images, etc.)
 *
 * \return void
 */
void begrun1(void)
{
  read_parameter_file(ParameterFile); /* ... read in parameters for this run */

  check_parameters(); /* consistency check of parameters */

#ifdef HAVE_HDF5
  H5Eset_auto(my_hdf5_error_handler, NULL);
#endif /* #ifdef HAVE_HDF5 */

  gsl_set_error_handler(my_gsl_error_handler);

#ifdef DEBUG
  enable_core_dumps_and_fpu_exceptions();
#endif /* #ifdef DEBUG */

  mpi_printf("BEGRUN: Size of particle structure       %3d  [bytes]\n", (int)sizeof(struct particle_data));
  mpi_printf("BEGRUN: Size of sph particle structure   %3d  [bytes]\n", (int)sizeof(struct sph_particle_data));
  mpi_printf("BEGRUN: Size of gravity tree node        %3d  [bytes]\n", (int)sizeof(struct NODE));
#ifdef MULTIPLE_NODE_SOFTENING
  mpi_printf("BEGRUN: Size of auxiliary gravity node   %3d  [bytes]\n", (int)sizeof(struct ExtNODE));
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

  set_units();

  if(RestartFlag == 1) /* this is needed here to allow domain decomposition right after restart */
    if(All.ComovingIntegrationOn)
      init_drift_table();

  init_io_fields();

  force_short_range_init();

#if defined(FORCETEST) && !defined(FORCETEST_TESTFORCELAW)
  forcetest_ewald_init();
#endif /* #if defined (FORCETEST) && !defined(FORCETEST_TESTFORCELAW) */

  /* set up random number generators */
  random_generator     = gsl_rng_alloc(gsl_rng_ranlxd1);
  random_generator_aux = gsl_rng_alloc(gsl_rng_ranlxd1);

  /* individual start-up seed */
  gsl_rng_set(random_generator, 42 + ThisTask);
  gsl_rng_set(random_generator_aux, 31452 + ThisTask);

  timebins_init(&TimeBinsHydro, "Hydro", &All.MaxPartSph);
  timebins_init(&TimeBinsGravity, "Gravity", &All.MaxPart);

#if defined(COOLING)
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();
  InitCool();
#endif /* #if defined(COOLING) */

#if !defined(PMGRID) && defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) && !defined(ONEDIMS_SPHERICAL)
  ewald_init();
#endif /* #if !defined(PMGRID) && defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) && !defined(ONEDIMS_SPHERICAL) */

#ifdef TILE_ICS
  All.BoxSize *= All.TileICsFactor;
#endif /* #ifdef TILE_ICS */

  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
#endif /* #ifdef LONG_X */
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
#endif /* #ifdef LONG_Y */
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
#endif /* #ifdef LONG_Z */

  EgyInjection = 0;

#ifdef PMGRID
  if((RestartFlag != 3) && (RestartFlag != 6))
    long_range_init();
#endif /* #ifdef PMGRID */

  if(RestartFlag <= 2)
    open_logfiles();

  All.TimeLastRestartFile = CPUThisRun;

#ifdef REDUCE_FLUSH
  All.FlushLast = CPUThisRun;
#endif /* #ifdef REDUCE_FLUSH */

  init_scalars();

  init_gradients();
}

/*! \brief Late setup, after the IC file has been loaded but before run() is
 *  called.
 *
 *  The output files are opened and various modules are initialized. The next
 *  output time is determined by find_next_outputtime() and various timers are
 *  set.
 *
 *  \return void
 */
void begrun2(void)
{
  char contfname[1000];
  sprintf(contfname, "%scont", All.OutputDir);
  unlink(contfname);

  delete_end_file();

  if(RestartFlag > 2)
    open_logfiles();

#if defined(USE_SFR)
  sfr_init();
#endif /* #if defined(USE_SFR) */

#ifdef PMGRID
  long_range_init_regionsize();
#endif /* #ifdef PMGRID */

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
  special_particle_create_list();
#endif /* #ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE */

  if(RestartFlag != 1) /* this needs to be done here because here All.TimeBegin has the correct value */
    if(All.ComovingIntegrationOn)
      init_drift_table();

  {
    if(RestartFlag == 2)
      All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 100);
    else
      All.Ti_nextoutput = find_next_outputtime(All.Ti_Current);
  }

  All.TimeLastRestartFile = CPUThisRun;

#ifdef REDUCE_FLUSH
  All.FlushLast = CPUThisRun;
#endif /* #ifdef REDUCE_FLUSH */

#if defined(FORCETEST) && defined(FORCETEST_TESTFORCELAW)
  gravity_forcetest_testforcelaw();
#endif /* #if defined(FORCETEST) && defined(FORCETEST_TESTFORCELAW) */
}

/*! \brief Computes conversion factors between internal code units and the
 *  cgs-system.
 *
 *  In addition constants like the gravitation constant are set.
 *
 *  \return void
 */
void set_units(void)
{
  double meanweight;

#ifdef STATICNFW
  double Mtot;
#endif /* #ifdef STATICNFW */

  All.UnitTime_in_s         = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;

  All.UnitDensity_in_cgs     = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs    = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs      = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;

  mpi_printf("BEGRUN: Hubble (internal units)   = %g\n", All.Hubble);
  mpi_printf("BEGRUN: G (internal units)        = %g\n", All.G);
  mpi_printf("BEGRUN: UnitMass_in_g             = %g\n", All.UnitMass_in_g);
  mpi_printf("BEGRUN: UnitTime_in_s             = %g\n", All.UnitTime_in_s);
  mpi_printf("BEGRUN: UnitVelocity_in_cm_per_s  = %g\n", All.UnitVelocity_in_cm_per_s);
  mpi_printf("BEGRUN: UnitDensity_in_cgs        = %g\n", All.UnitDensity_in_cgs);
  mpi_printf("BEGRUN: UnitEnergy_in_cgs         = %g\n", All.UnitEnergy_in_cgs);
  mpi_printf("\n");

  meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC); /* note: assuming NEUTRAL GAS */

  if(All.MinEgySpec == 0)
    {
      All.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
      All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

      mpi_printf("BEGRUN: MinEgySpec set to %g based on MinGasTemp=%g\n", All.MinEgySpec, All.MinGasTemp);
    }

#if defined(USE_SFR)
  set_units_sfr();
#endif /* #if defined(USE_SFR) */

#ifdef STATICNFW
  R200    = pow(NFW_M200 * All.G / (100 * All.Hubble * All.Hubble), 1.0 / 3);
  Rs      = R200 / NFW_C;
  Dc      = 200.0 / 3 * NFW_C * NFW_C * NFW_C / (log(1 + NFW_C) - NFW_C / (1 + NFW_C));
  RhoCrit = 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
  V200    = 10 * All.Hubble * R200;
  mpi_printf("V200= %g\n", V200);

  fac  = 1.0;
  Mtot = enclosed_mass(R200);
  mpi_printf("M200= %g\n", Mtot);
  fac  = V200 * V200 * V200 / (10 * All.G * All.Hubble) / Mtot;
  Mtot = enclosed_mass(R200);
  mpi_printf("M200= %g\n", Mtot);
#endif /* #ifdef STATICNFW */
}

/*! \brief deletes the end file if it exists.
 *
 *  This is needed in case a already completed simulation is extended or
 *  overwritten. Note that the end-file is completely passive.
 *
 *  \return void
 */
static void delete_end_file(void)
{
  if(RestartFlag > 2)  // no simulation happening
    {
      return;
    }

  char endfname[1000];
  sprintf(endfname, "%send", All.OutputDir);
  unlink(endfname);
  return;
}
