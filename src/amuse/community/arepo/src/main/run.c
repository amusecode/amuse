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
 * \file        src/main/run.c
 * \date        05/2018
 * \brief       The main simulation loop.
 * \details     contains functions:
 *                void run(void)
 *                void do_second_order_source_terms_first_half(void)
 *                void do_second_order_source_terms_second_half(void)
 *                void set_non_standard_physics_for_current_time(void)
 *                void calculate_non_standard_physics_with_valid_gravity_tree(void)
 *                void calculate_non_standard_physics_with_valid_gravity_tree_always(void)
 *                void calculate_non_standard_physics_prior_mesh_construction(void)
 *                void calculate_non_standard_physics_end_of_step(void)
 *                int check_for_interruption_of_run(void)
 *                int check_for_interruption_of_run(void)
 *                integertime find_next_outputtime(integertime ti_curr)
 *                void execute_resubmit_command(void)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 06.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <ctype.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"
#include "../mesh/voronoi/voronoi.h"

static void do_second_order_source_terms_first_half(void);
static void do_second_order_source_terms_second_half(void);
static void create_end_file(void);

/*! \brief Contains the main simulation loop that iterates over
 *  single timesteps.
 *
 *  The loop terminates when the cpu-time limit is
 *  reached, when a `stop' file is found in the output directory, or
 *  when the simulation ends because we arrived at TimeMax.
 *
 *  If the simulation is started from initial conditions, a domain
 *  decomposition performed, the gravitational forces are computed and the
 *  Voronoi mesh is constructed.
 *
 *  The main loop is structured as follow:
 *   - find new timesteps: find_timesteps()
 *   - first gravitational half kick: do_gravity_step_first_half()
 *   - gradients are calculated: calculate_gradients()
 *   - vertex velocities are assigned: set_vertex_velocities()
 *   - computation of the hydro flux: compute_interface_fluxes() (first half)
 *   - (de)refinement of hydro cells: do_derefinements_and_refinements()
 *   - drifting particles to next sync point: find_next_sync_point()
 *   (Afterwards the timebins are updated, so different particles might
 *   now be active then before)
 *   - (if needed) a new domain decomposition: domain_Decomposition()
 *   - construction of the Voronoi mesh: create_mesh()
 *   - computation of the hydro flux: compute_interface_fluxes() (second half)
 *   - update of primitive variables: update_primitive_variables()
 *   - computation of gravitational forces: in do_gravity_step_second_half()
 *   - second gravitational half kick: do_gravity_step_second_half()
 *
 *  \return void
 */
void run(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  if(RestartFlag != 1) /* if we have restarted from restart files, no need to do the setup sequence */
    {
      mark_active_timebins();

      output_log_messages();

      set_non_standard_physics_for_current_time();

      ngb_treefree();
      domain_free();
      domain_Decomposition(); /* do domain decomposition if needed */

      ngb_treeallocate();
      ngb_treebuild(NumGas);

      calculate_non_standard_physics_prior_mesh_construction();

      create_mesh();

      mesh_setup_exchange();

      update_primitive_variables();

      calculate_non_standard_physics_end_of_step();

      exchange_primitive_variables();

      calculate_gradients();

      set_vertex_velocities(); /* determine the speed of the mesh-generating vertices */

      ngb_update_velocities(); /* update the neighbor tree with the new vertex and cell velocities */

      do_second_order_source_terms_second_half();

      do_gravity_step_second_half();
    }

#if defined(VORONOI_STATIC_MESH)
  if(RestartFlag == 1)
    {
      int n_hydro_backup   = TimeBinsHydro.NActiveParticles;
      int *time_bin_hydro  = (int *)malloc(NumGas * sizeof(int));
      int *hydro_particles = (int *)malloc(n_hydro_backup * sizeof(int));
      for(int j = 0; j < TimeBinsHydro.NActiveParticles; j++)
        hydro_particles[j] = TimeBinsHydro.ActiveParticleList[j];

      for(int j = 0; j < NumGas; j++)
        {
          time_bin_hydro[j]                   = P[j].TimeBinHydro;
          P[j].TimeBinHydro                   = All.HighestActiveTimeBin;
          TimeBinsHydro.ActiveParticleList[j] = j;
        }
      TimeBinsHydro.NActiveParticles = NumGas;

      create_mesh();
      mesh_setup_exchange();

      for(int j = 0; j < NumGas; j++)
        P[j].TimeBinHydro = time_bin_hydro[j];

      TimeBinsHydro.NActiveParticles = n_hydro_backup;
      for(int j = 0; j < TimeBinsHydro.NActiveParticles; j++)
        TimeBinsHydro.ActiveParticleList[j] = hydro_particles[j];

      free(time_bin_hydro);
      free(hydro_particles);
    }
#endif /* #if defined(VORONOI_STATIC_MESH) */

  while(1) /* main loop */
    {
      if(RestartFlag !=
         1) /* if we are starting from restart files, skip in the first iteration the parts until the restart files were written  */
        {
          compute_statistics();

          flush_everything();

          create_snapshot_if_desired();

          if(All.Ti_Current >= TIMEBASE) /* we reached the final time */
            {
              mpi_printf("\nFinal time=%g reached. Simulation ends.\n", All.TimeMax);

              if(All.Ti_lastoutput != All.Ti_Current) /* make a snapshot at the final time in case none has produced at this time */
                produce_dump(); /* this will be overwritten if All.TimeMax is increased and the run is continued */

              create_end_file();  // create empty file called end in output directory

              break;
            }

          find_timesteps_without_gravity(); /* find-timesteps */

          find_gravity_timesteps_and_do_gravity_step_first_half(); /* gravity half-step for hydrodynamics */
                                                                   /* kicks collisionless particles by half a step */

#if(defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)) && !defined(MESHRELAX)
          update_timesteps_from_gravity();
#endif /* #if (defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)) && !defined(MESHRELAX) \
        */

          do_second_order_source_terms_first_half();

          exchange_primitive_variables();

          /* let's reconstruct gradients for every cell using Green-Gauss gradient estimation */
          calculate_gradients();

          /* determine the speed of the mesh-generating vertices */
          set_vertex_velocities();

          /* update the neighbor tree with the new vertex and cell velocities */
          ngb_update_velocities();

          exchange_primitive_variables_and_gradients();

          /* compute intercell flux with Riemann solver and update the cells with the fluxes */
          compute_interface_fluxes(&Mesh);

#ifdef OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT
#ifndef VORONOI_STATIC_MESH
          free_mesh_structures_not_needed_for_derefinement_refinement();
#endif /* #ifndef VORONOI_STATIC_MESH */
#endif /* #ifdef OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT */

#ifdef REFINEMENT
          do_derefinements_and_refinements();
#endif /* #ifdef REFINEMENT */

          write_cpu_log(); /* output some CPU usage log-info (accounts for everything needed up to completion of the current
                              sync-point) */

          find_next_sync_point(); /* find next synchronization time */

          make_list_of_active_particles();

          output_log_messages(); /* write some info to log-files */

#if !defined(VORONOI_STATIC_MESH)
#ifdef OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT
          free_all_remaining_mesh_structures();
#else  /* #ifdef OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT */
          free_mesh();
#endif /* #ifdef OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT #else */
#endif /* #if !defined(VORONOI_STATIC_MESH) */
          /* Check whether we should write a restart file.
           * Note that at this place we do not need to store the mesh, not the gravity tree.
           */
          if(check_for_interruption_of_run())
            return;
        }
      else
        RestartFlag = 0;

      set_non_standard_physics_for_current_time();

#if defined(VORONOI_STATIC_MESH) && !defined(VORONOI_STATIC_MESH_DO_DOMAIN_DECOMPOSITION) /* may only be used if there is no gravity \
                                                                                           */
#else /* #if defined(VORONOI_STATIC_MESH) && !defined(VORONOI_STATIC_MESH_DO_DOMAIN_DECOMPOSITION) */

      if(All.HighestActiveTimeBin >= All.SmallestTimeBinWithDomainDecomposition) /* only do this for sufficiently large steps */
        {
#ifdef VORONOI_STATIC_MESH
          free_mesh();
#endif /* #ifdef VORONOI_STATIC_MESH */

          ngb_treefree();
          domain_free();

          drift_all_particles();

          domain_Decomposition(); /* do new domain decomposition, will also make a new chained-list of synchronized particles */

          ngb_treeallocate();
          ngb_treebuild(NumGas);

#if defined(VORONOI_STATIC_MESH)
          create_mesh();
          mesh_setup_exchange();
#endif /* #if defined(VORONOI_STATIC_MESH) */
        }
#endif /* #if defined(VORONOI_STATIC_MESH) && !defined(VORONOI_STATIC_MESH_DO_DOMAIN_DECOMPOSITION) #else */

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
      special_particle_update_list();
#endif /* #ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE */

      calculate_non_standard_physics_prior_mesh_construction();

#if !defined(VORONOI_STATIC_MESH)
      create_mesh();
      mesh_setup_exchange();
#endif /* #if !defined(VORONOI_STATIC_MESH) */

      exchange_primitive_variables_and_gradients();

      compute_interface_fluxes(&Mesh);

      update_primitive_variables(); /* these effectively closes off the hydro step */

      /* the masses and positions are updated, let's get new forces and potentials */

      do_second_order_source_terms_second_half();

      do_gravity_step_second_half(); /* this closes off the gravity half-step */

      /* do any extra physics, Strang-split (update both primitive and conserved variables as needed ) */
      calculate_non_standard_physics_end_of_step();
    }

  restart(0); /* write a restart file at final time - can be used to continue simulation beyond final time */

  write_cpu_log(); /* output final cpu measurements */
}

/*! \brief Source terms before hydrodynamics timestep.
 *
 *  \return void
 */
void do_second_order_source_terms_first_half(void)
{
#ifdef MHD
  do_mhd_source_terms_first_half();
#endif /* #ifdef MHD */
}

/* \brief Source terms after hydrodynamics timestep.
 *
 *  If there are multiple source terms, the order of the second half source
 *  terms should be applied inverse to the order of the source terms in
 *  do_second_order_source_terms_first_half().
 *
 *  \return void
 */
void do_second_order_source_terms_second_half(void)
{
#ifdef MHD
  do_mhd_source_terms_second_half();
#endif /* #ifdef MHD */
}

/*! \brief Calls extra modules after drift operator.
 *
 *  This routine is called after the active particles are drifted
 *  to the next syncpoint, but before a new domain decomposition
 *  is performed.
 *
 *  \return void
 */
void set_non_standard_physics_for_current_time(void)
{
#if defined(COOLING)
  IonizeParams(); /* set UV background for the current time */
#endif            /* #if defined(COOLING) */
}

/*! \brief calls extra modules after the gravitational force is recomputed.
 *
 *  Only called if full gravity tree is present.
 *  *** NOTICE *** if HIERARCHICAL_GRAVITY is adopted, this function is carried
 *  out once per synchronization time, with in general only a partial tree that
 *  does not necessarily contain all particles. The latter is the case only for
 *   steps where the highest timesteps are active ("full timesteps").
 *
 *  \return void
 */
void calculate_non_standard_physics_with_valid_gravity_tree(void) {}

/*! \brief Calls extra modules after the gravitational force is recomputed
 *
 *  This is for runs which have the full tree at each time step;
 *  no HIERARCHICAL_GRAVITY
 *
 *  \return void
 */
void calculate_non_standard_physics_with_valid_gravity_tree_always(void) {}

/*! \brief Calls extra modules before the Voronoi mesh is built.
 *
 *  \return void
 */
void calculate_non_standard_physics_prior_mesh_construction(void)
{
#if defined(COOLING) && defined(USE_SFR)
  sfr_create_star_particles();
#endif /* #if defined(COOLING) && defined(USE_SFR) */
}

/*! \brief Calls extra modules at the end of the run loop.
 *
 *  The second gravitational half kick is already applied to the
 *  particles and the voronoi mesh is updated.
 *
 * \return void
 */
void calculate_non_standard_physics_end_of_step(void)
{
#ifdef COOLING
#ifdef USE_SFR
  cooling_and_starformation();
#else  /* #ifdef USE_SFR */
  cooling_only();
#endif /* #ifdef USE_SFR #else */
#endif /* #ifdef COOLING */
}

/*! \brief Checks whether the run must interrupted.
 *
 *  The run is interrupted either if the stop file is present or,
 *  if 85% of the CPU time are up. This routine also handles the
 *  regular writing of restart files. The restart file is also
 *  written if the restart file is present.
 *
 *  \return 1 if the run has to be interrupted, 0 otherwise.
 */
int check_for_interruption_of_run(void)
{
  /* Check whether we need to interrupt the run */
  int stopflag = 0;
  if(ThisTask == 0)
    {
      FILE *fd;
      char stopfname[MAXLEN_PATH];

      sprintf(stopfname, "%sstop", All.OutputDir);
      if((fd = fopen(stopfname, "r"))) /* Is the stop-file present? If yes, interrupt the run. */
        {
          fclose(fd);
          printf("stop-file detected. stopping.\n");
          stopflag = 1;
          unlink(stopfname);
        }

      sprintf(stopfname, "%srestart", All.OutputDir);
      if((fd = fopen(stopfname, "r"))) /* Is the restart-file present? If yes, write a user-requested restart file. */
        {
          fclose(fd);
          printf("restart-file detected. writing restart files.\n");
          stopflag = 3;
          unlink(stopfname);
        }

      if(CPUThisRun > 0.85 * All.TimeLimitCPU) /* are we running out of CPU-time ? If yes, interrupt run. */
        {
          printf("reaching time-limit. stopping.\n");
          stopflag = 2;
        }
    }

  MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(stopflag)
    {
      restart(0); /* write restart file */

      MPI_Barrier(MPI_COMM_WORLD);

      if(stopflag == 3)
        return 0;

      if(stopflag == 2 && ThisTask == 0)
        {
          FILE *fd;
          char contfname[MAXLEN_PATH];
          sprintf(contfname, "%scont", All.OutputDir);
          if((fd = fopen(contfname, "w")))
            fclose(fd);

          if(All.ResubmitOn)
            execute_resubmit_command();
        }
      return 1;
    }

  /* is it time to write a regular restart-file? (for security) */
  if(ThisTask == 0)
    {
      if((CPUThisRun - All.TimeLastRestartFile) >= All.CpuTimeBetRestartFile)
        {
          All.TimeLastRestartFile = CPUThisRun;
          stopflag                = 3;
        }
      else
        stopflag = 0;
    }

  MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(stopflag == 3)
    {
      restart(0); /* write an occasional restart file */
      stopflag = 0;
    }
  return 0;
}

/*! \brief Returns the next output time that is equal or larger than
 *         ti_curr.
 *
 *  \param[in] ti_curr Current simulation time.
 *
 *  \return Next output time.
 */
integertime find_next_outputtime(integertime ti_curr)
{
  int i, iter = 0;
  integertime ti, ti_next;
  double next, time;

  DumpFlagNextSnap = 1;
  ti_next          = -1;

  if(All.OutputListOn)
    {
      for(i = 0; i < All.OutputListLength; i++)
        {
          time = All.OutputListTimes[i];

          if(time >= All.TimeBegin && time <= All.TimeMax)
            {
              if(All.ComovingIntegrationOn)
                ti = (integertime)(log(time / All.TimeBegin) / All.Timebase_interval);
              else
                ti = (integertime)((time - All.TimeBegin) / All.Timebase_interval);

#ifdef PROCESS_TIMES_OF_OUTPUTLIST
              /* first, determine maximum output interval based on All.MaxSizeTimestep */
              integertime timax = (integertime)(All.MaxSizeTimestep / All.Timebase_interval);

              /* make it a power 2 subdivision */
              integertime ti_min = TIMEBASE;
              while(ti_min > timax)
                ti_min >>= 1;
              timax = ti_min;

              double multiplier = ti / ((double)timax);

              /* now round this to the nearest multiple of timax */
              ti = ((integertime)(multiplier + 0.5)) * timax;
#endif /* #ifdef PROCESS_TIMES_OF_OUTPUTLIST */
              if(ti >= ti_curr)
                {
                  if(ti_next == -1)
                    {
                      ti_next          = ti;
                      DumpFlagNextSnap = All.OutputListFlag[i];
                    }

                  if(ti_next > ti)
                    {
                      ti_next          = ti;
                      DumpFlagNextSnap = All.OutputListFlag[i];
                    }
                }
            }
        }
    }
  else
    {
      if(All.ComovingIntegrationOn)
        {
          if(All.TimeBetSnapshot <= 1.0)
            terminate("TimeBetSnapshot > 1.0 required for your simulation.\n");
        }
      else
        {
          if(All.TimeBetSnapshot <= 0.0)
            terminate("TimeBetSnapshot > 0.0 required for your simulation.\n");
        }

      time = All.TimeOfFirstSnapshot;
      iter = 0;

      while(time < All.TimeBegin)
        {
          if(All.ComovingIntegrationOn)
            time *= All.TimeBetSnapshot;
          else
            time += All.TimeBetSnapshot;

          iter++;

          if(iter > 1000000)
            terminate("Can't determine next output time.\n");
        }

      while(time <= All.TimeMax)
        {
          if(All.ComovingIntegrationOn)
            ti = (integertime)(log(time / All.TimeBegin) / All.Timebase_interval);
          else
            ti = (integertime)((time - All.TimeBegin) / All.Timebase_interval);

          if(ti >= ti_curr)
            {
              ti_next = ti;
              break;
            }

          if(All.ComovingIntegrationOn)
            time *= All.TimeBetSnapshot;
          else
            time += All.TimeBetSnapshot;

          iter++;

          if(iter > 1000000)
            terminate("Can't determine next output time.\n");
        }
    }

  if(ti_next == -1)
    {
      ti_next = 2 * TIMEBASE; /* this will prevent any further output */

      mpi_printf("\nRUN: There is no valid time for a further snapshot file.\n");
    }
  else
    {
      if(All.ComovingIntegrationOn)
        next = All.TimeBegin * exp(ti_next * All.Timebase_interval);
      else
        next = All.TimeBegin + ti_next * All.Timebase_interval;

#ifdef TIMESTEP_OUTPUT_LIMIT
      mpi_printf("\nRUN: Limiting timestep to %g to fulfill output frequency", 0.1 * (next - All.Time));
      All.TimestepOutputLimit = 0.1 * (next - All.Time);
#endif /* #ifdef TIMESTEP_OUTPUT_LIMIT */

      mpi_printf("\nRUN: Setting next time for snapshot file to Time_next= %g  (DumpFlag=%d)\n\n", next, DumpFlagNextSnap);
    }

  return ti_next;
}

/*! \brief Creates an empty file called 'end' in the output directory.
 *
 *  The existence of this file can be used e.g. for analysis scripts to
 *  verify that the simulation has run up to its final time and ended without
 *  error. Note that the end-file is completely passive.
 *
 *  \return void
 */
static void create_end_file(void)
{
  FILE *fd;
  char contfname[MAXLEN_PATH];
  sprintf(contfname, "%send", All.OutputDir);
  if((fd = fopen(contfname, "w")))
    fclose(fd);
}

/*! \brief Executes the resubmit command.
 *
 *  \return void
 */
void execute_resubmit_command(void)
{
  char buf[1000];
  sprintf(buf, "%s", All.ResubmitCommand);
#ifndef NOCALLSOFSYSTEM
  system(buf);
#endif /* #ifndef NOCALLSOFSYSTEM */
}
