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
 * \file        src/main/main.c
 * \date        05/2018
 * \brief       Start of the program.
 * \details     contains functions:
 *                int main(int argc, char **argv)
 *                void endrun()
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 06.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif /* #ifdef HAVE_HDF5 */

/*! \brief The entry point of the program.
 *
 *  This function initializes the MPI communication packages, and sets
 *  cpu-time counters to 0. Then begrun1() is called, which sets up
 *  the simulation. Then either IC's or restart files are loaded. In
 *  case of IC's init() is called which prepares the IC's for the run.
 *  A call to begrun2() finishes the initialization. Finally, run() is
 *  started, the main simulation loop, which iterates over the timesteps.
 *
 *  \param[in] argc Argument count from command line.
 *  \param[in] argv Argument vector from command line.
 *
 *  \return status of exit; 0 for normal exit.
 */
int main(int argc, char **argv)
{
#ifdef IMPOSE_PINNING
  detect_topology();
  get_core_set();
#endif /* #ifdef IMPOSE_PINNING */

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  /* output a welcome message */
  hello();

  /* initialize CPU-time/Wallclock-time measurement */
  init_cpu_log();

  determine_compute_nodes();

#ifdef IMPOSE_PINNING
  /* pin the MPI ranks to the available core set */
  pin_to_core_set();
  report_pinning();
#endif /* #ifdef IMPOSE_PINNING */

#ifdef HOST_MEMORY_REPORTING
  mpi_report_committable_memory();
#endif /* #ifdef HOST_MEMORY_REPORTING */

  Argc = argc;
  Argv = argv;

  for(PTask = 0; NTask > (1 << PTask); PTask++)
    ;

  begrun0();

  if(argc < 2)
    {
      if(ThisTask == 0)
        {
          printf("\nParameters are missing. \n");
          printf("Call with <ParameterFile> [<RestartFlag>] [<RestartSnapNum>] [<SpecialOptions>]\n");
          printf("\n");
          printf("   RestartFlag    Action\n");
          printf("       0          Read initial conditions and start simulation\n");
          printf("       1          Read restart files and resume simulation\n");
          printf("       2          Restart from specified snapshot dump and resume simulation\n");
          printf("       3          Run FOF and optionally SUBFIND: [<SubboxSnapNum> for SUBBOX_SNAPSHOTS]\n");
          printf(
              "       6          Convert snapshot file to different format [input=ICFormat  output=SnapFormat   NOTE: derived "
              "quantities have round-off errors!\n");
          printf("      14          Write out the Voronoi mesh: <SnapNum>\n");
          printf("      17          Write out snapshot dump with measured gradients\n");
          printf("      18          Recalculate gravitational potential values for specified snaphot dump: <snapnum>\n");
          printf("\n");
        }
      endrun();
    }

  strcpy(ParameterFile, argv[1]);

  if(argc >= 3)
    RestartFlag = atoi(argv[2]);
  else
    RestartFlag = 0;

  if(argc >= 4)
    RestartSnapNum = atoi(argv[3]);
  else
    RestartSnapNum = -1;

  // Do minimal validation of arguments here rather than in random places in the code
  if((RestartFlag == 3 || RestartFlag == 6 || RestartFlag == 14 || RestartFlag == 17 || RestartFlag == 18) && RestartSnapNum < 0)
    {
      mpi_printf("Need to give the snapshot number\n");
      return (0);
    }

#ifndef RECOMPUTE_POTENTIAL_IN_SNAPSHOT
  if(RestartFlag == 18)
    {
      mpi_printf("Need RECOMPUTE_POTENTIAL_IN_SNAPSHOT for this option\n");
      return (0);
    }
#endif /* #ifndef RECOMPUTE_POTENTIAL_IN_SNAPSHOT */

#ifdef RUNNING_SAFETY_FILE
  /* do not run if 'running' safety file exists */
  int runningflag = 0;
  if(ThisTask == 0)
    {
      FILE *fd;
      char runningfname[MAXLEN_PATH];

      sprintf(runningfname, "./running");
      if((fd = fopen(runningfname, "r"))) /* Is the running-file present? If yes, interrupt the run. */
        {
          fclose(fd);
          printf("running-file detected. stopping.\n");
          runningflag = 1;
        }
    }
  MPI_Bcast(&runningflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(runningflag)
    {
      MPI_Finalize(); /* do not call endrun() */
      return 0;
    }
  else
    {
      /* touch a running safety file */
      if(ThisTask == 0)
        {
          FILE *fd;
          char runningfname[MAXLEN_PATH];

          sprintf(runningfname, "./running");
          if((fd = fopen(runningfname, "w")))
            {
              fclose(fd);
              printf("touching a running-file: %s \n", runningfname);
            }
          else
            terminate("could not touch a running-file: %s\n", runningfname);
        }
    }
#endif /* #ifdef RUNNING_SAFETY_FILE */

  begrun1(); /* set-up run  */

  /* see if we are loading a restart file or an IC file */
  if(RestartFlag == 1)
    loadrestart();
  else
    {
      /* We're reading an IC file. Is it a snapshot or really an IC? */
      char fname[MAXLEN_PATH];

      if(RestartFlag >= 2 && RestartSnapNum >= 0)
        {
          if(All.NumFilesPerSnapshot > 1)
            sprintf(fname, "%s/snapdir_%03d/%s_%03d", All.OutputDir, RestartSnapNum, All.SnapshotFileBase, RestartSnapNum);
          else
            sprintf(fname, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, RestartSnapNum);
        }
      else
        strcpy(fname, All.InitCondFile);

        /* now we can load the file */

#ifdef READ_DM_AS_GAS
      read_ic(fname, (RestartFlag == 14) ? 0x02 : LOAD_TYPES);
#else  /* #ifdef READ_DM_AS_GAS */
      read_ic(fname, (RestartFlag == 14) ? 0x01 : LOAD_TYPES);
#endif /* #ifdef READ_DM_AS_GAS #else */

      /* If we are supposed to just convert the file, write and exit here. */
      if(RestartFlag == 6)
        {
          /* important for proper functioning of FOF+SUBFIND */
          if(All.ComovingIntegrationOn) /* change to new velocity variable */
            {
              int i, j;
              for(i = 0; i < NumPart; i++)
                for(j = 0; j < 3; j++)
                  P[i].Vel[j] *= sqrt(All.Time) * All.Time;
            }
          set_softenings();
          All.TopNodeAllocFactor = 0.08;
          All.TreeAllocFactor    = 0.7;
          All.NgbTreeAllocFactor = 0.7;

          sprintf(All.SnapshotFileBase, "%s_converted", All.SnapshotFileBase);
          mpi_printf("Start writing file %s\nRestartSnapNum %d\n", All.SnapshotFileBase, RestartSnapNum);
          savepositions(RestartSnapNum, 0);
          endrun();
        }

      /* init returns a status code, where a value of >=0 means that endrun() should be called. */
      int status = init();

      if(status >= 0)
        {
          if(status > 0)
            mpi_printf("init() returned with %d\n", status);

          endrun();
        }
    }

  begrun2();

  run(); /* main simulation loop */

  endrun(); /* clean up & finalize MPI */

  return 0;
}

/*! \brief This function ends the simulations in case of no error.
 *
 *  This method has to be called by all processes. It should be used only
 *  if the simulation ends without a errors.
 *  Otherwise terminate() should be used instead.
 *
 *  \return void
 */
void endrun()
{
  mpi_printf("Code run for %f seconds!\n", timediff(StartOfRun, second()));
  mpi_printf("endrun called, calling MPI_Finalize()\nbye!\n\n");
  fflush(stdout);

#ifdef HAVE_HDF5
  /*The hdf5 library will sometimes register an atexit() handler that calls its
   * error handler. In AREPO this is set to my_hdf_error_handler, which calls
   * MPI_Abort. Calling MPI_Abort after MPI_Finalize is not allowed.
   * Hence unset the HDF error handler here
   */
  H5Eset_auto(NULL, NULL);
#endif /* #ifdef HAVE_HDF5 */

#ifdef RUNNING_SAFETY_FILE
  if(All.Ti_Current < TIMEBASE) /* simulation has not reached the final time */
    {
      char running_fname[MAXLEN_PATH], running_done_fname[MAXLEN_PATH];
      sprintf(running_fname, "./running");
      sprintf(running_done_fname, "./running_done");
      rename(running_fname, running_done_fname);
      mpi_printf("moved ./running file to ./running_done, job can now restart.\n");
    }
  else
    mpi_printf("leaving ./running file in place since run is complete to prevent any restarts.\n");
#endif /* #ifdef RUNNING_SAFETY_FILE */

  MPI_Finalize();
  exit(0);
}
