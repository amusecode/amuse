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
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  /* output a welcome message */
  hello();

  /* initialize CPU-time/Wallclock-time measurement */
  init_cpu_log();

  determine_compute_nodes();

  for(PTask = 0; NTask > (1 << PTask); PTask++)
    ;

  begrun0();

  strcpy(ParameterFile, "param.txt");  /* Removing command line parsing. argv[1] replaced with "param.txt". */
  RestartFlag = 0;

  begrun1(); /* set-up run  */

  char fname[MAXLEN_PATH];
  strcpy(fname, All.InitCondFile);

  /* now we can load the file */

#ifdef READ_DM_AS_GAS
      read_ic(fname, (RestartFlag == 14) ? 0x02 : LOAD_TYPES);
#else  /* #ifdef READ_DM_AS_GAS */
      read_ic(fname, (RestartFlag == 14) ? 0x01 : LOAD_TYPES);
#endif /* #ifdef READ_DM_AS_GAS #else */

  /* init returns a status code, where a value of >=0 means that endrun() should be called. */
  int status = init();

  if(status >= 0)
    {
      if(status > 0)
        mpi_printf("init() returned with %d\n", status);

      endrun();
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

  MPI_Finalize();
  exit(0);
}
