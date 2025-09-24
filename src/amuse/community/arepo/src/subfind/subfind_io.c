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
 * \file        src/subfind/subfind_io.c
 * \date        05/2018
 * \brief       Main output routine for subfind.
 * \details     contains functions:
 *                void subfind_save_final(int num)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 14.05.2018 Prepared file for public release -- Rainer Weinberger
 */

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
#include "../fof/fof.h"

#ifdef SUBFIND
#include "subfind.h"

/*! \brief Saves subfind group catalogue to disk.
 *
 *  Note that this routine calls the FoF I/O routines.
 *
 *  \param[in] num Index of this snapshot output.
 *
 *  \return void
 */
void subfind_save_final(int num)
{
  int i, filenr, gr, ngrps, masterTask, lastTask, totsubs;
  char buf[1000];
  double t0, t1;

  /* prepare list of ids with assigned group numbers */
#ifdef FOF_STOREIDS
  fof_subfind_prepare_ID_list();
#endif /* #ifdef FOF_STOREIDS */

  t0 = second();

  /* fill in the FirstSub-values */
  for(i = 0, totsubs = 0; i < Ngroups; i++)
    {
      if(i > 0)
        Group[i].FirstSub = Group[i - 1].FirstSub + Group[i - 1].Nsubs;
      else
        Group[i].FirstSub = 0;
      totsubs += Group[i].Nsubs;
    }

  MPI_Allgather(&totsubs, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 1, Send_offset[0] = 0; i < NTask; i++)
    Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];

  for(i = 0; i < Ngroups; i++)
    {
      if(Group[i].Nsubs > 0)
        Group[i].FirstSub += Send_offset[ThisTask];
      else
        Group[i].FirstSub = -1;
    }

  CommBuffer = mymalloc("CommBuffer", COMMBUFFERSIZE);

  if(NTask < All.NumFilesPerSnapshot)
    {
      warn(
          "Number of processors must be larger or equal than All.NumFilesPerSnapshot! Reducing All.NumFilesPerSnapshot "
          "accordingly.\n");
      All.NumFilesPerSnapshot = NTask;
    }

  if(All.SnapFormat < 1 || All.SnapFormat > 3)
    mpi_printf("Unsupported File-Format All.SnapFormat=%d \n", All.SnapFormat);

#ifndef HAVE_HDF5
  if(All.SnapFormat == 3)
    {
      mpi_terminate("Code wasn't compiled with HDF5 support enabled!\n");
    }
#endif /* #ifndef  HAVE_HDF5 */

  /* assign processors to output files */
  distribute_file(All.NumFilesPerSnapshot, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);

  if(All.NumFilesPerSnapshot > 1)
    {
      if(ThisTask == 0)
        {
          sprintf(buf, "%s/groups_%03d", All.OutputDir, num);
          mkdir(buf, 02755);
        }
      MPI_Barrier(MPI_COMM_WORLD);
    }

  if(All.NumFilesPerSnapshot > 1)
    sprintf(buf, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "fof_subhalo_tab", num, filenr);
  else
    sprintf(buf, "%s%s_%03d", All.OutputDir, "fof_subhalo_tab", num);

  ngrps = All.NumFilesPerSnapshot / All.NumFilesWrittenInParallel;
  if((All.NumFilesPerSnapshot % All.NumFilesWrittenInParallel))
    ngrps++;

  for(gr = 0; gr < ngrps; gr++)
    {
      if((filenr / All.NumFilesWrittenInParallel) == gr) /* ok, it's this processor's turn */
        fof_subfind_write_file(buf, masterTask, lastTask);

      MPI_Barrier(MPI_COMM_WORLD);
    }

  myfree(CommBuffer);

#ifdef FOF_STOREIDS
  myfree(ID_list);
#endif /* #ifdef FOF_STOREIDS */

  t1 = second();

  mpi_printf("SUBFIND: Subgroup catalogues saved. took = %g sec\n", timediff(t0, t1));
}

#endif /* #ifdef SUBFIND */
