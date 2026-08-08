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
 * \file        src/fof/fof_io.c
 * \date        05/2018
 * \brief       Output functions for parallel FoF; also used by subfind.
 * \details     contains functions:
 *                void fof_save_groups(int num)
 *                void fof_subfind_prepare_ID_list(void)
 *                void fof_subfind_write_file(char *fname, int writeTask,
 *                  int lastTask)
 *                void fof_subfind_fill_write_buffer(enum fof_subfind_iofields
 *                  blocknr, int *startindex, int pc)
 *                void fof_subfind_get_dataset_name(enum fof_subfind_iofields
 *                  blocknr, char *label)
 *                int fof_subfind_get_dataset_group(enum fof_subfind_iofields
 *                  blocknr)
 *                int fof_subfind_get_particles_in_block(enum
 *                  fof_subfind_iofields blocknr)
 *                int fof_subfind_get_values_per_blockelement(enum
 *                  fof_subfind_iofields blocknr)
 *                int fof_subfind_get_bytes_per_blockelement(enum
 *                  fof_subfind_iofields blocknr)
 *                int fof_subfind_get_datatype(enum fof_subfind_iofields
 *                  blocknr)
 *                int fof_subfind_blockpresent(enum fof_subfind_iofields
 *                  blocknr)
 *                void fof_subfind_get_Tab_IO_Label(enum fof_subfind_iofields
 *                  blocknr, char *label)
 *                void fof_subfind_write_header_attributes_in_hdf5(hid_t
 *                  handle)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 24.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_math.h>
#include <inttypes.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"
#include "../gitversion/version.h"
#include "../subfind/subfind.h"
#include "fof.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
void fof_subfind_write_header_attributes_in_hdf5(hid_t handle);
void write_parameters_attributes_in_hdf5(hid_t handle);
void write_compile_time_options_in_hdf5(hid_t handle);
#endif /* #ifdef HAVE_HDF5 */

#ifdef FOF

/*! \brief Make sure a position lies in the box in case of periodic boundaries.
 *
 *  \param[in] pos Single coordinate in one dimension to be wrapped
 *  \param[in] dim Index of coordinate [0/1/2]
 *
 *  \return double: wrapped coordinate
 */
MyOutputFloat static wrap_position(MyOutputFloat pos, int dim)
{
#if defined(REFLECTIVE_X)
  if(dim == 0)
    return pos;
#endif

#if defined(REFLECTIVE_Y)
  if(dim == 1)
    return pos;
#endif

#if defined(REFLECTIVE_Z)
  if(dim == 2)
    return pos;
#endif

  double boxsize = All.BoxSize;

#ifdef LONG_X
  if(dim == 0)
    boxsize *= LONG_X;
#endif
#ifdef LONG_Y
  if(dim == 1)
    boxsize *= LONG_Y;
#endif
#ifdef LONG_Z
  if(dim == 2)
    boxsize *= LONG_Z;
#endif

  while(pos < 0)
    pos += boxsize;

  while(pos >= boxsize)
    pos -= boxsize;

  return pos;
}

/*! \brief Main routine for group output.
 *
 *  \param[in] num Index of group file (snapshot index for this output).
 *
 *  \return void
 */
void fof_save_groups(int num)
{
  int filenr, gr, ngrps, masterTask, lastTask;
  double t0, t1;
  char buf[500];

#ifdef FOF_STOREIDS
  fof_subfind_prepare_ID_list();
#endif /* #ifdef FOF_STOREIDS */

  t0 = second();

  CommBuffer = mymalloc("CommBuffer", COMMBUFFERSIZE);

  if(NTask < All.NumFilesPerSnapshot)
    {
      warn(
          "Number of processors must be larger or equal than All.NumFilesPerSnapshot! Reducing All.NumFilesPerSnapshot "
          "accordingly.\n");
      All.NumFilesPerSnapshot = NTask;
    }

  if(All.SnapFormat < 1 || All.SnapFormat > 3)
    mpi_printf("Unsupported File-Format. All.SnapFormat=%d\n", All.SnapFormat);

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
    sprintf(buf, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "fof_tab", num, filenr);
  else
    sprintf(buf, "%s%s_%03d", All.OutputDir, "fof_tab", num);

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

  mpi_printf("FOF: Group catalogues saved. took = %g sec\n", timediff(t0, t1));
}

/*! \brief Prepares ID list for option FOF_STOREIDS.
 *
 *  \return void
 */
void fof_subfind_prepare_ID_list(void)
{
  int i, nids;
  long long totNids;
  double t0, t1;

  t0 = second();

  ID_list = mymalloc("ID_list", sizeof(struct id_list) * Nids);

  for(i = 0, nids = 0; i < NumPart; i++)
    {
      if(PS[i].GrNr < TotNgroups)
        {
          if(nids >= Nids)
            terminate("nids >= Nids");

          ID_list[nids].GrNr = PS[i].GrNr;
          ID_list[nids].Type = P[i].Type;
          ID_list[nids].ID   = P[i].ID;
#ifdef SUBFIND
          ID_list[nids].SubNr      = PS[i].SubNr;
          ID_list[nids].BindingEgy = PS[i].BindingEnergy;
#endif /* #ifdef SUBFIND */
          nids++;
        }
    }

  sumup_large_ints(1, &nids, &totNids);
  if(totNids != TotNids)
    {
      char buf[1000];
      sprintf(buf, "Task=%d Nids=%d totNids=%lld TotNids=%lld\n", ThisTask, Nids, totNids, TotNids);
      terminate(buf);
    }

    /* sort the particle IDs according to group-number, and optionally subhalo number and binding energy  */
#ifdef SUBFIND
  parallel_sort(ID_list, Nids, sizeof(struct id_list), subfind_compare_ID_list);
#else  /* #ifdef SUBFIND */
  parallel_sort(ID_list, Nids, sizeof(struct id_list), fof_compare_ID_list_GrNrID);
#endif /* #ifdef SUBFIND #else */

  t1 = second();
  mpi_printf("FOF/SUBFIND: Particle/cell IDs in groups globally sorted. took = %g sec\n", timediff(t0, t1));
}

/*! \brief Writes a file with name fname containing data from writeTask to
 *         lastTask.
 *
 *  \param[in] fname Filename of the output file.
 *  \param[in] writeTask Task responsible for writing the file.
 *  \param[in] lastTask Last task whose data is still in this file.
 *
 *  \return void
 */
void fof_subfind_write_file(char *fname, int writeTask, int lastTask)
{
  int bytes_per_blockelement, npart, nextblock;
  int n_for_this_task, n, p, pc, offset = 0, task;
  int blockmaxlen, n_type[3], ntot_type[3], nn[3];
  enum fof_subfind_iofields blocknr;
  char label[8];
  int bnr;
  int blksize;
  MPI_Status status;
  FILE *fd = 0;
#ifdef HAVE_HDF5
  hid_t hdf5_file = 0, hdf5_grp[3], hdf5_headergrp = 0, hdf5_dataspace_memory;
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  hid_t hdf5_paramsgrp = 0, hdf5_configgrp = 0;
  herr_t hdf5_status;
  hsize_t dims[2], count[2], start[2];
  int rank = 0, pcsum = 0;
  char buf[1000];
#endif /* #ifdef HAVE_HDF5 */

#define SKIP                                 \
  {                                          \
    my_fwrite(&blksize, sizeof(int), 1, fd); \
  }

  /* determine group/id numbers of each type in file */
  n_type[0] = Ngroups;
  n_type[1] = Nsubgroups;
  n_type[2] = Nids;

  if(ThisTask == writeTask)
    {
      for(n = 0; n < 3; n++)
        ntot_type[n] = n_type[n];

      for(task = writeTask + 1; task <= lastTask; task++)
        {
          MPI_Recv(&nn[0], 3, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
          for(n = 0; n < 3; n++)
            ntot_type[n] += nn[n];
        }

      for(task = writeTask + 1; task <= lastTask; task++)
        MPI_Send(&ntot_type[0], 3, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Send(&n_type[0], 3, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
      MPI_Recv(&ntot_type[0], 3, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
    }

  /* fill file header */
  catalogue_header.Ngroups    = ntot_type[0];
  catalogue_header.Nsubgroups = ntot_type[1];
  catalogue_header.Nids       = ntot_type[2];

  catalogue_header.TotNgroups    = TotNgroups;
  catalogue_header.TotNsubgroups = TotNsubgroups;
  catalogue_header.TotNids       = TotNids;

  catalogue_header.num_files = All.NumFilesPerSnapshot;

  catalogue_header.time = All.Time;
  if(All.ComovingIntegrationOn)
    catalogue_header.redshift = 1.0 / All.Time - 1;
  else
    catalogue_header.redshift = 0;
  catalogue_header.HubbleParam = All.HubbleParam;
  catalogue_header.BoxSize     = All.BoxSize;
  catalogue_header.Omega0      = All.Omega0;
  catalogue_header.OmegaLambda = All.OmegaLambda;

#ifdef OUTPUT_IN_DOUBLEPRECISION
  catalogue_header.flag_doubleprecision = 1;
#else  /* #ifdef OUTPUT_IN_DOUBLEPRECISION */
  catalogue_header.flag_doubleprecision = 0;
#endif /* #ifdef OUTPUT_IN_DOUBLEPRECISION #else */

  /* open file and write header */

  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == 3)
        {
#ifdef HAVE_HDF5
          sprintf(buf, "%s.hdf5", fname);
          hdf5_file = my_H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
          mpi_printf("FOF/SUBFIND: writing group catalogue: '%s' (file 1 of %d)\n", fname, All.NumFilesPerSnapshot);
          hdf5_headergrp = my_H5Gcreate(hdf5_file, "/Header", 0);

          hdf5_grp[0] = my_H5Gcreate(hdf5_file, "/Group", 0);
          hdf5_grp[1] = my_H5Gcreate(hdf5_file, "/Subhalo", 0);
          hdf5_grp[2] = my_H5Gcreate(hdf5_file, "/IDs", 0);

          fof_subfind_write_header_attributes_in_hdf5(hdf5_headergrp);

          hdf5_paramsgrp = my_H5Gcreate(hdf5_file, "/Parameters", 0);
          write_parameters_attributes_in_hdf5(hdf5_paramsgrp);

          hdf5_configgrp = my_H5Gcreate(hdf5_file, "/Config", 0);
          write_compile_time_options_in_hdf5(hdf5_configgrp);

#endif /* #ifdef HAVE_HDF5 */
        }
      else
        {
          if(!(fd = fopen(fname, "w")))
            {
              printf("can't open file `%s' for writing snapshot.\n", fname);
              terminate("file open error");
            }

          mpi_printf("FOF/SUBFIND: writing group catalogue: '%s' (file 1 of %d)\n", fname, All.NumFilesPerSnapshot);

          if(All.SnapFormat == 2)
            {
              blksize = sizeof(int) + 4 * sizeof(char);
              SKIP;
              my_fwrite((void *)"HEAD", sizeof(char), 4, fd);
              nextblock = sizeof(catalogue_header) + 2 * sizeof(int);
              my_fwrite(&nextblock, sizeof(int), 1, fd);
              SKIP;
            }

          blksize = sizeof(catalogue_header);

          SKIP;
          my_fwrite(&catalogue_header, sizeof(catalogue_header), 1, fd);
          SKIP;
        }
    }

  for(bnr = 0; bnr < 1000; bnr++)
    {
      blocknr = (enum fof_subfind_iofields)bnr;

      if(blocknr == IO_FOF_LASTENTRY)
        break;

      if(fof_subfind_blockpresent(blocknr))
        {
          bytes_per_blockelement = fof_subfind_get_bytes_per_blockelement(blocknr);

          blockmaxlen = (int)(COMMBUFFERSIZE / bytes_per_blockelement);

          npart   = fof_subfind_get_particles_in_block(blocknr);
          int grp = fof_subfind_get_dataset_group(blocknr);

          if(npart > 0)
            {
              if(ThisTask == 0)
                {
                  char buf[1000];

                  fof_subfind_get_dataset_name(blocknr, buf);
                  printf("FOF/SUBFIND: writing block %d (%s)...\n", blocknr, buf);
                }

              if(ThisTask == writeTask)
                {
                  if(All.SnapFormat == 1 || All.SnapFormat == 2)
                    {
                      if(All.SnapFormat == 2)
                        {
                          blksize = sizeof(int) + 4 * sizeof(char);
                          SKIP;
                          fof_subfind_get_Tab_IO_Label(blocknr, label);
                          my_fwrite(label, sizeof(char), 4, fd);
                          nextblock = npart * bytes_per_blockelement + 2 * sizeof(int);
                          my_fwrite(&nextblock, sizeof(int), 1, fd);
                          SKIP;
                        }

                      blksize = npart * bytes_per_blockelement;
                      SKIP;
                    }
                  else if(All.SnapFormat == 3)
                    {
#ifdef HAVE_HDF5
                      switch(fof_subfind_get_datatype(blocknr))
                        {
                          case 0:
                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_INT);
                            break;
                          case 1:
#ifdef OUTPUT_IN_DOUBLEPRECISION
                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_DOUBLE);
#else  /* #ifdef OUTPUT_IN_DOUBLEPRECISION */
                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_FLOAT);
#endif /* #ifdef OUTPUT_IN_DOUBLEPRECISION #else */
                            break;
                          case 2:
                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_UINT64);
                            break;
                        }

                      dims[0] = ntot_type[grp];
                      dims[1] = fof_subfind_get_values_per_blockelement(blocknr);
                      if(dims[1] == 1)
                        rank = 1;
                      else
                        rank = 2;

                      fof_subfind_get_dataset_name(blocknr, buf);

                      hdf5_dataspace_in_file = my_H5Screate_simple(rank, dims, NULL);

                      hdf5_dataset = my_H5Dcreate(hdf5_grp[grp], buf, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);

                      pcsum = 0;
#endif /* #ifdef HAVE_HDF5 */
                    }
                }

              for(task = writeTask, offset = 0; task <= lastTask; task++)
                {
                  if(task == ThisTask)
                    {
                      n_for_this_task = n_type[grp];

                      for(p = writeTask; p <= lastTask; p++)
                        if(p != ThisTask)
                          MPI_Send(&n_for_this_task, 1, MPI_INT, p, TAG_NFORTHISTASK, MPI_COMM_WORLD);
                    }
                  else
                    MPI_Recv(&n_for_this_task, 1, MPI_INT, task, TAG_NFORTHISTASK, MPI_COMM_WORLD, &status);

                  while(n_for_this_task > 0)
                    {
                      pc = n_for_this_task;

                      if(pc > blockmaxlen)
                        pc = blockmaxlen;

                      if(ThisTask == task)
                        fof_subfind_fill_write_buffer(blocknr, &offset, pc);

                      if(ThisTask == writeTask && task != writeTask)
                        MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task, TAG_PDATA, MPI_COMM_WORLD, &status);

                      if(ThisTask != writeTask && task == ThisTask)
                        MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, writeTask, TAG_PDATA, MPI_COMM_WORLD);

                      if(ThisTask == writeTask)
                        {
                          if(All.SnapFormat == 3)
                            {
#ifdef HAVE_HDF5
                              start[0] = pcsum;
                              start[1] = 0;

                              count[0] = pc;
                              count[1] = fof_subfind_get_values_per_blockelement(blocknr);
                              pcsum += pc;

                              my_H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

                              dims[0]               = pc;
                              dims[1]               = fof_subfind_get_values_per_blockelement(blocknr);
                              hdf5_dataspace_memory = my_H5Screate_simple(rank, dims, NULL);

                              hdf5_status = my_H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory, hdf5_dataspace_in_file,
                                                        H5P_DEFAULT, CommBuffer, buf);

                              (void)hdf5_status;

                              my_H5Sclose(hdf5_dataspace_memory, H5S_SIMPLE);
#endif /* #ifdef HAVE_HDF5 */
                            }
                          else
                            {
                              my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
                            }
                        }

                      n_for_this_task -= pc;
                    }
                }

              if(ThisTask == writeTask)
                {
                  if(All.SnapFormat == 3)
                    {
#ifdef HAVE_HDF5
                      my_H5Dclose(hdf5_dataset, buf);
                      my_H5Sclose(hdf5_dataspace_in_file, H5S_SIMPLE);
                      my_H5Tclose(hdf5_datatype);
#endif /* #ifdef HAVE_HDF5 */
                    }
                  else
                    SKIP;
                }
            }
        }
    }

  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == 3)
        {
#ifdef HAVE_HDF5
          my_H5Gclose(hdf5_grp[0], "/Group");
          my_H5Gclose(hdf5_grp[1], "/Subhalo");
          my_H5Gclose(hdf5_grp[2], "/IDs");
          my_H5Gclose(hdf5_headergrp, "/Header");
          my_H5Gclose(hdf5_paramsgrp, "/Parameters");
          my_H5Gclose(hdf5_configgrp, "/Config");

          my_H5Fclose(hdf5_file, fname);
#endif /* #ifdef HAVE_HDF5 */
        }
      else
        fclose(fd);
    }
}

/*! \brief Copies data from global group array to appropriate output buffer.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *  \param[in] startindex First particle index to be included.
 *  \param[in] pc Particle count; number of particles to be written.
 *
 *  \return void
 */
void fof_subfind_fill_write_buffer(enum fof_subfind_iofields blocknr, int *startindex, int pc)
{
  int n, k, pindex, *ip;
  MyOutputFloat *fp;
  MyIDType *idp;

  fp  = (MyOutputFloat *)CommBuffer;
  ip  = (int *)CommBuffer;
  idp = (MyIDType *)CommBuffer;

  pindex = *startindex;

  for(n = 0; n < pc; pindex++, n++)
    {
      switch(blocknr)
        {
          case IO_FOF_LEN:
            *ip++ = Group[pindex].Len;
            break;
          case IO_FOF_MTOT:
            *fp++ = Group[pindex].Mass;
            break;
          case IO_FOF_POS:
            for(k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = wrap_position(Group[pindex].Pos[k] - All.GlobalDisplacementVector[k], k);
#else  /* #ifdef SUBFIND */
              *fp++ = wrap_position(Group[pindex].CM[k] - All.GlobalDisplacementVector[k], k);
#endif /* #ifdef SUBFIND #else */
            break;
          case IO_FOF_CM:
            for(k = 0; k < 3; k++)
              *fp++ = wrap_position(Group[pindex].CM[k] - All.GlobalDisplacementVector[k], k);
            break;
          case IO_FOF_VEL:
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].Vel[k];
            break;
          case IO_FOF_LENTYPE:
            for(k = 0; k < NTYPES; k++)
              *ip++ = Group[pindex].LenType[k];
            break;
          case IO_FOF_MASSTYPE:
            for(k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].MassType[k];
            break;
          case IO_FOF_SFR:
#ifdef USE_SFR
            *fp++ = Group[pindex].Sfr;
#endif /* #ifdef USE_SFR */
            break;
          case IO_FOF_M_MEAN200:
#ifdef SUBFIND
            *fp++ = Group[pindex].M_Mean200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_R_MEAN200:
#ifdef SUBFIND
            *fp++ = Group[pindex].R_Mean200;
#endif /* #ifdef SUBFIND */
            break;

#ifdef SUBFIND_EXTENDED_PROPERTIES
          case IO_FOF_J_MEAN200:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].J_Mean200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_JDM_MEAN200:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].JDM_Mean200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_JGAS_MEAN200:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].JGas_Mean200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_JSTARS_MEAN200:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].JStars_Mean200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_MASSTYPE_MEAN200:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].MassType_Mean200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_LENTYPE_MEAN200:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *ip++ = Group[pindex].LenType_Mean200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_CMFRAC_MEAN200:
#ifdef SUBFIND
            *fp++ = Group[pindex].CMFrac_Mean200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_CMFRACTYPE_MEAN200:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].CMFracType_Mean200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_J_CRIT200:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].J_Crit200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_JDM_CRIT200:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].JDM_Crit200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_JGAS_CRIT200:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].JGas_Crit200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_JSTARS_CRIT200:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].JStars_Crit200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_MASSTYPE_CRIT200:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].MassType_Crit200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_LENTYPE_CRIT200:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *ip++ = Group[pindex].LenType_Crit200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_CMFRAC_CRIT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].CMFrac_Crit200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_CMFRACTYPE_CRIT200:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].CMFracType_Crit200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_J_CRIT500:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].J_Crit500[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_JDM_CRIT500:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].JDM_Crit500[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_JGAS_CRIT500:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].JGas_Crit500[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_JSTARS_CRIT500:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].JStars_Crit500[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_MASSTYPE_CRIT500:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].MassType_Crit500[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_LENTYPE_CRIT500:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *ip++ = Group[pindex].LenType_Crit500[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_CMFRAC_CRIT500:
#ifdef SUBFIND
            *fp++ = Group[pindex].CMFrac_Crit500;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_CMFRACTYPE_CRIT500:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].CMFracType_Crit500[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_J_TOPHAT200:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].J_TopHat200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_JDM_TOPHAT200:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].JDM_TopHat200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_JGAS_TOPHAT200:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].JGas_TopHat200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_JSTARS_TOPHAT200:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].JStars_TopHat200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_MASSTYPE_TOPHAT200:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].MassType_TopHat200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_LENTYPE_TOPHAT200:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *ip++ = Group[pindex].LenType_TopHat200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_CMFRAC_TOPHAT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].CMFrac_TopHat200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_CMFRACTYPE_TOPHAT200:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].CMFracType_TopHat200[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_EPOT_CRIT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Epot_Crit200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_EKIN_CRIT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ekin_Crit200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_ETHR_CRIT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ethr_Crit200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_EPOT_MEAN200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Epot_Mean200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_EKIN_MEAN200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ekin_Mean200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_ETHR_MEAN200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ethr_Mean200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_EPOT_TOPHAT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Epot_TopHat200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_EKIN_TOPHAT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ekin_TopHat200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_ETHR_TOPHAT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ethr_TopHat200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_EPOT_CRIT500:
#ifdef SUBFIND
            *fp++ = Group[pindex].Epot_Crit500;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_EKIN_CRIT500:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ekin_Crit500;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_ETHR_CRIT500:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ethr_Crit500;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_J:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].J[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_JDM:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].JDM[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_JGAS:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].JGas[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_JSTARS:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = Group[pindex].JStars[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_CMFRAC:
#ifdef SUBFIND
            *fp++ = Group[pindex].CMFrac;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_CMFRACTYPE:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = Group[pindex].CMFracType[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_EKIN:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ekin;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_ETHR:
#ifdef SUBFIND
            *fp++ = Group[pindex].Ethr;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_EPOT:
#ifdef SUBFIND
            *fp++ = Group[pindex].Epot;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_EKIN:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].Ekin;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_ETHR:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].Ethr;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_EPOT:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].Epot;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_J:
            for(k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].J[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_JDM:
            for(k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jdm[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_JGAS:
            for(k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jgas[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_JSTARS:
            for(k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jstars[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_JINHALFRAD:
            for(k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].J_inHalfRad[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_JDMINHALFRAD:
            for(k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jdm_inHalfRad[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_JGASINHALFRAD:
            for(k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jgas_inHalfRad[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_JSTARSINHALFRAD:
            for(k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jstars_inHalfRad[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_JINRAD:
            for(k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].J_inRad[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_JDMINRAD:
            for(k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jdm_inRad[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_JGASINRAD:
            for(k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jgas_inRad[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_JSTARSINRAD:
            for(k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Jstars_inRad[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_CMFRAC:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].CMFrac;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_CMFRACTYPE:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = SubGroup[pindex].CMFracType[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_CMFRACINHALFRAD:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].CMFrac_inHalfRad;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_CMFRACTYPEINHALFRAD:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = SubGroup[pindex].CMFracType_inHalfRad[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_CMFRACINRAD:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].CMFrac_inRad;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_CMFRACTYPEINRAD:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = SubGroup[pindex].CMFracType_inRad[k];
#endif /* #ifdef SUBFIND */
            break;
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

            break;
          case IO_FOF_M_CRIT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].M_Crit200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_R_CRIT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].R_Crit200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_M_CRIT500:
#ifdef SUBFIND
            *fp++ = Group[pindex].M_Crit500;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_R_CRIT500:
#ifdef SUBFIND
            *fp++ = Group[pindex].R_Crit500;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_M_TOPHAT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].M_TopHat200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_R_TOPHAT200:
#ifdef SUBFIND
            *fp++ = Group[pindex].R_TopHat200;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_NSUBS:
#ifdef SUBFIND
            *ip++ = Group[pindex].Nsubs;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_FIRSTSUB:
#ifdef SUBFIND
            *ip++ = Group[pindex].FirstSub;
#endif /* #ifdef SUBFIND */
            break;
          case IO_FOF_FUZZOFFTYPE:
            break;
          case IO_SUB_LEN:
#ifdef SUBFIND
            *ip++ = SubGroup[pindex].Len;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_MTOT:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].Mass;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_POS:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = wrap_position(SubGroup[pindex].Pos[k] - All.GlobalDisplacementVector[k], k);
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_VEL:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = SubGroup[pindex].Vel[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_LENTYPE:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *ip++ = SubGroup[pindex].LenType[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_MASSTYPE:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = SubGroup[pindex].MassType[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_CM:
#ifdef SUBFIND
            for(k = 0; k < 3; k++)
              *fp++ = wrap_position(SubGroup[pindex].CM[k] - All.GlobalDisplacementVector[k], k);
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_SPIN:
            for(k = 0; k < 3; k++)
#ifdef SUBFIND
              *fp++ = SubGroup[pindex].Spin[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_VELDISP:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].SubVelDisp;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_VMAX:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].SubVmax;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_VMAXRAD:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].SubVmaxRad;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_HALFMASSRAD:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].SubHalfMassRad;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_HALFMASSRADTYPE:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = SubGroup[pindex].SubHalfMassRadType[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_MASSINRAD:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].SubMassInRad;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_MASSINRADTYPE:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = SubGroup[pindex].SubMassInRadType[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_MASSINHALFRAD:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].SubMassInHalfRad;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_MASSINHALFRADTYPE:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = SubGroup[pindex].SubMassInHalfRadType[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_MASSINMAXRAD:
#ifdef SUBFIND
            *fp++ = SubGroup[pindex].SubMassInMaxRad;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_MASSINMAXRADTYPE:
#ifdef SUBFIND
            for(k = 0; k < NTYPES; k++)
              *fp++ = SubGroup[pindex].SubMassInMaxRadType[k];
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_IDMOSTBOUND:
#ifdef SUBFIND
            *idp++ = SubGroup[pindex].SubMostBoundID;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_GRNR:
#ifdef SUBFIND
            *ip++ = SubGroup[pindex].GrNr;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_PARENT:
#ifdef SUBFIND
            *ip++ = SubGroup[pindex].SubParent;
#endif /* #ifdef SUBFIND */
            break;
          case IO_SUB_BFLD_HALO:
#if defined(MHD) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].Bfld_Halo * sqrt(4. * M_PI);
#endif /* #if defined(MHD) && defined(SUBFIND) */
            break;
          case IO_SUB_BFLD_DISK:
#if defined(MHD) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].Bfld_Disk * sqrt(4. * M_PI);
#endif /* #if defined(MHD) && defined(SUBFIND) */
            break;
          case IO_SUB_SFR:
#if defined(USE_SFR) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].Sfr;
#endif /* #if defined(USE_SFR) && defined(SUBFIND) */
            break;
          case IO_SUB_SFRINRAD:
#if defined(USE_SFR) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].SfrInRad;
#endif /* #if defined(USE_SFR) && defined(SUBFIND) */
            break;
          case IO_SUB_SFRINHALFRAD:
#if defined(USE_SFR) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].SfrInHalfRad;
#endif /* #if defined(USE_SFR) && defined(SUBFIND) */
            break;
          case IO_SUB_SFRINMAXRAD:
#if defined(USE_SFR) && defined(SUBFIND)
            *fp++ = SubGroup[pindex].SfrInMaxRad;
#endif /* #if defined(USE_SFR) && defined(SUBFIND) */
            break;
          case IO_FOFSUB_IDS:
#ifdef FOF_STOREIDS
            *idp++ = ID_list[pindex].ID;
#endif /* #ifdef FOF_STOREIDS */
            break;

          case IO_FOF_LASTENTRY:
            terminate("should not be reached");
            break;
        }
    }
}

/*! \brief Associates the output variable blocknumber with its name.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *  \param[out] label Name of field.
 *
 *  \return void
 */
void fof_subfind_get_dataset_name(enum fof_subfind_iofields blocknr, char *label)
{
  switch(blocknr)
    {
      case IO_FOF_LEN:
        strcpy(label, "GroupLen");
        break;
      case IO_FOF_MTOT:
        strcpy(label, "GroupMass");
        break;
      case IO_FOF_POS:
        strcpy(label, "GroupPos");
        break;
      case IO_FOF_CM:
        strcpy(label, "GroupCM");
        break;
      case IO_FOF_VEL:
        strcpy(label, "GroupVel");
        break;
      case IO_FOF_LENTYPE:
        strcpy(label, "GroupLenType");
        break;
      case IO_FOF_MASSTYPE:
        strcpy(label, "GroupMassType");
        break;
      case IO_FOF_SFR:
        strcpy(label, "GroupSFR");
        break;
      case IO_FOF_M_MEAN200:
        strcpy(label, "Group_M_Mean200");
        break;
      case IO_FOF_R_MEAN200:
        strcpy(label, "Group_R_Mean200");
        break;

#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_J_MEAN200:
        strcpy(label, "Group_J_Mean200");
        break;
      case IO_FOF_JDM_MEAN200:
        strcpy(label, "Group_Jdm_Mean200");
        break;
      case IO_FOF_JGAS_MEAN200:
        strcpy(label, "Group_Jgas_Mean200");
        break;
      case IO_FOF_JSTARS_MEAN200:
        strcpy(label, "Group_Jstars_Mean200");
        break;
      case IO_FOF_MASSTYPE_MEAN200:
        strcpy(label, "Group_MassType_Mean200");
        break;
      case IO_FOF_LENTYPE_MEAN200:
        strcpy(label, "Group_LenType_Mean200");
        break;
      case IO_FOF_CMFRAC_MEAN200:
        strcpy(label, "Group_CMFrac_Mean200");
        break;
      case IO_FOF_CMFRACTYPE_MEAN200:
        strcpy(label, "Group_CMFracType_Mean200");
        break;
      case IO_FOF_J_CRIT200:
        strcpy(label, "Group_J_Crit200");
        break;
      case IO_FOF_JDM_CRIT200:
        strcpy(label, "Group_Jdm_Crit200");
        break;
      case IO_FOF_JGAS_CRIT200:
        strcpy(label, "Group_Jgas_Crit200");
        break;
      case IO_FOF_JSTARS_CRIT200:
        strcpy(label, "Group_Jstars_Crit200");
        break;
      case IO_FOF_MASSTYPE_CRIT200:
        strcpy(label, "Group_MassType_Crit200");
        break;
      case IO_FOF_LENTYPE_CRIT200:
        strcpy(label, "Group_LenType_Crit200");
        break;
      case IO_FOF_CMFRAC_CRIT200:
        strcpy(label, "Group_CMFrac_Crit200");
        break;
      case IO_FOF_CMFRACTYPE_CRIT200:
        strcpy(label, "Group_CMFracType_Crit200");
        break;
      case IO_FOF_J_CRIT500:
        strcpy(label, "Group_J_Crit500");
        break;
      case IO_FOF_JDM_CRIT500:
        strcpy(label, "Group_Jdm_Crit500");
        break;
      case IO_FOF_JGAS_CRIT500:
        strcpy(label, "Group_Jgas_Crit500");
        break;
      case IO_FOF_JSTARS_CRIT500:
        strcpy(label, "Group_Jstars_Crit500");
        break;
      case IO_FOF_MASSTYPE_CRIT500:
        strcpy(label, "Group_MassType_Crit500");
        break;
      case IO_FOF_LENTYPE_CRIT500:
        strcpy(label, "Group_LenType_Crit500");
        break;
      case IO_FOF_CMFRAC_CRIT500:
        strcpy(label, "Group_CMFrac_Crit500");
        break;
      case IO_FOF_CMFRACTYPE_CRIT500:
        strcpy(label, "Group_CMFracType_Crit500");
        break;
      case IO_FOF_J_TOPHAT200:
        strcpy(label, "Group_J_TopHat200");
        break;
      case IO_FOF_JDM_TOPHAT200:
        strcpy(label, "Group_Jdm_TopHat200");
        break;
      case IO_FOF_JGAS_TOPHAT200:
        strcpy(label, "Group_Jgas_TopHat200");
        break;
      case IO_FOF_JSTARS_TOPHAT200:
        strcpy(label, "Group_Jstars_TopHat200");
        break;
      case IO_FOF_MASSTYPE_TOPHAT200:
        strcpy(label, "Group_MassType_TopHat200");
        break;
      case IO_FOF_LENTYPE_TOPHAT200:
        strcpy(label, "Group_LenType_TopHat200");
        break;
      case IO_FOF_CMFRAC_TOPHAT200:
        strcpy(label, "Group_CMFrac_TopHat200");
        break;
      case IO_FOF_CMFRACTYPE_TOPHAT200:
        strcpy(label, "Group_CMFracType_TopHat200");
        break;
      case IO_FOF_EPOT_CRIT200:
        strcpy(label, "Group_Epot_Crit200");
        break;
      case IO_FOF_EKIN_CRIT200:
        strcpy(label, "Group_Ekin_Crit200");
        break;
      case IO_FOF_ETHR_CRIT200:
        strcpy(label, "Group_Ethr_Crit200");
        break;
      case IO_FOF_EPOT_MEAN200:
        strcpy(label, "Group_Epot_Mean200");
        break;
      case IO_FOF_EKIN_MEAN200:
        strcpy(label, "Group_Ekin_Mean200");
        break;
      case IO_FOF_ETHR_MEAN200:
        strcpy(label, "Group_Ethr_Mean200");
        break;
      case IO_FOF_EPOT_TOPHAT200:
        strcpy(label, "Group_Epot_TopHat200");
        break;
      case IO_FOF_EKIN_TOPHAT200:
        strcpy(label, "Group_Ekin_TopHat200");
        break;
      case IO_FOF_ETHR_TOPHAT200:
        strcpy(label, "Group_Ethr_TopHat200");
        break;
      case IO_FOF_EPOT_CRIT500:
        strcpy(label, "Group_Epot_Crit500");
        break;
      case IO_FOF_EKIN_CRIT500:
        strcpy(label, "Group_Ekin_Crit500");
        break;
      case IO_FOF_ETHR_CRIT500:
        strcpy(label, "Group_Ethr_Crit500");
        break;
      case IO_FOF_J:
        strcpy(label, "Group_J");
        break;
      case IO_FOF_JDM:
        strcpy(label, "Group_Jdm");
        break;
      case IO_FOF_JGAS:
        strcpy(label, "Group_Jgas");
        break;
      case IO_FOF_JSTARS:
        strcpy(label, "Group_Jstars");
        break;
      case IO_FOF_CMFRAC:
        strcpy(label, "Group_CMFrac");
        break;
      case IO_FOF_CMFRACTYPE:
        strcpy(label, "Group_CMFracType");
        break;
      case IO_FOF_EKIN:
        strcpy(label, "GroupEkin");
        break;
      case IO_FOF_ETHR:
        strcpy(label, "GroupEthr");
        break;
      case IO_FOF_EPOT:
        strcpy(label, "GroupEpot");
        break;
      case IO_SUB_EKIN:
        strcpy(label, "SubhaloEkin");
        break;
      case IO_SUB_ETHR:
        strcpy(label, "SubhaloEthr");
        break;
      case IO_SUB_EPOT:
        strcpy(label, "SubhaloEpot");
        break;
      case IO_SUB_J:
        strcpy(label, "Subhalo_J");
        break;
      case IO_SUB_JDM:
        strcpy(label, "Subhalo_Jdm");
        break;
      case IO_SUB_JGAS:
        strcpy(label, "Subhalo_Jgas");
        break;
      case IO_SUB_JSTARS:
        strcpy(label, "Subhalo_Jstars");
        break;
      case IO_SUB_JINHALFRAD:
        strcpy(label, "Subhalo_JInHalfRad");
        break;
      case IO_SUB_JDMINHALFRAD:
        strcpy(label, "Subhalo_JdmInHalfRad");
        break;
      case IO_SUB_JGASINHALFRAD:
        strcpy(label, "Subhalo_JgasInHalfRad");
        break;
      case IO_SUB_JSTARSINHALFRAD:
        strcpy(label, "Subhalo_JstarsInHalfRad");
        break;
      case IO_SUB_JINRAD:
        strcpy(label, "Subhalo_JInRad");
        break;
      case IO_SUB_JDMINRAD:
        strcpy(label, "Subhalo_JdmInRad");
        break;
      case IO_SUB_JGASINRAD:
        strcpy(label, "Subhalo_JgasInRad");
        break;
      case IO_SUB_JSTARSINRAD:
        strcpy(label, "Subhalo_JstarsInRad");
        break;
      case IO_SUB_CMFRAC:
        strcpy(label, "Subhalo_CMFrac");
        break;
      case IO_SUB_CMFRACTYPE:
        strcpy(label, "Subhalo_CMFracType");
        break;
      case IO_SUB_CMFRACINHALFRAD:
        strcpy(label, "Subhalo_CMFracInHalfRad");
        break;
      case IO_SUB_CMFRACTYPEINHALFRAD:
        strcpy(label, "Subhalo_CMFracTypeInHalfRad");
        break;
      case IO_SUB_CMFRACINRAD:
        strcpy(label, "Subhalo_CMFracInRad");
        break;
      case IO_SUB_CMFRACTYPEINRAD:
        strcpy(label, "Subhalo_CMFracTypeInRad");
        break;
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

      case IO_FOF_M_CRIT200:
        strcpy(label, "Group_M_Crit200");
        break;
      case IO_FOF_R_CRIT200:
        strcpy(label, "Group_R_Crit200");
        break;
      case IO_FOF_M_CRIT500:
        strcpy(label, "Group_M_Crit500");
        break;
      case IO_FOF_R_CRIT500:
        strcpy(label, "Group_R_Crit500");
        break;
      case IO_FOF_M_TOPHAT200:
        strcpy(label, "Group_M_TopHat200");
        break;
      case IO_FOF_R_TOPHAT200:
        strcpy(label, "Group_R_TopHat200");
        break;
      case IO_FOF_NSUBS:
        strcpy(label, "GroupNsubs");
        break;
      case IO_FOF_FIRSTSUB:
        strcpy(label, "GroupFirstSub");
        break;
      case IO_FOF_FUZZOFFTYPE:
        strcpy(label, "GroupFuzzOffsetType");
        break;
      case IO_SUB_LEN:
        strcpy(label, "SubhaloLen");
        break;
      case IO_SUB_MTOT:
        strcpy(label, "SubhaloMass");
        break;
      case IO_SUB_POS:
        strcpy(label, "SubhaloPos");
        break;
      case IO_SUB_VEL:
        strcpy(label, "SubhaloVel");
        break;
      case IO_SUB_LENTYPE:
        strcpy(label, "SubhaloLenType");
        break;
      case IO_SUB_MASSTYPE:
        strcpy(label, "SubhaloMassType");
        break;
      case IO_SUB_CM:
        strcpy(label, "SubhaloCM");
        break;
      case IO_SUB_SPIN:
        strcpy(label, "SubhaloSpin");
        break;
      case IO_SUB_VELDISP:
        strcpy(label, "SubhaloVelDisp");
        break;
      case IO_SUB_VMAX:
        strcpy(label, "SubhaloVmax");
        break;
      case IO_SUB_VMAXRAD:
        strcpy(label, "SubhaloVmaxRad");
        break;
      case IO_SUB_HALFMASSRAD:
        strcpy(label, "SubhaloHalfmassRad");
        break;
      case IO_SUB_HALFMASSRADTYPE:
        strcpy(label, "SubhaloHalfmassRadType");
        break;
      case IO_SUB_MASSINRAD:
        strcpy(label, "SubhaloMassInRad");
        break;
      case IO_SUB_MASSINHALFRAD:
        strcpy(label, "SubhaloMassInHalfRad");
        break;
      case IO_SUB_MASSINMAXRAD:
        strcpy(label, "SubhaloMassInMaxRad");
        break;
      case IO_SUB_MASSINRADTYPE:
        strcpy(label, "SubhaloMassInRadType");
        break;
      case IO_SUB_MASSINHALFRADTYPE:
        strcpy(label, "SubhaloMassInHalfRadType");
        break;
      case IO_SUB_MASSINMAXRADTYPE:
        strcpy(label, "SubhaloMassInMaxRadType");
        break;
      case IO_SUB_IDMOSTBOUND:
        strcpy(label, "SubhaloIDMostbound");
        break;
      case IO_SUB_GRNR:
        strcpy(label, "SubhaloGrNr");
        break;
      case IO_SUB_PARENT:
        strcpy(label, "SubhaloParent");
        break;
      case IO_SUB_BFLD_HALO:
        strcpy(label, "SubhaloBfldHalo");
        break;
      case IO_SUB_BFLD_DISK:
        strcpy(label, "SubhaloBfldDisk");
        break;
      case IO_SUB_SFR:
        strcpy(label, "SubhaloSFR");
        break;
      case IO_SUB_SFRINRAD:
        strcpy(label, "SubhaloSFRinRad");
        break;
      case IO_SUB_SFRINHALFRAD:
        strcpy(label, "SubhaloSFRinHalfRad");
        break;
      case IO_SUB_SFRINMAXRAD:
        strcpy(label, "SubhaloSFRinMaxRad");
        break;
      case IO_FOFSUB_IDS:
        strcpy(label, "ID");
        break;

      case IO_FOF_LASTENTRY:
        terminate("should not be reached");
        break;
    }
}

/*! \brief Is this output field a group or subhalo property?
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return 0: group property; 1 subhalo property; 2: both (unused)
 */
int fof_subfind_get_dataset_group(enum fof_subfind_iofields blocknr)
{
  switch(blocknr)
    {
      case IO_FOF_LEN:
      case IO_FOF_MTOT:
      case IO_FOF_POS:
      case IO_FOF_CM:
      case IO_FOF_VEL:
      case IO_FOF_LENTYPE:
      case IO_FOF_MASSTYPE:
      case IO_FOF_SFR:
      case IO_FOF_M_MEAN200:
      case IO_FOF_R_MEAN200:
      case IO_FOF_M_CRIT200:
      case IO_FOF_R_CRIT200:
      case IO_FOF_M_TOPHAT200:
      case IO_FOF_R_TOPHAT200:
      case IO_FOF_M_CRIT500:
      case IO_FOF_R_CRIT500:
      case IO_FOF_NSUBS:
      case IO_FOF_FIRSTSUB:
      case IO_FOF_FUZZOFFTYPE:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_J_MEAN200:
      case IO_FOF_JDM_MEAN200:
      case IO_FOF_JGAS_MEAN200:
      case IO_FOF_JSTARS_MEAN200:
      case IO_FOF_MASSTYPE_MEAN200:
      case IO_FOF_LENTYPE_MEAN200:
      case IO_FOF_CMFRAC_MEAN200:
      case IO_FOF_CMFRACTYPE_MEAN200:
      case IO_FOF_J_CRIT200:
      case IO_FOF_JDM_CRIT200:
      case IO_FOF_JGAS_CRIT200:
      case IO_FOF_JSTARS_CRIT200:
      case IO_FOF_MASSTYPE_CRIT200:
      case IO_FOF_LENTYPE_CRIT200:
      case IO_FOF_CMFRAC_CRIT200:
      case IO_FOF_CMFRACTYPE_CRIT200:
      case IO_FOF_J_TOPHAT200:
      case IO_FOF_JDM_TOPHAT200:
      case IO_FOF_JGAS_TOPHAT200:
      case IO_FOF_JSTARS_TOPHAT200:
      case IO_FOF_MASSTYPE_TOPHAT200:
      case IO_FOF_LENTYPE_TOPHAT200:
      case IO_FOF_CMFRAC_TOPHAT200:
      case IO_FOF_CMFRACTYPE_TOPHAT200:
      case IO_FOF_J_CRIT500:
      case IO_FOF_JDM_CRIT500:
      case IO_FOF_JGAS_CRIT500:
      case IO_FOF_JSTARS_CRIT500:
      case IO_FOF_MASSTYPE_CRIT500:
      case IO_FOF_LENTYPE_CRIT500:
      case IO_FOF_CMFRAC_CRIT500:
      case IO_FOF_CMFRACTYPE_CRIT500:
      case IO_FOF_J:
      case IO_FOF_JDM:
      case IO_FOF_JGAS:
      case IO_FOF_JSTARS:
      case IO_FOF_CMFRAC:
      case IO_FOF_CMFRACTYPE:
      case IO_FOF_EKIN:
      case IO_FOF_ETHR:
      case IO_FOF_EPOT:
      case IO_FOF_EPOT_CRIT200:
      case IO_FOF_EKIN_CRIT200:
      case IO_FOF_ETHR_CRIT200:
      case IO_FOF_EPOT_MEAN200:
      case IO_FOF_EKIN_MEAN200:
      case IO_FOF_ETHR_MEAN200:
      case IO_FOF_EPOT_TOPHAT200:
      case IO_FOF_EKIN_TOPHAT200:
      case IO_FOF_ETHR_TOPHAT200:
      case IO_FOF_EPOT_CRIT500:
      case IO_FOF_EKIN_CRIT500:
      case IO_FOF_ETHR_CRIT500:
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

        return 0;

      case IO_SUB_LEN:
      case IO_SUB_MTOT:
      case IO_SUB_POS:
      case IO_SUB_VEL:
      case IO_SUB_LENTYPE:
      case IO_SUB_MASSTYPE:
      case IO_SUB_CM:
      case IO_SUB_SPIN:
      case IO_SUB_VELDISP:
      case IO_SUB_VMAX:
      case IO_SUB_VMAXRAD:
      case IO_SUB_HALFMASSRAD:
      case IO_SUB_HALFMASSRADTYPE:
      case IO_SUB_MASSINRAD:
      case IO_SUB_MASSINHALFRAD:
      case IO_SUB_MASSINMAXRAD:
      case IO_SUB_MASSINRADTYPE:
      case IO_SUB_MASSINHALFRADTYPE:
      case IO_SUB_MASSINMAXRADTYPE:
      case IO_SUB_IDMOSTBOUND:
      case IO_SUB_GRNR:
      case IO_SUB_PARENT:
      case IO_SUB_BFLD_HALO:
      case IO_SUB_BFLD_DISK:
      case IO_SUB_SFR:
      case IO_SUB_SFRINRAD:
      case IO_SUB_SFRINHALFRAD:
      case IO_SUB_SFRINMAXRAD:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_SUB_EKIN:
      case IO_SUB_ETHR:
      case IO_SUB_EPOT:
      case IO_SUB_J:
      case IO_SUB_JDM:
      case IO_SUB_JGAS:
      case IO_SUB_JSTARS:
      case IO_SUB_JINHALFRAD:
      case IO_SUB_JDMINHALFRAD:
      case IO_SUB_JGASINHALFRAD:
      case IO_SUB_JSTARSINHALFRAD:
      case IO_SUB_JINRAD:
      case IO_SUB_JDMINRAD:
      case IO_SUB_JGASINRAD:
      case IO_SUB_JSTARSINRAD:
      case IO_SUB_CMFRAC:
      case IO_SUB_CMFRACTYPE:
      case IO_SUB_CMFRACINHALFRAD:
      case IO_SUB_CMFRACTYPEINHALFRAD:
      case IO_SUB_CMFRACINRAD:
      case IO_SUB_CMFRACTYPEINRAD:
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
        return 1;

      case IO_FOFSUB_IDS:
        return 2;

      case IO_FOF_LASTENTRY:
        terminate("reached last entry in switch - strange.");
        break;
    }

  terminate("reached end of function - this should not happen");
  return 0;
}

/*! \brief Returns number of particles of specific field.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return Number of entries of this property.
 */
int fof_subfind_get_particles_in_block(enum fof_subfind_iofields blocknr)
{
  switch(blocknr)
    {
      case IO_FOF_LEN:
      case IO_FOF_MTOT:
      case IO_FOF_POS:
      case IO_FOF_CM:
      case IO_FOF_VEL:
      case IO_FOF_LENTYPE:
      case IO_FOF_MASSTYPE:
      case IO_FOF_SFR:
      case IO_FOF_FUZZOFFTYPE:
        return catalogue_header.Ngroups;

      case IO_FOF_M_MEAN200:
      case IO_FOF_R_MEAN200:
      case IO_FOF_M_CRIT200:
      case IO_FOF_R_CRIT200:
      case IO_FOF_M_TOPHAT200:
      case IO_FOF_R_TOPHAT200:
      case IO_FOF_M_CRIT500:
      case IO_FOF_R_CRIT500:
      case IO_FOF_NSUBS:
      case IO_FOF_FIRSTSUB:

#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_J_MEAN200:
      case IO_FOF_JDM_MEAN200:
      case IO_FOF_JGAS_MEAN200:
      case IO_FOF_JSTARS_MEAN200:
      case IO_FOF_MASSTYPE_MEAN200:
      case IO_FOF_LENTYPE_MEAN200:
      case IO_FOF_CMFRAC_MEAN200:
      case IO_FOF_CMFRACTYPE_MEAN200:
      case IO_FOF_J_CRIT200:
      case IO_FOF_JDM_CRIT200:
      case IO_FOF_JGAS_CRIT200:
      case IO_FOF_JSTARS_CRIT200:
      case IO_FOF_MASSTYPE_CRIT200:
      case IO_FOF_LENTYPE_CRIT200:
      case IO_FOF_CMFRAC_CRIT200:
      case IO_FOF_CMFRACTYPE_CRIT200:
      case IO_FOF_J_TOPHAT200:
      case IO_FOF_JDM_TOPHAT200:
      case IO_FOF_JGAS_TOPHAT200:
      case IO_FOF_JSTARS_TOPHAT200:
      case IO_FOF_MASSTYPE_TOPHAT200:
      case IO_FOF_LENTYPE_TOPHAT200:
      case IO_FOF_CMFRAC_TOPHAT200:
      case IO_FOF_CMFRACTYPE_TOPHAT200:
      case IO_FOF_J_CRIT500:
      case IO_FOF_JDM_CRIT500:
      case IO_FOF_JGAS_CRIT500:
      case IO_FOF_JSTARS_CRIT500:
      case IO_FOF_MASSTYPE_CRIT500:
      case IO_FOF_LENTYPE_CRIT500:
      case IO_FOF_CMFRAC_CRIT500:
      case IO_FOF_CMFRACTYPE_CRIT500:
      case IO_FOF_J:
      case IO_FOF_JDM:
      case IO_FOF_JGAS:
      case IO_FOF_JSTARS:
      case IO_FOF_CMFRAC:
      case IO_FOF_CMFRACTYPE:
      case IO_FOF_EKIN:
      case IO_FOF_ETHR:
      case IO_FOF_EPOT:
      case IO_FOF_EPOT_CRIT200:
      case IO_FOF_EKIN_CRIT200:
      case IO_FOF_ETHR_CRIT200:
      case IO_FOF_EPOT_MEAN200:
      case IO_FOF_EKIN_MEAN200:
      case IO_FOF_ETHR_MEAN200:
      case IO_FOF_EPOT_TOPHAT200:
      case IO_FOF_EKIN_TOPHAT200:
      case IO_FOF_ETHR_TOPHAT200:
      case IO_FOF_EPOT_CRIT500:
      case IO_FOF_EKIN_CRIT500:
      case IO_FOF_ETHR_CRIT500:
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

#ifdef SUBFIND
        return catalogue_header.Ngroups;
#else  /* #ifdef SUBFIND */
        return 0;
#endif /* #ifdef SUBFIND #else */

      case IO_SUB_LEN:
      case IO_SUB_MTOT:
      case IO_SUB_POS:
      case IO_SUB_VEL:
      case IO_SUB_LENTYPE:
      case IO_SUB_MASSTYPE:
      case IO_SUB_CM:
      case IO_SUB_SPIN:
      case IO_SUB_VELDISP:
      case IO_SUB_VMAX:
      case IO_SUB_VMAXRAD:
      case IO_SUB_HALFMASSRAD:
      case IO_SUB_HALFMASSRADTYPE:
      case IO_SUB_MASSINRAD:
      case IO_SUB_MASSINHALFRAD:
      case IO_SUB_MASSINMAXRAD:
      case IO_SUB_MASSINRADTYPE:
      case IO_SUB_MASSINHALFRADTYPE:
      case IO_SUB_MASSINMAXRADTYPE:
      case IO_SUB_IDMOSTBOUND:
      case IO_SUB_GRNR:
      case IO_SUB_PARENT:
      case IO_SUB_BFLD_HALO:
      case IO_SUB_BFLD_DISK:
      case IO_SUB_SFR:
      case IO_SUB_SFRINRAD:
      case IO_SUB_SFRINHALFRAD:
      case IO_SUB_SFRINMAXRAD:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_SUB_EKIN:
      case IO_SUB_ETHR:
      case IO_SUB_EPOT:
      case IO_SUB_J:
      case IO_SUB_JDM:
      case IO_SUB_JGAS:
      case IO_SUB_JSTARS:
      case IO_SUB_JINHALFRAD:
      case IO_SUB_JDMINHALFRAD:
      case IO_SUB_JGASINHALFRAD:
      case IO_SUB_JSTARSINHALFRAD:
      case IO_SUB_JINRAD:
      case IO_SUB_JDMINRAD:
      case IO_SUB_JGASINRAD:
      case IO_SUB_JSTARSINRAD:
      case IO_SUB_CMFRAC:
      case IO_SUB_CMFRACTYPE:
      case IO_SUB_CMFRACINHALFRAD:
      case IO_SUB_CMFRACTYPEINHALFRAD:
      case IO_SUB_CMFRACINRAD:
      case IO_SUB_CMFRACTYPEINRAD:
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

#ifdef SUBFIND
        return catalogue_header.Nsubgroups;
#else  /* #ifdef SUBFIND */
        return 0;
#endif /* #ifdef SUBFIND #else */

      case IO_FOFSUB_IDS:
        return catalogue_header.Nids;

      case IO_FOF_LASTENTRY:
        terminate("reached last entry in switch - strange.");
        break;
    }

  terminate("reached end of function - this should not happen");
  return 0;
}

/*! \brief Returns the number of elements per entry of a given property.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return Number of values per element of the specified property.
 */
int fof_subfind_get_values_per_blockelement(enum fof_subfind_iofields blocknr)
{
  int values = 0;

  switch(blocknr)
    {
      case IO_FOF_LEN:
      case IO_FOF_NSUBS:
      case IO_FOF_FIRSTSUB:
      case IO_SUB_LEN:
      case IO_SUB_GRNR:
      case IO_SUB_PARENT:
      case IO_FOF_MTOT:
      case IO_FOF_SFR:
      case IO_FOF_M_MEAN200:
      case IO_FOF_R_MEAN200:
      case IO_FOF_M_CRIT200:
      case IO_FOF_R_CRIT200:
      case IO_FOF_M_TOPHAT200:
      case IO_FOF_R_TOPHAT200:
      case IO_FOF_M_CRIT500:
      case IO_FOF_R_CRIT500:
      case IO_SUB_MTOT:
      case IO_SUB_VELDISP:
      case IO_SUB_VMAX:
      case IO_SUB_VMAXRAD:
      case IO_SUB_HALFMASSRAD:
      case IO_SUB_MASSINRAD:
      case IO_SUB_MASSINHALFRAD:
      case IO_SUB_MASSINMAXRAD:
      case IO_SUB_IDMOSTBOUND:
      case IO_SUB_BFLD_HALO:
      case IO_SUB_BFLD_DISK:
      case IO_SUB_SFR:
      case IO_SUB_SFRINRAD:
      case IO_SUB_SFRINHALFRAD:
      case IO_SUB_SFRINMAXRAD:
      case IO_FOFSUB_IDS:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_CMFRAC_MEAN200:
      case IO_FOF_CMFRAC_CRIT200:
      case IO_FOF_CMFRAC_TOPHAT200:
      case IO_FOF_CMFRAC_CRIT500:
      case IO_FOF_EPOT_CRIT200:
      case IO_FOF_EKIN_CRIT200:
      case IO_FOF_ETHR_CRIT200:
      case IO_FOF_EPOT_MEAN200:
      case IO_FOF_EKIN_MEAN200:
      case IO_FOF_ETHR_MEAN200:
      case IO_FOF_EPOT_TOPHAT200:
      case IO_FOF_EKIN_TOPHAT200:
      case IO_FOF_ETHR_TOPHAT200:
      case IO_FOF_EPOT_CRIT500:
      case IO_FOF_EKIN_CRIT500:
      case IO_FOF_ETHR_CRIT500:
      case IO_FOF_EKIN:
      case IO_FOF_ETHR:
      case IO_FOF_EPOT:
      case IO_SUB_EKIN:
      case IO_SUB_ETHR:
      case IO_SUB_EPOT:
      case IO_SUB_CMFRAC:
      case IO_SUB_CMFRACINHALFRAD:
      case IO_SUB_CMFRACINRAD:
      case IO_FOF_CMFRAC:
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
        values = 1;
        break;

      case IO_FOF_LENTYPE:
      case IO_SUB_LENTYPE:
      case IO_FOF_MASSTYPE:
      case IO_SUB_MASSTYPE:
      case IO_SUB_HALFMASSRADTYPE:
      case IO_SUB_MASSINRADTYPE:
      case IO_SUB_MASSINHALFRADTYPE:
      case IO_SUB_MASSINMAXRADTYPE:
      case IO_FOF_FUZZOFFTYPE:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_CMFRACTYPE:
      case IO_SUB_CMFRACTYPE:
      case IO_SUB_CMFRACTYPEINHALFRAD:
      case IO_SUB_CMFRACTYPEINRAD:
      case IO_FOF_LENTYPE_MEAN200:
      case IO_FOF_LENTYPE_CRIT200:
      case IO_FOF_LENTYPE_CRIT500:
      case IO_FOF_LENTYPE_TOPHAT200:
      case IO_FOF_MASSTYPE_MEAN200:
      case IO_FOF_MASSTYPE_CRIT200:
      case IO_FOF_MASSTYPE_CRIT500:
      case IO_FOF_MASSTYPE_TOPHAT200:
      case IO_FOF_CMFRACTYPE_MEAN200:
      case IO_FOF_CMFRACTYPE_CRIT200:
      case IO_FOF_CMFRACTYPE_CRIT500:
      case IO_FOF_CMFRACTYPE_TOPHAT200:
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
        values = NTYPES;
        break;

      case IO_FOF_POS:
      case IO_FOF_CM:
      case IO_FOF_VEL:
      case IO_SUB_POS:
      case IO_SUB_VEL:
      case IO_SUB_CM:
      case IO_SUB_SPIN:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_SUB_J:
      case IO_SUB_JDM:
      case IO_SUB_JGAS:
      case IO_SUB_JSTARS:
      case IO_SUB_JINHALFRAD:
      case IO_SUB_JDMINHALFRAD:
      case IO_SUB_JGASINHALFRAD:
      case IO_SUB_JSTARSINHALFRAD:
      case IO_SUB_JINRAD:
      case IO_SUB_JDMINRAD:
      case IO_SUB_JGASINRAD:
      case IO_SUB_JSTARSINRAD:
      case IO_FOF_J_MEAN200:
      case IO_FOF_JDM_MEAN200:
      case IO_FOF_JGAS_MEAN200:
      case IO_FOF_JSTARS_MEAN200:
      case IO_FOF_J_CRIT200:
      case IO_FOF_JDM_CRIT200:
      case IO_FOF_JGAS_CRIT200:
      case IO_FOF_JSTARS_CRIT200:
      case IO_FOF_J_TOPHAT200:
      case IO_FOF_JDM_TOPHAT200:
      case IO_FOF_JGAS_TOPHAT200:
      case IO_FOF_JSTARS_TOPHAT200:
      case IO_FOF_J_CRIT500:
      case IO_FOF_JDM_CRIT500:
      case IO_FOF_JGAS_CRIT500:
      case IO_FOF_JSTARS_CRIT500:
      case IO_FOF_J:
      case IO_FOF_JDM:
      case IO_FOF_JGAS:
      case IO_FOF_JSTARS:
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
        values = 3;
        break;

      case IO_FOF_LASTENTRY:
        terminate("reached last entry in switch - should not get here");
        break;
    }
  return values;
}

/*! \brief Returns the number of bytes per element of a given property.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return Number of bytes per element for this property.
 */
int fof_subfind_get_bytes_per_blockelement(enum fof_subfind_iofields blocknr)
{
  int bytes_per_blockelement = 0;

  switch(blocknr)
    {
      case IO_FOF_LEN:
      case IO_FOF_NSUBS:
      case IO_FOF_FIRSTSUB:
      case IO_SUB_LEN:
      case IO_SUB_GRNR:
      case IO_SUB_PARENT:
        bytes_per_blockelement = sizeof(int);
        break;

      case IO_FOF_LENTYPE:
      case IO_SUB_LENTYPE:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_LENTYPE_MEAN200:
      case IO_FOF_LENTYPE_CRIT200:
      case IO_FOF_LENTYPE_CRIT500:
      case IO_FOF_LENTYPE_TOPHAT200:
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
        bytes_per_blockelement = NTYPES * sizeof(int);
        break;

      case IO_FOF_MTOT:
      case IO_FOF_SFR:
      case IO_FOF_M_MEAN200:
      case IO_FOF_R_MEAN200:
      case IO_FOF_M_CRIT200:
      case IO_FOF_R_CRIT200:
      case IO_FOF_M_TOPHAT200:
      case IO_FOF_R_TOPHAT200:
      case IO_FOF_M_CRIT500:
      case IO_FOF_R_CRIT500:
      case IO_SUB_MTOT:
      case IO_SUB_VELDISP:
      case IO_SUB_VMAX:
      case IO_SUB_VMAXRAD:
      case IO_SUB_HALFMASSRAD:
      case IO_SUB_MASSINRAD:
      case IO_SUB_MASSINHALFRAD:
      case IO_SUB_MASSINMAXRAD:
      case IO_SUB_BFLD_HALO:
      case IO_SUB_BFLD_DISK:
      case IO_SUB_SFR:
      case IO_SUB_SFRINRAD:
      case IO_SUB_SFRINHALFRAD:
      case IO_SUB_SFRINMAXRAD:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_CMFRAC_MEAN200:
      case IO_FOF_CMFRAC_CRIT200:
      case IO_FOF_CMFRAC_TOPHAT200:
      case IO_FOF_CMFRAC_CRIT500:
      case IO_FOF_CMFRAC:
      case IO_FOF_EKIN:
      case IO_FOF_ETHR:
      case IO_FOF_EPOT:
      case IO_SUB_EKIN:
      case IO_SUB_ETHR:
      case IO_SUB_EPOT:
      case IO_SUB_CMFRAC:
      case IO_SUB_CMFRACINHALFRAD:
      case IO_SUB_CMFRACINRAD:
      case IO_FOF_EPOT_CRIT200:
      case IO_FOF_EKIN_CRIT200:
      case IO_FOF_ETHR_CRIT200:
      case IO_FOF_EPOT_MEAN200:
      case IO_FOF_EKIN_MEAN200:
      case IO_FOF_ETHR_MEAN200:
      case IO_FOF_EPOT_TOPHAT200:
      case IO_FOF_EKIN_TOPHAT200:
      case IO_FOF_ETHR_TOPHAT200:
      case IO_FOF_EPOT_CRIT500:
      case IO_FOF_EKIN_CRIT500:
      case IO_FOF_ETHR_CRIT500:
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
        bytes_per_blockelement = sizeof(MyOutputFloat);
        break;

      case IO_FOF_POS:
      case IO_FOF_CM:
      case IO_FOF_VEL:
      case IO_SUB_POS:
      case IO_SUB_VEL:
      case IO_SUB_CM:
      case IO_SUB_SPIN:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_SUB_J:
      case IO_SUB_JDM:
      case IO_SUB_JGAS:
      case IO_SUB_JSTARS:
      case IO_SUB_JINHALFRAD:
      case IO_SUB_JDMINHALFRAD:
      case IO_SUB_JGASINHALFRAD:
      case IO_SUB_JSTARSINHALFRAD:
      case IO_SUB_JINRAD:
      case IO_SUB_JDMINRAD:
      case IO_SUB_JGASINRAD:
      case IO_SUB_JSTARSINRAD:
      case IO_FOF_J_MEAN200:
      case IO_FOF_JDM_MEAN200:
      case IO_FOF_JGAS_MEAN200:
      case IO_FOF_JSTARS_MEAN200:
      case IO_FOF_J_CRIT200:
      case IO_FOF_JDM_CRIT200:
      case IO_FOF_JGAS_CRIT200:
      case IO_FOF_JSTARS_CRIT200:
      case IO_FOF_J_TOPHAT200:
      case IO_FOF_JDM_TOPHAT200:
      case IO_FOF_JGAS_TOPHAT200:
      case IO_FOF_JSTARS_TOPHAT200:
      case IO_FOF_J_CRIT500:
      case IO_FOF_JDM_CRIT500:
      case IO_FOF_JGAS_CRIT500:
      case IO_FOF_JSTARS_CRIT500:
      case IO_FOF_J:
      case IO_FOF_JDM:
      case IO_FOF_JGAS:
      case IO_FOF_JSTARS:
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
        bytes_per_blockelement = 3 * sizeof(MyOutputFloat);
        break;

      case IO_FOF_MASSTYPE:
      case IO_SUB_MASSTYPE:
      case IO_SUB_HALFMASSRADTYPE:
      case IO_SUB_MASSINRADTYPE:
      case IO_SUB_MASSINHALFRADTYPE:
      case IO_SUB_MASSINMAXRADTYPE:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_MASSTYPE_MEAN200:
      case IO_FOF_MASSTYPE_CRIT200:
      case IO_FOF_MASSTYPE_CRIT500:
      case IO_FOF_MASSTYPE_TOPHAT200:
      case IO_FOF_CMFRACTYPE_MEAN200:
      case IO_FOF_CMFRACTYPE_CRIT200:
      case IO_FOF_CMFRACTYPE_CRIT500:
      case IO_FOF_CMFRACTYPE_TOPHAT200:
      case IO_FOF_CMFRACTYPE:
      case IO_SUB_CMFRACTYPE:
      case IO_SUB_CMFRACTYPEINHALFRAD:
      case IO_SUB_CMFRACTYPEINRAD:
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
        bytes_per_blockelement = NTYPES * sizeof(MyOutputFloat);
        break;

      case IO_SUB_IDMOSTBOUND:
      case IO_FOFSUB_IDS:
        bytes_per_blockelement = sizeof(MyIDType);
        break;

      case IO_FOF_FUZZOFFTYPE:
        bytes_per_blockelement = NTYPES * sizeof(long long);
        break;

      case IO_FOF_LASTENTRY:
        terminate("reached last entry in switch - should not get here");
        break;
    }
  return bytes_per_blockelement;
}

/*! \brief Returns key for datatype of element of a given property.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return Key for datatype: 0: int, 1: (output)float, 2: long long.
 */
int fof_subfind_get_datatype(enum fof_subfind_iofields blocknr)
{
  int typekey = 0;

  switch(blocknr)
    {
      case IO_FOF_LEN:
      case IO_FOF_LENTYPE:
      case IO_FOF_NSUBS:
      case IO_FOF_FIRSTSUB:
      case IO_SUB_LEN:
      case IO_SUB_LENTYPE:
      case IO_SUB_GRNR:
      case IO_SUB_PARENT:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_LENTYPE_MEAN200:
      case IO_FOF_LENTYPE_CRIT200:
      case IO_FOF_LENTYPE_CRIT500:
      case IO_FOF_LENTYPE_TOPHAT200:
#endif               /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
        typekey = 0; /* native int */
        break;

      case IO_FOF_MTOT:
      case IO_FOF_POS:
      case IO_FOF_CM:
      case IO_FOF_VEL:
      case IO_FOF_MASSTYPE:
      case IO_FOF_SFR:
      case IO_FOF_M_MEAN200:
      case IO_FOF_R_MEAN200:
      case IO_FOF_M_CRIT200:
      case IO_FOF_R_CRIT200:
      case IO_FOF_M_TOPHAT200:
      case IO_FOF_R_TOPHAT200:
      case IO_FOF_M_CRIT500:
      case IO_FOF_R_CRIT500:
      case IO_SUB_MTOT:
      case IO_SUB_POS:
      case IO_SUB_VEL:
      case IO_SUB_MASSTYPE:
      case IO_SUB_CM:
      case IO_SUB_SPIN:
      case IO_SUB_VELDISP:
      case IO_SUB_VMAX:
      case IO_SUB_VMAXRAD:
      case IO_SUB_HALFMASSRAD:
      case IO_SUB_HALFMASSRADTYPE:
      case IO_SUB_MASSINRAD:
      case IO_SUB_MASSINHALFRAD:
      case IO_SUB_MASSINMAXRAD:
      case IO_SUB_MASSINRADTYPE:
      case IO_SUB_MASSINHALFRADTYPE:
      case IO_SUB_MASSINMAXRADTYPE:
      case IO_SUB_BFLD_HALO:
      case IO_SUB_BFLD_DISK:
      case IO_SUB_SFR:
      case IO_SUB_SFRINRAD:
      case IO_SUB_SFRINHALFRAD:
      case IO_SUB_SFRINMAXRAD:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_MASSTYPE_MEAN200:
      case IO_FOF_MASSTYPE_CRIT200:
      case IO_FOF_MASSTYPE_CRIT500:
      case IO_FOF_MASSTYPE_TOPHAT200:
      case IO_FOF_J_MEAN200:
      case IO_FOF_JDM_MEAN200:
      case IO_FOF_JGAS_MEAN200:
      case IO_FOF_JSTARS_MEAN200:
      case IO_FOF_CMFRAC_MEAN200:
      case IO_FOF_CMFRACTYPE_MEAN200:
      case IO_FOF_J_CRIT200:
      case IO_FOF_JDM_CRIT200:
      case IO_FOF_JGAS_CRIT200:
      case IO_FOF_JSTARS_CRIT200:
      case IO_FOF_CMFRAC_CRIT200:
      case IO_FOF_CMFRACTYPE_CRIT200:
      case IO_FOF_J_TOPHAT200:
      case IO_FOF_JDM_TOPHAT200:
      case IO_FOF_JGAS_TOPHAT200:
      case IO_FOF_JSTARS_TOPHAT200:
      case IO_FOF_CMFRAC_TOPHAT200:
      case IO_FOF_CMFRACTYPE_TOPHAT200:
      case IO_FOF_J_CRIT500:
      case IO_FOF_JDM_CRIT500:
      case IO_FOF_JGAS_CRIT500:
      case IO_FOF_JSTARS_CRIT500:
      case IO_FOF_CMFRAC_CRIT500:
      case IO_FOF_CMFRACTYPE_CRIT500:
      case IO_FOF_J:
      case IO_FOF_JDM:
      case IO_FOF_JGAS:
      case IO_FOF_JSTARS:
      case IO_FOF_CMFRAC:
      case IO_FOF_CMFRACTYPE:
      case IO_FOF_EKIN:
      case IO_FOF_ETHR:
      case IO_FOF_EPOT:
      case IO_FOF_EPOT_CRIT200:
      case IO_FOF_EKIN_CRIT200:
      case IO_FOF_ETHR_CRIT200:
      case IO_FOF_EPOT_MEAN200:
      case IO_FOF_EKIN_MEAN200:
      case IO_FOF_ETHR_MEAN200:
      case IO_FOF_EPOT_TOPHAT200:
      case IO_FOF_EKIN_TOPHAT200:
      case IO_FOF_ETHR_TOPHAT200:
      case IO_FOF_EPOT_CRIT500:
      case IO_FOF_EKIN_CRIT500:
      case IO_FOF_ETHR_CRIT500:
      case IO_SUB_EKIN:
      case IO_SUB_ETHR:
      case IO_SUB_EPOT:
      case IO_SUB_J:
      case IO_SUB_JDM:
      case IO_SUB_JGAS:
      case IO_SUB_JSTARS:
      case IO_SUB_JINHALFRAD:
      case IO_SUB_JDMINHALFRAD:
      case IO_SUB_JGASINHALFRAD:
      case IO_SUB_JSTARSINHALFRAD:
      case IO_SUB_JINRAD:
      case IO_SUB_JDMINRAD:
      case IO_SUB_JGASINRAD:
      case IO_SUB_JSTARSINRAD:
      case IO_SUB_CMFRAC:
      case IO_SUB_CMFRACTYPE:
      case IO_SUB_CMFRACINHALFRAD:
      case IO_SUB_CMFRACTYPEINHALFRAD:
      case IO_SUB_CMFRACINRAD:
      case IO_SUB_CMFRACTYPEINRAD:
#endif               /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
        typekey = 1; /* native MyOutputFloat */
        break;

      case IO_SUB_IDMOSTBOUND:
      case IO_FOFSUB_IDS:
#ifdef LONGIDS
        typekey = 2; /* native long long */
#else                /* #ifdef LONGIDS */
        typekey = 0; /* native int */
#endif               /* #ifdef LONGIDS #else */
        break;

      case IO_FOF_FUZZOFFTYPE:
        typekey = 2; /* native long long */
        break;

      case IO_FOF_LASTENTRY:
        terminate("should not be reached");
        break;
    }

  return typekey;
}

/*! \brief Determines if block is present in the current code configuration.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return 0: not present; 1: present.
 */
int fof_subfind_blockpresent(enum fof_subfind_iofields blocknr)
{
  int present = 0;

  switch(blocknr)
    {
      case IO_FOF_LEN:
      case IO_FOF_LENTYPE:
      case IO_FOF_MTOT:
      case IO_FOF_POS:
      case IO_FOF_CM:
      case IO_FOF_VEL:
      case IO_FOF_MASSTYPE:
        present = 1;
        break;

      case IO_FOF_SFR:
      case IO_SUB_SFR:
      case IO_SUB_SFRINRAD:
      case IO_SUB_SFRINHALFRAD:
      case IO_SUB_SFRINMAXRAD:
#ifdef USE_SFR
        present = 1;
#endif /* #ifdef USE_SFR */
        break;

      case IO_SUB_BFLD_HALO:
      case IO_SUB_BFLD_DISK:
#ifdef MHD
        present = 1;
#endif /* #ifdef MHD */
        break;

      case IO_FOF_FUZZOFFTYPE:
        break;

      case IO_FOF_M_MEAN200:
      case IO_FOF_R_MEAN200:
      case IO_FOF_M_CRIT200:
      case IO_FOF_R_CRIT200:
      case IO_FOF_M_TOPHAT200:
      case IO_FOF_R_TOPHAT200:
      case IO_FOF_M_CRIT500:
      case IO_FOF_R_CRIT500:
      case IO_FOF_NSUBS:
      case IO_FOF_FIRSTSUB:
      case IO_SUB_LEN:
      case IO_SUB_LENTYPE:
      case IO_SUB_MTOT:
      case IO_SUB_POS:
      case IO_SUB_VEL:
      case IO_SUB_MASSTYPE:
      case IO_SUB_CM:
      case IO_SUB_SPIN:
      case IO_SUB_VELDISP:
      case IO_SUB_VMAX:
      case IO_SUB_VMAXRAD:
      case IO_SUB_HALFMASSRAD:
      case IO_SUB_HALFMASSRADTYPE:
      case IO_SUB_MASSINRAD:
      case IO_SUB_MASSINHALFRAD:
      case IO_SUB_MASSINMAXRAD:
      case IO_SUB_MASSINRADTYPE:
      case IO_SUB_MASSINHALFRADTYPE:
      case IO_SUB_MASSINMAXRADTYPE:
      case IO_SUB_IDMOSTBOUND:
      case IO_SUB_GRNR:
      case IO_SUB_PARENT:
#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_J_MEAN200:
      case IO_FOF_JDM_MEAN200:
      case IO_FOF_JGAS_MEAN200:
      case IO_FOF_JSTARS_MEAN200:
      case IO_FOF_CMFRAC_MEAN200:
      case IO_FOF_CMFRACTYPE_MEAN200:
      case IO_FOF_J_CRIT200:
      case IO_FOF_JDM_CRIT200:
      case IO_FOF_JGAS_CRIT200:
      case IO_FOF_JSTARS_CRIT200:
      case IO_FOF_CMFRAC_CRIT200:
      case IO_FOF_CMFRACTYPE_CRIT200:
      case IO_FOF_J_TOPHAT200:
      case IO_FOF_JDM_TOPHAT200:
      case IO_FOF_JGAS_TOPHAT200:
      case IO_FOF_JSTARS_TOPHAT200:
      case IO_FOF_CMFRAC_TOPHAT200:
      case IO_FOF_CMFRACTYPE_TOPHAT200:
      case IO_FOF_J_CRIT500:
      case IO_FOF_JDM_CRIT500:
      case IO_FOF_JGAS_CRIT500:
      case IO_FOF_JSTARS_CRIT500:
      case IO_FOF_CMFRAC_CRIT500:
      case IO_FOF_CMFRACTYPE_CRIT500:
      case IO_FOF_J:
      case IO_FOF_JDM:
      case IO_FOF_JGAS:
      case IO_FOF_JSTARS:
      case IO_FOF_CMFRAC:
      case IO_FOF_CMFRACTYPE:
      case IO_FOF_EKIN:
      case IO_FOF_ETHR:
      case IO_FOF_EPOT:
      case IO_FOF_MASSTYPE_MEAN200:
      case IO_FOF_MASSTYPE_CRIT200:
      case IO_FOF_MASSTYPE_CRIT500:
      case IO_FOF_MASSTYPE_TOPHAT200:
      case IO_FOF_LENTYPE_MEAN200:
      case IO_FOF_LENTYPE_CRIT200:
      case IO_FOF_LENTYPE_CRIT500:
      case IO_FOF_LENTYPE_TOPHAT200:
      case IO_FOF_EPOT_CRIT200:
      case IO_FOF_EKIN_CRIT200:
      case IO_FOF_ETHR_CRIT200:
      case IO_FOF_EPOT_MEAN200:
      case IO_FOF_EKIN_MEAN200:
      case IO_FOF_ETHR_MEAN200:
      case IO_FOF_EPOT_TOPHAT200:
      case IO_FOF_EKIN_TOPHAT200:
      case IO_FOF_ETHR_TOPHAT200:
      case IO_FOF_EPOT_CRIT500:
      case IO_FOF_EKIN_CRIT500:
      case IO_FOF_ETHR_CRIT500:
      case IO_SUB_EKIN:
      case IO_SUB_ETHR:
      case IO_SUB_EPOT:
      case IO_SUB_J:
      case IO_SUB_JDM:
      case IO_SUB_JGAS:
      case IO_SUB_JSTARS:
      case IO_SUB_JINHALFRAD:
      case IO_SUB_JDMINHALFRAD:
      case IO_SUB_JGASINHALFRAD:
      case IO_SUB_JSTARSINHALFRAD:
      case IO_SUB_JINRAD:
      case IO_SUB_JDMINRAD:
      case IO_SUB_JGASINRAD:
      case IO_SUB_JSTARSINRAD:
      case IO_SUB_CMFRAC:
      case IO_SUB_CMFRACTYPE:
      case IO_SUB_CMFRACINHALFRAD:
      case IO_SUB_CMFRACTYPEINHALFRAD:
      case IO_SUB_CMFRACINRAD:
      case IO_SUB_CMFRACTYPEINRAD:
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
#ifdef SUBFIND
        present = 1;
#else  /* #ifdef SUBFIND */
        present = 0;
#endif /* #ifdef SUBFIND #else */
        break;

      case IO_FOFSUB_IDS:
#ifdef FOF_STOREIDS
        present = 1;
#else  /* #ifdef FOF_STOREIDS */
        present = 0;
#endif /* #ifdef FOF_STOREIDS #else */
        break;

      case IO_FOF_LASTENTRY:
        terminate("should not be reached");
        break;
    }
  return present;
}

/*! \brief Get the 4 letter IO label for a given output field.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *  \param[out] label String with the label.
 *
 *  \return void
 */
void fof_subfind_get_Tab_IO_Label(enum fof_subfind_iofields blocknr, char *label)
{
  switch(blocknr)
    {
      case IO_FOF_LEN:
        strncpy(label, "FLEN", 4);
        break;
      case IO_FOF_MTOT:
        strncpy(label, "FMAS", 4);
        break;
      case IO_FOF_POS:
        strncpy(label, "FPOS", 4);
        break;
      case IO_FOF_CM:
        strncpy(label, "FGCM", 4);
        break;
      case IO_FOF_VEL:
        strncpy(label, "FVEL", 4);
        break;
      case IO_FOF_LENTYPE:
        strncpy(label, "FLTY", 4);
        break;
      case IO_FOF_MASSTYPE:
        strncpy(label, "FMTY", 4);
        break;
      case IO_FOF_SFR:
        strncpy(label, "FSFR", 4);
        break;
      case IO_FOF_M_MEAN200:
        strncpy(label, "FMM2", 4);
        break;
      case IO_FOF_R_MEAN200:
        strncpy(label, "FRM2", 4);
        break;
      case IO_FOF_M_CRIT200:
        strncpy(label, "FMC2", 4);
        break;
      case IO_FOF_R_CRIT200:
        strncpy(label, "FRC2", 4);
        break;
      case IO_FOF_M_TOPHAT200:
        strncpy(label, "FMT2", 4);
        break;
      case IO_FOF_R_TOPHAT200:
        strncpy(label, "FRT2", 4);
        break;
      case IO_FOF_M_CRIT500:
        strncpy(label, "FMC5", 4);
        break;
      case IO_FOF_R_CRIT500:
        strncpy(label, "FRC5", 4);
        break;
      case IO_FOF_NSUBS:
        strncpy(label, "FNSH", 4);
        break;
      case IO_FOF_FIRSTSUB:
        strncpy(label, "FFSH", 4);
        break;
      case IO_FOF_FUZZOFFTYPE:
        strncpy(label, "FUOF", 4);
        break;

      case IO_SUB_LEN:
        strncpy(label, "SLEN", 4);
        break;
      case IO_SUB_MTOT:
        strncpy(label, "SMAS", 4);
        break;
      case IO_SUB_POS:
        strncpy(label, "SPOS", 4);
        break;
      case IO_SUB_VEL:
        strncpy(label, "SVEL", 4);
        break;
      case IO_SUB_LENTYPE:
        strncpy(label, "SLTY", 4);
        break;
      case IO_SUB_MASSTYPE:
        strncpy(label, "SMTY", 4);
        break;
      case IO_SUB_CM:
        strncpy(label, "SCMP", 4);
        break;
      case IO_SUB_SPIN:
        strncpy(label, "SSPI", 4);
        break;
      case IO_SUB_VELDISP:
        strncpy(label, "SVDI", 4);
        break;
      case IO_SUB_VMAX:
        strncpy(label, "SVMX", 4);
        break;
      case IO_SUB_VMAXRAD:
        strncpy(label, "SVRX", 4);
        break;
      case IO_SUB_HALFMASSRAD:
        strncpy(label, "SHMR", 4);
        break;
      case IO_SUB_HALFMASSRADTYPE:
        strncpy(label, "SHMT", 4);
        break;
      case IO_SUB_MASSINRAD:
        strncpy(label, "SMIR", 4);
        break;
      case IO_SUB_MASSINHALFRAD:
        strncpy(label, "SMIH", 4);
        break;
      case IO_SUB_MASSINMAXRAD:
        strncpy(label, "SMIM", 4);
        break;
      case IO_SUB_MASSINRADTYPE:
        strncpy(label, "SMIT", 4);
        break;
      case IO_SUB_MASSINHALFRADTYPE:
        strncpy(label, "SMHT", 4);
        break;
      case IO_SUB_MASSINMAXRADTYPE:
        strncpy(label, "SMMT", 4);
        break;
      case IO_SUB_IDMOSTBOUND:
        strncpy(label, "SIDM", 4);
        break;
      case IO_SUB_GRNR:
        strncpy(label, "SGNR", 4);
        break;
      case IO_SUB_PARENT:
        strncpy(label, "SPRT", 4);
        break;
      case IO_SUB_BFLD_HALO:
        strncpy(label, "BFDH", 4);
        break;
      case IO_SUB_BFLD_DISK:
        strncpy(label, "BFDD", 4);
        break;
      case IO_SUB_SFR:
        strncpy(label, "SSFR", 4);
        break;
      case IO_SUB_SFRINRAD:
        strncpy(label, "SSFI", 4);
        break;
      case IO_SUB_SFRINHALFRAD:
        strncpy(label, "SSFH", 4);
        break;
      case IO_SUB_SFRINMAXRAD:
        strncpy(label, "SSFM", 4);
        break;
      case IO_FOFSUB_IDS:
        strncpy(label, "PIDS", 4);
        break;

#ifdef SUBFIND_EXTENDED_PROPERTIES
      case IO_FOF_J_MEAN200:
        strncpy(label, "FJM2", 4);
        break;
      case IO_FOF_JDM_MEAN200:
        strncpy(label, "JDM2", 4);
        break;
      case IO_FOF_JGAS_MEAN200:
        strncpy(label, "JGM2", 4);
        break;
      case IO_FOF_JSTARS_MEAN200:
        strncpy(label, "JSM2", 4);
        break;
      case IO_FOF_MASSTYPE_MEAN200:
        strncpy(label, "MTM2", 4);
        break;
      case IO_FOF_LENTYPE_MEAN200:
        strncpy(label, "LTM2", 4);
        break;
      case IO_FOF_CMFRAC_MEAN200:
        strncpy(label, "CFM2", 4);
        break;
      case IO_FOF_CMFRACTYPE_MEAN200:
        strncpy(label, "FTM2", 4);
        break;
      case IO_FOF_J_CRIT200:
        strncpy(label, "FJC2", 4);
        break;
      case IO_FOF_JDM_CRIT200:
        strncpy(label, "JDC2", 4);
        break;
      case IO_FOF_JGAS_CRIT200:
        strncpy(label, "JGC2", 4);
        break;
      case IO_FOF_JSTARS_CRIT200:
        strncpy(label, "JSC2", 4);
        break;
      case IO_FOF_MASSTYPE_CRIT200:
        strncpy(label, "MTC2", 4);
        break;
      case IO_FOF_LENTYPE_CRIT200:
        strncpy(label, "LTC2", 4);
        break;
      case IO_FOF_CMFRAC_CRIT200:
        strncpy(label, "CFC2", 4);
        break;
      case IO_FOF_CMFRACTYPE_CRIT200:
        strncpy(label, "FTC2", 4);
        break;
      case IO_FOF_J_TOPHAT200:
        strncpy(label, "FJT2", 4);
        break;
      case IO_FOF_JDM_TOPHAT200:
        strncpy(label, "JDT2", 4);
        break;
      case IO_FOF_JGAS_TOPHAT200:
        strncpy(label, "JGT2", 4);
        break;
      case IO_FOF_JSTARS_TOPHAT200:
        strncpy(label, "JST2", 4);
        break;
      case IO_FOF_MASSTYPE_TOPHAT200:
        strncpy(label, "MTT2", 4);
        break;
      case IO_FOF_LENTYPE_TOPHAT200:
        strncpy(label, "LTT2", 4);
        break;
      case IO_FOF_CMFRAC_TOPHAT200:
        strncpy(label, "CFT2", 4);
        break;
      case IO_FOF_CMFRACTYPE_TOPHAT200:
        strncpy(label, "FTT2", 4);
        break;
      case IO_FOF_J_CRIT500:
        strncpy(label, "FJC5", 4);
        break;
      case IO_FOF_JDM_CRIT500:
        strncpy(label, "JDC5", 4);
        break;
      case IO_FOF_JGAS_CRIT500:
        strncpy(label, "JGC5", 4);
        break;
      case IO_FOF_JSTARS_CRIT500:
        strncpy(label, "JSC5", 4);
        break;
      case IO_FOF_MASSTYPE_CRIT500:
        strncpy(label, "MTC5", 4);
        break;
      case IO_FOF_LENTYPE_CRIT500:
        strncpy(label, "LTC5", 4);
        break;
      case IO_FOF_CMFRAC_CRIT500:
        strncpy(label, "CFC5", 4);
        break;
      case IO_FOF_CMFRACTYPE_CRIT500:
        strncpy(label, "FTC5", 4);
        break;
      case IO_FOF_J:
        strncpy(label, "FOFJ", 4);
        break;
      case IO_FOF_JDM:
        strncpy(label, "FOJD", 4);
        break;
      case IO_FOF_JGAS:
        strncpy(label, "FOJG", 4);
        break;
      case IO_FOF_JSTARS:
        strncpy(label, "FOJS", 4);
        break;
      case IO_FOF_CMFRAC:
        strncpy(label, "FOCF", 4);
        break;
      case IO_FOF_CMFRACTYPE:
        strncpy(label, "FOFT", 4);
        break;
      case IO_FOF_EKIN:
        strncpy(label, "EKIN", 4);
        break;
      case IO_FOF_ETHR:
        strncpy(label, "ETHR", 4);
        break;
      case IO_FOF_EPOT:
        strncpy(label, "EPOT", 4);
        break;

      case IO_FOF_EPOT_CRIT200:
        strncpy(label, "EPO1", 4);
        break;
      case IO_FOF_EKIN_CRIT200:
        strncpy(label, "EKI1", 4);
        break;
      case IO_FOF_ETHR_CRIT200:
        strncpy(label, "ETH1", 4);
        break;
      case IO_FOF_EPOT_MEAN200:
        strncpy(label, "EPO2", 4);
        break;
      case IO_FOF_EKIN_MEAN200:
        strncpy(label, "EKI2", 4);
        break;
      case IO_FOF_ETHR_MEAN200:
        strncpy(label, "ETH2", 4);
        break;
      case IO_FOF_EPOT_TOPHAT200:
        strncpy(label, "EPO3", 4);
        break;
      case IO_FOF_EKIN_TOPHAT200:
        strncpy(label, "EKI3", 4);
        break;
      case IO_FOF_ETHR_TOPHAT200:
        strncpy(label, "ETH3", 4);
        break;
      case IO_FOF_EPOT_CRIT500:
        strncpy(label, "EPO4", 4);
        break;
      case IO_FOF_EKIN_CRIT500:
        strncpy(label, "EKI4", 4);
        break;
      case IO_FOF_ETHR_CRIT500:
        strncpy(label, "ETH4", 4);
        break;

      case IO_SUB_EKIN:
        strncpy(label, "SEKN", 4);
        break;
      case IO_SUB_ETHR:
        strncpy(label, "SETH", 4);
        break;
      case IO_SUB_EPOT:
        strncpy(label, "SEPT", 4);
        break;
      case IO_SUB_J:
        strncpy(label, "SUBJ", 4);
        break;
      case IO_SUB_JDM:
        strncpy(label, "SJDM", 4);
        break;
      case IO_SUB_JGAS:
        strncpy(label, "SJGS", 4);
        break;
      case IO_SUB_JSTARS:
        strncpy(label, "SJST", 4);
        break;
      case IO_SUB_JINHALFRAD:
        strncpy(label, "SJHR", 4);
        break;
      case IO_SUB_JDMINHALFRAD:
        strncpy(label, "SJDH", 4);
        break;
      case IO_SUB_JGASINHALFRAD:
        strncpy(label, "SJGH", 4);
        break;
      case IO_SUB_JSTARSINHALFRAD:
        strncpy(label, "SJSH", 4);
        break;
      case IO_SUB_JINRAD:
        strncpy(label, "SJMR", 4);
        break;
      case IO_SUB_JDMINRAD:
        strncpy(label, "SJDR", 4);
        break;
      case IO_SUB_JGASINRAD:
        strncpy(label, "SJGR", 4);
        break;
      case IO_SUB_JSTARSINRAD:
        strncpy(label, "SJSR", 4);
        break;
      case IO_SUB_CMFRAC:
        strncpy(label, "SCMF", 4);
        break;
      case IO_SUB_CMFRACTYPE:
        strncpy(label, "SCMT", 4);
        break;
      case IO_SUB_CMFRACINHALFRAD:
        strncpy(label, "SCMH", 4);
        break;
      case IO_SUB_CMFRACTYPEINHALFRAD:
        strncpy(label, "SCTH", 4);
        break;
      case IO_SUB_CMFRACINRAD:
        strncpy(label, "SCMR", 4);
        break;
      case IO_SUB_CMFRACTYPEINRAD:
        strncpy(label, "SCTR", 4);
        break;
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

      case IO_FOF_LASTENTRY:
        terminate("should not be reached");
        break;
    }
}

#ifdef HAVE_HDF5
/*! \brief Function that handles writing hdf5 header.
 *
 *  \param[in] handle Handle for header hdf5 group.
 *
 *  \return void
 */
void fof_subfind_write_header_attributes_in_hdf5(hid_t handle)
{
  hid_t hdf5_dataspace, hdf5_attribute;

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Ngroups_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &catalogue_header.Ngroups, "Ngroups_ThisFile");
  my_H5Aclose(hdf5_attribute, "Ngroups_ThisFile");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Nsubgroups_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &catalogue_header.Nsubgroups, "Nsubgroups_ThisFile");
  my_H5Aclose(hdf5_attribute, "Nsubgroups_ThisFile");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Nids_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &catalogue_header.Nids, "Nids_ThisFile");
  my_H5Aclose(hdf5_attribute, "Nids_ThisFile");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Ngroups_Total", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &catalogue_header.TotNgroups, "Ngroups_Total");
  my_H5Aclose(hdf5_attribute, "Ngroups_Total");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Nsubgroups_Total", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &catalogue_header.TotNsubgroups, "Nsubgroups_Total");
  my_H5Aclose(hdf5_attribute, "Nsubgroups_Total");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Nids_Total", H5T_NATIVE_INT64, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT64, &catalogue_header.TotNids, "Nids_Total");
  my_H5Aclose(hdf5_attribute, "Nids_Total");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "NumFiles", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &catalogue_header.num_files, "NumFiles");
  my_H5Aclose(hdf5_attribute, "NumFiles");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &catalogue_header.time, "Time");
  my_H5Aclose(hdf5_attribute, "Time");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &catalogue_header.redshift, "Redshift");
  my_H5Aclose(hdf5_attribute, "Redshift");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &catalogue_header.HubbleParam, "HubbleParam");
  my_H5Aclose(hdf5_attribute, "HubbleParam");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &catalogue_header.BoxSize, "BoxSize");
  my_H5Aclose(hdf5_attribute, "BoxSize");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &catalogue_header.Omega0, "Omega0");
  my_H5Aclose(hdf5_attribute, "Omega0");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &catalogue_header.OmegaLambda, "OmegaLambda");
  my_H5Aclose(hdf5_attribute, "OmegaLambda");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "FlagDoubleprecision", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &catalogue_header.flag_doubleprecision, "FlagDoubleprecision");
  my_H5Aclose(hdf5_attribute, "FlagDoubleprecision");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hid_t atype = my_H5Tcopy(H5T_C_S1);

  my_H5Tset_size(atype, strlen(GIT_COMMIT));
  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Git_commit", atype, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, atype, GIT_COMMIT, "Git_commit");
  my_H5Aclose(hdf5_attribute, "Git_commit");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  my_H5Tset_size(atype, strlen(GIT_DATE));
  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Git_date", atype, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, atype, GIT_DATE, "Git_date");
  my_H5Aclose(hdf5_attribute, "Git_date");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
}
#endif /* #ifdef HAVE_HDF5 */

#endif /* #ifdef FOF */
