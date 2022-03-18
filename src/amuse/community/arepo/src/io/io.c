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
 * \file        src/io.c
 * \date        05/2018
 * \brief       Routines for input and output of snapshot files to disk.
 * \details     contains functions:
 *                void init_field
 *                void init_units
 *                void init_snapshot_type
 *                void write_error
 *                void create_snapshot_if_desired(void)
 *                void produce_dump(void)
 *                void savepositions(int num, int subbox_flag)
 *                void fill_write_buffer
 *                int get_bytes_per_blockelement
 *                int get_datatype_in_block(enum iofields blocknr, int mode)
 *                int get_values_per_blockelement(enum iofields blocknr)
 *                int get_particles_in_block(enum iofields blocknr, int
 *                  *typelist)
 *                int blockpresent(enum iofields blocknr, int write)
 *                void get_Tab_IO_Label(enum iofields blocknr, char *label)
 *                void get_dataset_name(enum iofields blocknr, char *buf)
 *                void write_file(char *fname, int writeTask, int lastTask,
 *                  int subbox_flag)
 *                void write_header_attributes_in_hdf5(hid_t handle)
 *                void write_parameters_attributes_in_hdf5(hid_t handle)
 *                herr_t my_hdf5_error_handler(void *unused)
 *                void write_dataset_attributes(hid_t hdf5_dataset, enum
 *                  iofields blocknr)
 *                void write_xdmf(char *fname)
 *                size_t my_fwrite(void *ptr, size_t size, size_t nmemb,
 *                  FILE * stream)
 *                size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *
 *                  stream)
 *                void mpi_printf(const char *fmt, ...)
 *                void mpi_fprintf(FILE * stream, const char *fmt, ...)
 *                void mpi_printf_each(const char *fmt, ...)
 *                FILE *open_file(char *fnam)
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

/* needs to be included after allvars.h */
#ifdef OUTPUT_XDMF
#include <libgen.h> /* for basename() function */
#endif /* #ifdef OUTPUT_XDMF */

#include "../fof/fof.h"
#include "../gitversion/version.h"
#include "../mesh/voronoi/voronoi.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
void write_header_attributes_in_hdf5(hid_t handle);
void write_parameters_attributes_in_hdf5(hid_t handle);
void write_compile_time_options_in_hdf5(hid_t handle);
void write_dataset_attributes(hid_t hdf5_dataset, enum iofields blocknr);
#endif /* #ifdef HAVE_HDF5 */

#ifdef TOLERATE_WRITE_ERROR
static char alternative_fname[MAXLEN_PATH];
#endif /* #ifdef TOLERATE_WRITE_ERROR */

#ifdef OUTPUT_XDMF
static void write_xdmf(char *fname);
#endif /* #ifdef OUTPUT_XDMF */

static int n_type[NTYPES]; /**< contains the local (for a single task) number of particles of each type in the snapshot file */
static long long ntot_type_all[NTYPES]; /**< contains the global number of particles of each type in the snapshot file */
static int subbox_dump = 0;

/*! \brief Function for registering an output field.
 *
 *  Don't forget to add the new IO_FLAG to allvars.h.
 *
 *  \param[in] field Specifies the field as an enumeration type iofields
 *             (allvars.h), e.g. IO_POS. Don't forget to insert new fields
 *             also in allvars.h.
 *  \param[in] label The label of the dataset (4 characters).
 *  \param[in] datasetname The name of the hdf5 dataset (maximum 256
 *             characters).
 *  \param[in] type_in_memory The type of the field in the memory (use
 *             MEM_NONE if specifying io_func).
 *  \param[in] type_in_file_output The output type in the hdf5 file.
 *  \param[in] type_in_file_input The input type in the hdf5 file (use
 *             FILE_MY_OUTPUT_TYPE for MyInputType, input is disabled with
 *             FILE_NONE).
 *  \param[in] values_per_block The number of values per field, e.g. 1 for
 *             mass, 3 for velocities.
 *  \param[in] array The array in which the value is stored. For an io_func
 *             this influences the particle index, the default (A_NONE) is an
 *             index into P/SphP, can be changed if required.
 *  \param[in] pointer_to_field A Pointer to the field in one of the global
 *             arrays, e.g. &SphP[0].Density, or &P[0].Vel[0].
 *  \param[in] io_func Alternatively, if the value to output/input is not a
 *             simple field, you can define a function which handles i/o.
 *  \param[in] typelist_bitmask Specifies for which particle type the field is
 *             present, e.g. 1+2+8 => field present for particle types 0,1,3
 *             (or use ALL_TYPES, GAS_ONLY,...).
 *
 *  \return void
 */
void init_field(enum iofields field, const char *label, const char *datasetname, enum types_in_memory type_in_memory,
                enum types_in_file type_in_file_output, enum types_in_file type_in_file_input, int values_per_block, enum arrays array,
                void *pointer_to_field, void (*io_func)(int, int, void *, int), int typelist_bitmask)
{
  int alloc_step = 5;

  if(Max_IO_Fields == 0)
    {
      IO_Fields     = (IO_Field *)mymalloc("IO_Fields", alloc_step * sizeof(IO_Field));
      Max_IO_Fields = alloc_step;
    }
  else if(Max_IO_Fields == N_IO_Fields)
    {
      Max_IO_Fields = ((Max_IO_Fields / alloc_step) + 1) * alloc_step;
      IO_Fields     = (IO_Field *)myrealloc(IO_Fields, Max_IO_Fields * sizeof(IO_Field));
    }

  IO_Fields[N_IO_Fields].field = field;
  strncpy(IO_Fields[N_IO_Fields].label, label, 4);
  strncpy(IO_Fields[N_IO_Fields].datasetname, datasetname, 256);
  IO_Fields[N_IO_Fields].type_in_memory      = type_in_memory;
  IO_Fields[N_IO_Fields].type_in_file_output = type_in_file_output;
  IO_Fields[N_IO_Fields].type_in_file_input  = type_in_file_input;
  IO_Fields[N_IO_Fields].values_per_block    = values_per_block;
  IO_Fields[N_IO_Fields].snap_type           = SN_FULL;
  IO_Fields[N_IO_Fields].typelist            = typelist_bitmask;

  IO_Fields[N_IO_Fields].array = array;

  if(array == A_NONE)
    {
      IO_Fields[N_IO_Fields].offset = 0;
    }
  else if(array == A_SPHP)
    {
      IO_Fields[N_IO_Fields].offset = (size_t)pointer_to_field - (size_t)SphP;
    }
  else if(array == A_P)
    {
      IO_Fields[N_IO_Fields].offset = (size_t)pointer_to_field - (size_t)P;
    }
  else if(array == A_PS)
    {
      IO_Fields[N_IO_Fields].offset = (size_t)pointer_to_field - (size_t)PS;
    }

  IO_Fields[N_IO_Fields].io_func = io_func;

  // validate types
  if(type_in_memory == MEM_INT &&
     ((type_in_file_input != FILE_NONE && type_in_file_input != FILE_INT) || type_in_file_output != FILE_INT))
    {
      terminate("combination of datatypes not supported (field %s)", datasetname);
    }

  if(type_in_memory == MEM_MY_ID_TYPE &&
     ((type_in_file_input != FILE_NONE && type_in_file_input != FILE_MY_ID_TYPE) || type_in_file_output != FILE_MY_ID_TYPE))
    {
      terminate("combination of datatypes not supported (field %s)", datasetname);
    }

  if((type_in_memory == MEM_FLOAT || type_in_memory == MEM_MY_SINGLE || type_in_memory == MEM_DOUBLE) &&
     ((type_in_file_input != FILE_NONE && (type_in_file_input == FILE_MY_ID_TYPE || type_in_file_input == FILE_INT)) ||
      type_in_file_output == FILE_INT || type_in_file_output == FILE_MY_ID_TYPE))
    {
      terminate("combination of datatypes not supported (field %s)", datasetname);
    }

  IO_Fields[N_IO_Fields].a       = 0.;
  IO_Fields[N_IO_Fields].h       = 0.;
  IO_Fields[N_IO_Fields].L       = 0.;
  IO_Fields[N_IO_Fields].M       = 0.;
  IO_Fields[N_IO_Fields].V       = 0.;
  IO_Fields[N_IO_Fields].c       = 0.;
  IO_Fields[N_IO_Fields].hasunit = 0;

  N_IO_Fields++;
}

/*! \brief Function for adding units to output field.
 *
 *  This only works for fields registered with init_field.
 *
 *  \param[in] field Specifies the field as an enumeration type iofields
 *             (allvars.h), e.g. IO_POS.
 *  \param[in] a the exponent of the cosmological a factor.
 *  \param[in] h the exponent of the hubble parameter.
 *  \param[in] L the length unit scaling.
 *  \param[in] M the mass unit scaling.
 *  \param[in] V the velocity unit scaling.
 *  \param[in] c conversion factor to cgs units (zero indicates dimensionless
 *             quantity, integer count, etc).
 *
 *  \return void
 */
void init_units(enum iofields field, double a, double h, double L, double M, double V, double c)
{
  for(int i = 0; i < N_IO_Fields; i++)
    {
      if(IO_Fields[i].field == field)
        {
          IO_Fields[i].hasunit = 1;
          IO_Fields[i].a       = a;
          IO_Fields[i].h       = h;
          IO_Fields[i].L       = L;
          IO_Fields[i].M       = M;
          IO_Fields[i].V       = V;
          IO_Fields[i].c       = c;
          break;
        }
    }
}

/*! \brief Function for determining whether a field is dumped in snapshot.
 *
 *  This only works for fields registered with init_field.
 *  The member snap_type is initialized to SN_FULL in init_field.
 *
 *  \param[in] field Specifies the field as an enumeration type iofields
 *             (allvars.h), e.g. IO_POS.
 *  \param[in] type In which snapshot types this field should be present
 *             (e.g. SN_FULL).
 *
 *  \return void
 */
void init_snapshot_type(enum iofields field, enum sn_type type)
{
  for(int i = 0; i < N_IO_Fields; i++)
    {
      if(IO_Fields[i].field == field)
        {
          IO_Fields[i].snap_type = type;
        }
    }
}

#ifdef TOLERATE_WRITE_ERROR
/*! \brief Print information about a write error.
 *
 *  If a write error occurs, this function prints some useful debug information
 *  and sets to 1 the variable WriteErrorFlag so that the write operation that
 *  caused the error can be performed again.
 *
 *  \param[in] check Flag that indicates where the function was called [0 and 1
 *             in my_fwrite(), 2 in my_hdf5_error_handler(), 3 in
 *             hdf5_header_error_handler()].
 *  \param[in] nwritten Number of elements actually written.
 *  \param[in] nmemb Number of elements that should be written.
 *
 *  \return void
 */
void write_error(int check, size_t nwritten, size_t nmemb)
{
  if(!WriteErrorFlag)
    {
      int len;
      char hostname[MPI_MAX_PROCESSOR_NAME];
      MPI_Get_processor_name(hostname, &len);

      printf("TOLERATE_WRITE_ERROR: write failed node=%s  nwritten=%lld  nmemb=%lld  errno=%s  task=%d  check=%d\n", hostname,
             (long long)nwritten, (long long)nmemb, strerror(errno), ThisTask, check);
      myflush(stdout);
      WriteErrorFlag = 1;
    }
}
#endif /* #ifdef TOLERATE_WRITE_ERROR */

/*! \brief Checks if a snapshot should be saved.
 *
 *  This function checks whether a snapshot file or other kinds of output
 *  files, such as a projection, should be saved at the current time-step.
 *  If that is the case, the appropriate functions to produce the desired
 *  file are called and the parameter controlling the output are updated
 *  accordingly.
 *
 *  \return void
 */
void create_snapshot_if_desired(void)
{
#ifdef OUTPUT_EVERY_STEP
  All.Ti_nextoutput = All.Ti_Current;
#endif /* #ifdef OUTPUT_EVERY_STEP */

  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) /* allow only top-level synchronization points */
    if(All.Ti_Current >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
      {
        DumpFlag = DumpFlagNextSnap;
        produce_dump();

        All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);
      }
}

/*! \brief A wrapper function used to create a snapshot.
 *
 *  This function wraps together savepositions(), the function that
 *  saves the snapshot file to the disk, with functions used for
 *  special output needs.
 *
 *  \return void
 */
void produce_dump(void)
{
#ifdef UPDATE_GRADIENTS_FOR_OUTPUT
  exchange_primitive_variables();
  calculate_gradients();
#endif /* #ifdef UPDATE_GRADIENTS_FOR_OUTPUT */

  savepositions(All.SnapshotFileCount++, 0); /* write snapshot file */
}

/*! \brief Saves snapshot to disk.
 *
 *  This function writes a snapshot of the particle distribution to one or
 *  several files. If NumFilesPerSnapshot>1, the snapshot is distributed
 *  into several files, which are written simultaneously. Each file contains
 *  data from a group of processors of size roughly NTask/NumFilesPerSnapshot.
 *
 *  \param[in] num The snapshot number.
 *  \param[in] subbox_flag If greater than 0 instructs the code to output only
 *             a subset of the whole domain.
 *
 *  \return void
 */
void savepositions(int num, int subbox_flag)
{
  char buf[500];
  int n, filenr, gr, ngroups, masterTask, lastTask;
  double t0, t1;

  t0 = second();
  CPU_Step[CPU_MISC] += measure_time();

  if(DumpFlag)
    {
      subbox_dump = 0;

      if(subbox_flag > 0)
        {
          mpi_printf("\nwriting small subbox #%d snapshot file #%d @ time %g ... \n", subbox_flag - 1, num, All.Time);
          subbox_dump = 1;
        }
      else
        mpi_printf("\nwriting snapshot file #%d @ time %g ... (DumpFlag=%d)\n", num, All.Time, DumpFlag);

#ifdef FOF
      if(RestartFlag != 3 && RestartFlag != 18 && subbox_flag == 0 && DumpFlag != 2)
        {
          {
            mpi_printf("\nWe shall first compute a group catalogue for this snapshot file\n");

            fof_fof(num);
          }
        }
#endif /* #ifdef FOF */

      if(DumpFlag != 4)
        {
          CommBuffer = mymalloc("CommBuffer", COMMBUFFERSIZE);

          if(NTask < All.NumFilesPerSnapshot)
            {
              warn(
                  "Number of processors must be larger or equal than All.NumFilesPerSnapshot! Reducing All.NumFilesPerSnapshot "
                  "accordingly.\n");
              All.NumFilesPerSnapshot = NTask;
            }

          if(All.SnapFormat < 1 || All.SnapFormat > 3)
            terminate("Unsupported File-Format.  All.SnapFormat=%d\n", All.SnapFormat);

#ifndef HAVE_HDF5
          if(All.SnapFormat == 3)
            {
              mpi_terminate("Code wasn't compiled with HDF5 support enabled!\n");
            }
#endif /* #ifndef  HAVE_HDF5 */

          /* determine global and local particle numbers */
          for(n = 0; n < NTYPES; n++)
            n_type[n] = 0;

          for(n = 0; n < NumPart; n++)
            {
              n_type[P[n].Type]++;
            }

          sumup_large_ints(NTYPES, n_type, ntot_type_all);

          /* assign processors to output files */
          distribute_file(All.NumFilesPerSnapshot, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);

          if(All.NumFilesPerSnapshot > 1)
            {
              if(ThisTask == 0)
                {
                  sprintf(buf, "%s/snapdir_%03d", All.OutputDir, num);
                  mkdir(buf, 02755);

#ifdef TOLERATE_WRITE_ERROR
                  sprintf(alternative_fname, "%s/snapdir_%03d", AlternativeOutputDir, num);
                  mkdir(alternative_fname, 02755);
#endif /* #ifdef TOLERATE_WRITE_ERROR */
                }

              MPI_Barrier(MPI_COMM_WORLD);
            }

          if(All.NumFilesPerSnapshot > 1)
            sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", All.OutputDir, num, All.SnapshotFileBase, num, filenr);
          else
            sprintf(buf, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, num);

#ifdef TOLERATE_WRITE_ERROR
          if(All.NumFilesPerSnapshot > 1)
            sprintf(alternative_fname, "%s/snapdir_%03d/%s_%03d.%d", AlternativeOutputDir, num, All.SnapshotFileBase, num, filenr);
          else
            sprintf(alternative_fname, "%s%s_%03d", AlternativeOutputDir, All.SnapshotFileBase, num);
#endif /* #ifdef TOLERATE_WRITE_ERROR */

          if(RestartFlag == 3)
            {
#ifndef FOF_STOREIDS
              if(All.NumFilesPerSnapshot > 1)
                sprintf(buf, "%s/snapdir_%03d/%s-groupordered_%03d.%d", All.OutputDir, num, All.SnapshotFileBase, num, filenr);
              else
                sprintf(buf, "%s%s-groupordered_%03d", All.OutputDir, All.SnapshotFileBase, num);
#else  /* #ifndef FOF_STOREIDS */
              if(All.NumFilesPerSnapshot > 1)
                sprintf(buf, "%s/snapdir_%03d/%s-storeids_%03d.%d", All.OutputDir, num, All.SnapshotFileBase, num, filenr);
              else
                sprintf(buf, "%s%s-storeids_%03d", All.OutputDir, All.SnapshotFileBase, num);
#endif /* #ifndef FOF_STOREIDS #else */
            }

#ifdef ADDBACKGROUNDGRID
          if(All.NumFilesPerSnapshot > 1)
            sprintf(buf, "%s-with-grid.%d", All.InitCondFile, filenr);
          else
            sprintf(buf, "%s-with-grid", All.InitCondFile);
#endif /* #ifdef ADDBACKGROUNDGRID */

          ngroups = All.NumFilesPerSnapshot / All.NumFilesWrittenInParallel;
          if((All.NumFilesPerSnapshot % All.NumFilesWrittenInParallel))
            ngroups++;

          for(gr = 0; gr < ngroups; gr++)
            {
              if((filenr / All.NumFilesWrittenInParallel) == gr) /* ok, it's this processor's turn */
                {
                  if(ThisTask == masterTask && (filenr % All.NumFilesWrittenInParallel) == 0)
                    printf("writing snapshot files group %d out of %d - files %d-%d (total of %d files): '%s'\n", gr + 1, ngroups,
                           filenr, filenr + All.NumFilesWrittenInParallel - 1, All.NumFilesPerSnapshot, buf);
                  write_file(buf, masterTask, lastTask, subbox_flag);
#ifdef OUTPUT_XDMF
                  if(All.SnapFormat == 3)
                    {
                      write_xdmf(buf);
                    }
#endif /* #ifdef OUTPUT_XDMF */
                }
              MPI_Barrier(MPI_COMM_WORLD);
            }

          myfree(CommBuffer);

          t1 = second();
          CPU_Step[CPU_SNAPSHOT] += measure_time();

          mpi_printf("done with writing snapshot (took %g sec).\n", timediff(t0, t1));
        }
      else
        {
          mpi_printf("done with writing files: no dump of snapshot (DumpFlag = %d).\n", DumpFlag);
        }  // if(DumpFlag !=4)

#ifdef FOF
      if(RestartFlag != 3 && RestartFlag != 6 && RestartFlag != 18 && subbox_flag == 0 && DumpFlag != 2)
        {
          {
#ifndef FOF_STOREIDS
            /* now revert from output order to the original order */
            for(n = 0; n < NumPart; n++)
              {
                PS[n].TargetTask  = PS[n].OriginTask;
                PS[n].TargetIndex = PS[n].OriginIndex;
              }

            fof_subfind_exchange(MPI_COMM_WORLD);

            myfree(PS);

            /* do resize because subfind may have increased these limits */
            if(All.MaxPart != fof_OldMaxPart)
              {
                All.MaxPart = fof_OldMaxPart;
                reallocate_memory_maxpart();
              }
            if(All.MaxPartSph != fof_OldMaxPartSph)
              {
                All.MaxPartSph = fof_OldMaxPartSph;
                reallocate_memory_maxpartsph();
              }

            CPU_Step[CPU_FOF] += measure_time();
#endif /* #ifndef FOF_STOREIDS */

            /* recreate the mesh that we had free to reduce peak memory usage */
            create_mesh();
            mesh_setup_exchange();
          }
        }
#endif /* #ifdef FOF */

      All.Ti_lastoutput = All.Ti_Current;

      CPU_Step[CPU_SNAPSHOT] += measure_time();
    }
}

/*! \brief This function fills the write buffer with particle data.
 *
 *  \param[out] buffer Buffer to be filled.
 *  \param[in] blocknr ID of the output block (i.e. position, velocities...).
 *  \param[in, out] startindex Pointer containing the offset in write buffer.
 *  \param[in] pc Number of particle to be put in the buffer.
 *  \param[in] type Particle type.
 *  \param[in] subbox_flag If greater than 0 instructs the code to output
 *             only a subset of the whole domain.
 *
 *  \return void
 */
void fill_write_buffer(void *buffer, enum iofields blocknr, int *startindex, int pc, int type, int subbox_flag)
{
  int n, k, pindex, f;
  MyOutputFloat *fp;
  MyIDType *ip;
  int *intp;

  /* determine which field we are working on */
  int field = -1;

  for(f = 0; f < N_IO_Fields; f++)
    {
      if(IO_Fields[f].field == blocknr)
        {
          field = f;
          break;
        }
    }

  if(field < 0)
    terminate("IO field=%d not registered with init_field()", (int)blocknr);

  set_cosmo_factors_for_current_time();

  fp              = (MyOutputFloat *)buffer;
  ip              = (MyIDType *)buffer;
  intp            = (int *)buffer;
  double *doublep = (double *)buffer;
  float *floatp   = (float *)buffer;

  pindex = *startindex;

  for(n = 0; n < pc; pindex++)
    {
      /* SUBBOX_SNAPSHOTS specialized output */

      /* normal particle output */
      if(P[pindex].Type == type)
        {
          if(IO_Fields[field].io_func)
            {
              int particle;
              switch(IO_Fields[field].array)
                {
                  case A_NONE:
                  case A_SPHP:
                  case A_P:
                    particle = pindex;
                    break;
                  case A_PS:
                    terminate("Not good, trying to read into PS[]?\n");
                    break;
                  default:
                    terminate("ERROR in fill_write_buffer: Array not found!\n");
                    break;
                }

              switch(IO_Fields[field].type_in_file_output)
                {
                  case FILE_NONE:
                    terminate("error");
                    break;
                  case FILE_INT:
                    IO_Fields[field].io_func(particle, IO_Fields[field].values_per_block, intp, 0);
                    intp += IO_Fields[field].values_per_block;
                    n++;
                    break;
                  case FILE_MY_ID_TYPE:
                    IO_Fields[field].io_func(particle, IO_Fields[field].values_per_block, ip, 0);
                    ip += IO_Fields[field].values_per_block;
                    n++;
                    break;
                  case FILE_MY_IO_FLOAT:
                    IO_Fields[field].io_func(particle, IO_Fields[field].values_per_block, fp, 0);
                    fp += IO_Fields[field].values_per_block;
                    n++;
                    break;
                  case FILE_DOUBLE:
                    IO_Fields[field].io_func(particle, IO_Fields[field].values_per_block, doublep, 0);
                    doublep += IO_Fields[field].values_per_block;
                    n++;
                    break;
                  case FILE_FLOAT:
                    IO_Fields[field].io_func(particle, IO_Fields[field].values_per_block, floatp, 0);
                    floatp += IO_Fields[field].values_per_block;
                    n++;
                    break;
                }
            }
          else
            {
              void *array_pos;

              switch(IO_Fields[field].array)
                {
                  case A_NONE:
                    array_pos = 0;
                    break;

                  case A_SPHP:
                    array_pos = SphP + pindex;
                    break;

                  case A_P:
                    array_pos = P + pindex;
                    break;
                  case A_PS:
                    array_pos = PS + pindex;
                    break;

                  default:
                    terminate("ERROR in fill_write_buffer: Array not found!\n");
                    break;
                }

              for(k = 0; k < IO_Fields[field].values_per_block; k++)
                {
                  double value = 0.;

                  switch(IO_Fields[field].type_in_memory)
                    {
                      case MEM_INT:
                        *intp = *((int *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(int)));
                        intp++;
                        break;

                      case MEM_MY_ID_TYPE:
                        *ip = *((MyIDType *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(MyIDType)));
                        ip++;
                        break;

                      case MEM_FLOAT:
                        value = *((float *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(float)));
                        break;

                      case MEM_DOUBLE:
                        value = *((double *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(double)));
                        break;

                      case MEM_MY_SINGLE:
                        value = *((MySingle *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(MySingle)));
                        break;

                      case MEM_MY_FLOAT:
                        value = *((MyFloat *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(MyFloat)));
                        break;

                      case MEM_MY_DOUBLE:
                        value = *((MyDouble *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(MyDouble)));
                        break;

                      case MEM_NONE:
                        terminate("ERROR in fill_write_buffer: reached MEM_NONE with no io_func specified!\n");
                        break;

                      default:
                        terminate("ERROR in fill_write_buffer: Type not found!\n");
                        break;
                    }

                  switch(IO_Fields[field].type_in_file_output)
                    {
                      case FILE_MY_IO_FLOAT:
                        *fp = value;
                        fp++;
                        break;

                      case FILE_DOUBLE:
                        *doublep = value;
                        doublep++;
                        break;

                      case FILE_FLOAT:
                        *floatp = value;
                        floatp++;
                        break;

                      default:
                        break;
                    }
                }

              n++;
            }  // end io_func/not
        }      // end type if
    }          // end particle loop

  *startindex = pindex;
}

/*! \brief This function tells the size in bytes of one data entry in each of
 *         the blocks defined for the output file.
 *
 *  \param[in] blocknr ID of the output block (i.e. position, velocities...).
 *  \param[in] mode Used to distinguish whether the function is called in input
 *             mode (mode > 0) or in output mode (mode = 0). The size of one
 *             data entry may vary depending on the mode.
 *
 *  \return Size of the data entry in bytes.
 */
int get_bytes_per_blockelement(enum iofields blocknr, int mode)
{
  int bytes_per_blockelement = 0;
  int f;

  for(f = 0; f < N_IO_Fields; f++)
    {
      if(IO_Fields[f].field == blocknr)
        {
          if(mode)
            {
              switch(IO_Fields[f].type_in_file_input)
                {
                  case FILE_NONE:
                    terminate("error");
                    break;
                  case FILE_INT:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(int);
                    break;
                  case FILE_MY_ID_TYPE:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(MyIDType);
                    break;
                  case FILE_MY_IO_FLOAT:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(MyInputFloat);
                    break;
                  case FILE_DOUBLE:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(double);
                    break;
                  case FILE_FLOAT:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(float);
                    break;
                }
            }
          else
            {
              switch(IO_Fields[f].type_in_file_output)
                {
                  case FILE_NONE:
                    terminate("error");
                    break;
                  case FILE_INT:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(int);
                    break;
                  case FILE_MY_ID_TYPE:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(MyIDType);
                    break;
                  case FILE_MY_IO_FLOAT:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(MyOutputFloat);
                    break;
                  case FILE_DOUBLE:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(double);
                    break;
                  case FILE_FLOAT:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(float);
                    break;
                }
            }
          break;
        }
    }

  return bytes_per_blockelement;
}

/*! \brief This function determines the type of one data entry in each of the
 *         blocks defined for the output file.
 *
 *  Used only if output in HDF5 format is enabled.
 *
 *  \param[in] blocknr ID of the output block (i.e. position, velocities...).
 *  \param[in] mode For input mode > 0, for output mode = 0.
 *
 *  \return typekey, a flag that indicates the type of the data entry.
 */
int get_datatype_in_block(enum iofields blocknr, int mode)
{
  int typekey, f;

  for(f = 0; f < N_IO_Fields; f++)
    {
      if(IO_Fields[f].field == blocknr)
        {
          if(mode)
            typekey = IO_Fields[f].type_in_file_input;
          else
            typekey = IO_Fields[f].type_in_file_output;

          return typekey;
        }
    }

  terminate("error invalid field");
  return typekey;
}

/*! \brief This function determines the number of elements composing one data
 *         entry in each of the blocks defined for the output file.
 *
 *  Used only if output in HDF5 format is enabled.
 *
 *  \param[in] blocknr ID of the output block (i.e. position, velocities...).
 *
 *  \return Number of elements of one data entry.
 */
int get_values_per_blockelement(enum iofields blocknr)
{
  int values = 0;
  int f;

  for(f = 0; f < N_IO_Fields; f++)
    {
      if(IO_Fields[f].field == blocknr)
        {
          values = IO_Fields[f].values_per_block;
          return values;
        }
    }

  terminate("reached last entry in switch - strange.");
  return values;
}

/*! \brief Gets particle number in an output block.
 *
 *  This function determines how many particles there are in a given block,
 *  based on the information in the header-structure.  It also flags particle
 *  types that are present in the block in the typelist array.
 *
 *  \param[in] blocknr ID of the output block (i.e. position, velocities...).
 *  \param[in] typelist Array that contains the number of particles of each
 *             type in the block.
 *
 *  \return The total number of particles in the block.
 */
int get_particles_in_block(enum iofields blocknr, int *typelist)
{
  int i, f;
  int npart = 0;

  switch(blocknr)
    {
      case IO_MASS:
        for(i = 0; i < NTYPES; i++)
          {
            typelist[i] = 0;
            if(All.MassTable[i] == 0)
              if(header.npart[i] > 0)
                {
                  typelist[i] = 1;
                  npart += header.npart[i];
                }
          }
        return npart; /* with masses */
        break;

      case IO_LASTENTRY:
        terminate("reached last entry in switch - strange.");
        break;

      default:
        for(f = 0; f < N_IO_Fields; f++)
          {
            if(IO_Fields[f].field == blocknr)
              {
                for(i = 0; i < NTYPES; i++)
                  {
                    if((IO_Fields[f].typelist & (1 << i)) && header.npart[i] > 0)
                      {
                        typelist[i] = 1;
                        npart += header.npart[i];
                      }
                    else
                      typelist[i] = 0;
                  }

                return npart;
              }
          }
        break;

    }  // end switch

  terminate("reached end of function - this should not happen");
  return 0;
}

/*! \brief Checks if a block is expected for file input or output.
 *
 *  This function tells whether a block in the input/output file is requested
 *  or not. Because the blocks processed in the two cases are different, the
 *  mode is indicated with the flag write (1=write, 0=read).
 *
 *  \param[in] blocknr ID of the output block (i.e. position, velocities...).
 *  \param[in] write If 0 the function is in read mode, if 1 the function is
 *             in write mode.
 *
 *  \return 0 if the block is not present, 1 otherwise.
 */
int blockpresent(enum iofields blocknr, int write)
{
  int f;

  if(!write)
    {
#ifdef PASSIVE_SCALARS
      if(RestartFlag == 0 && blocknr == IO_PASS)
        return 1;
#endif /* #ifdef PASSIVE_SCALARS */
#if defined(MHD) && !defined(MHD_SEEDFIELD)
      if(All.ICFormat != 3 && RestartFlag == 0 && (blocknr > IO_U && blocknr != IO_BFLD))
#else  /* #if defined(MHD) && !defined(MHD_SEEDFIELD) */
      if(All.ICFormat != 3 && RestartFlag == 0 && blocknr > IO_U)
#endif /* #if defined(MHD) && !defined(MHD_SEEDFIELD) #else */
#ifdef READ_LEGACY_ICS
        if(RestartFlag == 0 && blocknr > IO_U && blocknr != IO_BFLD)
#else               /* #ifdef  READ_LEGACY_ICS */
        if(RestartFlag == 0)
#endif              /* #ifdef  READ_LEGACY_ICS #else */
          return 0; /* ignore all other blocks in non-HDF5 initial conditions */
    }

  for(f = 0; f < N_IO_Fields; f++)
    {
      if(IO_Fields[f].field == blocknr)
        {
          if(!write)
            {
              if(IO_Fields[f].type_in_file_input != FILE_NONE)
                {
                  return 1;
                }
            }
          else
            {
              if(IO_Fields[f].type_in_file_output == FILE_NONE)
                return 0;

              /* subboxes: write all fields except those marked by SN_NO_SUBBOX or SN_MINI_ONLY
                 (must come first to ignore DumpFlag) */
              if(subbox_dump)
                {
                  if(IO_Fields[f].snap_type == SN_NO_SUBBOX || IO_Fields[f].snap_type == SN_MINI_ONLY)
                    return 0;

                  return 1;
                }

              /* normal full snapshot (with or without groupcat): only skip fields marked by SN_MINI_ONLY */
              if(DumpFlag == 1 || DumpFlag == 2)
                {
                  if(IO_Fields[f].snap_type == SN_MINI_ONLY)
                    return 0;

                  return 1;
                }

              /* mini-snaps: write only those fields marked by either SN_MINI or SN_MINI_ONLY */
              if(DumpFlag == 3)
                {
                  if(IO_Fields[f].snap_type == SN_MINI || IO_Fields[f].snap_type == SN_MINI_ONLY)
                    return 1;

                  if(IO_Fields[f].typelist == BHS_ONLY)
                    return 1;  // temporarily hard-coded that all BH fields are included in mini-snaps

                  return 0;  // specifically do not include any other fields in mini-snaps
                }
            }
          return 0;
        }
    }

  return 0; /* default: not present */
}

/*! \brief This function associates a short 4-character block name with each
 *         block number.
 *
 *   This is stored in front of each block for snapshot FileFormat=2.
 *
 *  \param[in] blocknr ID of the output block (i.e. position, velocities...).
 *  \param[in] label string containing the dataset name.
 *
 *  \return void
 */
void get_Tab_IO_Label(enum iofields blocknr, char *label)
{
  int f;
  for(f = 0; f < N_IO_Fields; f++)
    {
      if(IO_Fields[f].field == blocknr)
        {
          strncpy(label, IO_Fields[f].label, 4);
          return;
        }
    }

  terminate("error invalid field");
}

/*! \brief This function associates a dataset name with each block number.
 *
 *   This is needed to name the dataset if the output is written in HDF5
 *   format.
 *
 *  \param[in] blocknr ID of the output block (i.e. position, velocities...).
 *  \param[in] buf String containing the dataset name.
 *
 *  \return void
 */
void get_dataset_name(enum iofields blocknr, char *buf)
{
  int f;
  for(f = 0; f < N_IO_Fields; f++)
    {
      if(IO_Fields[f].field == blocknr)
        {
          strcpy(buf, IO_Fields[f].datasetname);
          return;
        }
    }

  terminate("error invalid field");
}

/*! \brief Actually write the snapshot file to the disk.
 *
 *  This function writes a snapshot file containing the data from processors
 *  'writeTask' to 'lastTask'. 'writeTask' is the one that actually writes.
 *  Each snapshot file contains a header and cell/particle details. The
 *  output fields for each particle type depend on included physics
 *  and compile-time flags.
 *
 *  \param[in] fname String containing the file name.
 *  \param[in] writeTask The rank of the task in a writing group that which
 *             is responsible for the output operations.
 *  \param[in] lastTask The rank of the last task in a writing group.
 *  \param[in] subbox_flag If greater than 0 instructs the code to output
 *             only a subset of the whole domain.
 *
 *  \return void
 */
void write_file(char *fname, int writeTask, int lastTask, int subbox_flag)
{
  int type, bytes_per_blockelement, npart, nextblock, typelist[NTYPES];
  int n_for_this_task, n, p, pc, offset = 0, task;
  int blockmaxlen, ntot_type[NTYPES], nn[NTYPES];
  enum iofields blocknr;
  char label[8];
  int bnr;
  int blksize;
  MPI_Status status;
  FILE *fd  = 0;
  int pcsum = 0;

#ifdef HAVE_HDF5
  hid_t hdf5_file = 0, hdf5_grp[NTYPES], hdf5_headergrp = 0, hdf5_dataspace_memory;
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  hsize_t dims[2], count[2], start[2];
  int rank = 0;
  char buf[500];
#ifdef HDF5_FILTERS
  hid_t hdf5_properties;
#endif /* #ifdef HDF5_FILTERS */
  hid_t hdf5_paramsgrp = 0;
  hid_t hdf5_configgrp = 0;
#endif /* #ifdef HAVE_HDF5 */

#define SKIP                                 \
  {                                          \
    my_fwrite(&blksize, sizeof(int), 1, fd); \
  }

#ifdef TOLERATE_WRITE_ERROR
  for(int try_io = 0; try_io < 2; try_io++)
    {
      WriteErrorFlag = 0;
#ifdef HAVE_HDF5
      H5Eget_current_stack(); /* clears current error stack */
#endif                        /* #ifdef HAVE_HDF5 */
#endif                        /* #ifdef TOLERATE_WRITE_ERROR */

      /* determine particle numbers of each type in file */
      if(ThisTask == writeTask)
        {
          for(n = 0; n < NTYPES; n++)
            ntot_type[n] = n_type[n];

          for(task = writeTask + 1; task <= lastTask; task++)
            {
              MPI_Recv(&nn[0], NTYPES, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
              for(n = 0; n < NTYPES; n++)
                ntot_type[n] += nn[n];
            }

          for(task = writeTask + 1; task <= lastTask; task++)
            MPI_Send(&ntot_type[0], NTYPES, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
        }
      else
        {
          MPI_Send(&n_type[0], NTYPES, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
          MPI_Recv(&ntot_type[0], NTYPES, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
        }

      /* fill file header */
      for(n = 0; n < NTYPES; n++)
        {
          header.npart[n]              = ntot_type[n];
          header.npartTotal[n]         = (unsigned int)ntot_type_all[n];
          header.npartTotalHighWord[n] = (unsigned int)(ntot_type_all[n] >> 32);
        }

      for(n = 0; n < NTYPES; n++)
        header.mass[n] = All.MassTable[n];

      header.time = All.Time;

      if(All.ComovingIntegrationOn)
        header.redshift = 1.0 / All.Time - 1;
      else
        header.redshift = 0;

      header.flag_sfr        = 0;
      header.flag_feedback   = 0;
      header.flag_cooling    = 0;
      header.flag_stellarage = 0;
      header.flag_metals     = 0;

      header.flag_tracer_field = 0;

#ifdef COOLING
      header.flag_cooling = 1;
#endif /* #ifdef COOLING */

#ifdef USE_SFR
      header.flag_sfr      = 1;
      header.flag_feedback = 1;
#endif /* #ifdef USE_SFR */

      header.num_files   = All.NumFilesPerSnapshot;
      header.BoxSize     = All.BoxSize;
      header.Omega0      = All.Omega0;
      header.OmegaLambda = All.OmegaLambda;
      header.HubbleParam = All.HubbleParam;

#ifdef OUTPUT_IN_DOUBLEPRECISION
      header.flag_doubleprecision = 1;
#else  /* #ifdef OUTPUT_IN_DOUBLEPRECISION */
  header.flag_doubleprecision = 0;
#endif /* #ifdef OUTPUT_IN_DOUBLEPRECISION #else */

      /* open file and write header */

      if(ThisTask == writeTask)
        {
          if(All.SnapFormat == 3)
            {
#ifdef HAVE_HDF5
              sprintf(buf, "%s.hdf5", fname);
              hdf5_file = my_H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

              hdf5_headergrp = my_H5Gcreate(hdf5_file, "/Header", 0);

              for(type = 0; type < NTYPES; type++)
                {
                  if(header.npart[type] > 0)
                    {
                      sprintf(buf, "/PartType%d", type);
                      hdf5_grp[type] = my_H5Gcreate(hdf5_file, buf, 0);
                    }
                }

              write_header_attributes_in_hdf5(hdf5_headergrp);

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

              if(All.SnapFormat == 2)
                {
                  blksize = sizeof(int) + 4 * sizeof(char);
                  SKIP;
                  my_fwrite((void *)"HEAD", sizeof(char), 4, fd);
                  nextblock = sizeof(header) + 2 * sizeof(int);
                  my_fwrite(&nextblock, sizeof(int), 1, fd);
                  SKIP;
                }

              blksize = sizeof(header);
              SKIP;
              my_fwrite(&header, sizeof(header), 1, fd);
              SKIP;
            }
        }

      for(bnr = 0; bnr < 1000; bnr++)
        {
          blocknr = (enum iofields)bnr;

          if(blocknr == IO_LASTENTRY)
            break;

          if(blockpresent(blocknr, 1))
            {
              bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 0);

              blockmaxlen = (int)(COMMBUFFERSIZE / bytes_per_blockelement);

              npart = get_particles_in_block(blocknr, &typelist[0]);

              if(npart > 0)
                {
                  if(ThisTask == 0)
                    {
                      char buf[1000];

                      get_dataset_name(blocknr, buf);
                      if(subbox_flag == 0)
                        printf("writing block %d (%s)...\n", blocknr, buf);
                    }

                  if(ThisTask == writeTask)
                    {
                      if(All.SnapFormat == 1 || All.SnapFormat == 2)
                        {
                          if(All.SnapFormat == 2)
                            {
                              blksize = sizeof(int) + 4 * sizeof(char);
                              SKIP;
                              get_Tab_IO_Label(blocknr, label);
                              my_fwrite(label, sizeof(char), 4, fd);
                              nextblock = npart * bytes_per_blockelement + 2 * sizeof(int);
                              my_fwrite(&nextblock, sizeof(int), 1, fd);
                              SKIP;
                            }

                          blksize = npart * bytes_per_blockelement;
                          SKIP;
                        }
                    }

                  for(type = 0; type < NTYPES; type++)
                    {
                      if(typelist[type])
                        {
#ifdef HAVE_HDF5
                          if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
                            {
                              switch(get_datatype_in_block(blocknr, 0))
                                {
                                  case FILE_INT:
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_UINT);
                                    break;
                                  case FILE_MY_IO_FLOAT:
#ifdef OUTPUT_IN_DOUBLEPRECISION
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_DOUBLE);
#else  /* #ifdef OUTPUT_IN_DOUBLEPRECISION */
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_FLOAT);
#endif /* #ifdef OUTPUT_IN_DOUBLEPRECISION #else */
                                    break;
                                  case FILE_MY_ID_TYPE:
#ifdef LONGIDS
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_UINT64);
#else  /* #ifdef LONGIDS */
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_UINT32);
#endif /* #ifdef LONGIDS #else */
                                    break;
                                  case FILE_DOUBLE:
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_DOUBLE);
                                    break;
                                  case FILE_FLOAT:
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_FLOAT);
                                    break;
                                }

                              dims[0] = header.npart[type];
                              dims[1] = get_values_per_blockelement(blocknr);
                              if(dims[1] == 1)
                                rank = 1;
                              else
                                rank = 2;

                              get_dataset_name(blocknr, buf);

                              hdf5_dataspace_in_file = my_H5Screate_simple(rank, dims, NULL);
#ifdef HDF5_FILTERS
                              hdf5_properties = my_H5Pcreate(H5P_DATASET_CREATE);
                              my_H5Pset_chunk(hdf5_properties, rank, dims); /* set chunk size */
                              my_H5Pset_shuffle(hdf5_properties);           /* reshuffle bytes to get better compression ratio */
                              my_H5Pset_deflate(hdf5_properties, 9);        /* gzip compression level 9 */
                              my_H5Pset_fletcher32(hdf5_properties);        /* Fletcher32 checksum on dataset */

                              if(my_H5Pall_filters_avail(hdf5_properties))
                                hdf5_dataset =
                                    my_H5Dcreate(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file, hdf5_properties);
                              else
                                {
                                  printf("HDF5_FILTERS: Warning selected filters not available! Writing data without filters! \n");
                                  myflush(stdout);
                                  hdf5_dataset = my_H5Dcreate(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);
                                }
#else  /* #ifdef HDF5_FILTERS */
                              hdf5_dataset = my_H5Dcreate(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);
#endif /* #ifdef HDF5_FILTERS #else */
                              write_dataset_attributes(hdf5_dataset, blocknr);
                            }
#endif /* #ifdef HAVE_HDF5 */

                          pcsum               = 0;
                          int remaining_space = blockmaxlen;
                          int bufferstart     = 0;

                          for(task = writeTask, offset = 0; task <= lastTask; task++)
                            {
                              if(task == ThisTask)
                                {
                                  n_for_this_task = n_type[type];

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

                                  if(pc > remaining_space)
                                    pc = remaining_space;

                                  void *buffer = (void *)((char *)CommBuffer + bufferstart * bytes_per_blockelement);

                                  if(ThisTask == task)
                                    fill_write_buffer(buffer, blocknr, &offset, pc, type, subbox_flag);

                                  if(ThisTask == writeTask && task != writeTask)
                                    MPI_Recv(buffer, bytes_per_blockelement * pc, MPI_BYTE, task, TAG_PDATA, MPI_COMM_WORLD, &status);

                                  if(ThisTask != writeTask && task == ThisTask)
                                    MPI_Ssend(buffer, bytes_per_blockelement * pc, MPI_BYTE, writeTask, TAG_PDATA, MPI_COMM_WORLD);

                                  remaining_space -= pc;
                                  bufferstart += pc;

                                  if(remaining_space == 0)
                                    {
                                      /* write stuff (number of elements equal to bufferstart) */
                                      if(ThisTask == writeTask)
                                        {
                                          if(All.SnapFormat == 3)
                                            {
#ifdef HAVE_HDF5
                                              start[0] = pcsum;
                                              start[1] = 0;

                                              count[0] = bufferstart;
                                              count[1] = get_values_per_blockelement(blocknr);

                                              my_H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

                                              dims[0]               = bufferstart;
                                              dims[1]               = get_values_per_blockelement(blocknr);
                                              hdf5_dataspace_memory = my_H5Screate_simple(rank, dims, NULL);

                                              my_H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory, hdf5_dataspace_in_file,
                                                          H5P_DEFAULT, CommBuffer, buf);

                                              my_H5Sclose(hdf5_dataspace_memory, H5S_SIMPLE);
#endif /* #ifdef HAVE_HDF5 */
                                            }
                                          else
                                            {
                                              my_fwrite(CommBuffer, bytes_per_blockelement, bufferstart, fd);
                                            }
                                        }

                                      pcsum += bufferstart;
                                      remaining_space = blockmaxlen;
                                      bufferstart     = 0;
                                    }

                                  n_for_this_task -= pc;
                                }
                            }

                          if(bufferstart > 0)
                            {
                              /* write remaining stuff (number of elements equal to bufferstart) */
                              if(ThisTask == writeTask)
                                {
                                  if(All.SnapFormat == 3)
                                    {
#ifdef HAVE_HDF5
                                      start[0] = pcsum;
                                      start[1] = 0;

                                      count[0] = bufferstart;
                                      count[1] = get_values_per_blockelement(blocknr);

                                      my_H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

                                      dims[0]               = bufferstart;
                                      dims[1]               = get_values_per_blockelement(blocknr);
                                      hdf5_dataspace_memory = my_H5Screate_simple(rank, dims, NULL);

                                      my_H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory, hdf5_dataspace_in_file,
                                                  H5P_DEFAULT, CommBuffer, buf);

                                      my_H5Sclose(hdf5_dataspace_memory, H5S_SIMPLE);
#endif /* #ifdef HAVE_HDF5 */
                                    }
                                  else
                                    {
                                      my_fwrite(CommBuffer, bytes_per_blockelement, bufferstart, fd);
                                    }
                                }

                              pcsum += bufferstart;
                              remaining_space = blockmaxlen;
                              bufferstart     = 0;
                            }

#ifdef HAVE_HDF5
                          if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
                            {
                              if(All.SnapFormat == 3)
                                {
                                  my_H5Dclose(hdf5_dataset, buf);
#ifdef HDF5_FILTERS
                                  my_H5Pclose(hdf5_properties);
#endif /* #ifdef HDF5_FILTERS */
                                  my_H5Sclose(hdf5_dataspace_in_file, H5S_SIMPLE);
                                  my_H5Tclose(hdf5_datatype);
                                }
                            }
#endif /* #ifdef HAVE_HDF5 */
                        }
                    }

                  if(ThisTask == writeTask)
                    {
                      if(All.SnapFormat == 1 || All.SnapFormat == 2)
                        SKIP;
                    }
                }

#ifdef TOLERATE_WRITE_ERROR
              if(ThisTask == writeTask)
                {
                  for(int p = writeTask; p <= lastTask; p++)
                    if(p != ThisTask)
                      MPI_Send(&WriteErrorFlag, 1, MPI_INT, p, TAG_KEY, MPI_COMM_WORLD);
                }
              else
                MPI_Recv(&WriteErrorFlag, 1, MPI_INT, writeTask, TAG_KEY, MPI_COMM_WORLD, &status);
#endif /* #ifdef TOLERATE_WRITE_ERROR */
            }

#ifdef TOLERATE_WRITE_ERROR
          if(WriteErrorFlag) /* don't write further blocks in this case */
            break;
#endif /* #ifdef TOLERATE_WRITE_ERROR */
        }

      if(ThisTask == writeTask)
        {
          if(All.SnapFormat == 3)
            {
#ifdef HAVE_HDF5
              for(type = NTYPES - 1; type >= 0; type--)
                if(header.npart[type] > 0)
                  my_H5Gclose(hdf5_grp[type], buf);
              my_H5Gclose(hdf5_headergrp, "/Header");
              my_H5Gclose(hdf5_paramsgrp, "/Parameters");
              my_H5Gclose(hdf5_configgrp, "/Config");

              sprintf(buf, "%s.hdf5", fname);
              my_H5Fclose(hdf5_file, buf);
#endif /* #ifdef HAVE_HDF5 */
            }
          else
            fclose(fd);
        }

#ifdef TOLERATE_WRITE_ERROR
      if(WriteErrorFlag == 0)
        break;

      if(try_io == 0)
        {
          if(ThisTask == writeTask)
            {
              printf(
                  "TOLERATE_WRITE_ERROR: Try to write to alternative file: masterTask=%d  lastTask=%d  try_io=%d "
                  "alternative-filename='%s'\n",
                  writeTask, lastTask, try_io, alternative_fname);
              myflush(stdout);
            }
          fname = alternative_fname; /* try on a different output directory */
        }
      else
        {
          terminate("TOLERATE_WRITE_ERROR: Second try with alternative file failed too.\n");
        }
    }
#endif /* #ifdef TOLERATE_WRITE_ERROR */
}

#ifdef HAVE_HDF5
/*! \brief Write the fields contained in the header group of the HDF5 snapshot
 *         file.
 *
 *  This function stores the fields of the structure io_header as attributes
 *  belonging to the header group of the HDF5 file.
 *
 *  \param[in] handle A handle for the header group.
 *
 *  \return void
 */
void write_header_attributes_in_hdf5(hid_t handle)
{
  hsize_t adim[1] = {NTYPES};
  hid_t hdf5_dataspace, hdf5_attribute;

  hdf5_dataspace = my_H5Screate(H5S_SIMPLE);
  my_H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL, "NumPart_ThisFile");
  hdf5_attribute = my_H5Acreate(handle, "NumPart_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, header.npart, "NumPart_ThisFile");
  my_H5Aclose(hdf5_attribute, "NumPart_ThisFile");
  my_H5Sclose(hdf5_dataspace, H5S_SIMPLE);

  hdf5_dataspace = my_H5Screate(H5S_SIMPLE);
  my_H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL, "NumPart_Total");
  hdf5_attribute = my_H5Acreate(handle, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal, "NumPart_Total");
  my_H5Aclose(hdf5_attribute, "NumPart_Total");
  my_H5Sclose(hdf5_dataspace, H5S_SIMPLE);

  hdf5_dataspace = my_H5Screate(H5S_SIMPLE);
  my_H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL, "NumPart_Total_HighWord");
  hdf5_attribute = my_H5Acreate(handle, "NumPart_Total_HighWord", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord, "NumPart_Total_HighWord");
  my_H5Aclose(hdf5_attribute, "NumPart_Total_HighWord");
  my_H5Sclose(hdf5_dataspace, H5S_SIMPLE);

  hdf5_dataspace = my_H5Screate(H5S_SIMPLE);
  my_H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL, "MassTable");
  hdf5_attribute = my_H5Acreate(handle, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass, "MassTable");
  my_H5Aclose(hdf5_attribute, "MassTable");
  my_H5Sclose(hdf5_dataspace, H5S_SIMPLE);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time, "Time");
  my_H5Aclose(hdf5_attribute, "Time");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift, "Redshift");
  my_H5Aclose(hdf5_attribute, "Redshift");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize, "BoxSize");
  my_H5Aclose(hdf5_attribute, "BoxSize");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.num_files, "NumFilesPerSnapshot");
  my_H5Aclose(hdf5_attribute, "NumFilesPerSnapshot");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0, "Omega0");
  my_H5Aclose(hdf5_attribute, "Omega0");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda, "OmegaLambda");
  my_H5Aclose(hdf5_attribute, "OmegaLambda");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "OmegaBaryon", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.OmegaBaryon, "OmegaBaryon");
  my_H5Aclose(hdf5_attribute, "OmegaBaryon");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam, "HubbleParam");
  my_H5Aclose(hdf5_attribute, "HubbleParam");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Flag_Sfr", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr, "Flag_Sfr");
  my_H5Aclose(hdf5_attribute, "Flag_Sfr");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Flag_Cooling", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling, "Flag_Cooling");
  my_H5Aclose(hdf5_attribute, "Flag_Cooling");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Flag_StellarAge", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage, "Flag_StellarAge");
  my_H5Aclose(hdf5_attribute, "Flag_StellarAge");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Flag_Metals", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals, "Flag_Metals");
  my_H5Aclose(hdf5_attribute, "Flag_Metals");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Flag_Feedback", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback, "Flag_Feedback");
  my_H5Aclose(hdf5_attribute, "Flag_Feedback");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Flag_DoublePrecision", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_doubleprecision, "Flag_DoublePrecision");
  my_H5Aclose(hdf5_attribute, "Flag_DoublePrecision");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Composition_vector_length", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.composition_vector_length, "Composition_vector_length");
  my_H5Aclose(hdf5_attribute, "Composition_vector_length");
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

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "UnitLength_in_cm", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.UnitLength_in_cm, "UnitLength_in_cm");
  my_H5Aclose(hdf5_attribute, "UnitLength_in_cm");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "UnitMass_in_g", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.UnitMass_in_g, "UnitMass_in_g");
  my_H5Aclose(hdf5_attribute, "UnitMass_in_g");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "UnitVelocity_in_cm_per_s", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.UnitVelocity_in_cm_per_s, "UnitVelocity_in_cm_per_s");
  my_H5Aclose(hdf5_attribute, "UnitVelocity_in_cm_per_s");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
}

/*! \brief Write the parameters read from the parameter file in the HDF5
 *         snapshot file.
 *
 *  This function stores the parameter io_header as attributes belonging
 *  to the parameter group of the HDF5 file.
 *
 *  \param[in] handle A handle for the parameter group.
 *
 *  \return void
 */
void write_parameters_attributes_in_hdf5(hid_t handle)
{
  hid_t hdf5_dataspace, hdf5_attribute, atype = my_H5Tcopy(H5T_C_S1);
  int i = 0;

  my_H5Tset_size(atype, MAXLEN_PARAM_VALUE);

  for(i = 0; i < All.NParameters; i++)
    {
      switch(ParametersType[i])
        {
          case 1:  // REAL
            hdf5_dataspace = my_H5Screate(H5S_SCALAR);
            hdf5_attribute = my_H5Acreate(handle, Parameters[i], H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
            my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, ParametersValue[i], Parameters[i]);
            my_H5Aclose(hdf5_attribute, Parameters[i]);
            my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
            break;
          case 2:  // STRING
            hdf5_dataspace = my_H5Screate(H5S_SCALAR);
            hdf5_attribute = my_H5Acreate(handle, Parameters[i], atype, hdf5_dataspace, H5P_DEFAULT);
            my_H5Awrite(hdf5_attribute, atype, ParametersValue[i], Parameters[i]);
            my_H5Aclose(hdf5_attribute, Parameters[i]);
            my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
            break;
          case 3:  // INT
            hdf5_dataspace = my_H5Screate(H5S_SCALAR);
            hdf5_attribute = my_H5Acreate(handle, Parameters[i], H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
            my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, ParametersValue[i], Parameters[i]);
            my_H5Aclose(hdf5_attribute, Parameters[i]);
            my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
            break;
        }
    }

  my_H5Tclose(atype);
}

/*! \brief A simple error handler for HDF5.
 *
 *  This function terminates the run or if write errors are tolerated, calls
 *  the write_error() function to print information about the error and returns
 *  a positive integer to allow the repetition of the write operation
 *  (see also the HDF5 documentation).
 *
 *  \param[in] unused The parameter is not used, but it is necessary for
 *             compatibility with the HDF5 library.
 *
 *  \return 1 if the write error is tolerated, otherwise the run is terminated.
 */
herr_t my_hdf5_error_handler(void *unused)
{
#ifdef TOLERATE_WRITE_ERROR
  if(FlagNyt == 0)
    write_error(2, 0, 0);
  return 1;
#else
  return 0;
#endif
}

/*! \brief Write attributes to dataset, scaling with a and h (cosmological)
 *         and units.
 *
 *  Only for hdf5 output.
 *
 *  \param[in] hdf5_dataset Dataset identifier.
 *  \param[in] blocknumber Number of field which is written.
 *
 *  \return void
 */
void write_dataset_attributes(hid_t hdf5_dataset, enum iofields blocknr)
{
  int ind = -1;

  for(int f = 0; f < N_IO_Fields; f++)
    {
      if(IO_Fields[f].field == blocknr)
        {
          ind = f;
          break;
        }
    }

  if(ind < 0)
    {
      return;
    }

  if(IO_Fields[ind].hasunit == 0)
    return;

  if(All.ComovingIntegrationOn)
    {
      hid_t hdf5_dataspace = my_H5Screate(H5S_SCALAR);
      hid_t hdf5_attribute = my_H5Acreate(hdf5_dataset, "a_scaling", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
      my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &IO_Fields[ind].a, "a_scaling");
      my_H5Aclose(hdf5_attribute, "a_scaling");
      my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

      hdf5_dataspace = my_H5Screate(H5S_SCALAR);
      hdf5_attribute = my_H5Acreate(hdf5_dataset, "h_scaling", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
      my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &IO_Fields[ind].h, "h_scaling");
      my_H5Aclose(hdf5_attribute, "h_scaling");
      my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
    }
  else
    {
      double zero          = 0;
      hid_t hdf5_dataspace = my_H5Screate(H5S_SCALAR);
      hid_t hdf5_attribute = my_H5Acreate(hdf5_dataset, "a_scaling", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
      my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &zero, "a_scaling");
      my_H5Aclose(hdf5_attribute, "a_scaling");
      my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

      hdf5_dataspace = my_H5Screate(H5S_SCALAR);
      hdf5_attribute = my_H5Acreate(hdf5_dataset, "h_scaling", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
      my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &zero, "h_scaling");
      my_H5Aclose(hdf5_attribute, "h_scaling");
      my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
    }

  hid_t hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hid_t hdf5_attribute = my_H5Acreate(hdf5_dataset, "length_scaling", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &IO_Fields[ind].L, "length_scaling");
  my_H5Aclose(hdf5_attribute, "length_scaling");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(hdf5_dataset, "mass_scaling", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &IO_Fields[ind].M, "mass_scaling");
  my_H5Aclose(hdf5_attribute, "mass_scaling");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(hdf5_dataset, "velocity_scaling", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &IO_Fields[ind].V, "velocity_scaling");
  my_H5Aclose(hdf5_attribute, "velocity_scaling");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(hdf5_dataset, "to_cgs", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &IO_Fields[ind].c, "to_cgs");
  my_H5Aclose(hdf5_attribute, "to_cgs");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
}
#endif /* #ifdef HAVE_HDF5 */

#ifdef OUTPUT_XDMF
/*! \brief Outputs a xdmf file corresponding to this snapshot.
 *
 *  This xdmf file can be used to load the snapshot into programs like visit.
 *  This option only works with output format 3 (hdf5).
 *
 *  \param[in] fname Name of the snapshot.
 *
 *  \return void
 */
static void write_xdmf(char *fname)
{
  FILE *f;
  char buf[256], buf2[256];
  int i;
  int npresent[NTYPES];

  for(i = 0; i < NTYPES; i++)
    npresent[i] = 0;

#ifdef OUTPUT_IN_DOUBLEPRECISION
  int prec = 8;
#else  /* #ifdef OUTPUT_IN_DOUBLEPRECISION */
  int prec = 4;
#endif /* #ifdef OUTPUT_IN_DOUBLEPRECISION */

  sprintf(buf, "%s.xmf", fname);
  f = fopen(buf, "w");

  fprintf(f, "<?xml version=\"1.0\" ?>\n");
  fprintf(f, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
  fprintf(f, "<Xdmf Version=\"2.0\">\n");
  fprintf(f, " <Domain>");

  /* hdf5 file path relative to xmf file, uses basename function of libgen.h,
   * i.e. POSIX version of basename() */
  sprintf(buf, "./%s.hdf5", basename(fname));
  int type = 0;
  for(; type < NTYPES; type++)
    {
      int bnr;

      for(bnr = 0; bnr < 1000; bnr++)
        {
          enum iofields i = (enum iofields)bnr;

          if(i == IO_LASTENTRY)
            break;

          if(blockpresent(i, 1))
            {
              // get_particles_in_block(i, ntypes);

              if(header.npart[type] > 0)
                {
                  if(i == IO_POS)
                    {
                      fprintf(f, "  <Grid Name=\"PartType%d\" GridType=\"Uniform\">\n", type);
                      fprintf(f, "   <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%d\"/>\n", header.npart[type]);
                      fprintf(f, "   <Geometry GeometryType=\"XYZ\">\n");
                      fprintf(f, "    <DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",
                              header.npart[type], prec);
                      fprintf(f, "     %s:/PartType0/Coordinates\n", buf);
                      fprintf(f, "    </DataItem>\n");
                      fprintf(f, "   </Geometry>\n");

                      npresent[type] = 1;
                    }
                  else
                    {
                      int dim   = get_values_per_blockelement(i);
                      int dtype = get_datatype_in_block(i, 0);
                      get_dataset_name(i, buf2);

                      if(dim == 1 || dim == 3)
                        {
                          if(dtype == 1)
                            {
                              if(dim == 1)
                                {
                                  fprintf(f, "   <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n", buf2);
                                  fprintf(f, "    <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",
                                          header.npart[type], prec);
                                }
                              else
                                {
                                  fprintf(f, "   <Attribute Name=\"%s\" AttributeType=\"Vector\" Center=\"Node\">\n", buf2);
                                  fprintf(f,
                                          "    <DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",
                                          header.npart[type], prec);
                                }

                              fprintf(f, "     %s:/PartType%d/%s\n", buf, type, buf2);
                              fprintf(f, "    </DataItem>\n");
                              fprintf(f, "   </Attribute>\n");
                            }
                        }
                    }
                }
            }
        }
      if(npresent[type] == 1)
        {
          fprintf(f, "  </Grid>\n");
        }
    }

  fprintf(f, " </Domain>\n");
  fprintf(f, "</Xdmf>");

  fclose(f);
}
#endif /* #ifdef OUTPUT_XDMF */

/*! \brief  A wrapper for the fwrite() function.
 *
 *  This catches I/O errors occuring for fwrite(). In this case we
 *  better stop. If stream is null, no attempt at writing is done.
 *
 *  \param[in] ptr Pointer to the beginning of data to write.
 *  \param[in] size Size in bytes of a single data element.
 *  \param[in] nmemb Number of elements to be written.
 *  \param[in] stream Pointer to the output stream.
 *
 *  \return Number of elements written to stream.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nwritten;

#ifdef TOLERATE_WRITE_ERROR
  if(WriteErrorFlag)
    return 0;
#endif /* #ifdef TOLERATE_WRITE_ERROR */

  if(!stream)
    return 0;

  if(size * nmemb > 0)
    {
      if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
        {
#ifdef TOLERATE_WRITE_ERROR
          write_error(0, nwritten, nmemb);
#else  /* #ifdef TOLERATE_WRITE_ERROR */
          printf("I/O error (fwrite) on task=%d has occured: %s\n", ThisTask, strerror(errno));
          myflush(stdout);
          terminate("write error");
#endif /* #ifdef TOLERATE_WRITE_ERROR #else */
        }
    }
  else
    nwritten = 0;

#ifdef TOLERATE_WRITE_ERROR
  if(ferror(stream))
    write_error(1, nwritten, nmemb);
#endif /* #ifdef TOLERATE_WRITE_ERROR */

  return nwritten;
}

/*! \brief  A wrapper for the fread() function.
 *
 *  This catches I/O errors occuring for fread(). In this case we
 *  better stop. If stream is null, no attempt at readingis done.
 *
 *  \param[out] ptr Pointer to the beginning of memory location where to
 *              store data.
 *  \param[in] size Size in bytes of a single data element.
 *  \param[in] nmemb Number of elements to be read.
 *  \param[in] stream Pointer to the input stream.
 *
 *  \return Number of elements read from stream.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nread;

  if(!stream)
    return 0;

  if(size * nmemb > 0)
    {
      if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
        {
          if(feof(stream))
            printf("I/O error (fread) on task=%d has occured: end of file\n", ThisTask);
          else
            printf("I/O error (fread) on task=%d has occured: %s\n", ThisTask, strerror(errno));
          myflush(stdout);
          terminate("read error");
        }
    }
  else
    nread = 0;

  return nread;
}

/*! \brief A wrapper for the printf() function.
 *
 *  This function has the same functionalities of the standard printf()
 *  function. However, data is written to the standard output only for
 *  the task with rank 0.
 *
 *  \param[in] fmt String that contains format arguments.
 *
 *  \return void
 */
void mpi_printf(const char *fmt, ...)
{
  if(ThisTask == 0)
    {
      va_list l;
      va_start(l, fmt);
      vprintf(fmt, l);
      myflush(stdout);
      va_end(l);
    }
}

/*! \brief A wrapper for the fprintf() function.
 *
 *  This function has the same functionalities of the standard fprintf()
 *  function. However, data is written to the standard output only for
 *  the task with rank 0.
 *
 *  \param[in] fmt String that contains format arguments.
 *
 *  \return void
 */
void mpi_fprintf(FILE *stream, const char *fmt, ...)
{
  if(ThisTask == 0)
    {
      va_list l;
      va_start(l, fmt);
      vfprintf(stream, fmt, l);
      myflush(stream);
      va_end(l);
    }
}

/*! \brief A function for printing debug information in parallel.
 *
 *  This function works like printf, however it takes care
 *  that the output is contigous in the stdout from task 0 to task NTask-1.
 *  Run this debug function only in code parts which all tasks reach.
 *
 *
 *  \param[in] fmt String that contains format arguments.
 *
 *  \return void
 */
void mpi_printf_each(const char *fmt, ...)
{
  char buffer[2048];

  va_list l;
  va_start(l, fmt);
  vsprintf(buffer, fmt, l);
  va_end(l);

  if(ThisTask == 0)
    {
      // print own message
      printf("%s", buffer);

      // print message from other tasks
      unsigned int i;

      for(i = 1; i < NTask; i++)
        {
          MPI_Recv(buffer, 2048, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          printf("%s", buffer);
        }
    }

  else
    {
      MPI_Send(buffer, strlen(buffer) + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
}

/*! \brief Opens the requested file name and returns the file descriptor.
 *
 *  If opening fails, an error is printed and the file descriptor is
 *  null.
 *
 *  \param[in] fnam The file name.
 *
 *  \return A file descriptor to the file.
 */
FILE *open_file(char *fnam)
{
  FILE *fd;

  if(!(fd = fopen(fnam, "w")))
    {
      printf("can't open file `%s' for writing.\n", fnam);
    }
  return fd;
}
