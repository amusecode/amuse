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
 * \file        src/io/restart.c
 * \date        05/2018
 * \brief       Handling of the loading/writing of restart files.
 * \details     contains functions:
 *                void loadrestart(void)
 *                void reread_params_after_loading_restart(void)
 *                static int compare_seq_data(const void *a, const void *b)
 *                static void create_restartfiles_dir()
 *                static void get_restart_filename(char *buf, int task,
 *                  int modus)
 *                static void backup_restartfiles(int task)
 *                static int get_file_to_check(int task)
 *                static void check_restart_files(char *buf, struct check *ch,
 *                  int *success)
 *                static void send_work_request(int modus, int i)
 *                static void polling(int modus)
 *                static void work_files(int modus)
 *                void restart(int modus)
 *                static void write_or_read_this_processors_restart_file(int
 *                  modus, char *buf, struct check *ch)
 *                static int execute_write_or_read(int modus, char *buf,
 *                  struct check *ch)
 *                static void contents_restart_file(int modus)
 *                void readjust_timebase(double TimeMax_old,
 *                  double TimeMax_new)
 *                void in(int *x, int modus)
 *                void byten(void *x, size_t n, int modus)
 *                void byten_nohash(void *x, size_t n, int modus)
 *                void byten_hash(void *x, size_t n, int modus, int hash)
 *                void allocate_iobuf(void)
 *                void deallocate_iobuf(int modus)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 21.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifndef __USE_GNU
#define _GNU_SOURCE /* needed for USE_DIRECT_IO_FOR_RESTARTS */
#endif              /* #ifndef __USE_GNU */

#include <fcntl.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../debug_md5/Md5.h"
#include "../domain/domain.h"
#include "../mesh/voronoi/voronoi.h"

#define MODUS_WRITE 0
#define MODUS_READ 1
#define MODUS_READCHECK 2
#define MODUS_CHECK 3

/*! \brief Data for scheduling restart file IO.
 */
static struct seq_data
{
  int thistask;
  int rankinnode;
  int thisnode;
} * seq;

/*! \brief Metadata of restart files to be compared to when checking.
 */
static struct check
{
  long long byte_count;
  unsigned char hash[16];
} * checks;

static char *write_success;

static int fdint;

static void in(int *x, int modus);
static void byten(void *x, size_t n, int modus);
static void byten_nohash(void *x, size_t n, int modus);
static void byten_hash(void *x, size_t n, int modus, int hash);
static void write_or_read_this_processors_restart_file(int modus, char *fname, struct check *ch);
static int execute_write_or_read(int modus, char *buf, struct check *ch);
static void contents_restart_file(int modus);

#define MAX_BLOCK_SIZE (32 * 1024 * 1024)

static int PageSize;
static char *iobuf_aligned, *io_buf;
static size_t fillp, iop;
void allocate_iobuf(void);
void deallocate_iobuf(int modus);

static long long byte_count;
static int files_started;
static int files_completed;
static int files_concurrent;
static int files_groups;

static MD5_CTX mysum;

static struct global_data_all_processes all;

/*! \brief This function loads the last restart file.
 *
 *  Some parameters of the parameter file might be changed between restarting.
 *  This function ensures that only the allowed parameters change,
 *  otherwise the old value from the restart file is taken.
 *  If the end time of the simulation changed readjust_timebase() is called in
 *  the end.
 *
 *  \return void
 */
void loadrestart(void)
{
  /* save global variables. (will be read from restart file) */
  all = All;

  /* Read restart files.
     Note: This also resets all variables in the struct `All'. */
  restart(MODUS_READ);

  /* However, during the run, some variables in the parameter
     file are allowed to be changed, if desired. These are copied here. */
  reread_params_after_loading_restart();
}

/*! \brief This function takes from the parameter file values that are allowed
 *         to change after restart.
 *
 *  \return void
 */
void reread_params_after_loading_restart(void)
{
  if(ThisTask == 0 && All.MinSizeTimestep != all.MinSizeTimestep)
    warn("MinSizeTimestep modified from %g to %g while restarting at Time=%g", All.MinSizeTimestep, all.MinSizeTimestep, All.Time);
  All.MinSizeTimestep = all.MinSizeTimestep;
  if(ThisTask == 0 && All.MaxSizeTimestep != all.MaxSizeTimestep)
    warn("MaxSizeTimestep modified from %g to %g while restarting at Time=%g", All.MaxSizeTimestep, all.MaxSizeTimestep, All.Time);
  All.MaxSizeTimestep = all.MaxSizeTimestep;
  if(ThisTask == 0 && All.TimeLimitCPU != all.TimeLimitCPU)
    warn("TimeLimitCPU modified from %g to %g while restarting at Time=%g", All.TimeLimitCPU, all.TimeLimitCPU, All.Time);
  All.TimeLimitCPU = all.TimeLimitCPU;
  if(ThisTask == 0 && All.ResubmitOn != all.ResubmitOn)
    warn("ResubmitOn modified from %d to %d while restarting at Time=%g", All.ResubmitOn, all.ResubmitOn, All.Time);
  All.ResubmitOn = all.ResubmitOn;
  if(ThisTask == 0 && All.TimeBetSnapshot != all.TimeBetSnapshot)
    warn("TimeBetSnapshot modified from %g to %g while restarting at Time=%g", All.TimeBetSnapshot, all.TimeBetSnapshot, All.Time);
  All.TimeBetSnapshot = all.TimeBetSnapshot;
  if(ThisTask == 0 && All.TimeBetStatistics != all.TimeBetStatistics)
    warn("TimeBetStatistics modified from %g to %g while restarting at Time=%g", All.TimeBetStatistics, all.TimeBetStatistics,
         All.Time);
  All.TimeBetStatistics = all.TimeBetStatistics;
  if(ThisTask == 0 && All.CpuTimeBetRestartFile != all.CpuTimeBetRestartFile)
    warn("CpuTimeBetRestartFile modified from %g to %g while restarting at Time=%g", All.CpuTimeBetRestartFile,
         all.CpuTimeBetRestartFile, All.Time);
  All.CpuTimeBetRestartFile = all.CpuTimeBetRestartFile;
  if(ThisTask == 0 && All.ErrTolIntAccuracy != all.ErrTolIntAccuracy)
    warn("ErrTolIntAccuracy modified from %g to %g while restarting at Time=%g", All.ErrTolIntAccuracy, all.ErrTolIntAccuracy,
         All.Time);
  All.ErrTolIntAccuracy = all.ErrTolIntAccuracy;
  if(ThisTask == 0 && All.SnapFormat != all.SnapFormat)
    warn("SnapFormat modified from %d to %d while restarting at Time=%g", All.SnapFormat, all.SnapFormat, All.Time);
  All.SnapFormat = all.SnapFormat;

  if(ThisTask == 0 && All.ErrTolForceAcc != all.ErrTolForceAcc)
    warn("ErrTolForceAcc modified from %g to %g while restarting at Time=%g", All.ErrTolForceAcc, all.ErrTolForceAcc, All.Time);
  All.ErrTolForceAcc = all.ErrTolForceAcc;
  if(ThisTask == 0 && All.TypeOfTimestepCriterion != all.TypeOfTimestepCriterion)
    warn("TypeOfTimestepCriterion modified from %d to %d while restarting at Time=%g", All.TypeOfTimestepCriterion,
         all.TypeOfTimestepCriterion, All.Time);
  All.TypeOfTimestepCriterion = all.TypeOfTimestepCriterion;
  if(ThisTask == 0 && All.TypeOfOpeningCriterion != all.TypeOfOpeningCriterion)
    warn("TypeOfOpeningCriterion modified from %d to %d while restarting at Time=%g", All.TypeOfOpeningCriterion,
         all.TypeOfOpeningCriterion, All.Time);
  All.TypeOfOpeningCriterion = all.TypeOfOpeningCriterion;
  if(ThisTask == 0 && All.NumFilesWrittenInParallel != all.NumFilesWrittenInParallel)
    warn("NumFilesWrittenInParallel modified from %d to %d while restarting at Time=%g", All.NumFilesWrittenInParallel,
         all.NumFilesWrittenInParallel, All.Time);
  All.NumFilesWrittenInParallel = all.NumFilesWrittenInParallel;
  if(ThisTask == 0 && All.NumFilesPerSnapshot != all.NumFilesPerSnapshot)
    warn("NumFilesPerSnapshot modified from %d to %d while restarting at Time=%g", All.NumFilesPerSnapshot, all.NumFilesPerSnapshot,
         All.Time);
  All.NumFilesPerSnapshot = all.NumFilesPerSnapshot;

  if(ThisTask == 0 && All.LimitUBelowThisDensity != all.LimitUBelowThisDensity)
    warn("LimitUBelowThisDensity modified from %g to %g while restarting at Time=%g", All.LimitUBelowThisDensity,
         all.LimitUBelowThisDensity, All.Time);
  All.LimitUBelowThisDensity = all.LimitUBelowThisDensity;
  if(ThisTask == 0 && All.LimitUBelowCertainDensityToThisValue != all.LimitUBelowCertainDensityToThisValue)
    warn("LimitUBelowCertainDensityToThisValue modified from %g to %g while restarting at Time=%g",
         All.LimitUBelowCertainDensityToThisValue, all.LimitUBelowCertainDensityToThisValue, All.Time);
  All.LimitUBelowCertainDensityToThisValue = all.LimitUBelowCertainDensityToThisValue;
  if(ThisTask == 0 && All.MinimumDensityOnStartUp != all.MinimumDensityOnStartUp)
    warn("MinimumDensityOnStartUp modified from %g to %g while restarting at Time=%g", All.MinimumDensityOnStartUp,
         all.MinimumDensityOnStartUp, All.Time);
  All.MinimumDensityOnStartUp = all.MinimumDensityOnStartUp;
  if(ThisTask == 0 && All.MultipleDomains != all.MultipleDomains)
    warn("MultipleDomains modified from %d to %d while restarting at Time=%g", All.MultipleDomains, all.MultipleDomains, All.Time);
  All.MultipleDomains = all.MultipleDomains;
  if(ThisTask == 0 && All.TopNodeFactor != all.TopNodeFactor)
    warn("TopNodeFactor modified from %g to %g while restarting at Time=%g", All.TopNodeFactor, all.TopNodeFactor, All.Time);
  All.TopNodeFactor = all.TopNodeFactor;
  if(ThisTask == 0 && All.ActivePartFracForNewDomainDecomp != all.ActivePartFracForNewDomainDecomp)
    warn("ActivePartFracForNewDomainDecomp modified from %g to %g while restarting at Time=%g", All.ActivePartFracForNewDomainDecomp,
         all.ActivePartFracForNewDomainDecomp, All.Time);
  All.ActivePartFracForNewDomainDecomp = all.ActivePartFracForNewDomainDecomp;
  if(ThisTask == 0 && All.OutputListOn != all.OutputListOn)
    warn("OutputListOn modified from %d to %d while restarting at Time=%g", All.OutputListOn, all.OutputListOn, All.Time);
  All.OutputListOn = all.OutputListOn;
  if(ThisTask == 0 && All.CourantFac != all.CourantFac)
    warn("CourantFac modified from %g to %g while restarting at Time=%g", All.CourantFac, all.CourantFac, All.Time);
  All.CourantFac = all.CourantFac;
#ifdef REGULARIZE_MESH_FACE_ANGLE
  if(ThisTask == 0 && All.CellMaxAngleFactor != all.CellMaxAngleFactor)
    warn("CellMaxAngleFactor modified from %g to %g while restarting at Time=%g", All.CellMaxAngleFactor, all.CellMaxAngleFactor,
         All.Time);
  All.CellMaxAngleFactor = all.CellMaxAngleFactor;
#else  /* #ifdef REGULARIZE_MESH_FACE_ANGLE */
  if(ThisTask == 0 && All.CellShapingFactor != all.CellShapingFactor)
    warn("CellShapingFactor modified from %g to %g while restarting at Time=%g", All.CellShapingFactor, all.CellShapingFactor,
         All.Time);
  All.CellShapingFactor = all.CellShapingFactor;
#endif /* #ifdef REGULARIZE_MESH_FACE_ANGLE #else */
  if(ThisTask == 0 && All.CellShapingSpeed != all.CellShapingSpeed)
    warn("CellShapingSpeed modified from %g to %g while restarting at Time=%g", All.CellShapingSpeed, all.CellShapingSpeed, All.Time);
  All.CellShapingSpeed = all.CellShapingSpeed;

  if(ThisTask == 0 && All.OutputListLength != all.OutputListLength)
    warn("OutputListLength modified from %d to %d while restarting at Time=%g", All.OutputListLength, all.OutputListLength, All.Time);
  All.OutputListLength = all.OutputListLength;
  if(ThisTask == 0 && memcmp(All.OutputListTimes, all.OutputListTimes, sizeof(double) * All.OutputListLength) != 0)
    warn("OutputListTimes modified while restarting at Time=%g", All.Time);
  memcpy(All.OutputListTimes, all.OutputListTimes, sizeof(double) * All.OutputListLength);
  if(ThisTask == 0 && memcmp(All.OutputListFlag, all.OutputListFlag, sizeof(char) * All.OutputListLength) != 0)
    warn("OutputListFlag modified while restarting at Time=%g", All.Time);
  memcpy(All.OutputListFlag, all.OutputListFlag, sizeof(char) * All.OutputListLength);

  if(ThisTask == 0 && strcmp(All.ResubmitCommand, all.ResubmitCommand) != 0)
    warn("ResubmitCommand modified from %s to %s while restarting at Time=%g", All.ResubmitCommand, all.ResubmitCommand, All.Time);
  strcpy(All.ResubmitCommand, all.ResubmitCommand);
  if(ThisTask == 0 && strcmp(All.OutputListFilename, all.OutputListFilename) != 0)
    warn("OutputListFilename modified from %s to %s while restarting at Time=%g", All.OutputListFilename, all.OutputListFilename,
         All.Time);
  strcpy(All.OutputListFilename, all.OutputListFilename);
  if(ThisTask == 0 && strcmp(All.OutputDir, all.OutputDir) != 0)
    warn("OutputDir modified from %s to %s while restarting at Time=%g", All.OutputDir, all.OutputDir, All.Time);
  strcpy(All.OutputDir, all.OutputDir);
  if(ThisTask == 0 && strcmp(All.SnapshotFileBase, all.SnapshotFileBase) != 0)
    warn("SnapshotFileBase modified from %s to %s while restarting at Time=%g", All.SnapshotFileBase, all.SnapshotFileBase, All.Time);
  strcpy(All.SnapshotFileBase, all.SnapshotFileBase);

#ifdef MHD_SEEDFIELD
  if(ThisTask == 0 && All.B_dir != all.B_dir)
    warn("B_dir modified from %d to %d while restarting at Time=%g", All.B_dir, all.B_dir, All.Time);
  All.B_dir = all.B_dir;
  if(ThisTask == 0 && All.B_value != all.B_value)
    warn("B_value modified from %g to %g while restarting at Time=%g", All.B_value, all.B_value, All.Time);
  All.B_value = all.B_value;
#endif /* #ifdef MHD_SEEDFIELD */

  if(All.TimeMax != all.TimeMax)
    {
      if(ThisTask == 0)
        warn("TimeMax modified from %g to %g while restarting at Time=%g", All.TimeMax, all.TimeMax, All.Time);
      readjust_timebase(All.TimeMax, all.TimeMax);
    }
}

/*! \brief Sorting kernel for seq_data strucutre.
 *
 *  Compares (top priority first)
 *    rankinnode
 *    thisnode
 *    thistask
 *
 *  \return (-1,0,1), -1 if a < b.
 */
static int compare_seq_data(const void *a, const void *b)
{
  if(((struct seq_data *)a)->rankinnode < ((struct seq_data *)b)->rankinnode)
    return -1;

  if(((struct seq_data *)a)->rankinnode > ((struct seq_data *)b)->rankinnode)
    return +1;

  if(((struct seq_data *)a)->thisnode < ((struct seq_data *)b)->thisnode)
    return -1;

  if(((struct seq_data *)a)->thisnode > ((struct seq_data *)b)->thisnode)
    return +1;

  if(((struct seq_data *)a)->thistask < ((struct seq_data *)b)->thistask)
    return -1;

  if(((struct seq_data *)a)->thistask > ((struct seq_data *)b)->thistask)
    return +1;

  return 0;
}

/*! \brief Creates the restart file directory with appropriate permissions.
 *
 *  \return void
 */
static void create_restartfiles_dir()
{
  char buf[MAXLEN_PATH];
#ifdef MULTIPLE_RESTARTS
  printf(", All.RestartFileCount=%03d", All.RestartFileCount);
#endif /* #ifdef MULTIPLE_RESTARTS */
  printf(".\n");
  sprintf(buf, "%s/restartfiles", All.OutputDir);
#ifdef MULTIPLE_RESTARTS
  sprintf(buf, "%s/restartfiles_%03d", All.OutputDir, All.RestartFileCount);
#endif /* #ifdef MULTIPLE_RESTARTS */
  mkdir(buf, 02755);

#ifdef TOLERATE_WRITE_ERROR
  sprintf(buf, "%s/restartfiles", AlternativeOutputDir);
  mkdir(buf, 02755);
#endif /* #ifdef TOLERATE_WRITE_ERROR */
}

/*! \brief Sets filename of restart file on local task.
 *
 *  \param[out] buf Buffer to which filename is written.
 *  \param[in] task Task for which restart file should be written.
 *  \param[in] modus Read or write mode flag.
 *
 *  \return void
 */
static void get_restart_filename(char *buf, int task, int modus)
{
  sprintf(buf, "%s/restartfiles/%s.%d", All.OutputDir, "restart", task);

#ifdef MULTIPLE_RESTARTS
  if(modus == MODUS_WRITE)
    sprintf(buf, "%s/restartfiles_%03d/%s.%d", All.OutputDir, All.RestartFileCount++, "restart", task);
  if((modus == MODUS_READ) || (modus == MODUS_READCHECK) || (modus == MODUS_CHECK))
    sprintf(buf, "%s/restartfiles_%03d/%s.%d", All.OutputDir, All.RestartFileCount - 1, "restart", task);
#endif /* #ifdef MULTIPLE_RESTARTS */
}

/*! \brief Renames existing restartfiles to backup-restartfiles.
 *
 *  This way the code ensures that there are two sets of restart-files per
 *  run.
 *
 *  \param[in] task Task for which restart file is renamed.
 *
 *  \return void
 */
static void backup_restartfiles(int task)
{
  char buf[MAXLEN_PATH];

  FILE *fcheck = NULL;
  char buf_bak[MAXLEN_PATH];

  int bak_files_status = 0;

  mpi_printf("RESTART: Backup restart files...\n");
  myflush(stdout);

  get_restart_filename(buf, task, MODUS_READ);

  sprintf(buf_bak, "%s/restartfiles/bak-%s.%d", All.OutputDir, "restart", ThisTask);
  if((fcheck = fopen(buf, "r")))
    {
      fclose(fcheck);

      rename(buf, buf_bak);
      bak_files_status = 1;
    }
#ifdef TOLERATE_WRITE_ERROR
  char alternative_fname[MAXLEN_PATH];
  sprintf(alternative_fname, "%s/restartfiles/%s.%d", AlternativeOutputDir, "restart", ThisTask);
  sprintf(buf_bak, "%s/restartfiles/bak-%s.%d", AlternativeOutputDir, "restart", ThisTask);

  if((fcheck = fopen(alternative_fname, "r")))
    {
      fclose(fcheck);

      rename(alternative_fname, buf_bak);
      bak_files_status = 1;
    }
#endif /* #ifdef TOLERATE_WRITE_ERROR */

  int bak_files_status_sum;
  MPI_Allreduce(&bak_files_status, &bak_files_status_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(bak_files_status_sum != NTask && bak_files_status_sum != 0)
    warn("RESTART: some (%d) restart files were renamed to bak, but some (%d) weren't - something is very possibly wrong!",
         bak_files_status, NTask - bak_files_status);
  if(bak_files_status_sum == NTask)
    mpi_printf("RESTART: done renaming pre-existing restart files to bak files.\n");
  else if(bak_files_status_sum == 0)
    mpi_printf("RESTART: no pre-existing restart files found.\n");

  myflush(stdout);
}

/*! \brief Returns the index of file which is to be checked by local task.
 *
 *  Task gets assigned a restart file to check which it has not written.
 *  This is to ensure that the check is actually read from disk (not from some
 *  buffer).
 *
 *  \param[in] task Local task.
 *
 *  \return File number.
 */
static int get_file_to_check(int task) { return (task + NTask / 2) % NTask; }

/*! \brief Checks restart files via an md5sum.
 *
 *  This is to ensure that they have been written correctly to the file system.
 *
 *  \param[in] buf Filename of restart file.
 *  \param[in] ch Metadata about data to be checked.
 *  \param[out] success Flag whether check was a success.
 *
 *  \return void
 */
static void check_restart_files(char *buf, struct check *ch, int *success)
{
#ifdef USE_DIRECT_IO_FOR_RESTARTS
  struct stat st;
  if(stat(buf, &st) == 0)
    {
      size_t size = st.st_size;
      if(size % PageSize > 0)
        {
          FILE *fd = fopen(buf, "a");
          if(fd)
            {
              size_t n = PageSize - (size % PageSize);
              char *p  = calloc(n, 1);
              if(p == NULL)
                terminate("p == NULL");
              printf("RESTART: Topping of restart file '%s' by %lld bytes\n", buf, (long long)n);
              fwrite(p, n, 1, fd);
              fclose(fd);
              free(p);
            }
          else
            terminate("can't increase length of restart file '%s'", buf);
        }
    }
  else
    terminate("Restart file '%s' not found.\n", buf);
#endif /* #ifdef USE_DIRECT_IO_FOR_RESTARTS */
  int oflag = O_RDONLY;
#ifdef USE_DIRECT_IO_FOR_RESTARTS
  oflag |= O_DIRECT;
#endif /* #ifdef USE_DIRECT_IO_FOR_RESTARTS */

  if((fdint = open(buf, oflag)) < 0)
    terminate("Restart file '%s' not found.\n", buf);

  allocate_iobuf();

  MD5Init(&mysum);

  long long readLen = ch->byte_count;
  while(readLen > 0)
    {
      int readChunk = 1024 * 1024 * 32;
      if(readChunk > readLen)
        readChunk = readLen;

      byten(NULL, readChunk, MODUS_CHECK);
      readLen -= readChunk;
    }

  MD5Final(&mysum);

  unsigned char has_hash[16], written_hash[16];

  for(int k = 0; k < 16; k++)
    has_hash[k] = mysum.digest[k];

  byten_nohash(written_hash, 16, MODUS_READ);

  if(memcmp(has_hash, ch->hash, 16) != 0 || memcmp(has_hash, written_hash, 16) != 0)
    {
      char str_has[48], str_expected[48], str_written[48];
      for(int i = 0; i < 16; i++)
        {
          sprintf(str_has + 2 * i, "%02X", has_hash[i]);
          sprintf(str_expected + 2 * i, "%02X", ch->hash[i]);
          sprintf(str_written + 2 * i, "%02X", written_hash[i]);
        }

      str_has[32] = str_expected[32] = str_written[32] = 0;

      char newname[10000];
      sprintf(newname, "%s-damaged", buf);
      rename(buf, newname);

      terminate("RESTART: file '%s' has MD5 hash of '%s', does not match expected hash '%s' or written hash '%s'.", newname, str_has,
                str_expected, str_written);
      *success = 0;
    }
  else
    {
#ifdef VERBOSE
      char str_has[48], str_expected[48], str_written[48];
      for(int i = 0; i < 16; i++)
        {
          sprintf(str_has + 2 * i, "%02X", has_hash[i]);
          sprintf(str_expected + 2 * i, "%02X", ch->hash[i]);
          sprintf(str_written + 2 * i, "%02X", written_hash[i]);
        }

      str_has[32] = str_expected[32] = str_written[32] = 0;

      printf("RESTART: Task %d: file '%s' has MD5 hash of '%s', does match expected hash '%s' and written hash '%s'.\n", ThisTask, buf,
             str_has, str_expected, str_written);
#endif /* #ifdef VERBOSE */
      *success = 1;
    }
  deallocate_iobuf(MODUS_CHECK);

  close(fdint);
}

/*! \brief Distributes information and meta-data to task that is supposed to
 *         check the restart file which has just been written.
 *
 *  \param[in] modus Write or check mode.
 *  \param[in] i Index in seq array.
 *
 *  \return void
 */
static void send_work_request(int modus, int i)
{
  int type = 0;

  if(modus == MODUS_WRITE)
    {
      if(write_success[seq[i].thistask])
        type = 1;
    }

  if(modus == MODUS_CHECK)
    {
      int task = get_file_to_check(seq[i].thistask);
      if(write_success[task])
        type = 1;
    }

  MPI_Ssend(&type, 1, MPI_INT, seq[i].thistask, TAG_N, MPI_COMM_WORLD);

  if(modus == MODUS_CHECK)
    {
      int task = get_file_to_check(seq[i].thistask);
      if(!write_success[task])
        MPI_Ssend(&checks[task], sizeof(struct check), MPI_BYTE, seq[i].thistask, TAG_N, MPI_COMM_WORLD);
    }
}

/*! \brief Gets work request.
 *
 *  \param[in] modus Write or check files.
 *
 *  \return void
 */
static void polling(int modus)
{
  if(ThisTask == 0)
    if(files_completed < NTask)
      {
        MPI_Status status;
        int flag;

        /* now check for a completion message  */
        MPI_Iprobe(MPI_ANY_SOURCE, TAG_KEY, MPI_COMM_WORLD, &flag, &status);

        if(flag)
          {
            int source = status.MPI_SOURCE;

            if(modus == MODUS_WRITE)
              {
                MPI_Recv(&checks[source], sizeof(struct check), MPI_BYTE, source, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
              }
            else if(modus == MODUS_CHECK)
              {
                int success;
                MPI_Recv(&success, 1, MPI_INT, source, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                write_success[get_file_to_check(source)] = success;
              }
            else
              {
                int dummy;
                MPI_Recv(&dummy, 1, MPI_INT, source, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
              }
            files_completed++;

            if(files_started < NTask)
              {
                if((files_started % files_concurrent) == 0)
                  {
                    if(modus == MODUS_READ)
                      mpi_printf("RESTART: Loading restart files group #%d out of %d...\n", (files_started / files_concurrent) + 1,
                                 files_groups);
                    else if(modus == MODUS_WRITE)
                      mpi_printf("RESTART: Writing restart files group #%d out of %d...\n", (files_started / files_concurrent) + 1,
                                 files_groups);
                    else
                      mpi_printf("RESTART: Checking restart files group #%d out of %d...\n", (files_started / files_concurrent) + 1,
                                 files_groups);
                  }

                send_work_request(modus, files_started++);
              }
          }
      }
}

/*! \brief Schedule the reading/writing/checking of restart files to ensure
 *         only NumFilesWrittenInParallel are written in parallel.
 *
 *  \param[in] modus Read, write or check files.
 *
 *  \return void
 */
static void work_files(int modus)
{
  if(ThisTask == 0)
    if(!(seq = malloc(NTask * sizeof(struct seq_data))))
      terminate("can't allocate seq_data");

  struct seq_data seq_loc;
  seq_loc.thistask   = ThisTask;
  seq_loc.rankinnode = RankInThisNode;
  seq_loc.thisnode   = ThisNode;

  MPI_Gather(&seq_loc, sizeof(struct seq_data), MPI_BYTE, seq, sizeof(struct seq_data), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      qsort(seq, NTask, sizeof(struct seq_data), compare_seq_data);
      if(seq[0].thistask != 0)
        terminate("unexpected");

      files_started   = 0;
      files_completed = 0;

      if((files_started % files_concurrent) == 0)
        {
          if(modus == MODUS_READ)
            mpi_printf("RESTART: Loading restart files group #%d out of %d...\n", (files_started / files_concurrent) + 1,
                       files_groups);
          else if(modus == MODUS_WRITE)
            mpi_printf("RESTART: Writing restart files group #%d out of %d...\n", (files_started / files_concurrent) + 1,
                       files_groups);
          else
            mpi_printf("RESTART: Checking restart files group #%d out of %d...\n", (files_started / files_concurrent) + 1,
                       files_groups);
        }

      for(int i = 1; i < All.NumFilesWrittenInParallel; i++)
        {
          files_started++;
          send_work_request(modus, i);
        }

      files_started++;
      if(!((modus == MODUS_WRITE && write_success[ThisTask]) || (modus == MODUS_CHECK && write_success[get_file_to_check(ThisTask)])))
        {
          if(modus == MODUS_CHECK)
            {
              char buf[MAXLEN_PATH];
              int task = get_file_to_check(ThisTask);
              get_restart_filename(buf, task, modus);

              int success;
              check_restart_files(buf, &checks[task], &success);
              write_success[task] = success;
            }
          else
            {
              char buf[MAXLEN_PATH];
              get_restart_filename(buf, ThisTask, modus);
              write_or_read_this_processors_restart_file(modus, buf, &checks[0]);
            }
        }
      files_completed++;

      if(files_started < NTask)
        {
          if((files_started % files_concurrent) == 0)
            {
              if(modus == MODUS_READ)
                mpi_printf("RESTART: Loading restart files group #%d out of %d...\n", (files_started / files_concurrent) + 1,
                           files_groups);
              else if(modus == MODUS_WRITE)
                mpi_printf("RESTART: Writing restart files group #%d out of %d...\n", (files_started / files_concurrent) + 1,
                           files_groups);
              else
                mpi_printf("RESTART: Checking restart files group #%d out of %d...\n", (files_started / files_concurrent) + 1,
                           files_groups);
            }

          send_work_request(modus, files_started++);
        }

      while(files_completed < NTask)
        polling(modus);

      free(seq);
    }
  else
    {
      int type;
      MPI_Recv(&type, 1, MPI_INT, 0, TAG_N, MPI_COMM_WORLD, MPI_STATUS_IGNORE); /* wait until we are told to start */

      if(type == 0)
        {
          if(modus == MODUS_CHECK)
            {
              struct check ch;
              MPI_Recv(&ch, sizeof(struct check), MPI_BYTE, 0, TAG_N, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              char buf[MAXLEN_PATH];
              get_restart_filename(buf, get_file_to_check(ThisTask), modus);

              int success;
              check_restart_files(buf, &ch, &success);
              MPI_Ssend(&success, 1, MPI_INT, 0, TAG_KEY, MPI_COMM_WORLD);
            }
          else
            {
              char buf[MAXLEN_PATH];
              get_restart_filename(buf, ThisTask, modus);
              struct check ch;
              write_or_read_this_processors_restart_file(modus, buf, &ch);

              if(modus == MODUS_WRITE)
                {
                  MPI_Ssend(&ch, sizeof(struct check), MPI_BYTE, 0, TAG_KEY, MPI_COMM_WORLD);
                }
              else
                {
                  int dummy = 0;
                  MPI_Ssend(&dummy, 1, MPI_INT, 0, TAG_KEY, MPI_COMM_WORLD);
                }
            }
        }
      else
        {
          int dummy = 1;
          MPI_Ssend(&dummy, 1, MPI_INT, 0, TAG_KEY, MPI_COMM_WORLD);
        }
    }
}

/*! \brief This function reads or writes the restart files.
 *
 *  Each processor writes its own restart file, with the
 *  I/O being done in parallel. To avoid congestion of the disks
 *  you can tell the program to restrict the number of files
 *  that are simultaneously written to NumFilesWrittenInParallel.
 *
 *  \param[in] modus if modus==MODUS_READ  the restart()-routine reads,
 *  if modus==MODUS_WRITE it writes a restart file.
 *
 *  \return void
 */
void restart(int modus)
{
  CPU_Step[CPU_MISC] += measure_time();
  double t0  = second();
  byte_count = 0;

  PageSize = getpagesize();
  mpi_printf("RESTART: PageSize = %d\n", PageSize);

  if(modus == MODUS_READ)
    mpi_printf("RESTART: Loading restart files...\n");

  if(ThisTask == 0 && modus == MODUS_WRITE)
    {
      printf("RESTART: Writing restart files");
      create_restartfiles_dir();
    }
  MPI_Barrier(MPI_COMM_WORLD);

  if(NTask < All.NumFilesWrittenInParallel)
    {
      warn("Number of processors should be a smaller or equal than `NumFilesWrittenInParallel'. We're adjusting the latter.\n");
      All.NumFilesWrittenInParallel = NTask;
    }

  if(All.NumFilesWrittenInParallel < 1)
    All.NumFilesWrittenInParallel = 1;

  files_concurrent = All.NumFilesWrittenInParallel;
  files_groups     = NTask / All.NumFilesWrittenInParallel;
  if(NTask % All.NumFilesWrittenInParallel)
    files_groups++;

#ifndef MULTIPLE_RESTARTS
  if(modus == MODUS_WRITE) /* write */
    backup_restartfiles(ThisTask);
#endif /* #ifndef MULTIPLE_RESTARTS */

  if(modus == MODUS_WRITE)
    if(ThisTask == 0)
      {
        if(!(checks = malloc(NTask * sizeof(struct check))))
          terminate("can't allocate checks");
        if(!(write_success = malloc(NTask)))
          terminate("can't allocate write_success");

        for(int i = 0; i < NTask; i++)
          {
            checks[i].byte_count = 0;
            write_success[i]     = 0;
          }
      }

  work_files(modus);

  MPI_Barrier(MPI_COMM_WORLD);

  if(modus == MODUS_WRITE)
    {
      int iter    = 0;
      int success = 0;
      while(!success)
        {
          work_files(MODUS_CHECK);

          if(ThisTask == 0)
            {
              int count = 0;
              for(int i = 0; i < NTask; i++)
                {
                  if(!write_success[i])
                    count++;
                }

              if(count == 0)
                {
                  printf("All restart files written successfully.\n");
                  success = 1;
                }
              else
                {
                  printf("Need to repeat writing for %d restartfiles.\n", count);
                }
            }

          MPI_Bcast(&success, 1, MPI_INT, 0, MPI_COMM_WORLD);

          if(success)
            break;

          iter++;
          if(iter > 4)
            terminate("Too many iterations, fix your file system.");

          work_files(MODUS_WRITE);
        };

      free(checks);
    }

  /* check whether the restarts are all at the same time */
  if(modus == MODUS_READ) /* read */
    {
      struct global_data_all_processes all_task0;

      if(ThisTask == 0)
        all_task0 = All;

      MPI_Bcast(&all_task0, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);

      if(all_task0.Time != All.Time)
        terminate("The restart file on task=%d is not consistent with the one on task=0\n", ThisTask);
    }

  long long byte_count_all;
  sumup_longs(1, &byte_count, &byte_count_all);

  double t1 = second();

  mpi_printf("RESTART: load/save took %g sec, corresponds to I/O rate of %g MB/sec\n", timediff(t0, t1),
             byte_count_all / (1024.0 * 1024.0) / timediff(t0, t1));

  CPU_Step[CPU_RESTART] += measure_time();
  mpi_printf("RESTART: done.\n");
}

/*! \brief Reads or writes restart file.
 *
 *  Try write until successful.
 *
 *  \param[in] modus Flag for write or read.
 *  \param[in] buf File name.
 *  \param[in] ch Check metadata.
 *
 *  \return void
 */
static void write_or_read_this_processors_restart_file(int modus, char *buf, struct check *ch)
{
  if(modus == MODUS_READ)
    {
      execute_write_or_read(MODUS_READ, buf, ch);
    }
  else
    {
      int failed = 0;

      do
        {
          execute_write_or_read(MODUS_WRITE, buf, ch);
        }
      while(failed > 0);
    }
}

/*! \brief Reads or writes a restart file.
 *
 *  A single attempt which either is successful or fails.
 *
 *  \param[in] modus Flag for write or read.
 *  \param[in] buf File name.
 *  \param[in] ch Check metadata.
 *
 *  \return 0: success, 1: failed.
 */
static int execute_write_or_read(int modus, char *buf, struct check *ch)
{
  if(modus == MODUS_WRITE)
    ch->byte_count = byte_count;

  int failed_flag = 0;

#ifdef TOLERATE_WRITE_ERROR
  for(int try_io = 0; try_io < 2; try_io++)
    {
      WriteErrorFlag = 0;
#endif /* #ifdef TOLERATE_WRITE_ERROR */
      if(modus == MODUS_READ || modus == MODUS_READCHECK)
        {
#ifdef USE_DIRECT_IO_FOR_RESTARTS
          struct stat st;
          if(stat(buf, &st) == 0)
            {
              size_t size = st.st_size;
              if(size % PageSize > 0)
                {
                  FILE *fd = fopen(buf, "a");
                  if(fd)
                    {
                      size_t n = PageSize - (size % PageSize);
                      char *p  = calloc(n, 1);
                      if(p == NULL)
                        terminate("p == NULL");
                      printf("RESTART: Topping of restart file '%s' by %lld bytes\n", buf, (long long)n);
                      fwrite(p, n, 1, fd);
                      fclose(fd);
                      free(p);
                    }
                  else
                    terminate("can't increase length of restart file '%s'", buf);
                }
            }
          else
            terminate("Restart file '%s' not found.\n", buf);
#endif /* #ifdef USE_DIRECT_IO_FOR_RESTARTS */
          int oflag = O_RDONLY;
#ifdef USE_DIRECT_IO_FOR_RESTARTS
          oflag |= O_DIRECT;
#endif /* #ifdef USE_DIRECT_IO_FOR_RESTARTS */
          if((fdint = open(buf, oflag)) < 0)
            terminate("Restart file '%s' not found.\n", buf);

          allocate_iobuf();
        }
      else
        {
#ifdef TOLERATE_WRITE_ERROR
          int try_open = 0;

          while(try_open < IO_TRIALS)
            {
              int oflag = O_WRONLY | O_CREAT | O_TRUNC;
#ifdef USE_DIRECT_IO_FOR_RESTARTS
              oflag |= O_DIRECT;
#endif /* #ifdef USE_DIRECT_IO_FOR_RESTARTS */
              if((fdint = open(buf, oflag, S_IRUSR | S_IWUSR | S_IRGRP)) < 0)
                {
                  printf("Restart file '%s' cannot be opened. Trying again...\n", buf);
                  myflush(stdout);

                  try_open++;

                  sleep(IO_SLEEP_TIME);
                }
              else
                break;
            }

          if(try_open == IO_TRIALS)
            terminate("Opening of restart file failed too often!");
#else /* #ifdef TOLERATE_WRITE_ERROR */
      int oflag = O_WRONLY | O_CREAT | O_TRUNC;
#ifdef USE_DIRECT_IO_FOR_RESTARTS
      oflag |= O_DIRECT;
#endif /* #ifdef USE_DIRECT_IO_FOR_RESTARTS */
      if((fdint = open(buf, oflag, S_IRUSR | S_IWUSR | S_IRGRP)) < 0)
        terminate("Restart file '%s' cannot be opened.\n", buf);
#endif /* #ifdef TOLERATE_WRITE_ERROR #else */
          allocate_iobuf();
        }

      MD5Init(&mysum);

      contents_restart_file(modus);

      MD5Final(&mysum);

      unsigned char has_hash[16];
      static unsigned char should_hash[16];

      for(int k = 0; k < 16; k++)
        has_hash[k] = mysum.digest[k];

      if(modus == MODUS_READ)
        {
          /* read */
          unsigned char written_hash[16];
          byten_nohash(written_hash, 16, modus);
          if(memcmp(has_hash, written_hash, 16) != 0)
            {
              char str_has[48], str_written[48];
              for(int i = 0; i < 16; i++)
                {
                  sprintf(str_has + 2 * i, "%02X", has_hash[i]);
                  sprintf(str_written + 2 * i, "%02X", written_hash[i]);
                }

              str_has[32] = str_written[32] = 0;

              terminate("RESTART: file '%s' does not match expected MD5 hash of '%s', found '%s' instead.", buf, str_has, str_written);
            }
        }
      else if(modus == MODUS_READCHECK)
        {
          if(memcmp(should_hash, has_hash, 16) != 0)
            {
              char str_should[48], str_has[48];
              for(int i = 0; i < 16; i++)
                {
                  sprintf(str_should + 2 * i, "%02X", should_hash[i]);
                  sprintf(str_has + 2 * i, "%02X", has_hash[i]);
                }

              str_should[32] = str_has[32] = 0;

              failed_flag = 1;

              terminate(
                  "RESTART-READCHECK: file '%s' does not match expected MD5 hash of '%s' after read-back check, has '%s' instead.",
                  buf, str_should, str_has);
            }
#ifdef VERBOSE
          else
            {
              char str_should[48], str_has[48];
              for(int i = 0; i < 16; i++)
                {
                  sprintf(str_should + 2 * i, "%02X", should_hash[i]);
                  sprintf(str_has + 2 * i, "%02X", has_hash[i]);
                }

              str_should[32] = str_has[32] = 0;

              printf("RESTART-READCHECK: Task %d: file '%s' does match expected MD5 hash of '%s' after read-back check, has '%s'.\n",
                     ThisTask, buf, str_should, str_has);
            }
#endif /* #ifdef VERBOSE */
        }
      else if(modus == MODUS_WRITE)
        {
          ch->byte_count = byte_count - ch->byte_count;
          for(int k = 0; k < 16; k++)
            ch->hash[k] = has_hash[k];

          /* write */
          byten_nohash(has_hash, 16, modus);

          for(int k = 0; k < 16; k++)
            should_hash[k] = has_hash[k];
        }
      else
        terminate("This should not happen - wrong modus!");

      deallocate_iobuf(modus);

      close(fdint);

#ifdef TOLERATE_WRITE_ERROR
      if(WriteErrorFlag == 0)
        break;

      if(try_io == 0)
        {
          char alternative_fname[MAXLEN_PATH];
          sprintf(alternative_fname, "%s/restartfiles/%s.%d", AlternativeOutputDir, "restart", ThisTask);

          printf("TOLERATE_WRITE_ERROR: Try to write to alternative file: Task=%d try_io=%d alternative-filename='%s'\n", ThisTask,
                 try_io, alternative_fname);
          myflush(stdout);
          strncpy(buf, alternative_fname, MAXLEN_PATH); /* try on a different output directory */
        }
      else
        {
          terminate("TOLERATE_WRITE_ERROR: Second try with alternative file failed too.\n");
        }
    }
#endif /* #ifdef TOLERATE_WRITE_ERROR */

  return failed_flag;
}

/*! \brief Defines contents of restart file.
 *
 *  \param[in] modus Read or write (0: write; >0 read).
 *
 *  \return void
 */
static void contents_restart_file(int modus)
{
  /* common data  */
  byten(&All, sizeof(struct global_data_all_processes), modus);

  /* individual allocation factors for meshes */
  byten(&Mesh.Indi, sizeof(struct individual_alloc_data), modus);
  byten(&DeRefMesh.Indi, sizeof(struct individual_alloc_data), modus);

  polling(modus);

  if(modus == MODUS_READ) /* read */
    allocate_memory();

  int ntask = NTask;
  in(&ntask, modus);

  if(modus == MODUS_READ)
    if(ntask != NTask)
      terminate("The restart files were written for ntask=%d while you're using now %d MPI ranks\n", ntask, NTask);

  in(&NumPart, modus);

  /* Particle data  */
  byten(&P[0], NumPart * sizeof(struct particle_data), modus);

  polling(modus);

  in(&NumGas, modus);

  if(NumGas > 0)
    {
      /* Sph-Particle data  */
      byten(&SphP[0], NumGas * sizeof(struct sph_particle_data), modus);
    }

  polling(modus);

  in(&Nvc, modus);
  in(&MaxNvc, modus);
  in(&FirstUnusedConnection, modus);

  if(modus == MODUS_READ) /* read */
    DC = mymalloc_movable(&DC, "DC", MaxNvc * sizeof(connection));

  byten(DC, MaxNvc * sizeof(connection), modus);

  polling(modus);

  /* write state of random number generators */
  byten(gsl_rng_state(random_generator), gsl_rng_size(random_generator), modus);
  byten(gsl_rng_state(random_generator_aux), gsl_rng_size(random_generator_aux), modus);

  /* now store variables for time integration bookkeeping */
  byten(TimeBinSynchronized, TIMEBINS * sizeof(int), modus);

  in(&TimeBinsHydro.NActiveParticles, modus);
  in(&TimeBinsGravity.NActiveParticles, modus);
  byten(&TimeBinsHydro.GlobalNActiveParticles, sizeof(long long), modus);
  byten(&TimeBinsGravity.GlobalNActiveParticles, sizeof(long long), modus);
  byten(TimeBinsHydro.ActiveParticleList, TimeBinsHydro.NActiveParticles * sizeof(int), modus);
  byten(TimeBinsGravity.ActiveParticleList, TimeBinsGravity.NActiveParticles * sizeof(int), modus);
  byten(TimeBinsHydro.NextInTimeBin, NumGas * sizeof(int), modus);
  byten(TimeBinsGravity.NextInTimeBin, NumPart * sizeof(int), modus);
  byten(TimeBinsHydro.PrevInTimeBin, NumGas * sizeof(int), modus);
  byten(TimeBinsGravity.PrevInTimeBin, NumPart * sizeof(int), modus);
  byten(TimeBinsHydro.TimeBinCount, TIMEBINS * sizeof(int), modus);
  byten(TimeBinsGravity.TimeBinCount, TIMEBINS * sizeof(int), modus);
  byten(TimeBinsHydro.FirstInTimeBin, TIMEBINS * sizeof(int), modus);
  byten(TimeBinsGravity.FirstInTimeBin, TIMEBINS * sizeof(int), modus);
  byten(TimeBinsHydro.LastInTimeBin, TIMEBINS * sizeof(int), modus);
  byten(TimeBinsGravity.LastInTimeBin, TIMEBINS * sizeof(int), modus);

#ifdef USE_SFR
  byten(TimeBinSfr, TIMEBINS * sizeof(double), modus);
#endif

  polling(modus);

  /* now store custom data for optional Config settings */
#ifdef USE_SFR
  in(&Stars_converted, modus);
#endif

  polling(modus);

  /* now store relevant data for tree */

  in(&NTopleaves, modus);
  in(&NTopnodes, modus);

  in(&Ngb_MaxPart, modus);
  in(&Ngb_MaxNodes, modus);
  in(&Ngb_NumNodes, modus);
  in(&Ngb_MarkerValue, modus);
  in(&Ngb_FirstNonTopLevelNode, modus);

  polling(modus);

  if(modus == MODUS_READ) /* read */
    {
      domain_allocate();
      ngb_treeallocate();
    }

  if(All.TotNumGas > 0)
    {
#ifdef TREE_BASED_TIMESTEPS
      byten(ExtNgb_Nodes + Ngb_MaxPart, Ngb_NumNodes * sizeof(struct ExtNgbNODE), modus);
#endif /* #ifdef TREE_BASED_TIMESTEPS */
      byten(Ngb_Nodes + Ngb_MaxPart, Ngb_NumNodes * sizeof(struct NgbNODE), modus);
      byten(Ngb_DomainNodeIndex, NTopleaves * sizeof(int), modus);
      byten(Ngb_Nextnode, (Ngb_MaxPart + NTopleaves) * sizeof(int), modus);
      byten(Ngb_Father, Ngb_MaxPart * sizeof(int), modus);
      byten(Ngb_Marker, (Ngb_MaxPart + NTopleaves) * sizeof(int), modus);
    }

  polling(modus);

  byten(TopNodes, NTopnodes * sizeof(struct topnode_data), modus);
  byten(DomainTask, NTopleaves * sizeof(int), modus);
  byten(DomainCorner, 3 * sizeof(double), modus);
  byten(DomainCenter, 3 * sizeof(double), modus);
  byten(&DomainLen, sizeof(double), modus);
  byten(&DomainFac, sizeof(double), modus);
  byten(&DomainInverseLen, sizeof(double), modus);
  byten(&DomainBigFac, sizeof(double), modus);
}

/*! \brief Adjusts the timeline if the TimeMax variable is
 *  increased between a restart.
 *
 *  The approach taken here is to reduce the resolution of the
 *  integer timeline by factors of 2 until the new final time
 *  can be reached within TIMEBASE.
 *
 *  \param[in] TimeMax_old old final time.
 *  \param[in] TimeMax_new new final time (must be larger than old one).
 *
 *  \return void
 */
void readjust_timebase(double TimeMax_old, double TimeMax_new)
{
  int i;
  long long ti_end;

  if(sizeof(long long) != 8)
    terminate("\nType 'long long' is not 64 bit on this platform\n\n");

  mpi_printf("\nRESTART: All.TimeMax has been changed in the parameterfile\nNeed to adjust integer timeline\n\n\n");

  if(TimeMax_new < TimeMax_old)
    terminate("\nIt is not allowed to reduce All.TimeMax\n\n");

  if(All.ComovingIntegrationOn)
    ti_end = (long long)(log(TimeMax_new / All.TimeBegin) / All.Timebase_interval);
  else
    ti_end = (long long)((TimeMax_new - All.TimeBegin) / All.Timebase_interval);

  while(ti_end > TIMEBASE)
    {
      All.Timebase_interval *= 2.0;

      ti_end /= 2;
      All.Ti_Current /= 2;
      All.Previous_Ti_Current /= 2;

#ifdef PMGRID
      All.PM_Ti_begstep /= 2;
      All.PM_Ti_endstep /= 2;
#endif /* #ifdef PMGRID */

      for(i = 0; i < NumPart; i++)
        {
          P[i].Ti_Current /= 2;

          if(P[i].TimeBinGrav > 0)
            {
              P[i].TimeBinGrav--;
              if(P[i].TimeBinGrav <= 0)
                {
                  char buf[1000];
                  sprintf(buf, "Error in readjust_timebase(). Minimum Timebin for particle %d reached.\n", i);
                  terminate(buf);
                }
            }

          if(P[i].Type == 0)
            if(P[i].TimeBinHydro > 0)
              {
                P[i].TimeBinHydro--;
                if(P[i].TimeBinHydro <= 0)
                  {
                    char buf[1000];
                    sprintf(buf, "Error in readjust_timebase(). Minimum Timebin for particle %d reached.\n", i);
                    terminate(buf);
                  }
              }
        }
    }

  All.TimeMax = TimeMax_new;
}

/*! \brief Reads/writes one integer to a restart file.
 *
 *  \param[in, out] x pointer to the integer.
 *  \param[in] modus if modus>0  the restart()-routine reads,
 *             if modus==0 it writes a restart file.
 *
 *  \return void
 */
void in(int *x, int modus) { byten(x, sizeof(int), modus); }

/*! \brief Reads/writes n bytes to restart file buffer.
 *
 * \param[in, out] x Pointer to the data.
 * \param[in] n Number of bytes.
 * \param[in] modus If modus>0  the restart()-routine reads,
 *            if modus==0 it writes a restart file.
 *
 *  \return void
 */
void byten(void *x, size_t n, int modus) { byten_hash(x, n, modus, 1); }

/*! \brief Wrapper for byten; called with hash=0.
 *
 *  I.e. writes something without including it in calculating the md5sum. This
 *  should only be done for the md5sum itself, but not for actual data.
 *
 *  \param[in, out] x Pointer to the data.
 *  \param[in] n Number of bytes.
 *  \param[in] modus If modus>0  the restart()-routine reads,
 *             if modus==0 it writes a restart file.
 *
 *  \return void
 */
void byten_nohash(void *x, size_t n, int modus) { byten_hash(x, n, modus, 0); }

/*! \brief Reads/writes n bytes to restart file buffer.
 *
 * \param[in, out] x pointer to the data
 * \param[in] n number of bytes
 * \param[in] modus if modus>0  the restart()-routine reads,
 *            if modus==0 it writes a restart file.
 * \param[in] hash If nonzero, this part is considered in md5sum.
 *
 *
 *  \return void
 */
void byten_hash(void *x, size_t n, int modus, int hash)
{
  byte_count += n;

  if(n > 0)
    {
      size_t nin = n;

      if(modus == MODUS_READ || modus == MODUS_READCHECK || modus == MODUS_CHECK) /* read */
        {
          if(modus == MODUS_READCHECK || modus == MODUS_CHECK)
            x = mymalloc("x", n);

          unsigned char *ptr = x;

          while(n > 0)
            {
              if(iop != fillp)
                {
                  size_t nn = n;
                  if(nn > (fillp - iop))
                    nn = fillp - iop;

                  memcpy(ptr, iobuf_aligned + iop, nn);

                  n -= nn;
                  ptr += nn;
                  iop += nn;
                }
              else
                {
                  if(iop == MAX_BLOCK_SIZE)
                    {
                      iop   = 0;
                      fillp = 0;
                    }

                  size_t nn = n;
                  if(nn % PageSize > 0)
                    nn = (nn / PageSize + 1) * PageSize;

                  if(nn > MAX_BLOCK_SIZE - fillp)
                    nn = MAX_BLOCK_SIZE - fillp;

                  if(read(fdint, iobuf_aligned + fillp, nn) != nn)
                    terminate("read error");

                  fillp += nn;
                }
            }

          if(hash) /* to prevent call if we write/load the checksum itself */
            MD5UpdateLong(&mysum, x, nin);

          if(modus == MODUS_READCHECK || modus == MODUS_CHECK)
            myfree(x);
        }
      else /* write */
        {
          unsigned char *ptr = x;

          while(n > 0)
            {
              if(iop < MAX_BLOCK_SIZE)
                {
                  size_t nn = n;
                  if(nn > MAX_BLOCK_SIZE - iop)
                    nn = MAX_BLOCK_SIZE - iop;
                  memcpy(iobuf_aligned + iop, ptr, nn);

                  n -= nn;
                  ptr += nn;
                  iop += nn;
                }
              else
                {
                  size_t nn = MAX_BLOCK_SIZE;
                  if(write(fdint, iobuf_aligned, nn) != nn)
                    terminate("write error");

                  iop = 0;
                }
            }

          if(hash) /* to prevent call if we write/load the checksum itself */
            MD5UpdateLong(&mysum, x, nin);
        }
    }
}

/*! \brief Allocates the IO buffer for reading/writing the restart-file buffer.
 *
 *  \return void
 */
void allocate_iobuf(void)
{
  if((MAX_BLOCK_SIZE % PageSize) > 0)
    terminate("MAX_BLOCK_SIZE must be a multiple of PageSize");

  if(!(io_buf = malloc(MAX_BLOCK_SIZE + PageSize)))
    terminate("cannot allocated IO buffer");

  iobuf_aligned = (char *)(((((size_t)io_buf) + (PageSize - 1)) / PageSize) * PageSize);

  fillp = 0;
  iop   = 0;
}

/*! \brief Frees the IO buffer for reading/writing the restart-files.
 *
 *  Writes buffer before freeing it if in MODUS_WRITE.
 *
 *  \param[in] modus Read or write.
 *
 *  \return void
 */
void deallocate_iobuf(int modus)
{
  if(modus == MODUS_WRITE) /* write */
    {
      if(iop > 0)
        {
          if(iop % PageSize > 0)
            iop = ((iop / PageSize) + 1) * PageSize;

          if(write(fdint, iobuf_aligned, iop) != iop)
            terminate("write error");
        }
    }

  free(io_buf);
}
