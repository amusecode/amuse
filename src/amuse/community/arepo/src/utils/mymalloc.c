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
 * \file        src/utils/mymalloc.c
 * \date        05/2018
 * \brief       Manager for dynamic memory allocation.
 * \details     This module handles the dynamic memory allocation for AREPO.
 *              To avoid memory allocation/dellocation overhead a big chunk of
 *              memory (which will be the maximum amount of dinamically
 *              allocatable memory) is allocated upon initialization. This
 *              chunk is then filled by the memory blocks as in a stack
 *              structure. The blocks are automatically aligned to a 64 bit
 *              boundary. Memory blocks come in two flavours: movable and
 *              non-movable. In non-movable blocks the starting address is
 *              fixed once the block is allocated and cannot be changed.
 *              Due to the stack structure of the dynamic memory, this implies
 *              that the last (non-movable) block allocated must be the first
 *              block to be deallocated. If this condition is not met, an abort
 *              condition is triggered. If more flexibility is needed, movable
 *              memory blocks can be used. In this case, the starting address
 *              of the block is again fixed upon allocation but the block can
 *              be shifted (therefore its initial address changes) according
 *              to needs. For a movable block to be successfully shifted it is
 *              required that all the subsequent allocated blocks are movable.
 *              Again, an abort condition is triggered if this condition is
 *              not met. Movable blocks can be deallocated in any order
 *              provided that the condition just described holds. The gap
 *              resulting form the deallocation of a block that is not in
 *              the last position will be automatically filled by shifting all
 *              the blocks coming after the deallocated block.
 *
 *              contains functions:
 *                static void *hmalloc(size_t size) (HUGEPAGES)
 *                static void *hmalloc(size_t size)
 *                void mymalloc_init(void)
 *                void report_memory_usage(int rank, char *tabbuf)
 *                void report_detailed_memory_usage_of_largest_task(void)
 *                void dump_memory_table(void)
 *                int dump_memory_table_buffer(char *p)
 *                void *mymalloc_fullinfo
 *                void *mymalloc_movable_fullinfo
 *                size_t roundup_to_multiple_of_cacheline_size(size_t n)
 *                void myfree_fullinfo(void *p, const char *func, const char
 *                  *file, int line)
 *                void *myfree_query_last_block(void)
 *                void myfree_movable_fullinfo(void *p, const char *func,
 *                  const char *file, int line)
 *                void *myrealloc_fullinfo(void *p, size_t n, const char
 *                  *func, const char *file, int line)
 *                void *myrealloc_movable_fullinfo(void *p, size_t n,
 *                  const char *func, const char *file, int line)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 07.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#define CACHELINESIZE 64

#define MAXBLOCKS 5000
#define MAXCHARS 40

static size_t AllocatedBytesGeneric;

static size_t HighMarkBytes;
static size_t HighMarkBytesWithoutGeneric;

static double OldGlobHighMarkMB;
static double OldGlobHighMarkMBWithoutGeneric;

static size_t TotBytes; /*!< The total dimension (in bytes) of dynamic memory available to the current task. */
static void *Base;      /*!< Base pointer (initial memory address) of the stack. */

static unsigned long Nblocks; /*!< The current number of allocated memory blocks. */

static void **Table;         /*!< Table containing the initial addresses of the allocated memory blocks. */
static size_t *BlockSize;    /*!< Array containing the size (in bytes) of all the allocated memory blocks. */
static char *MovableFlag;    /*!< Identifies whether a block is movable. */
static char *GenericFlag;    /*!< Identifies whether a block has been identified in the generic allocation routines. */
static void ***BasePointers; /*!< Base pointers containing the initial addresses of movable memory blocks */
static char *VarName;        /*!< The name of the variable with which the block has been allocated. */
static char *FunctionName;   /*!< The function name that has allocated the memory block. */
static char *ParentFileName; /*!< The location from which the generich routines were called */
static char *FileName;       /*!< The file name where the function that has allocated the block is called. */
static int *LineNumber;      /*!< The line number in FileName where the function that allocated the block has been called. */
static char *HighMarkTabBuf; /*!< This is a buffer that holds the log-file output corresponding to the largest memory use that has
                                occurred on this task */
static char *HighMarkTabBufWithoutGeneric; /*!< This is a buffer that holds the log-file output corresponding to the largest memory use
                                              that has occurred on this task */

#ifdef HUGEPAGES
#include <hugetlbfs.h>
/*! \brief Allocation function wrapper for hugepages usage.
 *
 *  \param[in] size Size of the allocated memory.
 *
 *  \return void pointer to address in memory.
 */
static void *hmalloc(size_t size)
{
  void *p = get_hugepage_region(size, GHR_STRICT);

  if(!p)
    {
      warn("Failed to get_hugepage_region of size %g\n", size / (1024.0 * 1024));

      p = malloc(size);

      if(!p)
        terminate("Failed to allocate memory of size %g\n", size / (1024.0 * 1024));
    }

  memset(p, 255, size);
  memset(p, 0, size);

  return p;
}
#else  /* #ifdef HUGEPAGES */

/*! \brief Allocation function wrapper without hugepages usage.
 *
 *  \param[in] size Size of the allocated memory.
 *
 *  \return void pointer to address in memory.
 */
static void *hmalloc(size_t size) { return malloc(size); }
#endif /* #ifdef HUGEPAGES #else */

/*! \brief Initializes memory manager.
 *
 *  This function initializes the memory manager. In particular, it sets
 *  the global variables of the module to their initial value and allocates
 *  the memory for the stack.
 *
 *  \return void
 */
void mymalloc_init(void)
{
  BlockSize                    = (size_t *)hmalloc(MAXBLOCKS * sizeof(size_t));
  Table                        = (void **)hmalloc(MAXBLOCKS * sizeof(void *));
  MovableFlag                  = (char *)hmalloc(MAXBLOCKS * sizeof(char));
  GenericFlag                  = (char *)hmalloc(MAXBLOCKS * sizeof(char));
  BasePointers                 = (void ***)hmalloc(MAXBLOCKS * sizeof(void **));
  VarName                      = (char *)hmalloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  FunctionName                 = (char *)hmalloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  ParentFileName               = (char *)hmalloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  FileName                     = (char *)hmalloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  LineNumber                   = (int *)hmalloc(MAXBLOCKS * sizeof(int));
  HighMarkTabBuf               = (char *)hmalloc((100 + 4 * MAXCHARS) * (MAXBLOCKS + 10));
  HighMarkTabBufWithoutGeneric = (char *)hmalloc((100 + 4 * MAXCHARS) * (MAXBLOCKS + 10));

  memset(VarName, 0, MAXBLOCKS * MAXCHARS);
  memset(FunctionName, 0, MAXBLOCKS * MAXCHARS);
  memset(ParentFileName, 0, MAXBLOCKS * MAXCHARS);
  memset(FileName, 0, MAXBLOCKS * MAXCHARS);

  size_t n = All.MaxMemSize * ((size_t)1024 * 1024);

  n = roundup_to_multiple_of_cacheline_size(n);

  if(!(Base = hmalloc(n)))
    terminate("Failed to allocate memory for `Base' (%d Mbytes).\n", All.MaxMemSize);

  TotBytes = FreeBytes = n;

  AllocatedBytes                  = 0;
  Nblocks                         = 0;
  HighMarkBytes                   = 0;
  HighMarkBytesWithoutGeneric     = 0;
  OldGlobHighMarkMB               = 0;
  OldGlobHighMarkMBWithoutGeneric = 0;
}

/*! \brief Writes memory usage in FdMemory.
 *
 *  \param[in] rank Number of tasks involved.
 *  \param[in] tabbuf Header message written in FdMemory.
 *
 *  \return void
 */
void report_memory_usage(int rank, char *tabbuf)
{
  if(ThisTask == rank)
    {
      char *buf = mymalloc("buf", (100 + 4 * MAXCHARS) * (Nblocks + 10));
      int cc    = 0;
      cc += sprintf(buf + cc, "\nMEMORY:  Largest Allocation = %g Mbyte  |  Largest Allocation Without Generic = %g Mbyte\n\n",
                    OldGlobHighMarkMB, OldGlobHighMarkMBWithoutGeneric);

      cc += sprintf(buf + cc, "%s", tabbuf);
      if(ThisTask == 0)
        {
          if(RestartFlag <= 2)
            {
              fprintf(FdMemory, "%s", buf);
              fflush(FdMemory);
            }
        }
      else
        {
          MPI_Send(&cc, 1, MPI_INT, 0, TAG_N, MPI_COMM_WORLD);
          MPI_Send(buf, cc + 1, MPI_BYTE, 0, TAG_PDATA, MPI_COMM_WORLD);
        }
      myfree(buf);
    }

  if(ThisTask == 0 && rank > 0)
    {
      int cc;
      MPI_Recv(&cc, 1, MPI_INT, rank, TAG_N, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      char *buf = mymalloc("buf", cc + 1);
      MPI_Recv(buf, cc + 1, MPI_BYTE, rank, TAG_PDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if(RestartFlag <= 2)
        {
          fprintf(FdMemory, "%s", buf);
          fflush(FdMemory);
        }
      myfree(buf);
    }
}

/*! \brief Writes memory useage of largest task in FdMemory.
 *
 *  \return void
 */
void report_detailed_memory_usage_of_largest_task(void)
{
  int flag = 0;

  struct
  {
    double mem;
    int rank;
  } local, global;

  local.mem  = HighMarkBytes / (1024.0 * 1024.0);
  local.rank = ThisTask;

  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  if(global.mem >= 1.05 * OldGlobHighMarkMB)
    {
      OldGlobHighMarkMB = global.mem;
      flag |= 1;
    }

  local.mem  = HighMarkBytesWithoutGeneric / (1024.0 * 1024.0);
  local.rank = ThisTask;

  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  if(global.mem >= 1.05 * OldGlobHighMarkMBWithoutGeneric)
    {
      OldGlobHighMarkMBWithoutGeneric = global.mem;
      flag |= 2;
    }

  if(flag & 2)
    report_memory_usage(global.rank, HighMarkTabBufWithoutGeneric);

  if(flag & 1)
    report_memory_usage(global.rank, HighMarkTabBuf);
}

/*! \brief Dumps the buffer where the memory information is stored to the
 *         standard output.
 *
 *  \return void
 */
void dump_memory_table(void)
{
  char *buf = malloc(200 * (Nblocks + 10));
  dump_memory_table_buffer(buf);
  printf("%s", buf);
  free(buf);
}

/*! \brief Fills the output buffer with the memory log.
 *
 *  \param[out] p Output buffer.
 *
 *  \return The number of charcter written to p.
 */
int dump_memory_table_buffer(char *p)
{
  int cc              = 0;
  size_t totBlocksize = 0;

  cc +=
      sprintf(p + cc, "-------------------------- Allocated Memory Blocks---- ( Step %8d )------------------\n", All.NumCurrentTiStep);
  cc += sprintf(p + cc, "Task    Nr F                  Variable      MBytes   Cumulative  Function|File|Linenumber\n");
  cc += sprintf(p + cc, "------------------------------------------------------------------------------------------\n");
  for(int i = 0; i < Nblocks; i++)
    {
      totBlocksize += BlockSize[i];

      cc += sprintf(p + cc, "%4d %5d %d %40s  %10.4f   %10.4f  %s%s()|%s|%d\n", ThisTask, i, MovableFlag[i], VarName + i * MAXCHARS,
                    BlockSize[i] / (1024.0 * 1024.0), totBlocksize / (1024.0 * 1024.0), ParentFileName + i * MAXCHARS,
                    FunctionName + i * MAXCHARS, FileName + i * MAXCHARS, LineNumber[i]);
    }
  cc += sprintf(p + cc, "------------------------------------------------------------------------------------------\n");

  return cc;
}

/*! \brief Allocates a non-movable memory block and store the relative
 *         information.
 *
 *  \param[in] varname Name of the variable to be stored in the allocated
 *             block.
 *  \param[in] n Size of the memory block in bytes.
 *  \param[in] func Name of function that has called the allocation routine
 *             (usually given by the __FUNCTION__ macro).
 *  \param[in] file File where the function that has called the allocation
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of file where the allocation routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return A pointer to the beginning of the allocated memory block.
 */
void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int line, int clear_flag, char *callorigin)
{
  if((n % CACHELINESIZE) > 0)
    n = (n / CACHELINESIZE + 1) * CACHELINESIZE;

  if(n < CACHELINESIZE)
    n = CACHELINESIZE;

  if(Nblocks >= MAXBLOCKS)
    terminate("Task=%d: No blocks left in mymalloc_fullinfo() at %s()/%s/line %d. MAXBLOCKS=%d\n", ThisTask, func, file, line,
              MAXBLOCKS);

  if(n > FreeBytes)
    {
      dump_memory_table();
      terminate(
          "\nTask=%d: Not enough memory in mymalloc_fullinfo() to allocate %g MB for variable '%s' at %s()/%s/line %d (FreeBytes=%g "
          "MB).\n",
          ThisTask, n / (1024.0 * 1024.0), varname, func, file, line, FreeBytes / (1024.0 * 1024.0));
    }
  Table[Nblocks] = Base + (TotBytes - FreeBytes);
  FreeBytes -= n;

  strncpy(VarName + Nblocks * MAXCHARS, varname, MAXCHARS - 1);
  if(callorigin)
    {
      strncpy(ParentFileName + Nblocks * MAXCHARS, callorigin, MAXCHARS - 1);
      GenericFlag[Nblocks] = 1;
      AllocatedBytesGeneric += n;
    }
  else
    {
      memset(ParentFileName + Nblocks * MAXCHARS, 0, MAXCHARS);
      GenericFlag[Nblocks] = 0;
    }
  strncpy(FunctionName + Nblocks * MAXCHARS, func, MAXCHARS - 1);
  strncpy(FileName + Nblocks * MAXCHARS, file, MAXCHARS - 1);
  LineNumber[Nblocks] = line;

  AllocatedBytes += n;
  BlockSize[Nblocks]   = n;
  MovableFlag[Nblocks] = 0;

  Nblocks += 1;

  if(AllocatedBytes - AllocatedBytesGeneric > HighMarkBytesWithoutGeneric)
    {
      HighMarkBytesWithoutGeneric = AllocatedBytes - AllocatedBytesGeneric;
      dump_memory_table_buffer(HighMarkTabBufWithoutGeneric);
    }

  if(AllocatedBytes > HighMarkBytes)
    {
      HighMarkBytes = AllocatedBytes;
      dump_memory_table_buffer(HighMarkTabBuf);
    }

  if(clear_flag)
    memset(Table[Nblocks - 1], 0, n);

  return Table[Nblocks - 1];
}

/*! \brief Allocates a movable memory block and store the relative information.
 *
 *  \param[in] ptr Pointer to the initial memory address of the block.
 *  \param[in] varname Name of the variable to be stored in the allocated block.
 *  \param[in] n Size of the memory block in bytes.
 *  \param[in] func Name of function that has called the allocation routine
 *             (usually given by the __FUNCTION__ macro).
 *  \param[in] file File where the function that has called the allocation
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of file where the allocation routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return A pointer to the beginning of the allocated memory block.
 */
void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file, int line,
                                char *callorigin)
{
  if((n % CACHELINESIZE) > 0)
    n = (n / CACHELINESIZE + 1) * CACHELINESIZE;

  if(n < CACHELINESIZE)
    n = CACHELINESIZE;

  if(Nblocks >= MAXBLOCKS)
    terminate("Task=%d: No blocks left in mymalloc_fullinfo() at %s()/%s/line %d. MAXBLOCKS=%d\n", ThisTask, func, file, line,
              MAXBLOCKS);

  if(n > FreeBytes)
    {
      dump_memory_table();
      terminate(
          "\nTask=%d: Not enough memory in mymalloc_fullinfo() to allocate %g MB for variable '%s' at %s()/%s/line %d (FreeBytes=%g "
          "MB).\n",
          ThisTask, n / (1024.0 * 1024.0), varname, func, file, line, FreeBytes / (1024.0 * 1024.0));
    }
  Table[Nblocks] = Base + (TotBytes - FreeBytes);
  FreeBytes -= n;

  strncpy(VarName + Nblocks * MAXCHARS, varname, MAXCHARS - 1);
  if(callorigin)
    {
      strncpy(ParentFileName + Nblocks * MAXCHARS, callorigin, MAXCHARS - 1);
      GenericFlag[Nblocks] = 1;
      AllocatedBytesGeneric += n;
    }
  else
    {
      memset(ParentFileName + Nblocks * MAXCHARS, 0, MAXCHARS);
      GenericFlag[Nblocks] = 0;
    }
  strncpy(FunctionName + Nblocks * MAXCHARS, func, MAXCHARS - 1);
  strncpy(FileName + Nblocks * MAXCHARS, file, MAXCHARS - 1);
  LineNumber[Nblocks] = line;

  AllocatedBytes += n;
  BlockSize[Nblocks]    = n;
  MovableFlag[Nblocks]  = 1;
  BasePointers[Nblocks] = ptr;

  Nblocks += 1;

  if(AllocatedBytes - AllocatedBytesGeneric > HighMarkBytesWithoutGeneric)
    {
      HighMarkBytesWithoutGeneric = AllocatedBytes - AllocatedBytesGeneric;
      dump_memory_table_buffer(HighMarkTabBufWithoutGeneric);
    }

  if(AllocatedBytes > HighMarkBytes)
    {
      HighMarkBytes = AllocatedBytes;
      dump_memory_table_buffer(HighMarkTabBuf);
    }

  return Table[Nblocks - 1];
}

/*! \brief Rounds up size to cachline size.
 *
 *  \param[in] n Size.
 *
 *  \return Rounded up size.
 */
size_t roundup_to_multiple_of_cacheline_size(size_t n)
{
  if((n % CACHELINESIZE) > 0)
    n = (n / CACHELINESIZE + 1) * CACHELINESIZE;

  return n;
}

/*! \brief Deallocates a non-movable memory block.
 *
 *  For this operation to be successful the block that has to be deallocated
 *  must be the last allocated one.
 *
 *  \param[in] p Pointer to the memory block to be deallocated.
 *  \param[in] func Name of function that has called the deallocation routine
 *             (usually given by the __FUNCTION__ macro).
 *  \param[in] file File where the function that has called the deallocation
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of file where the deallocation routine was
 *             called (usually given by the __LINE__ macro).
 */
void myfree_fullinfo(void *p, const char *func, const char *file, int line)
{
  if(Nblocks == 0)
    terminate("no allocated blocks that could be freed");

  if(p != Table[Nblocks - 1])
    {
      dump_memory_table();
      terminate("Task=%d: Wrong call of myfree() at %s()/%s/line %d: not the last allocated block!\n", ThisTask, func, file, line);
    }

  Nblocks -= 1;
  AllocatedBytes -= BlockSize[Nblocks];

  if(GenericFlag[Nblocks])
    AllocatedBytesGeneric -= BlockSize[Nblocks];

  FreeBytes += BlockSize[Nblocks];
}

/*! \brief Finds last allocated block.
 *
 *  \return void pointer to last allocated block.
 */
void *myfree_query_last_block(void)
{
  if(Nblocks == 0)
    terminate("no allocated blocks that could be returned");

  return Table[Nblocks - 1];
}

/*! \brief Deallocates a movable memory block.
 *
 *  For this operation to be successful all the blocks allocated after the
 *  block that has to be freed must be of movable type.
 *
 *  \param[in] p pointer to the memory block to be deallocated.
 *  \param[in] func name of function that has called the deallocation routine
 *             (usually given by the __FUNCTION__ macro).
 *  \param[in] file file where the function that has called the deallocation
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line line number of file where the deallocation routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return void
 */
void myfree_movable_fullinfo(void *p, const char *func, const char *file, int line)
{
  int i;

  if(Nblocks == 0)
    terminate("no allocated blocks that could be freed");

  /* first, let's find the block */
  int nr;

  for(nr = Nblocks - 1; nr >= 0; nr--)
    if(p == Table[nr])
      break;

  if(nr < 0)
    {
      dump_memory_table();
      terminate("Task=%d: Wrong call of myfree_movable() from %s()/%s/line %d - this block has not been allocated!\n", ThisTask, func,
                file, line);
    }

  if(nr < Nblocks - 1) /* the block is not the last allocated block */
    {
      /* check that all subsequent blocks are actually movable */
      for(i = nr + 1; i < Nblocks; i++)
        if(MovableFlag[i] == 0)
          {
            dump_memory_table();
            myflush(stdout);
            terminate(
                "Task=%d: Wrong call of myfree_movable() from %s()/%s/line %d - behind block=%d there are subsequent non-movable "
                "allocated blocks\n",
                ThisTask, func, file, line, nr);
          }
    }

  if(GenericFlag[nr])
    AllocatedBytesGeneric -= BlockSize[nr];

  AllocatedBytes -= BlockSize[nr];
  FreeBytes += BlockSize[nr];

  ptrdiff_t offset = -BlockSize[nr];
  size_t length    = 0;

  for(i = nr + 1; i < Nblocks; i++)
    length += BlockSize[i];

  if(nr < Nblocks - 1)
    memmove(Table[nr + 1] + offset, Table[nr + 1], length);

  for(i = nr + 1; i < Nblocks; i++)
    {
      Table[i] += offset;
      *BasePointers[i] = *BasePointers[i] + offset;
    }

  for(i = nr + 1; i < Nblocks; i++)
    {
      Table[i - 1]        = Table[i];
      BasePointers[i - 1] = BasePointers[i];
      BlockSize[i - 1]    = BlockSize[i];
      MovableFlag[i - 1]  = MovableFlag[i];
      GenericFlag[i - 1]  = GenericFlag[i];

      strncpy(VarName + (i - 1) * MAXCHARS, VarName + i * MAXCHARS, MAXCHARS - 1);
      strncpy(FunctionName + (i - 1) * MAXCHARS, FunctionName + i * MAXCHARS, MAXCHARS - 1);
      strncpy(ParentFileName + (i - 1) * MAXCHARS, ParentFileName + i * MAXCHARS, MAXCHARS - 1);
      strncpy(FileName + (i - 1) * MAXCHARS, FileName + i * MAXCHARS, MAXCHARS - 1);
      LineNumber[i - 1] = LineNumber[i];
    }

  Nblocks -= 1;
}

/*! \brief Reallocates an existing non-movable memory block.
 *
 *  For this operation to be successful this must be the last allocated block.
 *
 *  \param[in] p Pointer to the existing memory block to be reallocated.
 *  \param[in] n The new size of the memory block in bytes.
 *  \param[in] func Name of function that has called the reallocation routine
 *             (usually given by the __FUNCTION__ macro).
 *  \param[in] file File where the function that has called the reallocation
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of file where the reallocation routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return A pointer to the beginning of the newly allocated memory block.
 */
void *myrealloc_fullinfo(void *p, size_t n, const char *func, const char *file, int line)
{
  if((n % CACHELINESIZE) > 0)
    n = (n / CACHELINESIZE + 1) * CACHELINESIZE;

  if(n < CACHELINESIZE)
    n = CACHELINESIZE;

  if(Nblocks == 0)
    terminate("no allocated blocks that could be reallocated");

  if(p != Table[Nblocks - 1])
    {
      dump_memory_table();
      terminate("Task=%d: Wrong call of myrealloc() at %s()/%s/line %d - not the last allocated block!\n", ThisTask, func, file, line);
    }

  if(GenericFlag[Nblocks - 1])
    AllocatedBytesGeneric -= BlockSize[Nblocks - 1];

  AllocatedBytes -= BlockSize[Nblocks - 1];
  FreeBytes += BlockSize[Nblocks - 1];

  if(n > FreeBytes)
    {
      dump_memory_table();
      terminate("Task=%d: Not enough memory in myremalloc(n=%g MB) at %s()/%s/line %d. previous=%g FreeBytes=%g MB\n", ThisTask,
                n / (1024.0 * 1024.0), func, file, line, BlockSize[Nblocks - 1] / (1024.0 * 1024.0), FreeBytes / (1024.0 * 1024.0));
    }
  Table[Nblocks - 1] = Base + (TotBytes - FreeBytes);
  FreeBytes -= n;

  AllocatedBytes += n;
  BlockSize[Nblocks - 1] = n;

  if(AllocatedBytes > HighMarkBytes)
    {
      HighMarkBytes = AllocatedBytes;
      dump_memory_table_buffer(HighMarkTabBuf);
    }

  return Table[Nblocks - 1];
}

/*! \brief Reallocates an existing movable memory block.
 *
 *  For this operation to be successful all the blocks allocated after the
 *  block that has to be reallocated must be of movable type.
 *
 *  \param[in] p Pointer to the existing memory block to be reallocated.
 *  \param[in] n The new size of the memory block in bytes.
 *  \param[in] func Name of function that has called the reallocation routine
 *             (usually given by the __FUNCTION__ macro).
 *  \param[in] file File where the function that has called the reallocation
 *             routine resides (usually given by the __FILE__ macro).
 *  \param[in] line Line number of file where the reallocation routine was
 *             called (usually given by the __LINE__ macro).
 *
 *  \return A pointer to the beginning of the newly allocated memory block.
 */
void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line)
{
  int i;

  if((n % CACHELINESIZE) > 0)
    n = (n / CACHELINESIZE + 1) * CACHELINESIZE;

  if(n < CACHELINESIZE)
    n = CACHELINESIZE;

  if(Nblocks == 0)
    terminate("no allocated blocks that could be reallocated");

  /* first, let's find the block */
  int nr;

  for(nr = Nblocks - 1; nr >= 0; nr--)
    if(p == Table[nr])
      break;

  if(nr < 0)
    {
      dump_memory_table();
      terminate("Task=%d: Wrong call of myrealloc_movable() from %s()/%s/line %d - this block has not been allocated!\n", ThisTask,
                func, file, line);
    }

  if(nr < Nblocks - 1) /* the block is not the last allocated block */
    {
      /* check that all subsequent blocks are actually movable */
      for(i = nr + 1; i < Nblocks; i++)
        if(MovableFlag[i] == 0)
          {
            dump_memory_table();
            terminate(
                "Task=%d: Wrong call of myrealloc_movable() from %s()/%s/line %d - behind block=%d there are subsequent non-movable "
                "allocated blocks\n",
                ThisTask, func, file, line, nr);
          }
    }

  if(GenericFlag[nr])
    terminate("unexpected");

  AllocatedBytes -= BlockSize[nr];
  FreeBytes += BlockSize[nr];

  if(n > FreeBytes)
    {
      dump_memory_table();
      terminate("Task=%d: at %s()/%s/line %d: Not enough memory in myremalloc_movable(n=%g MB). previous=%g FreeBytes=%g MB\n",
                ThisTask, func, file, line, n / (1024.0 * 1024.0), BlockSize[nr] / (1024.0 * 1024.0), FreeBytes / (1024.0 * 1024.0));
    }

  ptrdiff_t offset = n - BlockSize[nr];
  size_t length    = 0;

  for(i = nr + 1; i < Nblocks; i++)
    length += BlockSize[i];

  if(nr < Nblocks - 1)
    memmove(Table[nr + 1] + offset, Table[nr + 1], length);

  for(i = nr + 1; i < Nblocks; i++)
    {
      Table[i] += offset;

      *BasePointers[i] = *BasePointers[i] + offset;
    }

  FreeBytes -= n;
  AllocatedBytes += n;
  BlockSize[nr] = n;

  if(AllocatedBytes > HighMarkBytes)
    {
      HighMarkBytes = AllocatedBytes;
      dump_memory_table_buffer(HighMarkTabBuf);
    }

  return Table[nr];
}
