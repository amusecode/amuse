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
 * \file        src/utils/generic_comm_helpers.h
 * \date        05/2018
 * \brief       Generic 'template' MPI communication structure used in many
 *              parts of the code.
 * \details     Usage:
 *                see e.g. src/init/density.c
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 04.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifdef USE_SUBCOMM_COMMUNICATOR
#define MYCOMMUNICATOR SubComm
#define MyThisTask SubThisTask
#define MyNTask SubNTask
#else /* #ifdef USE_SUBCOMM_COMMUNICATOR */
#define MYCOMMUNICATOR MPI_COMM_WORLD
#define MyThisTask ThisTask
#define MyNTask NTask
#endif /* #ifdef USE_SUBCOMM_COMMUNICATOR #else */

#define EXTRA_SPACE 16384

typedef struct datanodelist datanodelist;
typedef struct data_partlist data_partlist;

static size_t ExportSpace;
static size_t MinSpace;
static int NextParticle;
static int Nexport, Nimport;
static int NexportNodes, NimportNodes;
static long long SumNexport;
static int *NodeDataIn;
static int *NodeDataGet;

static char callorigin[1000];

#ifdef USE_DSDE
static void generic_prepare_import_counts_ibarrier(void);
#endif /* #ifdef USE_DSDE */

#ifdef USE_INLINED_IBARRIER
static void generic_prepare_import_counts_inlined_ibarrier(void);
#endif /* #ifdef USE_INLINED_IBARRIER */

#define generic_set_MaxNexport(...)                     \
  {                                                     \
    generic_set_info(__FUNCTION__, __FILE__, __LINE__); \
  }

/*! \brief This function determines how much buffer space we may use based on
 *         the memory that is locally still free, and it computes how much
 *         memory may at most be needed to process a single particle. We will
 *         only continue with a particle if this can still be safely processed.
 */
static void generic_set_info(const char *func, const char *file, int line)
{
  ExportSpace = 0.3 * (FreeBytes); /* we just grab at most 30% of the still available memory here */
  ExportSpace /= NUM_THREADS;
  ExportSpace -= NumPart * sizeof(int); /* to account for the neighbor list buffer that every thread allocated */

  /* make the size a multiple both of data_partlist and datanodelist */
  ExportSpace /= (sizeof(data_partlist) * sizeof(datanodelist));
  ExportSpace *= (sizeof(data_partlist) * sizeof(datanodelist));

  MinSpace =
      (MyNTask - 1) * (sizeof(data_partlist) + sizeof(data_in) + sizeof(data_out)) + NTopleaves * (sizeof(datanodelist) + sizeof(int));

  sprintf(callorigin, "%s|%d|", file, line);

#ifdef VERBOSE
  mpi_printf(
      "GENERIC: function %s(), file %s, line %d: MinSpace = %g MB  NTopleaves = %d  ExportSpace = %g MB sizeof(data_in)=%d "
      "sizeof(data_out)=%d\n",
      func, file, line, MinSpace / (1024.0 * 1024.0), NTopleaves, ExportSpace / (1024.0 * 1024.0)),
      (int)sizeof(data_in), (int)sizeof(data_out);
#endif /* #ifdef VERBOSE */

  if(ExportSpace < MinSpace)
    terminate(
        "Bummer. Can't even safely process a single particle for the available memory. FreeBytes=%lld  ExportSpace=%lld  "
        "MinSpace=%lld  MyNTask=%d  NTopleaves=%d",
        (long long)FreeBytes, (long long)ExportSpace, (long long)MinSpace, MyNTask, NTopleaves);
}

/*! \brief This function does the memory allocation at the beginning of a loop
 *         over the remaining local particles. The fields PartList[] and
 *         NodeList[] share the buffer space of size "ExportSpace" (in bytes).
 *         Here PartList will be filled in from the beginning, while NodeList
 *         will be filled in from the end. Since we do not know a priory the
 *         relative share of these two fields, we can make optimum use of
 *         the available space in this way.
 */
static void generic_alloc_partlist_nodelist_ngblist_threadbufs(void)
{
  for(int i = 0; i < NUM_THREADS; i++)
    {
      Thread[i].Nexport      = 0;
      Thread[i].NexportNodes = 0;
      Thread[i].ExportSpace  = ExportSpace;
      Thread[i].InitialSpace = ExportSpace;
      Thread[i].ItemSize     = (sizeof(data_partlist) + sizeof(data_in) + sizeof(data_out));

      Thread[i].PartList = (struct data_partlist *)mymalloc_movable_g(&Thread[i].PartList, "PartList", ExportSpace);
      /* note: the NodeList array will be attached to the end of this buffer, growing backwards */
      /* Thread[i].NodeList = (struct datanodelist *) (((char *) Thread[i].PartList) + InitialSpace);
       */
      Thread[i].Ngblist    = (int *)mymalloc_movable_g(&Thread[i].Ngblist, "Ngblist", NumPart * sizeof(int));
      Thread[i].R2list     = (double *)mymalloc_movable_g(&Thread[i].R2list, "R2list", NumPart * sizeof(double));
      Thread[i].Exportflag = Exportflag + i * ((((MyNTask - 1) / 16) + 1) * 16);
    }
}

/*! \brief The corresponding deallocation routine.
 */
static void generic_free_partlist_nodelist_ngblist_threadbufs(void)
{
  for(int i = NUM_THREADS - 1; i >= 0; i--)
    {
      myfree(Thread[i].R2list);
      myfree(Thread[i].Ngblist);
      myfree(Thread[i].PartList);
      Thread[i].R2list   = NULL;
      Thread[i].Ngblist  = NULL;
      Thread[i].PartList = NULL;
    }
}

static void generic_prepare_export_counts(void)
{
  for(int j = 0; j < MyNTask; j++)
    {
      Send[j].Count      = 0;
      Send[j].CountNodes = 0;
    }

  Nexport      = 0;
  NexportNodes = 0;

  for(int i = 0; i < NUM_THREADS; i++)
    {
      for(int j = 0; j < Thread[i].Nexport; j++)
        Send[Thread[i].PartList[j].Task].Count++;

      struct datanodelist *nodelist = (struct datanodelist *)(((char *)Thread[i].PartList) + Thread[i].InitialSpace);

      for(int j = 0; j < Thread[i].NexportNodes; j++)
        Send[nodelist[-1 - j].Task].CountNodes++;

      Nexport += Thread[i].Nexport;
      NexportNodes += Thread[i].NexportNodes;
    }

  SumNexport += Nexport;
}

/*! \brief Establishes the Recv counts from the Send counts (effectively a big
 *         transpose).
 */
static void generic_prepare_import_counts(void)
{
  /* our standard approach for this is to use an all-to-all communication. For very large processor counts,
   * this in principle becomes inefficient since mostly zeros need to be communicated.
   * we have also two option experimental communication routines that use a sparse=communication pattern instead.
   */
#ifdef USE_DSDE
  generic_prepare_import_counts_ibarrier();
#else /* #ifdef USE_DSDE */
#ifdef USE_INLINED_IBARRIER
  generic_prepare_import_counts_inlined_ibarrier();
#else  /* #ifdef USE_INLINED_IBARRIER */
  /* the default */
  MPI_Alltoall(Send, sizeof(struct send_recv_counts), MPI_BYTE, Recv, sizeof(struct send_recv_counts), MPI_BYTE, MYCOMMUNICATOR);
#endif /* #ifdef USE_INLINED_IBARRIER #else */
#endif /* #ifdef USE_DSDE #else */
}

/*! \brief Initializes offset tables that we need for the communication.
 */
static void generic_prepare_export_offsets(void)
{
  Send_offset[0]       = 0;
  Send_offset_nodes[0] = 0;

  for(int j = 1; j < MyNTask; j++)
    {
      Send_offset[j]       = Send_offset[j - 1] + Send[j - 1].Count;
      Send_offset_nodes[j] = Send_offset_nodes[j - 1] + Send[j - 1].CountNodes;
    }
}

/*! \brief Organizes the particle and node data for export in contiguous
 *         memory regions.
 */
static void generic_prepare_particle_data_for_export(void)
{
  int *rel_node_index = (int *)mymalloc_g("rel_node_index", MyNTask * sizeof(int));

  for(int j = 0; j < MyNTask; j++)
    {
      Send[j].Count      = 0;
      Send[j].CountNodes = 0;
      rel_node_index[j]  = 0;
    }

  for(int i = 0; i < NUM_THREADS; i++)
    {
      struct datanodelist *nodelist = (struct datanodelist *)(((char *)Thread[i].PartList) + Thread[i].InitialSpace);

      for(int j = 0, jj = 0; j < Thread[i].Nexport; j++)
        {
          int task = Thread[i].PartList[j].Task;
          int off  = Send_offset[task] + Send[task].Count++;

          int target = Thread[i].PartList[j].Index;

          particle2in(&DataIn[off], target, rel_node_index[task]);

          if(j < Thread[i].Nexport - 1)
            if(Thread[i].PartList[j].Index == Thread[i].PartList[j + 1].Index)
              continue;

          while(jj < Thread[i].NexportNodes && Thread[i].PartList[j].Index == nodelist[-1 - jj].Index)
            {
              int task = nodelist[-1 - jj].Task;
              int off  = Send_offset_nodes[task] + Send[task].CountNodes++;

              NodeDataIn[off] = nodelist[-1 - jj].Node;

              rel_node_index[task]++;
              jj++;
            }
        }
    }

  myfree(rel_node_index);
}

/*! \brief Driver routine to process the results that we obtained for a
 *         particle from a remote processor by working on it with the supplied
 *         out2particle() routine.
 */
static void generic_add_results_to_local(void)
{
  for(int j = 0; j < MyNTask; j++)
    Send[j].Count = 0;

  for(int i = 0; i < NUM_THREADS; i++)
    for(int j = 0; j < Thread[i].Nexport; j++)
      {
        int task = Thread[i].PartList[j].Task;
        int off  = Send_offset[task] + Send[task].Count++;

        int target = Thread[i].PartList[j].Index;

        out2particle(&DataOut[off], target, MODE_IMPORTED_PARTICLES);
      }
}

/*! \brief This function is called in the actual tree walk routine to find out
 *         how the number and starting index of the section in the node-list
 *         that needs to be processed for the imported particle.
 */
static void generic_get_numnodes(int target, int *numnodes, int **firstnode)
{
  if(target == Nimport - 1)
    *numnodes = NimportNodes - DataGet[target].Firstnode;
  else
    *numnodes = DataGet[target + 1].Firstnode - DataGet[target].Firstnode;

  *firstnode = &NodeDataGet[DataGet[target].Firstnode];
}

/*! \brief Calculates how many space we need to allocate to safely process a
 *         certain number of nodes and particles that are imported.
 */
static size_t generic_calc_import_storage(int nimport, int nimportnodes)
{
  size_t needed = nimport * sizeof(data_in) + nimportnodes * sizeof(int) + nimport * sizeof(data_out);

  /* add some extra space to not go to the last byte */
  needed += EXTRA_SPACE;

  return needed;
}

/*! \brief This routine carries out the communication step in several phases
 *         if needed.
 */
static void generic_multiple_phases(void (*kernel)(void))
{
  int ncycles;

  for(int ngrpstart = 1; ngrpstart < (1 << PTask); ngrpstart += ncycles)
    {
      /* now decide how many cycles we can process in this iteration */
      ncycles = (1 << PTask) - ngrpstart;

      do
        {
          Nimport      = 0;
          NimportNodes = 0;

          for(int ngrp = ngrpstart; ngrp < ngrpstart + ncycles; ngrp++)
            {
              int recvTask = MyThisTask ^ ngrp;

              if(recvTask < MyNTask)
                {
                  if(Recv[recvTask].Count > 0)
                    {
                      Nimport += Recv[recvTask].Count;
                      NimportNodes += Recv[recvTask].CountNodes;
                    }
                }
            }

          int flag = 0, flagall;

          if(generic_calc_import_storage(Nimport, NimportNodes) > FreeBytes)
            flag = 1;

          MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, MYCOMMUNICATOR);

          if(flagall)
            ncycles /= 2;
          else
            break;
        }
      while(ncycles > 0);

      if(ncycles == 0)
        terminate(
            "Seems like we can't even do one cycle: ncycles=%d  ngrpstart=%d  Nimport=%d  NimportNodes=%d  FreeBytes=%lld  needed "
            "storage=%lld",
            ncycles, ngrpstart, Nimport, NimportNodes, (long long)FreeBytes,
            (long long)generic_calc_import_storage(Nimport, NimportNodes));

      if(ngrpstart == 1 && ncycles != ((1 << PTask) - ngrpstart) && MyThisTask == 0)
        warn("need multiple import/export phases to avoid memory overflow");

      /* now allocated the import and results buffers */

      DataGet     = (data_in *)mymalloc_movable_g(&DataGet, "DataGet", Nimport * sizeof(data_in));
      NodeDataGet = (int *)mymalloc_movable_g(&NodeDataGet, "NodeDataGet", NimportNodes * sizeof(int));
      DataResult  = (data_out *)mymalloc_movable_g(&DataResult, "DataResult", Nimport * sizeof(data_out));

      Nimport      = 0;
      NimportNodes = 0;

      /* exchange particle data */
      for(int ngrp = ngrpstart; ngrp < ngrpstart + ncycles; ngrp++)
        {
          int recvTask = MyThisTask ^ ngrp;

          if(recvTask < MyNTask)
            {
              if(Send[recvTask].Count > 0 || Recv[recvTask].Count > 0)
                {
                  size_t len = sizeof(data_in);

                  /* get the particles */
                  MPI_Sendrecv(&DataIn[Send_offset[recvTask]], Send[recvTask].Count * len, MPI_BYTE, recvTask, TAG_HYDRO_A,
                               &DataGet[Nimport], Recv[recvTask].Count * len, MPI_BYTE, recvTask, TAG_HYDRO_A, MYCOMMUNICATOR,
                               MPI_STATUS_IGNORE);

                  /* get the nodes */
                  MPI_Sendrecv(&NodeDataIn[Send_offset_nodes[recvTask]], Send[recvTask].CountNodes, MPI_INT, recvTask, TAG_GRAV_B,
                               &NodeDataGet[NimportNodes], Recv[recvTask].CountNodes, MPI_INT, recvTask, TAG_GRAV_B, MYCOMMUNICATOR,
                               MPI_STATUS_IGNORE);

                  for(int k = 0; k < Recv[recvTask].Count; k++)
                    DataGet[Nimport + k].Firstnode += NimportNodes;

                  Nimport += Recv[recvTask].Count;
                  NimportNodes += Recv[recvTask].CountNodes;
                }
            }
        }

      /* now do the actual work for the imported points */
      kernel();

      /* send the results */
      Nimport      = 0;
      NimportNodes = 0;

      for(int ngrp = ngrpstart; ngrp < ngrpstart + ncycles; ngrp++)
        {
          int recvTask = MyThisTask ^ ngrp;
          if(recvTask < MyNTask)
            {
              if(Send[recvTask].Count > 0 || Recv[recvTask].Count > 0)
                {
                  size_t len = sizeof(data_out);

                  /* exchange the results */
                  MPI_Sendrecv(&DataResult[Nimport], Recv[recvTask].Count * len, MPI_BYTE, recvTask, TAG_HYDRO_B,
                               &DataOut[Send_offset[recvTask]], Send[recvTask].Count * len, MPI_BYTE, recvTask, TAG_HYDRO_B,
                               MYCOMMUNICATOR, MPI_STATUS_IGNORE);

                  Nimport += Recv[recvTask].Count;
                  NimportNodes += Recv[recvTask].CountNodes;
                }
            }
        }

      myfree(DataResult);
      myfree(NodeDataGet);
      myfree(DataGet);
    }
}

/*! \brief This function deals with the communication step, and then processes
 *         the imported particles, and finally computes the results back. If
 *         there is not enough memory available to hold all the data sent to
 *         us from other processors, we process the incoming data in multiple
 *         stages, which will always be possible.
 */
static void generic_exchange(void (*kernel)(void))
{
  /* set up Sendcount table */
  generic_prepare_export_counts();

  /* do the all-to-all exchange so that we have the Recvcount table as well */
  generic_prepare_import_counts();

  /* prepare offsets in export tables */
  generic_prepare_export_offsets();

  /* allocate particle data buffers */
  DataIn     = (data_in *)mymalloc_movable_g(&DataIn, "DataIn", Nexport * sizeof(data_in));
  NodeDataIn = (int *)mymalloc_movable_g(&NodeDataIn, "NodeDataIn", NexportNodes * sizeof(int));
  DataOut    = (data_out *)mymalloc_movable_g(&DataOut, "DataOut", Nexport * sizeof(data_out));

  /* prepare particle data for export */
  generic_prepare_particle_data_for_export();

  /* export particles and process them, if needed in several installments */
  generic_multiple_phases(kernel);

  /* add the results to the local particles */
  generic_add_results_to_local();

  myfree(DataOut);
  myfree(NodeDataIn);
  myfree(DataIn);
}

/* \brief Implements a repeated loop over the local particles in the list,
 *        processing them with the local kernel function, until we're done or
 *        the export buffer is full. Then we exchange the data, and process
 *        the imported ones with the provided kernel. We repeat if neeed until
 *        all processors are done.
 */
static int generic_comm_pattern(int nactive, void (*kernel_loc)(void), void (*kernel_imp)(void))
{
  int ndone_flag, ndone, iter = 0;

  SumNexport = 0; /* can be queried as a book-keeping variable */

  NextParticle = 0; /* first particle index for this task */

  do
    {
      iter++;

      /* allocate buffers to arrange communication */
      generic_alloc_partlist_nodelist_ngblist_threadbufs();

      /* do local particles */
      kernel_loc();

      /* do all necessary bookkeeping, data exchange, and processing of imported particles */
      generic_exchange(kernel_imp);

      /* free the rest of the buffers */
      generic_free_partlist_nodelist_ngblist_threadbufs();

      /* check whether we are done */
      if(NextParticle >= nactive)
        ndone_flag = 1;
      else
        ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MYCOMMUNICATOR);
    }
  while(ndone < MyNTask);

  return iter;
}

/*! \brief Same as generic_comm_pattern but you can pass the indices of the
 *         particles to be processed.
 */
static int generic_comm_pattern_for_given_particles(int nactive, int indices[], void (*kernel_loc)(int, int *),
                                                    void (*kernel_imp)(void))
{
  int ndone_flag, ndone, iter = 0;

  SumNexport = 0; /* can be queried as a book-keeping variable */

  NextParticle = 0; /* first particle index for this task */

  do
    {
      iter++;

      /* allocate buffers to arrange communication */
      generic_alloc_partlist_nodelist_ngblist_threadbufs();

      /* do local particles */
      kernel_loc(nactive, indices);

      /* do all necessary bookkeeping, data exchange, and processing of imported particles */
      generic_exchange(kernel_imp);

      /* free the rest of the buffers */
      generic_free_partlist_nodelist_ngblist_threadbufs();

      /* check whether we are done */
      if(NextParticle >= nactive)
        ndone_flag = 1;
      else
        ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MYCOMMUNICATOR);
    }
  while(ndone < MyNTask);

  return iter;
}

#ifdef USE_INLINED_IBARRIER
/*! \brief Can replace
 *         MPI_Alltoall(Send, sizeof(struct send_recv_counts), MPI_INT, Recv,
 *           sizeof(struct send_recv_counts), MPI_INT, MYCOMMUNICATOR);
 *         with a space communication pattern that effectively involves a
 *         home-grown non-blocking barrier to establish that we can stop
 *         listening.
 */
static void generic_prepare_import_counts_inlined_ibarrier(void)
{
  int nLevels         = my_fls(MyNTask - 1);
  int received_levels = 0, sent_levels = 0;

  int *stagelist = (int *)mymalloc("stagelist", nLevels * sizeof(int));
  for(int j = 0; j < nLevels; j++)
    stagelist[j] = j;

  MPI_Request *level_requests = (MPI_Request *)mymalloc("level_requests", nLevels * sizeof(MPI_Request));

  MPI_Request *requests = (MPI_Request *)mymalloc("requests", MyNTask * sizeof(MPI_Request));
  int n_requests        = 0;

  for(int j = 0; j < MyNTask; j++)
    {
      if(Send[j].Count > 0)
        MPI_Issend(&Send[j], sizeof(struct send_recv_counts), MPI_BYTE, j, TAG_N, MYCOMMUNICATOR, &requests[n_requests++]);

      Recv[j].Count      = 0;
      Recv[j].CountNodes = 0;
    }

  int barrier_active = 0;

  while(1)
    {
      int flag;
      MPI_Status status;

      MPI_Iprobe(MPI_ANY_SOURCE, TAG_N, MYCOMMUNICATOR, &flag, &status);

      if(flag)
        {
          int source = status.MPI_SOURCE;
          int tag    = status.MPI_TAG;

          MPI_Recv(&Recv[source], sizeof(struct send_recv_counts), MPI_BYTE, source, tag, MYCOMMUNICATOR, MPI_STATUS_IGNORE);
        }

      MPI_Iprobe(MPI_ANY_SOURCE, TAG_BARRIER, MYCOMMUNICATOR, &flag, &status);

      if(flag)
        {
          int source = status.MPI_SOURCE;
          int tag    = status.MPI_TAG;

          int stage;
          MPI_Recv(&stage, 1, MPI_INT, source, tag, MYCOMMUNICATOR, MPI_STATUS_IGNORE);
          received_levels |= (1 << stage);
        }

      if(barrier_active)
        {
          for(int stage = 0; stage < nLevels; stage++)
            if(!(sent_levels & (1 << stage)))
              {
                int mask = ((1 << stage) - 1);

                if((mask & received_levels) == mask)
                  {
                    sent_levels |= (1 << stage);

                    int target = (MyThisTask + (1 << stage)) % MyNTask;

                    MPI_Issend(&stagelist[stage], 1, MPI_INT, target, TAG_BARRIER, MYCOMMUNICATOR, &level_requests[stage]);
                  }
              }

          if(received_levels == ((1 << nLevels) - 1) && send_levels == ((1 << nLevels) - 1))
            break;
        }
      else
        {
          MPI_Testall(n_requests, requests, &flag, MPI_STATUSES_IGNORE);

          if(flag)
            barrier_active = 1;
        }
    }

  MPI_Waitall(nLevels, level_requests, MPI_STATUSES_IGNORE); /* as we are going to free stagelist */

  myfree(requests);
  myfree(level_requests);
  myfree(stagelist);
}
#endif /* #ifdef USE_INLINED_IBARRIER */

#ifdef USE_DSDE
/*! \brief Can replace
 *         MPI_Alltoall(Send, sizeof(struct send_recv_counts), MPI_INT, Recv,
 *           sizeof(struct send_recv_counts), MPI_INT, MYCOMMUNICATOR);
 *         with a space communication pattern that involves a non-blocking
 *         barrier (requires MPI-3.0).
 */
static int generic_prepare_import_counts_ibarrier(void)
{
  MPI_Request barrier_request;
  MPI_Request *requests = (MPI_Request *)mymalloc_movable(&requests, "requests", MyNTask * sizeof(MPI_Request));
  int n_requests        = 0;

  for(int j = 0; j < MyNTask; j++)
    {
      if(Send[j].Count > 0)
        MPI_Issend(&Send[j], sizeof(struct send_recv_counts), MPI_BYTE, j, TAG_N, MYCOMMUNICATOR, &requests[n_requests++]);

      Recv[j].Count      = 0;
      Recv[j].CountNodes = 0;
    }

  int barrier_active = 0;

  while(1)
    {
      int flag;
      MPI_Status status;

      MPI_Iprobe(MPI_ANY_SOURCE, TAG_N, MYCOMMUNICATOR, &flag, &status);

      if(flag)
        {
          int source = status.MPI_SOURCE;
          int tag    = status.MPI_TAG;

          int count;
          MPI_Get_count(&status, MPI_BYTE, &count);

          if(tag == TAG_N && source != MyThisTask)
            {
              if(count != 8)
                terminate("count=%d\n", count);

              MPI_Recv(&Recv[source], sizeof(struct send_recv_counts), MPI_BYTE, source, tag, MYCOMMUNICATOR, MPI_STATUS_IGNORE);
            }
        }

      if(barrier_active)
        {
          int flag2;

          MPI_Test(&barrier_request, &flag2, &status);

          if(flag2 != 0)
            break;
        }
      else
        {
          MPI_Testall(n_requests, requests, &flag, MPI_STATUSES_IGNORE);

          if(flag)
            {
              barrier_active = 1;

              MPI_Ibarrier(MYCOMMUNICATOR, &barrier_request);
            }
        }
    }

  myfree(requests);
}
#endif /* #ifdef USE_DSDE */
