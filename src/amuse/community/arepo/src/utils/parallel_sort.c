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
 * \file        src/utils/parallel_sort.c
 * \date        05/2018
 * \brief       MPI parallel sorting routine.
 * \details     contains functions:
 *                int parallel_sort_indirect_compare(const void *a,
 *                  const void *b)
 *                double parallel_sort(void *base, size_t nmemb, size_t size,
 *                  int (*compar) (const void *, const void *))
 *                double parallel_sort_comm(void *base, size_t nmemb, size_t
 *                  size, int (*compar) (const void *, const void *),
 *                  MPI_Comm comm)
 *                static void get_local_rank(char *element, size_t
 *                  tie_braking_rank, char *base, size_t nmemb, size_t size,
 *                  size_t noffs_thistask, long long left, long long right,
 *                  size_t * loc, int (*compar) (const void *, const void *))
 *                static void check_local_rank(char *element, size_t
 *                  tie_braking_rank, char *base, size_t nmemb, size_t size,
 *                  size_t noffs_thistask, long long left, long long right,
 *                  size_t loc, int (*compar) (const void *, const void *))
 *                static void serial_sort(char *base, size_t nmemb, size_t
 *                  size, int (*compar) (const void *, const void *))
 *                static void msort_serial_with_tmp(char *base, size_t n,
 *                  size_t s, int (*compar) (const void *, const void *),
 *                  char *t)
 *                void parallel_sort_test_order(char *base, size_t nmemb,
 *                  size_t size, int (*compar) (const void *, const void *))
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 21.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#define TRANSFER_SIZE_LIMIT 1000000000
#define MAX_ITER_PARALLEL_SORT 500

/* Note: For gcc-4.1.2, I found that the compiler produces incorrect code for this routune if optimization level O1 or higher is used.
 *       In  gcc-4.3.4, this problem is absent.
 */

#define TAG_TRANSFER 100

static void serial_sort(char *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *));
static void msort_serial_with_tmp(char *base, size_t n, size_t s, int (*compar)(const void *, const void *), char *t);
static void get_local_rank(char *element, size_t tie_braking_rank, char *base, size_t nmemb, size_t size, size_t noffs_thistask,
                           long long left, long long right, size_t *loc, int (*compar)(const void *, const void *));

static int (*comparfunc)(const void *, const void *);
static char *median_element_list;
static size_t element_size;

/*! \brief Wrapper for comparison of  two elements.
 *
 *  \param[in] a First element.
 *  \param[in] b Second element.
 *
 *  \return (-1,0,+1) -1 if a < b.
 */
int parallel_sort_indirect_compare(const void *a, const void *b)
{
  return (*comparfunc)(median_element_list + *((int *)a) * element_size, median_element_list + *((int *)b) * element_size);
}

/*! \brief Main function to perform a parallel sort.
 *
 *   Using MPI_COMM_WORLD as communicator.
 *
 *  \param[in, out] base Array to be sorted.
 *  \param nmemb Number of entries in array.
 *  \param[in] size Size of an element in array to be sorted.
 *  \param[in] compar Comparison function.
 *
 *  \return Time it took to sort array.
 */
double parallel_sort(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *))
{
  return parallel_sort_comm(base, nmemb, size, compar, MPI_COMM_WORLD);
}

/*! \brief Function to perform a parallel sort with specified MPI communicator.
 *
 *  \param[in, out] base Array to be sorted.
 *  \param[in] nmemb Number of entries in array.
 *  \param[in] size Size of an element in array to be sorted.
 *  \param[in] compar Comparison function.
 *  \param[in] comm MPI communicator.
 *
 *  \return Time it took to sort array.
 */
double parallel_sort_comm(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *), MPI_Comm comm)
{
  int i, j, ranks_not_found, Local_ThisTask, Local_NTask, Local_PTask, Color, new_max_loc;
  size_t tie_braking_rank, new_tie_braking_rank, rank;
  MPI_Comm MPI_CommLocal;

  double ta = second();

  /* do a serial sort of the local data up front */
  serial_sort((char *)base, nmemb, size, compar);

  /* we create a communicator that contains just those tasks with nmemb > 0. This makes
   *  it easier to deal with CPUs that do not hold any data.
   */
  if(nmemb)
    Color = 1;
  else
    Color = 0;

  MPI_Comm_split(comm, Color, ThisTask, &MPI_CommLocal);
  MPI_Comm_rank(MPI_CommLocal, &Local_ThisTask);
  MPI_Comm_size(MPI_CommLocal, &Local_NTask);

  if(Local_NTask > 1 && Color == 1)
    {
      for(Local_PTask = 0; Local_NTask > (1 << Local_PTask); Local_PTask++)
        ;

      size_t *nlist = (size_t *)mymalloc("nlist", Local_NTask * sizeof(size_t));
      size_t *noffs = (size_t *)mymalloc("noffs", Local_NTask * sizeof(size_t));

      MPI_Allgather(&nmemb, sizeof(size_t), MPI_BYTE, nlist, sizeof(size_t), MPI_BYTE, MPI_CommLocal);

      for(i = 1, noffs[0] = 0; i < Local_NTask; i++)
        noffs[i] = noffs[i - 1] + nlist[i - 1];

      char *element_guess              = mymalloc("element_guess", Local_NTask * size);
      size_t *element_tie_braking_rank = mymalloc("element_tie_braking_rank", Local_NTask * sizeof(size_t));
      size_t *desired_glob_rank        = mymalloc("desired_glob_rank", Local_NTask * sizeof(size_t));
      size_t *current_glob_rank        = mymalloc("current_glob_rank", Local_NTask * sizeof(size_t));
      size_t *current_loc_rank         = mymalloc("current_loc_rank", Local_NTask * sizeof(size_t));
      long long *range_left            = mymalloc("range_left", Local_NTask * sizeof(long long));
      long long *range_right           = mymalloc("range_right", Local_NTask * sizeof(long long));
      int *max_loc                     = mymalloc("max_loc", Local_NTask * sizeof(int));

      size_t *list                         = mymalloc("list", Local_NTask * sizeof(size_t));
      size_t *range_len_list               = mymalloc("range_len_list", Local_NTask * sizeof(long long));
      char *median_element                 = mymalloc("median_element", size);
      median_element_list                  = mymalloc("median_element_list", Local_NTask * size);
      size_t *tie_braking_rank_list        = mymalloc("tie_braking_rank_list", Local_NTask * sizeof(size_t));
      int *index_list                      = mymalloc("index_list", Local_NTask * sizeof(int));
      int *max_loc_list                    = mymalloc("max_loc_list", Local_NTask * sizeof(int));
      size_t *source_range_len_list        = mymalloc("source_range_len_list", Local_NTask * sizeof(long long));
      size_t *source_tie_braking_rank_list = mymalloc("source_tie_braking_rank_list", Local_NTask * sizeof(long long));
      char *source_median_element_list     = mymalloc("source_median_element_list", Local_NTask * size);
      char *new_element_guess              = mymalloc("new_element_guess", size);

      for(i = 0; i < Local_NTask - 1; i++)
        {
          desired_glob_rank[i] = noffs[i + 1];
          current_glob_rank[i] = 0;
          range_left[i]        = 0;     /* first element that it can be */
          range_right[i]       = nmemb; /* first element that it can not be */
        }

      /* now we determine the first split element guess, which is the same for all divisions in the first iteration */

      /* find the median of each processor, and then take the median among those values.
       * This should work reasonably well even for extremely skewed distributions
       */
      long long range_len = range_right[0] - range_left[0];

      if(range_len >= 1)
        {
          long long mid = (range_left[0] + range_right[0]) / 2;
          memcpy(median_element, (char *)base + mid * size, size);
          tie_braking_rank = mid + noffs[Local_ThisTask];
        }

      MPI_Gather(&range_len, sizeof(long long), MPI_BYTE, range_len_list, sizeof(long long), MPI_BYTE, 0, MPI_CommLocal);
      MPI_Gather(median_element, size, MPI_BYTE, median_element_list, size, MPI_BYTE, 0, MPI_CommLocal);
      MPI_Gather(&tie_braking_rank, sizeof(size_t), MPI_BYTE, tie_braking_rank_list, sizeof(size_t), MPI_BYTE, 0, MPI_CommLocal);

      if(Local_ThisTask == 0)
        {
          for(j = 0; j < Local_NTask; j++)
            max_loc_list[j] = j;

          /* eliminate the elements that are undefined because the corresponding CPU has zero range left */
          int nleft = Local_NTask;

          for(j = 0; j < nleft; j++)
            {
              if(range_len_list[j] < 1)
                {
                  range_len_list[j] = range_len_list[nleft - 1];
                  if(range_len_list[nleft - 1] >= 1 && j != (nleft - 1))
                    {
                      memcpy(median_element_list + j * size, median_element_list + (nleft - 1) * size, size);
                      memcpy(tie_braking_rank_list + j, tie_braking_rank_list + (nleft - 1), sizeof(size_t));
                      max_loc_list[j] = max_loc_list[nleft - 1];
                    }

                  nleft--;
                  j--;
                }
            }

          /* do a serial sort of the remaining elements (indirectly, so that we have the order of tie braking list as well) */
          comparfunc   = compar;
          element_size = size;
          for(j = 0; j < nleft; j++)
            index_list[j] = j;
          qsort(index_list, nleft, sizeof(int), parallel_sort_indirect_compare);

          /* now select the median of the medians */
          int mid = nleft / 2;
          memcpy(&element_guess[0], median_element_list + index_list[mid] * size, size);
          element_tie_braking_rank[0] = tie_braking_rank_list[index_list[mid]];
          max_loc[0]                  = max_loc_list[index_list[mid]];
        }

      MPI_Bcast(element_guess, size, MPI_BYTE, 0, MPI_CommLocal);
      MPI_Bcast(&element_tie_braking_rank[0], sizeof(size_t), MPI_BYTE, 0, MPI_CommLocal);
      MPI_Bcast(&max_loc[0], 1, MPI_INT, 0, MPI_CommLocal);

      for(i = 1; i < Local_NTask - 1; i++)
        {
          memcpy(element_guess + i * size, element_guess, size);
          element_tie_braking_rank[i] = element_tie_braking_rank[0];
          max_loc[i]                  = max_loc[0];
        }

      int iter = 0;

      do
        {
          for(i = 0; i < Local_NTask - 1; i++)
            {
              if(current_glob_rank[i] != desired_glob_rank[i])
                {
                  get_local_rank(element_guess + i * size, element_tie_braking_rank[i], (char *)base, nmemb, size,
                                 noffs[Local_ThisTask], range_left[i], range_right[i], &current_loc_rank[i], compar);
                }
            }

          /* now compute the global ranks by summing the local ranks */
          /* Note: the last element in current_loc_rank is not defined. It will be summed by the last processor, and stored in the last
           * element of current_glob_rank */
          MPI_Alltoall(current_loc_rank, sizeof(size_t), MPI_BYTE, list, sizeof(size_t), MPI_BYTE, MPI_CommLocal);
          for(j = 0, rank = 0; j < Local_NTask; j++)
            rank += list[j];
          MPI_Allgather(&rank, sizeof(size_t), MPI_BYTE, current_glob_rank, sizeof(size_t), MPI_BYTE, MPI_CommLocal);

          for(i = 0, ranks_not_found = 0; i < Local_NTask - 1; i++)
            {
              if(current_glob_rank[i] != desired_glob_rank[i]) /* here we're not yet done */
                {
                  ranks_not_found++;

                  if(current_glob_rank[i] < desired_glob_rank[i])
                    {
                      range_left[i] = current_loc_rank[i];

                      if(Local_ThisTask == max_loc[i])
                        range_left[i]++;
                    }

                  if(current_glob_rank[i] > desired_glob_rank[i])
                    range_right[i] = current_loc_rank[i];
                }
            }

          /* now we need to determine new element guesses */
          for(i = 0; i < Local_NTask - 1; i++)
            {
              if(current_glob_rank[i] != desired_glob_rank[i]) /* here we're not yet done */
                {
                  /* find the median of each processor, and then take the median among those values.
                   * This should work reasonably well even for extremely skewed distributions
                   */
                  source_range_len_list[i] = range_right[i] - range_left[i];

                  if(source_range_len_list[i] >= 1)
                    {
                      long long middle = (range_left[i] + range_right[i]) / 2;
                      memcpy(source_median_element_list + i * size, (char *)base + middle * size, size);
                      source_tie_braking_rank_list[i] = middle + noffs[Local_ThisTask];
                    }
                }
            }

          MPI_Alltoall(source_range_len_list, sizeof(long long), MPI_BYTE, range_len_list, sizeof(long long), MPI_BYTE, MPI_CommLocal);
          MPI_Alltoall(source_median_element_list, size, MPI_BYTE, median_element_list, size, MPI_BYTE, MPI_CommLocal);
          MPI_Alltoall(source_tie_braking_rank_list, sizeof(size_t), MPI_BYTE, tie_braking_rank_list, sizeof(size_t), MPI_BYTE,
                       MPI_CommLocal);

          if(Local_ThisTask < Local_NTask - 1)
            {
              if(current_glob_rank[Local_ThisTask] !=
                 desired_glob_rank[Local_ThisTask]) /* in this case we're not yet done for this split point */
                {
                  for(j = 0; j < Local_NTask; j++)
                    max_loc_list[j] = j;

                  /* eliminate the elements that are undefined because the corresponding CPU has zero range left */
                  int nleft = Local_NTask;

                  for(j = 0; j < nleft; j++)
                    {
                      if(range_len_list[j] < 1)
                        {
                          range_len_list[j] = range_len_list[nleft - 1];
                          if(range_len_list[nleft - 1] >= 1 && j != (nleft - 1))
                            {
                              memcpy(median_element_list + j * size, median_element_list + (nleft - 1) * size, size);
                              memcpy(tie_braking_rank_list + j, tie_braking_rank_list + (nleft - 1), sizeof(size_t));
                              max_loc_list[j] = max_loc_list[nleft - 1];
                            }

                          nleft--;
                          j--;
                        }
                    }

                  if((iter & 1))
                    {
                      int max_range, maxj;

                      for(j = 0, maxj = 0, max_range = 0; j < nleft; j++)
                        if(range_len_list[j] > max_range)
                          {
                            max_range = range_len_list[j];
                            maxj      = j;
                          }

                      /* now select the median element from the task which has the largest range */
                      memcpy(new_element_guess, median_element_list + maxj * size, size);
                      new_tie_braking_rank = tie_braking_rank_list[maxj];
                      new_max_loc          = max_loc_list[maxj];
                    }
                  else
                    {
                      /* do a serial sort of the remaining elements (indirectly, so that we have the order of tie braking list as well)
                       */
                      comparfunc   = compar;
                      element_size = size;
                      for(j = 0; j < nleft; j++)
                        index_list[j] = j;
                      qsort(index_list, nleft, sizeof(int), parallel_sort_indirect_compare);

                      /* now select the median of the medians */
                      int mid = nleft / 2;
                      memcpy(new_element_guess, median_element_list + index_list[mid] * size, size);
                      new_tie_braking_rank = tie_braking_rank_list[index_list[mid]];
                      new_max_loc          = max_loc_list[index_list[mid]];
                    }
                }
              else
                {
                  /* in order to preserve existing guesses */
                  memcpy(new_element_guess, element_guess + Local_ThisTask * size, size);
                  new_tie_braking_rank = element_tie_braking_rank[Local_ThisTask];
                  new_max_loc          = max_loc[Local_ThisTask];
                }
            }

          MPI_Allgather(new_element_guess, size, MPI_BYTE, element_guess, size, MPI_BYTE, MPI_CommLocal);
          MPI_Allgather(&new_tie_braking_rank, sizeof(size_t), MPI_BYTE, element_tie_braking_rank, sizeof(size_t), MPI_BYTE,
                        MPI_CommLocal);
          MPI_Allgather(&new_max_loc, 1, MPI_INT, max_loc, 1, MPI_INT, MPI_CommLocal);

          iter++;

          if(iter > (MAX_ITER_PARALLEL_SORT - 100) && Local_ThisTask == 0)
            {
              printf("PSORT: iter=%d: ranks_not_found=%d  Local_NTask=%d\n", iter, ranks_not_found, Local_NTask);
              myflush(stdout);
              if(iter > MAX_ITER_PARALLEL_SORT)
                terminate("can't find the split points. That's odd");
            }
        }
      while(ranks_not_found);

      myfree(new_element_guess);
      myfree(source_median_element_list);
      myfree(source_tie_braking_rank_list);
      myfree(source_range_len_list);
      myfree(max_loc_list);
      myfree(index_list);
      myfree(tie_braking_rank_list);
      myfree(median_element_list);
      myfree(median_element);

      /* At this point we have found all the elements corresponding to the desired split points */
      /* we can now go ahead and determine how many elements of the local CPU have to go to each other CPU */

      if(nmemb * size > (1LL << 31))
        terminate("currently, local data must be smaller than 2 GB");
      /* note: to restrict this limitation, the send/recv count arrays have to made 64-bit,
       * and the MPI data exchange though MPI_Alltoall has to be modified such that buffers > 2 GB become possible
       */

      int *send_count  = mymalloc("send_count", Local_NTask * sizeof(int));
      int *recv_count  = mymalloc("recv_count", Local_NTask * sizeof(int));
      int *send_offset = mymalloc("send_offset", Local_NTask * sizeof(int));
      int *recv_offset = mymalloc("recv_offset", Local_NTask * sizeof(int));

      for(i = 0; i < Local_NTask; i++)
        send_count[i] = 0;

      int target = 0;

      for(i = 0; i < nmemb; i++)
        {
          while(target < Local_NTask - 1)
            {
              int cmp = compar((char *)base + i * size, element_guess + target * size);
              if(cmp == 0)
                {
                  if(i + noffs[Local_ThisTask] < element_tie_braking_rank[target])
                    cmp = -1;
                  else if(i + noffs[Local_ThisTask] > element_tie_braking_rank[target])
                    cmp = +1;
                }
              if(cmp >= 0)
                target++;
              else
                break;
            }
          send_count[target]++;
        }

      MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, MPI_CommLocal);

      size_t nimport;

      for(j = 0, nimport = 0, recv_offset[0] = 0, send_offset[0] = 0; j < Local_NTask; j++)
        {
          nimport += recv_count[j];

          if(j > 0)
            {
              send_offset[j] = send_offset[j - 1] + send_count[j - 1];
              recv_offset[j] = recv_offset[j - 1] + recv_count[j - 1];
            }
        }

      if(nimport != nmemb)
        terminate("nimport != nmemb");

      for(j = 0; j < Local_NTask; j++)
        {
          send_count[j] *= size;
          recv_count[j] *= size;

          send_offset[j] *= size;
          recv_offset[j] *= size;
        }

      char *basetmp = mymalloc("basetmp", nmemb * size);

      /* exchange the data */
      MPI_Alltoallv(base, send_count, send_offset, MPI_BYTE, basetmp, recv_count, recv_offset, MPI_BYTE, MPI_CommLocal);

      memcpy(base, basetmp, nmemb * size);
      myfree(basetmp);

      serial_sort((char *)base, nmemb, size, compar);

      myfree(recv_offset);
      myfree(send_offset);
      myfree(recv_count);
      myfree(send_count);

      myfree(range_len_list);
      myfree(list);
      myfree(max_loc);
      myfree(range_right);
      myfree(range_left);
      myfree(current_loc_rank);
      myfree(current_glob_rank);
      myfree(desired_glob_rank);
      myfree(element_tie_braking_rank);
      myfree(element_guess);
      myfree(noffs);
      myfree(nlist);
    }

  MPI_Comm_free(&MPI_CommLocal);

  double tb = second();
  return timediff(ta, tb);
}

/*! \brief Get rank of an element.
 *
 *  \param[in] element Element of which we want the rank.
 *  \param[in] tie_braking_rank The inital global rank of this element (needed
 *             for braking ties).
 *  \param[in] base Base address of local data.
 *  \param[in] nmemb Number of elements in array.
 *  \param[in] size Size of local data.
 *  \param[in] noffs_thistask Cumulative length of data on lower tasks.
 *  \param[in] left Range of elements on local task that may hold the element.
 *  \param[in] right Range of elements on local task that may hold the element.
 *  \param[out] loc Local rank of the element.
 *  \param[in] compar User-specified  comparison function.
 *
 *  \return void
 */
static void get_local_rank(char *element, size_t tie_braking_rank, char *base, size_t nmemb, size_t size, size_t noffs_thistask,
                           long long left, long long right, size_t *loc, int (*compar)(const void *, const void *))
{
  if(right < left)
    terminate("right < left");

  if(left == 0 && right == nmemb + 1)
    {
      if(compar(base + (nmemb - 1) * size, element) < 0)
        {
          *loc = nmemb;
          return;
        }
      else if(compar(base, element) > 0)
        {
          *loc = 0;
          return;
        }
    }

  if(right == left) /* looks like we already converged to the proper rank */
    {
      *loc = left;
    }
  else
    {
      if(compar(base + (right - 1) * size, element) < 0) /* the last element is smaller, hence all elements are on the left */
        *loc = (right - 1) + 1;
      else if(compar(base + left * size, element) > 0) /* the first element is already larger, hence no element is on the left */
        *loc = left;
      else
        {
          while(right > left)
            {
              long long mid = ((right - 1) + left) / 2;

              int cmp = compar(base + mid * size, element);
              if(cmp == 0)
                {
                  if(mid + noffs_thistask < tie_braking_rank)
                    cmp = -1;
                  else if(mid + noffs_thistask > tie_braking_rank)
                    cmp = +1;
                }

              if(cmp == 0) /* element has exactly been found */
                {
                  *loc = mid;
                  break;
                }

              if((right - 1) == left) /* elements is not on this CPU */
                {
                  if(cmp < 0)
                    *loc = mid + 1;
                  else
                    *loc = mid;
                  break;
                }

              if(cmp < 0)
                {
                  left = mid + 1;
                }
              else
                {
                  if((right - 1) == left + 1)
                    {
                      if(mid != left)
                        terminate("Can't be: -->left=%lld  right=%lld\n", left, right);

                      *loc = left;
                      break;
                    }

                  right = mid;
                }
            }
        }
    }
}

/*! \brief Wrapper for serial sorting algorithm.
 *
 *  Calls a merge sort algorithm.
 *
 *  \param[in, out] base Array to be sorted.
 *  \param[in] nmemb Number of elements in array.
 *  \param[in] size Size of each element.
 *  \param[in] compar Comparison funciton.
 *
 *  \return void
 */
static void serial_sort(char *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *))
{
  size_t storage = nmemb * size;
  char *tmp      = (char *)mymalloc("tmp", storage);

  msort_serial_with_tmp(base, nmemb, size, compar, tmp);

  myfree(tmp);
}

/*! \brief Merge sort algorithm (serial).
 *
 *  \param[in, out] base Array to be sorted.
 *  \param[in] n Number of elements.
 *  \param[in] s Size of each element.
 *  \param[in] compar Comparison function.
 *  \param[in, out] t Array for temporary data storage.
 *
 *  \return void
 */
static void msort_serial_with_tmp(char *base, size_t n, size_t s, int (*compar)(const void *, const void *), char *t)
{
  char *tmp;
  char *b1, *b2;
  size_t n1, n2;

  if(n <= 1)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = base;
  b2 = base + n1 * s;

  msort_serial_with_tmp(b1, n1, s, compar, t);
  msort_serial_with_tmp(b2, n2, s, compar, t);

  tmp = t;

  while(n1 > 0 && n2 > 0)
    {
      if(compar(b1, b2) < 0)
        {
          --n1;
          memcpy(tmp, b1, s);
          tmp += s;
          b1 += s;
        }
      else
        {
          --n2;
          memcpy(tmp, b2, s);
          tmp += s;
          b2 += s;
        }
    }

  if(n1 > 0)
    memcpy(tmp, b1, n1 * s);

  memcpy(base, t, (n - n2) * s);
}

/*! \brief Test function for parallel sort.
 *
 *  \param[in] base Array to be checked.
 *  \param[in] nmemb Number of elements in array.
 *  \param[in] size Size of each element.
 *  \param[in] compar Comparison function.
 *
 *  \return void
 */
void parallel_sort_test_order(char *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *))
{
  int i, recv, send;
  size_t *nlist;

  nlist = (size_t *)mymalloc("nlist", NTask * sizeof(size_t));

  MPI_Allgather(&nmemb, sizeof(size_t), MPI_BYTE, nlist, sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0, recv = -1; i < ThisTask && nmemb > 0; i++)
    if(nlist[i] > 0)
      recv = i;

  for(i = ThisTask + 1, send = -1; nmemb > 0 && i < NTask; i++)
    if(nlist[i] > 0)
      {
        send = i;
        break;
      }

  char *element = mymalloc("element", size);

  MPI_Request requests[2];
  int nreq = 0;

  if(send >= 0)
    MPI_Isend(base + (nmemb - 1) * size, size, MPI_BYTE, send, TAG_TRANSFER, MPI_COMM_WORLD, &requests[nreq++]);

  if(recv >= 0)
    MPI_Irecv(element, size, MPI_BYTE, recv, TAG_TRANSFER, MPI_COMM_WORLD, &requests[nreq++]);

  MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE);

  if(recv >= 0)
    {
      for(i = 0; i < nmemb; i++)
        {
          if(compar(element, base + i * size) > 0)
            terminate("wrong order");
        }
    }

  myfree(element);
  myfree(nlist);
}
