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
 * \file        src/domain_sort_kernels.c
 * \date        05/2018
 * \brief       Comparison and sorting functions for Peano-Hilbert data.
 * \details     contains functions:
 *                int domain_compare_count(const void *a, const void *b)
 *                int domain_compare_key(const void *a, const void *b)
 *                static void msort_domain_with_tmp(struct
 *                  domain_peano_hilbert_data *b, size_t n, struct
 *                  domain_peano_hilbert_data *t)
 *                void mysort_domain(void *b, size_t n, size_t s)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 04.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../mesh/voronoi/voronoi.h"
#include "domain.h"

/*! \brief Comparison function for domain_count_data objects.
 *
 *  Compares the variable count.
 *
 *  \param[in] a Pointer to first domain_count_data object.
 *  \param[in] b Pointer to second domain_count_data object.
 *
 *  \return 1 if b>a; -1 if a>b; otherwise 0.
 */
int domain_compare_count(const void *a, const void *b)
{
  if(((struct domain_count_data *)a)->count > (((struct domain_count_data *)b)->count))
    return -1;

  if(((struct domain_count_data *)a)->count < (((struct domain_count_data *)b)->count))
    return +1;

  return 0;
}

/*! \brief Comparison function for domain_peano_hilbert_data objects.
 *
 *  Compares element key.
 *
 *  \param[in] a Pointer to first domain_peano_hilbert_data object.
 *  \param[in] b Pointer to second domain_peano_hilbert_data object.
 *
 *  \return 1 if b>a; -1 if a>b; otherwise 0.
 */
int domain_compare_key(const void *a, const void *b)
{
  if(((struct domain_peano_hilbert_data *)a)->key < (((struct domain_peano_hilbert_data *)b)->key))
    return -1;

  if(((struct domain_peano_hilbert_data *)a)->key > (((struct domain_peano_hilbert_data *)b)->key))
    return +1;

  return 0;
}

/*! \brief Customized mergesort sorting routine, requires temporary array.
 *
 *  \param[in, out] b domain_peano_hilbert data array that is to be sorted.
 *  \param[in] n Number of elements in array.
 *  \param[in, out] t Temporary domain_peano_hilbert data array.
 *
 *  \return void
 */
static void msort_domain_with_tmp(struct domain_peano_hilbert_data *b, size_t n, struct domain_peano_hilbert_data *t)
{
  struct domain_peano_hilbert_data *tmp;
  struct domain_peano_hilbert_data *b1, *b2;
  size_t n1, n2;

  if(n <= 1)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = b;
  b2 = b + n1;

  msort_domain_with_tmp(b1, n1, t);
  msort_domain_with_tmp(b2, n2, t);

  tmp = t;

  while(n1 > 0 && n2 > 0)
    {
      if(b1->key <= b2->key)
        {
          --n1;
          *tmp++ = *b1++;
        }
      else
        {
          --n2;
          *tmp++ = *b2++;
        }
    }

  if(n1 > 0)
    memcpy(tmp, b1, n1 * sizeof(struct domain_peano_hilbert_data));

  memcpy(b, t, (n - n2) * sizeof(struct domain_peano_hilbert_data));
}

/*! \brief Customized mergesort sorting routine.
 *
 *  This function tends to work slightly faster than a call of qsort() for
 *  this particular list, at least on most platforms.
 *
 *  \param[in, out] b domain_peano_hilbert data array that is to be sorted.
 *  \param[in] n Number of elements.
 *  \param[in] s Size of structure.
 *
 *  \return void
 */
void mysort_domain(void *b, size_t n, size_t s)
{
  const size_t size = n * s;
  struct domain_peano_hilbert_data *tmp;

  tmp = (struct domain_peano_hilbert_data *)mymalloc("tmp", size);

  msort_domain_with_tmp((struct domain_peano_hilbert_data *)b, n, tmp);

  myfree(tmp);
}
