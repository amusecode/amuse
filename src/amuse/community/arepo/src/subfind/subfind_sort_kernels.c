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
 * \file        src/subfind/subfind_sort_kernels.c
 * \date        05/2018
 * \brief       Comparison functions that serve as sorting kernels for various
 *              different structs used in subfind.
 * \details     contains functions:
 *                int subfind_compare_procassign_GrNr(const void *a,
 *                  const void *b)
 *                int subfind_compare_submp_GrNr_DM_Density(const void *a,
 *                  const void *b)
 *                int subfind_compare_submp_OldIndex(const void *a,
 *                  const void *b)
 *                int subfind_compare_ID_list(const void *a, const void *b)
 *                int subfind_compare_SubGroup_GrNr_SubNr(const void *a, const
 *                  void *b)
 *                int subfind_compare_dist_rotcurve(const void *a, const void
 *                  *b)
 *                int subfind_compare_rlist_mhd(const void *a, const void *b)
 *                int subfind_compare_binding_energy(const void *a, const void
 *                  *b)
 *                int subfind_compare_serial_candidates_boundlength(const void
 *                  *a, const void *b)
 *                int subfind_compare_serial_candidates_rank(const void *a,
 *                  const void *b)
 *                int subfind_compare_serial_candidates_subnr(const void *a,
 *                  const void *b)
 *                int subfind_compare_coll_candidates_subnr(const void *a,
 *                  const void *b)
 *                int subfind_compare_coll_candidates_nsubs(const void *a,
 *                  const void *b)
 *                int subfind_compare_coll_candidates_boundlength(const void
 *                  *a, const void *b)
 *                int subfind_compare_coll_candidates_rank(const void *a,
 *                  const void *b)
 *                int subfind_fof_compare_ID(const void *a, const void *b)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 11.05.2018 Prepared file for public release -- Rainer Weinberger
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
#include "subfind.h"

#ifdef SUBFIND

/*! \brief Comparison function for proc_assign_data objects.
 *
 *  Sorting kernel comparing element GrNr.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_procassign_GrNr(const void *a, const void *b)
{
  if(((struct proc_assign_data *)a)->GrNr < ((struct proc_assign_data *)b)->GrNr)
    return -1;

  if(((struct proc_assign_data *)a)->GrNr > ((struct proc_assign_data *)b)->GrNr)
    return +1;

  return 0;
}

/*! \brief Comparison function for submp_data objects.
 *
 *  Sorting kernel comparing element (most important first):
 *  GrNr, DM_Density.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b, except for DM density, where -1 if a > b
 */
int subfind_compare_submp_GrNr_DM_Density(const void *a, const void *b)
{
  if(((struct submp_data *)a)->GrNr < ((struct submp_data *)b)->GrNr)
    return -1;

  if(((struct submp_data *)a)->GrNr > ((struct submp_data *)b)->GrNr)
    return +1;

  if(((struct submp_data *)a)->DM_Density > ((struct submp_data *)b)->DM_Density)
    return -1;

  if(((struct submp_data *)a)->DM_Density < ((struct submp_data *)b)->DM_Density)
    return +1;

  return 0;
}

/*! \brief Comparison function for submp_data objects.
 *
 *  Sorting kernel comparing element OldIndex.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_submp_OldIndex(const void *a, const void *b)
{
  if(((struct submp_data *)a)->OldIndex < ((struct submp_data *)b)->OldIndex)
    return -1;

  if(((struct submp_data *)a)->OldIndex > ((struct submp_data *)b)->OldIndex)
    return +1;

  return 0;
}

/*! \brief Comparison function for id_list objects.
 *
 *  Sorting kernel comparing elements (most important first):
 *  GrNr, SubNr, Type, BindingEgy.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_ID_list(const void *a, const void *b)
{
  if(((struct id_list *)a)->GrNr < ((struct id_list *)b)->GrNr)
    return -1;

  if(((struct id_list *)a)->GrNr > ((struct id_list *)b)->GrNr)
    return +1;

  if(((struct id_list *)a)->SubNr < ((struct id_list *)b)->SubNr)
    return -1;

  if(((struct id_list *)a)->SubNr > ((struct id_list *)b)->SubNr)
    return +1;

  if(((struct id_list *)a)->Type < ((struct id_list *)b)->Type)
    return -1;

  if(((struct id_list *)a)->Type > ((struct id_list *)b)->Type)
    return +1;

  if(((struct id_list *)a)->BindingEgy < ((struct id_list *)b)->BindingEgy)
    return -1;

  if(((struct id_list *)a)->BindingEgy > ((struct id_list *)b)->BindingEgy)
    return +1;

  return 0;
}

/*! \brief Comparison function for subgroup_properties objects.
 *
 *  Sorting kernel comparing elements (most important first):
 *  GrNr and SubNr.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_SubGroup_GrNr_SubNr(const void *a, const void *b)
{
  if(((struct subgroup_properties *)a)->GrNr < ((struct subgroup_properties *)b)->GrNr)
    return -1;

  if(((struct subgroup_properties *)a)->GrNr > ((struct subgroup_properties *)b)->GrNr)
    return +1;

  if(((struct subgroup_properties *)a)->SubNr < ((struct subgroup_properties *)b)->SubNr)
    return -1;

  if(((struct subgroup_properties *)a)->SubNr > ((struct subgroup_properties *)b)->SubNr)
    return +1;

  return 0;
}

/*! \brief Comparison function for sort_r2list objects.
 *
 *  Sorting kernel comparing element r.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_dist_rotcurve(const void *a, const void *b)
{
  if(((sort_r2list *)a)->r < ((sort_r2list *)b)->r)
    return -1;

  if(((sort_r2list *)a)->r > ((sort_r2list *)b)->r)
    return +1;

  return 0;
}

/*! \brief Comparison function for variables of type double.
 *
 *  Sorting kernel.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_binding_energy(const void *a, const void *b)
{
  if(*((double *)a) > *((double *)b))
    return -1;

  if(*((double *)a) < *((double *)b))
    return +1;

  return 0;
}

/*! \brief Comparison function for cand_dat objects.
 *
 *  Sorting kernel comparing elements (most important first):
 *  bound_length and rank.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b, excpet bound length, where -1 if a > b.
 */
int subfind_compare_serial_candidates_boundlength(const void *a, const void *b)
{
  if(((struct cand_dat *)a)->bound_length > ((struct cand_dat *)b)->bound_length)
    return -1;

  if(((struct cand_dat *)a)->bound_length < ((struct cand_dat *)b)->bound_length)
    return +1;

  if(((struct cand_dat *)a)->rank < ((struct cand_dat *)b)->rank)
    return -1;

  if(((struct cand_dat *)a)->rank > ((struct cand_dat *)b)->rank)
    return +1;

  return 0;
}

/*! \brief Comparison function for cand_dat objects.
 *
 *  Sorting kernel comparing elements (most important first):
 *  rank and len.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b, except for len where -1 if a>b.
 */
int subfind_compare_serial_candidates_rank(const void *a, const void *b)
{
  if(((struct cand_dat *)a)->rank < ((struct cand_dat *)b)->rank)
    return -1;

  if(((struct cand_dat *)a)->rank > ((struct cand_dat *)b)->rank)
    return +1;

  if(((struct cand_dat *)a)->len > ((struct cand_dat *)b)->len)
    return -1;

  if(((struct cand_dat *)a)->len < ((struct cand_dat *)b)->len)
    return +1;

  return 0;
}

/*! \brief Comparison function for cand_dat objects.
 *
 *  Sorting kernel comparing element subnr.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_serial_candidates_subnr(const void *a, const void *b)
{
  if(((struct cand_dat *)a)->subnr < ((struct cand_dat *)b)->subnr)
    return -1;

  if(((struct cand_dat *)a)->subnr > ((struct cand_dat *)b)->subnr)
    return +1;

  return 0;
}

/*! \brief Comparison function for coll_cand_dat objects.
 *
 *  Sorting kernel comparing element subnr.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_coll_candidates_subnr(const void *a, const void *b)
{
  if(((struct coll_cand_dat *)a)->subnr < ((struct coll_cand_dat *)b)->subnr)
    return -1;

  if(((struct coll_cand_dat *)a)->subnr > ((struct coll_cand_dat *)b)->subnr)
    return +1;

  return 0;
}

/*! \brief Comparison function for coll_cand_dat objects.
 *
 *  Sorting kernel comparing element nsub.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_compare_coll_candidates_nsubs(const void *a, const void *b)
{
  if(((struct coll_cand_dat *)a)->nsub < ((struct coll_cand_dat *)b)->nsub)
    return -1;

  if(((struct coll_cand_dat *)a)->nsub > ((struct coll_cand_dat *)b)->nsub)
    return +1;

  return 0;
}

/*! \brief Comparison function for coll_cand_dat objects.
 *
 *  Sorting kernel comparing elements (most important first):
 *  bound_length, rank.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b, except for bound length where -1 if a > b.
 */
int subfind_compare_coll_candidates_boundlength(const void *a, const void *b)
{
  if(((struct coll_cand_dat *)a)->bound_length > ((struct coll_cand_dat *)b)->bound_length)
    return -1;

  if(((struct coll_cand_dat *)a)->bound_length < ((struct coll_cand_dat *)b)->bound_length)
    return +1;

  if(((struct coll_cand_dat *)a)->rank < ((struct coll_cand_dat *)b)->rank)
    return -1;

  if(((struct coll_cand_dat *)a)->rank > ((struct coll_cand_dat *)b)->rank)
    return +1;

  return 0;
}

/*! \brief Comparison function for coll_cand_dat objects.
 *
 *  Sorting kernel comparing elements (most important first):
 *  rank and len.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b, except for len, where -1 if a > b
 */
int subfind_compare_coll_candidates_rank(const void *a, const void *b)
{
  if(((struct coll_cand_dat *)a)->rank < ((struct coll_cand_dat *)b)->rank)
    return -1;

  if(((struct coll_cand_dat *)a)->rank > ((struct coll_cand_dat *)b)->rank)
    return +1;

  if(((struct coll_cand_dat *)a)->len > ((struct coll_cand_dat *)b)->len)
    return -1;

  if(((struct coll_cand_dat *)a)->len < ((struct coll_cand_dat *)b)->len)
    return +1;

  return 0;
}

/*! \brief Comparison function for variables of MyIDType.
 *
 *  Sorting kernel.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int subfind_fof_compare_ID(const void *a, const void *b)
{
  if(*((MyIDType *)a) < *((MyIDType *)b))
    return -1;

  if(*((MyIDType *)a) > *((MyIDType *)b))
    return +1;

  return 0;
}

#endif /* #ifdef SUBFIND */
