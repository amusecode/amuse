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
 * \file        src/mesh/criterion_derefinement.c
 * \date        05/2018
 * \brief       Criteria for the de-refinement of a cell.
 * \details     Routines which are checking whether a cell should be
 *              de-refined.
 *              contains functions:
 *                int derefine_should_this_cell_be_merged(int i, int flag)
 *                static int derefine_criterion_default(int i)
 *                static int derefine_criterion_jeans_ref(int i)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 04.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#if defined(REFINEMENT_MERGE_CELLS) && !defined(ONEDIMS)
static int derefine_criterion_jeans_ref(int i);
static int derefine_criterion_default(int i);
static int jeans_derefinement_criteria(int i);

/*! \brief Should this cell be dissolved?
 *
 *  This function signals whether a cell should be dissolved. This needs to be
 *  adjusted according to the needs of the simulation in question. One may also
 *  set the SphP[].Flag variable beforehand, these cells will also be
 *  dissolved.
 *
 *  \param[in] i Index of cell in P and SphP arrays.
 *  \param[in] flag If this is nonzero, flag is returned.
 *
 *  \return Flag if this cell should be dissolved.
 */
int derefine_should_this_cell_be_merged(int i, int flag)
{
#ifdef REFINEMENT_HIGH_RES_GAS
  if(SphP[i].AllowRefinement == 0)
    return 0;
#endif /* #ifdef REFINEMENT_HIGH_RES_GAS */

#ifdef NODEREFINE_BACKGROUND_GRID
  /* Keep in mind that this is used in cosmological zoom simulations.
   * I.e. this enforces no derefinement for cells in low-res region, while not
   * affecting the high-res region.
   */
  if(SphP[i].Volume > 0.1 * All.MeanVolume)
    return 0;
#endif /* #ifdef NODEREFINE_BACKGROUND_GRID */

#if defined(REFINEMENT_VOLUME_LIMIT)
  double maxvolume = All.MaxVolume;
  double minvolume = All.MinVolume;

  if(SphP[i].Volume > 0.5 * maxvolume)
    return 0;

  if(SphP[i].Volume < 0.5 * minvolume)
    return 1;

  if(All.MaxVolumeDiff > 0 && SphP[i].Volume > 0.3 * All.MaxVolumeDiff * SphP[i].MinNgbVolume)
    return 0;
#endif /* #if defined(REFINEMENT_VOLUME_LIMIT) */

  if(flag)
    return flag;

  switch(All.DerefinementCriterion)
    {
      case 0:
        return 0;
        break;

      case 1:
        return derefine_criterion_default(i);
        break;

      case 2:
        return derefine_criterion_jeans_ref(i);
        break;

      default:
        terminate("invalid derefinement criterion specified");
        break;
    }

  return 0;
}

/*
 * static functions; i.e. functions that are only called within this file
 */

/*! \brief Default de-refinement criterion.
 *
 *  Checks if cell is within a factor of 2 of the target gas mass.
 *
 *  \param[in] i Index of cell in P and SphP arrays.
 *
 *  \return Flag if this cell should be dissolved.
 */
static int derefine_criterion_default(int i)
{
#if defined(REFINEMENT_SPLIT_CELLS) && defined(REFINEMENT_MERGE_CELLS)

  if(P[i].Mass < 0.5 * All.TargetGasMass)
    return 1;
#endif /* #if defined(REFINEMENT_SPLIT_CELLS) && defined(REFINEMENT_MERGE_CELLS) */

  return 0;
}

/*! \brief Wrapper for Jeans de-refinement criterion.
 *
 *  \param[in] i Index of cell in P and SphP arrays.
 *
 *  \return Flag if this cell should be dissolved.
 */
static int derefine_criterion_jeans_ref(int i)
{
#ifdef JEANS_REFINEMENT
  return jeans_derefinement_criteria(i);
#endif /* #ifdef JEANS_REFINEMENT */
  return 0;
}

/*! \brief De-refinement criterion according to Jeans stability of a cell.
 *
 *  The cell can only be de-refined if the Jeans length is resolved by
 *  1.5 * JEANS_REFINEMENT cells. Otherwise, no de-refinement is possible even
 *  if the cell has a low mass.
 *
 *  \param[in] i Index of cell in P and SphP arrays.
 *
 *  \return Flag if this cell should be dissolved.
 */
static int jeans_derefinement_criteria(int i)
{
  if(P[i].Mass < 0.5 * All.TargetGasMass)
    return 1;

#ifdef JEANS_REFINEMENT
  double jeans_number, jeans_length, sound_speed, dx;
  sound_speed  = sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density);
  jeans_length = sqrt(M_PI / All.G / SphP[i].Density) * sound_speed;
  dx           = 2.0 * get_cell_radius(i);
  jeans_number = jeans_length / dx;

  if(jeans_number > 1.5 * JEANS_REFINEMENT && P[i].Mass < 0.5 * All.TargetGasMass)
    return 1;
#endif /* #ifdef JEANS_REFINEMENT */
  return 0;
}

#endif /* #if defined(REFINEMENT_MERGE_CELLS) && !defined(ONEDIMS) */
