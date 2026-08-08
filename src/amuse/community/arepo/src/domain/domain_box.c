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
 * \file        src/domain_box.c
 * \date        05/2018
 * \brief       Routines that determine domain box and do periodic wrapping.
 * \details     contains files:
 *                void domain_findExtent(void)
 *                void do_box_wrapping(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 05.05.2018 Prepared file for public release -- Rainer Weinberger
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

/*! \brief Move the coordinate in pos by the global displacement vector
 *
 *  \param[in] pos coordinate vector (3 entries).
 *  \param[in] mode displacement mode, either DISPLACE_POSITION_FORWARD or DISPLACE_POSITION_BACKWARD
 *
 *  \return void
 */
void domain_displacePosition(MyDouble *pos, enum domain_displace_mode mode)
{
  if(mode == DISPLACE_POSITION_FORWARD)
    {
      double xtmp, ytmp, ztmp;
      pos[0] = WRAP_X(pos[0] + All.GlobalDisplacementVector[0]);
      pos[1] = WRAP_Y(pos[1] + All.GlobalDisplacementVector[1]);
      pos[2] = WRAP_Z(pos[2] + All.GlobalDisplacementVector[2]);
    }
  else if(mode == DISPLACE_POSITION_BACKWARD)
    {
      double xtmp, ytmp, ztmp;
      pos[0] = WRAP_X(pos[0] - All.GlobalDisplacementVector[0]);
      pos[1] = WRAP_Y(pos[1] - All.GlobalDisplacementVector[1]);
      pos[2] = WRAP_Z(pos[2] - All.GlobalDisplacementVector[2]);
    }
  else
    terminate("Unkown mode %d.", mode);
}

/*! \brief Move the coordinate for all positions by the global displacement vector
 *
 *  \param[in] mode displacement mode, either DISPLACE_POSITION_FORWARD or DISPLACE_POSITION_BACKWARD
 *
 *  \return void
 */
static void domain_displacePositions(enum domain_displace_mode mode)
{
  for(int i = 0; i < NumPart; i++)
    {
      if(P[i].ID == 0 && P[i].Mass == 0) /* derefined */
        continue;

      domain_displacePosition(P[i].Pos, mode);

      if(i < NumGas)
        domain_displacePosition(SphP[i].Center, mode);
    }

#ifdef PLACEHIGHRESREGION
  domain_displacePosition(All.Xmintot[1], mode);
  domain_displacePosition(All.Xmaxtot[1], mode);
  domain_displacePosition(All.Corner[1], mode);
  domain_displacePosition(All.UpperCorner[1], mode);
#endif
}

/*! \brief Finds the extent of the global domain grid.
 *
 *  The minimum extent is the box size.
 *
 *  \return void
 */
void domain_findExtent(void)
{
  int i, j;
  double len, xmin[3], xmax[3], xmin_glob[3], xmax_glob[3];

  /* determine local extension */
  for(j = 0; j < 3; j++)
    {
      /* preset to simulation box */
      xmin[j] = 0;
      xmax[j] = boxSize;
    }
    // Take care of stretched box
#ifdef LONG_X
  xmax[0] = boxSize_X;
#endif /* #ifdef LONG_X */
#ifdef LONG_Y
  xmax[1] = boxSize_Y;
#endif /* #ifdef LONG_Y */
#ifdef LONG_Z
  xmax[2] = boxSize_Z;
#endif /* #ifdef LONG_Z */

  for(i = 0; i < NumPart; i++)
    {
#ifdef ADDBACKGROUNDGRID
      if(P[i].Type != 0)
        continue;
#endif /* #ifdef ADDBACKGROUNDGRID */
      for(j = 0; j < 3; j++)
        {
          if(xmin[j] > P[i].Pos[j])
            xmin[j] = P[i].Pos[j];

          if(xmax[j] < P[i].Pos[j])
            xmax[j] = P[i].Pos[j];
        }
    }

  MPI_Allreduce(xmin, xmin_glob, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(xmax, xmax_glob, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

#ifdef ADDBACKGROUNDGRID
  for(j = 0; j < 3; j++)
    if(xmax_glob[j] < All.BoxSize)
      xmax_glob[j] = All.BoxSize;

  for(j = 0; j < 3; j++)
    if(xmin_glob[j] > 0)
      xmin_glob[j] = 0;
#endif /* #ifdef ADDBACKGROUNDGRID */

  len = 0;
  for(j = 0; j < 3; j++)
    if(xmax_glob[j] - xmin_glob[j] > len)
      len = xmax_glob[j] - xmin_glob[j];

#if defined(GRAVITY_NOT_PERIODIC) && !defined(ADDBACKGROUNDGRID)
  len *= 1.2; /* enlarge box a bit to avoid triggering of an out of box recovery */
#else         /* #if defined(GRAVITY_NOT_PERIODIC) && !defined(ADDBACKGROUNDGRID) */
  len *= 1.00001;
#endif        /* #if defined(GRAVITY_NOT_PERIODIC) && !defined(ADDBACKGROUNDGRID) #else */

#if defined(DO_NOT_RANDOMIZE_DOMAINCENTER) || !defined(GRAVITY_NOT_PERIODIC) || defined(ONEDIMS) || defined(TWODIMS)
  for(j = 0; j < 3; j++)
    {
      DomainCenter[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]);
      DomainCorner[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]) - 0.5 * len;
    }
#else  /* #if defined(DO_NOT_RANDOMIZE_DOMAINCENTER) || !defined(GRAVITY_NOT_PERIODIC) || defined(ONEDIMS) || defined(TWODIMS) */
  for(j = 0; j < 3; j++)
    {
      DomainCenter[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]);
      DomainCenter[j] += (2. * get_random_number() - 1.) * 0.5 * len;
    }

  MPI_Bcast(DomainCenter, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  len *= 2;

  for(j = 0; j < 3; j++)
    DomainCorner[j] = DomainCenter[j] - 0.5 * len;
#endif /* #if defined(DO_NOT_RANDOMIZE_DOMAINCENTER) || !defined(GRAVITY_NOT_PERIODIC) || defined(ONEDIMS) || defined(TWODIMS) #else \
        */

  DomainLen = len;

  DomainInverseLen = 1.0 / DomainLen;
  DomainFac        = 1.0 / len * (((peanokey)1) << (BITS_PER_DIMENSION));
  DomainBigFac     = (DomainLen / (((long long)1) << 52));
}

/*! \brief Makes sure all particles are within box.
 *
 *  This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 *
 *  \return void
 */
void do_box_wrapping(void)
{
  int j;
  double boxsize[3];

#ifdef ADDBACKGROUNDGRID
  return;
#endif /* #ifdef ADDBACKGROUNDGRID */

  for(j = 0; j < 3; j++)
    boxsize[j] = All.BoxSize;

#ifdef LONG_X
  boxsize[0] *= LONG_X;
#endif /* #ifdef LONG_X */
#ifdef LONG_Y
  boxsize[1] *= LONG_Y;
#endif /* #ifdef LONG_Y */
#ifdef LONG_Z
  boxsize[2] *= LONG_Z;
#endif /* #ifdef LONG_Z */

#if !defined(GRAVITY_NOT_PERIODIC) && !defined(DO_NOT_RANDOMIZE_DOMAINCENTER) && defined(SELFGRAVITY) && (NUMDIMS > 2)
  domain_displacePositions(DISPLACE_POSITION_BACKWARD);

  if(ThisTask == 0)
    {
      double prefac = 1.;
#ifdef PLACEHIGHRESREGION
      prefac = 0.5;
#endif
      for(j = 0; j < 3; j++)
        All.GlobalDisplacementVector[j] = (get_random_number() - 0.5) * boxsize[j] * prefac;
    }

  mpi_printf("DOMAIN: New global displacement vector: %g, %g, %g\n", All.GlobalDisplacementVector[0], All.GlobalDisplacementVector[1],
             All.GlobalDisplacementVector[2]);
  MPI_Bcast(All.GlobalDisplacementVector, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  domain_displacePositions(DISPLACE_POSITION_FORWARD);
#endif /* #if !defined(GRAVITY_NOT_PERIODIC) && !defined(DO_NOT_RANDOMIZE_DOMAINCENTER) && defined(SELFGRAVITY) && (NUMDIMS > 2) */

  int i;
  for(i = 0; i < NumPart; i++)
    {
      if(i < NumGas)
        trans_table[i].wrapped = 0;

#if defined(GRAVITY_NOT_PERIODIC)
      if(P[i].Type != 0)
        continue;
#endif /* #if defined(GRAVITY_NOT_PERIODIC) */

#if !defined(REFLECTIVE_X)
      while(P[i].Pos[0] < 0)
        {
          P[i].Pos[0] += boxsize[0];
          if(i < NumGas)
            trans_table[i].wrapped |= 1;
        }

      while(P[i].Pos[0] >= boxsize[0])
        {
          P[i].Pos[0] -= boxsize[0];
          if(i < NumGas)
            trans_table[i].wrapped |= 2;
        }

#else  /* #if !defined(REFLECTIVE_X) */
      if(P[i].Pos[0] < 0 || P[i].Pos[0] >= boxsize[0])
        {
          char buf[1000];

          sprintf(buf, "i=%d ID=%d type=%d moved out of box. x=%g", i, P[i].ID, P[i].Type, P[i].Pos[0]);
          terminate(buf);
        }
#endif /* #if !defined(REFLECTIVE_X) #else */

#if !defined(REFLECTIVE_Y)
      while(P[i].Pos[1] < 0)
        {
          P[i].Pos[1] += boxsize[1];
          if(i < NumGas)
            trans_table[i].wrapped |= 4;
        }

      while(P[i].Pos[1] >= boxsize[1])
        {
          P[i].Pos[1] -= boxsize[1];
          if(i < NumGas)
            trans_table[i].wrapped |= 8;
        }

#else  /* #if !defined(REFLECTIVE_Y) */
      if(P[i].Pos[1] < 0 || P[i].Pos[1] >= boxsize[1])
        {
          char buf[1000];

          sprintf(buf, "i=%d ID=%d type=%d moved out of box. y=%g", i, P[i].ID, P[i].Type, P[i].Pos[1]);
          terminate(buf);
        }
#endif /* #if !defined(REFLECTIVE_Y) #else */

#if !defined(REFLECTIVE_Z)
      while(P[i].Pos[2] < 0)
        {
          P[i].Pos[2] += boxsize[2];
          if(i < NumGas)
            trans_table[i].wrapped |= 16;
        }

      while(P[i].Pos[2] >= boxsize[2])
        {
          P[i].Pos[2] -= boxsize[2];
          if(i < NumGas)
            trans_table[i].wrapped |= 32;
        }

#else  /* #if !defined(REFLECTIVE_Z) */
      if(P[i].Pos[2] < 0 || P[i].Pos[2] >= boxsize[2])
        {
          char buf[1000];

          sprintf(buf, "i=%d ID=%d type=%d moved out of box. z=%g", i, P[i].ID, P[i].Type, P[i].Pos[2]);
          terminate(buf);
        }
#endif /* #if !defined(REFLECTIVE_Z) #else */
    }
}
