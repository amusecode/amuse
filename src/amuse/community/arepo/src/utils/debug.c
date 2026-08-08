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
 * \file        src/utils/debug.c
 * \date        05/2018
 * \brief       Print relevant information about a particle / face for
 *              debugging.
 * \details     The functions contained in this file are mostly called when a
 *              condition, that causes the abort of the run, is met. In that
 *              case, the information about the state of the particle / face
 *              which triggered that condition is printed to the standard
 *              output.
 *              contains functions:
 *                void print_particle_info(int i)
 *                void print_particle_info_from_ID(MyIDType ID)
 *                void print_state_info(struct state *st)
 *                void print_state_face_info(struct state_face *st)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 03.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/*! \brief Prints particle / cell information to standard output.
 *
 *  \param[in] i Index of particle / cell.
 *
 *  \return void
 */
void print_particle_info(int i)
{
  printf("Task=%d, ID=%llu, Type=%d, TimeBinGrav=%d, TimeBinHydro=%d, Mass=%g, pos=%g|%g|%g, vel=%g|%g|%g\n", ThisTask,
         (unsigned long long)P[i].ID, P[i].Type, P[i].TimeBinGrav, P[i].TimeBinHydro, P[i].Mass, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],
         P[i].Vel[0], P[i].Vel[1], P[i].Vel[2]);
#ifdef PMGRID
  printf("GravAccel=%g|%g|%g, GravPM=%g|%g|%g, Soft=%g, SoftType=%d, OldAcc=%g\n", P[i].GravAccel[0], P[i].GravAccel[1],
         P[i].GravAccel[2], P[i].GravPM[0], P[i].GravPM[1], P[i].GravPM[2], All.ForceSoftening[P[i].SofteningType], P[i].SofteningType,
         P[i].OldAcc);
#else  /* #ifdef PMGRID */
  printf("GravAccel=%g|%g|%g, Soft=%g, SoftType=%d, OldAcc=%g\n", P[i].GravAccel[0], P[i].GravAccel[1], P[i].GravAccel[2],
         All.ForceSoftening[P[i].SofteningType], P[i].SofteningType, P[i].OldAcc);
#endif /* #ifdef PMGRID #else */

  if(P[i].Type == 0)
    {
      printf("Vol=%g, rad=%g, rho=%g, p=%g,u=%g, velVertex=%g|%g|%g, csnd=%g\n", SphP[i].Volume, get_cell_radius(i), SphP[i].Density,
             SphP[i].Pressure, SphP[i].Utherm, SphP[i].VelVertex[0], SphP[i].VelVertex[1], SphP[i].VelVertex[2], get_sound_speed(i));
      printf("Center-Pos=%g|%g|%g\n", SphP[i].Center[0] - P[i].Pos[0], SphP[i].Center[1] - P[i].Pos[1],
             SphP[i].Center[2] - P[i].Pos[2]);
#ifndef MHD
      printf("Mom=%g|%g|%g, Energy=%g, EInt=%g, EKin=%g\n", SphP[i].Momentum[0], SphP[i].Momentum[1], SphP[i].Momentum[2],
             SphP[i].Energy, SphP[i].Utherm * P[i].Mass,
             0.5 * P[i].Mass *
                 ((SphP[i].Momentum[0] / P[i].Mass) * (SphP[i].Momentum[0] / P[i].Mass) +
                  (SphP[i].Momentum[1] / P[i].Mass) * (SphP[i].Momentum[1] / P[i].Mass) +
                  (SphP[i].Momentum[2] / P[i].Mass) * (SphP[i].Momentum[2] / P[i].Mass)));
#else  /* #ifndef MHD */
      printf("Mom=%g|%g|%g, Energy=%g, EInt=%g, EKin=%g, EB=%g\n", SphP[i].Momentum[0], SphP[i].Momentum[1], SphP[i].Momentum[2],
             SphP[i].Energy, SphP[i].Utherm * P[i].Mass,
             0.5 * P[i].Mass *
                 ((SphP[i].Momentum[0] / P[i].Mass) * (SphP[i].Momentum[0] / P[i].Mass) +
                  (SphP[i].Momentum[1] / P[i].Mass) * (SphP[i].Momentum[1] / P[i].Mass) +
                  (SphP[i].Momentum[2] / P[i].Mass) * (SphP[i].Momentum[2] / P[i].Mass)),
             0.5 * SphP[i].Volume * (SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]));
#endif /* #ifndef MHD #else */

#ifdef MHD
      double err = pow(SphP[i].Volume, 1. / 3.) * fabs(SphP[i].DivB) /
                   sqrt(SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]);
      printf("B=%g|%g|%g, divb=%g, err=%g\n", SphP[i].B[0], SphP[i].B[1], SphP[i].B[2], SphP[i].DivB, err);
#endif /* #ifdef MHD */

#ifdef TREE_BASED_TIMESTEPS
      printf("ID=%llu SphP[p].CurrentMaxTiStep=%g\n", (unsigned long long)P[i].ID, SphP[i].CurrentMaxTiStep);
#endif /* #ifdef TREE_BASED_TIMESTEPS */
    }
}

/*! \brief Prints particle / cell information of the cell with a specific ID.
 *
 *  \param[in] ID particle / cell ID.
 *
 *  \return void
 */
void print_particle_info_from_ID(MyIDType ID)
{
  int i;
  for(i = 0; i < NumPart; i++)
    if(P[i].ID == ID)
      print_particle_info(i);
}

/*! \brief Prints information of the left or right state of a face to standard
 *         output.
 *
 *  \param[in] st Structure containing the left or right state of a face.
 *
 *  \return void
 */
void print_state_info(struct state *st)
{
  printf("Task=%d, ID=%llu rho=%g, p=%g, vel=%g|%g|%g, velVertex=%g|%g|%g\n", ThisTask, (unsigned long long)st->ID, st->rho, st->press,
         st->velx, st->vely, st->velz, st->velVertex[0], st->velVertex[1], st->velVertex[2]);
  printf("dx=%g, dy=%g, dz=%g, dt_half=%g\n", st->dx, st->dy, st->dz, st->dt_half);
  printf("timeBin=%d, volume=%g, activearea=%g, surfacearea=%g, csnd=%g\n", st->timeBin, st->volume, st->activearea, st->surfacearea,
         st->csnd);
#ifdef MHD
  printf("B=%g|%g|%g\n", st->Bx, st->By, st->Bz);
#endif /* #ifdef MHD */
}

/*! \brief Prints information of the state the of a face as determined by
 *         the Riemman solver to standard output.
 *
 *  \param[in] st Structure containing the state of a face after the solution
 *             of the Riemann problem.
 *
 *  \return void
 */
void print_state_face_info(struct state_face *st)
{
  printf("rho=%g, p=%g, vel=%g|%g|%g\n", st->rho, st->press, st->velx, st->vely, st->velz);
}
