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
 * \file        src/mesh/voronoi/voronoi_1d.c
 * \date        05/2018
 * \brief       Routines to build a 1d Voronoi mesh
 * \details     Note that some of these routines have the same name as the ones
 *              in voronoi_2d.c and voronoi_3d.c and just replace them in case
 *              the Config-option ONEDIMS is active. This is also the reason
 *              why some of these functions are empty but nonetheless have to
 *              exist in this file.
 *              contains functions:
 *                void write_voronoi_mesh(tessellation * T, char *fname,
 *                  int writeTask, int lastTask)
 *                void initialize_and_create_first_tetra(tessellation * T)
 *                void compute_circumcircles(tessellation * T)
 *                void set_integers_for_point(tessellation * T, int pp)
 *                int insert_point(tessellation * T, int pp, int ttstart)
 *                int voronoi_ghost_search(tessellation * T)
 *                int count_undecided_tetras(tessellation * T)
 *                int voronoi_ghost_search_alternative(tessellation * T)
 *                void compute_voronoi_faces_and_volumes(void)
 *                void voronoi_1D_order(void)
 *                int voronoi_1D_compare_key(const void *a, const void *b)
 *                void voronoi_1D_reorder_gas(void)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 21.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../main/allvars.h"
#include "../../main/proto.h"

#include "voronoi.h"

#if defined(ONEDIMS) && !defined(ONEDIMS_SPHERICAL) /* will only be compiled in 1D case */

/*! \brief Output of Voroioi mesh to file.
 *
 *  Not supported for 1d.
 *
 *  \return void
 */
void write_voronoi_mesh(tessellation *T, char *fname, int writeTask, int lastTask)
{
  terminate("write_voronoi_mesh not supported in 1d case!");
}

/*! \brief Initialises 1d tessellation and create all-enclosing segment.
 *
 *  \param[out] T Pointer to tessllation structure which is set and its arrays
 *              are allocated in this routine.
 *
 *  \return void
 */
void initialize_and_create_first_tetra(tessellation *T)
{
  char msg[200];

  if(NTask > 1)
    {
      mpi_printf("1D code works only for 1 CPU\n");
      endrun();
    }

  T->MaxNdp = NumGas + 4;
  T->MaxNdt = 4 + T->MaxNdp * 2;
  T->MaxNvf = T->MaxNdt;

  if(NumGas == 0)
    {
      sprintf(msg, "NumGas=%d on Task=%d, but need at least one particle!\n", NumGas, ThisTask);
      terminate(msg);
    }

  T->Ndp = 0;
  T->Nvf = 0;
  T->Ndt = 0;

  T->VF = mymalloc_movable(&T->VF, "VF", T->MaxNvf * sizeof(face));

  T->DP = mymalloc_movable(&T->DP, "DP", (T->MaxNdp + 5) * sizeof(point));
  T->DP += 5;

  T->DT = mymalloc_movable(&T->DT, "DT", T->MaxNdt * sizeof(tetra));
}

/*! \brief Computes circumcircles in 1d.
 *
 *  Not necessary in 1d. However, this function has to exist for the 1d code
 *  to work.
 *
 *  \param[in] T Pointer to tessllation structure.
 *
 *  \return void
 */
void compute_circumcircles(tessellation *T) {}

/*! \brief Empty funciton in 1d case.
 *
 *  Not necessary in 1d. However, this function has to exist for the 1d code
 *  to work.
 *
 * \return void
 */
void set_integers_for_point(tessellation *T, int pp) {}

/*! \brief Empty funciton in 1d case.
 *
 *  Not necessary in 1d. However, this function has to exist for the 1d code
 *  to work.
 *
 * \return 0
 */
int insert_point(tessellation *T, int pp, int ttstart) { return 0; }

/*! \brief Wrapper routine to search for ghost cells for boundary cells.
 *
 *  \param[out] T Pointer to tessellation.
 *
 *  \return 0
 */
int voronoi_ghost_search(tessellation *T) { return voronoi_ghost_search_alternative(T); }

/*! \brief Empty funciton in 1d case.
 *
 *  Not necessary in 1d. However, this function has to exist for the 1d code
 *  to work.
 *
 * \return 0
 */
int count_undecided_tetras(tessellation *T) { return 0; }

/*! \brief Searches for ghost cells in 1d Voronoi mesh.
 *
 *  This routine assumes an x ordered cell array.
 *
 *  \param[out] T pointer to tessellation.
 *
 *  \return 0
 */
int voronoi_ghost_search_alternative(tessellation *T)
{
  double xl, xr;
  int index_l, index_r;

#if defined(REFLECTIVE_X)
  xl      = -P[0].Pos[0];
  index_l = 0;

  xr      = boxSize_X + (boxSize_X - P[NumGas - 1].Pos[0]);
  index_r = NumGas - 1;
#else  /* #if defined(REFLECTIVE_X) */
  xl      = P[NumGas - 1].Pos[0] - boxSize_X;
  index_l = NumGas - 1;

  xr      = P[0].Pos[0] + boxSize_X;
  index_r = 0;
#endif /* #if defined(REFLECTIVE_X) #else */

  point *DP = T->DP;

  DP[-1].x     = xl;
  DP[-1].y     = 0;
  DP[-1].z     = 0;
  DP[-1].task  = ThisTask;
  DP[-1].ID    = P[index_l].ID;
  DP[-1].index = index_l + NumGas; /* this is a mirrored local point */
#if defined(REFLECTIVE_X)
  DP[-1].image_flags = REFL_X_FLAGS;
#if(REFLECTIVE_X == 2)
  DP[-1].image_flags |= OUTFLOW_X;
#endif /* #if (REFLECTIVE_X == 2) */
#endif /* #if defined(REFLECTIVE_X) */
  DP[NumGas].x     = xr;
  DP[NumGas].y     = 0;
  DP[NumGas].z     = 0;
  DP[NumGas].task  = ThisTask;
  DP[NumGas].ID    = P[index_r].ID;
  DP[NumGas].index = index_r + NumGas; /* this is a mirrored local point */
#if defined(REFLECTIVE_X)
  DP[NumGas].image_flags = REFL_X_FLAGS;
#if(REFLECTIVE_X == 2)
  DP[NumGas].image_flags |= OUTFLOW_X;
#endif /* #if (REFLECTIVE_X == 2) */
#endif /* #if defined(REFLECTIVE_X) */
  return 0;
}

/*! \brief Computes faces and volume of cells in 1d Voronoi mesh.
 *
 *  Also computes the center of mass.
 *
 *  \return void
 */
void compute_voronoi_faces_and_volumes(void)
{
  int i;

  tessellation *T = &Mesh;

  T->Nvf    = 0;
  point *DP = T->DP;
  face *VF  = T->VF;

  for(i = -1; i < NumGas; i++)
    {
      VF[T->Nvf].p1 = i;
      VF[T->Nvf].p2 = i + 1;

      VF[T->Nvf].cx = 0.5 * (DP[i].x + DP[i + 1].x);

      VF[T->Nvf].cy   = 0;
      VF[T->Nvf].cz   = 0;
      VF[T->Nvf].area = 1;

      T->Nvf++;
    }

  for(i = 0; i < NumGas; i++)
    {
      SphP[i].Volume    = VF[i + 1].cx - VF[i].cx;
      SphP[i].Center[0] = 0.5 * (VF[i + 1].cx + VF[i].cx);
      SphP[i].Center[1] = 0;
      SphP[i].Center[2] = 0;

      SphP[i].SurfaceArea = 2.;
    }
}

/*! \brief Data for 1d Voronoi mesh.
 */
static struct voronoi_1D_data
{
  double x;
  int index;
} * mp;

static int *Id;

/*! \brief Sort cells by their position and reorder in P and SphP array.
 *
 *  \return void
 */
void voronoi_1D_order(void)
{
  int i;

  mpi_printf("begin 1D order...\n");

  if(NumGas)
    {
      mp = (struct voronoi_1D_data *)mymalloc("mp", sizeof(struct voronoi_1D_data) * NumGas);
      Id = (int *)mymalloc("Id", sizeof(int) * NumGas);

      for(i = 0; i < NumGas; i++)
        {
          mp[i].index = i;
          mp[i].x     = P[i].Pos[0];
        }

      mysort(mp, NumGas, sizeof(struct voronoi_1D_data), voronoi_1D_compare_key);

      for(i = 0; i < NumGas; i++)
        Id[mp[i].index] = i;

      voronoi_1D_reorder_gas();

      myfree(Id);
      myfree(mp);
    }

  mpi_printf("1D order done.\n");
}

/*! \brief Compare x value of voronoi_1D_data objects.
 *
 *  \param[in] a Pointer to first voronoi_1D_data object.
 *  \param[in] b Pointer to second voronoi_1D_data object.
 *
 *  \return (-1,0,1) -1 if a->x < b->x.
 */
int voronoi_1D_compare_key(const void *a, const void *b)
{
  if(((struct voronoi_1D_data *)a)->x < (((struct voronoi_1D_data *)b)->x))
    return -1;

  if(((struct voronoi_1D_data *)a)->x > (((struct voronoi_1D_data *)b)->x))
    return +1;

  return 0;
}

/*! \brief Order the gas cells according to the index given in the ID array.
 *
 *  \return void
 */
void voronoi_1D_reorder_gas(void)
{
  int i;
  struct particle_data Psave, Psource;
  struct sph_particle_data SphPsave, SphPsource;
  int idsource, idsave, dest;

  for(i = 0; i < NumGas; i++)
    {
      if(Id[i] != i)
        {
          Psource    = P[i];
          SphPsource = SphP[i];

          idsource = Id[i];
          dest     = Id[i];

          do
            {
              Psave    = P[dest];
              SphPsave = SphP[dest];
              idsave   = Id[dest];

              P[dest]    = Psource;
              SphP[dest] = SphPsource;
              Id[dest]   = idsource;

              if(dest == i)
                break;

              Psource    = Psave;
              SphPsource = SphPsave;
              idsource   = idsave;

              dest = idsource;
            }
          while(1);
        }
    }
}

#endif /* #if defined (ONEDIMS) && !defined (ONEDIMS_SPHERICAL) */
