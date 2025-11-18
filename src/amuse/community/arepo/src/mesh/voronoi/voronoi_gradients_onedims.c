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
 * \file        src/mesh/voronoi/voronoi_gradients.c
 * \date        05/2018
 * \brief       Algorithms to calculate the gradients in 1d simulations.
 * \details     contains functions:
 *                double getValue(int i, int k)
 *                void calculate_gradients(void)
 *                void compute_divvel()
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 23.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include "../../main/allvars.h"
#include "../../main/proto.h"

#if defined(ONEDIMS)

#ifdef OUTPUT_DIVVEL
static void compute_divvel();
#endif /* #ifdef OUTPUT_DIVVEL */

/*! \brief Gets a value of a quantity.
 *
 *  \param[i] Index of cell in P and SphP array.
 *  \param[i] Index in grad_elements array (determines which quantity).
 *
 *  \return value
 */
double getValue(int i, int k)
{
  if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) ||
     (grad_elements[k].type == GRADIENT_TYPE_VELZ))
    return *(MyFloat *)(((char *)(&P[i])) + grad_elements[k].offset);
  else
    return *(MyFloat *)(((char *)(&SphP[i])) + grad_elements[k].offset);
}

/*! \brief Calculates gradients in a 1d simulation.
 *
 *  \return void
 */
void calculate_gradients(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  printf("Calculating 1D gradients...\n");

  int idx, i, k;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(k = 0; k < N_Grad; k++)
        {
          double Value = getValue(i, k);
          double Pos   = P[i].Pos[0];

#if defined(ONEDIMS_SPHERICAL) || defined(REFLECTIVE_X)
          if(i == 0 || i == NumGas - 1)
            {
              MySingle *data = (MySingle *)(((char *)(&(SphP[i].Grad))) + grad_elements[k].offset_grad);
              memset(data, 0, 3 * sizeof(MySingle));
              continue;
            }
#endif /* #if defined (ONEDIMS_SPHERICAL) || defined (REFLECTIVE_X) */
          /* if we get here, we have periodic boundary conditions or are not at the boundaries */
          double ValueL, ValueR;

          if(i == 0)
            ValueL = getValue(NumGas - 1, k);
          else
            ValueL = getValue(i - 1, k);

          if(i == NumGas - 1)
            ValueR = getValue(0, k);
          else
            ValueR = getValue(i + 1, k);

          double PosL = Mesh.DP[i - 1].x;
          double PosR = Mesh.DP[i + 1].x;

          double grad = (ValueL - ValueR) / (PosL - PosR);

          MySingle *data = (MySingle *)(((char *)(&(SphP[i].Grad))) + grad_elements[k].offset_grad);
          data[0]        = grad;
          data[1]        = 0;
          data[2]        = 0;

          double ValueMin = dmin(ValueL, ValueR);
          double ValueMax = dmax(ValueL, ValueR);

          if(Value + grad * (PosL - Pos) < ValueMin)
            {
              if(ValueMin < Value)
                grad = (ValueMin - Value) / (PosL - Pos);
              else
                grad = 0.;
            }

          if(Value + grad * (PosL - Pos) > ValueMax)
            {
              if(ValueMax > Value)
                grad = (ValueMax - Value) / (PosL - Pos);
              else
                grad = 0.;
            }

          if(Value + grad * (PosR - Pos) < ValueMin)
            {
              if(ValueMin < Value)
                grad = (ValueMin - Value) / (PosR - Pos);
              else
                grad = 0.;
            }

          if(Value + grad * (PosR - Pos) > ValueMax)
            {
              if(ValueMax > Value)
                grad = (ValueMax - Value) / (PosR - Pos);
              else
                grad = 0.;
            }

          data[0] = grad;
        }
    }

#ifdef OUTPUT_DIVVEL
  compute_divvel();
#endif /* #ifdef OUTPUT_DIVVEL */

  CPU_Step[CPU_GRADIENTS] += measure_time();
}

#ifdef OUTPUT_DIVVEL
/*! \brief Calculates velocity divergence in 1d simulation.
 *
 *  Using Gauss' theorem.
 *
 *  \return void
 */
void compute_divvel()
{
  face *VF = Mesh.VF;
  double VelxL, VelxR;

  int idx, i;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(i == 0)
        {
#if defined(ONEDIMS_SPHERICAL) || defined(REFLECTIVE_X)
          VelxL = P[i].Vel[0];
#else  /* #if defined (ONEDIMS_SPHERICAL) || defined (REFLECTIVE_X) */
          VelxL = P[NumGas - 1].Vel[0];
#endif /* #if defined (ONEDIMS_SPHERICAL) || defined (REFLECTIVE_X) #else */
        }
      else
        VelxL = P[i - 1].Vel[0];

      if(i == NumGas - 1)
        {
#if defined(ONEDIMS_SPHERICAL) || defined(REFLECTIVE_X)
          VelxR = P[i].Vel[0];
#else  /* #if defined (ONEDIMS_SPHERICAL) || defined (REFLECTIVE_X) */
          VelxR = P[0].Vel[0];
#endif /* #if defined (ONEDIMS_SPHERICAL) || defined (REFLECTIVE_X) #else */
        }
      else
        VelxR = P[i + 1].Vel[0];

      SphP[i].DivVel = 0.5 * (VF[i].area * VelxR - VF[i - 1].area * VelxL) / SphP[i].Volume;
    }
}
#endif /* #ifdef OUTPUT_DIVVEL */

#endif /* #if defined(ONEDIMS) */
