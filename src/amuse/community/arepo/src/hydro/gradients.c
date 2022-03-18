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
 * \file        src/gradients.c
 * \date        05/2018
 * \brief       Routines to initialize gradient data.
 * \details     contains functions:
 *                void init_gradients()
 *                void gradient_init(MyFloat * addr, MyFloat * addr_exch,
 *                  MySingle * addr_grad, int type)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 05.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../mesh/voronoi/voronoi.h"

int N_Grad = 0;

struct grad_elements grad_elements[MAXGRADIENTS], *GDensity, *GVelx, *GVely, *GVelz, *GPressure, *GUtherm;

/*! \brief Initializes all gradient fields.
 *
 *  Density, velocity, pressure and if needed magnetic fields and passive
 *  scalars.
 *
 *  \return void
 */
void init_gradients()
{
#if defined(MAXSCALARS)
  int k;
#endif /* #if defined(MAXSCALARS) */

  gradient_init(&SphP[0].Density, &PrimExch[0].Density, SphP[0].Grad.drho, GRADIENT_TYPE_DENSITY);

  gradient_init(&P[0].Vel[0], &PrimExch[0].VelGas[0], SphP[0].Grad.dvel[0], GRADIENT_TYPE_VELX);
  gradient_init(&P[0].Vel[1], &PrimExch[0].VelGas[1], SphP[0].Grad.dvel[1], GRADIENT_TYPE_VELY);
  gradient_init(&P[0].Vel[2], &PrimExch[0].VelGas[2], SphP[0].Grad.dvel[2], GRADIENT_TYPE_VELZ);

  gradient_init(&SphP[0].Pressure, &PrimExch[0].Pressure, SphP[0].Grad.dpress, GRADIENT_TYPE_PRESSURE);

#ifdef MHD
  gradient_init(&SphP[0].B[0], &PrimExch[0].B[0], SphP[0].Grad.dB[0], GRADIENT_TYPE_NORMAL);
  gradient_init(&SphP[0].B[1], &PrimExch[0].B[1], SphP[0].Grad.dB[1], GRADIENT_TYPE_NORMAL);
  gradient_init(&SphP[0].B[2], &PrimExch[0].B[2], SphP[0].Grad.dB[2], GRADIENT_TYPE_NORMAL);
#endif /* #ifdef MHD */

#ifdef MAXSCALARS
  MyFloat *addr;

  for(k = 0; k < N_Scalar; k++)
    {
      addr = (MyFloat *)(((char *)(&SphP[0])) + scalar_elements[k].offset);
      gradient_init(addr, &PrimExch[0].Scalars[k], SphP[0].Grad.dscalars[k], GRADIENT_TYPE_NORMAL);
    }
#endif /* #ifdef MAXSCALARS */

  mpi_printf("INIT: %d/%d Gradients used.\n", N_Grad, MAXGRADIENTS);
}

/*! \brief Initialize a gradient field.
 *
 *  Each time this initialization routine is called, the global variable
 *  NGrad is incremented by 1.
 *
 *  \param[in] addr Pointer to element in SphP[0] struct (for Vel in P[0])
 *  \param[in] addr_exch Pointer to element in PrimExch[0] struct
 *  \param[in] addr_grad Pointer to element in SphP[0].Grad struct
 *  \param[in] type Type of gradient
 *
 *  \return void
 */
void gradient_init(MyFloat *addr, MyFloat *addr_exch, MySingle *addr_grad, int type)
{
  if(N_Grad == MAXGRADIENTS)
    {
      mpi_printf("Failed to register gradient, maximum of %d already reached\n", MAXGRADIENTS);
      terminate("MAXGRADIENTS reached");
    }

  grad_elements[N_Grad].type = type;

  if((type == GRADIENT_TYPE_VELX) || (type == GRADIENT_TYPE_VELY) || (type == GRADIENT_TYPE_VELZ))
    {
      /* basic structure is P */
      grad_elements[N_Grad].offset = ((char *)addr) - ((char *)&P[0]);
    }
  else
    {
      /* basic structure is SphP */
      grad_elements[N_Grad].offset = ((char *)addr) - ((char *)&SphP[0]);
    }

  grad_elements[N_Grad].offset_exch = ((char *)addr_exch) - ((char *)&PrimExch[0]);
  grad_elements[N_Grad].offset_grad = ((char *)addr_grad) - ((char *)&(SphP[0].Grad));

  switch(type)
    {
      case GRADIENT_TYPE_VELX:
        GVelx = &grad_elements[N_Grad];
        break;
      case GRADIENT_TYPE_VELY:
        GVely = &grad_elements[N_Grad];
        break;
      case GRADIENT_TYPE_VELZ:
        GVelz = &grad_elements[N_Grad];
        break;
      case GRADIENT_TYPE_DENSITY:
        GDensity = &grad_elements[N_Grad];
        break;
      case GRADIENT_TYPE_PRESSURE:
        GPressure = &grad_elements[N_Grad];
        break;
      case GRADIENT_TYPE_UTHERM:
        GUtherm = &grad_elements[N_Grad];
        break;
      default:
        break;
    }

  N_Grad++;
}
