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
 * \file        src/gravity/gravtree.c
 * \date        05/2018
 * \brief       Special gravity routines for external forces.
 * \details     contains functions:
 *                void gravity_external(void)
 *                static void gravity_external_get_force( double pos[3],
 *                  int type, MyIDType ID, double acc[3], double *pot, int
 *                  *flag_set )
 *                void gravity_monopole_1d_spherical()
 *                double enclosed_mass(double R)
 *                void calc_exact_gravity_for_particle_type(void)
 *                void special_particle_create_list()
 *                void special_particle_update_list()
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 05.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

#ifdef EXTERNALGRAVITY
static void gravity_external_get_force(double pos[3], int type, MyIDType ID, double acc[3], double *pot, int *flag_set);

/*! \brief Main routine to add contribution of external gravitational potential
 *  to accelerations.
 *
 *  Function is called in gravity() (in accel.c). Function also evaluates
 *  the gradient of the accelerations which is needed for the timestep
 *  criterion due to the external potential.
 *
 *  \return void
 */
void gravity_external(void)
{
  mpi_printf("EXTERNALGRAVITY: execute\n");

  TIMER_START(CPU_TREE);

  for(int idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      int i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double *pos;

#ifdef CELL_CENTER_GRAVITY
      if(P[i].Type == 0)
        pos = SphP[i].Center;
      else
#endif /* #ifdef CELL_CENTER_GRAVITY */
        pos = P[i].Pos;

      double acc[3], pot;
      int flag_set = 0;
      gravity_external_get_force(pos, P[i].Type, P[i].ID, acc, &pot, &flag_set);

      if(flag_set)
        {
          for(int k = 0; k < NUMDIMS; k++)
            P[i].GravAccel[k] = acc[k];
          for(int k = NUMDIMS; k < 3; k++)
            P[i].GravAccel[k] = 0;
          P[i].ExtPotential = pot;
        }
      else
        {
          for(int k = 0; k < NUMDIMS; k++)
            P[i].GravAccel[k] += acc[k];
#ifdef EVALPOTENTIAL
          P[i].Potential += pot;
#endif
          P[i].ExtPotential += pot;
        }

      double dx;
      if(P[i].Type == 0)
        dx = 0.1 * get_cell_radius(i);
      else
        dx = 0.1 * All.ForceSoftening[P[i].SofteningType];

      P[i].dGravAccel = 0;
      for(int dim = 0; dim < NUMDIMS; dim++)
        {
          double accL[3], posL[3];
          for(int k = 0; k < 3; k++)
            posL[k] = pos[k];
          posL[dim] -= dx;
          gravity_external_get_force(posL, P[i].Type, P[i].ID, accL, &pot, &flag_set);

          double accR[3], posR[3];
          for(int k = 0; k < 3; k++)
            posR[k] = pos[k];
          posR[dim] += dx;
          gravity_external_get_force(posR, P[i].Type, P[i].ID, accR, &pot, &flag_set);

          for(int k = 0; k < NUMDIMS; k++)
            {
              double dGrav = accR[k] - accL[k];
              P[i].dGravAccel += dGrav * dGrav;
            }
        }
      P[i].dGravAccel = sqrt(P[i].dGravAccel) / (2. * dx);
    }

  TIMER_STOP(CPU_TREE);
}

/*! \brief Calculates the force from the external potential given a position.
 *
 *  \param[in] pos Position at which force is to be evaluated.
 *  \param[in] type (unused)
 *  \param[in] ID (unused)
 *  \param[in, out] acc Acceleration array.
 *  \param[in, out] pot Pointer to potential.
 *  \param[in] flag_set (unused)
 *
 *  \return void
 */
static void gravity_external_get_force(double pos[3], int type, MyIDType ID, double acc[3], double *pot, int *flag_set)
{
  for(int k = 0; k < 3; k++)
    acc[k] = 0;

  *pot = 0;

#ifdef EXTERNALGY
  acc[1] += EXTERNALGY;
  *pot = -(EXTERNALGY)*pos[1];
#endif /* #ifdef EXTERNALGY */

#ifdef STATICISO
  {
    double r, m;
    double dx, dy, dz;

    dx = pos[0] - boxHalf_X;
    dy = pos[1] - boxHalf_Y;
    dz = pos[2] - boxHalf_Z;

    r = sqrt(dx * dx + dy * dy + dz * dz);

    if(r > ISO_R200)
      m = ISO_M200;
    else
      m = ISO_M200 * r / ISO_R200;

#ifdef ISO_FRACTION
    m *= ISO_FRACTION;
#endif /* #ifdef ISO_FRACTION */

    if(r > 0)
      {
        acc[0] += -All.G * m * dx / r / (r * r + ISO_Eps * ISO_Eps);
        acc[1] += -All.G * m * dy / r / (r * r + ISO_Eps * ISO_Eps);
        acc[2] += -All.G * m * dz / r / (r * r + ISO_Eps * ISO_Eps);
      }
  }
#endif /* #ifdef STATICISO */

#ifdef STATICNFW
  {
    double r, m;
    double dx, dy, dz;

    dx = pos[0] - boxHalf_X;
    dy = pos[1] - boxHalf_Y;
    dz = pos[2] - boxHalf_Z;

    r = sqrt(dx * dx + dy * dy + dz * dz);
    m = enclosed_mass(r);
#ifdef NFW_DARKFRACTION
    m *= NFW_DARKFRACTION;
#endif /* #ifdef NFW_DARKFRACTION */
    if(r > 0)
      {
        acc[0] += -All.G * m * dx / (r * r * r);
        acc[1] += -All.G * m * dy / (r * r * r);
        acc[2] += -All.G * m * dz / (r * r * r);
      }
  }
#endif /* #ifdef STATICNFW */

#ifdef STATICHQ
  {
    double r, m, a;
    double dx, dy, dz;

    dx = pos[0] - boxHalf_X;
    dy = pos[1] - boxHalf_Y;
    dz = pos[2] - boxHalf_Z;

    r = sqrt(dx * dx + dy * dy + dz * dz);

    a = pow(All.G * HQ_M200 / (100 * All.Hubble * All.Hubble), 1.0 / 3) / HQ_C * sqrt(2 * (log(1 + HQ_C) - HQ_C / (1 + HQ_C)));

    m = HQ_M200 * pow(r / (r + a), 2);
#ifdef HQ_DARKFRACTION
    m *= HQ_DARKFRACTION;
#endif /* #ifdef HQ_DARKFRACTION */
    if(r > 0)
      {
        acc[0] += -All.G * m * dx / (r * r * r);
        acc[1] += -All.G * m * dy / (r * r * r);
        acc[2] += -All.G * m * dz / (r * r * r);
      }
  }
#endif /* #ifdef STATICHQ */
}
#endif /* #ifdef EXTERNALGRAVITY */

#ifdef ONEDIMS_SPHERICAL
/*! \brief One-dimensional gravity in the spherically symmetric case.
 *
 *  \return void
 */
void gravity_monopole_1d_spherical()
{
  printf("Doing 1D gravity...\n");

  int i;
  double msum = All.CoreMass;

  for(i = 0; i < NumGas; i++)
    {
      double r0;
      if(i > 0)
        r0 = 0.5 * (P[i].Pos[0] + P[i - 1].Pos[0]);
      else
        r0 = All.CoreRadius;
      double dm  = 4. / 3. * M_PI * (SphP[i].Center[0] * SphP[i].Center[0] * SphP[i].Center[0] - r0 * r0 * r0) * SphP[i].Density;
      double rad = SphP[i].Center[0];

      P[i].GravAccel[0] = -(msum + dm) * All.G / (rad * rad);

#ifdef EVALPOTENTIAL
      P[i].Potential = -(msum + dm) * All.G / rad;
#endif /* #ifdef EVALPOTENTIAL */

      msum += P[i].Mass;

      P[i].GravAccel[1] = 0;
      P[i].GravAccel[2] = 0;
    }

  printf("... 1D gravity done.\n");
}
#endif /* #ifdef ONEDIMS_SPHERICAL */

#ifdef STATICNFW
/*! \brief Auxiliary function for static NFW potential.
 *
 *  \param[in] R Radius from center of potential.
 *
 *  \return Enclosed mass (which causes the external potential).
 */
double enclosed_mass(double R)
{
  /* Eps is in units of Rs !!!! */

  if(R > Rs * NFW_C)
    R = Rs * NFW_C;

  return fac * 4 * M_PI * RhoCrit * Dc *
         (-(Rs * Rs * Rs * (1 - NFW_Eps + log(Rs) - 2 * NFW_Eps * log(Rs) + NFW_Eps * NFW_Eps * log(NFW_Eps * Rs))) /
              ((NFW_Eps - 1) * (NFW_Eps - 1)) +
          (Rs * Rs * Rs *
           (Rs - NFW_Eps * Rs - (2 * NFW_Eps - 1) * (R + Rs) * log(R + Rs) + NFW_Eps * NFW_Eps * (R + Rs) * log(R + NFW_Eps * Rs))) /
              ((NFW_Eps - 1) * (NFW_Eps - 1) * (R + Rs)));
}
#endif /* #ifdef STATICNFW */

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
/*! \brief Routine that computes gravitational force by direct summation.
 *
 *  Called by gravity() (in accel.c).
 *
 *  \return void
 */
void calc_exact_gravity_for_particle_type(void)
{
  int i, idx;
#ifdef EXACT_GRAVITY_REACTION
  double *accx, *accy, *accz;
  accx = (double *)mymalloc("accx", All.MaxPartSpecial * sizeof(double));
  accy = (double *)mymalloc("accy", All.MaxPartSpecial * sizeof(double));
  accz = (double *)mymalloc("accz", All.MaxPartSpecial * sizeof(double));
#ifdef EVALPOTENTIAL
  double *pot;
  pot = (double *)mymalloc("pot", All.MaxPartSpecial * sizeof(double));
#endif /* #ifdef EVALPOTENTIAL */
  int n;
  for(n = 0; n < All.MaxPartSpecial; n++)
    {
      accx[n] = accy[n] = accz[n] = 0.0;
#ifdef EVALPOTENTIAL
      pot[n] = 0.0;
#endif /* #ifdef EVALPOTENTIAL */
    }
#endif /* #ifdef EXACT_GRAVITY_REACTION */

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double fac, wp;
      double dx, dy, dz, r, r2;
      double h, h_inv, h3_inv, u;
      int k;

      /* set softening to corresponding particle's softening length */
      h = All.ForceSoftening[All.SofteningTypeOfPartType[EXACT_GRAVITY_FOR_PARTICLE_TYPE]];

      for(k = 0; k < All.MaxPartSpecial; k++)
        {
          if(PartSpecialListGlobal[k].ID == P[i].ID)
            continue;

          dx = P[i].Pos[0] - PartSpecialListGlobal[k].pos[0];
          dy = P[i].Pos[1] - PartSpecialListGlobal[k].pos[1];
          dz = P[i].Pos[2] - PartSpecialListGlobal[k].pos[2];

          r2 = dx * dx + dy * dy + dz * dz;
          r  = sqrt(r2);

          // using spline softening
          if(r >= h)
            {
              fac = 1 / (r2 * r);
              wp  = -1 / r;
            }
          else
            {
              h_inv  = 1.0 / h;
              h3_inv = h_inv * h_inv * h_inv;
              u      = r * h_inv;

              if(u < 0.5)
                {
                  fac = h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
                  wp  = h_inv * (-2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6)));
                }
              else
                {
                  fac = h3_inv *
                        (21.333333333333 - 48.0 * u + 38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
                  wp = h_inv * (-3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u))));
                }
            }

          P[i].GravAccel[0] -= All.G * PartSpecialListGlobal[k].mass * fac * dx;
          P[i].GravAccel[1] -= All.G * PartSpecialListGlobal[k].mass * fac * dy;
          P[i].GravAccel[2] -= All.G * PartSpecialListGlobal[k].mass * fac * dz;

#ifdef EVALPOTENTIAL
          P[i].Potential += All.G * PartSpecialListGlobal[k].mass * wp;
#endif /* #ifdef EVALPOTENTIAL */
#ifdef EXACT_GRAVITY_REACTION
          /* avoid double counting */
          if(P[i].Type != EXACT_GRAVITY_FOR_PARTICLE_TYPE)
            {
              accx[k] += All.G * P[i].Mass * fac * dx;
              accy[k] += All.G * P[i].Mass * fac * dy;
              accz[k] += All.G * P[i].Mass * fac * dz;
#ifdef EVALPOTENTIAL
              pot[k] += All.G * P[i].Mass * wp;
#endif /* #ifdef EVALPOTENTIAL */
            }
#endif /* #ifdef EXACT_GRAVITY_REACTION */
        }
    }
#ifdef EXACT_GRAVITY_REACTION
  double *buf = (double *)mymalloc("buf", All.MaxPartSpecial * sizeof(double));

  MPI_Allreduce(accx, buf, All.MaxPartSpecial, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(n = 0; n < All.MaxPartSpecial; n++)
    accx[n] = buf[n];
  MPI_Allreduce(accy, buf, All.MaxPartSpecial, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(n = 0; n < All.MaxPartSpecial; n++)
    accy[n] = buf[n];
  MPI_Allreduce(accz, buf, All.MaxPartSpecial, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(n = 0; n < All.MaxPartSpecial; n++)
    accz[n] = buf[n];
#ifdef EVALPOTENTIAL
  MPI_Allreduce(pot, buf, All.MaxPartSpecial, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(n = 0; n < All.MaxPartSpecial; n++)
    pot[n] = buf[n];
#endif /* #ifdef EVALPOTENTIAL */
  myfree(buf);

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;
      for(n = 0; n < All.MaxPartSpecial; n++)
        {
          if(PartSpecialListGlobal[n].ID == P[i].ID)
            {
              P[i].GravAccel[0] += accx[n];
              P[i].GravAccel[1] += accy[n];
              P[i].GravAccel[2] += accz[n];
#ifdef EVALPOTENTIAL
              P[i].Potential += pot[n];
#endif /* #ifdef EVALPOTENTIAL */
            }
        }
    }

#ifdef EVALPOTENTIAL
  myfree(pot);
#endif /* #ifdef EVALPOTENTIAL */
  myfree(accz);
  myfree(accy);
  myfree(accx);
#endif /* #ifdef EXACT_GRAVITY_REACTION */
}

/*! \brief Creates list of special particles, i.e. particles for which gravity
 *  is calculated by direct summation.
 *
 *  Called in begrund2() (begrun.c), i.e. only at startup of the simulation.
 *
 *  \return void
 */
void special_particle_create_list()
{
  struct special_particle_data *SpecialPartList;
  SpecialPartList =
      (struct special_particle_data *)mymalloc("SpecialPartList", All.MaxPartSpecial * sizeof(struct special_particle_data));

  int i, j, nsrc, nimport, ngrp;
  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if(P[i].Type == EXACT_GRAVITY_FOR_PARTICLE_TYPE)
        {
          SpecialPartList[nsrc].ID = P[i].ID;

          SpecialPartList[nsrc].pos[0] = P[i].Pos[0];
          SpecialPartList[nsrc].pos[1] = P[i].Pos[1];
          SpecialPartList[nsrc].pos[2] = P[i].Pos[2];

          SpecialPartList[nsrc++].mass = P[i].Mass;
        }
    }

  for(j = 0; j < NTask; j++)
    Send_count[j] = nsrc;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = 0;
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&SpecialPartList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct special_particle_data),
                           MPI_BYTE, recvTask, TAG_DENS_A, &PartSpecialListGlobal[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct special_particle_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(SpecialPartList);
}

/*! \brief Updates list of special particles, i.e. particles for which gravity
 *  is calculated by direct summation.
 *
 *  Called in run() (run.c).
 *
 *  \return void
 */
void special_particle_update_list()
{
  struct special_particle_data *SpecialPartList;
  SpecialPartList =
      (struct special_particle_data *)mymalloc("SpecialPartList", All.MaxPartSpecial * sizeof(struct special_particle_data));

  int i, j, nsrc, nimport, ngrp;
  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if(P[i].Type == EXACT_GRAVITY_FOR_PARTICLE_TYPE)
        {
          SpecialPartList[nsrc].ID = P[i].ID;

          SpecialPartList[nsrc].pos[0] = P[i].Pos[0];
          SpecialPartList[nsrc].pos[1] = P[i].Pos[1];
          SpecialPartList[nsrc].pos[2] = P[i].Pos[2];

          SpecialPartList[nsrc++].mass = P[i].Mass;
        }
    }

  for(j = 0; j < NTask; j++)
    Send_count[j] = nsrc;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = 0;
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&SpecialPartList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct special_particle_data),
                           MPI_BYTE, recvTask, TAG_DENS_A, &PartSpecialListGlobal[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct special_particle_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(SpecialPartList);
}
#endif /* #ifdef  EXACT_GRAVITY_FOR_PARTICLE_TYPE */
