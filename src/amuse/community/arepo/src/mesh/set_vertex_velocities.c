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
 * \file        src/mesh/set_vertex_velocities.c
 * \date        05/2018
 * \brief       Algorithms that decide how individual cells are moving.
 * \details     contains functions:
 *                void set_vertex_velocities(void)
 *                static void validate_vertex_velocities_1d()
 *                void validate_vertex_velocities(void)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 08.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../mesh/voronoi/voronoi.h"

#ifdef ONEDIMS_SPHERICAL
static void validate_vertex_velocities_1d();
#endif /* #ifdef ONEDIMS_SPHERICAL */

/*! \brief Sets velocities of individual mesh-generating points.
 *
 *  \retur void
 */
void set_vertex_velocities(void)
{
  TIMER_START(CPU_SET_VERTEXVELS);

  int idx, i, j;
  double dt;

#if defined(VORONOI_STATIC_MESH) || defined(NOHYDRO)
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(j = 0; j < 3; j++)
        SphP[i].VelVertex[j] = 0;
    }
  TIMER_STOP(CPU_SET_VERTEXVELS);
  return;
#endif /* #if defined (VORONOI_STATIC_MESH) || defined (NOHYDRO) */

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef MESHRELAX
      for(j = 0; j < 3; j++)
        SphP[i].VelVertex[j] = 0;
#else  /* #ifdef MESHRELAX */
      for(j = 0; j < 3; j++)
        SphP[i].VelVertex[j] = P[i].Vel[j]; /* make cell velocity equal to fluid's velocity */
#endif /* #ifdef MESHRELAX #else */

      double acc[3];

      /*  the actual time-step of particle */
      integertime ti_step = P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0;
      dt                  = ti_step * All.Timebase_interval;
      dt /= All.cf_hubble_a; /* this gives the actual timestep: dt = dloga/ (adot/a) */

      /* now let's add the gradient of the pressure force
       * note that the gravity half-step was already included in P[i].Vel[j]
       * prior to calling this function, thus it does not need to be accounted
       * here explicitly.
       */
      if(SphP[i].Density > 0)
        {
          acc[0] = -SphP[i].Grad.dpress[0] / SphP[i].Density;
          acc[1] = -SphP[i].Grad.dpress[1] / SphP[i].Density;
          acc[2] = -SphP[i].Grad.dpress[2] / SphP[i].Density;

#ifdef MHD
          /* we also add the acceleration due to the Lorentz force */
          acc[0] += (SphP[i].CurlB[1] * SphP[i].B[2] - SphP[i].CurlB[2] * SphP[i].B[1]) / SphP[i].Density;
          acc[1] += (SphP[i].CurlB[2] * SphP[i].B[0] - SphP[i].CurlB[0] * SphP[i].B[2]) / SphP[i].Density;
          acc[2] += (SphP[i].CurlB[0] * SphP[i].B[1] - SphP[i].CurlB[1] * SphP[i].B[0]) / SphP[i].Density;

#endif /* #ifdef MHD */

          SphP[i].VelVertex[0] += 0.5 * dt * acc[0];
          SphP[i].VelVertex[1] += 0.5 * dt * acc[1];
          SphP[i].VelVertex[2] += 0.5 * dt * acc[2];
        }
    }

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef REGULARIZE_MESH_CM_DRIFT

      double dx, dy, dz, d, fraction;

      dx = nearest_x(P[i].Pos[0] - SphP[i].Center[0]);
      dy = nearest_y(P[i].Pos[1] - SphP[i].Center[1]);
      dz = nearest_z(P[i].Pos[2] - SphP[i].Center[2]);

      /*  the actual time-step of particle */
      dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      dt /= All.cf_hubble_a; /* this is dt, the actual timestep  */

      double cellrad = get_cell_radius(i);

#if !defined(REGULARIZE_MESH_FACE_ANGLE)
      /* if there is a density gradient, use a center that is displaced slightly in the direction of the gradient.
       * This makes sure that the Lloyd scheme does not simply iterate towards cells of equal volume, instead
       * we keep cells of roughly equal mass.
       */
      double dgrad = sqrt(SphP[i].Grad.drho[0] * SphP[i].Grad.drho[0] + SphP[i].Grad.drho[1] * SphP[i].Grad.drho[1] +
                          SphP[i].Grad.drho[2] * SphP[i].Grad.drho[2]);

      if(dgrad > 0)
        {
          double scale = SphP[i].Density / dgrad;
          double tmp   = 3 * cellrad + scale;
          double x     = (tmp - sqrt(tmp * tmp - 8 * cellrad * cellrad)) / 4;

          if(x < 0.25 * cellrad)
            {
              dx = nearest_x(P[i].Pos[0] - (SphP[i].Center[0] + x * SphP[i].Grad.drho[0] / dgrad));
              dy = nearest_y(P[i].Pos[1] - (SphP[i].Center[1] + x * SphP[i].Grad.drho[1] / dgrad));
              dz = nearest_z(P[i].Pos[2] - (SphP[i].Center[2] + x * SphP[i].Grad.drho[2] / dgrad));
            }
        }
#endif /* #if !defined(REGULARIZE_MESH_FACE_ANGLE) */

      d = sqrt(dx * dx + dy * dy + dz * dz);

      fraction = 0;

#if !defined(REGULARIZE_MESH_FACE_ANGLE)
      if(d > 0.75 * All.CellShapingFactor * cellrad && dt > 0)
        {
          if(d > All.CellShapingFactor * cellrad)
            fraction = All.CellShapingSpeed;
          else
            fraction = All.CellShapingSpeed * (d - 0.75 * All.CellShapingFactor * cellrad) / (0.25 * All.CellShapingFactor * cellrad);
        }
#else  /* #if !defined(REGULARIZE_MESH_FACE_ANGLE) */
      if(SphP[i].MaxFaceAngle > 0.75 * All.CellMaxAngleFactor && dt > 0)
        {
          if(SphP[i].MaxFaceAngle > All.CellMaxAngleFactor)
            fraction = All.CellShapingSpeed;
          else
            fraction = All.CellShapingSpeed * (SphP[i].MaxFaceAngle - 0.75 * All.CellMaxAngleFactor) / (0.25 * All.CellMaxAngleFactor);
        }
#endif /* #if !defined(REGULARIZE_MESH_FACE_ANGLE) #else */

      if(d > 0 && fraction > 0)
        {
          double v;
#ifdef REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED

          v = All.cf_atime * get_sound_speed(i);

#if defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)
          /* calculate gravitational velocity scale */
          double ax, ay, az, ac, vgrav;
#ifdef HIERARCHICAL_GRAVITY
          ax = SphP[i].FullGravAccel[0];
          ay = SphP[i].FullGravAccel[1];
          az = SphP[i].FullGravAccel[2];
#else  /* #ifdef HIERARCHICAL_GRAVITY */
          ax = P[i].GravAccel[0];
          ay = P[i].GravAccel[1];
          az = P[i].GravAccel[2];
#endif /* #ifdef HIERARCHICAL_GRAVITY #else */
#ifdef PMGRID
          ax += P[i].GravPM[0];
          ay += P[i].GravPM[1];
          az += P[i].GravPM[2];
#endif /* #ifdef PMGRID */
          ac    = sqrt(ax * ax + ay * ay + az * az);
          vgrav = 4 * sqrt(All.cf_atime * cellrad * ac);
          if(v < vgrav)
            v = vgrav;
#endif /* #if defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE) */

          double vcurl = cellrad * SphP[i].CurlVel;
          if(v < vcurl)
            v = vcurl;

#else  /* #ifdef REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED */
          v = All.cf_atime * All.cf_atime * d / dt; /* use fiducial velocity */

          double vel  = sqrt(P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
          double vmax = dmax(All.cf_atime * get_sound_speed(i), vel);
          if(v > vmax)
            v = vmax;
#endif /* #ifdef REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED #else */

#ifdef REFINEMENT_SPLIT_CELLS
          double proj = SphP[i].SepVector[0] * dx + SphP[i].SepVector[1] * dy + SphP[i].SepVector[2] * dz;

          if(proj != 0)
            {
              dx = proj * SphP[i].SepVector[0];
              dy = proj * SphP[i].SepVector[1];
              dz = proj * SphP[i].SepVector[2];
            }

          SphP[i].SepVector[0] = 0;
          SphP[i].SepVector[1] = 0;
          SphP[i].SepVector[2] = 0;
#endif /* #ifdef REFINEMENT_SPLIT_CELLS */

          SphP[i].VelVertex[0] += fraction * v * (-dx / d);
          SphP[i].VelVertex[1] += fraction * v * (-dy / d);
          SphP[i].VelVertex[2] += fraction * v * (-dz / d);
        }
#endif /* #ifdef REGULARIZE_MESH_CM_DRIFT */

      for(j = NUMDIMS; j < 3; j++)
        SphP[i].VelVertex[j] = 0; /* vertex velocities for unused dimensions set to zero */
    }

#ifdef OUTPUT_VERTEX_VELOCITY_DIVERGENCE
  voronoi_exchange_primitive_variables();
  calculate_vertex_velocity_divergence();
#endif /* #ifdef OUTPUT_VERTEX_VELOCITY_DIVERGENCE */

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
  validate_vertex_velocities();
#endif /* #if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z) */

#ifdef ONEDIMS_SPHERICAL
  validate_vertex_velocities_1d();
#endif /* #ifdef ONEDIMS_SPHERICAL */

  TIMER_STOP(CPU_SET_VERTEXVELS);
}

#ifdef ONEDIMS_SPHERICAL
/*! \brief Handles inner boundary cells in 1d spherical case.
 *
 *  \return void
 */
static void validate_vertex_velocities_1d()
{
  double dt = (P[0].TimeBinHydro ? (((integertime)1) << P[0].TimeBinHydro) : 0) * All.Timebase_interval;
  if(P[0].Pos[0] + dt * SphP[0].VelVertex[0] < All.CoreRadius)
    SphP[0].VelVertex[0] = 0.;
}
#endif /* #ifdef ONEDIMS_SPHERICAL */

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
/*! \brief Checks validity of vertex velocities with boundary conditions.
 *
 *  In case we have reflecting boundaries, make sure that cell does not drift
 *  beyond boundary.
 *
 *  \return void
 */
void validate_vertex_velocities(void)
{
  int idx, i;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      integertime ti_step = P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0;
      double dt_drift;

      if(All.ComovingIntegrationOn)
        dt_drift = get_drift_factor(All.Ti_Current, All.Ti_Current + ti_step);
      else
        dt_drift = ti_step * All.Timebase_interval;

#if defined(REFLECTIVE_X)
      if((P[i].Pos[0] + dt_drift * SphP[i].VelVertex[0]) < 0 || (P[i].Pos[0] + dt_drift * SphP[i].VelVertex[0]) >= boxSize_X)
        SphP[i].VelVertex[0] = 0;
#endif /* #if defined(REFLECTIVE_X) */
#if defined(REFLECTIVE_Y)
      if((P[i].Pos[1] + dt_drift * SphP[i].VelVertex[1]) < 0 || (P[i].Pos[1] + dt_drift * SphP[i].VelVertex[1]) >= boxSize_Y)
        SphP[i].VelVertex[1] = 0;
#endif /* #if defined(REFLECTIVE_Y) */
#if defined(REFLECTIVE_Z)
      if((P[i].Pos[2] + dt_drift * SphP[i].VelVertex[2]) < 0 || (P[i].Pos[2] + dt_drift * SphP[i].VelVertex[2]) >= boxSize_Z)
        SphP[i].VelVertex[2] = 0;
#endif /* #if defined(REFLECTIVE_Z) */
    }
}
#endif /* #if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z) */
