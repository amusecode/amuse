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
 * \file        src/finite_volume_solver.c
 * \date        05/2018
 * \brief       Core algorithms of the finite-volume solver.
 * \details     contains functions:
 *                void compute_interface_fluxes(tessellation * T)
 *                void backup_face_areas(tessellation * T)
 *                void restore_face_areas(tessellation * T)
 *                int face_get_state(tessellation * T, int p, int i, struct
 *                  state *st)
 *                void face_boundary_check_vertex(tessellation * T, int p,
 *                  MyFloat * velx, MyFloat * vely, MyFloat * velz)
 *                void face_boundary_check(point * p, double *velx, double
 *                  *vely, double *velz)
 *                int face_check_responsibility_of_this_task(tessellation * T,
 *                  int p1, int p2, struct state *st_L, struct state *st_R)
 *                double face_timestep(struct state *state_L, struct state
 *                  *state_R, double *hubble_a, double *atime)
 *                void state_convert_to_local_frame(struct state *st, double
 *                  *vel_face, double hubble_a, double atime)
 *                void face_do_time_extrapolation(struct state *delta,
 *                  struct state *st, double atime)
 *                void face_do_spatial_extrapolation(struct state *delta,
 *                  struct state *st, struct state *st_other)
 *                void face_do_spatial_extrapolation_single_quantity(double
 *                  *delta, double st, double st_other, MySingle * grad,
 *                  double *dx, double *r)
 *                void face_add_extrapolations(struct state *st_face, struct
 *                  state *delta_time, struct state *delta_space, struct
 *                  fvs_stat *stat)
 *                void face_add_extrapolation(struct state *st_face, struct
 *                  state *delta, struct fvs_stat *stat)
 *                void face_add_extrapolation_with_check(struct state *st_face,
 *                  struct state *delta, struct fvs_stat *stat)
 *                void face_turn_velocities(struct state *st, struct geometry
 *                  *geom)
 *                void solve_advection(struct state *st_L, struct state *st_R,
 *                  struct state_face *st_face, struct geometry *geom,
 *                  double *vel_face)
 *                void face_turnback_velocities(struct state_face *st_face,
 *                  struct geometry *geom)
 *                void face_set_scalar_states_and_fluxes(struct state *st_L,
 *                  struct state *st_R, struct state_face *st_face, struct
 *                  fluxes *flux)
 *                void flux_convert_to_lab_frame(struct state *st_L, struct
 *                  state *st_R, double *vel_face, struct fluxes *flux)
 *                void face_turn_momentum_flux(struct fluxes *flux, struct
 *                  geometry *geom)
 *                void face_get_fluxes(struct state *st_L, struct state *st_R,
 *                  struct state_face *st_face, struct fluxes *flux, struct
 *                  geometry *geom, double *vel_face)
 *                void face_limit_fluxes(struct state *st_L, struct state
 *                  *st_R, struct state *st_center_L, struct state
 *                  *st_center_R, struct fluxes *flux, double dt, double
 *                  *count, double *count_reduced)
 *                void face_clear_fluxes(struct fluxes *flux)
 *                void face_add_fluxes_advection(struct state_face *st_face,
 *                  struct fluxes *flux, struct geometry *geom, double
 *                  *vel_face)
 *                int flux_list_data_compare(const void *a, const void *b)
 *                void apply_flux_list(void)
 *                void fvs_initialize_statistics(struct fvs_stat *stat)
 *                void fvs_evaluate_statistics(struct fvs_stat *stat)
 *                void apply_spherical_source_terms()
 *                void add_spin_source_term_from_grid_movement()
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 17.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../mesh/voronoi/voronoi.h"

/*! \brief Data needed for flux calculation.
 */
static struct flux_list_data
{
  int task, index;
  double dM, dP[3];
#ifdef MHD
  double dB[3];
#endif /* #ifdef MHD */

#ifndef ISOTHERM_EQS
  double dEnergy;
#endif /* #ifndef ISOTHERM_EQS */
#ifdef MAXSCALARS
  double dConservedScalars[MAXSCALARS];
#endif /* #ifdef MAXSCALARS */
} * FluxList;

static int Nflux, MaxNflux;

struct primexch *PrimExch;
struct grad_data *GradExch;

/*! state on a face determined by Riemann solver */
struct state_face state_face;

/*! flux through a face */
struct fluxes fluxes;

struct geometry geom;

#ifdef ONEDIMS_SPHERICAL
void apply_spherical_source_terms();
#endif /* #ifdef ONEDIMS_SPHERICAL */

static void face_add_extrapolation_with_check(struct state *st_face, struct state *delta, struct fvs_stat *stat);
static void fvs_initialize_statistics(struct fvs_stat *stat);
static void fvs_evaluate_statistics(struct fvs_stat *stat);

#ifdef VORONOI_BACKUP_RESTORE_FACE_AREAS
void backup_face_areas(tessellation *T);
void restore_face_areas(tessellation *T);
#endif /* #ifdef VORONOI_BACKUP_RESTORE_FACE_AREAS */

/*! \brief Main routine to compute fluxes across interfaces given am mesh T.
 *
 *  Adds these fluxes to conserved variables.
 *
 *  \param[in] T Pointer to tessellation.
 *
 *  \return void
 */
void compute_interface_fluxes(tessellation *T)
{
#ifdef NOHYDRO
  return;
#endif /* #ifdef NOHYDRO */
  TIMER_START(CPU_FLUXES);

  int i, j;
  double count = 0, count_reduced = 0, tot_count, tot_count_reduced;
  double face_dt, hubble_a, atime;
  struct fvs_stat stat;
#ifdef MHD
  double sqrtatime;
#endif /* #ifdef MHD */

#ifdef GODUNOV_STATS
  FILE *fdstats;
  char buf[1000];

  sprintf(buf, "%s/godunov_stats_%d.txt", All.OutputDir, ThisTask);
  if(!(fdstats = fopen(buf, "w")))
    terminate("error in opening file '%s'", buf);
#endif /* #ifdef GODUNOV_STATS */

#ifdef VORONOI_BACKUP_RESTORE_FACE_AREAS
  backup_face_areas(T);
#endif /* #ifdef VORONOI_BACKUP_RESTORE_FACE_AREAS */

  fvs_initialize_statistics(&stat);

  MaxNflux = T->Indi.AllocFacNflux;
  Nflux    = 0;
  FluxList = mymalloc_movable(&FluxList, "FluxList", MaxNflux * sizeof(struct flux_list_data));

  face *VF  = T->VF;
  point *DP = T->DP;

  for(i = 0; i < T->Nvf; i++)
    {
      struct state state_L, state_center_L, delta_time_L, delta_space_L;
      struct state state_R, state_center_R, delta_time_R, delta_space_R;

      face_dt = 0; /* the default is that this face is not active */

      /* calculate normal vectors */
      if(face_get_normals(T, i, &geom))
        continue;

      /* get the values of the states at the center of the cells */
      if(face_get_state(T, VF[i].p1, i, &state_center_L))
        continue;

      if(face_get_state(T, VF[i].p2, i, &state_center_R))
        continue;

      /* only treat faces where one of the two sides is active */
      if(!TimeBinSynchronized[state_center_L.timeBin] && !TimeBinSynchronized[state_center_R.timeBin])
        continue;

      /* clarify whether the face should be done by this task (it may be present also on another task) */
      if(face_check_responsibility_of_this_task(T, VF[i].p1, VF[i].p2, &state_center_L, &state_center_R))
        continue;

      /* calculate timestep of the face */
      face_dt = face_timestep(&state_center_L, &state_center_R, &hubble_a, &atime);
#ifdef MHD
      sqrtatime = sqrt(atime);
#endif /* #ifdef MHD */

      if(!(face_dt > 0))
        continue;

      /* now estimate the velocity of the midpoint of the face based on the velocities of the generators of the mesh. */
      double vel_face[3];

      if(All.ComovingIntegrationOn)
        for(j = 0; j < 3; j++)
          {
            state_center_L.velVertex[j] /= atime; /* convert vertex motion to peculiar velocity */
            state_center_R.velVertex[j] /= atime;
          }

      /* rough motion of mid-point of edge */
      vel_face[0] = 0.5 * (state_center_L.velVertex[0] + state_center_R.velVertex[0]);
      vel_face[1] = 0.5 * (state_center_L.velVertex[1] + state_center_R.velVertex[1]);
      vel_face[2] = 0.5 * (state_center_L.velVertex[2] + state_center_R.velVertex[2]);

      double cx, cy, cz, facv;

      cx = VF[i].cx - 0.5 * (DP[VF[i].p2].x + DP[VF[i].p1].x);
      cy = VF[i].cy - 0.5 * (DP[VF[i].p2].y + DP[VF[i].p1].y);
      cz = VF[i].cz - 0.5 * (DP[VF[i].p2].z + DP[VF[i].p1].z);

      facv = (cx * (state_center_L.velVertex[0] - state_center_R.velVertex[0]) +
              cy * (state_center_L.velVertex[1] - state_center_R.velVertex[1]) +
              cz * (state_center_L.velVertex[2] - state_center_R.velVertex[2])) /
             geom.nn;

      /* put in a limiter for highly distorted cells */
      double cc = sqrt(cx * cx + cy * cy + cz * cz);
      if(cc > 0.9 * geom.nn)
        facv *= (0.9 * geom.nn) / cc;

      vel_face[0] += facv * geom.nx;
      vel_face[1] += facv * geom.ny;
      vel_face[2] += facv * geom.nz;

#if defined(VORONOI_STATIC_MESH)
      vel_face[0] = 0;
      vel_face[1] = 0;
      vel_face[2] = 0;
#endif /* #if defined(VORONOI_STATIC_MESH) */

#if defined(RIEMANN_HLLC) || defined(RIEMANN_HLLD)
      double vel_face_turned[3];
      /* for these riemann solvers, the riemann problem is not solved in the
       * restframe of the face, instead the mesh motion is accounted for via
       * an advection step.
       */

      /* turn the face velocity */
      vel_face_turned[0] = vel_face[0] * geom.nx + vel_face[1] * geom.ny + vel_face[2] * geom.nz;
      vel_face_turned[1] = vel_face[0] * geom.mx + vel_face[1] * geom.my + vel_face[2] * geom.mz;
      vel_face_turned[2] = vel_face[0] * geom.px + vel_face[1] * geom.py + vel_face[2] * geom.pz;
#endif /* #if defined(RIEMANN_HLLC) || defined(RIEMANN_HLLD) */

      state_convert_to_local_frame(&state_center_L, vel_face, hubble_a, atime);
      state_convert_to_local_frame(&state_center_R, vel_face, hubble_a, atime);

      /* copy center state to state at interface, then add extrapolation terms */
      state_L = state_center_L;
      state_R = state_center_R;

      face_do_time_extrapolation(&delta_time_L, &state_center_L, atime);
      face_do_time_extrapolation(&delta_time_R, &state_center_R, atime);

      face_do_spatial_extrapolation(&delta_space_L, &state_center_L, &state_center_R);
      face_do_spatial_extrapolation(&delta_space_R, &state_center_R, &state_center_L);

      face_add_extrapolations(&state_L, &delta_time_L, &delta_space_L, &stat);
      face_add_extrapolations(&state_R, &delta_time_R, &delta_space_R, &stat);

#ifdef MHD
      if(All.ComovingIntegrationOn)
        {
          state_L.Bx /= sqrtatime;
          state_L.By /= sqrtatime;
          state_L.Bz /= sqrtatime;

          state_R.Bx /= sqrtatime;
          state_R.By /= sqrtatime;
          state_R.Bz /= sqrtatime;
        }
#endif /* #ifdef MHD */

#ifndef MESHRELAX
#ifndef ISOTHERM_EQS
      /* check for crazy values */
      if(state_L.press < 0 || state_R.press < 0 || state_L.rho < 0 || state_R.rho < 0)
        {
          printf("i=%d press_L=%g press_R=%g rho_L=%g rho_R=%g\n", i, state_L.press, state_R.press, state_L.rho, state_R.rho);
          printf("area=%g lx=%g ly=%g   rx=%g ry=%g\n", VF[i].area, state_L.dx, state_L.dy, state_R.dx, state_R.dy);
          terminate("found crazy values");
        }
#else  /* #ifndef ISOTHERM_EQS */
      if(state_L.press < 0 || state_R.press < 0 || state_L.rho < 0 || state_R.rho < 0)
        {
          printf("i=%d rho_L=%g rho_R=%g\n", i, state_L.rho, state_R.rho);
          printf("area=%g lx=%g ly=%g   rx=%g ry=%g\n", VF[i].area, state_L.dx, state_L.dy, state_R.dx, state_R.dy);
          terminate("found crazy values");
        }
#endif /* #ifndef ISOTHERM_EQS #else */
#endif /* #ifndef MESHRELAX */

      /* mirror velocity in case of reflecting boundaries */
      face_boundary_check(&T->DP[VF[i].p1], &state_L.velx, &state_L.vely, &state_L.velz);
      face_boundary_check(&T->DP[VF[i].p2], &state_R.velx, &state_R.vely, &state_R.velz);

#ifdef MHD
      /* mirror magnetic field in case of reflecting boundaries */
      face_boundary_check(&T->DP[VF[i].p1], &state_L.Bx, &state_L.By, &state_L.Bz);
      face_boundary_check(&T->DP[VF[i].p2], &state_R.Bx, &state_R.By, &state_R.Bz);
#endif /* #ifdef MHD */

      /* turn the velocities to get velx perpendicular and vely and velz in the plane of the face */
      face_turn_velocities(&state_L, &geom);
      face_turn_velocities(&state_R, &geom);

#ifndef MESHRELAX

      /* call Riemann solver */

      double press;
#ifdef RIEMANN_HLLC
      press = godunov_flux_3d_hllc(&state_L, &state_R, &state_face, &fluxes);
#else /* #ifdef RIEMANN_HLLC */
#ifdef RIEMANN_HLLD
      press = godunov_flux_3d_hlld(&state_L, &state_R, vel_face_turned, &state_face, &fluxes);
#else  /* #ifdef RIEMANN_HLLD */
      press = godunov_flux_3d(&state_L, &state_R, &state_face); /* exact ideal gas solver */
#endif /* #ifdef RIEMANN_HLLD #else */
#endif /* #ifdef RIEMANN_HLLC #else */

      if(press < 0)
        terminate("press < 0: ID_L: %d, ID_R: %d", VF[i].p1, VF[i].p2);

#ifdef GODUNOV_STATS
      get_mach_numbers(&state_L, &state_R, press);
      if(st_L.rho > 1.0e-6 && st_R.rho > 1.0e-6)
        fprintf(fdstats, "%g %g %g   %g %g %g  %g %g %g  %g %g %g\n", state_L.rho, state_L.velx, state_L.press, state_L.rho,
                state_L.velx, state_L.press, state_face.rho, state_face.velx, state_face.press, state_L.mach, state_R.mach,
                VF[i].area);
#endif /* GODUNOV_STATS */

#endif /* #ifndef MESHRELAX */

      /* turn the velocity field back */
      face_turnback_velocities(&state_face, &geom);

      /* add the face velocity again */
      state_face.velx += vel_face[0];
      state_face.vely += vel_face[1];
      state_face.velz += vel_face[2];

#ifndef MESHRELAX

#if defined(RIEMANN_HLLC) || defined(RIEMANN_HLLD)
      /* for non-exact Riemann solver, fluxes are already computed in the local frame, so convert to lab frame and turn momentum fluxes
       * to the lab orientation  */
      flux_convert_to_lab_frame(&state_L, &state_R, vel_face_turned, &fluxes);
      face_turn_momentum_flux(&fluxes, &geom);

#else /* #if defined(RIEMANN_HLLC) || defined(RIEMANN_HLLD) */

      /* calculate fluxes for exact Riemann problem */
      /* compute net flux with dot-product of outward normal and area of face */
      /* multiplication with area and time-step comes later */

      face_get_fluxes(&state_L, &state_R, &state_face, &fluxes, &geom, vel_face);

#endif /* #if defined(RIEMANN_HLLC) || defined(RIEMANN_HLLD)  #else */

      /* set the face states and fluxes of those quantities that are passively advected */
      face_set_scalar_states_and_fluxes(&state_L, &state_R, &state_face, &fluxes);

      face_limit_fluxes(&state_L, &state_R, &state_center_L, &state_center_R, &fluxes, face_dt, &count, &count_reduced);

      /* put in cosmological factors */
      if(All.ComovingIntegrationOn)
        {
          fluxes.momentum[0] *= atime;
          fluxes.momentum[1] *= atime;
          fluxes.momentum[2] *= atime;
          fluxes.energy *= atime * atime;
#ifdef MHD
          fluxes.B[0] *= sqrtatime;
          fluxes.B[1] *= sqrtatime;
          fluxes.B[2] *= sqrtatime;
#ifdef MHD_POWELL
          state_face.Bx *= sqrtatime;
#endif /* #ifdef MHD_POWELL */
#endif /* #ifdef MHD */
        }

#else /* #ifndef MESHRELAX */

      /* just solve the advection equation instead of Riemann problem */

      solve_advection(&state_L, &state_R, &state_face, &geom, vel_face);
      face_clear_fluxes(&fluxes);
      face_add_fluxes_advection(&state_face, &fluxes, &geom, vel_face);
      face_set_scalar_states_and_fluxes(&state_L, &state_R, &state_face, &fluxes);

#endif /* #ifndef MESHRELAX #else */

#ifndef ISOTHERM_EQS
      if(!gsl_finite(fluxes.energy))
        {
          printf("i=%d eFlux-Bummer: %g %g %g\n", i, fluxes.energy, state_face.press, state_face.rho);
          printf("rho_L=%g velx_L=%g vely_L=%g velz_L=%g press_L=%g\n", state_L.rho, state_L.velx, state_L.vely, state_L.velz,
                 state_L.press);
          printf("rho_R=%g velx_R=%g vely_R=%g velz_R=%g press_R=%g\n", state_R.rho, state_R.velx, state_R.vely, state_R.velz,
                 state_R.press);
          print_particle_info(i);
          terminate("infinity encountered");
        }
#endif /* #ifndef ISOTHERM_EQS */

      /* now apply the flux to update the conserved states of the cells */

      if(face_dt > 0) /* selects active faces */
        {
          int k, p, q;
          double dir;
          double fac = face_dt * VF[i].area;
#if defined(MAXSCALARS)
          int m;
#endif /* #if defined(MAXSCALARS) */

          fac *= 0.5;

#if defined(MHD_POWELL)
          struct state *state_center, *delta_time;
#endif /* #if defined(MHD_POWELL) */
          for(k = 0; k < 2; k++)
            {
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
              int qother;
#endif /* #if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z) */
              if(k == 0)
                {
                  q   = VF[i].p1;
                  p   = DP[q].index;
                  dir = -fac;
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
                  qother = VF[i].p2;
#endif /* #if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z) */
#if defined(MHD_POWELL)
                  state_center = &state_center_L;
                  delta_time   = &delta_time_L;
#endif /* #if defined(MHD_POWELL) */
                }
              else
                {
                  q   = VF[i].p2;
                  p   = DP[q].index;
                  dir = +fac;
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
                  qother = VF[i].p1;
#endif /* #if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z) */
#if defined(MHD_POWELL)
                  state_center = &state_center_R;
                  delta_time   = &delta_time_R;
#endif /* #if defined(MHD_POWELL) */
                }

              if(DP[q].task == ThisTask)
                {
                  if(DP[q].index >= NumGas) /* this is a local ghost point */
                    {
                      if(DP[VF[i].p1].ID == DP[VF[i].p2].ID) /* this may happen for reflective points */
                        continue;
                      p -= NumGas;
                    }

                  /* note: this will be executed if P[p] is a local point, independent of active or not */
                  P[p].Mass += dir * fluxes.mass;
                  SphP[p].Momentum[0] += dir * fluxes.momentum[0];
                  SphP[p].Momentum[1] += dir * fluxes.momentum[1];
                  SphP[p].Momentum[2] += dir * fluxes.momentum[2];

#ifdef MHD
                  SphP[p].BConserved[0] += dir * fluxes.B[0];
                  SphP[p].BConserved[1] += dir * fluxes.B[1];
                  SphP[p].BConserved[2] += dir * fluxes.B[2];
#if defined(MHD_POWELL)
                  double Velx = state_center->velx + delta_time->velx + vel_face[0];
                  double Vely = state_center->vely + delta_time->vely + vel_face[1];
                  double Velz = state_center->velz + delta_time->velz + vel_face[2];

                  if(All.ComovingIntegrationOn)
                    {
                      Velx += atime * hubble_a * state_center->dx;
                      Vely += atime * hubble_a * state_center->dy;
                      Velz += atime * hubble_a * state_center->dz;
                    }

                  double Bx = state_center->Bx + delta_time->Bx;
                  double By = state_center->By + delta_time->By;
                  double Bz = state_center->Bz + delta_time->Bz;

                  SphP[p].BConserved[0] += dir * Velx * state_face.Bx;
                  SphP[p].BConserved[1] += dir * Vely * state_face.Bx;
                  SphP[p].BConserved[2] += dir * Velz * state_face.Bx;

                  SphP[p].Momentum[0] += dir * Bx * state_face.Bx;
                  SphP[p].Momentum[1] += dir * By * state_face.Bx;
                  SphP[p].Momentum[2] += dir * Bz * state_face.Bx;

                  SphP[p].Energy += dir * (Bx * Velx + By * Vely + Bz * Velz) * state_face.Bx * atime;

                  {
                    double dMomX = dir * Bx * state_face.Bx;
                    double dMomY = dir * By * state_face.Bx;
                    double dMomZ = dir * Bz * state_face.Bx;

                    All.Powell_Momentum[0] += dMomX;
                    All.Powell_Momentum[1] += dMomY;
                    All.Powell_Momentum[2] += dMomZ;

                    double dx = SphP[p].Center[0] - 0.5 * All.BoxSize;
                    double dy = SphP[p].Center[1] - 0.5 * All.BoxSize;
                    double dz = SphP[p].Center[2] - 0.5 * All.BoxSize;

                    All.Powell_Angular_Momentum[0] += dy * dMomZ - dz * dMomY;
                    All.Powell_Angular_Momentum[1] += dz * dMomX - dx * dMomZ;
                    All.Powell_Angular_Momentum[2] += dx * dMomY - dy * dMomX;
                    All.Powell_Energy += dir * (Bx * Velx + By * Vely + Bz * Velz) * state_face.Bx * atime;
                  }
#endif /* #if defined(MHD_POWELL) */
#endif /* #ifdef MHD */

#ifdef MAXSCALARS
                  for(m = 0; m < N_Scalar; m++)
                    {
                      *(MyFloat *)(((char *)(&SphP[p])) + scalar_elements[m].offset_mass) += dir * fluxes.scalars[m];
                    }
#endif /* #ifdef MAXSCALARS */

#if !defined(ISOTHERM_EQS)
                  SphP[p].Energy += dir * fluxes.energy;
#endif /* #if !defined(ISOTHERM_EQS)  */
                }
              else
                {
                  /* here we have a foreign ghost point */
                  if(DP[q].originalindex < 0)
                    terminate("should not happen");

                  if(Nflux >= MaxNflux)
                    {
                      T->Indi.AllocFacNflux *= ALLOC_INCREASE_FACTOR;
                      MaxNflux = T->Indi.AllocFacNflux;
#ifdef VERBOSE
                      printf("Task=%d: increase memory allocation, MaxNflux=%d Indi.AllocFacNflux=%g\n", ThisTask, MaxNflux,
                             T->Indi.AllocFacNflux);
#endif /* #ifdef VERBOSE */
                      FluxList = myrealloc_movable(FluxList, MaxNflux * sizeof(struct flux_list_data));

                      if(Nflux >= MaxNflux)
                        terminate("Nflux >= MaxNflux");
                    }

                  FluxList[Nflux].task  = DP[q].task;
                  FluxList[Nflux].index = DP[q].originalindex;

                  FluxList[Nflux].dM = dir * fluxes.mass;

                  FluxList[Nflux].dP[0] = dir * fluxes.momentum[0];
                  FluxList[Nflux].dP[1] = dir * fluxes.momentum[1];
                  FluxList[Nflux].dP[2] = dir * fluxes.momentum[2];

#if !defined(ISOTHERM_EQS)
                  FluxList[Nflux].dEnergy = dir * fluxes.energy;
#endif /* #if !defined(ISOTHERM_EQS)  */

#ifdef MHD
                  FluxList[Nflux].dB[0] = dir * fluxes.B[0];
                  FluxList[Nflux].dB[1] = dir * fluxes.B[1];
                  FluxList[Nflux].dB[2] = dir * fluxes.B[2];
#if defined(MHD_POWELL)
                  double Velx = state_center->velx + delta_time->velx + vel_face[0];
                  double Vely = state_center->vely + delta_time->vely + vel_face[1];
                  double Velz = state_center->velz + delta_time->velz + vel_face[2];

                  if(All.ComovingIntegrationOn)
                    {
                      Velx += atime * hubble_a * state_center->dx;
                      Vely += atime * hubble_a * state_center->dy;
                      Velz += atime * hubble_a * state_center->dz;
                    }

                  double Bx = state_center->Bx + delta_time->Bx;
                  double By = state_center->By + delta_time->By;
                  double Bz = state_center->Bz + delta_time->Bz;

                  FluxList[Nflux].dB[0] += dir * Velx * state_face.Bx;
                  FluxList[Nflux].dB[1] += dir * Vely * state_face.Bx;
                  FluxList[Nflux].dB[2] += dir * Velz * state_face.Bx;

                  FluxList[Nflux].dP[0] += dir * Bx * state_face.Bx;
                  FluxList[Nflux].dP[1] += dir * By * state_face.Bx;
                  FluxList[Nflux].dP[2] += dir * Bz * state_face.Bx;
#ifndef ISOTHERM_EQS
                  FluxList[Nflux].dEnergy += dir * (Bx * Velx + By * Vely + Bz * Velz) * state_face.Bx * atime;
#endif /* #ifndef ISOTHERM_EQS */

                  {
                    double dMomX = dir * Bx * state_face.Bx;
                    double dMomY = dir * By * state_face.Bx;
                    double dMomZ = dir * Bz * state_face.Bx;

                    All.Powell_Momentum[0] += dMomX;
                    All.Powell_Momentum[1] += dMomY;
                    All.Powell_Momentum[2] += dMomZ;

                    double dx = PrimExch[p].Center[0] - 0.5 * All.BoxSize;
                    double dy = PrimExch[p].Center[1] - 0.5 * All.BoxSize;
                    double dz = PrimExch[p].Center[2] - 0.5 * All.BoxSize;

                    All.Powell_Angular_Momentum[0] += dy * dMomZ - dz * dMomY;
                    All.Powell_Angular_Momentum[1] += dz * dMomX - dx * dMomZ;
                    All.Powell_Angular_Momentum[2] += dx * dMomY - dy * dMomX;
                    All.Powell_Energy += dir * (Bx * Velx + By * Vely + Bz * Velz) * state_face.Bx * atime;
                  }
#endif /* #if defined(MHD_POWELL) */
#endif /* #ifdef MHD */

#ifdef MAXSCALARS
                  for(m = 0; m < N_Scalar; m++)
                    FluxList[Nflux].dConservedScalars[m] = dir * fluxes.scalars[m];
#endif /* #ifdef MAXSCALARS */

                  Nflux++;
                }
            }
        }
    }
  /* end of big loop over all faces */

  TIMER_STOPSTART(CPU_FLUXES, CPU_FLUXES_COMM);

  /* now exchange the flux-list and apply it when needed */
  apply_flux_list();

  TIMER_STOPSTART(CPU_FLUXES_COMM, CPU_FLUXES);

  myfree(FluxList);

  double in[2] = {count, count_reduced}, out[2];
  MPI_Reduce(in, out, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0)
    {
      tot_count         = out[0];
      tot_count_reduced = out[1];

      printf("FLUX: exchanged fluxes over %g faces, with %g reduced (fraction %g), cumulative fraction %g\n", tot_count,
             tot_count_reduced, tot_count_reduced / (tot_count + 1.0e-30), All.TotCountReducedFluxes / (All.TotCountFluxes + 1.0e-30));
      All.TotCountReducedFluxes += tot_count_reduced;
      All.TotCountFluxes += tot_count;
    }

  fvs_evaluate_statistics(&stat);

#ifdef MESHRELAX
  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Mass < 0)
        {
          terminate("negative mass reached for cell=%d mass=%g", P[i].ID, P[i].Mass);

          P[i].Mass           = 0;
          SphP[i].Energy      = 0;
          SphP[i].Momentum[0] = 0;
          SphP[i].Momentum[1] = 0;
          SphP[i].Momentum[2] = 0;
        }
    }
#endif /* #ifdef MESHRELAX */

#ifdef GODUNOV_STATS
  endrun();
#endif /* #ifdef GODUNOV_STATS */

#ifdef ONEDIMS_SPHERICAL
  apply_spherical_source_terms();
#endif /* #ifdef ONEDIMS_SPHERICAL */

#if defined(MHD_POWELL) && defined(VERBOSE)
  double Powell_Momentum[3];
  double Powell_Angular_Momentum[3];
  double Powell_Energy;

  MPI_Reduce(All.Powell_Momentum, Powell_Momentum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(All.Powell_Angular_Momentum, Powell_Angular_Momentum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&All.Powell_Energy, &Powell_Energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    printf("MHD_POWELL: Total ST contribution: Mom=%g,%g,%g   AngMom=%g,%g,%g   Energy=%g\n", Powell_Momentum[0], Powell_Momentum[1],
           Powell_Momentum[2], Powell_Angular_Momentum[0], Powell_Angular_Momentum[1], Powell_Angular_Momentum[2], Powell_Energy);
#endif /* #if defined(MHD_POWELL) && defined(VERBOSE) */

#ifdef VORONOI_BACKUP_RESTORE_FACE_AREAS
  restore_face_areas(T);
#endif /* #ifdef VORONOI_BACKUP_RESTORE_FACE_AREAS */

  TIMER_STOP(CPU_FLUXES);
}

#ifdef VORONOI_BACKUP_RESTORE_FACE_AREAS
/*! \brief Writes face areas to a backup variable.
 *
 *  \param[in, out] T Pointer to tessellation.
 *
 *  \return void
 */
void backup_face_areas(tessellation *T)
{
  for(int i = 0; i < T->Nvf; i++)
    T->VF[i].area_backup = T->VF[i].area;
}

/*! \brief Restores face areas from a backup variable.
 *
 *  \param[in, out] T Pointer to tessellation.
 *
 *  \return void
 */
void restore_face_areas(tessellation *T)
{
  for(int i = 0; i < T->Nvf; i++)
    T->VF[i].area = T->VF[i].area_backup;
}
#endif /* #ifdef VORONOI_BACKUP_RESTORE_FACE_AREAS */

/*! \brief Gets value of hydrodynamial quantities at face.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] p Index in DP array.
 *  \param[in] i Index in VF array.
 *  \param[out] st State at face.
 *
 *  \return 0
 */
int face_get_state(tessellation *T, int p, int i, struct state *st)
{
  int particle;
#if defined(MAXSCALARS)
  int j;
#endif /* #if defined(MAXSCALARS) */
  double aBegin;

  point *DP = T->DP;
  face *VF  = T->VF;

  particle = DP[p].index;

  if(particle < 0)
    return -1;

  if(particle >= NumGas && DP[p].task == ThisTask)
    particle -= NumGas;

  /* interpolation vector for the left state */
  if(DP[p].task == ThisTask)
    {
      st->dx = VF[i].cx - SphP[particle].Center[0];
      st->dy = VF[i].cy - SphP[particle].Center[1];
      st->dz = VF[i].cz - SphP[particle].Center[2];
    }
  else
    {
      st->dx = VF[i].cx - PrimExch[particle].Center[0];
      st->dy = VF[i].cy - PrimExch[particle].Center[1];
      st->dz = VF[i].cz - PrimExch[particle].Center[2];
    }

    /* correct for periodicity */
#if !defined(REFLECTIVE_X) && !defined(ONEDIMS_SPHERICAL)
  if(st->dx < -boxHalf_X)
    st->dx += boxSize_X;
  if(st->dx > boxHalf_X)
    st->dx -= boxSize_X;
#endif /* #if !defined(REFLECTIVE_X) && !defined(ONEDIMS_SPHERICAL) */
#if !defined(REFLECTIVE_Y)
  if(st->dy < -boxHalf_Y)
    st->dy += boxSize_Y;
  if(st->dy > boxHalf_Y)
    st->dy -= boxSize_Y;
#endif /* #if !defined(REFLECTIVE_Y) */
#if !defined(REFLECTIVE_Z)
  if(st->dz < -boxHalf_Z)
    st->dz += boxSize_Z;
  if(st->dz > boxHalf_Z)
    st->dz -= boxSize_Z;
#endif /* #if !defined(REFLECTIVE_Z) */

#ifdef ONEDIMS_SPHERICAL
  if(DP[p].task == ThisTask)
    st->radius = SphP[particle].Center[0];
  else
    st->radius = PrimExch[particle].Center[0];
#endif /* #ifdef ONEDIMS_SPHERICAL */

  if(DP[p].task == ThisTask)
    {
      st->velGas[0] = P[particle].Vel[0];
      st->velGas[1] = P[particle].Vel[1];
      st->velGas[2] = P[particle].Vel[2];

      st->velVertex[0] = SphP[particle].VelVertex[0];
      st->velVertex[1] = SphP[particle].VelVertex[1];
      st->velVertex[2] = SphP[particle].VelVertex[2];

      st->rho = SphP[particle].Density;

      st->press = SphP[particle].Pressure;

      st->grad = &SphP[particle].Grad;

      st->timeBin = P[particle].TimeBinHydro;

      st->volume = SphP[particle].Volume;

#ifdef MHD
      st->Bx = SphP[particle].B[0];
      st->By = SphP[particle].B[1];
      st->Bz = SphP[particle].B[2];
#ifdef MHD_POWELL
      st->divB = SphP[particle].DivB;
#endif /* #ifdef MHD_POWELL */
#endif /* #ifdef MHD */

#ifdef MAXSCALARS
      for(j = 0; j < N_Scalar; j++)
        st->scalars[j] = *(MyFloat *)(((char *)(&SphP[particle])) + scalar_elements[j].offset);
#endif /* #ifdef MAXSCALARS */

      aBegin = SphP[particle].TimeLastPrimUpdate;

      st->oldmass     = SphP[particle].OldMass;
      st->surfacearea = SphP[particle].SurfaceArea;
      st->activearea  = SphP[particle].ActiveArea;
      st->csnd        = get_sound_speed(particle);
      st->ID          = P[particle].ID;
    }
  else
    {
      st->velGas[0] = PrimExch[particle].VelGas[0];
      st->velGas[1] = PrimExch[particle].VelGas[1];
      st->velGas[2] = PrimExch[particle].VelGas[2];

      st->velVertex[0] = PrimExch[particle].VelVertex[0];
      st->velVertex[1] = PrimExch[particle].VelVertex[1];
      st->velVertex[2] = PrimExch[particle].VelVertex[2];

      st->rho = PrimExch[particle].Density;

      st->press = PrimExch[particle].Pressure;

      st->grad = &GradExch[particle];

      st->timeBin = PrimExch[particle].TimeBinHydro; /* This is the hydro timestep */

      st->volume = PrimExch[particle].Volume;

#ifdef MHD
      st->Bx = PrimExch[particle].B[0];
      st->By = PrimExch[particle].B[1];
      st->Bz = PrimExch[particle].B[2];
#ifdef MHD_POWELL
      st->divB = PrimExch[particle].DivB;
#endif /* #ifdef MHD_POWELL */
#endif /* #ifdef MHD */

#ifdef MAXSCALARS
      for(j = 0; j < N_Scalar; j++)
        st->scalars[j] = PrimExch[particle].Scalars[j];
#endif /* #ifdef MAXSCALARS */

      aBegin = PrimExch[particle].TimeLastPrimUpdate;

      st->oldmass     = PrimExch[particle].OldMass;
      st->surfacearea = PrimExch[particle].SurfaceArea;
      st->activearea  = PrimExch[particle].ActiveArea;
      st->csnd        = PrimExch[particle].Csnd;
      st->ID          = DP[p].ID;
    }

  st->dtExtrapolation = All.Time - aBegin;

  /* check for reflecting or outflowing boundaries */
  face_boundary_check_vertex(T, p, &st->velVertex[0], &st->velVertex[1], &st->velVertex[2]);

  return 0;
}

/*! \brief Checks for boundary cells with non-periodic boundary conditions.
 *
 *  Adjusts the velocities accordingly.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] p Index in DP array.
 *  \param[in, out] velx Velocity in x coordinate.
 *  \param[in, out] vely Velocity in y coordinate.
 *  \param[in, out] velz Velocity in z coordinate.
 *
 *  \return void
 */
void face_boundary_check_vertex(tessellation *T, int p, MyFloat *velx, MyFloat *vely, MyFloat *velz)
{
  /* check for reflecting or outflowing boundaries */
#if defined(REFLECTIVE_X)
  if((T->DP[p].image_flags & REFL_X_FLAGS))
    *velx *= -1;
#endif /* #if defined(REFLECTIVE_X) */
#if defined(REFLECTIVE_Y)
  if((T->DP[p].image_flags & REFL_Y_FLAGS))
    *vely *= -1;
#endif /* #if defined(REFLECTIVE_Y) */
#if defined(REFLECTIVE_Z)
  if((T->DP[p].image_flags & REFL_Z_FLAGS))
    *velz *= -1;
#endif /* #if defined(REFLECTIVE_Z) */

#ifdef ONEDIMS_SPHERICAL
  if(p == -1)
    *velx *= -1;
#endif /* #ifdef ONEDIMS_SPHERICAL */
}

/*! \brief Checks for boundary cells with non-periodic boundary conditions.
 *
 *  \param[in] p Pointer to point.
 *  \param[in, out] velx Velocity in x direction.
 *  \param[in, out] vely Velocity in y direction.
 *  \param[in, out] velz Velocity in z direction.
 *
 *  \return void
 */
void face_boundary_check(point *p, double *velx, double *vely, double *velz)
{
  /* check for reflecting or outflowing boundaries */
#if defined(REFLECTIVE_X)
  if((p->image_flags & REFL_X_FLAGS) && !(p->image_flags & OUTFLOW_X))
    *velx *= -1;
#endif /* #if defined(REFLECTIVE_X) */
#if defined(REFLECTIVE_Y)
  if((p->image_flags & REFL_Y_FLAGS) && !(p->image_flags & OUTFLOW_Y))
    *vely *= -1;
#endif /* #if defined(REFLECTIVE_Y) */
#if defined(REFLECTIVE_Z)
  if((p->image_flags & REFL_Z_FLAGS) && !(p->image_flags & OUTFLOW_Z))
    *velz *= -1;
#endif /* #if defined(REFLECTIVE_Z) */

#ifdef ONEDIMS_SPHERICAL
  if(p == &Mesh.DP[-1])
    *velx *= -1;
#endif /* #ifdef ONEDIMS_SPHERICAL */
}

/*! \brief Checks whether local task is responsible for a face.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] p1 Index in DP array of point1 making up the face.
 *  \param[in] p2 Index in DP array of point2 making up the face.
 *  \param[in] st_L Left hand side state of the face.
 *  \param[in] st_R Right hand side state of the face.
 *
 *  \return -1 if not local responsibility, 0 if it is.
 */
int face_check_responsibility_of_this_task(tessellation *T, int p1, int p2, struct state *st_L, struct state *st_R)
{
  int low_p, high_p;
  struct state *low_state, *high_state;

  point *DP = T->DP;

  if(DP[p1].ID < DP[p2].ID)
    {
      low_p      = p1;
      high_p     = p2;
      low_state  = st_L;
      high_state = st_R;
    }
  else if(DP[p1].ID > DP[p2].ID)
    {
      low_p      = p2;
      high_p     = p1;
      low_state  = st_R;
      high_state = st_L;
    }
  else
    {
      /* equality of the IDs should only occur for reflective boundaries */
      if(DP[p1].task == ThisTask && DP[p1].index < NumGas)
        {
          low_p      = p1;
          high_p     = p2;
          low_state  = st_L;
          high_state = st_R;
        }
      else
        {
          low_p      = p2;
          high_p     = p1;
          low_state  = st_R;
          high_state = st_L;
        }
    }

  if(TimeBinSynchronized[low_state->timeBin]) /* the one with the lower ID is active */
    {
      /* we need to check whether the one with the lower ID is a local particle */
      if(DP[low_p].task == ThisTask && DP[low_p].index < NumGas)
        return 0;
    }
  else if(TimeBinSynchronized[high_state->timeBin]) /* only the side with the higher ID is active */
    {
      /* we need to check whether we hold the one with the higher ID, if yes, we'll do it */
      if(DP[high_p].task == ThisTask && DP[high_p].index < NumGas)
        return 0;
    }

  return -1; /* we can skip this face on the local task */
}

/*! \brief Determines timestep of face.
 *
 *  \param[in] state_L Left hand side state of face.
 *  \param[in] state_R Right hand side state of face.
 *  \param[out] hubble_a Value of Hubble function at scalefactor
 *              a(cosmological).
 *  \param[out] atime Scalefactor (cosmological).
 *
 *  \return Face timestep.
 */
double face_timestep(struct state *state_L, struct state *state_R, double *hubble_a, double *atime)
{
  integertime ti_begin_L, ti_begin_R;
  short int timeBin;
  double face_dt;

  /* determine most recent start of the time bins */
  ti_begin_L = (All.Ti_Current >> state_L->timeBin) << state_L->timeBin;
  ti_begin_R = (All.Ti_Current >> state_R->timeBin) << state_R->timeBin;

  /* take the minimum of the two */
  timeBin = state_L->timeBin;
  if(timeBin > state_R->timeBin)
    timeBin = state_R->timeBin;

  /* compute the half-step prediction times */
  state_L->dt_half = (All.Ti_Current + (((integertime)1) << (timeBin - 1)) - ti_begin_L) * All.Timebase_interval;
  state_R->dt_half = (All.Ti_Current + (((integertime)1) << (timeBin - 1)) - ti_begin_R) * All.Timebase_interval;

  if(All.ComovingIntegrationOn)
    {
      /* calculate scale factor at middle of timestep */
      *atime    = All.TimeBegin * exp((All.Ti_Current + (((integertime)1) << (timeBin - 1))) * All.Timebase_interval);
      *hubble_a = hubble_function(*atime);
    }
  else
    *atime = *hubble_a = 1.0;

  /* set the actual time-step for the face */
  face_dt = (((integertime)1) << timeBin) * All.Timebase_interval;

  if(All.ComovingIntegrationOn)
    {
      /* converts to delta_t */
      state_L->dt_half /= *hubble_a;
      state_R->dt_half /= *hubble_a;
      face_dt /= *hubble_a;

      face_dt /= *atime; /* we need dt/a, the (1/a) takes care of the gradient in the cosmological euler equations */

      state_L->dtExtrapolation /= *hubble_a;
      state_L->dtExtrapolation /= *atime;
      state_R->dtExtrapolation /= *hubble_a;
      state_R->dtExtrapolation /= *atime;
    }

  return face_dt;
}

/*! \brief Converts the velocities to local frame, compensating for the
 *         movement of the face.
 *
 *  \param[in, out] st State to be converted to local frame.
 *  \param[in] vel_face Face velocity.
 *  \param[in] hubble_a Value of Hubble function at scalefactor
 *             a (cosmological).
 *  \param[in] atime Scalefactor (cosmological).
 *
 *  \return void
 */
void state_convert_to_local_frame(struct state *st, double *vel_face, double hubble_a, double atime)
{
  if(All.ComovingIntegrationOn)
    {
      st->velGas[0] /= atime; /* convert to peculiar velocity */
      st->velGas[1] /= atime;
      st->velGas[2] /= atime;
    }

  st->velx = st->velGas[0] - vel_face[0];
  st->vely = st->velGas[1] - vel_face[1];
  st->velz = st->velGas[2] - vel_face[2];

  if(All.ComovingIntegrationOn)
    {
      st->velx -= atime * hubble_a * st->dx; /* need to get the physical velocity relative to the face */
      st->vely -= atime * hubble_a * st->dy;
      st->velz -= atime * hubble_a * st->dz;
    }
}

/*! \brief Extrapolates the state in time.
 *
 *  \param[out] delta Change due to time extrapolation.
 *  \param[in] st State to be extrapolated.
 *  \param[in] atime Scalefactor at this time (cosmological).
 *
 *  \return void
 */
void face_do_time_extrapolation(struct state *delta, struct state *st, double atime)
{
  /* st is the state at the center of the cell */

  /* the code still allows for emtpy cells but we are going to divide
   * by rho, so ...
   */
  if(st->rho <= 0)
    return;

#if defined(MESHRELAX) || defined(DISABLE_TIME_EXTRAPOLATION)
  /* do not time extrapolation */
  (void)st;
  (void)atime;
  memset(delta, 0, sizeof(struct state));
  return;
#endif /* #if defined (MESHRELAX) || defined (DISABLE_TIME_EXTRAPOLATION) */

  struct grad_data *grad = st->grad;

  double dt_half = st->dtExtrapolation;

  if(All.ComovingIntegrationOn)
    dt_half /= atime;

  delta->rho = -dt_half * (st->velx * grad->drho[0] + st->rho * grad->dvel[0][0] + st->vely * grad->drho[1] +
                           st->rho * grad->dvel[1][1] + st->velz * grad->drho[2] + st->rho * grad->dvel[2][2]);

  delta->velx = -dt_half * (1.0 / st->rho * grad->dpress[0] + st->velx * grad->dvel[0][0] + st->vely * grad->dvel[0][1] +
                            st->velz * grad->dvel[0][2]);

  delta->vely = -dt_half * (1.0 / st->rho * grad->dpress[1] + st->velx * grad->dvel[1][0] + st->vely * grad->dvel[1][1] +
                            st->velz * grad->dvel[1][2]);

  delta->velz = -dt_half * (1.0 / st->rho * grad->dpress[2] + st->velx * grad->dvel[2][0] + st->vely * grad->dvel[2][1] +
                            st->velz * grad->dvel[2][2]);

  delta->press = -dt_half * (GAMMA * st->press * (grad->dvel[0][0] + grad->dvel[1][1] + grad->dvel[2][2]) +
                             st->velx * grad->dpress[0] + st->vely * grad->dpress[1] + st->velz * grad->dpress[2]);

#ifdef ONEDIMS_SPHERICAL
  delta->velx += dt_half * 2. * st->press / (st->rho * st->radius);
#endif /* #ifdef ONEDIMS_SPHERICAL */

#ifdef MHD
  delta->velx +=
      -dt_half * (1.0 / st->rho *
                  (st->By * grad->dB[1][0] + st->Bz * grad->dB[2][0] - st->By * grad->dB[0][1] - st->Bz * grad->dB[0][2]) / atime);

  delta->vely +=
      -dt_half * (1.0 / st->rho *
                  (st->Bx * grad->dB[0][1] + st->Bz * grad->dB[2][1] - st->Bx * grad->dB[1][0] - st->Bz * grad->dB[1][2]) / atime);

  delta->velz +=
      -dt_half * (1.0 / st->rho *
                  (st->Bx * grad->dB[0][2] + st->By * grad->dB[1][2] - st->Bx * grad->dB[2][0] - st->By * grad->dB[2][1]) / atime);

  delta->Bx =
      -dt_half * (-st->velx * grad->dB[1][1] - grad->dvel[0][1] * st->By + st->vely * grad->dB[0][1] + grad->dvel[1][1] * st->Bx +
                  st->velz * grad->dB[0][2] + grad->dvel[2][2] * st->Bx - st->velx * grad->dB[2][2] - grad->dvel[0][2] * st->Bz);

  delta->By =
      -dt_half * (+st->velx * grad->dB[1][0] + grad->dvel[0][0] * st->By - st->vely * grad->dB[0][0] - grad->dvel[1][0] * st->Bx -
                  st->vely * grad->dB[2][2] - grad->dvel[1][2] * st->Bz + st->velz * grad->dB[1][2] + grad->dvel[2][2] * st->By);

  delta->Bz =
      -dt_half * (-st->velz * grad->dB[0][0] - grad->dvel[2][0] * st->Bx + st->velx * grad->dB[2][0] + grad->dvel[0][0] * st->Bz +
                  st->vely * grad->dB[2][1] + grad->dvel[1][1] * st->Bz - st->velz * grad->dB[1][1] - grad->dvel[2][1] * st->By);
#endif /* #ifdef MHD */

#if defined(MAXSCALARS)
  int k;
  for(k = 0; k < N_Scalar; k++)
    {
      delta->scalars[k] =
          -dt_half * (st->velx * grad->dscalars[k][0] + st->vely * grad->dscalars[k][1] + st->velz * grad->dscalars[k][2]);
    }
#endif /* #if defined(MAXSCALARS) */
}

/*! \brief Extrapolates the state in space.
 *
 *  Linear extrapolation with neighbor cell to their common face.
 *
 *  \param[out] delta Change due to time extrapolation.
 *  \param[in] st State to be extrapolated.
 *  \param[in] st_other state of other cell.
 *
 *  \return void
 */
void face_do_spatial_extrapolation(struct state *delta, struct state *st, struct state *st_other)
{
#ifdef DISABLE_SPATIAL_RECONSTRUCTION
  memset(delta, 0, sizeof(struct state));
  return;
#endif /* #ifdef DISABLE_SPATIAL_RECONSTRUCTION */

#ifdef NO_RECONSTRUCTION_AT_STRONG_SHOCKS
  if(dmax(st->press, st_other->press) > 100. * dmin(st->press, st_other->press))
    {
      memset(delta, 0, sizeof(struct state));
      return;
    }
#endif /* #ifdef NO_RECONSTRUCTION_AT_STRONG_SHOCKS */

  struct grad_data *grad = st->grad;

  double dx[3];
  dx[0] = st->dx;
  dx[1] = st->dy;
  dx[2] = st->dz;

  double r[3];
  r[0] = -st_other->dx + st->dx;
  r[1] = -st_other->dy + st->dy;
  r[2] = -st_other->dz + st->dz;

  face_do_spatial_extrapolation_single_quantity(&delta->rho, st->rho, st_other->rho, grad->drho, dx, r);

  face_do_spatial_extrapolation_single_quantity(&delta->velx, st->velx, st_other->velx, grad->dvel[0], dx, r);
  face_do_spatial_extrapolation_single_quantity(&delta->vely, st->vely, st_other->vely, grad->dvel[1], dx, r);
  face_do_spatial_extrapolation_single_quantity(&delta->velz, st->velz, st_other->velz, grad->dvel[2], dx, r);

  face_do_spatial_extrapolation_single_quantity(&delta->press, st->press, st_other->press, grad->dpress, dx, r);

#ifdef MHD
  face_do_spatial_extrapolation_single_quantity(&delta->Bx, st->Bx, st_other->Bx, grad->dB[0], dx, r);
  face_do_spatial_extrapolation_single_quantity(&delta->By, st->By, st_other->By, grad->dB[1], dx, r);
  face_do_spatial_extrapolation_single_quantity(&delta->Bz, st->Bz, st_other->Bz, grad->dB[2], dx, r);
#endif /* #ifdef MHD */

#ifdef MAXSCALARS
  int k;
  for(k = 0; k < N_Scalar; k++)
    {
      face_do_spatial_extrapolation_single_quantity(&delta->scalars[k], st->scalars[k], st_other->scalars[k], grad->dscalars[k], dx,
                                                    r);
    }
#endif /* #ifdef MAXSCALARS */
}

/*! \brief Extrapolates a single quantity in space.
 *
 *  Linear interpolation with neighbor cell to their common face.
 *
 *  \param[out] delta Change due to time extrapolation.
 *  \param[in] st State to be extrapolated (unused).
 *  \param[in] st_other state of other cell (unused).
 *  \param[in] grad Gradient used for extrapolation.
 *  \param[in] dx normal vector.
 *  \param[in] r (unused).
 *
 *  \return void
 */
void face_do_spatial_extrapolation_single_quantity(double *delta, double st, double st_other, MySingle *grad, double *dx, double *r)
{
  (void)st;
  (void)st_other;
  (void)r;
  *delta = grad[0] * dx[0] + grad[1] * dx[1] + grad[2] * dx[2];
}

/*! \brief Adds space and time extrapolation to state.
 *
 *  \param[in, out] st_face State that is modified.
 *  \param[in] delta_time Change of state due to time extrapolation.
 *  \param[in] delta_space Change of state due to space extrapolation.
 *  \param[in, out] stat Structure that counts face value statistics.
 *
 *  \return void
 */
void face_add_extrapolations(struct state *st_face, struct state *delta_time, struct state *delta_space, struct fvs_stat *stat)
{
  stat->count_disable_extrapolation += 1;

  if(st_face->rho <= 0)
    return;

  if(st_face->rho + delta_time->rho + delta_space->rho < 0 || st_face->press + delta_time->press + delta_space->press < 0)
    return;

  stat->count_disable_extrapolation -= 1;

#if !defined(MESHRELAX) && !defined(DISABLE_TIME_EXTRAPOLATION)
  face_add_extrapolation(st_face, delta_time, stat);
#endif /* #if !defined(MESHRELAX) && !defined(DISABLE_TIME_EXTRAPOLATION)  */

#if !defined(DISABLE_SPATIAL_EXTRAPOLATION)
  face_add_extrapolation(st_face, delta_space, stat);
#endif /* #if !defined(DISABLE_SPATIAL_EXTRAPOLATION) */
}

/*! \brief Adds an extrapolation to state.
 *
 *  Called in face_add_extrapolations(..).
 *
 *  \param[in, out] st_face State that is modified.
 *  \param[in] delta Change of state due to extrapolation.
 *  \param[in] stat (unused)
 *
 *  \return void
 */
void face_add_extrapolation(struct state *st_face, struct state *delta, struct fvs_stat *stat)
{
  st_face->rho += delta->rho;
  st_face->velx += delta->velx;
  st_face->vely += delta->vely;
  st_face->velz += delta->velz;
  st_face->press += delta->press;

#ifdef MHD
#ifndef ONEDIMS
  /* in one dimension, Bx has to be constant! */
  st_face->Bx += delta->Bx;
#endif /* #ifndef ONEDIMS */
  st_face->By += delta->By;
  st_face->Bz += delta->Bz;
#endif /* #ifdef MHD */

#ifdef MAXSCALARS
  int k;
  for(k = 0; k < N_Scalar; k++)
    st_face->scalars[k] += delta->scalars[k];
#endif /* #ifdef MAXSCALARS */
}

/*! \brief Adds an extrapolation to state.
 *
 *  But checks for positivity of density.
 *
 *  \param[in, out] st_face State that is modified.
 *  \param[in] delta Change of state due to extrapolation.
 *  \param[in, out] stat Structure that counts face value statistics.
 *
 *  \return void
 */
void face_add_extrapolation_with_check(struct state *st_face, struct state *delta, struct fvs_stat *stat)
{
  stat->count_disable_extrapolation += 1;

  if(st_face->rho <= 0)
    return;

  if(st_face->rho + delta->rho < 0 || st_face->press + delta->press < 0)
    return;

  stat->count_disable_extrapolation -= 1;

  face_add_extrapolation(st_face, delta, stat);
}

/*! \brief Rotates velocities and magnetic field.
 *
 *  \param[in, out] st State that containes velocities to be rotated.
 *  \param[in] geom Geometry with a rotation matrix.
 *
 *  \return void
 */
void face_turn_velocities(struct state *st, struct geometry *geom)
{
  double velx, vely, velz;

  velx = st->velx;
  vely = st->vely;
  velz = st->velz;

  st->velx = velx * geom->nx + vely * geom->ny + velz * geom->nz;
  st->vely = velx * geom->mx + vely * geom->my + velz * geom->mz;
  st->velz = velx * geom->px + vely * geom->py + velz * geom->pz;

#ifdef MHD
  double Bx, By, Bz;

  Bx = st->Bx;
  By = st->By;
  Bz = st->Bz;

  st->Bx = Bx * geom->nx + By * geom->ny + Bz * geom->nz;
  st->By = Bx * geom->mx + By * geom->my + Bz * geom->mz;
  st->Bz = Bx * geom->px + By * geom->py + Bz * geom->pz;
#endif /* #ifdef MHD */
}

/*! \brief Sets the state at the face to its upwind value.
 *
 *  \param[in] st_L Left hand side hydrodynamical state.
 *  \param[in] st_R Right hand side hydrodynamical state.
 *  \param[out] st_face State at face.
 *  \param[in] geom Geometry structure that includes normal vector of face.
 *  \param[in] vel_face Velocity vector of face.
 *
 *  \return void
 */
void solve_advection(struct state *st_L, struct state *st_R, struct state_face *st_face, struct geometry *geom, double *vel_face)
{
  double ev = vel_face[0] * geom->nx + vel_face[1] * geom->ny + vel_face[2] * geom->nz;

  if(ev < 0)
    {
      st_face->rho   = st_L->rho;
      st_face->velx  = st_L->velx;
      st_face->vely  = st_L->vely;
      st_face->velz  = st_L->velz;
      st_face->press = st_L->press;
    }
  else
    {
      st_face->rho   = st_R->rho;
      st_face->velx  = st_R->velx;
      st_face->vely  = st_R->vely;
      st_face->velz  = st_R->velz;
      st_face->press = st_R->press;
    }
}

/*! \brief Rotates velocities backwards.
 *
 *  Inverse operation to face_turn_velocities(...).
 *
 *  \param[in, out] st State that containes velocities to be rotated.
 *  \param[in] geom Geometry with a rotation matrix.
 *
 *  \return void
 */
void face_turnback_velocities(struct state_face *st_face, struct geometry *geom)
{
  double velx, vely, velz;

  velx = st_face->velx;
  vely = st_face->vely;
  velz = st_face->velz;

  st_face->velx = velx * geom->nx + vely * geom->mx + velz * geom->px;
  st_face->vely = velx * geom->ny + vely * geom->my + velz * geom->py;
  st_face->velz = velx * geom->nz + vely * geom->mz + velz * geom->pz;
}

/*! \brief Sets the scalar states compute the scalar flux from mass flux.
 *
 *  \param[in] st_L Left hand side state.
 *  \param[in] st_R Right hand side state.
 *  \param[out] st_face Face state.
 *  \param[out] flux Flux over face.
 *
 *  \return void
 */
void face_set_scalar_states_and_fluxes(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux)
{
#if defined(MAXSCALARS)
  int i;

  double normfac, normifac;

  if(flux->mass > 0)
    st_face->scalars = st_L->scalars;
  else
    st_face->scalars = st_R->scalars;

  /* Normalize species here */
  normfac = 0;

  for(i = 0; i < N_Scalar; i++)
    {
      flux->scalars[i] = st_face->scalars[i] * flux->mass;

      if(scalar_elements[i].type == SCALAR_TYPE_SPECIES)
        normfac += st_face->scalars[i];
    }

  if(normfac != 0)
    {
      normifac = 1.0 / normfac;

      for(i = 0; i < N_Scalar; i++)
        if(scalar_elements[i].type == SCALAR_TYPE_SPECIES || scalar_elements[i].type == SCALAR_TYPE_NORMALIZE)
          flux->scalars[i] *= normifac;
    }

#endif /* #if defined(MAXSCALARS) */
}

#if defined(RIEMANN_HLLC) || defined(RIEMANN_HLLD)
/*! \brief Converts flux from face frame to simulation box frame.
 *
 *  \param[in] st_L Left hand side state.
 *  \param[in] st_R Right hand side state.
 *  \param[in] vel_face Velocity vector of face.
 *  \param[in, out] flux Flux vector accross face.
 *
 *  \return void
 */
void flux_convert_to_lab_frame(struct state *st_L, struct state *st_R, double *vel_face, struct fluxes *flux)
{
  double momx = flux->momentum[0];
  double momy = flux->momentum[1];
  double momz = flux->momentum[2];

  flux->momentum[0] += vel_face[0] * flux->mass;
  flux->momentum[1] += vel_face[1] * flux->mass;
  flux->momentum[2] += vel_face[2] * flux->mass;

  flux->energy += momx * vel_face[0] + momy * vel_face[1] + momz * vel_face[2] +
                  0.5 * flux->mass * (vel_face[0] * vel_face[0] + vel_face[1] * vel_face[1] + vel_face[2] * vel_face[2]);

#ifdef MHD
  double Bx;
  Bx = 0.5 * (st_L->Bx + st_R->Bx);

  flux->B[0] -= vel_face[0] * Bx;
  flux->B[1] -= vel_face[1] * Bx;
  flux->B[2] -= vel_face[2] * Bx;
#endif /* #ifdef MHD */
}
#endif /* #if defined(RIEMANN_HLLC) || defined(RIEMANN_HLLD) */

/*! \brief Rotates momenum flux and magnetic flux vector.
 *
 *  flux->momentum vector needs to be turned in case the HLLC or Rosunov
 *  Riemann solvers are used.
 *
 *  \param[in, out] flux Flux vector which is rotated.
 *  \param[in] geom Geometry structure that holds rotation matrix.
 *
 *  \return void
 */
void face_turn_momentum_flux(struct fluxes *flux, struct geometry *geom)
{
  double momx = flux->momentum[0];
  double momy = flux->momentum[1];
  double momz = flux->momentum[2];

  flux->momentum[0] = momx * geom->nx + momy * geom->mx + momz * geom->px;
  flux->momentum[1] = momx * geom->ny + momy * geom->my + momz * geom->py;
  flux->momentum[2] = momx * geom->nz + momy * geom->mz + momz * geom->pz;

#ifdef MHD
  double Bx = flux->B[0];
  double By = flux->B[1];
  double Bz = flux->B[2];

  flux->B[0] = Bx * geom->nx + By * geom->mx + Bz * geom->px;
  flux->B[1] = Bx * geom->ny + By * geom->my + Bz * geom->py;
  flux->B[2] = Bx * geom->nz + By * geom->mz + Bz * geom->pz;
#endif /* #ifdef MHD */
}

/*! \brief Calculates the flux from face states.
 *
 *  \param[in] st_L (unused)
 *  \param[in] st_R (unused)
 *  \param[in] st_face State at face.
 *  \param[out] flux Flux at face.
 *  \param[in] geom Geometry structure containing normal vector of face.
 *  \param[in] vel_face Velocity vector of face.
 *
 *  \return void
 */
void face_get_fluxes(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux, struct geometry *geom,
                     double *vel_face)
{
  double fac;

  /* calculate fluxes for ordinary Riemann solver */

  fac = (st_face->velx - vel_face[0]) * geom->nx + (st_face->vely - vel_face[1]) * geom->ny + (st_face->velz - vel_face[2]) * geom->nz;

  flux->mass = st_face->rho * fac;

  flux->momentum[0] = (st_face->rho * st_face->velx * fac + st_face->press * geom->nx);
  flux->momentum[1] = (st_face->rho * st_face->vely * fac + st_face->press * geom->ny);
  flux->momentum[2] = (st_face->rho * st_face->velz * fac + st_face->press * geom->nz);

#ifndef ISOTHERM_EQS
  flux->energy =
      (0.5 * st_face->rho * (st_face->velx * st_face->velx + st_face->vely * st_face->vely + st_face->velz * st_face->velz) +
       st_face->press / GAMMA_MINUS1) *
          fac +
      st_face->press * (st_face->velx * geom->nx + st_face->vely * geom->ny + st_face->velz * geom->nz);
#endif /* #ifndef ISOTHERM_EQS */
}

/*! \brief Flux limiter.
 *
 *  Make sure cell cannot loose more mass than it contains...
 *
 *  \param[in] st_L Left hand side hydrodynamical state.
 *  \param[in] st_R Right hand side hydrodynamical state.
 *  \param[in] st_center_L (unused)
 *  \param[in] st_center_R (unused)
 *  \param[in, out] fulx Flux vector.
 *  \param[in] dt Timestep.
 *  \param[in, out] count Number of calls of this function.
 *  \param[in, out] count_reduced Number if flux reductions caused by this
 *                  function.
 *
 *  \return void
 */
void face_limit_fluxes(struct state *st_L, struct state *st_R, struct state *st_center_L, struct state *st_center_R,
                       struct fluxes *flux, double dt, double *count, double *count_reduced)
{
  *count = *count + 1.0;

  /* choose upwind mass to determine a stability bound on the maximum allowed mass exchange,
     (we do this to prevent negative masses under all circumstances) */
  double upwind_mass, upwind_activearea, reduc_fac;
  integertime upwind_timebin, downstream_timebin;

  if(flux->mass > 0)
    {
      upwind_mass        = st_L->oldmass;
      upwind_activearea  = st_L->activearea;
      upwind_timebin     = st_L->timeBin;
      downstream_timebin = st_R->timeBin;
    }
  else
    {
      upwind_mass        = st_R->oldmass;
      upwind_activearea  = st_R->activearea;
      upwind_timebin     = st_R->timeBin;
      downstream_timebin = st_L->timeBin;
    }

  if(upwind_timebin > downstream_timebin)
    dt *= pow(2, upwind_timebin - downstream_timebin);

  if(fabs(flux->mass * dt * upwind_activearea) > 0.9 * upwind_mass)
    {
      reduc_fac = 0.9 * upwind_mass / fabs(flux->mass * dt * upwind_activearea);

      *count_reduced = *count_reduced + 1.0;

      flux->mass *= reduc_fac;
      flux->energy *= reduc_fac;
      flux->momentum[0] *= reduc_fac;
      flux->momentum[1] *= reduc_fac;
      flux->momentum[2] *= reduc_fac;

      /* remark: do not reduce the magnetic field flux, as it is not coupled to the mass flux */
#ifdef MAXSCALARS
      for(int i = 0; i < N_Scalar; i++)
        flux->scalars[i] *= reduc_fac;
#endif /* #ifdef MAXSCALARS */
    }
}

/*! \brief Set flux vector entries to zero.
 *
 *  \param[out] flux Flux vector.
 *
 *  \return void
 */
void face_clear_fluxes(struct fluxes *flux)
{
  flux->mass        = 0;
  flux->momentum[0] = 0;
  flux->momentum[1] = 0;
  flux->momentum[2] = 0;
  flux->energy      = 0;
#ifdef MHD
  flux->B[0] = 0;
  flux->B[1] = 0;
  flux->B[2] = 0;
#endif /* #ifdef MHD */
}

/*! \brief Adds flux due to advection to flux vector.
 *
 *  \param[in] st_face State at face.
 *  \param[in, out] flux Flux vector.
 *  \param[in] geom Geometry structure containing the face normal vector.
 *  \param[in] vel_face Velocity vector of the face.
 *
 *  \return void
 */
void face_add_fluxes_advection(struct state_face *st_face, struct fluxes *flux, struct geometry *geom, double *vel_face)
{
  double fac = -vel_face[0] * geom->nx - vel_face[1] * geom->ny - vel_face[2] * geom->nz;

  flux->mass += st_face->rho * fac;

  flux->momentum[0] += st_face->rho * st_face->velx * fac;
  flux->momentum[1] += st_face->rho * st_face->vely * fac;
  flux->momentum[2] += st_face->rho * st_face->velz * fac;

  flux->energy +=
      0.5 * st_face->rho * fac * (st_face->velx * st_face->velx + st_face->vely * st_face->vely + st_face->velz * st_face->velz) +
      st_face->press / GAMMA_MINUS1 * fac;
}

/*! \brief Compares tasks of flux list data.
 *
 *  Sort kernel for flux list data.
 *
 *  \param[in] a First flux list data object.
 *  \param[in] b Second flux list data object.
 *
 *  \return (-1,0,1) -1 if a->task < b->task.
 */
int flux_list_data_compare(const void *a, const void *b)
{
  if(((struct flux_list_data *)a)->task < (((struct flux_list_data *)b)->task))
    return -1;

  if(((struct flux_list_data *)a)->task > (((struct flux_list_data *)b)->task))
    return +1;

  return 0;
}

/*! \brief Communicates flux list and applies fluxes to conserved hydro
 *         variables.
 *
 *  \return void
 */
void apply_flux_list(void)
{
  int i, j, p, nimport, ngrp, recvTask;
#if defined(MAXSCALARS)
  int k;
#endif /* #if defined(MAXSCALARS) */

  /* now exchange the flux-list and apply it when needed */

  mysort(FluxList, Nflux, sizeof(struct flux_list_data), flux_list_data_compare);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < Nflux; i++)
    Send_count[FluxList[i].task]++;

  if(Send_count[ThisTask] > 0)
    terminate("Send_count[ThisTask]");

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct flux_list_data *FluxListGet = (struct flux_list_data *)mymalloc("FluxListGet", nimport * sizeof(struct flux_list_data));

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&FluxList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct flux_list_data), MPI_BYTE, recvTask,
                           TAG_DENS_A, &FluxListGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct flux_list_data),
                           MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  /* apply the fluxes */

  for(i = 0; i < nimport; i++)
    {
      p = FluxListGet[i].index;

      P[p].Mass += FluxListGet[i].dM;

      SphP[p].Momentum[0] += FluxListGet[i].dP[0];
      SphP[p].Momentum[1] += FluxListGet[i].dP[1];
      SphP[p].Momentum[2] += FluxListGet[i].dP[2];
#ifdef MHD
      SphP[p].BConserved[0] += FluxListGet[i].dB[0];
      SphP[p].BConserved[1] += FluxListGet[i].dB[1];
      SphP[p].BConserved[2] += FluxListGet[i].dB[2];
#endif /* #ifdef MHD */

#ifdef MAXSCALARS
      for(k = 0; k < N_Scalar; k++)
        *(MyFloat *)(((char *)(&SphP[p])) + scalar_elements[k].offset_mass) += FluxListGet[i].dConservedScalars[k];
#endif /* #ifdef MAXSCALARS */

#ifndef ISOTHERM_EQS
      SphP[p].Energy += FluxListGet[i].dEnergy;
#endif /* #ifndef ISOTHERM_EQS */
    }
  myfree(FluxListGet);
}

/*! \brief Initializes statistics of finite volume solver.
 *
 *  \param[out] stat Statistics structure.
 *
 *  \return void
 */
void fvs_initialize_statistics(struct fvs_stat *stat) { stat->count_disable_extrapolation = 0; }

/*! \brief Gathers statistics properties from all tasks and prints information.
 *
 *  \param[in] stat Finite volume solver statistics structure.
 *
 *  \return void
 */
void fvs_evaluate_statistics(struct fvs_stat *stat)
{
#ifdef VERBOSE
  int count_disable_extrapolation = 0;
  MPI_Reduce(&stat->count_disable_extrapolation, &count_disable_extrapolation, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_printf("FLUX: Disabled extrapolation for %d interfaces.\n", count_disable_extrapolation);
#endif /* #ifdef VERBOSE */
}

#ifdef ONEDIMS_SPHERICAL
/*! \brief Applies source terms that occur due to spherical symmetry.
 *
 *  \return void
 */
void apply_spherical_source_terms()
{
  int idx, i;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double Pressure         = SphP[i].Pressure;
      double dt_Extrapolation = All.Time - SphP[i].TimeLastPrimUpdate;
      struct grad_data *grad  = &SphP[i].Grad;

      Pressure += -dt_Extrapolation * (GAMMA * Pressure * (grad->dvel[0][0] + grad->dvel[1][1] + grad->dvel[2][2]) +
                                       P[i].Vel[0] * grad->dpress[0] + P[i].Vel[1] * grad->dpress[1] + P[i].Vel[2] * grad->dpress[2]);

      double dt = 0.5 * (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      SphP[i].Momentum[0] += dt * Pressure * (Mesh.VF[i + 1].area - Mesh.VF[i].area);
    }
}
#endif /* #ifdef ONEDIMS_SPHERICAL */
