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
 * \file        src/gravity/gravtree_forcetest.c
 * \date        05/2018
 * \brief       Test short range gravity evaluation.
 * \details     contains functions:
 *                static void particle2in(data_in * in, int i, int firstnode)
 *                static void out2particle(data_out * out, int i, int mode)
 *                static void kernel_local(void)
 *                static void kernel_imported(void)
 *                void gravity_forcetest(void)
 *                static void gravity_forcetest_evaluate(int target, int mode,
 *                  int threadid)
 *                void gravity_forcetest_testforcelaw(void)
 *                static void ewald_other_images(double x, double y, double z,
 *                  double alpha, double force[4])
 *                static void ewald_correction_force(double x, double y,
 *                  double z, double force[4])
 *                void forcetest_ewald_init(void)
 *                static void ewald_correction_force_table_lookup(double dx,
 *                  double dy, double dz, double force[4])
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 20.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

#ifdef FORCETEST

#if !defined(EVALPOTENTIAL) && defined(FORCETEST)
#error "When you enable FORCETEST you should also switch on EVALPOTENTIAL"
#endif /* #if !defined(EVALPOTENTIAL) && defined(FORCETEST) */

static void gravity_forcetest_evaluate(int target, int mode, int threadid);
static void ewald_correction_force(double x, double y, double z, double force[4]);
static void ewald_other_images(double x, double y, double z, double alpha, double force[4]);
static void ewald_correction_force_table_lookup(double x, double y, double z, double force[4]);

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  unsigned char Type;
  unsigned char SofteningType;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/*! \brief Routine that fills the relevant particle/cell data into the input
 *         structure defined above. Needed by generic_comm_helpers2.
 *
 *  \param[out] in Data structure to fill.
 *  \param[in] i Index of particle in P and SphP arrays.
 *  \param[in] firstnode First note of communication.
 *
 *  \return void
 */
static void particle2in(data_in *in, int i, int firstnode)
{
#ifdef CELL_CENTER_GRAVITY
  if(P[i].Type == 0)
    {
      for(int k = 0; k < 3; k++)
        in->Pos[k] = SphP[i].Center[k];
    }
  else
#endif /* #ifdef CELL_CENTER_GRAVITY */
    {
      for(int k = 0; k < 3; k++)
        in->Pos[k] = P[i].Pos[k];
    }

  in->Type          = P[i].Type;
  in->SofteningType = P[i].SofteningType;

  in->Firstnode = firstnode;
}

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
  MyFloat Acc[3];
  MyFloat Pot;
  MyFloat DistToID1;
#ifdef PMGRID
  MyFloat AccLongRange[3];
  MyFloat AccShortRange[3];
  MyFloat PotLongRange;
  MyFloat PotShortRange;
#endif /* #ifdef PMGRID */
} data_out;

static data_out *DataResult, *DataOut;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *  particle and cell data arrays (P, SphP,...)
 *  \param[in] i Index of particle in P and SphP arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *  communicated from other tasks and has to be added locally?
 *
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      P[i].GravAccelDirect[0] = out->Acc[0];
      P[i].GravAccelDirect[1] = out->Acc[1];
      P[i].GravAccelDirect[2] = out->Acc[2];
      P[i].PotentialDirect    = out->Pot;
      P[i].DistToID1          = out->DistToID1;
#ifdef PMGRID
      P[i].GravAccelLongRange[0]  = out->AccLongRange[0];
      P[i].GravAccelLongRange[1]  = out->AccLongRange[1];
      P[i].GravAccelLongRange[2]  = out->AccLongRange[2];
      P[i].GravAccelShortRange[0] = out->AccShortRange[0];
      P[i].GravAccelShortRange[1] = out->AccShortRange[1];
      P[i].GravAccelShortRange[2] = out->AccShortRange[2];
      P[i].PotentialLongRange     = out->PotLongRange;
      P[i].PotentialShortRange    = out->PotShortRange;
#endif /* #ifdef PMGRID */
    }
  else /* combine */
    {
      P[i].GravAccelDirect[0] += out->Acc[0];
      P[i].GravAccelDirect[1] += out->Acc[1];
      P[i].GravAccelDirect[2] += out->Acc[2];
      P[i].PotentialDirect += out->Pot;
      if(out->DistToID1 > 0)
        P[i].DistToID1 = out->DistToID1;
#ifdef PMGRID
      P[i].GravAccelLongRange[0] += out->AccLongRange[0];
      P[i].GravAccelLongRange[1] += out->AccLongRange[1];
      P[i].GravAccelLongRange[2] += out->AccLongRange[2];
      P[i].GravAccelShortRange[0] += out->AccShortRange[0];
      P[i].GravAccelShortRange[1] += out->AccShortRange[1];
      P[i].GravAccelShortRange[2] += out->AccShortRange[2];
      P[i].PotentialLongRange += out->PotLongRange;
      P[i].PotentialShortRange += out->PotShortRange;
#endif /* #ifdef PMGRID */
    }
}

#include "../utils/generic_comm_helpers2.h"

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int i;

  /* do local particles */
  {
    int j, threadid = get_thread_num();

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

        i = NextParticle++;

        if(i >= TimeBinsGravity.NActiveParticles)
          break;

        i = TimeBinsGravity.ActiveParticleList[i];
        if(i < 0)
          continue;

        if(P[i].TimeBinGrav < 0)
          gravity_forcetest_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

/*! \brief Routine that defines what to do with imported particles.
 *
 *  Calls the *_evaluate function in MODE_IMPORTED_PARTICLES.
 *
 *  \return void
 */
static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
  {
    int threadid = get_thread_num();

    while(1)
      {
        i = cnt++;

        if(i >= Nimport)
          break;

        gravity_forcetest_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*! \brief This function computes the gravitational forces for all active
 *  particles.
 *
 *  A new tree is constructed, if the number of force computations since
 *  it's last construction exceeds some fraction of the total
 *  particle number, otherwise tree nodes are dynamically updated if needed.
 *
 *  \return void
 */
void gravity_forcetest(void)
{
  int nthis, nloc, ntot;
  int idx, i, j;
  double fac1;
  char buf[200];

  nloc = 0;
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(get_random_number() < FORCETEST)
        {
          P[i].TimeBinGrav = -P[i].TimeBinGrav - 1; /* Mark as selected */
          nloc++;
        }
    }

  MPI_Allreduce(&nloc, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("FORCETEST: Testing forces of %d particles\n", ntot);

  double t0 = second();

  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsGravity.NActiveParticles, kernel_local, kernel_imported);

  double t1   = second();
  double maxt = timediff(t0, t1);

  /*  muliply by G */
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].TimeBinGrav < 0)
        {
          for(j = 0; j < 3; j++)
            {
              P[i].GravAccelDirect[j] *= All.G;
#ifdef PMGRID
              P[i].GravAccelLongRange[j] *= All.G;
              P[i].GravAccelShortRange[j] *= All.G;
#endif /* #ifdef PMGRID */
            }

          P[i].PotentialDirect *= All.G;
#ifdef PMGRID
          P[i].PotentialLongRange *= All.G;
          P[i].PotentialShortRange *= All.G;
#endif /* #ifdef PMGRID */
        }
    }

  /* Finally, the following factor allows a computation of cosmological simulation
     with vacuum energy in physical coordinates */

  if(All.ComovingIntegrationOn == 0)
    {
      fac1 = All.OmegaLambda * All.Hubble * All.Hubble;

      for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(P[i].TimeBinGrav < 0)
            for(j = 0; j < 3; j++)
              P[i].GravAccelDirect[j] += fac1 * P[i].Pos[j];
        }
    }

  /* now output the forces to a file */

  for(nthis = 0; nthis < NTask; nthis++)
    {
      if(nthis == ThisTask)
        {
          sprintf(buf, "%s%s", All.OutputDir, "forcetest.txt");

          if(!(FdForceTest = fopen(buf, "a")))
            terminate("error in opening file '%s'\n", buf);

          for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
            {
              i = TimeBinsGravity.ActiveParticleList[idx];
              if(i < 0)
                continue;

              if(P[i].TimeBinGrav < 0)
                {
#ifdef PMGRID
                  fprintf(FdForceTest,
                          "%d %d %lld  %g  %g %g %g  %g  %15.10g %15.10g %15.10g  %15.10g %15.10g %15.10g  %15.10g %15.10g %15.10g  "
                          "%15.10g %15.10g %15.10g  %15.10g %15.10g %15.10g  %15.10g %15.10g %15.10g %15.10g  %15.10g\n",
                          P[i].Type, ThisTask, (long long)P[i].ID, All.Time, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[i].DistToID1,
                          P[i].GravAccelDirect[0], P[i].GravAccelDirect[1], P[i].GravAccelDirect[2], P[i].GravAccelShortRange[0],
                          P[i].GravAccelShortRange[1], P[i].GravAccelShortRange[2], P[i].GravAccelLongRange[0],
                          P[i].GravAccelLongRange[1], P[i].GravAccelLongRange[2], P[i].GravAccel[0], P[i].GravAccel[1],
                          P[i].GravAccel[2], P[i].GravPM[0], P[i].GravPM[1], P[i].GravPM[2], P[i].PotentialDirect,
                          P[i].PotentialShortRange, P[i].PotentialLongRange, P[i].Potential, P[i].PM_Potential);
#else  /* #ifdef PMGRID */
                  fprintf(FdForceTest,
                          "%d %d %lld %g  %g %g %g %g  %15.10g %15.10g %15.10g  %15.10g %15.10g %15.10g  %15.10g %15.10g\n", P[i].Type,
                          ThisTask, (long long)P[i].ID, All.Time, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[i].DistToID1,
                          P[i].GravAccelDirect[0], P[i].GravAccelDirect[1], P[i].GravAccelDirect[2], P[i].GravAccel[0],
                          P[i].GravAccel[1], P[i].GravAccel[2], P[i].PotentialDirect, P[i].Potential);
#endif /* #ifdef PMGRID #else */
                }
            }

          fclose(FdForceTest);
        }

      MPI_Barrier(MPI_COMM_WORLD);
    }

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].TimeBinGrav < 0)
        P[i].TimeBinGrav = -P[i].TimeBinGrav - 1;
    }

  /* Now the force computation is finished */

  if(ThisTask == 0)
    {
      double costtotal = NumPart * ntot;

      fprintf(FdTimings, "DIRECT Nf= %d    part/sec=%g | %g  ia/part=%g\n\n", ntot, ((double)ntot) / (NTask * maxt + 1.0e-20),
              ntot / ((maxt + 1.0e-20) * NTask), ((double)(costtotal)) / (ntot + 1.0e-20));

      myflush(FdTimings);
    }
}

/*! \brief This function does the gravitational force computation with direct
 *  summation for the specified particle.
 *
 *  This can be useful for debugging purposes, in particular for explicit
 *  checks of the force accuracy reached with the tree. Depending on whether
 *  or not a PMGRID is used, the code does a short-range tree-walk or a full
 *  one.
 *
 *  \param i Index of the particle to be processed.
 *  \param mode 0: process local particle (phase 1), 1: process imported
 *         particle (phase 2).
 *  \param thread_id Id of this thread.
 *  \param measure_cost_flag Whether the cost of the tree walk should be
 *         measured.
 *
 *  \return Number of interactions processed for particle i.
 */
static void gravity_forcetest_evaluate(int target, int mode, int threadid)
{
  int j;
  double h_i, h_j, hmax, mass, dx, dy, dz, r, r2, fac, wp, fac_newton, wp_newton;
  double pos_x, pos_y, pos_z;
#ifdef PMGRID
  double asmth = All.Asmth[0];
#endif /* #ifdef PMGRID */
#if !defined(GRAVITY_NOT_PERIODIC)
  double xtmp, ytmp, ztmp;
#endif /* #if !defined(GRAVITY_NOT_PERIODIC) */

  double acc_x     = 0;
  double acc_y     = 0;
  double acc_z     = 0;
  double pot       = 0;
  double disttoid1 = 0;

  data_out out;
  data_in local, *target_data;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      /* make sure that the particle is exported to all other tasks */
      for(int task = 0; task < NTask; task++)
        if(task != ThisTask)
          {
            if(Thread[threadid].Exportflag[task] != target)
              {
                Thread[threadid].Exportflag[task]     = target;
                int nexp                              = Thread[threadid].Nexport++;
                Thread[threadid].PartList[nexp].Task  = task;
                Thread[threadid].PartList[nexp].Index = target;
                Thread[threadid].ExportSpace -= Thread[threadid].ItemSize;
              }

            int nexp = Thread[threadid].NexportNodes++;
            nexp     = -1 - nexp;
            struct datanodelist *nodelist =
                (struct datanodelist *)(((char *)Thread[threadid].PartList) + Thread[threadid].InitialSpace);
            nodelist[nexp].Task  = task;
            nodelist[nexp].Index = target;
            nodelist[nexp].Node  = 0; /* the node doesn't matter here */
            Thread[threadid].ExportSpace -= sizeof(struct datanodelist) + sizeof(int);
          }
    }
  else
    {
      target_data = &DataGet[target];
    }

  pos_x = target_data->Pos[0];
  pos_y = target_data->Pos[1];
  pos_z = target_data->Pos[2];
  h_i   = All.ForceSoftening[target_data->SofteningType];

#ifdef PLACEHIGHRESREGION
  if(pmforce_is_particle_high_res(target_data->Type, target_data->Pos))
    asmth = All.Asmth[1];
#endif /* #ifdef PLACEHIGHRESREGION */

  out.Pot = 0;
#ifdef PMGRID
  out.PotShortRange = 0;
  out.PotLongRange  = 0;
#endif /* #ifdef PMGRID */

  for(int i = 0; i < 3; i++)
    {
      out.Acc[i] = 0;
#ifdef PMGRID
      out.AccShortRange[i] = 0;
      out.AccLongRange[i]  = 0;
#endif /* #ifdef PMGRID */
    }

  for(j = 0; j < NumPart; j++)
    {
      h_j = All.ForceSoftening[P[j].SofteningType];

      if(h_j > h_i)
        hmax = h_j;
      else
        hmax = h_i;

#ifdef CELL_CENTER_GRAVITY
      if(P[j].Type == 0)
        {
          dx = GRAVITY_NEAREST_X(SphP[j].Center[0] - pos_x);
          dy = GRAVITY_NEAREST_Y(SphP[j].Center[1] - pos_y);
          dz = GRAVITY_NEAREST_Z(SphP[j].Center[2] - pos_z);
        }
      else
#endif /* #ifdef CELL_CENTER_GRAVITY */
        {
          dx = GRAVITY_NEAREST_X(P[j].Pos[0] - pos_x);
          dy = GRAVITY_NEAREST_Y(P[j].Pos[1] - pos_y);
          dz = GRAVITY_NEAREST_Z(P[j].Pos[2] - pos_z);
        }

      r2 = dx * dx + dy * dy + dz * dz;

      mass = P[j].Mass;

      /* now evaluate the multipole moment */

      r = sqrt(r2);

      if(P[j].ID == 1)
        disttoid1 = r;

      /* we compute 3 different forces:
       * (1) The correct direct summation force, if needed with Ewald correction: ftrue
       * In the case of PM:
       * (2) The short range direct summation force with only the erfc cut-off (this is what the tree can at best deliver): fsr
       * (3) The expected PM force based on the long-range part of the Ewald sum. This is equal to ftrue - fsr - fsfr_periodic_images
       * */

      if(r > 0)
        {
          fac_newton = mass / (r2 * r);
          wp_newton  = -mass / r;
        }
      else
        {
          fac_newton = 0;
          wp_newton  = 0;
        }

      if(r >= hmax)
        {
          fac = fac_newton;
          wp  = wp_newton;
        }
      else
        {
          double h_inv  = 1.0 / hmax;
          double h3_inv = h_inv * h_inv * h_inv;
          double u      = r * h_inv;

          if(u < 0.5)
            {
              double u2 = u * u;
              fac       = mass * h3_inv * (SOFTFAC1 + u2 * (SOFTFAC2 * u + SOFTFAC3));
              wp        = mass * h_inv * (SOFTFAC4 + u2 * (SOFTFAC5 + u2 * (SOFTFAC6 * u + SOFTFAC7)));
            }
          else
            {
              double u2 = u * u, u3 = u2 * u;
              fac = mass * h3_inv * (SOFTFAC8 + SOFTFAC9 * u + SOFTFAC10 * u2 + SOFTFAC11 * u3 + SOFTFAC12 / u3);
              wp  = mass * h_inv * (SOFTFAC13 + SOFTFAC14 / u + u2 * (SOFTFAC1 + u * (SOFTFAC15 + u * (SOFTFAC16 + SOFTFAC17 * u))));
            }
        }

      double acc_newton_x = dx * fac;
      double acc_newton_y = dy * fac;
      double acc_newton_z = dz * fac;
      double pot_newton   = wp;

#ifdef PMGRID
      double u = 0.5 / asmth * r;

      double factor_force = (erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u) - 1.0);
      double factor_pot   = erfc(u);

      fac += fac_newton * factor_force;
      wp += wp_newton * (factor_pot - 1.0);

      double acc_short_x = dx * fac;
      double acc_short_y = dy * fac;
      double acc_short_z = dz * fac;
      double pot_short   = wp + mass * M_PI / (asmth * asmth * boxSize_X * boxSize_Y * boxSize_Z);

      out.AccShortRange[0] += acc_short_x;
      out.AccShortRange[1] += acc_short_y;
      out.AccShortRange[2] += acc_short_z;
      out.PotShortRange += pot_short;
#endif /* #ifdef PMGRID */

#if defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) && !defined(ONEDIMS_SPHERICAL)
      double fcorr[4];

#if !defined(FORCETEST_TESTFORCELAW)
      ewald_correction_force_table_lookup(dx, dy, dz, fcorr);
#else  /* #if !defined(FORCETEST_TESTFORCELAW) */
      ewald_correction_force(dx, dy, dz, fcorr);
#endif /* #if !defined(FORCETEST_TESTFORCELAW) #else */

      acc_x = acc_newton_x + mass * fcorr[0];
      acc_y = acc_newton_y + mass * fcorr[1];
      acc_z = acc_newton_z + mass * fcorr[2];

      pot = pot_newton + mass * fcorr[3];
#else  /* #if defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) && !defined(ONEDIMS_SPHERICAL) */
      acc_x = acc_newton_x;
      acc_y = acc_newton_y;
      acc_z = acc_newton_z;
      pot = pot_newton;
#endif /* #if defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) && !defined(ONEDIMS_SPHERICAL) #else */

      out.Acc[0] += acc_x;
      out.Acc[1] += acc_y;
      out.Acc[2] += acc_z;
      out.Pot += pot;

#ifdef PMGRID
      double fimages[4] = {0, 0, 0, 0};
#ifdef FORCETEST_TESTFORCELAW
      ewald_other_images(dx, dy, dz, 0.5 / asmth, fimages);
#endif /* #ifdef FORCETEST_TESTFORCELAW */
      out.AccLongRange[0] += acc_x - acc_short_x - mass * fimages[0];
      out.AccLongRange[1] += acc_y - acc_short_y - mass * fimages[1];
      out.AccLongRange[2] += acc_z - acc_short_z - mass * fimages[2];
      out.PotLongRange += pot - pot_short - mass * fimages[3];
#endif /* #ifdef PMGRID */
    }

  out.DistToID1 = disttoid1;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;
}

#ifdef FORCETEST_TESTFORCELAW
/*! \brief Places particle with ID 1 radomly in box and calculates force on it.
 *
 *  \return void
 */
void gravity_forcetest_testforcelaw(void)
{
  int Ncycles = 40;
  double xyz[3], eps;

  ngb_treefree();
  mark_active_timebins();

  for(int cycle = 0; cycle < Ncycles; cycle++)
    {
      mpi_printf("\nTEST-FORCE-LAW: cycle=%d|%d ----------------------------------\n\n", cycle, Ncycles);

      double epsloc = 0, xyzloc[3] = {0, 0, 0};

      /* set particle with ID=1 to new random coordinate in box */
      for(int n = 0; n < NumPart; n++)
        {
          P[n].Type = 1;

          if(P[n].ID == 1)
            {
              xyzloc[0] = All.BoxSize * STRETCHX * get_random_number();
              xyzloc[1] = All.BoxSize * STRETCHY * get_random_number();
              xyzloc[2] = All.BoxSize * STRETCHZ * get_random_number();

#if defined(PLACEHIGHRESREGION) && (FORCETEST_TESTFORCELAW == 1)
              for(int j = 0; j < 3; j++)
                xyzloc[j] = 0.5 * (All.Xmintot[1][j] + All.Xmaxtot[1][j]);
#endif /* #if defined(PLACEHIGHRESREGION) && (FORCETEST_TESTFORCELAW == 1) */

#if defined(PLACEHIGHRESREGION) && (FORCETEST_TESTFORCELAW == 2)
              if(get_random_number() < 0.5)
                {
                  for(int j = 0; j < 3; j++)
                    xyzloc[j] = All.Xmintot[1][j] + get_random_number() * (All.Xmaxtot[1][j] - All.Xmintot[1][j]);
                }
#endif /* #if defined(PLACEHIGHRESREGION) && (FORCETEST_TESTFORCELAW == 2) */

              for(int i = 0; i < 3; i++)
                P[n].Pos[i] = xyzloc[i];

              epsloc = All.ForceSoftening[P[n].SofteningType];
            }
        }

      MPI_Allreduce(xyzloc, xyz, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&epsloc, &eps, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      double rmin = 0.01 * eps;
      double rmax =
          sqrt(pow(0.5 * All.BoxSize * STRETCHX, 2) + pow(0.5 * All.BoxSize * STRETCHY, 2) + pow(0.5 * All.BoxSize * STRETCHZ, 2));

      for(int n = 0; n < NumPart; n++)
        {
          if(P[n].ID != 1)
            {
              double r     = exp(log(rmin) + (log(rmax) - log(rmin)) * get_random_number());
              double theta = acos(2 * get_random_number() - 1);
              double phi   = 2 * M_PI * get_random_number();

              double dx = r * sin(theta) * cos(phi);
              double dy = r * sin(theta) * sin(phi);
              double dz = r * cos(theta);

              double xtmp, ytmp, ztmp;
              P[n].Pos[0] = WRAP_X(xyz[0] + dx);
              P[n].Pos[1] = WRAP_Y(xyz[1] + dy);
              P[n].Pos[2] = WRAP_Z(xyz[2] + dz);
            }
        }

      domain_free();
      domain_Decomposition(); /* do domain decomposition if needed */

#ifdef PMGRID
      long_range_force();
#endif /* #ifdef PMGRID */

      compute_grav_accelerations(All.HighestActiveTimeBin, FLAG_FULL_TREE);
    }

  endrun();
}
#endif /* #ifdef FORCETEST_TESTFORCELAW */

/*! \brief Periodicity effects in gravity.
 *
 *  \param[in] x X coordinate of point.
 *  \param[in] y Y coordinate of point.
 *  \param[in] z Z coordinate of point.
 *  \param[in] alpha Cutoff for tree-PM.
 *  \param[out] force Force vector.
 */
static void ewald_other_images(double x, double y, double z, double alpha, double force[4])
{
  double signx, signy, signz;

  for(int i = 0; i < 4; i++)
    force[i] = 0;

  double r2 = x * x + y * y + z * z;

  if(r2 == 0)
    return;

  if(x < 0)
    {
      x     = -x;
      signx = +1;
    }
  else
    signx = -1;
  if(y < 0)
    {
      y     = -y;
      signy = +1;
    }
  else
    signy = -1;
  if(z < 0)
    {
      z     = -z;
      signz = +1;
    }
  else
    signz = -1;

  double alpha2 = alpha * alpha;

  const int nmax = 4;

  for(int nx = -nmax; nx <= nmax; nx++)
    for(int ny = -nmax; ny <= nmax; ny++)
      for(int nz = -nmax; nz <= nmax; nz++)
        {
          if(nx != 0 || ny != 0 || nz != 0)
            {
              double dx   = x - nx * STRETCHX * All.BoxSize;
              double dy   = y - ny * STRETCHY * All.BoxSize;
              double dz   = z - nz * STRETCHZ * All.BoxSize;
              double r2   = dx * dx + dy * dy + dz * dz;
              double r    = sqrt(r2);
              double val  = erfc(alpha * r) + 2.0 * alpha * r / sqrt(M_PI) * exp(-alpha2 * r2);
              double val2 = val / (r2 * r);
              double val3 = erfc(alpha * r) / r;

              force[0] -= dx * val2;
              force[1] -= dy * val2;
              force[2] -= dz * val2;
              force[3] -= val3;
            }
        }

  force[0] *= signx;
  force[1] *= signy;
  force[2] *= signz;
}

/*! \brief Force due to periodic boundary conditions.
 *
 *  \param[in] x X coordinate of point.
 *  \param[in] y Y coordinate of point.
 *  \param[in] z Z coordinate of point.
 *  \param[out] force Force vector.
 */
static void ewald_correction_force(double x, double y, double z, double force[4])
{
  double signx, signy, signz;

  for(int i = 0; i < 4; i++)
    force[i] = 0;

  double r2 = x * x + y * y + z * z;

  if(r2 == 0)
    return;

  if(x < 0)
    {
      x     = -x;
      signx = +1;
    }
  else
    signx = -1;
  if(y < 0)
    {
      y     = -y;
      signy = +1;
    }
  else
    signy = -1;
  if(z < 0)
    {
      z     = -z;
      signz = +1;
    }
  else
    signz = -1;

  double lmin   = imin(imin(STRETCHX, STRETCHY), STRETCHZ);
  double alpha  = 2.0 / lmin / All.BoxSize;
  double alpha2 = alpha * alpha;
  double r      = sqrt(r2);
  double r3inv  = 1.0 / (r2 * r);

  force[0] += r3inv * x;
  force[1] += r3inv * y;
  force[2] += r3inv * z;

  const int nmax = 6;

  for(int nx = -nmax; nx <= nmax; nx++)
    for(int ny = -nmax; ny <= nmax; ny++)
      for(int nz = -nmax; nz <= nmax; nz++)
        {
          double dx   = x - nx * STRETCHX * All.BoxSize;
          double dy   = y - ny * STRETCHY * All.BoxSize;
          double dz   = z - nz * STRETCHZ * All.BoxSize;
          double r2   = dx * dx + dy * dy + dz * dz;
          double r    = sqrt(r2);
          double val  = erfc(alpha * r) + 2.0 * alpha * r / sqrt(M_PI) * exp(-alpha2 * r2);
          double val2 = val / (r2 * r);
          double val3 = erfc(alpha * r) / r; /* for potential */

          force[0] -= dx * val2;
          force[1] -= dy * val2;
          force[2] -= dz * val2;
          force[3] -= val3;
        }

  int nxmax = (int)(4 * alpha * All.BoxSize * (STRETCHX / lmin) + 0.5);
  int nymax = (int)(4 * alpha * All.BoxSize * (STRETCHY / lmin) + 0.5);
  int nzmax = (int)(4 * alpha * All.BoxSize * (STRETCHZ / lmin) + 0.5);

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      for(int nz = -nzmax; nz <= nzmax; nz++)
        {
          double kx = (2.0 * M_PI / (All.BoxSize * STRETCHX)) * nx;
          double ky = (2.0 * M_PI / (All.BoxSize * STRETCHY)) * ny;
          double kz = (2.0 * M_PI / (All.BoxSize * STRETCHZ)) * nz;
          double k2 = kx * kx + ky * ky + kz * kz;

          if(k2 > 0)
            {
              double kdotx = (x * kx + y * ky + z * kz);
              double vv    = 4.0 * M_PI / (k2 * pow(All.BoxSize, 3) * STRETCHX * STRETCHY * STRETCHZ) * exp(-k2 / (4.0 * alpha2));
              double val   = vv * sin(kdotx);
              double val2  = vv * cos(kdotx);
              force[0] -= kx * val;
              force[1] -= ky * val;
              force[2] -= kz * val;
              force[3] -= val2;
            }
        }

  force[3] += M_PI / (alpha2 * pow(All.BoxSize, 3) * STRETCHX * STRETCHY * STRETCHZ) + 1.0 / r;

  force[0] *= signx;
  force[1] *= signy;
  force[2] *= signz;
}

#if !defined(FORCETEST_TESTFORCELAW)

#define TEW_N 128

#define TEW_NX (DBX * STRETCHX * TEW_N)
#define TEW_NY (DBY * STRETCHY * TEW_N)
#define TEW_NZ (DBZ * STRETCHZ * TEW_N)

static double Ewd_table[4][TEW_NX + 1][TEW_NY + 1][TEW_NZ + 1];
static double Ewd_table_intp;

/*! \brief Initializes Ewald correction force test.
 *
 *  \return void
 */
void forcetest_ewald_init(void)
{
  double t0 = second();

  mpi_printf("FORCETEST: initialize high-res Ewald lookup table...\n");

#ifdef LONG_X
  if(LONG_X != (int)(LONG_X))
    terminate("LONG_X must be an integer");
#endif /* #ifdef LONG_X */

#ifdef LONG_Y
  if(LONG_Y != (int)(LONG_Y))
    terminate("LONG_Y must be an integer");
#endif /* #ifdef LONG_Y */

#ifdef LONG_Z
  if(LONG_Z != (int)(LONG_Z))
    terminate("LONG_Z must be an integer");
#endif /* #ifdef LONG_Z */

  /* ok, let's compute things. Actually, we do that in parallel. */
  int size = (TEW_NX + 1) * (TEW_NY + 1) * (TEW_NZ + 1);
  int first, count;

  subdivide_evenly(size, NTask, ThisTask, &first, &count);

  for(int n = first; n < first + count; n++)
    {
      int i = n / ((TEW_NY + 1) * (TEW_NZ + 1));
      int j = (n - i * (TEW_NY + 1) * (TEW_NZ + 1)) / (TEW_NZ + 1);
      int k = (n - i * (TEW_NY + 1) * (TEW_NZ + 1) - j * (TEW_NZ + 1));

      if(ThisTask == 0)
        {
          if(((n - first) % (count / 20)) == 0)
            {
              printf("%4.1f percent done\n", (n - first) / (count / 100.0));
              myflush(stdout);
            }
        }

      double xx = 0.5 * DBX * STRETCHX * ((double)i) / TEW_NX * All.BoxSize;
      double yy = 0.5 * DBY * STRETCHY * ((double)j) / TEW_NY * All.BoxSize;
      double zz = 0.5 * DBZ * STRETCHZ * ((double)k) / TEW_NZ * All.BoxSize;

      double fcorr[4];
      ewald_correction_force(xx, yy, zz, fcorr);

      for(int rep = 0; rep < 4; rep++)
        Ewd_table[rep][i][j][k] = fcorr[rep];
    }

  int *recvcnts = (int *)mymalloc("recvcnts", NTask * sizeof(int));
  int *recvoffs = (int *)mymalloc("recvoffs", NTask * sizeof(int));

  for(int i = 0; i < NTask; i++)
    {
      int off, cnt;
      subdivide_evenly(size, NTask, i, &off, &cnt);
      recvcnts[i] = cnt * sizeof(double);
      recvoffs[i] = off * sizeof(double);
    }

  for(int rep = 0; rep < 4; rep++)
    MPI_Allgatherv(MPI_IN_PLACE, size * sizeof(double), MPI_BYTE, Ewd_table[rep], recvcnts, recvoffs, MPI_BYTE, MPI_COMM_WORLD);

  myfree(recvoffs);
  myfree(recvcnts);

  /* now scale things to the boxsize that is actually used */
  Ewd_table_intp = 2 * TEW_N / All.BoxSize;

  double t1 = second();
  mpi_printf("FORCETEST: Initialization of high-res Ewald table finished, took %g sec.\n", timediff(t0, t1));
}

/*! \brief Looks up Ewald force from tabulated values.
 *
 *  \param[in] dx X position.
 *  \param[in] dy Y position.
 *  \param[in] dz Z position.
 *  \param[out] force Ewald force correction.
 *
 *  \return void
 */
static void ewald_correction_force_table_lookup(double dx, double dy, double dz, double force[4])
{
  int signx, signy, signz;
  int i, j, k;
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;

  if(dx < 0)
    {
      dx    = -dx;
      signx = -1;
    }
  else
    signx = +1;

  if(dy < 0)
    {
      dy    = -dy;
      signy = -1;
    }
  else
    signy = +1;

  if(dz < 0)
    {
      dz    = -dz;
      signz = -1;
    }
  else
    signz = +1;

  u = dx * Ewd_table_intp;
  i = (int)u;
  if(i >= TEW_NX)
    i = TEW_NX - 1;
  u -= i;
  v = dy * Ewd_table_intp;
  j = (int)v;
  if(j >= TEW_NY)
    j = TEW_NY - 1;
  v -= j;
  w = dz * Ewd_table_intp;
  k = (int)w;
  if(k >= TEW_NZ)
    k = TEW_NZ - 1;
  w -= k;

  f1 = (1 - u) * (1 - v) * (1 - w);
  f2 = (1 - u) * (1 - v) * (w);
  f3 = (1 - u) * (v) * (1 - w);
  f4 = (1 - u) * (v) * (w);
  f5 = (u) * (1 - v) * (1 - w);
  f6 = (u) * (1 - v) * (w);
  f7 = (u) * (v) * (1 - w);
  f8 = (u) * (v) * (w);

  for(int rep = 0; rep < 4; rep++)
    {
      force[rep] = Ewd_table[rep][i][j][k] * f1 + Ewd_table[rep][i][j][k + 1] * f2 + Ewd_table[rep][i][j + 1][k] * f3 +
                   Ewd_table[rep][i][j + 1][k + 1] * f4 + Ewd_table[rep][i + 1][j][k] * f5 + Ewd_table[rep][i + 1][j][k + 1] * f6 +
                   Ewd_table[rep][i + 1][j + 1][k] * f7 + Ewd_table[rep][i + 1][j + 1][k + 1] * f8;
    }

  force[0] *= signx;
  force[1] *= signy;
  force[2] *= signz;
}

#endif /* #if !defined(FORCETEST_TESTFORCELAW) */

#endif /* #ifdef FORCETEST */
