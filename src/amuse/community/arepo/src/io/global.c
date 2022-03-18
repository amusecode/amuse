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
 * \file        src/global.c
 * \date        05/2018
 * \brief       Routines to compute statistics of the global state of the
 *              code.
 * \details     contains functions:
 *                void compute_statistics(void)
 *                void energy_statistics(void)
 *                void compute_global_quantities_of_system(void)
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

#include "../main/allvars.h"
#include "../main/proto.h"

/*! \brief Computes new global statistics if needed (call of
 *         energy_statistics()).
 *
 *  \return void
 */
void compute_statistics(void)
{
  /* check whether we want a full energy statistics */
  if((All.Time - All.TimeLastStatistics) >= All.TimeBetStatistics &&
     All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) /* allow only top-level synchronization points */
    {
      TIMER_START(CPU_LOGS);

      energy_statistics(); /* compute and output energy statistics */

      All.TimeLastStatistics += All.TimeBetStatistics;

      TIMER_STOP(CPU_LOGS);
    }
}

/*! \brief Compute global statistics of the system.
 *
 *  This function first calls a computation of various global
 *  quantities of the particle distribution
 *  (compute_global_quantities_of_system() ), and then writes some statistics
 *  about the energies of the various particle types to the file FdEnergy
 *  (energy.txt).
 *
 *  \return void
 */
void energy_statistics(void)
{
  double egyinj_tot;

  compute_global_quantities_of_system();

  MPI_Reduce(&EgyInjection, &egyinj_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdEnergy, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", All.Time,
              SysState.EnergyInt, SysState.EnergyPot, SysState.EnergyKin, SysState.EnergyIntComp[0], SysState.EnergyPotComp[0],
              SysState.EnergyKinComp[0], SysState.EnergyIntComp[1], SysState.EnergyPotComp[1], SysState.EnergyKinComp[1],
              SysState.EnergyIntComp[2], SysState.EnergyPotComp[2], SysState.EnergyKinComp[2], SysState.EnergyIntComp[3],
              SysState.EnergyPotComp[3], SysState.EnergyKinComp[3], SysState.EnergyIntComp[4], SysState.EnergyPotComp[4],
              SysState.EnergyKinComp[4], SysState.EnergyIntComp[5], SysState.EnergyPotComp[5], SysState.EnergyKinComp[5],
              SysState.MassComp[0], SysState.MassComp[1], SysState.MassComp[2], SysState.MassComp[3], SysState.MassComp[4],
              SysState.MassComp[5], egyinj_tot);

      myflush(FdEnergy);
    }
}

/*! \brief This routine computes various global properties of the particle
 *         distribution and stores the result in the struct `SysState'.
 *
 *  Currently, not all the information that's computed here is
 *  actually used (e.g. momentum is not really used anywhere),
 *  just the energies are written to a log-file every once in a while.
 *
 *  \return void
 */
void compute_global_quantities_of_system(void)
{
  int i, j, n;
  struct state_of_system sys;
  double egyspec, vel[3];

  for(n = 0; n < NTYPES; n++)
    {
      sys.MassComp[n] = sys.EnergyKinComp[n] = sys.EnergyPotComp[n] = sys.EnergyIntComp[n] = 0;

      for(j = 0; j < 4; j++)
        sys.CenterOfMassComp[n][j] = sys.MomentumComp[n][j] = sys.AngMomentumComp[n][j] = 0;
    }

  for(i = 0; i < NumPart; i++)
    {
      sys.MassComp[P[i].Type] += P[i].Mass;

#if defined(SELFGRAVITY)
#ifdef EVALPOTENTIAL
#ifndef EXACT_GRAVITY_FOR_PARTICLE_TYPE
      sys.EnergyPotComp[P[i].Type] +=
          0.5 * P[i].Mass * (P[i].Potential + All.G * P[i].Mass / (All.ForceSoftening[P[i].SofteningType] / 2.8)) / All.cf_atime;
#else  /* #ifndef EXACT_GRAVITY_FOR_PARTICLE_TYPE */
      /* ignore self-contribution from gravity if exact gravity is used */
      if(P[i].Type == EXACT_GRAVITY_FOR_PARTICLE_TYPE)
        sys.EnergyPotComp[P[i].Type] += 0.5 * P[i].Mass * P[i].Potential / All.cf_atime;
      else
        sys.EnergyPotComp[P[i].Type] +=
            0.5 * P[i].Mass * (P[i].Potential + All.G * P[i].Mass / (All.ForceSoftening[P[i].SofteningType] / 2.8)) / All.cf_atime;
#endif /* #ifndef EXACT_GRAVITY_FOR_PARTICLE_TYPE #else */
#endif /* #ifdef EVALPOTENTIAL */
#endif /* #if defined(SELFGRAVITY) */

#if defined(EXTERNALGRAVITY)
#if defined(SELFGRAVITY)
      sys.EnergyPotComp[P[i].Type] += 0.5 * P[i].Mass * P[i].ExtPotential; /* note: ExtPotential already included on P[].p.Potential,
                                                                              that's why only 0.5 is needed here to recover the rest */
#else                                                                      /* #if defined(SELFGRAVITY) */
      sys.EnergyPotComp[P[i].Type] += 1.0 * P[i].Mass * P[i].ExtPotential;
#endif                                                                     /* #if defined(SELFGRAVITY) #else */
#endif                                                                     /* #if defined(EXTERNALGRAVITY) */

      if(P[i].Type == 0)
        {
          for(j = 0; j < 3; j++)
            {
              vel[j] = P[i].Vel[j];
            }

          sys.EnergyKinComp[0] += 0.5 * P[i].Mass * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);

          egyspec = SphP[i].Utherm;

          sys.EnergyIntComp[0] += P[i].Mass * egyspec;
        }
      else
        {
          for(j = 0; j < 3; j++)
            {
              vel[j] = P[i].Vel[j];
            }
          sys.EnergyKinComp[P[i].Type] += 0.5 * P[i].Mass * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) * All.cf_a2inv;
        }

      for(j = 0; j < 3; j++)
        {
          sys.MomentumComp[P[i].Type][j] += P[i].Mass * vel[j];
          sys.CenterOfMassComp[P[i].Type][j] += P[i].Mass * P[i].Pos[j];
        }

      sys.AngMomentumComp[P[i].Type][0] += P[i].Mass * (P[i].Pos[1] * vel[2] - P[i].Pos[2] * vel[1]);
      sys.AngMomentumComp[P[i].Type][1] += P[i].Mass * (P[i].Pos[2] * vel[0] - P[i].Pos[0] * vel[2]);
      sys.AngMomentumComp[P[i].Type][2] += P[i].Mass * (P[i].Pos[0] * vel[1] - P[i].Pos[1] * vel[0]);
    }

  /* some the stuff over all processors */
  MPI_Reduce(&sys.MassComp[0], &SysState.MassComp[0], NTYPES, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.EnergyPotComp[0], &SysState.EnergyPotComp[0], NTYPES, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.EnergyIntComp[0], &SysState.EnergyIntComp[0], NTYPES, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.EnergyKinComp[0], &SysState.EnergyKinComp[0], NTYPES, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.MomentumComp[0][0], &SysState.MomentumComp[0][0], NTYPES * 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.AngMomentumComp[0][0], &SysState.AngMomentumComp[0][0], NTYPES * 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.CenterOfMassComp[0][0], &SysState.CenterOfMassComp[0][0], NTYPES * 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      for(i = 0; i < NTYPES; i++)
        SysState.EnergyTotComp[i] = SysState.EnergyKinComp[i] + SysState.EnergyPotComp[i] + SysState.EnergyIntComp[i];

      SysState.Mass = SysState.EnergyKin = SysState.EnergyPot = SysState.EnergyInt = SysState.EnergyTot = 0;

      for(j = 0; j < 3; j++)
        SysState.Momentum[j] = SysState.AngMomentum[j] = SysState.CenterOfMass[j] = 0;

      for(i = 0; i < NTYPES; i++)
        {
          SysState.Mass += SysState.MassComp[i];
          SysState.EnergyKin += SysState.EnergyKinComp[i];
          SysState.EnergyPot += SysState.EnergyPotComp[i];
          SysState.EnergyInt += SysState.EnergyIntComp[i];
          SysState.EnergyTot += SysState.EnergyTotComp[i];

          for(j = 0; j < 3; j++)
            {
              SysState.Momentum[j] += SysState.MomentumComp[i][j];
              SysState.AngMomentum[j] += SysState.AngMomentumComp[i][j];
              SysState.CenterOfMass[j] += SysState.CenterOfMassComp[i][j];
            }
        }

      for(i = 0; i < NTYPES; i++)
        for(j = 0; j < 3; j++)
          if(SysState.MassComp[i] > 0)
            SysState.CenterOfMassComp[i][j] /= SysState.MassComp[i];

      for(j = 0; j < 3; j++)
        if(SysState.Mass > 0)
          SysState.CenterOfMass[j] /= SysState.Mass;

      for(i = 0; i < NTYPES; i++)
        {
          SysState.CenterOfMassComp[i][3] = SysState.MomentumComp[i][3] = SysState.AngMomentumComp[i][3] = 0;
          for(j = 0; j < 3; j++)
            {
              SysState.CenterOfMassComp[i][3] += SysState.CenterOfMassComp[i][j] * SysState.CenterOfMassComp[i][j];
              SysState.MomentumComp[i][3] += SysState.MomentumComp[i][j] * SysState.MomentumComp[i][j];
              SysState.AngMomentumComp[i][3] += SysState.AngMomentumComp[i][j] * SysState.AngMomentumComp[i][j];
            }
          SysState.CenterOfMassComp[i][3] = sqrt(SysState.CenterOfMassComp[i][3]);
          SysState.MomentumComp[i][3]     = sqrt(SysState.MomentumComp[i][3]);
          SysState.AngMomentumComp[i][3]  = sqrt(SysState.AngMomentumComp[i][3]);
        }

      SysState.CenterOfMass[3] = SysState.Momentum[3] = SysState.AngMomentum[3] = 0;

      for(j = 0; j < 3; j++)
        {
          SysState.CenterOfMass[3] += SysState.CenterOfMass[j] * SysState.CenterOfMass[j];
          SysState.Momentum[3] += SysState.Momentum[j] * SysState.Momentum[j];
          SysState.AngMomentum[3] += SysState.AngMomentum[j] * SysState.AngMomentum[j];
        }

      SysState.CenterOfMass[3] = sqrt(SysState.CenterOfMass[3]);
      SysState.Momentum[3]     = sqrt(SysState.Momentum[3]);
      SysState.AngMomentum[3]  = sqrt(SysState.AngMomentum[3]);
    }

  /* give everyone the result, maybe the want to do something with it */
  MPI_Bcast(&SysState, sizeof(struct state_of_system), MPI_BYTE, 0, MPI_COMM_WORLD);
}
