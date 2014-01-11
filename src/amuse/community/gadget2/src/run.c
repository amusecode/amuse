#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef NOMPI
#include <mpi.h>
#endif
#include <unistd.h>

#include "allvars.h"
#include "proto.h"

/*! \file run.c
 *  \brief  iterates over timesteps, main loop
 */

/*! This routine contains the main simulation loop that iterates over single
 *  timesteps. The loop terminates when the cpu-time limit is reached, when a
 *  `stop' file is found in the output directory, or when the simulation ends
 *  because we arrived at TimeMax.
 */
void run(void)
{
  FILE *fd;
  int stopflag = 0;
  char stopfname[200], contfname[200];
  double t0, t1;


  sprintf(stopfname, "%sstop", All.OutputDir);
  sprintf(contfname, "%scont", All.OutputDir);
  unlink(contfname);

  do				/* main loop */
    {
      t0 = second();

      find_next_sync_point_and_drift();	/* find next synchronization point and drift particles to this time.
					 * If needed, this function will also write an output file
					 * at the desired time.
					 */

      every_timestep_stuff();	/* write some info to log-files */


      domain_Decomposition();	/* do domain decomposition if needed */


      compute_accelerations(0);	/* compute accelerations for
				 * the particles that are to be advanced
				 */

      /* check whether we want a full energy statistics */
      if((All.Time - All.TimeLastStatistics) >= All.TimeBetStatistics)
	{
#ifdef COMPUTE_POTENTIAL_ENERGY
	  compute_potential();
#endif
	  energy_statistics();	/* compute and output energy statistics */
	  All.TimeLastStatistics += All.TimeBetStatistics;
	}

      advance_and_find_timesteps();	/* 'kick' active particles in
					 * momentum space and compute new
					 * timesteps for them
					 */
      All.NumCurrentTiStep++;

      /* Check whether we need to interrupt the run */
      if(ThisTask == 0)
	{
	  /* Is the stop-file present? If yes, interrupt the run. */
	  if((fd = fopen(stopfname, "r")))
	    {
	      fclose(fd);
	      stopflag = 1;
	      unlink(stopfname);
	    }

	  /* are we running out of CPU-time ? If yes, interrupt run. */
	  if(CPUThisRun > 0.85 * All.TimeLimitCPU)
	    {
	      printf("reaching time-limit. stopping.\n");
	      stopflag = 2;
	    }
	}

#ifndef NOMPI
      MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
      if(stopflag)
	{
	  restart(0);		/* write restart file */
#ifndef NOMPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif

	  if(stopflag == 2 && ThisTask == 0)
	    {
	      if((fd = fopen(contfname, "w")))
		fclose(fd);
	    }

	  if(stopflag == 2 && All.ResubmitOn && ThisTask == 0)
	    {
	      close_outputfiles();
	      system(All.ResubmitCommand);
	    }
	  return;
	}

      /* is it time to write a regular restart-file? (for security) */
      if(ThisTask == 0)
	{
	  if((CPUThisRun - All.TimeLastRestartFile) >= All.CpuTimeBetRestartFile)
	    {
	      All.TimeLastRestartFile = CPUThisRun;
	      stopflag = 3;
	    }
	  else
	    stopflag = 0;
	}

#ifndef NOMPI
      MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
      if(stopflag == 3)
	{
	  restart(0);		/* write an occasional restart file */
	  stopflag = 0;
	}

      t1 = second();

      All.CPU_Total += timediff(t0, t1);
      CPUThisRun += timediff(t0, t1);
    }
  while(All.Ti_Current < TIMEBASE && All.Time <= All.TimeMax);

  restart(0);

  savepositions(All.SnapshotFileCount++);	/* write a last snapshot
						 * file at final time (will
						 * be overwritten if
						 * All.TimeMax is increased
						 * and the run is continued)
						 */
}


/*! This function finds the next synchronization point of the system (i.e. the
 *  earliest point of time any of the particles needs a force computation),
 *  and drifts the system to this point of time.  If the system drifts over
 *  the desired time of a snapshot file, the function will drift to this
 *  moment, generate an output, and then resume the drift.
 */
void find_next_sync_point_and_drift(void)
{
  int n, min, min_glob, flag, *temp;
  double timeold;
  double t0, t1;

  t0 = second();

  timeold = All.Time;

  for(n = 1, min = P[0].Ti_endstep; n < NumPart; n++)
    if(min > P[n].Ti_endstep)
      min = P[n].Ti_endstep;

#ifndef NOMPI
  MPI_Allreduce(&min, &min_glob, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#else
    min_glob = min;
#endif
  /* We check whether this is a full step where all particles are synchronized */
  flag = 1;
  for(n = 0; n < NumPart; n++)
    if(P[n].Ti_endstep > min_glob)
      flag = 0;

#ifndef NOMPI
  MPI_Allreduce(&flag, &Flag_FullStep, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#else
    Flag_FullStep = flag;
#endif
#ifdef PMGRID
  if(min_glob >= All.PM_Ti_endstep)
    {
      min_glob = All.PM_Ti_endstep;
      Flag_FullStep = 1;
    }
#endif

  /* Determine 'NumForceUpdate', i.e. the number of particles on this processor that are going to be active */
  for(n = 0, NumForceUpdate = 0; n < NumPart; n++)
    {
      if(P[n].Ti_endstep == min_glob)
#ifdef SELECTIVE_NO_GRAVITY
        if(!((1 << P[n].Type) & (SELECTIVE_NO_GRAVITY)))
#endif
          NumForceUpdate++;
    }

  /* note: NumForcesSinceLastDomainDecomp has type "long long" */
  temp = malloc(NTask * sizeof(int));
#ifndef NOMPI
  MPI_Allgather(&NumForceUpdate, 1, MPI_INT, temp, 1, MPI_INT, MPI_COMM_WORLD);
#else
    temp[0] = NumForceUpdate;
#endif
  for(n = 0; n < NTask; n++)
    All.NumForcesSinceLastDomainDecomp += temp[n];
  free(temp);



  t1 = second();

  All.CPU_Predict += timediff(t0, t1);

  while(min_glob >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
    {
      move_particles(All.Ti_Current, All.Ti_nextoutput);

      All.Ti_Current = All.Ti_nextoutput;

      if(All.ComovingIntegrationOn)
	All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
      else
	All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;

#ifdef OUTPUTPOTENTIAL
      All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;
      domain_Decomposition();
      compute_potential();
#endif
      savepositions(All.SnapshotFileCount++);	/* write snapshot file */

      All.Ti_nextoutput = find_next_outputtime(All.Ti_nextoutput + 1);
    }

  move_particles(All.Ti_Current, min_glob);

  All.Ti_Current = min_glob;

  if(All.ComovingIntegrationOn)
    All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
  else
    All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;

  All.TimeStep = All.Time - timeold;
}



/*! this function returns the next output time that is equal or larger to
 *  ti_curr
 */
int find_next_outputtime(int ti_curr)
{
  int i, ti, ti_next, iter = 0;
  double next, time;

  ti_next = -1;

  if(All.OutputListOn)
    {
      for(i = 0; i < All.OutputListLength; i++)
	{
	  time = All.OutputListTimes[i];

	  if(time >= All.TimeBegin && time <= All.TimeMax)
	    {
	      if(All.ComovingIntegrationOn)
		ti = log(time / All.TimeBegin) / All.Timebase_interval;
	      else
		ti = (time - All.TimeBegin) / All.Timebase_interval;

	      if(ti >= ti_curr)
		{
		  if(ti_next == -1)
		    ti_next = ti;

		  if(ti_next > ti)
		    ti_next = ti;
		}
	    }
	}
    }
  else
    {
      if(All.ComovingIntegrationOn)
	{
	  if(All.TimeBetSnapshot <= 1.0)
	    {
	      printf("TimeBetSnapshot > 1.0 required for your simulation.\n");
	      endrun(13123);
	    }
	}
      else
	{
	  if(All.TimeBetSnapshot <= 0.0)
	    {
	      printf("TimeBetSnapshot > 0.0 required for your simulation.\n");
	      endrun(13123);
	    }
	}

      time = All.TimeOfFirstSnapshot;

      iter = 0;

      while(time < All.TimeBegin)
	{
	  if(All.ComovingIntegrationOn)
	    time *= All.TimeBetSnapshot;
	  else
	    time += All.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun(110);
	    }
	}

      while(time <= All.TimeMax)
	{
	  if(All.ComovingIntegrationOn)
	    ti = log(time / All.TimeBegin) / All.Timebase_interval;
	  else
	    ti = (time - All.TimeBegin) / All.Timebase_interval;

	  if(ti >= ti_curr)
	    {
	      ti_next = ti;
	      break;
	    }

	  if(All.ComovingIntegrationOn)
	    time *= All.TimeBetSnapshot;
	  else
	    time += All.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun(111);
	    }
	}
    }

  if(ti_next == -1)
    {
      ti_next = 2 * TIMEBASE;	/* this will prevent any further output */

      if(ThisTask == 0)
	printf("\nThere is no valid time for a further snapshot file.\n");
    }
  else
    {
      if(All.ComovingIntegrationOn)
	next = All.TimeBegin * exp(ti_next * All.Timebase_interval);
      else
	next = All.TimeBegin + ti_next * All.Timebase_interval;

      if(ThisTask == 0)
	printf("\nSetting next time for snapshot file to Time_next= %g\n\n", next);
    }

  return ti_next;
}




/*! This routine writes one line for every timestep to two log-files.  In
 *  FdInfo, we just list the timesteps that have been done, while in FdCPU the
 *  cumulative cpu-time consumption in various parts of the code is stored.
 */
void every_timestep_stuff(void)
{
  double z;

  if(ThisTask == 0)
    {
      if(All.ComovingIntegrationOn)
	{
	  z = 1.0 / (All.Time) - 1;
	  fprintf(FdInfo, "\nBegin Step %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n",
		  All.NumCurrentTiStep, All.Time, z, All.TimeStep,
		  log(All.Time) - log(All.Time - All.TimeStep));
	  printf("\nBegin Step %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n", All.NumCurrentTiStep,
		 All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
	  fflush(FdInfo);
	}
      else
	{
	  fprintf(FdInfo, "\nBegin Step %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time,
		  All.TimeStep);
	  printf("\nBegin Step %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time, All.TimeStep);
	  fflush(FdInfo);
	}

      fprintf(FdCPU, "Step %d, Time: %g, CPUs: %d\n", All.NumCurrentTiStep, All.Time, NTask);

      fprintf(FdCPU,
	      "%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
	      All.CPU_Total, All.CPU_Gravity, All.CPU_Hydro, All.CPU_Domain, All.CPU_Potential,
	      All.CPU_Predict, All.CPU_TimeLine, All.CPU_Snapshot, All.CPU_TreeWalk, All.CPU_TreeConstruction,
	      All.CPU_CommSum, All.CPU_Imbalance, All.CPU_HydCompWalk, All.CPU_HydCommSumm,
	      All.CPU_HydImbalance, All.CPU_EnsureNgb, All.CPU_PM, All.CPU_Peano);
      fflush(FdCPU);
    }

  set_random_numbers();
}


/*! This routine first calls a computation of various global quantities of the
 *  particle distribution, and then writes some statistics about the energies
 *  in the various particle components to the file FdEnergy.
 */
void energy_statistics(void)
{
  compute_global_quantities_of_system();

  if(ThisTask == 0)
    {
      fprintf(FdEnergy,
	      "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      All.Time, SysState.EnergyInt, SysState.EnergyPot, SysState.EnergyKin, SysState.EnergyIntComp[0],
	      SysState.EnergyPotComp[0], SysState.EnergyKinComp[0], SysState.EnergyIntComp[1],
	      SysState.EnergyPotComp[1], SysState.EnergyKinComp[1], SysState.EnergyIntComp[2],
	      SysState.EnergyPotComp[2], SysState.EnergyKinComp[2], SysState.EnergyIntComp[3],
	      SysState.EnergyPotComp[3], SysState.EnergyKinComp[3], SysState.EnergyIntComp[4],
	      SysState.EnergyPotComp[4], SysState.EnergyKinComp[4], SysState.EnergyIntComp[5],
	      SysState.EnergyPotComp[5], SysState.EnergyKinComp[5], SysState.MassComp[0],
	      SysState.MassComp[1], SysState.MassComp[2], SysState.MassComp[3], SysState.MassComp[4],
	      SysState.MassComp[5]);

      fflush(FdEnergy);
    }
}
