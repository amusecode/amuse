#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef NOMPI
#include <mpi.h>
#endif
#include <unistd.h>

//#include "allvars.hpp"
#include "proto.hpp"

/*! \file run.c
 *  \brief  iterates over timesteps, main loop
 */






/*! This routine writes one line for every timestep to two log-files.  In
 *  FdInfo, we just list the timesteps that have been done, while in FdCPU the
 *  cumulative cpu-time consumption in various parts of the code is stored.
 */
void gadgetmg2::every_timestep_stuff(void)
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
void gadgetmg2::energy_statistics(void)
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
