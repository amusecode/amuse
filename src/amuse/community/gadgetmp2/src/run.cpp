/* ################################################################################## */
/* ###                                                                            ### */
/* ###                                 Gadgetmp2                                  ### */
/* ###                                                                            ### */
/* ###   Original: Gadget2 in the version used in Amuse                           ### */
/* ###   Author: Gadget2 and Amuse contributors                                   ### */
/* ###                                                                            ### */
/* ###   Modified: July 2020                                                      ### */
/* ###   Author: Thomas Schano                                                    ### */
/* ###                                                                            ### */
/* ###   Changes are intended to enable precise calculations in                   ### */
/* ###   non periodic small domain simulations in which comoving parts            ### */
/* ###   are simulated in std precision                                           ### */
/* ###                                                                            ### */
/* ################################################################################## */
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
void gadgetmp2::every_timestep_stuff(void)
{
  my_float z;

  if(ThisTask == 0)
    {
      if(All.ComovingIntegrationOn)
	{
	  z = 1.0 / (All.Time) - 1;
	  fprintf(FdInfo, "\nBegin Step %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n",
		  All.NumCurrentTiStep, All.Time, z.toDouble(), All.TimeStep,
		  (log(All.Time) - log(All.Time - All.TimeStep)));
	  printf("\nBegin Step %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n", All.NumCurrentTiStep,
		 All.Time, z.toDouble(), All.TimeStep, (log(All.Time) - log(All.Time - All.TimeStep)));
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
void gadgetmp2::energy_statistics(void)
{
  compute_global_quantities_of_system();

  if(ThisTask == 0)
    {
      fprintf(FdEnergy,
	      "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      All.Time, SysState.EnergyInt.toDouble(), SysState.EnergyPot.toDouble(), SysState.EnergyKin.toDouble(), SysState.EnergyIntComp[0].toDouble(),
	      SysState.EnergyPotComp[0].toDouble(), SysState.EnergyKinComp[0].toDouble(), SysState.EnergyIntComp[1].toDouble(),
	      SysState.EnergyPotComp[1].toDouble(), SysState.EnergyKinComp[1].toDouble(), SysState.EnergyIntComp[2].toDouble(),
	      SysState.EnergyPotComp[2].toDouble(), SysState.EnergyKinComp[2].toDouble(), SysState.EnergyIntComp[3].toDouble(),
	      SysState.EnergyPotComp[3].toDouble(), SysState.EnergyKinComp[3].toDouble(), SysState.EnergyIntComp[4].toDouble(),
	      SysState.EnergyPotComp[4].toDouble(), SysState.EnergyKinComp[4].toDouble(), SysState.EnergyIntComp[5].toDouble(),
	      SysState.EnergyPotComp[5].toDouble(), SysState.EnergyKinComp[5].toDouble(), SysState.MassComp[0].toDouble(),
	      SysState.MassComp[1].toDouble(), SysState.MassComp[2].toDouble(), SysState.MassComp[3].toDouble(), SysState.MassComp[4].toDouble(),
	      SysState.MassComp[5].toDouble());

      fflush(FdEnergy);
    }
}
