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
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

//#include "allvars.hpp"
#include "proto.hpp"


/*! \file begrun.c
 *  \brief initial set-up of a simulation run
 *
 *  This file contains various functions to initialize a simulation run. In
 *  particular, the parameterfile is read in and parsed, the initial
 *  conditions or restart files are read, and global variables are
 *  initialized to their proper values.
 */





/*! Computes conversion factors between internal code units and the
 *  cgs-system.
 */
void gadgetmp2::set_units(void)
{
  my_float meanweight;

  All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / std::pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * std::pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;

  All.UnitDensity_in_cgs = All.UnitMass_in_g / std::pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / std::pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * std::pow(All.UnitLength_in_cm, 2) / std::pow(All.UnitTime_in_s, 2);

  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;

  if(ThisTask == 0)
    {
      printf("\nHubble (internal units) = %g\n", All.Hubble);
      printf("G (internal units) = %g\n", All.G);
      printf("UnitMass_in_g = %g \n", All.UnitMass_in_g);
      printf("UnitTime_in_s = %g \n", All.UnitTime_in_s);
      printf("UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
      printf("UnitDensity_in_cgs = %g \n", All.UnitDensity_in_cgs);
      printf("UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
      printf("\n");
    }

  meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: we assume neutral gas here */

#ifdef ISOTHERM_EQS
  All.MinEgySpec = 0;
#else
  All.MinEgySpec = (1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp).toDouble();
  All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
#endif

}



/*!  This function opens various log-files that report on the status and
 *   performance of the simulstion. On restart from restart-files
 *   (start-option 1), the code will append to these files.
 */
void gadgetmp2::open_outputfiles(void)
{
  char mode[2], buf[200];

    sprintf(buf, "%s%s%d%s", All.OutputDir, "debug",(int)ThisTask,".txt");
    DEBUG.open (buf);
    DEBUG << "debugged task " << ThisTask <<"\n";
    DEBUG.flush();

  if(ThisTask != 0)		/* only the root processor writes to the log files */
    return;

  if(RestartFlag == 0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");


  sprintf(buf, "%s%s", All.OutputDir, All.CpuFile);
  if(!(FdCPU = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.InfoFile);
  if(!(FdInfo = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.EnergyFile);
  if(!(FdEnergy = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.TimingsFile);
  if(!(FdTimings = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
/*  sprintf(buf, "%s%s", All.OutputDir, "debug.txt");
    DEBUG.open (buf);
    DEBUG << "test  te" << 2 <<"\n";
    DEBUG.flush();*/

#ifdef FORCETEST
  if(RestartFlag == 0)
    {
      sprintf(buf, "%s%s", All.OutputDir, "forcetest.txt");
      if(!(FdForceTest = fopen(buf, "w")))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}
      fclose(FdForceTest);
    }
#endif
}


/*!  This function closes the global log-files.
 */
void gadgetmp2::close_outputfiles(void)
{
  if(ThisTask != 0)		/* only the root processor writes to the log files */
    return;

  fclose(FdCPU);
  fclose(FdInfo);
  fclose(FdEnergy);
  fclose(FdTimings);
#ifdef FORCETEST
  fclose(FdForceTest);
#endif
}




