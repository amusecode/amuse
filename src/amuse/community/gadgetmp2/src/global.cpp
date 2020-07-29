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

//#include "allvars.hpp"
#include "proto.hpp"


/*! \file global.c
 *  \brief Computes global physical properties of the system
 */


/*! This routine computes various global properties of the particle
 *  distribution and stores the result in the struct `SysState'.
 *  Currently, not all the information that's computed here is actually
 *  used (e.g. momentum is not really used anywhere), just the energies are
 *  written to a log-file every once in a while.
 */
void gadgetmp2::compute_global_quantities_of_system(void)
{
  int i, j, n;
  struct state_of_system sys;
  my_float a1, a2, a3;
  my_float entr = 0, egyspec, vel[3];
  my_float dt_entr, dt_gravkick, dt_hydrokick;
  my_float_buff *gather_buff;
  my_float_buff *sys_buff=my_float_buff::generate_buffer(state_of_system::max_cnt1*state_of_system::max_cnt2);
  size_t buff_size=sys_buff->get_transfer_buff_size();

  if(buff_size*NTask > All.BufferSize * 1024 * 1024)
    exit(40);

  if(All.ComovingIntegrationOn)
    {
      a1 = All.Time;
      a2 = All.Time * All.Time;
      a3 = All.Time * All.Time * All.Time;
    }
  else
    {
      a1 = a2 = a3 = 1;
    }


  for(n = 0; n < 6; n++)
    {
      sys.MassComp[n] = sys.EnergyKinComp[n] = sys.EnergyPotComp[n] = sys.EnergyIntComp[n] = 0;

      for(j = 0; j < 4; j++)
	sys.CenterOfMassComp[n][j] = sys.MomentumComp[n][j] = sys.AngMomentumComp[n][j] = 0;
    }

  for(i = 0; i < NumPart; i++)
    {
      sys.MassComp[P[i].Type] += P[i].Mass;

      sys.EnergyPotComp[P[i].Type] += 0.5 * P[i].Mass * P[i].Potential / a1;

      if(All.ComovingIntegrationOn)
	{
	  dt_entr = (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;
	  dt_gravkick = get_gravkick_factor(P[i].Ti_begstep, All.Ti_Current) -
	    get_gravkick_factor(P[i].Ti_begstep, (P[i].Ti_begstep + P[i].Ti_endstep) / 2);
	  dt_hydrokick = get_hydrokick_factor(P[i].Ti_begstep, All.Ti_Current) -
	    get_hydrokick_factor(P[i].Ti_begstep, (P[i].Ti_begstep + P[i].Ti_endstep) / 2);
	}
      else
	dt_entr = dt_gravkick = dt_hydrokick =
	  (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;

      for(j = 0; j < 3; j++)
	{
	  vel[j] = P[i].Vel[j] + P[i].GravAccel[j] * dt_gravkick;
	  if(P[i].Type == 0)
	    vel[j] += SphP[i].HydroAccel[j] * dt_hydrokick;
	}
      if(P[i].Type == 0)
	entr = SphP[i].Entropy + SphP[i].DtEntropy * dt_entr;


      sys.EnergyKinComp[P[i].Type] +=
	0.5 * P[i].Mass * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) / a2;

      if(P[i].Type == 0)
	{
#ifdef ISOTHERM_EQS
          egyspec = entr;
#else
	  egyspec = entr / (GAMMA_MINUS1) * pow(SphP[i].Density / a3, GAMMA_MINUS1);
#endif
	  sys.EnergyIntComp[0] += P[i].Mass * egyspec;
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
#ifndef NOMPI
    for(i=0;i<state_of_system::max_cnt1;i++)
    {
        sys_buff[i]=sys.MassComp[i];
        sys_buff[i+state_of_system::max_cnt1]=sys.EnergyPotComp[i];
        sys_buff[i+state_of_system::max_cnt1*2]=sys.EnergyIntComp[i];
        sys_buff[i+state_of_system::max_cnt1*3]=sys.EnergyKinComp[i];
        SysState.MassComp[i]=0;
        SysState.EnergyPotComp[i]=0;
        SysState.EnergyIntComp[i]=0;
        SysState.EnergyKinComp[i]=0;
    }
    MPI_Gather(sys_buff->get_transfer_buff(),buff_size,MPI_BYTE,CommBuffer,buff_size,MPI_BYTE,0,GADGET_WORLD);
    if(ThisTask == 0)
    {
        gather_buff=my_float_buff::re_org_transfer_buff(CommBuffer);
        for(n=0;n<NTask;n++)
        {
            for(i=0;i<state_of_system::max_cnt1;i++)
            {
                SysState.MassComp[i]+=gather_buff[i];
                SysState.EnergyPotComp[i]+=gather_buff[i+state_of_system::max_cnt1];
                SysState.EnergyIntComp[i]+=gather_buff[i+state_of_system::max_cnt1*2];
                SysState.EnergyKinComp[i]+=gather_buff[i+state_of_system::max_cnt1*3];
            }
            gather_buff=my_float_buff::re_org_transfer_buff(gather_buff->next_pointer());
        }
    }
  //MPI_Reduce(&sys.MassComp[0], &SysState.MassComp[0], 6, MPI_DOUBLE, MPI_SUM, 0, GADGET_WORLD);
  //MPI_Reduce(&sys.EnergyPotComp[0], &SysState.EnergyPotComp[0], 6, MPI_DOUBLE, MPI_SUM, 0, GADGET_WORLD);
  //MPI_Reduce(&sys.EnergyIntComp[0], &SysState.EnergyIntComp[0], 6, MPI_DOUBLE, MPI_SUM, 0, GADGET_WORLD);
  //MPI_Reduce(&sys.EnergyKinComp[0], &SysState.EnergyKinComp[0], 6, MPI_DOUBLE, MPI_SUM, 0, GADGET_WORLD);
  for(j=0;j<state_of_system::max_cnt2;j++)
    for(i=0;i<state_of_system::max_cnt1;i++)
    {
        sys_buff[i+state_of_system::max_cnt1*j]=sys.MomentumComp[i][j];
        SysState.MomentumComp[i][j]=0;
    }
    MPI_Gather(sys_buff->get_transfer_buff(),buff_size,MPI_BYTE,CommBuffer,buff_size,MPI_BYTE,0,GADGET_WORLD);
    if(ThisTask == 0)
    {
        gather_buff=my_float_buff::re_org_transfer_buff(CommBuffer);
        for(n=0;n<NTask;n++)
        {
          for(j=0;j<state_of_system::max_cnt2;j++)
            for(i=0;i<state_of_system::max_cnt1;i++)
            {
                SysState.MomentumComp[i][j]+=gather_buff[i+state_of_system::max_cnt1*j];
            }
            gather_buff=my_float_buff::re_org_transfer_buff(gather_buff->next_pointer());
        }
    }
  //MPI_Reduce(&sys.MomentumComp[0][0], &SysState.MomentumComp[0][0], 6 * 4, MPI_DOUBLE, MPI_SUM, 0, GADGET_WORLD);

  for(j=0;j<state_of_system::max_cnt2;j++)
    for(i=0;i<state_of_system::max_cnt1;i++)
    {
        sys_buff[i+state_of_system::max_cnt1*j]=sys.AngMomentumComp[i][j];
        SysState.AngMomentumComp[i][j]=0;
    }
    MPI_Gather(sys_buff->get_transfer_buff(),buff_size,MPI_BYTE,CommBuffer,buff_size,MPI_BYTE,0,GADGET_WORLD);
    if(ThisTask == 0)
    {
        gather_buff=my_float_buff::re_org_transfer_buff(CommBuffer);
        for(n=0;n<NTask;n++)
        {
          for(j=0;j<state_of_system::max_cnt2;j++)
            for(i=0;i<state_of_system::max_cnt1;i++)
            {
                SysState.AngMomentumComp[i][j]+=gather_buff[i+state_of_system::max_cnt1*j];
            }
            gather_buff=my_float_buff::re_org_transfer_buff(gather_buff->next_pointer());
        }
    }
  //MPI_Reduce(&sys.AngMomentumComp[0][0], &SysState.AngMomentumComp[0][0], 6 * 4, MPI_DOUBLE, MPI_SUM, 0, GADGET_WORLD);
  for(j=0;j<state_of_system::max_cnt2;j++)
    for(i=0;i<state_of_system::max_cnt1;i++)
    {
        sys_buff[i+state_of_system::max_cnt1*j]=sys.CenterOfMassComp[i][j];
        SysState.CenterOfMassComp[i][j]=0;
    }
    MPI_Gather(sys_buff->get_transfer_buff(),buff_size,MPI_BYTE,CommBuffer,buff_size,MPI_BYTE,0,GADGET_WORLD);
    if(ThisTask == 0)
    {
        gather_buff=my_float_buff::re_org_transfer_buff(CommBuffer);
        for(n=0;n<NTask;n++)
        {
          for(j=0;j<state_of_system::max_cnt2;j++)
            for(i=0;i<state_of_system::max_cnt1;i++)
            {
                SysState.CenterOfMassComp[i][j]+=gather_buff[i+state_of_system::max_cnt1*j];
            }
            gather_buff=my_float_buff::re_org_transfer_buff(gather_buff->next_pointer());
        }
    }
  //MPI_Reduce(&sys.CenterOfMassComp[0][0], &SysState.CenterOfMassComp[0][0], 6 * 4, MPI_DOUBLE, MPI_SUM, 0, GADGET_WORLD);
#else
    for(i = 0; i < 6; i++){
        SysState.MassComp[i] = sys.MassComp[i];
        SysState.EnergyPotComp[i] = sys.EnergyPotComp[i];
        SysState.EnergyIntComp[i] = sys.EnergyIntComp[i];
        SysState.EnergyKinComp[i] = sys.EnergyKinComp[i];
        SysState.MassComp[i] = sys.MassComp[i];
        for(j = 0; j < 4; j++){
            SysState.MomentumComp[i][j] = sys.MomentumComp[i][j];
            SysState.AngMomentumComp[i][j] = sys.AngMomentumComp[i][j];
            SysState.CenterOfMassComp[i][j] = sys.CenterOfMassComp[i][j];
        }
    }
#endif

  if(ThisTask == 0)
    {
      for(i = 0; i < 6; i++)
	SysState.EnergyTotComp[i] = SysState.EnergyKinComp[i] +
	  SysState.EnergyPotComp[i] + SysState.EnergyIntComp[i];

      SysState.Mass = SysState.EnergyKin = SysState.EnergyPot = SysState.EnergyInt = SysState.EnergyTot = 0;

      for(j = 0; j < 3; j++)
	SysState.Momentum[j] = SysState.AngMomentum[j] = SysState.CenterOfMass[j] = 0;

      for(i = 0; i < 6; i++)
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

      for(i = 0; i < 6; i++)
	for(j = 0; j < 3; j++)
	  if(SysState.MassComp[i] > 0)
	    SysState.CenterOfMassComp[i][j] /= SysState.MassComp[i];

      for(j = 0; j < 3; j++)
	if(SysState.Mass > 0)
	  SysState.CenterOfMass[j] /= SysState.Mass;

      for(i = 0; i < 6; i++)
	{
	  SysState.CenterOfMassComp[i][3] = SysState.MomentumComp[i][3] = SysState.AngMomentumComp[i][3] = 0;
	  for(j = 0; j < 3; j++)
	    {
	      SysState.CenterOfMassComp[i][3] +=
		SysState.CenterOfMassComp[i][j] * SysState.CenterOfMassComp[i][j];
	      SysState.MomentumComp[i][3] += SysState.MomentumComp[i][j] * SysState.MomentumComp[i][j];
	      SysState.AngMomentumComp[i][3] +=
		SysState.AngMomentumComp[i][j] * SysState.AngMomentumComp[i][j];
	    }
	  SysState.CenterOfMassComp[i][3] = sqrt(SysState.CenterOfMassComp[i][3]);
	  SysState.MomentumComp[i][3] = sqrt(SysState.MomentumComp[i][3]);
	  SysState.AngMomentumComp[i][3] = sqrt(SysState.AngMomentumComp[i][3]);
	}

      SysState.CenterOfMass[3] = SysState.Momentum[3] = SysState.AngMomentum[3] = 0;

      for(j = 0; j < 3; j++)
	{
	  SysState.CenterOfMass[3] += SysState.CenterOfMass[j] * SysState.CenterOfMass[j];
	  SysState.Momentum[3] += SysState.Momentum[j] * SysState.Momentum[j];
	  SysState.AngMomentum[3] += SysState.AngMomentum[j] * SysState.AngMomentum[j];
	}

      SysState.CenterOfMass[3] = sqrt(SysState.CenterOfMass[3]);
      SysState.Momentum[3] = sqrt(SysState.Momentum[3]);
      SysState.AngMomentum[3] = sqrt(SysState.AngMomentum[3]);
    }

  /* give everyone the result, maybe the want to do something with it */
#ifndef NOMPI
//  MPI_Bcast(&SysState, sizeof(struct state_of_system), MPI_BYTE, 0, GADGET_WORLD);
#endif // NOMPI
}
