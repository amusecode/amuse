/* ################################################################################## */
/* ###                                                                            ### */
/* ###                              timestep_limit.c                              ### */
/* ###                                                                            ### */
/* ###   Original: timestep.c (public version of Gadget 2)                        ### */
/* ###   Author: Volker Springel                                                  ### */
/* ###                                                                            ### */
/* ###   Modified: February 2011                                                  ### */
/* ###   Author: Fabrice Durier                                                   ### */
/* ###                                                                            ### */
/* ################################################################################## */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <assert.h>

#include "allvars.h"
#include "proto.h"

/*! \file timestep_limit.c 
 *
 *  This file is a modified version of timestep.c, which include the routines 
 *  for 'kicking' particles in momentum space and assigning new timesteps.
 *
 *  This file now includes the Timestep Limiter routine for communication of 
 *  the new timestep between neighbouring SPH particles
 */

#ifdef PERIODIC
static double boxSize, boxHalf;

#ifdef LONG_X
static double boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef LONG_Y
static double boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef LONG_Z
static double boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif
#endif

static double fac1, fac2, fac3, hubble_a, atime, a3inv;
static double dt_displacement = 0;

static int NDone;
static long long NTotDone;


/*! This function advances the system in momentum space, i.e. it does apply
 *  the 'kick' operation after the forces have been computed. Additionally, it
 *  assigns new timesteps to particles following the Limiter criterion. 
 *  At start-up, a half-timestep is carried out, as well as at the end of the simulation. 
 *  In between, the half-step kick that ends the previous timestep and 
 *  the half-step kick for the new timestep are combined into one operation.
 */

void advance_and_find_timesteps(void)
{
  long long ntot, ntotleft;
  int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
  int n, ndone, maxfill, source;
  int level, ngrp, sendTask, recvTask, place, nexport;
  double t0, t1;
  double sumt = 0, sumimbalance = 0, sumcomm = 0;
  double timecomp = 0, timeimbalance = 0, timecommsumm = 0;
  MPI_Status status;

  int CptLimit = 0;
  int shrinkcount = 0, shrinktot = 0;

  int i, j, no, ti_step, ti_min, tend, tstart;
  double dt_entr, dt_entr2, dt_gravkick, dt_hydrokick, dt_gravkick2, dt_hydrokick2;
  double minentropy, aphys;
  FLOAT dv[3];

#ifdef FLEXSTEPS
  int ti_grp;
#endif
#if defined(PSEUDOSYMMETRIC) && !defined(FLEXSTEPS)
  double apred, prob;
  int ti_step2;
#endif
#ifdef PMGRID
  double dt_gravkickA, dt_gravkickB;
#endif
#ifdef MAKEGLASS
  double disp, dispmax, globmax, dmean, fac, disp2sum, globdisp2sum;
#endif

#ifdef PERIODIC
  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
#endif
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
#endif
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
#endif
#endif

  noffset = (int*)malloc(sizeof(int) * NTask);  /* offsets of bunches in common list */
  nbuffer = (int*)malloc(sizeof(int) * NTask);
  nsend_local = (int*)malloc(sizeof(int) * NTask);
  nsend = (int*)malloc(sizeof(int) * NTask * NTask);
  ndonelist = (int*)malloc(sizeof(int) * NTask);

  if(ThisTask == 0)
    {
      fprintf(stdout,"Start time-stepping evaluation ...\n");
      fflush(stdout);
    }

  if(All.ComovingIntegrationOn)
    {
      fac1 = 1 / (All.Time * All.Time);
      fac2 = 1 / pow(All.Time, 3 * GAMMA - 2);
      fac3 = pow(All.Time, 3 * (1 - GAMMA) / 2.0);
      hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
	+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) + All.OmegaLambda;

      hubble_a = All.Hubble * sqrt(hubble_a);
      a3inv = 1 / (All.Time * All.Time * All.Time);
      atime = All.Time;
    }
  else
    fac1 = fac2 = fac3 = hubble_a = a3inv = atime = 1;

  if(Flag_FullStep || dt_displacement == 0)
    find_dt_displacement_constraint(hubble_a * atime * atime);

#ifdef PMGRID
  if(All.ComovingIntegrationOn)
    dt_gravkickB = get_gravkick_factor(All.PM_Ti_begstep, All.Ti_Current) -
      get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
  else
    dt_gravkickB = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;

  if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
      /* make sure that we reconstruct the domain/tree next time because we don't kick the tree nodes in this case */
      All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;
    }
#endif


#ifdef MAKEGLASS
  for(i = 0, dispmax = 0, disp2sum = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].GravPM[j] *= -1;
	  P[i].GravAccel[j] *= -1;
	  P[i].GravAccel[j] += P[i].GravPM[j];
	  P[i].GravPM[j] = 0;
	}

      disp = sqrt(P[i].GravAccel[0] * P[i].GravAccel[0] +
		  P[i].GravAccel[1] * P[i].GravAccel[1] + P[i].GravAccel[2] * P[i].GravAccel[2]);

      disp *= 2.0 / (3 * All.Hubble * All.Hubble);

      disp2sum += disp * disp;

      if(disp > dispmax)
	dispmax = disp;
    }

  MPI_Allreduce(&dispmax, &globmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&disp2sum, &globdisp2sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  dmean = pow(P[0].Mass / (All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)), 1.0 / 3);

  if(globmax > dmean)
    fac = dmean / globmax;
  else
    fac = 1.0;

  if(ThisTask == 0)
    {
      printf("\nglass-making:  dmean= %g  global disp-maximum= %g  rms= %g\n\n",
	     dmean, globmax, sqrt(globdisp2sum / All.TotNumPart));
      fflush(stdout);
    }

  for(i = 0, dispmax = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].Vel[j] = 0;
	  P[i].Pos[j] += fac * P[i].GravAccel[j] * 2.0 / (3 * All.Hubble * All.Hubble);
	  P[i].GravAccel[j] = 0;
	}
    }
#endif



  /* Now assign new timesteps and kick */

#ifdef FLEXSTEPS
  if((All.Ti_Current % (4 * All.PresentMinStep)) == 0)
    if(All.PresentMinStep < TIMEBASE)
      All.PresentMinStep *= 2;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Ti_endstep == All.Ti_Current)
	{
	  ti_step = get_timestep(i, &aphys, 0);

	  /* make it a power 2 subdivision */
	  ti_min = TIMEBASE;
	  while(ti_min > ti_step)
	    ti_min >>= 1;
	  ti_step = ti_min;

	  if(ti_step < All.PresentMinStep)
	    All.PresentMinStep = ti_step;
	}
    }

  ti_step = All.PresentMinStep;
  MPI_Allreduce(&ti_step, &All.PresentMinStep, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif


  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Ti_endstep == All.Ti_Current)
	{
	  ti_step = get_timestep(i, &aphys, 0);

	  /* make it a power 2 subdivision */
	  ti_min = TIMEBASE;
	  while(ti_min > ti_step)
	    ti_min >>= 1;
	  ti_step = ti_min;

#ifdef FLEXSTEPS
	  ti_grp = P[i].FlexStepGrp % All.PresentMaxStep;
	  ti_grp = (ti_grp / All.PresentMinStep) * All.PresentMinStep;
	  ti_step = ((P[i].Ti_endstep + ti_grp + ti_step) / ti_step) * ti_step - (P[i].Ti_endstep + ti_grp);
#else

#ifdef PSEUDOSYMMETRIC
	  if(P[i].Type != 0)
	    {
	      if(P[i].Ti_endstep > P[i].Ti_begstep)
		{
		  apred = aphys + ((aphys - P[i].AphysOld) / (P[i].Ti_endstep - P[i].Ti_begstep)) * ti_step;
		  if(fabs(apred - aphys) < 0.5 * aphys)
		    {
		      ti_step2 = get_timestep(i, &apred, -1);
		      ti_min = TIMEBASE;
		      while(ti_min > ti_step2)
			ti_min >>= 1;
		      ti_step2 = ti_min;

		      if(ti_step2 < ti_step)
			{
			  get_timestep(i, &apred, ti_step);
			  prob =
			    ((apred - aphys) / (aphys - P[i].AphysOld) * (P[i].Ti_endstep -
									  P[i].Ti_begstep)) / ti_step;
			  if(prob < get_random_number(P[i].ID))
			    ti_step /= 2;
			}
		      else if(ti_step2 > ti_step)
			{
			  get_timestep(i, &apred, 2 * ti_step);
			  prob =
			    ((apred - aphys) / (aphys - P[i].AphysOld) * (P[i].Ti_endstep -
									  P[i].Ti_begstep)) / ti_step;
			  if(prob < get_random_number(P[i].ID + 1))
			    ti_step *= 2;
			}
		    }
		}
	      P[i].AphysOld = aphys;
	    }
#endif

#ifdef SYNCHRONIZATION
	  if(ti_step > (P[i].Ti_endstep - P[i].Ti_begstep))	/* timestep wants to increase */
	    {
	      if(((TIMEBASE - P[i].Ti_endstep) % ti_step) > 0)
		ti_step = P[i].Ti_endstep - P[i].Ti_begstep;	/* leave at old step */
	    }
#endif
#endif /* end of FLEXSTEPS */

	  if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	    ti_step = 0;

	  if((TIMEBASE - All.Ti_Current) < ti_step)	/* check that we don't run beyond the end */
	    ti_step = TIMEBASE - All.Ti_Current;

#ifdef TIMESTEP_LIMITER
/*  ------------------------------------------------------------------------------------------------------------------  */
/* 	  Main change in the time intregration module:  */
/* 	  loop over all active gas particles to communicate  */
/* 	  their NEXT timestep to the neighbours */
/*  ------------------------------------------------------------------------------------------------------------------  */

	  P[i].Ti_sizestep = ti_step;
	}
    }
  t0 = second();
  NumSphUpdate = 0;
  
  if(All.NumCurrentTiStep > 0)
    {
      for(n = 0; n < N_gas; n++)
	if(P[n].Type == 0)
	  {
	    if(P[n].Ti_endstep == All.Ti_Current)
	      NumSphUpdate++;
	  }
    }
  t1 = second();
  timecomp += timediff(t0, t1);
  
  t0 = second();
  numlist = (int*)malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);
  t1 = second();
  timeimbalance += timediff(t0, t1);
  
  if(ThisTask == 0)
    {
      fprintf(stdout,"\n      ---> number of timestep active gas particles:    %6lld \n", ntot);
      fflush(stdout);
    }
  
  if(ntot)
    NTotDone = 1;
  else
    NTotDone = 0;
  
  /* loop for time-step limiter */
  while(NTotDone > 0)
    {
      NDone = 0;
      
      i = 0;                        /* first particle for this task */
      ntotleft = ntot;              /* particles left for all tasks together */
      
      while(ntotleft > 0)
	{
	  for(j = 0; j < NTask; j++)
	    nsend_local[j] = 0;
	  
	  /* do local particles and prepare export list */
	  t0 = second();
	  for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeTime - NTask; i++)
	    if(P[i].Type == 0)
	      {
		if(P[i].Ti_endstep == All.Ti_Current)
		  {
		    ndone++;
		    
		    for(j = 0; j < NTask; j++)
		      Exportflag[j] = 0;
		    
		    /* -------------------------------- */
		    NDone += time_limiter_evaluate(i, 0);
		    /* -------------------------------- */
		    
		    for(j = 0; j < NTask; j++)
		      {
			if(Exportflag[j])
			  {
			    TimeDataIn[nexport].Pos[0] =    P[i].Pos[0];
			    TimeDataIn[nexport].Pos[1] =    P[i].Pos[1];
			    TimeDataIn[nexport].Pos[2] =    P[i].Pos[2];
			    TimeDataIn[nexport].Hsml   = SphP[i].Hsml;
			    TimeDataIn[nexport].Size   =    P[i].Ti_sizestep;
			    TimeDataIn[nexport].Begin  =    P[i].Ti_endstep;
			    TimeDataIn[nexport].Index  = i;
			    TimeDataIn[nexport].Task   = j;
			    nexport++;
			    nsend_local[j]++;
			  }
		      }
		  }
	      }
	  t1 = second();
	  timecomp += timediff(t0, t1);
	  
	  qsort(TimeDataIn, nexport, sizeof(struct timedata_in), time_compare_key);
	  
	  for(j = 1, noffset[0] = 0; j < NTask; j++)
	    noffset[j] = noffset[j - 1] + nsend_local[j - 1];
	  
	  t0 = second();
	  MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);
	  t1 = second();
	  timeimbalance += timediff(t0, t1);
	  
	  /* now do the particles that need to be exported */
	  
	  for(level = 1; level < (1 << PTask); level++)
	    {
	      t0 = second();
	      for(j = 0; j < NTask; j++)
		nbuffer[j] = 0;
	      for(ngrp = level; ngrp < (1 << PTask); ngrp++)
		{
		  maxfill = 0;
		  for(j = 0; j < NTask; j++)
		    {
		      if((j ^ ngrp) < NTask)
			if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
			  maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		    }
		  if(maxfill >= All.BunchSizeTime)
		    break;
		  
		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;
		  
		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* get the particles */
			  MPI_Sendrecv(&TimeDataIn[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct timedata_in), MPI_BYTE,
				       recvTask, 0,
				       &TimeDataGet[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct timedata_in),
				       MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, &status);
			}
		    }
		  
		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		}
	      t1 = second();
	      timecommsumm += timediff(t0, t1);
	      
	      t0 = second();
	      /* --------------------------------- */
	      for(j = 0; j < nbuffer[ThisTask]; j++)
		NDone += time_limiter_evaluate(j, 1);
	      /* --------------------------------- */
	      t1 = second();
	      timecomp += timediff(t0, t1);
	      
	      /* do a block to explicitly measure imbalance */
	      t0 = second();
	      MPI_Barrier(MPI_COMM_WORLD);
	      t1 = second();
	      timeimbalance += timediff(t0, t1);
	      
	      level = ngrp - 1;
	    }
	  
	  t0 = second();
	  MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
	  for(j = 0; j < NTask; j++)
	    ntotleft -= ndonelist[j];
	  t1 = second();
	  timeimbalance += timediff(t0, t1);
	}
      
      t0 = second();
      numlist = (int*)malloc(NTask * sizeof(int) * NTask);
      MPI_Allgather(&NDone, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
      for(i = 0, NTotDone = 0; i < NTask; i++)
        NTotDone += numlist[i];
      free(numlist);
      t1 = second();
      timeimbalance += timediff(t0, t1);
      
      if(ThisTask == 0)
	{
	  fprintf(stdout," %3d) ---> number of timestep shrinked gas neighbors:  %6lld \n", CptLimit++, NTotDone);
	  fflush(stdout);
	}
    }
  
  /*   do final operations on results */
  t0 = second();
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Ti_endstep == All.Ti_Current)
	{
	  ti_step =  P[i].Ti_sizestep;
	  
	  /*  ------------------------------------------------------------------------------------------------------------------  */
#endif
	  
	  tstart  = (P[i].Ti_begstep + P[i].Ti_endstep) / 2;	/* midpoint of old step */
	  tend    =  P[i].Ti_endstep + ti_step / 2;	        /* midpoint of new step */
	  
	  if(All.ComovingIntegrationOn)
	    {
	      dt_entr = (tend - tstart) * All.Timebase_interval;
	      dt_entr2 = (tend - P[i].Ti_endstep) * All.Timebase_interval;
	      dt_gravkick = get_gravkick_factor(tstart, tend);
	      dt_hydrokick = get_hydrokick_factor(tstart, tend);
	      dt_gravkick2 = get_gravkick_factor(P[i].Ti_endstep, tend);
	      dt_hydrokick2 = get_hydrokick_factor(P[i].Ti_endstep, tend);
	    }
	  else
	    {
	      dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;
	      dt_gravkick2 = dt_hydrokick2 = dt_entr2 = (tend - P[i].Ti_endstep) * All.Timebase_interval;
	    }

	  P[i].Ti_begstep = P[i].Ti_endstep;
	  P[i].Ti_endstep = P[i].Ti_begstep + ti_step;

	  /* do the kick */

	  for(j = 0; j < 3; j++)
	    {
	      dv[j] = P[i].GravAccel[j] * dt_gravkick;
	      P[i].Vel[j] += dv[j];
	    }

	  if(P[i].Type == 0)	/* SPH stuff */
	    {
	      for(j = 0; j < 3; j++)
		{
		  dv[j] += SphP[i].HydroAccel[j] * dt_hydrokick;
		  P[i].Vel[j] += SphP[i].HydroAccel[j] * dt_hydrokick;

		  SphP[i].VelPred[j] =
		    P[i].Vel[j] - dt_gravkick2 * P[i].GravAccel[j] - dt_hydrokick2 * SphP[i].HydroAccel[j];
#ifdef PMGRID
		  SphP[i].VelPred[j] += P[i].GravPM[j] * dt_gravkickB;
#endif
		}

	      /* In case of cooling, we prevent that the entropy (and
	         hence temperature decreases by more than a factor 0.5 */

	      if(SphP[i].DtEntropy * dt_entr > -0.5 * SphP[i].Entropy)
		SphP[i].Entropy += SphP[i].DtEntropy * dt_entr;
	      else
		SphP[i].Entropy *= 0.5;

	      if(All.MinEgySpec)
		{
		  minentropy = All.MinEgySpec * GAMMA_MINUS1 / pow(SphP[i].Density * a3inv, GAMMA_MINUS1);
		  if(SphP[i].Entropy < minentropy)
		    {
		      SphP[i].Entropy = minentropy;
		      SphP[i].DtEntropy = 0;
		    }
		}

	      /* In case the timestep increases in the new step, we
	         make sure that we do not 'overcool' when deriving
	         predicted temperatures. The maximum timespan over
	         which prediction can occur is ti_step/2, i.e. from
	         the middle to the end of the current step */

	      dt_entr = ti_step / 2 * All.Timebase_interval;
	      if(SphP[i].Entropy + SphP[i].DtEntropy * dt_entr < 0.5 * SphP[i].Entropy)
		SphP[i].DtEntropy = -0.5 * SphP[i].Entropy / dt_entr;
	    }


	  /* if tree is not going to be reconstructed, kick parent nodes dynamically.
	   */
	  if(All.NumForcesSinceLastDomainDecomp < All.TotNumPart * All.TreeDomainUpdateFrequency)
	    {
	      no = Father[i];
	      while(no >= 0)
		{
		  for(j = 0; j < 3; j++)
		    Extnodes[no].vs[j] += dv[j] * P[i].Mass / Nodes[no].u.d.mass;

		  no = Nodes[no].u.d.father;
		}
	    }
	}
    }


  t1 = second();
  timecomp += timediff(t0, t1);

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);


#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
      ti_step = TIMEBASE;
      while(ti_step > (dt_displacement / All.Timebase_interval))
	ti_step >>= 1;

      if(ti_step > (All.PM_Ti_endstep - All.PM_Ti_begstep))	/* PM-timestep wants to increase */
	{
	  /* we only increase if an integer number of steps will bring us to the end */
	  if(((TIMEBASE - All.PM_Ti_endstep) % ti_step) > 0)
	    ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;	/* leave at old step */
	}

      if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	ti_step = 0;

      tstart = (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2;
      tend = All.PM_Ti_endstep + ti_step / 2;

      if(All.ComovingIntegrationOn)
	dt_gravkick = get_gravkick_factor(tstart, tend);
      else
	dt_gravkick = (tend - tstart) * All.Timebase_interval;

      All.PM_Ti_begstep = All.PM_Ti_endstep;
      All.PM_Ti_endstep = All.PM_Ti_begstep + ti_step;

      if(All.ComovingIntegrationOn)
	dt_gravkickB = -get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
      else
	dt_gravkickB =
	  -((All.PM_Ti_begstep + All.PM_Ti_endstep) / 2 - All.PM_Ti_begstep) * All.Timebase_interval;

      for(i = 0; i < NumPart; i++)
	{
	  for(j = 0; j < 3; j++)	/* do the kick */
	    P[i].Vel[j] += P[i].GravPM[j] * dt_gravkick;

	  if(P[i].Type == 0)
	    {
	      if(All.ComovingIntegrationOn)
		{
		  dt_gravkickA = get_gravkick_factor(P[i].Ti_begstep, All.Ti_Current) -
		    get_gravkick_factor(P[i].Ti_begstep, (P[i].Ti_begstep + P[i].Ti_endstep) / 2);
		  dt_hydrokick = get_hydrokick_factor(P[i].Ti_begstep, All.Ti_Current) -
		    get_hydrokick_factor(P[i].Ti_begstep, (P[i].Ti_begstep + P[i].Ti_endstep) / 2);
		}
	      else
		dt_gravkickA = dt_hydrokick =
		  (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;

	      for(j = 0; j < 3; j++)
		SphP[i].VelPred[j] = P[i].Vel[j]
		  + P[i].GravAccel[j] * dt_gravkickA
		  + SphP[i].HydroAccel[j] * dt_hydrokick + P[i].GravPM[j] * dt_gravkickB;
	      
	    }
	}
    }
#endif
  
}



/*! This function normally (for flag==0) returns the maximum allowed timestep
 *  of a particle, expressed in terms of the integer mapping that is used to
 *  represent the total simulated timespan. The physical acceleration is
 *  returned in `aphys'. The latter is used in conjunction with the
 *  PSEUDOSYMMETRIC integration option, which also makes of the second
 *  function of get_timestep. When it is called with a finite timestep for
 *  flag, it returns the physical acceleration that would lead to this
 *  timestep, assuming timestep criterion 0.
 */
int get_timestep(int p,         /*!< particle index */
                 double *aphys, /*!< acceleration (physical units) */
                 int flag       /*!< either 0 for normal operation, or finite timestep to get corresponding
				  aphys */ )
{
  double ax, ay, az, ac, csnd;
  double dt = 0, dt_courant = 0, dt_accel;
  int ti_step;
  
#ifdef CONDUCTION
  double dt_cond;
#endif
  
  if(flag == 0)
    {
      ax = fac1 * P[p].GravAccel[0] * 0;
      ay = fac1 * P[p].GravAccel[1] * 0;
      az = fac1 * P[p].GravAccel[2] * 0;
      
#ifdef PMGRID
      ax += fac1 * P[p].GravPM[0] * 0;
      ay += fac1 * P[p].GravPM[1] * 0;
      az += fac1 * P[p].GravPM[2] * 0;
#endif
      
      if(P[p].Type == 0)
        {
#ifndef TIMESTEP_UPDATE
          ax += fac2 * SphP[p].HydroAccel[0];
          ay += fac2 * SphP[p].HydroAccel[1];
          az += fac2 * SphP[p].HydroAccel[2];
#else
          if(SphP[p].FeedbackFlag)
            {
              ax += fac2 * SphP[p].FeedAccel[0];
              ay += fac2 * SphP[p].FeedAccel[1];
              az += fac2 * SphP[p].FeedAccel[2];
            }
          else
            {
              ax += fac2 * SphP[p].HydroAccel[0];
              ay += fac2 * SphP[p].HydroAccel[1];
              az += fac2 * SphP[p].HydroAccel[2];
            }
#endif
        }
      
      ac = sqrt(ax * ax + ay * ay + az * az);   /* this is now the physical acceleration */
      *aphys = ac;
    }
  else
    ac = *aphys;
  
  if(ac == 0)
    ac = 1.0e-30;
  
  switch (All.TypeOfTimestepCriterion)
    {
    case 0:
      if(flag > 0)
	{
          dt = flag * All.Timebase_interval;
          dt /= hubble_a;       /* convert dloga to physical timestep  */
          ac = 2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / (dt * dt);
          *aphys = ac;
          return flag;
	}
      
      if(P[p].Type == 0)
        dt = dt_accel = sqrt(2 * All.ErrTolIntAccuracy * atime * dmin(SphP[p].Hsml, All.SofteningTable[P[p].Type]) / ac);
      else
        dt = dt_accel = sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac);
      
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(P[p].Type == 0)
        dt = dt_accel = sqrt(2 * All.ErrTolIntAccuracy * atime * SphP[p].Hsml / 2.8 / ac);
#endif
      break;
    default:
      endrun(888);
      break;
    }
  
  if(P[p].Type == 0)
    {
      csnd = sqrt(GAMMA * SphP[p].Pressure / SphP[p].Density);
      
      if(All.ComovingIntegrationOn)
        dt_courant = 2 * All.CourantFac * All.Time * SphP[p].Hsml / (fac3 * SphP[p].MaxSignalVel);
      else
        dt_courant = 2 * All.CourantFac * SphP[p].Hsml / SphP[p].MaxSignalVel;
      
      if(dt_courant < dt)
	dt = dt_courant;
    }
  
  dt *= hubble_a;
  
  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;
  
  if(dt >= dt_displacement)
    dt = dt_displacement;
  
  
  if(dt < All.MinSizeTimestep)
    {
#ifndef NOSTOP_WHEN_BELOW_MINTIMESTEP
      printf("warning: Timestep wants to be below the limit `MinSizeTimestep'\n");
      
      if(P[p].Type == 0)
        {
          printf
            ("Part-ID=%d  dt=%g dtc=%g ac=%g xyz=(%g|%g|%g)  hsml=%g  maxsignalvel=%g dt0=%g eps=%g\n",
             (int) P[p].ID, dt, dt_courant * hubble_a, ac, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2],
             SphP[p].Hsml, SphP[p].MaxSignalVel,
             sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac) * hubble_a,
             All.SofteningTable[P[p].Type]);
        }
      else
        {
          printf("Part-ID=%d  dt=%g ac=%g xyz=(%g|%g|%g)\n", (int) P[p].ID, dt, ac, P[p].Pos[0], P[p].Pos[1],
                 P[p].Pos[2]);
	}
      fflush(stdout);
      endrun(888);
#endif
      dt = All.MinSizeTimestep;
    }
  
  ti_step = dt / All.Timebase_interval;
  
  if(!(ti_step > 0 && ti_step < TIMEBASE) && !ZeroTimestepEncountered)
    {
      printf("\nError: A timestep of size zero was assigned on the integer timeline!\n"
             "We better stop.\n"
             "Task=%d Part-ID=%d dt=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g) tree=(%g|%g%g)\n\n",
             ThisTask, (int) P[p].ID, dt, All.Timebase_interval, ti_step, ac,
             P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], P[p].GravAccel[0], P[p].GravAccel[1], P[p].GravAccel[2]);
#ifdef PMGRID
      printf("pm_force=(%g|%g|%g)\n", P[p].GravPM[0], P[p].GravPM[1], P[p].GravPM[2]);
#endif
      if(P[p].Type == 0)
	printf("hydro-frc=(%g|%g|%g)\n", SphP[p].HydroAccel[0], SphP[p].HydroAccel[1], SphP[p].HydroAccel[2]);
      
      fflush(stdout);
      //endrun(818);
      ZeroTimestepEncountered = 1; // Do not terminate the run, but let AMUSE handle it
    }
  
  return ti_step;
}



/*! This function computes an upper limit ('dt_displacement') to the global
 *  timestep of the system based on the rms velocities of particles. For
 *  cosmological simulations, the criterion used is that the rms displacement
 *  should be at most a fraction MaxRMSDisplacementFac of the mean particle
 *  separation. Note that the latter is estimated using the assigned particle
 *  masses, separately for each particle type. If comoving integration is not
 *  used, the function imposes no constraint on the timestep.
 */
void find_dt_displacement_constraint(double hfac /*!<  should be  a^2*H(a)  */ )
{
  int i, j, type, *temp;
  int count[6];
  long long count_sum[6];
  double v[6], v_sum[6], mim[6], min_mass[6];
  double dt, dmean, asmth = 0;

  dt_displacement = All.MaxSizeTimestep;

  if(All.ComovingIntegrationOn)
    {
      for(type = 0; type < 6; type++)
	{
	  count[type] = 0;
	  v[type] = 0;
	  mim[type] = 1.0e30;
	}

      for(i = 0; i < NumPart; i++)
	{
	  v[P[i].Type] += P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2];
	  if(mim[P[i].Type] > P[i].Mass)
	    mim[P[i].Type] = P[i].Mass;
	  count[P[i].Type]++;
	}

      MPI_Allreduce(v, v_sum, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(mim, min_mass, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      temp = malloc(NTask * 6 * sizeof(int));
      MPI_Allgather(count, 6, MPI_INT, temp, 6, MPI_INT, MPI_COMM_WORLD);
      for(i = 0; i < 6; i++)
	{
	  count_sum[i] = 0;
	  for(j = 0; j < NTask; j++)
	    count_sum[i] += temp[j * 6 + i];
	}
      free(temp);

      for(type = 0; type < 6; type++)
	{
	  if(count_sum[type] > 0)
	    {
	      if(type == 0)
		dmean =
		  pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);
	      else
		dmean =
		  pow(min_mass[type] /
		      ((All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);

	      dt = All.MaxRMSDisplacementFac * hfac * dmean / sqrt(v_sum[type] / count_sum[type]);

#ifdef PMGRID
	      asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
	      if(((1 << type) & (PLACEHIGHRESREGION)))
		asmth = All.Asmth[1];
#endif
	      if(asmth < dmean)
		dt = All.MaxRMSDisplacementFac * hfac * asmth / sqrt(v_sum[type] / count_sum[type]);
#endif

	      if(ThisTask == 0)
		printf("type=%d  dmean=%g asmth=%g minmass=%g a=%g  sqrt(<p^2>)=%g  dlogmax=%g\n",
		       type, dmean, asmth, min_mass[type], All.Time, sqrt(v_sum[type] / count_sum[type]), dt);

	      if(dt < dt_displacement)
		dt_displacement = dt;
	    }
	}

      if(ThisTask == 0)
	printf("displacement time constraint: %g  (%g)\n", dt_displacement, All.MaxSizeTimestep);
    }
}


#ifdef TIMESTEP_LIMITER

int time_limiter_evaluate(int target, int mode)
{
  int j, k, n, startnode, numngb_inbox;
  double h, h2;
  double r2, dx, dy, dz;

  FLOAT *pos;

  int  sizestep_i, begstep_i;
  int  sizestep_j, begstep_j, endstep_j;
  int  sizestep_j_old, begstep_j_old;

  int  tstart, tend;

  int  CptShrink = 0;

  double dt_entr;
  double dt_gravkick;
  double dt_hydrokick;

  FLOAT dv[3];

  if(mode == 0)
    {
      pos        =    P[target].Pos;
      h          = SphP[target].Hsml;
      sizestep_i =    P[target].Ti_sizestep;
      begstep_i  =    P[target].Ti_endstep;
    }
  else
    {
      pos        =  TimeDataGet[target].Pos;
      h          =  TimeDataGet[target].Hsml;
      sizestep_i =  TimeDataGet[target].Size;
      begstep_i  =  TimeDataGet[target].Begin;
    }

  h2   = h * h;

  startnode = All.MaxPart;

  do
    {
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode);

      for(n = 0; n < numngb_inbox; n++)
        {
          j = Ngblist[n];

          dx = P[j].Pos[0] - pos[0];
          dy = P[j].Pos[1] - pos[1];
          dz = P[j].Pos[2] - pos[2];

#ifdef PERIODIC                 /*  now find the closest image in the given box size  */
          if(dx > boxHalf_X)
            dx -= boxSize_X;
          if(dx < -boxHalf_X)
            dx += boxSize_X;
          if(dy > boxHalf_Y)
            dy -= boxSize_Y;
          if(dy < -boxHalf_Y)
            dy += boxSize_Y;
          if(dz > boxHalf_Z)
            dz -= boxSize_Z;
          if(dz < -boxHalf_Z)
            dz += boxSize_Z;
#endif
          r2 = dx * dx + dy * dy + dz * dz;
          if(r2 < h2)
            {
              /* ################################ */
              sizestep_j_old = P[j].Ti_sizestep;
              sizestep_j     = P[j].Ti_sizestep;
              begstep_j_old  = P[j].Ti_begstep;
              begstep_j      = P[j].Ti_begstep;
              endstep_j      = P[j].Ti_endstep;
              /* ################################ */

              if( (endstep_j - All.Ti_Current > TIMESTEP_LIMITER * sizestep_i)
		  || ( (sizestep_j > TIMESTEP_LIMITER * sizestep_i) && endstep_j == All.Ti_Current ) )
                {
                  if( endstep_j == All.Ti_Current ) /* the particle is synchronized */
                    {
                      /* just change the particle time-step. note: this may not work correctly since the
                         particle is active at this time, and may not be able to communicate the new
                         time-step to neighbours. in the end that is enforced in Saitoh's scheme, so we
                         may do better here */

		      CptShrink++;

                      P[j].Ti_sizestep = TIMESTEP_LIMITER * sizestep_i;
                    }
                  else
                    {
		      CptShrink++;

#ifdef SYNCHRONIZATION
                      for(k = 0; begstep_j + k * TIMESTEP_LIMITER * sizestep_i <= All.Ti_Current; k++);
                      
		      endstep_j  = begstep_j + k * TIMESTEP_LIMITER * sizestep_i;
                      sizestep_j = TIMESTEP_LIMITER * sizestep_i;
                      begstep_j  = endstep_j - sizestep_j;
#else
		      endstep_j  = begstep_i + TIMESTEP_LIMITER * sizestep_i;
                      sizestep_j = endstep_j - begstep_j;
#endif

		      /* ###################################################### */
                      /* it's needed to synchronize mid-step quantities to the 
			 new fictious time-step. */
                      tstart = begstep_j_old + sizestep_j_old / 2;
                      tend = (begstep_j + endstep_j) / 2;
		      /* ###################################################### */

                      if(All.ComovingIntegrationOn)
                        {
                          dt_entr = (tend - tstart) * All.Timebase_interval;
                          dt_gravkick = get_gravkick_factor(tstart, tend);
                          dt_hydrokick = get_hydrokick_factor(tstart, tend);
                        }
                      else
			dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;

		      for(k = 0; k < 3; k++)
                        {
                          dv[k] = P[j].GravAccel[k] * dt_gravkick;
                          P[j].Vel[k] += dv[k];
			}

                      if(P[j].Type == 0)        /* SPH stuff */
                        {
                          /* note that VelPred is already at the right time */
                          for(k = 0; k < 3; k++)
                            P[j].Vel[k] += SphP[j].HydroAccel[k] * dt_hydrokick;

                          SphP[j].Entropy += SphP[j].DtEntropy * dt_entr;
                        }
		      /* ######################### */
		      P[j].Ti_begstep  = begstep_j;
		      P[j].Ti_endstep  = endstep_j;
		      P[j].Ti_sizestep = sizestep_j;
		      /* ######################### */
		    }
		}
            }
        }
    }
  while(startnode >= 0);

  return CptShrink;
}





/*! This is a comparison kernel for a sort routine, which is used to group
 *  particles that are going to be exported to the same CPU.
 */
int time_compare_key(const void *a, const void *b)
{
  if(((struct timedata_in *) a)->Task < (((struct timedata_in *) b)->Task))
    return -1;
  if(((struct timedata_in *) a)->Task > (((struct timedata_in *) b)->Task))
    return +1;
  return 0;
}




void make_it_active(int target)
{
  int k;
  int tstart, tend;

  double dt_entr;
  double dt_gravkick;
  double dt_hydrokick;

  FLOAT dv[3];

#ifdef PMGRID
  double dt_gravkick_pm;

  if(All.ComovingIntegrationOn)
    dt_gravkick_pm =
      get_gravkick_factor(All.PM_Ti_begstep, All.Ti_Current) -
      get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
  else
    dt_gravkick_pm = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;
#endif

  tstart = (P[target].Ti_begstep + All.Ti_Current) / 2;       /* midpoint of shrinked step */
  tend = (P[target].Ti_begstep + P[target].Ti_endstep) / 2;   /* midpoint of old step */

  if(All.ComovingIntegrationOn)
    {
      dt_entr = (tend - tstart) * All.Timebase_interval;
      dt_gravkick = get_gravkick_factor(tstart, tend);
      dt_hydrokick = get_hydrokick_factor(tstart, tend);
    }
  else
    dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;

  for(k = 0; k < 3; k++)
    {
      dv[k] = P[target].GravAccel[k] * dt_gravkick +
        SphP[target].HydroAccel[k] * dt_hydrokick;
#ifdef PMGRID
      dv[k] += P[target].GravPM[k] * dt_gravkick_pm;
#endif
      P[target].Vel[k] -= dv[k];
    }

  SphP[target].Entropy -= SphP[target].DtEntropy * dt_entr;

  P[target].Ti_endstep = All.Ti_Current;
  NumForceUpdate++;
}
 
#endif
