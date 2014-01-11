#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#ifndef NOMPI
#include <mpi.h>
#endif

#include "allvars.h"
#include "proto.h"


/*! \file gravtree_forcetest.c
 *  \brief routines for direct summation forces
 *
 *  The code in this file allows to compute checks of the force accuracy by
 *  an independent direct summation computation.  To this end, one can
 *  instruct GADGET in the Makefile to coompute direct summation forces for
 *  a certain random subfraction of particles.
 */



#ifdef FORCETEST

/*! This routine does the test of the gravitational tree force by computing
 *  the force for a random subset of particles with direct summation.
 */
void gravity_forcetest(void)
{
  int ntot, iter = 0, ntotleft, nthis;
  double tstart, tend, timetree = 0;
  int i, j, ndone, ngrp, maxfill, place, ndonetot;

#ifndef NOGRAVITY
  int *noffset, *nbuffer, *nsend, *nsend_local;
  int k, nexport;
  int level, sendTask, recvTask;
  double fac1;
#ifndef NOMPI
  MPI_Status status;
#endif
#endif
  double costtotal, *costtreelist;
  double maxt, sumt, *timetreelist;
  double fac;
  char buf[200];

#ifdef PMGRID
  if(All.PM_Ti_endstep != All.Ti_Current)
    return;
#endif

  if(All.ComovingIntegrationOn)
    set_softenings();		/* set new softening lengths */

  for(i = 0, NumForceUpdate = 0; i < NumPart; i++)
    {
      if(P[i].Ti_endstep == All.Ti_Current)
	{
	  if(get_random_number(P[i].ID) < FORCETEST)
	    {
	      P[i].Ti_endstep = -P[i].Ti_endstep - 1;
	      NumForceUpdate++;
	    }
	}
    }

  /* NumForceUpdate is the number of particles on this processor that want a force update */
#ifndef NOMPI
  MPI_Allreduce(&NumForceUpdate, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  ntot = NumForceUpdate;
#endif

  costtotal = 0;

  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);

  i = 0;			/* beginn with this index */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      iter++;

      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */
      tstart = second();
      for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeForce - NTask; i++)
	if(P[i].Ti_endstep < 0)
	  {
	    ndone++;

	    for(j = 0; j < NTask; j++)
	      Exportflag[j] = 1;
	    Exportflag[ThisTask] = 0;

	    costtotal += force_treeevaluate_direct(i, 0);

	    for(j = 0; j < NTask; j++)
	      {
		if(Exportflag[j])
		  {
		    for(k = 0; k < 3; k++)
		      GravDataGet[nexport].u.Pos[k] = P[i].Pos[k];

#ifdef UNEQUALSOFTENINGS
		    GravDataGet[nexport].Type = P[i].Type;
#endif
		    GravDataGet[nexport].w.OldAcc = P[i].OldAcc;

		    GravDataIndexTable[nexport].Task = j;
		    GravDataIndexTable[nexport].Index = i;
		    GravDataIndexTable[nexport].SortIndex = nexport;

		    nexport++;
		    nsend_local[j]++;
		  }
	      }
	  }
      tend = second();
      timetree += timediff(tstart, tend);

      qsort(GravDataIndexTable, nexport, sizeof(struct gravdata_index), grav_tree_compare_key);

      for(j = 0; j < nexport; j++)
	GravDataIn[j] = GravDataGet[GravDataIndexTable[j].SortIndex];

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

#ifndef NOMPI
      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);
#else
    nsend[0] = nsend_local[i];
#endif
      /* now do the particles that need to be exported */

      for(level = 1; level < (1 << PTask); level++)
	{
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
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&GravDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_DIRECT_A,
				   &GravDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_DIRECT_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }

	  tstart = second();
	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    {
	      costtotal += force_treeevaluate_direct(j, 1);
	    }
	  tend = second();
	  timetree += timediff(tstart, tend);


	  /* get the result */
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
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&GravDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_DIRECT_B,
				   &GravDataOut[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_DIRECT_B, MPI_COMM_WORLD, &status);

		      /* add the result to the particles */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  place = GravDataIndexTable[noffset[recvTask] + j].Index;

			  for(k = 0; k < 3; k++)
			    P[place].GravAccelDirect[k] += GravDataOut[j + noffset[recvTask]].u.Acc[k];
			}
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }

	  level = ngrp - 1;
	}

      MPI_Allreduce(&ndone, &ndonetot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      ntotleft -= ndonetot;
    }

  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);


  /* now add things for comoving integration */

  if(All.ComovingIntegrationOn)
    {
#ifndef PERIODIC
      fac1 = 0.5 * All.Hubble * All.Hubble * All.Omega0 / All.G;

      for(i = 0; i < NumPart; i++)
	if(P[i].Ti_endstep < 0)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccelDirect[j] += fac1 * P[i].Pos[j];
#endif
    }



  /*  muliply by G */

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep < 0)
      for(j = 0; j < 3; j++)
	P[i].GravAccelDirect[j] *= All.G;



  /* Finally, the following factor allows a computation of cosmological simulation
     with vacuum energy in physical coordinates */

  if(All.ComovingIntegrationOn == 0)
    {
      fac1 = All.OmegaLambda * All.Hubble * All.Hubble;

      for(i = 0; i < NumPart; i++)
	if(P[i].Ti_endstep < 0)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccelDirect[j] += fac1 * P[i].Pos[j];
    }

  /* now output the forces to a file */

  for(nthis = 0; nthis < NTask; nthis++)
    {
      if(nthis == ThisTask)
	{
	  sprintf(buf, "%s%s", All.OutputDir, "forcetest.txt");
	  if(!(FdForceTest = fopen(buf, "a")))
	    {
	      printf("error in opening file '%s'\n", buf);
	      endrun(17);
	    }
	  for(i = 0; i < NumPart; i++)
	    if(P[i].Ti_endstep < 0)
	      {
#ifndef PMGRID
		fprintf(FdForceTest, "%d %g %g %g %g %g %g %g %g %g %g %g\n",
			P[i].Type, All.Time, All.Time - TimeOfLastTreeConstruction,
			P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],
			P[i].GravAccelDirect[0], P[i].GravAccelDirect[1], P[i].GravAccelDirect[2],
			P[i].GravAccel[0], P[i].GravAccel[1], P[i].GravAccel[2]);
#else
		fprintf(FdForceTest, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
			P[i].Type, All.Time, All.Time - TimeOfLastTreeConstruction,
			P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],
			P[i].GravAccelDirect[0], P[i].GravAccelDirect[1], P[i].GravAccelDirect[2],
			P[i].GravAccel[0], P[i].GravAccel[1], P[i].GravAccel[2],
			P[i].GravPM[0] + P[i].GravAccel[0],
			P[i].GravPM[1] + P[i].GravAccel[1], P[i].GravPM[2] + P[i].GravAccel[2]);
#endif
	      }
	  fclose(FdForceTest);
	}

#ifndef NOMPI
	      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep < 0)
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;

  /* Now the force computation is finished */



  timetreelist = malloc(sizeof(double) * NTask);
  costtreelist = malloc(sizeof(double) * NTask);

  MPI_Gather(&costtotal, 1, MPI_DOUBLE, costtreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&timetree, 1, MPI_DOUBLE, timetreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fac = NTask / ((double) All.TotNumPart);

      for(i = 0, maxt = timetreelist[0], sumt = 0, costtotal = 0; i < NTask; i++)
	{
	  costtotal += costtreelist[i];

	  if(maxt < timetreelist[i])
	    maxt = timetreelist[i];
	  sumt += timetreelist[i];
	}

      fprintf(FdTimings, "DIRECT Nf= %d    part/sec=%g | %g  ia/part=%g \n", ntot, ntot / (sumt + 1.0e-20),
	      ntot / (maxt * NTask), ((double) (costtotal)) / ntot);
      fprintf(FdTimings, "\n");

      fflush(FdTimings);
    }

  free(costtreelist);
  free(timetreelist);
}

#endif
