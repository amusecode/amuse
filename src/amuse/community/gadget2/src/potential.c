#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>


#include "allvars.h"
#include "proto.h"


/*! \file potential.c 
 *  \brief Computation of the gravitational potential of particles
 */


/*! This function computes the gravitational potential for ALL the particles.
 *  First, the (short-range) tree potential is computed, and then, if needed,
 *  the long range PM potential is added.
 */
void compute_potential(void)
{
  int i;

#ifndef NOGRAVITY
  long long ntot, ntotleft;
  int j, k, level, sendTask, recvTask;
  int ndone;
  int maxfill, ngrp, place, nexport;
  int *nsend, *noffset, *nsend_local, *nbuffer, *ndonelist, *numlist;
  double fac;
  double t0, t1, tstart, tend;
  MPI_Status status;
  double r2;

  t0 = second();

  if(All.ComovingIntegrationOn)
    set_softenings();

  if(ThisTask == 0)
    {
      printf("Start computation of potential for all particles...\n");
      fflush(stdout);
    }


  tstart = second();
  if(TreeReconstructFlag)
    {
      if(ThisTask == 0)
	printf("Tree construction.\n");

      force_treebuild(NumPart);

      TreeReconstructFlag = 0;

      if(ThisTask == 0)
	printf("Tree construction done.\n");
    }
  tend = second();
  All.CPU_TreeConstruction += timediff(tstart, tend);

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumPart, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);

  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  i = 0;			/* beginn with this index */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */
      for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeForce - NTask; i++)
	{
	  ndone++;

	  for(j = 0; j < NTask; j++)
	    Exportflag[j] = 0;

#ifndef PMGRID
	  force_treeevaluate_potential(i, 0);
#else
	  force_treeevaluate_potential_shortrange(i, 0);
#endif

	  for(j = 0; j < NTask; j++)
	    {
	      if(Exportflag[j])
		{
		  for(k = 0; k < 3; k++)
		    GravDataGet[nexport].u.Pos[k] = P[i].Pos[k];
#ifdef UNEQUALSOFTENINGS
		  GravDataGet[nexport].Type = P[i].Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
                  if(P[i].Type == 0)
                    GravDataGet[nexport].Soft = SphP[i].Hsml;
#endif
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

      qsort(GravDataIndexTable, nexport, sizeof(struct gravdata_index), grav_tree_compare_key);

      for(j = 0; j < nexport; j++)
	GravDataIn[j] = GravDataGet[GravDataIndexTable[j].SortIndex];

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

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
				   recvTask, TAG_POTENTIAL_A,
				   &GravDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_POTENTIAL_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }

	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    {
#ifndef PMGRID
	      force_treeevaluate_potential(j, 1);
#else
	      force_treeevaluate_potential_shortrange(j, 1);
#endif
	    }


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
				   MPI_BYTE, recvTask, TAG_POTENTIAL_B,
				   &GravDataOut[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_POTENTIAL_B, MPI_COMM_WORLD, &status);

		      /* add the result to the particles */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  place = GravDataIndexTable[noffset[recvTask] + j].Index;

			  P[place].Potential += GravDataOut[j + noffset[recvTask]].u.Potential;
			}
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }

	  level = ngrp - 1;
	}

      MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
      for(j = 0; j < NTask; j++)
	ntotleft -= ndonelist[j];
    }

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);


  /* add correction to exclude self-potential */

  for(i = 0; i < NumPart; i++)
    {
      /* remove self-potential */
      P[i].Potential += P[i].Mass / All.SofteningTable[P[i].Type];

      if(All.ComovingIntegrationOn)
	if(All.PeriodicBoundariesOn)
	  P[i].Potential -= 2.8372975 * pow(P[i].Mass, 2.0 / 3) *
	    pow(All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G), 1.0 / 3);
    }


  /* multiply with the gravitational constant */

  for(i = 0; i < NumPart; i++)
    P[i].Potential *= All.G;


#ifdef PMGRID

#ifdef PERIODIC
  pmpotential_periodic();
#ifdef PLACEHIGHRESREGION
  i = pmpotential_nonperiodic(1);
  if(i == 1)  /* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmpotential_nonperiodic(1);	/* try again */
    }
  if(i == 1)
    endrun(88686);
#endif
#else
  i = pmpotential_nonperiodic(0);
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmpotential_nonperiodic(0);	/* try again */
    }
  if(i == 1)
    endrun(88687);
#ifdef PLACEHIGHRESREGION
  i = pmpotential_nonperiodic(1);
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();

      i = pmpotential_nonperiodic(1);
    }
  if(i != 0)
    endrun(88688);
#endif
#endif

#endif



  if(All.ComovingIntegrationOn)
    {
#ifndef PERIODIC
      fac = -0.5 * All.Omega0 * All.Hubble * All.Hubble;

      for(i = 0; i < NumPart; i++)
	{
	  for(k = 0, r2 = 0; k < 3; k++)
	    r2 += P[i].Pos[k] * P[i].Pos[k];

	  P[i].Potential += fac * r2;
	}
#endif
    }
  else
    {
      fac = -0.5 * All.OmegaLambda * All.Hubble * All.Hubble;
      if(fac != 0)
	{
	  for(i = 0; i < NumPart; i++)
	    {
	      for(k = 0, r2 = 0; k < 3; k++)
		r2 += P[i].Pos[k] * P[i].Pos[k];

	      P[i].Potential += fac * r2;
	    }
	}
    }


  if(ThisTask == 0)
    {
      printf("potential done.\n");
      fflush(stdout);
    }

  t1 = second();

  All.CPU_Potential += timediff(t0, t1);

#else
  for(i = 0; i < NumPart; i++)
    P[i].Potential = 0;
#endif
}
