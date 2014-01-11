#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef NOMPI
#include <mpi.h>
#endif

#include "allvars.h"
#include "proto.h"



/*! \file domain.c
 *  \brief code for domain decomposition
 *
 *  This file contains the code for the domain decomposition of the
 *  simulation volume.  The domains are constructed from disjoint subsets
 *  of the leaves of a fiducial top-level tree that covers the full
 *  simulation volume. Domain boundaries hence run along tree-node
 *  divisions of a fiducial global BH tree. As a result of this method, the
 *  tree force are in principle strictly independent of the way the domains
 *  are cut. The domain decomposition can be carried out for an arbitrary
 *  number of CPUs. Individual domains are not cubical, but spatially
 *  coherent since the leaves are traversed in a Peano-Hilbert order and
 *  individual domains form segments along this order.  This also ensures
 *  that each domain has a small surface to volume ratio, which minimizes
 *  communication.
 */

#define TOPNODEFACTOR  20.0

#define REDUC_FAC      0.98


/*! toGo[task*NTask + partner] gives the number of particles in task 'task'
 *  that have to go to task 'partner'
 */
static int *toGo, *toGoSph;
static int *local_toGo, *local_toGoSph;
static int *list_NumPart;
static int *list_N_gas;
static int *list_load;
static int *list_loadsph;
static double *list_work;

static long long maxload, maxloadsph;

static struct topnode_exchange
{
  peanokey Startkey;
  int Count;
}
 *toplist, *toplist_local;



/*! This is the main routine for the domain decomposition.  It acts as a
 *  driver routine that allocates various temporary buffers, maps the
 *  particles back onto the periodic box if needed, and then does the
 *  domain decomposition, and a final Peano-Hilbert order of all particles
 *  as a tuning measure.
 */
void domain_Decomposition(void)
{
  double t0, t1;

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)
    {
      All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;
      /* to make sure that we do a domain decomposition before the PM-force is evaluated.
         this is needed to make sure that the particles are wrapped into the box */
    }
#endif

  /* Check whether it is really time for a new domain decomposition */
  if(All.NumForcesSinceLastDomainDecomp > All.TotNumPart * All.TreeDomainUpdateFrequency)
    {
      t0 = second();

#ifdef PERIODIC
      do_box_wrapping();	/* map the particles back onto the box */
#endif
      All.NumForcesSinceLastDomainDecomp = 0;
      TreeReconstructFlag = 1;	/* ensures that new tree will be constructed */

      if(ThisTask == 0)
	{
	  printf("domain decomposition... \n");
	  fflush(stdout);
	}

      Key = malloc(sizeof(peanokey) * All.MaxPart);
      KeySorted = malloc(sizeof(peanokey) * All.MaxPart);

      toGo = malloc(sizeof(int) * NTask * NTask);
      toGoSph = malloc(sizeof(int) * NTask * NTask);
      local_toGo = malloc(sizeof(int) * NTask);
      local_toGoSph = malloc(sizeof(int) * NTask);
      list_NumPart = malloc(sizeof(int) * NTask);
      list_N_gas = malloc(sizeof(int) * NTask);
      list_load = malloc(sizeof(int) * NTask);
      list_loadsph = malloc(sizeof(int) * NTask);
      list_work = malloc(sizeof(double) * NTask);

#ifndef NOMPI
      MPI_Allgather(&NumPart, 1, MPI_INT, list_NumPart, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Allgather(&N_gas, 1, MPI_INT, list_N_gas, 1, MPI_INT, MPI_COMM_WORLD);
#else
      list_NumPart[0] = NumPart;
      list_N_gas[0] = N_gas;
#endif
      maxload = All.MaxPart * REDUC_FAC;
      maxloadsph = All.MaxPartSph * REDUC_FAC;

      domain_decompose();

      free(list_work);
      free(list_loadsph);
      free(list_load);
      free(list_N_gas);
      free(list_NumPart);
      free(local_toGoSph);
      free(local_toGo);
      free(toGoSph);
      free(toGo);


      if(ThisTask == 0)
	{
	  printf("domain decomposition done. \n");
	  fflush(stdout);
	}

      t1 = second();
      All.CPU_Domain += timediff(t0, t1);

#ifdef PEANOHILBERT
      t0 = second();
      peano_hilbert_order();
      t1 = second();
      All.CPU_Peano += timediff(t0, t1);
#endif

      free(KeySorted);
      free(Key);
    }

}



/*! This function carries out the actual domain decomposition for all
 *  particle types. It will try to balance the work-load for each domain,
 *  as estimated based on the P[i]-GravCost values.  The decomposition will
 *  respect the maximum allowed memory-imbalance given by the value of
 *  PartAllocFactor.
 */
void domain_decompose(void)
{
  int i, j, status;
  int ngrp, task, partner, sendcount, recvcount;
  long long sumtogo, sumload;
  int maxload, *temp;
  double sumwork, maxwork;

  for(i = 0; i < 6; i++)
    NtypeLocal[i] = 0;

  for(i = 0; i < NumPart; i++)
    NtypeLocal[P[i].Type]++;

  /* because Ntype[] is of type `long long', we cannot do a simple
   * MPI_Allreduce() to sum the total particle numbers
   */
  temp = malloc(NTask * 6 * sizeof(int));
#ifndef NOMPI
  MPI_Allgather(NtypeLocal, 6, MPI_INT, temp, 6, MPI_INT, MPI_COMM_WORLD);
#else
    for(i = 0; i < 6; i++){temp[i] = NtypeLocal[i];}
#endif
  for(i = 0; i < 6; i++)
    {
      Ntype[i] = 0;
      for(j = 0; j < NTask; j++)
	Ntype[i] += temp[j * 6 + i];
    }
  free(temp);

#ifndef UNEQUALSOFTENINGS
  for(i = 0; i < 6; i++)
    if(Ntype[i] > 0)
      break;

  for(ngrp = i + 1; ngrp < 6; ngrp++)
    {
      if(Ntype[ngrp] > 0)
	if(All.SofteningTable[ngrp] != All.SofteningTable[i])
	  {
	    if(ThisTask == 0)
	      {
		fprintf(stdout, "Code was not compiled with UNEQUALSOFTENINGS, but some of the\n");
		fprintf(stdout, "softening lengths are unequal nevertheless.\n");
		fprintf(stdout, "This is not allowed.\n");
	      }
	    endrun(0);
	  }
    }
#endif


  /* determine global dimensions of domain grid */
  domain_findExtent();

  domain_determineTopTree();

  /* determine cost distribution in domain grid */
  domain_sumCost();

  /* find the split of the domain grid recursively */
  status = domain_findSplit(0, NTask, 0, NTopleaves - 1);
  if(status != 0)
    {
      if(ThisTask == 0)
	printf("\nNo domain decomposition that stays within memory bounds is possible.\n");
      endrun(0);
    }

  /* now try to improve the work-load balance of the split */
  domain_shiftSplit();

  DomainMyStart = DomainStartList[ThisTask];
  DomainMyLast = DomainEndList[ThisTask];

  if(ThisTask == 0)
    {
      sumload = maxload = 0;
      sumwork = maxwork = 0;
      for(i = 0; i < NTask; i++)
	{
	  sumload += list_load[i];
	  sumwork += list_work[i];

	  if(list_load[i] > maxload)
	    maxload = list_load[i];

	  if(list_work[i] > maxwork)
	    maxwork = list_work[i];
	}

      printf("work-load balance=%g   memory-balance=%g\n",
	     maxwork / (sumwork / NTask), maxload / (((double) sumload) / NTask));
    }


  /* determine for each cpu how many particles have to be shifted to other cpus */
  domain_countToGo();

  for(i = 0, sumtogo = 0; i < NTask * NTask; i++)
    sumtogo += toGo[i];

  while(sumtogo > 0)
    {
      if(ThisTask == 0)
	{
	  printf("exchange of %d%09d particles\n", (int) (sumtogo / 1000000000),
		 (int) (sumtogo % 1000000000));
	  fflush(stdout);
	}

      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  for(task = 0; task < NTask; task++)
	    {
	      partner = task ^ ngrp;

	      if(partner < NTask && task < partner)
		{
		  /* treat SPH separately */
		  if(All.TotN_gas > 0)
		    {
		      domain_findExchangeNumbers(task, partner, 1, &sendcount, &recvcount);

		      list_NumPart[task] += recvcount - sendcount;
		      list_NumPart[partner] -= recvcount - sendcount;
		      list_N_gas[task] += recvcount - sendcount;
		      list_N_gas[partner] -= recvcount - sendcount;

		      toGo[task * NTask + partner] -= sendcount;
		      toGo[partner * NTask + task] -= recvcount;
		      toGoSph[task * NTask + partner] -= sendcount;
		      toGoSph[partner * NTask + task] -= recvcount;

		      if(task == ThisTask)	/* actually carry out the exchange */
			domain_exchangeParticles(partner, 1, sendcount, recvcount);
		      if(partner == ThisTask)
			domain_exchangeParticles(task, 1, recvcount, sendcount);
		    }

		  domain_findExchangeNumbers(task, partner, 0, &sendcount, &recvcount);

		  list_NumPart[task] += recvcount - sendcount;
		  list_NumPart[partner] -= recvcount - sendcount;

		  toGo[task * NTask + partner] -= sendcount;
		  toGo[partner * NTask + task] -= recvcount;

		  if(task == ThisTask)	/* actually carry out the exchange */
		    domain_exchangeParticles(partner, 0, sendcount, recvcount);
		  if(partner == ThisTask)
		    domain_exchangeParticles(task, 0, recvcount, sendcount);
		}
	    }
	}

      for(i = 0, sumtogo = 0; i < NTask * NTask; i++)
	sumtogo += toGo[i];
    }
}

/*! This function tries to find a split point in a range of cells in the
 *  domain-grid.  The range of cells starts at 'first', and ends at 'last'
 *  (inclusively). The number of cpus that holds the range is 'ncpu', with
 *  the first cpu given by 'cpustart'. If more than 2 cpus are to be split,
 *  the function calls itself recursively. The division tries to achieve a
 *  best particle-load balance under the constraint that 'maxload' and
 *  'maxloadsph' may not be exceeded, and that each cpu holds at least one
 *  cell from the domaingrid. If such a decomposition cannot be achieved, a
 *  non-zero error code is returned.
 *
 *  After successful completion, DomainMyStart[] and DomainMyLast[] contain
 *  the first and last cell of the domaingrid assigned to the local task
 *  for the given type. Also, DomainTask[] contains for each cell the task
 *  it was assigned to.
 */
int domain_findSplit(int cpustart, int ncpu, int first, int last)
{
  int i, split, ok_left, ok_right;
  long long load, sphload, load_leftOfSplit, sphload_leftOfSplit;
  int ncpu_leftOfSplit;
  double maxAvgLoad_CurrentSplit, maxAvgLoad_NewSplit;


  ncpu_leftOfSplit = ncpu / 2;

  for(i = first, load = 0, sphload = 0; i <= last; i++)
    {
      load += DomainCount[i];
      sphload += DomainCountSph[i];
    }

  split = first + ncpu_leftOfSplit;

  for(i = first, load_leftOfSplit = sphload_leftOfSplit = 0; i < split; i++)
    {
      load_leftOfSplit += DomainCount[i];
      sphload_leftOfSplit += DomainCountSph[i];
    }

  /* find the best split point in terms of work-load balance */

  while(split < last - (ncpu - ncpu_leftOfSplit - 1) && split > 0)
    {
      maxAvgLoad_CurrentSplit =
	dmax(load_leftOfSplit / ncpu_leftOfSplit, (load - load_leftOfSplit) / (ncpu - ncpu_leftOfSplit));

      maxAvgLoad_NewSplit =
	dmax((load_leftOfSplit + DomainCount[split]) / ncpu_leftOfSplit,
	     (load - load_leftOfSplit - DomainCount[split]) / (ncpu - ncpu_leftOfSplit));

      if(maxAvgLoad_NewSplit <= maxAvgLoad_CurrentSplit)
	{
	  load_leftOfSplit += DomainCount[split];
	  sphload_leftOfSplit += DomainCountSph[split];
	  split++;
	}
      else
	break;
    }


  /* we will now have to check whether this solution is possible given the restrictions on the maximum load */

  for(i = first, load_leftOfSplit = 0, sphload_leftOfSplit = 0; i < split; i++)
    {
      load_leftOfSplit += DomainCount[i];
      sphload_leftOfSplit += DomainCountSph[i];
    }

  if(load_leftOfSplit > maxload * ncpu_leftOfSplit ||
     (load - load_leftOfSplit) > maxload * (ncpu - ncpu_leftOfSplit))
    {
      /* we did not find a viable split */
      return -1;
    }

  if(sphload_leftOfSplit > maxloadsph * ncpu_leftOfSplit ||
     (sphload - sphload_leftOfSplit) > maxloadsph * (ncpu - ncpu_leftOfSplit))
    {
      /* we did not find a viable split */
      return -1;
    }

  if(ncpu_leftOfSplit >= 2)
    ok_left = domain_findSplit(cpustart, ncpu_leftOfSplit, first, split - 1);
  else
    ok_left = 0;

  if((ncpu - ncpu_leftOfSplit) >= 2)
    ok_right = domain_findSplit(cpustart + ncpu_leftOfSplit, ncpu - ncpu_leftOfSplit, split, last);
  else
    ok_right = 0;

  if(ok_left == 0 && ok_right == 0)
    {
      /* found a viable split */

      if(ncpu_leftOfSplit == 1)
	{
	  for(i = first; i < split; i++)
	    DomainTask[i] = cpustart;

	  list_load[cpustart] = load_leftOfSplit;
	  list_loadsph[cpustart] = sphload_leftOfSplit;
	  DomainStartList[cpustart] = first;
	  DomainEndList[cpustart] = split - 1;
	}

      if((ncpu - ncpu_leftOfSplit) == 1)
	{
	  for(i = split; i <= last; i++)
	    DomainTask[i] = cpustart + ncpu_leftOfSplit;

	  list_load[cpustart + ncpu_leftOfSplit] = load - load_leftOfSplit;
	  list_loadsph[cpustart + ncpu_leftOfSplit] = sphload - sphload_leftOfSplit;
	  DomainStartList[cpustart + ncpu_leftOfSplit] = split;
	  DomainEndList[cpustart + ncpu_leftOfSplit] = last;
	}

      return 0;
    }

  /* we did not find a viable split */
  return -1;
}



/*! This function tries to improve the domain decomposition found by
 *  domain_findSplit() with respect to work-load balance.  To this end, the
 *  boundaries in the existing domain-split solution (which was found by
 *  trying to balance the particle load) are shifted as long as this leads
 *  to better work-load while still remaining within the allowed
 *  memory-imbalance constraints.
 */
void domain_shiftSplit(void)
{
  int i, task, iter = 0, moved;
  double maxw, newmaxw;

  for(task = 0; task < NTask; task++)
    list_work[task] = 0;

  for(i = 0; i < NTopleaves; i++)
    list_work[DomainTask[i]] += DomainWork[i];

  do
    {
      for(task = 0, moved = 0; task < NTask - 1; task++)
	{
	  maxw = dmax(list_work[task], list_work[task + 1]);

	  if(list_work[task] < list_work[task + 1])
	    {
	      newmaxw = dmax(list_work[task] + DomainWork[DomainStartList[task + 1]],
			     list_work[task + 1] - DomainWork[DomainStartList[task + 1]]);
	      if(newmaxw <= maxw)
		{
		  if(list_load[task] + DomainCount[DomainStartList[task + 1]] <= maxload)
		    {
		      if(list_loadsph[task] + DomainCountSph[DomainStartList[task + 1]] > maxloadsph)
			continue;

		      /* ok, we can move one domain cell from right to left */
		      list_work[task] += DomainWork[DomainStartList[task + 1]];
		      list_load[task] += DomainCount[DomainStartList[task + 1]];
		      list_loadsph[task] += DomainCountSph[DomainStartList[task + 1]];
		      list_work[task + 1] -= DomainWork[DomainStartList[task + 1]];
		      list_load[task + 1] -= DomainCount[DomainStartList[task + 1]];
		      list_loadsph[task + 1] -= DomainCountSph[DomainStartList[task + 1]];

		      DomainTask[DomainStartList[task + 1]] = task;
		      DomainStartList[task + 1] += 1;
		      DomainEndList[task] += 1;

		      moved++;
		    }
		}
	    }
	  else
	    {
	      newmaxw = dmax(list_work[task] - DomainWork[DomainEndList[task]],
			     list_work[task + 1] + DomainWork[DomainEndList[task]]);
	      if(newmaxw <= maxw)
		{
		  if(list_load[task + 1] + DomainCount[DomainEndList[task]] <= maxload)
		    {
		      if(list_loadsph[task + 1] + DomainCountSph[DomainEndList[task]] > maxloadsph)
			continue;

		      /* ok, we can move one domain cell from left to right */
		      list_work[task] -= DomainWork[DomainEndList[task]];
		      list_load[task] -= DomainCount[DomainEndList[task]];
		      list_loadsph[task] -= DomainCountSph[DomainEndList[task]];
		      list_work[task + 1] += DomainWork[DomainEndList[task]];
		      list_load[task + 1] += DomainCount[DomainEndList[task]];
		      list_loadsph[task + 1] += DomainCountSph[DomainEndList[task]];

		      DomainTask[DomainEndList[task]] = task + 1;
		      DomainEndList[task] -= 1;
		      DomainStartList[task + 1] -= 1;

		      moved++;
		    }
		}

	    }
	}

      iter++;
    }
  while(moved > 0 && iter < 10 * NTopleaves);
}


/*! This function counts how many particles have to be exchanged between
 *  two CPUs according to the domain split. If the CPUs are already quite
 *  full and hold data from other CPUs as well, not all the particles may
 *  be exchanged at once. In this case the communication phase has to be
 *  repeated, until enough of the third-party particles have been moved
 *  away such that the decomposition can be completed.
 */
void domain_findExchangeNumbers(int task, int partner, int sphflag, int *send, int *recv)
{
  int numpartA, numpartsphA, ntobesentA, maxsendA, maxsendA_old;
  int numpartB, numpartsphB, ntobesentB, maxsendB, maxsendB_old;

  numpartA = list_NumPart[task];
  numpartsphA = list_N_gas[task];

  numpartB = list_NumPart[partner];
  numpartsphB = list_N_gas[partner];

  if(sphflag == 1)
    {
      ntobesentA = toGoSph[task * NTask + partner];
      ntobesentB = toGoSph[partner * NTask + task];
    }
  else
    {
      ntobesentA = toGo[task * NTask + partner] - toGoSph[task * NTask + partner];
      ntobesentB = toGo[partner * NTask + task] - toGoSph[partner * NTask + task];
    }

  maxsendA = imin(ntobesentA, All.BunchSizeDomain);
  maxsendB = imin(ntobesentB, All.BunchSizeDomain);

  do
    {
      maxsendA_old = maxsendA;
      maxsendB_old = maxsendB;

      maxsendA = imin(All.MaxPart - numpartB + maxsendB, maxsendA);
      maxsendB = imin(All.MaxPart - numpartA + maxsendA, maxsendB);
    }
  while((maxsendA != maxsendA_old) || (maxsendB != maxsendB_old));


  /* now make also sure that there is enough space for SPH particeles */
  if(sphflag == 1)
    {
      do
	{
	  maxsendA_old = maxsendA;
	  maxsendB_old = maxsendB;

	  maxsendA = imin(All.MaxPartSph - numpartsphB + maxsendB, maxsendA);
	  maxsendB = imin(All.MaxPartSph - numpartsphA + maxsendA, maxsendB);
	}
      while((maxsendA != maxsendA_old) || (maxsendB != maxsendB_old));
    }

  *send = maxsendA;
  *recv = maxsendB;
}




/*! This function exchanges particles between two CPUs according to the
 *  domain split. In doing this, the memory boundaries which may restrict
 *  the exhange process are observed.
 */
void domain_exchangeParticles(int partner, int sphflag, int send_count, int recv_count)
{
  int i, no, n, count, rep;
#ifndef NOMPI
  MPI_Status status;
#endif

  for(n = 0, count = 0; count < send_count && n < NumPart; n++)
    {
      if(sphflag)
	{
	  if(P[n].Type != 0)
	    continue;
	}
      else
	{
	  if(P[n].Type == 0)
	    continue;
	}

      no = 0;

      while(TopNodes[no].Daughter >= 0)
	no = TopNodes[no].Daughter + (Key[n] - TopNodes[no].StartKey) / (TopNodes[no].Size / 8);

      no = TopNodes[no].Leaf;

      if(DomainTask[no] == partner)
	{
	  if(sphflag)		/* special reorder routine for SPH particles (need to stay at beginning) */
	    {
	      DomainPartBuf[count] = P[n];	/* copy particle and collect in contiguous memory */
	      DomainKeyBuf[count] = Key[n];
	      DomainSphBuf[count] = SphP[n];

	      P[n] = P[N_gas - 1];
	      P[N_gas - 1] = P[NumPart - 1];

	      Key[n] = Key[N_gas - 1];
	      Key[N_gas - 1] = Key[NumPart - 1];

	      SphP[n] = SphP[N_gas - 1];

	      N_gas--;
	    }
	  else
	    {
	      DomainPartBuf[count] = P[n];	/* copy particle and collect in contiguous memory */
	      DomainKeyBuf[count] = Key[n];
	      P[n] = P[NumPart - 1];
	      Key[n] = Key[NumPart - 1];
	    }

	  count++;
	  NumPart--;
	  n--;
	}
    }

  if(count != send_count)
    {
      printf("Houston, we got a problem...\n");
      printf("ThisTask=%d count=%d send_count=%d\n", ThisTask, count, send_count);
      endrun(88);
    }

  /* transmit */

#ifndef NOMPI
  for(rep = 0; rep < 2; rep++)
    {
      if((rep == 0 && ThisTask < partner) || (rep == 1 && ThisTask > partner))
	{
	  if(send_count > 0)
	    {
	      MPI_Ssend(&DomainPartBuf[0], send_count * sizeof(struct particle_data), MPI_BYTE, partner,
			TAG_PDATA, MPI_COMM_WORLD);

	      MPI_Ssend(&DomainKeyBuf[0], send_count * sizeof(peanokey), MPI_BYTE, partner, TAG_KEY,
			MPI_COMM_WORLD);

	      if(sphflag)
		MPI_Ssend(&DomainSphBuf[0], send_count * sizeof(struct sph_particle_data), MPI_BYTE, partner,
			  TAG_SPHDATA, MPI_COMM_WORLD);
	    }
	}

      if((rep == 1 && ThisTask < partner) || (rep == 0 && ThisTask > partner))
	{
	  if(recv_count > 0)
	    {
	      if(sphflag)
		{
		  if((NumPart - N_gas) > recv_count)
		    {
		      for(i = 0; i < recv_count; i++)
			{
			  P[NumPart + i] = P[N_gas + i];
			  Key[NumPart + i] = Key[N_gas + i];
			}
		    }
		  else
		    {
		      for(i = NumPart - 1; i >= N_gas; i--)
			{
			  P[i + recv_count] = P[i];
			  Key[i + recv_count] = Key[i];
			}
		    }

		  MPI_Recv(&P[N_gas], recv_count * sizeof(struct particle_data), MPI_BYTE, partner, TAG_PDATA,
			   MPI_COMM_WORLD, &status);
		  MPI_Recv(&Key[N_gas], recv_count * sizeof(peanokey), MPI_BYTE, partner, TAG_KEY,
			   MPI_COMM_WORLD, &status);
		  MPI_Recv(&SphP[N_gas], recv_count * sizeof(struct sph_particle_data), MPI_BYTE, partner,
			   TAG_SPHDATA, MPI_COMM_WORLD, &status);

		  N_gas += recv_count;
		}
	      else
		{
		  MPI_Recv(&P[NumPart], recv_count * sizeof(struct particle_data), MPI_BYTE, partner,
			   TAG_PDATA, MPI_COMM_WORLD, &status);
		  MPI_Recv(&Key[NumPart], recv_count * sizeof(peanokey), MPI_BYTE, partner,
			   TAG_KEY, MPI_COMM_WORLD, &status);
		}

	      NumPart += recv_count;
	    }
	}
    }
#endif
}

/*! This function determines how many particles that are currently stored
 *  on the local CPU have to be moved off according to the domain
 *  decomposition.
 */
void domain_countToGo(void)
{
  int n, no;

  for(n = 0; n < NTask; n++)
    {
      local_toGo[n] = 0;
      local_toGoSph[n] = 0;
    }

  for(n = 0; n < NumPart; n++)
    {
      no = 0;

      while(TopNodes[no].Daughter >= 0)
	no = TopNodes[no].Daughter + (Key[n] - TopNodes[no].StartKey) / (TopNodes[no].Size / 8);

      no = TopNodes[no].Leaf;

      if(DomainTask[no] != ThisTask)
	{
	  local_toGo[DomainTask[no]] += 1;
	  if(P[n].Type == 0)
	    local_toGoSph[DomainTask[no]] += 1;
	}
    }
#ifndef NOMPI
  MPI_Allgather(local_toGo, NTask, MPI_INT, toGo, NTask, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(local_toGoSph, NTask, MPI_INT, toGoSph, NTask, MPI_INT, MPI_COMM_WORLD);
#else
    toGo[0] = local_toGo[0];
    toGoSph[0] = local_toGoSph[0];
#endif
}


/*! This function walks the global top tree in order to establish the
 *  number of leaves it has. These leaves are distributed to different
 *  processors.
 */
void domain_walktoptree(int no)
{
  int i;

  if(TopNodes[no].Daughter == -1)
    {
      TopNodes[no].Leaf = NTopleaves;
      NTopleaves++;
    }
  else
    {
      for(i = 0; i < 8; i++)
	domain_walktoptree(TopNodes[no].Daughter + i);
    }
}

/*! This routine bins the particles onto the domain-grid, i.e. it sums up the
 *  total number of particles and the total amount of work in each of the
 *  domain-cells. This information forms the basis for the actual decision on
 *  the adopted domain decomposition.
 */
void domain_sumCost(void)
{
  int i, n, no;
  double *local_DomainWork;
  int *local_DomainCount;
  int *local_DomainCountSph;

  local_DomainWork = malloc(NTopnodes * sizeof(double));
  local_DomainCount = malloc(NTopnodes * sizeof(int));
  local_DomainCountSph = malloc(NTopnodes * sizeof(int));



  NTopleaves = 0;

  domain_walktoptree(0);

  for(i = 0; i < NTopleaves; i++)
    {
      local_DomainWork[i] = 0;
      local_DomainCount[i] = 0;
      local_DomainCountSph[i] = 0;
    }

  if(ThisTask == 0)
    printf("NTopleaves= %d\n", NTopleaves);

  for(n = 0; n < NumPart; n++)
    {
      no = 0;

      while(TopNodes[no].Daughter >= 0)
	no = TopNodes[no].Daughter + (Key[n] - TopNodes[no].StartKey) / (TopNodes[no].Size / 8);

      no = TopNodes[no].Leaf;

      if(P[n].Ti_endstep > P[n].Ti_begstep)
	local_DomainWork[no] += (1.0 + P[n].GravCost) / (P[n].Ti_endstep - P[n].Ti_begstep);
      else
	local_DomainWork[no] += (1.0 + P[n].GravCost);

      local_DomainCount[no] += 1;
      if(P[n].Type == 0)
	local_DomainCountSph[no] += 1;
    }

#ifndef NOMPI
  MPI_Allreduce(local_DomainWork, DomainWork, NTopleaves, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(local_DomainCount, DomainCount, NTopleaves, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(local_DomainCountSph, DomainCountSph, NTopleaves, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
    for(i = 0; i < NTopleaves; i++){
        DomainWork[i] = local_DomainWork[i];
        DomainCount[i] = local_DomainCount[i];
        DomainCountSph[i] = local_DomainCountSph[i];
    }
#endif

  free(local_DomainCountSph);
  free(local_DomainCount);
  free(local_DomainWork);
}


/*! This routine finds the extent of the global domain grid.
 */
void domain_findExtent(void)
{
  int i, j;
  double len, xmin[3], xmax[3], xmin_glob[3], xmax_glob[3];

  /* determine local extension */
  for(j = 0; j < 3; j++)
    {
      xmin[j] = MAX_REAL_NUMBER;
      xmax[j] = -MAX_REAL_NUMBER;
    }

  for(i = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  if(xmin[j] > P[i].Pos[j])
	    xmin[j] = P[i].Pos[j];

	  if(xmax[j] < P[i].Pos[j])
	    xmax[j] = P[i].Pos[j];
	}
    }
#ifndef NOMPI
  MPI_Allreduce(xmin, xmin_glob, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(xmax, xmax_glob, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
    for(j = 0; j < 3; j++){
        xmin_glob[j] = xmin[j];
        xmax_glob[j] = xmax[j];
    }
#endif

  len = 0;
  for(j = 0; j < 3; j++)
    if(xmax_glob[j] - xmin_glob[j] > len)
      len = xmax_glob[j] - xmin_glob[j];

  len *= 1.001;

  for(j = 0; j < 3; j++)
    {
      DomainCenter[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]);
      DomainCorner[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]) - 0.5 * len;
    }

  DomainLen = len;
  DomainFac = 1.0 / len * (((peanokey) 1) << (BITS_PER_DIMENSION));
}


/*! This function constructs the global top-level tree node that is used
 *  for the domain decomposition. This is done by considering the string of
 *  Peano-Hilbert keys for all particles, which is recursively chopped off
 *  in pieces of eight segments until each segment holds at most a certain
 *  number of particles.
 */
void domain_determineTopTree(void)
{
  int i, ntop_local, ntop;
  int *ntopnodelist, *ntopoffset;

  for(i = 0; i < NumPart; i++)
    {
      KeySorted[i] = Key[i] = peano_hilbert_key((P[i].Pos[0] - DomainCorner[0]) * DomainFac,
						(P[i].Pos[1] - DomainCorner[1]) * DomainFac,
						(P[i].Pos[2] - DomainCorner[2]) * DomainFac,
						BITS_PER_DIMENSION);
    }

  qsort(KeySorted, NumPart, sizeof(peanokey), domain_compare_key);

  NTopnodes = 1;
  TopNodes[0].Daughter = -1;
  TopNodes[0].Size = PEANOCELLS;
  TopNodes[0].StartKey = 0;
  TopNodes[0].Count = NumPart;
  TopNodes[0].Pstart = 0;

  domain_topsplit_local(0, 0);

  toplist_local = malloc(NTopnodes * sizeof(struct topnode_exchange));

  for(i = 0, ntop_local = 0; i < NTopnodes; i++)
    {
      if(TopNodes[i].Daughter == -1)	/* only use leaves */
	{
	  toplist_local[ntop_local].Startkey = TopNodes[i].StartKey;
	  toplist_local[ntop_local].Count = TopNodes[i].Count;
	  ntop_local++;
	}
    }

  ntopnodelist = malloc(sizeof(int) * NTask);
  ntopoffset = malloc(sizeof(int) * NTask);
#ifndef NOMPI
  MPI_Allgather(&ntop_local, 1, MPI_INT, ntopnodelist, 1, MPI_INT, MPI_COMM_WORLD);
#else
   ntopnodelist[0] = ntop_local;
#endif
  for(i = 0, ntop = 0, ntopoffset[0] = 0; i < NTask; i++)
    {
      ntop += ntopnodelist[i];
      if(i > 0)
	ntopoffset[i] = ntopoffset[i - 1] + ntopnodelist[i - 1];
    }


  toplist = malloc(ntop * sizeof(struct topnode_exchange));

  for(i = 0; i < NTask; i++)
    {
      ntopnodelist[i] *= sizeof(struct topnode_exchange);
      ntopoffset[i] *= sizeof(struct topnode_exchange);
    }

#ifndef NOMPI
  MPI_Allgatherv(toplist_local, ntop_local * sizeof(struct topnode_exchange), MPI_BYTE,
		 toplist, ntopnodelist, ntopoffset, MPI_BYTE, MPI_COMM_WORLD);
#else
   for(i = 0; i < ntop; i++) {
    toplist[i].Count = toplist_local[i].Count;
    toplist[i].Startkey = toplist_local[i].Startkey;
    }
#endif

  qsort(toplist, ntop, sizeof(struct topnode_exchange), domain_compare_toplist);

  NTopnodes = 1;
  TopNodes[0].Daughter = -1;
  TopNodes[0].Size = PEANOCELLS;
  TopNodes[0].StartKey = 0;
  TopNodes[0].Count = All.TotNumPart;
  TopNodes[0].Pstart = 0;
  TopNodes[0].Blocks = ntop;

  domain_topsplit(0, 0);

  free(toplist);
  free(ntopoffset);
  free(ntopnodelist);
  free(toplist_local);

}



/*! This function is responsible for constructing the local top-level
 *  Peano-Hilbert segments. A segment is cut into 8 pieces recursively
 *  until the number of particles in the segment has fallen below
 *  All.TotNumPart / (TOPNODEFACTOR * NTask * NTask).
 */
void domain_topsplit_local(int node, peanokey startkey)
{
  int i, p, sub, bin;

  if(TopNodes[node].Size >= 8)
    {
      TopNodes[node].Daughter = NTopnodes;

      for(i = 0; i < 8; i++)
	{
	  if(NTopnodes < MAXTOPNODES)
	    {
	      sub = TopNodes[node].Daughter + i;
	      TopNodes[sub].Size = TopNodes[node].Size / 8;
	      TopNodes[sub].Count = 0;
	      TopNodes[sub].Daughter = -1;
	      TopNodes[sub].StartKey = startkey + i * TopNodes[sub].Size;
	      TopNodes[sub].Pstart = TopNodes[node].Pstart;

	      NTopnodes++;
	    }
	  else
	    {
	      printf("task=%d: We are out of Topnodes. Increasing the constant MAXTOPNODES might help.\n",
		     ThisTask);
	      fflush(stdout);
	      endrun(13213);
	    }
	}

      for(p = TopNodes[node].Pstart; p < TopNodes[node].Pstart + TopNodes[node].Count; p++)
	{
	  bin = (KeySorted[p] - startkey) / (TopNodes[node].Size / 8);

	  if(bin < 0 || bin > 7)
	    {
	      printf("task=%d: something odd has happened here. bin=%d\n", ThisTask, bin);
	      fflush(stdout);
	      endrun(13123123);
	    }

	  sub = TopNodes[node].Daughter + bin;

	  if(TopNodes[sub].Count == 0)
	    TopNodes[sub].Pstart = p;

	  TopNodes[sub].Count++;
	}

      for(i = 0; i < 8; i++)
	{
	  sub = TopNodes[node].Daughter + i;
	  if(TopNodes[sub].Count > All.TotNumPart / (TOPNODEFACTOR * NTask * NTask))
	    domain_topsplit_local(sub, TopNodes[sub].StartKey);
	}
    }
}



/*! This function is responsible for constructing the global top-level tree
 *  segments. Starting from a joint list of all local top-level segments,
 *  in which mulitple occurences of the same spatial segment have been
 *  combined, a segment is subdivided into 8 pieces recursively until the
 *  number of particles in each segment has fallen below All.TotNumPart /
 *  (TOPNODEFACTOR * NTask).
 */
void domain_topsplit(int node, peanokey startkey)
{
  int i, p, sub, bin;

  if(TopNodes[node].Size >= 8)
    {
      TopNodes[node].Daughter = NTopnodes;

      for(i = 0; i < 8; i++)
	{
	  if(NTopnodes < MAXTOPNODES)
	    {
	      sub = TopNodes[node].Daughter + i;
	      TopNodes[sub].Size = TopNodes[node].Size / 8;
	      TopNodes[sub].Count = 0;
	      TopNodes[sub].Blocks = 0;
	      TopNodes[sub].Daughter = -1;
	      TopNodes[sub].StartKey = startkey + i * TopNodes[sub].Size;
	      TopNodes[sub].Pstart = TopNodes[node].Pstart;
	      NTopnodes++;
	    }
	  else
	    {
	      printf("Task=%d: We are out of Topnodes. Increasing the constant MAXTOPNODES might help.\n",
		     ThisTask);
	      fflush(stdout);
	      endrun(137213);
	    }
	}

      for(p = TopNodes[node].Pstart; p < TopNodes[node].Pstart + TopNodes[node].Blocks; p++)
	{
	  bin = (toplist[p].Startkey - startkey) / (TopNodes[node].Size / 8);
	  sub = TopNodes[node].Daughter + bin;

	  if(bin < 0 || bin > 7)
	    endrun(77);

	  if(TopNodes[sub].Blocks == 0)
	    TopNodes[sub].Pstart = p;

	  TopNodes[sub].Count += toplist[p].Count;
	  TopNodes[sub].Blocks++;
	}

      for(i = 0; i < 8; i++)
	{
	  sub = TopNodes[node].Daughter + i;
	  if(TopNodes[sub].Count > All.TotNumPart / (TOPNODEFACTOR * NTask))
	    domain_topsplit(sub, TopNodes[sub].StartKey);
	}
    }
}


/*! This is a comparison kernel used in a sort routine.
 */
int domain_compare_toplist(const void *a, const void *b)
{
  if(((struct topnode_exchange *) a)->Startkey < (((struct topnode_exchange *) b)->Startkey))
    return -1;

  if(((struct topnode_exchange *) a)->Startkey > (((struct topnode_exchange *) b)->Startkey))
    return +1;

  return 0;
}

/*! This is a comparison kernel used in a sort routine.
 */
int domain_compare_key(const void *a, const void *b)
{
  if(*(peanokey *) a < *(peanokey *) b)
    return -1;

  if(*(peanokey *) a > *(peanokey *) b)
    return +1;

  return 0;
}
