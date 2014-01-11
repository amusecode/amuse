#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef NOMPI
#include <mpi.h>
#endif

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/file.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"

/*! \file restart.c
 *  \brief Code for reading and writing restart files
 */

static FILE *fd;

static void in(int *x, int modus);
static void byten(void *x, size_t n, int modus);


/*! This function reads or writes the restart files.  Each processor writes
 *  its own restart file, with the I/O being done in parallel. To avoid
 *  congestion of the disks you can tell the program to restrict the number of
 *  files that are simultaneously written to NumFilesWrittenInParallel.
 *
 *  If modus>0 the restart()-routine reads, if modus==0 it writes a restart
 *  file.
 */
void restart(int modus)
{
  char buf[200], buf_bak[200], buf_mv[500];
  double save_PartAllocFactor, save_TreeAllocFactor;
  int i, nprocgroup, masterTask, groupTask, old_MaxPart, old_MaxNodes;
  struct global_data_all_processes all_task0;


  sprintf(buf, "%s%s.%d", All.OutputDir, All.RestartFile, ThisTask);
  sprintf(buf_bak, "%s%s.%d.bak", All.OutputDir, All.RestartFile, ThisTask);
  sprintf(buf_mv, "mv %s %s", buf, buf_bak);


  if((NTask < All.NumFilesWrittenInParallel))
    {
      printf
	("Fatal error.\nNumber of processors must be a smaller or equal than `NumFilesWrittenInParallel'.\n");
      endrun(2131);
    }

  nprocgroup = NTask / All.NumFilesWrittenInParallel;

  if((NTask % All.NumFilesWrittenInParallel))
    {
      nprocgroup++;
    }

  masterTask = (ThisTask / nprocgroup) * nprocgroup;

  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))	/* ok, it's this processor's turn */
	{
	  if(modus)
	    {
	      if(!(fd = fopen(buf, "r")))
		{
		  printf("Restart file '%s' not found.\n", buf);
		  endrun(7870);
		}
	    }
	  else
	    {
	      system(buf_mv);	/* move old restart files to .bak files */

	      if(!(fd = fopen(buf, "w")))
		{
		  printf("Restart file '%s' cannot be opened.\n", buf);
		  endrun(7878);
		}
	    }


	  save_PartAllocFactor = All.PartAllocFactor;
	  save_TreeAllocFactor = All.TreeAllocFactor;

	  /* common data  */
	  byten(&All, sizeof(struct global_data_all_processes), modus);

	  if(ThisTask == 0 && modus > 0)
	    all_task0 = All;

#ifndef NOMPI
	  if(modus > 0 && groupTask == 0)	/* read */
	    {
	      MPI_Bcast(&all_task0, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);
	    }
#endif
	  old_MaxPart = All.MaxPart;
	  old_MaxNodes = All.TreeAllocFactor * All.MaxPart;

	  if(modus)		/* read */
	    {
	      if(All.PartAllocFactor != save_PartAllocFactor)
		{
		  All.PartAllocFactor = save_PartAllocFactor;
		  All.MaxPart = All.PartAllocFactor * (All.TotNumPart / NTask);
		  All.MaxPartSph = All.PartAllocFactor * (All.TotN_gas / NTask);
		  save_PartAllocFactor = -1;
		}

	      if(All.TreeAllocFactor != save_TreeAllocFactor)
		{
		  All.TreeAllocFactor = save_TreeAllocFactor;
		  save_TreeAllocFactor = -1;
		}

	      if(all_task0.Time != All.Time)
		{
		  printf("The restart file on task=%d is not consistent with the one on task=0\n", ThisTask);
		  fflush(stdout);
		  endrun(16);
		}

	      allocate_memory();
	    }

	  in(&NumPart, modus);

	  if(NumPart > All.MaxPart)
	    {
	      printf
		("it seems you have reduced(!) 'PartAllocFactor' below the value of %g needed to load the restart file.\n",
		 NumPart / (((double) All.TotNumPart) / NTask));
	      printf("fatal error\n");
	      endrun(22);
	    }

	  /* Particle data  */
	  byten(&P[0], NumPart * sizeof(struct particle_data), modus);

	  in(&N_gas, modus);

	  if(N_gas > 0)
	    {
	      if(N_gas > All.MaxPartSph)
		{
		  printf
		    ("SPH: it seems you have reduced(!) 'PartAllocFactor' below the value of %g needed to load the restart file.\n",
		     N_gas / (((double) All.TotN_gas) / NTask));
		  printf("fatal error\n");
		  endrun(222);
		}
	      /* Sph-Particle data  */
	      byten(&SphP[0], N_gas * sizeof(struct sph_particle_data), modus);
	    }

	  /* write state of random number generator */
	  byten(gsl_rng_state(random_generator), gsl_rng_size(random_generator), modus);


	  /* now store relevant data for tree */

	  if(modus)		/* read */
	    {
	      ngb_treeallocate(MAX_NGB);

	      force_treeallocate(All.TreeAllocFactor * All.MaxPart, All.MaxPart);
	    }


	  in(&Numnodestree, modus);

	  if(Numnodestree > MaxNodes)
	    {
	      printf
		("Tree storage: it seems you have reduced(!) 'PartAllocFactor' below the value needed to load the restart file (task=%d). "
		 "Numnodestree=%d  MaxNodes=%d\n", ThisTask, Numnodestree, MaxNodes);
	      endrun(221);
	    }

	  byten(Nodes_base, Numnodestree * sizeof(struct NODE), modus);
	  byten(Extnodes_base, Numnodestree * sizeof(struct extNODE), modus);

	  byten(Father, NumPart * sizeof(int), modus);

	  byten(Nextnode, NumPart * sizeof(int), modus);
	  byten(Nextnode + All.MaxPart, MAXTOPNODES * sizeof(int), modus);

	  byten(DomainStartList, NTask * sizeof(int), modus);
	  byten(DomainEndList, NTask * sizeof(int), modus);
	  byten(DomainTask, MAXTOPNODES * sizeof(int), modus);
	  byten(DomainNodeIndex, MAXTOPNODES * sizeof(int), modus);
	  byten(DomainTreeNodeLen, MAXTOPNODES * sizeof(FLOAT), modus);
	  byten(DomainHmax, MAXTOPNODES * sizeof(FLOAT), modus);
	  byten(DomainMoment, MAXTOPNODES * sizeof(struct DomainNODE), modus);

	  byten(DomainCorner, 3 * sizeof(double), modus);
	  byten(DomainCenter, 3 * sizeof(double), modus);
	  byten(&DomainLen, sizeof(double), modus);
	  byten(&DomainFac, sizeof(double), modus);
	  byten(&DomainMyStart, sizeof(int), modus);
	  byten(&DomainMyLast, sizeof(int), modus);

	  if(modus)		/* read */
	    if(All.PartAllocFactor != save_PartAllocFactor || All.TreeAllocFactor != save_TreeAllocFactor)
	      {
		for(i = 0; i < NumPart; i++)
		  Father[i] += (All.MaxPart - old_MaxPart);

		for(i = 0; i < NumPart; i++)
		  if(Nextnode[i] >= old_MaxPart)
		    {
		      if(Nextnode[i] >= old_MaxPart + old_MaxNodes)
			Nextnode[i] += (All.MaxPart - old_MaxPart) + (MaxNodes - old_MaxPart);
		      else
			Nextnode[i] += (All.MaxPart - old_MaxPart);
		    }

		for(i = 0; i < Numnodestree; i++)
		  {
		    if(Nodes_base[i].u.d.sibling >= old_MaxPart)
		      {
			if(Nodes_base[i].u.d.sibling >= old_MaxPart + old_MaxNodes)
			  Nodes_base[i].u.d.sibling +=
			    (All.MaxPart - old_MaxPart) + (MaxNodes - old_MaxNodes);
			else
			  Nodes_base[i].u.d.sibling += (All.MaxPart - old_MaxPart);
		      }

		    if(Nodes_base[i].u.d.father >= old_MaxPart)
		      {
			if(Nodes_base[i].u.d.father >= old_MaxPart + old_MaxNodes)
			  Nodes_base[i].u.d.father += (All.MaxPart - old_MaxPart) + (MaxNodes - old_MaxNodes);
			else
			  Nodes_base[i].u.d.father += (All.MaxPart - old_MaxPart);
		      }

		    if(Nodes_base[i].u.d.nextnode >= old_MaxPart)
		      {
			if(Nodes_base[i].u.d.nextnode >= old_MaxPart + old_MaxNodes)
			  Nodes_base[i].u.d.nextnode +=
			    (All.MaxPart - old_MaxPart) + (MaxNodes - old_MaxNodes);
			else
			  Nodes_base[i].u.d.nextnode += (All.MaxPart - old_MaxPart);
		      }
		  }

		for(i = 0; i < MAXTOPNODES; i++)
		  if(Nextnode[i + All.MaxPart] >= old_MaxPart)
		    {
		      if(Nextnode[i + All.MaxPart] >= old_MaxPart + old_MaxNodes)
			Nextnode[i + All.MaxPart] += (All.MaxPart - old_MaxPart) + (MaxNodes - old_MaxNodes);
		      else
			Nextnode[i + All.MaxPart] += (All.MaxPart - old_MaxPart);
		    }

		for(i = 0; i < MAXTOPNODES; i++)
		  if(DomainNodeIndex[i] >= old_MaxPart)
		    {
		      if(DomainNodeIndex[i] >= old_MaxPart + old_MaxNodes)
			DomainNodeIndex[i] += (All.MaxPart - old_MaxPart) + (MaxNodes - old_MaxNodes);
		      else
			DomainNodeIndex[i] += (All.MaxPart - old_MaxPart);
		    }
	      }

	  fclose(fd);
	}
      else			/* wait inside the group */
	{
#ifndef NOMPI
	  if(modus > 0 && groupTask == 0)	/* read */
	    {
	      MPI_Bcast(&all_task0, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);
	    }
#endif
	}

#ifndef NOMPI
	      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
}



/*! reads/writes n bytes in restart routine
 */
void byten(void *x, size_t n, int modus)
{
  if(modus)
    my_fread(x, n, 1, fd);
  else
    my_fwrite(x, n, 1, fd);
}


/*! reads/writes one `int' variable in restart routine
 */
void in(int *x, int modus)
{
  if(modus)
    my_fread(x, 1, sizeof(int), fd);
  else
    my_fwrite(x, 1, sizeof(int), fd);
}
