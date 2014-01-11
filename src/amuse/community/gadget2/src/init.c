#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef NOMPI
#include <mpi.h>
#endif

#include "allvars.h"
#include "proto.h"


/*! \file init.c
 *  \brief Code for initialisation of a simulation from initial conditions
 */


/*! This function reads the initial conditions, and allocates storage for the
 *  tree. Various variables of the particle data are initialised and An intial
 *  domain decomposition is performed. If SPH particles are present, the inial
 *  SPH smoothing lengths are determined.
 */
void init(void)
{
  int i, j;
  double a3;

  All.Time = All.TimeBegin;

  switch (All.ICFormat)
    {
    case 1:
#if (MAKEGLASS > 1)
      seed_glass();
#else
      read_ic(All.InitCondFile);
#endif
      break;
    case 2:
    case 3:
      read_ic(All.InitCondFile);
      break;
    default:
      if(ThisTask == 0)
	printf("ICFormat=%d not supported.\n", All.ICFormat);
      endrun(0);
    }

  All.Time = All.TimeBegin;
  All.Ti_Current = 0;

  if(All.ComovingIntegrationOn)
    {
      All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
      a3 = All.Time * All.Time * All.Time;
    }
  else
    {
      All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
      a3 = 1;
    }

  set_softenings();

  All.NumCurrentTiStep = 0;	/* setup some counters */
  All.SnapshotFileCount = 0;
  if(RestartFlag == 2)
    All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 3) + 1;

  All.TotNumOfForces = 0;
  All.NumForcesSinceLastDomainDecomp = 0;

  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn == 1)
      check_omega();

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;

  if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
    {
      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].Vel[j] *= sqrt(All.Time) * All.Time;
    }

  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] = 0;
#ifdef PMGRID
      for(j = 0; j < 3; j++)
	P[i].GravPM[j] = 0;
#endif
      P[i].Ti_endstep = 0;
      P[i].Ti_begstep = 0;

      P[i].OldAcc = 0;
      P[i].GravCost = 1;
      P[i].Potential = 0;
    }

#ifdef PMGRID
  All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif

#ifdef FLEXSTEPS
  All.PresentMinStep = TIMEBASE;
  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      P[i].FlexStepGrp = (int) (TIMEBASE * get_random_number(P[i].ID));
    }
#endif


  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      for(j = 0; j < 3; j++)
	{
	  SphP[i].VelPred[j] = P[i].Vel[j];
	  SphP[i].HydroAccel[j] = 0;
	}

#ifdef TIMESTEP_UPDATE
  SphP[i].FeedbackFlag = 0;

  for(j = 0; j < 3; j++)
      SphP[i].FeedAccel[j] = 0;
#endif

      SphP[i].DtEntropy = 0;

      if(RestartFlag == 0)
	{
	  SphP[i].Hsml = 0;
	  SphP[i].Density = -1;
	}
    }

  ngb_treeallocate(MAX_NGB);

  force_treeallocate(All.TreeAllocFactor * All.MaxPart, All.MaxPart);

  All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;

  Flag_FullStep = 1;		/* to ensure that Peano-Hilber order is done */

  domain_Decomposition();	/* do initial domain decomposition (gives equal numbers of particles) */

  ngb_treebuild();		/* will build tree */

  setup_smoothinglengths();

  TreeReconstructFlag = 1;

  /* at this point, the entropy variable normally contains the
   * internal energy, read in from the initial conditions file, unless the file
   * explicitly signals that the initial conditions contain the entropy directly.
   * Once the density has been computed, we can convert thermal energy to entropy.
   */
#ifndef ISOTHERM_EQS
  if(header.flag_entropy_instead_u == 0)
    for(i = 0; i < N_gas; i++)
      SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].Density / a3, GAMMA_MINUS1);
#endif
}


/*! This routine computes the mass content of the box and compares it to the
 *  specified value of Omega-matter.  If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass = 0, masstot, omega;
  int i;

  for(i = 0; i < NumPart; i++)
    mass += P[i].Mass;

#ifndef NOMPI
  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  masstot = mass;
#endif
  omega =
    masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));

  if(fabs(omega - All.Omega0) > 1.0e-3)
    {
      if(ThisTask == 0)
	{
	  printf("\n\nI've found something odd!\n");
	  printf
	    ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
	     omega, All.Omega0);
	  printf("\nI better stop.\n");

	  fflush(stdout);
	}
      endrun(1);
    }
}



/*! This function is used to find an initial smoothing length for each SPH
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
void setup_smoothinglengths(void)
{
  int i, no, p;

  if(RestartFlag == 0)
    {

      for(i = 0; i < N_gas; i++)
	{
	  no = Father[i];

	  while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }
#ifndef TWODIMS
	  SphP[i].Hsml =
	    pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#else
	  SphP[i].Hsml =
	    pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif
	}
    }

  density();
}


/*! If the code is run in glass-making mode, this function populates the
 *  simulation box with a Poisson sample of particles.
 */
#if (MAKEGLASS > 1)
void seed_glass(void)
{
  int i, k, n_for_this_task;
  double Range[3], LowerBound[3];
  double drandom, partmass;
  long long IDstart;

  All.TotNumPart = MAKEGLASS;
  partmass = All.Omega0 * (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G))
    * (All.BoxSize * All.BoxSize * All.BoxSize) / All.TotNumPart;

  All.MaxPart = All.PartAllocFactor * (All.TotNumPart / NTask);	/* sets the maximum number of particles that may */

  allocate_memory();

  header.npartTotal[1] = All.TotNumPart;
  header.mass[1] = partmass;

  if(ThisTask == 0)
    {
      printf("\nGlass initialising\nPartMass= %g\n", partmass);
      printf("TotNumPart= %d%09d\n\n",
	     (int) (All.TotNumPart / 1000000000), (int) (All.TotNumPart % 1000000000));
    }

  /* set the number of particles assigned locally to this task */
  n_for_this_task = All.TotNumPart / NTask;

  if(ThisTask == NTask - 1)
    n_for_this_task = All.TotNumPart - (NTask - 1) * n_for_this_task;

  NumPart = 0;
  IDstart = 1 + (All.TotNumPart / NTask) * ThisTask;

  /* split the temporal domain into Ntask slabs in z-direction */

  Range[0] = Range[1] = All.BoxSize;
  Range[2] = All.BoxSize / NTask;
  LowerBound[0] = LowerBound[1] = 0;
  LowerBound[2] = ThisTask * Range[2];

  srand48(ThisTask);

  for(i = 0; i < n_for_this_task; i++)
    {
      for(k = 0; k < 3; k++)
	{
	  drandom = drand48();

	  P[i].Pos[k] = LowerBound[k] + Range[k] * drandom;
	  P[i].Vel[k] = 0;
	}

      P[i].Mass = partmass;
      P[i].Type = 1;
      P[i].ID = IDstart + i;

      NumPart++;
    }
}
#endif
