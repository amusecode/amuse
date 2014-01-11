#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#ifndef NOMPI
#include <mpi.h>
#endif

/*! \file pm_periodic.c
 *  \brief routines for periodic PM-force computation
 */

#ifdef PMGRID
#ifdef PERIODIC

#ifdef NOTYPEPREFIX_FFTW
#include        <rfftw_mpi.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h>	/* double precision FFTW */
#else
#include     <srfftw_mpi.h>
#endif
#endif


#include "allvars.h"
#include "proto.h"

#define  PMGRID2 (2*(PMGRID/2 + 1))




static rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;

static int slab_to_task[PMGRID];
static int *slabs_per_task;
static int *first_slab_of_task;
static int *meshmin_list, *meshmax_list;

static int slabstart_x, nslab_x, slabstart_y, nslab_y, smallest_slab;

static int fftsize, maxfftsize;

static fftw_real *rhogrid, *forcegrid, *workspace;
static fftw_complex *fft_of_rhogrid;


static FLOAT to_slab_fac;


/*! This routines generates the FFTW-plans to carry out the parallel FFTs
 *  later on. Some auxiliary variables are also initialized.
 */
void pm_init_periodic(void)
{
  int i;
  int slab_to_task_local[PMGRID];

  All.Asmth[0] = ASMTH * All.BoxSize / PMGRID;
  All.Rcut[0] = RCUT * All.Asmth[0];

  /* Set up the FFTW plan files. */

#ifndef NOMPI
  fft_forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, PMGRID, PMGRID, PMGRID,
					     FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fft_inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, PMGRID, PMGRID, PMGRID,
					     FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
#else
  fft_forward_plan = rfftw3d_create_plan(PMGRID, PMGRID, PMGRID,
					     FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fft_inverse_plan = rfftw3d_create_plan(PMGRID, PMGRID, PMGRID,
					     FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
#endif
  /* Workspace out the ranges on each processor. */

  rfftwnd_mpi_local_sizes(fft_forward_plan, &nslab_x, &slabstart_x, &nslab_y, &slabstart_y, &fftsize);

  for(i = 0; i < PMGRID; i++)
    slab_to_task_local[i] = 0;

  for(i = 0; i < nslab_x; i++)
    slab_to_task_local[slabstart_x + i] = ThisTask;

  slabs_per_task = malloc(NTask * sizeof(int));
#ifndef NOMPI
  MPI_Allreduce(slab_to_task_local, slab_to_task, PMGRID, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&nslab_x, &smallest_slab, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  MPI_Allgather(&nslab_x, 1, MPI_INT, slabs_per_task, 1, MPI_INT, MPI_COMM_WORLD);
#else
    slab_to_task = slab_to_task_local;
    smallest_slab = nslab_x;
    slabs_per_task[0] = nslab_x;
#endif
  if(ThisTask == 0)
    {
      for(i = 0; i < NTask; i++)
	printf("Task=%d  FFT-Slabs=%d\n", i, slabs_per_task[i]);
    }

  first_slab_of_task = malloc(NTask * sizeof(int));
  MPI_Allgather(&slabstart_x, 1, MPI_INT, first_slab_of_task, 1, MPI_INT, MPI_COMM_WORLD);

  meshmin_list = malloc(3 * NTask * sizeof(int));
  meshmax_list = malloc(3 * NTask * sizeof(int));


  to_slab_fac = PMGRID / All.BoxSize;

#ifndef NOMPI
  MPI_Allreduce(&fftsize, &maxfftsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else
    maxfftsize = fftsize;
#endif
}


/*! This function allocates the memory neeed to compute the long-range PM
 *  force. Three fields are used, one to hold the density (and its FFT, and
 *  then the real-space potential), one to hold the force field obtained by
 *  finite differencing, and finally a workspace field, which is used both as
 *  workspace for the parallel FFT, and as buffer for the communication
 *  algorithm used in the force computation.
 */
void pm_init_periodic_allocate(int dimprod)
{
  static int first_alloc = 1;
  int dimprodmax;
  double bytes_tot = 0;
  size_t bytes;

  MPI_Allreduce(&dimprod, &dimprodmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  /* allocate the memory to hold the FFT fields */

  if(!(rhogrid = (fftw_real *) malloc(bytes = fftsize * sizeof(fftw_real))))
    {
      printf("failed to allocate memory for `FFT-rhogrid' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;


  if(!(forcegrid = (fftw_real *) malloc(bytes = imax(fftsize, dimprodmax) * sizeof(fftw_real))))
    {
      printf("failed to allocate memory for `FFT-forcegrid' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  if(!(workspace = (fftw_real *) malloc(bytes = imax(maxfftsize, dimprodmax) * sizeof(fftw_real))))
    {
      printf("failed to allocate memory for `FFT-workspace' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  if(first_alloc == 1)
    {
      first_alloc = 0;
      if(ThisTask == 0)
	printf("\nAllocated %g MByte for FFT data.\n\n", bytes_tot / (1024.0 * 1024.0));
    }

  fft_of_rhogrid = (fftw_complex *) & rhogrid[0];
}



/*! This routine frees the space allocated for the parallel FFT algorithm.
 */
void pm_init_periodic_free(void)
{
  /* allocate the memory to hold the FFT fields */
  free(workspace);
  free(forcegrid);
  free(rhogrid);
}



/*! Calculates the long-range periodic force given the particle positions
 *  using the PM method.  The force is Gaussian filtered with Asmth, given in
 *  mesh-cell units. We carry out a CIC charge assignment, and compute the
 *  potenial by Fourier transform methods. The potential is finite differenced
 *  using a 4-point finite differencing formula, and the forces are
 *  interpolated tri-linearly to the particle positions. The CIC kernel is
 *  deconvolved. Note that the particle distribution is not in the slab
 *  decomposition that is used for the FFT. Instead, overlapping patches
 *  between local domains and FFT slabs are communicated as needed.
 */
void pmforce_periodic(void)
{
  double k2, kx, ky, kz, smth;
  double dx, dy, dz;
  double fx, fy, fz, ff;
  double asmth2, fac, acc_dim;
  int i, j, slab, level, sendTask, recvTask;
  int x, y, z, xl, yl, zl, xr, yr, zr, xll, yll, zll, xrr, yrr, zrr, ip, dim;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int meshmin[3], meshmax[3], sendmin, sendmax, recvmin, recvmax;
  int rep, ncont, cont_sendmin[2], cont_sendmax[2], cont_recvmin[2], cont_recvmax[2];
  int dimx, dimy, dimz, recv_dimx, recv_dimy, recv_dimz;
#ifndef NOMPI
  MPI_Status status;
#endif


  if(ThisTask == 0)
    {
      printf("Starting periodic PM calculation.\n");
      fflush(stdout);
    }


  force_treefree();


  asmth2 = (2 * M_PI) * All.Asmth[0] / All.BoxSize;
  asmth2 *= asmth2;

  fac = All.G / (M_PI * All.BoxSize);	/* to get potential */
  fac *= 1 / (2 * All.BoxSize / PMGRID);	/* for finite differencing */

  /* first, establish the extension of the local patch in the PMGRID  */

  for(j = 0; j < 3; j++)
    {
      meshmin[j] = PMGRID;
      meshmax[j] = 0;
    }

  for(i = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  slab = to_slab_fac * P[i].Pos[j];
	  if(slab >= PMGRID)
	    slab = PMGRID - 1;

	  if(slab < meshmin[j])
	    meshmin[j] = slab;

	  if(slab > meshmax[j])
	    meshmax[j] = slab;
	}
    }

#ifndef NOMPI
  MPI_Allgather(meshmin, 3, MPI_INT, meshmin_list, 3, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(meshmax, 3, MPI_INT, meshmax_list, 3, MPI_INT, MPI_COMM_WORLD);
#else
    for( i = 0; i < 3; i++){
        mesmmin_list[0+i] = meshmin[i];
        meshmax_list[0+i] = meshmax[i];
    }
#endif
  dimx = meshmax[0] - meshmin[0] + 2;
  dimy = meshmax[1] - meshmin[1] + 2;
  dimz = meshmax[2] - meshmin[2] + 2;

  pm_init_periodic_allocate((dimx + 4) * (dimy + 4) * (dimz + 4));

  for(i = 0; i < dimx * dimy * dimz; i++)
    workspace[i] = 0;

  for(i = 0; i < NumPart; i++)
    {
      slab_x = to_slab_fac * P[i].Pos[0];
      if(slab_x >= PMGRID)
	slab_x = PMGRID - 1;
      dx = to_slab_fac * P[i].Pos[0] - slab_x;
      slab_x -= meshmin[0];
      slab_xx = slab_x + 1;

      slab_y = to_slab_fac * P[i].Pos[1];
      if(slab_y >= PMGRID)
	slab_y = PMGRID - 1;
      dy = to_slab_fac * P[i].Pos[1] - slab_y;
      slab_y -= meshmin[1];
      slab_yy = slab_y + 1;

      slab_z = to_slab_fac * P[i].Pos[2];
      if(slab_z >= PMGRID)
	slab_z = PMGRID - 1;
      dz = to_slab_fac * P[i].Pos[2] - slab_z;
      slab_z -= meshmin[2];
      slab_zz = slab_z + 1;

      workspace[(slab_x * dimy + slab_y) * dimz + slab_z] += P[i].Mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
      workspace[(slab_x * dimy + slab_yy) * dimz + slab_z] += P[i].Mass * (1.0 - dx) * dy * (1.0 - dz);
      workspace[(slab_x * dimy + slab_y) * dimz + slab_zz] += P[i].Mass * (1.0 - dx) * (1.0 - dy) * dz;
      workspace[(slab_x * dimy + slab_yy) * dimz + slab_zz] += P[i].Mass * (1.0 - dx) * dy * dz;

      workspace[(slab_xx * dimy + slab_y) * dimz + slab_z] += P[i].Mass * (dx) * (1.0 - dy) * (1.0 - dz);
      workspace[(slab_xx * dimy + slab_yy) * dimz + slab_z] += P[i].Mass * (dx) * dy * (1.0 - dz);
      workspace[(slab_xx * dimy + slab_y) * dimz + slab_zz] += P[i].Mass * (dx) * (1.0 - dy) * dz;
      workspace[(slab_xx * dimy + slab_yy) * dimz + slab_zz] += P[i].Mass * (dx) * dy * dz;
    }


  for(i = 0; i < fftsize; i++)	/* clear local density field */
    rhogrid[i] = 0;

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;
      if(recvTask < NTask)
	{
	  /* check how much we have to send */
	  sendmin = 2 * PMGRID;
	  sendmax = -1;
	  for(slab_x = meshmin[0]; slab_x < meshmax[0] + 2; slab_x++)
	    if(slab_to_task[slab_x % PMGRID] == recvTask)
	      {
		if(slab_x < sendmin)
		  sendmin = slab_x;
		if(slab_x > sendmax)
		  sendmax = slab_x;
	      }
	  if(sendmax == -1)
	    sendmin = 0;

	  /* check how much we have to receive */
	  recvmin = 2 * PMGRID;
	  recvmax = -1;
	  for(slab_x = meshmin_list[3 * recvTask]; slab_x < meshmax_list[3 * recvTask] + 2; slab_x++)
	    if(slab_to_task[slab_x % PMGRID] == sendTask)
	      {
		if(slab_x < recvmin)
		  recvmin = slab_x;
		if(slab_x > recvmax)
		  recvmax = slab_x;
	      }
	  if(recvmax == -1)
	    recvmin = 0;


	  if((recvmax - recvmin) >= 0 || (sendmax - sendmin) >= 0)	/* ok, we have a contribution to the slab */
	    {
	      recv_dimx = meshmax_list[3 * recvTask + 0] - meshmin_list[3 * recvTask + 0] + 2;
	      recv_dimy = meshmax_list[3 * recvTask + 1] - meshmin_list[3 * recvTask + 1] + 2;
	      recv_dimz = meshmax_list[3 * recvTask + 2] - meshmin_list[3 * recvTask + 2] + 2;

	      if(level > 0)
		{
#ifndef NOMPI
		  MPI_Sendrecv(workspace + (sendmin - meshmin[0]) * dimy * dimz,
			       (sendmax - sendmin + 1) * dimy * dimz * sizeof(fftw_real), MPI_BYTE, recvTask,
			       TAG_PERIODIC_A, forcegrid,
			       (recvmax - recvmin + 1) * recv_dimy * recv_dimz * sizeof(fftw_real), MPI_BYTE,
			       recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);
#else
		  memcpy(forcegrid, workspace + (sendmin - meshmin[0]) * dimy * dimz,
			 (sendmax - sendmin + 1) * dimy * dimz * sizeof(fftw_real));
#endif
		}
	      else
		{
		  memcpy(forcegrid, workspace + (sendmin - meshmin[0]) * dimy * dimz,
			 (sendmax - sendmin + 1) * dimy * dimz * sizeof(fftw_real));
		}

	      for(slab_x = recvmin; slab_x <= recvmax; slab_x++)
		{
		  slab_xx = (slab_x % PMGRID) - first_slab_of_task[ThisTask];

		  if(slab_xx >= 0 && slab_xx < slabs_per_task[ThisTask])
		    {
		      for(slab_y = meshmin_list[3 * recvTask + 1];
			  slab_y <= meshmax_list[3 * recvTask + 1] + 1; slab_y++)
			{
			  slab_yy = slab_y;
			  if(slab_yy >= PMGRID)
			    slab_yy -= PMGRID;

			  for(slab_z = meshmin_list[3 * recvTask + 2];
			      slab_z <= meshmax_list[3 * recvTask + 2] + 1; slab_z++)
			    {
			      slab_zz = slab_z;
			      if(slab_zz >= PMGRID)
				slab_zz -= PMGRID;

			      rhogrid[PMGRID * PMGRID2 * slab_xx + PMGRID2 * slab_yy + slab_zz] +=
				forcegrid[((slab_x - recvmin) * recv_dimy +
					   (slab_y - meshmin_list[3 * recvTask + 1])) * recv_dimz +
					  (slab_z - meshmin_list[3 * recvTask + 2])];
			    }
			}
		    }
		}
	    }
	}
    }

  /* Do the FFT of the density field */

  rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

  /* multiply with Green's function for the potential */

  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < PMGRID; x++)
      for(z = 0; z < PMGRID / 2 + 1; z++)
	{
	  if(x > PMGRID / 2)
	    kx = x - PMGRID;
	  else
	    kx = x;
	  if(y > PMGRID / 2)
	    ky = y - PMGRID;
	  else
	    ky = y;
	  if(z > PMGRID / 2)
	    kz = z - PMGRID;
	  else
	    kz = z;

	  k2 = kx * kx + ky * ky + kz * kz;

	  if(k2 > 0)
	    {
	      smth = -exp(-k2 * asmth2) / k2;

	      /* do deconvolution */

	      fx = fy = fz = 1;
	      if(kx != 0)
		{
		  fx = (M_PI * kx) / PMGRID;
		  fx = sin(fx) / fx;
		}
	      if(ky != 0)
		{
		  fy = (M_PI * ky) / PMGRID;
		  fy = sin(fy) / fy;
		}
	      if(kz != 0)
		{
		  fz = (M_PI * kz) / PMGRID;
		  fz = sin(fz) / fz;
		}
	      ff = 1 / (fx * fy * fz);
	      smth *= ff * ff * ff * ff;

	      /* end deconvolution */

	      ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;
	      fft_of_rhogrid[ip].re *= smth;
	      fft_of_rhogrid[ip].im *= smth;
	    }
	}

  if(slabstart_y == 0)
    fft_of_rhogrid[0].re = fft_of_rhogrid[0].im = 0.0;

  /* Do the FFT to get the potential */

  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

  /* Now rhogrid holds the potential */
  /* construct the potential for the local patch */


  dimx = meshmax[0] - meshmin[0] + 6;
  dimy = meshmax[1] - meshmin[1] + 6;
  dimz = meshmax[2] - meshmin[2] + 6;

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	{

	  /* check how much we have to send */
	  sendmin = 2 * PMGRID;
	  sendmax = -PMGRID;
	  for(slab_x = meshmin_list[3 * recvTask] - 2; slab_x < meshmax_list[3 * recvTask] + 4; slab_x++)
	    if(slab_to_task[(slab_x + PMGRID) % PMGRID] == sendTask)
	      {
		if(slab_x < sendmin)
		  sendmin = slab_x;
		if(slab_x > sendmax)
		  sendmax = slab_x;
	      }
	  if(sendmax == -PMGRID)
	    sendmin = sendmax + 1;


	  /* check how much we have to receive */
	  recvmin = 2 * PMGRID;
	  recvmax = -PMGRID;
	  for(slab_x = meshmin[0] - 2; slab_x < meshmax[0] + 4; slab_x++)
	    if(slab_to_task[(slab_x + PMGRID) % PMGRID] == recvTask)
	      {
		if(slab_x < recvmin)
		  recvmin = slab_x;
		if(slab_x > recvmax)
		  recvmax = slab_x;
	      }
	  if(recvmax == -PMGRID)
	    recvmin = recvmax + 1;

	  if((recvmax - recvmin) >= 0 || (sendmax - sendmin) >= 0)	/* ok, we have a contribution to the slab */
	    {
	      recv_dimx = meshmax_list[3 * recvTask + 0] - meshmin_list[3 * recvTask + 0] + 6;
	      recv_dimy = meshmax_list[3 * recvTask + 1] - meshmin_list[3 * recvTask + 1] + 6;
	      recv_dimz = meshmax_list[3 * recvTask + 2] - meshmin_list[3 * recvTask + 2] + 6;

	      ncont = 1;
	      cont_sendmin[0] = sendmin;
	      cont_sendmax[0] = sendmax;
	      cont_sendmin[1] = sendmax + 1;
	      cont_sendmax[1] = sendmax;

	      cont_recvmin[0] = recvmin;
	      cont_recvmax[0] = recvmax;
	      cont_recvmin[1] = recvmax + 1;
	      cont_recvmax[1] = recvmax;

	      for(slab_x = sendmin; slab_x <= sendmax; slab_x++)
		{
		  if(slab_to_task[(slab_x + PMGRID) % PMGRID] != ThisTask)
		    {
		      /* non-contiguous */
		      cont_sendmax[0] = slab_x - 1;
		      while(slab_to_task[(slab_x + PMGRID) % PMGRID] != ThisTask)
			slab_x++;
		      cont_sendmin[1] = slab_x;
		      ncont++;
		    }
		}

	      for(slab_x = recvmin; slab_x <= recvmax; slab_x++)
		{
		  if(slab_to_task[(slab_x + PMGRID) % PMGRID] != recvTask)
		    {
		      /* non-contiguous */
		      cont_recvmax[0] = slab_x - 1;
		      while(slab_to_task[(slab_x + PMGRID) % PMGRID] != recvTask)
			slab_x++;
		      cont_recvmin[1] = slab_x;
		      if(ncont == 1)
			ncont++;
		    }
		}


	      for(rep = 0; rep < ncont; rep++)
		{
		  sendmin = cont_sendmin[rep];
		  sendmax = cont_sendmax[rep];
		  recvmin = cont_recvmin[rep];
		  recvmax = cont_recvmax[rep];

		  /* prepare what we want to send */
		  if(sendmax - sendmin >= 0)
		    {
		      for(slab_x = sendmin; slab_x <= sendmax; slab_x++)
			{
			  slab_xx = ((slab_x + PMGRID) % PMGRID) - first_slab_of_task[ThisTask];

			  for(slab_y = meshmin_list[3 * recvTask + 1] - 2;
			      slab_y < meshmax_list[3 * recvTask + 1] + 4; slab_y++)
			    {
			      slab_yy = (slab_y + PMGRID) % PMGRID;

			      for(slab_z = meshmin_list[3 * recvTask + 2] - 2;
				  slab_z < meshmax_list[3 * recvTask + 2] + 4; slab_z++)
				{
				  slab_zz = (slab_z + PMGRID) % PMGRID;

				  forcegrid[((slab_x - sendmin) * recv_dimy +
					     (slab_y - (meshmin_list[3 * recvTask + 1] - 2))) * recv_dimz +
					    slab_z - (meshmin_list[3 * recvTask + 2] - 2)] =
				    rhogrid[PMGRID * PMGRID2 * slab_xx + PMGRID2 * slab_yy + slab_zz];
				}
			    }
			}
		    }

		  if(level > 0)
		    {

#ifndef NOMPI
		      MPI_Sendrecv(forcegrid,
				   (sendmax - sendmin + 1) * recv_dimy * recv_dimz * sizeof(fftw_real),
				   MPI_BYTE, recvTask, TAG_PERIODIC_B,
				   workspace + (recvmin - (meshmin[0] - 2)) * dimy * dimz,
				   (recvmax - recvmin + 1) * dimy * dimz * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_PERIODIC_B, MPI_COMM_WORLD, &status);
#else
		      memcpy(workspace + (recvmin - (meshmin[0] - 2)) * dimy * dimz,
			     forcegrid, (recvmax - recvmin + 1) * dimy * dimz * sizeof(fftw_real));
#endif
		    }
		  else
		    {
		      memcpy(workspace + (recvmin - (meshmin[0] - 2)) * dimy * dimz,
			     forcegrid, (recvmax - recvmin + 1) * dimy * dimz * sizeof(fftw_real));
		    }
		}
	    }
	}
    }


  dimx = meshmax[0] - meshmin[0] + 2;
  dimy = meshmax[1] - meshmin[1] + 2;
  dimz = meshmax[2] - meshmin[2] + 2;

  recv_dimx = meshmax[0] - meshmin[0] + 6;
  recv_dimy = meshmax[1] - meshmin[1] + 6;
  recv_dimz = meshmax[2] - meshmin[2] + 6;


  for(dim = 0; dim < 3; dim++)	/* Calculate each component of the force. */
    {
      /* get the force component by finite differencing the potential */
      /* note: "workspace" now contains the potential for the local patch, plus a suffiently large buffer region */

      for(x = 0; x < meshmax[0] - meshmin[0] + 2; x++)
	for(y = 0; y < meshmax[1] - meshmin[1] + 2; y++)
	  for(z = 0; z < meshmax[2] - meshmin[2] + 2; z++)
	    {
	      xrr = xll = xr = xl = x;
	      yrr = yll = yr = yl = y;
	      zrr = zll = zr = zl = z;

	      switch (dim)
		{
		case 0:
		  xr = x + 1;
		  xrr = x + 2;
		  xl = x - 1;
		  xll = x - 2;
		  break;
		case 1:
		  yr = y + 1;
		  yl = y - 1;
		  yrr = y + 2;
		  yll = y - 2;
		  break;
		case 2:
		  zr = z + 1;
		  zl = z - 1;
		  zrr = z + 2;
		  zll = z - 2;
		  break;
		}

	      forcegrid[(x * dimy + y) * dimz + z]
		=
		fac * ((4.0 / 3) *
		       (workspace[((xl + 2) * recv_dimy + (yl + 2)) * recv_dimz + (zl + 2)]
			- workspace[((xr + 2) * recv_dimy + (yr + 2)) * recv_dimz + (zr + 2)]) -
		       (1.0 / 6) *
		       (workspace[((xll + 2) * recv_dimy + (yll + 2)) * recv_dimz + (zll + 2)] -
			workspace[((xrr + 2) * recv_dimy + (yrr + 2)) * recv_dimz + (zrr + 2)]));
	    }

      /* read out the forces */

      for(i = 0; i < NumPart; i++)
	{
	  slab_x = to_slab_fac * P[i].Pos[0];
	  if(slab_x >= PMGRID)
	    slab_x = PMGRID - 1;
	  dx = to_slab_fac * P[i].Pos[0] - slab_x;
	  slab_x -= meshmin[0];
	  slab_xx = slab_x + 1;

	  slab_y = to_slab_fac * P[i].Pos[1];
	  if(slab_y >= PMGRID)
	    slab_y = PMGRID - 1;
	  dy = to_slab_fac * P[i].Pos[1] - slab_y;
	  slab_y -= meshmin[1];
	  slab_yy = slab_y + 1;

	  slab_z = to_slab_fac * P[i].Pos[2];
	  if(slab_z >= PMGRID)
	    slab_z = PMGRID - 1;
	  dz = to_slab_fac * P[i].Pos[2] - slab_z;
	  slab_z -= meshmin[2];
	  slab_zz = slab_z + 1;

	  acc_dim =
	    forcegrid[(slab_x * dimy + slab_y) * dimz + slab_z] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
	  acc_dim += forcegrid[(slab_x * dimy + slab_yy) * dimz + slab_z] * (1.0 - dx) * dy * (1.0 - dz);
	  acc_dim += forcegrid[(slab_x * dimy + slab_y) * dimz + slab_zz] * (1.0 - dx) * (1.0 - dy) * dz;
	  acc_dim += forcegrid[(slab_x * dimy + slab_yy) * dimz + slab_zz] * (1.0 - dx) * dy * dz;

	  acc_dim += forcegrid[(slab_xx * dimy + slab_y) * dimz + slab_z] * (dx) * (1.0 - dy) * (1.0 - dz);
	  acc_dim += forcegrid[(slab_xx * dimy + slab_yy) * dimz + slab_z] * (dx) * dy * (1.0 - dz);
	  acc_dim += forcegrid[(slab_xx * dimy + slab_y) * dimz + slab_zz] * (dx) * (1.0 - dy) * dz;
	  acc_dim += forcegrid[(slab_xx * dimy + slab_yy) * dimz + slab_zz] * (dx) * dy * dz;

	  P[i].GravPM[dim] = acc_dim;
	}
    }

  pm_init_periodic_free();
  force_treeallocate(All.TreeAllocFactor * All.MaxPart, All.MaxPart);

  All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;

  if(ThisTask == 0)
    {
      printf("done PM.\n");
      fflush(stdout);
    }
}


/*! Calculates the long-range potential using the PM method.  The potential is
 *  Gaussian filtered with Asmth, given in mesh-cell units. We carry out a CIC
 *  charge assignment, and compute the potenial by Fourier transform
 *  methods. The CIC kernel is deconvolved.
 */
void pmpotential_periodic(void)
{
  double k2, kx, ky, kz, smth;
  double dx, dy, dz;
  double fx, fy, fz, ff;
  double asmth2, fac;
  int i, j, slab, level, sendTask, recvTask;
  int x, y, z, ip;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int meshmin[3], meshmax[3], sendmin, sendmax, recvmin, recvmax;
  int rep, ncont, cont_sendmin[2], cont_sendmax[2], cont_recvmin[2], cont_recvmax[2];
  int dimx, dimy, dimz, recv_dimx, recv_dimy, recv_dimz;

#ifndef NOMPI
  MPI_Status status;
#endif

  if(ThisTask == 0)
    {
      printf("Starting periodic PM calculation.\n");
      fflush(stdout);
    }

  asmth2 = (2 * M_PI) * All.Asmth[0] / All.BoxSize;
  asmth2 *= asmth2;

  fac = All.G / (M_PI * All.BoxSize);	/* to get potential */

  force_treefree();

  /* first, establish the extension of the local patch in the PMGRID  */

  for(j = 0; j < 3; j++)
    {
      meshmin[j] = PMGRID;
      meshmax[j] = 0;
    }

  for(i = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  slab = to_slab_fac * P[i].Pos[j];
	  if(slab >= PMGRID)
	    slab = PMGRID - 1;

	  if(slab < meshmin[j])
	    meshmin[j] = slab;

	  if(slab > meshmax[j])
	    meshmax[j] = slab;
	}
    }

#ifndef NOMPI
  MPI_Allgather(meshmin, 3, MPI_INT, meshmin_list, 3, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(meshmax, 3, MPI_INT, meshmax_list, 3, MPI_INT, MPI_COMM_WORLD);
#else
    for( i = 0; i < 3; i++){
        mesmmin_list[0+i] = meshmin[i];
        meshmax_list[0+i] = meshmax[i];
    }
#endif
  dimx = meshmax[0] - meshmin[0] + 2;
  dimy = meshmax[1] - meshmin[1] + 2;
  dimz = meshmax[2] - meshmin[2] + 2;

  pm_init_periodic_allocate((dimx + 4) * (dimy + 4) * (dimz + 4));

  for(i = 0; i < dimx * dimy * dimz; i++)
    workspace[i] = 0;

  for(i = 0; i < NumPart; i++)
    {
      slab_x = to_slab_fac * P[i].Pos[0];
      if(slab_x >= PMGRID)
	slab_x = PMGRID - 1;
      dx = to_slab_fac * P[i].Pos[0] - slab_x;
      slab_x -= meshmin[0];
      slab_xx = slab_x + 1;

      slab_y = to_slab_fac * P[i].Pos[1];
      if(slab_y >= PMGRID)
	slab_y = PMGRID - 1;
      dy = to_slab_fac * P[i].Pos[1] - slab_y;
      slab_y -= meshmin[1];
      slab_yy = slab_y + 1;

      slab_z = to_slab_fac * P[i].Pos[2];
      if(slab_z >= PMGRID)
	slab_z = PMGRID - 1;
      dz = to_slab_fac * P[i].Pos[2] - slab_z;
      slab_z -= meshmin[2];
      slab_zz = slab_z + 1;

      workspace[(slab_x * dimy + slab_y) * dimz + slab_z] += P[i].Mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
      workspace[(slab_x * dimy + slab_yy) * dimz + slab_z] += P[i].Mass * (1.0 - dx) * dy * (1.0 - dz);
      workspace[(slab_x * dimy + slab_y) * dimz + slab_zz] += P[i].Mass * (1.0 - dx) * (1.0 - dy) * dz;
      workspace[(slab_x * dimy + slab_yy) * dimz + slab_zz] += P[i].Mass * (1.0 - dx) * dy * dz;

      workspace[(slab_xx * dimy + slab_y) * dimz + slab_z] += P[i].Mass * (dx) * (1.0 - dy) * (1.0 - dz);
      workspace[(slab_xx * dimy + slab_yy) * dimz + slab_z] += P[i].Mass * (dx) * dy * (1.0 - dz);
      workspace[(slab_xx * dimy + slab_y) * dimz + slab_zz] += P[i].Mass * (dx) * (1.0 - dy) * dz;
      workspace[(slab_xx * dimy + slab_yy) * dimz + slab_zz] += P[i].Mass * (dx) * dy * dz;
    }


  for(i = 0; i < fftsize; i++)	/* clear local density field */
    rhogrid[i] = 0;

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;
      if(recvTask < NTask)
	{
	  /* check how much we have to send */
	  sendmin = 2 * PMGRID;
	  sendmax = -1;
	  for(slab_x = meshmin[0]; slab_x < meshmax[0] + 2; slab_x++)
	    if(slab_to_task[slab_x % PMGRID] == recvTask)
	      {
		if(slab_x < sendmin)
		  sendmin = slab_x;
		if(slab_x > sendmax)
		  sendmax = slab_x;
	      }
	  if(sendmax == -1)
	    sendmin = 0;

	  /* check how much we have to receive */
	  recvmin = 2 * PMGRID;
	  recvmax = -1;
	  for(slab_x = meshmin_list[3 * recvTask]; slab_x < meshmax_list[3 * recvTask] + 2; slab_x++)
	    if(slab_to_task[slab_x % PMGRID] == sendTask)
	      {
		if(slab_x < recvmin)
		  recvmin = slab_x;
		if(slab_x > recvmax)
		  recvmax = slab_x;
	      }
	  if(recvmax == -1)
	    recvmin = 0;


	  if((recvmax - recvmin) >= 0 || (sendmax - sendmin) >= 0)	/* ok, we have a contribution to the slab */
	    {
	      recv_dimx = meshmax_list[3 * recvTask + 0] - meshmin_list[3 * recvTask + 0] + 2;
	      recv_dimy = meshmax_list[3 * recvTask + 1] - meshmin_list[3 * recvTask + 1] + 2;
	      recv_dimz = meshmax_list[3 * recvTask + 2] - meshmin_list[3 * recvTask + 2] + 2;

	      if(level > 0)
		{
#ifndef NOMPI
		  MPI_Sendrecv(workspace + (sendmin - meshmin[0]) * dimy * dimz,
			       (sendmax - sendmin + 1) * dimy * dimz * sizeof(fftw_real), MPI_BYTE, recvTask,
			       TAG_PERIODIC_C, forcegrid,
			       (recvmax - recvmin + 1) * recv_dimy * recv_dimz * sizeof(fftw_real), MPI_BYTE,
			       recvTask, TAG_PERIODIC_C, MPI_COMM_WORLD, &status);
#else
		  memcpy(forcegrid, workspace + (sendmin - meshmin[0]) * dimy * dimz,
			 (sendmax - sendmin + 1) * dimy * dimz * sizeof(fftw_real))
#endif
		}
	      else
		{
		  memcpy(forcegrid, workspace + (sendmin - meshmin[0]) * dimy * dimz,
			 (sendmax - sendmin + 1) * dimy * dimz * sizeof(fftw_real));
		}

	      for(slab_x = recvmin; slab_x <= recvmax; slab_x++)
		{
		  slab_xx = (slab_x % PMGRID) - first_slab_of_task[ThisTask];

		  if(slab_xx >= 0 && slab_xx < slabs_per_task[ThisTask])
		    {
		      for(slab_y = meshmin_list[3 * recvTask + 1];
			  slab_y <= meshmax_list[3 * recvTask + 1] + 1; slab_y++)
			{
			  slab_yy = slab_y;
			  if(slab_yy >= PMGRID)
			    slab_yy -= PMGRID;

			  for(slab_z = meshmin_list[3 * recvTask + 2];
			      slab_z <= meshmax_list[3 * recvTask + 2] + 1; slab_z++)
			    {
			      slab_zz = slab_z;
			      if(slab_zz >= PMGRID)
				slab_zz -= PMGRID;

			      rhogrid[PMGRID * PMGRID2 * slab_xx + PMGRID2 * slab_yy + slab_zz] +=
				forcegrid[((slab_x - recvmin) * recv_dimy +
					   (slab_y - meshmin_list[3 * recvTask + 1])) * recv_dimz +
					  (slab_z - meshmin_list[3 * recvTask + 2])];
			    }
			}
		    }
		}
	    }
	}
    }



  /* Do the FFT of the density field */

  rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

  /* multiply with Green's function for the potential */

  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < PMGRID; x++)
      for(z = 0; z < PMGRID / 2 + 1; z++)
	{
	  if(x > PMGRID / 2)
	    kx = x - PMGRID;
	  else
	    kx = x;
	  if(y > PMGRID / 2)
	    ky = y - PMGRID;
	  else
	    ky = y;
	  if(z > PMGRID / 2)
	    kz = z - PMGRID;
	  else
	    kz = z;

	  k2 = kx * kx + ky * ky + kz * kz;

	  if(k2 > 0)
	    {
	      smth = -exp(-k2 * asmth2) / k2 * fac;
	      /* do deconvolution */
	      fx = fy = fz = 1;
	      if(kx != 0)
		{
		  fx = (M_PI * kx) / PMGRID;
		  fx = sin(fx) / fx;
		}
	      if(ky != 0)
		{
		  fy = (M_PI * ky) / PMGRID;
		  fy = sin(fy) / fy;
		}
	      if(kz != 0)
		{
		  fz = (M_PI * kz) / PMGRID;
		  fz = sin(fz) / fz;
		}
	      ff = 1 / (fx * fy * fz);
	      smth *= ff * ff * ff * ff;
	      /* end deconvolution */

	      ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;
	      fft_of_rhogrid[ip].re *= smth;
	      fft_of_rhogrid[ip].im *= smth;
	    }
	}

  if(slabstart_y == 0)
    fft_of_rhogrid[0].re = fft_of_rhogrid[0].im = 0.0;

  /* Do the FFT to get the potential */

  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

  /* note: "rhogrid" now contains the potential */



  dimx = meshmax[0] - meshmin[0] + 6;
  dimy = meshmax[1] - meshmin[1] + 6;
  dimz = meshmax[2] - meshmin[2] + 6;

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	{

	  /* check how much we have to send */
	  sendmin = 2 * PMGRID;
	  sendmax = -PMGRID;
	  for(slab_x = meshmin_list[3 * recvTask] - 2; slab_x < meshmax_list[3 * recvTask] + 4; slab_x++)
	    if(slab_to_task[(slab_x + PMGRID) % PMGRID] == sendTask)
	      {
		if(slab_x < sendmin)
		  sendmin = slab_x;
		if(slab_x > sendmax)
		  sendmax = slab_x;
	      }
	  if(sendmax == -PMGRID)
	    sendmin = sendmax + 1;


	  /* check how much we have to receive */
	  recvmin = 2 * PMGRID;
	  recvmax = -PMGRID;
	  for(slab_x = meshmin[0] - 2; slab_x < meshmax[0] + 4; slab_x++)
	    if(slab_to_task[(slab_x + PMGRID) % PMGRID] == recvTask)
	      {
		if(slab_x < recvmin)
		  recvmin = slab_x;
		if(slab_x > recvmax)
		  recvmax = slab_x;
	      }
	  if(recvmax == -PMGRID)
	    recvmin = recvmax + 1;

	  if((recvmax - recvmin) >= 0 || (sendmax - sendmin) >= 0)	/* ok, we have a contribution to the slab */
	    {
	      recv_dimx = meshmax_list[3 * recvTask + 0] - meshmin_list[3 * recvTask + 0] + 6;
	      recv_dimy = meshmax_list[3 * recvTask + 1] - meshmin_list[3 * recvTask + 1] + 6;
	      recv_dimz = meshmax_list[3 * recvTask + 2] - meshmin_list[3 * recvTask + 2] + 6;

	      ncont = 1;
	      cont_sendmin[0] = sendmin;
	      cont_sendmax[0] = sendmax;
	      cont_sendmin[1] = sendmax + 1;
	      cont_sendmax[1] = sendmax;

	      cont_recvmin[0] = recvmin;
	      cont_recvmax[0] = recvmax;
	      cont_recvmin[1] = recvmax + 1;
	      cont_recvmax[1] = recvmax;

	      for(slab_x = sendmin; slab_x <= sendmax; slab_x++)
		{
		  if(slab_to_task[(slab_x + PMGRID) % PMGRID] != ThisTask)
		    {
		      /* non-contiguous */
		      cont_sendmax[0] = slab_x - 1;
		      while(slab_to_task[(slab_x + PMGRID) % PMGRID] != ThisTask)
			slab_x++;
		      cont_sendmin[1] = slab_x;
		      ncont++;
		    }
		}

	      for(slab_x = recvmin; slab_x <= recvmax; slab_x++)
		{
		  if(slab_to_task[(slab_x + PMGRID) % PMGRID] != recvTask)
		    {
		      /* non-contiguous */
		      cont_recvmax[0] = slab_x - 1;
		      while(slab_to_task[(slab_x + PMGRID) % PMGRID] != recvTask)
			slab_x++;
		      cont_recvmin[1] = slab_x;
		      if(ncont == 1)
			ncont++;
		    }
		}


	      for(rep = 0; rep < ncont; rep++)
		{
		  sendmin = cont_sendmin[rep];
		  sendmax = cont_sendmax[rep];
		  recvmin = cont_recvmin[rep];
		  recvmax = cont_recvmax[rep];

		  /* prepare what we want to send */
		  if(sendmax - sendmin >= 0)
		    {
		      for(slab_x = sendmin; slab_x <= sendmax; slab_x++)
			{
			  slab_xx = ((slab_x + PMGRID) % PMGRID) - first_slab_of_task[ThisTask];

			  for(slab_y = meshmin_list[3 * recvTask + 1] - 2;
			      slab_y < meshmax_list[3 * recvTask + 1] + 4; slab_y++)
			    {
			      slab_yy = (slab_y + PMGRID) % PMGRID;

			      for(slab_z = meshmin_list[3 * recvTask + 2] - 2;
				  slab_z < meshmax_list[3 * recvTask + 2] + 4; slab_z++)
				{
				  slab_zz = (slab_z + PMGRID) % PMGRID;

				  forcegrid[((slab_x - sendmin) * recv_dimy +
					     (slab_y - (meshmin_list[3 * recvTask + 1] - 2))) * recv_dimz +
					    slab_z - (meshmin_list[3 * recvTask + 2] - 2)] =
				    rhogrid[PMGRID * PMGRID2 * slab_xx + PMGRID2 * slab_yy + slab_zz];
				}
			    }
			}
		    }

		  if(level > 0)
		    {
#ifndef NOMPI
		      MPI_Sendrecv(forcegrid,
				   (sendmax - sendmin + 1) * recv_dimy * recv_dimz * sizeof(fftw_real),
				   MPI_BYTE, recvTask, TAG_PERIODIC_D,
				   workspace + (recvmin - (meshmin[0] - 2)) * dimy * dimz,
				   (recvmax - recvmin + 1) * dimy * dimz * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_PERIODIC_D, MPI_COMM_WORLD, &status);
#else
		      memcpy(workspace + (recvmin - (meshmin[0] - 2)) * dimy * dimz,
			     forcegrid, (recvmax - recvmin + 1) * dimy * dimz * sizeof(fftw_real));

#endif
		    }
		  else
		    {
		      memcpy(workspace + (recvmin - (meshmin[0] - 2)) * dimy * dimz,
			     forcegrid, (recvmax - recvmin + 1) * dimy * dimz * sizeof(fftw_real));
		    }
		}
	    }
	}
    }


  dimx = meshmax[0] - meshmin[0] + 2;
  dimy = meshmax[1] - meshmin[1] + 2;
  dimz = meshmax[2] - meshmin[2] + 2;

  recv_dimx = meshmax[0] - meshmin[0] + 6;
  recv_dimy = meshmax[1] - meshmin[1] + 6;
  recv_dimz = meshmax[2] - meshmin[2] + 6;



  for(x = 0; x < meshmax[0] - meshmin[0] + 2; x++)
    for(y = 0; y < meshmax[1] - meshmin[1] + 2; y++)
      for(z = 0; z < meshmax[2] - meshmin[2] + 2; z++)
	{
	  forcegrid[(x * dimy + y) * dimz + z] =
	    workspace[((x + 2) * recv_dimy + (y + 2)) * recv_dimz + (z + 2)];
	}


  /* read out the potential */

  for(i = 0; i < NumPart; i++)
    {
      slab_x = to_slab_fac * P[i].Pos[0];
      if(slab_x >= PMGRID)
	slab_x = PMGRID - 1;
      dx = to_slab_fac * P[i].Pos[0] - slab_x;
      slab_x -= meshmin[0];
      slab_xx = slab_x + 1;

      slab_y = to_slab_fac * P[i].Pos[1];
      if(slab_y >= PMGRID)
	slab_y = PMGRID - 1;
      dy = to_slab_fac * P[i].Pos[1] - slab_y;
      slab_y -= meshmin[1];
      slab_yy = slab_y + 1;

      slab_z = to_slab_fac * P[i].Pos[2];
      if(slab_z >= PMGRID)
	slab_z = PMGRID - 1;
      dz = to_slab_fac * P[i].Pos[2] - slab_z;
      slab_z -= meshmin[2];
      slab_zz = slab_z + 1;

      P[i].Potential +=
	forcegrid[(slab_x * dimy + slab_y) * dimz + slab_z] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
      P[i].Potential += forcegrid[(slab_x * dimy + slab_yy) * dimz + slab_z] * (1.0 - dx) * dy * (1.0 - dz);
      P[i].Potential += forcegrid[(slab_x * dimy + slab_y) * dimz + slab_zz] * (1.0 - dx) * (1.0 - dy) * dz;
      P[i].Potential += forcegrid[(slab_x * dimy + slab_yy) * dimz + slab_zz] * (1.0 - dx) * dy * dz;

      P[i].Potential += forcegrid[(slab_xx * dimy + slab_y) * dimz + slab_z] * (dx) * (1.0 - dy) * (1.0 - dz);
      P[i].Potential += forcegrid[(slab_xx * dimy + slab_yy) * dimz + slab_z] * (dx) * dy * (1.0 - dz);
      P[i].Potential += forcegrid[(slab_xx * dimy + slab_y) * dimz + slab_zz] * (dx) * (1.0 - dy) * dz;
      P[i].Potential += forcegrid[(slab_xx * dimy + slab_yy) * dimz + slab_zz] * (dx) * dy * dz;
    }

  pm_init_periodic_free();
  force_treeallocate(All.TreeAllocFactor * All.MaxPart, All.MaxPart);

  All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;

  if(ThisTask == 0)
    {
      printf("done PM-Potential.\n");
      fflush(stdout);
    }
}

#endif
#endif
