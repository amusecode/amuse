#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

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

static void hydro_state_evaluate(FLOAT h, FLOAT pos[3], FLOAT vel[3], 
  FLOAT *dhsml_out, FLOAT *rho_out, FLOAT *rhov_out, FLOAT *rhov2_out, FLOAT
  *rhoe_out);

void hydro_state_at_point(FLOAT pos[3], FLOAT vel[3], FLOAT *h_out, FLOAT *ngb_out, 
  FLOAT *dhsml_out, FLOAT *rho_out, FLOAT *rhov_out, FLOAT *rhov2_out, FLOAT *rhoe_out)
{
  double low,up,low_ngb,up_ngb,h,ngb,dhsml;
  double rho,rhov[3],rhov2,rhoe;
  int iter;

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

  if(h<=All.MinGasHsml) h=All.MinGasHsml;

  low=h;
  up=h;
  low_ngb=All.DesNumNgb*2;
  up_ngb=All.DesNumNgb/2;
  iter=0;
  do
  {
    if(low_ngb>All.DesNumNgb) low/=1.26;
    if(up_ngb<All.DesNumNgb) up*=1.26;
    hydro_state_evaluate(low, pos, vel, &low_ngb, &dhsml,&rho,rhov,&rhov2,&rhoe);
    MPI_AllReduce(MPI_IN_PLACE, &low_ngb, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    hydro_state_evaluate(up, pos, vel, &up_ngb, &dhsml,&rho,rhov,&rhov2,&rhoe);
    MPI_AllReduce(MPI_IN_PLACE, &up_ngb, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if( iter > MAXITER) enrun(3211);
  } while(low_ngb>All.DesNumNgb || up_ngb<All.DesNumNgb)

  iter=0;
  while(up-low > 1.e-3 * up)
  {
    h=(low+up)/2
    if(h <=All.MinGasHsml ) break;
    hydro_state_evaluate(h, pos, vel, &ngb, &dhsml,&rho,rhov,&rhov2,&rhoe)
    MPI_AllReduce(MPI_IN_PLACE, &ngb, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    
    if(ngb>All.DesNumNgb)
      up=h;
    else
      low=h;
    
    if( iter > MAXITER) endrun(3212);
  }

  h=(low+up)/2
  if(h <=All.MinGasHsml ) h=All.MinGasHsml;

  ngb=dhsml=rho=rhov2=rhoe=rhov[0]=rhov[1]=rhov[2]=0;
  hydro_state_evaluate(h, pos, vel, &ngb, &dhsml,&rho,rhov,&rhov2,&rhoe)
  MPI_AllReduce(MPI_IN_PLACE, &ngb, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_AllReduce(MPI_IN_PLACE, &dhsml, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_AllReduce(MPI_IN_PLACE, &rho, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_AllReduce(MPI_IN_PLACE, &rhov, 3, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_AllReduce(MPI_IN_PLACE, &rhov2, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_AllReduce(MPI_IN_PLACE, &rhoe, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  *h_out=h;
  *ngb_out=ngb;
  *dhsml_out=dhsml;
  *rho_out=rho;
  *rhov2_out=rhov2;
  *rhoe_out=rhoe;
  rhov_out[0]=rhov[0];
  rhov_out[1]=rhov[1];
  rhov_out[2]=rhov[2];
}

static void hydro_state_evaluate(FLOAT h, FLOAT pos[3], FLOAT vel[3], 
  FLOAT *numngb, FLOAT *dhsml_out, FLOAT *rho_out, FLOAT *rhov_out, 
  FLOAT *rhov2_out, FLOAT *rhoe_out)
{
  int j, n, startnode, numngb, numngb_inbox;
  double h, h2, fac, hinv, hinv3, hinv4;
  double rho, rhov[3], rhov2, rhoe, wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  double dvx, dvy, dvz;
  double weighted_numngb, dhsmlrho;

  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;

  rho = rhov2 = rhoe = rhov[0] = rhov[1] = rhov[2] = 0;
  weighted_numngb = 0;
  dhsmlrho = 0;

  startnode = All.MaxPart;
  numngb = 0;
  do
    {
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode);

      for(n = 0; n < numngb_inbox; n++)
	{
	  j = Ngblist[n];

	  dx = pos[0] - P[j].Pos[0];
	  dy = pos[1] - P[j].Pos[1];
	  dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
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
	      numngb++;

	      r = sqrt(r2);

	      u = r * hinv;

	      if(u < 0.5)
		{
		  wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		  dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		}
	      else
		{
		  wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		  dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		}

	      mass_j = P[j].Mass;

              fac = mass_j * wk;
	      rho += fac;

	      weighted_numngb += NORM_COEFF * wk / hinv3;

	      dhsmlrho += -mass_j * (NUMDIMS * hinv * wk + u * dwk);

              dvx = vel[0] - SphP[j].VelPred[0];
              dvy = vel[1] - SphP[j].VelPred[1];
	      dvz = vel[2] - SphP[j].VelPred[2];
 
              rhovx[0]-= fac*dvx;
              rhovx[1]-= fac*dvy;
              rhovx[2]-= fac*dvz;

              rhoe+= fac*SphP[j].Entropy;
                  
              rhov2+= fac * (dvx * dvx + dvy * dvy + dvz * dvz);                   
	    }
	}
    }
  while(startnode >= 0);

  *rho_out+=rho;
  *rhoe_out+=rhoe;
  *rhov2_out+=rhov2;
  *dhsml_out+=dhsmlrho;
  *numngb+=weighted_numngb;
  rhov_out[0]+=rhov[0]
  rhov_out[1]+=rhov[1]
  rhov_out[2]+=rhov[2]
}
