// CUDA implementation of ATON
//
// If you make a change to this file, remember to make the same change to
// cpu_numerical.cc.
//
// Dominique Aubert <dominique.aubert@unistra.fr>
// Timothy Stranex <timothy@stranex.com>
//
// This file uses only 'physical' quantities, never comoving quantities.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gpu.h"
#include "constants.h"

#define BLOCKCOOL 128 // Must be <128
#define GRIDCOOLX 4   // Can be increased (4 for 256^3) but make sure that GRIDCOOLY > 0.
#define GRIDCOOLY ((NCELLX*NCELLY*NCELLZ)/GRIDCOOLX/BLOCKCOOL)

#define ONE 0.99999f
#define ZERO 0.00001f
#define NCELL4 ((NCELLZ+NBOUND2)*(NCELLY+NBOUND2)*(NCELLX+NBOUND2))
#define TINY 1e-26
#define FLAGSOURCE 5.
#define KBOLTZ 1.3806e-23
#define EPSCOOL 0.0001


#ifdef TEST7_RAYTRACING_TRICKS

// Maximum eigenvalue in the y and z directions. This would normally be 1 for
// the GLF scheme but we can reduce it to 0.1 in the case of Test 7 because the
// photon flux is mostly directed in the x direction. The effect is to reduce
// diffusion in the y and z directions so that the clump casts a sharp shadow.
#define LAMBDA_YZ 0.1

#define OTSA

#else

#define LAMBDA_YZ 1.0

#endif


// Bounds on physical quantities
#define MIN_TEMP 1.0e-20
#define MIN_EGY 0
#define MIN_X 0
#define MAX_X 0.99999999999999f
#define MIN_EINT 0



// Add a discrete set of point sources.
//
// source_N: photon number density source [1/m^3/s]. Notice that this is in units of density, not photons per second.
// source_pos: grid positions of the point sources.
// source_count: number of point sources.
//
// dt: time step [s]
//
// rad_N: photon number density [1/m^3]
//
__global__ void cuAddPointSources(
	const double* source_N, const int* source_pos, int source_count,
	double dt,
	double* rad_N)
{
  int source_i = threadIdx.x + blockIdx.x*blockDim.x;
  int i = source_pos[source_i + 0*source_count];
  int j = source_pos[source_i + 1*source_count];
  int k = source_pos[source_i + 2*source_count];
  double s = source_N[source_i];

  int rad_i =
    (i+NBOUND) +
    (j+NBOUND)*(NCELLX+NBOUND2) +
    (k+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  rad_N[rad_i] += s*dt;
}

// Add a photon source field.
//
// photon_source: photon number density source [1/m^3/s]
// dt: time step [s]
// rad_N: photon number density [1/m^3]
//
__global__ void cuAddSourceField(
	const double *photon_source,
	double dt,
	double *rad_N)
{
  int tx = threadIdx.x;
  int bx = blockIdx.x;
  int by = blockIdx.y;
  int idx1 = tx + bx*blockDim.x + by*gridDim.x*blockDim.x;
  int k = idx1/(NCELLX*NCELLY);
  int j = (idx1-k*(NCELLX*NCELLY))/NCELLX;
  int i = idx1-k*(NCELLX*NCELLY)-j*(NCELLX);

  int rad_i =
    (i+NBOUND) +
    (j+NBOUND)*(NCELLX+NBOUND2) +
    (k+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  rad_N[rad_i] += photon_source[rad_i] * dt;
}


// Update the photon source density using the GLF scheme.
//
// rad_N: photon number density [1/m^3]
// rad_F: photon number flux [1/m^2 1/s]
//
// c: speed of light [m/s]
// dx: grid spacing [m]
// dt: time step [s]
//
// rad_N_new: updated photon number density [1/m^3]
//
__global__ void cuStepRadN(
  const double *rad_N, const double *rad_F,
  double c, double dx, double dt,
  double *rad_N_new)
{
  int tx = threadIdx.x + NBOUND;
  int bx = blockIdx.x + NBOUND;
  int by = blockIdx.y + NBOUND;

  double res;
  double um1,up1,fm1,fp1,u0;

  // Divergence along Z

  int baseidu = (tx)+bx*(NCELLX+NBOUND2);
  um1 = rad_N[baseidu+(by-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)];
  u0 = rad_N[baseidu+(by  )*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)];
  up1 = rad_N[baseidu+(by+1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)];
  
  int baseidf = 2*NCELL4+(tx)+bx*(NCELLX+NBOUND2);
  fm1 = rad_F[baseidf+(by-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)];
  fp1 = rad_F[baseidf+(by+1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)];

  res = u0 - 0.5*((fp1-fm1)+LAMBDA_YZ*c*(2*u0-um1-up1))/dx*dt;

  // Divergence along Y

  baseidu = (tx)+by*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  um1 = rad_N[baseidu+(bx-1)*(NCELLX+NBOUND2)];
  u0 = rad_N[baseidu+(bx  )*(NCELLX+NBOUND2)];
  up1 = rad_N[baseidu+(bx+1)*(NCELLX+NBOUND2)];
  
  baseidf = 1*NCELL4+(tx)+by*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  fm1 = rad_F[baseidf+(bx-1)*(NCELLX+NBOUND2)];
  fp1 = rad_F[baseidf+(bx+1)*(NCELLX+NBOUND2)];

  res = res - 0.5*((fp1-fm1)+LAMBDA_YZ*c*(2*u0-um1-up1))/dx*dt;
  
  // Divergence along X

  baseidu=bx*(NCELLX+NBOUND2)+by*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  baseidf=0*NCELL4+bx*(NCELLX+NBOUND2)+by*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  __shared__ double u[NCELLX+NBOUND2];
  __shared__ double f[NCELLX+NBOUND2];
  u[tx] = rad_N[baseidu+tx];
  f[tx] = rad_F[baseidf+tx];
  __syncthreads();

  if (tx-NBOUND == 0) {
    u[NBOUND-1] = rad_N[baseidu+tx-1];
    f[NBOUND-1] = rad_F[baseidf+tx-1];
  }
  if (tx-NBOUND == blockDim.x-1) {
     u[NCELLX+NBOUND] = rad_N[baseidu+tx+1];
     f[NCELLX+NBOUND] = rad_F[baseidf+tx+1];
  }

  res = res - 0.5*((f[tx+1]-f[tx-1])+c*(2*u[tx]-u[tx+1]-u[tx-1]))/dx*dt;

  rad_N_new[baseidu+tx] = res;
}

// Return the ij-th component of the Eddington tensor.
//
// fx, fy, fz: photon number flux [1/m^3]
// ee: photon number density [1/m^3]
// c: speed of light [m/s]
// i, j: desired tensor component
//
__device__ double Eddington(double fx, double fy, double fz, double ee, double c, int i, int j)
{
  double n[3];
  n[0] = 0.0;
  n[1] = 0.0;
  n[2] = 0.0;

  double ff = 0.0;
  if(ee > 0.0) {
    ff = sqrtf(fx*fx + fy*fy + fz*fz);
    if(ff > 0.0) {
      n[0] = fx/ff;
      n[1] = fy/ff;
      n[2] = fz/ff;
    }
    ff = ff/(c*ee);
  }
  
  double arg = fmaxf(4.0 - 3.0*ff*ff, 0.0);
  double chi = (3.0 + 4.0*ff*ff) / (5.0 + 2.0*sqrtf(arg));

  double c2e = ee*c*c;
  double res = 0.0;
  if (i == j) {
    res = (1.0 - chi) / 2.0 * c2e;
  }
  arg = (3.0 * chi - 1.0) / 2.0 * c2e;
  res += arg*n[i]*n[j];
  return res;
}

// Update the photon number flux using the GLF scheme.
//
// rad_N: photon number density [1/m^3]
// rad_F: photon number flux [1/m^2 1/s]
//
// c: speed of light [m/s]
// dx: grid spacing [m]
// dt: time step [s]
//
// rad_F_new: updated photon number flux [1/m^2 1/s]
//
__global__ void cuStepRadF(
	const double* rad_N, const double* rad_F,
        double c, double dx, double dt,
	double* rad_F_new)
{
  int tx=threadIdx.x+NBOUND;
  int bx=blockIdx.x +NBOUND;
  int by=blockIdx.y +NBOUND;

  double fm1,fp1;

  // REMINDER LF flux : (fl+fr-ur+ul)*0.5f;
  //  f_icjcks_p =cuCompute_FaceLF(f[2+idx*3],f[2+idxp*3],c*e[idx],c*e[idxp]);

  double resfx, resfy, resfz;

  __shared__ double u[(NCELLX+NBOUND2)*3],fp[(NCELLX+NBOUND2)*3],fm[(NCELLX+NBOUND2)*3],ep[(NCELLX+NBOUND2)],em[(NCELLX+NBOUND2)];

  //================================================ Z DIRECTION =============================================
  
  int baseidu=0*NCELL4+(tx)+bx*(NCELLX+NBOUND2);

  u[0*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)]; // FX local cell
  u[1*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)+NCELL4]; // FX local cell
  u[2*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(by)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)+2*NCELL4]; // FX local cell

  ep[tx]=rad_N[baseidu+(by+1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)]; // E Cell+1
  em[tx]=rad_N[baseidu+(by-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)]; // E Cell+1

  fm[0*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(by-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)]; // E Cell+1
  fm[1*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(by-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)+NCELL4]; // E Cell+1
  fm[2*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(by-1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)+NCELL4*2]; // E Cell+1


  fp[0*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(by+1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)]; // E Cell+1
  fp[1*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(by+1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)+NCELL4]; // E Cell+1
  fp[2*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(by+1)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2)+NCELL4*2]; // E Cell+1


  __syncthreads();

  // FX Divergence along Z
  
  fp1=Eddington(fp[0*(NCELLX+NBOUND2)+tx],fp[1*(NCELLX+NBOUND2)+tx],fp[2*(NCELLX+NBOUND2)+tx],ep[tx],c,0,2);
  fm1=Eddington(fm[0*(NCELLX+NBOUND2)+tx],fm[1*(NCELLX+NBOUND2)+tx],fm[2*(NCELLX+NBOUND2)+tx],em[tx],c,0,2);

  resfx=u[tx+0*(NCELLX+NBOUND2)]-0.5*((fp1-fm1)+LAMBDA_YZ*c*(2*u[tx+0*(NCELLX+NBOUND2)]-fm[tx+0*(NCELLX+NBOUND2)]-fp[tx+0*(NCELLX+NBOUND2)]))/dx*dt; 
  

 // FY Divergence along Z


  fp1=Eddington(fp[0*(NCELLX+NBOUND2)+tx],fp[1*(NCELLX+NBOUND2)+tx],fp[2*(NCELLX+NBOUND2)+tx],ep[tx],c,1,2);
  fm1=Eddington(fm[0*(NCELLX+NBOUND2)+tx],fm[1*(NCELLX+NBOUND2)+tx],fm[2*(NCELLX+NBOUND2)+tx],em[tx],c,1,2);
  resfy=u[tx+1*(NCELLX+NBOUND2)]-0.5*((fp1-fm1)+LAMBDA_YZ*c*(2*u[tx+1*(NCELLX+NBOUND2)]-fm[tx+1*(NCELLX+NBOUND2)]-fp[tx+1*(NCELLX+NBOUND2)]))/dx*dt; 
  

  // FZ Divergence along Z

  fp1=Eddington(fp[0*(NCELLX+NBOUND2)+tx],fp[1*(NCELLX+NBOUND2)+tx],fp[2*(NCELLX+NBOUND2)+tx],ep[tx],c,2,2);
  fm1=Eddington(fm[0*(NCELLX+NBOUND2)+tx],fm[1*(NCELLX+NBOUND2)+tx],fm[2*(NCELLX+NBOUND2)+tx],em[tx],c,2,2);
  resfz=u[tx+2*(NCELLX+NBOUND2)]-0.5*((fp1-fm1)+LAMBDA_YZ*c*(2*u[tx+2*(NCELLX+NBOUND2)]-fm[tx+2*(NCELLX+NBOUND2)]-fp[tx+2*(NCELLX+NBOUND2)]))/dx*dt; 
  
  __syncthreads();


  //================================================ Y DIRECTION =============================================
  
  baseidu=0*NCELL4+(tx)+by*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);

  u[0*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(bx)*(NCELLX+NBOUND2)]; // FX local cell
  u[1*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(bx)*(NCELLX+NBOUND2)+NCELL4]; // FX local cell
  u[2*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(bx)*(NCELLX+NBOUND2)+2*NCELL4]; // FX local cell

  ep[tx]=rad_N[baseidu+(bx+1)*(NCELLX+NBOUND2)]; // E Cell+1
  em[tx]=rad_N[baseidu+(bx-1)*(NCELLX+NBOUND2)]; // E Cell+1

  fm[0*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(bx-1)*(NCELLX+NBOUND2)]; // E Cell+1
  fm[1*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(bx-1)*(NCELLX+NBOUND2)+NCELL4]; // E Cell+1
  fm[2*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(bx-1)*(NCELLX+NBOUND2)+NCELL4*2]; // E Cell+1

  fp[0*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(bx+1)*(NCELLX+NBOUND2)]; // E Cell+1
  fp[1*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(bx+1)*(NCELLX+NBOUND2)+NCELL4]; // E Cell+1
  fp[2*(NCELLX+NBOUND2)+tx]=rad_F[baseidu+(bx+1)*(NCELLX+NBOUND2)+NCELL4*2]; // E Cell+1

  __syncthreads();

  // FX Divergence along Y
  
  fp1=Eddington(fp[0*(NCELLX+NBOUND2)+tx],fp[1*(NCELLX+NBOUND2)+tx],fp[2*(NCELLX+NBOUND2)+tx],ep[tx],c,0,1);
  fm1=Eddington(fm[0*(NCELLX+NBOUND2)+tx],fm[1*(NCELLX+NBOUND2)+tx],fm[2*(NCELLX+NBOUND2)+tx],em[tx],c,0,1);
  resfx=resfx-0.5*((fp1-fm1)+LAMBDA_YZ*c*(2*u[tx+0*(NCELLX+NBOUND2)]-fm[tx+0*(NCELLX+NBOUND2)]-fp[tx+0*(NCELLX+NBOUND2)]))/dx*dt; 

 // FY Divergence along Y

  fp1=Eddington(fp[0*(NCELLX+NBOUND2)+tx],fp[1*(NCELLX+NBOUND2)+tx],fp[2*(NCELLX+NBOUND2)+tx],ep[tx],c,1,1);
  fm1=Eddington(fm[0*(NCELLX+NBOUND2)+tx],fm[1*(NCELLX+NBOUND2)+tx],fm[2*(NCELLX+NBOUND2)+tx],em[tx],c,1,1);
  resfy=resfy-0.5*((fp1-fm1)+LAMBDA_YZ*c*(2*u[tx+1*(NCELLX+NBOUND2)]-fm[tx+1*(NCELLX+NBOUND2)]-fp[tx+1*(NCELLX+NBOUND2)]))/dx*dt; 
  
  // FZ Divergence along Y

  fp1=Eddington(fp[0*(NCELLX+NBOUND2)+tx],fp[1*(NCELLX+NBOUND2)+tx],fp[2*(NCELLX+NBOUND2)+tx],ep[tx],c,2,1);
  fm1=Eddington(fm[0*(NCELLX+NBOUND2)+tx],fm[1*(NCELLX+NBOUND2)+tx],fm[2*(NCELLX+NBOUND2)+tx],em[tx],c,2,1);
  resfz=resfz-0.5*((fp1-fm1)+LAMBDA_YZ*c*(2*u[tx+2*(NCELLX+NBOUND2)]-fm[tx+2*(NCELLX+NBOUND2)]-fp[tx+2*(NCELLX+NBOUND2)]))/dx*dt; 
  


  __syncthreads();


  //================================================ X DIRECTION =============================================

  baseidu=0*NCELL4+bx*(NCELLX+NBOUND2)+by*(NCELLX+NBOUND2)*(NCELLY+NBOUND2);
  
  u[0*(NCELLX+NBOUND2)+tx]=rad_F[tx+baseidu]; // FX local cell
  u[1*(NCELLX+NBOUND2)+tx]=rad_F[tx+baseidu+NCELL4]; // FX local cell
  u[2*(NCELLX+NBOUND2)+tx]=rad_F[tx+baseidu+2*NCELL4]; // FX local cell
  ep[tx]=rad_N[tx+baseidu]; // E Cell+1

  if(tx-NBOUND==0) 
    {
      u[NBOUND-1+0*(NCELLX+NBOUND2)]=rad_F[baseidu+tx-1];
      u[NBOUND-1+1*(NCELLX+NBOUND2)]=rad_F[baseidu+tx-1+NCELL4];
      u[NBOUND-1+2*(NCELLX+NBOUND2)]=rad_F[baseidu+tx-1+2*NCELL4];
      ep[NBOUND-1]=rad_N[tx-1+baseidu]; 
    }

  if(tx-NBOUND==blockDim.x-1) 
    {
      u[NCELLX+NBOUND+0*(NCELLX+NBOUND2)]=rad_F[baseidu+tx+1];
      u[NCELLX+NBOUND+1*(NCELLX+NBOUND2)]=rad_F[baseidu+tx+1+NCELL4];
      u[NCELLX+NBOUND+2*(NCELLX+NBOUND2)]=rad_F[baseidu+tx+1+2*NCELL4];
      ep[NCELLX+NBOUND]=rad_N[tx+1+baseidu]; 
    }

  __syncthreads();


  // FX Divergence along X

  fp1=Eddington(u[0*(NCELLX+NBOUND2)+tx+1],u[1*(NCELLX+NBOUND2)+tx+1],u[2*(NCELLX+NBOUND2)+tx+1],ep[tx+1],c,0,0);
  fm1=Eddington(u[0*(NCELLX+NBOUND2)+tx-1],u[1*(NCELLX+NBOUND2)+tx-1],u[2*(NCELLX+NBOUND2)+tx-1],ep[tx-1],c,0,0);
  resfx=resfx-0.5*((fp1-fm1)+c*(2*u[tx+0*(NCELLX+NBOUND2)]-u[tx+1+0*(NCELLX+NBOUND2)]-u[tx-1+0*(NCELLX+NBOUND2)]))/dx*dt;
  

  // FY Divergence along X

  fp1=Eddington(u[0*(NCELLX+NBOUND2)+tx+1],u[1*(NCELLX+NBOUND2)+tx+1],u[2*(NCELLX+NBOUND2)+tx+1],ep[tx+1],c,1,0);
  fm1=Eddington(u[0*(NCELLX+NBOUND2)+tx-1],u[1*(NCELLX+NBOUND2)+tx-1],u[2*(NCELLX+NBOUND2)+tx-1],ep[tx-1],c,1,0);
  resfy=resfy-0.5*((fp1-fm1)+c*(2*u[tx+1*(NCELLX+NBOUND2)]-u[tx+1+1*(NCELLX+NBOUND2)]-u[tx-1+1*(NCELLX+NBOUND2)]))/dx*dt;
  

  // FZ Divergence along X

  fp1=Eddington(u[0*(NCELLX+NBOUND2)+tx+1],u[1*(NCELLX+NBOUND2)+tx+1],u[2*(NCELLX+NBOUND2)+tx+1],ep[tx+1],c,2,0);
  fm1=Eddington(u[0*(NCELLX+NBOUND2)+tx-1],u[1*(NCELLX+NBOUND2)+tx-1],u[2*(NCELLX+NBOUND2)+tx-1],ep[tx-1],c,2,0);
  resfz=resfz-0.5*((fp1-fm1)+c*(2*u[tx+2*(NCELLX+NBOUND2)]-u[tx+1+2*(NCELLX+NBOUND2)]-u[tx-1+2*(NCELLX+NBOUND2)]))/dx*dt;
  

  rad_F_new[baseidu+tx]=resfx;
  rad_F_new[baseidu+tx+NCELL4]=resfy;
  rad_F_new[baseidu+tx+2*NCELL4]=resfz;

}

// Return the case B recombination rate, alpha_B [m^3/s].
// T: temperature [K]
__device__ double cuCalcAlphaB(double T) {
  T = max(T, MIN_TEMP);  // Protect against divide-by-zero errors.
  double lambda = 2e0 * 157807e0 / T;
  double alpha_b = 2.753e-14 * powf(lambda, 1.5) / powf(1e0 + powf(lambda/2.740, 0.407), 2.242);  // [cm^3/s]
  return alpha_b * 1e-6; // [m^3/s]
}

// Return the case A recombination rate, alpha_A [m^3/s].
// T: temperature [K]
__device__ double cuCalcAlphaA(double T) {
  T = max(T, MIN_TEMP);  // Protect against divide-by-zero errors.
  double lambda = 2e0 * 157807e0 / T;
  double alpha_a = 1.269e-13 * powf(lambda, 1.503) / powf(1e0 + powf(lambda/0.522, 0.470), 1.923);  // [cm^3/s]
  return alpha_a * 1e-6;  // [m^3/s]
}

// Return the ionization rate, beta [m^3/s].
// T: temperature [K]
__device__ double cuCalcBeta(double T) {
  T = max(T, MIN_TEMP);  // Protect against divide-by-zero errors.
  double T5 = T/1e5;
  double beta = 5.85e-11 * sqrtf(T) / (1 + sqrtf(T5)) * expf(-(157809e0/T));  // [cm^3/s]
  return beta * 1e-6;  // [m^3/s]
}


// Calculate the net cooling rate and time.
//
// T: temperature [K]
// x: ionized fraction [1]
// nH: hydrogen number density [1/m^3]
// aexp: cosmological scale factor (aexp=1 today).
//
// lambda: net cooling rate [J/m^3/s]
// tcool: cooling time [s]
//
__device__ void cuCalcCoolingRate(
	double T, double x, double nH,
	double aexp,
	double *lambda, double *tcool) {  
  T = max(T, MIN_TEMP);  // Protect against divide-by-zero errors. 
  double nh2 = nH*1e-6;  // [cm^-3]


  // 1. Collisional Ionization cooling [erg/cm^3/s]
  double c1 = expf(-157809.1e0/T)*1.27e-21*sqrtf(T)/(1e0+sqrtf(T/1e5))*x*(1-x)*nh2*nh2;
  
  // 2. Case A Recombination cooling [erg/cm^3/s]
  double c2 = 1.778e-29*T*powf(2e0*157807e0/T,1.965e0)/powf(1e0+powf(2e0*157807e0/T/0.541e0,0.502e0),2.697e0)*x*x*nh2*nh2;
    
  // 3. Disabled: Case B Recombination cooling [erg/cm^3/s]
  // TODO: Should we enable or disable this depending on OTSA?
  // c3 = 3.435e-30*T*powf(2e0*157807e0/T,1.970e0)/powf(1e0+(powf(2e0*157807e0/T/2.250e0,0.376e0)),3.720e0)*x*x*nh2*nh2;
  double c3 = 0.0;

  // 4. Collisional excitation cooling [erg/cm^3/s]
  double c4 = expf(-118348e0/T)*7.5e-19/(1+sqrtf(T/1e5))*x*(1-x)*nh2*nh2;  
  
  // 5. Bremmsstrahlung [erg/cm^3/s]
  double c5 = 1.42e-27*1.5e0*sqrtf(T)*x*x*nh2*nh2;
  
  // 6. Compton Cooling and Heating [erg/cm^3/s]
  double c6 = 1.017e-37*powf(2.727/aexp,4)*(T-2.727/aexp)*nh2*x;


  // Net cooling rate
  *lambda = c1+c2+c3+c4+c5+c6;   // [erg/cm^3/s]
  *lambda = (*lambda)*1e-7*1e6;  // [J/m^3/s]

  // Cooling time
  double maxc = fmaxf(c1, c2);
  maxc = fmaxf(maxc, c3);
  maxc = fmaxf(maxc, c4);
  maxc = fmaxf(maxc, c5);
  maxc = fmaxf(maxc, fabs(c6));
  maxc = maxc * 1e-7;  // [J/cm^3/s]
  *tcool = 1.5 * nh2*(1+x) * 1.3806e-23 * T / maxc;  // [s]
}


// Update the temperature, ionized fraction, photon density and flux due to the cooling and chemistry.
//
// cudensity: hydrogen number density [1/m^3]
// c: speed of light [m/s]
// dt: time step [s]
// fudgecool, aexp, hubblet: See aton_gpu_step_.
//
// cuxion: ionized fraction
// cutemperature: temperature [K]
// rad_N: photon number density [1/m^3]
// rad_F: photon number flux [1/m^2 1/s]
//
__global__ void cuStepCooling(
	const double* cudensity,
        double fudgecool, double c, double dt, double aexp, double hubblet,
	double* cuxion, double* cutemperature, double* rad_N, double* rad_F) {
  double hnu = AVG_EGY * 1.6022e-19;  // Average photon energy [J]
  double hnu0 = 13.6 * 1.6022e-19;    // Hydrogen ionization energy [J]

  // Notice that the following rates will be lower if the speed of light is reduced.
  // This is intended (see the paper).
  double alphai = AVG_CSN * c;
  double alphae = AVG_CSE * c;

  int tx=threadIdx.x;
  int bx=blockIdx.x;
  int by=blockIdx.y;
  int idx1 = tx+bx*blockDim.x+by*gridDim.x*blockDim.x;
  int k = idx1/(NCELLX*NCELLY);
  int j = (idx1-k*(NCELLX*NCELLY))/NCELLX;
  int i = idx1-k*(NCELLX*NCELLY)-j*(NCELLX);
  int idx = (i+NBOUND)+(j+NBOUND)*(NCELLX+NBOUND2)+(k+NBOUND)*(NCELLX+NBOUND2)*(NCELLY+NBOUND2); // following a convention a[k,j,i] where i varies the first
    
  int idloc = tx;
  int idloc3 = 3*idloc;

  // TODO: As far as I can tell, there's no need to use these shared arrays instead of local variables.
  // They should be removed.
  __shared__ double egyloc[BLOCKCOOL];
  __shared__ double floc[3*BLOCKCOOL];
  __shared__ double x0[BLOCKCOOL];
  __shared__ double nH[BLOCKCOOL];
  __shared__ double tloc[BLOCKCOOL];
  
  x0[idloc] = cuxion[idx];
  nH[idloc] = cudensity[idx];
  egyloc[idloc] = rad_N[idx];
  tloc[idloc] = cutemperature[idx]; 
  floc[0+idloc3] = rad_F[0*NCELL4+idx];
  floc[1+idloc3] = rad_F[1*NCELL4+idx];
  floc[2+idloc3] = rad_F[2*NCELL4+idx];
  __syncthreads();

  double currentcool_t = 0.f;
  int nitcool = 0;
  double dtlimit = dt;
  while (currentcool_t < dt) {
    nitcool++;
    double eint = 1.5*nH[idloc]*KBOLTZ*(1.+x0[idloc])*tloc[idloc];

    //== Getting a timestep
    // 1. Don't evolve past the full time step.
    double dtcool = dt - currentcool_t;
    // 2. Maybe subcycle based on the cooling time.
    double Cool, tcool1;
    cuCalcCoolingRate(tloc[idloc],x0[idloc],nH[idloc],aexp,&Cool,&tcool1);
    const double Heat = egyloc[idloc]*nH[idloc]*(1-x0[idloc]) *
      (alphae*hnu-alphai*hnu0);
    const double eint_rate = Heat - Cool;
    if (fabsf(eint_rate) > 1.0e-30) {
      double tcool = fabsf(eint/eint_rate);
      dtcool = min(dtcool, tcool);
    }
    dtcool = fmin(dtcool, dtlimit);

    //== Cross sections
#ifndef OTSA
    double alpha = cuCalcAlphaA(tloc[idloc]);
#else
    double alpha = cuCalcAlphaB(tloc[idloc]);
#endif
    double alphab = cuCalcAlphaB(tloc[idloc]);
    double beta = cuCalcBeta(tloc[idloc]);

    //== Update
    double q = (alpha-alphab)*x0[idloc]*x0[idloc]*nH[idloc]*nH[idloc];
    double p = alphai*nH[idloc]*(1-x0[idloc]);
    double r = 1.0;

    //== Update photon density
    double new_egy = (q*dtcool + egyloc[idloc]) / (1 + dtcool*(p + 3*hubblet));

    //== Update ionized fraction
    double new_x = 1 - (alpha*x0[idloc]*x0[idloc]*nH[idloc]*dtcool+1.f -x0[idloc])/(1.f+dtcool*(beta*x0[idloc]*nH[idloc]+alphai*new_egy+3*hubblet));

    //== Update internal energy
    // The heating and cooling rates are computed using the new density values.
    cuCalcCoolingRate(tloc[idloc],new_x,nH[idloc],aexp,&Cool,&tcool1);
    eint = (eint + dtcool*(new_egy*nH[idloc]*(1.f-new_x)*(alphae*hnu-alphai*hnu0)-Cool)) / (1 + 3*hubblet*dtcool);

    if (new_egy < MIN_EGY) {
      dtlimit = fudgecool * dtcool;
      continue;
    }
    new_egy = fmaxf(new_egy, MIN_EGY);
    if (new_x < MIN_X || new_x > MAX_X) {
      dtlimit = fudgecool * dtcool;
      continue;
    }
    new_x = fminf(fmaxf(new_x, MIN_X), MAX_X);
    if (eint < MIN_EINT) {
      dtlimit = fudgecool * dtcool;
      continue;
    }
    eint = fmaxf(eint, MIN_EINT);

    //== Save the new values
    egyloc[idloc] = new_egy;
    floc[0+idloc3] = r*floc[0+idloc3] / (1.f + dtcool*(p + 3*hubblet));
    floc[1+idloc3] = r*floc[1+idloc3] / (1.f + dtcool*(p + 3*hubblet));
    floc[2+idloc3] = r*floc[2+idloc3] / (1.f + dtcool*(p + 3*hubblet));
    x0[idloc] = new_x;
    tloc[idloc] = eint / (1.5f*nH[idloc]*KBOLTZ*(1 + x0[idloc]));

    currentcool_t += dtcool;
  }

  cutemperature[idx]=tloc[idloc];
  cuxion[idx]=x0[idloc];
  rad_N[idx]=egyloc[idloc];
  rad_F[0*NCELL4+idx]=floc[0+idloc3];
  rad_F[1*NCELL4+idx]=floc[1+idloc3];
  rad_F[2*NCELL4+idx]=floc[2+idloc3];
  __syncthreads();
}

// Run the transport part of the ATON step. The photon density and flux is updated.
// The function operates on the ATON global arrays.
//
// c_light: speed of light [m/s]
// dx: grid spacing [m]
// dt: time step [s]
//
void gpu_rad_transport(double c_light, double dx, double dt) {
  dim3 blocksimple(NCELLX);
  dim3 gridsimple(NCELLY,NCELLZ);

  cuStepRadN<<<gridsimple,blocksimple>>>(
	cuegy, cuflx, c_light, dx, dt, cuegy_new);
  cudaThreadSynchronize();
  cuStepRadF<<<gridsimple,blocksimple>>>(
	cuegy, cuflx, c_light, dx, dt, cuflx_new);
  cudaThreadSynchronize();
  
  cudaMemcpy(cuegy, cuegy_new, NCELL4*sizeof(double),
	     cudaMemcpyDeviceToDevice);
  cudaMemcpy(cuflx, cuflx_new, NCELL4*sizeof(double)*3,
	     cudaMemcpyDeviceToDevice);
  cudaThreadSynchronize();
}

// Add sources to the photon density.
// The function operates on the ATON global arrays.
//
// dx: grid spacing [m]
// dt: time step [s]
// source_count: number of point sources
//
void gpu_rad_add_sources(double dx, double dt, int source_count) {
  dim3 bcool(BLOCKCOOL);
  dim3 gcool(GRIDCOOLX,GRIDCOOLY);
  cuAddSourceField<<<gcool,bcool>>>(cu_photon_source, dt, cuegy);
  cudaThreadSynchronize();

  int nthreadsource = min(source_count, 128);
  dim3 gridsource((int)(round((double)(source_count)/double(nthreadsource))));
  dim3 blocksource(nthreadsource);
  cuAddPointSources<<<gridsource,blocksource>>>(cusrc0, cusrc0pos, source_count, dt, cuegy);
  cudaThreadSynchronize();
}

// Run the chemistry and cooling part of the ATON step.
// The temperature, ionized fraction, photon density and flux are updated.
// The function operates on the ATON global arrays.
//
// See cuStepCooling for details about the arguments.
//
void gpu_rad_cooling(double c_light, double dx, double dt,
		 double aexp, double hubblet, double fudgecool) {
  dim3 bcool(BLOCKCOOL);
  dim3 gcool(GRIDCOOLX,GRIDCOOLY);
  cuStepCooling<<<gcool,bcool>>>(
	cudensity,
	fudgecool, c_light, dt, aexp, hubblet,
	cuxion, cutemperature, cuegy, cuflx);
  cudaThreadSynchronize();
}

extern "C" void aton_gpu_step_(
	const double* c, const double* dx, const double* dt,
	const int* source_count, const double* fudgecool,
	const double* aexp, const double* hubblet) {
  if (GRIDCOOLY <= 0) {
    printf("Error: GRIDCOOLY <= 0. Adjust the grid size and defines in gpu_numerical.cu.\n");
    return;
  }

  gpu_rad_transport(*c, *dx, *dt);
  gpu_rad_add_sources(*dx, *dt, *source_count);
  gpu_rad_cooling(*c, *dx, *dt, *aexp, *hubblet, *fudgecool);
}
