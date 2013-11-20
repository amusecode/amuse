// CPU implementation of ATON
// This is just a copy of the GPU version, no special CPU optimizations.
//
// If you make a change to this file, remember to make the same change to
// gpu_numerical.cc.
//
// Dominique Aubert <dominique.aubert@unistra.fr>
// Timothy Stranex <timothy@stranex.com>
//
// This file uses only 'physical' quantities, never comoving quantities.

#include <algorithm>
#include <cstring>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "aton_cpp.h"
#include "aton_fortran.h"
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

using std::max;
using std::min;

static int CellIndex(int i, int j, int k) {
  return ((i + NBOUND) +
          (j + NBOUND)*(NCELLX + NBOUND2) +
          (k + NBOUND)*(NCELLX + NBOUND2)*(NCELLY + NBOUND2));
}

static int CellIndex4(int i, int j, int k, int component) {
  return CellIndex(i, j, k) + component * NCELL4;
}



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
static void AddPointSources(
	const double* source_N, const int* source_pos, int source_count,
	double dt,
	double* rad_N)
{
  for (int source_i = 0; source_i < source_count; ++source_i) {
    int i = source_pos[source_i + 0*source_count];
    int j = source_pos[source_i + 1*source_count];
    int k = source_pos[source_i + 2*source_count];
    double s = source_N[source_i];
    rad_N[CellIndex(i, j, k)] += s*dt;
  }
}

// Add a photon source field.
//
// photon_source: photon number density source [1/m^3/s]
// dt: time step [s]
// rad_N: photon number density [1/m^3]
//
static void AddSourceField(
	const double photon_source,
	double dt,
	double *rad_N)
{
  *rad_N += photon_source * dt;
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
static void StepRadN(
  int i, int j, int k,
  const double *rad_N, const double *rad_F,
  double c, double dx, double dt,
  double *rad_N_new)
{
  double res;
  double um1,up1,fm1,fp1,u0;

  // Divergence along Z

  um1 = rad_N[CellIndex(i, j, k-1)];
  u0  = rad_N[CellIndex(i, j, k  )];
  up1 = rad_N[CellIndex(i, j, k+1)];
  
  fm1 = rad_F[CellIndex4(i, j, k-1, 2)];
  fp1 = rad_F[CellIndex4(i, j, k+1, 2)];

  res = u0 - 0.5*((fp1-fm1)+LAMBDA_YZ*c*(2*u0-um1-up1))/dx*dt;

  // Divergence along Y

  um1 = rad_N[CellIndex(i, j-1, k)];
  u0  = rad_N[CellIndex(i, j,   k)];
  up1 = rad_N[CellIndex(i, j+1, k)];
  
  fm1 = rad_F[CellIndex4(i, j-1, k, 1)];
  fp1 = rad_F[CellIndex4(i, j+1, k, 1)];

  res = res - 0.5*((fp1-fm1)+LAMBDA_YZ*c*(2*u0-um1-up1))/dx*dt;
  
  // Divergence along X

  um1 = rad_N[CellIndex(i-1, j, k)];
  u0  = rad_N[CellIndex(i,   j, k)];
  up1 = rad_N[CellIndex(i+1, j, k)];

  fm1 = rad_F[CellIndex4(i-1, j, k, 0)];
  fp1 = rad_F[CellIndex4(i+1, j, k, 0)];

  res = res - 0.5*((fp1-fm1)+LAMBDA_YZ*c*(2*u0-um1-up1))/dx*dt;

  // Result

  rad_N_new[CellIndex(i, j, k)] = res;
}

// Return the ij-th component of the Eddington tensor.
//
// fx, fy, fz: photon number flux [1/m^3]
// ee: photon number density [1/m^3]
// c: speed of light [m/s]
// i, j: desired tensor component
//
static double Eddington(double fx, double fy, double fz, double ee, double c, int i, int j)
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
static void StepRadF(
	int i, int j, int k,
	const double* rad_N, const double* rad_F,
        double c, double dx, double dt,
	double* rad_F_new)
{
  double fm1,fp1;

  // REMINDER LF flux : (fl+fr-ur+ul)*0.5f;
  //  f_icjcks_p =cuCompute_FaceLF(f[2+idx*3],f[2+idxp*3],c*e[idx],c*e[idxp]);

  double resfx, resfy, resfz;

  double u[3], fp[3], fm[3], ep, em;

  //================================================ Z DIRECTION =============================================
  
  u[0] = rad_F[CellIndex4(i, j, k, 0)];
  u[1] = rad_F[CellIndex4(i, j, k, 1)];
  u[2] = rad_F[CellIndex4(i, j, k, 2)];

  ep = rad_N[CellIndex(i, j, k+1)];
  em = rad_N[CellIndex(i, j, k-1)];

  fm[0] = rad_F[CellIndex4(i, j, k-1, 0)];
  fm[1] = rad_F[CellIndex4(i, j, k-1, 1)];
  fm[2] = rad_F[CellIndex4(i, j, k-1, 2)];

  fp[0] = rad_F[CellIndex4(i, j, k+1, 0)];
  fp[1] = rad_F[CellIndex4(i, j, k+1, 1)];
  fp[2] = rad_F[CellIndex4(i, j, k+1, 2)];

  // FX Divergence along Z
  fp1 = Eddington(fp[0], fp[1], fp[2], ep, c, 0, 2);
  fm1 = Eddington(fm[0], fm[1], fm[2], em, c, 0, 2);
  resfx = u[0] - 0.5*((fp1-fm1) + LAMBDA_YZ*c*(2*u[0]-fm[0]-fp[0]))/dx*dt;
  
  // FY Divergence along Z
  fp1 = Eddington(fp[0], fp[1], fp[2], ep, c, 1, 2);
  fm1 = Eddington(fm[0], fm[1], fm[2], em, c, 1, 2);
  resfy = u[1] - 0.5*((fp1-fm1) + LAMBDA_YZ*c*(2*u[1] - fm[1] - fp[1]))/dx*dt; 
  
  // FZ Divergence along Z
  fp1 = Eddington(fp[0], fp[1], fp[2], ep, c, 2, 2);
  fm1 = Eddington(fm[0], fm[1], fm[2], em, c, 2, 2);
  resfz=u[2] - 0.5*((fp1-fm1) + LAMBDA_YZ*c*(2*u[2] - fm[2] - fp[2]))/dx*dt; 
  

  //================================================ Y DIRECTION =============================================
  
  u[0] = rad_F[CellIndex4(i, j, k, 0)];
  u[1] = rad_F[CellIndex4(i, j, k, 1)];
  u[2] = rad_F[CellIndex4(i, j, k, 2)];

  ep = rad_N[CellIndex(i, j+1, k)];
  em = rad_N[CellIndex(i, j-1, k)];

  fm[0] = rad_F[CellIndex4(i, j-1, k, 0)];
  fm[1] = rad_F[CellIndex4(i, j-1, k, 1)];
  fm[2] = rad_F[CellIndex4(i, j-1, k, 2)];

  fp[0] = rad_F[CellIndex4(i, j+1, k, 0)];
  fp[1] = rad_F[CellIndex4(i, j+1, k, 1)];
  fp[2] = rad_F[CellIndex4(i, j+1, k, 2)];

  // FX Divergence along Y
  fp1 = Eddington(fp[0], fp[1], fp[2], ep, c, 0, 1);
  fm1 = Eddington(fm[0], fm[1], fm[2], em, c, 0, 1);
  resfx = resfx - 0.5*((fp1-fm1) + LAMBDA_YZ*c*(2*u[0] - fm[0] - fp[0]))/dx*dt; 

  // FY Divergence along Y
  fp1 = Eddington(fp[0], fp[1], fp[2], ep, c, 1, 1);
  fm1 = Eddington(fm[0], fm[1], fm[2], em, c, 1, 1);
  resfy = resfy - 0.5*((fp1-fm1) + LAMBDA_YZ*c*(2*u[1] - fm[1] - fp[1]))/dx*dt; 
  
  // FZ Divergence along Y
  fp1 = Eddington(fp[0], fp[1], fp[2], ep, c, 2, 1);
  fm1 = Eddington(fm[0], fm[1], fm[2], em, c, 2, 1);
  resfz = resfz - 0.5*((fp1-fm1) + LAMBDA_YZ*c*(2*u[2] - fm[2] - fp[2]))/dx*dt; 


  //================================================ X DIRECTION =============================================

  u[0] = rad_F[CellIndex4(i, j, k, 0)];
  u[1] = rad_F[CellIndex4(i, j, k, 1)];
  u[2] = rad_F[CellIndex4(i, j, k, 2)];

  ep = rad_N[CellIndex(i+1, j, k)];
  em = rad_N[CellIndex(i-1, j, k)];

  fm[0] = rad_F[CellIndex4(i-1, j, k, 0)];
  fm[1] = rad_F[CellIndex4(i-1, j, k, 1)];
  fm[2] = rad_F[CellIndex4(i-1, j, k, 2)];

  fp[0] = rad_F[CellIndex4(i+1, j, k, 0)];
  fp[1] = rad_F[CellIndex4(i+1, j, k, 1)];
  fp[2] = rad_F[CellIndex4(i+1, j, k, 2)];

  // FX Divergence along X
  fp1 = Eddington(fp[0], fp[1], fp[2], ep, c, 0, 0);
  fm1 = Eddington(fm[0], fm[1], fm[2], em, c, 0, 0);
  resfx = resfx - 0.5*((fp1-fm1) + c*(2*u[0] - fm[0] - fp[0]))/dx*dt;

  // FY Divergence along X
  fp1 = Eddington(fp[0], fp[1], fp[2], ep, c, 1, 0);
  fm1 = Eddington(fm[0], fm[1], fm[2], em, c, 1, 0);
  resfy = resfy - 0.5*((fp1-fm1) + c*(2*u[1] - fm[1] - fp[1]))/dx*dt;

  // FZ Divergence along X
  fp1 = Eddington(fp[0], fp[1], fp[2], ep, c, 2, 0);
  fm1 = Eddington(fm[0], fm[1], fm[2], em, c, 2, 0);
  resfz = resfz - 0.5*((fp1-fm1) + c*(2*u[2] - fm[2] - fp[2]))/dx*dt;


  //====Result
  rad_F_new[CellIndex4(i, j, k, 0)] = resfx;
  rad_F_new[CellIndex4(i, j, k, 1)] = resfy;
  rad_F_new[CellIndex4(i, j, k, 2)] = resfz;
}

// Return the case B recombination rate, alpha_B [m^3/s].
// T: temperature [K]
static double CalcAlphaB(double T) {
  T = max(T, MIN_TEMP);  // Protect against divide-by-zero errors.
  double lambda = 2e0 * 157807e0 / T;
  double alpha_b = 2.753e-14 * powf(lambda, 1.5) / powf(1e0 + powf(lambda/2.740, 0.407), 2.242);  // [cm^3/s]
  return alpha_b * 1e-6; // [m^3/s]
}

// Return the case A recombination rate, alpha_A [m^3/s].
// T: temperature [K]
static double CalcAlphaA(double T) {
  T = max(T, MIN_TEMP);  // Protect against divide-by-zero errors.
  double lambda = 2e0 * 157807e0 / T;
  double alpha_a = 1.269e-13 * powf(lambda, 1.503) / powf(1e0 + powf(lambda/0.522, 0.470), 1.923);  // [cm^3/s]
  return alpha_a * 1e-6;  // [m^3/s]
}

// Return the ionization rate, beta [m^3/s].
// T: temperature [K]
static double CalcBeta(double T) {
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
static void CalcCoolingRate(
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
static void StepCooling(
	const double cudensity,
        const double fudgecool, const double c, const double dt, const double aexp, const double hubblet,
	double* cuxion, double* cutemperature, double* rad_N, double* rad_F) {
  double hnu = AVG_EGY * 1.6022e-19;  // Average photon energy [J]
  double hnu0 = 13.6 * 1.6022e-19;    // Hydrogen ionization energy [J]

  // Notice that the following rates will be lower if the speed of light is reduced.
  // This is intended (see the paper).
  double alphai = AVG_CSN * c;
  double alphae = AVG_CSE * c;

  // TODO: As far as I can tell, there's no need to use these shared arrays instead of local variables.
  // They should be removed.
  double egyloc;
  double floc[3*BLOCKCOOL];
  double x0;
  double nH;
  double tloc;
  
  x0 = *cuxion;
  nH = cudensity;
  egyloc = *rad_N;
  tloc = *cutemperature;
  floc[0] = rad_F[0*NCELL4];
  floc[1] = rad_F[1*NCELL4];
  floc[2] = rad_F[2*NCELL4];

  double currentcool_t = 0.f;
  int nitcool = 0;
  double dtlimit = dt;
  while (currentcool_t < dt) {
    nitcool++;
    double eint = 1.5*nH*KBOLTZ*(1.+x0)*tloc;

    //== Getting a timestep
    // 1. Don't evolve past the full time step.
    double dtcool = dt - currentcool_t;
    // 2. Maybe subcycle based on the cooling time.
    double Cool, tcool1;
    CalcCoolingRate(tloc,x0,nH,aexp,&Cool,&tcool1);
    const double Heat = egyloc*nH*(1-x0) *
      (alphae*hnu-alphai*hnu0);
    const double eint_rate = Heat - Cool;
    if (fabsf(eint_rate) > 1.0e-30) {
      double tcool = fabsf(eint/eint_rate);
      dtcool = min(dtcool, tcool);
    }
    dtcool = fmin(dtcool, dtlimit);

    //== Cross sections
#ifndef OTSA
    double alpha = CalcAlphaA(tloc);
#else
    double alpha = CalcAlphaB(tloc);
#endif
    double alphab = CalcAlphaB(tloc);
    double beta = CalcBeta(tloc);

    //== Update
    double q = (alpha-alphab)*x0*x0*nH*nH;
    double p = alphai*nH*(1-x0);
    double r = 1.0;

    //== Update photon density
    double new_egy = (q*dtcool + egyloc) / (1 + dtcool*(p + 3*hubblet));

    //== Update ionized fraction
    double new_x = 1 - (alpha*x0*x0*nH*dtcool+1.f -x0)/(1.f+dtcool*(beta*x0*nH+alphai*new_egy+3*hubblet));

    //== Update internal energy
    // The heating and cooling rates are computed using the new density values.
    CalcCoolingRate(tloc,new_x,nH,aexp,&Cool,&tcool1);
    eint = (eint + dtcool*(new_egy*nH*(1.f-new_x)*(alphae*hnu-alphai*hnu0)-Cool)) / (1 + 3*hubblet*dtcool);

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
    egyloc = new_egy;
    floc[0] = r*floc[0] / (1.f + dtcool*(p + 3*hubblet));
    floc[1] = r*floc[1] / (1.f + dtcool*(p + 3*hubblet));
    floc[2] = r*floc[2] / (1.f + dtcool*(p + 3*hubblet));
    x0 = new_x;
    tloc = eint / (1.5f*nH*KBOLTZ*(1 + x0));

    currentcool_t += dtcool;
  }

  *cutemperature = tloc;
  *cuxion = x0;
  *rad_N = egyloc;
  rad_F[0*NCELL4] = floc[0];
  rad_F[1*NCELL4] = floc[1];
  rad_F[2*NCELL4] = floc[2];
}

namespace aton {
  
// Run the transport part of the ATON step. The photon density and flux is updated.
// The function operates on the ATON global arrays.
//
// c_light: speed of light [m/s]
// dx: grid spacing [m]
// dt: time step [s]
//
void cpu_rad_transport(State state, double c_light, double dx, double dt) {
  // FIXME: move these to state and allocate only once
  double* E_new = (double*) malloc(NCELL4*sizeof(double));
  double* F_new = (double*) malloc(NCELL4*sizeof(double)*3);

  for (int k = 0; k < NCELLZ; k++) {
    for (int j = 0; j < NCELLY; j++) {
      for (int i = 0; i < NCELLX; i++) {
        StepRadN(i, j, k, state.E, state.F, c_light, dx, dt, E_new);
        StepRadF(i, j, k, state.E, state.F, c_light, dx, dt, F_new);
      }
    }
  }

  memcpy(state.E, E_new, NCELL4*sizeof(double));
  memcpy(state.F, F_new, NCELL4*sizeof(double)*3);

  free(E_new);
  free(F_new);
}

// Add sources to the photon density.
// The function operates on the ATON global arrays.
//
// dx: grid spacing [m]
// dt: time step [s]
// source_count: number of point sources
//
void cpu_rad_add_sources(State state, double dx, double dt) {
  for (int k = 0; k < NCELLZ; k++) {
    for (int j = 0; j < NCELLY; j++) {
      for (int i = 0; i < NCELLX; i++) {
        int cell = CellIndex(i, j, k);
        AddSourceField(state.photon_source[cell], dt, &state.E[cell]);
      }
    }
  }
  AddPointSources(state.point_source, state.point_source_pos, state.point_source_count, dt, state.E);
}

// Run the chemistry and cooling part of the ATON step.
// The temperature, ionized fraction, photon density and flux are updated.
// The function operates on the ATON global arrays.
//
// See cuStepCooling for details about the arguments.
//
void cpu_rad_cooling(State state, double c_light, double dx, double dt,
		 double aexp, double hubblet, double fudgecool) {
  for (int k = 0; k < NCELLZ; k++) {
    for (int j = 0; j < NCELLY; j++) {
      for (int i = 0; i < NCELLX; i++) {
        int cell = CellIndex(i, j, k);
        StepCooling(
	  state.nH[cell],
	  fudgecool, c_light, dt, aexp, hubblet,
	  &state.xHII[cell], &state.T[cell], &state.E[cell], &state.F[cell]);
      }
    }
  }
}

}  // namespace aton

extern "C" void aton_cpu_rad_(
    int *myid,
    double *c, double *dx, double *dt, int *source_count,
    double *fudgecool, double *aexp, double *hubblet,
    double *cpu_e, double* cpu_d, double* cpu_t, double* cpu_x,
    double *cpu_photon_source, double *cpu_f,
    const int* point_source_pos, const double* point_source) {
  if (GRIDCOOLY <= 0) {
    printf("Error: GRIDCOOLY <= 0. Adjust the grid size and defines in gpu_numerical.cu.\n");
    return;
  }

  int ncellx, ncelly, ncellz, nbnd;
  aton_get_grid_size_(&ncellx, &ncelly, &ncellz, &nbnd);
  int n = (ncellx + 2*nbnd)*(ncelly + 2*nbnd)*(ncellz + 2*nbnd);

  aton::State state;
  state.size = n;
  state.E = cpu_e;
  state.nH = cpu_d;
  state.T = cpu_t;
  state.xHII = cpu_x;
  state.photon_source = cpu_photon_source;
  state.F = cpu_f;

  state.point_source_count = *source_count;
  state.point_source_pos = point_source_pos;
  state.point_source = point_source;

  cpu_rad_transport(state, *c, *dx, *dt);
  cpu_rad_add_sources(state, *dx, *dt);
  cpu_rad_cooling(state, *c, *dx, *dt, *aexp, *hubblet, *fudgecool);
}
