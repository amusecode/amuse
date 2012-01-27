//#include"sapporo2.h"
//#include "/disks/botlek1/iwasawa/work/code/others/sapporo2/libheaders/6thorder.h"
#include "6thorder.h"
#include<iostream>
#include<cmath>


static const int JMEMSIZE = 131072;
static double POS[JMEMSIZE][3];
static double VEL[JMEMSIZE][3];
static double ACC[JMEMSIZE][3];
static double JRK[JMEMSIZE][3];
static double SNP[JMEMSIZE][3];
static double CRK[JMEMSIZE][3];
static double MASS[JMEMSIZE];
static double TIME[JMEMSIZE];
static int ID[JMEMSIZE];
static double EPS_SQ[JMEMSIZE];


static double POS_PRE[JMEMSIZE][3];
static double VEL_PRE[JMEMSIZE][3];
static double ACC_PRE[JMEMSIZE][3];
static double TIME_PRE;

static double OVER3 = 1.0/3.0;
static double OVER4 = 1.0/4.0;
static double OVER5 = 1.0/5.0;

using namespace std;

/*
Initialize the GPU library

SHOULD BE CALLED FIRST! AND ONLY ONCE!

*/
void initialize(){
  std::cerr<<std::endl;
  std::cerr<<std::endl;
  std::cerr<<"------------------------"<<std::endl;
  std::cerr<<"sapporo2 dummy"<<std::endl;
  std::cerr<<"JMEMSIZE="<<JMEMSIZE<<std::endl;
  std::cerr<<std::endl;
  std::cerr<<std::endl;
}



/*
add = address
pos, vel, acc, jrk, snp, crk,
mass,
time = current particle time
id = unique particle id
eps2 = softening of j-particle
*/ 
void set_j_particle(int add, double pos[3], double vel[3], double acc[3],
		    double jrk[3], double snp[3], double crk[3], double mass, double time, int id, double eps2){

  for(int k=0; k<3; k++){
    POS[add][k] = pos[k];
    VEL[add][k] = vel[k];
    ACC[add][k] = acc[k];
    JRK[add][k] = jrk[k];
    SNP[add][k] = snp[k];
    CRK[add][k] = crk[k];
  }

  MASS[add] = mass;
  TIME[add] = time;
  ID[add] = id;
  EPS_SQ[add] = eps2;
}


/*
Set time of the prediction

time = time to which particles are predicted
nj = amount of particles that are predicted

*/
void predict_all(double time, int nj){
  TIME_PRE = time;
  for(int i=0; i<nj; i++){
    double dt = TIME_PRE - TIME[i];
    for(int k=0; k<3; k++){
      POS_PRE[i][k] = ((((CRK[i][k]*dt*OVER5 + SNP[i][k])*dt*OVER4 + JRK[i][k])*dt*OVER3 + ACC[i][k])*dt*0.5 + VEL[i][k])*dt + POS[i][k];
      VEL_PRE[i][k] =  (((CRK[i][k]*dt*OVER4 + SNP[i][k])*dt*OVER3 + JRK[i][k])*dt*0.5   + ACC[i][k])*dt     + VEL[i][k];
      ACC_PRE[i][k] =   ((CRK[i][k]*dt*OVER3 + SNP[i][k])*dt*0.5   + JRK[i][k])*dt       + ACC[i][k];
    }

  }
}


/*

Do not execute prediction, but only copy the particles
into the predicted buffers.

*/
void no_predict_all(double time, int nj){
  TIME_PRE = time;
  for(int i=0; i<nj; i++){
    for(int k=0; k<3; k++){
      POS_PRE[i][k] = POS[i][k];
      VEL_PRE[i][k] = VEL[i][k];
      ACC_PRE[i][k] = ACC[i][k];
    }
  }
}




/*

Return the predicted values for a particle at an address

addr = address of the particle

pos = buffer to store predicted position
vel = buffer to store predicted velocity
acc = buffer to store predicted acceleration

*/
void pick_up_predictor_2(int addr,  int &id, double &mass, double &eps2,  double pos[3],  double vel[3],  double acc[3]){
  id = ID[addr];
  mass = MASS[addr];
  eps2 = EPS_SQ[addr];
  for(int k=0; k<3; k++){
    pos[k] = POS_PRE[addr][k];
    vel[k] = VEL_PRE[addr][k];
    acc[k] = ACC_PRE[addr][k];
  }
}




/*
Calculate the gravity on the i-particles

//Input
ni = number of particles to be integrated
nj = number of sources
pos, vel, acc, mass, eps2

//Output
acc, jrk, snp, potential (phi)
nnb = nearest neighbour ID
nnb_r2 = distance to the nearest neighbour. (Squared distance + softening)
nnb_r2 =   double r2 = EPS2 + dx*dx + dy*dy + dz*dz;


*/

/*
void calc_force_on_predictors(int ni, int nj, int ids[], double pos[][3], double vel[][3],double acc[][3],
                              double mass[], double eps2[],
                              double accNew[][3], double jrkNew[][3],
                              double snpNew[][3], double phi[], int nnb[], double nnb_r2[]){
  double x, y, z;
  double vx, vy, vz;
  double ax_i, ay_i, az_i;
  double ax_o, ay_o, az_o;
  double jx, jy, jz;
  double sx, sy, sz;
  double r_sq;
  double R, R2;
  double mjR3;
  double A1, A2;
  double A, J;
  double rmin_sq = 9999999999999.9;
  int id_rmin = -1;
  double pot;
  for(int i=0; i<ni; i++){
    pot = ax_o = ay_o = az_o = jx = jy = jz = sx = sy = sz = 0.0;
    for(int j=0; j<nj; j++){
      if( ids[i] == ID[j] ){continue;}
      x = pos[i][0] - POS_PRE[j][0];
      y = pos[i][1] - POS_PRE[j][1];
      z = pos[i][2] - POS_PRE[j][2];

      vx = vel[i][0] - VEL_PRE[j][0];
      vy = vel[i][1] - VEL_PRE[j][1];
      vz = vel[i][2] - VEL_PRE[j][2];

      ax_i = acc[i][0] - ACC_PRE[j][0];
      ay_i = acc[i][1] - ACC_PRE[j][1];
      az_i = acc[i][2] - ACC_PRE[j][2];

      r_sq = eps2[i] + x*x + y*y + z*z;

      R2 = 1.0/r_sq;
      R = sqrt(R2); 
      mjR3 = MASS[j]*R2*R;

      A1 = 3.0*R2*(x*vx + y*vy + z*vz);
      A2 = 3.0*( ((vx*vx + vy*vy + vz*vz)+ (x*ax_i + y*ay_i + z*az_i)) *R2 + A1*A1);

      A = -mjR3*x;
      ax_o += A;
      J = -mjR3*vx - A1*A;
      jx += J;
      sx += -mjR3*ax_i - 2.0*A1*J - A2*A;

      A = -mjR3*y;
      ay_o += A;
      J = -mjR3*vy - A1*A;
      jy += J;
      sy += -mjR3*ay_i - 2.0*A1*J - A2*A;

      A = -mjR3*z;
      az_o += A;
      J = -mjR3*vz - A1*A;
      jz += J;
      sz += -mjR3*az_i - 2.0*A1*J - A2*A;

      pot -=  MASS[j]*R;

      if(r_sq < rmin_sq){
	rmin_sq = r_sq;
	id_rmin = ids[i];
      }
    }
    accNew[i][0] = ax_o;
    accNew[i][1] = ay_o;
    accNew[i][2] = az_o;
    jrkNew[i][0] = jx;
    jrkNew[i][1] = jy;
    jrkNew[i][2] = jz;
    snpNew[i][0] = sx;
    snpNew[i][1] = sy;
    snpNew[i][2] = sz;
    phi[i] = pot;
    nnb[i] = id_rmin;
    nnb_r2[i] = rmin_sq;
  }
}
*/

void calc_force_on_predictors(int ni, int nj, int ids[], double pos[][3], double vel[][3],double acc[][3],
                              double mass[], double eps2[],
                              double accNew[][3], double jrkNew[][3],
                              double snpNew[][3], double crkNew[][3], 
			      double phi[], int nnb[], double nnb_r2[]){
  double rij[3], vij[3], aij[3];
  double r_sq, v_sq, rv, ra;
  double R, R2, R3;
  double mjR3;
  double A1, A2;
  double F0[3], F1[3], F2[3];
  double rmin_sq = 9999999999999.9;
  int id_rmin = -1;
  for(int i=0; i<ni; i++){
    rmin_sq = 9999999999999.9;
    id_rmin = -1;
    for(int j=0; j<nj; j++){
      if( ids[i] == ID[j] ){continue;}
      //std::cerr<<"ID[j]="<<ID[j]<<std::endl;
      //std::cerr<<"eps2[j]="<<eps2[j])<<std::endl;
      r_sq = eps2[i];
      v_sq = 0.0;
      rv = 0.0;
      ra = 0.0;
      for(int k=0; k<3; k++){
	rij[k] = pos[i][k] - POS_PRE[j][k];
	vij[k] = vel[i][k] - VEL_PRE[j][k];
	aij[k] = acc[i][k] - ACC_PRE[j][k];
	r_sq += rij[k]*rij[k];
	v_sq += vij[k]*vij[k];
	rv += rij[k]*vij[k];
	ra += rij[k]*aij[k];
      }
      R2 = 1.0/r_sq;
      R = sqrt(R2); 
      R3 = R2*R;
      mjR3 = MASS[j]*R3;
      A1 = rv*R2;
      A2 = (v_sq +  ra)*R2 + A1*A1;
      for(int k=0; k<3; k++){
	F0[k] = -mjR3*rij[k];
	F1[k] = -mjR3*vij[k] - 3.0*A1*F0[k];
	F2[k] = -mjR3*aij[k] - 6.0*A1*F1[k] - 3.0*A2*F0[k];
	accNew[i][k] += F0[k];
	jrkNew[i][k] += F1[k];
	snpNew[i][k] += F2[k];
      }
      phi[i] -= MASS[j]*R;
      if(r_sq < rmin_sq){
	rmin_sq = r_sq;
	id_rmin = ID[j];
      }
    }
    nnb[i] = id_rmin;
    nnb_r2[i] = rmin_sq;
    /*
    if(id_rmin == -1){
      cout<<"CHECK"<<endl;
      cout<<"nj="<<nj<<endl;
      cout<<"ids[i]="<<ids[i]<<endl;
      cout<<"pos[i][0]="<<pos[i][0]<<endl;
      cout<<"POS_PRE[0][0]="<<POS_PRE[0][0]<<endl;
    }
    */
  }
}


/*
void calc_force_on_predictors(int ni, int nj, int ids[], double pos[][3], double vel[][3],double acc[][3],
                              double mass[], double eps2[],
                              double accNew[][3], double jrkNew[][3],
                              double snpNew[][3], double phi[], int nnb[], double nnb_r2[]){
  double x, y, z;
  double vx, vy, vz;
  double ax, ay, az; 
  double jx, jy, jz; 
  double pot;
  double r_sq, v_sq, rv;
  double R, R2;
  double mjR3;
  double A1;
  double F0x, F0y, F0z; 
  double F1x, F1y, F1z; 
  double rmin_sq = 9999999999999.9;
  int id_rmin = -1;
  for(int i=0; i<ni; i++){
    pot = ax = ay = az = jx = jy = jz = 0.0;
    for(int j=0; j<nj; j++){
      if( ids[i] == ID[j] ){continue;}

      x = pos[i][0] - POS_PRE[j][0];
      y = pos[i][1] - POS_PRE[j][1];
      z = pos[i][2] - POS_PRE[j][2];

      vx = vel[i][0] - VEL_PRE[j][0];
      vy = vel[i][1] - VEL_PRE[j][1];
      vz = vel[i][2] - VEL_PRE[j][2];

      r_sq = eps2[i] + x*x + y*y + z*z;
      rv = x*vx + y*vy + z*vz;

      R2 = 1.0/r_sq;
      R = sqrt(R2); 
      mjR3 = MASS[j]*R2*R;
      A1 = 3.0*rv*R2;

      F0x = -mjR3*x;
      F0y = -mjR3*y;
      F0z = -mjR3*z;

      F1x = -mjR3*vx - A1*F0x;
      F1y = -mjR3*vy - A1*F0y;
      F1z = -mjR3*vz - A1*F0z;

      ax += F0x;
      ay += F0y;
      az += F0z;

      jx += F1x;
      jy += F1y;
      jz += F1z;

      pot -= MASS[j]*R;
      if(r_sq < rmin_sq){
	rmin_sq = r_sq;
	id_rmin = ids[i];
      }
    }
    accNew[i][0] = ax;
    accNew[i][1] = ay;
    accNew[i][2] = az;
    jrkNew[i][0] = jx;
    jrkNew[i][1] = jy;
    jrkNew[i][2] = jz;
    phi[i] = pot;
    nnb[i] = id_rmin;
    nnb_r2[i] = rmin_sq;
  }
}
*/


/*
void calc_force_on_predictors(int ni, int nj, int ids[], double pos[][3], double vel[][3],double acc[][3],
                              double mass[], double eps2[],
                              double accNew[][3], double jrkNew[][3],
                              double snpNew[][3], double phi[], int nnb[], double nnb_r2[]){
  double x, y, z;
  double vx, vy, vz;
  double ax, ay, az;
  double jx, jy, jz;
  double pot;
  double mjR, R2, mjR3;
  double rvR2, r_sq;
  double rmin_sq = 9999999999999.9;
  int id_rmin = -1;
  for(int i=0; i<ni; i++){
    pot = ax = ay = az = jx = jy = jz = 0.0;
    for(int j=0; j<nj; j++){
      if( ids[i] == ID[j] ){continue;}
      x = pos[i][0] - POS_PRE[j][0];
      y = pos[i][1] - POS_PRE[j][1];
      z = pos[i][2] - POS_PRE[j][2];

      vx = vel[i][0] - VEL_PRE[j][0];
      vy = vel[i][1] - VEL_PRE[j][1];
      vz = vel[i][2] - VEL_PRE[j][2];

      r_sq = eps2[i] + x*x + y*y + z*z;
      R2 = 1.0/r_sq;

      rvR2 = vx*x + vy*y + vz*z;
      rvR2 *= R2*3.0;

      mjR = sqrt(R2)*MASS[j]; 
      mjR3 = R2*mjR;

      pot -= mjR;

      x *= mjR3;     ax -= x;
      y *= mjR3;     ay -= y;
      z *= mjR3;     az -= z;

      vx *= mjR3;    jx -= vx;
      vy *= mjR3;    jy -= vy;
      vz *= mjR3;    jz -= vz;

      x *= rvR2;     jx += x;
      y *= rvR2;     jy += y;
      z *= rvR2;     jz += z;

      if(r_sq < rmin_sq){
	rmin_sq = r_sq;
	id_rmin = ids[i];
      }
    }
    accNew[i][0] = ax;
    accNew[i][1] = ay;
    accNew[i][2] = az;
    jrkNew[i][0] = jx;
    jrkNew[i][1] = jy;
    jrkNew[i][2] = jz;
    phi[i] = pot;
    nnb[i] = id_rmin;
    nnb_r2[i] = rmin_sq;
  }
}
*/

/*
void calc_force_on_predictors(int ni, int nj, int ids[], double pos[][3], double vel[][3],double acc[][3],
                              double mass[], double eps2[],
                              double accNew[][3], double jrkNew[][3],
                              double snpNew[][3], double phi[], int nnb[], double nnb_r2[]){
  double x, y, z;
  double vx, vy, vz;
  double ax_i, ay_i, az_i;
  double ax_o, ay_o, az_o;
  double jx, jy, jz;
  double sx, sy, sz;
  double pot;
  double mjR, R2, mjR3;
  double rvR2, r_sq;
  double beta;
  double rmin_sq = 9999999999999.9;
  int id_rmin = -1;
  for(int i=0; i<ni; i++){
    pot = ax_o = ay_o = az_o = jx = jy = jz = sx = sy = sz = 0.0;
    for(int j=0; j<nj; j++){
      if( ids[i] == ID[j] ){continue;}
      x = pos[i][0] - POS_PRE[j][0];
      y = pos[i][1] - POS_PRE[j][1];
      z = pos[i][2] - POS_PRE[j][2];

      vx = vel[i][0] - VEL_PRE[j][0];
      vy = vel[i][1] - VEL_PRE[j][1];
      vz = vel[i][2] - VEL_PRE[j][2];

      ax_i = acc[i][0] - ACC_PRE[j][0];
      ay_i = acc[i][1] - ACC_PRE[j][1];
      az_i = acc[i][2] - ACC_PRE[j][2];

      r_sq = eps2[i] + x*x + y*y + z*z;
      R2 = 1.0/r_sq;

      rvR2 = vx*x + vy*y + vz*z;
      beta = 3.0*((vx*vx + vy*vy + vz*vz + x*ax_i + y*ay_i + z*az_i)*R2 + rvR2*rvR2);
      rvR2 *= R2*3.0;

      mjR = sqrt(R2)*MASS[j]; 
      mjR3 = R2*mjR;

      pot -= mjR;

      x *= mjR3;     ax_o -= x;
      y *= mjR3;     ay_o -= y;
      z *= mjR3;     az_o -= z;

      vx *= mjR3;    jx -= vx;
      vy *= mjR3;    jy -= vy;
      vz *= mjR3;    jz -= vz;

      x *= rvR2;     jx += x;
      y *= rvR2;     jy += y;
      z *= rvR2;     jz += z;

      ax_i *= mjR3;    sx -= ax_i;
      ay_i *= mjR3;    sy -= ay_i;
      az_i *= mjR3;    sz -= az_i;

      vx *= rvR2;


      if(r_sq < rmin_sq){
	rmin_sq = r_sq;
	id_rmin = ids[i];
      }
    }
    accNew[i][0] = ax_o;
    accNew[i][1] = ay_o;
    accNew[i][2] = az_o;
    jrkNew[i][0] = jx;
    jrkNew[i][1] = jy;
    jrkNew[i][2] = jz;
    snpNew[i][0] = sx;
    snpNew[i][1] = sy;
    snpNew[i][2] = sz;
    phi[i] = pot;
    nnb[i] = id_rmin;
    nnb_r2[i] = rmin_sq;
  }
}
*/
