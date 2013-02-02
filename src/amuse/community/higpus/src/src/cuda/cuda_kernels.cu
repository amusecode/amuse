#include <kernel.h>

__global__ void update_local_time(int *next, double *local_time, double GTIME){

	unsigned int gtid = blockIdx.x*blockDim.x + threadIdx.x;
	int who = next[gtid];

	if(who < 0)
		return;

	local_time[who] = GTIME;

}

__global__ void sum_partial(double4 *a, double4 *b, unsigned int nextsize){

	unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

	if(i >= nextsize)
      return;

	extern __shared__ double4 shaccelerations[];
   double4 *shacc = (double4*) shaccelerations;
	double4 myacc;

   myacc = b[i];
   shacc[threadIdx.x] = a[i];

   myacc.x += shacc[threadIdx.x].x;
	myacc.y += shacc[threadIdx.x].y;
   myacc.z += shacc[threadIdx.x].z;

   b[i] = myacc;

}


__global__ void Corrector_gpu(double GTIME, 
										double *local_time, 
										double *step, 
										int *next, 
										unsigned long nextsize, 
										double4 *pos_CH, 
										double4 *vel_CH, 
										double4 *a_tot_D, 
										double4 *a1_tot_D, 
										double4 *a2_tot_D,
										double4 *a_H0, 
										double4 *a3_H, 
										double ETA6, 
										double ETA4, 
										double DTMAX, 
										double DTMIN, 
										unsigned int N){

	unsigned int gtid = blockIdx.x*blockDim.x + threadIdx.x;

		double dt;
      int who = next[gtid];
      int who1 = gtid + nextsize;
      int who2 = who1 + nextsize;

		if(gtid >= nextsize )
			return;

      a_H0[gtid].w = a_H0[gtid].x * a_H0[gtid].x +
                    a_H0[gtid].y * a_H0[gtid].y +
                    a_H0[gtid].z * a_H0[gtid].z ;

      a_H0[who1].w = a_H0[who1].x * a_H0[who1].x +
                       a_H0[who1].y * a_H0[who1].y +
                       a_H0[who1].z * a_H0[who1].z ;

      a_H0[who2].w =  a_H0[who2].x * a_H0[who2].x +
                          a_H0[who2].y * a_H0[who2].y +
                          a_H0[who2].z * a_H0[who2].z ;

      double h = GTIME-local_time[who];
      local_time[who] = GTIME;

		double h1 = 0.5*h;
      double h2 = h1*h1;
      double h3 = 0.75/(h1*h1*h1);
		double h4 = 1.5/(h2*h2);
		double h5 = 7.5/(h2*h2*h1);

		double Amin = a_H0[gtid].x - a_tot_D[who].x;
		double Aplu = a_H0[gtid].x + a_tot_D[who].x;
		double Jmin = h1 * (a_H0[who1].x - a1_tot_D[who].x);
	   double Jplu = h1 * (a_H0[who1].x + a1_tot_D[who].x);
    	double Smin = h1 * h1 * (a_H0[who2].x - a2_tot_D[who].x);
    	double Splu = h1 * h1 * (a_H0[who2].x + a2_tot_D[who].x);

		double over= 1.0/15.0;
		
		pos_CH[who].x = pos_CH[who].x + h1*vel_CH[who].x - 0.4*h2*Amin + over*h2*Jplu;
		vel_CH[who].x = vel_CH[who].x + h1*Aplu          - 0.4*h1*Jmin + over*h1*Splu;
		pos_CH[who].x += h1*vel_CH[who].x;

		a3_H[who].x = h3*(-5.0*Amin + 5.0*Jplu - Smin);
    	double a4halfx = h4*(-Jmin + Splu);
    	double a5halfx = h5*(3.0*Amin - 3.0*Jplu + Smin);
    	a3_H[who].x += h1*a4halfx + 0.5*h2*a5halfx;
    	a4halfx += h1*a5halfx;

		Amin = a_H0[gtid].y - a_tot_D[who].y;
		Aplu = a_H0[gtid].y + a_tot_D[who].y;
		Jmin = h1 * (a_H0[who1].y - a1_tot_D[who].y);
		Jplu = h1 * (a_H0[who1].y + a1_tot_D[who].y);
		Smin = h1 * h1 * (a_H0[who2].y - a2_tot_D[who].y);
		Splu = h1 * h1 * (a_H0[who2].y + a2_tot_D[who].y);

		pos_CH[who].y = pos_CH[who].y + h1*vel_CH[who].y - 0.4*h2*Amin + over*h2*Jplu;
      vel_CH[who].y = vel_CH[who].y + h1*Aplu          - 0.4*h1*Jmin + over*h1*Splu;
      pos_CH[who].y += h1*vel_CH[who].y;

      a3_H[who].y = h3*(-5.0*Amin + 5.0*Jplu - Smin);
      double a4halfy = h4*(-Jmin + Splu);
      double a5halfy = h5*(3.0*Amin - 3.0*Jplu + Smin);
      a3_H[who].y += h1*a4halfy + 0.5*h2*a5halfy;
      a4halfy += h1*a5halfy;

		Amin = a_H0[gtid].z - a_tot_D[who].z;
      Aplu = a_H0[gtid].z + a_tot_D[who].z;
      Jmin = h1 * (a_H0[who1].z - a1_tot_D[who].z);
      Jplu = h1 * (a_H0[who1].z + a1_tot_D[who].z);
      Smin = h1 * h1 * (a_H0[who2].z - a2_tot_D[who].z);
      Splu = h1 * h1 * (a_H0[who2].z + a2_tot_D[who].z);

      pos_CH[who].z = pos_CH[who].z + h1*vel_CH[who].z - 0.4*h2*Amin + over*h2*Jplu;
      vel_CH[who].z = vel_CH[who].z + h1*Aplu          - 0.4*h1*Jmin + over*h1*Splu;
      pos_CH[who].z += h1*vel_CH[who].z;

      a3_H[who].z = h3*(-5.0*Amin + 5.0*Jplu - Smin);
      double a4halfz = h4*(-Jmin + Splu);
      double a5halfz = h5*(3.0*Amin - 3.0*Jplu + Smin);
      a3_H[who].z += h1*a4halfz + 0.5*h2*a5halfz;
      a4halfz += h1*a5halfz;

    a3_H[who].w = sqrt(a3_H[who].x*a3_H[who].x + a3_H[who].y*a3_H[who].y + a3_H[who].z*a3_H[who].z);
    double a4mod = sqrt(a4halfx*a4halfx + a4halfy*a4halfy + a4halfz*a4halfz);
    double a5mod = sqrt(a5halfx*a5halfx + a5halfy*a5halfy + a5halfz*a5halfz);

	 double    dt6 = (sqrt(a_H0[gtid].w*a_H0[who2].w) + a_H0[who1].w) / (a5mod*a3_H[who].w + a4mod*a4mod);
    dt6 = ETA6 * pow(dt6,1.0/6.0);
   
	 double stp = h;
      double overh3 = 1.0/(stp*stp*stp);
      double overh2 = 1.0/(stp*stp);

    double a2dx = overh2 * (-6.0 * (a_tot_D[who].x - a_H0[gtid].x) -
             stp * (4.0 * a_H0[who1].x + 2.0 * a1_tot_D[who].x));
    double a2dy = overh2 * (-6.0 * (a_tot_D[who].y - a_H0[gtid].y) -
             stp * (4.0 * a_H0[who1].y + 2.0 * a1_tot_D[who].y));
    double a2dz = overh2 * (-6.0 * (a_tot_D[who].z - a_H0[gtid].z) -
             stp * (4.0 * a_H0[who1].z + 2.0 * a1_tot_D[who].z));

	 double a3dx = overh3 * (12.0 * (a_tot_D[who].x - a_H0[gtid].x) +
             6.0 * stp * (a_H0[who1].x + a1_tot_D[who].x));
    double a3dy = overh3 * (12.0 * (a_tot_D[who].y - a_H0[gtid].y) +
             6.0 * stp * (a_H0[who1].y + a1_tot_D[who].y));
    double a3dz = overh3 * (12.0 * (a_tot_D[who].z - a_H0[gtid].z) +
             6.0 * stp * (a_H0[who1].z + a1_tot_D[who].z));

    a2dx += h*a3dx;
    a2dy += h*a3dy;
    a2dx += h*a3dz;

    a_H0[who2].w =  a2dx*a2dx + a2dy*a2dy + a2dz*a2dz;
    a3_H[who].w = a3dx*a3dx + a3dy*a3dy + a3dz*a3dz;

    double dt4 = sqrt(ETA4*(sqrt(a_H0[gtid].w*a_H0[who2].w) + a_H0[who1].w) / (sqrt(a_H0[who1].w*a3_H[who].w) + a_H0[who2].w));

    dt = 0.5*dt4+0.5*dt6;

    double rest = GTIME / (2.0 * step[who]);
   rest = (double)((int)(rest)) - rest;

//	return;
//	pos_CH[who].x = step[who];
//	return;

   if(dt > 2.0*step[who] && rest == 0.0 && 2.0*step[who] <= DTMAX)
     step[who] *= 2.0;
   else if (dt < 0.5*step[who])
     step[who] *= 0.25;
   else if (dt < step[who])
      step[who]*=0.5;
   
	if(step[who] < DTMIN)
      step[who] = DTMIN;

	a_tot_D[who] = a_H0[gtid];
	a1_tot_D[who] = a_H0[who1];
	a2_tot_D[who] = a_H0[who2];

}


__global__ void Reconstruct(int *nex,
                            unsigned long nextsize,
									 double4 *pc,
									 double4 *vc,
									 double4 *a3,
									 double4 *a,
									 double4 *a1,
									 double4 *a2,
									 double4 *pva3,
									 double4 *aaa) {


unsigned int gtid = blockIdx.x*blockDim.x + threadIdx.x;

int k = gtid/nextsize;
int who = nex[gtid - k*nextsize];

if(gtid<nextsize){
   pc[who] = pva3[gtid];
}
else if(gtid >= nextsize && gtid < 2*nextsize){
   vc[who] = pva3[gtid];
}
else if(gtid >= 2*nextsize && gtid < 3*nextsize){
   a3[who] = pva3[gtid];
}
else if(gtid >= 3*nextsize && gtid < 4*nextsize){
   a[who] = aaa[gtid - 3*nextsize];
}
else if(gtid>= 4*nextsize && gtid < 5*nextsize){
   a1[who] = aaa[gtid - 3*nextsize];
}
else if(gtid>= 5*nextsize && gtid < 6*nextsize){
   a2[who] = aaa[gtid - 3*nextsize];
}


}


__global__ void initvectors(double4 *acc3, float4 *apred){
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	acc3[i].x = acc3[i].y = acc3[i].z = 0.0;
	apred[i].x = apred[i].y = apred[i].z = 0.0f;
}

__global__ void Predictor (const double TIME,
                           double4 *p_pred,
                           float4  *v_pred,
                           float4  *a_pred,
                           double4 *p_corr,
                           double4 *v_corr,
                           double  *loc_time,
                           double4 *acc,
									double4 *acc1,
									double4 *acc2,
                           double4 *acc3,
                				int istart,
                				int* nvec,
                				int ppgpus,
									unsigned int N){



int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
int cost = ppgpus+istart;

if(i>=cost){
   i = nvec[i - cost];
   if(i>=istart && i < cost)
      i=-1;
}
if(i<0)
   return;

 double timestep = TIME - loc_time[i];
  double t2 = timestep * timestep;
  double t3 = t2 * timestep;
  double t4 = t2 * t2;
  double t5 = t4 * timestep;

  t2 *= 0.5;
  t3 *= 0.1666666666666666666666;
  t4 *= 0.0416666666666666666666;
  t5 *= 0.0083333333333333333333;

  double4 myppred;
  myppred.x = p_pred[i].x;
  myppred.y = p_pred[i].y;
  myppred.z = p_pred[i].z;

  float4  mypred;
  mypred.x = v_pred[i].x;
  mypred.y = v_pred[i].y;
  mypred.z = v_pred[i].z;

  double4 mypcorr;
  mypcorr.x = p_corr[i].x;
  mypcorr.y = p_corr[i].y;
  mypcorr.z = p_corr[i].z;

  double4 myvcorr;
  myvcorr.x = v_corr[i].x;
  myvcorr.y = v_corr[i].y;
  myvcorr.z = v_corr[i].z;

  double4 myacc;
  myacc.x = acc[i].x;
  myacc.y = acc[i].y;
  myacc.z = acc[i].z;

  double4 myacc1;
  myacc1.x = acc1[i].x;
  myacc1.y = acc1[i].y;
  myacc1.z = acc1[i].z;

  double4 myacc2;
  myacc2.x = acc2[i].x;
  myacc2.y = acc2[i].y;
  myacc2.z = acc2[i].z;

  double4 myacc3;
  myacc3.x = acc3[i].x;
  myacc3.y = acc3[i].y;
  myacc3.z = acc3[i].z;


  myppred.x = mypcorr.x + timestep * myvcorr.x +
    t2 * myacc.x  +
    t3 * myacc1.x +
    t4 * myacc2.x +
    t5 * myacc3.x ;

  myppred.y = mypcorr.y + timestep * myvcorr.y +
     t2 * myacc.y  +
     t3 * myacc1.y +
     t4 * myacc2.y +
     t5 * myacc3.y ;

  myppred.z = mypcorr.z + timestep * myvcorr.z +
     t2 * myacc.z  +
     t3 * myacc1.z +
     t4 * myacc2.z +
     t5 * myacc3.z ;

  p_pred[i].x = myppred.x;
  p_pred[i].y = myppred.y;
  p_pred[i].z = myppred.z;

  mypred.x = myvcorr.x + timestep * myacc.x +
    t2 * myacc1.x +
    t3 * myacc2.x +
    t4 * myacc3.x ;

  mypred.y = myvcorr.y + timestep * myacc.y +
    t2 * myacc1.y +
    t3 * myacc2.y +
    t4 * myacc3.y ;

  mypred.z = myvcorr.z + timestep * myacc.z +
    t2 * myacc1.z +
    t3 * myacc2.z +
    t4 * myacc3.z ;

  v_pred[i].x = mypred.x;
  v_pred[i].y = mypred.y;
  v_pred[i].z = mypred.z;

  mypred.x = myacc.x + timestep * myacc1.x +
  t2 * myacc2.x +
    t3 * myacc3.x ;

  mypred.y = myacc.y + timestep * myacc1.y +
    t2 * myacc2.y +
    t3 * myacc3.y ;

  mypred.z = myacc.z + timestep * myacc1.z +
    t2 * myacc2.z +
    t3 * myacc3.z ;

  a_pred[i].x = mypred.x;
  a_pred[i].y = mypred.y;
  a_pred[i].z = mypred.z;
}


__global__ void reduce(double4 *ac, double4 *ac1, double4 *ac2, unsigned int bf_real, unsigned int dimension){

   unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int k = dimension*bf_real;
   double4 myacc;

	extern __shared__ double4 shaccelerations[];
   double4 *shacc = (double4*) shaccelerations;

	if(i < k){
		myacc = ac[i];

      shacc[threadIdx.x] = ac[i + k];

      myacc.x += shacc[threadIdx.x].x;
      myacc.y += shacc[threadIdx.x].y;
      myacc.z += shacc[threadIdx.x].z;

      ac[i] = myacc;
   }
	else if (i >= k && i < 2*k){
      myacc = ac1[i - k];

      shacc[threadIdx.x] = ac1[i];

      myacc.x += shacc[threadIdx.x].x;
      myacc.y += shacc[threadIdx.x].y;
      myacc.z += shacc[threadIdx.x].z;

      ac1[i - k] = myacc;
   }
   else {
      myacc = ac2[i - 2*k];

      shacc[threadIdx.x] = ac2[i - k];

      myacc.x += shacc[threadIdx.x].x;
      myacc.y += shacc[threadIdx.x].y;
      myacc.z += shacc[threadIdx.x].z;

      ac2[i - 2*k] = myacc;
   }
}


__global__ void reposition (double4 *ac, double4 *ac1, double4 *ac2, double4 *af, unsigned long nextsize)
{
   unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if(i < nextsize){
      af[i]              = ac[i];
      af[i + nextsize]   = ac1[i];
      af[i + 2*nextsize] = ac2[i];
   }


}

__forceinline__ __device__ void AddPlummerEnergy(double4 dr, double *pot, double *b, double *M, double mul){

#ifdef PLUMMER
      double distance = (dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
      *pot -= dr.w * (*M)/mul * rsqrt(distance + (*b)*(*b));
#endif

}

__forceinline__ __device__ void AddGalacticEnergy(double4 dr, double *pot, double Mscale, double Rscale, double mul){

#ifdef GALAXY
	double b1 = 387.3 / Rscale;
	double M1 = 1.40592e10 / Mscale;

	double a2 = 5317.8 / Rscale;
	double M2 = 8.56080e10 / Mscale;
	double b2 = 250.0 / Rscale;

	double a3 = 12000.0 / Rscale;
	double M3 = 10.70680e10 / Mscale;

	double distance = (dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
	double x = sqrt(distance)/a3;
	double r2d_2 = dr.x * dr.x + dr.y * dr.y;
	double disk = a2 + sqrt(dr.z*dr.z + b2*b2);

   *pot -= dr.w * M1/mul * rsqrt(distance + b1*b1 );
	*pot += dr.w * M3/(1.02*mul*a3) * (log(1.0 + pow(x, 1.02)) - 3.1863227746391254);
	*pot -= dr.w * M2/mul * rsqrt(r2d_2 + disk*disk);

#endif

}


__global__ void energy(double4 *posD, double4 *velD, double *K, double *P, unsigned int N, double EPS, unsigned int istart, unsigned int ppG, double plummer_core, double plummer_mass, double rscale, double mscale){

	unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

	extern __shared__ double4 shPosition[];
	double4 *shPos = (double4*)shPosition;

	double4 myposition = posD[i];
	double4 myvelocity = velD[i];

	double kin = 0.0, pE = 0.0;
	int k = 0, j;
	int totgpus = N/ppG;

      kin = 0.5/totgpus * myposition.w *(myvelocity.x*myvelocity.x +
                               myvelocity.y*myvelocity.y +
                               myvelocity.z*myvelocity.z);

	while(k < ppG){
      j = 0;
      shPos[threadIdx.x] = posD[k + threadIdx.x + istart];
      __syncthreads();

      while(j < blockDim.x){
         double4 dr = {shPos[j].x - myposition.x, shPos[j].y - myposition.y, shPos[j].z - myposition.z, 0.0};
         double distance = (dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
         if(distance != 0.0)
            pE -= 0.5 * myposition.w*shPos[j].w * rsqrt(distance + EPS*EPS);
      j++;
      }
   __syncthreads();
   k += blockDim.x;
   }

	AddPlummerEnergy(myposition, &pE, &plummer_core, &plummer_mass, totgpus);
	AddGalacticEnergy(myposition, &pE, mscale, rscale, totgpus);


	K[i] = kin;
	P[i] = pE;
}

__forceinline__ __device__ void AddGalacticBulge(double4 dr, float4 dv, float4 da, double4 *acc, double4 *jrk, double4 *snp, float M, float b, double mul){

#ifdef GALAXY

  float distance = (dr.x * dr.x + dr.y * dr.y + dr.z * dr.z) + b*b;
  float sqrdist = M * rsqrtf(distance*distance*distance);
  distance = 1.0f/distance;

  float alpha = (dv.x * dr.x + dv.y * dr.y + dv.z * dr.z)*distance;
  float beta = (dv.x*dv.x + dv.y*dv.y + dv.z*dv.z + dr.x*da.x + dr.y*da.y + dr.z*da.z)*distance + alpha*alpha;
  beta *= -3.0f * sqrdist;

  float scalar3 = -3.0f * alpha * sqrdist;

  float3 au = {sqrdist*dv.x + scalar3 * dr.x, sqrdist*dv.y + scalar3 * dr.y, sqrdist*dv.z + scalar3 * dr.z};

  sqrdist *= mul;
  acc->x -= dr.x*sqrdist;
  acc->y -= dr.y*sqrdist;
  acc->z -= dr.z*sqrdist;

  jrk->x -= au.x*mul;
  jrk->y -= au.y*mul;
  jrk->z -= au.z*mul;

  alpha *= -6.0f*mul;
  beta *= mul;

  snp->x -= sqrdist * da.x + beta * dr.x + au.x * alpha;
  snp->y -= sqrdist * da.y + beta * dr.y + au.y * alpha;
  snp->z -= sqrdist * da.z + beta * dr.z + au.z * alpha;

#endif

}

__forceinline__ __device__ void AddGalacticHalo(double4 dr, float4 dv, float4 da, double4 *acc, double4 *jrk, double4 *snp, float M3, float a3, double mul){

#ifdef GALAXY

	float distance = sqrtf(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
	float tmp1, tmp2, tmp3;

   float A2 = 1.0f + powf(distance/a3, 1.02f);
   float cost = - M3 * powf(a3,-2.02f)/A2*powf(distance,-0.98f);

   acc->x += mul*cost*dr.x;
   acc->y += mul*cost*dr.y;
   acc->z += mul*cost*dr.z;

	float A1 = dv.x * dr.x + dv.y * dr.y + dv.z * dr.z;
   float A3 = A1/(distance*distance);
   float alp = 1.02f*(A2-1.0f)*A3/A2;

   alp = - 0.98f*A3 - alp;

   tmp1 = cost*dr.x*(dv.x/dr.x + alp);
   tmp2 = cost*dr.y*(dv.y/dr.y + alp);
   tmp3 = cost*dr.z*(dv.z/dr.z + alp);

   jrk->x += mul*tmp1;
   jrk->y += mul*tmp2;
   jrk->z += mul*tmp3;

	float A1d = (dv.x*dv.x + dv.y*dv.y + dv.z*dv.z + dr.x*da.x + dr.y*da.y + dr.z*da.z);
   M3 = A1d/(distance*distance);
   alp = M3 - 2.0f*A3*A3;
   a3 = 1.02f*(A2-1.0f)*(M3 - 0.98f*A3*A3);
   alp = - 0.98f*alp - a3;

	float derivx = (dr.x*da.x - dv.x*dv.x)/(dr.x*dr.x);
	float derivy = (dr.y*da.y - dv.y*dv.y)/(dr.y*dr.y);
	float derivz = (dr.z*da.z - dv.z*dv.z)/(dr.z*dr.z);

   snp->x += mul*tmp1*tmp1/(cost*dr.x) + mul*cost*dr.x*(derivx + alp);
   snp->y += mul*tmp2*tmp2/(cost*dr.y) + mul*cost*dr.y*(derivy + alp);
   snp->z += mul*tmp3*tmp3/(cost*dr.z) + mul*cost*dr.z*(derivz + alp);

#endif
}

__forceinline__ __device__ void AddGalacticDisk(double4 dr, float4 dv, float4 da, double4 *acc, double4 *jrk, double4 *snp, float M2, float a2, float b2, double mul){

#ifdef GALAXY
   float alp = sqrtf(b2*b2 + dr.z*dr.z);

	float distance = sqrtf(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
   float A3 = distance*distance - dr.z*dr.z + (a2 + alp)*(a2 + alp);

   M2 = - M2/powf(A3,1.5f);

   acc->x += mul*M2*dr.x;
   acc->y += mul*M2*dr.y;
   acc->z += mul*M2*dr.z*(alp+a2)/alp;

	float A1 = dv.x * dr.x + dv.y * dr.y + dv.z * dr.z;
   float cost = A1 + a2/alp*dr.z*dv.z;
   cost = cost/A3;

   float spunto = dr.z*dv.z/alp;

   float tmp1 = M2*dr.x*(dv.x/dr.x - 3.0f * cost);
   float tmp2 = M2*dr.y*(dv.y/dr.y - 3.0f * cost);
   float tmp3 = M2*dr.z*(dv.z/dr.z - 3.0f * cost);

   jrk->x += mul*tmp1;
   jrk->y += mul*tmp2;
   jrk->z += mul*tmp3*(alp+a2)/alp - mul*M2*dr.z*a2/(alp*alp)*spunto;

   float A2 = dv.z*dv.z + dr.z*da.z;

   float vz2 = dv.z*dv.z*dr.z*dr.z;
   float s2punti = (b2*alp*alp - vz2)/(alp*alp*alp);
	float A1d = (dv.x*dv.x + dv.y*dv.y + dv.z*dv.z + dr.x*da.x + dr.y*da.y + dr.z*da.z);

   A2 = A1d + b2/alp*A2 - b2/(alp*alp*alp)*vz2;

   cost = -3.0f*A2/(A3*A3) + 3.0f*cost*cost;

	float derivx = (dr.x*da.x - dv.x*dv.x)/(dr.x*dr.x);
   float derivy = (dr.y*da.y - dv.y*dv.y)/(dr.y*dr.y);
   float derivz = (dr.z*da.z - dv.z*dv.z)/(dr.z*dr.z);

   snp->x += mul*tmp1*tmp1/(M2*dr.x) + mul*M2*dr.x*(derivx + cost);
   snp->y += mul*tmp2*tmp2/(M2*dr.y) + mul*M2*dr.y*(derivy + cost);
   snp->z += mul*tmp3*tmp3/(M2*dr.z) + mul*M2*dr.z*(derivz + cost)*(alp+a2)/alp - mul*2.0f*tmp3*a2/(alp*alp)*spunto +
     mul*M2*dr.z*a2*(s2punti*alp - 2.0f*spunto*spunto)/(alp*alp*alp);

#endif

}

__forceinline__ __device__ void AddAllenSantillan(double4 dr, float4 dv, float4 da, float Mscale, float Rscale, double4 *acc, double4 *jrk, double4 *snp, double mul){

#ifdef GALAXY
	
	AddGalacticBulge(dr, dv, da, acc, jrk, snp, 1.40592e10f/Mscale, 387.3f/Rscale, mul);

	AddGalacticHalo(dr, dv, da, acc, jrk, snp, 1.07068e11f/Mscale, 12000.0f/Rscale, mul);

	AddGalacticDisk(dr, dv, da, acc, jrk, snp, 8.5608e10f/Mscale, 5317.8f/Rscale, 250.0f/Rscale, mul);



#endif

}

__forceinline__ __device__ void AddPlummer(double4 dr, float4 dv, float4 da, double4 *acc, double4 *jrk, double4 *snp, double *b, double *M, double mul){

#ifdef PLUMMER

  float distance = (dr.x * dr.x + dr.y * dr.y + dr.z * dr.z) + (*b)*(*b);
  float sqrdist = *M * rsqrtf(distance*distance*distance);
  distance = 1.0f/distance;

  float alpha = (dv.x * dr.x + dv.y * dr.y + dv.z * dr.z)*distance;
  float beta = (dv.x*dv.x + dv.y*dv.y + dv.z*dv.z + dr.x*da.x + dr.y*da.y + dr.z*da.z)*distance + alpha*alpha;
  beta *= -3.0f * sqrdist;

  float scalar3 = -3.0f * alpha * sqrdist;

    float3 au = {sqrdist*dv.x + scalar3 * dr.x, sqrdist*dv.y + scalar3 * dr.y, sqrdist*dv.z + scalar3 * dr.z};

    sqrdist *= mul;
  acc->x -= dr.x*sqrdist;
  acc->y -= dr.y*sqrdist;
  acc->z -= dr.z*sqrdist;

  jrk->x -= au.x*mul;
  jrk->y -= au.y*mul;
  jrk->z -= au.z*mul;

  alpha *= -6.0f*mul;
  beta *= mul;

  snp->x -= sqrdist * da.x + beta * dr.x + au.x * alpha;
  snp->y -= sqrdist * da.y + beta * dr.y + au.y * alpha;
  snp->z -= sqrdist * da.z + beta * dr.z + au.z * alpha;

#endif

}


__global__  void evaluation(unsigned int N,
									 const double4 *const globalX,
                            const float4 *const globalV,
                            const float4 *const globalA,
                            double4 *aC,
									 double4 *aC1,
									 double4 *aC2,
                            unsigned int istart,
                            int ppGpus,
                            int Bfactor,
                            int dimension,
                            int *next,
                            double *loc_time,
                            double TIME,
									 double EPS,
									 double pl_core,
									 double pl_mass,
									 double rscale,
									 double mscale)
{

    extern __shared__ double4 shPosition[];

    double4 *shPos = (double4*)shPosition;
    float4  *shVel = (float4*)&shPos[blockDim.x];
    float4  *shAcc = (float4*)&shVel[blockDim.x];

    double4 myPosition = {0.0, 0.0, 0.0, 0.0};
    float4 myVelocity = {0.0f, 0.0f, 0.0f, 0.0f};
    float4 myAccelera = {0.0f, 0.0f, 0.0f, 0.0f};

    double4 acc = {0.0, 0.0, 0.0, 0.0};
    double4 jrk = {0.0, 0.0, 0.0, 0.0};
    double4 snp = {0.0, 0.0, 0.0, 0.0};

  int i = 0, j = 0;

  int gtid = blockIdx.x*blockDim.x + threadIdx.x;

int whoami = gtid/dimension;
int particles = ppGpus/Bfactor;
int costante2 = whoami * particles;
double tot_mul = 1.0/(N/ppGpus * Bfactor);

int whop = next[gtid - whoami*dimension];
int address;

if(whop != -1){
  myPosition = globalX[whop];
  myVelocity = globalV[whop];
  myAccelera = globalA[whop];
}

 while(i < particles){
    j = 0;
		 address = i + threadIdx.x + costante2 + istart;
       shPos[threadIdx.x] = globalX[address];
       shVel[threadIdx.x] = globalV[address];
       shAcc[threadIdx.x] = globalA[address];

     __syncthreads();

    while(j < blockDim.x){

float3 dr = {shPos[j].x - myPosition.x, shPos[j].y - myPosition.y, shPos[j].z - myPosition.z};

  float distance = (dr.x * dr.x + dr.y * dr.y + dr.z * dr.z) + EPS*EPS;
  if(distance != 0.0){
     float sqrdist = shPos[j].w * rsqrtf(distance*distance*distance);
     distance = 1.0f/distance;

     float3 dv = {shVel[j].x - myVelocity.x, shVel[j].y - myVelocity.y, shVel[j].z - myVelocity.z};
     float3 da = {shAcc[j].x - myAccelera.x, shAcc[j].y - myAccelera.y, shAcc[j].z - myAccelera.z};

     float alpha = (dv.x * dr.x + dv.y * dr.y + dv.z * dr.z)*distance;
     float beta = (dv.x*dv.x + dv.y*dv.y + dv.z*dv.z + dr.x*da.x + dr.y*da.y + dr.z*da.z)*distance + alpha*alpha;
     beta *= -3.0f * sqrdist;

     float scalar3 = -3.0f * alpha * sqrdist;

     float3 au = {sqrdist*dv.x + scalar3 * dr.x, sqrdist*dv.y + scalar3 * dr.y, sqrdist*dv.z + scalar3 * dr.z};

     acc.x = acc.x + dr.x*sqrdist;
     acc.y = acc.y + dr.y*sqrdist;
     acc.z = acc.z + dr.z*sqrdist;

     jrk.x += au.x;
     jrk.y += au.y;
     jrk.z += au.z;

     alpha *= -6.0f;

     snp.x += sqrdist * da.x + beta * dr.x + au.x * alpha;
     snp.y += sqrdist * da.y + beta * dr.y + au.y * alpha;
     snp.z += sqrdist * da.z + beta * dr.z + au.z * alpha;
  }  
  j++;
}

     __syncthreads();
     i += blockDim.x;
  }

		
	AddPlummer(myPosition, myVelocity, myAccelera, &acc, &jrk, &snp, &pl_core, &pl_mass, tot_mul);
	AddAllenSantillan(myPosition, myVelocity, myAccelera, mscale, rscale, &acc, &jrk, &snp, tot_mul);


   aC[gtid] = acc;
	aC1[gtid] = jrk;
	aC2[gtid] = snp;

}


