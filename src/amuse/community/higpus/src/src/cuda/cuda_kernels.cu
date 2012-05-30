
__global__ void Reconstruct(int *nex,
                            unsigned long nextsize,
                            double4 *A,
                            double4 *B,
									 unsigned int offset) {


unsigned int gtid = blockIdx.x*blockDim.x + threadIdx.x;

int who = nex[gtid];
int p = gtid + offset*nextsize;

if(gtid<nextsize){
  A[who].x = B[p].x;
  A[who].y = B[p].y;
  A[who].z = B[p].z;
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


__global__ void reduce(double4 *ac, unsigned int bf_real, unsigned int dimension){

   unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int k = dimension*bf_real/2;

   double4 myacc = ac[i];

   extern __shared__ double4 shaccelerations[];
   double4 *shacc = (double4*) shaccelerations;

   shacc[threadIdx.x] = ac[i + k];

   myacc.x += shacc[threadIdx.x].x;
   myacc.y += shacc[threadIdx.x].y;
   myacc.z += shacc[threadIdx.x].z;

   ac[i] = myacc;
}


__global__ void reposition (double4 *ac, double4 *af, unsigned int offset, unsigned long nextsize)
{
   unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i < nextsize)
		af[i + offset*nextsize] = ac[i];

}

__forceinline__ __device__ void AddPlummerEnergy(double4 dr, double *pot, double *b, double *M, double mul){

#ifdef PLUMMER
	double distance = (dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
	*pot -= dr.w * (*M)/mul * rsqrt(distance + (*b)*(*b));
#endif

}


__global__ void energy(double4 *posCD, double4 *velCD, float4 *velPD, double *E, unsigned int N, unsigned int istart, unsigned int ppG, double plummer_core, double plummer_mass){

	unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

	extern __shared__ double4 shPosition[];
	double4 *shPos = (double4*)shPosition;
   float4 *shVel = (float4*)&shPos[blockDim.x];
	
	double4 myposition = posCD[i];
	double4 myvelocity = velCD[i];
   float myEpsilon = velPD[i].w;

	double kin = 0.0, pE = 0.0;
	int k = 0, j;
	int totgpus = N/ppG;

   kin = 0.5/totgpus * myposition.w *(myvelocity.x*myvelocity.x +
                               myvelocity.y*myvelocity.y +
                               myvelocity.z*myvelocity.z);

	while(k < ppG){
      j = 0;
      shPos[threadIdx.x] = posCD[k + threadIdx.x + istart];
      shVel[threadIdx.x] = velPD[k + threadIdx.x + istart];
		__syncthreads();

      while(j < blockDim.x){
         double4 dr = {shPos[j].x - myposition.x, shPos[j].y - myposition.y, shPos[j].z - myposition.z, 0.0};
         double distance = (dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
         if(distance != 0.0)
            pE -= 0.5 * myposition.w * shPos[j].w * rsqrt(distance + myEpsilon * shVel[j].w);
      j++;
      }
   __syncthreads();
   k += blockDim.x;
   }

	AddPlummerEnergy(myposition, &pE, &plummer_core, &plummer_mass, totgpus);

	E[i] = kin + pE;
}

__global__ void potential_energy(double4 *posCD, float4 *velPD, double *E, unsigned int N, unsigned int istart, unsigned int ppG, double plummer_core, double plummer_mass){

	unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

	extern __shared__ double4 shPosition[];
	double4 *shPos = (double4*)shPosition;
   float4 *shVel = (float4*)&shPos[blockDim.x];

	double4 myposition = posCD[i];
   float myEpsilon = velPD[i].w;

	double pE = 0.0;
	int k = 0, j;
	int totgpus = N/ppG;

	while(k < ppG){
		j = 0;
		shPos[threadIdx.x] = posCD[k + threadIdx.x + istart];
		shVel[threadIdx.x] = velPD[k + threadIdx.x + istart];
		__syncthreads();

		while(j < blockDim.x){
			double4 dr = {shPos[j].x - myposition.x, shPos[j].y - myposition.y, shPos[j].z - myposition.z, 0.0};
			double distance = (dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
			if(distance != 0.0)
				pE -= 0.5 * myposition.w * shPos[j].w * rsqrt(distance + myEpsilon * shVel[j].w);
			j++;	
		}
		__syncthreads();
		k += blockDim.x;
	}

	AddPlummerEnergy(myposition, &pE, &plummer_core, &plummer_mass, totgpus);

	E[i] = pE;
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
									 double pl_core,
									 double pl_mass)
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
		if(whoami == 0)
			loc_time[whop] = TIME;

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
        float distance = (dr.x * dr.x + dr.y * dr.y + dr.z * dr.z) + myVelocity.w * shVel[j].w;
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

        j++;
	  }

	  __syncthreads();
     i += blockDim.x;
	}

		
	AddPlummer(myPosition, myVelocity, myAccelera, &acc, &jrk, &snp, &pl_core, &pl_mass, tot_mul);


   aC[gtid] = acc;
	aC1[gtid] = jrk;
	aC2[gtid] = snp;

}


