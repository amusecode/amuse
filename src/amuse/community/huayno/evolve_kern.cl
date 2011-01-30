#define RVTIMESTEP
#define RATIMESTEP
#define RVRARATIO   1.

#define FLOAT double
#define FLOAT4 double4
#define BIGNUM HUGE_VAL

#pragma OPENCL EXTENSION cl_amd_fp64 : enable 

__kernel void kick_kernel(
  uint nj, 
  int blocksize,  
  FLOAT eps2,
  __global FLOAT4* ipos,
  __global FLOAT4* jpos,
  __global FLOAT4* acc,
  __local FLOAT4* pblock)
{
  int global_id = get_global_id(0);
  int local_id = get_local_id(0);
  int global_n = get_global_size(0);
  int local_n = get_local_size(0);
  int nblocks = (nj-1)/blocksize+1; 
  int nb;
  
  FLOAT4 p = ipos[global_id];   
  FLOAT4 a = (FLOAT4) (0.0,0.0,0.0,0.0);
    
  for(int jb=0; jb < nblocks; jb++) 
  { 
    nb=blocksize; if((jb+1)*blocksize > nj ) nb=nj-jb*blocksize;
    event_t e=async_work_group_copy(pblock,jpos+jb*blocksize, nb,0);
    wait_group_events(1,&e);
    
    for(int j=0; j<nb; j++) 
    {
      FLOAT4 p2 = pblock[j];
      FLOAT4 d = p2-p;
      FLOAT dr2=d.x*d.x + d.y*d.y + d.z*d.z + eps2;
      if(dr2 > 0)
      {
        FLOAT invr = rsqrt(dr2);
        FLOAT invr3 = p2.w*invr*invr*invr;
        FLOAT4 f = (FLOAT4) (invr3,invr3,invr3,0.0);
        a += f*d; 
      }  
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }    
  acc[global_id] = a;
}

__kernel void timestep_kernel(
  uint nj,
  int blocksize,
  FLOAT eps2,
  FLOAT dt_param,
  __global FLOAT4* ipos,
  __global FLOAT4* ivel,
  __global FLOAT4* jpos,
  __global FLOAT4* jvel,
  __global FLOAT* timestep,
  __local FLOAT4* pblock,
  __local FLOAT4* vblock)
{
  int global_id = get_global_id(0);
  int local_id = get_local_id(0);
  int global_n = get_global_size(0);
  int local_n = get_local_size(0);
  int nblocks = (nj-1)/blocksize+1; 
  int nb;
  
  FLOAT4 p = ipos[global_id];   
  FLOAT4 v = ivel[global_id];   
  FLOAT ts = BIGNUM;
  
  for(int jb=0; jb < nblocks; jb++) 
  { 
    event_t e[2];
    nb=blocksize; if((jb+1)*blocksize > nj ) nb=nj-jb*blocksize;
    e[0]=async_work_group_copy(pblock,jpos+jb*blocksize, nb,0);
    e[1]=async_work_group_copy(vblock,jvel+jb*blocksize, nb,0);
    wait_group_events(2,e);
    
    for(int j=0; j<nb; j++) 
    {
      FLOAT4 p2 = pblock[j];
      FLOAT4 v2 = vblock[j];
      FLOAT4 d = p2-p;
      FLOAT dr2=d.x*d.x + d.y*d.y + d.z*d.z + eps2;
      if(dr2 > 0)
      {
        FLOAT dr = sqrt(dr2);
        FLOAT dr3 = dr*dr2;
        FLOAT4 dv = v2-v;
        FLOAT vdotdr2=(dv.x*d.x + dv.y*d.y + dv.z*d.z)/dr2;
        FLOAT dv2=dv.x*dv.x + dv.y*dv.y + dv.z*dv.z;
        FLOAT mu=p.w+p2.w;
        
#ifdef RATIMESTEP
        FLOAT tau=RVRARATIO*dt_param/M_SQRT2*sqrt(dr3/mu);
        FLOAT dtau=3/2.*tau*vdotdr2;
        if(dtau>1.) dtau=1.;
        tau/=(1-dtau/2);
        ts=fmin(ts,tau);
#endif

        if(dv2>0)
        {
#ifdef RVTIMESTEP
          tau=dt_param*dr/sqrt(dv2);
          dtau=tau*vdotdr2*(1+mu/(dv2*dr));
          if(dtau>1.) dtau=1.;
          tau/=(1-dtau/2);
          ts=fmin(ts,tau);
#endif
        }
      }  
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }    
  timestep[global_id] = ts;
}


__kernel void potential_kernel(
  uint nj, 
  int blocksize,  
  FLOAT eps2,
  __global FLOAT4* ipos,
  __global FLOAT4* jpos,
  __global FLOAT* potential,
  __local FLOAT4* pblock)
{
  int global_id = get_global_id(0);
  int local_id = get_local_id(0);
  int global_n = get_global_size(0);
  int local_n = get_local_size(0);
  int nblocks = (nj-1)/blocksize+1; 
  int nb;
  
  FLOAT4 p = ipos[global_id];   
  FLOAT pot = 0.0;
    
  for(int jb=0; jb < nblocks; jb++) 
  { 
    nb=blocksize; if((jb+1)*blocksize > nj ) nb=nj-jb*blocksize;
    event_t e=async_work_group_copy(pblock,jpos+jb*blocksize, nb,0);
    wait_group_events(1,&e);
    
    for(int j=0; j<nb; j++) 
    {
      FLOAT4 p2 = pblock[j];
      FLOAT4 d = p2-p;
      FLOAT dr2=d.x*d.x + d.y*d.y + d.z*d.z + eps2;
      if(dr2 > 0)
      {
        FLOAT invr = rsqrt(dr2);
        pot-= p2.w*invr;
      }  
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }    
  potential[global_id] = pot;
}

