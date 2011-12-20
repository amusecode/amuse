#include "support_kernels.cu"

#include <stdio.h>


__device__ int isinbox(real4 pos, double4 xlow, double4 xhigh)
{  
    if((pos.x < xlow.x)||(pos.x > xhigh.x))          
      return 0;
    if((pos.y < xlow.y)||(pos.y > xhigh.y))          
      return 0;
    if((pos.z < xlow.z)||(pos.z > xhigh.z))          
      return 0;
    
    return 1;
}


extern "C" __global__ void doDomainCheck(int    n_bodies,
                                           double4  xlow,
                                           double4  xhigh,
                                           real4  *body_pos,
                                           int    *validList    //Valid is 1 if particle is outside domain
){
  
  uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  uint tid = threadIdx.x;
  uint id  = bid * blockDim.x + tid;
  
  if (id >= n_bodies) return;

  real4 pos = body_pos[id];

  int valid      = isinbox(pos, xlow, xhigh);
  valid = !valid;
  validList[id] = id | ((valid) << 31);
}
  

//Checks the domain and computes the key list
//if a particle is outside the domain it gets a special key
//otherwise the normal key is used
extern "C" __global__ void doDomainCheckAdvanced(int    n_bodies,
                                           double4  xlow,
                                           double4  xhigh,
                                           real4  *body_pos,
                                           int    *validList    //Valid is 1 if particle is outside domain
){
  
  uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  uint tid = threadIdx.x;
  uint id  = bid * blockDim.x + tid;
  
  if (id >= n_bodies) return;

  real4 pos = body_pos[id];

  int valid      = isinbox(pos, xlow, xhigh);
  valid = !valid;
  validList[id] = id | ((valid) << 31);
}
  

extern "C" __global__ void extractSampleParticles(int    n_bodies,
                                                  int    sample_freq,
                                                  real4  *body_pos,
                                                  real4  *samplePosition
){
  uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  uint tid = threadIdx.x;
  uint id  = bid * blockDim.x + tid;
  

  int idx  = id*sample_freq;
  if  (idx >= n_bodies) return;

  samplePosition[id] =  body_pos[idx];
}

extern "C" __global__ void extractOutOfDomainParticlesR4(int n_extract,
                                                       int *extractList,
                                                       real4 *source,
                                                       real4 *destination)
{
  uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  uint tid = threadIdx.x;
  uint id  = bid * blockDim.x + tid;

  if(id >= n_extract) return;

  destination[id] = source[extractList[id]];

}



typedef struct bodyStruct
{
  real4 pos;
  real4 vel;
  real4 acc0;
  real4 acc1;
  real4 Ppos;
  real4 Pvel;
  float2 time;
  int   id;
  int   temp;
} bodyStruct;


extern "C" __global__ void extractOutOfDomainParticlesAdvanced(int n_extract,
                                                       int *extractList,
                                                       real4 *Ppos,
                                                       real4 *Pvel,
                                                       real4 *pos,
                                                       real4 *vel,
                                                       real4 *acc0,
                                                       real4 *acc1,
                                                       float2 *time,
                                                       int   *body_id,
                                                       bodyStruct *destination)
{
  uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  uint tid = threadIdx.x;
  uint id  = bid * blockDim.x + tid;

  if(id >= n_extract) return;

  //copy the data from a struct of arrays into a array of structs
  destination[id].Ppos = Ppos[extractList[id]];
  destination[id].Pvel = Pvel[extractList[id]];
  destination[id].pos  = pos[extractList[id]];
  destination[id].vel  = vel[extractList[id]];
  destination[id].acc0  = acc0[extractList[id]];
  destination[id].acc1  = acc1[extractList[id]];
  destination[id].time  = time[extractList[id]];
  destination[id].id    = body_id[extractList[id]];

}


extern "C" __global__ void internalMove(int       n_extract,                                       
                                        int       n_bodies,
                                        double4  xlow,
                                        double4  xhigh,
                                        int       *extractList,
                                        int       *indexList,
                                        real4     *Ppos,
                                        real4     *Pvel,
                                        real4     *pos,
                                        real4     *vel,
                                        real4     *acc0,
                                        real4     *acc1,
                                        float2    *time,
                                        int       *body_id)
{
  uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  uint tid = threadIdx.x;
  uint id  = bid * blockDim.x + tid;

  if(id >= n_extract) return;

  int srcIdx     = (n_bodies-n_extract) + id;
  real4 testpos  = Ppos[srcIdx];

  if(isinbox(testpos, xlow, xhigh))
  {
    int dstIdx = atomicAdd(indexList, 1);    
    dstIdx     = extractList[dstIdx];

    //Move!
    Ppos[dstIdx] = Ppos[srcIdx];
    Pvel[dstIdx] = Pvel[srcIdx];
    pos[dstIdx]  = pos[srcIdx];
    vel[dstIdx]  = vel[srcIdx];
    acc0[dstIdx] = acc0[srcIdx];
    acc1[dstIdx] = acc1[srcIdx];
    time[dstIdx] = time[srcIdx];
    body_id[dstIdx] = body_id[srcIdx];
  }//if isinbox

}

extern "C" __global__ void insertNewParticles(int       n_extract,
                                              int       n_insert,
                                              int       n_oldbodies,
                                              int       offset,
                                              real4     *Ppos,
                                              real4     *Pvel,
                                              real4     *pos,
                                              real4     *vel,
                                              real4     *acc0,
                                              real4     *acc1,
                                              float2    *time,
                                              int       *body_id,
                                              bodyStruct *source)
{
  uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  uint tid = threadIdx.x;
  uint id  = bid * blockDim.x + tid;

  if(id >= n_insert) return;

  //The newly added particles are added at the end of the array
  int idx = (n_oldbodies-n_extract) + id + offset;

  //copy the data from a struct of arrays into a array of structs
  Ppos[idx]     = source[id].Ppos;
  Pvel[idx]     = source[id].Pvel;
  pos[idx]      = source[id].pos;
  vel[idx]      = source[id].vel;
  acc0[idx]     = source[id].acc0;
  acc1[idx]     = source[id].acc1;
  time[idx]     = source[id].time;
  body_id[idx]  = source[id].id;
}





// extern "C" __global__ void insertNewParticles(int       n_extract,
//                                               int       n_insert,
//                                               int       n_oldbodies,
//                                               int       *extractList,
//                                               real4     *Ppos,
//                                               real4     *Pvel,
//                                               real4     *pos,
//                                               real4     *vel,
//                                               real4     *acc0,
//                                               real4     *acc1,
//                                               float2    *time,
//                                               int       *body_id,
//                                               bodyStruct *source)
// {
//   uint bid = blockIdx.y * gridDim.x + blockIdx.x;
//   uint tid = threadIdx.x;
//   uint id  = bid * blockDim.x + tid;
// 
//   int idx, srcidx = -1; 
// /*
// 
// //Situaties:
// - n_insert > n_extract -> particles moeten aan einde worden toegevoegd (meer toevoegen dan weggehaald)
//     id < n_extract -> idx = extractList[id]  ; uit source[id]
//     id >= n_extract & id < n_insert  --> idx = n_oldbodies + (id-n_extract); uit source[id]
//   
// - n_insert <= n_exract -> particles moeten van het einde naar het begin (meer verwijderd dan toegevoegd)
//     id < n_extract -> idx = extractList[id] ; uit source[id]
//     id >= n_extract & id < n_insert -> idx = extractList[id] ; uit dest[n_bodies-(n_extract-n_insert) + (id - n_insert)]
// 
//   */
// 
//   if(n_insert > n_extract)
//   {
//     if(id < n_extract)
//     {
//        idx = extractList[id];
//     }
//     else if(id >= n_extract && id < n_insert)
//     {
//       //Insert particles at the end of the array
//       idx = n_oldbodies + (id-n_extract);
//     }
//     else
//     {
//       return;
//     }
//   }
//   else
//   {
//     //n_insert <= n_extract
// 
//     if(id < n_insert)
//     {
//        idx = extractList[id];
//     }
//     else if(id >= n_insert && id < n_extract)
//     {
//       //Move particles from the back of the array to the empty spots
//       idx    = extractList[id];
//       srcidx = extractList[n_oldbodies-(n_extract-n_insert) + (id - n_insert)];
//     //  srcidx = n_oldbodies-(n_extract-n_insert) + (id - n_insert);
//     }
//     else
//     {
//       return;
//     }
//   }
// /*
// Gaat niet goed als n_insert < n_extract
// omdat we als we gaan moven we ook kans hebben dat we iets moven
// van het begin naar het eind als daar iets is uitgehaald
// we zouden dus de laatste verwijderde moeten vinden en zorgen dat er neits achter komt ofzo
// 
// 
// 
// */
// 
// /*
//   if(id < n_extract)
//   {
//     idx = extractList[id];
//   }
//   else if(id >= n_extract && id < n_insert)
//   {
//     if(n_insert > n_extract)
//     {
//       //Insert particles at the end of the array
//       idx = n_oldbodies + (id-n_extract);
//     }
//     else
//     {
//       //Move particles from the back of the array to the empty spots
//       idx    = extractList[id];
//       srcidx = n_oldbodies-(n_extract-n_insert) + (id - n_insert);
//     }
//   }
//   else
//   {
//     //Outside all array ranges
//     return;
//   }*/
// 
// 
//   if(srcidx < 0)
//   {
//     //copy the data from a struct of arrays into a array of structs
//     Ppos[idx] = source[id].Ppos;
//     Pvel[idx] = source[id].Pvel;
//     pos[idx]  = source[id].pos;
//     vel[idx]  = source[id].vel;
//     acc0[idx] = source[id].acc0;
//     acc1[idx] = source[id].acc1;
//     time[idx] = source[id].time;
//     body_id[idx] = source[id].id;
// 
// printf("%d  (CMOVE external %d) goes to: %d \n", source[id].id,n_insert, idx);
// 
// 
//   }
//   else
//   {
//     Ppos[idx] = Ppos[srcidx];
//     Pvel[idx] = Pvel[srcidx];
//     pos[idx]  = pos[srcidx];
//     vel[idx]  = vel[srcidx];
//     acc0[idx] = acc0[srcidx];
//     acc1[idx] = acc1[srcidx];
//     time[idx] = time[srcidx];
//   int temp = body_id[idx];
//     body_id[idx] = body_id[srcidx];
// 
// printf("%d stored at: %d (CMOVE internal %d) goes to: %d  overwr: %d \n", body_id[srcidx],srcidx, n_insert, idx, temp);
// 
// 
//   }//if srcidx < 0
// 
// 
// }


