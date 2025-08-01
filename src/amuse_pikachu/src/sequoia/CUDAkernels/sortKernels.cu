#include "scanKernels.cu"
#include "support_kernels.cu"

//Helper functions

//Reorders data 
extern "C" __global__ void dataReorderR4(const int n_particles,
                                         real4 *source,
                                         real4 *destination,
                                         uint  *permutation) {
  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;
  const int dim =  blockDim.x * blockDim.y;
  
  int idx = bid * dim + tid;
  if (idx >= n_particles) return;

   int newIndex = permutation[idx];
   destination[idx] = source[newIndex];  
}

extern "C" __global__ void dataReorderF2(const int n_particles,
                                         float2 *source,
                                         float2 *destination,
                                         uint  *permutation) {
  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;
  const int dim =  blockDim.x * blockDim.y;
  
  int idx = bid * dim + tid;
  if (idx >= n_particles) return;

  int newIndex = permutation[idx];
  destination[idx] = source[newIndex];  
}

extern "C" __global__ void dataReorderI1(const int n_particles,
                                         int *source,
                                         int *destination,
                                         uint  *permutation) {
  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;
  const int dim =  blockDim.x * blockDim.y;
  
  int idx = bid * dim + tid;
  if (idx >= n_particles) return;

   int newIndex = permutation[idx];
   destination[idx] = source[newIndex];  
}


//Convert a 64bit key uint2 key into a 96key with a permutation value build in
extern "C" __global__ void convertKey64to96(uint4 *keys,  uint4 *newKeys, const int N)
{
  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;
  const int dim =  blockDim.x * blockDim.y;

  int idx = bid * dim + tid;

  if (idx >= N) return;

  uint4 temp = keys[idx];
  newKeys[idx] = (uint4){temp.x, temp.y, temp.z, idx};
}

extern "C" __global__ void extractKeyAndPerm(uint4 *newKeys, uint4 *keys, uint *permutation, const int N)
{
  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;
  const int dim =  blockDim.x * blockDim.y;

  int idx = bid * dim + tid;

  if (idx >= N) return;

  uint4 temp = newKeys[idx];
  
  keys[idx]        = (uint4){temp.x, temp.y, temp.z, temp.w};
  permutation[idx] = temp.w;
}


//Extract 1 of the 4 items of an uint4 key and move it into a 32bit array
extern "C" __global__ void extractInt(uint4 *keys,  uint *simpleKeys, const int N, int keyIdx)
{
  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;
  const int dim =  blockDim.x * blockDim.y;

  int idx = bid * dim + tid;

  if (idx >= N) return;

  uint4 temp = keys[idx];
  int  simpleTemp;

  if(keyIdx == 0)
      simpleTemp = temp.x;
  else if(keyIdx == 1)
      simpleTemp = temp.y;
  else if(keyIdx == 2)
      simpleTemp = temp.z;

  simpleKeys[idx] = simpleTemp;
}

//Create range of 0 to N
extern "C" __global__ void fillSequence(uint *sequence, const int N)
{
  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;
  const int dim =  blockDim.x * blockDim.y;

  int idx = bid * dim + tid;

  if (idx >= N) return;

  sequence[idx] = idx;
}

//Reorder the data in the arrays according to a given permutation
extern "C" __global__ void reOrderKeysValues(uint4 *keysSrc, uint4 *keysDest, uint *permutation, const int N)
{
  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;
  const int dim =  blockDim.x * blockDim.y;

  int idx = bid * dim + tid;

  if (idx >= N) return;

  int newIndex = permutation[idx];
  keysDest[idx] = keysSrc[newIndex];
}

extern "C" __global__ void sort_count(volatile uint2 *valid, int *counts, const int N, setupParams sParam, int bitIdx/*, int2 *blaat*/)
{
  const int tid    =  threadIdx.x;
  const int bid    =  blockDim.y *  blockIdx.x + threadIdx.y;

  int totalNumThreads = gridDim.x*blockDim.y*blockDim.x; //120*4*32 // gridDim.x * blockDim.y; //2D !!!!

  volatile __shared__ int shmemSC[128];
  volatile __shared__ int shmemSCTEST[128];

  //Determine the parameters and loop over the particles
  int jobSize = (N / 2) / totalNumThreads;
  int offSet  = jobSize * bid;
  int count   = 0;

  jobSize = sParam.jobs;
  if(bid < sParam.blocksWithExtraJobs)
    jobSize++;

  if(bid <= sParam.blocksWithExtraJobs)
    offSet = (sParam.jobs+1)*64*bid;
  else
  {
    offSet = sParam.blocksWithExtraJobs*(sParam.jobs+1)*64;
    offSet += (bid-sParam.blocksWithExtraJobs)*(sParam.jobs)*64;
  }

  offSet /= 2;  //Divide by two since we do double loads (uint2)

  for(int i=0; i < jobSize; i++)
  {   
    count  += !(valid[offSet + tid].x & (1u<<bitIdx));
    count  += !(valid[offSet + tid].y & (1u<<bitIdx));
    offSet += blockDim.x;
  }

  //Reduce to get the count of this block
  shmemSC[32*threadIdx.y + tid] = count;
  reduce_block2(tid, &shmemSC[32*threadIdx.y], count);

  //Save the values / count of the current block
  if(threadIdx.x == 0)
    counts[bid] = shmemSC[32*threadIdx.y];

  //Block 0 handles any extra elements that couldn't be divided equally
  if(bid == 0)
  {
    //Here i use single element reads for ease of boundary conditions and steps
    count   = 0;
    offSet  = sParam.extraOffset;

    uint* valid2 = (uint*) valid;

    for(int i=0 ; i < sParam.extraElements;  i += blockDim.x)
    {
      if((offSet + i +  tid) < (N))  //Make sure we dont read more than there are items
      {
        count += !(valid2[offSet + i +  tid] & (1u<<bitIdx));
      }
    }

    //Reduce
    shmemSCTEST[tid] = count;

    __syncthreads();

    if(tid < 16){
      shmemSCTEST[tid] = count = count + shmemSCTEST[tid+16];
      shmemSCTEST[tid] = count = count + shmemSCTEST[tid+8];
      shmemSCTEST[tid] = count = count + shmemSCTEST[tid+4];
      shmemSCTEST[tid] = count = count + shmemSCTEST[tid+2];
      shmemSCTEST[tid] = count = count + shmemSCTEST[tid+1]; 
    }

    //Save the count
    if(tid == 0)
    {
      counts[gridDim.x*blockDim.y] = shmemSCTEST[0];
    }

    __syncthreads();

  }//end if  bid==0 
}//end compact_count


// __device__  __forceinline__ int testTest(volatile unsigned int tmp[], uint val, const int idx, long test)
// {
//   tmp[idx-16] = 0; tmp[idx] = val;
// 
//   // Since we set half the array to 0 we don't need ifs!
//   tmp[idx] = val = tmp[idx -  1]  + val;
//   tmp[idx] = val = tmp[idx -  2]  + val;
//   tmp[idx] = val = tmp[idx -  4]  + val;
//   tmp[idx] = val = tmp[idx -  8]  + val;
//   tmp[idx] = val = tmp[idx -  16] + val;
// 
//   return (idx > 0) ? tmp[idx-1] : 0;
// }


/*
For sorting it turns out that the stage kernels works faster than the non-staged
Might depend on how much has to be sorted/moved, have to do timings in the actual code
*/
extern "C" __global__ void sort_move_stage_key_value(uint2 *valid, int *output,
                                          uint2 *srcValues, uint *valuesOut,
                                          int *counts,
                                          const int N, setupParams sParam, int bitIdx)
{
  //Walk the values of this block
  const int tid    =  threadIdx.x;
  const int bid    =  blockDim.y *  blockIdx.x + threadIdx.y;

  volatile __shared__ unsigned int shmemSMSKV[192];
  volatile __shared__ int stage[64*4];
  volatile __shared__ int stage_values[64*4];

  //Determine the parameters and loop over the particles
  int jobSize, offSet;


  jobSize = sParam.jobs;
  if(bid < sParam.blocksWithExtraJobs)
    jobSize++;

  if(bid <= sParam.blocksWithExtraJobs)
    offSet = (sParam.jobs+1)*64*bid;
  else
  {
    offSet = sParam.blocksWithExtraJobs*(sParam.jobs+1)*64;
    offSet += (bid-sParam.blocksWithExtraJobs)*(sParam.jobs)*64;
  }

  int outputOffset = counts[bid];

  //Get the start of the output offset of the invalid items
  //this is calculated as follows:
  //totalValidItems + startReadOffset - startOutputOffset
  //startReadOffset - startOutputOffset <- is the total number of invalid items from any blocks
  //before the current block
  int rightOutputOffset = counts[gridDim.x*blockDim.y+1];
  rightOutputOffset     = rightOutputOffset + offSet - outputOffset;

  offSet /= 2;  //Divide by two since we do double loads (uint2) TODO what happens if offSet is uneven...?

  int curCount;
  int idx, ridx;

  outputOffset      += threadIdx.x;
  rightOutputOffset += threadIdx.x;

  //Do per step the prefix scan to determine the output locations
  for(int i=0; i < jobSize; i++)
  {
    uint2  validBase  = valid[offSet + tid];
    uint2  valuesBase = srcValues[offSet + tid];
    int value         = !(validBase.x  & (1u<<bitIdx));
    value            += !(validBase.y  & (1u<<bitIdx));

    idx  = hillisSteele5(&shmemSMSKV[48*threadIdx.y+16], curCount, value, threadIdx.x);

    ridx = curCount + threadIdx.x*2 - idx; //lane*2 - idx , *2 since we read 2 items a time

    if(!(validBase.x  & (1u<<bitIdx)))
    {
      stage[idx + threadIdx.y*64]          = validBase.x;
      stage_values[idx++ + threadIdx.y*64] = valuesBase.x;
    }
    else
    {
      stage[ridx + threadIdx.y*64]          = validBase.x;
      stage_values[ridx++ + threadIdx.y*64] = valuesBase.x;
    }

    if(!(validBase.y  & (1u<<bitIdx)))
    {
      stage[idx + threadIdx.y*64]        = validBase.y;
      stage_values[idx + threadIdx.y*64] = valuesBase.y;
    }
    else
    {
      stage[ridx + threadIdx.y*64]        = validBase.y;
      stage_values[ridx + threadIdx.y*64] = valuesBase.y;
    }

    //Reuse value as index
    value = outputOffset;
    //Flush output, first 32
    if(threadIdx.x >= curCount)
      value = rightOutputOffset-curCount;
    output[value]    = stage       [threadIdx.x + threadIdx.y*64];
    valuesOut[value] = stage_values[threadIdx.x + threadIdx.y*64];

    //2nd 32
    value = outputOffset + blockDim.x;
    if(threadIdx.x + blockDim.x >= curCount)
      value = rightOutputOffset + blockDim.x - curCount;

    output[value]    = stage       [threadIdx.x + blockDim.x + threadIdx.y*64];
    valuesOut[value] = stage_values[threadIdx.x + blockDim.x + threadIdx.y*64];

    outputOffset      += curCount;      //Increase the output offset
    rightOutputOffset += 64 - curCount; //64 (32*2) since we do 2 items a time
    offSet            += blockDim.x;    //Step to the next N threads
  }

  //Block 0 handles any extra elements that couldn't be divided equally
  if(bid == 0)
  {
    //Here i use single element reads for ease of boundary conditions and steps
    offSet              = sParam.extraOffset;
    outputOffset        = counts[gridDim.x*blockDim.y];
    rightOutputOffset   = counts[gridDim.x*blockDim.y+1];
    rightOutputOffset   = rightOutputOffset + offSet - outputOffset;

    uint* valid2 = (uint*) valid;
    uint* srcValues2 = (uint*) srcValues;

    for(int i=0; i < sParam.extraElements;  i += blockDim.x)
    {
      uint value = 0;
      uint srcValueItem = 0;
    
      if((offSet + i +  tid) < (N)){  //Make sure we dont read more than there are items
        value        = valid2[offSet + i +  tid];
        srcValueItem = srcValues2[offSet + i +  tid];
      }

      idx  = hillisSteele5(&shmemSMSKV[48*threadIdx.y+16], curCount, !(value & (1u<<bitIdx)), threadIdx.x);
      ridx = threadIdx.x - idx;

      if((offSet + i +  tid) < N)
        if(!(value & (1u<<bitIdx)))
        {
          output[idx + outputOffset]    = value;
          valuesOut[idx + outputOffset] = srcValueItem;
        }
        else
        {
          output[ridx + rightOutputOffset]     = value;
          valuesOut[ridx + rightOutputOffset]  = srcValueItem;
        }

      outputOffset      += curCount;       //Increase the output offset
      rightOutputOffset += 32-curCount;    //32 since we do only 1 at a time
    }
  }//end if bid==0 
}//end sort_move_stage_key_value

