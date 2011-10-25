typedef unsigned int uint;


//Warp based summation
__device__ int inexclusive_scan_warp(volatile int *ptr,bool inclusive, const unsigned int idx, int value) {
  const unsigned int lane = idx & 31;
  
  if (lane >=  1) ptr[idx] = value = ptr[idx -  1]   + value;
  if (lane >=  2) ptr[idx] = value = ptr[idx -  2]   + value;
  if (lane >=  4) ptr[idx] = value = ptr[idx -  4]   + value;
  if (lane >=  8) ptr[idx] = value = ptr[idx -  8]   + value;
  if (lane >= 16) ptr[idx] = value = ptr[idx -  16]  + value;

  if(inclusive)  
    return value;    //Inclusive
  else
    return(lane > 0) ? ptr[idx-1] : 0;  //Exclusive
}



//N is number of previous blocks in the count call
__device__ void exclusive_scan_blockD(int *ptr, const int N, int *count, volatile int *shmemESB) 
{
  const unsigned int idx = threadIdx.x;
  const unsigned int lane   = idx & 31;
  const unsigned int warpid = idx >> 5;

  int value;

  //Read the data in shmem
  if(idx < N + 1)
    shmemESB[idx] = value = ptr[idx];
  else
    shmemESB[idx] = value = 0;
  
  __syncthreads();

  // step 1: Intra-warp scan in each warp
  int val = inexclusive_scan_warp(&shmemESB[0], false, idx, value);
  __syncthreads();

  // step 2: Collect per-warp particle results
  if (lane == 31) shmemESB[warpid] = shmemESB[idx];
  __syncthreads();
  value = shmemESB[idx];
  // step 3: Use 1st warp to scan per-warp results
  if (warpid == 0) inexclusive_scan_warp(&shmemESB[0], false, idx, value);
  __syncthreads();
  
  // step 4: Accumulate results from Steps 1 and 3;
  if (warpid > 0) val = shmemESB[warpid - 1] + val;
  __syncthreads();

  // Step 5: Write and return the final result
  ptr[idx] = val;
  __syncthreads();

  // ptr[blockDim.x - 1] + lastValue; //count
  if(idx == 0)//Thread 0 saves the total count value
    *count = ptr[blockDim.x - 1];
}

//N is number of previous blocks in the count call
extern "C"  __global__ void exclusive_scan_block(int *ptr, const int N, int *count) 
{
  extern __shared__ int shmemESB[];
  exclusive_scan_blockD(ptr, N, count, shmemESB);
}


typedef struct setupParams
{
  int jobs;                     //Minimal number of jobs for each 'processor'
  int blocksWithExtraJobs;      //Some ' processors'  do one extra job all with bid < bWEJ
  int extraElements;            //The elements that didn't fit completely
  int extraOffset;              //Start of the extra elements
}setupParams;


//Warp based prefix sum, using extra buffer space to remove the need for if statements
// __device__  int hillisSteele4(volatile int *ptr, int *count, uint val, const unsigned int idx)
__device__  int hillisSteele4(volatile int *ptr, int &count, uint val, const unsigned int idx)
{
  //  const unsigned int lane   = idx & 31;  
  //We don't require lane here since idx is always < 32 in the way we start the blocks/threads

  //volatile int* tmp = ptr + (32 / 2);
  volatile int* tmp = &ptr[16];
  ptr[idx] = 0; tmp[idx] = val;

  //Since we set half the array to 0 we don't need ifs!
  tmp[idx] = val = tmp[idx -  1]  + val;
  tmp[idx] = val = tmp[idx -  2]  + val;
  tmp[idx] = val = tmp[idx -  4]  + val;
  tmp[idx] = val = tmp[idx -  8]  + val;
  tmp[idx] = val = tmp[idx -  16] + val;

  //Inclusive sum/count
  count = tmp[blockDim.x-1];

  //Exclusive index
  return (idx > 0) ? tmp[idx-1] : 0;
}

//Warp based prefix sum, using extra buffer space to remove the need for if statements
// __device__  int hillisSteele4(volatile int *ptr, int *count, uint val, const unsigned int idx)
__device__  __forceinline__ int hillisSteele5(volatile unsigned int tmp[], int &count, uint val, const int idx)
{
  //  const unsigned int lane   = idx & 31;  
  //We don't require lane here since idx is always < 32 in the way we start the blocks/threads

  //volatile int* tmp = ptr + (32 / 2);
//   volatile int* tmp = &ptr[16];
 
  tmp[idx-16] = 0; tmp[idx] = val;
 

  //Since we set half the array to 0 we don't need ifs!
  tmp[idx] = val = tmp[idx -  1]  + val;
  tmp[idx] = val = tmp[idx -  2]  + val;
  tmp[idx] = val = tmp[idx -  4]  + val;
  tmp[idx] = val = tmp[idx -  8]  + val;
  tmp[idx] = val = tmp[idx -  16] + val;

  //Inclusive sum/count
  count = tmp[blockDim.x-1];

  //Exclusive index
  return (idx > 0) ? tmp[idx-1] : 0;
}


__device__ void reduce_block2(int tid, volatile int *shmem, int val)
{
  //Reduce the 32 block
  if(tid < 16){
    shmem[tid] = val = val + shmem[tid+16];
    shmem[tid] = val = val + shmem[tid+8];
    shmem[tid] = val = val + shmem[tid+4];
    shmem[tid] = val = val + shmem[tid+2];
    shmem[tid] = val = val + shmem[tid+1];  
  }
}


//Count the number of valid elements in this BLOCK
__device__ void compact_countD(volatile uint2 *values,
                              uint *counts,  
                              const int N,                             
                              setupParams sParam, volatile int *shmemCC2) {
  const int tid    =  threadIdx.x;
  const int bid    =  blockDim.y *  blockIdx.x + threadIdx.y;


  volatile __shared__ int shmemCC[128];

  //Determine the parameters and loop over the particles
  int jobSize, offSet, count = 0;

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
    count  += (values[offSet + tid].x >> 31);
    count  += (values[offSet + tid].y >> 31);
    offSet += blockDim.x;
  }


  //Reduce to get the count of this block
  shmemCC[32*threadIdx.y + tid] = count;
  __syncthreads();
  reduce_block2(tid, &shmemCC[32*threadIdx.y], count);
  
  //Save the values / count of the current block
  if(threadIdx.x == 0)  
    counts[bid] = shmemCC[32*threadIdx.y];

  //Block 0 handles any extra elements that couldn't be divided equally
  if(bid == 0)
  {
    //Here i use single element reads for ease of boundary conditions and steps
    count   = 0;
    offSet  = sParam.extraOffset;   

    uint* value2 = (uint*) values;

    for(int i=0 ; i < sParam.extraElements;  i += blockDim.x)
    {
      if((offSet + i +  tid) < (N))  //Make sure we dont read more than there are items
      {
        count += (value2[offSet + i +  tid] >> 31);
      }
    }

    //Reduce
    shmemCC[tid] = count;
    __syncthreads();
    reduce_block2(tid, &shmemCC[0], count);
  
    //Save the count
    if(tid == 0)  
      counts[gridDim.x*blockDim.y] = shmemCC[0];

  }//end if  bid==0 
}//end compact_count

//Count the number of valid elements in this BLOCK
extern "C" __global__ void compact_count(volatile uint2 *values,
                              uint *counts,  
                              const int N,                             
                              setupParams sParam) {


  extern __shared__  int shmemCC[];
  compact_countD(values, counts, N, sParam,shmemCC) ;
}


//The kernel that actually moves the data
__device__ void compact_moveD( uint2 *values,
                             uint *output, 
                             uint *counts,  
                             const int N,
                            setupParams sParam, volatile unsigned int *shmemCM)
{
  //Walk the values of this block
  const int tid    =  threadIdx.x;
  const int bid    =  blockDim.y *  blockIdx.x + threadIdx.y;


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

  offSet /= 2;  //Divide by two since we do double loads (uint2) TODO what happens if offSet is uneven...?
  
  int outputOffset = counts[bid];
  int curCount;

  //Do per step the prefix scan to determine the output locations
  for(int i=0; i < jobSize; i++)
  {     
    uint2  validBase = values[offSet + tid];
    int value        = (validBase.x >> 31);  
    value           += (validBase.y >> 31);             

    int idx  = hillisSteele5(&shmemCM[48*threadIdx.y+16], curCount, value, threadIdx.x);

    if((validBase.x >> 31))
    {
      output[idx + outputOffset] = validBase.x & 0x7FFFFFFF;
      idx++;
    }
    if((validBase.y >> 31))
    {
      output[idx + outputOffset] = validBase.y & 0x7FFFFFFF;
    }

    outputOffset += curCount;       //Increase the output offset
    offSet       += blockDim.x;     //Step to the next N threads
  }
  

  //Block 0 handles any extra elements that couldn't be divided equally
  if(bid == 0)
  {
    //Here i use single element reads for ease of boundary conditions and steps
    offSet       = sParam.extraOffset;   
    outputOffset = counts[gridDim.x*blockDim.y]; 

    uint* value2 = (uint*) values;

    for(int i=0; i < sParam.extraElements;  i += blockDim.x)
    {
      uint value = 0;     
      if((offSet + i +  tid) < (N))  //Make sure we dont read more than there are items      
        value = value2[offSet + i +  tid];     
 
        int idx  = hillisSteele5(&shmemCM[48*threadIdx.y+16], curCount, value >> 31, tid);

      if((offSet + i +  tid) < N)  
        if(value >> 31)
          output[idx + outputOffset] = value & 0x7FFFFFFF; 

      outputOffset += curCount;       //Increase the output offset
    }
  }//end if bid==0 
}//end compact_move

//The kernel that actually moves the data
extern "C"  __global__ void compact_move( uint2 *values,
                             uint *output, 
                             uint *counts,  
                             const int N,
                            setupParams sParam)
{
  extern __shared__ unsigned int shmemCM[];
  compact_moveD(values, output, counts,N,sParam,shmemCM);
}


//The kernel that actually moves/splits the data
__device__ void split_moveD( uint2 *valid,
                                        uint *output, 
                                        uint *counts,  
                                        const int N,                        
                                        setupParams sParam, volatile unsigned int *shmemSM)
{
  //Walk the values of this block
  const int tid    =  threadIdx.x;
  const int bid    =  blockDim.y *  blockIdx.x + threadIdx.y;


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

  //Do per step the prefix scan to determine the output locations
  for(int i=0; i < jobSize; i++)
  {     
    uint2  validBase = valid[offSet + tid];
    int value    = (validBase.x >> 31);  
    value       += (validBase.y >> 31);             

    idx  = hillisSteele5(&shmemSM[48*threadIdx.y+16], curCount, value, tid);
    ridx = threadIdx.x*2 - idx; //lane*2 - idx , *2 since we read 2 items a time
    
    if((validBase.x >> 31))
    {
      output[idx + outputOffset] = validBase.x & 0x7FFFFFFF;
      idx++;
    }
    else
    {
      output[ridx + rightOutputOffset] = validBase.x & 0x7FFFFFFF;      
      ridx++;
    }

    if((validBase.y >> 31))
      output[idx + outputOffset] = validBase.y & 0x7FFFFFFF;
    else
      output[ridx + rightOutputOffset] = validBase.y & 0x7FFFFFFF;

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
    
    for(int i=0; i < sParam.extraElements;  i += blockDim.x)
    {
      uint value = 0;    
      if((offSet + i +  tid) < (N))  //Make sure we dont read more than there are items      
        value = valid2[offSet + i +  tid];     

      idx  = hillisSteele5(&shmemSM[48*threadIdx.y+16], curCount, value >> 31, tid);
      ridx = threadIdx.x - idx;

      if((offSet + i +  tid) < N)  
        if(value >> 31)
          output[idx + outputOffset]       = value & 0x7FFFFFFF; 
        else
          output[ridx + rightOutputOffset] = value & 0x7FFFFFFF;

      outputOffset      += curCount;       //Increase the output offset
      rightOutputOffset += 32-curCount;    //32 since we do only 1 at a time
    }
  }//end if bid==0 
}//end split_move


//The kernel that actually moves/splits the data
extern "C"  __global__ void split_move( uint2 *valid,
                                        uint *output, 
                                        uint *counts,  
                                        const int N,                        
                                        setupParams sParam)
{
  extern __shared__ int unsigned shmemSM[];
  split_moveD(valid, output, counts, N, sParam, shmemSM);

}

