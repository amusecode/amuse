#include "octree.h"


//Compacts an array of integers, the values in srcValid indicate if a
//value is valid (1 == valid anything else is UNvalid) returns the 
//compacted values in the output array and the total 
//number of valid items is stored in 'count' 
void octree::gpuCompact(my_dev::context &devContext, 
                        my_dev::dev_mem<uint> &srcValues,
                        my_dev::dev_mem<uint> &output,                        
                        int N, int *validCount)
{

  // In the next step we associate the GPU memory with the Kernel arguments
  
//   my_dev::dev_mem<uint> counts(devContext, 512), countx(devContext, 512);
  //Memory that should be alloced outside the function:
  //devMemCounts and devMemCountsx 

  //Kernel configuration parameters
  setupParams sParam;
  sParam.jobs = (N / 64) / 480  ; //64=32*2 2 items per look, 480 is 120*4, number of procs
  sParam.blocksWithExtraJobs = (N / 64) % 480; 
  sParam.extraElements = N % 64;
  sParam.extraOffset = N - sParam.extraElements;

  compactCount.set_arg<cl_mem>(0, srcValues.p());
  compactCount.set_arg<cl_mem>(1, this->devMemCounts.p());
  compactCount.set_arg<uint>(2, &N);
  compactCount.set_arg<int>(3, NULL, 128);
  compactCount.set_arg<setupParams>(4, &sParam);

  vector<size_t> localWork(2), globalWork(2);
  globalWork[0] = 32*120;   globalWork[1] = 4;
  localWork [0] = 32;       localWork[1] = 4;   
  compactCount.setWork(globalWork, localWork);

  ///////////////

  exScanBlock.set_arg<cl_mem>(0, this->devMemCounts.p());  
  int blocks = 120*4;
  exScanBlock.set_arg<int>(1, &blocks);
  exScanBlock.set_arg<cl_mem>(2, this->devMemCountsx.p());
  exScanBlock.set_arg<int>(3, NULL, 512); //shared memory allocation

  globalWork[0] = 512; globalWork[1] = 1;
  localWork [0] = 512; localWork [1] = 1;

  exScanBlock.setWork(globalWork, localWork);

  //////////////

  compactMove.set_arg<cl_mem>(0, srcValues.p());
  compactMove.set_arg<cl_mem>(1, output.p());
  compactMove.set_arg<cl_mem>(2, this->devMemCounts.p());
  compactMove.set_arg<uint>(3, &N);
  compactMove.set_arg<uint>(4, NULL, 192); //Dynamic shared memory
  compactMove.set_arg<setupParams>(5, &sParam);

  globalWork[0] = 120*32;  globalWork[1] = 4;
  localWork [0] = 32;      localWork [1] = 4;

  compactMove.setWork(globalWork, localWork);

  ////////////////////

  compactCount.execute();
  
  exScanBlock.execute();
  
//   
  compactMove.execute();
  
  #ifdef USE_CUDA
    cuCtxSynchronize();
  #else
    clFinish(devContext.get_command_queue());
  #endif
  this->devMemCountsx.d2h();
  *validCount = this->devMemCountsx[0];
  //printf("Total number of valid items: %d \n", countx[0]);
 

}

//Splits an array of integers, the values in srcValid indicate if a
//value is valid (1 == valid anything else is UNvalid) returns the 
//splitted values in the output array (first all valid 
//number and then the invalid ones) and the total
//number of valid items is stored in 'count' 
void octree::gpuSplit(my_dev::context &devContext, 
                        my_dev::dev_mem<uint> &srcValues,
                        my_dev::dev_mem<uint> &output,                        
                        int N, int *validCount)
{

  // In the next step we associate the GPU memory with the Kernel arguments
//   my_dev::dev_mem<uint> counts(devContext, 512), countx(devContext, 512);
  //Memory that should be alloced outside the function:
  //devMemCounts and devMemCountsx 
  

  //Kernel configuration parameters
  setupParams sParam;
  sParam.jobs = (N / 64) / 480  ; //64=32*2 2 items per look, 480 is 120*4, number of procs
  sParam.blocksWithExtraJobs = (N / 64) % 480; 
  sParam.extraElements = N % 64;
  sParam.extraOffset = N - sParam.extraElements;
  
  
//   printf("Param info: %d %d %d %d \n", sParam.jobs, sParam.blocksWithExtraJobs, sParam.extraElements, sParam.extraOffset);

  compactCount.set_arg<cl_mem>(0, srcValues.p());
  compactCount.set_arg<cl_mem>(1, this->devMemCounts.p());
  compactCount.set_arg<uint>(2, &N);
  compactCount.set_arg<int>(3, NULL, 128);
  compactCount.set_arg<setupParams>(4, &sParam);

  vector<size_t> localWork(2), globalWork(2);
  globalWork[0] = 32*120;   globalWork[1] = 4;
  localWork [0] = 32;       localWork[1] = 4;   
  compactCount.setWork(globalWork, localWork);

  ///////////////

  exScanBlock.set_arg<cl_mem>(0, this->devMemCounts.p());  
  int blocks = 120*4;
  exScanBlock.set_arg<int>(1, &blocks);
  exScanBlock.set_arg<cl_mem>(2, this->devMemCountsx.p());
  exScanBlock.set_arg<int>(3, NULL, 512); //shared memory allocation

  globalWork[0] = 512; globalWork[1] = 1;
  localWork [0] = 512; localWork [1] = 1;

  exScanBlock.setWork(globalWork, localWork);

  //////////////

  splitMove.set_arg<cl_mem>(0, srcValues.p());
  splitMove.set_arg<cl_mem>(1, output.p());
  splitMove.set_arg<cl_mem>(2, this->devMemCounts.p());
  splitMove.set_arg<uint>(3, &N);
  splitMove.set_arg<uint>(4, NULL, 192); //Dynamic shared memory
  splitMove.set_arg<setupParams>(5, &sParam);

  globalWork[0] = 120*32;  globalWork[1] = 4;
  localWork [0] = 32;      localWork [1] = 4;

  splitMove.setWork(globalWork, localWork);

  ////////////////////
  compactCount.execute();

//   exit(0);
//   counts.d2h();
//   for(int i=0; i < 482; i++)
//   {
//     printf("%d\t%d\n", i, counts[i]);
//   }
//   

  exScanBlock.execute();
  
  splitMove.execute();

  //TODO fix the damn clFinish function
  #ifdef USE_CUDA
    cuCtxSynchronize();
  #else
    clFinish(devContext.get_command_queue());
  #endif
  this->devMemCountsx.d2h();
  *validCount = this->devMemCountsx[0];
  //printf("Total number of valid items: %d \n", countx[0]);
}





/*
Sort an array of int4, the idea is that the key is somehow moved into x/y/z and the
value is put in w... 
Sorts values based on the last item so order becomes something like:
z y x
2 2 1
2 1 2
2 3 3
2 5 3

*/

// If srcValues and buffer are different, then the original values
// are preserved, if they are the same srcValues will be overwritten
void  octree::gpuSort(my_dev::context &devContext,
                      my_dev::dev_mem<uint4> &srcValues,
                      my_dev::dev_mem<uint4> &output,
                      my_dev::dev_mem<uint4> &buffer,
                      int N, int numberOfBits, int subItems,
                      tree_structure &tree) {

  //Extra buffer values

//   my_dev::dev_mem<uint> simpleKeys(devContext, N);    //Int keys,
//   my_dev::dev_mem<uint> permutation(devContext, N);   //Permutation values, for sorting the int4 data
//   my_dev::dev_mem<int> output32b(devContext, N); //Permutation values, for sorting the int4 data
//   my_dev::dev_mem<uint> valuesOutput(devContext, N);  //Buffers for the values which are the indexes
  
  
  my_dev::dev_mem<uint> simpleKeys(devContext);    //Int keys,
  my_dev::dev_mem<uint> permutation(devContext);   //Permutation values, for sorting the int4 data
  my_dev::dev_mem<int>  output32b(devContext);       //Permutation values, for sorting the int4 data
  my_dev::dev_mem<uint> valuesOutput(devContext);  //Buffers for the values which are the indexes
  
  int prevOffsetSum = getAllignmentOffset(4*N); //The offset of output

  
  simpleKeys.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                          tree.generalBuffer1.get_flags(), 
                          tree.generalBuffer1.get_devMem(),
                          &tree.generalBuffer1[8*N], 8*N,
                          N, prevOffsetSum + getAllignmentOffset(8*N + prevOffsetSum));    //Ofset 8 since we have 2 uint4 before
  
  prevOffsetSum += getAllignmentOffset(8*N + prevOffsetSum);
  
  permutation.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                          tree.generalBuffer1.get_flags(), 
                          tree.generalBuffer1.get_devMem(),
                          &tree.generalBuffer1[9*N], 9*N,
                          N, prevOffsetSum + getAllignmentOffset(9*N + prevOffsetSum));  //N elements after simpleKeys    

  prevOffsetSum += getAllignmentOffset(9*N + prevOffsetSum);
  

  output32b.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                          tree.generalBuffer1.get_flags(), 
                          tree.generalBuffer1.get_devMem(),
                          &tree.generalBuffer1[10*N], 10*N,
                          N, prevOffsetSum + getAllignmentOffset(10*N + prevOffsetSum));  //N elements after permutation      
  
  prevOffsetSum += getAllignmentOffset(10*N + prevOffsetSum);
  
  valuesOutput.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                          tree.generalBuffer1.get_flags(), 
                          tree.generalBuffer1.get_devMem(),
                          &tree.generalBuffer1[11*N], 11*N,
                          N, prevOffsetSum + getAllignmentOffset(11*N + prevOffsetSum));  //N elements after output32b        

    
  //Dimensions for the kernels that shuffle and extract data
  const int blockSize = 256;
  int ng = (N)/blockSize + 1;
  int nx = (int)sqrt(ng);
  int ny = (ng-1)/nx + 1;

  vector<size_t> localWork(2), globalWork(2);
  globalWork[0] = nx*blockSize;   globalWork[1] = ny;
  localWork [0] = blockSize;       localWork[1] = 1;

  extractInt.setWork(globalWork, localWork);
  fillSequence.setWork(globalWork, localWork);
  reOrderKeysValues.setWork(globalWork, localWork);
  
  //Idx depends on subitems, z goes first, x last if subitems = 3
  //subitems = 3, than idx=2
  //subitems = 2, than idx=1
  //subitems = 1, than idx=0
  //intIdx = subItems-1   
  int intIdx = subItems-1;

  extractInt.set_arg<cl_mem>(0, srcValues.p());
  extractInt.set_arg<cl_mem>(1, simpleKeys.p());
  extractInt.set_arg<uint>(2, &N);
  extractInt.set_arg<int>(3, &intIdx);//bit idx

  fillSequence.set_arg<cl_mem>(0, permutation.p());
  fillSequence.set_arg<uint>(1, &N);

  reOrderKeysValues.set_arg<cl_mem>(0, srcValues.p());
  reOrderKeysValues.set_arg<cl_mem>(1, output.p());
  reOrderKeysValues.set_arg<cl_mem>(2, valuesOutput.p());
  reOrderKeysValues.set_arg<uint>(3, &N);

  extractInt.execute();
  fillSequence.execute();

  //Now sort the first 32bit keys
  //Using 32bit sort with key and value seperated    
  gpuSort_32b(devContext, 
                   simpleKeys, permutation,
//                     output32b, aPing32b,
                   output32b, simpleKeys,
//                    valuesOutput,valuesAPing,
                   valuesOutput,permutation,
//                   count,
                   N, 32);


  //Now reorder the main keys
  //Use output as the new output/src value thing buffer
  reOrderKeysValues.execute();
  
  if(subItems == 1)
  {
    //Only doing one 32bit sort. Data is already in output so done
    return;
  }


  //2nd set of 32bit keys
  //Idx depends on subitems, z goes first, x last if subitems = 3  
  //subitems = 3, than idx=1
  //subitems = 2, than idx=0
  //subitems = 1, completed previous round
  //intIdx = subItems-2   
  intIdx = subItems-2;
  
  extractInt.set_arg<cl_mem>(0, output.p());
  extractInt.set_arg<int>(3, &intIdx);//smem size

  reOrderKeysValues.set_arg<cl_mem>(0, output.p());
  reOrderKeysValues.set_arg<cl_mem>(1, buffer.p());
 
  extractInt.execute();
  
  fillSequence.execute();

  //Now sort the 2nd 32bit keys
  //Using 32bit sort with key and value seperated    
  gpuSort_32b(devContext, 
                   simpleKeys, permutation,
                   output32b, simpleKeys,
//                    output32b, aPing32b,
//                   valuesOutput,valuesAPing,
                   valuesOutput,permutation,
                   //count,
                   N, 32);
                   
  reOrderKeysValues.execute();
  

  if(subItems == 2)
  {
    //Doing two 32bit sorts. Data is in buffer
    //so move the data from buffer to output    
    output.copy(buffer, buffer.get_size());    
    return;
  }

  //3th set of 32bit keys
  //Idx depends on subitems, z goes first, x last if subitems = 3  
  //subitems = 3, than idx=0
  //subitems = 2, completed previous round
  //subitems = 1, completed previous round
  //intIdx = subItems-2     
  intIdx = 0;
  
  extractInt.set_arg<cl_mem>(0, buffer.p());
  extractInt.set_arg<int>(3, &intIdx);//integer idx

  reOrderKeysValues.set_arg<cl_mem>(0, buffer.p());
  reOrderKeysValues.set_arg<cl_mem>(1, output.p());

  extractInt.execute();
  fillSequence.execute();
  //Now sort the 32bit keys
  //Using int2 with key and value combined
  //See sortArray4
  //Using key and value in a seperate array
  //Now sort the 2nd 32bit keys
  //Using 32bit sort with key and value seperated    
  gpuSort_32b(devContext, 
              simpleKeys, permutation,
              output32b, simpleKeys,
//               output32b, aPing32b,
//               valuesOutput,valuesAPing,
              valuesOutput,permutation,
              //count,
              N, 32);  

  reOrderKeysValues.execute();

  clFinish(devContext.get_command_queue());

//   fprintf(stderr, "sortArray2 done in %g sec (Without memory alloc & compilation) \n", get_time() - t0);
}


void octree::gpuSort_32b(my_dev::context &devContext, 
                    my_dev::dev_mem<uint> &srcKeys,     my_dev::dev_mem<uint> &srcValues,
                    my_dev::dev_mem<int>  &keysOutput,  my_dev::dev_mem<uint> &keysAPing,
                    my_dev::dev_mem<uint> &valuesOutput,my_dev::dev_mem<uint> &valuesAPing,
                    int N, int numberOfBits)
{

  int bitIdx = 0;

  //Step 1, do the count
  //Memory that should be alloced outside the function:

  setupParams sParam;
  sParam.jobs = (N / 64) / 480  ; //64=32*2 2 items per look, 480 is 120*4, number of procs
  sParam.blocksWithExtraJobs = (N / 64) % 480;
  sParam.extraElements = N % 64;
  sParam.extraOffset = N - sParam.extraElements;

  sortCount.set_arg<cl_mem>(0, srcKeys.p());
  sortCount.set_arg<cl_mem>(1, this->devMemCounts.p());
  sortCount.set_arg<uint>(2, &N);
  sortCount.set_arg<int>(3, NULL, 128);//smem size
  sortCount.set_arg<setupParams>(4, &sParam);
  sortCount.set_arg<int>(5, &bitIdx);


  vector<size_t> localWork(2), globalWork(2);
  globalWork[0] = 32*120;   globalWork[1] = 4;
  localWork [0] = 32;       localWork[1] = 4;
  sortCount.setWork(globalWork, localWork);

  ///////////////

  exScanBlock.set_arg<cl_mem>(0, this->devMemCounts.p());
  int blocks = 120*4;
  exScanBlock.set_arg<int>(1, &blocks);
  exScanBlock.set_arg<cl_mem>(2, this->devMemCountsx.p());
  exScanBlock.set_arg<int>(3, NULL, 512); //shared memory allocation

  globalWork[0] = 512; globalWork[1] = 1;
  localWork [0] = 512; localWork [1] = 1;

  exScanBlock.setWork(globalWork, localWork);

  //////////////

  sortMove.set_arg<cl_mem>(0, srcKeys.p());
  sortMove.set_arg<cl_mem>(1, keysOutput.p());
  sortMove.set_arg<cl_mem>(2, srcValues.p());
  sortMove.set_arg<cl_mem>(3, valuesOutput.p());
  sortMove.set_arg<cl_mem>(4, this->devMemCounts.p());
  sortMove.set_arg<uint>(5, &N);
  sortMove.set_arg<uint>(6, NULL, 192); //Dynamic shared memory 128+64 , prefux sum buffer
  sortMove.set_arg<uint>(7, NULL, 64*4); //Dynamic shared memory stage buffer
  sortMove.set_arg<uint>(8, NULL, 64*4); //Dynamic shared memory stage_values buffer
  sortMove.set_arg<setupParams>(9, &sParam);
  sortMove.set_arg<int>(10, &bitIdx);

  globalWork[0] = 120*32;  globalWork[1] = 4;
  localWork [0] = 32;      localWork [1] = 4;

  sortMove.setWork(globalWork, localWork);

  bool pingPong = false;

  //Execute bitIdx 0


  sortCount.execute();
  
  exScanBlock.execute();
     
  sortMove.execute();  

  //Swap buffers
  sortCount.set_arg<cl_mem>(0, keysOutput.p());
  sortMove.set_arg<cl_mem>(0, keysOutput.p());
  sortMove.set_arg<cl_mem>(1, keysAPing.p());
  sortMove.set_arg<cl_mem>(2, valuesOutput.p());
  sortMove.set_arg<cl_mem>(3, valuesAPing.p());

  //Remaining bits, ping ponging buffers
  for(int i=1; i < numberOfBits; i++)
  {
    bitIdx = i;
    sortCount.set_arg<int>(5, &bitIdx);
    sortMove.set_arg<int>(10, &bitIdx);

    sortCount.execute();
    exScanBlock.execute(); 
    
    sortMove.execute();

    //Switch buffers
    if(pingPong)
    {
      sortCount.set_arg<cl_mem>(0, keysOutput.p());

      sortMove.set_arg<cl_mem>(0, keysOutput.p());
      sortMove.set_arg<cl_mem>(1, keysAPing.p());

      sortMove.set_arg<cl_mem>(2, valuesOutput.p());
      sortMove.set_arg<cl_mem>(3, valuesAPing.p());

      pingPong = false;
    }
    else
    {
      sortCount.set_arg<cl_mem>(0, keysAPing.p());

      sortMove.set_arg<cl_mem>(0, keysAPing.p());
      sortMove.set_arg<cl_mem>(1, keysOutput.p());

      sortMove.set_arg<cl_mem>(2, valuesAPing.p());
      sortMove.set_arg<cl_mem>(3, valuesOutput.p());

      pingPong = true;
    }
  }
  

  #ifdef USE_CUDA
    cuCtxSynchronize();
  #else
    clFinish(devContext.get_command_queue());
  #endif
  
//   fprintf(stderr, "sortArray2_32b done in %g sec\n", get_time() - t0);
}

