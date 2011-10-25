#include "octree.h"
void octree::sort_keys() {

  my_ocl::kernel sortCountIntr4(oclContext);
  my_ocl::kernel exScanBlock(oclContext);
  my_ocl::kernel sortMoveInt4(oclContext);

  sortCountInt4.load_source("sortKernels.cl", "");
  sortCountInt4.create("sort_count_int4");
  
  exScanBlock.load_source("sortKernels.cl", "");
  exScanBlock.create("exclusive_scan_block");
  
  sortMoveInt4.load_source("sortKernels.cl", "");
  sortMoveInt4.create("sort_move_stage4");

  // In the next step we associate the GPU memory with the Kernel arguments

  my_ocl::ocl_mem<uint> counts(oclContext, 512), countx(oclContext, 512);
  
  //Kernel configuration parameters
  
  setupParams sParam;
  sParam.jobs = (n_b / 32) / 480  ; //32=32*1 1 items per look, 480 is 120*4, number of procs
  sParam.blocksWithExtraJobs = (n_b / 32) % 480; 
  sParam.extraElements = n_b % 32;
  sParam.extraOffset = n_b - sParam.extraElements;

  int bitIdx = 0;
  
  sortCountInt4.set_arg<cl_mem>     (0, srcValues.p());
  sortCountInt4.set_arg<cl_mem>     (1, counts.p());
  sortCountInt4.set_arg<uint>       (2, &n_b);
  sortCountInt4.set_arg<int>        (3, NULL, 128);
  sortCountInt4.set_arg<setupParams>(4, &sParam);
  sortCountInt4.set_arg<int>        (5, &bitIdx);

  vector<size_t> localWork(2), globalWork(2);
  globalWork[0] = 32*120;   globalWork[1] = 4;
  localWork [0] = 32;       localWork [1] = 4;   
  sortCountInt4.setWork(globalWork, localWork);
  
  ///////////////

  exScanBlock.set_arg<cl_mem>(0, counts.p());  
  int blocks = 120*4;
  exScanBlock.set_arg<int>   (1, &blocks);
  exScanBlock.set_arg<cl_mem>(2, countx.p());
  exScanBlock.set_arg<int>   (3, NULL, 512); //shared memory allocation

  globalWork[0] = 512; globalWork[1] = 1;
  localWork [0] = 512; localWork [1] = 1;
  
  exScanBlock.setWork(globalWork, localWork);

  //////////////

  sortMoveInt4.set_arg<cl_mem>     (0, srcValues.p());
  sortMoveInt4.set_arg<cl_mem>     (1, output.p());
  sortMoveInt4.set_arg<cl_mem>     (2, counts.p());
  sortMoveInt4.set_arg<uint>       (3, &N);
  sortMoveInt4.set_arg<uint>       (4, NULL, 192); //Dynamic shared memory 128+64 , prefux sum buffer
  sortMoveInt4.set_arg<uint4>      (5, NULL, 32*4); //Dynamic shared memory stage buffer
  sortMoveInt4.set_arg<setupParams>(6, &sParam);
  sortMoveInt4.set_arg<int>        (7, &bitIdx);

  globalWork[0] = 120*32;  globalWork[1] = 4;
  localWork [0] = 32;      localWork [1] = 4;

  sortMoveInt4.setWork(globalWork, localWork);

  ////////////////////
  

  bool pingPong = false;

  double t0 = get_time();

  //Execute bitIdx 0
  sortCountInt4.execute();
  exScanBlock.execute();
  sortMoveInt4.execute();

  //Swap buffers
  sortCountInt4.set_arg<cl_mem>(0, output.p());
  sortMoveInt4.set_arg<cl_mem>(0, output.p());
  sortMoveInt4.set_arg<cl_mem>(1, buffer.p());


  //Remaining bits, ping ponging buffers
  for(int i=1; i < numberOfBits; i++)
  {
    bitIdx = i;
    sortCountInt4.set_arg<int>(5, &bitIdx);
    sortMoveInt4.set_arg<int>(7, &bitIdx);

    sortCountInt4.execute();
    exScanBlock.execute();
    sortMoveInt4.execute();

    //Switch buffers
    if(pingPong)
    {
      sortCountInt4.set_arg<cl_mem>(0, output.p());
      sortMoveInt4.set_arg<cl_mem>(0, output.p());
      sortMoveInt4.set_arg<cl_mem>(1, buffer.p());
      pingPong = false;
    }
    else
    {
      sortCountInt4.set_arg<cl_mem>(0, buffer.p());
      sortMoveInt4.set_arg<cl_mem>(0, buffer.p());
      sortMoveInt4.set_arg<cl_mem>(1, output.p());
      pingPong = true;
    }
  }
  clFinish(oclContext.get_command_queue());

  fprintf(stderr, "done in %g sec\n", get_time() - t0);


}
