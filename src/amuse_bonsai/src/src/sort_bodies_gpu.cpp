#include "octree.h"

void octree::getBoundaries(tree_structure &tree, real4 &r_min, real4 &r_max)
{

  //Start reduction to get the boundary's of the system
  boundaryReduction.set_arg<int>(0, &tree.n);
//   boundaryReduction.set_arg<cl_mem>(1, tree.bodies_pos.p());
  boundaryReduction.set_arg<cl_mem>(1, tree.bodies_Ppos.p());
  boundaryReduction.set_arg<cl_mem>(2, devMemRMIN.p());
  boundaryReduction.set_arg<cl_mem>(3, devMemRMAX.p());

  boundaryReduction.setWork(tree.n, NTHREAD_BOUNDARY, NBLOCK_BOUNDARY);  //256 threads and 120 blocks in total
  boundaryReduction.execute();
  
   
  devMemRMIN.d2h();     //Need to be defined and initialized somewhere outside this function
  devMemRMAX.d2h();     //Need to be defined and initialized somewhere outside this function
  r_min = (real4){+1e10, +1e10, +1e10, +1e10}; 
  r_max = (real4){-1e10, -1e10, -1e10, -1e10};   
  
  //Reduce the blocks, done on host since its
  //A faster and B we need the results anyway
  for (int i = 0; i < 120; i++) {    
    r_min.x = fmin(r_min.x, devMemRMIN[i].x);
    r_min.y = fmin(r_min.y, devMemRMIN[i].y);
    r_min.z = fmin(r_min.z, devMemRMIN[i].z);
    
    r_max.x = fmax(r_max.x, devMemRMAX[i].x);
    r_max.y = fmax(r_max.y, devMemRMAX[i].y);
    r_max.z = fmax(r_max.z, devMemRMAX[i].z);    
//     printf("%f\t%f\t%f\t || \t%f\t%f\t%f\n", rMIN[i].x,rMIN[i].y,rMIN[i].z,rMAX[i].x,rMAX[i].y,rMAX[i].z);    
  }
  
  rMinLocalTree = r_min;
  rMaxLocalTree = r_max;
  
  printf("Found boundarys, number of particles %d : \n", tree.n);
  printf("min: %f\t%f\t%f\tmax: %f\t%f\t%f \n", r_min.x,r_min.y,r_min.z,r_max.x,r_max.y,r_max.z);
}

void octree::getBoundariesGroups(tree_structure &tree, real4 &r_min, real4 &r_max)
{

  //Start reduction to get the boundary's of the system
  boundaryReductionGroups.set_arg<int>(0, &tree.n_groups);
  boundaryReductionGroups.set_arg<cl_mem>(1, tree.groupCenterInfo.p());
  boundaryReductionGroups.set_arg<cl_mem>(2, tree.groupSizeInfo.p());
  boundaryReductionGroups.set_arg<cl_mem>(3, devMemRMIN.p());
  boundaryReductionGroups.set_arg<cl_mem>(4, devMemRMAX.p());

  boundaryReductionGroups.setWork(tree.n_groups, NTHREAD_BOUNDARY, NBLOCK_BOUNDARY);  //256 threads and 120 blocks in total
  boundaryReductionGroups.execute();

   
  devMemRMIN.d2h();     //Need to be defined and initialized somewhere outside this function
  devMemRMAX.d2h();     //Need to be defined and initialized somewhere outside this function
  r_min = (real4){+1e10, +1e10, +1e10, +1e10}; 
  r_max = (real4){-1e10, -1e10, -1e10, -1e10};   
  
  //Reduce the blocks, done on host since its
  //A faster and B we need the results anyway
  for (int i = 0; i < 120; i++) {    
    r_min.x = fmin(r_min.x, devMemRMIN[i].x);
    r_min.y = fmin(r_min.y, devMemRMIN[i].y);
    r_min.z = fmin(r_min.z, devMemRMIN[i].z);
    
    r_max.x = fmax(r_max.x, devMemRMAX[i].x);
    r_max.y = fmax(r_max.y, devMemRMAX[i].y);
    r_max.z = fmax(r_max.z, devMemRMAX[i].z);    
  }
  
  printf("Found group boundarys before increase, number of groups %d : \n", tree.n_groups);
  printf("min: %f\t%f\t%f\tmax: %f\t%f\t%f \n", r_min.x,r_min.y,r_min.z,r_max.x,r_max.y,r_max.z);
  
  //Prevent small-numerical differences by making the group/box slightly bigger
  
  double smallFac1 = 0.999;
  double smallFac2 = 1.001;
  
  //Note that we have to check the sign to move the border in the right
  //direction
  r_min.x = (r_min.x < 0) ? r_min.x * smallFac2 : r_min.x * smallFac1;
  r_min.y = (r_min.y < 0) ? r_min.y * smallFac2 : r_min.y * smallFac1;
  r_min.z = (r_min.z < 0) ? r_min.z * smallFac2 : r_min.z * smallFac1;

  r_max.x = (r_max.x < 0) ? r_max.x * smallFac1 : r_max.x * smallFac2;
  r_max.y = (r_max.y < 0) ? r_max.y * smallFac1 : r_max.y * smallFac2;
  r_max.z = (r_max.z < 0) ? r_max.z * smallFac1 : r_max.z * smallFac2;
  

  rMinLocalTreeGroups = r_min;
  rMaxLocalTreeGroups = r_max;
  
  
  printf("Found group boundarys after increase, number of groups %d : \n", tree.n_groups);
  printf("min: %f\t%f\t%f\tmax: %f\t%f\t%f \n", r_min.x,r_min.y,r_min.z,r_max.x,r_max.y,r_max.z);
}



void octree::sort_bodies(tree_structure &tree, bool doDomainUpdate) {

  //We assume the bodies are already onthe GPU
  devContext.startTiming();
  real4 r_min = {+1e10, +1e10, +1e10, +1e10}; 
  real4 r_max = {-1e10, -1e10, -1e10, -1e10};   
  
  if(doDomainUpdate)
  {
    getBoundaries(tree, r_min, r_max);  
    //Sync the boundary over the various processes
    if(this->mpiGetNProcs() > 1)
    {
      this->sendCurrentRadiusInfo(r_min, r_max);
    }
    rMinGlobal = r_min;    rMaxGlobal = r_max;
  }
  
  r_min = rMinGlobal;
  r_max = rMaxGlobal;
  
  //Compute the boundarys of the tree  
  real size     = 1.001*fmax(r_max.z - r_min.z,
                        fmax(r_max.y - r_min.y, r_max.x - r_min.x));
  
  tree.corner   = (real4){0.5*(r_min.x + r_max.x) - 0.5*size,
                          0.5*(r_min.y + r_max.y) - 0.5*size,
                          0.5*(r_min.z + r_max.z) - 0.5*size, size}; 
       
  tree.domain_fac   = size/(1 << MAXLEVELS);
  
  
  float idomain_fac = 1.0/tree.domain_fac;
  float domain_fac  = tree.domain_fac;
  
  tree.corner.w = domain_fac;  
  
  printf("Corner: %f %f %f idomain fac: %f domain_fac: %f\n", 
         tree.corner.x, tree.corner.y, tree.corner.z, idomain_fac, domain_fac);
  printf("domain fac: %f idomain_fac: %f size: %f MAXLEVELS: %d \n", domain_fac, idomain_fac, size, MAXLEVELS);
  //printf("Test theta: %f  \n", inv_theta);
  
  //Compute the keys
  build_key_list.set_arg<cl_mem>(0,   tree.bodies_key.p());
  build_key_list.set_arg<cl_mem>(1,   tree.bodies_Ppos.p());
  build_key_list.set_arg<int>(2,      &tree.n);
  build_key_list.set_arg<real4>(3,    &tree.corner);
  
  build_key_list.setWork(tree.n, 128); //128 threads per block
  build_key_list.execute();  
 
 
  //Call the GPUSort function, since we made it general 
  //into a uint4 so we can extend the tree to 96bit key
  //we have to convert to 64bit key to a 96bit for sorting
  //and back from 96 to 64    
  my_dev::dev_mem<uint4>  srcValues(devContext);
  my_dev::dev_mem<uint4>  output(devContext);
  
  //The generalBuffer1 has size uint*4*N*3
  //this buffer gets part: 0-uint*4*N
  srcValues.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                         tree.generalBuffer1.get_flags(), 
                         tree.generalBuffer1.get_devMem(),
                         &tree.generalBuffer1[0], 0,  
                         tree.n, getAllignmentOffset(0));  
 
  //this buffer gets part: uint*4*N-uint*4*N*2
  output.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                         tree.generalBuffer1.get_flags(), 
                         tree.generalBuffer1.get_devMem(),
                         &tree.generalBuffer1[4*tree.n], 4*tree.n,
                         tree.n, getAllignmentOffset(4*tree.n));
   //Extract the keys
  convertKey64to96.set_arg<cl_mem>(0,   tree.bodies_key.p());
  convertKey64to96.set_arg<cl_mem>(1,   srcValues.p());
  convertKey64to96.set_arg<int>(2,      &tree.n);  
  convertKey64to96.setWork(tree.n, 256);
  convertKey64to96.execute();
 
  //Sort the keys  
  
  // If srcValues and buffer are different, then the original values
  // are preserved, if they are the same srcValues will be overwritten  
  gpuSort(devContext, srcValues, output, srcValues, tree.n, 32, 3, tree);
 
  my_dev::dev_mem<uint>   sortPermutation(devContext);
  sortPermutation.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                         tree.generalBuffer1.get_flags(), 
                         tree.generalBuffer1.get_devMem(),
                         &tree.generalBuffer1[0], 0,  
                         tree.n, getAllignmentOffset(0));    
  

  //Extract the keys and get the permuation required to sort the other
  //properties of the particles
  //Extract the sorted keys
  extractKeyAndPerm.set_arg<cl_mem>(0,   output.p());
  extractKeyAndPerm.set_arg<cl_mem>(1,   tree.bodies_key.p());
  extractKeyAndPerm.set_arg<cl_mem>(2,   sortPermutation.p());  
  extractKeyAndPerm.set_arg<int>(3,      &tree.n);
  extractKeyAndPerm.setWork(tree.n, 256);
  extractKeyAndPerm.execute();  

  //Use calloc here so valgrind does not complain about uninitialized values
  my_dev::dev_mem<real4>  real4Buffer(devContext);

  real4Buffer.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                         tree.generalBuffer1.get_flags(), 
                         tree.generalBuffer1.get_devMem(),
                         &tree.generalBuffer1[4*tree.n], 4*tree.n, 
                         tree.n, getAllignmentOffset(4*tree.n));      
  
  devContext.stopTiming("Sorting", 0);  

  //Call the reorder data function

  //For the position
  dataReorderR4.set_arg<int>(0,      &tree.n);
  dataReorderR4.set_arg<cl_mem>(1,   tree.bodies_pos.p());
  dataReorderR4.set_arg<cl_mem>(2,   real4Buffer.p());
  dataReorderR4.set_arg<cl_mem>(3,   sortPermutation.p()); 
  
  dataReorderR4.setWork(tree.n, 256);  
  devContext.startTiming();

  dataReorderR4.execute();  
  //Copy the data back to the source and sync
  tree.bodies_pos.copy(real4Buffer, real4Buffer.get_size()); 
  clFinish(devContext.get_command_queue());

  if(tree.needToReorder)
  {
    //Fill the oriParticleOrder with a continues increasing sequence

    //Dimensions for the kernels that shuffle and extract data
    fillSequence.set_arg<cl_mem>(0, tree.oriParticleOrder.p());
    fillSequence.set_arg<uint>(1, &tree.n);
    fillSequence.setWork(tree.n, 256);
    fillSequence.execute();
    tree.needToReorder = false;
  }

  //Velocity
  dataReorderR4.set_arg<cl_mem>(1,   tree.bodies_vel.p());
  dataReorderR4.execute();  
  tree.bodies_vel.copy(real4Buffer, real4Buffer.get_size()); 
  clFinish(devContext.get_command_queue());  
     
  //Acceleration 0
  dataReorderR4.set_arg<cl_mem>(1,   tree.bodies_acc0.p());
  dataReorderR4.execute();  
  tree.bodies_acc0.copy(real4Buffer, real4Buffer.get_size()); 
  clFinish(devContext.get_command_queue());    

  //Acceleration 1
  dataReorderR4.set_arg<cl_mem>(1,   tree.bodies_acc1.p());
  dataReorderR4.execute();  
  tree.bodies_acc1.copy(real4Buffer, real4Buffer.get_size()); 
  clFinish(devContext.get_command_queue()); 
  
  //Predicted positions 
  dataReorderR4.set_arg<cl_mem>(1,   tree.bodies_Ppos.p());
  dataReorderR4.execute();  
  tree.bodies_Ppos.copy(real4Buffer, real4Buffer.get_size()); 
  clFinish(devContext.get_command_queue()); 
  
  //Predicted velocities
  dataReorderR4.set_arg<cl_mem>(1,   tree.bodies_Pvel.p());
  dataReorderR4.execute();  
  tree.bodies_Pvel.copy(real4Buffer, real4Buffer.get_size()); 
  clFinish(devContext.get_command_queue());     
  
  my_dev::dev_mem<double2> float2Buffer(devContext);
  float2Buffer.cmalloc_copy(tree.generalBuffer1.get_pinned(),   
                         tree.generalBuffer1.get_flags(), 
                         tree.generalBuffer1.get_devMem(),
                         &tree.generalBuffer1[4*tree.n], 4*tree.n,
                         tree.n, getAllignmentOffset(4*tree.n));   

  
  //Integration time
  dataReorderF2.set_arg<int>(0,      &tree.n);
  dataReorderF2.set_arg<cl_mem>(1,   tree.bodies_time.p());
  dataReorderF2.set_arg<cl_mem>(2,   float2Buffer.p());
  dataReorderF2.set_arg<cl_mem>(3,   sortPermutation.p());   
  dataReorderF2.setWork(tree.n, 256);  
  dataReorderF2.execute();  
  tree.bodies_time.copy(float2Buffer, float2Buffer.get_size()); 

  my_dev::dev_mem<int>  intBuffer(devContext);
  intBuffer.cmalloc_copy(tree.generalBuffer1.get_pinned(),   
                         tree.generalBuffer1.get_flags(), 
                         tree.generalBuffer1.get_devMem(),
                         &tree.generalBuffer1[4*tree.n], 4*tree.n,
                         tree.n, getAllignmentOffset(4*tree.n));  

  //Particle ids
  dataReorderI1.set_arg<int>(0,      &tree.n);
  dataReorderI1.set_arg<cl_mem>(1,   tree.bodies_ids.p());
  dataReorderI1.set_arg<cl_mem>(2,   intBuffer.p());
  dataReorderI1.set_arg<cl_mem>(3,   sortPermutation.p());   
  dataReorderI1.setWork(tree.n, 256);  
  dataReorderI1.execute();  
  tree.bodies_ids.copy(intBuffer, intBuffer.get_size());   

  //Reorder the particle order indixes
  dataReorderI1.set_arg<int>(0,      &tree.n);
  dataReorderI1.set_arg<cl_mem>(1,   tree.oriParticleOrder.p());
  dataReorderI1.set_arg<cl_mem>(2,   intBuffer.p());
  dataReorderI1.set_arg<cl_mem>(3,   sortPermutation.p());   
  dataReorderI1.setWork(tree.n, 256);  
  dataReorderI1.execute();  
  tree.oriParticleOrder.copy(intBuffer, intBuffer.get_size());   

  devContext.stopTiming("Data-reordering", 1);    
}

void octree::desort_bodies(tree_structure &tree) {

  my_dev::dev_mem<uint>   sortPermutation(devContext);
  sortPermutation.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                         tree.generalBuffer1.get_flags(), 
                         tree.generalBuffer1.get_devMem(),
                         &tree.generalBuffer1[0], 0,  
                         tree.n, getAllignmentOffset(0));    

  tree.oriParticleOrder.d2h();
  //Copy the original order to the permutation so next sort will
  //put everything in the original location
  for(int i=0; i < tree.n; i++){
    sortPermutation[tree.oriParticleOrder[i]] = i;
  }    

  sortPermutation.h2d();
  
  
  //Use calloc here so valgrind does not complain about uninitialized values
  my_dev::dev_mem<real4>  real4Buffer(devContext);

  real4Buffer.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                         tree.generalBuffer1.get_flags(), 
                         tree.generalBuffer1.get_devMem(),
                         &tree.generalBuffer1[4*tree.n], 4*tree.n, 
                         tree.n, getAllignmentOffset(4*tree.n));      
  
  //Call the reorder data function


  //For the position
  dataReorderR4.set_arg<int>(0,      &tree.n);
  dataReorderR4.set_arg<cl_mem>(1,   tree.bodies_pos.p());
  dataReorderR4.set_arg<cl_mem>(2,   real4Buffer.p());
  dataReorderR4.set_arg<cl_mem>(3,   sortPermutation.p()); 
  
  dataReorderR4.setWork(tree.n, 256);  
  devContext.startTiming();

  dataReorderR4.execute();  
  //Copy the data back to the source and sync
  tree.bodies_pos.copy(real4Buffer, real4Buffer.get_size()); 
  clFinish(devContext.get_command_queue());
  
  //Velocity
  dataReorderR4.set_arg<cl_mem>(1,   tree.bodies_vel.p());
  dataReorderR4.execute();  
  tree.bodies_vel.copy(real4Buffer, real4Buffer.get_size()); 
  clFinish(devContext.get_command_queue());  
     
  //Acceleration 0
  dataReorderR4.set_arg<cl_mem>(1,   tree.bodies_acc0.p());
  dataReorderR4.execute();  
  tree.bodies_acc0.copy(real4Buffer, real4Buffer.get_size()); 
  clFinish(devContext.get_command_queue());    

  //Acceleration 1
  dataReorderR4.set_arg<cl_mem>(1,   tree.bodies_acc1.p());
  dataReorderR4.execute();  
  tree.bodies_acc1.copy(real4Buffer, real4Buffer.get_size()); 
  clFinish(devContext.get_command_queue()); 
  
  //Predicted positions 
  dataReorderR4.set_arg<cl_mem>(1,   tree.bodies_Ppos.p());
  dataReorderR4.execute();  
  tree.bodies_Ppos.copy(real4Buffer, real4Buffer.get_size()); 
  clFinish(devContext.get_command_queue()); 
  
  //Predicted velocities
  dataReorderR4.set_arg<cl_mem>(1,   tree.bodies_Pvel.p());
  dataReorderR4.execute();  
  tree.bodies_Pvel.copy(real4Buffer, real4Buffer.get_size()); 
  clFinish(devContext.get_command_queue());     
  
  my_dev::dev_mem<double2> float2Buffer(devContext);
  float2Buffer.cmalloc_copy(tree.generalBuffer1.get_pinned(),   
                         tree.generalBuffer1.get_flags(), 
                         tree.generalBuffer1.get_devMem(),
                         &tree.generalBuffer1[4*tree.n], 4*tree.n,
                         tree.n, getAllignmentOffset(4*tree.n));   
  
  //Integration time
  dataReorderF2.set_arg<int>(0,      &tree.n);
  dataReorderF2.set_arg<cl_mem>(1,   tree.bodies_time.p());
  dataReorderF2.set_arg<cl_mem>(2,   float2Buffer.p());
  dataReorderF2.set_arg<cl_mem>(3,   sortPermutation.p());   
  dataReorderF2.setWork(tree.n, 256);  
  dataReorderF2.execute();  
  tree.bodies_time.copy(float2Buffer, float2Buffer.get_size()); 

  my_dev::dev_mem<int>  intBuffer(devContext);
  intBuffer.cmalloc_copy(tree.generalBuffer1.get_pinned(),   
                         tree.generalBuffer1.get_flags(), 
                         tree.generalBuffer1.get_devMem(),
                         &tree.generalBuffer1[4*tree.n], 4*tree.n,
                         tree.n, getAllignmentOffset(4*tree.n));  

  //Particle ids
  dataReorderI1.set_arg<int>(0,      &tree.n);
  dataReorderI1.set_arg<cl_mem>(1,   tree.bodies_ids.p());
  dataReorderI1.set_arg<cl_mem>(2,   intBuffer.p());
  dataReorderI1.set_arg<cl_mem>(3,   sortPermutation.p());   
  dataReorderI1.setWork(tree.n, 256);  
  dataReorderI1.execute();  
  tree.bodies_ids.copy(intBuffer, intBuffer.get_size());   
  
  
  devContext.stopTiming("Data reordering", 1);    
}
