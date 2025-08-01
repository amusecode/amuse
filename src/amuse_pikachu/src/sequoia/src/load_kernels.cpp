#include "octree.h"

// void octree::set_context( bool disable_timing) {
//   
//   devContext.create(disable_timing);
//   set_context2();
// }
// 
// 
// void octree::set_context(std::ofstream &log, bool disable_timing) {  
//   
//   devContext.create(log, disable_timing);
//   set_context2();
// }

void octree::set_context()
{
  devContext.create(true); //Create context without timing
  createDeviceQueue();
}

void octree::set_context(std::ofstream &log, bool disable_timing)
{ 
  devContext.create(log, disable_timing);
  createDeviceQueue();  
}

void octree::createDeviceQueue()
{  
  devContext.createQueue(devID);
  devContext_flag = true;
  
  //Assign the context to the local tree
  this->localTree.setContext(devContext);

  
  //Get the number of multiprocessors and compute number of 
  //blocks to be used during the tree-walk
  nMultiProcessors   = devContext.multiProcessorCount;
  nBlocksForTreeWalk = nMultiProcessors*TREE_WALK_BLOCKS_PER_SM;
  
  allocateInternalBuffers();
}


void octree::allocateInternalBuffers()
{
  //General memory buffers
    
  //Set the context for the memory
  this->tnext.setContext(devContext);  
  this->tnext.ccalloc(NBLOCK_REDUCE,false);  
  
  this->nactive.setContext(devContext);  
  this->nactive.ccalloc(NBLOCK_REDUCE,false);    

  this->devMemRMIN.setContext(devContext);
  this->devMemRMIN.cmalloc(NBLOCK_BOUNDARY, false);  

  this->devMemRMAX.setContext(devContext);
  this->devMemRMAX.cmalloc(NBLOCK_BOUNDARY, false);  
  
  this->devMemCounts.setContext(devContext);
  this->devMemCounts.cmalloc(NBLOCK_PREFIX, false);  

  this->devMemCountsx.setContext(devContext);
  this->devMemCountsx.cmalloc(NBLOCK_PREFIX, false);  
  
  this->localTree.level_list.cmalloc(MAXLEVELS);  
  this->localTree.node_level_list.cmalloc(MAXLEVELS*2 , false);

  
  
   //General buffer is used at multiple locations and reused in different functions
  //The generalBuffer is also used during the tree-walk, so the size has to be at least
  //large enough to store the tree-walk stack. Add 4096 for extra memory allignment space
  //Times 2 since first half is used for regular walks, 2nd half for walks that go outside 
  //the memory stack and require another walk with 'unlimited' stack size
  int treeWalkStackSize = (2*LMEM_STACK_SIZE*NTHREAD*nBlocksForTreeWalk) + 4096;
   //General buffer is used at multiple locations and reused in different functions
  localTree.generalBuffer1.cmalloc(treeWalkStackSize, true);
  
  int n_bodies = 1024;  //Some standard number so we can setup the buffers
  int n_groups = 128;
  
  localTree.n_children.cmalloc(n_bodies, false);   //No transfer to host so no need to pin
  localTree.node_bodies.cmalloc(n_bodies, false);  //No transfer to host so no need to pin
    
  //Initialize buffers related to tree properties
  localTree.leafNodeIdx.cmalloc(n_bodies, false);
  localTree.multipole.cmalloc(n_bodies*3, false);
  localTree.boxSizeInfo.cmalloc(n_bodies, false);
  localTree.boxCenterInfo.cmalloc(n_bodies, false);
  
  //Initialize buffers related to groups
  localTree.body2group_list.cmalloc(n_bodies, false);
  localTree.group_list.cmalloc(n_groups,      false);
  
  localTree.activePartList.cmalloc(n_bodies+5, false);   //Active particles, TODO for later
//   localTree.ngb.cmalloc(n_bodies, false);   //Nearest neighbour, TODO for later
  localTree.interactions.cmalloc(n_bodies, false);   //Interactions, TODO for later
  
  
  
   execStream = new my_dev::dev_stream(0);    
  
}

void octree::allocateParticleSpecificBuffers(int n_bodies)
{  
  bool reduce = false;
    
  //General buffer is used at multiple locations and reused in different functions
  int treeWalkStackSize = (2*LMEM_STACK_SIZE*NTHREAD*nBlocksForTreeWalk) + 4096;
  
  int tempSize   = max(n_bodies, 4096);   //Use minium of 4096 to prevent offsets mess up with small N
  tempSize       = 3*tempSize*4 + 4096;  //Add 4096 to give some space for memory allignment  
  //TempSize is large enough to store 3x n_bodies*float4 particles, used for sorting
  tempSize = max(tempSize, treeWalkStackSize);
  
  //General buffer is used at multiple locations and reused in different functions
  localTree.generalBuffer1.cresize(tempSize, reduce);  
  
  
  //Other particle related buffers
  
 // localTree.bodies_key.cresize(n_bodies+1, reduce); // +1 since we need the last index
  
  //n_children
  //node_bodies
  
  localTree.n_children.cresize(n_bodies,  reduce);  //used in build and compute properties    
  localTree.node_bodies.cresize(n_bodies, reduce);  //used in build and compute properties 

}

void octree::allocateNodeSpecificBuffers(int n_nodes)
{
  bool reduce = false;
  
  localTree.leafNodeIdx.cresize(n_nodes, reduce);
  localTree.multipole.cresize(n_nodes*3, reduce);
  localTree.boxSizeInfo.cresize(n_nodes, reduce);
  localTree.boxCenterInfo.cresize(n_nodes, reduce);
}

void octree::allocateGroupSpecificBuffers(int n_bodies, int n_groups)
{
  bool reduce = false;
  
  localTree.body2group_list.cresize(n_bodies, reduce);
  localTree.group_list.cresize(n_groups, reduce);
  
  localTree.activePartList.cresize(n_bodies+5, reduce);   //Active particles, TODO for later, +5 since we use counters in that buffer
//   localTree.ngb.cresize(n_bodies, reduce);   //Nearest neighbour, TODO for later
  localTree.interactions.cresize(n_bodies, reduce);   //Interactions, TODO for later
  
}






void octree::load_kernels() {

  if (!devContext_flag) set_context();
  
  //If we arive here we have aquired a device, configure parts of the code

  

  std::string pathName;


  //AMUSE specific
  if(this->src_directory != NULL)
  {
    pathName.assign(this->src_directory);
  }
  else
  {  
    //Strip the executable name, to get the path name
    std::string temp(execPath);
    int idx = temp.find_last_of("/\\");
    pathName.assign(temp.substr(0, idx+1));
  }
  

  // load scan & sort kernels

  compactCount.setContext(devContext);
  exScanBlock.setContext(devContext);
  compactMove.setContext(devContext);
  splitMove.setContext(devContext);
  sortCount.setContext(devContext);
  sortMove.setContext(devContext);
  extractInt.setContext(devContext);
  fillSequence.setContext(devContext);
  reOrderKeysValues.setContext(devContext);
  convertKey64to96.setContext(devContext);
  extractKeyAndPerm.setContext(devContext);
  dataReorderR4.setContext(devContext);
  dataReorderF2.setContext(devContext);
  dataReorderI1.setContext(devContext);
  
  
#ifdef USE_CUDA
  compactCount.load_source("./scanKernels.ptx", pathName.c_str());
  compactCount.create("compact_count");
  
  exScanBlock.load_source("./scanKernels.ptx", pathName.c_str());
  exScanBlock.create("exclusive_scan_block");
  
  compactMove.load_source("./scanKernels.ptx", pathName.c_str());
  compactMove.create("compact_move");
  
  splitMove.load_source("./scanKernels.ptx", pathName.c_str());
  splitMove.create("split_move");
  
  sortCount.load_source("./sortKernels.ptx", pathName.c_str());
  sortCount.create("sort_count");

  sortMove.load_source("./sortKernels.ptx", pathName.c_str());
  sortMove.create("sort_move_stage_key_value");  

  extractInt.load_source("./sortKernels.ptx", pathName.c_str());
  extractInt.create("extractInt");  
  
  fillSequence.load_source("./sortKernels.ptx", pathName.c_str());
  fillSequence.create("fillSequence");  
  
  reOrderKeysValues.load_source("./sortKernels.ptx", pathName.c_str());
  reOrderKeysValues.create("reOrderKeysValues");    
  
  extractKeyAndPerm.load_source("./sortKernels.ptx", pathName.c_str());
  extractKeyAndPerm.create("extractKeyAndPerm");  
  
  convertKey64to96.load_source("./sortKernels.ptx", pathName.c_str());
  convertKey64to96.create("convertKey64to96");
  
  dataReorderR4.load_source("./sortKernels.ptx", pathName.c_str());
  dataReorderR4.create("dataReorderR4");  
  
  dataReorderF2.load_source("./sortKernels.ptx", pathName.c_str());
  dataReorderF2.create("dataReorderF2");  

  dataReorderI1.load_source("./sortKernels.ptx", pathName.c_str());
  dataReorderI1.create("dataReorderI1");        
  
#else
  compactCount.load_source("scanKernels.cl", "OpenCLKernels");
  compactCount.create("compact_count");
  
  exScanBlock.load_source("scanKernels.cl", "OpenCLKernels");
  exScanBlock.create("exclusive_scan_block");
  
  compactMove.load_source("scanKernels.cl", "OpenCLKernels");
  compactMove.create("compact_move");
  
  splitMove.load_source("scanKernels.cl", "OpenCLKernels");
  splitMove.create("split_move");
#endif

  // load tree-build kernels
  
  /* set context */
  
  build_key_list.setContext(devContext);
  build_valid_list.setContext(devContext);
  build_nodes.setContext(devContext);
  link_tree.setContext(devContext);
  define_groups.setContext(devContext);
  build_level_list.setContext(devContext);
  boundaryReduction.setContext(devContext);
  boundaryReductionGroups.setContext(devContext);  
  build_body2group_list.setContext(devContext);
  store_groups.setContext(devContext);
  expand_leaflist.setContext(devContext);
  
//   build_phkey_list.setContext(devContext);
  
  /* load kernels tree properties */
  
#ifdef USE_CUDA
  build_key_list.load_source("./build_tree.ptx", pathName.c_str());
  build_valid_list.load_source("./build_tree.ptx", pathName.c_str());
  build_nodes.load_source("./build_tree.ptx", pathName.c_str());
  link_tree.load_source("./build_tree.ptx", pathName.c_str());
  define_groups.load_source("./build_tree.ptx", pathName.c_str());
  build_level_list.load_source("./build_tree.ptx", pathName.c_str());
  boundaryReduction.load_source("./build_tree.ptx", pathName.c_str());
  boundaryReductionGroups.load_source("./build_tree.ptx", pathName.c_str());
  build_body2group_list.load_source("./build_tree.ptx", pathName.c_str());
  store_groups.load_source("./build_tree.ptx", pathName.c_str());
  expand_leaflist.load_source("./build_tree.ptx", pathName.c_str());
  
//   build_phkey_list.load_source("./build_tree.ptx", pathName.c_str());
  
  /* create kernels */

  build_key_list.create("cl_build_key_list");  
  build_valid_list.create("cl_build_valid_list");
  build_nodes.create("cl_build_nodes");
  link_tree.create("cl_link_tree");
  define_groups.create("build_group_list_new");
  build_level_list.create("build_level_list");
  boundaryReduction.create("boundaryReduction");
  boundaryReductionGroups.create("boundaryReductionGroups");
  build_body2group_list.create("build_body2group_list");
  store_groups.create("store_group_list");
  expand_leaflist.create("expandLeafList");
  
//   build_phkey_list.create("build_phkey_list");
  
#else
  build_key_list.load_source("build_tree.cl", "");
  build_valid_list.load_source("build_tree.cl", "");
  build_nodes.load_source("build_tree.cl", "");
  link_tree.load_source("build_tree.cl", "");
  
  /* create kernels */

  build_key_list.create("cl_build_key_list");
  build_valid_list.create("cl_build_valid_list");
  build_nodes.create("cl_build_nodes");
  link_tree.create("cl_link_tree");
#endif

  

  // load tree-props kernels

//   propsNonLeaf.setContext(devContext);
//   propsLeaf.setContext(devContext);
//   propsScaling.setContext(devContext);

  propsNonLeafD.setContext(devContext);
  propsLeafD.setContext(devContext);
  propsScalingD.setContext(devContext);
  
  copyNodeDataToGroupData.setContext(devContext);
 
  /* load kernels */
  
#ifdef USE_CUDA
//   propsNonLeaf.load_source("./compute_properties.ptx", pathName.c_str());
//   propsLeaf.load_source("./compute_properties.ptx", pathName.c_str());
//   propsScaling.load_source("./compute_properties.ptx", pathName.c_str());

  propsNonLeafD.load_source("./compute_propertiesD.ptx", pathName.c_str(), "", -1, CU_TARGET_COMPUTE_20);
  propsLeafD.load_source("./compute_propertiesD.ptx", pathName.c_str(), "", -1, CU_TARGET_COMPUTE_20);
  propsScalingD.load_source("./compute_propertiesD.ptx", pathName.c_str(), "",-1, CU_TARGET_COMPUTE_20);
  
  copyNodeDataToGroupData.load_source("./compute_propertiesD.ptx", pathName.c_str());

  /* create kernels */

//   propsNonLeaf.create("compute_non_leaf"); 
//   propsLeaf.create("compute_leaf");
//   propsScaling.create("compute_scaling");
  
  propsNonLeafD.create("compute_non_leaf"); 
  propsLeafD.create("compute_leaf");
  propsScalingD.create("compute_scaling");

  copyNodeDataToGroupData.create("setPHGroupData");
  
#else
//   propsNonLeaf.load_source("compProps.cl", "");
//   propsLeaf.load_source("compProps.cl", "");
//   propsScaling.load_source("compProps.cl", ""); 
  
  /* create kernels */

  propsNonLeaf.create("compute_non_leaf");
  propsLeaf.create("compute_leaf");
  propsScaling.create("compute_scaling");  
#endif

  
  /* Tree iteration */
  /*
  getTNext.setContext(devContext);
  predictParticles.setContext(devContext);
  getNActive.setContext(devContext);
  correctParticles.setContext(devContext);
  computeDt.setContext(devContext);
  computeEnergy.setContext(devContext);
  setActiveGrps.setContext(devContext);
  distanceCheck.setContext(devContext);  
  approxGravLET.setContext(devContext);
  */
  
  approxGrav.setContext(devContext);

#ifdef USE_CUDA
  /*
  getTNext.load_source("./timestep.ptx", pathName.c_str(), "", -1, CU_TARGET_COMPUTE_20);
  predictParticles.load_source("./timestep.ptx", pathName.c_str(), "", -1, CU_TARGET_COMPUTE_20);
  getNActive.load_source("./timestep.ptx", pathName.c_str(), "", -1, CU_TARGET_COMPUTE_20);
//   approxGrav.load_source("./dev_approximate_gravity.cubin", "", "", 64);
  correctParticles.load_source("./timestep.ptx", pathName.c_str(), "", -1, CU_TARGET_COMPUTE_20);
  computeDt.load_source("./timestep.ptx", pathName.c_str(), "", -1, CU_TARGET_COMPUTE_20);
  computeEnergy.load_source("./timestep.ptx", pathName.c_str(), "", -1, CU_TARGET_COMPUTE_20);
  setActiveGrps.load_source("./timestep.ptx", pathName.c_str(), "", -1, CU_TARGET_COMPUTE_20);
  distanceCheck.load_source("./timestep.ptx", pathName.c_str(), "", -1, CU_TARGET_COMPUTE_20);  
  
  approxGravLET.load_source("./dev_approximate_gravity_let.ptx", pathName.c_str(), "", 64);  
  */
  
  approxGrav.load_source("./dev_approximate_gravity.ptx", pathName.c_str(), "", 64);
  //approxGrav.load_source("./dev_approximate_gravity.ptx", pathName.c_str(), "", 68); // modified by M.I. ( 64+sizeof(float) )
  /* create kernels */
/*
  getTNext.create("get_Tnext"); 
  predictParticles.create("predict_particles");   
  getNActive.create("get_nactive");
  correctParticles.create("correct_particles");
  computeDt.create("compute_dt");  
  setActiveGrps.create("setActiveGroups");

  computeEnergy.create("compute_energy_double");  
  distanceCheck.create("distanceCheck");  
  
  approxGravLET.create("dev_approximate_gravity");
  */
  approxGrav.create("dev_approximate_gravity");

#else
  getTNext.load_source("", "");
  
  /* create kernels */

  getTNext.create("");
  
#endif

#if 0
  //Parallel kernels
  domainCheck.setContext(devContext);  
  extractSampleParticles.setContext(devContext);  
  extractOutOfDomainR4.setContext(devContext);  
  extractOutOfDomainBody.setContext(devContext);  
  insertNewParticles.setContext(devContext);  
  internalMove.setContext(devContext);  

#ifdef USE_CUDA
  domainCheck.load_source("./parallel.ptx", pathName.c_str());
  extractSampleParticles.load_source("./parallel.ptx", pathName.c_str());
  extractOutOfDomainR4.load_source("./parallel.ptx", pathName.c_str());
  extractOutOfDomainBody.load_source("./parallel.ptx", pathName.c_str());
  insertNewParticles.load_source("./parallel.ptx", pathName.c_str());
  internalMove.load_source("./parallel.ptx", pathName.c_str());
  
  domainCheck.create("doDomainCheck");
  extractSampleParticles.create("extractSampleParticles");
  extractOutOfDomainR4.create("extractOutOfDomainParticlesR4");
  extractOutOfDomainBody.create("extractOutOfDomainParticlesAdvanced");
  insertNewParticles.create("insertNewParticles");
  internalMove.create("internalMove");

#else


#endif

#endif

}

