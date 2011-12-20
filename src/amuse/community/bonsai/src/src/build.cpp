#include "octree.h"

void octree::allocateParticleMemory(tree_structure &tree)
{
  //Allocates the memory to hold the particles data
  //and the arrays that have the same size as there are
  //particles. Eg valid arrays used in tree construction
  int n_bodies = tree.n;
  
  if(nProcs > 1)                //10% extra space, only in parallel when
    n_bodies = n_bodies*1.1;    //number of particles can fluctuate
  
  //Particle properties
  tree.bodies_pos.cmalloc(n_bodies+1, true);   //+1 to set end pos, host mapped? TODO not needed right since we use Ppos
  tree.bodies_key.cmalloc(n_bodies+1, false);   //+1 to set end key
  tree.bodies_ids.cmalloc(n_bodies+1, false);   //+1 to set end key
  
  tree.bodies_Ppos.cmalloc(n_bodies+1, true);   //Memory to store predicted positions, host mapped
  tree.bodies_Pvel.cmalloc(n_bodies+1, true);   //Memory to store predicted velocities, host mapped
  
  tree.bodies_vel.cmalloc(n_bodies, false);     
  tree.bodies_acc0.ccalloc(n_bodies, false);    //ccalloc -> init to 0
  tree.bodies_acc1.ccalloc(n_bodies, false);    //ccalloc -> init to 0
  tree.bodies_time.ccalloc(n_bodies, false);    //ccalloc -> init to 0

  tree.oriParticleOrder.cmalloc(n_bodies, false);      //To desort the bodies tree later on
  //iteration properties / information
  tree.activePartlist.ccalloc(n_bodies+2, false);   //+2 since we use the last two values as a atomicCounter (for grp count and semaphore access)
  tree.ngb.ccalloc(n_bodies, false);  
  tree.interactions.cmalloc(n_bodies, false);
  
  tree.body2group_list.cmalloc(n_bodies, false);
  
  tree.level_list.cmalloc(MAXLEVELS);  
  
  //The generalBuffer is also used during the tree-walk, so the size has to be at least
  //large enough to store the tree-walk stack. Add 4096 for extra memory allignment space
  //Times 2 since first half is used for regular walks, 2nd half for walks that go outside 
  //the memory stack and require another walk with 'unlimited' stack size
  int treeWalkStackSize = (2*LMEM_STACK_SIZE*NTHREAD*nBlocksForTreeWalk) + 4096;
    

  int tempSize   = max(n_bodies, 4096);   //Use minium of 4096 to prevent offsets mess up with small N
  tempSize       = 3*tempSize *4 + 4096;  //Add 4096 to give some space for memory allignment  
  tempSize = max(tempSize, treeWalkStackSize);
  
  //General buffer is used at multiple locations and reused in different functions
  tree.generalBuffer1.cmalloc(tempSize, true);  
  

  //Interactions needs; uint2   n 
  //ngb          needs:  uint   n
  //activePartList needs: uint  n  (to be used at same moment as interactions and ngb)
  
  //So we need at least (uint4+uint+uint)*n is uint4*n of memory
  //to be used at the same time, and to be used at the same time as 
  //multipole data


  #if 0
  Vanuit sort bodies

  Keys worden berekend voor de sort bodies, moet het dan nog wel
  in de build.cpp ? Roepen we die niet altijd in combinatie aan?
  

  Could combine it with multiple moments
  none of these are used at the same time
  if we sort, we compute new multiple moments
 
  Have to increase the size of multiple moments

  #endif


  //Tree properties, tree size is not known at forehand so
  //allocate worst possible outcome  
  n_bodies = n_bodies / 1;
  tree.node_key.cmalloc(n_bodies, false);
  tree.n_children.cmalloc(n_bodies, false);
  tree.node_bodies.cmalloc(n_bodies, false);  
  
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
  
  if(mpiGetNProcs() > 1)
  {
//    int remoteSize = (n_bodies*0.1) +  (n_bodies*0.1); //TODO some more realistic number
    int remoteSize = (n_bodies*0.5); //TODO some more realistic number
    this->remoteTree.fullRemoteTree.cmalloc(remoteSize, true);    
  }
  
}


void octree::reallocateParticleMemory(tree_structure &tree)
{
  //Realloc the memory to hold the particles data
  //and the arrays that have the same size as there are
  //particles. Eg valid arrays used in tree construction
  int n_bodies = tree.n;
  
  bool reduce = false;  //Set this to true to limit memory usage by only allocating what
                        //is required. If its false, then memory is not reduced and a larger
                        //buffer is kept
  
  //Particle properties
  tree.bodies_pos.cresize(n_bodies+1, reduce);   //+1 to set boundary condition
  tree.bodies_key.cresize(n_bodies+1, reduce);   //+1 to set boundary condition
  tree.bodies_ids.cresize(n_bodies+1, reduce);   //

  tree.bodies_Ppos.cresize(n_bodies+1, reduce);   //Memory to store predicted positions
  tree.bodies_Pvel.cresize(n_bodies+1, reduce);   //Memory to store predicted velocities  
  
  tree.bodies_vel.cresize (n_bodies, reduce);     
  tree.bodies_acc0.cresize(n_bodies, reduce);    //ccalloc -> init to 0
  tree.bodies_acc1.cresize(n_bodies, reduce);    //ccalloc -> init to 0
  tree.bodies_time.cresize(n_bodies, reduce);    //ccalloc -> init to 0
  
  tree.oriParticleOrder.cresize(n_bodies, reduce);     //To desort the bodies tree later on         
  //iteration properties / information
  tree.activePartlist.cresize(n_bodies+1, reduce);      //+1 since we use the last value as a atomicCounter
  tree.ngb.cresize(n_bodies, reduce);  
  tree.interactions.cresize(n_bodies, reduce);
  
  tree.body2group_list.cresize(n_bodies, reduce);  

  
  //Tree properties, tree size is not known at forehand so
  //allocate worst possible outcome  
  n_bodies = n_bodies / 1;
  tree.node_key.cresize(n_bodies, reduce);
  tree.n_children.cresize(n_bodies, reduce);
  tree.node_bodies.cresize(n_bodies, reduce);  
  
  my_dev::base_mem::printMemUsage();
}

void octree::allocateTreePropMemory(tree_structure &tree)
{ 
  int n_nodes = tree.n_nodes;

  //Allocate memory
  if(tree.groupCenterInfo.get_size() > 0)
  {
    //Resize, so we dont alloc if we already have mem alloced
    tree.multipole.cresize(3*n_nodes,     false);
    
    tree.boxSizeInfo.cresize(n_nodes,     false);  //host alloced
    tree.groupSizeInfo.cresize(n_nodes,   false);     

    tree.boxCenterInfo.cresize(n_nodes,   false); //host alloced
    tree.groupCenterInfo.cresize(n_nodes, false);    
  }
  else
  {    
    n_nodes = n_nodes * 1.1; 
    tree.multipole.cmalloc(3*n_nodes, true); //host alloced
        
    tree.boxSizeInfo.cmalloc(n_nodes, true);     //host alloced
    tree.groupSizeInfo.cmalloc(n_nodes, false);     

    tree.boxCenterInfo.cmalloc(n_nodes, true); //host alloced
    tree.groupCenterInfo.cmalloc(n_nodes,false);
  }
}



void octree::build (tree_structure &tree) {

  int level      = 0;
  int validCount = 0;
  int offset     = 0;


  /******** load kernels **********/

  /******** create memory buffers **********/


  my_dev::dev_mem<uint>  validList(devContext);
  my_dev::dev_mem<uint>  compactList(devContext);
  
  validList.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                                    tree.generalBuffer1.get_flags(), 
                                    tree.generalBuffer1.get_devMem(),
                                    &tree.generalBuffer1[0], 0,
                                    tree.n*2, getAllignmentOffset(0));
  validList.zeroMem(); 
                                    
  compactList.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                                    tree.generalBuffer1.get_flags(), 
                                    tree.generalBuffer1.get_devMem(),
                                    &tree.generalBuffer1[tree.n*2], tree.n*2,
                                    tree.n*2, getAllignmentOffset(tree.n*2));
                                                                        

  
  /******** set kernels parameters **********/
  

  build_key_list.set_arg<cl_mem>(0,   tree.bodies_key.p());
  build_key_list.set_arg<cl_mem>(1,   tree.bodies_Ppos.p());
  build_key_list.set_arg<int>(2,      &tree.n);
  build_key_list.set_arg<real4>(3,    &tree.corner);
  build_key_list.setWork(tree.n, 128);
  
  build_valid_list.set_arg<int>(0, &tree.n);
  build_valid_list.set_arg<int>(1, &level);
  build_valid_list.set_arg<cl_mem>(2,  tree.bodies_key.p());
  build_valid_list.set_arg<cl_mem>(3,  validList.p());  
  build_valid_list.setWork(tree.n, 128);
  


  build_nodes.set_arg<int>(0,     &level);
  build_nodes.set_arg<int>(1,     &validCount);
  build_nodes.set_arg<int>(2,     &offset);
  build_nodes.set_arg<cl_mem>(3,  compactList.p());
  build_nodes.set_arg<cl_mem>(4,  tree.bodies_key.p());
  build_nodes.set_arg<cl_mem>(5,  tree.node_key.p());
  build_nodes.set_arg<cl_mem>(6,  tree.n_children.p());
  build_nodes.set_arg<cl_mem>(7,  tree.node_bodies.p());

  link_tree.set_arg<int>(0,     &offset);
  link_tree.set_arg<cl_mem>(1,  tree.n_children.p());
  link_tree.set_arg<cl_mem>(2,  tree.node_bodies.p());
  link_tree.set_arg<cl_mem>(3,  tree.bodies_Ppos.p());
  link_tree.set_arg<real4>(4,   &tree.corner);
  link_tree.set_arg<cl_mem>(5,  tree.level_list.p());
  link_tree.set_arg<cl_mem>(6,  validList.p()); 
  link_tree.set_arg<cl_mem>(7,  tree.node_key.p());
  link_tree.set_arg<cl_mem>(8,  tree.bodies_key.p());
  link_tree.set_arg<int>(9,     &level);


  /********** build  list of keys ********/
  
  build_key_list.execute();  
  
  /******  build the levels *********/
  
  int nodeSum = 0;
  for (level = 0; level < MAXLEVELS; level++) {
    // mark bodies to be combined into nodes
    build_valid_list.set_arg<int>(1, &level);
    build_valid_list.execute();
      
    //gpuCompact to get number of created nodes    
    gpuCompact(devContext, validList, compactList, tree.n*2, &validCount);
                 
    nodeSum += validCount / 2;
    printf("ValidCount (%d): %d \tSum: %d Offset: %d\n", mpiGetRank(), validCount, nodeSum, offset);
    
    validCount /= 2;     
                  
    if (validCount == 0) break;                 
      
    // asssemble nodes           
    build_nodes.setWork(validCount, 128);
    build_nodes.set_arg<int>(0, &level);
    build_nodes.set_arg<int>(1, &validCount);
    build_nodes.set_arg<int>(2, &offset);    
    build_nodes.execute();
                 
    tree.level_list[level] = (uint2){offset, offset + validCount};
    offset += validCount;

  } //end for lvl
  

  //Put the last level + 1 index to 0,0 
  //so we dont need an extra if statement in the linking phase
  tree.level_list[level] = (uint2){0, 0};
  tree.level_list.h2d();
    
  int n_nodes  = offset;
  tree.n_nodes = n_nodes;
  
 
  /***** Link the tree ******/
  
  link_tree.set_arg<int>(0, &offset);   //Offset=number of nodes
  link_tree.set_arg<int>(9, &level);   //level=highest number of levels
  
  //The maximum number of levels that can be used is MAXLEVEl 
  //if max level is larger than that the program will exit
  printf("Max level : %d \n", level);
  if(level >= MAXLEVELS)
  {
    cerr << "The tree has become too deep, the program will exit. \n";
    cerr << "Consider the removal of far away particles to prevent a too large box. \n";
    exit(0);
  }
  
  link_tree.setWork(n_nodes, 128);
  printf("Link_tree: "); link_tree.printWorkSize();
  
  tree.n_levels = level-1;

  for(int i=0; i < level; i++)
    printf("%d\t%d\t%d\n", i, tree.level_list[i].x, tree.level_list[i].y);
 
  //Link the tree      
  link_tree.execute();
  

  //After executing link_tree, the id_list contains for each node
  //the ID of its parent.
  //Valid_list contains for each node if its a leaf (valid) or a normal
  //node -> non_valid
  //Execute a split on the validList to get seperate id lists 
  //for the leafs and nodes. Used when computing multipoles
    
  tree.leafNodeIdx.cmalloc(tree.n_nodes , false);
  
  //Split the leaf ids and non-leaf node ids
  gpuSplit(devContext, validList, tree.leafNodeIdx, tree.n_nodes, &tree.n_leafs);     
                 
  printf("Total nodes: %d N_leafs: %d  non-leafs: %d \n", tree.n_nodes, tree.n_leafs, tree.n_nodes - tree.n_leafs);
  

  build_level_list.set_arg<int>(0, &tree.n_nodes);
  build_level_list.set_arg<int>(1, &tree.n_leafs);
  build_level_list.set_arg<cl_mem>(2, tree.leafNodeIdx.p());
  build_level_list.set_arg<cl_mem>(3, tree.node_bodies.p());
  build_level_list.set_arg<cl_mem>(4, validList.p());  
  build_level_list.setWork(tree.n_nodes-tree.n_leafs, 128);
  
  validList.zeroMem();  

  //Build the level list based on the leafIdx list
  //required for easy access in the compute node properties
  build_level_list.execute();  

  tree.node_level_list.cmalloc(level*2 , false);

  int levelThing;
  
  gpuCompact(devContext, validList, tree.node_level_list, 
             2*(tree.n_nodes-tree.n_leafs), &levelThing);             
  
  tree.node_level_list.d2h();
  
  //We only care about end positions, so compress the list:
  int j=0;
  for(int i=0; i < levelThing; i+=2, j++)
    tree.node_level_list[j] = tree.node_level_list[i];
  
  tree.node_level_list[j] =tree.node_level_list[levelThing-1]+1; //Add 1 to make it the end position
  levelThing = j+1;
  tree.node_level_list.h2d();
  
  printf("Finished level list \n");
  
  for(int i=0; i < levelThing; i++)
  {
    printf("node_level_list: %d \t%d\n", i, tree.node_level_list[i]);
  }
  
  ///******   Start building the particle groups *******///////

  //Compute the box size, the max length of one of the sides of the rectangle
  real size     = fmax(fabs(rMaxLocalTree.z - rMinLocalTree.z), 
                  fmax(fabs(rMaxLocalTree.y - rMinLocalTree.y),
                       fabs(rMaxLocalTree.x - rMinLocalTree.x)));
  real dist     = ((rMaxLocalTree.z - rMinLocalTree.z) * (rMaxLocalTree.z - rMinLocalTree.z) + 
                   (rMaxLocalTree.y - rMinLocalTree.y) * (rMaxLocalTree.y - rMinLocalTree.y) +
                   (rMaxLocalTree.x - rMinLocalTree.x) * (rMaxLocalTree.x - rMinLocalTree.x));      
                   
  float maxDist = sqrt(dist) / 10;
  maxDist *= maxDist; //Square since we dont do sqrt on device
                       
  fprintf(stderr,"Box max size: %f en max dist: %f \t %f en %f  \n", size, dist, sqrt(dist), maxDist);
  
  //maxDist = 50;
  
  validList.zeroMem();
  //The newest group creation method!
  define_groups.set_arg<int>(0, &tree.n);  
  define_groups.set_arg<cl_mem>(1, validList.p());    
  define_groups.set_arg<cl_mem>(2, tree.bodies_Ppos.p());
  define_groups.set_arg<float>(3, &maxDist);     
  define_groups.setWork(tree.n, 128);  
  define_groups.execute();
  
  //gpuCompact    
  gpuCompact(devContext, validList, compactList, tree.n*2, &validCount);
  
  printf("Found number of groups: %d \n", validCount/2);

  tree.n_groups = validCount/2;
  //Now compact validList to get the list of group ids
  tree.group_list_test.cmalloc(tree.n_groups , false);  
  
  store_groups.set_arg<int>(0, &tree.n);  
  store_groups.set_arg<int>(1, &tree.n_groups);  
  store_groups.set_arg<cl_mem>(2, compactList.p());    
  store_groups.set_arg<cl_mem>(3, tree.body2group_list.p());     
  store_groups.set_arg<cl_mem>(4, tree.group_list_test.p());     
  store_groups.setWork(-1, NCRIT, tree.n_groups);  
  store_groups.execute();  

  //Memory allocation for the valid group lists
  if(tree.active_group_list.get_size() > 0)
  {
    tree.active_group_list.cresize(tree.n_groups, false);
    tree.activeGrpList.cresize(tree.n_groups, false);     
    
  }
  else
  {
    tree.active_group_list.cmalloc(tree.n_groups, false);
    tree.activeGrpList.cmalloc(tree.n_groups, false);     
  }


  printf("Tree built complete!\n");

  /*************************/

}
