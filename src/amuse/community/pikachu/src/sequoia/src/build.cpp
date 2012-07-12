#include "octree.h"


void octree::build(tree_structure &tree, 
		   my_dev::dev_mem<float4>  &bodies_pos, 
		   int n_bodies) {

    int level      = 0;
    int validCount = 0;
    int offset     = 0;

    /******** create memory buffers **********/
  
    //Make sure that buffers that rely on the n_bodies count have the minimal size
    this->allocateParticleSpecificBuffers(n_bodies);

    my_dev::dev_mem<uint>  validList(devContext);
    my_dev::dev_mem<uint>  compactList(devContext);
    my_dev::dev_mem<uint4>  bodies_key(devContext);
    my_dev::dev_mem<uint4>  node_key(devContext);
  
    validList.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
			   tree.generalBuffer1.get_flags(), 
			   tree.generalBuffer1.get_devMem(),
			   &tree.generalBuffer1[0], 0,
			   n_bodies*2, getAllignmentOffset(0));
    validList.zeroMem(); 
                                    
  compactList.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                                    tree.generalBuffer1.get_flags(), 
                                    tree.generalBuffer1.get_devMem(),
                                    &tree.generalBuffer1[n_bodies*2], n_bodies*2,
                                    n_bodies*2, getAllignmentOffset(n_bodies*2));  
  int prevOffset = getAllignmentOffset(n_bodies*2);
  
  
  node_key.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                                    tree.generalBuffer1.get_flags(), 
                                    tree.generalBuffer1.get_devMem(),
                                    &tree.generalBuffer1[n_bodies*4], n_bodies*4,
                                    n_bodies, prevOffset +                          //n_bodies items of uint4
                                    getAllignmentOffset(n_bodies*4 + prevOffset));
  
  prevOffset += getAllignmentOffset(n_bodies*4 + prevOffset);


  bodies_key.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                                    tree.generalBuffer1.get_flags(), 
                                    tree.generalBuffer1.get_devMem(),
                                    &tree.generalBuffer1[n_bodies*8], n_bodies*8,
                                    n_bodies, prevOffset +                          //n_bodies items of uint4
                                    getAllignmentOffset(n_bodies*8  + prevOffset));
  
  
  
  //Total amount of space in the generalbuffer: at least: 
  //- 3x n_bodies*uint4/float4 , so 12*int size
  // validList          = 0-2*n_bodies (int)             Remaining: 10*n_bodies
  // compactList        = 2*n_bodies - 4*n_bodies (int)  Remaining:  8*n_bodies
  // node_key           = 4*n_bodies - 8*n_bodies (int)  Remaining:  4*n_bodies
  // bodies_key         = 8*n_bodies = 12*n_bodies (int) Remaining:  0*n_bodies
  
  
  /********** Reserve memory for the build in buffers ****************/

  
  /******** set kernels parameters **********/
  
 
  build_valid_list.set_arg<int>(0,     &n_bodies);
  build_valid_list.set_arg<int>(1,     &level);
  build_valid_list.set_arg<cl_mem>(2,  bodies_key.p());
  build_valid_list.set_arg<cl_mem>(3,  validList.p());  
  build_valid_list.setWork(n_bodies, 128);

  build_nodes.set_arg<int>(0,     &level);
  build_nodes.set_arg<int>(1,     &validCount);
  build_nodes.set_arg<int>(2,     &offset);
  build_nodes.set_arg<cl_mem>(3,  compactList.p());
  build_nodes.set_arg<cl_mem>(4,  bodies_key.p());
  build_nodes.set_arg<cl_mem>(5,  node_key.p());
  build_nodes.set_arg<cl_mem>(6,  tree.n_children.p());
  build_nodes.set_arg<cl_mem>(7,  tree.node_bodies.p());

  link_tree.set_arg<int>(0,     &offset);
  link_tree.set_arg<cl_mem>(1,  tree.n_children.p());
  link_tree.set_arg<cl_mem>(2,  tree.node_bodies.p());
  link_tree.set_arg<cl_mem>(3,  bodies_pos.p());
  link_tree.set_arg<real4>(4,   &tree.corner);
  link_tree.set_arg<cl_mem>(5,  tree.level_list.p());
  link_tree.set_arg<cl_mem>(6,  validList.p()); 
  link_tree.set_arg<cl_mem>(7,  node_key.p());
  link_tree.set_arg<cl_mem>(8,  bodies_key.p());
  link_tree.set_arg<int>(9,     &level);


  /********** build  list of keys ********/
  
  this->compute_keys(bodies_pos, bodies_key, n_bodies, tree.corner, tree.domain_fac);
  
  
  /******  build the levels *********/
  
  int nodeSum = 0;
  for (level = 0; level < MAXLEVELS; level++) {
    // mark bodies to be combined into nodes
    build_valid_list.set_arg<int>(1, &level);
    build_valid_list.execute();
      
    //gpuCompact to get number of created nodes    
    gpuCompact(devContext, validList, compactList, n_bodies*2, &validCount);
                 
    nodeSum += validCount / 2;
    printf("ValidCount (%d): %d \tSum: %d Offset: %d\n", 0, validCount, nodeSum, offset);
    
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
  
  //Allocate the buffers that store node properties
  this->allocateNodeSpecificBuffers(tree.n_nodes);
    
  
  //Split the leaf ids and non-leaf node ids
  gpuSplit(devContext, validList, tree.leafNodeIdx, tree.n_nodes, &tree.n_leafs);     

  printf("Total nodes: %d N_leafs: %d  non-leafs: %d \n", 
	 tree.n_nodes, tree.n_leafs, tree.n_nodes - tree.n_leafs);

  build_level_list.set_arg<int>(0, &tree.n_nodes);
  build_level_list.set_arg<int>(1, &tree.n_leafs);
  build_level_list.set_arg<cl_mem>(2, tree.leafNodeIdx.p());
  build_level_list.set_arg<cl_mem>(3, tree.node_bodies.p());
  build_level_list.set_arg<cl_mem>(4, validList.p());  
  build_level_list.setWork(tree.n_nodes-tree.n_leafs, 128);
  
  validList.zeroMem();  

  //Build the level list based on the leafIdx list
  //required for easy access in the compute node properties
  build_level_list.execute();  // build_level_list (in build_tree.cu)

  int levelThing;
  
  gpuCompact(devContext, validList, tree.node_level_list, 
             2*(tree.n_nodes-tree.n_leafs), &levelThing);             
  
  tree.node_level_list.d2h();
  
  //We only care about end positions, so compress the list:
  int j=0;
  for(int i=0; i < levelThing; i+=2, j++)
    tree.node_level_list[j] = tree.node_level_list[i];
  
  tree.node_level_list[j] = tree.node_level_list[levelThing-1]+1; //Add 1 to make it the end position
  levelThing = j+1;
  tree.node_level_list.h2d();
  
  printf("Finished level list \n");
  
  for(int i=0; i < levelThing; i++){
      printf("node_level_list: %d \t%d\n", i, tree.node_level_list[i]);
  }

  /*************************/

}


void octree::createGroups(tree_structure &tree, my_dev::dev_mem<float4>  &bodies_pos, int n_bodies)
{
  //Tree-walks are only efficient if particles in similar space locations are group together
  //since there are sort routines available we assume here the particles are already sorted
  
  //It is possible that more particles are used for the groups than for the tree
  this->allocateParticleSpecificBuffers(n_bodies);
  
  my_dev::dev_mem<uint>  validList(devContext);
  my_dev::dev_mem<uint>  compactList(devContext);
  my_dev::dev_mem<uint4>  bodies_key(devContext);

  
  validList.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                                    tree.generalBuffer1.get_flags(), 
                                    tree.generalBuffer1.get_devMem(),
                                    &tree.generalBuffer1[0], 0,
                                    n_bodies*2, getAllignmentOffset(0));
  validList.zeroMem(); 
                                    
  compactList.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                                    tree.generalBuffer1.get_flags(), 
                                    tree.generalBuffer1.get_devMem(),
                                    &tree.generalBuffer1[n_bodies*2], n_bodies*2,
                                    n_bodies*2, getAllignmentOffset(n_bodies*2));  
  int prevOffset = getAllignmentOffset(n_bodies*2);
  
  
  bodies_key.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                                    tree.generalBuffer1.get_flags(), 
                                    tree.generalBuffer1.get_devMem(),
                                    &tree.generalBuffer1[n_bodies*4], n_bodies*4,
                                    n_bodies, prevOffset +                          //n_bodies items of uint4
                                    getAllignmentOffset(n_bodies*4 + prevOffset));  
  
  //Compute boundary and the keys
  float4 corner;
  float domain_fac;
  compute_keys(bodies_pos, bodies_key, n_bodies, corner, domain_fac);
  
  //Now build the groups, but first we have to determine the splitting level
  int level = 0;
  int validCount;
  
  
  build_valid_list.set_arg<int>(0,     &n_bodies);
  build_valid_list.set_arg<int>(1,     &level);
  build_valid_list.set_arg<cl_mem>(2,  bodies_key.p());
  build_valid_list.set_arg<cl_mem>(3,  validList.p());  
  build_valid_list.setWork(n_bodies, 128);  
  
  for (level = 0; level < MAXLEVELS; level++)
  {
    // mark bodies to be combined into nodes
    build_valid_list.set_arg<int>(1, &level);
    build_valid_list.execute();
      
    //gpuCompact to get number of created nodes    
    gpuCompact(devContext, validList, compactList, n_bodies*2, &validCount);
    validCount /= 2;     
    
    //This is a sort of arbirary number
    if(validCount > 50) break;    
  }

  
  //We can reuse the results of the build_valid_list
  //it already has the breaks in the right place 
  //Just add breaks every NGROUP items
  
  //The newest group creation method!
  define_groups.set_arg<int>(0, &n_bodies);  
  define_groups.set_arg<cl_mem>(1, validList.p());    
  define_groups.setWork(n_bodies, 256);  
  define_groups.execute();
  

  
  gpuCompact(devContext, validList, compactList, n_bodies*2, &validCount);
  int n_groups = validCount / 2;
  fprintf(stderr, "Ngroups: %d \n", n_groups);
  

  this->allocateGroupSpecificBuffers(n_bodies, n_groups);
   
  store_groups.set_arg<int>(0,     &n_groups);  
  store_groups.set_arg<cl_mem>(1,  compactList.p());    
  store_groups.set_arg<cl_mem>(2,  tree.body2group_list.p());     
  store_groups.set_arg<cl_mem>(3,  tree.group_list.p());     
  store_groups.setWork(-1, NCRIT,  n_groups);  
  store_groups.execute();

  tree.n_groups = n_groups;
  

  //Outcome is a list in group_list in which the x-element is the 
  //startparticle and y-element is the end-particle INCLUSIVE
}
#if 0

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
  
#endif  

