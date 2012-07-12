#include "octree.h"

void octree::compute_properties(tree_structure &tree,  
				my_dev::dev_mem<float4>  &bodies_pos, 
				int n_bodies) {

  /*****************************************************          
    Assign the memory buffers, note that we check the size first
    and if needed we increase the size of the generalBuffer1
    Size required:
      - multipoleD -> double4*3_n_nodes -> 2*3*n_nodes*uint4 
      - lower/upperbounds ->               2*n_nodes*uint4
      - node lower/upper  ->               2*n_nodes*uint4
      - SUM:                               10*n_nodes*uint4 
      - generalBuffer1 has default size: 3*N*uint4
      
    check if 10*n_nodes < 3*N if so realloc
    (Note that generalBuffer might be larger because of tree-walk stack)
   *****************************************************/
  
  if(10*tree.n_nodes > 3*tree.n)
  {
    fprintf(stderr, "Resizeing the generalBuffer1 \n");
    tree.generalBuffer1.cresize(10*tree.n_nodes*4, false);
  }
  
  my_dev::dev_mem<double4> multipoleD(devContext);
  my_dev::dev_mem<real4>  nodeLowerBounds(devContext); //Lower bounds used for scaling? TODO
  my_dev::dev_mem<real4>  nodeUpperBounds(devContext); //Upper bounds used for scaling? TODO  
  
  multipoleD.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                          tree.generalBuffer1.get_flags(), 
                          tree.generalBuffer1.get_devMem(),
                          &tree.generalBuffer1[0], 0, 
                          3*tree.n_nodes, getAllignmentOffset(0));

  //Offset is in uint, so: double4 = 8uint*3*n_nodes
  nodeLowerBounds.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                          tree.generalBuffer1.get_flags(), 
                          tree.generalBuffer1.get_devMem(),
                          &tree.generalBuffer1[8*3*tree.n_nodes],  8*3*tree.n_nodes,
                          tree.n_nodes, getAllignmentOffset(8*3*tree.n_nodes));
                         
  int prevOffsetSum = getAllignmentOffset(8*3*tree.n_nodes); //The offset of output
                          
  nodeUpperBounds.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                          tree.generalBuffer1.get_flags(), 
                          tree.generalBuffer1.get_devMem(),
                          &tree.generalBuffer1[8*3*tree.n_nodes + 4*tree.n_nodes], 
                          8*3*tree.n_nodes + 4*tree.n_nodes, 
                          tree.n_nodes, 
                          prevOffsetSum + getAllignmentOffset(8*3*tree.n_nodes + 4*tree.n_nodes + prevOffsetSum));     
       
  //Computes the tree-properties (size, cm, monopole, quadropole, etc)
  //start the kernel for the leaf-type nodes
  propsLeafD.set_arg<int>(0,    &tree.n_leafs);
  propsLeafD.set_arg<cl_mem>(1, tree.leafNodeIdx.p());
  propsLeafD.set_arg<cl_mem>(2, tree.node_bodies.p());
  propsLeafD.set_arg<cl_mem>(3, bodies_pos.p());
  propsLeafD.set_arg<cl_mem>(4, multipoleD.p());
  propsLeafD.set_arg<cl_mem>(5, nodeLowerBounds.p());
  propsLeafD.set_arg<cl_mem>(6, nodeUpperBounds.p());
//   propsLeafD.set_arg<cl_mem>(7, tree.bodies_Pvel.p());  //Velocity to get max eps
  
  //printf("tree.node_bodies.p()=%p, \n", tree.node_bodies.p());

  propsLeafD.setWork(tree.n_leafs, 128);
  printf("PropsLeaf: "); propsLeafD.printWorkSize();
  propsLeafD.execute(); 
   
  
  int temp = tree.n_nodes-tree.n_leafs;
  propsNonLeafD.set_arg<int>(0,    &temp);
  propsNonLeafD.set_arg<cl_mem>(1, tree.leafNodeIdx.p());
  propsNonLeafD.set_arg<cl_mem>(2, tree.node_level_list.p());
  propsNonLeafD.set_arg<cl_mem>(3, tree.n_children.p());  
  propsNonLeafD.set_arg<cl_mem>(4, multipoleD.p());
  propsNonLeafD.set_arg<cl_mem>(5, nodeLowerBounds.p());
  propsNonLeafD.set_arg<cl_mem>(6, nodeUpperBounds.p());

  //Work from the bottom up
  for(int i=tree.n_levels; i >= 1; i--)
  {   
      propsNonLeafD.set_arg<int>(0,    &i);  
      {    
        vector<size_t> localWork(2), globalWork(2);
        int totalOnThisLevel;
      
        totalOnThisLevel = tree.node_level_list[i]-tree.node_level_list[i-1];
        
        propsNonLeafD.setWork(totalOnThisLevel, 128);
        
        printf("PropsNonLeaf, nodes on level %d : %d (start: %d end: %d) , config: \t", i, totalOnThisLevel,
               tree.node_level_list[i-1], tree.node_level_list[i]); 
        propsNonLeafD.printWorkSize();
      }      
      propsNonLeafD.set_arg<int>(0,    &i); //set the level
      propsNonLeafD.execute();     
  }
  
  propsScalingD.set_arg<int>(0,    &tree.n_nodes);
  propsScalingD.set_arg<cl_mem>(1, multipoleD.p());
  propsScalingD.set_arg<cl_mem>(2, nodeLowerBounds.p());
  propsScalingD.set_arg<cl_mem>(3, nodeUpperBounds.p());
  propsScalingD.set_arg<cl_mem>(4, tree.n_children.p());  
  propsScalingD.set_arg<cl_mem>(5, tree.multipole.p());
  propsScalingD.set_arg<float >(6, &theta);
  propsScalingD.set_arg<cl_mem>(7, tree.boxSizeInfo.p());
  propsScalingD.set_arg<cl_mem>(8, tree.boxCenterInfo.p());
  propsScalingD.set_arg<cl_mem>(9, tree.node_bodies.p());
  
  propsScalingD.setWork(tree.n_nodes, 128);
  printf("propsScaling: \t "); propsScalingD.printWorkSize();
  propsScalingD.execute();   


#if 0
  #ifdef INDSOFT
    //If we use individual softening we need to get the max softening value
    //to be broadcasted during the exchange of the LET boundaries.
    //Only copy the root node that contains the max value
    my_dev::dev_stream memCpyStream;
    tree.multipole.d2h(3, false, memCpyStream.s());
  #endif

        
  //Set the group properties, note that it is not based on the nodes anymore
  //but on self created groups based on particle order setPHGroupData    
  copyNodeDataToGroupData.set_arg<int>(0,    &tree.n_groups);
  copyNodeDataToGroupData.set_arg<int>(1,    &tree.n);
  copyNodeDataToGroupData.set_arg<cl_mem>(2, tree.bodies_Ppos.p());  
  copyNodeDataToGroupData.set_arg<cl_mem>(3, tree.group_list_test.p());
  copyNodeDataToGroupData.set_arg<cl_mem>(4, tree.groupCenterInfo.p());  
  copyNodeDataToGroupData.set_arg<cl_mem>(5, tree.groupSizeInfo.p());
  copyNodeDataToGroupData.setWork(-1, NCRIT, tree.n_groups);    
  copyNodeDataToGroupData.printWorkSize();
  copyNodeDataToGroupData.execute();
  
  #ifdef INDSOFT  
    memCpyStream.sync();  
    this->maxLocalEps = tree.multipole[0*3 + 1].w; //Softening value
  #else
  #endif
  
  //Get the local domain boundary based on group positions and sizes
  real4 r_min, r_max;
  getBoundariesGroups(tree, r_min, r_max); 
  
#endif  
  
}
