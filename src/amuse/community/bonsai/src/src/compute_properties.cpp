#include "octree.h"

void octree::compute_properties (tree_structure &tree) {
  
  #if 0
      fprintf(stderr,"This file is not up to date anymore! %s\n", __FILE__);
    exit(0);
  
  //Computes the tree-properties (size, cm, monopole, quadropole, etc)
  //start the kernel for the leaf-type nodes
  propsLeaf.set_arg<int>(0,    &tree.n_leafs);
  propsLeaf.set_arg<cl_mem>(1, tree.leafNodeIdx.p());
  propsLeaf.set_arg<cl_mem>(2, tree.node_bodies.p());
  propsLeaf.set_arg<cl_mem>(3, tree.bodies_Ppos.p());
//   propsLeaf.set_arg<cl_mem>(3, tree.bodies_pos.p());  
  propsLeaf.set_arg<cl_mem>(4, tree.multipole.p());
  propsLeaf.set_arg<cl_mem>(5, tree.nodeLowerBounds.p());
  propsLeaf.set_arg<cl_mem>(6, tree.nodeUpperBounds.p());
  propsLeaf.set_arg<cl_mem>(7, tree.lowerBounds.p());
  propsLeaf.set_arg<cl_mem>(8, tree.upperBounds.p());  
  propsLeaf.set_arg<cl_mem>(9, tree.bodies_Pvel.p());  //Velocity to get max eps  
  
  
  propsLeaf.setWork(tree.n_leafs, 128);
  printf("PropsLeaf: "); propsLeaf.printWorkSize();

  propsLeaf.execute();  


  int temp = tree.n_nodes-tree.n_leafs;
  propsNonLeaf.set_arg<int>(0,    &temp);
  propsNonLeaf.set_arg<cl_mem>(1, tree.leafNodeIdx.p());
  propsNonLeaf.set_arg<cl_mem>(2, tree.node_level_list.p());
  propsNonLeaf.set_arg<cl_mem>(3, tree.n_children.p());  
  propsNonLeaf.set_arg<cl_mem>(4, tree.multipole.p());
  propsNonLeaf.set_arg<cl_mem>(5, tree.nodeLowerBounds.p());
  propsNonLeaf.set_arg<cl_mem>(6, tree.nodeUpperBounds.p());


  for(int i=tree.n_levels; i >= 1; i--)
  {   
      propsNonLeaf.set_arg<int>(0,    &i);  
      {    
        vector<size_t> localWork(2), globalWork(2);
        int totalOnThisLevel;
      
        totalOnThisLevel = tree.node_level_list[i]-tree.node_level_list[i-1];
        
        propsNonLeaf.setWork(totalOnThisLevel, 128);
        
        printf("PropsNonLeaf, nodes on level %d : %d (start: %d end: %d) , config: \t", i, totalOnThisLevel,
               tree.node_level_list[i-1], tree.node_level_list[i]); 
        propsNonLeaf.printWorkSize();
      }      
      propsNonLeaf.set_arg<int>(0,    &i); //set the level
      propsNonLeaf.execute();     
  }
  
  float theta2 = theta;
  
  propsScaling.set_arg<int>(0,    &tree.n_nodes);
  propsScaling.set_arg<real4>(1,  &tree.corner);
  propsScaling.set_arg<cl_mem>(2, tree.multipole.p());
  propsScaling.set_arg<cl_mem>(3, tree.nodeLowerBounds.p());
  propsScaling.set_arg<cl_mem>(4, tree.nodeUpperBounds.p());
  propsScaling.set_arg<cl_mem>(5, tree.n_children.p());  
  propsScaling.set_arg<cl_mem>(6, tree.node_data.p());
  propsScaling.set_arg<float >(7, &theta2);
  propsScaling.set_arg<cl_mem>(8, tree.boxSizeInfo.p());
  propsScaling.set_arg<cl_mem>(9, tree.boxCenterInfo.p());
  
  propsScaling.setWork(tree.n_nodes, 128);
  printf("propsScaling: \t "); propsScaling.printWorkSize();
  propsScaling.execute();     


  //tree.multipole.d2h();
  //printf("COM: %f %f %f %f \n",tree.multipole[0].x, tree.multipole[0].y, tree.multipole[0].z, tree.multipole[0].w);


  #ifdef USE_CUDA
    cuCtxSynchronize();
  #else
    clFinish(devContext.get_command_queue());
  #endif
  
  tree.nodeLowerBounds.d2h();
  tree.nodeUpperBounds.d2h();
  
  
  
  copyNodeDataToGroupData.set_arg<int>(0,    &tree.n_groups);
  copyNodeDataToGroupData.set_arg<int>(1,    &tree.n_nodes);
  copyNodeDataToGroupData.set_arg<cl_mem>(2, tree.node_data.p());
  copyNodeDataToGroupData.set_arg<cl_mem>(3, tree.group_data.p());
  copyNodeDataToGroupData.set_arg<cl_mem>(4, tree.node_bodies.p());
  copyNodeDataToGroupData.set_arg<cl_mem>(5, tree.group_list.p());
  copyNodeDataToGroupData.set_arg<cl_mem>(6, tree.boxCenterInfo.p());
  copyNodeDataToGroupData.set_arg<cl_mem>(7, tree.boxSizeInfo.p());
  copyNodeDataToGroupData.set_arg<cl_mem>(8, tree.groupCenterInfo.p());
  copyNodeDataToGroupData.set_arg<cl_mem>(9, tree.groupSizeInfo.p());
  copyNodeDataToGroupData.setWork(tree.n_nodes, 128);

  printf("copyNodeDataToGroupData: \t "); copyNodeDataToGroupData.printWorkSize();
  copyNodeDataToGroupData.execute();  
 
//    tree.multipole.d2h();
//   testRes.d2h();
  
//   for(int i=0; i < tree.n_nodes; i++)
//   for(int i=tree.n_nodes-10; i < tree.n_nodes; i++)
/*   for(int i=0; i < 10; i++)
   {
     fprintf(stderr,"%d\t%f\t%f\t%f\t%f\n", i, tree.multipole[i*3+0].x,tree.multipole[i*3+0].y,tree.multipole[i*3+0].z, tree.multipole[i*3+0].w);
//     fprintf(stderr,"%d\t%f\t%f\t%f\t%f\t%f\n", i, tree.multipole[i*3+1].x,tree.multipole[i*3+1].y,tree.multipole[i*3+1].z, tree.multipole[i*3+1].w, testRes[i]);
    fprintf(stderr,"%d\t%f\t%f\t%f\t%f\t%f\n", i, tree.multipole[i*3+1].x,tree.multipole[i*3+1].y,tree.multipole[i*3+1].z, tree.multipole[i*3+1].w, 0);
    fprintf(stderr,"%d\t%f\t%f\t%f\t%f\n", i, tree.multipole[i*3+2].x,tree.multipole[i*3+2].y,tree.multipole[i*3+2].z, tree.multipole[i*3+2].w);
   }

exit(0);
*/
  #else

  compute_properties_double(tree);
  
  #endif
  
}

void octree::compute_properties_double(tree_structure &tree) {

  /*****************************************************          
    Assign the memory buffers, note that we check the size first
    and if needed we increase the size of the generalBuffer1
    Size required:
      - multipoleD -> double4*3_n_nodes -> 6*n_nodes*uint4 
      - lower/upperbounds ->               2*n_nodes*uint4
      - node lower/upper  ->               2*n_nodes*uint4
      - SUM: 10*n_nodes*uint4 
      - generalBuffer1 has default size: 3*N*uint4
      
    check if 10*n_nodes < 3*N if so realloc
    
   *****************************************************/
  
  if(10*tree.n_nodes > 3*tree.n)
  {
    fprintf(stderr, "Resizeing the generalBuffer1 \n");
    tree.generalBuffer1.cresize(8*tree.n_nodes*4, false);
  }
  
  my_dev::dev_mem<double4> multipoleD(devContext);
  
  multipoleD.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                          tree.generalBuffer1.get_flags(), 
                          tree.generalBuffer1.get_devMem(),
                          &tree.generalBuffer1[0], 0, 
                          3*tree.n_nodes, getAllignmentOffset(0));

  //Offset is in uint, so: double4 = 8uint*3*n_nodes
  tree.nodeLowerBounds.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                          tree.generalBuffer1.get_flags(), 
                          tree.generalBuffer1.get_devMem(),
                          &tree.generalBuffer1[8*3*tree.n_nodes],  8*3*tree.n_nodes,
                          tree.n_nodes, getAllignmentOffset(8*3*tree.n_nodes));
                         
  int prevOffsetSum = getAllignmentOffset(8*3*tree.n_nodes); //The offset of output
                          
  tree.nodeUpperBounds.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
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
  propsLeafD.set_arg<cl_mem>(3, tree.bodies_Ppos.p());  
  propsLeafD.set_arg<cl_mem>(4, multipoleD.p());
  propsLeafD.set_arg<cl_mem>(5, tree.nodeLowerBounds.p());
  propsLeafD.set_arg<cl_mem>(6, tree.nodeUpperBounds.p());
  propsLeafD.set_arg<cl_mem>(7, tree.bodies_Pvel.p());  //Velocity to get max eps
  
  propsLeafD.setWork(tree.n_leafs, 128);
  printf("PropsLeaf: "); propsLeafD.printWorkSize();
  propsLeafD.execute(); 
   
  
  int temp = tree.n_nodes-tree.n_leafs;
  propsNonLeafD.set_arg<int>(0,    &temp);
  propsNonLeafD.set_arg<cl_mem>(1, tree.leafNodeIdx.p());
  propsNonLeafD.set_arg<cl_mem>(2, tree.node_level_list.p());
  propsNonLeafD.set_arg<cl_mem>(3, tree.n_children.p());  
  propsNonLeafD.set_arg<cl_mem>(4, multipoleD.p());
  propsNonLeafD.set_arg<cl_mem>(5, tree.nodeLowerBounds.p());
  propsNonLeafD.set_arg<cl_mem>(6, tree.nodeUpperBounds.p());

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
  

  float theta2 = theta;
  
  propsScalingD.set_arg<int>(0,    &tree.n_nodes);
  propsScalingD.set_arg<real4>(1,  &tree.corner);
  propsScalingD.set_arg<cl_mem>(2, multipoleD.p());
  propsScalingD.set_arg<cl_mem>(3, tree.nodeLowerBounds.p());
  propsScalingD.set_arg<cl_mem>(4, tree.nodeUpperBounds.p());
  propsScalingD.set_arg<cl_mem>(5, tree.n_children.p());  
  propsScalingD.set_arg<cl_mem>(6, tree.multipole.p());
  propsScalingD.set_arg<float >(7, &theta2);
  propsScalingD.set_arg<cl_mem>(8, tree.boxSizeInfo.p());
  propsScalingD.set_arg<cl_mem>(9, tree.boxCenterInfo.p());
  propsScalingD.set_arg<cl_mem>(10, tree.node_bodies.p());
  
  propsScalingD.setWork(tree.n_nodes, 128);
  printf("propsScaling: \t "); propsScalingD.printWorkSize();
  propsScalingD.execute();   


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
  
  #if 0
    //Write the tree structure to file

    string nodeFileName = "fullTreeStructure.txt";
    char fileName[256];
    sprintf(fileName, "fullTreeStructure-%d.txt", mpiGetRank());
    ofstream nodeFile;
    //nodeFile.open(nodeFileName.c_str());
    nodeFile.open(fileName);
    tree.multipole.d2h();
    tree.boxSizeInfo.d2h();
    tree.boxCenterInfo.d2h();
    
    for(int i=0; i < tree.n_nodes; i++)
    {
      //nodeFile << i << "\t" << tree.boxCenterInfo[i].x << "\t" << tree.boxCenterInfo[i].y;
      //nodeFile << "\t" << 2*tree.boxSizeInfo[i].x << "\t" << 2*tree.boxSizeInfo[i].y << "\t";
      
      nodeFile << i << "\t" << tree.boxCenterInfo[i].x-tree.boxSizeInfo[i].x << "\t" << tree.boxCenterInfo[i].y-tree.boxSizeInfo[i].y;
      nodeFile << "\t"      << tree.boxCenterInfo[i].x+tree.boxSizeInfo[i].x << "\t" << tree.boxCenterInfo[i].y+tree.boxSizeInfo[i].y << "\t";
      
      
      nodeFile << tree.multipole[i*3+0].x << "\t" << tree.multipole[i*3+0].w << "\n";
    }
    
    nodeFile.close();
    

    sprintf(fileName, "fullTreeStructureParticles-%d.txt", mpiGetRank());
    ofstream partFile;
    partFile.open(fileName);
    tree.bodies_pos.d2h();
                                                                                                                    
    for(int i=0; i < tree.n; i++)                                                                                     
    {                                                                                                                 
      float4  pos =  tree.bodies_pos[i];                                                                              
      partFile << i << "\t" << pos.x << "\t" << pos.y << "\t" << pos.z << endl;                                                        
    }                                                                                                                 
    partFile.close(); 
    

  
  #endif
 
}
