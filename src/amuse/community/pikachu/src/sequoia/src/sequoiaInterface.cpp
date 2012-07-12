#include "../include/sequoiaInterface.h"
#include "../include/octree.h"


//The memory counters
long long my_dev::base_mem::currentMemUsage;
long long my_dev::base_mem::maxMemUsage;

//The Bonsai Octree class

octree *sequoia;

bool initFlag = false;
ofstream logFile;
  
  
//extern "C" {
my_dev::context & sequoia_init(char** argv, 
			       int device, 
			       const float _theta, 
			       const float eps){  
    sequoia = new octree(argv, device, _theta, eps);
   
    //Used for profiler
    char *gpu_prof_log;
    gpu_prof_log=getenv("CUDA_PROFILE_LOG");
    if(gpu_prof_log){
	char tmp[50];
	sprintf(tmp,"process%d_%s",0,gpu_prof_log);
	setenv("CUDA_PROFILE_LOG",tmp,1);
    }

    string logFileName    = "gpuLog.log";
    logFile.open(logFileName.c_str());
    
    sequoia->set_context(logFile, false); //Do logging to file and enable timing (false = enabled)
  
    //Load the device kernels
    sequoia->load_kernels();
 
    initFlag = true;
  
    return sequoia->getDevContext();
}


// by M.I.
my_dev::context & sequoia_init(char** argv, 
			       int device, 
			       const float _theta, 
			       const float eps,
			       const float rsearch2){  
    sequoia = new octree(argv, device, _theta, eps, rsearch2);
   
    //Used for profiler
    char *gpu_prof_log;
    gpu_prof_log=getenv("CUDA_PROFILE_LOG");
    if(gpu_prof_log){
	char tmp[50];
	sprintf(tmp,"process%d_%s",0,gpu_prof_log);
	setenv("CUDA_PROFILE_LOG",tmp,1);
    }
    
    string logFileName    = "gpuLog.log";
    logFile.open(logFileName.c_str());
    
    sequoia->set_context(logFile, false); //Do logging to file and enable timing (false = enabled)
  
    //Load the device kernels
    sequoia->load_kernels();
    
    initFlag = true;
  
    return sequoia->getDevContext();
}


int sequoia_cleanup()
{
  assert(initFlag);
  
  
  delete sequoia;
  sequoia = NULL;
  
  
  
  initFlag = false;
  return 0;
}


int sequoia_sortBodies(my_dev::dev_mem<real4>  &bodies_pos, my_dev::dev_mem<uint> &permutation, int n_bodies)
{
  sequoia->sort_bodies(sequoia->localTree, bodies_pos, permutation, n_bodies);  
  return 0; 
}

int sequoia_reorderReal4(my_dev::dev_mem<real4>  &data, my_dev::dev_mem<uint> &permutation, int n_items)
{
  sequoia->reorder_dataR4(data, permutation, n_items);  
  return 0; 
}

int sequoia_reorderReal2(my_dev::dev_mem<real2>  &data, my_dev::dev_mem<uint> &permutation, int n_items)
{
  sequoia->reorder_dataR2(data, permutation, n_items);  
  return 0; 
}

int sequoia_reorderInt1(my_dev::dev_mem<int>  &data, my_dev::dev_mem<uint> &permutation, int n_items)
{
  sequoia->reorder_dataI1(data, permutation, n_items);  
  return 0; 
}

int sequoia_buildTreeStructure(my_dev::dev_mem<real4>  &bodies_pos, int n_bodies)
{
  sequoia->build(sequoia->localTree, bodies_pos, n_bodies);
  sequoia->compute_properties(sequoia->localTree, bodies_pos, n_bodies);
  return 0;
}

int sequoia_computeTreeProperties(my_dev::dev_mem<real4>  &bodies_pos, int n_bodies)
{
    sequoia->compute_properties(sequoia->localTree, bodies_pos, n_bodies);
    return 0;
}

int sequoia_createGroups(my_dev::dev_mem<real4> &bodies_pos, int n_bodies)
{
  sequoia->createGroups(sequoia->localTree, bodies_pos, n_bodies);
  return 0;
}

// by M.I.
int sequoia_computeGravity(my_dev::dev_mem<real4> &j_bodies_pos, 
			   my_dev::dev_mem<real4> &i_bodies_pos, 
			   my_dev::dev_mem<real4> &i_bodies_acc, 
			   my_dev::dev_mem<real>  &i_bodies_ds2, 
			   my_dev::dev_mem<int>   &i_bodies_ngb, 
			   my_dev::dev_mem<int>   &i_bodies_Nngb, 
			   int                     n_i_bodies)
{
  
  sequoia->approximate_gravity(sequoia->localTree, j_bodies_pos, i_bodies_pos,
			       n_i_bodies, i_bodies_acc, i_bodies_ds2,
			       i_bodies_ngb, i_bodies_Nngb);
  
  return 0;
}

/*
int get_nodeID(uint* leafsIdxs, uint id){
    return leafsIdxs[id];
}
*/

// by M.I.
int sequoia_setParticlesAndGetGravity_firsthalf(my_dev::dev_mem<real4> &j_bodies_pos,
						my_dev::dev_mem<int>   &j_bodies_ids,
						int                    n_j_bodies,
						bool                   sortJBodies,
						uint *& out_leafNodeIdx,
						uint2 *& out_node_bodies,
						uint *& out_n_children,
						uint2 *& out_level_list,
						uint *& out_node_level_list,
						real4 *& out_multipole,
						float4 *& out_boxSizeInfo,
						float4 *& out_boxCenterInfo,
						int &out_n_leafs,
						int &out_n_nodes,
						int &out_n_levels){
    assert(initFlag);
  
    //If required sort J particles
    if(sortJBodies){
	my_dev::dev_mem<uint> permutation(sequoia->getDevContext(), n_j_bodies);
	sequoia_sortBodies  (j_bodies_pos, permutation, n_j_bodies);
	sequoia_reorderReal4(j_bodies_pos, permutation, n_j_bodies);
	sequoia_reorderInt1 (j_bodies_ids, permutation, n_j_bodies);
    }

    //Build the tree-structure
    sequoia_buildTreeStructure(j_bodies_pos, n_j_bodies);

    sequoia->localTree.leafNodeIdx.d2h();
    sequoia->localTree.node_bodies.d2h();
    sequoia->localTree.n_children.d2h();
    sequoia->localTree.level_list.d2h();
    sequoia->localTree.node_level_list.d2h();
    sequoia->localTree.multipole.d2h();
    sequoia->localTree.boxSizeInfo.d2h();
    sequoia->localTree.boxCenterInfo.d2h();
    j_bodies_pos.d2h();
    j_bodies_ids.d2h();


    out_leafNodeIdx = sequoia->localTree.leafNodeIdx.h();
    out_node_bodies = sequoia->localTree.node_bodies.h();
    out_n_children = sequoia->localTree.n_children.h();
    out_level_list = sequoia->localTree.level_list.h();
    out_node_level_list = sequoia->localTree.node_level_list.h();
    out_multipole = sequoia->localTree.multipole.h();
    out_boxSizeInfo = sequoia->localTree.boxSizeInfo.h();
    out_boxCenterInfo = sequoia->localTree.boxCenterInfo.h();
    out_n_leafs = sequoia->localTree.n_leafs;
    out_n_nodes = sequoia->localTree.n_nodes;
    out_n_levels = sequoia->localTree.n_levels;
  
  return 0;
}

int sequoia_setParticlesAndGetGravity_firsthalf_for_neighbour_search(my_dev::dev_mem<real4> &j_bodies_pos,
								     my_dev::dev_mem<int>   &j_bodies_ids,
								     int                    n_j_bodies,
								     bool                   sortJBodies,
								     uint *& out_leafNodeIdx,
								     uint2 *& out_node_bodies,
								     uint *& out_n_children,
								     float4 *& out_boxSizeInfo,
								     float4 *& out_boxCenterInfo,
								     int &out_n_leafs,
								     int &out_n_nodes){

    assert(initFlag);
  
    //If required sort J particles
    if(sortJBodies){
	my_dev::dev_mem<uint> permutation(sequoia->getDevContext(), n_j_bodies);
	sequoia_sortBodies  (j_bodies_pos, permutation, n_j_bodies);
	sequoia_reorderReal4(j_bodies_pos, permutation, n_j_bodies);
	sequoia_reorderInt1 (j_bodies_ids, permutation, n_j_bodies);
    }

    //Build the tree-structure
    sequoia_buildTreeStructure(j_bodies_pos, n_j_bodies);


    sequoia->localTree.leafNodeIdx.d2h();
    sequoia->localTree.node_bodies.d2h();
    sequoia->localTree.n_children.d2h();
    sequoia->localTree.boxSizeInfo.d2h();
    sequoia->localTree.boxCenterInfo.d2h();
    //j_bodies_pos.d2h();
    //j_bodies_ids.d2h();
    j_bodies_pos.d2h(n_j_bodies);
    j_bodies_ids.d2h(n_j_bodies);

    out_leafNodeIdx = sequoia->localTree.leafNodeIdx.h();
    out_node_bodies = sequoia->localTree.node_bodies.h();
    out_n_children = sequoia->localTree.n_children.h();
    out_boxSizeInfo = sequoia->localTree.boxSizeInfo.h();
    out_boxCenterInfo = sequoia->localTree.boxCenterInfo.h();
    out_n_leafs = sequoia->localTree.n_leafs;
    out_n_nodes = sequoia->localTree.n_nodes;
  
  return 0;
}

int sequoia_setParticlesAndGetGravity_lasthalf(my_dev::dev_mem<real4> &j_bodies_pos,
					       my_dev::dev_mem<real4> &i_bodies_pos,
					       my_dev::dev_mem<int>   &i_bodies_ids,
					       int                    n_i_bodies, 
					       bool                   sortIBodies,
					       my_dev::dev_mem<real4> &i_bodies_acc,
					       my_dev::dev_mem<real>  &i_bodies_ds2,
					       my_dev::dev_mem<int>   &i_bodies_ngb,
					       my_dev::dev_mem<int>   &i_bodies_Nngb){
    if(sortIBodies){
	my_dev::dev_mem<uint> permutation(sequoia->getDevContext(), n_i_bodies);
        sequoia_sortBodies  (i_bodies_pos, permutation, n_i_bodies);
	sequoia_reorderReal4(i_bodies_pos, permutation, n_i_bodies);
	sequoia_reorderInt1 (i_bodies_ids, permutation, n_i_bodies);
    }

    //Create the groups that walk the tree
    sequoia_createGroups(i_bodies_pos, n_i_bodies);
   
    //Finally compute the gravity and get the nearest neighbour + distance 
    sequoia_computeGravity(j_bodies_pos, i_bodies_pos, i_bodies_acc, 
			   i_bodies_ds2, i_bodies_ngb, i_bodies_Nngb, 
			   n_i_bodies);  

}


int sequoia_setParticlesAndGetGravity(my_dev::dev_mem<real4> &j_bodies_pos,      //Positions J-particles
				      my_dev::dev_mem<int>   &j_bodies_ids,      //Particle IDs J-particles
				      int                    n_j_bodies,         //Number of J-particles
				      my_dev::dev_mem<real4> &i_bodies_pos,      //Positions I-particles 
				      my_dev::dev_mem<int>   &i_bodies_ids,      //Particle IDs J-particles
				      int                    n_i_bodies,         //Number of I-particles
				      bool                   sortJBodies,        //Do we need to sort J-particles?
				      bool                   sortIBodies,        //Do we need to sort I-particles?           
				      my_dev::dev_mem<real4> &i_bodies_acc,      //OUT  Accelerations for I-particles
				      my_dev::dev_mem<real>  &i_bodies_ds2,      //OUT  min distance squared for I-particles
				      my_dev::dev_mem<int>   &i_bodies_ngb,      //OUT  J-ID of the nearest neighbour for I-particles
				      my_dev::dev_mem<int>   &i_bodies_Nngb)     //OUT  the number of the nearest neighbour for I-particles
{
    assert(initFlag);
  
    //If required sort J particles
    if(sortJBodies){
	my_dev::dev_mem<uint> permutation(sequoia->getDevContext(), n_j_bodies);
	sequoia_sortBodies  (j_bodies_pos, permutation, n_j_bodies);
	sequoia_reorderReal4(j_bodies_pos, permutation, n_j_bodies);
	sequoia_reorderInt1 (j_bodies_ids, permutation, n_j_bodies);
    }
    //If required sort I particles
    if(sortIBodies){
	my_dev::dev_mem<uint> permutation(sequoia->getDevContext(), n_i_bodies);
        sequoia_sortBodies  (i_bodies_pos, permutation, n_i_bodies);
	sequoia_reorderReal4(i_bodies_pos, permutation, n_i_bodies);
	sequoia_reorderInt1 (i_bodies_ids, permutation, n_i_bodies);
    }

    //Build the tree-structure
    sequoia_buildTreeStructure(j_bodies_pos, n_j_bodies); // split 2 functions

    //Create the groups that walk the tree
    sequoia_createGroups(i_bodies_pos, n_i_bodies);
   
    //Finally compute the gravity and get the nearest neighbour + distance 
    sequoia_computeGravity(j_bodies_pos, i_bodies_pos, i_bodies_acc, 
			   i_bodies_ds2, i_bodies_ngb, i_bodies_Nngb, 
			   n_i_bodies);  
    return 0;
}



