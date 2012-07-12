#ifndef _OCTREE_H_
#define _OCTREE_H_


#define USE_CUDA

#ifdef USE_CUDA
  #include "my_cuda.h"
#else
  #include "my_ocl.h"
#endif


#include "node_specs.h"
#include <cmath>
#include <algorithm>
#include <sys/time.h>
#include <iostream>
#include <fstream>


using namespace std;

typedef float real;
typedef unsigned int uint;

#define NBLOCK_REDUCE     256
#define NBLOCK_BOUNDARY   120
#define NTHREAD_BOUNDARY  256
#define NBLOCK_PREFIX     512           //At the moment only used during memory alloc


inline int host_float_as_int(float val)
{
      union{float f; int i;} u; //__float_as_int
      u.f           = val;
      return u.i;
}

inline float host_int_as_float(int val)
{
      union{int i; float f;} itof; //__int_as_float
      itof.i           = val;
      return itof.f;
}


//For the sorting routines
typedef struct setupParams {
  int jobs;                     //Minimal number of jobs for each 'processor'
  int blocksWithExtraJobs;      //Some ' processors'  do one extra job all with bid < bWEJ
  int extraElements;            //The elements that didn't fit completely
  int extraOffset;              //Start of the extra elements

} setupParams;


//Structure and properties of a tree
class tree_structure
{
  private:
    my_dev::context *devContext;        //Pointer so destructor is only called once  

  public:
    int n;                                //Number of particles in the tree
    int n_leafs;                          //Number of leafs in the tree
    int n_nodes;                          //Total number of nodes in the tree (including leafs)
    int n_groups;                         //Number of groups
    int n_levels;                         //Depth of the tree
    bool needToReorder;			//Set to true if SetN is called so we know we need to change particle order

   
    my_dev::dev_mem<int>   activePartList;
    
    my_dev::dev_mem<int>   oriParticleOrder;         //Used to restore original particle order
    
    my_dev::dev_mem<uint2> level_list;    //List containing the start and end positions of each level
    my_dev::dev_mem<uint>  n_children;
    my_dev::dev_mem<uint2> node_bodies;
    my_dev::dev_mem<uint>  leafNodeIdx;    //First n_leaf items represent indices of leafs
                                           //remaining (n_nodes-n_leafs) are indices of non-leafs
    my_dev::dev_mem<uint2>  group_list;     //The group to particle relation
    my_dev::dev_mem<uint>  node_level_list; //List containing start and end idxs in (leafNode idx) for each level
    my_dev::dev_mem<uint>  body2group_list; //Contains per particle to which group it belongs TESTING TODO

    //Variables used for properties
    my_dev::dev_mem<real4>  multipole;      //Array storing the properties for each node (mass, mono, quad pole)
//     my_dev::dev_mem<real4>  nodeLowerBounds; //Lower bounds used for scaling? TODO
//     my_dev::dev_mem<real4>  nodeUpperBounds; //Upper bounds used for scaling? TODO

//     my_dev::dev_mem<uint>  ngb;                 //Compacted list of active groups
    my_dev::dev_mem<int2>  interactions;        //Counts the number of interactions, mainly for debugging and performance

    //BH Opening criteria
    my_dev::dev_mem<float4> boxSizeInfo;
    my_dev::dev_mem<float4> boxCenterInfo;

    //     my_dev::dev_mem<float4> groupSizeInfo;
    //     my_dev::dev_mem<float4> groupCenterInfo;

    //Combined buffers:
    /*
     * TODO this list is not complete anymore
    Buffer1 used during: Sorting, Tree-construction, and Tree-traverse:
    Sorting:
      - SrcValues, Output, simpleKeys, permutation, output32b, valuesOutput
    Tree-construction:
      - ValidList, compactList
    Tree-traverse:
      - Interactions, NGB, Active_partliist
    */

    my_dev::dev_mem<uint> generalBuffer1;

    real4 corner;                         //Corner of tree-structure
    real  domain_fac;                     //Domain_fac of tree-structure

  tree_structure(){ n = 0;}

  tree_structure(my_dev::context &context)
  {
    n = 0;
    devContext = &context;
    setMemoryContexts();
    needToReorder = true;
  }
  void setContext(my_dev::context &context)
  {
    devContext = &context;
    setMemoryContexts();
    needToReorder = true;
  }

  void setN(int particles)
  {
    n = particles;
    needToReorder = true;
  }

  void setMemoryContexts()
  {
    n_children.setContext(*devContext);
    node_bodies.setContext(*devContext);
    leafNodeIdx.setContext(*devContext);
    group_list.setContext(*devContext);
    node_level_list.setContext(*devContext);
    multipole.setContext(*devContext);

    body2group_list.setContext(*devContext);

    activePartList.setContext(*devContext);

    oriParticleOrder.setContext(*devContext);
    level_list.setContext(*devContext);
//     ngb.setContext(*devContext);
    interactions.setContext(*devContext);

    //BH Opening
    boxSizeInfo.setContext(*devContext);
    boxCenterInfo.setContext(*devContext);

    //     groupSizeInfo.setContext(*devContext);
    //     groupCenterInfo.setContext(*devContext);

    //General buffer
    generalBuffer1.setContext(*devContext);
  }
  

  my_dev::context getContext()
  {
    return *devContext;
  }
};



class octree {
protected:
  int devID;
  
  char *execPath;
  char *src_directory;
  
  //Device configuration
  int nMultiProcessors;
  int nBlocksForTreeWalk;

   //Simulation properties
    float eps2;
    float inv_theta;
    int   dt_limit;
    float eta;
    float theta;
    float rsearch_sq;    // modified by M.I.

    //Sim stats
    double Ekin, Ekin0, Ekin1;
    double Epot, Epot0, Epot1;
    double Etot, Etot0, Etot1;

 
  my_dev::context devContext;
  bool devContext_flag;
  
  //Streams
  my_dev::dev_stream *execStream;

  // scan & sort algorithm

  my_dev::kernel  compactCount, exScanBlock, compactMove, splitMove;
  my_dev::kernel  sortCount, sortMove;
  my_dev::kernel  extractInt, fillSequence, reOrderKeysValues;
  my_dev::kernel  convertKey64to96, extractKeyAndPerm;
  my_dev::kernel  dataReorderR4;
  my_dev::kernel  dataReorderF2;
  my_dev::kernel  dataReorderI1;

  // tree construction kernels
  my_dev::kernel  build_key_list;
  my_dev::kernel  build_valid_list;
  my_dev::kernel  build_nodes;
  my_dev::kernel  link_tree;
  my_dev::kernel  define_groups;
  my_dev::kernel  build_level_list;
  my_dev::kernel  store_groups;
  my_dev::kernel  expand_leaflist;

  my_dev::kernel  boundaryReduction;
  my_dev::kernel  boundaryReductionGroups;
  my_dev::kernel  build_body2group_list;
  my_dev::kernel  build_phkey_list;


  // tree properties kernels
  my_dev::kernel  propsNonLeaf, propsLeaf, propsScaling;
  my_dev::kernel  propsNonLeafD, propsLeafD, propsScalingD;

  my_dev::kernel  copyNodeDataToGroupData;
  my_dev::kernel  setActiveGrps;

  //Iteraction kernels
  my_dev::kernel approxGrav;

  ///////////
  
private:
  void createDeviceQueue();  

public:
  
  
   my_dev::context & getDevContext()
   {
     return devContext;
   }
   
   void set_src_directory(string src_dir);

   //Memory used in the whole system, not depending on a certain number of particles
   my_dev::dev_mem<float> tnext;
   my_dev::dev_mem<uint>  nactive;
   //General memory buffers
   my_dev::dev_mem<float3>  devMemRMIN;
   my_dev::dev_mem<float3>  devMemRMAX;

   my_dev::dev_mem<uint> devMemCounts;
   my_dev::dev_mem<uint> devMemCountsx;

   tree_structure localTree;

//   void set_context(bool disable_timing = false);
  void set_context();   
  void set_context(std::ofstream &log, bool disable_timing);
  

//   void set_context2();

  int getAllignmentOffset(int n);
  int getTextureAllignmentOffset(int n, int size);
  
  

  //GPU kernels and functions
  void load_kernels();
  
  //Sorting, compacting, splitting functions
  void gpuCompact(my_dev::context&, my_dev::dev_mem<uint> &srcValues,
                  my_dev::dev_mem<uint> &output, int N, int *validCount);
  void gpuSplit(my_dev::context&, my_dev::dev_mem<uint> &srcValues,
                  my_dev::dev_mem<uint> &output, int N, int *validCount);
  void gpuSort(my_dev::context &devContext,
                       my_dev::dev_mem<uint4> &srcValues,
                       my_dev::dev_mem<uint4> &output,
                       my_dev::dev_mem<uint4> &buffer,
                       int N, int numberOfBits, int subItems,
                       tree_structure &tree);

  void gpuSort_32b(my_dev::context &devContext,
                    my_dev::dev_mem<uint> &srcKeys,     my_dev::dev_mem<uint> &srcValues,
                    my_dev::dev_mem<int>  &keysOutput,  my_dev::dev_mem<uint> &keysAPing,
                    my_dev::dev_mem<uint> &valuesOutput,my_dev::dev_mem<uint> &valuesAPing,
                   int N, int numberOfBits);


  void desort_bodies(tree_structure &tree);
  void sort_bodies(tree_structure &tree, my_dev::dev_mem<float4>  &bodies_pos,
                   my_dev::dev_mem<uint>  &sortPermutation, int n_bodies);
  
  void getCorner(my_dev::dev_mem<float4>  &bodies_pos, int n_bodies, float4 &corner, float &domain_fac);
  
  void compute_keys(my_dev::dev_mem<float4>  &bodies_pos, 
                    my_dev::dev_mem<uint4>  &bodies_key,
                    int n_bodies, float4 &corner, float &domain_fac);

  
  void reorder_dataR4(my_dev::dev_mem<real4>  &data, my_dev::dev_mem<uint>  &sortPermutation, int n_items);
  void reorder_dataR2(my_dev::dev_mem<real2>  &data, my_dev::dev_mem<uint>  &sortPermutation, int n_items);
  void reorder_dataI1(my_dev::dev_mem<int>    &data, my_dev::dev_mem<uint>  &sortPermutation, int n_items);  
  
  
  void getBoundaries(my_dev::dev_mem<float4>  &bodies_pos, int n_bodies, real4 &r_min, real4 &r_max);
  
  void getBoundariesGroups(tree_structure &tree, real4 &r_min, real4 &r_max);  

  //ALlocate memory
  void allocateParticleMemory(tree_structure &tree);
  void allocateTreePropMemory(tree_structure &tree);
  void reallocateParticleMemory(tree_structure &tree);
  
  void allocateInternalBuffers();
  void allocateParticleSpecificBuffers(int n_bodies);
  void allocateNodeSpecificBuffers(int n_nodes);
  void allocateGroupSpecificBuffers(int n_bodies, int n_groups);

  void build (tree_structure &tree, my_dev::dev_mem<float4>  &bodies_pos, int n_bodies);
  void compute_properties (tree_structure &tree,  my_dev::dev_mem<float4>  &bodies_pos, int n_bodies);
  
  void createGroups(tree_structure &tree, my_dev::dev_mem<float4>  &bodies_pos, int n_bodies);
    
  
  void setActiveGrpsFunc(tree_structure &tree);

  void iterate();

  //Subfunctions of iterate, should probally be private TODO
//   void approximate_gravity(tree_structure &tree);
  /*void approximate_gravity(tree_structure &tree,
                                 my_dev::dev_mem<float4>  &bodies_pos,  //Bodies that are part of the tree-structure
                                 int n_bodies,                          //Number of bodies that are part of tree-structure
                                 my_dev::dev_mem<float4>  &group_bodies_pos,  //Bodies that are part of the groups
                                 int n_groupBodies,                     //Number of bodies that are part of the groups
                                 my_dev::dev_mem<float4>  &bodies_acc);  
  */
    /*  
  void  approximate_gravity(tree_structure &tree,
                                 my_dev::dev_mem<float4>  &j_bodies_pos,  //Bodies that are part of the tree-structure
                                 my_dev::dev_mem<float4>  &i_bodies_pos,  //Bodies that are part of the groups
                                 int n_groupBodies,                     //Number of bodies that are part of the groups
                                 my_dev::dev_mem<real4>  &i_bodies_acc,
                                 my_dev::dev_mem<real>  &i_bodies_ds2,
                                 my_dev::dev_mem<int>    &i_bodies_ngb);
    */
    // by M.I.
    void  approximate_gravity(tree_structure &tree,
			      my_dev::dev_mem<float4>  &j_bodies_pos,  //Bodies that are part of the tree-structure
			      my_dev::dev_mem<float4>  &i_bodies_pos,  //Bodies that are part of the groups
			      int n_groupBodies,                     //Number of bodies that are part of the groups
			      my_dev::dev_mem<real4>  &i_bodies_acc,
			      my_dev::dev_mem<real>  &i_bodies_ds2,
			      my_dev::dev_mem<int>    &i_bodies_ngb,
			      my_dev::dev_mem<int>    &i_bodies_Nngb);
  
  
  
  
  double compute_energies(tree_structure &tree);

  //Library interface functions  
  void  setEps(float eps);
  float getEps();
  void  setDt(float dt);
  float getDt();
  void  setTheta(float theta);
  float getTheta();
  void  setTEnd(float tEnd);
  float getTEnd();
  void  setTime(float);
  float getTime();
  float getPot();
  float getKin();

  //End library functions

    octree(char** argv, 
	   const int device = 0, 
	   const float _theta = 0.75, 
	   const float eps = 0.05) {
	//octree(char **argv, const int device = 0, const float _theta = 0.75, const float eps = 0.05) {

	devContext_flag = false;

	src_directory = NULL;

	if(argv != NULL)
	    execPath = argv[0];

	//First init mpi
	//     int argc = 0;
	//     mpiInit(argc, argv, procId, nProcs);

	devID = device;
	//devID = procId % nProcs;
	//     devID = 1;
	//    devID = procId % getNumberOfCUDADevices();
	//      devID = 1;
	//     cerr << "Preset device : "  << devID << "\t" << device<<endl;

	//TODO!
	inv_theta   = 1.0/_theta;
	eps2        = eps*eps;
	eta         = 0.02;
	theta       = _theta;
	rsearch_sq = 99999999.9; // by M.I.

	execStream = NULL;    
    
	getNumberOfCUDADevices();  //Stop compiler warnings
    }

    octree(char** argv, 
	   const int device, 
	   const float _theta, 
	   const float eps,
	   const float _rsearch_sq) {

	devContext_flag = false;

	src_directory = NULL;

	if(argv != NULL)
	    execPath = argv[0];

	//First init mpi
	//     int argc = 0;
	//     mpiInit(argc, argv, procId, nProcs);

	devID = device;
	//devID = procId % nProcs;
	//     devID = 1;
	//    devID = procId % getNumberOfCUDADevices();
	//      devID = 1;
	//     cerr << "Preset device : "  << devID << "\t" << device<<endl;
	
	//TODO!
	inv_theta   = 1.0/_theta;
	eps2        = eps*eps;
	eta         = 0.02;
	theta       = _theta;
	rsearch_sq = _rsearch_sq;

	execStream = NULL;    
	
	getNumberOfCUDADevices();  //Stop compiler warnings
    }
    
  ~octree() {    
  };
};


#endif // _OCTREE_H_
