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

#include "mpi.h"

#define PRINT_MPI_DEBUG


using namespace std;

typedef float real;
typedef unsigned int uint;

#define NBLOCK_REDUCE     256

#define NBLOCK_BOUNDARY   120
#define NTHREAD_BOUNDARY  256

#define NBLOCK_PREFIX     512           //At the moment only used during memory alloc


struct morton_struct {
  uint2 key;
  int   value;
};


inline int cmp_uint2(const uint2 a, const uint2 b) {
  if      (a.x < b.x) return -1;
  else if (a.x > b.x) return +1;
  else {
    if       (a.y < b.y) return -1;
    else  if (a.y > b.y) return +1;
    return 0;
  }
}

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

typedef struct setupParams {
  int jobs;                     //Minimal number of jobs for each 'processor'
  int blocksWithExtraJobs;      //Some ' processors'  do one extra job all with bid < bWEJ
  int extraElements;            //The elements that didn't fit completely
  int extraOffset;              //Start of the extra elements

} setupParams;

typedef struct bodyStruct
{
  real4 pos;
  real4 vel;
  real4 acc0;
  real4 acc1;
  float2 time;
  int   id;
} bodyStruct;


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
    my_dev::dev_mem<real4> bodies_pos;    //The particles positions
    my_dev::dev_mem<uint4> bodies_key;    //The particles keys
    my_dev::dev_mem<real4> bodies_vel;    //Velocities
    my_dev::dev_mem<real4> bodies_acc0;    //Acceleration
    my_dev::dev_mem<real4> bodies_acc1;    //Acceleration
    my_dev::dev_mem<float2> bodies_time;  //The timestep details (.x=tb, .y=te
    my_dev::dev_mem<int>   bodies_ids;
    my_dev::dev_mem<int>   oriParticleOrder;         //Used to restore original particle order
    
    my_dev::dev_mem<real4> bodies_Ppos;    //Predicted position
    my_dev::dev_mem<real4> bodies_Pvel;    //Predicted velocity

    my_dev::dev_mem<uint2> level_list;    //List containing the start and end positions of each level
    my_dev::dev_mem<uint4> node_key;
    my_dev::dev_mem<uint>  n_children;
    my_dev::dev_mem<uint2> node_bodies;
    my_dev::dev_mem<uint>  leafNodeIdx;    //First n_leaf items represent indices of leafs
                                           //remaining (n_nodes-n_leafs) are indices of non-leafs
    my_dev::dev_mem<uint>  group_list;     //The id's of nodes that form a group
    my_dev::dev_mem<uint>  node_level_list; //List containing start and end idxs in (leafNode idx) for each level
    my_dev::dev_mem<uint>  body2group_list; //Contains per particle to which group it belongs TESTING TODO
//     my_dev::dev_mem<uint4> group_data;      //Contains copies of the node_data of nodes that form the groups, has to
                                            //be copies since for parallel version the node_data can contain data
                                            //from the remove tree instead of the local tree
    my_dev::dev_mem<uint2>  group_list_test;     //The group to particle relation

    //Variables used for properties
    my_dev::dev_mem<real4>  multipole;      //Array storing the properties for each node (mass, mono, quad pole)
    my_dev::dev_mem<real4>  nodeLowerBounds; //Lower bounds used for scaling? TODO
    my_dev::dev_mem<real4>  nodeUpperBounds; //Upper bounds used for scaling? TODO

    //Variables used for iteration
    int n_active_groups;
    int n_active_particles;

    my_dev::dev_mem<uint>  activeGrpList;       //Non-compacted list of active grps
    my_dev::dev_mem<uint>  active_group_list;   //Compacted list of active groups
    my_dev::dev_mem<uint>  activePartlist;      //Compacted list of active groups
    my_dev::dev_mem<uint>  ngb;                 //Compacted list of active groups

    my_dev::dev_mem<int2>  interactions;        //Counts the number of interactions, mainly for debugging and performance

    //BH Opening criteria
    my_dev::dev_mem<float4> boxSizeInfo;
    my_dev::dev_mem<float4> groupSizeInfo;

    my_dev::dev_mem<float4> boxCenterInfo;
    my_dev::dev_mem<float4> groupCenterInfo;

    //Combined buffers:
    /*
    Buffer1 used during: Sorting, Tree-construction, and Tree-traverse:
    Sorting:
      - SrcValues, Output, simpleKeys, permutation, output32b, valuesOutput
    Tree-construction:
      - ValidList, compactList
    Tree-traverse:
      - Interactions, NGB, Active_partliist
    */

    my_dev::dev_mem<uint> generalBuffer1;


    my_dev::dev_mem<float4> fullRemoteTest;

    uint4 remoteTreeStruct;

    real4 corner;                         //Corner of tree-structure
    real  domain_fac;                     //Domain_fac of tree-structure

  tree_structure(){}

  tree_structure(my_dev::context &context)
  {
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
    bodies_pos.setContext(*devContext);
    bodies_key.setContext(*devContext);
    node_key.setContext(*devContext);
    n_children.setContext(*devContext);
    node_bodies.setContext(*devContext);
    leafNodeIdx.setContext(*devContext);
    group_list.setContext(*devContext);
    node_level_list.setContext(*devContext);
    multipole.setContext(*devContext);
    nodeLowerBounds.setContext(*devContext);
    nodeUpperBounds.setContext(*devContext);

    body2group_list.setContext(*devContext);

    bodies_vel.setContext(*devContext);
    bodies_acc0.setContext(*devContext);
    bodies_acc1.setContext(*devContext);
    bodies_time.setContext(*devContext);
    bodies_ids.setContext(*devContext);

    bodies_Ppos.setContext(*devContext);
    bodies_Pvel.setContext(*devContext);

    oriParticleOrder.setContext(*devContext);
//     group_data.setContext(*devContext);
    activeGrpList.setContext(*devContext);
    active_group_list.setContext(*devContext);
    level_list.setContext(*devContext);
    activePartlist.setContext(*devContext);
    ngb.setContext(*devContext);
    interactions.setContext(*devContext);

    group_list_test.setContext(*devContext);

    //BH Opening
    boxSizeInfo.setContext(*devContext);
    groupSizeInfo.setContext(*devContext);
    boxCenterInfo.setContext(*devContext);
    groupCenterInfo.setContext(*devContext);

    //General buffer
    generalBuffer1.setContext(*devContext);

    fullRemoteTest.setContext(*devContext);
  }

  my_dev::context getContext()
  {
    return *devContext;
  }
};



class octree {
protected:
  int devID;

  char *src_directory;
  
   //Simulation properties
  int   iter;
  float t_current, t_previous;
  int   snapshotIter;
  string snapshotFile;
  int nextSnapTime;

  int   NTotal, NFirst, NSecond, NThird, snapShotAdd;
  float killDistance;
  float removeDistance;

  float eps2;
  float inv_theta;
  int   dt_limit;
  float eta;
  float timeStep;
  int   tEnd;
  float theta;

  //Sim stats
  double Ekin, Ekin0, Ekin1;
  double Epot, Epot0, Epot1;
  double Etot, Etot0, Etot1;

  bool store_energy_flag;
  double tinit;

  // OpenCL context

  my_dev::context devContext;
  bool devContext_flag;

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
  my_dev::kernel  build_body2group_list;

  my_dev::kernel  build_phkey_list;


  // tree properties kernels
  my_dev::kernel  propsNonLeaf, propsLeaf, propsScaling;

  my_dev::kernel  propsNonLeafD, propsLeafD, propsScalingD;

  my_dev::kernel  copyNodeDataToGroupData;
  my_dev::kernel  setActiveGrps;

  //Iteraction kernels
  my_dev::kernel getTNext;
  my_dev::kernel predictParticles;
  my_dev::kernel getNActive;
  my_dev::kernel approxGrav;
  my_dev::kernel correctParticles;
  my_dev::kernel computeDt;
  my_dev::kernel computeEnergy;

  my_dev::kernel distanceCheck;
  my_dev::kernel approxGravLET;

  ///////////////////////

  /////////////////
  void write_dumbp_snapshot(real4 *bodyPositions, real4 *bodyVelocities, int *ids, int n, string fileName);

  void   to_binary(int);
  void   to_binary(uint2);
  uint2  dilate3(int);
  int    undilate3(uint2);

  uint2  get_key(int3);
  int3   get_crd(uint2);
  real4  get_pos(uint2, real, tree_structure&);

  uint2  get_mask(int);
  uint2  get_imask(uint2);

  int find_key(uint2, vector<uint2>&,         int, int);
  int find_key(uint2, vector<morton_struct>&, int, int);

  ///////////

public:
   void set_src_directory(string src_dir);
   
   double get_time();

   void write_dumbp_snapshot_parallel(real4 *bodyPositions, real4 *bodyVelocities, int* bodyIds, int n, string fileName) ;
   void write_dumbp_snapshot_parallel_tipsy(real4 *bodyPositions, real4 *bodyVelocities, int* bodyIds, int n, string fileName,
                                            int NCombTotal, int NCombFirst, int NCombSecond, int NCombThird);

   //Memory used in the whole system, not depending on a certain number of particles
   my_dev::dev_mem<float> tnext;
   my_dev::dev_mem<uint>  nactive;
   //General memory buffers
   my_dev::dev_mem<float3>  devMemRMIN;
   my_dev::dev_mem<float3>  devMemRMAX;

   my_dev::dev_mem<uint> devMemCounts;
   my_dev::dev_mem<uint> devMemCountsx;

   tree_structure localTree;
   tree_structure remoteTree;


  //TODO jb made these functions public for testing
  void set_context(bool disable_timing = false);
  void set_context(std::ofstream &log, bool disable_timing = false);
  void set_context2();

  int getAllignmentOffset(int n);

  //GPU kernels and functions
  void load_kernels();
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
//                     my_dev::dev_mem<int>  &keysOutput,  my_dev::dev_mem<int> &keysAPing,
                    my_dev::dev_mem<int>  &keysOutput,  my_dev::dev_mem<uint> &keysAPing,
                    my_dev::dev_mem<uint> &valuesOutput,my_dev::dev_mem<uint> &valuesAPing,
                    //my_dev::dev_mem<uint> &count,
                   int N, int numberOfBits);
  //end TODO

//  void compute_keys();
//  void sort_keys();
  void desort_bodies(tree_structure &tree);
  void sort_bodies(tree_structure &tree);
  void getBoundaries(tree_structure &tree, real4 &r_min, real4 &r_max);

  void allocateParticleMemory(tree_structure &tree);
  void allocateTreePropMemory(tree_structure &tree);
  void reallocateParticleMemory(tree_structure &tree);

  void build(tree_structure &tree);
  void compute_properties (tree_structure &tree);
  void compute_properties_double(tree_structure &tree);
  void setActiveGrpsFunc(tree_structure &tree);

  void iterate();

  //Subfunctions of iterate, should probally be private TODO
  void predict(tree_structure &tree);
  void approximate_gravity(tree_structure &tree);
  void correct(tree_structure &tree);
  double compute_energies(tree_structure &tree);

  int  checkMergingDistance(tree_structure &tree, int iter, double dE);
  void checkRemovalDistance(tree_structure &tree);


  //Parallel version functions
  //Approximate for LET
  void approximate_gravity_let(tree_structure &tree, tree_structure &remoteTree);

  //Parallel version functions
  int procId, nProcs;   //Proccess ID in the mpi stack, number of processors in the commm world
  int nx, ny, nz;       //The division parameters, number of procs in the x,y and z axis
  int sampleFreq;       //Sample frequency for the division

  double4 *xlow, *xhigh;          //Domain boundaries for each processor
  double4 *cur_xlow, *cur_xhigh;  //The current boundaries for each processor, this can
                                  //be different from xlow and xhigh since the domain
                                  //is not updated during each step
  double4 *cur_xlow_xhigh;        //Combines the values cur_xlow and cur_xhigh in one array to reduce communication

  double4 *let_xlow, *let_xhigh;  //The current boundaries for each processor, this can
                                  //be different from xlow and xhigh since the domain
                                  //is not updated during each step

  real maxLocalEps;               //Contains the maximum local eps/softening value
                                  //will be stored in cur_xlow[i].w after exchange
  real4 rMinLocalTree;
  real4 rMaxLocalTree;



  //Functions
  void mpiInit(int argc,char *argv[], int &procId, int &nProcs);

  //Utility
  void mpiSync();
  int mpiGetRank();
  int mpiGetNProcs();
  void mpiRadiusFind(real4 &rmin, real4 &rmax);
  void AllSum(double &value);
  int  SumOnRootRank(int &value);

  //Main MPI functions

  //Functions for domain division
  void createORB();
  void determine_sample_freq(int numberOfParticles);
  void createDistribution(real4 *bodies, int n_bodies);
  void collect_sample_particles(real4 *bodies, int nbody,
                              int sample_freq, vector<real4> &sampleArray,
                              int &nsample, double &rmax);
  void determine_division(int np, vector<real4> &pos, int nx, int ny, int nz,
                          double rmax, double4 xlow[], double4 xhigh[]);
  void sortCoordinates(real4 *r, int lo, int up, int cid );
  void calculate_boxdim(int np, real4 pos[], int cid, int istart, int iend,
                      double rmax, double & xlow, double & xhigh);
  void updateDistribution(real4 *bodies, int n_bodies);

  int  exchange_particles_with_overflow_check(tree_structure &localTree);

  template<class T>
  int MP_exchange_particle_with_overflow_check(int ibox,
                                              T *source_buffer,
                                              vector<T> &recv_buffer,
                                              int firstloc,
                                              int nparticles,
                                              int isource,
                                              int &nsend,
                                              unsigned int &recvCount);

   //Local Essential Tree related functions
  void essential_tree_exchange(vector<real4> &LETParticles,  tree_structure &tree, tree_structure &remoteTree);
  void create_local_essential_tree(real4* bodies, real4* multipole, real4* nodeSizeInfo, real4* nodeCenterInfo,
                                    double4 boxCenter, double4 boxSize, float group_eps, int start, int end,
                                    vector<real4> &particles,    vector<real4> &multipoleData,
                                    vector<real4> &nodeSizeData, vector<real4> &nodeCenterData);
  void create_local_essential_tree_count(real4* bodies, real4* multipole, real4* nodeSizeInfo, real4* nodeCenterInfo,
                                    double4 boxCenter, double4 boxSize, float group_eps, int start, int end,
                                    int &particleCount, int &nodeCount);
  void create_local_essential_tree_fill(real4* bodies, real4* velocities, real4* multipole, real4* nodeSizeInfo, real4* nodeCenterInfo,
                                    double4 boxCenter, double4 boxSize,  float group_eps, int start, int end,
                                    int particleCount, int nodeCount, real4 *dataBuffer);

  real4* MP_exchange_bhlist(int ibox, int isource,
                                int bufferSize, real4 *letDataBuffer);


  void ICRecv(int procId, vector<real4> &bodyPositions, vector<real4> &bodyVelocities,  vector<int> &bodiesIDs);
  void ICSend(int destination, real4 *bodyPositions, real4 *bodyVelocities,  int *bodiesIDs, int size);

//   void MP_exchange_bhlist(int ibox, int isource,
//                                 vector<real4> &particles,         vector<real4> &multipoleData,
//                                 vector<real4> &nodeSizeData,      vector<real4> &nodeCenterData,
//                                 vector<real4> &recv_particles,    vector<real4> &recv_multipoleData,
//                                 vector<real4> &recv_nodeSizeData, vector<real4> &recv_nodeCenterData);

  void makeLET();

  void getAllLETBoxes(real4 bMin, real4 bMax);

  //End functions for parallel code
  //

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
/*
  
  void calcGravityOnParticles(real4 *bodyPositions, real4 *bodyVelocities, int *bodyIDs);


  void calcGravityOnParticles(real4 *bodyPositions, real4 *bodyVelocities, int *bodyIDs, int n)
  {
	  tree->setN(n);

	  tree->allocateParticleMemory();

	  memcpy(&tree->bodies_pos[0], bodyPositions,  sizeof(real4)*n);
	  memcpy(&tree->bodies_vel[0], bodyVelocities, sizeof(real4)*n);
	  memcpy(&tree->bodies_idsi[0], bodyIDs,  sizeof(int)*n);

	  memcpy(&tree->bodies_Ppos[0], bodyPositions,  sizeof(real4)*n);
	  memcpy(&tree->bodies_Pvel[0], bodyVelocities, sizeof(real4)*n);

  tree->localTree.bodies_pos.h2d();
  tree->localTree.bodies_vel.h2d();
  tree->localTree.bodies_Ppos.h2d();
  tree->localTree.bodies_Pvel.h2d();
  tree->localTree.bodies_ids.h2d();
  
          

  //Start construction of the tree
  tree->sort_bodies(tree->localTree);
  tree->build(tree->localTree);
  tree->allocateTreePropMemory(tree->localTree);
  tree->compute_properties(tree->localTree);
  tree->setActiveGrpsFunc(this->localTree);
  tree->approximate_gravity(this->localTree);
	  //sort
	  //build
	  //properties
	  //integrate, need i-particles!!, they should be sorted...??!!!!


  }*/

  //End library functions


  void setDataSetProperties(int NTotalT = -1, int NFirstT = -1, int NSecondT = -1, int NThirdT = -1)
  {
    NTotal   = NTotalT;
    NFirst   = NFirstT;
    NSecond  = NSecondT;
    NThird   = NThirdT;
  }

//   octree(const int device = 0, const float _theta = 0.75, const float eps = 0.05,
//          string snapF = "", int snapI = -1, float tempTimeStep = 1.0 / 16.0, int tempTend = 10000) {
  octree(const int device = 0, const float _theta = 0.75, const float eps = 0.05,
         string snapF = "", int snapI = -1,  float tempTimeStep = 1.0 / 16.0, int tempTend = 1000,
         float killDistanceT = -1, int maxDistT = -1, int snapAdd = 0) {

    devContext_flag = false;
    iter = 0;
    t_current = t_previous = 0;

    //First init mpi
    int argc = 0;
    char **argv = NULL;
    mpiInit(argc, argv, procId, nProcs);

//     devID = device;
    //devID = procId % nProcs;
//     devID = 1;
     devID = procId % 2;
    //devID = 1;
    cerr << "Preset device : "  << devID << "\t" << device << "\t" << nProcs <<endl;

    snapshotIter = snapI;
    snapshotFile = snapF;

    timeStep = tempTimeStep;
    tEnd     = tempTend;

    killDistance   = killDistanceT;
    removeDistance = maxDistT;
    snapShotAdd    = snapAdd;

    //TODO!
    inv_theta   = 1.0/_theta;
    eps2        = eps*eps;
    eta         = 0.02;
    theta       = _theta;

    nextSnapTime = 0;
    //Calc dt_limit
    const float dt_max = 1.0f / (1 << 4);

    dt_limit = int(-log(dt_max)/log(2.0f));

    store_energy_flag = true;
  }
  ~octree() {
  };

};


#endif // _OCTREE_H_
