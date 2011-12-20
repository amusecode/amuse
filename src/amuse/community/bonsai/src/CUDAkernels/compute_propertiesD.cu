#include "support_kernels.cu"
#include <stdio.h>


//////////////////////////////
//////////////////////////////
//////////////////////////////

//Helper functions for leaf-nodes
__device__ void compute_monopole(double &mass, double &posx,
                                 double &posy, double &posz,
                                 float4 pos)
{
  mass += pos.w;
  posx += pos.w*pos.x;
  posy += pos.w*pos.y;
  posz += pos.w*pos.z;
}

__device__ void compute_quadropole(double &oct_q11, double &oct_q22, double &oct_q33,
                                   double &oct_q12, double &oct_q13, double &oct_q23,
                                   float4 pos)
{
  oct_q11 += pos.w * pos.x*pos.x;
  oct_q22 += pos.w * pos.y*pos.y;
  oct_q33 += pos.w * pos.z*pos.z;
  oct_q12 += pos.w * pos.x*pos.y;
  oct_q13 += pos.w * pos.y*pos.z;
  oct_q23 += pos.w * pos.z*pos.x;
}

__device__ void compute_bounds(float3 &r_min, float3 &r_max,
                               float4 pos)
{
  r_min.x = fminf(r_min.x, pos.x);
  r_min.y = fminf(r_min.y, pos.y);
  r_min.z = fminf(r_min.z, pos.z);

  r_max.x = fmaxf(r_max.x, pos.x);
  r_max.y = fmaxf(r_max.y, pos.y);
  r_max.z = fmaxf(r_max.z, pos.z);
}

//Non-leaf node helper functions
__device__ void compute_monopole_node(double &mass, double &posx,
                                 double &posy, double &posz,
                                 double4  pos)
{
  mass += pos.w;
  posx += pos.w*pos.x;
  posy += pos.w*pos.y;
  posz += pos.w*pos.z;
}


__device__ void compute_quadropole_node(double &oct_q11, double &oct_q22, double &oct_q33,
                                        double &oct_q12, double &oct_q13, double &oct_q23,
                                        double4 Q0, double4 Q1)
{
  oct_q11 += Q0.x;
  oct_q22 += Q0.y;
  oct_q33 += Q0.z;
  oct_q12 += Q1.x;
  oct_q13 += Q1.y;
  oct_q23 += Q1.z;
}

__device__ void compute_bounds_node(float3 &r_min, float3 &r_max,
                                    float4 node_min, float4 node_max)
{
  r_min.x = fminf(r_min.x, node_min.x);
  r_min.y = fminf(r_min.y, node_min.y);
  r_min.z = fminf(r_min.z, node_min.z);

  r_max.x = fmaxf(r_max.x, node_max.x);
  r_max.y = fmaxf(r_max.y, node_max.y);
  r_max.z = fmaxf(r_max.z, node_max.z);
}


extern "C" __global__ void compute_leaf(const int n_leafs,
                                            uint *leafsIdxs,
                                            uint2 *node_bodies,
                                            real4 *body_pos,
                                            double4 *multipole,
                                            real4 *nodeLowerBounds,
                                            real4 *nodeUpperBounds,
//                                             float3 *lowerBounds,
//                                             float3 *upperBounds,
                                            real4  *body_vel) {

  const uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  const uint tid = threadIdx.x;
  const uint id  = bid * blockDim.x + tid;


  volatile __shared__ float3 shmem[256];
  volatile float3 *sh_rmin = (float3*)&shmem [ 0];
  volatile float3 *sh_rmax = (float3*)&shmem[128];

  //Set the shared memory for these threads and exit the thread
  if (id >= n_leafs)
  {
    sh_rmin[tid].x = +1e10f; sh_rmin[tid].y = +1e10f; sh_rmin[tid].z = +1e10f;
    sh_rmax[tid].x = -1e10f; sh_rmax[tid].y = -1e10f; sh_rmax[tid].z = -1e10f;
    return;
  }


  //Since nodes are intermixes with non-leafs in the node_bodies array
  //we get a leaf-id from the leafsIdxs array
  int nodeID = leafsIdxs[id];

  const uint2 bij          =  node_bodies[nodeID];
  const uint firstChild    =  bij.x & ILEVELMASK;
  const uint lastChild     =  bij.y;  //TODO maybe have to increase it by 1

  //Variables holding properties and intermediate answers
  float4 p;

  double mass, posx, posy, posz;
  mass = posx = posy = posz = 0.0;

  double oct_q11, oct_q22, oct_q33;
  double oct_q12, oct_q13, oct_q23;

  oct_q11 = oct_q22 = oct_q33 = 0.0;
  oct_q12 = oct_q13 = oct_q23 = 0.0;
  float3 r_min, r_max;
  r_min = (float3){+1e10f, +1e10f, +1e10f};
  r_max = (float3){-1e10f, -1e10f, -1e10f};

  //Loop over the children=>particles=>bodys
  //unroll increases register usage #pragma unroll 16
  float maxEps = -100.0f;
  int count=0;
  for(int i=firstChild; i < lastChild; i++)
  {
    p      = body_pos[i];
    maxEps = fmaxf(body_vel[i].w, maxEps);      //Determine the max softening within this leaf
    count++;
    compute_monopole(mass, posx, posy, posz, p);
    compute_quadropole(oct_q11, oct_q22, oct_q33, oct_q12, oct_q13, oct_q23, p);
    compute_bounds(r_min, r_max, p);
  }

  double4 mon = {posx, posy, posz, mass};

  double im = 1.0/mon.w;
  if(mon.w == 0) im = 0;        //Allow tracer/massless particles
  mon.x *= im;
  mon.y *= im;
  mon.z *= im;

  double4 Q0, Q1;
  Q0 = (double4){oct_q11, oct_q22, oct_q33, maxEps}; //Store max softening
  Q1 = (double4){oct_q12, oct_q13, oct_q23, 0};

  //Store the leaf properties
  multipole[3*nodeID + 0] = mon;       //Monopole
  multipole[3*nodeID + 1] = Q0;        //Quadropole
  multipole[3*nodeID + 2] = Q1;        //Quadropole

  //Store the node boundaries
  nodeLowerBounds[nodeID] = (float4){r_min.x, r_min.y, r_min.z, 0.0f};
  nodeUpperBounds[nodeID] = (float4){r_max.x, r_max.y, r_max.z, 1.0f};  //4th parameter is set to 1 to indicate this is a leaf

  //Global domain boundaries using reduction
/*  sh_rmin[tid].x = r_min.x; sh_rmin[tid].y = r_min.y; sh_rmin[tid].z = r_min.z;
  sh_rmax[tid].x = r_max.x; sh_rmax[tid].y = r_max.y; sh_rmax[tid].z = r_max.z;
  __syncthreads();

  //TODO what happens here if some of the threads have been cancelled because of
  // id >= nleaf
  //Than we run the risk of comparing datas with garbade in shem

  //Reduction of the global boundaries of the system
  if(blockDim.x >= 128) if (tid < 64) {sh_MinMax(tid, tid + 64, &r_min, &r_max, sh_rmin, sh_rmax);} __syncthreads();
  if (tid < 32)
  {
    sh_MinMax(tid, tid + 32, &r_min, &r_max, sh_rmin,sh_rmax);
    sh_MinMax(tid, tid + 16, &r_min, &r_max, sh_rmin,sh_rmax);
    sh_MinMax(tid, tid +  8, &r_min, &r_max, sh_rmin,sh_rmax);
    sh_MinMax(tid, tid +  4, &r_min, &r_max, sh_rmin,sh_rmax);
    sh_MinMax(tid, tid +  2, &r_min, &r_max, sh_rmin,sh_rmax);
    sh_MinMax(tid, tid +  1, &r_min, &r_max, sh_rmin,sh_rmax);
  }
  __syncthreads();

  //Store the results
  if(tid == 0)
  {
    //Compiler doesnt allow: volatile float3 = float3
    lowerBounds[bid].x = sh_rmin[0].x; lowerBounds[bid].y = sh_rmin[0].y; lowerBounds[bid].z = sh_rmin[0].z;
    upperBounds[bid].x = sh_rmax[0].x; upperBounds[bid].y = sh_rmax[0].y; upperBounds[bid].z = sh_rmax[0].z;
  }
*/
  return;
}


//Function goes level by level (starting from deepest) and computes
//the properties of the non-leaf nodes
extern "C" __global__ void compute_non_leaf(const int curLevel,         //Level for which we calc
                                            uint  *leafsIdxs,           //Conversion of ids
                                            uint  *node_level_list,     //Contains the start nodes of each lvl
                                            uint  *n_children,          //Reference from node to first child and number of childs
                                            double4 *multipole,
                                            real4 *nodeLowerBounds,
                                            real4 *nodeUpperBounds){

  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;

  const int idx = bid * (blockDim.x * blockDim.y) + tid;

  const int endNode   = node_level_list[curLevel];
  const int startNode = node_level_list[curLevel-1];


  if(idx >= (endNode-startNode))     return;

  const int nodeID = leafsIdxs[idx + startNode];

  //Get the children info
  const uint firstChild = n_children[nodeID] & 0x0FFFFFFF;                  //TODO make this name/define?
  const uint nChildren  = ((n_children[nodeID]  & 0xF0000000) >> 28); //TODO make this name/define?

  //Variables
  double mass, posx, posy, posz;
  mass = posx = posy = posz = 0.0;

  double oct_q11, oct_q22, oct_q33;
  double oct_q12, oct_q13, oct_q23;

  oct_q11 = oct_q22 = oct_q33 = 0.0;
  oct_q12 = oct_q13 = oct_q23 = 0.0;

  float3 r_min, r_max;
  r_min = (float3){+1e10f, +1e10f, +1e10f};
  r_max = (float3){-1e10f, -1e10f, -1e10f};

  //Process the children (1 to 8)
  float maxEps = -100.0f;
  for(int i=firstChild; i < firstChild+nChildren; i++)
  {
    //Gogo process this data!
    double4 tmon = multipole[3*i + 0];

    maxEps = max(multipole[3*i + 1].w, maxEps);

    compute_monopole_node(mass, posx, posy, posz, tmon);
    compute_quadropole_node(oct_q11, oct_q22, oct_q33, oct_q12, oct_q13, oct_q23,
                            multipole[3*i + 1], multipole[3*i + 2]);
    compute_bounds_node(r_min, r_max, nodeLowerBounds[i], nodeUpperBounds[i]);
  }

  //Save the bounds
  nodeLowerBounds[nodeID] = (float4){r_min.x, r_min.y, r_min.z, 0.0f};
  nodeUpperBounds[nodeID] = (float4){r_max.x, r_max.y, r_max.z, 0.0f}; //4th is set to 0 to indicate a non-leaf

  //Regularize and store the results
  double4 mon = {posx, posy, posz, mass};
  double im = 1.0/mon.w;
  if(mon.w == 0) im = 0; //Allow tracer/massless particles

  mon.x *= im;
  mon.y *= im;
  mon.z *= im;

  double4 Q0, Q1;
  Q0 = (double4){oct_q11, oct_q22, oct_q33, maxEps}; //store max Eps
  Q1 = (double4){oct_q12, oct_q13, oct_q23, 0};

  multipole[3*nodeID + 0] = mon;        //Monopole
  multipole[3*nodeID + 1] = Q0;         //Quadropole1
  multipole[3*nodeID + 2] = Q1;         //Quadropole2

  return;
}
extern "C" __global__ void compute_scaling(const int node_count,
                                           real4 corner2, //TODO can be removed
                                           double4 *multipole,
                                           real4 *nodeLowerBounds,
                                           real4 *nodeUpperBounds,
                                           uint  *n_children,
//                                            uint4 *node_data,
                                           real4 *multipoleF,
                                           float theta,
                                           real4 *boxSizeInfo,
					   real4 *boxCenterInfo,
                                           uint2 *node_bodies){

  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;

  const int idx = bid * (blockDim.x * blockDim.y) + tid;

  if(idx >= node_count)     return;

  double4 monD, Q0, Q1;

  monD = multipole[3*idx + 0];        //Monopole
  Q0   = multipole[3*idx + 1];        //Quadropole1
  Q1   = multipole[3*idx + 2];        //Quadropole2

  //Scale the quadropole
  double im = 1.0 / monD.w;
  if(monD.w == 0) im = 0;               //Allow tracer/massless particles
  Q0.x = Q0.x*im - monD.x*monD.x;
  Q0.y = Q0.y*im - monD.y*monD.y;
  Q0.z = Q0.z*im - monD.z*monD.z;
  Q1.x = Q1.x*im - monD.x*monD.y;
  Q1.y = Q1.y*im - monD.y*monD.z;
  Q1.z = Q1.z*im - monD.x*monD.z;

  //Switch the y and z parameter
  double temp = Q1.y;
  Q1.y = Q1.z; Q1.z = temp;

  //Convert the doubles to floats
  float4 mon            = (float4){monD.x, monD.y, monD.z, monD.w};
  multipoleF[3*idx + 0] = mon;
  multipoleF[3*idx + 1] = (float4){Q0.x, Q0.y, Q0.z, Q0.w};        //Quadropole1
  multipoleF[3*idx + 2] = (float4){Q1.x, Q1.y, Q1.z, Q1.w};        //Quadropole2

  float4 r_min, r_max;
  r_min = nodeLowerBounds[idx];
  r_max = nodeUpperBounds[idx];


  //Determine the size of the node based on the center of mass and the bounds of the node
  #if 0
  float size =  powf(fmaxf(fabs(mon.x-r_min.x), fabs(mon.x-r_max.x)), 2) +
                powf(fmaxf(fabs(mon.y-r_min.y), fabs(mon.y-r_max.y)), 2) +
                powf(fmaxf(fabs(mon.z-r_min.z), fabs(mon.z-r_max.z)), 2);
  size = sqrtf(size)*0.5;
  #else
//  float3 size3  = (float3){fmaxf(fabs(mon.x-r_min.x), fabs(mon.x-r_max.x)),
//                          fmaxf(fabs(mon.y-r_min.y), fabs(mon.y-r_max.y)),
//                          fmaxf(fabs(mon.z-r_min.z), fabs(mon.z-r_max.z))};
//  float size    = fmaxf(size3.x, fmaxf(size3.y, size3.z));
//   size *= 0.5f;
  #endif

/*
  #ifdef IMPBH
  //Calculations for improved barnes hut opening
  float3 boxCenter;
  boxCenter.x = 0.5*(r_min.x + r_max.x);
  boxCenter.y = 0.5*(r_min.y + r_max.y);
  boxCenter.z = 0.5*(r_min.z + r_max.z);

  float3 boxSize = (float3){fmaxf(fabs(boxCenter.x-r_min.x), fabs(boxCenter.x-r_max.x)),
                          fmaxf(fabs(boxCenter.y-r_min.y), fabs(boxCenter.y-r_max.y)),
                          fmaxf(fabs(boxCenter.z-r_min.z), fabs(boxCenter.z-r_max.z))};

//   float3 size3  = (float3){fmaxf(fabs(mon.x-r_min.x), fabs(mon.x-r_max.x)),
//                           fmaxf(fabs(mon.y-r_min.y), fabs(mon.y-r_max.y)),
//                           fmaxf(fabs(mon.z-r_min.z), fabs(mon.z-r_max.z))};

  //Calculate distance between center of the box and the center of mass
  float3 s3     = (float3){(boxCenter.x - mon.x), (boxCenter.y - mon.y), (boxCenter.z - mon.z)};
  double s      = sqrt((s3.x*s3.x) + (s3.y*s3.y) + (s3.z*s3.z));

//   if(idx < -10)
//   {
//     printf("%d\tcom:\t (%f, %f, %f) \tbox:  (%f, %f, %f) \n", idx, mon.x, mon.y, mon.z, boxCenter.x, boxCenter.y, boxCenter.z);
//     printf("%d\tmin:\t (%f, %f, %f) \tmax:  (%f, %f, %f) \n", idx, r_min.x, r_min.y, r_min.z, r_max.x, r_max.y, r_max.z);
//     printf("%d\tsize:\t (%f, %f, %f) \tbox:  (%f, %f, %f) \n", idx, boxSize.x, boxSize.y, boxSize.z, size3.x, size3.y, size3.z);
//     printf("%d\ts:\t (%f, %f, %f) \tbox:  (%f, %f, %f) \n", idx, s3.x, s3.y, s3.z, s, size3.y, size3.z);
//   }

  //BH: l/theta + s < d
  float l = 2*fmaxf(boxSize.x, fmaxf(boxSize.y, boxSize.z));

  float cellOp = (l/theta) + s;
  cellOp = cellOp*cellOp;

  //Store the box size and opening criteria
  cellOpening[idx].x = boxSize.x;
  cellOpening[idx].y = boxSize.y;
  cellOpening[idx].z = boxSize.z;
  cellOpening[idx].w = cellOp;

  #endif
*/


  float3 boxCenter;
  boxCenter.x = 0.5*(r_min.x + r_max.x);
  boxCenter.y = 0.5*(r_min.y + r_max.y);
  boxCenter.z = 0.5*(r_min.z + r_max.z);

  float3 boxSize = (float3){fmaxf(fabs(boxCenter.x-r_min.x), fabs(boxCenter.x-r_max.x)),
                          fmaxf(fabs(boxCenter.y-r_min.y), fabs(boxCenter.y-r_max.y)),
                          fmaxf(fabs(boxCenter.z-r_min.z), fabs(boxCenter.z-r_max.z))};

  //Calculate distance between center of the box and the center of mass
  float3 s3     = (float3){(boxCenter.x - mon.x), (boxCenter.y - mon.y), (boxCenter.z - mon.z)};
  double s      = sqrt((s3.x*s3.x) + (s3.y*s3.y) + (s3.z*s3.z));

  //Length of the box, note times 2 since we only computed half the distance before
  float l = 2*fmaxf(boxSize.x, fmaxf(boxSize.y, boxSize.z));

  //Extra check, shouldnt be necessary
//  if(l < 0.000001)
 //   l = 0.000001;

  //Store the box size and opening criteria
  boxSizeInfo[idx].x = boxSize.x;
  boxSizeInfo[idx].y = boxSize.y;
  boxSizeInfo[idx].z = boxSize.z;
  boxSizeInfo[idx].w = __int_as_float(n_children[idx]);

  boxCenterInfo[idx].x = boxCenter.x;
  boxCenterInfo[idx].y = boxCenter.y;
  boxCenterInfo[idx].z = boxCenter.z;

  //Extra check, shouldnt be necessary, probably it is otherwise the test for leaf can fail
  //So it IS important Otherwise 0.0 < 0 can fail, now it will be: -1e-12 < 0 
  if(l < 0.000001)
    l = 0.000001;

  #ifdef IMPBH
    float cellOp = (l/theta) + s;
  #else
    //Minimum distance method
    float cellOp = (l/theta); 
  #endif
    
  cellOp = cellOp*cellOp;

  if(r_max.w > 0)
  {
    cellOp = -cellOp;       //This is a leaf node
  }
  

  boxCenterInfo[idx].w = cellOp;
 
  //Change the indirections of the leaf nodes so they point to
  //the particle data
  bool leaf = (r_max.w > 0);
  if(leaf)
  {
    uint2 bij     = node_bodies[idx];
    uint pfirst   = bij.x & ILEVELMASK;
    uint nchild   = bij.y - pfirst;

    pfirst = pfirst | ((nchild-1) << LEAFBIT);
    boxSizeInfo[idx].w = __int_as_float(pfirst);
  }


  //float size = s; //We store distance com-center in size


  //Calculate the key based on geometrical center
 /* int4 crd;

  float domain_fac = corner.w;

  float idomain_fac = 1.0f / domain_fac;
  crd.x = (int)((boxCenter.x - corner.x) * idomain_fac + 0.5);
  crd.y = (int)((boxCenter.y - corner.y) * idomain_fac + 0.5);
  crd.z = (int)((boxCenter.z - corner.z) * idomain_fac + 0.5);
  uint2 key = get_key(crd);

*/
/*
  //Calculate the key
  int4 crd;

  float domain_fac = corner.w;

  float idomain_fac = 1.0f / domain_fac;
  crd.x = (int)((mon.x - corner.x) * idomain_fac + 0.5);
  crd.y = (int)((mon.y - corner.y) * idomain_fac + 0.5);
  crd.z = (int)((mon.z - corner.z) * idomain_fac + 0.5);
  uint2 key = get_key(crd);

  //Use the key to calculate back the position
  float3 pos;
  pos.x = crd.x*domain_fac + corner.x;
  pos.y = crd.y*domain_fac + corner.y;
  pos.z = crd.z*domain_fac + corner.z;

  //Adjust size based on the  key-based position of the node
  float ds = fmax(fabs(pos.x - mon.x), max(fabs(pos.y - mon.y), fabs(pos.z - mon.z)));
  temp     = size;
  size    += ds;

  #ifdef IMPBH
  //Box size, max size for now
  size = l;
  if(l < 0.000001)
    size = 0.000001;

  #endif
*/

/*
  if(r_max.w > 0)
  {
    size = -size;       //This is a leaf node
  }

  //nchildren contains the node to node references
  //we also need to use node_bodies to get the
  //leaf-particle references
  node_data[idx] = (uint4){key.x, key.y,
                          float_as_int(size),
                          n_children[idx]};

*/
//   r_min.w = size;
//   nodeLowerBounds[idx] = r_min;


  return;
}


//Modify the references to the fist body and the number of bodys
//for the leafs
//Also copy the node_data to the group data
extern "C" __global__ void copyNodeDataToGroupData(const int n_groups,
                                                   const int n_nodes,
                                                   uint4 *node_data,
//                                                    uint4 *group_data,
                                                   uint2 *node_bodies,
                                                   int   *group_list,                                                
						   real4 *boxCenterInfo,
                                                   real4 *boxSizeInfo,
						   real4 *groupCenterInfo,
                                                   real4 *groupSizeInfo){
  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;

  const int idx = bid * (blockDim.x * blockDim.y) + tid;

  if(idx >= n_nodes)     return;

  //Copy the data and change the children data
  //Instead of pointing to child nodes we want it to point to
  //particles
//   uint4 nodeData = node_data[idx];
//   bool leaf =  __int_as_float(nodeData.z) <= 0;

  float temp = boxCenterInfo[idx].w;
  bool leaf = temp <= 0;

  //uint2 bij2     = node_bodies[idx];
  //uint pfirst2   = bij2.x & ILEVELMASK;
  //uint nchild2   = bij2.y - pfirst2;

  
  //Change the indirections of the leaf nodes so they point to
  //the particle data
  if(leaf)
  {
    uint2 bij     = node_bodies[idx];
    uint pfirst   = bij.x & ILEVELMASK;
    uint nchild   = bij.y - pfirst;

    pfirst = pfirst | ((nchild-1) << LEAFBIT);
    boxSizeInfo[idx].w = __int_as_float(pfirst);
  }

  //Now fill in the group data
  if(idx >= n_groups)     return;

  int nodeID         = group_list[idx];
  real4 nodeData     = boxSizeInfo[nodeID];

  uint2 bij     = node_bodies[nodeID];
  int pfirst    = bij.x & ILEVELMASK;
  int nchild    = bij.y - pfirst;

  pfirst = pfirst | (nchild-1) << CRITBIT;
  nodeData.w = __int_as_float(pfirst);

  groupSizeInfo[idx]   = nodeData;  
  groupCenterInfo[idx] = boxCenterInfo[nodeID];
}



//Compute the properties for the groups
extern "C" __global__ void setPHGroupData(const int n_groups,
                                          const int n_particles,   
                                          real4 *bodies_pos,
                                          int2  *group_list,                                                
                                          real4 *groupCenterInfo,
                                          real4 *groupSizeInfo){
  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;

  if(bid >= n_groups)     return;

  //Do a reduction on the particles assigned to this group

  volatile __shared__ float3 shmem[2*NCRIT];
  volatile float3 *sh_rmin = (float3*)&shmem [ 0];
  volatile float3 *sh_rmax = (float3*)&shmem[NCRIT];

  float3 r_min = (float3){+1e10f, +1e10f, +1e10f};
  float3 r_max = (float3){-1e10f, -1e10f, -1e10f};

  int start = group_list[bid].x;
  int end   = group_list[bid].y;
  
  int partIdx = start + threadIdx.x;

  //Set the shared memory with the data
  if (partIdx >= end)
  {
    sh_rmin[tid].x = r_min.x; sh_rmin[tid].y = r_min.y; sh_rmin[tid].z = r_min.z;
    sh_rmax[tid].x = r_max.x; sh_rmax[tid].y = r_max.y; sh_rmax[tid].z = r_max.z;
  }
  else
  {
    sh_rmin[tid].x = r_min.x = bodies_pos[partIdx].x; sh_rmin[tid].y = r_min.y = bodies_pos[partIdx].y; sh_rmin[tid].z = r_min.z = bodies_pos[partIdx].z;
    sh_rmax[tid].x = r_max.x = bodies_pos[partIdx].x; sh_rmax[tid].y = r_max.y = bodies_pos[partIdx].y; sh_rmax[tid].z = r_max.z = bodies_pos[partIdx].z;
  }


  __syncthreads();
  // do reduction in shared mem  
  if(blockDim.x >= 512) if (tid < 256) {sh_MinMax(tid, tid + 256, &r_min, &r_max, sh_rmin, sh_rmax);} __syncthreads();
  if(blockDim.x >= 256) if (tid < 128) {sh_MinMax(tid, tid + 128, &r_min, &r_max, sh_rmin, sh_rmax);} __syncthreads();
  if(blockDim.x >= 128) if (tid < 64)  {sh_MinMax(tid, tid + 64,  &r_min, &r_max, sh_rmin, sh_rmax);} __syncthreads();

  if(blockDim.x >= 64) if (tid < 32)  {sh_MinMax(tid, tid + 32, &r_min, &r_max, sh_rmin, sh_rmax); }
  if(blockDim.x >= 32) if (tid < 16) { sh_MinMax(tid, tid + 16, &r_min, &r_max, sh_rmin, sh_rmax); }

  if(tid < 8)
  {
    sh_MinMax(tid, tid +  8, &r_min, &r_max, sh_rmin, sh_rmax);
    sh_MinMax(tid, tid +  4, &r_min, &r_max, sh_rmin, sh_rmax);
    sh_MinMax(tid, tid +  2, &r_min, &r_max, sh_rmin, sh_rmax);
    sh_MinMax(tid, tid +  1, &r_min, &r_max, sh_rmin, sh_rmax);
  }
  
  // write result for this block to global mem
  if (tid == 0)
  {

    //Compute the group center and size
    float3 grpCenter;
    grpCenter.x = 0.5*(r_min.x + r_max.x);
    grpCenter.y = 0.5*(r_min.y + r_max.y);
    grpCenter.z = 0.5*(r_min.z + r_max.z);

    float3 grpSize = (float3){fmaxf(fabs(grpCenter.x-r_min.x), fabs(grpCenter.x-r_max.x)),
                              fmaxf(fabs(grpCenter.y-r_min.y), fabs(grpCenter.y-r_max.y)),
                              fmaxf(fabs(grpCenter.z-r_min.z), fabs(grpCenter.z-r_max.z))};

    //Store the box size and opening criteria
    groupSizeInfo[bid].x = grpSize.x;
    groupSizeInfo[bid].y = grpSize.y;
    groupSizeInfo[bid].z = grpSize.z;

    int nchild             = end-start;
    start                  = start | (nchild-1) << CRITBIT;
    groupSizeInfo[bid].w   = __int_as_float(start);  

    float l = max(grpSize.x, max(grpSize.y, grpSize.z));

    groupCenterInfo[bid].x = grpCenter.x;
    groupCenterInfo[bid].y = grpCenter.y;
    groupCenterInfo[bid].z = grpCenter.z;

    //Test stats for physical group size
    groupCenterInfo[bid].w = l;

    //groupCenterInfo[idx].w = free variable
  } //end tid == 0
}//end copyNode2grp

