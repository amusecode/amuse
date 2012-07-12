
#ifndef _BONSAI_LIB_
#define _BONSAI_LIB_

#define USE_CUDA

#ifdef USE_CUDA
  #include "my_cuda.h"
#else
  #include "my_ocl.h"
#endif


//extern "C" {
  
  my_dev::context &  sequoia_init(char** argv, int device = 0, const float _theta = 0.75, const float eps = 0.05);

my_dev::context &  sequoia_init(char** argv, 
				int device, 
				const float _theta, 
				const float eps, 
				const float rsearch2);

  int sequoia_cleanup();
  
  int sequoia_sortBodies(my_dev::dev_mem<float4>  &bodies_pos, my_dev::dev_mem<uint> &permutation, int n_bodies);
  int sequoia_reorderReal4(my_dev::dev_mem<real4>  &data, my_dev::dev_mem<uint> &permutation, int n_items);
  int sequoia_reorderInt1(my_dev::dev_mem<int>  &data, my_dev::dev_mem<uint> &permutation, int n_items);
  
  int sequoia_buildTreeStructure(my_dev::dev_mem<real4>  &bodies_pos, int n_bodies);
  int sequoia_computeTreeProperties(my_dev::dev_mem<real4>  &bodies_pos, int n_bodies);
  int sequoia_createGroups(my_dev::dev_mem<real4> &bodies_pos, int n_bodies);
  
  int sequoia_computeGravity(my_dev::dev_mem<real4> &j_bodies_pos, 
                          my_dev::dev_mem<real4> &i_bodies_pos, 
                          my_dev::dev_mem<real4> &i_bodies_acc, 
                          my_dev::dev_mem<real>  &i_bodies_ds2, 
                          my_dev::dev_mem<int>   &i_bodies_ngb, 
                          int                     n_i_bodies);
/*
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
                                     my_dev::dev_mem<int>   &i_bodies_ngb);      //OUT  J-ID of the nearest neighbour for I-particles
*/
// by M.I.
int sequoia_setParticlesAndGetGravity(my_dev::dev_mem<real4> &j_bodies_pos,         //Positions J-particles
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
				      my_dev::dev_mem<int>   &i_bodies_Nngb);      //OUT  the nmuber of the nearest neighbour for I-particles

  
//}

int sequoia_setParticlesAndGetGravity_firsthalf(my_dev::dev_mem<real4> &j_bodies_pos,
						 my_dev::dev_mem<int>   &j_bodies_ids,
						 int                    n_j_bodies,	
						 bool                   sortJBodies);

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
						int &out_n_levels);

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
								     int &out_n_nodes);

int sequoia_setParticlesAndGetGravity_lasthalf(my_dev::dev_mem<real4> &j_bodies_pos,
					       my_dev::dev_mem<real4> &i_bodies_pos,
					       my_dev::dev_mem<int>   &i_bodies_ids,
					       int                    n_i_bodies, 
					       bool                   sortIBodies,
					       my_dev::dev_mem<real4> &i_bodies_acc,
					       my_dev::dev_mem<real>  &i_bodies_ds2,
					       my_dev::dev_mem<int>   &i_bodies_ngb,
					       my_dev::dev_mem<int>   &i_bodies_Nngb);

#endif
