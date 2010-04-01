#include <cutil.h>
#include <stdio.h>
#include "dev_evaluate_gravity.cu"

#include "n_per_cell.h"

#ifdef __DEVICE_EMULATION__
#endif

inline int n_norm(int n, int j) {
  n = ((n-1)/j) * j + j;
  if (n == 0) n = j;
  return n;
}

double get_time();

extern "C"
{
  int n_alloc;
  int    cuda_interaction_list_len;
  int    cuda_interaction_node_len;
  int    cuda_interaction_leaf_len;
  int    cuda_interaction_node_list;
  int    cuda_interaction_leaf_list;
  int    cuda_n_node;
  int    cuda_n_leaf;

  int3 *dev_interaction_list_len;
  int2 *dev_interaction_node_len;
  int2 *dev_interaction_leaf_len;
  int  *dev_interaction_node_list;
  int  *dev_interaction_leaf_list;
  int  *dev_n_node;
  int  *dev_n_leaf;


  void initCUDA() {   
    // CUT_DEVICE_INIT();
  }
  
  void allocateCUDAarray(void** pos, int n) {
    CUDA_SAFE_CALL(cudaMalloc((void**)pos, n));
    CUT_CHECK_ERROR("cudaMalloc failed!\n");
  }
  void deleteCUDAarray(void* pos) {
    CUDA_SAFE_CALL(cudaFree((void*)pos));
    CUT_CHECK_ERROR("cudaFree failed!\n");
  }
  void copyArrayToDevice(void* device, const void* host, int n) {
    CUDA_SAFE_CALL(cudaMemcpy(device, host, n, cudaMemcpyHostToDevice));
    CUT_CHECK_ERROR("cudaMemcpy (host->device) failed!\n");
  }
  void copyArrayFromDevice(void* host, const void* device, int n) {   
    CUDA_SAFE_CALL(cudaMemcpy(host, device, n, cudaMemcpyDeviceToHost));
    CUT_CHECK_ERROR("cudaMemcpy (device->host) failed!\n");
  }
  void threadSync() { cudaThreadSynchronize(); }

#define SAFE_ALLOC(what, oldv, newv) {\
    if ((oldv) < (newv)) { \
      if ((oldv) > 0) { \
	deleteCUDAarray((void*)(what)); \
      } \
      allocateCUDAarray((void**)&(what),  max(int(1.3*(oldv)), (int)(newv)) ); \
      n_alloc++; \
      oldv = max(int(1.3*(oldv)), (int)(newv)); \
    } }
 
  double host_evaluate_gravity(float  inv_opening_angle,
			       float  softening_squared,
			       
			       int    n_bodies,
			       float4 *bodies_pos,
			       float4 *bodies_grav,

			       int    n_children,
			       int4   *children,

			       int    n_nodes,
			       float4 root_pos,
			       float4 root_com,
			       float4 *node_pos,
			       float4 *node_com,
			       float4 *node_Qu,
			       float4 *node_Qd,
			       float4 *Oct1,
			       float4 *Oct2,
			       float2 *Oct3,
			       int    *n_in_node,
			       int    *node_bodies_offset,
			       
			       int    n_cells,
			       float4 *cell_pos,
			       float4 *cell_com,
			       int    *n_in_cell,
			       int    *cell_bodies_offset) {

    double t_begin = get_time();
    
    bodies_pos_tex.addressMode[0] = cudaAddressModeWrap;
    bodies_pos_tex.addressMode[1] = cudaAddressModeWrap;
    bodies_pos_tex.filterMode     = cudaFilterModePoint;
    bodies_pos_tex.normalized     = false;
    CUDA_SAFE_CALL(cudaBindTexture(0, bodies_pos_tex, bodies_pos, n_bodies * sizeof(float4)));

    /***************************************************/

    children_tex.addressMode[0] = cudaAddressModeWrap;
    children_tex.addressMode[1] = cudaAddressModeWrap;
    children_tex.filterMode     = cudaFilterModePoint;
    children_tex.normalized     = false;
    CUDA_SAFE_CALL(cudaBindTexture(0, children_tex, children, n_children * sizeof(int4)));
    
    /***************************************************/

    node_bodies_offset_tex.addressMode[0] = cudaAddressModeWrap;
    node_bodies_offset_tex.addressMode[1] = cudaAddressModeWrap;
    node_bodies_offset_tex.filterMode     = cudaFilterModePoint;
    node_bodies_offset_tex.normalized     = false;
    CUDA_SAFE_CALL(cudaBindTexture(0, node_bodies_offset_tex, node_bodies_offset, n_nodes * sizeof(int)));

    cell_bodies_offset_tex.addressMode[0] = cudaAddressModeWrap;
    cell_bodies_offset_tex.addressMode[1] = cudaAddressModeWrap;
    cell_bodies_offset_tex.filterMode     = cudaFilterModePoint;
    cell_bodies_offset_tex.normalized     = false;
    CUDA_SAFE_CALL(cudaBindTexture(0, cell_bodies_offset_tex, cell_bodies_offset, n_cells * sizeof(int)));
    
    /***************************************************/    

    node_pos_tex.addressMode[0] = cudaAddressModeWrap;
    node_pos_tex.addressMode[1] = cudaAddressModeWrap;
    node_pos_tex.filterMode     = cudaFilterModePoint;
    node_pos_tex.normalized     = false;
    CUDA_SAFE_CALL(cudaBindTexture(0, node_pos_tex, node_pos, n_nodes * sizeof(float4)));
    
    node_com_tex.addressMode[0] = cudaAddressModeWrap;
    node_com_tex.addressMode[1] = cudaAddressModeWrap;
    node_com_tex.filterMode     = cudaFilterModePoint;
    node_com_tex.normalized     = false;
    CUDA_SAFE_CALL(cudaBindTexture(0, node_com_tex, node_com, n_nodes * sizeof(float4)));

    node_Qu_tex.addressMode[0]  = cudaAddressModeWrap;
    node_Qu_tex.addressMode[1]  = cudaAddressModeWrap;
    node_Qu_tex.filterMode      = cudaFilterModePoint;
    node_Qu_tex.normalized      = false;
    CUDA_SAFE_CALL(cudaBindTexture(0, node_Qu_tex,  node_Qu,  n_nodes * sizeof(float4)));

    node_Qd_tex.addressMode[0]  = cudaAddressModeWrap;
    node_Qd_tex.addressMode[1]  = cudaAddressModeWrap;
    node_Qd_tex.filterMode      = cudaFilterModePoint;
    node_Qd_tex.normalized      = false;
    CUDA_SAFE_CALL(cudaBindTexture(0, node_Qd_tex,  node_Qd,  n_nodes * sizeof(float4)));

    Oct1_tex.addressMode[0]  = cudaAddressModeWrap;
    Oct1_tex.addressMode[1]  = cudaAddressModeWrap;
    Oct1_tex.filterMode      = cudaFilterModePoint;
    Oct1_tex.normalized      = false;
    CUDA_SAFE_CALL(cudaBindTexture(0, Oct1_tex,  Oct1,  n_nodes * sizeof(float4)));

    Oct2_tex.addressMode[0]  = cudaAddressModeWrap;
    Oct2_tex.addressMode[1]  = cudaAddressModeWrap;
    Oct2_tex.filterMode      = cudaFilterModePoint;
    Oct2_tex.normalized      = false;
    CUDA_SAFE_CALL(cudaBindTexture(0, Oct2_tex,  Oct2,  n_nodes * sizeof(float4)));

    Oct3_tex.addressMode[0]  = cudaAddressModeWrap;
    Oct3_tex.addressMode[1]  = cudaAddressModeWrap;
    Oct3_tex.filterMode      = cudaFilterModePoint;
    Oct3_tex.normalized      = false;
    CUDA_SAFE_CALL(cudaBindTexture(0, Oct3_tex,  Oct3,  n_nodes * sizeof(float2)));

    n_in_node_tex.addressMode[0] = cudaAddressModeWrap;
    n_in_node_tex.addressMode[1] = cudaAddressModeWrap;
    n_in_node_tex.filterMode     = cudaFilterModePoint;
    n_in_node_tex.normalized     = false;
    CUDA_SAFE_CALL(cudaBindTexture(0, n_in_node_tex, n_in_node, n_nodes * sizeof(int)));

    /***************************************************/

    cell_pos_tex.addressMode[0] = cudaAddressModeWrap;
    cell_pos_tex.addressMode[1] = cudaAddressModeWrap;
    cell_pos_tex.filterMode     = cudaFilterModePoint;
    cell_pos_tex.normalized     = false;
    CUDA_SAFE_CALL(cudaBindTexture(0, cell_pos_tex, cell_pos, n_cells * sizeof(float4)));

    cell_com_tex.addressMode[0] = cudaAddressModeWrap;
    cell_com_tex.addressMode[1] = cudaAddressModeWrap;
    cell_com_tex.filterMode     = cudaFilterModePoint;
    cell_com_tex.normalized     = false;
    CUDA_SAFE_CALL(cudaBindTexture(0, cell_com_tex, cell_com, n_cells * sizeof(float4)));

    n_in_cell_tex.addressMode[0] = cudaAddressModeWrap;
    n_in_cell_tex.addressMode[1] = cudaAddressModeWrap;
    n_in_cell_tex.filterMode     = cudaFilterModePoint;
    n_in_cell_tex.normalized     = false;
    CUDA_SAFE_CALL(cudaBindTexture(0, n_in_cell_tex, n_in_cell, n_cells * sizeof(int)));

    /***************************************************/
    
    CUDA_SAFE_CALL(cudaMemcpyToSymbol("inv_opening_angle", &inv_opening_angle, 
				      sizeof(float), 0, 
				      cudaMemcpyHostToDevice));

    CUDA_SAFE_CALL(cudaMemcpyToSymbol("softening_squared", &softening_squared, 
				      sizeof(float), 0, 
				      cudaMemcpyHostToDevice));

    CUDA_SAFE_CALL(cudaMemcpyToSymbol("root_pos", &root_pos, 
				      sizeof(float4), 0, 
				      cudaMemcpyHostToDevice));

    CUDA_SAFE_CALL(cudaMemcpyToSymbol("root_com", &root_com, 
				      sizeof(float4), 0, 
				      cudaMemcpyHostToDevice));

    CUDA_SAFE_CALL(cudaMemcpyToSymbol("n_nodes", &n_nodes, 
				      sizeof(int), 0, 
				      cudaMemcpyHostToDevice));

    CUDA_SAFE_CALL(cudaMemcpyToSymbol("n_cells", &n_cells, 
				      sizeof(int), 0, 
				      cudaMemcpyHostToDevice));

    CUDA_SAFE_CALL(cudaMemcpyToSymbol("n_bodies", &n_bodies, 
				      sizeof(int), 0, 
				      cudaMemcpyHostToDevice));
    
    /***************************************************/
  
    /*****   building interaction list   *****/
    
    int p = 128;
    int n_cells_dev = n_cells;
    if (n_cells_dev < p) 
      p = n_cells_dev;
    else
      n_cells_dev = n_norm(n_cells_dev, p);
    
    dim3 threads(p,1,1);
    dim3 grid(n_cells_dev/p, 1, 1);

    /*** computing interaction list length ***/
    
    int3 *hst_interaction_list_len;

//     int3 *dev_interaction_list_len = NULL
    SAFE_ALLOC(dev_interaction_list_len, cuda_interaction_list_len, n_cells * sizeof(int3));
//     allocateCUDAarray((void**)&dev_interaction_list_len,  n_cells * sizeof(int3));

    hst_interaction_list_len = (int3*)malloc(n_cells * sizeof(int3));
    
    double t1 = get_time();
    double dt_ilen = 0;
    fprintf(stderr, "   compute_interaction_list_len ... ");
    dev_compute_interaction_list_len<<<grid, threads>>>(dev_interaction_list_len);
    threadSync();
    fprintf(stderr, "   done in %lf seconds \n", (dt_ilen = get_time() - t1));

    /************************/
    /*** computing offset ***/
    /************************/
    
    copyArrayFromDevice(hst_interaction_list_len,
			dev_interaction_list_len, n_cells * sizeof(int3));

    int* hst_n_in_cell = (int*)malloc(n_cells * sizeof(int));
    copyArrayFromDevice(hst_n_in_cell, n_in_cell, n_cells * sizeof(int));


    int2 *hst_interaction_node_len; //, *dev_interaction_node_len;
    int2 *hst_interaction_leaf_len; //, *dev_interaction_leaf_len;
    hst_interaction_node_len = (int2*)malloc((n_cells + 1) * sizeof(int2));
    hst_interaction_leaf_len = (int2*)malloc((n_cells + 1) * sizeof(int2));

    SAFE_ALLOC(dev_interaction_node_len, cuda_interaction_node_len, n_cells * sizeof(int2));
    SAFE_ALLOC(dev_interaction_leaf_len, cuda_interaction_leaf_len, n_cells * sizeof(int2));
//     allocateCUDAarray((void**)&dev_interaction_node_len,  n_cells * sizeof(int2));
//     allocateCUDAarray((void**)&dev_interaction_leaf_len,  n_cells * sizeof(int2));

    long long int n_io = 0;
    int n_interacting_nodes_total = 0, n_interacting_leaves_total = 0;

    int n_interacting_nodes_max = 0;
    int n_interacting_leaves_max = 0;

    int n_blocks = 0;
    for (int i = 0; i < n_cells; i += p*NBLOCKS) {
      n_blocks += 1;
      int n_in_block = min(n_cells - i, p*NBLOCKS);
//       fprintf(stderr, "     block %d  n_in_block= %d cell_offset= %d\n",
// 	      n_blocks, n_in_block, i);
      
      hst_interaction_node_len[i].x = 0;
      hst_interaction_leaf_len[i].x = 0;
      int n_interacting_nodes = 0, n_interacting_leaves = 0;
      for (int j = 0; j < n_in_block; j++) {
	int3 val = hst_interaction_list_len[i+j];
	n_interacting_nodes  += val.x;
	n_interacting_leaves += val.y;
	n_io                 += val.z;
	
	hst_interaction_node_len[i+j+1].x = hst_interaction_node_len[i+j].x + val.x;
	hst_interaction_leaf_len[i+j+1].x = hst_interaction_leaf_len[i+j].x + val.y;
	
	hst_interaction_node_len[i+j].y = val.x;
	hst_interaction_leaf_len[i+j].y = val.y;
      }
      n_interacting_nodes_max     = max(n_interacting_nodes_max,  n_interacting_nodes);
      n_interacting_leaves_max    = max(n_interacting_leaves_max, n_interacting_leaves);
      n_interacting_nodes_total  += n_interacting_nodes;
      n_interacting_leaves_total += n_interacting_leaves;
    }
    
    copyArrayToDevice(dev_interaction_node_len,
		      hst_interaction_node_len, n_cells * sizeof(int2));
    copyArrayToDevice(dev_interaction_leaf_len,
		      hst_interaction_leaf_len, n_cells * sizeof(int2));
    free(hst_n_in_cell);

    fprintf(stderr, " *****************************************************\n");
    fprintf(stderr, "       n_blocks= %d  n_in_block= %d\n", n_blocks, p*NBLOCKS);
    fprintf(stderr, "   #interacting nodes=  %d [max in block= %d]\n", n_interacting_nodes_total, n_interacting_nodes_max);
    fprintf(stderr, "   #interacting leaves= %d [max in block= %d]\n", n_interacting_leaves_total, n_interacting_leaves_max);
    fprintf(stderr, "   read+write_len= %lg GB  (%lg GB/s)\n", 
	    n_io*4.0/pow(1024.0,3.0),
	    n_io*4.0/pow(1024.0,3.0)/dt_ilen);
    fprintf(stderr, " *****************************************************\n");
    
    
    /*********************************/
    /*** building interaction list ***/
    /*********************************/

//     int *dev_interaction_node_list, *dev_interaction_leaf_list;
//     allocateCUDAarray((void**)&dev_interaction_node_list,  n_interacting_nodes_max  * sizeof(int));
//     allocateCUDAarray((void**)&dev_interaction_leaf_list,  n_interacting_leaves_max * sizeof(int));
    SAFE_ALLOC(dev_interaction_node_list, cuda_interaction_node_list, n_interacting_nodes_max  * sizeof(int));
    SAFE_ALLOC(dev_interaction_leaf_list, cuda_interaction_leaf_list, n_interacting_leaves_max * sizeof(int));
    
    int *hst_n_node;
    int *hst_n_leaf;
    SAFE_ALLOC(dev_n_node, cuda_n_node, n_bodies * sizeof(int));
    SAFE_ALLOC(dev_n_leaf, cuda_n_leaf, n_bodies * sizeof(int));

//     allocateCUDAarray((void**)&dev_n_node,  n_bodies * sizeof(int));
//     allocateCUDAarray((void**)&dev_n_leaf,  n_bodies * sizeof(int));
    hst_n_node = (int*)malloc(n_bodies * sizeof(int));
    hst_n_leaf = (int*)malloc(n_bodies * sizeof(int));

    double dt_ibuild = 0, dt_node = 0, dt_leaf = 0;
    int cur_block = 0;
    for (int i = 0; i < n_cells; i += p*NBLOCKS) {
      cur_block++;
      int n_in_block = min(n_cells - i, p*NBLOCKS);
      dim3 threads(p, 1, 1);
      dim3 grid(n_norm(n_in_block, p)/p, 1, 1);
      if (n_in_block < p) {
	threads.x = n_in_block;
	grid.x    = 1;
      }
      fprintf(stderr, "   block %d out of %d\n", cur_block, n_blocks);
//       fprintf(stderr, "     n_in_block= %d cell_offset= %d\n",
// 	      n_in_block, i);
//        fprintf(stderr, "threads= [%d %d %d]\n", threads.x, threads.y, threads.z);
//        fprintf(stderr, "   grid= [%d %d %d]\n", grid.x, grid.y, grid.z);
      
      double dt_ibuild0 = 0;
      t1 = get_time();
      fprintf(stderr, "   build_interaction_list ... ");
      dev_build_interaction_list<<<grid, threads>>>(i, 
						    dev_interaction_node_list,
						    dev_interaction_node_len,
						    dev_interaction_leaf_list,
						    dev_interaction_leaf_len);
      threadSync();
      dt_ibuild0 = get_time() - t1;
      dt_ibuild += dt_ibuild0;
      fprintf(stderr, "   done in %lf seconds [%lf] \n", dt_ibuild0, dt_ibuild);
      
    
      /***************************************************/
      
      interaction_node_tex.addressMode[0] = cudaAddressModeWrap;
      interaction_node_tex.addressMode[1] = cudaAddressModeWrap;
      interaction_node_tex.filterMode     = cudaFilterModePoint;
      interaction_node_tex.normalized     = false;
      CUDA_SAFE_CALL(cudaBindTexture(0, interaction_node_tex, dev_interaction_node_list, n_interacting_nodes_max * sizeof(int)));
      
      interaction_leaf_tex.addressMode[0] = cudaAddressModeWrap;
      interaction_leaf_tex.addressMode[1] = cudaAddressModeWrap;
      interaction_leaf_tex.filterMode     = cudaFilterModePoint;
      interaction_leaf_tex.normalized     = false;
      CUDA_SAFE_CALL(cudaBindTexture(0, interaction_leaf_tex, dev_interaction_leaf_list, n_interacting_leaves_max * sizeof(int)));
      
      /***************************************************/

      int p_node = NTHREADS;
      int p_leaf = NCRIT;
      
      dim3 thread_node(p_node, 1, 1);
      dim3 thread_leaf(p_leaf, 1, 1);
      
      dim3 grid_node(n_in_block, 1, 1);
      dim3 grid_leaf(n_in_block, 1, 1);
      
      int shared_mem_size_leaf = 3 * p_leaf * sizeof(float4);
      int shared_mem_size_node = p_node * (3*sizeof(float4) + sizeof(float2));
      
      double dt_node0 = 0;
      fprintf(stderr, "   evaluate_gravity_node  ... "); t1 = get_time();
      dev_evaluate_gravity_node<<<grid_node, thread_node, shared_mem_size_node>>>(i,
										  bodies_grav,
										  dev_n_node,
										  dev_interaction_node_len);
      threadSync();
      dt_node0 = get_time() - t1;
      dt_node += dt_node0;
      fprintf(stderr, "   done in %lf seconds [%lf] \n", dt_node0, dt_node);
      
      double dt_leaf0 = 0;
      fprintf(stderr, "   evaluate_gravity_leaf  ... "); t1 = get_time();
      dev_evaluate_gravity_leaf<<<grid_leaf, thread_leaf, shared_mem_size_leaf>>>(i,
										  bodies_grav,
										  dev_n_leaf,
										  dev_interaction_leaf_len);
      threadSync();
      dt_leaf0 = get_time() - t1;
      dt_leaf += dt_leaf0;
      fprintf(stderr, "   done in %lf seconds [%lf] \n", dt_leaf0, dt_leaf);
      

      CUDA_SAFE_CALL(cudaUnbindTexture(interaction_leaf_tex));
      CUDA_SAFE_CALL(cudaUnbindTexture(interaction_node_tex));
      
    }      
    
    copyArrayFromDevice(hst_n_node, dev_n_node, n_bodies * sizeof(int));
    copyArrayFromDevice(hst_n_leaf, dev_n_leaf, n_bodies * sizeof(int));
    long long n_leaf = 0, n_node = 0;
    for (int i = 0; i < n_bodies; i++) {
      n_node += hst_n_node[i];
      n_leaf += hst_n_leaf[i];
    }

    double flops_per_node = 57;
    double flops_per_leaf = 21;
    

    double flops_per_node1 = 57 + 16;
    double flops_per_leaf1 = 21 + 16;
#ifdef OCTUPOLE
    flops_per_node  += 108;
    flops_per_node  += 108;
    flops_per_node1 += 108;
    flops_per_node1 += 108;
#endif
 
    fprintf(stderr, " *****************************************************\n");
    fprintf(stderr, "   #interacting nodes=  %d [max in block= %d]\n", n_interacting_nodes_total, n_interacting_nodes_max);
    fprintf(stderr, "   #interacting leaves= %d [max in block= %d]\n", n_interacting_leaves_total, n_interacting_leaves_max);
    fprintf(stderr, "   read+write_len= %lg GB  (%lg GB/s)\n", 
	    n_io*4.0/pow(1024.0,3.0),
	    n_io*4.0/pow(1024.0,3.0)/dt_ilen);
    n_io += n_interacting_nodes_total + n_interacting_leaves_total;
    fprintf(stderr, "   read+write_bld= %lg GB  (%lg GB/s)\n", 
	    n_io*4.0/pow(1024.0,3.0),
	    n_io*4.0/pow(1024.0,3.0)/dt_ibuild);
   
    fprintf(stderr, "    interaction statistics: \n");
    fprintf(stderr, "      n_nodes=  %ld (%lg [%lg] GLFOPS)\n", n_node, n_node*flops_per_node/dt_node/1e9, n_node*flops_per_node1/dt_node/1e9);
    fprintf(stderr, "      n_leaves= %ld (%lg [%lg] GFLOPS)\n", n_leaf, n_leaf*flops_per_leaf/dt_leaf/1e9,  n_leaf*flops_per_leaf1/dt_leaf/1e9);
    fprintf(stderr, " *****************************************************\n");


    /***************************************************/

    CUDA_SAFE_CALL(cudaUnbindTexture(n_in_cell_tex));
    CUDA_SAFE_CALL(cudaUnbindTexture(cell_pos_tex));
    CUDA_SAFE_CALL(cudaUnbindTexture(cell_com_tex));

    CUDA_SAFE_CALL(cudaUnbindTexture(n_in_node_tex));
    CUDA_SAFE_CALL(cudaUnbindTexture(node_Qu_tex));
    CUDA_SAFE_CALL(cudaUnbindTexture(node_Qd_tex));

    CUDA_SAFE_CALL(cudaUnbindTexture(Oct1_tex));
    CUDA_SAFE_CALL(cudaUnbindTexture(Oct2_tex));
    CUDA_SAFE_CALL(cudaUnbindTexture(Oct3_tex));

    CUDA_SAFE_CALL(cudaUnbindTexture(node_com_tex));
    CUDA_SAFE_CALL(cudaUnbindTexture(node_pos_tex));

    CUDA_SAFE_CALL(cudaUnbindTexture(cell_bodies_offset_tex));
    CUDA_SAFE_CALL(cudaUnbindTexture(node_bodies_offset_tex));

    CUDA_SAFE_CALL(cudaUnbindTexture(children_tex));
    
    CUDA_SAFE_CALL(cudaUnbindTexture(bodies_pos_tex));

    /***************************************************/


    free(hst_n_node);
    free(hst_n_leaf);
//     deleteCUDAarray((void*)dev_n_node);
//     deleteCUDAarray((void*)dev_n_leaf);

    free(hst_interaction_list_len);
//     deleteCUDAarray((void*)dev_interaction_list_len);

    free(hst_interaction_node_len);
    free(hst_interaction_leaf_len);
//     deleteCUDAarray((void*)dev_interaction_node_len);
//     deleteCUDAarray((void*)dev_interaction_leaf_len);

//     deleteCUDAarray((void*)dev_interaction_node_list);
//     deleteCUDAarray((void*)dev_interaction_leaf_list);
    
    return get_time() - t_begin;
  }

};
