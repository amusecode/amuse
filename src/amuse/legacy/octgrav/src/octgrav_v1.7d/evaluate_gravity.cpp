#include "octgrav.h"

double get_time() {
  struct timeval Tvalue;
  struct timezone dummy;
  
  gettimeofday(&Tvalue,&dummy);
  return ((double) Tvalue.tv_sec +
	  1.e-6*((double) Tvalue.tv_usec));
}

void octgrav::evaluate_gravity(vector<float4> &bodies_pos_in,
			       vector<float4> &bodies_grav) {
  
  if (first_run) {
    initCUDA();
    first_run = false;
  }
  
  double t1   = get_time();
  fprintf(stderr, "Load data ... ");
  float4 com = load_data(bodies_pos_in);
  fprintf(stderr, "done in %lf sec \n", get_time() - t1);
  fprintf(stderr, "system_com= [%g %g %g %g] \n\n",
	  com.x, com.y, com.z, com.w);
  
  fprintf(stderr, "\n\n old_tree= %d \n rebuild_tree= %d\n\n", old_tree, rebuild_tree);

  rebuild_tree = true;
  if (old_tree == false || rebuild_tree == true) {
    t1   = get_time();
    fprintf(stderr, "Building gravity tree ... \n");
    tree.destroy_tree();
    tree.build_tree(bodies_pos, bodies_id);
    old_tree = false;
    rebuild_tree = false;
    fprintf(stderr, "done in %lf sec \n", get_time() - t1);
  }
  
  bool repeat = true;

  while(repeat) {
    
    t1 = get_time();
    fprintf(stderr, "Compute multipole moments ... ");
    rebuild_tree = tree.compute_multipole_moments(bodies_pos, 
						  old_tree,
						  rebuild_tree);
    fprintf(stderr, "\n\nrebuild_tree= %d\n\n", rebuild_tree);
    old_tree = true;
    fprintf(stderr, "done in %lf sec \n", get_time() - t1);
    if (rebuild_tree) {
      t1   = get_time();
      fprintf(stderr, "Building gravity tree ... \n");
      tree.destroy_tree();
      tree.build_tree(bodies_pos, bodies_id);
      old_tree = false;
      rebuild_tree = false;
      fprintf(stderr, "done in %lf sec \n", get_time() - t1);
    } else {
      repeat = false;
    }

  }    

  float4 root_com = tree.get_com();
  fprintf(stderr, "root_com= [%g %g %g %g] \n",
	  root_com.x, root_com.y, root_com.z, root_com.w);
  
  t1 = get_time();
  fprintf(stderr, "Generating oct_tree lists ... \n");
  child_struct child;
  children_list.resize(2);
  node_list.resize(8);
  cell_list.clear();
  leaf_list.clear();
  n_crit = NCRIT;
  tree.generate_node_list(child, 0,
			  children_list,
			  node_list);
  tree.generate_cell_list(n_crit, cell_list);
  fprintf(stderr, "done in %lf sec \n", get_time() - t1);
  
  t1 = get_time();
  fprintf(stderr, "Reorder bodies ... ");

  hst.bodies_pos  = (float4*)malloc(n_bodies * sizeof(float4));
  hst.bodies_grav = (float4*)malloc(n_bodies * sizeof(float4));

  reorder_bodies();
  fprintf(stderr, "done in %lf sec \n", get_time() - t1);
  
  fprintf(stderr, "  children_list.size= %d \n", (int)children_list.size());
  fprintf(stderr, "  node_list.size=     %d \n", (int)node_list.size());
  fprintf(stderr, "  leaf_list.size=     %d \n", (int)leaf_list.size());
  fprintf(stderr, "  cell_list.size=     %d \n", (int)cell_list.size());

  t1 = get_time();
  fprintf(stderr, "Allocating host memory ... ");
  allocate_host_memory();
  fprintf(stderr, "done in %lf sec \n", get_time() - t1);
  t1 = get_time();
  fprintf(stderr, "Allocating device memory ... ");
  allocate_cuda_memory();
  fprintf(stderr, "done in %lf sec \n", get_time() - t1);
  

  children_list[0] = children_list[children_list.size() - 2];
  children_list[1] = children_list[children_list.size() - 1];

  for (int i = 0; i < 8; i++)
    node_list[i] = node_list[node_list.size() - 8 + i];

  t1 = get_time();
  fprintf(stderr, "Preparing & sending data to the device ... ");
  prepare_data_for_device();
  copy_data_to_device();
  fprintf(stderr, "done in %lf sec \n", get_time() - t1);

  double dt_cuda = CUDA_evaluate_gravity();
  fprintf(stderr, "dt_cuda= %lf\n", dt_cuda);
  
  t1 = get_time();
  fprintf(stderr, "Copying data form the device ... ");
  copy_data_from_device();
  fprintf(stderr, "done in %lf sec \n", get_time() - t1);

  for (int i = 0; i < n_bodies; i++) {
    int j = bodies_id_orig[i];
//     fprintf(stderr, "i= %d, j= %d: %f %f \n", i, j, hst.bodies_grav[i].x, hst.bodies_grav[i].w);
    bodies_grav[j] = hst.bodies_grav[i];
//      bodies_grav[i] = hst.bodies_grav[i];
  }
  
  
  t1 = get_time();
  fprintf(stderr, "Freeing memory ... ");
  free_cuda_memory();
  free_host_memory();
  fprintf(stderr, "done in %lf sec \n", get_time() - t1);

}
