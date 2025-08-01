#include "octree.h"



struct cmp_morton_key{
  bool operator () (const morton_struct &a, const morton_struct &b){
    return (cmp_uint2(a.key, b.key) < 1);
  }
};

void octree::sort_bodies(tree_structure &tree) {
  
  int n_b = tree.n;
  
  real4 r_min = {+1e10, +1e10, +1e10, +1e10}; 
  real4 r_max = {-1e10, -1e10, -1e10, -1e10}; 
  double  t0 = get_time();
  for (int i = 0; i < n_b; i++) {
    real4 pos = tree.bodies_pos[i];
    r_min.x = fmin(r_min.x, pos.x);
    r_min.y = fmin(r_min.y, pos.y);
    r_min.z = fmin(r_min.z, pos.z);
    
    r_max.x = fmax(r_max.x, pos.x);
    r_max.y = fmax(r_max.y, pos.y);
    r_max.z = fmax(r_max.z, pos.z);
  }
  printf("Determine min/max took: %lg \n", get_time()-t0);
  printf("Found boundarys: \n");
  printf("min: %f\t%f\t%f\tmax: %f\t%f\t%f \n", r_min.x,r_min.y,r_min.z,r_max.x,r_max.y,r_max.z);    
  
  real size = 1.001*fmax(r_max.z - r_min.z,
			 fmax(r_max.y - r_min.y, r_max.x - r_min.x));
  
  tree.corner = (real4){0.5*(r_min.x + r_max.x) - 0.5*size,
		   0.5*(r_min.y + r_max.y) - 0.5*size,
		   0.5*(r_min.z + r_max.z) - 0.5*size, size}; 
  
  tree.domain_fac   = size/(1 << MAXLEVELS);
  
  tree.corner.w = tree.domain_fac;
  
  
  float idomain_fac = 1.0/tree.domain_fac;
  float domain_fac = tree.domain_fac;
  
  printf("domain fac: %f idomain_fac: %f size: %f MAXLEVELS: %d \n", domain_fac, idomain_fac, size, MAXLEVELS);
  
  std::vector<morton_struct> keys(n_b);

  for (int i = 0; i < n_b; i++) {
    real4 pos = tree.bodies_pos[i];
    int3 crd;
    crd.x = (int)((pos.x - tree.corner.x) * idomain_fac + 0.5);
    crd.y = (int)((pos.y - tree.corner.y) * idomain_fac + 0.5);
    crd.z = (int)((pos.z - tree.corner.z) * idomain_fac + 0.5);
//     crd.x = (int)((pos.x - tree.corner.x) * domain_fac + 0.5);
//     crd.y = (int)((pos.y - tree.corner.y) * domain_fac + 0.5);
//     crd.z = (int)((pos.z - tree.corner.z) * domain_fac + 0.5);

    if(i < 10)
      printf("%d\t%d\t%d \n", crd.x, crd.y, crd.z);
//     
    keys[i].key   = get_key(crd);
    keys[i].value = i;    
    
    tree.bodies_key[i] = keys[i].key;
  }
  
  
  for(int i=0; i < 10; i++)
  {
    printf("%d\t%d\t%d\t%d\n", i,  0, tree.bodies_key[i].y, tree.bodies_key[i].x);
  }    
  for(int i=0; i < 10; i++)
  {
    printf("%d\t%f\t%f\t%f\n", i,  tree.bodies_pos[i].z, tree.bodies_pos[i].y, tree.bodies_pos[i].x);
  }  

  std::sort(keys.begin(), keys.end(), cmp_morton_key());

  std::vector<real4> bodies_pos_tmp(n_b);
  std::vector<uint2> bodies_key_tmp(n_b);
  for (int i = 0; i < n_b; i++) {
    bodies_pos_tmp[i] = tree.bodies_pos[i];
    bodies_key_tmp[i] = tree.bodies_key[i];
  }
  for (int i = 0; i < n_b; i++)  {
    int j = keys[i].value;
    tree.bodies_pos[i] = bodies_pos_tmp[j];
    tree.bodies_key[i] = bodies_key_tmp[j];
  }
  
//   tree.bodies_pos.d2h();
  for(int i=0; i < 10; i++)
  {
    printf("%d\t%d\t%d\t%d\n", i,  0, tree.bodies_key[i].y, tree.bodies_key[i].x);
  }

  for(int i=0; i < 10; i++)
  {
    printf("%d\t%f\t%f\t%f\n", i,  tree.bodies_pos[i].z, tree.bodies_pos[i].y, tree.bodies_pos[i].x);
  }
/*    
exit(0);*/
}
