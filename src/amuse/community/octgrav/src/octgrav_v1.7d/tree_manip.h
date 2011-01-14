#ifndef _TREE_MANIP_
#define _TREE_MANIP_

#include "oct_tree_base.h"

/* --- get list of particles in a block --- */
// template<int N_BODIES>
// oct_tree<N_BODIES>* oct_tree<N_BODIES>::get_node(vec& body_pos) {
//   oct_tree<N_BODIES>* node;
//   if (n_bodies < 0) {
//     node = children[octant(body_pos)].get_node(body_pos);
//   } else if (n_bodies > 0) {
//     return this;
//   } else {
//     return NULL;
//   }
//   return node;
// }

// template<int N_BODIES>
// void oct_tree<N_BODIES>::get_pblocks(int node_id, real node_range,
// 				     vec& body_pos, 
// 				     vector<vec>& bodies_pos,
// 				     int pblocks[], int n[]) {
//   if (n_bodies < 0) {
    
//     real s = node_range + 0.5 * size;
//     for (int i = 0; i < N_CHILDREN; i++) {
//       if (children[i].n_bodies != 0) {
// 	vec dr = body_pos - children[i].posc;
// 	if (abs(dr.x)<=s && abs(dr.y)<=s && abs(dr.z)<=s)
// 	  children[i].get_pblocks(node_id, node_range, body_pos, 
// 				  bodies_pos, pblocks, n);
//       }
//     }
    
//   } else {
    
//     if (id == node_id) {
//       n[1] = n[0];
//       for (int i = 0; i < n_bodies; i++) 
// 	pblocks[n[0]++] = bodies_id[i];
//       n[2] = n[0];
//     } else {
//       real s = node_range;
//       for (int i = 0; i < n_bodies; i++) {
//   	vec dr = body_pos - bodies_pos[bodies_id[i]];
//        	if (abs(dr.x)<=s && abs(dr.y)<=s && abs(dr.z)<=s)
// 	  pblocks[n[0]++] = bodies_id[i];
//       }
//     }
    
//   }
// }


// /* --- remove a body from the tree --- */
// template<int N_BODIES>
// void oct_tree<N_BODIES>::remove_body(vec& body_pos, int body_id) {
  
//   if (n_bodies < 0) {
//     children[octant(body_pos)].remove_body(body_pos, body_id);
//     --total_bodies;
    
//     /* 
//        if total_bodies < N_BODIES, 
//        copy bodies from children to parent
//     */
    
//     if (total_bodies <= N_BODIES) {
//       for (int i = 0; i < N_CHILDREN; i++) 
// 	for (int j = 0; j < children[i].n_bodies; j++) {
// 	  bodies_id.push_back(children[i].bodies_id[j]);
// 	}
//       delete[] children;
//       children = NULL;
//       n_bodies = total_bodies;
//     }
    
//   } else if (n_bodies <= N_BODIES) {
    
//     for (int i = 0; i < n_bodies; i++)
//       if (bodies_id[i] == body_id) {
// 	bodies_id.erase(bodies_id.begin() + i);
// 	n_bodies = bodies_id.size();
// 	total_bodies = n_bodies;
// 	return;
//       }
//     cerr << "oct_tree::remove(body) -- body not found " << endl;
//     sqrt(-1.0);
//     exit(-1);
//   } else {
//     fprintf(stderr, "n_bodies= %d  N_BODIES= %d \n", n_bodies, N_BODIES);
//     cerr << "oct_tree::remove(body) --- strange --- " << endl;
//     sqrt(-1.0);
//     exit(-1);
//   }
// }


/* ---  insert a body into the tree --- */

template<int N_BODIES>
void oct_tree<N_BODIES>::insert_body(register float4 body_pos,
				     register int body_id,
				     vector<float4> &bodies_pos) {
  
  /* --- n_bodies < 0, travel down the tree --- */
  
  if (n_bodies < 0) {
    children[octant(body_pos)].insert_body(body_pos, body_id, bodies_pos);
    total_bodies++;
    
    /* --- n_bodies < N_BODIES, insert body here --- */
    
  } else if (n_bodies < N_BODIES)  {

    bodies_id[n_bodies] = body_id; //.push_back(body_id);
    n_bodies++;
    total_bodies++;
    
    /* --- otherwise, leaf should be splitted further --- */
    
  } else {
    
    /* --- create new children --- */
    children = new oct_tree<N_BODIES>[N_CHILDREN];
    
    for (int i = 0; i < N_CHILDREN; i++) {
      children[i].pos = child_pos(i);
      children[i].level = level + 1;
      level_counter() = max(level_counter(), (long int)level + 1);
    }
    
    /* --- copy bodies to children --- */
    
    n_bodies = -1;
    
    for (int i = 0; i < N_BODIES; i++) 
      insert_body(bodies_pos[bodies_id[i]], bodies_id[i], bodies_pos);
    bodies_id.clear();
    total_bodies -= N_BODIES;
    
    insert_body(body_pos, body_id, bodies_pos);
  }
}

/* --- this routine builds the tree  --- */

template<int N_BODIES>
void oct_tree<N_BODIES>::build_tree(vector<float4>& bodies_pos_in,
				    vector<int>&    bodies_id_in) {
  int n_bodies_in = bodies_pos_in.size();
  fprintf(stderr, "n_bodies_in= %d\n", n_bodies_in);
  
  float4 r_min = bodies_pos_in[0];
  float4 r_max = r_min;

  for (int i = 1; i < n_bodies_in; ++i) {
    float4 pos = bodies_pos_in[i];
    r_max.x = max(r_max.x, pos.x);
    r_max.y = max(r_max.y, pos.y);
    r_max.z = max(r_max.z, pos.z);


    r_min.x = min(r_min.x, pos.x);
    r_min.y = min(r_min.y, pos.y);
    r_min.z = min(r_min.z, pos.z);
  }

  pos.x = 0.5 * (r_min.x + r_max.x);
  pos.y = 0.5 * (r_min.y + r_max.y);
  pos.z = 0.5 * (r_min.z + r_max.z);
  pos.w = 0.5 * max(r_max.z - r_min.z, max(r_max.y - r_min.y, r_max.x - r_min.x));

  level = 0;
  
  id_counter() = 0;
  n_leaves()   = 0;

  for(int i = 0; i < n_bodies_in; ++i)  {
//     fprintf(stderr, "insterting %d-th body [%f %f %f %f]\n", i,
// 	    bodies_pos_in[i].x,
// 	    bodies_pos_in[i].y,
// 	    bodies_pos_in[i].z,
// 	    bodies_pos_in[i].w
// 	    );
    insert_body(bodies_pos_in[i], bodies_id_in[i], bodies_pos_in);
  }
  fprintf(stderr, " # of nodes= %ld \n", id_counter());
  fprintf(stderr, " # of level= %ld \n", level_counter());
  
  id = 0;
}

template<int N_BODIES>
void oct_tree<N_BODIES>::destroy_tree() {
  if (n_bodies < 0) {
    n_bodies = 0;
    total_bodies = 0;
  } else {
    n_bodies = 0;
    total_bodies = 0;
    bodies_id.clear();
  }
  if (children != NULL) {
    delete[] children;
    children = NULL;
  }
}

/* --- computes multipole-moments --- */

template<int N_BODIES>
bool oct_tree<N_BODIES>::compute_multipole_moments(vector<float4> &bodies_pos_in, 
						   bool old_tree,
						   bool rebuild_tree) {
  
  com   = (double4){0,0,0,0};
  Qu    = (double4){0,0,0,0};
  Qd    = (double4){0,0,0,0};
  Oct1  = (double4){0,0,0,0};
  Oct2  = (double4){0,0,0,0};
  Oct3  = (double2){0,0};
  r_min = (float3){+1e30, +1e30, +1e30};
  r_max = (float3){-1e30, -1e30, -1e30};

#ifdef OCTUPOLE
  double del[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
#endif
  double S[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};
  double S123 = 0;
  
  /* --- n_bodies < 0, compute com's data from children nodes --- */
  
  if (n_bodies < 0) {
    
    for (int i = 0; i < N_CHILDREN; i++) {
      if (children[i].n_bodies != 0) {
	rebuild_tree = rebuild_tree | children[i].compute_multipole_moments(bodies_pos_in,
									    old_tree,
									    rebuild_tree);

	r_min.x = min(r_min.x, children[i].r_min.x);
	r_min.y = min(r_min.y, children[i].r_min.y);
	r_min.z = min(r_min.z, children[i].r_min.z);

	r_max.x = max(r_max.x, children[i].r_max.x);
	r_max.y = max(r_max.y, children[i].r_max.y);
	r_max.z = max(r_max.z, children[i].r_max.z);

	com.x += children[i].com.w * children[i].com.x;
	com.y += children[i].com.w * children[i].com.y;
	com.z += children[i].com.w * children[i].com.z;
	com.w += children[i].com.w;
      }
    }

    if (com.w > 0) {
      com.x *= 1.0/com.w;
      com.y *= 1.0/com.w;
      com.z *= 1.0/com.w;
    }


    for (int k = 0; k < N_CHILDREN; k++) {
      if (children[k].n_bodies != 0) {

	double4 Ccom = children[k].com;

	Ccom.x -= com.x;
	Ccom.y -= com.y;
	Ccom.z -= com.z;
	
	float ds2 = Ccom.x*Ccom.x + Ccom.y*Ccom.y + Ccom.z*Ccom.z;

#ifdef QUADRUPOLE
	Qu.x += Ccom.w *  3*Ccom.x*Ccom.y        + children[k].Qu.x;
	Qu.y += Ccom.w *  3*Ccom.x*Ccom.z        + children[k].Qu.y;
	Qu.z += Ccom.w *  3*Ccom.y*Ccom.z        + children[k].Qu.z;
	Qd.x += Ccom.w * (3*Ccom.x*Ccom.x - ds2) + children[k].Qd.x;
	Qd.y += Ccom.w * (3*Ccom.y*Ccom.y - ds2) + children[k].Qd.y;
	Qd.z += Ccom.w * (3*Ccom.z*Ccom.z - ds2) + children[k].Qd.z;
#endif

#ifdef OCTUPOLE
	double pos[4] = {Ccom.x, Ccom.y, Ccom.z, Ccom.w};
	double Q[3][3] = {{children[k].Qd.x, children[k].Qu.x, children[k].Qu.y}, 
			  {children[k].Qu.x, children[k].Qd.y, children[k].Qu.z}, 
			  {children[k].Qu.y, children[k].Qu.z, children[k].Qd.z}};
	
	for (int i = 0; i < 3; i++) {
	  for (int j = 0; j < 3; j++) {
	    S[i][j] += pos[3] * (5*(3 - 2*del[i][j])*pos[i]*pos[i] - 3*ds2)*pos[j] +
	      5 * (1-del[i][j]) * pos[i]*Q[i][j] + 2.5*pos[j]*Q[i][i] -
	      (pos[0]*Q[j][0] + pos[1]*Q[j][1] +  pos[2]*Q[j][2]);
	  }
	}
	S[0][0] += children[k].Oct1.x;
	S[0][1] += children[k].Oct1.y;
	S[0][2] += children[k].Oct1.z;
	S[1][0] += children[k].Oct2.x;
	S[1][1] += children[k].Oct2.y;
	S[1][2] += children[k].Oct2.z;
	S[2][0] += children[k].Oct1.w;
	S[2][1] += children[k].Oct2.w;
	S[2][2] += children[k].Oct3.x;
	
	S123 += 15 * (pos[3]*pos[0]*pos[1]*pos[2] + 
		     5.0/3*(pos[0]*Q[1][2] + pos[1]*Q[2][0] * pos[2]*Q[0][1]) + children[k].Oct3.y);
#endif
      }
    }

    Oct1.x = S[0][0];
    Oct1.y = S[0][1];
    Oct1.z = S[0][2];
    Oct2.x = S[1][0];
    Oct2.y = S[1][1];
    Oct2.z = S[1][2];
    Oct1.w = S[2][0];
    Oct2.w = S[2][1];
    Oct3.x = S[2][2];
    Oct3.y = S123;

    /* --- n_bodies > 0, the compute com's data from bodies --- */
  } else {
    
    for (int i = 0; i < n_bodies; i++) {
      float4 body_pos = bodies_pos_in[bodies_id[i]];
      
      r_min.x = min(r_min.x, body_pos.x);
      r_min.y = min(r_min.y, body_pos.y);
      r_min.z = min(r_min.z, body_pos.z);

      r_max.x = max(r_max.x, body_pos.x);
      r_max.y = max(r_max.y, body_pos.y);
      r_max.z = max(r_max.z, body_pos.z);

      com.x += body_pos.w * body_pos.x;
      com.y += body_pos.w * body_pos.y;
      com.z += body_pos.w * body_pos.z;
      com.w += body_pos.w;
    }

    if (com.w > 0) {
      com.x *= 1.0/com.w;
      com.y *= 1.0/com.w;
      com.z *= 1.0/com.w;
    }

    for (int i = 0; i < n_bodies; i++) {
      float4 Bcom = bodies_pos_in[bodies_id[i]];

      Bcom.x -= com.x;
      Bcom.y -= com.y;
      Bcom.z -= com.z;
      
      float ds2 = Bcom.x*Bcom.x + Bcom.y*Bcom.y + Bcom.z*Bcom.z;

#ifdef QUADRUPOLE
      Qu.x += Bcom.w *  3*Bcom.x*Bcom.y;
      Qu.y += Bcom.w *  3*Bcom.x*Bcom.z;
      Qu.z += Bcom.w *  3*Bcom.y*Bcom.z;
      Qd.x += Bcom.w * (3*Bcom.x*Bcom.x - ds2);
      Qd.y += Bcom.w * (3*Bcom.y*Bcom.y - ds2);
      Qd.z += Bcom.w * (3*Bcom.z*Bcom.z - ds2);
#endif

#ifdef OCTUPOLE
      double pos[4] = {Bcom.x, Bcom.y, Bcom.z, Bcom.w};
      
      for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	  S[i][j] += pos[3] * (5*(3 - 2*del[i][j])*pos[i]*pos[i] - 3*ds2)*pos[j];
	}
      }
      S123 += 15 *(pos[3]*pos[0]*pos[1]*pos[2]);
#endif
    }
  }
  
  Oct1.x = S[0][0];
  Oct1.y = S[0][1];
  Oct1.z = S[0][2];
  Oct2.x = S[1][0];
  Oct2.y = S[1][1];
  Oct2.z = S[1][2];
  Oct1.w = S[2][0];
  Oct2.w = S[2][1];
  Oct3.x = S[2][2];
  Oct3.y = S123;
  
  boundary.x = 0.5 * (r_min.x + r_max.x);
  boundary.y = 0.5 * (r_min.y + r_max.y);
  boundary.z = 0.5 * (r_min.z + r_max.z);
  boundary.w = 0.5 * max(r_max.z - r_min.z,
			 max(r_max.y - r_min.y, r_max.x - r_min.x));

  if (old_tree == true) {
    
    float f = fabs((boundary.w - boundary_old.w)/boundary_old.w);
    if (f > 0.5 ) rebuild_tree = true;

//     f = fabs((boundary.x - boundary_old.x)/boundary_old.x);
//     if (f > 0.5) rebuild_tree = true;
    
//     f = fabs((boundary.y - boundary_old.y)/boundary_old.y);
//     if (f > 0.5) rebuild_tree = true;

//     f = fabs((boundary.z - boundary_old.z)/boundary_old.z);
//     if (f > 0.5) rebuild_tree = true;
    
  } else {

    boundary_old = boundary;
    rebuild_tree = false;

  }
    
  return rebuild_tree;
}

/*
  this routine builds a list
   of childrens, leafs * cells 
*/
template<int N_BODIES>
child_struct oct_tree<N_BODIES>::generate_node_list(child_struct child, int octant,
						    vector<int4>                 &children_list,
						    vector<oct_tree<N_BODIES>*>  &node_list) {
  
  /* a parent, so has children then */
  int4 up, dn;
  if (n_bodies < 0) {
    oct_tree<N_BODIES>* child_ptr[N_CHILDREN];
    for (int i = 0; i < N_CHILDREN; i++) {
      child = children[i].generate_node_list(child, i,
					     children_list,
					     node_list);
      switch(i) {
      case 0:  up.x = child.num;  break;
      case 1:  up.y = child.num;  break;
      case 2:  up.z = child.num;  break;
      case 3:  up.w = child.num;  break;
      case 4:  dn.x = child.num;  break;
      case 5:  dn.y = child.num;  break;
      case 6:  dn.z = child.num;  break;
      case 7:  dn.w = child.num;  break;
      }
    
      child_ptr[i] = (oct_tree<N_BODIES>*)child.node_ptr;
    }
    
    children_list.push_back(up);
    children_list.push_back(dn);
    for (int i = 0; i < N_CHILDREN; i++)
      node_list.push_back(child_ptr[i]);
    
    child.num      = children_list.size() - 2;
    child.node_ptr = (void*)this;
    
  } else if (n_bodies > 0) {
    
    child.num      = (1 << (24));
    child.node_ptr = (void*)this;
    
  } else {

    child.num      = 0;
    child.node_ptr = NULL;
    
  }
  
  return child;
}

template<int N_BODIES>
void oct_tree<N_BODIES>::generate_cell_list(int n_crit,
					    vector<oct_tree<N_BODIES>*>  &cell_list) {
  
  /* a parent, so has children then */
  if (total_bodies > n_crit) {
    
    for (int i = 0; i < N_CHILDREN; i++) 
      children[i].generate_cell_list(n_crit, cell_list);
    
  } else if (total_bodies > 0) {
    
    cell_list.push_back(this);
    
  }
}

template<int N_BODIES>
void oct_tree<N_BODIES>::generate_leaf_list(vector<oct_tree<N_BODIES>*>  &leaf_list) {
  
  /* a parent, so has children then */
  if (n_bodies < 0) {
    for (int i = 0; i < N_CHILDREN; i++) 
      children[i].generate_leaf_list(leaf_list);
    
    
  } else if (n_bodies > 0) {
    
    leaf_list.push_back(this);
    
  }
}


#endif // _TREE_MANIP_
