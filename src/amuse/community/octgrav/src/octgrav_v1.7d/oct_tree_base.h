#ifndef _OCT_TREE_BASE_
#define _OCT_TREE_BASE_

#define N_CHILDREN 8

#include "octgravdefs.h"

struct child_struct {
  int   num;
  void* node_ptr;
};


/* ----------------------------------------------------------------------- */

template<int N_BODIES = 32>
class oct_tree {
protected:

  /* ---  array with children, allocated only for a cell --- */
  
  oct_tree<N_BODIES> *children;
  
  /* --- pointers to bodies, allocated only for a leaf --- */
  
  vector<int> bodies_id;
  int n_bodies;                   // < 0, if cell, > 0 if leaf & = 0 if empty
  
  /* --- total number of bodies contained in all descending cell ---- */
  
  int total_bodies;     
  
  /* ---- node id --- */
  
  int  level;
  long int id;
  
  /* --- half-length of the cell's edge and location of the geometrical centre --- */

  float4 pos;

  float3 r_min, r_max;
  float4 boundary;
  float4 boundary_old;
  
  /* --- multipole moments of the cell --- */

  double4 com;
  double4 Qu;     /* 12, 13, 23 */
  double4 Qd;     /* 11, 22, 33 */


  double4 Oct1;   /* S11, S12, S13, S31 */
  double4 Oct2;   /* S21, S22, S23, S32 */
  double2 Oct3;   /* S33, S123 */

  
  long int& n_leaves() {
    static long int __n_leaves;
    return __n_leaves;
  }

  long int& id_counter() {
    static long int __id_counter;
    return __id_counter;
  }

  long int& level_counter() {
    static long int __level;
    return __level;
  }

public:
  
  int offset;
  
  /* --- constructor/destructor --- */
  
  explicit oct_tree() : children(NULL),
			n_bodies(0),
			total_bodies(0),
			offset(-1)
  {
    id = ++id_counter();
    bodies_id.resize(N_BODIES);
  }
  
  ~oct_tree() {
    if (children != NULL) {
      delete[] children;
      children = NULL;
    }
    n_bodies = 0;
    total_bodies = 0;
    bodies_id.clear();
  }

  /* ----- accessor methods ----- */

  int          get_n_nodes()    {return id_counter();}
  int          get_n_levels()   {return level_counter();}
  int          get_n_leaves()   {return n_leaves();}
  oct_tree&    get_child(int i) {return children[i];}
  vector<int>& get_bodies_id()  {return bodies_id;}
  
  int          get_n_bodies()     {return n_bodies;}
  unsigned int get_total_bodies() {return total_bodies;}
  int          get_id()           {return id;}
 
  float4 get_boundary() {return boundary;}
  float4 get_pos() {return pos;}
  float4 get_com() {return (float4){com.x, com.y, com.z, com.w};}
  float4 get_Qu()  {return (float4){Qu.x, Qu.y, Qu.z, Qu.w};} 
  float4 get_Qd()  {return (float4){Qd.x, Qd.y, Qd.z, Qd.w};} 
  float4 get_Oct1() {return (float4){Oct1.x, Oct1.y, Oct1.z, Oct1.w};}
  float4 get_Oct2() {return (float4){Oct2.x, Oct2.y, Oct2.z, Oct2.w};}
  float2 get_Oct3() {return (float2){Oct3.x, Oct3.y};}
 
  /* --- tree management methods --- */
  
//   oct_tree<N_BODIES>* get_node(float4);
  // void get_pblocks(int  node_id,  real node_range,
// 		   vec& body_pos, 
// 		   vector<vec>& bodies_pos,
// 		   int pblocks[], int n[]);

  void insert_body(float4, int, vector<float4>&);
  void remove_body(float4, int);
  void build_tree(vector<float4>&, vector<int>&);
  void destroy_tree();
  bool compute_multipole_moments(vector<float4>&, bool, bool);
  child_struct generate_node_list(child_struct, int,
				  vector<int4>&,
				  vector<oct_tree<N_BODIES>*>&);
  void generate_cell_list(int, vector<oct_tree<N_BODIES>*>&);
  void generate_leaf_list(vector<oct_tree<N_BODIES>*>&);
  
  
  /* --- obtain octant of a given locaiton wrt geometrical cell's centre --- */
  
  inline int octant(float4 r) {
    return int(((r.x < pos.x) << (0)) + 
	       ((r.y < pos.y) << (1)) +
	       ((r.z < pos.z) << (2)));
  }
  
  /* --- computes child's geometrical centre --- */
  
  inline float4 child_pos(int oct) { 
    return (float4){pos.x - ( ((oct & 1) << (1)) - 1) * 0.5 * pos.w,
    	            pos.y - ( ((oct & 2) << (0)) - 1) * 0.5 * pos.w,
                    pos.z - ( ((oct & 4) >> (1)) - 1) * 0.5 * pos.w,
	            0.5 * pos.w};
  }
  
};

#endif // _OCT_TREE_BASE_
