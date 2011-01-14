#ifndef  BHTREE_H
#  define   BHTREE_H
/*-----------------------------------------------------------------------------
 *  BHtree : basic class for C++ implementation of BH treecode
 *  J. Makino 1998/12/14
 *-----------------------------------------------------------------------------
 */

#include "stdinc.h"
#include "vec.h"
#include "nbody_particle.h"




// const int default_key_length = 15;	// 64-bit machines
const int default_key_length = 10;	// 32-bit machines





const int default_ix_offset = 1<<(default_key_length-1);

typedef nbody_particle real_particle;
typedef nbody_system real_system;

typedef long BHlong;

class bhparticle
{
private:
    real_particle * rp;
    BHlong key;
    
public:
    bhparticle(){
	rp = NULL;
	key = 0;
    }
    void set_key(real rscale, int ioffset, int keybits);
    int friend compare_key( bhparticle * p1,  bhparticle * p2);
    BHlong get_key(){return key;}
    void set_rp(real_particle * p){rp = p;}
    real_particle * get_rp(){return rp ;}
    void friend sort_bh_array( bhparticle * r, int lo, int up );
    
};


class bhnode
{
private:
    vec pos;
    real l;
    bhnode * child[8];
    bhparticle * bpfirst;
    int nparticle;
    int isleaf;
#ifdef SPH    
    real hmax_for_sph;
#endif    
    vec cmpos;
    real cmmass;
    
public:
    bhnode(){
	pos = 0.0;
	l = 0.0;
	for(int i = 0; i<8;i++)child[i] = NULL;
	bpfirst = NULL;
	nparticle = 0;
	isleaf = 1;
#ifdef SPH	
	hmax_for_sph = 0;
#endif	
	cmpos = 0.0;
	cmmass = 0.0;
    }
    void clear(){
	pos = 0.0;
	l = 0.0;
	for(int i = 0; i<8;i++)child[i] = NULL;
	bpfirst = NULL;
	nparticle = 0;
	isleaf = 1;
#ifdef SPH	
	hmax_for_sph = 0;
#endif	
	cmpos = 0.0;
	cmmass = 0.0;
    }
    void set_pos(vec newpos){pos = newpos;}
    vec get_pos(){return pos;}
    void set_length(real newl){l = newl;}
    real get_length(){return l;}
    void create_tree_recursive(bhnode * & heap_top, int & heap_remainder,
			       BHlong current_key,
			       int current_level,
			       int n_critical);
    void assign_root(vec root_pos, real length, bhparticle * bp, int nparticle);
    void dump(int indent);
    int sanity_check();
#ifdef SPH    
    void set_hmax_for_sph();
    real get_hmax_for_sph(){return hmax_for_sph;}
#endif    
    int friend check_and_set_nbl(bhnode * p1,bhnode * p2);
    void set_cm_quantities();
    void accumulate_force_from_tree(vec & ipos, real eps2, real theta2,
				   vec & acc,
				   real & phi);
    void add_to_interaction_list(bhnode & dest_node, real theta2,
				 vec * pos_list,
				 real * mass_list,
				 int & nlist,
				 int list_max,
				 int & first_leaf);
    void evaluate_gravity_using_tree_and_list(bhnode & source_node,
					      real theta2,
					      real eps2,
					      int ncrit);
};


void clear_tree_counters();
void print_tree_counters();
#endif
