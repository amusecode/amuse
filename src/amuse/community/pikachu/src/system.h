#ifndef SYSTEM_H
#define SYSTEM_H


#ifdef HASH_MAP
#include<unordered_map>
typedef std::unordered_map<int, int> mymap;
#elif TREE_MAP
#include<map>
typedef std::map<int, int> mymap;
#endif //HASH_MAP

#include<iostream>
#include<fstream>
#include<cstdio>
#include<cstdlib>
#include<cstring>

#include"mpi.h"
#include"Matrix3.h"
#include"particle.h"
#include"BHtree.h"
#include"schedule.h"
#include"IO.h"
#include"mpi_interface.h"
#include"distribution.h"

#include"const.h"
#include"global.h"

#define dump_cout(x) cout<<#x<<" = "<<x<<endl;
#define dump_cerr(x) cerr<<#x<<" = "<<x<<endl;

#ifdef SEQUOIA
#include"my_cuda.h"
#endif

#include<omp.h>


class Soft_System{

public:
    Particle* prt; // first address of prt_loc in nbody_system
    int first;
    my_dev::context* devContext;
    //Input
    my_dev::dev_mem<real4>* j_pos_buff;  //Bodies positions
    my_dev::dev_mem<int>*   j_adr_buff;  //Bodies ids
    // if NPROC > 1
    my_dev::dev_mem<real4>* i_pos_buff;  //Bodies positions
    my_dev::dev_mem<int>*   i_adr_buff;  //Bodies ids

    //Output
    my_dev::dev_mem<real4>* i_acc_buff; //Bodies Accelerations
    my_dev::dev_mem<real>*  i_ngh_ds2_buff; //Bodies distance to nearest neighbour squared
    my_dev::dev_mem<int>*   i_ngh_adr_buff;  //Bodies nearest neighbour
    my_dev::dev_mem<int>*   i_Nngh_buff;  //Bodies nearest neighbour
    
    Soft_System();
    
    ~Soft_System(){
        delete[] cell_tree;
        delete[] pos_list_send;
        delete[] pos_list_recv;
        delete[] mass_list_send;
        delete[] mass_list_recv;
        delete[] idx_list_send;
        delete[] idx_list_recv;
    }
    
    void init(Particle prt_loc[]);
    void clear();
    void exchange_LET(Vector3 xlow[],
		    Vector3 xhigh[],
		    int dst_rank[],
		    int src_rank[],
		    Particle BH_glb[]);

    //void setup_tree_adding_LET();
    //void setup_tree_adding_LET_for_sequoia();

    void check_tree(ostream& fout);
    
    void kick(const double &dt);

    void drift(const double &dt);

    void set(const double &_eps2_FS_FS,
	     const double &_eps2_FS_BH,
	     const double &_eps2_BH_BH,
	     const double &_rcut_out_FS_FS,
	     const double &_rcut_out_FS_BH,
	     const double &_rcut_out_BH_BH,
	     const double &_rsearch_FS_FS,
	     const double &_rsearch_FS_BH,  
	     const double &_rsearch_BH_BH, 
	     const double &_theta2,
	     const int &_Ncrit,
	     const int &_Nleaf,
	     const int &_quad_flag);

    void evaluate_gravity_of_FSloc_BHglb_jpara(Particle BH_glb[]);

    void search_neighbour_index(Particle BH_glb[],
				int &ngh_list_len_loc,
				int &ngh_list_len_glb,
				Neighbour_List_Index ngh_list_idx[],
				int &Nshort_loc,
				int &Nshort_glb,
				const int &flag_min_box);

    void evaluate_acc_short_ipara(Neighbour_List_Index ngh_list_idx[],
				  const int &ngh_list_len_loc,
				  mymap &adr_from_idx);

    void dump_cell(ofstream &fout_node);

    void dump(ostream &fout=cout);

#ifdef SEQUOIA
    void build_tree_and_get_force_using_sequoia(char *argv[],
						Vector3 xlow[], 
						Vector3 xhigh[],
						int dst_rank[],
						int src_rank[],
						Particle BH_glb[]);

    void build_tree_on_host_using_sequoia(uint leafNodeIdx[], 
					  uint2 node_bodies[],
					  uint n_children[],
					  float4 boxCenterInfo[],
					  float4 boxSizeInfo[],
					  int n_leafs, 
					  int n_nodes,
					  my_dev::dev_mem<int>  &j_adr_buff);

#endif

    double eps2_FS_FS, rcut_in_FS_FS, rcut_out_FS_FS;

private:
    double rcut_out_FS_BH, rcut_out_BH_BH;
    double rcut_in_FS_BH, rcut_in_BH_BH;
    double eps2_FS_BH, eps2_BH_BH;
    double rsearch2_FS_FS, rsearch2_FS_BH, rsearch2_BH_BH;

    // for tree
    double theta2;
    int Ncrit;
    int Nleaf;
    int quad_flag;

    int heap_remainder;
    Cell_Tree *cell_tree; // array
    Cell_Tree *heap_top; // pointer


    int list_len_send, list_len_recv;
    Vector3 *pos_list_send, *pos_list_recv;
    double *mass_list_send, *mass_list_recv;
    int *idx_list_send, *idx_list_recv;

    void mpi_exchange_interaction_list(const int &box_dest,
				     const int &box_source);

};


class dr2_rank{
public:
  double dr2;
  int rank;
};

class Particle_Comm{
    int index;
    double mass;
    Vector3 pos;
    Vector3 vel;
    int Nj;
public:
    int adr_org;
    void set(Particle prt, const int &_adr){
	index = prt.index;
	mass = prt.mass;
	pos = prt.pos;
	vel = prt.vel;
	Nj = prt.Nj;
	adr_org = _adr;
    }
    void give(Particle_Short &prt){
	prt.index = index;
	prt.mass = mass;
	prt.pos = pos;
	prt.vel = vel;
	prt.Nj = Nj;
    }

#ifdef SEQUOIA
    void calc_soft_forces_using_sequoia(char *argv[]);
#endif

};


class Hard_System{

public:
    Particle_Short *prt_short;
    Hard_System();
    
    ~Hard_System(){
        delete[] prt_short;
        delete[] Tnext;
        delete[] adr_time;
        delete[] prt_comm_send;
        delete[] prt_comm_recv;
    }
    
    void init(); // generate new comunicator for hermite
    void set(const double &_eps2_FS_FS,  
	     const double &_eps2_FS_BH,  
	     const double &_eps2_BH_BH,  
	     const double &_rcut_out_FS_FS,  
	     const double &_rcut_out_FS_BH,  
	     const double &_rcut_out_BH_BH, 
	     const double &_eta_s,
	     const double &_eta_FS,
	     const double &_eta_BH,
	     const double &_mass_min);

    void merge_ngh_list(Particle prt_loc[],
			Neighbour_List_Index ngh_list_idx[],
			mymap adr_from_idx,
			const int &Nshort_loc,
			int &Nshort_glb,
			const int &ngh_list_len_loc,
			int &ngh_list_len_glb,
			Neighbour_List ngh_list[]);

    int hermite4_para(const double &Tsys,
		      const double &dt_glb,
		      const int &Nshort_glb);

    void copy_prt_short_to_prt_loc(Particle prt_loc[],
				   Particle BH_glb[]);

    void dump_prt(const int &id);

    void dump(ostream &fout=cout);

private:
    double eps2_FS_FS, eps2_FS_BH, eps2_BH_BH;
    double rcut_out_FS_FS, rcut_out_FS_BH, rcut_out_BH_BH;
    double rcut_in_FS_FS, rcut_in_FS_BH, rcut_in_BH_BH;
    double eta_s, eta_FS, eta_BH;

    // for mpi
    MPI_Comm comm_new_array[NPROC_MAX+1]; // comm_new_array[0] is null
    int myrank_new_array[NPROC_MAX+1]; // myrank_new_array[0] is null

    double *Tnext;
    int *adr_time;

    // for communication
    Particle_Comm *prt_comm_send;
    Particle_Comm *prt_comm_recv;
    int Nshort_disp[NPROC_MAX+1];

    double acc_offset_sq; // for timestep

    void generate_new_communicator();

    int evaluate_acc_jrk_short_jpara(const int &Ni,
				     const double &Tnext,
				     const int &mode);

    int evaluate_acc_jrk_short_jpara_sse(const int &Ni,
					 const double &Tnext,
					 const int &mode);

    int evaluate_acc_jrk_short_sirial(const int &Ni,
				      const double &Tnext,
				      const int &mode);

    int hermite4_2body(Particle_Short *prt_i,
		       Particle_Short *prt_j,
		       const double &Tsys,
		       const double &dt_glb);



    void calc_acc_jrk_duncan4_from_jarray_sse_mix(Particle_Short* i_prt,
						  const int& j_head,
						  const int& j_tale,
						  const double& Tnext,
						  const int& mode);

    void calc_acc_jrk_duncan4_from_jarray_sse_fulld(Particle_Short* i_prt,
						    const int& j_head,
						    const int& j_tale,
						    const double& Tnext,
						    const int& mode);

    void calc_acc_jrk_duncan4_from_jarray_scholar(Particle_Short* i_prt,
						  const int& j_head,
						  const int& j_tale,
						  const double& Tnext,
						  const int& mode);

};


class Nbody_System{
public:
    double Tsys;
    double Tend;
    double dt_glb;

    Nbody_System();
    
    ~Nbody_System(){
        delete[] prt_loc;
        delete[] BH_glb;
        delete[] prt_dead_glb;
        delete[] ngh_list_idx;
        delete[] ngh_list;
    }

    void read_file(char param_file[],
		   char output_dir[]);
    
    void divide_particles();

    void calc_soft_forces(char *argv[]);

    void create_ngh_list_glb();

    void initialize_division();

    void calc_energy(double &Ek,
		     double &Ep,
		     double &E);

    void calc_energy(const int &mode);
    void kick_half();
    void kick_full();
    void merge_ngh_list();
    void write_file(char output_dir[]);
    int evolve(char output_dir[],
	       char *argv[],
	       int &after_stellar_evolution);
    void dump_prt();
    void dump_ngh_list();
    void dump_prt_short();
    void dump_tree_strcture();


    Soft_System soft_system;
    Hard_System hard_system;

    Particle *prt_loc;
    Particle *BH_glb;
    mymap adr_from_idx;


    Particle *prt_dead_glb;
    int snp_id;
    double snp_interval;
    double Ek0, Ep0, E0, Emerge0, Egr0;
    double Ek1, Ep1, E1, Emerge1, Egr1;

    double pos_max;
    double mass_min;

private:

    int Nshort_loc;
    int Nshort_glb;
    int ngh_list_len_loc;
    int ngh_list_len_glb;

    Neighbour_List *ngh_list;
    Neighbour_List_Index *ngh_list_idx;


    // for MPI
    int dst_rank[NPROC_MAX];
    int src_rank[NPROC_MAX];
    int npdim[3];
    int sample_freq;
    Vector3 xlow[NPROC_MAX]; //front, left and bottome vertex of domain
    Vector3 xhigh[NPROC_MAX]; //back, right and top vertex of domain

    // IO stream
    ofstream fout_tcal;
    ofstream fout_energy;
    ofstream fout_BH;

    void evolve_hard_system(char *argv[]);
    void evolve_soft_system(const double &dt_kick);
};




#endif //SYSTEM_H

