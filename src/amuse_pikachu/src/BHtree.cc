#include<iostream>
#include"BHtree.h"

static const int LETINDEXBASE = (int)2e7; // for local essential tree
static const int SPINDEXBASE = (int)3e7; // for super particle
static const int BPWORKMAX = 16384; // more than Nleaf


static const int INT_LIST_MAX = (int)1e6;
static Vector3 pos_list[INT_LIST_MAX];
static double mass_list[INT_LIST_MAX];
static Particle *(prt_list[INT_LIST_MAX]);
static int index_list[INT_LIST_MAX];


void dev_open(){};
void dev_close(){};

//////////////////////////////////////////
//
//       make interaction list
//
//////////////////////////////////////////

// if MPI is used, 
// the particles in other node possibly is included in this box.
// fuction "pack_real_particles()" exclude these particles.
// these particles must be integrated in other node.
inline void pack_real_particles(int &Ni,  const int &first_leaf){
    int Ninew = 0;
    for(int i=first_leaf; i<first_leaf+Ni; i++){
	if(prt_list[i]->node_org == MYRANK){
	    if(i > Ninew+first_leaf){
		// swap location i and ninew
		Particle *tmp = prt_list[i];
		prt_list[i] = prt_list[Ninew+first_leaf];
		prt_list[Ninew+first_leaf] = tmp;
	    }
	    Ninew++;
	}
    }
    if(Ninew < Ni){
	for(int i=first_leaf; i<first_leaf+Ni; i++){
	    pos_list[i] = prt_list[i]->pos;
	    mass_list[i] = prt_list[i]->mass;
	    index_list[i] = prt_list[i]->index;
	}
	Ni = Ninew;
    }
}

//////////////////////////////////////////
//
//       add particles or cells to LET
//
//////////////////////////////////////////

void Cell_Tree::add_particles_to_essential_tree(Vector3 pos_list_send[], 
                                                double mass_list_send[], 
                                                int &list_len, 
                                                const int &LIST_COMM_MAX,
                                                int index_list_send[]){
    if (isleaf){
        Particle *prt_tmp = prt_first;
        for(int i=0; i<Nprt; i++){
            pos_list_send[list_len] = prt_tmp->pos;
            mass_list_send[list_len] = prt_tmp->mass;
            index_list_send[list_len] = prt_tmp->index;
            list_len ++;

            if(list_len > INT_LIST_MAX){
                cerr <<"myrank="<<MYRANK<<endl;
                cerr <<"List len exceeded"<<endl;
                cerr <<"list_len="<<list_len<<endl;
                cerr <<"INT_LIST_MAX="<<INT_LIST_MAX<<endl;
                halt_program();
            }
            prt_tmp = prt_tmp->prt_next;
        }
    }
    else{
        cerr << "add_particles_to_essential_tree intenal error: non-leaf\n";
        halt_program();
    }
}


void Cell_Tree::add_to_essential_tree(const Vector3 &box_dest_pos, // center of pos in box_dest
				      const Vector3 &box_dest_length, // len of box_dest
				      const double &theta2,
				      Vector3 pos_list_send[],
				      double mass_list_send[],
				      int &list_len,
				      const int &LIST_COMM_MAX,
				      const double &rsearch2_FS_FS,
				      int index_list_send[],
				      const int &quad_flag,
				      Particle BH_dst[],
				      const int &NBH_dst, 
				      const double &rsearch2_FS_BH){ 
    // j particles will be scattered
    // own node is considerd as j particles
    double rmin2_box_box = shortest_separation_squared_box_box(pos_center, size, 
							       box_dest_pos, 
							       box_dest_length);
    double rmin2_box_BH = LARGEDOUBLE;
    for(int dst=0; dst<NBH_dst; dst++){
        double rtmp2 = shortest_separation_squared_point_box(pos_center, size,
							     BH_dst[dst].pos);
        if(rtmp2 < rmin2_box_BH){
            rmin2_box_BH = rtmp2;
        }
    }
    if( rmin2_box_box*theta2 > size*size
	&& ! are_overlapped_box_box(pos_center, size, box_dest_pos, box_dest_length) 
	&& rmin2_box_box > rsearch2_FS_FS
        && rmin2_box_BH > rsearch2_FS_BH){
        // node and box are well separated;
        if(!quad_flag){
            // dipole approximation
            pos_list_send[list_len] = cmpos;
            mass_list_send[list_len] = cmmass;
            index_list_send[list_len] = list_len + LETINDEXBASE;
            list_len ++;
        }
        else{
            // quadrupole approximation
            for(int k=0; k<3; k++){
                if(p2mass[k] == 0.0){break;} // don't forget
                pos_list_send[list_len] = p2pos[k];
                mass_list_send[list_len] = p2mass[k];
                index_list_send[list_len] = list_len + LETINDEXBASE;
                list_len ++;
            }
        }
        if(list_len >= LIST_COMM_MAX){
            cerr <<"myrank="<<MYRANK<<endl;
            cerr <<"LIST_COMM_MAX is too small (in function add_to_essential_tree()))"<<endl;
            cerr <<"list_len="<<list_len<<endl;
            cerr <<"LIST_COMM_MAX="<<LIST_COMM_MAX<<endl;
            halt_program();
        }
    }
    else{
        if(isleaf){
            add_particles_to_essential_tree(pos_list_send,  mass_list_send,  
                                            list_len,  LIST_COMM_MAX,  
                                            index_list_send);
        }
        else{
            for(int i=0;i<8;i++){
                if (child[i] != NULL){
                    child[i]->add_to_essential_tree(box_dest_pos,  box_dest_length,
						    theta2,  
						    pos_list_send,  mass_list_send,
						    list_len,  LIST_COMM_MAX,
						    rsearch2_FS_FS,  
						    index_list_send,  
						    quad_flag,
						    BH_dst,  NBH_dst,  rsearch2_FS_BH);
                }
            }
        }
    }
}






////////////////////////////////////
//
//      search neighbours
//
////////////////////////////////////

void Cell_Tree::search_neighbour_index(Particle &prti,
                                       const double &r2_ngh, 
                                       Neighbour_List_Index ngh_list_idx[],
                                       int &ngh_list_len,
                                       const int& ngh_list_max,
				       const int& flag_min_box = 0){
    double separation_prti_node_squared;
    if(flag_min_box){
	separation_prti_node_squared = shortest_separation_squared_point_box(pos_center, size, prti.pos);
    }
    else{
	separation_prti_node_squared = shortest_separation_squared_point_box(pos_center, size, prti.pos);
    }
    if(isleaf){
        if( separation_prti_node_squared < r2_ngh
            || separation_prti_node_squared == 0.0 ){
            Particle *prt_tree_j = prt_first;
            for(int j=0; j<Nprt; j++){
                Particle *prtj = prt_tree_j;
                Vector3 rijtmp = prti.pos - prtj->pos;
                double separation_prti_prtj_squared=rijtmp*rijtmp;
                if( separation_prti_prtj_squared < r2_ngh 
                    && prti.index != prtj->index
                    && prtj->index < NFS_GLB_ORG+NBH_GLB_ORG){
                    if(ngh_list_len >= ngh_list_max){
                        cerr<<"ERROR: neighbour list len is too small in function search_neighbour"<<endl;
                        cerr<<"ngh_list_len="<<ngh_list_len<<endl;
                        cerr<<"ngh_list_max="<<ngh_list_max<<endl;
                        halt_program();
                    }
                    ngh_list_idx[ngh_list_len].idx_i = prti.index;
                    ngh_list_idx[ngh_list_len].idx_j = prtj->index;
                    prti.have_ngh = 1;
                    prtj->have_ngh = 1;
                    prti.Nj++;
                    ngh_list_len++;
                }
                prt_tree_j = prt_tree_j->prt_next;
            }
        }
    }
    else if ( separation_prti_node_squared == 0.0 
              || separation_prti_node_squared < r2_ngh){
        for(int i=0; i<8; i++){
            if(child[i] != NULL){
                child[i]->search_neighbour_index(prti,  r2_ngh,  ngh_list_idx,  ngh_list_len, ngh_list_max, flag_min_box);
            }
        }
    }
}


