#include"system.h"
#include"const.h"



Soft_System::Soft_System(){
    cerr<<"Soft_System constructor"<<endl;

    cell_tree = new Cell_Tree[NCELL_TREE_LOC_MAX];

    pos_list_send = new Vector3[LIST_COMM_MAX];
    pos_list_recv = new Vector3[LIST_COMM_MAX];

    mass_list_send = new double[LIST_COMM_MAX];
    mass_list_recv = new double[LIST_COMM_MAX];

    idx_list_send = new int[LIST_COMM_MAX];
    idx_list_recv = new int[LIST_COMM_MAX];

    eps2_FS_FS = eps2_FS_BH = eps2_BH_BH = 0.0;
    rcut_out_FS_FS = rcut_out_FS_BH = rcut_out_BH_BH = 0.0;
    rcut_in_FS_FS = rcut_in_FS_BH = rcut_in_BH_BH = 0.0;
    rsearch2_FS_FS = rsearch2_FS_BH = rsearch2_BH_BH = 0.0;
    theta2 = 0.0;
    Ncrit = 8;
    Nleaf = 8;
    quad_flag = 0;
    heap_top = cell_tree;
    heap_remainder = NCELL_TREE_LOC_MAX;
    list_len_send = list_len_recv = 0;
    cerr<<"Soft_System constructed"<<endl;
}

void Soft_System::init(Particle prt_loc[]){
    prt = prt_loc;
}

void Soft_System::clear(){
    for(int i=0; i<NCELL_TREE_LOC_MAX - heap_remainder; i++){
        cell_tree[i].clear();
    }
    heap_top = cell_tree;
    heap_remainder = NCELL_TREE_LOC_MAX;
    //NALL_LOC = 0;
}

void Soft_System::set(const double &_eps2_FS_FS,  
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
		      const int &_quad_flag){
    eps2_FS_FS = _eps2_FS_FS; 
    eps2_FS_BH = _eps2_FS_BH; 
    eps2_BH_BH = _eps2_BH_BH; 
    rcut_out_FS_FS = _rcut_out_FS_FS;
    rcut_out_FS_BH = _rcut_out_FS_BH;
    rcut_out_BH_BH = _rcut_out_BH_BH;
    rcut_in_FS_FS = rcut_out_FS_FS * RCUT_IN_FACTOR;
    rcut_in_FS_BH = rcut_out_FS_BH * RCUT_IN_FACTOR;
    rcut_in_BH_BH = rcut_out_BH_BH * RCUT_IN_FACTOR;

    theta2 = _theta2;
    Ncrit = _Ncrit;
    Nleaf = _Nleaf;
    quad_flag = _quad_flag;

    rsearch2_FS_FS = _rsearch_FS_FS * _rsearch_FS_FS;
    rsearch2_FS_BH = _rsearch_FS_BH * _rsearch_FS_BH;
    rsearch2_BH_BH = _rsearch_BH_BH * _rsearch_BH_BH;
}

/////////////////////////////
//                         //  
//           DUMP          //
//                         //
/////////////////////////////

void Soft_System::dump(ostream &fout){

    fout<<"dump soft system"<<endl;
    fout<<"rcut_out_FS_FS="<<rcut_out_FS_FS<<endl;
    fout<<"rcut_out_FS_BH="<<rcut_out_FS_BH<<endl;
    fout<<"rcut_out_BH_BH="<<rcut_out_BH_BH<<endl;
    fout<<"rcut_in_FS_FS="<<rcut_in_FS_FS<<endl;
    fout<<"rcut_in_FS_BH="<<rcut_in_FS_BH<<endl;
    fout<<"rcut_in_BH_BH="<<rcut_in_BH_BH<<endl;
    fout<<"eps2_FS_FS="<<eps2_FS_FS<<endl;
    fout<<"eps2_FS_BH="<<eps2_FS_BH<<endl;
    fout<<"eps2_BH_BH="<<eps2_BH_BH<<endl;
    fout<<"rsearch2_FS_FS="<<rsearch2_FS_FS<<endl;
    fout<<"rsearch2_FS_BH="<<rsearch2_FS_BH<<endl;
    fout<<"rsearch2_BH_BH="<<rsearch2_BH_BH<<endl;
    fout<<"theta2="<<theta2<<endl;
    fout<<"Ncrit="<<Ncrit<<endl;
    fout<<"Nleaf="<<Nleaf<<endl;
    fout<<"quad_flag="<<quad_flag<<endl;
}


/////////////////////////////////////
//                                 //    
//           KICK & DRIFT          //
//                                 //    
/////////////////////////////////////

void Soft_System::kick(const double &dt){
    for(int i=0; i<NFS_LOC+NBH_LOC; i++){
        prt[i].kick_for_tree(dt);
    }
}

void Soft_System::drift(const double &dt){
    for(int i=0; i<NFS_LOC+NBH_LOC; i++){
	if(!prt[i].have_ngh){
	    prt[i].drift(dt);
	}
    }
}


/////////////////////////////////////
//                                 //    
//            SETUP TREE           //
//                                 //    
/////////////////////////////////////


void Soft_System::mpi_exchange_interaction_list(const int &box_dest,
						const int &box_source){
    //first send and get the number of particles to send and get
    MPI_Status status;
    MPI_Sendrecv(&list_len_send, 1, MPI_INT, box_dest, MYRANK*10,
                 &list_len_recv, 1, MPI_INT, box_source, box_source*10, 
                 MPI_COMM_WORLD, &status);

    if (list_len_recv >= LIST_COMM_MAX){
        cerr<<"myrank = "<<MYRANK<<": mpi_exchange_interaction_list: LIST_COMM_MAX is too small"<<endl;
        cerr<<"LIST_COMM_MAX="<<LIST_COMM_MAX<<", list_len_recv="<<list_len_recv<<endl;
        halt_program();
    }

    MPI_Sendrecv(pos_list_send, list_len_send*3, MPI_DOUBLE, box_dest, MYRANK*10+1,
                 pos_list_recv, list_len_recv*3, MPI_DOUBLE, box_source, box_source*10+1,
                 MPI_COMM_WORLD, &status);

    MPI_Sendrecv(mass_list_send, list_len_send, MPI_DOUBLE, box_dest, MYRANK*10+2,
                 mass_list_recv, list_len_recv, MPI_DOUBLE, box_source, box_source*10+2,
                 MPI_COMM_WORLD, &status);

    MPI_Sendrecv(idx_list_send, list_len_send, MPI_INT, box_dest, MYRANK*10+3,
                 idx_list_recv, list_len_recv, MPI_INT, box_source, box_source*10+3,
                 MPI_COMM_WORLD, &status);
} 


// set prt_loc
void Soft_System::exchange_LET(Vector3 xlow[],
			       Vector3 xhigh[],
			       int dst_rank[],
			       int src_rank[],
			       Particle BH_glb[]){
    static Particle BH_dest[NBH_GLB_MAX];
    int totalsent = 0;
    int i_loc = NFS_LOC + NBH_LOC;
    for(int ib=0; ib<NPROC-1; ib++){
        int box_dest = dst_rank[ib];
        int box_source = src_rank[ib];

        // construct LET
        list_len_send = 0;
        Vector3 box_dest_pos = 0.5*(xhigh[box_dest] + xlow[box_dest]);
        Vector3 box_dest_length = xhigh[box_dest] - xlow[box_dest];

        for(int ibh=0; ibh<NBH_LOC; ibh++){
            if(shortest_separation_squared_point_box(box_dest_pos, box_dest_length, prt[ibh].pos) <= rsearch2_FS_BH){
                // for BHs (in thisnode) - boundary of dest node
                pos_list_send[list_len_send] = prt[ibh].pos;
                mass_list_send[list_len_send] = prt[ibh].mass;
                idx_list_send[list_len_send] = prt[ibh].index;
                list_len_send++;
            }
            else{
                for(int jbh=0; jbh<NBH_GLB; jbh++){
                    if(BH_glb[jbh].node_org == box_dest 
                       && (BH_glb[jbh].pos - prt[ibh].pos)*(BH_glb[jbh].pos - prt[ibh].pos) <= rsearch2_BH_BH){
                        // for BHs (in thisnode) - BH of dest node (using BH_glb)
                        pos_list_send[list_len_send] = prt[ibh].pos;
                        mass_list_send[list_len_send] = prt[ibh].mass;
                        idx_list_send[list_len_send] = prt[ibh].index;
                        list_len_send++;
                    }
                }
            }
        }

        // for FSs (in thisnode) - boundary of dest node and BHs
        int NBH_dest = 0;
        for(int jbh=0; jbh<NBH_GLB; jbh++){
            if(BH_glb[jbh].node_org == box_dest){
                BH_dest[NBH_dest] = BH_glb[jbh];
                NBH_dest++;
            }
        }
        cell_tree->add_to_essential_tree(box_dest_pos,  box_dest_length,
					 theta2,  
					 pos_list_send,  mass_list_send,  
					 list_len_send,  LIST_COMM_MAX,
					 rsearch2_FS_FS,
					 idx_list_send,
					 quad_flag,   
					 BH_dest,  NBH_dest,  rsearch2_FS_BH);

        list_len_recv=0;
        totalsent += list_len_send;
	
        mpi_exchange_interaction_list(box_dest, box_source);

        if(i_loc + list_len_recv >= NPRT_LOC_MAX){
            cerr<<"NPR_LOC_MAX is too small (in function exchange_loc_essential_trees)"<<endl;
            cerr<<"i_loc="<<i_loc<<",  list_len_recv="<<list_len_recv<<endl;
            cerr<<"NPRT_LOC_MAX="<<NPRT_LOC_MAX<<endl;
            halt_program();
        }

        for(int i = 0; i<list_len_recv; i++){
            prt[i_loc].pos = pos_list_recv[i];
            prt[i_loc].mass = mass_list_recv[i];
            prt[i_loc].index = idx_list_recv[i];
            prt[i_loc].node_org = box_source;
            prt[i_loc].Nj = 0;
            prt[i_loc].vel = 0.0;
            i_loc++;
        }
    }

    for(int i=0; i<NBH_GLB_ORG; i++){
        if(BH_glb[i].node_org == MYRANK || isdead(BH_glb[i]) ){continue;}
        prt[i_loc].pos = BH_glb[i].pos;
        prt[i_loc].mass = BH_glb[i].mass;
        prt[i_loc].index = BH_glb[i].index;
        prt[i_loc].node_org = BH_glb[i].node_org;
        prt[i_loc].Nj = 0;
        i_loc++;
    }
    NALL_LOC = i_loc;
    int sent_loc = totalsent;
    MPI_Reduce(&sent_loc, &totalsent, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
    if (MYRANK == 0){
        cerr << "Exchanged treenodes = " << totalsent << endl;
    }
}


void Soft_System::check_tree(ostream& fout = cout){
    cell_tree[0].sanity_check(fout);
}


///////////////////////////////////////////////////////
//                                                   //    
//            LONG RANGE FORCE CALCULATION           //
//                                                   //    
///////////////////////////////////////////////////////

void Soft_System::evaluate_gravity_of_FSloc_BHglb_jpara(Particle BH_glb[]){

    // i: all BH particles
    // j: all local FS & BH particles
    // 1) evaluate forces between all BH <-> local FS and count up BH.Nj  (FS only)
    // 2) using all reduce, evaluate forces (all BH <- all FS) and BH.Nj (j para)
    // 3) evaluate all BH <- all BH
    // NOTE, prt.Nj is not changed. Nj is still be the number of prt
    // in this function BH_glb[] is used insted of prt_loc, finally the imformation of BH_glb are sent to prt_loc

    static Vector3 acc_list_send[NBH_GLB_MAX];
    static Vector3 acc_list_recv[NBH_GLB_MAX];
    static double pot_list_send[NBH_GLB_MAX];
    static double pot_list_recv[NBH_GLB_MAX];
    static int Nj_list_send[NBH_GLB_MAX];
    static int Nj_list_recv[NBH_GLB_MAX];
    static dr2_rank dr2_rank_FS_send[NBH_GLB_MAX];
    static dr2_rank dr2_rank_FS_recv[NBH_GLB_MAX];
    double r2_tmp;
    for(int i=0; i<NBH_GLB_ORG; i++){
	acc_list_send[i] = 0.0;
	pot_list_send[i] = 0.0;
	Nj_list_send[i] = 0;
	dr2_rank_FS_send[i].rank = MYRANK;
	dr2_rank_FS_send[i].dr2 = 0.0;
	BH_glb[i].r2_ngh_FS = LARGEDOUBLE;
	BH_glb[i].idx_ngh_FS = -1;
	BH_glb[i].r2_ngh_BH = LARGEDOUBLE;
	BH_glb[i].idx_ngh_BH = -1;
	BH_glb[i].r2_ngh_BH = LARGEDOUBLE;
	BH_glb[i].Nj = 0;
    }
    for(int i=NBH_LOC; i<NFS_LOC+NBH_LOC; i++){
	prt[i].r2_ngh_BH = LARGEDOUBLE;
	prt[i].idx_ngh_BH = LARGEDOUBLE;
    }

    // 1) evaluate forces between all BH <-> local FS
    for(int i=0; i<NBH_GLB_ORG; i++){
	Particle *prti = BH_glb+i;
	if(isdead(*prti)){continue;}
	for(int j=NBH_LOC; j<NFS_LOC+NBH_LOC; j++){
	    Particle *prtj = prt+j;
	    if(prti->index == prtj->index){continue;}
	    pairwise_acc(prti->mass,  prti->pos,  acc_list_send[i],  pot_list_send[i],
			 prtj->mass,  prtj->pos,  prtj->acc,  prtj->pot,
			 eps2_FS_BH,  r2_tmp);
	    if(r2_tmp < prti->r2_ngh_FS){
		prti->r2_ngh_FS = r2_tmp;
		prti->idx_ngh_FS = prtj->index;
		dr2_rank_FS_send[i].dr2 = r2_tmp;
	    }
	    if(r2_tmp < prtj->r2_ngh_BH){
		prtj->r2_ngh_BH = r2_tmp;
		prtj->idx_ngh_BH = prti->index;
	    }
	    if(r2_tmp < rsearch2_FS_BH)
		Nj_list_send[i]++;
	}
    }
    MPI_Allreduce(acc_list_send, acc_list_recv, NBH_GLB*3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(pot_list_send, pot_list_recv, NBH_GLB, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(Nj_list_send, Nj_list_recv, NBH_GLB, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(dr2_rank_FS_send, dr2_rank_FS_recv, NBH_GLB,  MPI_DOUBLE_INT,  MPI_MINLOC, MPI_COMM_WORLD);
    for(int i=0; i<NBH_GLB_ORG; i++){
	BH_glb[i].acc += acc_list_recv[i];
	BH_glb[i].pot += pot_list_recv[i];
	BH_glb[i].Nj += Nj_list_recv[i];
	BH_glb[i].r2_ngh_FS = dr2_rank_FS_recv[i].dr2;
	MPI_Bcast(&(BH_glb[i].idx_ngh_FS),  1,  MPI_INT,  dr2_rank_FS_recv[i].rank,  MPI_COMM_WORLD);
    }


    // 3) evaluate all BH <- all BH
    for(int i=0; i<NBH_GLB_ORG; i++){
	Particle *prti = BH_glb+i;
	for(int j=i+1; j<NBH_GLB_ORG; j++){
	    Particle *prtj = BH_glb+j;
	    if(prti->index == prtj->index){continue;}
	    pairwise_acc(prti->mass,  prti->pos,  prti->acc,  prti->pot,
			 prtj->mass,  prtj->pos,  prtj->acc,  prtj->pot,
			 eps2_BH_BH,  r2_tmp);
	    if(r2_tmp < prti->r2_ngh_BH){ 
		prti->r2_ngh_BH = r2_tmp;
		prti->idx_ngh_BH = prtj->index;
	    }
	    if(r2_tmp < prtj->r2_ngh_BH){ 
		prtj->r2_ngh_BH = r2_tmp;
		prtj->idx_ngh_BH = prti->index;
	    }
	}
    }

    for(int i=0; i<NBH_LOC; i++){
	int adr_tmp = prt[i].index;
	prt[i].acc = BH_glb[adr_tmp].acc;
	prt[i].pot = BH_glb[adr_tmp].pot;
	prt[i].Nj = BH_glb[adr_tmp].Nj;
	prt[i].r2_ngh_FS = BH_glb[adr_tmp].r2_ngh_FS;
	prt[i].r2_ngh_BH = BH_glb[adr_tmp].r2_ngh_BH;
	prt[i].idx_ngh_FS = BH_glb[adr_tmp].idx_ngh_FS;
	prt[i].idx_ngh_BH = BH_glb[adr_tmp].idx_ngh_BH;
    }
}



///////////////////////////////////////////
//                                       //
//            SEARCH NEIGHBOUR           //
//                                       //
///////////////////////////////////////////
void Soft_System::search_neighbour_index(Particle BH_glb[],
					 int &ngh_list_len_loc,
					 int &ngh_list_len_glb,
					 Neighbour_List_Index ngh_list_idx[],
					 int &Nshort_loc,
					 int &Nshort_glb, 
					 const int &flag_min_box = 0){
    ngh_list_len_loc = 0;
    double rsearch2_tmp = 0.0;
    for(int i=0; i<NFS_LOC+NBH_LOC; i++){
	Particle *prti = prt+i;
	int njfs_tmp = prti->Nj;
	prti->Nj = 0;

	/////////////////////////
	// search FS(iprt) - BH(jprt) and BH - BH, directly
	if(prti->index < NBH_GLB){
	    rsearch2_tmp = rsearch2_BH_BH;
	}
	else if(NBH_GLB <= prti->index){
	    rsearch2_tmp = rsearch2_FS_BH;
	}
	else{
	    prti->dump();
	    cerr<<"error (search_neighbour)"<<endl;
	    halt_program();
	}
	for(int j=0; j<NBH_GLB; j++){
	    Particle *BHj = BH_glb+j;
	    if(prti->index == BHj->index || rsearch2_tmp < prti->r2_ngh_BH){continue;}
	    double r2_tmp = (prti->pos - BHj->pos) * (prti->pos - BHj->pos);
	    if(r2_tmp > rsearch2_tmp){continue;}
	    ngh_list_idx[ngh_list_len_loc].idx_i = prti->index;
	    ngh_list_idx[ngh_list_len_loc].idx_j = BHj->index;
	    prti->Nj++;
	    ngh_list_len_loc++;
	    prti->have_ngh = 1;
	    BHj->have_ngh = 1;
	}

	/////////////////////////
	// search FS as j-partcile using tree
	if(prti->index < NBH_GLB){
	    rsearch2_tmp = rsearch2_FS_BH;
	}
	else if(NBH_GLB <= prti->index){
	    rsearch2_tmp = rsearch2_FS_FS;
	}
#ifdef SEQUOIA
	if( njfs_tmp == 1 
	    && prti->idx_ngh_FS >= 0 
	    && prti->idx_ngh_FS < NFS_GLB_ORG+NBH_GLB_ORG ){
	    ngh_list_idx[ngh_list_len_loc].idx_i = prti->index;
	    ngh_list_idx[ngh_list_len_loc].idx_j =  prti->idx_ngh_FS;
	    ngh_list_len_loc++;
	    prti->Nj += njfs_tmp;
	}
	else if( njfs_tmp > 1){
	    cell_tree[0].search_neighbour_index(prt[i], rsearch2_tmp, 
						ngh_list_idx,  ngh_list_len_loc, 
						NGH_LIST_MAX, flag_min_box);
	}
#else
	if(rsearch2_tmp < prti->r2_ngh_FS 
	   && prti->idx_ngh_FS >= 0 
	   && prti->idx_ngh_FS < NFS_GLB_ORG+NBH_GLB_ORG ){
	    continue;
	}
	else{
	    cell_tree[0].search_neighbour_index(prt[i], rsearch2_tmp, 
						ngh_list_idx,  ngh_list_len_loc, 
						NGH_LIST_MAX, flag_min_box);
	}
#endif // SEQUOIA
    }

    Nshort_loc = 0;
    for(int i=0; i<NFS_LOC+NBH_LOC; i++){
	Particle *prti = prt+i;
	if(prti->Nj > 0){
	    Nshort_loc++;
	}
    }
    Nshort_glb = mpi_sum_int_np(Nshort_loc);
    ngh_list_len_glb = mpi_sum_int_np(ngh_list_len_loc);
    for(int i=NALL_LOC-1; i >NALL_LOC-1-(NBH_GLB-NBH_LOC); i--){
	int adr = prt[i].index;
	prt[i].have_ngh = BH_glb[adr].have_ngh;
    }
}


/////////////////////////////////////////////////////
//                                                 //    
//     SHORT RANGE FORCE CALCULATION (acc only)    //
//                                                 //    
/////////////////////////////////////////////////////

// probablly update of acc is not need
void Soft_System::evaluate_acc_short_ipara(Neighbour_List_Index ngh_list_idx[],
					   const int &ngh_list_len_loc,
					   mymap &adr_from_idx){
    double eps2_tmp = 0.0;
    double r2_tmp = 0.0;
    double rcut_in_tmp = 0.0;
    double rcut_out_tmp = 0.0;
    for(int i=0; i<ngh_list_len_loc; ){
	int adr_i = adr_from_idx[ngh_list_idx[i].idx_i];
	Particle *prti = prt+adr_i;
	prti->acc_short = 0.0;
	prti->pot_short = 0.0;
	int Nj = prti->Nj;
	for(int j=0; j<Nj; j++){
	    int adr_j = adr_from_idx[ngh_list_idx[i+j].idx_j];
	    Particle *prtj = prt+adr_j;
	    if(prti->index == prtj->index){continue;}
	    if(prti->index < NBH_GLB_ORG && prtj->index < NBH_GLB_ORG){
		rcut_out_tmp = rcut_out_BH_BH; 
		rcut_in_tmp = rcut_in_BH_BH; 
		eps2_tmp = eps2_BH_BH;
	    }
	    else if(prti->index >= NBH_GLB_ORG && prtj->index >= NBH_GLB_ORG){
		rcut_out_tmp = rcut_out_FS_FS; 
		rcut_in_tmp = rcut_in_FS_FS;
		eps2_tmp = eps2_FS_FS;
	    }
	    else if(prti->index < NBH_GLB_ORG || prtj->index < NBH_GLB_ORG){
		rcut_out_tmp = rcut_out_FS_BH;
		rcut_in_tmp = rcut_in_FS_BH;
		eps2_tmp = eps2_FS_BH;
	    }
	    calc_acc_duncan4(prti->pos,  
			     prti->acc_short, prti->pot_short, 
			     prtj->pos,
			     prtj->mass,
			     eps2_tmp,  rcut_in_tmp,  rcut_out_tmp, 
			     r2_tmp);
	    if(prtj->index < NBH_GLB_ORG && r2_tmp < prti->r2_ngh_BH){
		prti->r2_ngh_BH = r2_tmp;
		prti->idx_ngh_BH = prtj->index;
	    }
	    if(prtj->index >= NBH_GLB_ORG && r2_tmp < prti->r2_ngh_FS){
		prti->r2_ngh_FS = r2_tmp;
		prti->idx_ngh_FS = prtj->index;
	    }
	}
	i += Nj;
    }
}

void Soft_System::dump_cell(ofstream &fout_cell){
    if(MYRANK!=0){return;}
    if(MYRANK == 0){
	fout_cell<<"XXXXX"<<endl;
	for(int i=0; i<NCELL_TREE_LOC_MAX-heap_remainder; i++){
	    fout_cell<<cell_tree[i].pos_center<<"   "<<cell_tree[i].size<<endl;
	}
    }
}



#ifdef SEQUOIA

#include <stdlib.h>
#include <vector>
#include <fstream>
#include"sequoiaInterface.h"
#include"my_cuda.h"
#include"node_specs.h"
#include"tipsydefs.h"

// 1) construct tree & evaluate gravity through SEQUOIA
// 2) neighbour search (construct ngh_list_idx[])
// 3) calc acc from close particles
// set ngh_list_len_loc, 
// set ngh_list_len_glb, 
// set Nshort_loc,  
// set Nshort_glb, 
// set ngh_list_idx
// set adr_from_idx
class HOST_RANK_DEVID{
public:
    char hostname[STRINGSIZE];
    int rank;
    int devid;
    HOST_RANK_DEVID(){
	rank = -1;
	devid = 0;
    }
};

void Soft_System::build_tree_on_host_using_sequoia(uint leafNodeIdx[], 
						   uint2 node_bodies[],
						   uint n_children[],
						   float4 boxCenterInfo[],
						   float4 boxSizeInfo[],
						   int n_leafs, 
						   int n_nodes,
						   my_dev::dev_mem<int>  &j_adr_buff){
    this->clear();
    for(int i_node = 0; i_node < n_nodes; i_node++){
	int nodeID = leafNodeIdx[i_node];
	Cell_Tree* current = this->cell_tree+nodeID;
	current->pos_center.set((double)boxCenterInfo[nodeID].x, (double)boxCenterInfo[nodeID].y, (double)boxCenterInfo[nodeID].z);
	current->size.set((double)boxSizeInfo[nodeID].x, (double)boxSizeInfo[nodeID].y, (double)boxSizeInfo[nodeID].z);
	current->size *= 2.0001;

	if(NPROC > 1){
	    // do something
	}
	current->Nchild = 0;
	uint2 bij          =  node_bodies[nodeID];
	uint firstChild    =  bij.x & ILEVELMASK;
	uint lastChild     =  bij.y;
	current->Nprt = lastChild - firstChild;
	int level = (bij.x&LEVELMASK)>>BITLEVELS;
	current->level  = level;
	if(i_node < n_leafs){
	    current->isleaf = 1;
	    Particle* prt_tmp = prt + j_adr_buff[firstChild];
	    current->prt_first = prt_tmp;
	    for(int idx = firstChild+1; idx < lastChild; idx++){
		int adr = j_adr_buff[idx];
		prt_tmp->prt_next = prt + adr;
		prt_tmp = prt_tmp->prt_next;
	    }
	    prt_tmp->prt_next = NULL;
	}
	else{
	    current->isleaf = 0;
	}
	this->heap_remainder--;
	this->heap_top++;
    }
    for(int i_node = n_leafs; i_node < n_nodes; i_node++){
	int nodeID = leafNodeIdx[i_node];
	uint firstChild = n_children[nodeID]&BODYMASK;
	uint nChildren = (n_children[nodeID]&INVBMASK)>>LEAFBIT;
	Cell_Tree* current = this->cell_tree + nodeID;
	for(int i=0; i<8; i++){
	    current->child[i] = NULL;
	}
	for(int i_child = firstChild; i_child < firstChild + nChildren; i_child++){
	    Cell_Tree* child = this->cell_tree + i_child;
	    int idxchild = current->Nchild;
	    if(current->child[idxchild] != NULL){
		fout_debug<<"someone open the child"<<endl;
		halt_program();
	    }
	    current->child[idxchild] = child;
	    current->Nchild++;
	}
    }
}

void Soft_System::build_tree_and_get_force_using_sequoia(char *argv[],
							 Vector3 xlow[], 
							 Vector3 xhigh[],
							 int dst_rank[],
							 int src_rank[],
							 Particle BH_glb[]){
    
    double tcal_offset, tcal_offset_long;
    float eps2 = (float)(eps2_FS_FS);
    float eps = sqrt(eps2);
    float theta = sqrt(theta2);

    /////////////////////////////////////////
    ////////// initialize sequoia ///////////

    static int first = 1;

    static my_dev::context devContext;
    //Input
    static my_dev::dev_mem<real4> j_pos_buff;  //Bodies positions
    static my_dev::dev_mem<int>   j_adr_buff;  //Bodies ids
    // if NPROC > 1
    static my_dev::dev_mem<real4> i_pos_buff;  //Bodies positions
    static my_dev::dev_mem<int>   i_adr_buff;  //Bodies ids

    //Output
    static my_dev::dev_mem<real4> i_acc_buff; //Bodies Accelerations
    static my_dev::dev_mem<real>  i_ngh_ds2_buff; //Bodies distance to nearest neighbour squared
    static my_dev::dev_mem<int>   i_ngh_adr_buff;  //Bodies nearest neighbour
    static my_dev::dev_mem<int>   i_Nngh_buff;  //Bodies nearest neighbour

    if(first){
	char hostname[STRINGSIZE];
	gethostname(hostname, STRINGSIZE);
	HOST_RANK_DEVID host_rank_devid_send;
	strcpy(host_rank_devid_send.hostname, hostname);
	host_rank_devid_send.rank = MYRANK;
	host_rank_devid_send.devid = 0;
	HOST_RANK_DEVID* host_rank_devid_recv = new HOST_RANK_DEVID[NPROC];
	mpi_allgatherv_T(&host_rank_devid_send, 1, host_rank_devid_recv, NPROC, NPROC);
	int devID = 0;
	for(int i=0; i<MYRANK; i++){
	    if( strcmp(host_rank_devid_recv[i].hostname, hostname) == 0){
		devID++;
	    }
	}
	delete [] host_rank_devid_recv;
	cerr<<"hostname="<<hostname<<endl;
	cerr<<"MYRANK="<<MYRANK<<endl;
	cerr<<"devID="<<devID<<endl;

	devContext = sequoia_init(argv, devID, theta, eps, rsearch2_FS_FS);

        j_pos_buff.setContext(devContext);
        j_pos_buff.cmalloc(NPRT_LOC_MAX);
        j_adr_buff.setContext(devContext);
        j_adr_buff.cmalloc(NPRT_LOC_MAX);
        i_acc_buff.setContext(devContext);
        i_acc_buff.cmalloc(NPRT_LOC_MAX);
        i_ngh_ds2_buff.setContext(devContext);
        i_ngh_ds2_buff.cmalloc(NPRT_LOC_MAX);
        i_ngh_adr_buff.setContext(devContext);
        i_ngh_adr_buff.cmalloc(NPRT_LOC_MAX);
        i_Nngh_buff.setContext(devContext);
        i_Nngh_buff.cmalloc(NPRT_LOC_MAX);

	if(NPROC > 1){
	    // do something
	}

        first = 0;

    }

    ///////////////////////////////////////////////////
    ////////// send particle data to sequoia //////////
    tcal_offset_long = mpi_wtime();
    tcal_offset = mpi_wtime();
    int NFS_add = 0;
    for(int j=0; j<NALL_LOC; j++){
        if(prt[j].index < NBH_GLB_ORG){continue;}
        j_pos_buff[NFS_add].w = prt[j].mass;
        j_pos_buff[NFS_add].x = prt[j].pos[0];
        j_pos_buff[NFS_add].y = prt[j].pos[1];
        j_pos_buff[NFS_add].z = prt[j].pos[2];
  	j_adr_buff[NFS_add] = j;
        NFS_add++;
    }
    j_pos_buff.h2d(NFS_add);
    j_adr_buff.h2d(NFS_add);

    ///////////////////////////////////////////////////////////
    ////////// build tree on GPU and get tree info ////////////
    ////////// first tree construction  

    // IN
    int Nj = NFS_LOC;
    bool j_sort = true;

    // OUT
    uint* leafNodeIdx;
    uint2* node_bodies;
    uint* n_children;
    float4* boxSizeInfo;
    float4* boxCenterInfo;
    int n_leafs;
    int n_nodes;

    // OUT (for MPI)
    uint2* level_list;
    uint* node_level_list;
    real4* multipole;
    int n_levels;

    cerr<<"Nj="<<Nj<<endl;
    sequoia_setParticlesAndGetGravity_firsthalf_for_neighbour_search(j_pos_buff, j_adr_buff, Nj, j_sort,
								     leafNodeIdx, node_bodies, n_children,  
								     boxSizeInfo, boxCenterInfo,
								     n_leafs, n_nodes);
    TCAL_FIRST_HALF = mpi_wtime() - tcal_offset;


    //////////////////////////////////////////
    ////////// build tree on HOST ////////////
    tcal_offset = mpi_wtime();
    build_tree_on_host_using_sequoia(leafNodeIdx, node_bodies, n_children, 
				     boxCenterInfo, boxSizeInfo, 
				     n_leafs, n_nodes, j_adr_buff);
    TCAL_TREE_SETUP_LOC = mpi_wtime() - tcal_offset;

    tcal_offset = mpi_wtime();
    if(NPROC > 1){
	// do something
    }
    TCAL_TREE_SET_CM_LOC = mpi_wtime() - tcal_offset;

    tcal_offset = mpi_wtime();
    exchange_LET(xlow, xhigh, dst_rank, src_rank, BH_glb); // set prt & prt_tree (include BH)
    TCAL_TREE_LET_EX = mpi_wtime() - tcal_offset;


    tcal_offset = mpi_wtime();
    if(NPROC > 1){
	// do something
    }
    TCAL_TREE_SETUP = mpi_wtime() - tcal_offset_long;

    //////////////////////////////////////////////
    ////////// evaluate gravity on GPU ///////////
    tcal_offset = mpi_wtime();

    for(int i=0; i<NFS_LOC+NBH_LOC; i++){
        prt[i].acc = 0.0;
        prt[i].pot = 0.0;
    }
    for(int i=0; i<NBH_GLB_ORG; i++){
        BH_glb[i].acc = 0.0;
        BH_glb[i].pot = 0.0;
    }
    for(int i=0; i<NALL_LOC; i++){
	prt[i].have_ngh = 0;
	prt[i].Nj = 0;
    }

    int Ni = NFS_add;
    if(NPROC > 1){
	// do something
    }
    else{
	Ni = Nj;
	bool i_sort = false;
	sequoia_setParticlesAndGetGravity_lasthalf(j_pos_buff, j_pos_buff, 
						   j_adr_buff, Ni,
						   i_sort, i_acc_buff, i_ngh_ds2_buff,
						   i_ngh_adr_buff, i_Nngh_buff);
    }
    i_acc_buff.d2h(Ni);
    i_ngh_ds2_buff.d2h(Ni);
    i_ngh_adr_buff.d2h(Ni);
    i_Nngh_buff.d2h(Ni);
    if(NPROC > 1){
	// do something
    }
    else{
	j_adr_buff.d2h(Ni);
	for(int i=0; i<Ni; i++){
	    int adr = j_adr_buff[i];
	    prt[adr].acc[0] = i_acc_buff[i].x;
	    prt[adr].acc[1] = i_acc_buff[i].y;
	    prt[adr].acc[2] = i_acc_buff[i].z;
	    prt[adr].pot = i_acc_buff[i].w;
	    prt[adr].r2_ngh_FS = i_ngh_ds2_buff[i] - eps2;
	    if(i_ngh_adr_buff[i] < 0){
		prt[adr].idx_ngh_FS = -1;
	    }
	    else{
		prt[adr].idx_ngh_FS = prt[j_adr_buff[i_ngh_adr_buff[i]]].index;
	    }
	    prt[adr].Nj = i_Nngh_buff[i];
	    if(prt[adr].Nj > 0){
		prt[adr].have_ngh = 1;
		prt[j_adr_buff[i_ngh_adr_buff[i]]].have_ngh = 1;
	    }
	}
    }
    TCAL_TREE_EVALUATE = mpi_wtime() - tcal_offset;
}

#endif
