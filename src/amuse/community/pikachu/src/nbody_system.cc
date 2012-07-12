#include"system.h"
#include"const.h"
#include"misc.h"

Nbody_System::Nbody_System(){
    cerr<<"Nbosy_System constructor"<<endl;

    prt_loc = new Particle[NPRT_LOC_MAX];
    BH_glb = new Particle[NBH_GLB_MAX];
    prt_dead_glb  = new Particle[NPRT_LOC_MAX / 1000];

    Tsys = Tend = pos_max = mass_min = 0.0;
    ngh_list_len_loc = ngh_list_len_glb = Nshort_loc = Nshort_glb = 0;

    ngh_list_idx = new Neighbour_List_Index[NGH_LIST_MAX];
    ngh_list = new Neighbour_List[NGH_LIST_MAX];

    soft_system.init(prt_loc);

    int I, ID, IS;
    for(I=0, ID=NPROC-1, IS=1; I<NPROC-1; I++, ID--, IS++){
	dst_rank[I] = (ID+MYRANK)%NPROC;
	src_rank[I] = (IS+MYRANK)%NPROC;
	if(MYRANK == 0){
	    cerr<<"dst_rank[I]="<<dst_rank[I]<<endl;
	    cerr<<"src_rank[I]="<<src_rank[I]<<endl;
	}
    }
    npdim[0] = npdim[1] = npdim[2] = sample_freq = 0;
}

void Nbody_System::calc_energy(double& Ek,
			       double& Ep,
			       double& E){ 
    Ek = Ep = E = 0.0;
    double Ek_loc = 0.0;
    double Ep_loc = 0.0;
    double E_loc = 0.0;
    for(int i=0; i<NFS_LOC+NBH_LOC; i++){
	Ek_loc += prt_loc[i].mass*(prt_loc[i].vel*prt_loc[i].vel);
	Ep_loc += prt_loc[i].mass*prt_loc[i].pot;
    }
    Ek_loc *= 0.5;
    Ep_loc *= 0.5;
    E_loc = Ek_loc + Ep_loc;
    Ek = mpi_sum_double_np(Ek_loc);
    Ep = mpi_sum_double_np(Ep_loc);
    E = mpi_sum_double_np(E_loc);
}
void Nbody_System::calc_energy(const int &mode){
    if(mode == 0)
	calc_energy(Ek0, Ep0, E0);
    else if(mode == 1)
	calc_energy(Ek1, Ep1, E1);
}

void Nbody_System::kick_half(){
    soft_system.kick(dt_glb*0.5);
}

void Nbody_System::kick_full(){
    soft_system.kick(dt_glb);
}

void Nbody_System::merge_ngh_list(){
    hard_system.merge_ngh_list(prt_loc,  ngh_list_idx, adr_from_idx,
			       Nshort_loc,  Nshort_glb,
			       ngh_list_len_loc,  ngh_list_len_glb,
			       ngh_list);
}



//////////////////////////////////
//                              //
//      HIGH LEVEL FUNCTION     //
//                              //
//////////////////////////////////
// call once at the begening
// set sample_freq
// set npdim[3]
void Nbody_System::initialize_division(){
    sample_freq = determine_sample_freq(NFS_GLB+NBH_GLB);
    create_division(NPROC, npdim[0], npdim[1], npdim[2]);
}

// assigne particles to each node
// set NFS_LOC
// set NBH_LOC
// set xlow[], xhigh[]
// set pos_max
void Nbody_System::divide_particles(){
    distribute_particle(prt_loc,  BH_glb,
                        NBH_GLB_ORG, 
                        NFS_LOC,  NBH_LOC,  NPRT_LOC_MAX,
                        npdim,  sample_freq,  xlow,  xhigh,
                        dst_rank, src_rank,
                        pos_max);
    NALL_LOC = NFS_LOC+NBH_LOC;
}


// 1) construct tree
// 2) evaluate gravity
// 3) neighbour search (construct ngh_list_idx[])
// 4) calc acc from close particles
// set ngh_list_len_loc, 
// set ngh_list_len_glb, 
// set Nshort_loc,  
// set Nshort_glb, 
// set ngh_list_idx
// set adr_from_idx
void Nbody_System::calc_soft_forces(char *argv[]){
    double tcal_offset, tcal_offset_long;
    int flag_min_box = 1;
#ifdef SEQUOIA
    flag_min_box = 1;
    // 1 FS get forces from all FSs (exclude the contribution for BHs)
    // 2 FS get r2_ngh_FS and idx_ngh_FS
    // 3 FS get Nj (NOTE BHs is not included)
    // Nj means the number of FSs before making neighbour list,
    // because in function "search_neighbour_index", if prt.Nj == 1, neighbour search can be skipped.
    soft_system.build_tree_and_get_force_using_sequoia(argv, xlow, xhigh,
						       dst_rank, src_rank, BH_glb);
#else
    flag_min_box = 0;
    tcal_offset_long = mpi_wtime();
    soft_system.clear();

    tcal_offset = mpi_wtime();
    soft_system.setup_tree_loc(pos_max);
    TCAL_TREE_SETUP_LOC = mpi_wtime() - tcal_offset;

    tcal_offset = mpi_wtime();
    soft_system.exchange_LET(xlow, xhigh, dst_rank, src_rank, BH_glb); // set prt_loc & prt_tree (include BH)
    TCAL_TREE_LET_EX = mpi_wtime() - tcal_offset;

    tcal_offset = mpi_wtime();
    soft_system.setup_tree_adding_LET();
    TCAL_TREE_SETUP_LET = mpi_wtime() - tcal_offset;

    TCAL_TREE_SETUP = mpi_wtime() - tcal_offset_long;

    ////////////////////////////////
    ////////// initialize //////////
    for(int i=0; i<NFS_LOC+NBH_LOC; i++){
        prt_loc[i].acc = 0.0;
        prt_loc[i].pot = 0.0;
    }
    for(int i=0; i<NBH_GLB_ORG; i++){
        BH_glb[i].acc = 0.0;
        BH_glb[i].pot = 0.0;
    }
    for(int i=0; i<NALL_LOC; i++){
	prt_loc[i].have_ngh = 0;
	prt_loc[i].Nj = 0;
    }


    tcal_offset = mpi_wtime();
    soft_system.evaluate_gravity(pos_max);
    TCAL_TREE_EVALUATE = mpi_wtime() - tcal_offset;
#endif // SEQUOIA
    // now, prt.Nj is the number of FSs (w/o BHs).


    ////////////////////////////////////////////////////////////////
    ////////// direct evaluate gravity from all particles //////////
    // 1 FS get forces from all BHs 
    // 2 BH get forces from all particles (all FSs and all BHs)
    // 3 if iprt is FS, nearest BH index and its distance.
    //   if iprt is BH, nearest BH and FS index and its distance.
    // NOTE FS.Nj and BH.Nj are the number of neighbour FSs only (not BHs)
    tcal_offset = mpi_wtime();
    soft_system.evaluate_gravity_of_FSloc_BHglb_jpara(BH_glb);
    TCAL_DIRECT_EVALUATE = mpi_wtime() - tcal_offset;

    //////////////////////////////////////
    ////////// search neighbour //////////
    tcal_offset_long = mpi_wtime();
    soft_system.search_neighbour_index(BH_glb,  ngh_list_len_loc,  ngh_list_len_glb,
				       ngh_list_idx, Nshort_loc, Nshort_glb,
				       flag_min_box);

    ////////////////////////////////////
    ////////// make dictionary /////////
    tcal_offset = mpi_wtime();
    adr_from_idx.clear();
    for(int i=0; i<NALL_LOC; i++){
        Particle *p = prt_loc+i;
	if(!(p->have_ngh)){continue;}
        adr_from_idx.insert(mymap::value_type(p->index, i));
    }

    TCAL_NGH_SEARCH_MAKE_DIC = mpi_wtime() - tcal_offset;
    TCAL_NGH_SEARCH = mpi_wtime() - tcal_offset_long;

    ////////////////////////////////////
    ////////// calc short acc //////////
    tcal_offset = mpi_wtime();
    for(int i=0; i<NFS_LOC+NBH_LOC; i++){
        prt_loc[i].acc_short = 0.0;
        prt_loc[i].pot_short = 0.0;
    }
    soft_system.evaluate_acc_short_ipara(ngh_list_idx,  ngh_list_len_loc, adr_from_idx);
    for(int i=0; i<NFS_LOC+NBH_LOC; i++){
        // Both prti->acc and prti->acc_short include PN term.
        // Thus prti->acc_long has Newton term only
        prt_loc[i].acc_long = prt_loc[i].acc - prt_loc[i].acc_short;
    }
    TCAL_ACC_ONLY= mpi_wtime() - tcal_offset;
    for(int i=0; i<NALL_LOC; i++){
        prt_loc[i].prt_next = NULL;
    }
}

////////////////////////////////////////////
/// evolve through the inner hamiltonian ///
////////////////////////////////////////////
void Nbody_System::evolve_hard_system(char *argv[]){

    ///////////////////////////
    ////////// drift //////////
    double tcal_offset = mpi_wtime();
    soft_system.drift(dt_glb);
    TCAL_DRIFT = mpi_wtime() - tcal_offset;

    //////////////////////////////////
    ////////// hermite integ /////////
    tcal_offset = mpi_wtime();
    TCAL_HERMITE4_COMM = 0.0;
    TCAL_HERMITE4_PRE = 0.0;
    hard_system.hermite4_para(Tsys, dt_glb, Nshort_glb);
    NGH_LIST_LEN_MULBODY = ngh_list_len_glb - NGH_LIST_LEN_2BODY;
    TCAL_HERMITE4= mpi_wtime() - tcal_offset;
    tcal_offset = mpi_wtime();
    hard_system.copy_prt_short_to_prt_loc(prt_loc,  BH_glb);
    TCAL_COPY_PRT_SHORT = mpi_wtime() - tcal_offset;
}


int Nbody_System::evolve(char output_dir[],
			 char *argv[],
			 int &after_stellar_evolution){

    double tcal_offset_loop = mpi_wtime();
    double tcal_offset = mpi_wtime();

    if(after_stellar_evolution){
	// new forces are required
	calc_soft_forces(argv);
	after_stellar_evolution = 0;
    }

    cerr<<"kick 1"<<endl;
    //////////////////////////////////////////////////////////
    ///////// evolve under OUTER Hamiltonian (KICK) //////////
    tcal_offset = mpi_wtime();
    kick_half();
    TCAL_KICK = mpi_wtime() - tcal_offset;


    cerr<<"copy"<<endl;
    ////////////////////////////////////////////////////////////////////
    ///////// copy and transfered the imformation of prt_short /////////
    tcal_offset = mpi_wtime();
    merge_ngh_list();
    TCAL_MERGE_LIST = mpi_wtime() - tcal_offset;

    //////////////////////////////////
    ///// dump calculation time /////
    if(MYRANK == 0)
	dump_calc_time(fout_tcal, Tsys, Nshort_glb, ngh_list_len_glb, dt_glb);

    static int first = 1;
    if(first){
        E1 = E0;
        Ek1 = Ek0;
        Ep1 = Ep0;
	if(MYRANK == 0){
	    char sout[STRINGSIZE];
	    sprintf(sout,"%s/energy.dat", output_dir);
	    fout_energy.open(sout);
	    sprintf(sout,"%s/BH.dat", output_dir);
	    fout_BH.open(sout);
	    sprintf(sout,"%s/tcal.dat", output_dir);
	    fout_tcal.open(sout);

	    fout_energy<<setprecision(15);
	    fout_BH<<setprecision(15);
	    fout_tcal<<setprecision(15);

	    fout_energy<<Tsys<<"   "<<(E1-E0)/E0<<"  "<<E1<<"   "<<Ek1<<"   "<<Ep1<<"   "<<E0<<"   "<<Ek0<<"   "<<Ep0<<endl;
	    cerr<<Tsys<<"   "<<(E1-E0)/E0<<"   "<<E1<<"   "<<Ek1<<"   "<<Ep1<<endl;
	}
        first = 0;
    }



    cerr<<"DRIFT"<<endl;
    ///////////////////////////////////////////////////////////
    ///////// evolve under INNER Hamiltonian (DRIFT) //////////
    evolve_hard_system(argv);
    Tsys += dt_glb;
    TCAL_EVOLVE_HARD = mpi_wtime() - tcal_offset;

    cerr<<"divide particles"<<endl;
    //////////////////////////////////////
    ///////// divide particles ///////////
    tcal_offset = mpi_wtime();
    if(NPROC > 1){
	divide_particles();
    }
    TCAL_DIVIDE_PRT = mpi_wtime() - tcal_offset;

    cerr<<"evaluate soft forces"<<endl;
    /////////////////////////////////////////
    ///////// evaluate soft forces //////////
    tcal_offset = mpi_wtime();
    calc_soft_forces(argv);
    TCAL_SOFT_FORCES = mpi_wtime() - tcal_offset;

    cerr<<"kick2"<<endl;
    //////////////////////////////////////////////////////////
    ///////// evolve under OUTER Hamiltonian (KICK) //////////
    tcal_offset = mpi_wtime();
    kick_half();
    TCAL_KICK += mpi_wtime() - tcal_offset;

    calc_energy(1);
    if(fmod(Tsys, snp_interval) == 0.0){
	write_file(output_dir);
    }
    if(MYRANK == 0){ 
	fout_energy<<Tsys<<"   "<<(E1-E0)/E0<<"  "<<E1<<"   "<<Ek1<<"   "<<Ep1<<"   "<<E0<<"   "<<Ek0<<"   "<<Ep0<<endl;
	if(NBH_GLB_ORG > 1){
	    fout_BH<<Tsys;
	    for(int i=0; i<NBH_GLB_ORG; i++){
		fout_BH<<"   "<<BH_glb[i].mass<<"   "<<BH_glb[i].pos<<"   "<<BH_glb[i].vel<<"   "<<BH_glb[i].pot;
	    }
	    fout_BH<<endl;
	}
    }


    TCAL_LOOP = mpi_wtime() - tcal_offset_loop;

    return 0;
}







//////////////////////////
//                      //
//      IO FUNCTION     //
//                      //
//////////////////////////
void Nbody_System::write_file(char output_dir[]){
    write0_gather(prt_loc, prt_dead_glb, Tsys, Egr0, output_dir, snp_id);
    snp_id++;
}

void Nbody_System::read_file(char param_file[],
			     char output_dir[]){
    char input_file[STRINGSIZE];
    double eps2_FS_FS_tmp, eps2_FS_BH_tmp, eps2_BH_BH_tmp;
    double rcut_out_FS_FS_tmp, rcut_out_FS_BH_tmp, rcut_out_BH_BH_tmp;
    double rsearch_FS_FS_tmp, rsearch_FS_BH_tmp, rsearch_BH_BH_tmp;
    double Tend_tmp, dt_glb_tmp;
    int snp_id_tmp, read_flag;
    double snp_interval_tmp;
    double theta2_tmp;
    int Nleaf_tmp, Ncrit_tmp, quad_flag_tmp;
    double eta_s_tmp, eta_FS_tmp, eta_BH_tmp;
    double vel_light_tmp;

    readparam(param_file,  input_file,  output_dir,
              BH_glb,   NBH_GLB,
              eps2_FS_FS_tmp,  eps2_FS_BH_tmp,  eps2_BH_BH_tmp,
              rcut_out_FS_FS_tmp,  rcut_out_FS_BH_tmp,  rcut_out_BH_BH_tmp,
	      rsearch_FS_FS_tmp,  rsearch_FS_BH_tmp,  rsearch_BH_BH_tmp,
              Tend_tmp,  dt_glb_tmp,
              snp_id_tmp,  snp_interval_tmp,  read_flag,
              theta2_tmp,  Nleaf_tmp,  Ncrit_tmp, quad_flag_tmp,
              eta_s_tmp,  eta_FS_tmp,  eta_BH_tmp,
              vel_light_tmp);

    Tend = Tend_tmp;
    snp_id = snp_id_tmp;
    snp_interval = snp_interval_tmp;
    Egr0 = 0.0;

    double Tsys_tmp;
    if(MYRANK == 0){
        for(int i=0; i<NBH_GLB; i++){
            prt_loc[i] = BH_glb[i];
            prt_loc[i].type = blackhole;
        }
        NBH_LOC = NBH_GLB;
    }
    else{
        NBH_LOC = 0;
    }

    if(read_flag == 0){
        read0_scatter(prt_loc,  NFS_GLB,  NFS_LOC,
                      Tsys_tmp,  input_file, 
                      NBH_GLB,  NBH_LOC);
    }
    else if(read_flag == 1){
        read1_scatter(prt_loc, BH_glb, prt_dead_glb,
                      Tsys_tmp, input_file);
    }
    cerr<<"finish reading snap"<<endl;

    dump_cerr(prt_loc[0].pos);
    dump_cerr(soft_system.prt[0].pos);
    dump_cerr(prt_loc[1].pos);
    dump_cerr(soft_system.prt[1].pos);

    int NDEAD_LOC = 0;
    for(int i=0; i<NFS_LOC+NBH_LOC; i++){
        if(prt_loc[i].mass == 0.0){
            NDEAD_LOC++;
        }
    }

    Tsys = Tsys_tmp;
    dt_glb = dt_glb_tmp;
    NDEAD_GLB = mpi_sum_int_np(NDEAD_LOC);
    NFS_GLB_ORG = NFS_GLB + NDEAD_GLB;
    NBH_GLB_ORG = NBH_GLB + NDEAD_GLB;

    dump_cout(eps2_FS_FS_tmp);
    dump_cout(eps2_FS_BH_tmp);
    dump_cout(eps2_BH_BH_tmp);
    dump_cout(rcut_out_FS_FS_tmp);
    dump_cout(rcut_out_FS_BH_tmp);
    dump_cout(rcut_out_BH_BH_tmp);
    dump_cout(rsearch_FS_FS_tmp);
    dump_cout(rsearch_FS_BH_tmp);
    dump_cout(rsearch_BH_BH_tmp);
    dump_cout(NBH_GLB_ORG);
    dump_cout(NBH_GLB);
    dump_cout(NBH_LOC);
    dump_cout(NFS_GLB_ORG);
    dump_cout(NFS_GLB);
    dump_cout(NFS_LOC);
    dump_cout(NDEAD_GLB);

    double mass_min_loc = LARGEDOUBLE;
    for(int i=0; i<NFS_LOC+NBH_LOC; i++){
        if(prt_loc[i].mass < mass_min_loc){
            mass_min_loc = prt_loc[i].mass;
        }
    }
    mass_min = mpi_min_double(&mass_min_loc);

    soft_system.set(eps2_FS_FS_tmp, eps2_FS_BH_tmp, eps2_BH_BH_tmp,
                    rcut_out_FS_FS_tmp, rcut_out_FS_BH_tmp, rcut_out_BH_BH_tmp,
                    rsearch_FS_FS_tmp, rsearch_FS_BH_tmp, rsearch_BH_BH_tmp,
                    theta2_tmp, Ncrit_tmp, Nleaf_tmp, quad_flag_tmp);

    hard_system.set(eps2_FS_FS_tmp, eps2_FS_BH_tmp, eps2_BH_BH_tmp,
                    rcut_out_FS_FS_tmp, rcut_out_FS_BH_tmp, rcut_out_BH_BH_tmp,
                    eta_s_tmp, eta_FS_tmp, eta_BH_tmp, mass_min);


    soft_system.dump();
    hard_system.dump();
}
