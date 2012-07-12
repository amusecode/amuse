#include"system.h"
#include"const.h"


Hard_System::Hard_System(){
    cerr<<"Hard_System constructor"<<endl;
    prt_short = new Particle_Short[NSHORT_MAX];
    Tnext = new double[NSHORT_MAX];
    adr_time = new int[NSHORT_MAX];
    prt_comm_send = new Particle_Comm[NSHORT_MAX];
    prt_comm_recv = new Particle_Comm[NSHORT_MAX];
    cerr<<"Hard_System constructed"<<endl;
}

void Hard_System::init(){
    generate_new_communicator();
}

void Hard_System::set(const double &_eps2_FS_FS,  
		      const double &_eps2_FS_BH,  
		      const double &_eps2_BH_BH,  
		      const double &_rcut_out_FS_FS,  
		      const double &_rcut_out_FS_BH,  
		      const double &_rcut_out_BH_BH, 
		      const double &_eta_s,
		      const double &_eta_FS,
		      const double &_eta_BH,
		      const double &_mass_min){
    eps2_FS_FS = _eps2_FS_FS; 
    eps2_FS_BH = _eps2_FS_BH; 
    eps2_BH_BH = _eps2_BH_BH; 
    rcut_out_FS_FS = _rcut_out_FS_FS;
    rcut_out_FS_BH = _rcut_out_FS_BH;
    rcut_out_BH_BH = _rcut_out_BH_BH;
    rcut_in_FS_FS = rcut_out_FS_FS * RCUT_IN_FACTOR;
    rcut_in_FS_BH = rcut_out_FS_BH * RCUT_IN_FACTOR;
    rcut_in_BH_BH = rcut_out_BH_BH * RCUT_IN_FACTOR;
    eta_s = _eta_s;
    eta_FS = _eta_FS;
    eta_BH = _eta_BH;

    acc_offset_sq = _mass_min*_mass_min/(rcut_out_FS_FS*rcut_out_FS_FS*rcut_out_FS_FS*rcut_out_FS_FS);

}


/////////////////////////////
//                         //  
//           DUMP          //
//                         //
/////////////////////////////

void Hard_System::dump(ostream &fout){

    fout<<"dump hard system"<<endl;
    fout<<"rcut_out_FS_FS="<<rcut_out_FS_FS<<endl;
    fout<<"rcut_out_FS_BH="<<rcut_out_FS_BH<<endl;
    fout<<"rcut_out_BH_BH="<<rcut_out_BH_BH<<endl;
    fout<<"rcut_in_FS_FS="<<rcut_in_FS_FS<<endl;
    fout<<"rcut_in_FS_BH="<<rcut_in_FS_BH<<endl;
    fout<<"rcut_in_BH_BH="<<rcut_in_BH_BH<<endl;
    fout<<"eps2_FS_FS="<<eps2_FS_FS<<endl;
    fout<<"eps2_FS_BH="<<eps2_FS_BH<<endl;
    fout<<"eps2_BH_BH="<<eps2_BH_BH<<endl;
    fout<<"eta_s="<<eta_s<<endl;
    fout<<"eta_FS="<<eta_FS<<endl;
    fout<<"eta_BH="<<eta_BH<<endl;
}

///////////////////////////////////////////////////////
//                                                   //    
//            SHORT RANGE FORCE CALCULATION          //
//                                                   //    
///////////////////////////////////////////////////////

void Hard_System::generate_new_communicator(){
    for(int np=1; np<NPROC_MAX+1; np++){
	MPI_Comm_split(MPI_COMM_WORLD, MYRANK/np, MYRANK, comm_new_array+np);
	MPI_Comm_rank(comm_new_array[np], myrank_new_array+np);
    }
}

// NOTE: Order of calculations of predictors.
// Initailly, all Ni particles are predicted .
// During force loop on an i particle, 
// if a j particle has been alreaddy predicted, 
// this j particle must not predicted, because
// if this j particle is already integrated, 
// acc and jrk of this particle have been new values in the future.
// Thus each particle must be predicted once.
int Hard_System::evaluate_acc_jrk_short_jpara(const int &Ni,
					      const double &Tnext,
					      const int &mode){
    static Vector3 acc_short_send;
    static Vector3 acc_short_recv;
    static Vector3 jrk_short_send;
    static Vector3 jrk_short_recv;
    static int idx_ngh_short_FS;
    static int idx_ngh_short_BH;
    static dr2_rank dr2_rank_short_FS_send;
    static dr2_rank dr2_rank_short_FS_recv;
    static dr2_rank dr2_rank_short_BH_send;
    static dr2_rank dr2_rank_short_BH_recv;
    double eps2_tmp = 0.0;
    double r2_tmp = 0.0;
    double rcut_in_tmp = 0.0;
    double rcut_out_tmp = 0.0;


    if(mode == 1){
	for(int i=0; i<Ni; i++){
	    Particle_Short *prti = prt_short+adr_time[i];
	    double dt = Tnext - prti->time;
	    prti->predict_short_h4(dt);
	    prti->time_pre = Tnext;
	}
    }

    for(int i=0; i<Ni; i++){
	Particle_Short *prti = prt_short+adr_time[i];
	int Nj = prti->Nj;
	if(mode == 1){
	    if(prti->time_pre != Tnext){
		double dt = Tnext - prti->time;
		prti->predict_short_h4(dt);
		prti->time_pre = Tnext;
	    }
	}
	prti->acc_short_old = prti->acc_short;
	prti->jrk_short_old = prti->jrk_short;
	prti->acc_short = prti->jrk_short = 0.0;
	prti->pot_short = 0.0;
	prti->r2_ngh_FS = prti->r2_ngh_BH = LARGEDOUBLE;


	if(Nj < NJCRIT){
	    for(int j=0; j<Nj; j++){
		Particle_Short *prtj = ((prti->ngh_list_first)+j)->prtj;
		if(prti->index == prtj->index){continue;}
		if(prti->index < NBH_GLB && prtj->index < NBH_GLB){
		    rcut_out_tmp = rcut_out_BH_BH; 
		    rcut_in_tmp = rcut_in_BH_BH; 
		    eps2_tmp = eps2_BH_BH;
		}
		else if(prti->index >= NBH_GLB && prtj->index >= NBH_GLB){
		    rcut_out_tmp = rcut_out_FS_FS; 
		    rcut_in_tmp = rcut_in_FS_FS; 
		    eps2_tmp = eps2_FS_FS;
		}
		else if(prti->index < NBH_GLB || prtj->index < NBH_GLB){
		    rcut_out_tmp = rcut_out_FS_BH;
		    rcut_in_tmp = rcut_in_FS_BH; 
		    eps2_tmp = eps2_FS_BH;
		}
		if(mode == 0){
		    calc_acc_jrk_duncan4(prti->pos,  prti->vel, 
					 prti->acc_short,  prti->jrk_short,
					 prti->pot_short, 
					 prtj->pos,  prtj->vel, 
					 prtj->mass,
					 eps2_tmp,  rcut_in_tmp,  rcut_out_tmp, 
					 r2_tmp);
		}
		else if(mode == 1){
		    if(prtj->time_pre != Tnext){
			double dt = Tnext - prtj->time;
			prtj->predict_short_h4(dt);
			prtj->time_pre = Tnext;
		    }
		    calc_acc_jrk_duncan4(prti->pos_pre,  prti->vel_pre, 
					 prti->acc_short,  prti->jrk_short,
					 prti->pot_short, 
					 prtj->pos_pre,  prtj->vel_pre, 
					 prtj->mass,
					 eps2_tmp,  rcut_in_tmp,  rcut_out_tmp, 
					 r2_tmp);
		}
		if(prtj->index < NBH_GLB && r2_tmp < prti->r2_ngh_BH){
		    prti->r2_ngh_BH = r2_tmp;
		    prti->idx_ngh_BH = prtj->index;
		}
		if(prtj->index >= NBH_GLB && r2_tmp < prti->r2_ngh_FS){
		    prti->r2_ngh_FS = r2_tmp;
		    prti->idx_ngh_FS = prtj->index;
		}
	    }
	}
	else{
	    acc_short_send = 0.0;
	    jrk_short_send = 0.0;

	    dr2_rank_short_FS_send.dr2 = dr2_rank_short_BH_send.dr2 = LARGEDOUBLE;
	    dr2_rank_short_FS_send.rank = dr2_rank_short_BH_send.rank = MYRANK;

	    int Nproc_new = (Nj/NJCRIT+1) >= NPROC ? NPROC : (Nj/NJCRIT+1);
	    int Nj_loc = Nj/Nproc_new;
	    int j_head_loc = Nj_loc*MYRANK;
	    if(MYRANK < Nproc_new){
		if(MYRANK < Nj%Nproc_new){
		    Nj_loc++;
		    j_head_loc += MYRANK;
		}
		else{
		    j_head_loc += Nj%Nproc_new;
		}
	    }
	    else{
		Nj_loc = 0;
		j_head_loc = 0;
	    }

	    for(int j=j_head_loc; j<j_head_loc+Nj_loc; j++){
		Particle_Short *prtj = ((prti->ngh_list_first)+j)->prtj;
		if(prti->index == prtj->index){continue;}
		if(prti->index < NBH_GLB && prtj->index < NBH_GLB){
		    rcut_out_tmp = rcut_out_BH_BH; 
		    rcut_in_tmp = rcut_in_BH_BH; 
		    eps2_tmp = eps2_BH_BH;
		}
		else if(prti->index >= NBH_GLB && prtj->index >= NBH_GLB){
		    rcut_out_tmp = rcut_out_FS_FS; 
		    rcut_in_tmp = rcut_in_FS_FS; 
		    eps2_tmp = eps2_FS_FS;
		}
		else if(prti->index < NBH_GLB || prtj->index < NBH_GLB){
		    rcut_out_tmp = rcut_out_FS_BH;
		    rcut_in_tmp = rcut_in_FS_BH; 
		    eps2_tmp = eps2_FS_BH;
		}
		//rcut_in_tmp = rcut_out_tmp*0.1;
		if(mode == 0){
		    calc_acc_jrk_duncan4(prti->pos,  prti->vel, 
					 acc_short_send,  jrk_short_send,
					 prti->pot_short, 
					 prtj->pos,  prtj->vel, 
					 prtj->mass,
					 eps2_tmp,  rcut_in_tmp,  rcut_out_tmp, 
					 r2_tmp);
		}
		else if(mode == 1){
		    if(prtj->time_pre != Tnext){
			double dt = Tnext - prtj->time;
			prtj->predict_short_h4(dt);
			prtj->time_pre = Tnext;
		    }
		    calc_acc_jrk_duncan4(prti->pos_pre,  prti->vel_pre, 
					 acc_short_send,  jrk_short_send,
					 prti->pot_short, 
					 prtj->pos_pre,  prtj->vel_pre, 
					 prtj->mass,
					 eps2_tmp,  rcut_in_tmp,  rcut_out_tmp, 
					 r2_tmp);
		}
		if(prtj->index < NBH_GLB && r2_tmp < prti->r2_ngh_BH){
		    dr2_rank_short_BH_send.dr2 = r2_tmp;
		    idx_ngh_short_BH = prtj->index;
		}
		if(prtj->index >= NBH_GLB && r2_tmp < prti->r2_ngh_FS){
		    dr2_rank_short_FS_send.dr2 = r2_tmp;
		    idx_ngh_short_FS = prtj->index;
		}
	    }
	
	    MPI_Allreduce(&acc_short_send, &acc_short_recv, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(&jrk_short_send, &jrk_short_recv, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(&dr2_rank_short_FS_send, &dr2_rank_short_FS_recv, 1,  MPI_DOUBLE_INT,  MPI_MINLOC, MPI_COMM_WORLD);
	    MPI_Allreduce(&dr2_rank_short_BH_send, &dr2_rank_short_BH_recv, 1,  MPI_DOUBLE_INT,  MPI_MINLOC, MPI_COMM_WORLD);

	    prti->acc_short = acc_short_recv;
	    prti->jrk_short = jrk_short_recv;
	    prti->r2_ngh_BH = r2_tmp;

	    MPI_Bcast(&dr2_rank_short_FS_recv,  1,  MPI_DOUBLE_INT,  dr2_rank_short_FS_recv.rank,  MPI_COMM_WORLD);
	    MPI_Bcast(&dr2_rank_short_BH_recv,  1,  MPI_DOUBLE_INT,  dr2_rank_short_BH_recv.rank,  MPI_COMM_WORLD);
	    MPI_Bcast(&idx_ngh_short_FS,  1,  MPI_INT,  dr2_rank_short_FS_recv.rank,  MPI_COMM_WORLD);
	    MPI_Bcast(&idx_ngh_short_BH,  1,  MPI_INT,  dr2_rank_short_BH_recv.rank,  MPI_COMM_WORLD);
	}

    }
    return 0;
}



int Hard_System::evaluate_acc_jrk_short_sirial(const int &Ni,
					       const double &Tnext,
					       const int &mode){
    double eps2_tmp = 0.0;
    double r2_tmp = 0.0;
    double rcut_in_tmp = 0.0;
    double rcut_out_tmp = 0.0;

    if(mode == 1){
	for(int i=0; i<Ni; i++){
	    Particle_Short *prti = prt_short+adr_time[i];
	    double dt = Tnext - prti->time;
	    prti->predict_short_h4(dt);
	    prti->time_pre = Tnext;
	}
    }

    for(int i=0; i<Ni; i++){
	Particle_Short *prti = prt_short+adr_time[i];
	int Nj = prti->Nj;

	if(mode == 1){
	    if(prti->time_pre != Tnext){
		double dt = Tnext - prti->time;
		prti->predict_short_h4(dt);
		prti->time_pre = Tnext;
	    }
	}
	prti->acc_short_old = prti->acc_short;
	prti->jrk_short_old = prti->jrk_short;
	prti->acc_short = prti->jrk_short = 0.0;
	prti->pot_short = 0.0;
	prti->r2_ngh_FS = prti->r2_ngh_BH = LARGEDOUBLE;
	for(int j=0; j<Nj; j++){
	    Particle_Short *prtj = ((prti->ngh_list_first)+j)->prtj;
	    if(prti->index == prtj->index){continue;}
	    if(prti->index < NBH_GLB && prtj->index < NBH_GLB){
		rcut_out_tmp = rcut_out_BH_BH; 
		rcut_in_tmp = rcut_in_BH_BH; 
		eps2_tmp = eps2_BH_BH;
	    }
	    else if(prti->index >= NBH_GLB && prtj->index >= NBH_GLB){
		rcut_out_tmp = rcut_out_FS_FS; 
		rcut_in_tmp = rcut_in_FS_FS; 
		eps2_tmp = eps2_FS_FS;
	    }
	    else if(prti->index < NBH_GLB || prtj->index < NBH_GLB){
		rcut_out_tmp = rcut_out_FS_BH;
		rcut_in_tmp = rcut_in_FS_BH; 
		eps2_tmp = eps2_FS_BH;
	    }
	    if(mode == 0){
		calc_acc_jrk_duncan4(prti->pos,  prti->vel, 
				     prti->acc_short,  prti->jrk_short,
				     prti->pot_short, 
				     prtj->pos,  prtj->vel, 
				     prtj->mass,
				     eps2_tmp,  rcut_in_tmp,  rcut_out_tmp, 
				     r2_tmp);
	    }
	    else if(mode == 1){
		if(prtj->time_pre != Tnext){
		    double dt = Tnext - prtj->time;
		    prtj->predict_short_h4(dt);
		    prtj->time_pre = Tnext;
		}
		calc_acc_jrk_duncan4(prti->pos_pre,  prti->vel_pre, 
				     prti->acc_short,  prti->jrk_short,
				     prti->pot_short, 
				     prtj->pos_pre,  prtj->vel_pre, 
				     prtj->mass,
				     eps2_tmp,  rcut_in_tmp,  rcut_out_tmp, 
				     r2_tmp);
	    }
	    if(prtj->index < NBH_GLB && r2_tmp < prti->r2_ngh_BH){
		prti->r2_ngh_BH = r2_tmp;
		prti->idx_ngh_BH = prtj->index;
	    }
	    if(prtj->index >= NBH_GLB && r2_tmp < prti->r2_ngh_FS){
		prti->r2_ngh_FS = r2_tmp;
		prti->idx_ngh_FS = prtj->index;
	    }
	}
    }
    return 0;
}


//////////////////////////////////////////////////////
//                                                  //
//            4th-order HERMITE INTEGRATION         //
//                                                  //
//////////////////////////////////////////////////////

static double maxdt_factor = 1.0/4.0;
int Hard_System::hermite4_2body(Particle_Short* prt_i,
				Particle_Short* prt_j,
				const double& Tsys,
				const double& dt_glb){

    prt_i->time = prt_j->time = Tsys;

    double Tend_tmp = Tsys+dt_glb;
    double dt_max = dt_glb*maxdt_factor;
    double Toffset = Tsys;
    double Tsync = Tend_tmp;
    int loop = 0;

    double eta_tmp = 0.0;
    double eta_s_tmp = eta_s;
    if(prt_i->index < NBH_GLB || prt_j->index < NBH_GLB){
	eta_tmp = eta_BH;
	//eta_s_tmp = eta_s*0.25;
    }
    else{
	eta_tmp = eta_FS;
    }

    double rcut_out_tmp = 0.0;
    double rcut_in_tmp = 0.0;
    double eps2_tmp = 0.0;
    if(prt_i->index < NBH_GLB && prt_j->index < NBH_GLB){
	rcut_out_tmp = rcut_out_BH_BH; 
	rcut_in_tmp = rcut_in_BH_BH; 
	eps2_tmp = eps2_BH_BH;
    }
    else if(prt_i->index >= NBH_GLB && prt_j->index >= NBH_GLB){
	rcut_out_tmp = rcut_out_FS_FS; 
	rcut_in_tmp = rcut_in_FS_FS; 
	eps2_tmp = eps2_FS_FS;
    }
    else if(prt_i->index < NBH_GLB || prt_j->index < NBH_GLB){
	rcut_out_tmp = rcut_out_FS_BH;
	rcut_in_tmp = rcut_in_FS_BH; 
	eps2_tmp = eps2_FS_BH;
    }

    double r2_tmp = 0.0; 
    while(1){
	if(loop == 0){
	    prt_i->set_dt_2nd(dt_max, eta_s_tmp, Toffset, Tsync, acc_offset_sq);
	}
	else{
	    prt_i->set_dt_4th(dt_max, eta_tmp, Toffset, Tsync, acc_offset_sq);
	}
	prt_j->delta_time = prt_i->delta_time;
	Tnext[0] = prt_i->time + prt_i->delta_time;

	prt_i->predict_short_h4(prt_i->delta_time);
	prt_j->predict_short_h4(prt_j->delta_time);
	prt_i->time_pre = prt_j->time_pre = Tnext[0];

	prt_i->acc_short_old = prt_i->acc_short;
	prt_i->jrk_short_old = prt_i->jrk_short;

	prt_j->acc_short_old = prt_j->acc_short;
	prt_j->jrk_short_old = prt_j->jrk_short;

	prt_i->acc_short = prt_i->jrk_short = prt_j->acc_short = prt_j->jrk_short = 0.0;
	prt_i->pot_short = prt_j->pot_short = 0.0;

	prt_i->r2_ngh_FS = prt_i->r2_ngh_BH = prt_j->r2_ngh_FS = prt_j->r2_ngh_BH = LARGEDOUBLE;

	calc_acc_jrk_duncan4(prt_i->pos_pre,  prt_i->vel_pre, 
			     prt_i->acc_short,  prt_i->jrk_short,
			     prt_i->pot_short, 
			     prt_j->pos_pre,  prt_j->vel_pre, 
			     prt_j->mass,
			     eps2_tmp,  rcut_in_tmp,  rcut_out_tmp, 
			     r2_tmp);

	prt_j->acc_short = -1.0*prt_i->acc_short;
	prt_j->jrk_short = -1.0*prt_i->jrk_short;
	prt_j->pot_short = prt_i->pot_short;

	prt_i->step++;
	prt_j->step++;
	STEPS_HARD += 2.0;

	if(prt_i->index < NBH_GLB){
	    prt_j->r2_ngh_BH = r2_tmp;
	}
	else{
	    prt_j->r2_ngh_FS = r2_tmp;
	}

	if(prt_j->index < NBH_GLB){
	    prt_i->r2_ngh_BH = r2_tmp;
	}
	else{
	    prt_i->r2_ngh_FS = r2_tmp;
	}

	prt_i->pos_old = prt_i->pos;
	prt_i->vel_old = prt_i->vel;
	prt_i->correct_h4();
	prt_i->interpolate_h4();
	prt_i->time += prt_i->delta_time;

	prt_j->pos_old = prt_j->pos;
	prt_j->vel_old = prt_j->vel;
	prt_j->correct_h4();
	prt_j->interpolate_h4();
	prt_j->time += prt_j->delta_time;

	loop++;
	if(Tnext[0] == Tend_tmp){
	    break;
	}
    }
}

class Particle_Short_Comm{
    Vector3 pos;
    Vector3 vel;
public:
    int adr;
    void set(const Particle_Short &prt, const int &_adr){
	pos = prt.pos;
	vel = prt.vel;
	adr = _adr;
    }
    void give(Particle_Short &prt){
	prt.pos = pos;
	prt.vel = vel;
    }
};

int Hard_System::hermite4_para(const double &Tsys,
			       const double &dt_glb,
			       const int &Nshort_glb){
    if(Nshort_glb == 0){ return 0;}
    double Tend_tmp = Tsys+dt_glb;
    int loop = 0;
    double eta_tmp = 0.0;
    double eta_s_tmp = eta_s;
    int Ni = Nshort_glb;
    double dt_max = dt_glb*maxdt_factor;
    double Toffset = Tsys;
    double Tsync = Tend_tmp;
    for(int i=0; i<Ni; i++){
	prt_short[i].time = Tsys;
	prt_short[i].step = 0;
	adr_time[i] = i;
	Tnext[i] = LARGEDOUBLE;
    }
    STEPS_HARD = 0.0;

#ifdef USING_SSE
    evaluate_acc_jrk_short_jpara_sse(Ni, Tsys, 0);
#else
    //evaluate_acc_jrk_short_sirial(Ni, Tsys, 0);
    evaluate_acc_jrk_short_jpara(Ni, Tsys, 0);
#endif //USING_SSE

    double toffset = mpi_wtime();
    /////////////////////////////
    // ipara for 2body system only
    int Npair_glb = 0;
    int Npair_loc = 0;
    static int adr_2body_loc[NSHORT_MAX];
    for(int i=0; i<Ni; i++){
	Particle_Short* prt_i = prt_short+i;
	if(prt_i->Nj == 1 && prt_i->time != Tend_tmp){
	    Particle_Short* prt_j = prt_i->ngh_list_first->prtj;
	    if( prt_j->Nj == 1 && prt_j->time != Tend_tmp){
		if(Npair_glb % NPROC == MYRANK){
		    adr_2body_loc[Npair_loc] = i;
		    Npair_loc ++;
		}
		prt_i->time = prt_j->time = Tend_tmp;
		Npair_glb ++;
	    }
	}
    }

    /////////////////////////////
    // exchange pos & vel of particles
    NGH_LIST_LEN_2BODY = Npair_glb*2;
    static Particle_Short_Comm prt_short_comm_send[NSHORT_MAX];
    static Particle_Short_Comm prt_short_comm_recv[NSHORT_MAX];
    for(int i=0; i<Npair_loc; i++){
	Particle_Short* prt_i = prt_short+adr_2body_loc[i];
	Particle_Short* prt_j = prt_i->ngh_list_first->prtj;
	prt_i->time = Tsys;
	prt_j->time = Tsys;
	hermite4_2body(prt_i, prt_j, Tsys, dt_glb);
	prt_short_comm_send[2*i].set(*prt_i, adr_2body_loc[i]);
	prt_short_comm_send[2*i+1].set(*prt_j, -1);
    }
    int Ni_loc = Npair_loc*2;
    int Ni_glb = Npair_glb*2;
    mpi_allgatherv_T(prt_short_comm_send, Ni_loc, prt_short_comm_recv, Ni_glb, NPROC);
    for(int i=0; i<Npair_glb; i++){
	int adr_i = prt_short_comm_recv[2*i].adr;
	Particle_Short* prt_i = prt_short+adr_i;
	prt_short_comm_recv[2*i].give(*prt_i);
	Particle_Short* prt_j = prt_i->ngh_list_first->prtj;
	prt_short_comm_recv[2*i+1].give(*prt_j);
	prt_i->time = prt_j->time = Tend_tmp;
    }
    double steps_send = STEPS_HARD;
    MPI_Allreduce(&steps_send, &STEPS_HARD, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dump_cout(STEPS_HARD);
    TCAL_HERMITE4_2BODY = mpi_wtime() - toffset;
    // exchange pos & vel of particles 
    /////////////////////////////


    ////////////////////////////////
    // j-para for multi body system
    toffset = mpi_wtime();
    while(1){
	for(int i=0; i<Ni; i++){
	    int adr = adr_time[i] ;
	    Particle_Short *prti = prt_short + adr;
	    if(prti->time == Tend_tmp){
		Tnext[adr] = prti->time + 0.125; // to skip integration for pair paritcles whose neighbours are each other
		continue;
	    }
	    if(prti->index < NBH_GLB){
		eta_tmp = eta_BH;
		//eta_s_tmp = eta_s*0.25;
	    }
	    else{
		eta_tmp = eta_FS;
	    }
	    if(loop == 0){
		prti->set_dt_2nd(dt_max, eta_s_tmp, Toffset, Tsync, acc_offset_sq);
	    }
	    else{
		prti->set_dt_4th(dt_max, eta_tmp, Toffset, Tsync, acc_offset_sq);
	    }
	    Tnext[adr] = prti->time + prti->delta_time;
	}
	sort_select_iparticle(Ni,  Nshort_glb,  adr_time,  Tnext);

#ifdef USING_SSE
	evaluate_acc_jrk_short_jpara_sse(Ni, Tnext[adr_time[0]], 1);
#else
	//evaluate_acc_jrk_short_sirial(Ni, Tnext[adr_time[0]], 1);
	evaluate_acc_jrk_short_jpara(Ni, Tnext[adr_time[0]], 1);
#endif //USING_SSE

	// all nodes do the same calculation...
	for(int i=0; i<Ni; i++){
	    Particle_Short *prti = prt_short+adr_time[i];
	    prti->pos_old = prti->pos;
	    prti->vel_old = prti->vel;
	    prti->correct_h4();
	    prti->interpolate_h4();
	    prti->time += prti->delta_time;
	    prti->step++;
	}
	STEPS_HARD += Ni;
	loop++;
	if(Tnext[adr_time[0]] == Tend_tmp){
	    cerr<<"Ni="<<Ni<<endl;
	    cerr<<"Nshort_glb="<<Nshort_glb<<endl;
	    cerr<<"loop="<<loop<<endl;
	    break;
	}
    }
    TCAL_HERMITE4_MULBODY = mpi_wtime() - toffset;


    //STEPS_HARD = 0.0;
    for(int i=0; i<Nshort_glb; i++){
	prt_short[i].acc_short = prt_short[i].jrk_short = 0.0;
	//STEPS_HARD += prt_short[i].step;
	prt_short[i].step = 0;
    }
    STEPS_HARD /= Nshort_glb;
    return 0;
}




// set prt_short
// index, mass, pos, vel, Nj are transfered
// but, imformations of neighbours (index, r2 etc...) are not transfered
// these values are meanless.
static int adr_from_idx_BH[NBH_GLB_MAX];
void Hard_System::merge_ngh_list(Particle prt_loc[],
				 Neighbour_List_Index ngh_list_idx[],
				 mymap adr_from_idx,
				 const int &Nshort_loc,
				 int &Nshort_glb,
				 const int &ngh_list_len_loc,
				 int &ngh_list_len_glb,
				 Neighbour_List ngh_list[]){
    static Neighbour_List_Index ngh_list_idx_recv[NGH_LIST_MAX];
    static mymap adr_from_idx_short;
    ////////////////////////////////
    /// allgathe  short particle 
    int Ncnt = 0;
    for(int i=0; i<Nshort_loc; i++){
	int adr = adr_from_idx[ngh_list_idx[Ncnt].idx_i];
	prt_comm_send[i].set(prt_loc[adr], adr);
	Ncnt += prt_loc[adr].Nj;
    }

    mpi_allgatherv_T_2(prt_comm_send, Nshort_loc, prt_comm_recv, Nshort_glb, NPROC, Nshort_disp);

    adr_from_idx_short.clear();
    for(int i=0; i<NBH_GLB_ORG; i++){
	adr_from_idx_BH[i] = -1;
    }
    for(int i=0; i<Nshort_glb; i++){
	prt_comm_recv[i].give(prt_short[i]);
	prt_short[i].ngh_list_first = NULL;
	prt_short[i].Nj = 0;
	adr_from_idx_short.insert(mymap::value_type(prt_short[i].index, i));
	if(prt_short[i].index < NBH_GLB_ORG){
	    adr_from_idx_BH[prt_short[i].index] = i;
	}
    }
    /// allgathe  short particle 
    ////////////////////////////////


    ////////////////////////////////
    /// allgathe ngh list 
    mpi_allgatherv_T(ngh_list_idx, ngh_list_len_loc, ngh_list_idx_recv, ngh_list_len_glb, NPROC);
    int idx_i_old = 0;
    for(int i=0; i<ngh_list_len_glb; i++){
	int idx_i_tmp = ngh_list_idx_recv[i].idx_i;
	int idx_j_tmp = ngh_list_idx_recv[i].idx_j;
	Particle_Short *prti_tmp = prt_short + adr_from_idx_short[idx_i_tmp];
	Particle_Short *prtj_tmp = prt_short + adr_from_idx_short[idx_j_tmp];
	prti_tmp->Nj++;
	if(i == 0 || idx_i_old != idx_i_tmp){
	    prti_tmp->ngh_list_first = ngh_list + i;
	}
	ngh_list[i].prti = prti_tmp;
	ngh_list[i].prtj = prtj_tmp;
	idx_i_old = idx_i_tmp;

	// initialize (in most cases, not needed)
	prti_tmp->r2_ngh_FS = prti_tmp->r2_ngh_BH = prtj_tmp->r2_ngh_FS = prtj_tmp->r2_ngh_BH = 0.0;
	prti_tmp->idx_ngh_FS = prti_tmp->idx_ngh_BH = prtj_tmp->idx_ngh_FS = prtj_tmp->idx_ngh_BH = 0;

    }
    /// allgathe ngh list 
    ////////////////////////////////
}


void Hard_System::copy_prt_short_to_prt_loc(Particle prt_loc[],
					    Particle BH_glb[]){
    for(int i=Nshort_disp[MYRANK]; i<Nshort_disp[MYRANK+1]; i++){
	int adr = prt_comm_recv[i].adr_org;
	prt_loc[adr].mass = prt_short[i].mass;
	prt_loc[adr].pos = prt_short[i].pos;
	prt_loc[adr].vel = prt_short[i].vel;
    }

    for(int i=0; i<NBH_GLB_ORG; i++){
	if(adr_from_idx_BH[i] == -1){continue;}
	BH_glb[i].mass = prt_short[adr_from_idx_BH[i]].mass;
	BH_glb[i].pos = prt_short[adr_from_idx_BH[i]].pos;
	BH_glb[i].vel = prt_short[adr_from_idx_BH[i]].vel;
    }
}




#ifdef USING_SSE

//////////////////////////////////////////////
//                                          //
//            CALC FORCES USIGN SSE         //
//                                          //
//////////////////////////////////////////////

typedef int v4si __attribute__ ((vector_size(16)));
typedef float v4sf __attribute__ ((vector_size(16)));
typedef double v2df __attribute__ ((vector_size(16)));

v4sf inv(const v4sf &val){
    v4sf tmp = __builtin_ia32_rcpps(val);
    return tmp*( (v4sf){2.0, 2.0, 2.0, 2.0} - tmp*val);
    //return v4sf(1.0f)/(*this);
}
v2df inv(const v2df &val){
    //v4sf x = __builtin_ia32_cvtpd2ps(val);
    //return __builtin_ia32_cvtps2pd(inv(x));
    return (v2df){1.0, 1.0}/val;
}

v4sf rsqrt(const v4sf &val){
    v4sf x = val;
    v4sf y = __builtin_ia32_rsqrtps(x);
    return ( (v4sf){-0.5, -0.5, -0.5, -0.5} * y) * (x*y*y - (v4sf){3.0, 3.0, 3.0, 3.0});
}

v2df rsqrt(const v2df &val){
    v4sf x = __builtin_ia32_cvtpd2ps(val);
    return __builtin_ia32_cvtps2pd(rsqrt(x));
}

float sum(const v4sf &val){
    v4sf x = __builtin_ia32_haddps(val, val);
    x = __builtin_ia32_haddps(x, x);
    return __builtin_ia32_vec_ext_v4sf(x, 0); 
}

double sum(const v2df &val){
    v2df x = __builtin_ia32_haddpd(val, val);
    return __builtin_ia32_vec_ext_v2df(x, 0); 
}


v2df pre_pos(const v2df &x, const v2df &v, const v2df &a, const v2df &j, const v2df &dt){
    static v2df HALF = (v2df){0.5, 0.5};
    static v2df OVER3 = (v2df){1.0/3.0, 1.0/3.0};
    return ( (j*dt*OVER3 + a)*dt*HALF + v)*dt + x; 
}
v2df pre_vel(const v2df &v, const v2df &a, const v2df &j, const v2df &dt){
    static v2df HALF = (v2df){0.5, 0.5};
    return (j*dt*HALF + a)*dt + v;
}

// NOTE: Order of predictor calculation.
// Initailly, all Ni particles are predicted .
// During force loop on an i particle, 
// if a j particle has been alreaddy predicted, 
// this j particle must not be predicted, because
// if this j particle is already integrated, 
// acc and jrk of this particle have been new values at t+dt.
// Thus each particle must be predicted once.
int Hard_System::evaluate_acc_jrk_short_jpara_sse(const int &Ni,
						  const double &Tnext,
						  const int &mode){
    static Vector3 acc_short_send;
    static Vector3 acc_short_recv;
    static Vector3 jrk_short_send;
    static Vector3 jrk_short_recv;
    static int idx_ngh_short_FS;
    static int idx_ngh_short_BH;
    static dr2_rank dr2_rank_short_FS_send;
    static dr2_rank dr2_rank_short_FS_recv;
    static dr2_rank dr2_rank_short_BH_send;
    static dr2_rank dr2_rank_short_BH_recv;
    double eps2_tmp = 0.0;
    double r2_tmp = 0.0;
    double rcut_in_tmp = 0.0;
    double rcut_out_tmp = 0.0;

    if(mode == 1){
	// all i particles are predicted with double precision
	v2df xp, vxp, yp, vyp, zp, vzp;
	v2df xold, vold, aold, jold, dt;
	for(int i=0; i<Ni; i+=2){
	    int ni_tmp = i+2 < Ni ? 2 : Ni-i;
	    Particle_Short *i_prt0 = prt_short+adr_time[i];
	    if(ni_tmp == 2){
		Particle_Short *i_prt1 = prt_short+adr_time[i+1];
		dt = (v2df){Tnext-i_prt0->time, Tnext-i_prt1->time}; 
		xold = (v2df){i_prt0->pos[0], i_prt1->pos[0]};
		vold = (v2df){i_prt0->vel[0], i_prt1->vel[0]};
		aold = (v2df){i_prt0->acc_short[0], i_prt1->acc_short[0]};
		jold = (v2df){i_prt0->jrk_short[0], i_prt1->jrk_short[0]};
		xp = pre_pos(xold, vold, aold, jold, dt);
		vxp = pre_vel(vold, aold, jold, dt);

		xold = (v2df){i_prt0->pos[1], i_prt1->pos[1]};
		vold = (v2df){i_prt0->vel[1], i_prt1->vel[1]};
		aold = (v2df){i_prt0->acc_short[1], i_prt1->acc_short[1]};
		jold = (v2df){i_prt0->jrk_short[1], i_prt1->jrk_short[1]};
		yp = pre_pos(xold, vold, aold, jold, dt);
		vyp = pre_vel(vold, aold, jold, dt);

		xold = (v2df){i_prt0->pos[2], i_prt1->pos[2]};
		vold = (v2df){i_prt0->vel[2], i_prt1->vel[2]};
		aold = (v2df){i_prt0->acc_short[2], i_prt1->acc_short[2]};
		jold = (v2df){i_prt0->jrk_short[2], i_prt1->jrk_short[2]};
		zp = pre_pos(xold, vold, aold, jold, dt);
		vzp = pre_vel(vold, aold, jold, dt);

		i_prt0->pos_pre[0] = __builtin_ia32_vec_ext_v2df(xp, 0);
		i_prt1->pos_pre[0] = __builtin_ia32_vec_ext_v2df(xp, 1);
		i_prt0->vel_pre[0] = __builtin_ia32_vec_ext_v2df(vxp, 0);
		i_prt1->vel_pre[0] = __builtin_ia32_vec_ext_v2df(vxp, 1);

		i_prt0->pos_pre[1] = __builtin_ia32_vec_ext_v2df(yp, 0);
		i_prt1->pos_pre[1] = __builtin_ia32_vec_ext_v2df(yp, 1);
		i_prt0->vel_pre[1] = __builtin_ia32_vec_ext_v2df(vyp, 0);
		i_prt1->vel_pre[1] = __builtin_ia32_vec_ext_v2df(vyp, 1);

		i_prt0->pos_pre[2] = __builtin_ia32_vec_ext_v2df(zp, 0);
		i_prt1->pos_pre[2] = __builtin_ia32_vec_ext_v2df(zp, 1);
		i_prt0->vel_pre[2] = __builtin_ia32_vec_ext_v2df(vzp, 0);
		i_prt1->vel_pre[2] = __builtin_ia32_vec_ext_v2df(vzp, 1);

		i_prt0->time_pre = Tnext;
		i_prt1->time_pre = Tnext;
	    }
	    else if(ni_tmp == 1){
		double dt = Tnext - i_prt0->time;
		i_prt0->predict_short_h4(dt);
		i_prt0->time_pre = Tnext;
	    }
	}
    }

    for(int i=0; i<Ni; i++){
	Particle_Short* i_prt = prt_short+adr_time[i];
	int Nj = i_prt->Nj;
	i_prt->acc_short_old = i_prt->acc_short;
	i_prt->jrk_short_old = i_prt->jrk_short;
	i_prt->acc_short = i_prt->jrk_short = 0.0;
	i_prt->pot_short = 0.0;
	i_prt->r2_ngh_FS = i_prt->r2_ngh_BH = LARGEDOUBLE;

	if(Nj < NJCRIT){
	    // all nodes execute the same calculations
	    int j_head = 0;
	    int j_tale = Nj;
	    //calc_acc_jrk_duncan4_from_jarray_scholar(i_prt, j_head, j_tale, Tnext, mode);
	    calc_acc_jrk_duncan4_from_jarray_sse_fulld(i_prt, j_head, j_tale, Tnext, mode);
	    //calc_acc_jrk_duncan4_from_jarray_sse_mix(i_prt, j_head, j_tale, Tnext, mode);
	}
	else{
	    acc_short_send = 0.0;
	    jrk_short_send = 0.0;

	    dr2_rank_short_FS_send.dr2 = dr2_rank_short_BH_send.dr2 = LARGEDOUBLE;
	    dr2_rank_short_FS_send.rank = dr2_rank_short_BH_send.rank = MYRANK;

	    int Nproc_new = (Nj/NJCRIT+1) >= NPROC ? NPROC : (Nj/NJCRIT+1);
	    int Nj_loc = Nj/Nproc_new;
	    int j_head_loc = Nj_loc*MYRANK;
	    if(MYRANK < Nproc_new){
		if(MYRANK < Nj%Nproc_new){
		    Nj_loc++;
		    j_head_loc += MYRANK;
		}
		else{
		    j_head_loc += Nj%Nproc_new;
		}
	    }
	    else{
		Nj_loc = 0;
		j_head_loc = 0;
	    }
	    int j_tale_loc = j_head_loc + Nj_loc; 
	    //calc_acc_jrk_duncan4_from_jarray_scholar(i_prt, j_head_loc, j_tale_loc, Tnext, mode);
	    calc_acc_jrk_duncan4_from_jarray_sse_fulld(i_prt, j_head_loc, j_tale_loc, Tnext, mode);
	    //calc_acc_jrk_duncan4_from_jarray_sse_mix(i_prt, j_head_loc, j_tale_loc, Tnext, mode);

	    MPI_Allreduce(&acc_short_send, &acc_short_recv, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(&jrk_short_send, &jrk_short_recv, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(&dr2_rank_short_FS_send, &dr2_rank_short_FS_recv, 1,  MPI_DOUBLE_INT,  MPI_MINLOC, MPI_COMM_WORLD);
	    MPI_Allreduce(&dr2_rank_short_BH_send, &dr2_rank_short_BH_recv, 1,  MPI_DOUBLE_INT,  MPI_MINLOC, MPI_COMM_WORLD);

	    i_prt->acc_short = acc_short_recv;
	    i_prt->jrk_short = jrk_short_recv;
	    i_prt->r2_ngh_BH = r2_tmp;

	    MPI_Bcast(&dr2_rank_short_FS_recv,  1,  MPI_DOUBLE_INT,  dr2_rank_short_FS_recv.rank,  MPI_COMM_WORLD);
	    MPI_Bcast(&dr2_rank_short_BH_recv,  1,  MPI_DOUBLE_INT,  dr2_rank_short_BH_recv.rank,  MPI_COMM_WORLD);
	    MPI_Bcast(&idx_ngh_short_FS,  1,  MPI_INT,  dr2_rank_short_FS_recv.rank,  MPI_COMM_WORLD);
	    MPI_Bcast(&idx_ngh_short_BH,  1,  MPI_INT,  dr2_rank_short_BH_recv.rank,  MPI_COMM_WORLD);
	}

    }
    return 0;
}


void Hard_System::calc_acc_jrk_duncan4_from_jarray_sse_mix(Particle_Short* i_prt,
							   const int& j_head,
							   const int& j_tale,
							   const double& Tnext,
							   const int& mode){

    static v4sf ZERO = (v4sf){0.0, 0.0, 0.0, 0.0};
    static v4sf ONE = (v4sf){1.0, 1.0, 1.0, 1.0};
    static v4sf THREE = (v4sf){3.0, 3.0, 3.0, 3.0};
    static v4sf v4sf_20 = (v4sf){20.0, 20.0, 20.0, 20.0};
    static v4sf v4sf_70 = (v4sf){70.0, 70.0, 70.0, 70.0};
    static v4sf v4sf_84 = (v4sf){84.0, 84.0, 84.0, 84.0};
    static v4sf v4sf_35 = (v4sf){35.0, 35.0, 35.0, 35.0};
    static v4sf v4sf_140 = (v4sf){140.0, 140.0, 140.0, 140.0};
    static v4sf v4sf_420 = (v4sf){420.0, 420.0, 420.0, 420.0};
    static v4sf v4sf_NBH_GLB_ORG = (v4sf){(float)NBH_GLB_ORG*0.9999999999, (float)NBH_GLB_ORG*0.9999999999, (float)NBH_GLB_ORG*0.9999999999, (float)NBH_GLB_ORG*0.9999999999};

    v4sf eps_sq_FS_FS = (v4sf){eps2_FS_FS, eps2_FS_FS};
    v4sf eps_sq_FS_BH = (v4sf){eps2_FS_BH, eps2_FS_BH};
    v4sf eps_sq_BH_BH = (v4sf){eps2_BH_BH, eps2_BH_BH};
    v4sf rout_FS_FS = (v4sf){rcut_out_FS_FS, rcut_out_FS_FS};
    v4sf rout_FS_BH = (v4sf){rcut_out_FS_BH, rcut_out_FS_BH};
    v4sf rout_BH_BH = (v4sf){rcut_out_BH_BH, rcut_out_BH_BH};
    v4sf rin_FS_FS = (v4sf){rcut_in_FS_FS, rcut_in_FS_FS};
    v4sf rin_FS_BH = (v4sf){rcut_in_FS_BH, rcut_in_FS_BH};
    v4sf rin_BH_BH = (v4sf){rcut_in_BH_BH, rcut_in_BH_BH};

    v4sf eps_sq, rout, rin;
    v2df tnext = (v2df){Tnext, Tnext};

    v2df x_pi, y_pi, z_pi;
    v4sf x_vi, y_vi, z_vi;
    if(mode == 0){
	double pi = i_prt->pos[0];
	x_pi = (v2df){pi, pi};
	pi = i_prt->pos[1];
	y_pi = (v2df){pi, pi};
	pi = i_prt->pos[2];
	z_pi = (v2df){pi, pi};

	float vi = i_prt->vel[0];
	x_vi = (v4sf){vi, vi, vi, vi};
	vi = i_prt->vel[1];
	y_vi = (v4sf){vi, vi, vi, vi};
	vi = i_prt->vel[2];
	z_vi = (v4sf){vi, vi, vi, vi};
    }
    else if(mode == 1){
	double pi = i_prt->pos_pre[0];
	x_pi = (v2df){pi, pi};
	pi = i_prt->pos_pre[1];
	y_pi = (v2df){pi, pi};
	pi = i_prt->pos_pre[2];
	z_pi = (v2df){pi, pi};

	float vi = i_prt->vel_pre[0];
	x_vi = (v4sf){vi, vi, vi, vi};
	vi = i_prt->vel_pre[1];
	y_vi = (v4sf){vi, vi, vi, vi};
	vi = i_prt->vel_pre[2];
	z_vi = (v4sf){vi, vi, vi, vi};
    }

    v4sf x_ai, y_ai, z_ai, x_ji, y_ji, z_ji, pot;
    x_ai = y_ai = z_ai = x_ji = y_ji = z_ji = pot = ZERO;

    v4sf idx_pj;
    v2df x_pj0, x_pj1, y_pj0, y_pj1, z_pj0, z_pj1;
    v4sf x_vj, y_vj, z_vj, m_j;
    Neighbour_List* j_first = i_prt->ngh_list_first;
    for(int j = j_head; j < j_tale; j += 4){
	int nj_tmp = j+4 < j_tale ? 4 : j_tale - j; 
	// set pos and vel of j paritcles
	if(nj_tmp == 4){
	    m_j = (v4sf){(float)(j_first+j)->prtj->mass, (float)(j_first+j+1)->prtj->mass, (float)(j_first+j+2)->prtj->mass, (float)(j_first+j+3)->prtj->mass};
	    idx_pj = (v4sf){(j_first+j)->prtj->index, (j_first+j+1)->prtj->index, (j_first+j)->prtj->index, (j_first+j+1)->prtj->index};
	}
	else if(nj_tmp == 3){
	    m_j = (v4sf){(float)(j_first+j)->prtj->mass, (float)(j_first+j+1)->prtj->mass, (float)(j_first+j+2)->prtj->mass, 0.0};
	    idx_pj = (v4sf){(j_first+j)->prtj->index, (j_first+j+1)->prtj->index, (j_first+j)->prtj->index, -1.0};
	}
	else if(nj_tmp == 2){
	    m_j = (v4sf){(float)(j_first+j)->prtj->mass, (float)(j_first+j+1)->prtj->mass, 0.0, 0.0};
	    idx_pj = (v4sf){(j_first+j)->prtj->index, (j_first+j+1)->prtj->index, -1.0, -1.0};
	}
	else if(nj_tmp == 1){
	    m_j = (v4sf){(float)(j_first+j)->prtj->mass, 0.0, 0.0, 0.0};
	    idx_pj = (v4sf){(j_first+j)->prtj->index, -1.0, -1.0, -1.0};
	}
	if(mode == 0){
	    if(nj_tmp == 4){
		x_pj0 = (v2df){(j_first+j)->prtj->pos[0],   (j_first+j+1)->prtj->pos[0]}; 
		x_pj1 = (v2df){(j_first+j+2)->prtj->pos[0], (j_first+j+3)->prtj->pos[0]};
		x_vj = (v4sf){(float)(j_first+j)->prtj->vel[0],   (float)(j_first+j+1)->prtj->vel[0],
			      (float)(j_first+j+2)->prtj->vel[0], (float)(j_first+j+3)->prtj->vel[0]};

		y_pj0 = (v2df){(j_first+j)->prtj->pos[1],   (j_first+j+1)->prtj->pos[1]}; 
		y_pj1 = (v2df){(j_first+j+2)->prtj->pos[1], (j_first+j+3)->prtj->pos[1]};
		y_vj = (v4sf){(float)(j_first+j)->prtj->vel[1],   (float)(j_first+j+1)->prtj->vel[1],
			      (float)(j_first+j+2)->prtj->vel[1], (float)(j_first+j+3)->prtj->vel[1]};

		z_pj0 = (v2df){(j_first+j)->prtj->pos[2],   (j_first+j+1)->prtj->pos[2]}; 
		z_pj1 = (v2df){(j_first+j+2)->prtj->pos[2], (j_first+j+3)->prtj->pos[2]};
		z_vj = (v4sf){(float)(j_first+j)->prtj->vel[2],   (float)(j_first+j+1)->prtj->vel[2],
			      (float)(j_first+j+2)->prtj->vel[2], (float)(j_first+j+3)->prtj->vel[2]};

	    }
	    else if(nj_tmp == 3){
		x_pj0 = (v2df){(j_first+j)->prtj->pos[0],    (j_first+j+1)->prtj->pos[0]}; 
		x_pj1 = (v2df){(j_first+j+2)->prtj->pos[0],  LARGEDOUBLE};
		x_vj = (v4sf){(float)(j_first+j)->prtj->vel[0],    (float)(j_first+j+1)->prtj->vel[0],
			      (float)(j_first+j+2)->prtj->vel[0],  LARGESINGLE};

		y_pj0 = (v2df){(j_first+j)->prtj->pos[1],    (j_first+j+1)->prtj->pos[1]}; 
		y_pj1 = (v2df){(j_first+j+2)->prtj->pos[1],  LARGEDOUBLE};
		y_vj = (v4sf){(float)(j_first+j)->prtj->vel[1],    (float)(j_first+j+1)->prtj->vel[1],
			      (float)(j_first+j+2)->prtj->vel[1],  LARGESINGLE};

		z_pj0 = (v2df){(j_first+j)->prtj->pos[2],    (j_first+j+1)->prtj->pos[2]}; 
		z_pj1 = (v2df){(j_first+j+2)->prtj->pos[2],  LARGEDOUBLE};
		z_vj = (v4sf){(float)(j_first+j)->prtj->vel[2],    (float)(j_first+j+1)->prtj->vel[2],
			      (float)(j_first+j+2)->prtj->vel[2],  LARGESINGLE};
	    }
	    else if(nj_tmp == 2){
		x_pj0 = (v2df){(j_first+j)->prtj->pos[0],  (j_first+j+1)->prtj->pos[0]}; 
		x_pj1 = (v2df){LARGEDOUBLE,  LARGEDOUBLE};
		x_vj = (v4sf){(float)(j_first+j)->prtj->vel[0],  (float)(j_first+j+1)->prtj->vel[0],  LARGESINGLE,  LARGESINGLE};

		y_pj0 = (v2df){(j_first+j)->prtj->pos[1],  (j_first+j+1)->prtj->pos[1]}; 
		y_pj1 = (v2df){LARGEDOUBLE,  LARGEDOUBLE};
		y_vj = (v4sf){(float)(j_first+j)->prtj->vel[1],  (float)(j_first+j+1)->prtj->vel[1],  LARGESINGLE,  LARGESINGLE};

		z_pj0 = (v2df){(j_first+j)->prtj->pos[2],  (j_first+j+1)->prtj->pos[2]}; 
		z_pj1 = (v2df){LARGEDOUBLE,  LARGEDOUBLE};
		z_vj = (v4sf){(float)(j_first+j)->prtj->vel[2],  (float)(j_first+j+1)->prtj->vel[2],  LARGESINGLE,  LARGESINGLE};
	    }
	    else if(nj_tmp == 1){
		x_pj0 = (v2df){(j_first+j)->prtj->pos[0],  LARGEDOUBLE};
		x_pj1 = (v2df){LARGEDOUBLE,  LARGEDOUBLE};
		x_vj = (v4sf){(float)(j_first+j)->prtj->vel[0],  LARGESINGLE,  LARGESINGLE,  LARGESINGLE};

		y_pj0 = (v2df){(j_first+j)->prtj->pos[1],  LARGEDOUBLE};
		y_pj1 = (v2df){LARGEDOUBLE,  LARGEDOUBLE};
		y_vj = (v4sf){(float)(j_first+j)->prtj->vel[1],  LARGESINGLE,  LARGESINGLE,  LARGESINGLE};

		z_pj0 = (v2df){(j_first+j)->prtj->pos[2], LARGEDOUBLE};
		z_pj1 = (v2df){LARGEDOUBLE,  LARGEDOUBLE};
		z_vj = (v4sf){(float)(j_first+j)->prtj->vel[2],  LARGESINGLE,  LARGESINGLE,  LARGESINGLE};
	    }
	}
	else if(mode == 1){
	    v2df t_pre0, t_pre1;
	    if(nj_tmp == 4){
		t_pre0 = (v2df){(j_first+j)->prtj->time_pre,   (j_first+j+1)->prtj->time_pre};
		t_pre1 = (v2df){(j_first+j+2)->prtj->time_pre, (j_first+j+3)->prtj->time_pre};
	    }
	    else if(nj_tmp == 3){
		t_pre0 = (v2df){(j_first+j)->prtj->time_pre,   (j_first+j+1)->prtj->time_pre};
		t_pre1 = (v2df){(j_first+j+2)->prtj->time_pre, Tnext};
	    }
	    else if(nj_tmp == 2){
		t_pre0 = (v2df){(j_first+j)->prtj->time_pre,   (j_first+j+1)->prtj->time_pre};
		t_pre1 = (v2df){Tnext,   Tnext};
	    }
	    else if(nj_tmp == 1){
		t_pre0 = (v2df){(j_first+j)->prtj->time_pre,    Tnext};
		t_pre1 = (v2df){Tnext,   Tnext};
	    }
	    v2df mask0 = (v2df)__builtin_ia32_cmpneqpd(t_pre0, tnext); // t_pre != tnext ? 0xffffffff : 0
	    v2df mask1 = (v2df)__builtin_ia32_cmpneqpd(t_pre1, tnext); // t_pre != tnext ? 0xffffffff : 0
	    int bits0 = __builtin_ia32_movmskpd(mask0);
	    int bits1 = __builtin_ia32_movmskpd(mask1);
	    if(bits0){
		if(bits0 & 1){
		    double dt = Tnext - (j_first+j)->prtj->time;
		    (j_first+j)->prtj->predict_short_h4(dt);
		    (j_first+j)->prtj->time_pre = Tnext;
		}
		if(bits0 & 2){
		    double dt = Tnext - (j_first+j+1)->prtj->time;
		    (j_first+j+1)->prtj->predict_short_h4(dt);
		    (j_first+j+1)->prtj->time_pre = Tnext;
		}
	    }
	    if(bits1){
		if(bits1 & 1){
		    double dt = Tnext - (j_first+j+2)->prtj->time;
		    (j_first+j+2)->prtj->predict_short_h4(dt);
		    (j_first+j+2)->prtj->time_pre = Tnext;
		}
		if(bits1 & 2){
		    double dt = Tnext - (j_first+j+3)->prtj->time;
		    (j_first+j+3)->prtj->predict_short_h4(dt);
		    (j_first+j+3)->prtj->time_pre = Tnext;
		}
	    }

	    if(nj_tmp == 4){
		x_pj0 = (v2df){(j_first+j)->prtj->pos_pre[0],    (j_first+j+1)->prtj->pos_pre[0]}; 
		x_pj1 = (v2df){(j_first+j+2)->prtj->pos_pre[0],  (j_first+j+3)->prtj->pos_pre[0]};
		x_vj = (v4sf){(float)(j_first+j)->prtj->vel_pre[0],   (float)(j_first+j+1)->prtj->vel_pre[0],
			      (float)(j_first+j+2)->prtj->vel_pre[0], (float)(j_first+j+3)->prtj->vel_pre[0]};

		y_pj0 = (v2df){(j_first+j)->prtj->pos_pre[1],    (j_first+j+1)->prtj->pos_pre[1]}; 
		y_pj1 = (v2df){(j_first+j+2)->prtj->pos_pre[1],  (j_first+j+3)->prtj->pos_pre[1]};
		y_vj = (v4sf){(float)(j_first+j)->prtj->vel_pre[1],   (float)(j_first+j+1)->prtj->vel_pre[1],
			      (float)(j_first+j+2)->prtj->vel_pre[1], (float)(j_first+j+3)->prtj->vel_pre[1]};

		z_pj0 = (v2df){(j_first+j)->prtj->pos_pre[2],    (j_first+j+1)->prtj->pos_pre[2]}; 
		z_pj1 = (v2df){(j_first+j+2)->prtj->pos_pre[2],  (j_first+j+3)->prtj->pos_pre[2]};
		z_vj = (v4sf){(float)(j_first+j)->prtj->vel_pre[2],   (float)(j_first+j+1)->prtj->vel_pre[2],
			      (float)(j_first+j+2)->prtj->vel_pre[2], (float)(j_first+j+3)->prtj->vel_pre[2]};
	    }
	    else if(nj_tmp == 3){
		x_pj0 = (v2df){(j_first+j)->prtj->pos_pre[0],    (j_first+j+1)->prtj->pos_pre[0]}; 
		x_pj1 = (v2df){(j_first+j+2)->prtj->pos_pre[0],  LARGEDOUBLE};
		x_vj = (v4sf){(float)(j_first+j)->prtj->vel_pre[0],    (float)(j_first+j+1)->prtj->vel_pre[0],
			      (float)(j_first+j+2)->prtj->vel_pre[0],  LARGESINGLE};

		y_pj0 = (v2df){(j_first+j)->prtj->pos_pre[1],    (j_first+j+1)->prtj->pos_pre[1]}; 
		y_pj1 = (v2df){(j_first+j+2)->prtj->pos_pre[1],  LARGEDOUBLE};
		y_vj = (v4sf){(float)(j_first+j)->prtj->vel_pre[1],    (float)(j_first+j+1)->prtj->vel_pre[1],
			      (float)(j_first+j+2)->prtj->vel_pre[1],  LARGESINGLE};

		z_pj0 = (v2df){(j_first+j)->prtj->pos_pre[2],    (j_first+j+1)->prtj->pos_pre[2]}; 
		z_pj1 = (v2df){(j_first+j+2)->prtj->pos_pre[2],  LARGEDOUBLE};
		z_vj = (v4sf){(float)(j_first+j)->prtj->vel_pre[2],    (float)(j_first+j+1)->prtj->vel_pre[2],
			      (float)(j_first+j+2)->prtj->vel_pre[2],  LARGESINGLE};
	    }
	    else if(nj_tmp == 2){
		x_pj0 = (v2df){(j_first+j)->prtj->pos_pre[0],  (j_first+j+1)->prtj->pos_pre[0]}; 
		x_pj1 = (v2df){LARGEDOUBLE,   LARGEDOUBLE};
		x_vj = (v4sf){(float)(j_first+j)->prtj->vel_pre[0],  (float)(j_first+j+1)->prtj->vel_pre[0],  LARGESINGLE,  LARGESINGLE};

		y_pj0 = (v2df){(j_first+j)->prtj->pos_pre[1],  (j_first+j+1)->prtj->pos_pre[1]}; 
		y_pj1 = (v2df){LARGEDOUBLE, LARGEDOUBLE};
		y_vj = (v4sf){(float)(j_first+j)->prtj->vel_pre[1],  (float)(j_first+j+1)->prtj->vel_pre[1],  LARGESINGLE,  LARGESINGLE};

		z_pj0 = (v2df){(j_first+j)->prtj->pos_pre[2],  (j_first+j+1)->prtj->pos_pre[2]}; 
		z_pj1 = (v2df){LARGEDOUBLE, LARGEDOUBLE};
		z_vj = (v4sf){(float)(j_first+j)->prtj->vel_pre[2],  (float)(j_first+j+1)->prtj->vel_pre[2],  LARGESINGLE,  LARGESINGLE};
	    }
	    else if(nj_tmp == 1){
		x_pj0 = (v2df){(j_first+j)->prtj->pos_pre[0],  LARGEDOUBLE};
		x_pj1 = (v2df){LARGEDOUBLE,  LARGEDOUBLE};
		x_vj = (v4sf){(float)(j_first+j)->prtj->vel_pre[0],  LARGESINGLE,  LARGESINGLE,  LARGESINGLE};

		y_pj0 = (v2df){(j_first+j)->prtj->pos_pre[1],  LARGEDOUBLE};
		y_pj1 = (v2df){LARGEDOUBLE,  LARGEDOUBLE};
		y_vj = (v4sf){(float)(j_first+j)->prtj->vel_pre[1],  LARGESINGLE,  LARGESINGLE,  LARGESINGLE};

		z_pj0 = (v2df){(j_first+j)->prtj->pos_pre[2],  LARGEDOUBLE};
		z_pj1 = (v2df){LARGEDOUBLE,  LARGEDOUBLE};
		z_vj = (v4sf){(float)(j_first+j)->prtj->vel_pre[2],  LARGESINGLE,  LARGESINGLE,  LARGESINGLE};
	    }
	}

	v4sf mask = (v4sf)__builtin_ia32_cmpltps(idx_pj, v4sf_NBH_GLB_ORG); // idx_pj < NBH_GLB_ORG ? 0xffffffff : 0
	if(i_prt->index < NBH_GLB_ORG){
	    eps_sq = __builtin_ia32_andps(mask, eps_sq_BH_BH) + __builtin_ia32_andnps(mask, eps_sq_FS_BH);
	    rout = __builtin_ia32_andps(mask, rout_BH_BH) + __builtin_ia32_andnps(mask, rout_FS_BH);
	    rin = __builtin_ia32_andps(mask, rin_BH_BH) + __builtin_ia32_andnps(mask, rin_FS_BH);
	}
	else{
	    eps_sq = __builtin_ia32_andps(mask, eps_sq_FS_BH) + __builtin_ia32_andnps(mask, eps_sq_FS_FS);
	    rout = __builtin_ia32_andps(mask, rout_FS_BH) + __builtin_ia32_andnps(mask, rout_FS_FS);
	    rin = __builtin_ia32_andps(mask, rin_FS_BH) + __builtin_ia32_andnps(mask, rin_FS_FS);
	}

	v4sf dx = __builtin_ia32_cvtpd2ps(x_pi - x_pj0);
	v4sf ftmp = __builtin_ia32_cvtpd2ps(x_pi - x_pj1);
	dx = __builtin_ia32_movlhps(dx, ftmp);

        v4sf dy = __builtin_ia32_cvtpd2ps(y_pi - y_pj0);
        ftmp = __builtin_ia32_cvtpd2ps(y_pi - y_pj1);
        dy = __builtin_ia32_movlhps(dy, ftmp);

        v4sf dz = __builtin_ia32_cvtpd2ps(z_pi - z_pj0);
        ftmp = __builtin_ia32_cvtpd2ps(z_pi - z_pj1);
        dz = __builtin_ia32_movlhps(dz, ftmp);

        v4sf r2 = dx*dx + dy*dy + dz*dz;
	v4sf rout_sq = rout*rout;
	v4sf mask_sf = (v4sf)__builtin_ia32_cmpltps(r2, rout_sq); // r2 < rout2 ? 0xffffffff : 0
	int bit = __builtin_ia32_movmskps(mask_sf);
	if(!bit){ continue; } // all j particles are outside the neighbour sphere

	v4sf r_wo_eps_inv = rsqrt(r2);
	v4sf r_wo_eps = inv(r_wo_eps_inv);

	r2 = r2 + eps_sq;
	v4sf r_inv = rsqrt(r2);
	v4sf r2_inv = r_inv*r_inv;
	v4sf r3_inv = r2_inv*r_inv;

        v4sf dvx = x_vi - x_vj;
        v4sf dvy = y_vi - y_vj;
        v4sf dvz = z_vi - z_vj;

        v4sf rv = dx*dvx + dy*dvy + dz*dvz;

	// evaluate change over function by Duncan for 4th-order
	v4sf dr_cut_inv = inv(rout-rin);
        v4sf rtmp = (r_wo_eps-rin)*dr_cut_inv;
        v4sf rtmpdot = rv*r_wo_eps_inv*dr_cut_inv;

        rtmp = __builtin_ia32_andps(rtmp, (v4sf)__builtin_ia32_cmpltps(ZERO, rtmp)); // rtmp = 0.0 < rtmp ?  rtmp : 0.0
	rtmp = __builtin_ia32_andps((rtmp-ONE), (v4sf)__builtin_ia32_cmpltps(rtmp, ONE)) + ONE; // rtmp = rtmp < 1.0 ?  rtmp : 1.0
	rtmpdot = __builtin_ia32_andps(rtmpdot, (v4sf)__builtin_ia32_cmpltps(ZERO, rtmp));  // rtmpdot = (0.0 < rtmp ?  rtmpdot : 0.0)
	rtmpdot = __builtin_ia32_andps(rtmpdot, (v4sf)__builtin_ia32_cmpltps(rtmp, ONE)); // rtmpdot = (rtmp < 1.0 ?  rtmpdot : 0.0)
        v4sf rtmp2 = rtmp*rtmp;
        v4sf rtmp3 = rtmp2*rtmp;
	v4sf rtmp4 = rtmp2*rtmp2;
        v4sf K = ( ( (-v4sf_20*rtmp + v4sf_70)*rtmp - v4sf_84)*rtmp + v4sf_35)*rtmp4;
        v4sf Kdot = ( ((( -v4sf_140*rtmp + v4sf_420)*rtmp - v4sf_420)*rtmp + v4sf_140)*rtmp3) * rtmpdot;


        v4sf A0 = m_j * r3_inv * (ONE-K);
        v4sf f0x = A0 * -dx;
        v4sf f0y = A0 * -dy;
        v4sf f0z = A0 * -dz;

        v4sf A1 = THREE*rv*r2_inv;
	v4sf A2 = m_j * r3_inv * Kdot;
        v4sf f1x = A0*-dvx - A1*f0x + A2*dx;
        v4sf f1y = A0*-dvy - A1*f0y + A2*dy;
	v4sf f1z = A0*-dvz - A1*f0z + A2*dz;

        x_ai += f0x;
        y_ai += f0y;
        z_ai += f0z;

        x_ji += f1x;
        y_ji += f1y;
        z_ji += f1z;

    }
    i_prt->acc_short.set(sum(x_ai), sum(y_ai), sum(z_ai));
    i_prt->jrk_short.set(sum(x_ji), sum(y_ji), sum(z_ji));
}

void Hard_System::calc_acc_jrk_duncan4_from_jarray_sse_fulld(Particle_Short* i_prt,
							     const int& j_head,
							     const int& j_tale,
							     const double& Tnext,
							     const int& mode){
    static v2df ZERO = (v2df){0.0, 0.0};
    static v2df ONE = (v2df){1.0, 1.0};
    static v2df THREE = (v2df){3.0, 3.0};
    static v2df v2df_20 = (v2df){20.0, 20.0};
    static v2df v2df_70 = (v2df){70.0, 70.0};
    static v2df v2df_84 = (v2df){84.0, 84.0};
    static v2df v2df_35 = (v2df){35.0, 35.0};
    static v2df v2df_140 = (v2df){140.0, 140.0};
    static v2df v2df_420 = (v2df){420.0, 420.0};
    static v2df LARGE = (v2df){999999999.9, 999999999.9};
    //static int v4si v4si_NBH_ORG_GLB = (v4si){NBH_ORG_GLB, NBH_ORG_GLB, NBH_ORG_GLB, NBH_ORG_GLB};
    static v2df v2df_NBH_GLB_ORG = (v2df){(double)NBH_GLB_ORG*0.9999999999, (double)NBH_GLB_ORG*0.9999999999};

    v2df eps_sq_FS_FS = (v2df){eps2_FS_FS, eps2_FS_FS};
    v2df eps_sq_FS_BH = (v2df){eps2_FS_BH, eps2_FS_BH};
    v2df eps_sq_BH_BH = (v2df){eps2_BH_BH, eps2_BH_BH};
    v2df rout_FS_FS = (v2df){rcut_out_FS_FS, rcut_out_FS_FS};
    v2df rout_FS_BH = (v2df){rcut_out_FS_BH, rcut_out_FS_BH};
    v2df rout_BH_BH = (v2df){rcut_out_BH_BH, rcut_out_BH_BH};
    v2df rin_FS_FS = (v2df){rcut_in_FS_FS, rcut_in_FS_FS};
    v2df rin_FS_BH = (v2df){rcut_in_FS_BH, rcut_in_FS_BH};
    v2df rin_BH_BH = (v2df){rcut_in_BH_BH, rcut_in_BH_BH};

    v2df eps_sq, rout, rin;
    v2df tnext = (v2df){Tnext, Tnext};

    v2df x_pi, y_pi, z_pi, x_vi, y_vi, z_vi;
    if(mode == 0){
	double pvi = i_prt->pos[0];
	x_pi = (v2df){pvi, pvi};
	pvi = i_prt->pos[1];
	y_pi = (v2df){pvi, pvi};
	pvi = i_prt->pos[2];
	z_pi = (v2df){pvi, pvi};
	pvi = i_prt->vel[0];
	x_vi = (v2df){pvi, pvi};
	pvi = i_prt->vel[1];
	y_vi = (v2df){pvi, pvi};
	pvi = i_prt->vel[2];
	z_vi = (v2df){pvi, pvi};
    }
    else if(mode == 1){
	double pvi = i_prt->pos_pre[0];
	x_pi = (v2df){pvi, pvi};
	pvi = i_prt->pos_pre[1];
	y_pi = (v2df){pvi, pvi};
	pvi = i_prt->pos_pre[2];
	z_pi = (v2df){pvi, pvi};
	pvi = i_prt->vel_pre[0];
	x_vi = (v2df){pvi, pvi};
	pvi = i_prt->vel_pre[1];
	y_vi = (v2df){pvi, pvi};
	pvi = i_prt->vel_pre[2];
	z_vi = (v2df){pvi, pvi};
    }

    v2df x_ai, y_ai, z_ai, x_ji, y_ji, z_ji, pot;
    x_ai = y_ai = z_ai = x_ji = y_ji = z_ji = pot = ZERO;

    v2df idx_pj;
    v2df x_pj, y_pj, z_pj, x_vj, y_vj, z_vj, m_j;
    Neighbour_List* j_first = i_prt->ngh_list_first;
    for(int j = j_head; j < j_tale; j += 2){
	int nj_tmp = j+2 < j_tale ? 2 : j_tale - j; 
	// set pos and vel of j paritcles
	if(nj_tmp == 2){
	    m_j = (v2df){(j_first+j)->prtj->mass, (j_first+j+1)->prtj->mass};
	    idx_pj = (v2df){(j_first+j)->prtj->index, (j_first+j+1)->prtj->index};
	}
	else if(nj_tmp == 1){
	    m_j = (v2df){(j_first+j)->prtj->mass, 0.0};
	    idx_pj = (v2df){(j_first+j)->prtj->index, -1.0};
	}

	if(mode == 0){
	    if(nj_tmp == 2){
		x_pj = (v2df){(j_first+j)->prtj->pos[0], (j_first+j+1)->prtj->pos[0]};
		x_vj = (v2df){(j_first+j)->prtj->vel[0], (j_first+j+1)->prtj->vel[0]};

		y_pj = (v2df){(j_first+j)->prtj->pos[1], (j_first+j+1)->prtj->pos[1]}; 
		y_vj = (v2df){(j_first+j)->prtj->vel[1], (j_first+j+1)->prtj->vel[1]};

		z_pj = (v2df){(j_first+j)->prtj->pos[2], (j_first+j+1)->prtj->pos[2]}; 
		z_vj = (v2df){(j_first+j)->prtj->vel[2], (j_first+j+1)->prtj->vel[2]};
	    }
	    else if(nj_tmp == 1){
		// the second element must be larger than rcut.
		x_pj = (v2df){(j_first+j)->prtj->pos[0], LARGEDOUBLE};
		x_vj = (v2df){(j_first+j)->prtj->vel[0], LARGEDOUBLE};

		y_pj = (v2df){(j_first+j)->prtj->pos[1], LARGEDOUBLE};
		y_vj = (v2df){(j_first+j)->prtj->vel[1], LARGEDOUBLE};

		z_pj = (v2df){(j_first+j)->prtj->pos[2], LARGEDOUBLE};
		z_vj = (v2df){(j_first+j)->prtj->vel[2], LARGEDOUBLE};
	    }
	}
	else if(mode == 1){
	    v2df t_pre;
	    if(nj_tmp == 2){
		t_pre = (v2df){(j_first+j)->prtj->time_pre, (j_first+j+1)->prtj->time_pre};
	    }
	    else if(nj_tmp == 1){
		t_pre = (v2df){(j_first+j)->prtj->time_pre, Tnext};
	    }
	    v2df mask = (v2df)__builtin_ia32_cmpneqpd(t_pre, tnext); // t_pre != tnext ? 0xffffffff : 0
	    int bits = __builtin_ia32_movmskpd(mask);
	    if(bits){
		if(bits & 1){
		    double dt = Tnext - (j_first+j)->prtj->time;
		    (j_first+j)->prtj->predict_short_h4(dt);
		    (j_first+j)->prtj->time_pre = Tnext;
		}
		if(bits & 2){
		    double dt = Tnext - (j_first+j+1)->prtj->time;
		    (j_first+j+1)->prtj->predict_short_h4(dt);
		    (j_first+j+1)->prtj->time_pre = Tnext;
		}
	    }
	    if(nj_tmp == 2){
		x_pj = (v2df){(j_first+j)->prtj->pos_pre[0], (j_first+j+1)->prtj->pos_pre[0]}; 
		x_vj = (v2df){(j_first+j)->prtj->vel_pre[0], (j_first+j+1)->prtj->vel_pre[0]};

		y_pj = (v2df){(j_first+j)->prtj->pos_pre[1], (j_first+j+1)->prtj->pos_pre[1]}; 
		y_vj = (v2df){(j_first+j)->prtj->vel_pre[1], (j_first+j+1)->prtj->vel_pre[1]};

		z_pj = (v2df){(j_first+j)->prtj->pos_pre[2], (j_first+j+1)->prtj->pos_pre[2]}; 
		z_vj = (v2df){(j_first+j)->prtj->vel_pre[2], (j_first+j+1)->prtj->vel_pre[2]};
	    }
	    else if(nj_tmp == 1){
		x_pj = (v2df){(j_first+j)->prtj->pos_pre[0], LARGEDOUBLE};
		x_vj = (v2df){(j_first+j)->prtj->vel_pre[0], LARGEDOUBLE};

		y_pj = (v2df){(j_first+j)->prtj->pos_pre[1], LARGEDOUBLE};
		y_vj = (v2df){(j_first+j)->prtj->vel_pre[1], LARGEDOUBLE};

		z_pj = (v2df){(j_first+j)->prtj->pos_pre[2], LARGEDOUBLE};
		z_vj = (v2df){(j_first+j)->prtj->vel_pre[2], LARGEDOUBLE};
	    }
	}
	v2df mask = (v2df)__builtin_ia32_cmpltpd(idx_pj, v2df_NBH_GLB_ORG); // idx_pj < NBH_GLB_ORG ? 0xffffffff : 0
	if(i_prt->index < NBH_GLB_ORG){
	    eps_sq = __builtin_ia32_andpd(mask, eps_sq_BH_BH) + __builtin_ia32_andnpd(mask, eps_sq_FS_BH);
	    rout = __builtin_ia32_andpd(mask, rout_BH_BH) + __builtin_ia32_andnpd(mask, rout_FS_BH);
	    rin = __builtin_ia32_andpd(mask, rin_BH_BH) + __builtin_ia32_andnpd(mask, rin_FS_BH);
	}
	else{
	    eps_sq = __builtin_ia32_andpd(mask, eps_sq_FS_BH) + __builtin_ia32_andnpd(mask, eps_sq_FS_FS);
	    rout = __builtin_ia32_andpd(mask, rout_FS_BH) + __builtin_ia32_andnpd(mask, rout_FS_FS);
	    rin = __builtin_ia32_andpd(mask, rin_FS_BH) + __builtin_ia32_andnpd(mask, rin_FS_FS);
	}
	v2df dx = x_pi - x_pj;
        v2df dy = y_pi - y_pj;
	v2df dz = z_pi - z_pj;

        v2df r2 = dx*dx + dy*dy + dz*dz;
	v2df rout_sq = rout*rout;
	v2df mask_df = (v2df)__builtin_ia32_cmpltpd(r2, rout_sq); // r2 < rout_sq ? 0xffffffff : 0
	int bit = __builtin_ia32_movmskpd(mask_df);
	if(!bit){ continue; } // all j particles are outside the neighbour sphere

	v2df r_wo_eps = __builtin_ia32_sqrtpd(r2);
	v2df r_wo_eps_inv = inv(r_wo_eps);

	r2 += eps_sq;
        v2df r2_inv = inv(r2);
	v2df r_inv = __builtin_ia32_sqrtpd(r2_inv);
        v2df r3_inv = r2_inv*r_inv;

	v2df dvx = x_vi - x_vj;
	v2df dvy = y_vi - y_vj;
	v2df dvz = z_vi - z_vj;

        v2df rv = dx*dvx + dy*dvy + dz*dvz;

	// evaluate change over function by Duncan for 4th-order
	v2df dr_cut_inv = inv(rout-rin);
        v2df rtmp = (r_wo_eps-rin)*dr_cut_inv;
        v2df rtmpdot = rv*r_wo_eps_inv*dr_cut_inv;

        rtmp = __builtin_ia32_andpd(rtmp, (v2df)__builtin_ia32_cmpltpd(ZERO, rtmp)); // rtmp = (0.0 < rtmp ?  rtmp : 0.0)
	rtmp = __builtin_ia32_andpd((rtmp - ONE), (v2df)__builtin_ia32_cmpltpd(rtmp, ONE)) + ONE; // rtmp = (rtmp < 1.0 ?  rtmp : 1.0)
	rtmpdot = __builtin_ia32_andpd(rtmpdot, (v2df)__builtin_ia32_cmpltpd(ZERO, rtmp));  // rtmpdot = (0.0 < rtmp ?  rtmpdot : 0.0)
	rtmpdot = __builtin_ia32_andpd(rtmpdot, (v2df)__builtin_ia32_cmpltpd(rtmp, ONE)); // rtmpdot = (rtmp < 1.0 ?  rtmpdot : 0.0)
        v2df rtmp2 = rtmp*rtmp;
        v2df rtmp3 = rtmp2*rtmp;
	v2df rtmp4 = rtmp2*rtmp2;
        v2df K = ( ( (-v2df_20*rtmp + v2df_70)*rtmp - v2df_84)*rtmp + v2df_35)*rtmp4;
        v2df Kdot = ( ((( -v2df_140*rtmp + v2df_420)*rtmp - v2df_420)*rtmp + v2df_140)*rtmp3) * rtmpdot;

        v2df A0 = m_j * r3_inv * (ONE-K);
        v2df f0x = A0 * -dx;
        v2df f0y = A0 * -dy;
        v2df f0z = A0 * -dz;

        v2df A1 = THREE*rv*r2_inv;
	v2df A2 = m_j * r3_inv * Kdot;
        v2df f1x = A0*-dvx - A1*f0x + A2*dx;
        v2df f1y = A0*-dvy - A1*f0y + A2*dy;
	v2df f1z = A0*-dvz - A1*f0z + A2*dz;

        x_ai += f0x;
        y_ai += f0y;
        z_ai += f0z;

        x_ji += f1x;
        y_ji += f1y;
        z_ji += f1z;

    }
    i_prt->acc_short.set(sum(x_ai), sum(y_ai), sum(z_ai));
    i_prt->jrk_short.set(sum(x_ji), sum(y_ji), sum(z_ji));
}


void Hard_System::calc_acc_jrk_duncan4_from_jarray_scholar(Particle_Short* i_prt,
							   const int& j_head,
							   const int& j_tale,
							   const double& Tnext,
							   const int& mode){

    double rout, rin, eps_sq; 
    double x_pi, y_pi, z_pi;
    double x_vi, y_vi, z_vi;
    if(mode == 0){
	double pvi = i_prt->pos[0];
	x_pi = pvi;
	pvi = i_prt->pos[1];
	y_pi = pvi;
	pvi = i_prt->pos[2];
	z_pi = pvi;

	pvi = i_prt->vel[0];
	x_vi = pvi;
	pvi = i_prt->vel[1];
	y_vi = pvi;
	pvi = i_prt->vel[2];
	z_vi = pvi;
    }
    else if(mode == 1){
	double pvi = i_prt->pos_pre[0];
	x_pi = pvi;
	pvi = i_prt->pos_pre[1];
	y_pi = pvi;
	pvi = i_prt->pos_pre[2];
	z_pi = pvi;

	pvi = i_prt->vel_pre[0];
	x_vi = pvi;
	pvi = i_prt->vel_pre[1];
	y_vi = pvi;
	pvi = i_prt->vel_pre[2];
	z_vi = pvi;
    }

    double x_ai, y_ai, z_ai, x_ji, y_ji, z_ji, pot;
    x_ai = y_ai = z_ai = x_ji = y_ji = z_ji = pot = 0.0;

    double x_pj, y_pj, z_pj, x_vj, y_vj, z_vj, m_j;
    double r2min = 99999999.9;
    int adrmin = -1;
    Neighbour_List* j_first = i_prt->ngh_list_first;
    for(int j = j_head; j < j_tale; j++){
	// set pos and vel of j paritcles
	m_j = (j_first+j)->prtj->mass; // modify
	if(mode == 0){
	    x_pj = (j_first+j)->prtj->pos[0];
	    x_vj = (j_first+j)->prtj->vel[0];
	    
	    y_pj = (j_first+j)->prtj->pos[1];
	    y_vj = (j_first+j)->prtj->vel[1];

	    z_pj = (j_first+j)->prtj->pos[2];
	    z_vj = (j_first+j)->prtj->vel[2];
	}
	else if(mode == 1){
	    double t_pre;
	    t_pre = (j_first+j)->prtj->time_pre;
	    if(t_pre != Tnext){
		double dt = Tnext - (j_first+j)->prtj->time;
		(j_first+j)->prtj->predict_short_h4(dt);
		(j_first+j)->prtj->time_pre = Tnext;
	    }
	    x_pj = (j_first+j)->prtj->pos_pre[0];
	    x_vj = (j_first+j)->prtj->vel_pre[0];

	    y_pj = (j_first+j)->prtj->pos_pre[1];
	    y_vj = (j_first+j)->prtj->vel_pre[1];

	    z_pj = (j_first+j)->prtj->pos_pre[2];
	    z_vj = (j_first+j)->prtj->vel_pre[2];

	}

	if(i_prt->index < NBH_GLB_ORG && (j_first+j)->prtj->index < NBH_GLB_ORG){
	    eps_sq = eps2_BH_BH;
	    rout = rcut_out_BH_BH;
	    rin = rcut_in_BH_BH;
	}
	else if(i_prt->index >= NBH_GLB_ORG && (j_first+j)->prtj->index >= NBH_GLB_ORG){
	    eps_sq = eps2_FS_FS;
	    rout = rcut_out_FS_FS;
	    rin = rcut_in_FS_FS;
	}
	else{
	    eps_sq = eps2_FS_BH;
	    rout = rcut_out_FS_BH;
	    rin = rcut_in_FS_BH;
	}

	double dx = x_pi - x_pj;
        double dy = y_pi - y_pj;
	double dz = z_pi - z_pj;
	double rout2 = rout*rout;
	double rin =  rcut_in_FS_FS;
	double dr_cut = rout - rin;
	double dr_cut_inv = 1.0/dr_cut;

        double r2 = dx*dx + dy*dy + dz*dz;

	if(rout2 < r2){continue;}

	double r_wo_eps = sqrt(r2);
	double r_wo_eps_inv = 1.0/r_wo_eps;

	r2 += eps_sq;

	double dvx = x_vi - x_vj;
	double dvy = y_vi - y_vj;
	double dvz = z_vi - z_vj;

	double r2_inv = 1.0/r2;
	double r_inv = sqrt(r2_inv);
        double r3_inv = r2_inv*r_inv;

        double rv = dx*dvx + dy*dvy + dz*dvz;

	// evaluate change over function by Duncan for 4th-order
        double rtmp = (r_wo_eps - rin)*dr_cut_inv;
        double rtmpdot = rv*r_wo_eps_inv*dr_cut_inv;
	double K, Kdot;
	if(rtmp < 0.0){
	    K = 0.0;
	    Kdot = 0.0;
	}
	else if(rtmp > 1.0){
	    K = 1.0;
	    Kdot = 0.0;
	}
	else{
	    double rtmp2 = rtmp*rtmp;
	    double rtmp3 = rtmp2*rtmp;
	    double rtmp4 = rtmp2*rtmp2;
	    K = ((( -20.0*rtmp + 70.0)*rtmp - 84.0)*rtmp + 35.0)*rtmp4;
	    Kdot = ( ((( -140.0*rtmp + 420.0)*rtmp - 420.0)*rtmp + 140.0)*rtmp3) * rtmpdot;
	}
        double A0 = m_j * r3_inv * (1.0-K);
        double f0x = A0 * -dx;
        double f0y = A0 * -dy;
        double f0z = A0 * -dz;

        double A1 = 3.0*rv*r2_inv;
        double A2 = m_j * r3_inv * Kdot;
        double f1x = A0*-dvx - A1*f0x + A2*dx;
        double f1y = A0*-dvy - A1*f0y + A2*dy;
	double f1z = A0*-dvz - A1*f0z + A2*dz;

        x_ai += f0x;
        y_ai += f0y;
        z_ai += f0z;

        x_ji += f1x;
        y_ji += f1y;
        z_ji += f1z;
    }
    i_prt->acc_short.set(x_ai, y_ai, z_ai);
    i_prt->jrk_short.set(x_ji, y_ji, z_ji);
}




#endif
