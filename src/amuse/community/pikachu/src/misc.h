#ifndef __MISC_H__
#define __MISC_H__

#include<fstream>

void dump_calc_time(ostream& fout,
		    const double& Tsys,
		    const int& Nshort_glb,
		    const int& ngh_list_len_glb, 
		    const double dt_glb){

    const double inv_dt_glb = 1.0/dt_glb;
    fout<<endl;
    fout<<"------------------------------"<<endl;
    fout<<"Tsys(after drift)="<<Tsys<<endl;
    fout<<"TCAL_LOOP="<<TCAL_LOOP*inv_dt_glb<<endl;
    fout<<endl;
    fout<<"NFS_LOC+NBH_LOC(after drift)="<<NFS_LOC+NBH_LOC<<endl;
    fout<<"NALL_LOC(after drift)="<<NALL_LOC<<endl;
    fout<<"--- HERMITE INFO ---"<<endl;
    fout<<"Nshort_glb(after drift)="<<Nshort_glb<<endl;
    fout<<"ngh_list_len_glb(after drift)="<<ngh_list_len_glb<<endl;
    fout<<"NGH_LIST_LEN_2BODY+NGH_LIST_LEN_MULBODY="<<NGH_LIST_LEN_2BODY+NGH_LIST_LEN_MULBODY<<endl;
    fout<<"NGH_LIST_LEN_2BODY="<<NGH_LIST_LEN_2BODY<<endl;
    fout<<"NGH_LIST_LEN_MULBODY="<<NGH_LIST_LEN_MULBODY<<endl;
    //fout<<"STEPS_HARD(par body)="<<STEPS_HARD*inv_dt_glb<<endl;
    fout<<"STEPS_HARD(par body)="<<STEPS_HARD<<endl;
    fout<<"--- TREE INFO ---"<<endl;
    fout<<"LIST_LEN_TOT="<<LIST_LEN_TOT<<endl;
    fout<<"NPRT_TOT="<<NPRT_TOT<<endl;
    fout<<"LOOP_WALK="<<LOOP_WALK<<endl;
    if(LOOP_WALK != 0){
	fout<<"Nlist="<<LIST_LEN_TOT/LOOP_WALK<<endl;
	fout<<"Ng="<<NPRT_TOT/LOOP_WALK<<endl;
    }
    else{
	fout<<"Nlist=0"<<endl;
	fout<<"Ng=0"<<endl;
    }
    fout<<endl;
    fout<<"TCAL_EVOLVE_HARD="<<TCAL_EVOLVE_HARD*inv_dt_glb<<endl;
    fout<<"    TCAL_DRIFT="<<TCAL_DRIFT*inv_dt_glb<<endl;
    fout<<"    TCAL_HERMITE4="<<TCAL_HERMITE4*inv_dt_glb<<endl;
    fout<<"        TCAL_HERMITE4_2BODY="<<TCAL_HERMITE4_2BODY*inv_dt_glb<<endl;
    fout<<"        TCAL_HERMITE4_MULBODY="<<TCAL_HERMITE4_MULBODY*inv_dt_glb<<endl;
    fout<<"        TCAL_HERMITE4_PRE="<<TCAL_HERMITE4_PRE*inv_dt_glb<<endl;
    fout<<"        TCAL_HERMITE4_COMM="<<TCAL_HERMITE4_COMM*inv_dt_glb<<endl;
    fout<<"    TCAL_COPY_PRT_SHORT="<<TCAL_COPY_PRT_SHORT*inv_dt_glb<<endl;
    fout<<endl;
    fout<<"TCAL_DIVIDE_PRT="<<TCAL_DIVIDE_PRT*inv_dt_glb<<endl;
    fout<<endl;
    fout<<"TCAL_EVOLVE_SOFT="<<(TCAL_SOFT_FORCES+TCAL_KICK)*inv_dt_glb<<endl;
    fout<<"    TCAL_SOFT_FORCES="<<TCAL_SOFT_FORCES*inv_dt_glb<<endl;
    fout<<"        TCAL_TREE_SETUP="<<TCAL_TREE_SETUP*inv_dt_glb<<endl;
    fout<<"            TCAL_FIRST_HALF="<<TCAL_FIRST_HALF*inv_dt_glb<<endl;
    fout<<"            TCAL_TREE_SETUP_LOC="<<TCAL_TREE_SETUP_LOC*inv_dt_glb<<endl;
    fout<<"            TCAL_TREE_INSERT_PRT_LOC="<<TCAL_TREE_INSERT_PRT_LOC*inv_dt_glb<<endl;
    fout<<"            TCAL_TREE_SET_CM_LOC="<<TCAL_TREE_SET_CM_LOC*inv_dt_glb<<endl;
    fout<<"            TCAL_TREE_LET_EX="<<TCAL_TREE_LET_EX*inv_dt_glb<<endl;
    fout<<"            TCAL_TREE_SETUP_LET="<<TCAL_TREE_SETUP_LET*inv_dt_glb<<endl;
    fout<<"        TCAL_TREE_EVALUATE="<<TCAL_TREE_EVALUATE*inv_dt_glb<<endl;
    fout<<"            TCAL_TREE_WALK="<<TCAL_TREE_WALK*inv_dt_glb<<endl;
    fout<<"            TCAL_PACK_PRT="<<TCAL_PACK_PRT*inv_dt_glb<<endl;
    fout<<"            TCAL_TREE_FORCE="<<TCAL_TREE_FORCE*inv_dt_glb<<endl;
    fout<<"                TCAL_DEV_COMM="<<TCAL_DEV_COMM*inv_dt_glb<<endl;
    fout<<"        TCAL_DIRECT_EVALUATE="<<TCAL_DIRECT_EVALUATE*inv_dt_glb<<endl;
    fout<<"        TCAL_NGH_SEARCH="<<TCAL_NGH_SEARCH*inv_dt_glb<<endl;
    fout<<"            TCAL_NGH_SEARCH_MAKE_DIC="<<TCAL_NGH_SEARCH_MAKE_DIC*inv_dt_glb<<endl;
    fout<<"        TCAL_ACC_ONLY="<<TCAL_ACC_ONLY*inv_dt_glb<<endl;
    fout<<"    TCAL_KICK="<<TCAL_KICK*inv_dt_glb<<endl;
    fout<<endl;
    fout<<"TCAL_MERGE_LIST="<<TCAL_MERGE_LIST*inv_dt_glb<<endl;
    fout<<"------------------------------"<<endl;
    fout<<endl;
    fout<<endl;
}

#endif
