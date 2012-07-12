#ifdef GLOBAL_VALUE_DEFINE
#define GLOBAL
#else
#define GLOBAL extern
#endif

GLOBAL int MYRANK;
GLOBAL int NPROC;

GLOBAL int NFS_GLB;
GLOBAL int NFS_LOC;
GLOBAL int NFS_GLB_ORG;

GLOBAL int NBH_GLB;
GLOBAL int NBH_LOC;
GLOBAL int NBH_GLB_ORG;

GLOBAL int NDEAD_GLB;
GLOBAL int NALL_LOC; // NBH_LOC + NFS_LOC + NTREE

GLOBAL int LIST_LEN_TOT;
GLOBAL int NPRT_TOT;
GLOBAL int LOOP_WALK;


// for total
GLOBAL double TCAL_LOOP;

// for hard system
GLOBAL double TCAL_HERMITE4_COMM;
GLOBAL double TCAL_HERMITE4_PRE;
GLOBAL double TCAL_HERMITE4_2BODY;
GLOBAL double TCAL_HERMITE4_MULBODY;

GLOBAL int NGH_LIST_LEN_2BODY;
GLOBAL int NGH_LIST_LEN_MULBODY;


// for soft system
GLOBAL double TCAL_TREE_WALK;
GLOBAL double TCAL_PACK_PRT;
GLOBAL double TCAL_TREE_FORCE;
GLOBAL double TCAL_DEV_COMM;


// for nbody_system.h 
GLOBAL double TCAL_TREE_SETUP;
GLOBAL double TCAL_TREE_SETUP_LOC;
GLOBAL double TCAL_TREE_INSERT_PRT_LOC;
GLOBAL double TCAL_TREE_SET_CM_LOC;
GLOBAL double TCAL_TREE_LET_EX;
GLOBAL double TCAL_TREE_SETUP_LET;
GLOBAL double TCAL_TREE_EVALUATE;
GLOBAL double TCAL_DIRECT_EVALUATE;
GLOBAL double TCAL_NGH_SEARCH;
GLOBAL double TCAL_NGH_SEARCH_MAKE_DIC;
GLOBAL double TCAL_ACC_ONLY;

GLOBAL double TCAL_DRIFT;
GLOBAL double TCAL_HERMITE4;
GLOBAL double TCAL_COPY_PRT_SHORT;

GLOBAL double TCAL_EVOLVE_HARD;
GLOBAL double TCAL_DIVIDE_PRT;
GLOBAL double TCAL_SOFT_FORCES;
GLOBAL double TCAL_KICK;
GLOBAL double TCAL_MERGE_LIST;

GLOBAL double TCAL_FIRST_HALF;
GLOBAL double TCAL_LAST_HALF;

// for debug
#include<fstream>
GLOBAL ofstream fout_debug;

// for measuring the performance
GLOBAL double STEPS_HARD;


GLOBAL void dev_open();
GLOBAL void dev_close();
GLOBAL void halt_program();



#ifndef GLORBAL_H
#define GLORBAL_H
template<class T> inline int isdead(const T &prt){
    return (prt.mass <= 0.0) ? 1 : 0 ;
}


/*
inline int isdead(const Particle &prt){
    return (prt.mass <= 0.0) ? 1 : 0 ;
}
*/


#endif //GLORBAL_H
