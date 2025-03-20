#ifndef CONST_H
#define CONST_H

#include<sys/time.h>

/////////// astrophisical constant (MKS) //////////////
const double Gsi = 6.67259E-11; // m^3/kg/sec^2
const double csi = 2.99792458E08;  // m/sec
const double Pc = 3.08567802E16;   //1pc
const double AUsi = 1.49597892E11;  // 1AU
const double Rsolsi = 6.9599E08;  //1R_sol
const double Msolsi = 1.989E30; //1M_sol
const double Yearsi = 3.1558149984E07;  //1year

const double G = 1.0;

const double n2 = 0.5;
const double n3 = 1.0/3.0;
const double n4 = 1.0/4.0;
const double n5 = 1.0/5.0;
const double n6 = 1.0/6.0;
const double n7 = 1.0/7.0;
const double n8 = 1.0/8.0;
const double n9 = 1.0/9.0;
const double n10 = 1.0/10.0;
const double n12 = 1.0/12.0;
const double n24 = 1.0/24.0;
const double n30 = 1.0/30.0;
const double n32 = 1.0/32.0;
const double n60 = 1.0/60.0;
const double n40 = 1.0/40.0;
const double n42 = 1.0/42.0;
const double n84 = 1.0/84.0;
const double n120 = 1.0/120.0;
const double n280 = 1.0/280.0;
const double n720 = 1.0/720.0;
const double n1024 = 1.0/1024.0;
const double n1680 = 1.0/1680.0;
const double n3_28 = 3.0/28.0;
const double n5_6 = 5.0/6.0;
const double n25_6 = 25.0/6.0;
const double n4_9 = 4.0/9.0;

const int STRINGSIZE = 1024;

const double LARGEDOUBLE = 99999999999999.9;
const float LARGESINGLE = 99999999999999.9;
const int LARGEINT = 99999999;
const int MEGAINT = 1048576;
const int NSAMPLE_MAX = 10000;
const int NPROC_MAX = 1024;

static int NJCRIT = 2000; // j paralellization factor

static double RCUT_IN_FACTOR = 0.1; // rcut_in = rcut_out * RCUT_IN_FACTOR

const int NBH_GLB_MAX = 100;

#ifdef SMALL
const int NFS_LOC_MAX = 100000;
#else
const int NFS_LOC_MAX = (int)(8.4e6);
#endif

const int NPRT_LOC_MAX = NFS_LOC_MAX + NBH_GLB_MAX;
const int NGH_LIST_MAX = NPRT_LOC_MAX*2.0; // maxmum neighbour list length
const int NSHORT_MAX = NGH_LIST_MAX * 0.5; //neighbour list max
const int NCELL_TREE_LOC_MAX = NPRT_LOC_MAX / 4; //  node
const int LIST_COMM_MAX = (int)(1e6); // communication list

#endif
