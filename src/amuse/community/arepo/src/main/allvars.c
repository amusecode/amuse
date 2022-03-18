/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/main/allvars.c
 * \date        05/2018
 * \brief       Contains all global variables.
 * \details     This file contains the global variables used in Arepo.
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 21.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include "../main/allvars.h"

struct data_nodelist *DataNodeList; /* to be deleted */

MyDouble boxSize, boxHalf;

#ifdef LONG_X
MyDouble boxSize_X, boxHalf_X;
#else  /* #ifdef LONG_X */
#endif /* #ifdef LONG_X #else */
#ifdef LONG_Y
MyDouble boxSize_Y, boxHalf_Y;
#else  /* #ifdef LONG_Y */
#endif /* #ifdef LONG_Y #else */
#ifdef LONG_Z
MyDouble boxSize_Z, boxHalf_Z;
#else  /* #ifdef LONG_Z */
#endif /* #ifdef LONG_Z #else */

#ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
MPI_Status mpistat;
#endif /* #ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG */

/*********************************************************/
/*  Global variables                                     */
/*********************************************************/

int ThisTask; /*!< the number of the local processor  */
int NTask;    /*!< number of processors */
int PTask;    /*!< note: NTask = 2^PTask */

int ThisNode;        /*!< the rank of the current compute node  */
int NumNodes;        /*!< the number of compute nodes used  */
int MinTasksPerNode; /*!< the minimum number of MPI tasks that is found on any of the nodes  */
int MaxTasksPerNode; /*!< the maximum number of MPI tasks that is found on any of the nodes  */
int TasksInThisNode; /*!< number of MPI tasks on  current compute node */
int RankInThisNode;  /*!< rank of the MPI task on the current compute node */
long long MemoryOnNode;
double CPUThisRun; /*!< Sums CPU time of current process */
int MaxTopNodes;   /*!< Maximum number of nodes in the top-level tree used for domain decomposition */
int RestartFlag;   /*!< taken from command line used to start code. 0 is normal start-up from
                      initial conditions, 1 is resuming a run from a set of restart files, while 2
                      marks a restart from a snapshot file. */
int RestartSnapNum;
int Argc;
char **Argv;

size_t AllocatedBytes;
size_t FreeBytes;

int Nforces;
int *TargetList;
struct thread_data Thread[NUM_THREADS];

#ifdef IMPOSE_PINNING
hwloc_cpuset_t cpuset_thread[NUM_THREADS];
#endif /* #ifdef IMPOSE_PINNING */

int *Exportflag,
    *ThreadsExportflag[NUM_THREADS]; /*!< Buffer used for flagging whether a particle needs to be exported to another process */
int *Exportnodecount;
int *Exportindex;

int *Send_offset, *Send_count, *Recv_count, *Recv_offset;
int *Send_offset_nodes, *Send_count_nodes, *Recv_count_nodes, *Recv_offset_nodes;
int *TasksThatSend, *TasksThatRecv, NSendTasks, NRecvTasks;
struct send_recv_counts *Send, *Recv;

int Mesh_nimport, Mesh_nexport, *Mesh_Send_offset, *Mesh_Send_count, *Mesh_Recv_count, *Mesh_Recv_offset;
int Force_nimport, Force_nexport, *Force_Send_offset, *Force_Send_count, *Force_Recv_count, *Force_Recv_offset;

int TakeLevel;
int TagOffset;

int TimeBinSynchronized[TIMEBINS];
struct TimeBinData TimeBinsHydro, TimeBinsGravity;

#ifdef USE_SFR
double TimeBinSfr[TIMEBINS];
#endif

#ifdef SUBFIND
int GrNr;
int NumPartGroup;
#endif /* #ifdef SUBFIND */

char DumpFlag         = 1;
char DumpFlagNextSnap = 1;

int FlagNyt = 0;

double CPU_Step[CPU_LAST];
double CPU_Step_Stored[CPU_LAST];

double WallclockTime; /*!< This holds the last wallclock time measurement for timings measurements */
double StartOfRun;    /*!< This stores the time of the start of the run for evaluating the elapsed time */

double EgyInjection;

int NumPart; /*!< number of particles on the LOCAL processor */
int NumGas;  /*!< number of gas particles on the LOCAL processor  */

gsl_rng *random_generator;     /*!< a random number generator  */
gsl_rng *random_generator_aux; /*!< an auxialiary random number generator for use if one doesn't want to influence the main code's
                                  random numbers  */

#ifdef USE_SFR
int Stars_converted; /*!< current number of star particles in gas particle block */
#endif

#ifdef TOLERATE_WRITE_ERROR
int WriteErrorFlag;
char AlternativeOutputDir[MAXLEN_PATH];
#endif /* #ifdef TOLERATE_WRITE_ERROR */

double TimeOfLastDomainConstruction; /*!< holds what it says */

int *Ngblist; /*!< Buffer to hold indices of neighbours retrieved by the neighbour search
                 routines */

double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
double DomainInverseLen, DomainBigFac;
int *DomainStartList, *DomainEndList;
double *DomainCost, *TaskCost;
int *DomainCount, *TaskCount;
struct no_list_data *ListNoData;

int domain_bintolevel[TIMEBINS];
int domain_refbin[TIMEBINS];
int domain_grav_weight[TIMEBINS];
int domain_hydro_weight[TIMEBINS];
int domain_to_be_balanced[TIMEBINS];

int *DomainTask;
int *DomainNewTask;
int *DomainNodeIndex;

peanokey *Key, *KeySorted;

struct topnode_data *TopNodes;

int NTopnodes, NTopleaves;

/* variables for input/output , usually only used on process 0 */

char ParameterFile[MAXLEN_PATH]; /*!< file name of parameterfile used for starting the simulation */

FILE *FdInfo,   /*!< file handle for info.txt log-file. */
    *FdEnergy,  /*!< file handle for energy.txt log-file. */
    *FdTimings, /*!< file handle for timings.txt log-file. */
    *FdDomain,  /*!< file handle for domain.txt log-file. */
    *FdBalance, /*!< file handle for balance.txt log-file. */
    *FdMemory,  /*!< file handle for memory.txt log-file. */
    *FdTimebin, /*!< file handle for timebins.txt log-file. */
    *FdCPU;     /*!< file handle for cpu.txt log-file. */

#ifdef DETAILEDTIMINGS
FILE *FdDetailed;
#endif /* #ifdef DETAILEDTIMINGS */

#ifdef OUTPUT_CPU_CSV
FILE *FdCPUCSV;
#endif /* #ifdef OUTPUT_CPU_CSV */

#ifdef RESTART_DEBUG
FILE *FdRestartTest;
#endif /* #ifdef RESTART_DEBUG */

#ifdef USE_SFR
FILE *FdSfr; /*!< file handle for sfr.txt log-file. */
#endif

struct pair_data *Pairlist;

#ifdef FORCETEST
FILE *FdForceTest; /*!< file handle for forcetest.txt log-file. */
#endif             /* #ifdef FORCETEST */

int WriteMiscFiles = 1;

void *CommBuffer; /*!< points to communication buffer, which is used at a few places */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
struct global_data_all_processes All;

/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
struct particle_data *P, /*!< holds particle data on local processor */
    *DomainPartBuf;      /*!< buffer for particle data used in domain decomposition */

struct subfind_data *PS;

/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
struct sph_particle_data *SphP, /*!< holds SPH particle data on local processor */
    *DomainSphBuf;              /*!< buffer for SPH particle data in domain decomposition */

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
struct special_particle_data *PartSpecialListGlobal;
#endif /* #ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE */

peanokey *DomainKeyBuf;

/*! global state of system
 */
struct state_of_system SysState, SysStateAtStart, SysStateAtEnd;

/*! Various structures for communication during the gravity computation.
 */
struct directdata *DirectDataIn, *DirectDataAll;
struct accdata *DirectAccOut, *DirectAccIn;
int ThreadsNexport[NUM_THREADS], ThreadsNexportNodes[NUM_THREADS];
struct data_partlist *PartList, *ThreadsPartList[NUM_THREADS];
struct datanodelist *NodeList, *ThreadsNodeList[NUM_THREADS];
struct potdata_out *PotDataResult, /*!< holds the partial results computed for imported particles. Note: We use GravDataResult =
                                      GravDataGet, such that the result replaces the imported data */
    *PotDataOut; /*!< holds partial results received from other processors. This will overwrite the GravDataIn array */

/*! Header for the standard file format.
 */
struct io_header header; /*!< holds header for snapshot files */
#ifdef NTYPES_ICS
struct io_header_ICs header_ICs; /*!< holds header for IC files */
#endif                           /* #ifdef NTYPES_ICS */
char (*Parameters)[MAXLEN_PARAM_TAG];
char (*ParametersValue)[MAXLEN_PARAM_VALUE];
char *ParametersType;

/*! Variables for gravitational tree
 * ------------------
 */
int Tree_MaxPart;
int Tree_NumNodes;
int Tree_MaxNodes;
int Tree_FirstNonTopLevelNode;
int Tree_NumPartImported;
int Tree_NumPartExported;
int Tree_ImportedNodeOffset;
int Tree_NextFreeNode;
MyDouble *Tree_Pos_list;
unsigned long long *Tree_IntPos_list;
int *Tree_Task_list;
int *Tree_ResultIndexList;

struct treepoint_data *Tree_Points;
struct resultsactiveimported_data *Tree_ResultsActiveImported;

int *Nextnode; /*!< gives next node in tree walk  (nodes array) */
int *Father;   /*!< gives parent node in tree (Prenodes array) */

struct NODE *Nodes; /*!< points to the actual memory allocted for the nodes */
                    /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart]
                       gives the first allocated node */

#ifdef MULTIPLE_NODE_SOFTENING
struct ExtNODE *ExtNodes;
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

float *Nodes_GravCost;

/*! Variables for neighbor tree
 * -----------------
 */
int Ngb_MaxPart;
int Ngb_NumNodes;
int Ngb_MaxNodes;
int Ngb_FirstNonTopLevelNode;
int Ngb_NextFreeNode;
int *Ngb_Father;
int *Ngb_Marker;
int Ngb_MarkerValue;

int *Ngb_DomainNodeIndex;
int *DomainListOfLocalTopleaves;
int *DomainNLocalTopleave;
int *DomainFirstLocTopleave;
int *Ngb_Nextnode;

/*! The ngb-tree data structure
 */
struct NgbNODE *Ngb_Nodes;
struct ExtNgbNODE *ExtNgb_Nodes;

#ifdef STATICNFW
double Rs, R200;
double Dc;
double RhoCrit, V200;
double fac;
#endif /* #ifdef STATICNFW */

int MaxThreads = 1;

IO_Field *IO_Fields;
int N_IO_Fields   = 0;
int Max_IO_Fields = 0;
