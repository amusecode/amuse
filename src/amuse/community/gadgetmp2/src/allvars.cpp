/* ################################################################################## */
/* ###                                                                            ### */
/* ###                                 Gadgetmp2                                  ### */
/* ###                                                                            ### */
/* ###   Original: Gadget2 in the version used in Amuse                           ### */
/* ###   Author: Gadget2 and Amuse contributors                                   ### */
/* ###                                                                            ### */
/* ###   Modified: July 2020                                                      ### */
/* ###   Author: Thomas Schano                                                    ### */
/* ###                                                                            ### */
/* ###   Changes are intended to enable precise calculations in                   ### */
/* ###   non periodic small domain simulations in which comoving parts            ### */
/* ###   are simulated in std precision                                           ### */
/* ###                                                                            ### */
/* ################################################################################## */
/*! \file allvars.c
 *  \brief provides instances of all global variables.
 */

#ifndef NOMPI
#include <mpi.h>
#endif

#include "proto.hpp"


int gadgetmp2::ThisTask;		/*!< the rank of the local processor */
int gadgetmp2::NTask;               /*!< number of processors */
int gadgetmp2::PTask;	        /*!< smallest integer such that NTask <= 2^PTask */

int gadgetmp2::NumPart;		/*!< number of particles on the LOCAL processor */
int gadgetmp2::N_gas;		/*!< number of gas particles on the LOCAL processor  */
long long gadgetmp2::Ntype[6];      /*!< total number of particles of each type */
int gadgetmp2::NtypeLocal[6];       /*!< local number of particles of each type */

int gadgetmp2::NumForceUpdate;      /*!< number of active particles on local processor in current timestep  */
int gadgetmp2::NumSphUpdate;        /*!< number of active SPH particles on local processor in current timestep  */

double gadgetmp2::CPUThisRun;	/*!< Sums the CPU time for the process (current submission only) */


int gadgetmp2::RestartFlag;         /*!< taken from command line used to start code. 0 is normal start-up from
                                     initial conditions, 1 is resuming a run from a set of restart files, while 2
                                     marks a restart from a snapshot file. */

char *gadgetmp2::Exportflag=nullptr;        /*!< Buffer used for flagging whether a particle needs to be exported to another process */

int  *gadgetmp2::Ngblist;           /*!< Buffer to hold indices of neighbours retrieved by the neighbour search routines */

int gadgetmp2::TreeReconstructFlag; /*!< Signals that a new tree needs to be constructed */

int gadgetmp2::Flag_FullStep;       /*!< This flag signals that the current step involves all particles */

int gadgetmp2::ZeroTimestepEncountered;  /*!< Flag used by AMUSE. When a particle is assigned a timestep of zero, an
                                     exception is raised instead of forcing the application to exit. */


gsl_rng *gadgetmp2::random_generator; /*!< the employed random number generator of the GSL library */

my_float gadgetmp2::RndTable[RNDTABLE]; /*!< Hold a table with random numbers, refreshed every timestep */


my_float gadgetmp2::DomainCorner[3];    /*!< gives the lower left corner of simulation volume */
my_float gadgetmp2::DomainCenter[3];    /*!< gives the center of simulation volume */
my_float gadgetmp2::DomainLen;          /*!< gives the (maximum) side-length of simulation volume */
my_float gadgetmp2::DomainFac;          /*!< factor used for converting particle coordinates to a Peano-Hilbert mesh covering the simulation volume */
int    gadgetmp2::DomainMyStart;      /*!< first domain mesh cell that resides on the local processor */
int    gadgetmp2::DomainMyLast;       /*!< last domain mesh cell that resides on the local processor */
int    *gadgetmp2::DomainStartList=nullptr;   /*!< a table that lists the first domain mesh cell for all processors */
int    *gadgetmp2::DomainEndList=nullptr;     /*!< a table that lists the last domain mesh cell for all processors */
my_float *gadgetmp2::DomainWork=nullptr;        /*!< a table that gives the total "work" due to the particles stored by each processor */
int    *gadgetmp2::DomainCount=nullptr;       /*!< a table that gives the total number of particles held by each processor */
int    *gadgetmp2::DomainCountSph=nullptr;    /*!< a table that gives the total number of SPH particles held by each processor */

int    *gadgetmp2::DomainTask=nullptr;        /*!< this table gives for each leaf of the top-level tree the processor it was assigned to */
int    *gadgetmp2::DomainNodeIndex=nullptr;   /*!< this table gives for each leaf of the top-level tree the corresponding node of the gravitational tree */
my_float  *gadgetmp2::DomainTreeNodeLen=nullptr; /*!< this table gives for each leaf of the top-level tree the side-length of the corresponding node of the gravitational tree */
All_Reduce_buff  *gadgetmp2::DomainHmax=nullptr;        /*!< this table gives for each leaf of the top-level tree the maximum SPH smoothing length among the particles of the corresponding node of the gravitational tree */

DomainNODE *gadgetmp2::DomainMoment=nullptr;                    /*!< this table stores for each node of the top-level tree corresponding node data from the gravitational tree */

size_t DomainNODE::s0_off;
size_t DomainNODE::s1_off;
size_t DomainNODE::s2_off;
size_t DomainNODE::vs0_off;
size_t DomainNODE::vs1_off;
size_t DomainNODE::vs2_off;
size_t DomainNODE::mass_off;
    #ifdef UNEQUALSOFTENINGS
    #ifndef ADAPTIVE_GRAVSOFT_FORGAS
size_t DomainNODE::bitflags_off;
    #else
size_t DomainNODE::maxsoft_off;
    #endif
    #endif
size_t DomainNODE::tot_size;
mpfr_prec_t DomainNODE::prec;

peanokey *gadgetmp2::DomainKeyBuf;    /*!< this points to a buffer used during the exchange of particle data */

peanokey *gadgetmp2::Key;             /*!< a table used for storing Peano-Hilbert keys for particles */
peanokey *gadgetmp2::KeySorted;       /*!< holds a sorted table of Peano-Hilbert keys for all particles, used to construct top-level tree */


int gadgetmp2::NTopnodes;             /*!< total number of nodes in top-level tree */
int gadgetmp2::NTopleaves;            /*!< number of leaves in top-level tree. Each leaf can be assigned to a different processor */

struct topnode_data
 *gadgetmp2::TopNodes;                      /*!< points to the root node of the top-level tree */


my_float gadgetmp2::TimeOfLastTreeConstruction; /*!< holds what it says, only used in connection with FORCETEST */



/* variables for input/output, usually only used on process 0 */

char gadgetmp2::ParameterFile[MAXLEN_FILENAME];  /*!< file name of parameterfile used for starting the simulation */

FILE *gadgetmp2::FdInfo;       /*!< file handle for info.txt log-file. */
FILE *gadgetmp2::FdEnergy;     /*!< file handle for energy.txt log-file. */
FILE *gadgetmp2::FdTimings;    /*!< file handle for timings.txt log-file. */
FILE *gadgetmp2::FdCPU;        /*!< file handle for cpu.txt log-file. */
ofstream gadgetmp2::DEBUG;

#ifdef FORCETEST
FILE *gadgetmp2::FdForceTest;  /*!< file handle for forcetest.txt log-file. */
#endif


my_float gadgetmp2::DriftTable[DRIFT_TABLE_LENGTH];      /*!< table for the cosmological drift factors */
my_float gadgetmp2::GravKickTable[DRIFT_TABLE_LENGTH];   /*!< table for the cosmological kick factor for gravitational forces */
my_float gadgetmp2::HydroKickTable[DRIFT_TABLE_LENGTH];  /*!< table for the cosmological kick factor for hydrodynmical forces */

void *gadgetmp2::CommBuffer=nullptr;   /*!< points to communication buffer, which is used in the domain decomposition, the
                                parallel tree-force computation, the SPH routines, etc. */



/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
struct global_data_all_processes
 gadgetmp2::All;



/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
particle_data *gadgetmp2::P=nullptr;              /*!< holds particle data on local processor */
particle_data_buff *gadgetmp2::DomainPartBuf_s,
                   *gadgetmp2::DomainPartBuf_r;  /*!< buffer for particle data used in domain decomposition */

size_t particle_data_buff::Pos0_off;
size_t particle_data_buff::Pos1_off;
size_t particle_data_buff::Pos2_off;
size_t particle_data_buff::Mass_off;
size_t particle_data_buff::Vel0_off;
size_t particle_data_buff::Vel1_off;
size_t particle_data_buff::Vel2_off;
size_t particle_data_buff::Radius_off;
size_t particle_data_buff::GravAccel0_off;
size_t particle_data_buff::GravAccel1_off;
size_t particle_data_buff::GravAccel2_off;
        #ifdef FORCETEST
size_t particle_data_buff::GravAccelDirect0_off;
size_t particle_data_buff::GravAccelDirect1_off;
size_t particle_data_buff::GravAccelDirect2_off;
        #endif
size_t particle_data_buff::Potential_off;
size_t particle_data_buff::OldAcc_off;
size_t particle_data_buff::ID_off;
size_t particle_data_buff::Type_off;
size_t particle_data_buff::Ti_endstep_off;
size_t particle_data_buff::Ti_begstep_off;
        #ifdef TIMESTEP_LIMITER
size_t particle_data_buff::Ti_sizestep_off;
        #endif
        #ifdef FLEXSTEPS
size_t particle_data_buff::FlexStepGrp_off;
        #endif
size_t particle_data_buff::GravCost_off;
        #ifdef PSEUDOSYMMETRIC
size_t particle_data_buff::AphysOld_off;
        #endif
size_t particle_data_buff::tot_size;
mpfr_prec_t particle_data_buff::prec;

/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
struct sph_particle_data *gadgetmp2::SphP=nullptr;;                        	/*!< holds SPH particle data on local processor */
sph_particle_data_buff *gadgetmp2::DomainSphBuf_s,                 /*!< buffer for SPH particle data in domain decomposition */
                       *gadgetmp2::DomainSphBuf_r;


size_t sph_particle_data_buff::Entropy_off;
size_t sph_particle_data_buff::Density_off;
size_t sph_particle_data_buff::Hsml_off;
size_t sph_particle_data_buff::Left_off;
size_t sph_particle_data_buff::Right_off;
size_t sph_particle_data_buff::NumNgb_off;
size_t sph_particle_data_buff::Pressure_off;
size_t sph_particle_data_buff::DtEntropy_off;
size_t sph_particle_data_buff::HydroAccel0_off;
size_t sph_particle_data_buff::HydroAccel1_off;
size_t sph_particle_data_buff::HydroAccel2_off;
size_t sph_particle_data_buff::VelPred0_off;
size_t sph_particle_data_buff::VelPred1_off;
size_t sph_particle_data_buff::VelPred2_off;
size_t sph_particle_data_buff::DivVel_off;
size_t sph_particle_data_buff::CurlVel_off;
size_t sph_particle_data_buff::Rot0_off;
size_t sph_particle_data_buff::Rot1_off;
size_t sph_particle_data_buff::Rot2_off;
size_t sph_particle_data_buff::DhsmlDensityFactor_off;
size_t sph_particle_data_buff::MaxSignalVel_off;
        #ifdef TIMESTEP_UPDATE
size_t sph_particle_data_buff::FeedbackFlag_off;
size_t sph_particle_data_buff::FeedAccel0_off;
size_t sph_particle_data_buff::FeedAccel1_off;
size_t sph_particle_data_buff::FeedAccel2_off;
        #endif
        #ifdef MORRIS97VISC
size_t sph_particle_data_buff::Alpha_off;
size_t sph_particle_data_buff::DAlphaDt_off;
        #endif
size_t sph_particle_data_buff::tot_size;
mpfr_prec_t sph_particle_data_buff::prec;





/*  Variables for Tree
 */

int gadgetmp2::MaxNodes;		/*!< maximum allowed number of internal nodes */
int gadgetmp2::Numnodestree;	/*!< number of (internal) nodes in each tree */

NODE
 *gadgetmp2::Nodes_base,                   /*!< points to the actual memory allocted for the nodes */
 *gadgetmp2::Nodes;                        /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart]
 				     gives the first allocated node */


int *gadgetmp2::Nextnode;	        /*!< gives next node in tree walk */
int *gadgetmp2::Father;	        /*!< gives parent node in tree    */


extNODE           /*!< this structure holds additional tree-node information which is not needed in the actual gravity computation */
 *gadgetmp2::Extnodes_base,                /*!< points to the actual memory allocted for the extended node information */
 *gadgetmp2::Extnodes;                     /*!< provides shifted access to extended node information, parallel to Nodes/Nodes_base */





/*! Header for the standard file format.
 */
struct io_header
 gadgetmp2::header;  /*!< holds header for snapshot files */



char gadgetmp2::Tab_IO_Labels[IO_NBLOCKS][4];   /*<! This table holds four-byte character tags used for fileformat 2 */



/* global state of system, used for global statistics
 */
struct state_of_system
 gadgetmp2::SysState;      /*<! Structure for storing some global statistics about the simulation. */



/* Various structures for communication
 */
gravdata_in
 *gadgetmp2::GravDataIn,                   /*!< holds particle data to be exported to other processors */
 *gadgetmp2::GravDataGet,                  /*!< holds particle data imported from other processors */
 *gadgetmp2::GravDataResult,               /*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *gadgetmp2::GravDataOut;                  /*!< holds partial results received from other processors. This will overwrite the GravDataIn array */

size_t gravdata_in::u0_off=0;
size_t gravdata_in::u1_off=0;
size_t gravdata_in::u2_off=0;
size_t gravdata_in::Ninteractions_off=0;
size_t gravdata_in::OldAcc_off=0;
#ifdef UNEQUALSOFTENINGS
size_t gravdata_in::Type_off=0;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
size_t gravdata_in::Soft_off=0;
#endif
#endif
size_t gravdata_in::tot_size=0;
mpfr_prec_t gravdata_in::prec=0;

struct gravdata_index
 *gadgetmp2::GravDataIndexTable;           /*!< the particles to be exported are grouped by task-number. This table allows the results to be disentangled again and to be assigned to the correct particle */



struct densdata_in
 *gadgetmp2::DensDataIn,                   /*!< holds particle data for SPH density computation to be exported to other processors */
 *gadgetmp2::DensDataGet;                  /*!< holds imported particle data for SPH density computation */

size_t densdata_in::Pos0_off;
size_t densdata_in::Pos1_off;
size_t densdata_in::Pos2_off;
size_t densdata_in::Vel0_off;
size_t densdata_in::Vel1_off;
size_t densdata_in::Vel2_off;
size_t densdata_in::Hsml_off;
size_t densdata_in::Index_off;
size_t densdata_in::Task_off;
size_t densdata_in::tot_size;
mpfr_prec_t densdata_in::prec;

struct densdata_out
 *gadgetmp2::DensDataResult,               /*!< stores the locally computed SPH density results for imported particles */
 *gadgetmp2::DensDataPartialResult;        /*!< imported partial SPH density results from other processors */


size_t densdata_out::Rho_off;
size_t densdata_out::Div_off;
size_t densdata_out::Rot0_off;
size_t densdata_out::Rot1_off;
size_t densdata_out::Rot2_off;
size_t densdata_out::DhsmlDensity_off;
size_t densdata_out::Ngb_off;
size_t densdata_out::tot_size;
mpfr_prec_t densdata_out::prec;


struct hydrodata_in
 *gadgetmp2::HydroDataIn,                  /*!< holds particle data for SPH hydro-force computation to be exported to other processors */
 *gadgetmp2::HydroDataGet;                 /*!< holds imported particle data for SPH hydro-force computation */


size_t hydrodata_in::Pos0_off;
size_t hydrodata_in::Pos1_off;
size_t hydrodata_in::Pos2_off;
size_t hydrodata_in::Vel0_off;
size_t hydrodata_in::Vel1_off;
size_t hydrodata_in::Vel2_off;
size_t hydrodata_in::Hsml_off;
size_t hydrodata_in::Mass_off;
size_t hydrodata_in::Density_off;
size_t hydrodata_in::Pressure_off;
size_t hydrodata_in::F1_off;
size_t hydrodata_in::DhsmlDensityFactor_off;
size_t hydrodata_in::Timestep_off;
size_t hydrodata_in::Task_off;
size_t hydrodata_in::Index_off;
        #ifdef MORRIS97VISC
size_t hydrodata_in::Alpha_off;
        #endif
size_t hydrodata_in::tot_size;
mpfr_prec_t hydrodata_in::prec;

struct hydrodata_out
 *gadgetmp2::HydroDataResult,              /*!< stores the locally computed SPH hydro results for imported particles */
 *gadgetmp2::HydroDataPartialResult;       /*!< imported partial SPH hydro-force results from other processors */

size_t hydrodata_out::Acc0_off;
size_t hydrodata_out::Acc1_off;
size_t hydrodata_out::Acc2_off;
size_t hydrodata_out::DtEntropy_off;
size_t hydrodata_out::MaxSignalVel_off;
size_t hydrodata_out::tot_size;
mpfr_prec_t hydrodata_out::prec;

#ifdef TIMESTEP_LIMITER
timedata_in *gadgetmp2::TimeDataIn, *gadgetmp2::TimeDataGet;

size_t timedata_in::Pos0_off=0;
size_t timedata_in::Pos1_off=0;
size_t timedata_in::Pos2_off=0;
size_t timedata_in::Hsml_off=0;
size_t timedata_in::Size_off=0;
size_t timedata_in::Begin_off=0;
size_t timedata_in::Index_off=0;
size_t timedata_in::Task_off=0;
size_t timedata_in::tot_size=0;
mpfr_prec_t timedata_in::prec=0;

All_Reduce_buff *gadgetmp2::all_reduce_buff;

size_t All_Reduce_buff::tot_size=0;
mpfr_prec_t All_Reduce_buff::prec=0;


#endif

#ifndef NOMPI
#include <mpi.h>
MPI_Comm gadgetmp2::GADGET_WORLD;
#endif


