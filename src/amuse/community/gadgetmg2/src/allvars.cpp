/*! \file allvars.c
 *  \brief provides instances of all global variables.
 */

#ifndef NOMPI
#include <mpi.h>
#endif

#include "proto.hpp"


int gadgetmg2::ThisTask;		/*!< the rank of the local processor */
int gadgetmg2::NTask;               /*!< number of processors */
int gadgetmg2::PTask;	        /*!< smallest integer such that NTask <= 2^PTask */

int gadgetmg2::NumPart;		/*!< number of particles on the LOCAL processor */
int gadgetmg2::N_gas;		/*!< number of gas particles on the LOCAL processor  */
long long gadgetmg2::Ntype[6];      /*!< total number of particles of each type */
int gadgetmg2::NtypeLocal[6];       /*!< local number of particles of each type */

int gadgetmg2::NumForceUpdate;      /*!< number of active particles on local processor in current timestep  */
int gadgetmg2::NumSphUpdate;        /*!< number of active SPH particles on local processor in current timestep  */

double gadgetmg2::CPUThisRun;	/*!< Sums the CPU time for the process (current submission only) */


int gadgetmg2::RestartFlag;         /*!< taken from command line used to start code. 0 is normal start-up from
                                     initial conditions, 1 is resuming a run from a set of restart files, while 2
                                     marks a restart from a snapshot file. */

char *gadgetmg2::Exportflag;        /*!< Buffer used for flagging whether a particle needs to be exported to another process */

int  *gadgetmg2::Ngblist;           /*!< Buffer to hold indices of neighbours retrieved by the neighbour search routines */

int gadgetmg2::TreeReconstructFlag; /*!< Signals that a new tree needs to be constructed */

int gadgetmg2::Flag_FullStep;       /*!< This flag signals that the current step involves all particles */

int gadgetmg2::ZeroTimestepEncountered;  /*!< Flag used by AMUSE. When a particle is assigned a timestep of zero, an
                                     exception is raised instead of forcing the application to exit. */


gsl_rng *gadgetmg2::random_generator; /*!< the employed random number generator of the GSL library */

double gadgetmg2::RndTable[RNDTABLE]; /*!< Hold a table with random numbers, refreshed every timestep */


double gadgetmg2::DomainCorner[3];    /*!< gives the lower left corner of simulation volume */
double gadgetmg2::DomainCenter[3];    /*!< gives the center of simulation volume */
double gadgetmg2::DomainLen;          /*!< gives the (maximum) side-length of simulation volume */
double gadgetmg2::DomainFac;          /*!< factor used for converting particle coordinates to a Peano-Hilbert mesh covering the simulation volume */
int    gadgetmg2::DomainMyStart;      /*!< first domain mesh cell that resides on the local processor */
int    gadgetmg2::DomainMyLast;       /*!< last domain mesh cell that resides on the local processor */
int    *gadgetmg2::DomainStartList;   /*!< a table that lists the first domain mesh cell for all processors */
int    *gadgetmg2::DomainEndList;     /*!< a table that lists the last domain mesh cell for all processors */
double *gadgetmg2::DomainWork;        /*!< a table that gives the total "work" due to the particles stored by each processor */
int    *gadgetmg2::DomainCount;       /*!< a table that gives the total number of particles held by each processor */
int    *gadgetmg2::DomainCountSph;    /*!< a table that gives the total number of SPH particles held by each processor */

int    *gadgetmg2::DomainTask;        /*!< this table gives for each leaf of the top-level tree the processor it was assigned to */
int    *gadgetmg2::DomainNodeIndex;   /*!< this table gives for each leaf of the top-level tree the corresponding node of the gravitational tree */
double  *gadgetmg2::DomainTreeNodeLen; /*!< this table gives for each leaf of the top-level tree the side-length of the corresponding node of the gravitational tree */
double  *gadgetmg2::DomainHmax;        /*!< this table gives for each leaf of the top-level tree the maximum SPH smoothing length among the particles of the corresponding node of the gravitational tree */

struct DomainNODE
 *gadgetmg2::DomainMoment;                    /*!< this table stores for each node of the top-level tree corresponding node data from the gravitational tree */

peanokey *gadgetmg2::DomainKeyBuf;    /*!< this points to a buffer used during the exchange of particle data */

peanokey *gadgetmg2::Key;             /*!< a table used for storing Peano-Hilbert keys for particles */
peanokey *gadgetmg2::KeySorted;       /*!< holds a sorted table of Peano-Hilbert keys for all particles, used to construct top-level tree */


int gadgetmg2::NTopnodes;             /*!< total number of nodes in top-level tree */
int gadgetmg2::NTopleaves;            /*!< number of leaves in top-level tree. Each leaf can be assigned to a different processor */

struct topnode_data
 *gadgetmg2::TopNodes;                      /*!< points to the root node of the top-level tree */


double gadgetmg2::TimeOfLastTreeConstruction; /*!< holds what it says, only used in connection with FORCETEST */



/* variables for input/output, usually only used on process 0 */

char gadgetmg2::ParameterFile[MAXLEN_FILENAME];  /*!< file name of parameterfile used for starting the simulation */

FILE *gadgetmg2::FdInfo;       /*!< file handle for info.txt log-file. */
FILE *gadgetmg2::FdEnergy;     /*!< file handle for energy.txt log-file. */
FILE *gadgetmg2::FdTimings;    /*!< file handle for timings.txt log-file. */
FILE *gadgetmg2::FdCPU;        /*!< file handle for cpu.txt log-file. */
ofstream gadgetmg2::DEBUG;

#ifdef FORCETEST
FILE *gadgetmg2::FdForceTest;  /*!< file handle for forcetest.txt log-file. */
#endif


double gadgetmg2::DriftTable[DRIFT_TABLE_LENGTH];      /*!< table for the cosmological drift factors */
double gadgetmg2::GravKickTable[DRIFT_TABLE_LENGTH];   /*!< table for the cosmological kick factor for gravitational forces */
double gadgetmg2::HydroKickTable[DRIFT_TABLE_LENGTH];  /*!< table for the cosmological kick factor for hydrodynmical forces */

void *gadgetmg2::CommBuffer;   /*!< points to communication buffer, which is used in the domain decomposition, the
                                parallel tree-force computation, the SPH routines, etc. */



/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
struct global_data_all_processes
 gadgetmg2::All;



/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
struct particle_data
 *gadgetmg2::P,              /*!< holds particle data on local processor */
 *gadgetmg2::DomainPartBuf;  /*!< buffer for particle data used in domain decomposition */


/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
struct sph_particle_data
 *gadgetmg2::SphP,                        	/*!< holds SPH particle data on local processor */
 *gadgetmg2::DomainSphBuf;                 /*!< buffer for SPH particle data in domain decomposition */





/*  Variables for Tree
 */

int gadgetmg2::MaxNodes;		/*!< maximum allowed number of internal nodes */
int gadgetmg2::Numnodestree;	/*!< number of (internal) nodes in each tree */

struct NODE
 *gadgetmg2::Nodes_base,                   /*!< points to the actual memory allocted for the nodes */
 *gadgetmg2::Nodes;                        /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart]
 				     gives the first allocated node */


int *gadgetmg2::Nextnode;	        /*!< gives next node in tree walk */
int *gadgetmg2::Father;	        /*!< gives parent node in tree    */


struct extNODE           /*!< this structure holds additional tree-node information which is not needed in the actual gravity computation */
 *gadgetmg2::Extnodes_base,                /*!< points to the actual memory allocted for the extended node information */
 *gadgetmg2::Extnodes;                     /*!< provides shifted access to extended node information, parallel to Nodes/Nodes_base */





/*! Header for the standard file format.
 */
struct io_header
 gadgetmg2::header;  /*!< holds header for snapshot files */



char gadgetmg2::Tab_IO_Labels[IO_NBLOCKS][4];   /*<! This table holds four-byte character tags used for fileformat 2 */



/* global state of system, used for global statistics
 */
struct state_of_system
 gadgetmg2::SysState;      /*<! Structure for storing some global statistics about the simulation. */



/* Various structures for communication
 */
struct gravdata_in
 *gadgetmg2::GravDataIn,                   /*!< holds particle data to be exported to other processors */
 *gadgetmg2::GravDataGet,                  /*!< holds particle data imported from other processors */
 *gadgetmg2::GravDataResult,               /*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *gadgetmg2::GravDataOut;                  /*!< holds partial results received from other processors. This will overwrite the GravDataIn array */

struct gravdata_index
 *gadgetmg2::GravDataIndexTable;           /*!< the particles to be exported are grouped by task-number. This table allows the results to be disentangled again and to be assigned to the correct particle */



struct densdata_in
 *gadgetmg2::DensDataIn,                   /*!< holds particle data for SPH density computation to be exported to other processors */
 *gadgetmg2::DensDataGet;                  /*!< holds imported particle data for SPH density computation */

struct densdata_out
 *gadgetmg2::DensDataResult,               /*!< stores the locally computed SPH density results for imported particles */
 *gadgetmg2::DensDataPartialResult;        /*!< imported partial SPH density results from other processors */



struct hydrodata_in
 *gadgetmg2::HydroDataIn,                  /*!< holds particle data for SPH hydro-force computation to be exported to other processors */
 *gadgetmg2::HydroDataGet;                 /*!< holds imported particle data for SPH hydro-force computation */

struct hydrodata_out
 *gadgetmg2::HydroDataResult,              /*!< stores the locally computed SPH hydro results for imported particles */
 *gadgetmg2::HydroDataPartialResult;       /*!< imported partial SPH hydro-force results from other processors */

#ifdef TIMESTEP_LIMITER
struct timedata_in *gadgetmg2::TimeDataIn, *gadgetmg2::TimeDataGet;
#endif

#ifndef NOMPI
#include <mpi.h>
MPI_Comm gadgetmg2::GADGET_WORLD;
#endif


