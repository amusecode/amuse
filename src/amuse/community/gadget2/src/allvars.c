/*! \file allvars.c
 *  \brief provides instances of all global variables.
 */

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "tags.h"
#include "allvars.h"


int ThisTask;		/*!< the rank of the local processor */
int NTask;               /*!< number of processors */
int PTask;	        /*!< smallest integer such that NTask <= 2^PTask */

int NumPart;		/*!< number of particles on the LOCAL processor */
int N_gas;		/*!< number of gas particles on the LOCAL processor  */
long long Ntype[6];      /*!< total number of particles of each type */
int NtypeLocal[6];       /*!< local number of particles of each type */

int NumForceUpdate;      /*!< number of active particles on local processor in current timestep  */
int NumSphUpdate;        /*!< number of active SPH particles on local processor in current timestep  */

double CPUThisRun;	/*!< Sums the CPU time for the process (current submission only) */


int RestartFlag;         /*!< taken from command line used to start code. 0 is normal start-up from
                                     initial conditions, 1 is resuming a run from a set of restart files, while 2
                                     marks a restart from a snapshot file. */

char *Exportflag;        /*!< Buffer used for flagging whether a particle needs to be exported to another process */

int  *Ngblist;           /*!< Buffer to hold indices of neighbours retrieved by the neighbour search routines */

int TreeReconstructFlag; /*!< Signals that a new tree needs to be constructed */

int Flag_FullStep;       /*!< This flag signals that the current step involves all particles */

int ZeroTimestepEncountered;  /*!< Flag used by AMUSE. When a particle is assigned a timestep of zero, an
                                     exception is raised instead of forcing the application to exit. */


gsl_rng *random_generator; /*!< the employed random number generator of the GSL library */

double RndTable[RNDTABLE]; /*!< Hold a table with random numbers, refreshed every timestep */


double DomainCorner[3];    /*!< gives the lower left corner of simulation volume */
double DomainCenter[3];    /*!< gives the center of simulation volume */
double DomainLen;          /*!< gives the (maximum) side-length of simulation volume */
double DomainFac;          /*!< factor used for converting particle coordinates to a Peano-Hilbert mesh covering the simulation volume */
int    DomainMyStart;      /*!< first domain mesh cell that resides on the local processor */
int    DomainMyLast;       /*!< last domain mesh cell that resides on the local processor */
int    *DomainStartList;   /*!< a table that lists the first domain mesh cell for all processors */
int    *DomainEndList;     /*!< a table that lists the last domain mesh cell for all processors */
double *DomainWork;        /*!< a table that gives the total "work" due to the particles stored by each processor */
int    *DomainCount;       /*!< a table that gives the total number of particles held by each processor */
int    *DomainCountSph;    /*!< a table that gives the total number of SPH particles held by each processor */

int    *DomainTask;        /*!< this table gives for each leaf of the top-level tree the processor it was assigned to */
int    *DomainNodeIndex;   /*!< this table gives for each leaf of the top-level tree the corresponding node of the gravitational tree */
FLOAT  *DomainTreeNodeLen; /*!< this table gives for each leaf of the top-level tree the side-length of the corresponding node of the gravitational tree */
FLOAT  *DomainHmax;        /*!< this table gives for each leaf of the top-level tree the maximum SPH smoothing length among the particles of the corresponding node of the gravitational tree */

struct DomainNODE
 *DomainMoment;                    /*!< this table stores for each node of the top-level tree corresponding node data from the gravitational tree */

peanokey *DomainKeyBuf;    /*!< this points to a buffer used during the exchange of particle data */

peanokey *Key;             /*!< a table used for storing Peano-Hilbert keys for particles */
peanokey *KeySorted;       /*!< holds a sorted table of Peano-Hilbert keys for all particles, used to construct top-level tree */


int NTopnodes;             /*!< total number of nodes in top-level tree */
int NTopleaves;            /*!< number of leaves in top-level tree. Each leaf can be assigned to a different processor */

struct topnode_data
 *TopNodes;                      /*!< points to the root node of the top-level tree */


double TimeOfLastTreeConstruction; /*!< holds what it says, only used in connection with FORCETEST */



/* variables for input/output, usually only used on process 0 */

char ParameterFile[MAXLEN_FILENAME];  /*!< file name of parameterfile used for starting the simulation */

FILE *FdInfo;       /*!< file handle for info.txt log-file. */
FILE *FdEnergy;     /*!< file handle for energy.txt log-file. */
FILE *FdTimings;    /*!< file handle for timings.txt log-file. */
FILE *FdCPU;        /*!< file handle for cpu.txt log-file. */

#ifdef FORCETEST
FILE *FdForceTest;  /*!< file handle for forcetest.txt log-file. */
#endif


double DriftTable[DRIFT_TABLE_LENGTH];      /*!< table for the cosmological drift factors */
double GravKickTable[DRIFT_TABLE_LENGTH];   /*!< table for the cosmological kick factor for gravitational forces */
double HydroKickTable[DRIFT_TABLE_LENGTH];  /*!< table for the cosmological kick factor for hydrodynmical forces */

void *CommBuffer;   /*!< points to communication buffer, which is used in the domain decomposition, the
                                parallel tree-force computation, the SPH routines, etc. */



/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
struct global_data_all_processes
 All;



/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
struct particle_data
 *P,              /*!< holds particle data on local processor */
 *DomainPartBuf;  /*!< buffer for particle data used in domain decomposition */


/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
struct sph_particle_data
 *SphP,                        	/*!< holds SPH particle data on local processor */
 *DomainSphBuf;                 /*!< buffer for SPH particle data in domain decomposition */





/*  Variables for Tree
 */

int MaxNodes;		/*!< maximum allowed number of internal nodes */
int Numnodestree;	/*!< number of (internal) nodes in each tree */

struct NODE
 *Nodes_base,                   /*!< points to the actual memory allocted for the nodes */
 *Nodes;                        /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 
 				     gives the first allocated node */


int *Nextnode;	        /*!< gives next node in tree walk */
int *Father;	        /*!< gives parent node in tree    */


struct extNODE           /*!< this structure holds additional tree-node information which is not needed in the actual gravity computation */
 *Extnodes_base,                /*!< points to the actual memory allocted for the extended node information */
 *Extnodes;                     /*!< provides shifted access to extended node information, parallel to Nodes/Nodes_base */





/*! Header for the standard file format.
 */
struct io_header
 header;  /*!< holds header for snapshot files */



char Tab_IO_Labels[IO_NBLOCKS][4];   /*<! This table holds four-byte character tags used for fileformat 2 */



/* global state of system, used for global statistics
 */
struct state_of_system
 SysState;      /*<! Structure for storing some global statistics about the simulation. */
 


/* Various structures for communication
 */
struct gravdata_in
 *GravDataIn,                   /*!< holds particle data to be exported to other processors */
 *GravDataGet,                  /*!< holds particle data imported from other processors */
 *GravDataResult,               /*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;                  /*!< holds partial results received from other processors. This will overwrite the GravDataIn array */

struct gravdata_index
 *GravDataIndexTable;           /*!< the particles to be exported are grouped by task-number. This table allows the results to be disentangled again and to be assigned to the correct particle */



struct densdata_in
 *DensDataIn,                   /*!< holds particle data for SPH density computation to be exported to other processors */
 *DensDataGet;                  /*!< holds imported particle data for SPH density computation */

struct densdata_out
 *DensDataResult,               /*!< stores the locally computed SPH density results for imported particles */
 *DensDataPartialResult;        /*!< imported partial SPH density results from other processors */



struct hydrodata_in
 *HydroDataIn,                  /*!< holds particle data for SPH hydro-force computation to be exported to other processors */
 *HydroDataGet;                 /*!< holds imported particle data for SPH hydro-force computation */

struct hydrodata_out
 *HydroDataResult,              /*!< stores the locally computed SPH hydro results for imported particles */
 *HydroDataPartialResult;       /*!< imported partial SPH hydro-force results from other processors */


