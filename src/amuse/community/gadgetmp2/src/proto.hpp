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
/*! \file proto.hpp
 *  \brief this file contains all function prototypes of the code
 */


#ifndef ALLVARS_H
//#include "allvars.hpp"
#endif

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "tags.hpp"
#include "allvars.hpp"
#include <iostream>
#include <fstream>
using namespace std;
class gadgetmp2{
public:
static int ThisTask;		/*!< the rank of the local processor */
static int NTask;               /*!< number of processors */
static int PTask;	        /*!< smallest integer such that NTask <= 2^PTask */

static int NumPart;		/*!< number of particles on the LOCAL processor */
static int N_gas;		/*!< number of gas particles on the LOCAL processor  */
static long long Ntype[6];      /*!< total number of particles of each type */
static int NtypeLocal[6];       /*!< local number of particles of each type */

static int NumForceUpdate;      /*!< number of active particles on local processor in current timestep  */
static int NumSphUpdate;        /*!< number of active SPH particles on local processor in current timestep  */

static double CPUThisRun;	/*!< Sums the CPU time for the process (current submission only) */


static int RestartFlag;         /*!< taken from command line used to start code. 0 is normal start-up from
                                     initial conditions, 1 is resuming a run from a set of restart files, while 2
                                     marks a restart from a snapshot file. */

static char *Exportflag;        /*!< Buffer used for flagging whether a particle needs to be exported to another process */

static int  *Ngblist;           /*!< Buffer to hold indices of neighbours retrieved by the neighbour search routines */

static int TreeReconstructFlag; /*!< Signals that a new tree needs to be constructed */

static int Flag_FullStep;       /*!< This flag signals that the current step involves all particles */

static int ZeroTimestepEncountered;  /*!< Flag used by AMUSE. When a particle is assigned a timestep of zero, an
                                     exception is raised instead of forcing the application to exit. */


static gsl_rng *random_generator; /*!< the employed random number generator of the GSL library */

static my_float RndTable[RNDTABLE]; /*!< Hold a table with random numbers, refreshed every timestep */


static my_float DomainCorner[3];    /*!< gives the lower left corner of simulation volume */
static my_float DomainCenter[3];    /*!< gives the center of simulation volume */
static my_float DomainLen;          /*!< gives the (maximum) side-length of simulation volume */
static my_float DomainFac;          /*!< factor used for converting particle coordinates to a Peano-Hilbert mesh covering the simulation volume */
static int    DomainMyStart;      /*!< first domain mesh cell that resides on the local processor */
static int    DomainMyLast;       /*!< last domain mesh cell that resides on the local processor */
static int    *DomainStartList;   /*!< a table that lists the first domain mesh cell for all processors */
static int    *DomainEndList;     /*!< a table that lists the last domain mesh cell for all processors */
static my_float *DomainWork;        /*!< a table that gives the total "work" due to the particles stored by each processor */
static int    *DomainCount;       /*!< a table that gives the total number of particles held by each processor */
static int    *DomainCountSph;    /*!< a table that gives the total number of SPH particles held by each processor */

static int    *DomainTask;        /*!< this table gives for each leaf of the top-level tree the processor it was assigned to */
static int    *DomainNodeIndex;   /*!< this table gives for each leaf of the top-level tree the corresponding node of the gravitational tree */
static my_float  *DomainTreeNodeLen; /*!< this table gives for each leaf of the top-level tree the side-length of the corresponding node of the gravitational tree */
static All_Reduce_buff  *DomainHmax;        /*!< this table gives for each leaf of the top-level tree the maximum SPH smoothing length among the particles of the corresponding node of the gravitational tree */

static DomainNODE
 *DomainMoment;                    /*!< this table stores for each node of the top-level tree corresponding node data from the gravitational tree */

static peanokey *DomainKeyBuf;    /*!< this points to a buffer used during the exchange of particle data */

static peanokey *Key;             /*!< a table used for storing Peano-Hilbert keys for particles */
static peanokey *KeySorted;       /*!< holds a sorted table of Peano-Hilbert keys for all particles, used to construct top-level tree */


static int NTopnodes;             /*!< total number of nodes in top-level tree */
static int NTopleaves;            /*!< number of leaves in top-level tree. Each leaf can be assigned to a different processor */

static struct topnode_data
 *TopNodes;                      /*!< points to the root node of the top-level tree */


static my_float TimeOfLastTreeConstruction; /*!< holds what it says, only used in connection with FORCETEST */



/* variables for input/output, usually only used on process 0 */

static char ParameterFile[MAXLEN_FILENAME];  /*!< file name of parameterfile used for starting the simulation */

static FILE *FdInfo;       /*!< file handle for info.txt log-file. */
static FILE *FdEnergy;     /*!< file handle for energy.txt log-file. */
static FILE *FdTimings;    /*!< file handle for timings.txt log-file. */
static FILE *FdCPU;        /*!< file handle for cpu.txt log-file. */
static ofstream DEBUG;

#ifdef FORCETEST
static FILE *FdForceTest;  /*!< file handle for forcetest.txt log-file. */
#endif


static my_float DriftTable[DRIFT_TABLE_LENGTH];      /*!< table for the cosmological drift factors */
static my_float GravKickTable[DRIFT_TABLE_LENGTH];   /*!< table for the cosmological kick factor for gravitational forces */
static my_float HydroKickTable[DRIFT_TABLE_LENGTH];  /*!< table for the cosmological kick factor for hydrodynmical forces */

static void *CommBuffer;   /*!< points to communication buffer, which is used in the domain decomposition, the
                                parallel tree-force computation, the SPH routines, etc. */



/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
static struct global_data_all_processes
 All;



/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
static struct particle_data *P;              /*!< holds particle data on local processor */
static particle_data_buff *DomainPartBuf_s,  /*!< buffer for particle data used in domain decomposition */
                          *DomainPartBuf_r;

/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
static struct sph_particle_data *SphP;                     	/*!< holds SPH particle data on local processor */
static sph_particle_data_buff *DomainSphBuf_s,
                              *DomainSphBuf_r;            /*!< buffer for SPH particle data in domain decomposition */





/*  Variables for Tree
 */

static int MaxNodes;		/*!< maximum allowed number of internal nodes */
static int Numnodestree;	/*!< number of (internal) nodes in each tree */

static struct NODE
 *Nodes_base,                   /*!< points to the actual memory allocted for the nodes */
 *Nodes;                        /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart]
 				     gives the first allocated node */


static int *Nextnode;	        /*!< gives next node in tree walk */
static int *Father;	        /*!< gives parent node in tree    */


static struct extNODE           /*!< this structure holds additional tree-node information which is not needed in the actual gravity computation */
 *Extnodes_base,                /*!< points to the actual memory allocted for the extended node information */
 *Extnodes;                     /*!< provides shifted access to extended node information, parallel to Nodes/Nodes_base */





/*! Header for the standard file format.
 */
static struct io_header
 header;  /*!< holds header for snapshot files */



static char Tab_IO_Labels[IO_NBLOCKS][4];   /*<! This table holds four-byte character tags used for fileformat 2 */



/* global state of system, used for global statistics
 */
static struct state_of_system
 SysState;      /*<! Structure for storing some global statistics about the simulation. */



/* Various structures for communication
 */
static gravdata_in
 *GravDataIn,                   /*!< holds particle data to be exported to other processors */
 *GravDataGet,                  /*!< holds particle data imported from other processors */
 *GravDataResult,               /*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;                  /*!< holds partial results received from other processors. This will overwrite the GravDataIn array */

static struct gravdata_index
 *GravDataIndexTable;           /*!< the particles to be exported are grouped by task-number. This table allows the results to be disentangled again and to be assigned to the correct particle */



static struct densdata_in
 *DensDataIn,                   /*!< holds particle data for SPH density computation to be exported to other processors */
 *DensDataGet;                  /*!< holds imported particle data for SPH density computation */

static struct densdata_out
 *DensDataResult,               /*!< stores the locally computed SPH density results for imported particles */
 *DensDataPartialResult;        /*!< imported partial SPH density results from other processors */



static struct hydrodata_in
 *HydroDataIn,                  /*!< holds particle data for SPH hydro-force computation to be exported to other processors */
 *HydroDataGet;                 /*!< holds imported particle data for SPH hydro-force computation */

static struct hydrodata_out
 *HydroDataResult,              /*!< stores the locally computed SPH hydro results for imported particles */
 *HydroDataPartialResult;       /*!< imported partial SPH hydro-force results from other processors */

#ifdef TIMESTEP_LIMITER
static struct timedata_in *TimeDataIn, *TimeDataGet;
#endif

static All_Reduce_buff *all_reduce_buff;

#ifndef NOMPI
#include <mpi.h>
static MPI_Comm GADGET_WORLD;
#endif

static void   advance_and_find_timesteps(void);
static void   allocate_commbuffers(void);
static void   allocate_memory(void);
//static void   begrun(void);
//static int    blockpresent(enum iofields blocknr);
static void   catch_abort(int sig);
static void   catch_fatal(int sig);
//static void   check_omega(void);
static void   close_outputfiles(void);
static int    compare_key(const void *a, const void *b);
static void   compute_accelerations(int mode);
static void   compute_global_quantities_of_system(void);
static void   compute_potential(void);
static int    dens_compare_key(const void *a, const void *b);
static void   density(void);
static void   density_decouple(void);
static void   density_evaluate(int i, int mode);

//static void   distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master, int *last);
static my_float dmax(my_float, my_float);
static my_float dmin(my_float, my_float);
static void   do_box_wrapping(void);

static void   domain_Decomposition(void);
static int    domain_compare_key(const void *a, const void *b);
static int    domain_compare_toplist(const void *a, const void *b);
static void   domain_countToGo(void);
static void   domain_decompose(void);
static void   domain_determineTopTree(void);
static void   domain_exchangeParticles(int partner, int sphflag, int send_count, int recv_count);
static void   domain_findExchangeNumbers(int task, int partner, int sphflag, int *send, int *recv);
static void   domain_findExtent(void);
static int    domain_findSplit(int cpustart, int ncpu, int first, int last);
static void   domain_shiftSplit(void);
static void   domain_sumCost(void);
static void   domain_topsplit(int node, peanokey startkey);
static void   domain_topsplit_local(int node, peanokey startkey);

static double drift_integ(double a, void *param);
static void   dump_particles(void);
//static void   empty_read_buffer(enum iofields blocknr, int offset, int pc, int type);
static void   endrun(int);
static void   energy_statistics(void);
static void   every_timestep_stuff(void);

static void   ewald_corr(my_float dx, my_float dy, my_float dz, my_float *fper);
static void   ewald_force(int ii, int jj, int kk, my_float x[3], my_float force[3]);
static void   ewald_init(void);
static my_float ewald_pot_corr(my_float dx, my_float dy, my_float dz);
static my_float ewald_psi(my_float x[3]);

//static void   fill_Tab_IO_Labels(void);
//static void   fill_write_buffer(enum iofields blocknr, int *pindex, int pc, int type);
static void   find_dt_displacement_constraint(my_float hfac);
//static int    find_files(char *fname);
//static int    find_next_outputtime(int time);
static void   find_next_sync_point_and_drift(void);

static void   force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z, int *nodecount, int *nextfree);
static void   force_exchange_pseudodata(void);
static void   force_flag_localnodes(void);
static void   force_insert_pseudo_particles(void);
static void   force_setupnonrecursive(int no);
static void   force_treeallocate(int maxnodes, int maxpart);
static int    force_treebuild(int npart);
static int    force_treebuild_single(int npart);
static int    force_treeevaluate(int target, int mode, double *ewaldcountsum);
static int    force_treeevaluate_direct(int target, int mode);
static int    force_treeevaluate_ewald_correction(int target, int mode, my_float pos_x, my_float pos_y, my_float pos_z, my_float aold);
static void   force_treeevaluate_potential(int target, int type);
static void   force_treeevaluate_potential_shortrange(int target, int mode);
static int    force_treeevaluate_shortrange(int target, int mode);
static void   force_treefree(void);
static void   force_treeupdate_pseudos(void);
static void   force_update_hmax(void);
static void   force_update_len(void);
static void   force_update_node(int no, int flag);
static void   force_update_node_hmax_local(void);
static void   force_update_node_hmax_toptree(void);
static void   force_update_node_len_local(void);
static void   force_update_node_len_toptree(void);
static void   force_update_node_recursive(int no, int sib, int father);
static void   force_update_pseudoparticles(void);
static void   force_update_size_of_parent_node(int no);

static void   free_memory(void);

static int    get_bytes_per_blockelement(enum iofields blocknr);
//static void   get_dataset_name(enum iofields blocknr, char *buf);
//static int    get_datatype_in_block(enum iofields blocknr);
static my_float get_drift_factor(int time0, int time1);
static my_float get_gravkick_factor(int time0, int time1);
static my_float get_hydrokick_factor(int time0, int time1);
//static int    get_particles_in_block(enum iofields blocknr, int *typelist);
static my_float get_random_number(int id);
static int    get_timestep(int p, my_float *a, int flag);
//static int    get_values_per_blockelement(enum iofields blocknr);

static int    grav_tree_compare_key(const void *a, const void *b);
static void   gravity_forcetest(void);
static void   gravity_tree(void);
static void   gravity_tree_shortrange(void);
static double gravkick_integ(double a, void *param);

static int    hydro_compare_key(const void *a, const void *b);
static void   hydro_evaluate(int target, int mode);
static void   hydro_force(void);
static double hydrokick_integ(double a, void *param);

static int    imax(int, int);
static int    imin(int, int);

//static void   init(void);
static void   init_drift_table(void);
static void   init_peano_map(void);

static void   long_range_force(void);
static void   long_range_init(void);
static void   long_range_init_regionsize(void);
static void   move_particles(int time0, int time1);
//static size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
static size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);

static int    ngb_clear_buf(my_float searchcenter[3], my_float hguess, int numngb);
static void   ngb_treeallocate(int npart);
static void   ngb_treebuild(void);
static int    ngb_treefind_pairs(my_float searchcenter[3], my_float hsml, int *startnode);
static int    ngb_treefind_variable(my_float searchcenter[3], my_float hguess, int *startnode);
static void   ngb_treefree(void);
static void   ngb_treesearch(int);
static void   ngb_treesearch_pairs(int);
static void   ngb_update_nodes(void);

static void   open_outputfiles(void);

static peanokey peano_hilbert_key(int x, int y, int z, int bits);
static void   peano_hilbert_order(void);

//static void   pm_init_nonperiodic(void);
//static void   pm_init_nonperiodic_allocate(int dimprod);
//static void   pm_init_nonperiodic_free(void);
//static void   pm_init_periodic(void);
//static void   pm_init_periodic_allocate(int dimprod);
//static void   pm_init_periodic_free(void);
//static void   pm_init_regionsize(void);
//static void   pm_setup_nonperiodic_kernel(void);
//static int    pmforce_nonperiodic(int grnr);
//static void   pmforce_periodic(void);
//static int    pmpotential_nonperiodic(int grnr);
//static void   pmpotential_periodic(void);

#ifdef NOPOW
static my_float pow(my_float, my_float);  /* on some old DEC Alphas, the correct prototype for pow() is missing, even when math.h is included */
#endif

//static void   read_file(char *fname, int readTask, int lastTask);
//static void   read_header_attributes_in_hdf5(char *fname);
//static void   read_ic(char *fname);
//static int    read_outputlist(char *fname);
//static void   read_parameter_file(char *fname);
//static void   readjust_timebase(my_float TimeMax_old, my_float TimeMax_new);

static void   reorder_gas(void);
static void   reorder_particles(void);
//static void   restart(int mod);
static void   run(void);
//static void   savepositions(int num);

static double second(void);

//void   seed_glass(void);
static void   set_random_numbers(void);
static void   set_softenings(void);
static void   set_units(void);

static void   setup_smoothinglengths(void);
static void   statistics(void);
static void   terminate_processes(void);
static double timediff(double t0, double t1);

#ifdef HAVE_HDF5
static void   write_header_attributes_in_hdf5(hid_t handle);
#endif
//static void   write_file(char *fname, int readTask, int lastTask);
static void   write_pid_file(void);


#ifdef TIMESTEP_LIMITER
static   int  time_compare_key(const void *a, const void *b);
static   int  time_limiter_evaluate(int target, int mode);
#endif

#ifdef TIMESTEP_UPDATE
static   void get_sigvel(void);
static   void get_sigvel_evaluate(int target, int mode);
#endif
static void hydro_state_evaluate(my_float h, my_float pos[3], my_float vel[3], my_float *numngb,
    my_float *dhsml_out, my_float *rho_out, my_float *rhov_out, my_float *rhov2_out, my_float
    *rhoe_out);
static double growthfactor_integ(double a, void *param);
static void domain_walktoptree(int no);
//static peanokey peano_hilbert_key(int x, int y, int z, int bits);
static void peano_hilbert_key_inverse(peanokey key, int bits, int *x, int *y, int *z);
static void hydro_state_at_point(my_float pos[3], my_float vel[3], my_float *h_out, my_float *ngb_out,
    my_float *dhsml_out, my_float *rho_out, my_float *rhov_out, my_float *rhov2_out, my_float *rhoe_out);
static void make_it_active(int target);
};

static inline void show_buffer(void * buff, size_t size_in_byte)
{
    size_t *tmp= (size_t*)buff;
    void *hex;
    size_t raw;
    for (size_t i=0;i<=size_in_byte/sizeof(size_t);i++)
    {
        raw=tmp[i];
        hex=(void*)raw;
        gadgetmp2::DEBUG << "Pos " << &tmp[i] <<" : "<<  hex <<" : "<<  tmp[i] <<"\n";
        gadgetmp2::DEBUG.flush();
 }
 gadgetmp2::DEBUG <<"\n";
        gadgetmp2::DEBUG.flush();
}

#ifndef NOMPI
static inline my_float mpi_all_sum(my_float val)
{
  my_float retval=0;
  All_Reduce_buff *buff =  gadgetmp2::all_reduce_buff;
  size_t tt= gadgetmp2::ThisTask;
  size_t nt= gadgetmp2::NTask;
  size_t s= buff->get_size();

  buff->set_init(val,tt);
  MPI_Allgather(buff->get_buff_start(tt),s,MPI_BYTE, buff,s,MPI_BYTE,gadgetmp2::GADGET_WORLD);
  for(size_t i=0;i<nt;i++)
   retval += buff->read_re_init(i);
  return retval;
}

static inline my_float mpi_all_min(my_float val)
{
  my_float retval=val;
  All_Reduce_buff *buff =  gadgetmp2::all_reduce_buff;
  size_t tt= gadgetmp2::ThisTask;
  size_t nt= gadgetmp2::NTask;
  size_t s= buff->get_size();
  buff->set_init(val,tt);
  MPI_Allgather(buff->get_buff_start(tt),s,MPI_BYTE, buff,s,MPI_BYTE,gadgetmp2::GADGET_WORLD);
  for(size_t i=0;i<nt;i++)
  {
    if(buff->read_re_init(i)<retval)
    {
        retval = buff->read(i);
    }
  }
  return retval;
}

static inline my_float mpi_all_max(my_float val)
{
  my_float retval=val;
  All_Reduce_buff *buff =  gadgetmp2::all_reduce_buff;
  size_t tt= gadgetmp2::ThisTask;
  size_t nt= gadgetmp2::NTask;
  size_t s= buff->get_size();

  buff->set_init(val,tt);
  MPI_Allgather(buff->get_buff_start(tt),s,MPI_BYTE, buff,s,MPI_BYTE,gadgetmp2::GADGET_WORLD);
  for(size_t i=0;i<nt;i++)
   if(buff->read_re_init(i)>retval)
      retval = buff->read(i);
  return retval;
}

#endif

