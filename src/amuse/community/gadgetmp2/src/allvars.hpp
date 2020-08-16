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
/*! \file allvars.hpp
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should be added here, and declared as
 *  'extern'. The actual existence of these variables is provided by the file 'allvars.c'. To produce
 *  'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define's, typedef's, and enum's
 *     - add #include "allvars.h", delete the #ifndef ALLVARS_H conditional
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include "tags.hpp"

#include <pmpreal.h> //usage of debian packet libmpfrc++-dev
using namespace mpfr;
typedef mpreal  my_float;
typedef pmpreal  my_float_buff;

#define  GADGETVERSION   "2.0"   /*!< code version string */
#ifndef TIMESTEP_LIMITER
#define TIMESTEP_LIMITER
#endif // TIMESTEP_LIMITER

#define  TIMEBASE        (1<<28) /*!< The simulated timespan is mapped onto the integer interval [0,TIMESPAN],
*   where TIMESPAN needs to be a power of 2. Note that (1<<28) corresponds to 2^29
*/

#define  MAXTOPNODES     200000   /*!< Maximum number of nodes in the top-level tree used for domain decomposition */


typedef  long long  peanokey;    /*!< defines the variable type used for Peano-Hilbert keys */

#define  BITS_PER_DIMENSION 18	 /*!< Bits per dimension available for Peano-Hilbert order.
Note: If peanokey is defined as type int, the allowed maximum is 10.
If 64-bit integers are used, the maximum is 21 */

#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))  /*!< The number of different Peano-Hilbert cells */


#define  RNDTABLE         3000   /*!< gives the length of a table with random numbers, refreshed at every timestep.
This is used to allow application of random numbers to a specific particle
in a way that is independent of the number of processors used. */
#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37

#define  MAXLEN_FILENAME  4096    /*!< Maximum number of characters for filenames (including the full path) */

#ifdef   ISOTHERM_EQS
#define  GAMMA         (1.0)     /*!< index for isothermal gas */
#else
#define  GAMMA         (5.0/3)   /*!< adiabatic index of simulated gas */
#endif

#define  GAMMA_MINUS1  (GAMMA-1)

#define  HYDROGEN_MASSFRAC 0.76  /*!< mass fraction of hydrogen, relevant only for radiative cooling */

/* Some physical constants in cgs units */

#define  GRAVITY           6.672e-8   /*!< Gravitational constant (in cgs units) */
#define  SOLAR_MASS        1.989e33
#define  SOLAR_LUM         3.826e33
#define  RAD_CONST         7.565e-15
#define  AVOGADRO          6.0222e23
#define  BOLTZMANN         1.3806e-16
#define  GAS_CONST         8.31425e7
//#define  C                 2.9979e10
#define  PLANCK            6.6262e-27
#define  CM_PER_MPC        3.085678e24
#define  PROTONMASS        1.6726e-24
#define  ELECTRONMASS      9.10953e-28
#define  THOMPSON          6.65245e-25
#define  ELECTRONCHARGE    4.8032e-10
#define  HUBBLE            3.2407789e-18	/* in h/sec */

/* Some conversion factors */

#define  SEC_PER_MEGAYEAR  3.155e13
#define  SEC_PER_YEAR      3.155e7

#ifndef ASMTH
#define ASMTH 1.25  /*!< ASMTH gives the scale of the short-range/long-range force split in units of FFT-mesh cells */
#endif

#ifndef RCUT
#define RCUT  4.5   /*!< RCUT gives the maximum distance (in units of the scale used for the force split) out to
which short-range forces are evaluated in the short-range tree walk. */
#endif

#define MAX_NGB             20000  /*!< defines maximum length of neighbour list */

#define MAXLEN_OUTPUTLIST   500	   /*!< maxmimum number of entries in list of snapshot output times */

#define DRIFT_TABLE_LENGTH  1000   /*!< length of the lookup table used to hold the drift and kick factors */

#define MAXITER             150    /*!< maxmimum number of steps for SPH neighbour iteration */

//#ifdef DOUBLEPRECISION             /*!< If defined, the variable type FLOAT is set to "my_float", otherwise to FLOAT */
//#define FLOAT my_float
//#else
//#define FLOAT my_float
//#endif


#ifndef  TWODIMS
#define  NUMDIMS 3                                      /*!< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  2.546479089470                 /*!< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#define  KERNEL_COEFF_4  30.557749073644
#define  KERNEL_COEFF_5  5.092958178941
#define  KERNEL_COEFF_6  (-15.278874536822)
#define  NORM_COEFF      4.188790204786                 /*!< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */
#else
#define  NUMDIMS 2                                      /*!< For 2D-normalized kernel */
#define  KERNEL_COEFF_1  (5.0/7*2.546479089470)         /*!< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  (5.0/7*15.278874536822)
#define  KERNEL_COEFF_3  (5.0/7*45.836623610466)
#define  KERNEL_COEFF_4  (5.0/7*30.557749073644)
#define  KERNEL_COEFF_5  (5.0/7*5.092958178941)
#define  KERNEL_COEFF_6  (5.0/7*(-15.278874536822))
#define  NORM_COEFF      M_PI                           /*!< Coefficient for kernel normalization. */
#endif


    class DomainNODE
    {
    private:
        static size_t s0_off;
        static size_t s1_off;
        static size_t s2_off;
        static size_t vs0_off;
        static size_t vs1_off;
        static size_t vs2_off;
        static size_t mass_off;
    #ifdef UNEQUALSOFTENINGS
    #ifndef ADAPTIVE_GRAVSOFT_FORGAS
        static size_t bitflags_off;
    #else
        static size_t maxsoft_off;
    #endif
    #endif
        static size_t tot_size;
        static mpfr_prec_t prec;
    public:
//    my_float s[3];                     /*!< center-of-mass coordinates */
//    my_float vs[3];                    /*!< center-of-mass velocities */
//    my_float mass;                     /*!< mass of node */
//    #ifdef UNEQUALSOFTENINGS
//    #ifndef ADAPTIVE_GRAVSOFT_FORGAS
//    int   bitflags;                 /*!< this bit-field encodes the particle type with the largest softening among the particles of the nodes, and whether there are particles with different softening in the node */
//    #else
//    my_float maxsoft;                  /*!< hold the maximum gravitational softening of particles in the  node if the ADAPTIVE_GRAVSOFT_FORGAS option is selected */
//    #endif
//    #endif
        static inline size_t get_size()
        {
            return tot_size;
        };
        static inline size_t gen_size()
        {
            prec = my_float_buff::get_default_prec();
            size_t my_float_buff_size = my_float_buff::get_needed_mem_single(prec);
            size_t offset=0;
            s0_off = offset;
            offset += my_float_buff_size;
            s1_off = offset;
            offset += my_float_buff_size;
            s2_off = offset;
            offset += my_float_buff_size;
            vs0_off = offset;
            offset += my_float_buff_size;
            vs1_off = offset;
            offset += my_float_buff_size;
            vs2_off = offset;
            offset += my_float_buff_size;
            mass_off = offset;
            offset += my_float_buff_size;
    #ifdef UNEQUALSOFTENINGS
    #ifndef ADAPTIVE_GRAVSOFT_FORGAS
            bitflags_off = offset;
            offset += sizeof(int);
    #else
            maxsoft_off = offset;
            offset += my_float_buff_size;
    #endif
    #endif
            tot_size=offset;
            return tot_size;
        };
        inline void* get_buff_start(size_t index=0)
        {
            return (void*)((size_t)this + index * tot_size);
        }
        static inline void lcopy(DomainNODE* to, size_t tindex, DomainNODE* from, size_t findex)
        {
            void* l_from;
            void* l_to;
            l_from = (void*)((size_t)from + findex * tot_size);
            l_to = (void*)((size_t)to + tindex * tot_size);
            memcpy ( l_to, l_from, tot_size );
        };
        inline void set_init_s0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + s0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_s0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + s0_off);
            *to_store =value;
        };
        inline my_float read_re_init_s0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + s0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_s0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + s0_off);
            return *to_read;
        };
        inline void set_init_s1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + s1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_s1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + s1_off);
            *to_store =value;
        };
        inline my_float read_re_init_s1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + s1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_s1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + s1_off);
            return *to_read;
        };
        inline void set_init_s2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + s2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_s2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + s2_off);
            *to_store =value;
        };
        inline my_float read_re_init_s2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + s2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_s2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + s2_off);
            return *to_read;
        };
        inline void set_init_vs0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + vs0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_vs0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + vs0_off);
            *to_store =value;
        };
        inline my_float read_re_init_vs0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + vs0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_vs0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + vs0_off);
            return *to_read;
        };
        inline void set_init_vs1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + vs1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_vs1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + vs1_off);
            *to_store =value;
        };
        inline my_float read_re_init_vs1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + vs1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_vs1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + vs1_off);
            return *to_read;
        };
        inline void set_init_vs2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + vs2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_vs2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + vs2_off);
            *to_store =value;
        };
        inline my_float read_re_init_vs2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + vs2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_vs2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + vs2_off);
            return *to_read;
        };
        inline void set_init_mass(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + mass_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_mass(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + mass_off);
            *to_store =value;
        };
        inline my_float read_re_init_mass(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + mass_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_mass(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + mass_off);
            return *to_read;
        };
    #ifdef UNEQUALSOFTENINGS
    #ifndef ADAPTIVE_GRAVSOFT_FORGAS
        inline void set_bitflags(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + bitflags_off);
            *to_store =value;
        };
        inline int read_bitflags(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + bitflags_off);
            return *to_read;
        };
    #else
        inline void set_init_maxsoft(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + maxsoft_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_maxsoft(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + maxsoft_off);
            *to_store =value;
        };
        inline my_float read_re_init_maxsoft(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + maxsoft_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_maxsoft(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + maxsoft_off);
            return *to_read;
        };
    #endif
    #endif
    }
    ;

struct topnode_data
{
    int Daughter;                   /*!< index of first daughter cell (out of 8) of top-level node */
    int Pstart;                     /*!< for the present top-level node, this gives the index of the first node in the concatenated list of topnodes collected from all processors */
    int Blocks;                     /*!< for the present top-level node, this gives the number of corresponding nodes in the concatenated list of topnodes collected from all processors */
    int Leaf;                       /*!< if the node is a leaf, this gives its number when all leaves are traversed in Peano-Hilbert order */
    peanokey Size;                  /*!< number of Peano-Hilbert mesh-cells represented by top-level node */
    peanokey StartKey;              /*!< first Peano-Hilbert key in top-level node */
    long long Count;                /*!< counts the number of particles in this top-level node */
}
;

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
struct global_data_all_processes
{
    long long TotNumPart;		/*!< total particle numbers (global value) */
    long long TotN_gas;		/*!< total gas particle number (global value) */

    int MaxPart;			/*!< This gives the maxmimum number of particles that can be stored on one processor. */
    int MaxPartSph;		/*!< This gives the maxmimum number of SPH particles that can be stored on one processor. */

    double BoxSize;               /*!< Boxsize in case periodic boundary conditions are used */

    int ICFormat;			/*!< selects different versions of IC file-format */

    int SnapFormat;		/*!< selects different versions of snapshot file-formats */

    int NumFilesPerSnapshot;      /*!< number of files in multi-file snapshot dumps */
    int NumFilesWrittenInParallel;/*!< maximum number of files that may be written simultaneously when
    writing/reading restart-files, or when writing snapshot files */

    int BufferSize;		/*!< size of communication buffer in MB */
    int BunchSizeForce;		/*!< number of particles fitting into the buffer in the parallel tree-force algorithm  */
    int BunchSizeDensity;         /*!< number of particles fitting into the communication buffer in the density computation */
    int BunchSizeHydro;           /*!< number of particles fitting into the communication buffer in the SPH hydrodynamical force computation */
    int BunchSizeDomain;          /*!< number of particles fitting into the communication buffer in the domain decomposition */

    double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually
    NOT be balanced.  Each processor allocates memory for PartAllocFactor times
    the average number of particles to allow for that */

    double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
    the maximum(!) number of particles.  Note: A typical local tree for N
    particles needs usually about ~0.65*N nodes. */

    /* some SPH parameters */

    double DesNumNgb;             /*!< Desired number of SPH neighbours */
    double MaxNumNgbDeviation;    /*!< Maximum allowed deviation neighbour number */

    double ArtBulkViscConst;      /*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
    double ArtBulkViscBeta;      /*!< Sets the parameter \f$\beta\f$ of the artificial viscosity */
    double InitGasTemp;		/*!< may be used to set the temperature in the IC's */
    double MinGasTemp;		/*!< may be used to set a floor for the gas temperature */
    double MinEgySpec;            /*!< the minimum allowed temperature expressed as energy per unit mass */


    /* some force counters  */

    long long TotNumOfForces;	             /*!< counts total number of force computations  */
    long long NumForcesSinceLastDomainDecomp;  /*!< count particle updates since last domain decomposition */


    /* system of units  */

    double G;                        /*!< Gravity-constant in internal units */
    double UnitTime_in_s;   	   /*!< factor to convert internal time unit to seconds/h */
    double UnitMass_in_g;            /*!< factor to convert internal mass unit to grams/h */
    double UnitVelocity_in_cm_per_s; /*!< factor to convert intqernal velocity unit to cm/sec */
    double UnitLength_in_cm;         /*!< factor to convert internal length unit to cm/h */
    double UnitPressure_in_cgs;      /*!< factor to convert internal pressure unit to cgs units (little 'h' still around!) */
    double UnitDensity_in_cgs;       /*!< factor to convert internal length unit to g/cm^3*h^2 */
    double UnitCoolingRate_in_cgs;   /*!< factor to convert internal cooling rate to cgs units */
    double UnitEnergy_in_cgs;        /*!< factor to convert internal energy to cgs units */
    double UnitTime_in_Megayears;    /*!< factor to convert internal time to megayears/h */
    double GravityConstantInternal;  /*!< If set to zero in the parameterfile, the internal value of the
    gravitational constant is set to the Newtonian value based on the system of
    units specified. Otherwise the value provided is taken as internal gravity constant G. */


    /* Cosmological parameters */

    double Hubble;       /*!< Hubble-constant in internal units */
    double Omega0;       /*!< matter density in units of the critical density (at z=0)*/
    double OmegaLambda;  /*!< vaccum energy density relative to crictical density (at z=0) */
    double OmegaBaryon;  /*!< baryon density in units of the critical density (at z=0)*/
    double HubbleParam;  /*!< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute physical values for cooling physics */


    /* Code options */

    int ComovingIntegrationOn;	/*!< flags that comoving integration is enabled */
    int PeriodicBoundariesOn;     /*!< flags that periodic boundaries are enabled */
    int ResubmitOn;               /*!< flags that automatic resubmission of job to queue system is enabled */
    int TypeOfOpeningCriterion;   /*!< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative criterion */
    int TypeOfTimestepCriterion;  /*!< gives type of timestep criterion (only 0 supported right now - unlike gadget-1.1) */
    int OutputListOn;             /*!< flags that output times are listed in a specified file */


    /* Parameters determining output frequency */

    int SnapshotFileCount;        /*!< number of snapshot that is written next */
    double TimeBetSnapshot;       /*!< simulation time interval between snapshot files */
    double TimeOfFirstSnapshot;   /*!< simulation time of first snapshot files */
    double CpuTimeBetRestartFile; /*!< cpu-time between regularly generated restart files */
    double TimeLastRestartFile;   /*!< cpu-time when last restart-file was written */
    double TimeBetStatistics;     /*!< simulation time interval between computations of energy statistics */
    double TimeLastStatistics;    /*!< simulation time when the energy statistics was computed the last time */
    int NumCurrentTiStep;         /*!< counts the number of system steps taken up to this point */


    /* Current time of the simulation, global step, and end of simulation */

    double Time;                  /*!< current time of the simulation */
    double TimeBegin;             /*!< time of initial conditions of the simulation */
    double TimeStep;              /*!< difference between current times of previous and current timestep */
    double TimeMax;	        /*!< marks the point of time until the simulation is to be evolved */


    /* variables for organizing discrete timeline */

    double Timebase_interval;     /*!< factor to convert from floating point time interval to integer timeline */
    int Ti_Current;               /*!< current time on integer timeline */
    int Ti_nextoutput;            /*!< next output time on integer timeline */
    #ifdef FLEXSTEPS
    int PresentMinStep;           /*!< If FLEXSTEPS is used, particle timesteps are chosen as multiples of the present minimum timestep. */
    int PresentMaxStep;		/*!< If FLEXSTEPS is used, this is the maximum timestep in timeline units, rounded down to the next power 2 division */
    #endif

    /* Variables that keep track of cumulative CPU consumption */

    double TimeLimitCPU;          /*!< CPU time limit as defined in parameterfile */
    double CPU_TreeConstruction;  /*!< time spent for constructing the gravitational tree */
    double CPU_TreeWalk;          /*!< actual time spent for pure tree-walks */
    double CPU_Gravity;           /*!< cumulative time used for gravity computation (tree-algorithm only) */
    double CPU_Potential;         /*!< time used for computing gravitational potentials */
    double CPU_Domain;            /*!< cumulative time spent for domain decomposition */
    double CPU_Snapshot;          /*!< time used for writing snapshot files */
    double CPU_Total;             /*!< cumulative time spent for domain decomposition */
    double CPU_CommSum;           /*!< accumulated time used for communication, and for collecting partial results, in tree-gravity */
    double CPU_Imbalance;         /*!< cumulative time lost accross all processors as work-load imbalance in gravitational tree */
    double CPU_HydCompWalk;       /*!< time used for actual SPH computations, including neighbour search */
    double CPU_HydCommSumm;       /*!< cumulative time used for communication in SPH, and for collecting partial results */
    double CPU_HydImbalance;      /*!< cumulative time lost due to work-load imbalance in SPH */
    double CPU_Hydro;             /*!< cumulative time spent for SPH related computations */
    double CPU_EnsureNgb;         /*!< time needed to iterate on correct neighbour numbers */
    double CPU_Predict;           /*!< cumulative time to drift the system forward in time, including dynamic tree updates */
    double CPU_TimeLine;          /*!< time used for determining new timesteps, and for organizing the timestepping, including kicks of active particles */
    double CPU_PM;                /*!< time used for long-range gravitational force */
    double CPU_Peano;             /*!< time required to establish Peano-Hilbert order */

    /* tree code opening criterion */

    double ErrTolTheta;		/*!< BH tree opening angle */
    double ErrTolForceAcc;	/*!< parameter for relative opening criterion in tree walk */


    /* adjusts accuracy of time-integration */

    double ErrTolIntAccuracy;	/*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
    timestep is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

    double MinSizeTimestep;       /*!< minimum allowed timestep. Normally, the simulation terminates if the
    timestep determined by the timestep criteria falls below this limit. */
    double MaxSizeTimestep;       /*!< maximum allowed timestep */

    double MaxRMSDisplacementFac; /*!< this determines a global timestep criterion for cosmological simulations
    in comoving coordinates.  To this end, the code computes the rms velocity
    of all particles, and limits the timestep such that the rms displacement
    is a fraction of the mean particle separation (determined from the
    particle mass and the cosmological parameters). This parameter specifies
    this fraction. */

    double CourantFac;		/*!< SPH-Courant factor */

    int BunchSizeTime;		/*!< number of particles fitting into the communication buffer in the timestep communication */


    /* frequency of tree reconstruction/domain decomposition */

    double TreeDomainUpdateFrequency; /*!< controls frequency of domain decompositions  */


    /* Gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening length).
     * Five groups of particles are supported 0="gas", 1="halo", 2="disk", 3="bulge", 4="stars", 5="bndry"
     */

    double MinGasHsmlFractional;  /*!< minimum allowed SPH smoothing length in units of SPH gravitational softening length */
    double MinGasHsml;            /*!< minimum allowed SPH smoothing length */


    double SofteningGas;          /*!< comoving gravitational softening lengths for type 0 */
    double SofteningHalo;         /*!< comoving gravitational softening lengths for type 1 */
    double SofteningDisk;         /*!< comoving gravitational softening lengths for type 2 */
    double SofteningBulge;        /*!< comoving gravitational softening lengths for type 3 */
    double SofteningStars;        /*!< comoving gravitational softening lengths for type 4 */
    double SofteningBndry;        /*!< comoving gravitational softening lengths for type 5 */

    double SofteningGasMaxPhys;   /*!< maximum physical softening length for type 0 */
    double SofteningHaloMaxPhys;  /*!< maximum physical softening length for type 1 */
    double SofteningDiskMaxPhys;  /*!< maximum physical softening length for type 2 */
    double SofteningBulgeMaxPhys; /*!< maximum physical softening length for type 3 */
    double SofteningStarsMaxPhys; /*!< maximum physical softening length for type 4 */
    double SofteningBndryMaxPhys; /*!< maximum physical softening length for type 5 */

    double SofteningTable[6];     /*!< current (comoving) gravitational softening lengths for each particle type */
    double ForceSoftening[6];     /*!< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */


    double MassTable[6];          /*!< Table with particle masses for particle types with equal mass.
    If particle masses are all equal for one type, the corresponding entry in MassTable
    is set to this value, allowing the size of the snapshot files to be reduced. */



    /* some filenames */

    char InitCondFile[MAXLEN_FILENAME];          /*!< filename of initial conditions */
    char OutputDir[MAXLEN_FILENAME];             /*!< output directory of the code */
    char SnapshotFileBase[MAXLEN_FILENAME];      /*!< basename to construct the names of snapshotf files */
    char EnergyFile[MAXLEN_FILENAME];            /*!< name of file with energy statistics */
    char CpuFile[MAXLEN_FILENAME];               /*!< name of file with cpu-time statistics */
    char InfoFile[MAXLEN_FILENAME];              /*!< name of log-file with a list of the timesteps taken */
    char TimingsFile[MAXLEN_FILENAME];           /*!< name of file with performance metrics of gravitational tree algorithm */
    char RestartFile[MAXLEN_FILENAME];           /*!< basename of restart-files */
    char ResubmitCommand[MAXLEN_FILENAME];       /*!< name of script-file that will be executed for automatic restart */
    char OutputListFilename[MAXLEN_FILENAME];    /*!< name of file with list of desired output times */

    double OutputListTimes[MAXLEN_OUTPUTLIST];   /*!< table with desired output times */
    int OutputListLength;                        /*!< number of output times stored in the table of desired output times */

    }
    ;


    /*! This structure holds all the information that is
     * stored for each particle of the simulation.
     */
    class particle_data
    {
        public:
        my_float Pos[3];			/*!< particle position at its current time */
        my_float Mass;			/*!< particle mass */
        my_float Vel[3];			/*!< particle velocity at its current time */
        my_float radius;
        my_float GravAccel[3];		/*!< particle acceleration due to gravity */

        #ifdef FORCETEST
        my_float GravAccelDirect[3];	/*!< particle acceleration when computed with direct summation */
        #endif
        my_float Potential;		/*!< gravitational potential */
        my_float OldAcc;			/*!< magnitude of old gravitational force. Used in relative opening criterion */
        #ifndef LONGIDS
        unsigned int ID;		/*!< particle identifier */
        #else
        unsigned long long ID;        /*!< particle identifier */
        #endif

        int Type;		        /*!< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
        int Ti_endstep;               /*!< marks start of current timestep of particle on integer timeline */
        int Ti_begstep;               /*!< marks end of current timestep of particle on integer timeline */
        #ifdef TIMESTEP_LIMITER
        int Ti_sizestep;
        #endif
        #ifdef FLEXSTEPS
        int FlexStepGrp;		/*!< a random 'offset' on the timeline to create a smooth groouping of particles */
        #endif
        my_float GravCost;		/*!< weight factor used for balancing the work-load */
        #ifdef PSEUDOSYMMETRIC
        my_float AphysOld;               /*!< magnitude of acceleration in last timestep. Used to make a first order
        prediction of the change of acceleration expected in the future, thereby
        allowing to guess whether a decrease/increase of the timestep should occur
        in the timestep that is started. */
        #endif
        static inline size_t get_size_()
        {
            size_t retval=0;
            size_t my_float_size = my_float_buff::get_needed_mem_single();
            retval += my_float_size*3; //Pos[3];
            retval += my_float_size; //Mass;
            retval += my_float_size*3; //Vel[3];
            retval += my_float_size; //radius;
            retval += my_float_size*3; //GravAccel[3];
            #ifdef FORCETEST
            retval += my_float_size; //GravAccelDirect[3];
            #endif
            retval += my_float_size; //Potential;
            retval += my_float_size; //OldAcc;
            #ifndef LONGIDS
            retval += sizeof(unsigned int); //ID;
            #else
            retval += sizeof(unsigned long long); //ID;
            #endif
            retval += sizeof(int); //Type;
            retval += sizeof(int); //Ti_endstep;
            retval += sizeof(int); //Ti_begstep;
            #ifdef TIMESTEP_LIMITER
            retval += sizeof(int); //Ti_sizestep;
            #endif
            #ifdef FLEXSTEPS
            retval += sizeof(int); //FlexStepGrp;
            #endif
            retval += my_float_size; //GravCost;
            #ifdef PSEUDOSYMMETRIC
            retval += my_float_size; //AphysOld;
            #endif
            return retval;
        };

        inline void change_prec ()
        {
            mpfr_prec_t prec = my_float::get_default_prec();
            Pos[0].setPrecision(prec);
            Pos[1].setPrecision(prec);
            Pos[2].setPrecision(prec);
            Mass.setPrecision(prec);
            Vel[0].setPrecision(prec);
            Vel[1].setPrecision(prec);
            Vel[2].setPrecision(prec);
            radius.setPrecision(prec);
            GravAccel[0].setPrecision(prec);
            GravAccel[1].setPrecision(prec);
            GravAccel[2].setPrecision(prec);
            #ifdef FORCETEST
            GravAccelDirect[3].setPrecision(prec);
            #endif
            Potential.setPrecision(prec);
            OldAcc.setPrecision(prec);
            GravCost.setPrecision(prec);
            #ifdef PSEUDOSYMMETRIC
            AphysOld.setPrecision(prec);
            #endif
        };
    }
    ;

    class particle_data_buff
    {
    private:
        static size_t Pos0_off;
        static size_t Pos1_off;
        static size_t Pos2_off;
        static size_t Mass_off;
        static size_t Vel0_off;
        static size_t Vel1_off;
        static size_t Vel2_off;
        static size_t Radius_off;
        static size_t GravAccel0_off;
        static size_t GravAccel1_off;
        static size_t GravAccel2_off;
        #ifdef FORCETEST
        static size_t GravAccelDirect0_off;
        static size_t GravAccelDirect1_off;
        static size_t GravAccelDirect2_off;
        #endif
        static size_t Potential_off;
        static size_t OldAcc_off;
        static size_t ID_off;
        static size_t Type_off;
        static size_t Ti_endstep_off;
        static size_t Ti_begstep_off;
        #ifdef TIMESTEP_LIMITER
        static size_t Ti_sizestep_off;
        #endif
        #ifdef FLEXSTEPS
        static size_t FlexStepGrp_off;
        #endif
        static size_t GravCost_off;
        #ifdef PSEUDOSYMMETRIC
        static size_t AphysOld_off;
        #endif
        static size_t tot_size;
        static mpfr_prec_t prec;

    public:
        //my_float Pos[3];
        //my_float Mass;
        //my_float Vel[3];
        //my_float GravAccel[3];
        //#ifdef FORCETEST
        //my_float GravAccelDirect[3];
        //#endif
        //my_float Potential;
        //my_float OldAcc;
        //#ifndef LONGIDS
        //unsigned int ID;
        //#else
        //unsigned long long ID;
        //#endif
        //int Type;
        //int Ti_endstep;
        //int Ti_begstep;
        //#ifdef TIMESTEP_LIMITER
        //int Ti_sizestep;
        //#endif
        //#ifdef FLEXSTEPS
        //int FlexStepGrp;
        //#endif
        //my_float GravCost;
        //#ifdef PSEUDOSYMMETRIC
        //my_float AphysOld;
        //#endif

        static inline size_t get_size()
        {
            return tot_size;
        };
        static inline size_t gen_size()
        {
            prec = my_float_buff::get_default_prec();
            size_t my_float_buff_size = my_float_buff::get_needed_mem_single(prec);
            size_t offset=0;
            Pos0_off = offset;
            offset += my_float_buff_size;
            Pos1_off = offset;
            offset += my_float_buff_size;
            Pos2_off = offset;
            offset += my_float_buff_size;
            Mass_off = offset;
            offset += my_float_buff_size;
            Vel0_off = offset;
            offset += my_float_buff_size;
            Vel1_off = offset;
            offset += my_float_buff_size;
            Vel2_off = offset;
            offset += my_float_buff_size;
            Radius_off = offset;
            offset += my_float_buff_size;
            GravAccel0_off = offset;
            offset += my_float_buff_size;
            GravAccel1_off = offset;
            offset += my_float_buff_size;
            GravAccel2_off = offset;
            offset += my_float_buff_size;
            #ifdef FORCETEST
            GravAccelDirect0_off = offset;
            offset += my_float_buff_size;
            GravAccelDirect1_off = offset;
            offset += my_float_buff_size;
            GravAccelDirect2_off = offset;
            offset += my_float_buff_size;
            #endif
            Potential_off = offset;
            offset += my_float_buff_size;
            OldAcc_off = offset;
            offset += my_float_buff_size;
            #ifndef LONGIDS
            ID_off = offset;
            offset += sizeof(unsigned int);
            #else
            ID_off = offset;
            offset += sizeof(unsigned long long);
            #endif
            Type_off = offset;
            offset += sizeof(int);
            Ti_endstep_off = offset;
            offset += sizeof(int);
            Ti_begstep_off = offset;
            offset += sizeof(int);
            #ifdef TIMESTEP_LIMITER
            Ti_sizestep_off = offset;
            offset += sizeof(int);
            #endif
            #ifdef FLEXSTEPS
            FlexStepGrp_off = offset;
            offset += sizeof(int);
            #endif
            GravCost_off = offset;
            offset += my_float_buff_size;
            #ifdef PSEUDOSYMMETRIC
            AphysOld_off = offset;
            offset += my_float_buff_size;
            #endif
            tot_size=offset;
            return tot_size;
        };
        static inline void l_copy(particle_data_buff* to, size_t tindex, particle_data_buff* from, size_t findex)
        {
            void* l_from;
            void* l_to;
            l_from = (void*)((size_t)from + findex * tot_size);
            l_to = (void*)((size_t)to + tindex * tot_size);
            memcpy ( l_to, l_from, tot_size );
        };
        inline void set_init_Pos0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Pos0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Pos0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Pos0_off);
            *to_store =value;
        };
        inline my_float read_re_init_Pos0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Pos0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos0_off);
            return *to_read;
        };
        inline void set_init_Pos1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Pos1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Pos1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Pos1_off);
            *to_store =value;
        };
        inline my_float read_re_init_Pos1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Pos1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos1_off);
            return *to_read;
        };
        inline void set_init_Pos2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Pos2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Pos2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Pos2_off);
            *to_store =value;
        };
        inline my_float read_re_init_Pos2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Pos2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos2_off);
            return *to_read;
        };
        inline void set_init_Mass(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Mass_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Mass(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Mass_off);
            *to_store =value;
        };
        inline my_float read_re_init_Mass(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Mass_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Mass(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Mass_off);
            return *to_read;
        };
        inline void set_init_Vel0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Vel0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Vel0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Vel0_off);
            *to_store =value;
        };
        inline my_float read_re_init_Vel0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Vel0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel0_off);
            return *to_read;
        };
        inline void set_init_Vel1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Vel1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Vel1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Vel1_off);
            *to_store =value;
        };
        inline my_float read_re_init_Vel1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Vel1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel1_off);
            return *to_read;
        };
        inline void set_init_Vel2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Vel2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Vel2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Vel2_off);
            *to_store =value;
        };
        inline my_float read_re_init_Vel2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Vel2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel2_off);
            return *to_read;
        };
        inline void set_init_Radius(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Radius_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Radius(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Radius_off);
            *to_store =value;
        };
        inline my_float read_re_init_Radius(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Radius_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Radius(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Radius_off);
            return *to_read;
        };
        inline void set_init_GravAccel0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + GravAccel0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_GravAccel0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + GravAccel0_off);
            *to_store =value;
        };
        inline my_float read_re_init_GravAccel0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + GravAccel0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_GravAccel0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + GravAccel0_off);
            return *to_read;
        };
        inline void set_init_GravAccel1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + GravAccel1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_GravAccel1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + GravAccel1_off);
            *to_store =value;
        };
        inline my_float read_re_init_GravAccel1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + GravAccel1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_GravAccel1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + GravAccel1_off);
            return *to_read;
        };
        inline void set_init_GravAccel2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + GravAccel2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_GravAccel2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + GravAccel2_off);
            *to_store =value;
        };
        inline my_float read_re_init_GravAccel2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + GravAccel2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_GravAccel2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + GravAccel2_off);
            return *to_read;
        };
        #ifdef FORCETEST
        inline void set_init_GravAccelDirect0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + GravAccelDirect0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_GravAccelDirect0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + GravAccelDirect0_off);
            *to_store =value;
        };
        inline my_float read_re_init_GravAccelDirect0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + GravAccelDirect0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_GravAccelDirect0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + GravAccelDirect0_off);
            return *to_read;
        };
        inline void set_init_GravAccelDirect1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + GravAccelDirect1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_GravAccelDirect1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + GravAccelDirect1_off);
            *to_store =value;
        };
        inline my_float read_re_init_GravAccelDirect1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + GravAccelDirect1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_GravAccelDirect1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + GravAccelDirect1_off);
            return *to_read;
        };
        inline void set_init_GravAccelDirect2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + GravAccelDirect2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_GravAccelDirect2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + GravAccelDirect2_off);
            *to_store =value;
        };
        inline my_float read_re_init_GravAccelDirect2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + GravAccelDirect2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_GravAccelDirect2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + GravAccelDirect2_off);
            return *to_read;
        };
        #endif
        inline void set_init_Potential(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Potential_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Potential(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Potential_off);
            *to_store =value;
        };
        inline my_float read_re_init_Potential(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Potential_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Potential(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Potential_off);
            return *to_read;
        };
        inline void set_init_OldAcc(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + OldAcc_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_OldAcc(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + OldAcc_off);
            *to_store =value;
        };
        inline my_float read_re_init_OldAcc(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + OldAcc_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_OldAcc(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + OldAcc_off);
            return *to_read;
        };
        #ifndef LONGIDS
        inline void set_ID(unsigned int value, size_t index=0)
        {
            unsigned int* to_store;
            to_store = (unsigned int*)((size_t)this + index * tot_size + ID_off);
            *to_store =value;
        };
        inline unsigned int read_ID(size_t index=0)
        {
            unsigned int* to_read;
            to_read = (unsigned int*)((size_t)this + index * tot_size + ID_off);
            return *to_read;
        };
        #else
        inline void set_ID(unsigned long long value, size_t index=0)
        {
            unsigned long long* to_store;
            to_store = (unsigned long long*)((size_t)this + index * tot_size + ID_off);
            *to_store =value;
        };
        inline unsigned long long read_ID(size_t index=0)
        {
            unsigned long long* to_read;
            to_read = (unsigned long long*)((size_t)this + index * tot_size + ID_off);
            return *to_read;
        };
        #endif
        inline void set_Type(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + Type_off);
            *to_store =value;
        };
        inline int read_Type(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + Type_off);
            return *to_read;
        };
        inline void set_Ti_endstep(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + Ti_endstep_off);
            *to_store =value;
        };
        inline int read_Ti_endstep(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + Ti_endstep_off);
            return *to_read;
        };
        inline void set_Ti_begstep(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + Ti_begstep_off);
            *to_store =value;
        };
        inline int read_Ti_begstep(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + Ti_begstep_off);
            return *to_read;
        };
        #ifdef TIMESTEP_LIMITER
        inline void set_Ti_sizestep(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + Ti_sizestep_off);
            *to_store =value;
        };
        inline int read_Ti_sizestep(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + Ti_sizestep_off);
            return *to_read;
        };
        #endif
        #ifdef FLEXSTEPS
        inline void set_FlexStepGrp(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + FlexStepGrp_off);
            *to_store =value;
        };
        inline int read_FlexStepGrp(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + FlexStepGrp_off);
            return *to_read;
        };
        #endif
        inline void set_init_GravCost(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + GravCost_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_GravCost(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + GravCost_off);
            *to_store =value;
        };
        inline my_float read_re_init_GravCost(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + GravCost_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_GravCost(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + GravCost_off);
            return *to_read;
        };
        #ifdef PSEUDOSYMMETRIC
        inline void set_init_AphysOld(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + AphysOld_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_AphysOld(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + AphysOld_off);
            *to_store =value;
        };
        inline my_float read_re_init_AphysOld(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + AphysOld_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_AphysOld(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + AphysOld_off);
            return *to_read;
        };
        #endif
        inline void l_fill(struct particle_data in, size_t index=0)
        {
        set_init_Pos0(in.Pos[0],index);
        set_init_Pos1(in.Pos[1],index);
        set_init_Pos2(in.Pos[2],index);
        set_init_Mass(in.Mass,index);
        set_init_Vel0(in.Vel[0],index);
        set_init_Vel1(in.Vel[1],index);
        set_init_Vel2(in.Vel[2],index);
        set_init_Radius(in.radius,index);
        set_init_GravAccel0(in.GravAccel[0],index);
        set_init_GravAccel1(in.GravAccel[1],index);
        set_init_GravAccel2(in.GravAccel[2],index);
        #ifdef FORCETEST
        set_init_GravAccelDirect0(in.GravAccelDirect[0],index);
        set_init_GravAccelDirect1(in.GravAccelDirect[1],index);
        set_init_GravAccelDirect2(in.GravAccelDirect[2],index);
        #endif
        set_init_Potential(in.Potential,index);
        set_init_OldAcc(in.OldAcc,index);
        set_ID(in.ID,index);
        set_Type(in.Type,index);
        set_Ti_endstep(in.Ti_endstep,index);
        set_Ti_begstep(in.Ti_begstep,index);
        #ifdef TIMESTEP_LIMITER
        set_Ti_sizestep(in.Ti_sizestep,index);
        #endif
        #ifdef FLEXSTEPS
        set_FlexStepGrp(in.FlexStepGrp,index);
        #endif
        set_init_GravCost(in.GravCost,index);
        #ifdef PSEUDOSYMMETRIC
        set_init_AphysOld(in.AphysOld,index);
        #endif
        }
        inline particle_data unpack_particle(size_t index=0)
        {
            particle_data retval;
        retval.Pos[0] = read_re_init_Pos0(index);
        retval.Pos[1] = read_re_init_Pos1(index);
        retval.Pos[2] = read_re_init_Pos2(index);
        retval.Mass = read_re_init_Mass(index);
        retval.Vel[0] = read_re_init_Vel0(index);
        retval.Vel[1] = read_re_init_Vel1(index);
        retval.Vel[2] = read_re_init_Vel2(index);
        retval.radius = read_re_init_Radius(index);
        retval.GravAccel[0] = read_re_init_GravAccel0(index);
        retval.GravAccel[1] = read_re_init_GravAccel1(index);
        retval.GravAccel[2] = read_re_init_GravAccel2(index);
        #ifdef FORCETEST
        retval.GravAccelDirect[0] = read_re_init_GravAccelDirect0(index);
        retval.GravAccelDirect[1] = read_re_init_GravAccelDirect1(index);
        retval.GravAccelDirect[2] = read_re_init_GravAccelDirect2(index);
        #endif
        retval.Potential = read_re_init_Potential(index);
        retval.OldAcc = read_re_init_OldAcc(index);
        retval.ID = read_ID(index);
        retval.Type=  read_Type(index);
        retval.Ti_endstep = read_Ti_endstep(index);
        retval.Ti_begstep = read_Ti_begstep(index);
        #ifdef TIMESTEP_LIMITER
        retval.Ti_sizestep = read_Ti_sizestep(index);
        #endif
        #ifdef FLEXSTEPS
        retval.FlexStepGrp = read_FlexStepGrp(index);
        #endif
        retval.GravCost = read_re_init_GravCost(index);
        #ifdef PSEUDOSYMMETRIC
        retval.AphysOld = read_re_init_AphysOld(index);
        #endif
        return retval;
        }
    }
    ;

    /* the following struture holds data that is stored for each SPH particle in addition to the collisionless
     * variables.
     */
    struct sph_particle_data
    {
        my_float Entropy;                /*!< current value of entropy (actually entropic function) of particle */
        my_float Density;		/*!< current baryonic mass density of particle */
        my_float Hsml;			/*!< current smoothing length */
        my_float Left;                   /*!< lower bound in iterative smoothing length search */
        my_float Right;                  /*!< upper bound in iterative smoothing length search */
        my_float NumNgb;                 /*!< weighted number of neighbours found */
        my_float Pressure;		/*!< current pressure */
        my_float DtEntropy;              /*!< rate of change of entropy */
        my_float HydroAccel[3];		/*!< acceleration due to hydrodynamical force */
        my_float VelPred[3];		/*!< predicted SPH particle velocity at the current time */
        my_float DivVel;			/*!< local velocity divergence */
        my_float CurlVel;		/*!< local velocity curl */
        my_float Rot[3];		        /*!< local velocity curl */
        my_float DhsmlDensityFactor;     /*!< correction factor needed in the equation of motion of the conservative entropy formulation of SPH */
        my_float MaxSignalVel;           /*!< maximum "signal velocity" occuring for this particle */
        #ifdef TIMESTEP_UPDATE
        int   FeedbackFlag;
        my_float FeedAccel[3];  /*!< acceleration due to feedback force */
        #endif
        #ifdef MORRIS97VISC
        my_float Alpha;		        /*!< viscosity coefficient */
        my_float DAlphaDt;       		/*!< time rate of change of viscosity coefficient */
        #endif
        static inline size_t get_size_()
        {
            size_t retval=0;
            size_t my_float_size = my_float_buff::get_needed_mem_single();
            retval += my_float_size; //Entropy
            retval += my_float_size; //Density
            retval += my_float_size; //Hsml
            retval += my_float_size; //Left
            retval += my_float_size; //Right
            retval += my_float_size; //NumNgb
            retval += my_float_size; //Pressure
            retval += my_float_size; //DtEntropy
            retval += my_float_size*3; //HydroAccel[3]
            retval += my_float_size*3; //VelPred[3]
            retval += my_float_size; //DivVel
            retval += my_float_size; //CurlVel
            retval += my_float_size*3; //Rot[3]
            retval += my_float_size; //DhsmlDensityFactor
            retval += my_float_size; //MaxSignalVel
            #ifdef TIMESTEP_UPDATE
            retval += sizeof(int); //FeedbackFlag;
            retval += my_float_size*3; //FeedAccel[3]
            #endif
            #ifdef MORRIS97VISC
            retval += my_float_size; //Alpha;
            retval += my_float_size; //DAlphaDt;
            #endif
            return retval;
        };

        inline void change_prec ()
        {
            mpfr_prec_t prec = my_float::get_default_prec();
            Entropy.setPrecision(prec);
            Density.setPrecision(prec);
            Hsml.setPrecision(prec);
            Left.setPrecision(prec);
            Right.setPrecision(prec);
            NumNgb.setPrecision(prec);
            Pressure.setPrecision(prec);
            DtEntropy.setPrecision(prec);
            HydroAccel[0].setPrecision(prec);
            HydroAccel[1].setPrecision(prec);
            HydroAccel[2].setPrecision(prec);
            VelPred[0].setPrecision(prec);
            VelPred[1].setPrecision(prec);
            VelPred[2].setPrecision(prec);
            DivVel.setPrecision(prec);
            CurlVel.setPrecision(prec);
            Rot[0].setPrecision(prec);
            Rot[1].setPrecision(prec);
            Rot[2].setPrecision(prec);
            DhsmlDensityFactor.setPrecision(prec);
            MaxSignalVel.setPrecision(prec);
            #ifdef TIMESTEP_UPDATE
            FeedAccel[3].setPrecision(prec);
            #endif
            #ifdef MORRIS97VISC
            Alpha.setPrecision(prec);
            DAlphaDt.setPrecision(prec);
            #endif
        };
    }
    ;

    class sph_particle_data_buff
    {
    private:
        static size_t Entropy_off;
        static size_t Density_off;
        static size_t Hsml_off;
        static size_t Left_off;
        static size_t Right_off;
        static size_t NumNgb_off;
        static size_t Pressure_off;
        static size_t DtEntropy_off;
        static size_t HydroAccel0_off;
        static size_t HydroAccel1_off;
        static size_t HydroAccel2_off;
        static size_t VelPred0_off;
        static size_t VelPred1_off;
        static size_t VelPred2_off;
        static size_t DivVel_off;
        static size_t CurlVel_off;
        static size_t Rot0_off;
        static size_t Rot1_off;
        static size_t Rot2_off;
        static size_t DhsmlDensityFactor_off;
        static size_t MaxSignalVel_off;
        #ifdef TIMESTEP_UPDATE
        static size_t FeedbackFlag_off;
        static size_t FeedAccel0_off;
        static size_t FeedAccel1_off;
        static size_t FeedAccel2_off;
        #endif
        #ifdef MORRIS97VISC
        static size_t Alpha_off;
        static size_t DAlphaDt_off;
        #endif
        static size_t tot_size;
        static mpfr_prec_t prec;

    public:
        static inline size_t get_size()
        {
            return tot_size;
        };
        static inline size_t gen_size()
        {
            prec = my_float_buff::get_default_prec();
            size_t my_float_buff_size = my_float_buff::get_needed_mem_single(prec);
            size_t offset=0;
            Entropy_off = offset;
            offset += my_float_buff_size;
            Density_off = offset;
            offset += my_float_buff_size;
            Hsml_off = offset;
            offset += my_float_buff_size;
            Left_off = offset;
            offset += my_float_buff_size;
            Right_off = offset;
            offset += my_float_buff_size;
            NumNgb_off = offset;
            offset += my_float_buff_size;
            Pressure_off = offset;
            offset += my_float_buff_size;
            DtEntropy_off = offset;
            offset += my_float_buff_size;
            HydroAccel0_off = offset;
            offset += my_float_buff_size;
            HydroAccel1_off = offset;
            offset += my_float_buff_size;
            #ifdef FORCETEST
            HydroAccel2_off = offset;
            offset += my_float_buff_size;
            VelPred0_off = offset;
            offset += my_float_buff_size;
            VelPred1_off = offset;
            offset += my_float_buff_size;
            #endif
            VelPred2_off = offset;
            offset += my_float_buff_size;
            DivVel_off = offset;
            offset += my_float_buff_size;
            CurlVel_off = offset;
            offset += my_float_buff_size;
            Rot0_off = offset;
            offset += my_float_buff_size;
            Rot1_off = offset;
            offset += my_float_buff_size;
            Rot2_off = offset;
            offset += my_float_buff_size;
            DhsmlDensityFactor_off = offset;
            offset += my_float_buff_size;
            MaxSignalVel_off = offset;
            offset += my_float_buff_size;
        #ifdef TIMESTEP_UPDATE
            FeedbackFlag_off = offset;
            offset += sizeof(int);
            FeedAccel0_off = offset;
            offset += my_float_buff_size;
            FeedAccel1_off = offset;
            offset += my_float_buff_size;
            FeedAccel2_off = offset;
            offset += my_float_buff_size;
        #endif
        #ifdef MORRIS97VISC
            Alpha_off = offset;
            offset += my_float_buff_size;
            DAlphaDt_off = offset;
            offset += my_float_buff_size;
        #endif
            tot_size=offset;
            return tot_size;
        };
        static inline void l_copy(sph_particle_data_buff* to, size_t tindex, sph_particle_data_buff* from, size_t findex)
        {
            void* l_from;
            void* l_to;
            l_from = (void*)((size_t)from + findex * tot_size);
            l_to = (void*)((size_t)to + tindex * tot_size);
            memcpy ( l_to, l_from, tot_size );
        };
        inline void set_init_Entropy(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Entropy_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Entropy(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Entropy_off);
            *to_store =value;
        };
        inline my_float read_re_init_Entropy(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Entropy_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Entropy(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Entropy_off);
            return *to_read;
        };
        inline void set_init_Density(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Density_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Density(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Density_off);
            *to_store =value;
        };
        inline my_float read_re_init_Density(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Density_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Density(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Density_off);
            return *to_read;
        };
        inline void set_init_Hsml(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Hsml_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Hsml(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Hsml_off);
            *to_store =value;
        };
        inline my_float read_re_init_Hsml(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Hsml_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Hsml(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Hsml_off);
            return *to_read;
        };
        inline void set_init_Left(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Left_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Left(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Left_off);
            *to_store =value;
        };
        inline my_float read_re_init_Left(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Left_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Left(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Left_off);
            return *to_read;
        };
        inline void set_init_Right(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Right_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Right(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Right_off);
            *to_store =value;
        };
        inline my_float read_re_init_Right(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Right_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Right(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Right_off);
            return *to_read;
        };
        inline void set_init_NumNgb(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + NumNgb_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_NumNgb(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + NumNgb_off);
            *to_store =value;
        };
        inline my_float read_re_init_NumNgb(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + NumNgb_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_NumNgb(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + NumNgb_off);
            return *to_read;
        };
        inline void set_init_Pressure(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Pressure_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Pressure(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Pressure_off);
            *to_store =value;
        };
        inline my_float read_re_init_Pressure(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pressure_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Pressure(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pressure_off);
            return *to_read;
        };
        inline void set_init_DtEntropy(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + DtEntropy_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_DtEntropy(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + DtEntropy_off);
            *to_store =value;
        };
        inline my_float read_re_init_DtEntropy(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + DtEntropy_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_DtEntropy(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + DtEntropy_off);
            return *to_read;
        };
        inline void set_init_HydroAccel0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + HydroAccel0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_HydroAccel0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + HydroAccel0_off);
            *to_store =value;
        };
        inline my_float read_re_init_HydroAccel0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + HydroAccel0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_HydroAccel0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + HydroAccel0_off);
            return *to_read;
        };
        inline void set_init_HydroAccel1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + HydroAccel1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_HydroAccel1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + HydroAccel1_off);
            *to_store =value;
        };
        inline my_float read_re_init_HydroAccel1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + HydroAccel1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_HydroAccel1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + HydroAccel1_off);
            return *to_read;
        };
        inline void set_init_HydroAccel2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + HydroAccel2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_HydroAccel2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + HydroAccel2_off);
            *to_store =value;
        };
        inline my_float read_re_init_HydroAccel2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + HydroAccel2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_HydroAccel2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + HydroAccel2_off);
            return *to_read;
        };
        inline void set_init_VelPred0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + VelPred0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_VelPred0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + VelPred0_off);
            *to_store =value;
        };
        inline my_float read_re_init_VelPred0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + VelPred0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_VelPred0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + VelPred0_off);
            return *to_read;
        };
        inline void set_init_VelPred1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + VelPred1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_VelPred1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + VelPred1_off);
            *to_store =value;
        };
        inline my_float read_re_init_VelPred1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + VelPred1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_VelPred1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + VelPred1_off);
            return *to_read;
        };
        inline void set_init_VelPred2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + VelPred2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_VelPred2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + VelPred2_off);
            *to_store =value;
        };
        inline my_float read_re_init_VelPred2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + VelPred2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_VelPred2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + VelPred2_off);
            return *to_read;
        };
        inline void set_init_DivVel(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + DivVel_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_DivVel(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + DivVel_off);
            *to_store =value;
        };
        inline my_float read_re_init_DivVel(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + DivVel_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_DivVel(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + DivVel_off);
            return *to_read;
        };
        inline void set_init_CurlVel(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + CurlVel_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_CurlVel(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + CurlVel_off);
            *to_store =value;
        };
        inline my_float read_re_init_CurlVel(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + CurlVel_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_CurlVel(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + CurlVel_off);
            return *to_read;
        };
        inline void set_init_Rot0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Rot0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Rot0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Rot0_off);
            *to_store =value;
        };
        inline my_float read_re_init_Rot0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Rot0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Rot0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Rot0_off);
            return *to_read;
        };
        inline void set_init_Rot1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Rot1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Rot1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Rot1_off);
            *to_store =value;
        };
        inline my_float read_re_init_Rot1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Rot1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Rot1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Rot1_off);
            return *to_read;
        };
        inline void set_init_Rot2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Rot2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Rot2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Rot2_off);
            *to_store =value;
        };
        inline my_float read_re_init_Rot2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Rot2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Rot2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Rot2_off);
            return *to_read;
        };
        inline void set_init_DhsmlDensityFactor(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + DhsmlDensityFactor_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_DhsmlDensityFactor(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + DhsmlDensityFactor_off);
            *to_store =value;
        };
        inline my_float read_re_init_DhsmlDensityFactor(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + DhsmlDensityFactor_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_DhsmlDensityFactor(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + DhsmlDensityFactor_off);
            return *to_read;
        };
        inline void set_init_MaxSignalVel(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + MaxSignalVel_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_MaxSignalVel(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + MaxSignalVel_off);
            *to_store =value;
        };
        inline my_float read_re_init_MaxSignalVel(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + MaxSignalVel_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_MaxSignalVel(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + MaxSignalVel_off);
            return *to_read;
        };
        #ifdef TIMESTEP_UPDATE
        inline void set_FeedbackFlag(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + FeedbackFlag_off);
            *to_store =value;
        };
        inline int read_FeedbackFlag(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + FeedbackFlag_off);
            return *to_read;
        };
        inline void set_init_FeedAccel0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + FeedAccel0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_FeedAccel0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + FeedAccel0_off);
            *to_store =value;
        };
        inline my_float read_re_init_FeedAccel0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + FeedAccel0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_FeedAccel0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + FeedAccel0_off);
            return *to_read;
        };
        inline void set_init_FeedAccel1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + FeedAccel1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_FeedAccel1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + FeedAccel1_off);
            *to_store =value;
        };
        inline my_float read_re_init_FeedAccel1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + FeedAccel1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_FeedAccel1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + FeedAccel1_off);
            return *to_read;
        };
        inline void set_init_FeedAccel2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + FeedAccel2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_FeedAccel2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + FeedAccel2_off);
            *to_store =value;
        };
        inline my_float read_re_init_FeedAccel2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + FeedAccel2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_FeedAccel2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + FeedAccel2_off);
            return *to_read;
        };
        #endif
        #ifdef MORRIS97VISC
        inline void set_init_Alpha(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Alpha_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Alpha(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Alpha_off);
            *to_store =value;
        };
        inline my_float read_re_init_Alpha(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Alpha_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Alpha(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Alpha_off);
            return *to_read;
        };
        inline void set_init_DAlphaDt(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + DAlphaDt_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_DAlphaDt(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + DAlphaDt_off);
            *to_store =value;
        };
        inline my_float read_re_init_DAlphaDt(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + DAlphaDt_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_DAlphaDt(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + DAlphaDt_off);
            return *to_read;
        };
        #endif
        inline void l_fill(struct sph_particle_data in, size_t index=0)
        {
        set_init_Entropy(in.Entropy,index);
        set_init_Density(in.Density,index);
        set_init_Hsml(in.Hsml,index);
        set_init_Left(in.Left,index);
        set_init_Right(in.Right,index);
        set_init_NumNgb(in.NumNgb,index);
        set_init_Pressure(in.Pressure,index);
        set_init_DtEntropy(in.DtEntropy,index);
        set_init_HydroAccel0(in.HydroAccel[0],index);
        set_init_HydroAccel1(in.HydroAccel[1],index);
        set_init_HydroAccel2(in.HydroAccel[2],index);
        set_init_VelPred0(in.VelPred[0],index);
        set_init_VelPred1(in.VelPred[1],index);
        set_init_VelPred2(in.VelPred[2],index);
        set_init_DivVel(in.DivVel,index);
        set_init_CurlVel(in.CurlVel,index);
        set_init_Rot0(in.Rot[0],index);
        set_init_Rot1(in.Rot[1],index);
        set_init_Rot2(in.Rot[2],index);
        set_init_DhsmlDensityFactor(in.DhsmlDensityFactor,index);
        set_init_MaxSignalVel(in.MaxSignalVel,index);
        #ifdef TIMESTEP_UPDATE
        set_FeedbackFlag(in.FeedbackFlag,index);
        set_init_FeedAccel0(in.FeedAccel[0],index);
        set_init_FeedAccel1(in.FeedAccel[1],index);
        set_init_FeedAccel2(in.FeedAccel[2],index);
        #endif
        #ifdef MORRIS97VISC
        set_init_Alpha(in.Alpha,index);
        set_init_DAlphaDt(in.DAlphaDt,index);
        #endif
        }
        inline sph_particle_data unpack_particle(size_t index=0)
        {
            sph_particle_data retval;
        retval.Entropy = read_re_init_Entropy(index);
        retval.Density = read_re_init_Density(index);
        retval.Hsml = read_re_init_Hsml(index);
        retval.Left = read_re_init_Left(index);
        retval.Right = read_re_init_Right(index);
        retval.NumNgb = read_re_init_NumNgb(index);
        retval.Pressure = read_re_init_Pressure(index);
        retval.DtEntropy = read_re_init_DtEntropy(index);
        retval.HydroAccel[0] = read_re_init_HydroAccel0(index);
        retval.HydroAccel[1] = read_re_init_HydroAccel1(index);
        retval.HydroAccel[2] = read_re_init_HydroAccel2(index);
        retval.VelPred[0] = read_re_init_VelPred0(index);
        retval.VelPred[1] = read_re_init_VelPred1(index);
        retval.VelPred[2] = read_re_init_VelPred2(index);
        retval.DivVel = read_re_init_DivVel(index);
        retval.CurlVel = read_re_init_CurlVel(index);
        retval.Rot[0] = read_re_init_Rot0(index);
        retval.Rot[1] = read_re_init_Rot1(index);
        retval.Rot[2] = read_re_init_Rot2(index);
        retval.DhsmlDensityFactor = read_re_init_DhsmlDensityFactor(index);
        retval.MaxSignalVel = read_re_init_MaxSignalVel(index);
        #ifdef TIMESTEP_UPDATE
        retval.FeedbackFlag = read_FeedbackFlag(index);
        retval.FeedAccel[0] = read_re_init_FeedAccel0(index);
        retval.FeedAccel[1] = read_re_init_FeedAccel1(index);
        retval.FeedAccel[2] = read_re_init_FeedAccel2(index);
        #endif
        #ifdef MORRIS97VISC
        retval.GravAccel[2] = read_re_init_Alpha(index);
        retval.DAlphaDt = read_re_init_DAlphaDt(index);
        #endif
        return retval;
        }
    }
    ;

    class NODE
    {
    public:
        my_float len;			/*!< sidelength of treenode */
        my_float center[3];		/*!< geometrical center of node */
        #ifdef ADAPTIVE_GRAVSOFT_FORGAS
        my_float maxsoft;                /*!< hold the maximum gravitational softening of particles in the
        node if the ADAPTIVE_GRAVSOFT_FORGAS option is selected */
        #endif
        my_float u_d_s[3];               /*!< center of mass of node */
        my_float u_d_mass;               /*!< mass of node */
        union
        {
            int suns[8];		/*!< temporary pointers to daughter nodes */
            struct
            {
                int bitflags;             /*!< a bit-field with various information on the node */
                int sibling;              /*!< this gives the next node in the walk in case the current node can be used */
                int nextnode;             /*!< this gives the next node in case the current node needs to be opened */
                int father;               /*!< this gives the parent node of each node (or -1 if we have the root node) */
            }
            d;
        }
        u;
    }
    ;

    class extNODE           /*!< this structure holds additional tree-node information which is not needed in the actual gravity computation */
    {
    public:
        my_float hmax;			/*!< maximum SPH smoothing length in node. Only used for gas particles */
        my_float vs[3];			/*!< center-of-mass velocity */
    }
    ;

    /*! Header for the standard file format.
     */
    struct io_header
    {
        int npart[6];                        /*!< number of particles of each type in this file */
        my_float mass[6];                      /*!< mass of particles of each type. If 0, then the masses are explicitly
        stored in the mass-block of the snapshot file, otherwise they are omitted */
        my_float time;                         /*!< time of snapshot file */
        my_float redshift;                     /*!< redshift of snapshot file */
        int flag_sfr;                        /*!< flags whether the simulation was including star formation */
        int flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
        unsigned int npartTotal[6];          /*!< total number of particles of each type in this snapshot. This can be
        different from npart if one is dealing with a multi-file snapshot. */
        int flag_cooling;                    /*!< flags whether cooling was included  */
        int num_files;                       /*!< number of files in multi-file snapshot */
        my_float BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
        my_float Omega0;                       /*!< matter density in units of critical density */
        my_float OmegaLambda;                  /*!< cosmological constant parameter */
        my_float HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
        int flag_stellarage;                 /*!< flags whether the file contains formation times of star particles */
        int flag_metals;                     /*!< flags whether the file contains metallicity values for gas and star particles */
        unsigned int npartTotalHighWord[6];  /*!< High word of the total number of particles of each type */
        int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
        char fill[60];	               /*!< fills to 256 Bytes */
    }
    ;

    #define IO_NBLOCKS 11   /*!< total number of defined information blocks for snapshot files.
    Must be equal to the number of entries in "enum iofields" */

    enum iofields           /*!< this enumeration lists the defined output blocks in snapshot files. Not all of them need to be present. */
    {
        IO_POS,
        IO_VEL,
        IO_ID,
        IO_MASS,
        IO_U,
        IO_RHO,
        IO_HSML,
        IO_POT,
        IO_ACCEL,
        IO_DTENTR,
        IO_TSTP,
    };


    /* global state of system, used for global statistics
     */
    class state_of_system
    {
     public:
        static const size_t max_cnt1=6;
        static const size_t max_cnt2=4;
        my_float Mass;
        my_float EnergyKin;
        my_float EnergyPot;
        my_float EnergyInt;
        my_float EnergyTot;
        my_float Momentum[4];
        my_float AngMomentum[4];
        my_float CenterOfMass[4];
        my_float MassComp[6];
        my_float EnergyKinComp[6];
        my_float EnergyPotComp[6];
        my_float EnergyIntComp[6];
        my_float EnergyTotComp[6];
        my_float MomentumComp[6][4];
        my_float AngMomentumComp[6][4];
        my_float CenterOfMassComp[6][4];
    }
    ;


    /* Various structures for communication
     */
    class gravdata_in
    {
    private:
        static size_t u0_off;
        static size_t u1_off;
        static size_t u2_off;
        static size_t Ninteractions_off;
        static size_t OldAcc_off;
        #ifdef UNEQUALSOFTENINGS
        static size_t Type_off;
        #ifdef ADAPTIVE_GRAVSOFT_FORGAS
        static size_t Soft_off;
        #endif
        #endif
        static size_t tot_size;
        static mpfr_prec_t prec;
    public:
        static inline size_t get_size()
        {
            return tot_size;
        };
        static inline size_t gen_size()
        {
            prec = my_float_buff::get_default_prec();
            size_t my_float_buff_size = my_float_buff::get_needed_mem_single(prec);
            size_t offset=0;
            u0_off = offset;
            offset += my_float_buff_size;
            u1_off = offset;
            offset += my_float_buff_size;
            u2_off = offset;
            offset += my_float_buff_size;
            OldAcc_off = offset;
            offset += my_float_buff_size;
            Ninteractions_off = offset;
            offset += sizeof(int);
            #ifdef UNEQUALSOFTENINGS
            Type_off = offset;
            offset += sizeof(int);
            #ifdef ADAPTIVE_GRAVSOFT_FORGAS
            Soft_off = offset;
            offset += my_float_buff_size;
            #endif
            #endif
            tot_size=offset;
            return tot_size;
        };
        inline void* get_buff_start(size_t index=0)
        {
            return (void*)((size_t)this + index * tot_size);
        }
        static inline void lcopy(gravdata_in* to, size_t tindex, gravdata_in* from, size_t findex)
        {
            void* l_from;
            void* l_to;
            l_from = (void*)((size_t)from + findex * tot_size);
            l_to = (void*)((size_t)to + tindex * tot_size);
            memcpy ( l_to, l_from, tot_size );
        };
        inline void set_init_u0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + u0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_u0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + u0_off);
            *to_store =value;
        };
        inline my_float read_re_init_u0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + u0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_u0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + u0_off);
            return *to_read;
        };
        inline void set_init_u1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + u1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_u1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + u1_off);
            *to_store =value;
        };
        inline my_float read_re_init_u1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + u1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_u1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + u1_off);
            return *to_read;
        };
        inline void set_init_u2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + u2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_u2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + u2_off);
            *to_store =value;
        };
        inline my_float read_re_init_u2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + u2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_u2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + u2_off);
            return *to_read;
        };
        inline void set_init_OldAcc(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + OldAcc_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_OldAcc(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + OldAcc_off);
            *to_store =value;
        };
        inline my_float read_re_init_OldAcc(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + OldAcc_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_OldAcc(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + OldAcc_off);
            return *to_read;
        };
        inline void set_Ninteractions(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + Ninteractions_off);
            *to_store =value;
        };
        inline int read_Ninteractions(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + Ninteractions_off);
            return *to_read;
        };
        #ifdef UNEQUALSOFTENINGS
        inline void set_Type(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + Type_off);
            *to_store =value;
        };
        inline int read_Type(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + Type_off);
            return *to_read;
        };
        #ifdef ADAPTIVE_GRAVSOFT_FORGAS
        inline void set_init_Soft(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Soft_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Soft(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Soft_off);
            *to_store =value;
        };
        inline my_float read_re_init_Soft(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Soft_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Soft(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Soft_off);
            return *to_read;
        };
        #endif
        #endif
    }
    ;


    struct gravdata_index
    {
        int Task;
        int Index;
        int SortIndex;
    }
    ;


    class densdata_in
    {
    private:
        static size_t Pos0_off;
        static size_t Pos1_off;
        static size_t Pos2_off;
        static size_t Vel0_off;
        static size_t Vel1_off;
        static size_t Vel2_off;
        static size_t Hsml_off;
        static size_t Index_off;
        static size_t Task_off;
        static size_t tot_size;
        static mpfr_prec_t prec;
    public:
        //        my_float_buff Pos[3];
        //        my_float_buff Vel[3];
        //        my_float_buff Hsml;
        //        int Index;
        //        int Task;
        static inline size_t get_size()
        {
            return tot_size;
        };
        static inline size_t gen_size()
        {
            prec = my_float_buff::get_default_prec();
            size_t my_float_buff_size = my_float_buff::get_needed_mem_single(prec);
            size_t offset=0;
            Pos0_off = offset;
            offset += my_float_buff_size;
            Pos1_off = offset;
            offset += my_float_buff_size;
            Pos2_off = offset;
            offset += my_float_buff_size;
            Vel0_off = offset;
            offset += my_float_buff_size;
            Vel1_off = offset;
            offset += my_float_buff_size;
            Vel2_off = offset;
            offset += my_float_buff_size;
            Hsml_off = offset;
            offset += my_float_buff_size;
            Index_off = offset;
            offset += sizeof(int);
            Task_off = offset;
            offset += sizeof(int);
            tot_size=offset;
            return tot_size;
        };
        inline void* get_buff_start(size_t index=0)
        {
            return (void*)((size_t)this + index * tot_size);
        }
        static inline void lcopy(densdata_in* to, size_t tindex, densdata_in* from, size_t findex)
        {
            void* l_from;
            void* l_to;
            l_from = (void*)((size_t)from + findex * tot_size);
            l_to = (void*)((size_t)to + tindex * tot_size);
            memcpy ( l_to, l_from, tot_size );
        };
        inline void set_init_Pos0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Pos0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Pos0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Pos0_off);
            *to_store =value;
        };
        inline my_float read_re_init_Pos0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Pos0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos0_off);
            return *to_read;
        };
        inline void set_init_Pos1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Pos1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Pos1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Pos1_off);
            *to_store =value;
        };
        inline my_float read_re_init_Pos1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Pos1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos1_off);
            return *to_read;
        };
        inline void set_init_Pos2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Pos2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Pos2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Pos2_off);
            *to_store =value;
        };
        inline my_float read_re_init_Pos2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Pos2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos2_off);
            return *to_read;
        };
        inline void set_init_Vel0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Vel0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Vel0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Vel0_off);
            *to_store =value;
        };
        inline my_float read_re_init_Vel0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Vel0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel0_off);
            return *to_read;
        };
        inline void set_init_Vel1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Vel1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Vel1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Vel1_off);
            *to_store =value;
        };
        inline my_float read_re_init_Vel1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Vel1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel1_off);
            return *to_read;
        };
        inline void set_init_Vel2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Vel2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Vel2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Vel2_off);
            *to_store =value;
        };
        inline my_float read_re_init_Vel2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Vel2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel2_off);
            return *to_read;
        };
        inline void set_init_Hsml(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Hsml_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Hsml(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Hsml_off);
            *to_store =value;
        };
        inline my_float read_re_init_Hsml(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Hsml_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Hsml(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Hsml_off);
            return *to_read;
        };
        inline void set_Index(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + Index_off);
            *to_store =value;
        };
        inline int read_Index(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + Index_off);
            return *to_read;
        };
        inline void set_Task(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + Task_off);
            *to_store =value;
        };
        inline int read_Task(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + Task_off);
            return *to_read;
        };
    }
    ;

    class densdata_out
    {
    private:
        static size_t Rho_off;
        static size_t Div_off;
        static size_t Rot0_off;
        static size_t Rot1_off;
        static size_t Rot2_off;
        static size_t DhsmlDensity_off;
        static size_t Ngb_off;
        static size_t tot_size;
        static mpfr_prec_t prec;
    public:
        //my_float_buff Rho;
        //my_float_buff Div, Rot[3];
        //my_float_buff DhsmlDensity;
        //my_float_buff Ngb;
        static inline size_t get_size()
        {
            return tot_size;
        };
        static inline size_t gen_size()
        {
            prec = my_float_buff::get_default_prec();
            size_t my_float_buff_size = my_float_buff::get_needed_mem_single(prec);
            size_t offset=0;
            Rho_off = offset;
            offset += my_float_buff_size;
            Div_off = offset;
            offset += my_float_buff_size;
            Rot0_off = offset;
            offset += my_float_buff_size;
            Rot1_off = offset;
            offset += my_float_buff_size;
            Rot2_off = offset;
            offset += my_float_buff_size;
            DhsmlDensity_off = offset;
            offset += my_float_buff_size;
            Ngb_off = offset;
            offset += my_float_buff_size;
            tot_size=offset;
            return tot_size;
        };
        inline void* get_buff_start(size_t index=0)
        {
            return (void*)((size_t)this + index * tot_size);
        }
        static inline void lcopy(densdata_out* to, size_t tindex, densdata_out* from, size_t findex)
        {
            void* l_from;
            void* l_to;
            l_from = (void*)((size_t)from + findex * tot_size);
            l_to = (void*)((size_t)to + tindex * tot_size);
            memcpy ( l_to, l_from, tot_size );
        };
        inline void set_init_Rho(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Rho_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Rho(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Rho_off);
            *to_store =value;
        };
        inline my_float read_re_init_Rho(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Rho_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Rho(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Rho_off);
            return *to_read;
        };
        inline void set_init_Div(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Div_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Div(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Div_off);
            *to_store =value;
        };
        inline my_float read_re_init_Div(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Div_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Div(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Div_off);
            return *to_read;
        };
        inline void set_init_Rot0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Rot0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Rot0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Rot0_off);
            *to_store =value;
        };
        inline my_float read_re_init_Rot0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Rot0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Rot0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Rot0_off);
            return *to_read;
        };
        inline void set_init_Rot1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Rot1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Rot1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Rot1_off);
            *to_store =value;
        };
        inline my_float read_re_init_Rot1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Rot1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Rot1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Rot1_off);
            return *to_read;
        };
        inline void set_init_Rot2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Rot2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Rot2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Rot2_off);
            *to_store =value;
        };
        inline my_float read_re_init_Rot2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Rot2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Rot2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Rot2_off);
            return *to_read;
        };
        inline void set_init_DhsmlDensity(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + DhsmlDensity_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_DhsmlDensity(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + DhsmlDensity_off);
            *to_store =value;
        };
        inline my_float read_re_init_DhsmlDensity(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + DhsmlDensity_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_DhsmlDensity(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + DhsmlDensity_off);
            return *to_read;
        };
        inline void set_init_Ngb(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Ngb_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Ngb(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Ngb_off);
            *to_store =value;
        };
        inline my_float read_re_init_Ngb(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Ngb_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Ngb(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Ngb_off);
            return *to_read;
        };
    }
    ;


    class hydrodata_in
    {
    private:
        static size_t Pos0_off;
        static size_t Pos1_off;
        static size_t Pos2_off;
        static size_t Vel0_off;
        static size_t Vel1_off;
        static size_t Vel2_off;
        static size_t Hsml_off;
        static size_t Mass_off;
        static size_t Density_off;
        static size_t Pressure_off;
        static size_t F1_off;
        static size_t DhsmlDensityFactor_off;
        static size_t Timestep_off;
        static size_t Task_off;
        static size_t Index_off;
        #ifdef MORRIS97VISC
        static size_t Alpha_off;
        #endif
        static size_t tot_size;
        static mpfr_prec_t prec;
    public:
        //my_float_buff Pos[3];
        //my_float_buff Vel[3];
        // my_float_buff Hsml;
        //my_float_buff Mass;
        //my_float_buff Density;
        //my_float_buff Pressure;
        //my_float_buff F1;
        //my_float_buff DhsmlDensityFactor;
        //int   Timestep;
        //int   Task;
        //int   Index;
        //#ifdef MORRIS97VISC
        //my_float Alpha;
        //#endif

        static inline size_t get_size()
        {
            return tot_size;
        };
        static inline size_t gen_size()
        {
            prec = my_float_buff::get_default_prec();
            size_t my_float_buff_size = my_float_buff::get_needed_mem_single(prec);
            size_t offset=0;
            Pos0_off = offset;
            offset += my_float_buff_size;
            Pos1_off = offset;
            offset += my_float_buff_size;
            Pos2_off = offset;
            offset += my_float_buff_size;
            Vel0_off = offset;
            offset += my_float_buff_size;
            Vel1_off = offset;
            offset += my_float_buff_size;
            Vel2_off = offset;
            offset += my_float_buff_size;
            Hsml_off = offset;
            offset += my_float_buff_size;
            Mass_off = offset;
            offset += my_float_buff_size;
            Density_off = offset;
            offset += my_float_buff_size;
            Pressure_off = offset;
            offset += my_float_buff_size;
            F1_off = offset;
            offset += my_float_buff_size;
            DhsmlDensityFactor_off = offset;
            offset += my_float_buff_size;
            Timestep_off = offset;
            offset += sizeof(int);
            Task_off = offset;
            offset += sizeof(int);
            Index_off = offset;
            offset += sizeof(int);
            #ifdef MORRIS97VISC
            DhsmlDensityFactor_off = offset;
            offset += my_float_buff_size;
            #endif
            tot_size=offset;
            return tot_size;
        };
        inline void* get_buff_start(size_t index=0)
        {
            return (void*)((size_t)this + index * tot_size);
        }
        static inline void lcopy(hydrodata_in* to, size_t tindex, hydrodata_in* from, size_t findex)
        {
            void* l_from;
            void* l_to;
            l_from = (void*)((size_t)from + findex * tot_size);
            l_to = (void*)((size_t)to + tindex * tot_size);
            memcpy ( l_to, l_from, tot_size );
        };
        inline void set_init_Pos0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Pos0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Pos0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Pos0_off);
            *to_store =value;
        };
        inline my_float read_re_init_Pos0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Pos0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos0_off);
            return *to_read;
        };
        inline void set_init_Pos1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Pos1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Pos1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Pos1_off);
            *to_store =value;
        };
        inline my_float read_re_init_Pos1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Pos1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos1_off);
            return *to_read;
        };
        inline void set_init_Pos2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Pos2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Pos2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Pos2_off);
            *to_store =value;
        };
        inline my_float read_re_init_Pos2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Pos2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos2_off);
            return *to_read;
        };
        inline void set_init_Vel0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Vel0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Vel0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Vel0_off);
            *to_store =value;
        };
        inline my_float read_re_init_Vel0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Vel0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel0_off);
            return *to_read;
        };
        inline void set_init_Vel1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Vel1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Vel1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Vel1_off);
            *to_store =value;
        };
        inline my_float read_re_init_Vel1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Vel1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel1_off);
            return *to_read;
        };
        inline void set_init_Vel2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Vel2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Vel2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Vel2_off);
            *to_store =value;
        };
        inline my_float read_re_init_Vel2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Vel2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Vel2_off);
            return *to_read;
        };
        inline void set_init_Hsml(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Hsml_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Hsml(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Hsml_off);
            *to_store =value;
        };
        inline my_float read_re_init_Hsml(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Hsml_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Hsml(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Hsml_off);
            return *to_read;
        };
        inline void set_init_Mass(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Mass_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Mass(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Mass_off);
            *to_store =value;
        };
        inline my_float read_re_init_Mass(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Mass_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Mass(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Mass_off);
            return *to_read;
        };
        inline void set_init_Density(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Density_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Density(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Density_off);
            *to_store =value;
        };
        inline my_float read_re_init_Density(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Density_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Density(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Density_off);
            return *to_read;
        };
        inline void set_init_Pressure(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Pressure_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Pressure(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Pressure_off);
            *to_store =value;
        };
        inline my_float read_re_init_Pressure(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pressure_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Pressure(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pressure_off);
            return *to_read;
        };
        inline void set_init_F1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + F1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_F1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + F1_off);
            *to_store =value;
        };
        inline my_float read_re_init_F1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + F1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_F1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + F1_off);
            return *to_read;
        };
        inline void set_init_DhsmlDensityFactor(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + DhsmlDensityFactor_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_DhsmlDensityFactor(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + DhsmlDensityFactor_off);
            *to_store =value;
        };
        inline my_float read_re_init_DhsmlDensityFactor(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + DhsmlDensityFactor_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_DhsmlDensityFactor(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + DhsmlDensityFactor_off);
            return *to_read;
        };
        inline void set_Timestep(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + Timestep_off);
            *to_store =value;
        };
        inline int read_Timestep(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + Timestep_off);
            return *to_read;
        };
        inline void set_Task(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + Task_off);
            *to_store =value;
        };
        inline int read_Task(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + Task_off);
            return *to_read;
        };
        inline void set_Index(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + Index_off);
            *to_store =value;
        };
        inline int read_Index(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + Index_off);
            return *to_read;
        };
        #ifdef MORRIS97VISC
        inline void set_init_Alpha(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Alpha_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Alpha(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Alpha_off);
            *to_store =value;
        };
        inline my_float read_re_init_Alpha(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Alpha_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Alpha(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Alpha_off);
            return *to_read;
        };
        #endif
    }
    ;

    class hydrodata_out
    {
    private:
        static size_t Acc0_off;
        static size_t Acc1_off;
        static size_t Acc2_off;
        static size_t DtEntropy_off;
        static size_t MaxSignalVel_off;
        static size_t tot_size;
        static mpfr_prec_t prec;
    public:
        static inline size_t get_size()
        {
            return tot_size;
        };
        static inline size_t gen_size()
        {
            prec = my_float_buff::get_default_prec();
            size_t my_float_buff_size = my_float_buff::get_needed_mem_single(prec);
            size_t offset=0;
            Acc0_off = offset;
            offset += my_float_buff_size;
            Acc1_off = offset;
            offset += my_float_buff_size;
            Acc2_off = offset;
            offset += my_float_buff_size;
            DtEntropy_off = offset;
            offset += my_float_buff_size;
            MaxSignalVel_off = offset;
            offset += my_float_buff_size;
            tot_size=offset;
            return tot_size;
        };
        inline void* get_buff_start(size_t index=0)
        {
            return (void*)((size_t)this + index * tot_size);
        }
        static inline void lcopy(hydrodata_out* to, size_t tindex, hydrodata_out* from, size_t findex)
        {
            void* l_from;
            void* l_to;
            l_from = (void*)((size_t)from + findex * tot_size);
            l_to = (void*)((size_t)to + tindex * tot_size);
            memcpy ( l_to, l_from, tot_size );
        };
        inline void set_init_Acc0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Acc0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Acc0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Acc0_off);
            *to_store =value;
        };
        inline my_float read_re_init_Acc0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Acc0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Acc0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Acc0_off);
            return *to_read;
        };
        inline void set_init_Acc1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Acc1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Acc1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Acc1_off);
            *to_store =value;
        };
        inline my_float read_re_init_Acc1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Acc1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Acc1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Acc1_off);
            return *to_read;
        };
        inline void set_init_Acc2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Acc2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Acc2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Acc2_off);
            *to_store =value;
        };
        inline my_float read_re_init_Acc2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Acc2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Acc2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Acc2_off);
            return *to_read;
        };
        inline void set_init_DtEntropy(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + DtEntropy_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_DtEntropy(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + DtEntropy_off);
            *to_store =value;
        };
        inline my_float read_re_init_DtEntropy(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + DtEntropy_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_DtEntropy(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + DtEntropy_off);
            return *to_read;
        };
        inline void set_init_MaxSignalVel(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + MaxSignalVel_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_MaxSignalVel(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + MaxSignalVel_off);
            *to_store =value;
        };
        inline my_float read_re_init_MaxSignalVel(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + MaxSignalVel_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_MaxSignalVel(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + MaxSignalVel_off);
            return *to_read;
        };
        //my_float_buff Acc[3];
        //my_float_buff DtEntropy;
        //my_float_buff MaxSignalVel;
    }
    ;

    #ifdef TIMESTEP_LIMITER
    class timedata_in
    {
    private:
        static size_t Pos0_off;
        static size_t Pos1_off;
        static size_t Pos2_off;
        static size_t Hsml_off;
        static size_t Size_off;
        static size_t Begin_off;
        static size_t Index_off;
        static size_t Task_off;
        static size_t tot_size;
        static mpfr_prec_t prec;
    public:
        static inline size_t get_size()
        {
            return tot_size;
        };
        static inline size_t gen_size()
        {
            prec = my_float_buff::get_default_prec();
            size_t my_float_buff_size = my_float_buff::get_needed_mem_single(prec);
            size_t offset=0;
            Pos0_off = offset;
            offset += my_float_buff_size;
            Pos1_off = offset;
            offset += my_float_buff_size;
            Pos2_off = offset;
            offset += my_float_buff_size;
            Hsml_off = offset;
            offset += my_float_buff_size;
            Size_off = offset;
            offset += sizeof(int);
            Begin_off = offset;
            offset += sizeof(int);
            Index_off = offset;
            offset += sizeof(int);
            Task_off = offset;
            offset += sizeof(int);
            tot_size=offset;
            return tot_size;
        };
        inline void* get_buff_start(size_t index=0)
        {
            return (void*)((size_t)this + index * tot_size);
        }
        static inline void lcopy(timedata_in* to, size_t tindex, timedata_in* from, size_t findex)
        {
            void* l_from;
            void* l_to;
            l_from = (void*)((size_t)from + findex * tot_size);
            l_to = (void*)((size_t)to + tindex * tot_size);
            memcpy ( l_to, l_from, tot_size );
        };
        inline void set_init_Pos0(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Pos0_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Pos0(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Pos0_off);
            *to_store =value;
        };
        inline my_float read_re_init_Pos0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos0_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Pos0(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos0_off);
            return *to_read;
        };
        inline void set_init_Pos1(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Pos1_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Pos1(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Pos1_off);
            *to_store =value;
        };
        inline my_float read_re_init_Pos1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos1_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Pos1(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos1_off);
            return *to_read;
        };
        inline void set_init_Pos2(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Pos2_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Pos2(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Pos2_off);
            *to_store =value;
        };
        inline my_float read_re_init_Pos2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos2_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Pos2(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Pos2_off);
            return *to_read;
        };
        inline void set_init_Hsml(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size + Hsml_off);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set_Hsml(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size + Hsml_off);
            *to_store =value;
        };
        inline my_float read_re_init_Hsml(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Hsml_off);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read_Hsml(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size + Hsml_off);
            return *to_read;
        };
        inline void set_Size(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + Size_off);
            *to_store =value;
        };
        inline int read_Size(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + Size_off);
            return *to_read;
        };
        inline void set_Begin(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + Begin_off);
            *to_store =value;
        };
        inline int read_Begin(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + Begin_off);
            return *to_read;
        };
        inline void set_Index(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + Index_off);
            *to_store =value;
        };
        inline int read_Index(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + Index_off);
            return *to_read;
        };
        inline void set_Task(int value, size_t index=0)
        {
            int* to_store;
            to_store = (int*)((size_t)this + index * tot_size + Task_off);
            *to_store =value;
        };
        inline int read_Task(size_t index=0)
        {
            int* to_read;
            to_read = (int*)((size_t)this + index * tot_size + Task_off);
            return *to_read;
        };
    }
    ;
    #endif

    class All_Reduce_buff
    {
//    private:
public:
        static size_t tot_size;
        static mpfr_prec_t prec;
//    public:
        static inline size_t get_size()
        {
            return tot_size;
        };
        static inline size_t gen_size()
        {
            prec = my_float_buff::get_default_prec();
            size_t my_float_buff_size = my_float_buff::get_needed_mem_single(prec);
            tot_size=my_float_buff_size;
            return tot_size;
        };
        inline void* get_buff_start(size_t index=0)
        {
            return (void*)((size_t)this + index * tot_size);
        }
        static inline void lcopy(All_Reduce_buff* to, size_t tindex, All_Reduce_buff* from, size_t findex)
        {
            void* l_from;
            void* l_to;
            l_from = (void*)((size_t)from + findex * tot_size);
            l_to = (void*)((size_t)to + tindex * tot_size);
            memcpy ( l_to, l_from, tot_size );
        };
        inline void set_init(my_float value, size_t index=0)
        {
            my_float_buff* to_prep;
            to_prep = (my_float_buff*)((size_t)this + index * tot_size);
            my_float_buff::place_pmpreal((void*)to_prep, prec);
            *to_prep= value;
        };
        inline void set(my_float value, size_t index=0)
        {
            my_float_buff* to_store;
            to_store = (my_float_buff*)((size_t)this + index * tot_size);
            *to_store =value;
        };
        inline my_float read_re_init(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size);
            my_float_buff::re_org_pmpreal(to_read);
            return *to_read;
        };
        inline my_float read(size_t index=0)
        {
            my_float_buff* to_read;
            to_read = (my_float_buff*)((size_t)this + index * tot_size);
            return *to_read;
        };
    }
    ;
    #ifndef NOMPI

    #endif

    #endif


