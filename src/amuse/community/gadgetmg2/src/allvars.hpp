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
#include <gsl/gsl_rng.h>
#include "tags.hpp"

//#include <mpreal.h>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost::multiprecision;

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

//#ifdef DOUBLEPRECISION             /*!< If defined, the variable type FLOAT is set to "double", otherwise to FLOAT */
//#define FLOAT double
//#else
//#define FLOAT double
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


struct DomainNODE
{
  double s[3];                     /*!< center-of-mass coordinates */
  double vs[3];                    /*!< center-of-mass velocities */
  double mass;                     /*!< mass of node */
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
  int   bitflags;                 /*!< this bit-field encodes the particle type with the largest softening among the particles of the nodes, and whether there are particles with different softening in the node */
#else
  double maxsoft;                  /*!< hold the maximum gravitational softening of particles in the
                                       node if the ADAPTIVE_GRAVSOFT_FORGAS option is selected */
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
struct particle_data
{
  double Pos[3];			/*!< particle position at its current time */
  double Mass;			/*!< particle mass */
  double Vel[3];			/*!< particle velocity at its current time */
  double GravAccel[3];		/*!< particle acceleration due to gravity */

#ifdef FORCETEST
  double GravAccelDirect[3];	/*!< particle acceleration when computed with direct summation */
#endif
  double Potential;		/*!< gravitational potential */
  double OldAcc;			/*!< magnitude of old gravitational force. Used in relative opening criterion */
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
  double GravCost;		/*!< weight factor used for balancing the work-load */
#ifdef PSEUDOSYMMETRIC
  double AphysOld;               /*!< magnitude of acceleration in last timestep. Used to make a first order
                                     prediction of the change of acceleration expected in the future, thereby
                                     allowing to guess whether a decrease/increase of the timestep should occur
                                     in the timestep that is started. */
#endif
}
 ;

/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
struct sph_particle_data
{
  double Entropy;                /*!< current value of entropy (actually entropic function) of particle */
  double Density;		/*!< current baryonic mass density of particle */
  double Hsml;			/*!< current smoothing length */
  double Left;                   /*!< lower bound in iterative smoothing length search */
  double Right;                  /*!< upper bound in iterative smoothing length search */
  double NumNgb;                 /*!< weighted number of neighbours found */
  double Pressure;		/*!< current pressure */
  double DtEntropy;              /*!< rate of change of entropy */
  double HydroAccel[3];		/*!< acceleration due to hydrodynamical force */
  double VelPred[3];		/*!< predicted SPH particle velocity at the current time */
  double DivVel;			/*!< local velocity divergence */
  double CurlVel;		/*!< local velocity curl */
  double Rot[3];		        /*!< local velocity curl */
  double DhsmlDensityFactor;     /*!< correction factor needed in the equation of motion of the conservative entropy formulation of SPH */
  double MaxSignalVel;           /*!< maximum "signal velocity" occuring for this particle */
#ifdef TIMESTEP_UPDATE
  int   FeedbackFlag;
  double FeedAccel[3];  /*!< acceleration due to feedback force */
#endif
#ifdef MORRIS97VISC
  double Alpha;		        /*!< viscosity coefficient */
  double DAlphaDt;       		/*!< time rate of change of viscosity coefficient */
#endif
}
 ;


struct NODE
{
  double len;			/*!< sidelength of treenode */
  double center[3];		/*!< geometrical center of node */
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double maxsoft;                /*!< hold the maximum gravitational softening of particles in the
                                     node if the ADAPTIVE_GRAVSOFT_FORGAS option is selected */
#endif
  union
  {
    int suns[8];		/*!< temporary pointers to daughter nodes */
    struct
    {
      double s[3];               /*!< center of mass of node */
      double mass;               /*!< mass of node */
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

struct extNODE           /*!< this structure holds additional tree-node information which is not needed in the actual gravity computation */
{
  double hmax;			/*!< maximum SPH smoothing length in node. Only used for gas particles */
  double vs[3];			/*!< center-of-mass velocity */
}
 ;

/*! Header for the standard file format.
 */
struct io_header
{
  int npart[6];                        /*!< number of particles of each type in this file */
  double mass[6];                      /*!< mass of particles of each type. If 0, then the masses are explicitly
                                            stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                         /*!< time of snapshot file */
  double redshift;                     /*!< redshift of snapshot file */
  int flag_sfr;                        /*!< flags whether the simulation was including star formation */
  int flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];          /*!< total number of particles of each type in this snapshot. This can be
                                            different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;                    /*!< flags whether cooling was included  */
  int num_files;                       /*!< number of files in multi-file snapshot */
  double BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;                       /*!< matter density in units of critical density */
  double OmegaLambda;                  /*!< cosmological constant parameter */
  double HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
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
struct state_of_system
{
  double Mass;
  double EnergyKin;
  double EnergyPot;
  double EnergyInt;
  double EnergyTot;
  double Momentum[4];
  double AngMomentum[4];
  double CenterOfMass[4];
  double MassComp[6];
  double EnergyKinComp[6];
  double EnergyPotComp[6];
  double EnergyIntComp[6];
  double EnergyTotComp[6];
  double MomentumComp[6][4];
  double AngMomentumComp[6][4];
  double CenterOfMassComp[6][4];
}
 ;


/* Various structures for communication
 */
struct gravdata_in
{
  union
  {
    double Pos[3];
    double Acc[3];
    double Potential;
  }
  u;
#ifdef UNEQUALSOFTENINGS
  int Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double Soft;
#endif
#endif
  union
  {
    double OldAcc;
    int Ninteractions;
  }
  w;
}
 ;

struct gravdata_index
{
  int Task;
  int Index;
  int SortIndex;
}
 ;


struct densdata_in
{
  double Pos[3];
  double Vel[3];
  double Hsml;
  int Index;
  int Task;
}
 ;

struct densdata_out
{
  double Rho;
  double Div, Rot[3];
  double DhsmlDensity;
  double Ngb;
}
 ;


struct hydrodata_in
{
  double Pos[3];
  double Vel[3];
  double Hsml;
  double Mass;
  double Density;
  double Pressure;
  double F1;
  double DhsmlDensityFactor;
  int   Timestep;
  int   Task;
  int   Index;
#ifdef MORRIS97VISC
  double Alpha;
#endif
}
 ;

struct hydrodata_out
{
  double Acc[3];
  double DtEntropy;
  double MaxSignalVel;
}
 ;

#ifdef TIMESTEP_LIMITER
  struct timedata_in
  {
    double  Pos[3];
    double  Hsml;
    int    Size;
    int    Begin;
    int    Index;
    int    Task;
  }
  ;
#endif

#ifndef NOMPI

#endif

#endif
