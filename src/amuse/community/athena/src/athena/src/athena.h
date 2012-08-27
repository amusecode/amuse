#ifndef ATHENA_H
#define ATHENA_H 
/*============================================================================*/
/*! \file athena.h
 *  \brief Contains definitions of many data types and structures.
 *
 * PURPOSE: Contains definitions of the following data types and structures:
 * - Real    - either float or double, depending on configure option
 * - ConsS   - cell-centered conserved variables
 * - PrimS   - cell-centered primitive variables
 * - Cons1DS - conserved variables in 1D: same as ConsS minus Bx
 * - Prim1DS - primitive variables in 1D: same as PrimS minus Bx
 * - GrainS  - basic properties of particles
 * - GridS   - everything in a single Grid
 * - DomainS - everything in a single Domain (potentially many Grids)
 * - MeshS   - everything across whole Mesh (potentially many Domains)
 * - OutputS - everything associated with an individual output		      */
/*============================================================================*/
#include "defs.h"

#ifdef MPI_PARALLEL
#include "mpi.h"
#endif

/*! \typedef Real
 *  \brief Variable precision float, depends on macro set by configure.
 */ 
#if defined(SINGLE_PREC)
#ifdef MPI_PARALLEL
#error: MPI requires double precision
#endif /*MPI_PARALLEL */
typedef float  Real;
#elif defined(DOUBLE_PREC)
typedef double Real;
#else
# error "Not a valid precision flag"
#endif

/*! \struct Real3Vect
 *  \brief General 3-vectors of Reals.
 */
typedef struct Real3Vect_s{
  Real x, y, z;
}Real3Vect;
/*! \struct Int3Vect
 *  \brief General 3-vectors of ints.
 */
typedef struct Int3Vect_s{
  int i, j, k;
}Int3Vect;

/*! \struct SideS
 *  \brief Sides of a cube, used to find overlaps between Grids 
 *   at different levels.
 */
typedef struct Side_s{
  int ijkl[3];    /*!< indices of left-sides  in each dir [0,1,2]=[i,j,k] */ 
  int ijkr[3];    /*!< indices of right-sides in each dir [0,1,2]=[i,j,k] */ 
}SideS;

/*! \struct GridsDataS
 *  \brief Number of zones in, and identifying information about, a Grid.
 */
typedef struct GridsData_s{
  int Nx[3];       /*!< number of zones in each dir [0,1,2]=[x1,x2,x3] */
  int Disp[3];   /*!< i,j,k displacements from origin of root [0,1,2]=[i,j,k] */
  int ID_Comm_world;      /*!< ID of process for this Grid in MPI_COMM_WORLD */
  int ID_Comm_Domain;     /*!< ID of process for this Grid in Comm_Domain    */
#ifdef STATIC_MESH_REFINEMENT
  int ID_Comm_Children;     /*!< ID updating this Grid in Comm_Domain    */
  int ID_Comm_Parent;     /*!< ID updating this Grid in Comm_Domain    */
#endif
}GridsDataS;



/*----------------------------------------------------------------------------*/
/*! \struct ConsS
 *  \brief Conserved variables.
 *  IMPORTANT!! The order of the elements in ConsS CANNOT be changed.
 */
typedef struct Cons_s{
  Real d;			/*!< density */
  Real M1;			/*!< momentum density in 1-direction*/
  Real M2;			/*!< momentum density in 2-direction*/
  Real M3;			/*!< momentum density in 3-direction*/
#ifndef BAROTROPIC
  Real E;			/*!< total energy density */
#endif /* BAROTROPIC */
#ifdef MHD
  Real B1c;			/*!< cell centered magnetic fields in 1-dir*/
  Real B2c;			/*!< cell centered magnetic fields in 2-dir*/
  Real B3c;			/*!< cell centered magnetic fields in 3-dir*/
#endif /* MHD */
#if (NSCALARS > 0)
  Real s[NSCALARS];             /*!< passively advected scalars */
#endif
#ifdef CYLINDRICAL
  Real Pflux;	 		/*!< pressure component of flux */
#endif
}ConsS;

/*----------------------------------------------------------------------------*/
/*! \struct PrimS
 *  \brief Primitive variables.
 *  IMPORTANT!! The order of the elements in PrimS CANNOT be changed.
 */
typedef struct Prim_s{
  Real d;			/*!< density */
  Real V1;			/*!< velocity in 1-direction */
  Real V2;			/*!< velocity in 2-direction */
  Real V3;			/*!< velocity in 3-direction */
#ifndef BAROTROPIC
  Real P;			/*!< pressure */
#endif /* BAROTROPIC */
#ifdef MHD
  Real B1c;                     /*!< cell centered magnetic fields in 1-dir */
  Real B2c;                     /*!< cell centered magnetic fields in 2-dir */
  Real B3c;                     /*!< cell centered magnetic fields in 3-dir */
#endif /* MHD */
#if (NSCALARS > 0)
  Real r[NSCALARS];             /*!< density-normalized advected scalars */
#endif
}PrimS;

/*----------------------------------------------------------------------------*/
/*! \struct Cons1DS
 *  \brief Conserved variables in 1D (does not contain Bx).
 *  IMPORTANT!! The order of the elements in Cons1DS CANNOT be changed.
 */
typedef struct Cons1D_s{
  Real d;			/*!< density */
  Real Mx;			/*!< momentum density in X,Y,Z; where X is    */
  Real My;                      /*!< direction longitudinal to 1D slice; which*/
  Real Mz;                      /*!< can be in any dimension: 1,2,or 3        */
#ifndef BAROTROPIC
  Real E;			/*!< total energy density */
#endif /* BAROTROPIC */
#ifdef MHD
  Real By;			/*!< cell centered magnetic fields in Y */
  Real Bz;			/*!< cell centered magnetic fields in Z */
#endif /* MHD */
#if (NSCALARS > 0)
  Real s[NSCALARS];             /*!< passively advected scalars */
#endif
#ifdef CYLINDRICAL
  Real Pflux;	 		/*!< pressure component of flux */
#endif
}Cons1DS;

/*----------------------------------------------------------------------------*/
/*! \struct Prim1DS
 *  \brief Primitive variables in 1D (does not contain Bx).
 *  IMPORTANT!! The order of the elements in Prim1DS CANNOT be changed.
 */
typedef struct Prim1D_s{
  Real d;			/*!< density */
  Real Vx;			/*!< velocity in X-direction */
  Real Vy;			/*!< velocity in Y-direction */
  Real Vz;			/*!< velocity in Z-direction */
#ifndef BAROTROPIC
  Real P;			/*!< pressure */
#endif /* BAROTROPIC */
#ifdef MHD
  Real By;			/*!< cell centered magnetic fields in Y-dir */
  Real Bz;			/*!< cell centered magnetic fields in Z-dir */
#endif /* MHD */
#if (NSCALARS > 0)
  Real r[NSCALARS];             /*!< density-normalized advected scalars */
#endif
}Prim1DS;

/*----------------------------------------------------------------------------*/
/*! \struct GrainS
 *  \brief Basic quantities for one pseudo-particle.
 */
#ifdef PARTICLES

/* Physical quantities of a dust particle */
typedef struct Grain_s{
  Real x1,x2,x3;	/*!< coordinate in X,Y,Z */
  Real v1,v2,v3;	/*!< velocity in X,Y,Z */
  int property;		/*!< index of grain properties */
  short pos;		/*!< position: 0: ghost; 1: grid; >=10: cross out/in; */
  long my_id;		/*!< particle id */
#ifdef MPI_PARALLEL
  int init_id;          /*!< particle's initial host processor id */
#endif
#ifdef FARGO
  Real shift;           /*!< amount of shift in x2 direction */
#endif
}GrainS;

/*! \struct GrainAux
 *  \brief Auxilary quantities for a dust particle. */
typedef struct GrainAux_s{
  Real dpar;            /*!< local particle density */
#ifdef FARGO
  Real shift;           /*!< amount of shift in x2 direction */
#endif
}GrainAux;

/*! \struct Grain_Property
 *  \brief List of physical grain properties. */
typedef struct Grain_Property_s{
#ifdef FEEDBACK
  Real m;		/*!< mass of this type of particle */
#endif
  Real rad;		/*!< radius of this type of particle (cm) */
  Real rho;		/*!< solid density of this type of particle (g/cm^3) */
  long num;		/*!< number of particles with this property */
  short integrator;	/*!< integrator type: exp (1), semi (2) or full (3) */
}Grain_Property;

/*! \struct GPCouple
 *  \brief Grid elements for gas-particle coupling. */
typedef struct GPCouple_s{
  Real grid_d;		/*!< gas density (at 1/2 step) */
  Real grid_v1;		/*!< gas 1-velocity (at 1/2 step) */
  Real grid_v2;		/*!< gas 2-velocity (at 1/2 step) */
  Real grid_v3;		/*!< gas 3-velocity (at 1/2 step) */
#ifndef BAROTROPIC
  Real grid_cs;		/*!< gas sound speed */
#endif
#ifdef FEEDBACK
  Real fb1;             /*!< 1-momentum feedback to the grid */
  Real fb2;             /*!< 2-momentum feedback to the grid */
  Real fb3;             /*!< 3-momentum feedback to the grid */
  Real FBstiff;         /*!< stiffness of the feedback term */
  Real Eloss;           /*!< energy dissipation */
#endif
}GPCouple;

#endif /* PARTICLES */

/*----------------------------------------------------------------------------*/
/*! \struct GridOvrlpS
 *  \brief Contains information about Grid overlaps, used for SMR.
 */
#ifdef STATIC_MESH_REFINEMENT
typedef struct GridOvrlp_s{
  int ijks[3];         /*!< start ijk on this Grid of overlap [0,1,2]=[i,j,k] */
  int ijke[3];         /*!< end   ijk on this Grid of overlap [0,1,2]=[i,j,k] */
  int ID, DomN;        /*!< processor ID, and Domain #, of OVERLAP Grid */
  int nWordsRC, nWordsP; /*!< # of words communicated for Rest/Corr and Prol */
  ConsS **myFlx[6];   /*!< fluxes of conserved variables at 6 boundaries */
#ifdef MHD
  Real **myEMF1[6];      /*!< fluxes of magnetic field (EMF1) at 6 boundaries */
  Real **myEMF2[6];      /*!< fluxes of magnetic field (EMF2) at 6 boundaries */
  Real **myEMF3[6];      /*!< fluxes of magnetic field (EMF3) at 6 boundaries */
#endif
}GridOvrlpS;
#endif /* STATIC_MESH_REFINEMENT */


#ifdef AMUSE
typedef struct BoundaryCell_s {
    ConsS ***LeftX1;
    ConsS ***RightX1;
    ConsS ***LeftX2;
    ConsS ***RightX2;
    ConsS ***LeftX3;
    ConsS ***RightX3;
} BoundaryCellS;
#endif /* AMUSE */

/*----------------------------------------------------------------------------*/
/*! \struct GridS
 *  \brief 3D arrays of dependent variables, plus grid data, plus particle data,
 *   plus data about child and parent Grids, plus MPI rank information for a
 *   Grid.
 *
 *   Remember a Grid is defined as the region of a Domain at some
 *   refinement level being updated by a single processor.  Uses an array of
 *   ConsS, rather than arrays of each variable, to increase locality of data
 *   for a given cell in memory.  */

typedef struct Grid_s{
  ConsS ***U;                /*!< conserved variables */
#ifdef MHD
  Real ***B1i,***B2i,***B3i;    /*!< interface magnetic fields */
#ifdef RESISTIVITY
  Real ***eta_Ohm,***eta_Hall,***eta_AD; /*!< magnetic diffusivities */ 
#endif
#endif /* MHD */
#ifdef SELF_GRAVITY
  Real ***Phi, ***Phi_old;      /*!< gravitational potential */
  Real ***x1MassFlux;           /*!< x1 mass flux for source term correction */
  Real ***x2MassFlux;           /*!< x2 mass flux for source term correction */
  Real ***x3MassFlux;           /*!< x3 mass flux for source term correction */
#endif /* GRAVITY */
  Real MinX[3];       /*!< min(x) in each dir on this Grid [0,1,2]=[x1,x2,x3] */
  Real MaxX[3];       /*!< max(x) in each dir on this Grid [0,1,2]=[x1,x2,x3] */
  Real dx1,dx2,dx3;   /*!< cell size on this Grid */
  Real time, dt;           /*!< current time and timestep  */
  int is,ie;		   /*!< start/end cell index in x1 direction */
  int js,je;		   /*!< start/end cell index in x2 direction */
  int ks,ke;		   /*!< start/end cell index in x3 direction */
  int Nx[3];     /*!< # of zones in each dir on Grid [0,1,2]=[x1,x2,x3] */
  int Disp[3];   /*!< i,j,k displacements of Grid from origin [0,1,2]=[i,j,k] */

  int rx1_id, lx1_id;  /*!< ID of Grid to R/L in x1-dir (default=-1; no Grid) */
  int rx2_id, lx2_id;  /*!< ID of Grid to R/L in x2-dir (default=-1; no Grid) */
  int rx3_id, lx3_id;  /*!< ID of Grid to R/L in x3-dir (default=-1; no Grid) */

#ifdef PARTICLES
  int partypes;              /*!< number of particle types */
  Grain_Property *grproperty;/*!< array of particle properties of all types */
  long nparticle;            /*!< number of particles */
  long arrsize;              /*!< size of the particle array */
  Grain *particle;           /*!< array of all particles */
  GrainAux *parsub;          /*!< supplemental particle information */
  GPCouple ***Coup;          /*!< array of gas-particle coupling */
#endif /* PARTICLES */

#ifdef STATIC_MESH_REFINEMENT
  int NCGrid;        /*!< # of child  Grids that overlap this Grid */
  int NPGrid;        /*!< # of parent Grids that this Grid overlaps */
  int NmyCGrid;      /*!< # of child  Grids on same processor as this Grid */
  int NmyPGrid;      /*!< # of parent Grids on same processor (either 0 or 1) */

  GridOvrlpS *CGrid;  /*!< 1D array of data for NCGrid child  overlap regions */
  GridOvrlpS *PGrid;  /*!< 1D array of data for NPGrid parent overlap regions */
/* NB: The first NmyCGrid[NmyPGrid] elements of these arrays contain overlap
 * regions being updated by the same processor as this Grid */
#endif /* STATIC_MESH_REFINEMENT */

#ifdef CYLINDRICAL
  Real *r,*ri;                  /*!< cylindrical scaling factors */ 
#endif /* CYLINDRICAL */

#ifdef AMUSE
  BoundaryCellS * boundary;
#endif
}GridS;

/*! \fn void (*VGFun_t)(GridS *pG)
 *  \brief Generic void function of Grid. */
typedef void (*VGFun_t)(GridS *pG);    /* generic void function of Grid */

/*----------------------------------------------------------------------------*/
/*! \struct DomainS
 *  \brief Information about one region of Mesh at some particular level.
 *
 * Contains pointer to a single Grid, even though the Domain may contain many
 * Grids, because for any general parallelization mode, no more than one Grid
 * can exist per Domain per processor.
 *
 * The i,j,k displacements are measured in units of grid cells on this Domain
 */
typedef struct Domain_s{
  Real RootMinX[3]; /*!< min(x) in each dir on root Domain [0,1,2]=[x1,x2,x3] */
  Real RootMaxX[3]; /*!< max(x) in each dir on root Domain [0,1,2]=[x1,x2,x3] */
  Real MinX[3];     /*!< min(x) in each dir on this Domain [0,1,2]=[x1,x2,x3] */
  Real MaxX[3];     /*!< max(x) in each dir on this Domain [0,1,2]=[x1,x2,x3] */
  Real dx[3];                /*!< cell size in this Domain [0,1,2]=[x1,x2,x3] */
  int Nx[3];    /*!< # of zones in each dir in this Domain [0,1,2]=[x1,x2,x3] */
  int NGrid[3]; /*!< # of Grids in each dir in this Domain [0,1,2]=[x1,x2,x3] */
  int Disp[3]; /*!< i,j,k displacements of Domain from origin [0,1,2]=[i,j,k] */
  int Level,DomNumber;   /*!< level and ID number of this Domain */
  int InputBlock;      /*!< # of <domain> block in input file for this Domain */
  GridS *Grid;     /*!< pointer to Grid in this Dom updated on this processor */

  GridsDataS ***GData;/*!< size,location,& processor IDs of Grids in this Dom */

  VGFun_t ix1_BCFun, ox1_BCFun;/*!< ix1/ox1 BC function pointers for this Dom */
  VGFun_t ix2_BCFun, ox2_BCFun;/*!< ix1/ox1 BC function pointers for this Dom */
  VGFun_t ix3_BCFun, ox3_BCFun;/*!< ix1/ox1 BC function pointers for this Dom */

#ifdef MPI_PARALLEL
  MPI_Comm Comm_Domain;      /*!< MPI communicator between Grids on this Dom */
  MPI_Group Group_Domain;    /*!< MPI group for Domain communicator */
#ifdef STATIC_MESH_REFINEMENT
  MPI_Comm Comm_Parent;      /*!< MPI communicator to Grids in parent Domain  */
  MPI_Comm Comm_Children;    /*!< MPI communicator to Grids in  child Domains */
  MPI_Group Group_Children;  /*!< MPI group for Children communicator */
#endif /* STATIC_MESH_REFINEMENT */
#endif /* MPI_PARALLEL */
}DomainS;

/*! \fn void (*VDFun_t)(DomainS *pD)
 *  \brief Generic void function of Domain. */
typedef void (*VDFun_t)(DomainS *pD);  /* generic void function of Domain */

/*----------------------------------------------------------------------------*/
/*! \struct MeshS
 *  \brief Information about entire mesh hierarchy, including array of Domains.
 */

typedef struct Mesh_s{
  Real RootMinX[3]; /*!< min(x) in each dir on root Domain [0,1,2]=[x1,x2,x3] */
  Real RootMaxX[3]; /*!< max(x) in each dir on root Domain [0,1,2]=[x1,x2,x3] */
  Real dx[3];     /*!< cell size on root Domain [0,1,2]=[x1,x2,x3] */
  Real time, dt;  /*!< current time and timestep for entire Mesh */
  int Nx[3];    /*!< # of zones in each dir on root Domain [0,1,2]=[x1,x2,x3] */
  int nstep;                 /*!< number of integration steps taken */
  int BCFlag_ix1, BCFlag_ox1;  /*!< BC flag on root domain for inner/outer x1 */
  int BCFlag_ix2, BCFlag_ox2;  /*!< BC flag on root domain for inner/outer x2 */
  int BCFlag_ix3, BCFlag_ox3;  /*!< BC flag on root domain for inner/outer x3 */
  int NLevels;               /*!< overall number of refinement levels in mesh */
  int *DomainsPerLevel;      /*!< number of Domains per level (DPL) */
  DomainS **Domain;        /*!< array of Domains, indexed over levels and DPL */
  char *outfilename;         /*!< basename for output files containing -id#  */
}MeshS;

/*----------------------------------------------------------------------------*/
/* OutputS: Everything for outputs. */

struct Output_s;
/*! \fn void (*VOutFun_t)(MeshS *pM, struct Output_s *pout)
 *  \brief Output function pointer. */
typedef void (*VOutFun_t)(MeshS *pM, struct Output_s *pout);
/*! \fn void (*VResFun_t)(MeshS *pM, struct Output_s *pout)
 *  \brief Restart function pointer. */
typedef void (*VResFun_t)(MeshS *pM, struct Output_s *pout);
/*! \fn Real (*ConsFun_t)(const GridS *pG, const int i,const int j,const int k) 
 *  \brief Pointer to expression that computes quant for output.*/
typedef Real (*ConsFun_t)(const GridS *pG, const int i,const int j,const int k);
#ifdef PARTICLES
/*! \fn int (*PropFun_t)(const Grain *gr, const GrainAux *grsub)
 *  \brief Particle property selection function */
typedef int (*PropFun_t)(const Grain *gr, const GrainAux *grsub);
typedef Real (*Parfun_t)(const Grid *pG, const Grain *gr);
typedef Real (*Parfun_t)(const Grid *pG, const Grain *gr);
#endif

/*! \struct OutputS
 *  \brief Everything for outputs. */
typedef struct Output_s{
  int n;          /*!< the N from the <outputN> block of this output */
  Real dt;        /*!< time interval between outputs  */
  Real t;         /*!< next time to output */
  int num;        /*!< dump number (0=first) */
  char *out;      /*!< variable (or user fun) to be output */
  char *id;       /*!< filename is of the form <basename>[.idump][.id].<ext> */
#ifdef PARTICLES
  int out_pargrid;    /*!< output grid binned particles (=1) or not (=0) */
  PropFun_t par_prop; /*!< particle property selection function */
#endif

/* level and domain number of output (default = [-1,-1] = output all levels) */

  int nlevel, ndomain;

/* variables which describe data min/max */
  Real dmin,dmax;   /*!< user defined min/max for scaling data */
  Real gmin,gmax;   /*!< computed global min/max (over all output data) */
  int sdmin,sdmax;  /*!< 0 = auto scale, otherwise use dmin/dmax */

/* variables which describe coordinates of output data volume */
  int ndim;       /*!< 3=cube 2=slice 1=vector 0=scalar */
  int reduce_x1;  /*!< flag to denote reduction in x1 (0=no reduction) */
  int reduce_x2;  /*!< flag to denote reduction in x2 (0=no reduction) */
  int reduce_x3;  /*!< flag to denote reduction in x3 (0=no reduction) */
  Real x1l, x1u;  /*!< lower/upper x1 range for data slice  */
  Real x2l, x2u;  /*!< lower/upper x2 range for data slice  */
  Real x3l, x3u;  /*!< lower/upper x3 range for data slice  */

/* variables which describe output format */
  char *out_fmt;  /*!< output format = {bin, tab, hdf, hst, pgm, ppm, ...} */
  char *dat_fmt;  /*!< format string for tabular type output, e.g. "%10.5e" */
  char *palette;  /*!< name of palette for RGB conversions */
  float *rgb;     /*!< array of RGB[256*3] values derived from palette */
  float *der;     /*!< helper array of derivatives for interpolation into RGB */

/* pointers to output functions; data expressions */
  VOutFun_t out_fun; /*!< output function pointer */
  VResFun_t res_fun; /*!< restart function pointer */
  ConsFun_t expr;   /*!< pointer to expression that computes quant for output */

}OutputS;


/*----------------------------------------------------------------------------*/
/* typedefs for functions:
 */
/* for static gravitational potential and cooling, set in problem generator,
 * and used by integrators */

/*! \fn Real (*GravPotFun_t)(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential function. */
typedef Real (*GravPotFun_t)(const Real x1, const Real x2, const Real x3);
#ifdef CYLINDRICAL
/*! \fn Real (*StaticGravAcc_t)(const Real x1, const Real x2, const Real x3)
 *  \brief Static gravitational acceleration. */
typedef Real (*StaticGravAcc_t)(const Real x1, const Real x2, const Real x3);
#ifdef FARGO
/*! \fn Real (*OrbitalFun_t)(const Real x1)
 *  \brief Orbital function for FARGO */
typedef Real (*OrbitalFun_t)(const Real x1);
/*! \fn Real (*ShearFun_t)(const Real x1)
 *  \brief Shear function for FARGO */
typedef Real (*ShearFun_t)(const Real x1);
#endif
#endif /* Cylindrical */
/*! \fn Real (*CoolingFun_t)(const Real d, const Real p, const Real dt);
 *  \brief Cooling function. */
typedef Real (*CoolingFun_t)(const Real d, const Real p, const Real dt);
#ifdef RESISTIVITY
/*! \fn void (*EtaFun_t)(GridS *pG, int i, int j, int k,
                         Real *eta_O, Real *eta_H, Real *eta_A)
 *  \brief Resistivity Eta Function. */
typedef void (*EtaFun_t)(GridS *pG, int i, int j, int k,
                         Real *eta_O, Real *eta_H, Real *eta_A);
#endif /* RESISTIVITY */

#ifdef PARTICLES
/* function types for interpolation schemes and stopping time */
/*! \fn void (*WeightFun_t)(GridS *pG, Real x1, Real x2, Real x3,
  Real3Vector cell1, Real weight[3][3][3], int *is, int *js, int *ks);
 *  \brief Interpolation scheme for particles. */
typedef void (*WeightFun_t)(GridS *pG, Real x1, Real x2, Real x3,
  Real3Vector cell1, Real weight[3][3][3], int *is, int *js, int *ks);
/*! \fn Real (*TSFun_t)(GridS *pG, int type, Real rho, Real cs, Real vd)
 *  \brief Stopping time function for particles. */
typedef Real (*TSFun_t)(GridS *pG, int type, Real rho, Real cs, Real vd);
#endif /* PARTICLES */

/*----------------------------------------------------------------------------*/
/*! \enum BCDirection
 *  \brief Directions for the set_bvals_fun() function */
enum BCDirection {left_x1, right_x1, left_x2, right_x2, left_x3, right_x3};

#endif /* ATHENA_H */
