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
 * \file        src/main/allvars.h
 * \date        05/2018
 * \brief       All (global) variables.
 * \details
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 30.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifndef ALLVARS_H
#define ALLVARS_H

#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>
#include <mpi.h>
#include <stddef.h>
#include <stdio.h>

#include "./arepoconfig.h"

#ifdef IMPOSE_PINNING
#include <hwloc.h>
#endif /* #ifdef IMPOSE_PINNING */

#include "../time_integration/timestep.h"
#include "../utils/dtypes.h"
#include "../utils/tags.h"

#define AREPO_VERSION "Arepo public 1.0" /* code version string */

/* default values for unspecified config options */

#if defined(__linux__) && !defined(HOST_MEMORY_REPORTING)
#define HOST_MEMORY_REPORTING
#endif /* #if defined(__linux__) && !defined(HOST_MEMORY_REPORTING) */

#ifndef LOAD_TYPES
#define LOAD_TYPES 0xff
#endif /* #ifndef LOAD_TYPES */

#if defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS)
#define REFINEMENT
#else /* #if defined (REFINEMENT_SPLIT_CELLS) || defined (REFINEMENT_MERGE_CELLS) */
#undef REFINEMENT
#endif /* #if defined (REFINEMENT_SPLIT_CELLS) || defined (REFINEMENT_MERGE_CELLS) #else */

#ifndef NTYPES
#define NTYPES 6
#endif /* #ifndef NTYPES */

#ifndef NSOFTTYPES
#define NSOFTTYPES NTYPES
#endif /* #ifndef NSOFTTYPES */

#if !defined(OUTPUT_PRESSURE_GRADIENT) && !defined(OUTPUT_DENSITY_GRADIENT) && !defined(OUTPUT_VELOCITY_GRADIENT) && \
    !defined(OUTPUT_BFIELD_GRADIENT) && !defined(OUTPUT_DIVVEL) && !defined(OUTPUT_CURLVEL) && !defined(OUTPUT_VORTICITY)
// only if no gradient output defined, no need to update them directly before output.
#else /* #if !defined(OUTPUT_PRESSURE_GRADIENT) && !defined(OUTPUT_DENSITY_GRADIENT) && !defined(OUTPUT_VELOCITY_GRADIENT) && \
         !defined(OUTPUT_BFIELD_GRADIENT) && !defined(OUTPUT_DIVVEL) && !defined(OUTPUT_CURLVEL) && !defined(OUTPUT_VORTICITY) */
#define UPDATE_GRADIENTS_FOR_OUTPUT
#endif /* #if !defined(OUTPUT_PRESSURE_GRADIENT) && !defined(OUTPUT_DENSITY_GRADIENT) && !defined(OUTPUT_VELOCITY_GRADIENT) &&        \
          !defined(OUTPUT_BFIELD_GRADIENT) && !defined(OUTPUT_DIVVEL) && !defined(OUTPUT_CURLVEL) && !defined(OUTPUT_VORTICITY) #else \
        */

#ifdef ADAPTIVE_HYDRO_SOFTENING
#ifndef NSOFTTYPES_HYDRO
#define NSOFTTYPES_HYDRO 64
#endif /* #ifndef NSOFTTYPES_HYDRO */
#else  /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#undef NSOFTTYPES_HYDRO
#define NSOFTTYPES_HYDRO 0
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING #else */

#if defined(SAVE_HSML_IN_SNAPSHOT)
#define SUBFIND_CALC_MORE
#endif /* #if defined(SAVE_HSML_IN_SNAPSHOT) */

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
#define NO_SELFGRAVITY_TYPE \
  EXACT_GRAVITY_FOR_PARTICLE_TYPE                        // exclude particle type from self-gravity (can be used with exact gravity)
#define NO_GRAVITY_TYPE EXACT_GRAVITY_FOR_PARTICLE_TYPE  // disable computation of gravity on particle type
#define EXACT_GRAVITY_REACTION                           // include reaction to other particle types when using exact gravity
#endif                                                   /* #ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE */

/* restrictions on config option combinations */
#if(NSOFTTYPES + NSOFTTYPES_HYDRO) >= 254
#error "(NSOFTTYPES + NSOFTTYPES_HYDRO) >= 254"
#endif /* #if (NSOFTTYPES + NSOFTTYPES_HYDRO) >= 254 */

#if NSOFTTYPES < 2
#error "NSOFTTYPES < 2"
#endif /* #if NSOFTTYPES < 2 */

#if defined(HOST_MEMORY_REPORTING) && !defined(__linux__)
#error "HOST_MEMORY_REPORTING only works under Linux."
#endif /* #if defined(HOST_MEMORY_REPORTING) && !defined(__linux__) */

#if defined(USE_DIRECT_IO_FOR_RESTARTS) && !defined(__linux__)
#error "USE_DIRECT_IO_FOR_RESTARTS only works under Linux."
#endif /* #if defined(USE_DIRECT_IO_FOR_RESTARTS) && !defined(__linux__) */

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
#if !((INDIVIDUAL_GRAVITY_SOFTENING + 0) >= 1)
#error "set INDIVIDUAL_GRAVITY_SOFTENING to a bitmask of particle types"
#endif /* #if !((INDIVIDUAL_GRAVITY_SOFTENING+0) >= 1) */
#endif /* #ifdef INDIVIDUAL_GRAVITY_SOFTENING */

#ifdef OUTPUTPOTENTIAL
#ifndef EVALPOTENTIAL
#error "the option OUTPUTPOTENTIAL requires EVALPOTENTIAL"
#endif /* #ifndef EVALPOTENTIAL */
#endif /* #ifdef OUTPUTPOTENTIAL */

#if defined(CELL_CENTER_GRAVITY) && defined(SELFGRAVITY)
#ifndef HIERARCHICAL_GRAVITY
#error "the of option CELL_CENTER_GRAVITY requires HIERARCHICAL_GRAVITY"
#endif /* #ifndef HIERARCHICAL_GRAVITY */
#endif /* #if defined(CELL_CENTER_GRAVITY) && defined(SELFGRAVITY) */

#ifdef MHD
#ifndef RIEMANN_HLLD
#error "the of option MHD requires RIEMANN_HLLD"
#endif /* #ifndef RIEMANN_HLLD */
#endif /* #ifdef MHD */

/* optional additional headers based on config options */

#include "../utils/timer.h"

#if defined(COOLING)
#include "../cooling/cooling_vars.h"
#endif /* #if defined(COOLING) */

#ifdef ADDBACKGROUNDGRID
#include "../add_backgroundgrid/add_bggrid.h"
#endif /* #ifdef ADDBACKGROUNDGRID */

/* function mappings and macros */

#ifdef MPI_HYPERCUBE_ALLGATHERV
#define MPI_Allgatherv MPI_hypercube_Allgatherv
#endif /* #ifdef MPI_HYPERCUBE_ALLGATHERV */

#ifdef MPISENDRECV_CHECKSUM
#define MPI_Sendrecv MPI_Check_Sendrecv
#endif /* #ifdef MPISENDRECV_CHECKSUM */

#define terminate(...)                                                                                                            \
  {                                                                                                                               \
    if(FlagNyt == 0)                                                                                                              \
      {                                                                                                                           \
        char termbuf1[1000], termbuf2[1000];                                                                                      \
        sprintf(termbuf1, "TERMINATE: ******!!!!!******  Code termination on task=%d, function %s(), file %s, line %d", ThisTask, \
                __FUNCTION__, __FILE__, __LINE__);                                                                                \
        sprintf(termbuf2, __VA_ARGS__);                                                                                           \
        printf("%s: %s\n", termbuf1, termbuf2);                                                                                   \
        fflush(stdout);                                                                                                           \
        FlagNyt = 1;                                                                                                              \
        MPI_Abort(MPI_COMM_WORLD, 1);                                                                                             \
      }                                                                                                                           \
    exit(1);                                                                                                                      \
  }
#define mpi_terminate(...)    \
  {                           \
    if(ThisTask == 0)         \
      terminate(__VA_ARGS__); \
  }
#define warn(...)                                                                                                            \
  {                                                                                                                          \
    char termbuf1[1000], termbuf2[1000];                                                                                     \
    sprintf(termbuf1, "WARNING: Code warning on task=%d, function %s(), file %s, line %d", ThisTask, __FUNCTION__, __FILE__, \
            __LINE__);                                                                                                       \
    sprintf(termbuf2, __VA_ARGS__);                                                                                          \
    printf("%s: %s\n", termbuf1, termbuf2);                                                                                  \
    myflush(stdout);                                                                                                         \
    FILE *fd = fopen("WARNINGS", "a");                                                                                       \
    fprintf(fd, "%s: %s\n", termbuf1, termbuf2);                                                                             \
    fclose(fd);                                                                                                              \
  }

/* define an "assert" macro which outputs MPI task (we do NOT want to
   call MPI_Abort, because then the assertion failure isn't caught in
   the debugger) */
#define myassert(cond)                                                                                                              \
  if(!(cond))                                                                                                                       \
    {                                                                                                                               \
      char termbuf[1000];                                                                                                           \
      sprintf(termbuf, "Assertion failure!\n\ttask=%d, function %s(), file %s, line %d:\n\t%s\n", ThisTask, __FUNCTION__, __FILE__, \
              __LINE__, #cond);                                                                                                     \
      printf("%s", termbuf);                                                                                                        \
      myflush(stdout);                                                                                                              \
      assert(0);                                                                                                                    \
    }

/* memory manager */
#define mymalloc(x, y) mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__, 0, NULL)
#define mymalloc_g(x, y) mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__, 0, callorigin)
#define mymalloc_clear(x, y) mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__, 1, NULL)
#define mymalloc_movable(x, y, z) mymalloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__, NULL)
#define mymalloc_movable_g(x, y, z) mymalloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__, callorigin)
#define myrealloc(x, y) myrealloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define myrealloc_movable(x, y) myrealloc_movable_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define myfree(x) myfree_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)
#define myfree_movable(x) myfree_movable_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)

#define MAX_FIRST_ELEMENTS_CONSIDERED \
  5 /* This sets the number of lowest loaded tasks to be considered for assignment of next domain patch */

#define NUMBER_OF_MEASUREMENTS_TO_RECORD 6

#ifndef GRAVCOSTLEVELS
#define GRAVCOSTLEVELS 6
#endif /* #ifndef  GRAVCOSTLEVELS */

#define MODE_LOCAL_NO_EXPORT -1
#define MODE_LOCAL_PARTICLES 0
#define MODE_IMPORTED_PARTICLES 1
#define MODE_FINISHED 2

#ifndef DIRECT_SUMMATION_THRESHOLD
#define DIRECT_SUMMATION_THRESHOLD 3000
#endif /* #ifndef DIRECT_SUMMATION_THRESHOLD */

#define MODE_FIRST_HALFSTEP 0
#define MODE_SECOND_HALFSTEP 1

#define FLAG_PARTIAL_TREE 0
#define FLAG_FULL_TREE 1

#ifndef MPI_MESSAGE_SIZELIMIT_IN_MB
#define MPI_MESSAGE_SIZELIMIT_IN_MB 200
#endif /* #ifndef MPI_MESSAGE_SIZELIMIT_IN_MB */

#define MPI_MESSAGE_SIZELIMIT_IN_BYTES ((MPI_MESSAGE_SIZELIMIT_IN_MB)*1024LL * 1024LL)

#define COMMBUFFERSIZE (32 * 1024LL * 1024LL)

#define NUM_THREADS 1 /* no OpenMP support in this code! */

extern int Nforces;
extern int *TargetList;

extern struct thread_data
{
  int Nexport __attribute__((__aligned__(64))); /* to align on different cache lines */
  int NexportNodes;
  int Interactions;
  int dummy;
  double Cost;

  double Costtotal;  /*!< The total cost of the particles/nodes processed by each thread */
  double Ewaldcount; /*!< The total cost for the Ewald correction per thread */
  int FirstExec;     /*!< Keeps track, if a given thread executes the gravity_primary_loop() for the first time */

  size_t ExportSpace;
  size_t InitialSpace;
  size_t ItemSize;

  int *P_CostCount;
  int *TreePoints_CostCount;
  int *Node_CostCount;

  struct data_partlist *PartList;

  int *Ngblist;
  double *R2list;
  int *Exportflag;
  int *toGoDM;
  int *toGoSph;

} Thread[NUM_THREADS];

/* If we use a static Voronoi mesh with local timestepping and no rebuild of
 * the static mesh, then we need to backup the face areas before calling
 * compute_interface_fluxes(), because this function calls face_get_normals()
 * which sets some face area to 0 under some circumstances */
#if defined(VORONOI_STATIC_MESH) && !defined(FORCE_EQUAL_TIMESTEPS) && !defined(VORONOI_STATIC_MESH_DO_DOMAIN_DECOMPOSITION)
#define VORONOI_BACKUP_RESTORE_FACE_AREAS
#else /* #if defined(VORONOI_STATIC_MESH) && !defined(FORCE_EQUAL_TIMESTEPS) && !defined(VORONOI_STATIC_MESH_DO_DOMAIN_DECOMPOSITION) \
       */
#undef VORONOI_BACKUP_RESTORE_FACE_AREAS
#endif /* #if defined(VORONOI_STATIC_MESH) && !defined(FORCE_EQUAL_TIMESTEPS) && \
          !defined(VORONOI_STATIC_MESH_DO_DOMAIN_DECOMPOSITION) #else */

#ifdef IMPOSE_PINNING
extern hwloc_cpuset_t cpuset_thread[NUM_THREADS];
#endif /* #ifdef IMPOSE_PINNING */

#ifdef ONEDIMS
#define ALLOC_TOLERANCE 0.3
#else /* #ifdef ONEDIMS */
#define ALLOC_TOLERANCE 0.1
#endif /* #ifdef ONEDIMS #else */
#define ALLOC_STARBH_ROOM 0.02

#ifdef TOLERATE_WRITE_ERROR
#define IO_TRIALS 20
#define IO_SLEEP_TIME 10
#endif /* #ifdef TOLERATE_WRITE_ERROR */

/* calculate appropriate value of MAXSCALARS */

#if defined(REFINEMENT_HIGH_RES_GAS) || defined(PASSIVE_SCALARS)

#ifdef REFINEMENT_HIGH_RES_GAS
#define COUNT_REFINE 1
#else /* #ifdef  REFINEMENT_HIGH_RES_GAS */
#define COUNT_REFINE 0
#endif /* #ifdef  REFINEMENT_HIGH_RES_GAS #else */

#ifdef PASSIVE_SCALARS
#define COUNT_PASSIVE_SCALARS PASSIVE_SCALARS
#else /* #ifdef PASSIVE_SCALARS */
#define COUNT_PASSIVE_SCALARS 0
#endif /* #ifdef PASSIVE_SCALARS #else */

#define MAXSCALARS (COUNT_REFINE + COUNT_PASSIVE_SCALARS)
#endif /* #if defined(REFINEMENT_HIGH_RES_GAS) ||  defined(PASSIVE_SCALARS)*/

/* calculate appropriate value of MAXGRADIENTS */

#define COUNT_GRAD_DEFAULT 5

#ifdef MHD
#define COUNT_GRAD_MHD 3
#else /* #ifdef MHD */
#define COUNT_GRAD_MHD 0
#endif /* #ifdef MHD #else */

#ifdef MAXSCALARS
#define COUNT_GRAD_SCALARS MAXSCALARS
#else /* #ifdef MAXSCALARS */
#define COUNT_GRAD_SCALARS 0
#endif /* #ifdef MAXSCALARS #else*/

#define MAXGRADIENTS (COUNT_GRAD_DEFAULT + COUNT_GRAD_MHD + COUNT_GRAD_SCALARS)

/*************************************/

/*! For Peano-Hilbert order.
 *  Note: Maximum is 10 to fit in 32-bit integer,
 *  maximum is 21 to fit into 64-bit integer,
 *  and 42 is the absolute maximum, for which 128-bit integers are needed
 */
#ifndef BITS_PER_DIMENSION
#define BITS_PER_DIMENSION 42
#endif /* #ifndef  BITS_PER_DIMENSION */
#if(BITS_PER_DIMENSION <= 21)
typedef unsigned long long peanokey;
#else  /* #if (BITS_PER_DIMENSION <= 21) */
typedef __int128 peanokey;
#endif /* #if (BITS_PER_DIMENSION <= 21) #else */
#if(BITS_PER_DIMENSION <= 31)
typedef unsigned int peano1D;
#else /* #if (BITS_PER_DIMENSION <= 31) */
#if(BITS_PER_DIMENSION <= 42)
typedef unsigned long long peano1D;
#else /* #if (BITS_PER_DIMENSION <= 42) */
#error "BITS_PER_DIMENSION can be at most 42"
#endif /* #if (BITS_PER_DIMENSION <= 42) #else */
#endif /* #if (BITS_PER_DIMENSION <= 31) #else */

#define PEANOCELLS (((peanokey)1) << (3 * BITS_PER_DIMENSION))

#define MAX_FLOAT_NUMBER 1e37
#define MIN_FLOAT_NUMBER 1e-37
#define MAX_DOUBLE_NUMBER 1e306
#define MIN_DOUBLE_NUMBER 1e-306

#ifdef DOUBLEPRECISION
#if(DOUBLEPRECISION == 2)
#define MAX_REAL_NUMBER MAX_FLOAT_NUMBER
#define MIN_REAL_NUMBER MIN_FLOAT_NUMBER
#else /* #if (DOUBLEPRECISION==2) */
#define MAX_REAL_NUMBER MAX_DOUBLE_NUMBER
#define MIN_REAL_NUMBER MIN_DOUBLE_NUMBER
#endif /* #if (DOUBLEPRECISION==2) #else */
#else  /* #ifdef DOUBLEPRECISION */
#define MAX_REAL_NUMBER MAX_FLOAT_NUMBER
#define MIN_REAL_NUMBER MIN_FLOAT_NUMBER
#endif /* #ifdef DOUBLEPRECISION #else */

#ifndef GAMMA
#define GAMMA (5. / 3.) /*!< adiabatic index of simulated gas */
#endif                  /* #ifndef  GAMMA */
#define GAMMA_MINUS1 (GAMMA - 1.)
#define GAMMA_PLUS1 (GAMMA + 1.)

#define HYDROGEN_MASSFRAC 0.76 /*!< mass fraction of hydrogen, relevant only for radiative cooling */
#define HE_ABUND ((1. / HYDROGEN_MASSFRAC - 1.) / 4.)

/* ... often used physical constants (cgs units; NIST 2010) */

#define GRAVITY 6.6738e-8
#define SOLAR_MASS 1.989e33
#define SOLAR_LUM 3.826e33
#define SOLAR_EFF_TEMP 5.780e3
#define RAD_CONST 7.5657e-15
#define AVOGADRO 6.02214e23
#define BOLTZMANN 1.38065e-16
#define GAS_CONST 8.31446e7
#define CLIGHT 2.99792458e10

#define PLANCK 6.6260695e-27
#define PARSEC 3.085678e18
#define KILOPARSEC 3.085678e21
#define MEGAPARSEC 3.085678e24
#define ASTRONOMICAL_UNIT 1.49598e13
#define PROTONMASS 1.67262178e-24
#define ELECTRONMASS 9.1093829e-28
#define THOMPSON 6.65245873e-25
#define ELECTRONCHARGE 4.8032042e-10
#define HUBBLE 3.2407789e-18      /* in h/sec */
#define LYMAN_ALPHA 1215.6e-8     /* 1215.6 Angstroem */
#define LYMAN_ALPHA_HeII 303.8e-8 /* 303.8 Angstroem */
#define OSCILLATOR_STRENGTH 0.41615
#define OSCILLATOR_STRENGTH_HeII 0.41615
#define ELECTRONVOLT_IN_ERGS 1.60217656e-12

#define SEC_PER_GIGAYEAR 3.15576e16
#define SEC_PER_MEGAYEAR 3.15576e13
#define SEC_PER_YEAR 3.15576e7

#ifndef FOF_PRIMARY_LINK_TYPES
#define FOF_PRIMARY_LINK_TYPES 2
#endif /* #ifndef FOF_PRIMARY_LINK_TYPES */

#ifndef FOF_SECONDARY_LINK_TYPES
#define FOF_SECONDARY_LINK_TYPES 0
#endif /* #ifndef FOF_SECONDARY_LINK_TYPES */

#ifndef ASMTH
/*! ASMTH gives the scale of the short-range/long-range force split in units
 * of FFT-mesh cells
 */
#define ASMTH 1.25
#endif /* #ifndef ASMTH */

#ifndef RCUT
/*! RCUT gives the maximum distance (in units of the scale used for the force
 * split) out to which short-range forces are evaluated in the short-range
 * tree walk.
 */
#define RCUT 4.5
#endif /* #ifndef RCUT */

#define MAXLEN_OUTPUTLIST 1100  /*!< maxmimum number of entries in output list */
#define MAXLEN_PATH 256         /*!< maximum length of various filenames (full path) */
#define MAXLEN_PARAM_TAG 50     /*!< maximum length of the tag of a parameter in the parameter file */
#define MAXLEN_PARAM_VALUE 200  /*!< maximum length of the value of a parameter in the parameter file */
#define MAX_PARAMETERS 300      /*!< maximum number of parameters in the parameter file */
#define DRIFT_TABLE_LENGTH 1000 /*!< length of the lookup table used to hold the drift and kick factors */

#define BASENUMBER 100
#define HIGHRESMASSFAC 0.5

#define MAXITER 300000 /*! Maximum number of iterations before process is terminated */

#ifndef FOF_LINKLENGTH
#define FOF_LINKLENGTH 0.2
#endif /* #ifndef FOF_LINKLENGTH */

#ifndef FOF_GROUP_MIN_LEN
#define FOF_GROUP_MIN_LEN 32
#endif /* #ifndef FOF_GROUP_MIN_LEN */

typedef struct
{
  double r;
  double mass;
} sort_r2list;

typedef struct
{
  MyFloat r2;
  int index;
} r2type;

#include "../mesh/mesh.h"
#include "../mesh/voronoi/voronoi.h"

struct unbind_data
{
  int index;
};

#ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
extern MPI_Status mpistat;
#undef MPI_STATUS_IGNORE
#define MPI_STATUS_IGNORE &mpistat
#endif /* #ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG */

#define FLT(x) (x)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif /* #ifndef M_PI */

#define TO_MBYTE_FAC (1.0 / (1024.0 * 1024.0))

#ifdef ONEDIMS
#define NUMDIMS 1
#define KERNEL_COEFF_1 (4.0 / 3)
#define KERNEL_COEFF_2 (8.0)
#define KERNEL_COEFF_3 (24.0)
#define KERNEL_COEFF_4 (16.0)
#define KERNEL_COEFF_5 (8.0 / 3)
#define KERNEL_COEFF_6 (-8.0)
#define NORM_COEFF 2.0
#else /* #ifdef   ONEDIMS */
#ifndef TWODIMS
#define NUMDIMS 3                     /*!< For 3D-normalized kernel */
#define KERNEL_COEFF_1 2.546479089470 /*!< Coefficients for SPH spline kernel and its derivative */
#define KERNEL_COEFF_2 15.278874536822
#define KERNEL_COEFF_3 45.836623610466
#define KERNEL_COEFF_4 30.557749073644
#define KERNEL_COEFF_5 5.092958178941
#define KERNEL_COEFF_6 (-15.278874536822)
#define NORM_COEFF 4.188790204786                 /*!< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */
#else                                             /* #ifndef  TWODIMS */
#define NUMDIMS 2                                 /*!< For 2D-normalized kernel */
#define KERNEL_COEFF_1 (5.0 / 7 * 2.546479089470) /*!< Coefficients for SPH spline kernel and its derivative */
#define KERNEL_COEFF_2 (5.0 / 7 * 15.278874536822)
#define KERNEL_COEFF_3 (5.0 / 7 * 45.836623610466)
#define KERNEL_COEFF_4 (5.0 / 7 * 30.557749073644)
#define KERNEL_COEFF_5 (5.0 / 7 * 5.092958178941)
#define KERNEL_COEFF_6 (5.0 / 7 * (-15.278874536822))
#define NORM_COEFF M_PI /*!< Coefficient for kernel normalization. */
#endif                  /* #ifndef  TWODIMS #else */
#endif                  /* #ifdef   ONEDIMS #else*/

#define SOFTFAC1 10.666666666667 /*!< Coefficients for gravitational softening */
#define SOFTFAC2 32.0
#define SOFTFAC3 (-38.4)
#define SOFTFAC4 (-2.8)
#define SOFTFAC5 5.333333333333
#define SOFTFAC6 6.4
#define SOFTFAC7 (-9.6)
#define SOFTFAC8 21.333333333333
#define SOFTFAC9 (-48.0)
#define SOFTFAC10 38.4
#define SOFTFAC11 (-10.666666666667)
#define SOFTFAC12 (-0.066666666667)
#define SOFTFAC13 (-3.2)
#define SOFTFAC14 0.066666666667
#define SOFTFAC15 (-16.0)
#define SOFTFAC16 9.6
#define SOFTFAC17 (-2.133333333333)
#define SOFTFAC18 128.0
#define SOFTFAC19 (-115.2)
#define SOFTFAC20 21.333333333333
#define SOFTFAC21 (-96.0)
#define SOFTFAC22 115.2
#define SOFTFAC23 (-42.666666666667)
#define SOFTFAC24 0.1333333333333

extern MyDouble boxSize, boxHalf;
#ifdef LONG_X
extern MyDouble boxSize_X, boxHalf_X;
#else /* #ifdef LONG_X */
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif /* #ifdef LONG_X #else */
#ifdef LONG_Y
extern MyDouble boxSize_Y, boxHalf_Y;
#else /* #ifdef LONG_Y */
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif /* #ifdef LONG_Y #else */
#ifdef LONG_Z
extern MyDouble boxSize_Z, boxHalf_Z;
#else /* #ifdef LONG_Z */
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif /* #ifdef LONG_Z #else */

#if !defined(GRAVITY_NOT_PERIODIC)
#define GRAVITY_NEAREST_X(x) \
  (xtmp = (x), (xtmp > boxHalf_X) ? (xtmp - boxSize_X) : ((xtmp < -boxHalf_X) ? (xtmp + boxSize_X) : (xtmp)))
#define GRAVITY_NEAREST_Y(x) \
  (ytmp = (x), (ytmp > boxHalf_Y) ? (ytmp - boxSize_Y) : ((ytmp < -boxHalf_Y) ? (ytmp + boxSize_Y) : (ytmp)))
#define GRAVITY_NEAREST_Z(x) \
  (ztmp = (x), (ztmp > boxHalf_Z) ? (ztmp - boxSize_Z) : ((ztmp < -boxHalf_Z) ? (ztmp + boxSize_Z) : (ztmp)))
#else /* #if !defined(GRAVITY_NOT_PERIODIC) */
#define GRAVITY_NEAREST_X(x) (x)
#define GRAVITY_NEAREST_Y(x) (x)
#define GRAVITY_NEAREST_Z(x) (x)
#endif /* #if !defined(GRAVITY_NOT_PERIODIC) #else */

#if !defined(GRAVITY_NOT_PERIODIC)
#define FOF_NEAREST_LONG_X(x) (xtmp = fabs(x), (xtmp > boxHalf_X) ? (boxSize_X - xtmp) : xtmp)
#define FOF_NEAREST_LONG_Y(x) (ytmp = fabs(x), (ytmp > boxHalf_Y) ? (boxSize_Y - ytmp) : ytmp)
#define FOF_NEAREST_LONG_Z(x) (ztmp = fabs(x), (ztmp > boxHalf_Z) ? (boxSize_Z - ztmp) : ztmp)
#else /* #if !defined(GRAVITY_NOT_PERIODIC) */
#define FOF_NEAREST_LONG_X(x) fabs(x)
#define FOF_NEAREST_LONG_Y(x) fabs(x)
#define FOF_NEAREST_LONG_Z(x) fabs(x)
#endif /* #if !defined(GRAVITY_NOT_PERIODIC) #else */

/* periodicity of gas */
#ifndef REFLECTIVE_X
#define NGB_PERIODIC_LONG_X(x) (xtmp = fabs(x), (xtmp > boxHalf_X) ? (boxSize_X - xtmp) : xtmp)
#define NEAREST_X(x) (xtmp = (x), (xtmp > boxHalf_X) ? (xtmp - boxSize_X) : ((xtmp < -boxHalf_X) ? (xtmp + boxSize_X) : (xtmp)))
#define WRAP_X(x) (xtmp = (x), (xtmp > boxSize_X) ? (xtmp - boxSize_X) : ((xtmp < 0) ? (xtmp + boxSize_X) : (xtmp)))
#else /* #ifndef REFLECTIVE_X */
#define NGB_PERIODIC_LONG_X(x) fabs(x)
#define NEAREST_X(x) (x)
#define WRAP_X(x) (x)
#endif /* #ifndef REFLECTIVE_X #else */

#ifndef REFLECTIVE_Y
#define NGB_PERIODIC_LONG_Y(x) (ytmp = fabs(x), (ytmp > boxHalf_Y) ? (boxSize_Y - ytmp) : ytmp)
#define NEAREST_Y(x) (ytmp = (x), (ytmp > boxHalf_Y) ? (ytmp - boxSize_Y) : ((ytmp < -boxHalf_Y) ? (ytmp + boxSize_Y) : (ytmp)))
#define WRAP_Y(x) (ytmp = (x), (ytmp > boxSize_Y) ? (ytmp - boxSize_Y) : ((ytmp < 0) ? (ytmp + boxSize_Y) : (ytmp)))
#else /* #ifndef REFLECTIVE_Y */
#define NGB_PERIODIC_LONG_Y(x) fabs(x)
#define NEAREST_Y(x) (x)
#define WRAP_Y(x) (x)
#endif /* #ifndef REFLECTIVE_Y #else */

#ifndef REFLECTIVE_Z
#define NGB_PERIODIC_LONG_Z(x) (ztmp = fabs(x), (ztmp > boxHalf_Z) ? (boxSize_Z - ztmp) : ztmp)
#define NEAREST_Z(x) (ztmp = (x), (ztmp > boxHalf_Z) ? (ztmp - boxSize_Z) : ((ztmp < -boxHalf_Z) ? (ztmp + boxSize_Z) : (ztmp)))
#define WRAP_Z(x) (ztmp = (x), (ztmp > boxSize_Z) ? (ztmp - boxSize_Z) : ((ztmp < 0) ? (ztmp + boxSize_Z) : (ztmp)))
#else /* #ifndef REFLECTIVE_Z */
#define NGB_PERIODIC_LONG_Z(x) fabs(x)
#define NEAREST_Z(x) (x)
#define WRAP_Z(x) (x)
#endif /* #ifndef REFLECTIVE_Z #else */

#define FACT1 0.366025403785 /* FACT1 = 0.5 * (sqrt(3)-1) */
#define FAC_TWO_TO_TWO_THIRDS 1.5874011

/*********************************************************/
/*  Global variables                                     */
/*********************************************************/

extern int TimeBinSynchronized[TIMEBINS];
extern struct TimeBinData TimeBinsHydro, TimeBinsGravity;

#ifdef USE_SFR
extern double TimeBinSfr[TIMEBINS];
#endif /* #ifdef USE_SFR */

extern int ThisTask; /*!< the number of the local processor  */
extern int NTask;    /*!< number of processors */
extern int PTask;    /*!< note: NTask = 2^PTask */

extern int ThisNode;        /*!< the rank of the current compute node  */
extern int NumNodes;        /*!< the number of compute nodes used  */
extern int MinTasksPerNode; /*!< the minimum number of MPI tasks that is found on any of the nodes  */
extern int MaxTasksPerNode; /*!< the maximum number of MPI tasks that is found on any of the nodes  */
extern int TasksInThisNode; /*!< number of MPI tasks on  current compute node */
extern int RankInThisNode;  /*!< rank of the MPI task on the current compute node */
extern long long MemoryOnNode;

extern double CPUThisRun; /*!< Sums CPU time of current process */

extern int MaxTopNodes; /*!< Maximum number of nodes in the top-level tree used for domain decomposition */

extern int RestartFlag; /*!< taken from command line used to start code. 0 is normal start-up from
                             initial conditions, 1 is resuming a run from a set of restart files, while 2
                             marks a restart from a snapshot file. */
extern int RestartSnapNum;
extern int TakeLevel;
extern int TagOffset;

extern int Argc;
extern char **Argv;

extern double CPU_Step[CPU_LAST];
extern double CPU_Step_Stored[CPU_LAST];

extern double WallclockTime; /*!< This holds the last wallclock time measurement for timings measurements */
extern double StartOfRun;    /*!< This stores the time of the start of the run for evaluating the elapsed time */

extern size_t AllocatedBytes;
extern size_t FreeBytes;

extern char DumpFlag;
extern char DumpFlagNextSnap;

extern int FlagNyt;

extern int NumPart; /*!< number of particles on the LOCAL processor */
extern int NumGas;  /*!< number of gas particles on the LOCAL processor  */

extern gsl_rng *random_generator;     /*!< a random number generator  */
extern gsl_rng *random_generator_aux; /*!< an auxialiary random number generator for use if one doesn't want to influence the main
                                         code's random numbers  */

#ifdef USE_SFR
extern int Stars_converted; /*!< current number of star particles in gas particle block */
#endif                      /* #ifdef USE_SFR */

#ifdef TOLERATE_WRITE_ERROR
extern int WriteErrorFlag;
extern char AlternativeOutputDir[MAXLEN_PATH];
#endif /* #ifdef TOLERATE_WRITE_ERROR */

extern double EgyInjection;

extern double TimeOfLastDomainConstruction; /*!< holds what it says */

extern double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
extern double DomainInverseLen, DomainBigFac;
extern int *DomainStartList, *DomainEndList;
extern double *DomainCost, *TaskCost;
extern int *DomainCount, *TaskCount;
extern struct no_list_data
{
  int task;
  int no;
  int domainCount;
  double domainCost;
} * ListNoData;

extern int domain_bintolevel[TIMEBINS];
extern int domain_refbin[TIMEBINS];
extern int domain_grav_weight[TIMEBINS];
extern int domain_hydro_weight[TIMEBINS];
extern int domain_to_be_balanced[TIMEBINS];

/*! Array of task numbers holding the respective top-level nodes. For
    the topnodes entries, it is indexed by the Leaf member, for
    pseudoparticles it is indexed by the node
    number-MaxPart-MaxNodes.  */
extern int *DomainTask;
extern int *DomainNewTask;

/*! Array of indices of the main tree nodes that are identical to the
 *  top-level nodes. For the topnodes entries, it is indexed by the
 *  Leaf member, for pseudoparticles it is indexed by the node
 *  number-MaxPart-MaxNodes.
 */
extern int *DomainNodeIndex;

extern peanokey *Key, *KeySorted;

/*! The top node structure is an octree used for encoding the domain
 *  decomposition. Its leaf nodes are the units into which the domain
 *  is decomposed.
 */
extern struct topnode_data
{
  peanokey Size;
  peanokey StartKey;
  long long Count;
  /*! The index of the first daughter node. The remaining 7 follow
      sequentially, I think. */
  int Daughter;
  /*! The index of this topnode in the DomainTask etc arrays. Is this
      only valid for topnodes that have daughter=-1, i.e. the actual
      leaves? */
  int Leaf;
  unsigned char MortonToPeanoSubnode[8];
} * TopNodes;

extern int NTopnodes, NTopleaves;

/*! Variables for gravitational tree */
extern int Tree_MaxPart;
extern int Tree_NumNodes;
extern int Tree_MaxNodes;
extern int Tree_FirstNonTopLevelNode;
extern int Tree_NumPartImported;
extern int Tree_NumPartExported;
extern int Tree_ImportedNodeOffset;
extern int Tree_NextFreeNode;

extern int *Tree_ResultIndexList;
extern int *Tree_Task_list;
extern MyDouble *Tree_Pos_list;
extern unsigned long long *Tree_IntPos_list;

extern struct treepoint_data
{
  MyDouble Pos[3];
  unsigned long long IntPos[3];
  MyDouble Mass;
  float OldAcc;
  int index;
  int th;
  unsigned char level;
  unsigned char Type;
  unsigned char SofteningType : 7;
#ifndef HIERARCHICAL_GRAVITY
  unsigned char ActiveFlag : 1;
#endif /* #ifndef HIERARCHICAL_GRAVITY */

#if defined(SUBFIND) && defined(SUBFIND_EXTENDED_PROPERTIES)
  MyFloat GroupRad;
  int GrNr;
#endif /* #if defined(SUBFIND) && defined(SUBFIND_EXTENDED_PROPERTIES) */
} * Tree_Points;

extern struct resultsactiveimported_data
{
  MyFloat GravAccel[3];
#ifdef EVALPOTENTIAL
  MyFloat Potential;
#endif /* #ifdef EVALPOTENTIAL */
  int index;
} * Tree_ResultsActiveImported;

extern char ParameterFile[MAXLEN_PATH]; /*!< file name of parameterfile used for starting the simulation */

extern FILE *FdInfo, /*!< file handle for info.txt log-file. */
    *FdEnergy,       /*!< file handle for energy.txt log-file. */
    *FdTimings,      /*!< file handle for timings.txt log-file. */
    *FdBalance,      /*!< file handle for balance.txt log-file. */
    *FdTimebin,      /*!< file handle for timebins.txt log-file. */
    *FdDomain,       /*!< file handle for domain.txt log-file. */
    *FdMemory,       /*!< file handle for memory.txt log-file. */
    *FdCPU;          /*!< file handle for cpu.txt log-file. */

#ifdef DETAILEDTIMINGS
extern FILE *FdDetailed;
#endif /* #ifdef DETAILEDTIMINGS */

#ifdef OUTPUT_CPU_CSV
extern FILE *FdCPUCSV; /**< file handle for cpu.csv log-file. Used if the cpu log is printed in csv format as well. */
#endif                 /* #ifdef OUTPUT_CPU_CSV */

#ifdef RESTART_DEBUG
extern FILE *FdRestartTest;
#endif /* #ifdef RESTART_DEBUG */

#ifdef USE_SFR
extern FILE *FdSfr; /**< file handle for sfr.txt log-file. */
#endif              /* #ifdef USE_SFR */

#ifdef FORCETEST
extern FILE *FdForceTest; /*!< file handle for forcetest.txt log-file. */
#endif                    /* #ifdef FORCETEST */

/*! Determines whether various dump files are written. Normally true,
    set to false by Sunrise to avoid creating them. */
extern int WriteMiscFiles;

extern void *CommBuffer; /*!< points to communication buffer, which is used at a few places */

/*! \brief Global simulation data.
 *
 *  Data which is the SAME for all tasks (mostly code parameters read
 *  from the parameter file).  Holding this data in a structure is
 *  convenient for writing/reading the restart file, and it allows the
 *  introduction of new global variables in a simple way. The only
 *  thing to do is to introduce them into this structure.
 */
extern struct global_data_all_processes
{
  long long TotNumPart; /*!<  total particle numbers (global value) */
  long long TotNumGas;  /*!<  total gas particle number (global value) */

  int MaxPart;    /*!< This gives the maxmimum number of particles that can be stored on one
                     processor. */
  int MaxPartSph; /*!< This gives the maxmimum number of SPH particles that can be stored on one
                     processor. */

#if defined(COOLING)
  char TreecoolFile[MAXLEN_PATH];
#endif /* #if defined(COOLING) */

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
  int TotPartSpecial, MaxPartSpecial;
#endif /* #ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE */

#if defined(REFINEMENT)
  double ReferenceGasPartMass;
#endif /* #if defined(REFINEMENT) */

#ifdef REFINEMENT
  double TargetGasMass;
  double TargetGasMassFactor;
  int RefinementCriterion;
  int DerefinementCriterion;
#endif /* #ifdef REFINEMENT */

  double TotGravCost;

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  double AvgType1Mass;
#endif /* #ifdef INDIVIDUAL_GRAVITY_SOFTENING */

  double MeanVolume;

  int MultipleDomains;
  double TopNodeFactor;

  int ICFormat; /*!< selects different versions of IC file-format */

  int SnapFormat; /*!< selects different versions of snapshot file-formats */

  int NumFilesPerSnapshot;       /*!< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel; /*!< maximum number of files that may be written/read simultaneously when
                                    writing/reading restart-files, or when writing snapshot files */

  double TreeAllocFactor; /*!< Each processor allocates a number of nodes which is TreeAllocFactor times
                             the maximum(!) number of particles.  Note: A typical local tree for N
                             particles needs usually about ~0.65*N nodes. */

  double TopNodeAllocFactor; /*!< Each processor allocates a number of nodes which is TreeAllocFactor times
                                the maximum(!) number of particles.  Note: A typical local tree for N
                                particles needs usually about ~0.65*N nodes. */

  double NgbTreeAllocFactor; /*!< Each processor allocates a number of nodes for the neighbor search which is NgbTreeAllocFactor times
                                 the maximum(!) number of gas particles.  Note: A typical local tree for N
                                 particles needs usually about ~0.65*N nodes. */

  int MaxMemSize; /*!< size of maximum memory consumption in MB */

  /* some SPH parameters */

  int DesNumNgb; /*!< Desired number of SPH neighbours */

#ifdef SUBFIND
  int DesLinkNgb;
  double ErrTolThetaSubfind;
#endif /* #ifdef SUBFIND */

  double TotCountReducedFluxes;
  double TotCountFluxes;

  double DtDisplacement;

  double MaxNumNgbDeviation; /*!< Maximum allowed deviation neighbour number */

  double InitGasTemp; /*!< may be used to set the temperature in the IC's */
  double InitGasU;    /*!< the same, but converted to thermal energy per unit mass */
  double MinGasTemp;  /*!< may be used to set a floor for the gas temperature */
  double MinEgySpec;  /*!< the minimum allowed temperature expressed as energy per unit mass; code will inject energy if a cell falls
                         below this limit */

  double MinimumDensityOnStartUp;

  double GasSoftFactor;

  double LimitUBelowThisDensity;
  double LimitUBelowCertainDensityToThisValue;

  /* some force counters  */
  long long TotNumOfForces; /*!< counts total number of force computations  */

#ifdef MULTIPLE_RESTARTS
  int RestartFileCount;
#endif /* #ifdef MULTIPLE_RESTARTS */

  /* various cosmological factors that are only a function of the current scale factor, and in non-comoving runs are set to 1 */
  double cf_atime, cf_a2inv, cf_a3inv, cf_afac1, cf_afac2, cf_afac3, cf_hubble_a, cf_time_hubble_a, cf_redshift;
  /* Hubble rate at the current time, valid both for comoving and non-comoving integration */
  double cf_H;
  /* Hubble expansion rate, but in non-comoving integration set to zero */
  double cf_Hrate;

  /* system of units  */
  double UnitTime_in_s,         /*!< factor to convert internal time unit to seconds/h */
      UnitMass_in_g,            /*!< factor to convert internal mass unit to grams/h */
      UnitVelocity_in_cm_per_s, /*!< factor to convert internal velocity unit to cm/sec */
      UnitLength_in_cm,         /*!< factor to convert internal length unit to cm/h */
      UnitPressure_in_cgs,      /*!< factor to convert internal pressure unit to cgs units (little 'h' still
                                   around!) */
      UnitDensity_in_cgs,       /*!< factor to convert internal mass density unit to g/cm^3*h^2 */
      UnitCoolingRate_in_cgs,   /*!< factor to convert internal cooling rate to cgs units */
      UnitEnergy_in_cgs,        /*!< factor to convert internal energy to cgs units */
      UnitTime_in_Megayears,    /*!< factor to convert internal time to megayears/h */
      GravityConstantInternal,  /*!< If set to zero in the parameterfile, the internal value of the
                                   gravitational constant is set to the Newtonian value based on the system of
                                   units specified. Otherwise the value provided is taken as internal gravity
                                   constant G. */
      G;                        /*!< Gravity-constant in internal units */

  /* Cosmology */

  double Hubble;   /*!< Hubble-constant in internal units */
  double Omega0,   /*!< matter density in units of the critical density (at z=0) */
      OmegaLambda, /*!< vaccum energy density relative to crictical density (at z=0) */
      OmegaBaryon, /*!< baryon density in units of the critical density (at z=0) */
      HubbleParam; /*!< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute
                    * physical values for cooling physics
                    */

  double BoxSize; /*!< Boxsize in case periodic boundary conditions are used */

  /* Code options */

  int ComovingIntegrationOn;   /*!< flags that comoving integration is enabled */
  int PeriodicBoundariesOn;    /*!< flags that periodic boundaries are enabled for gravity */
  int ResubmitOn;              /*!< flags that automatic resubmission of job to queue system is enabled */
  int TypeOfOpeningCriterion;  /*!< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative
                                  criterion */
  int TypeOfTimestepCriterion; /*!< gives type of timestep criterion (only 0 supported right now - unlike
                                  gadget-1.1) */
  int OutputListOn;            /*!< flags that output times are listed in a specified file */
  int CoolingOn;               /*!< flags that cooling is enabled */
  int StarformationOn;         /*!< flags that star formation is enabled */

  int NParameters;

  int LowestActiveTimeBin;
  int HighestActiveTimeBin;
  int LowestOccupiedTimeBin;
  int HighestOccupiedTimeBin;
  int LowestOccupiedGravTimeBin;
  int HighestOccupiedGravTimeBin;
  int HighestSynchronizedTimeBin;
  int SmallestTimeBinWithDomainDecomposition;
  double ActivePartFracForNewDomainDecomp;

  /* parameters determining output frequency */

  int SnapshotFileCount;     /*!< number of snapshot that is written next */
  double TimeBetSnapshot,    /*!< simulation time interval between snapshot files */
      TimeOfFirstSnapshot,   /*!< simulation time of first snapshot files */
      CpuTimeBetRestartFile, /*!< cpu-time between regularly generated restart files */
      TimeLastRestartFile,   /*!< cpu-time when last restart-file was written */
      TimeBetStatistics,     /*!< simulation time interval between computations of energy statistics */
      TimeLastStatistics;    /*!< simulation time when the energy statistics was computed the last time */
  int NumCurrentTiStep;      /*!< counts the number of system steps taken up to this point */

  /* Current time of the simulation, global step, and end of simulation */

  double Time,   /*!< current time of the simulation */
      TimeBegin, /*!< time of initial conditions of the simulation */
      TimeStep,  /*!< difference between current times of previous and current timestep */
      TimeMax;   /*!< marks the point of time until the simulation is to be evolved */

  /* variables for organizing discrete timeline */

  double Timebase_interval; /*!< factor to convert from floating point time interval to integer timeline */
  integertime Ti_Current;   /*!< current time on integer timeline */
  integertime Previous_Ti_Current;
  integertime Ti_nextoutput; /*!< next output time on integer timeline */
  integertime Ti_lastoutput;

  integertime Ti_begstep[TIMEBINS]; /*!< marks start of current step of each timebin on integer timeline */

#ifdef PMGRID
  integertime PM_Ti_endstep, PM_Ti_begstep;
  double Asmth[2], Rcut[2];
  double Corner[2][3], UpperCorner[2][3], Xmintot[2][3], Xmaxtot[2][3];
  double TotalMeshSize[2];
#if defined(EVALPOTENTIAL) && defined(PMGRID) && !defined(GRAVITY_NOT_PERIODIC)
  double MassPMregions[2];
#endif /* #if defined(EVALPOTENTIAL) && defined(PMGRID) && !defined(GRAVITY_NOT_PERIODIC) */
#endif /* #ifdef PMGRID */

  long long GlobalNSynchronizedHydro;
  long long GlobalNSynchronizedGravity;

  int LevelToTimeBin[GRAVCOSTLEVELS];
  int LevelHasBeenMeasured[GRAVCOSTLEVELS];

  /* variables that keep track of cumulative CPU consumption */

  double TimeLimitCPU;
  double CPU_Sum[CPU_LAST]; /*!< sums wallclock time/CPU consumption in whole run */

  /* tree code opening criterion */

  double ErrTolTheta;    /*!< BH tree opening angle */
  double ErrTolForceAcc; /*!< parameter for relative opening criterion in tree walk */

  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy; /*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
                               timesteps is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep, /*!< minimum allowed timestep. Normally, the simulation terminates if the
                             timestep determined by the timestep criteria falls below this limit. */
      MaxSizeTimestep;    /*!< maximum allowed timestep */

#ifdef TIMESTEP_OUTPUT_LIMIT
  double TimestepOutputLimit;
#endif /* #ifdef TIMESTEP_OUTPUT_LIMIT */

#ifdef FORCE_EQUAL_TIMESTEPS
  integertime GlobalTimeStep;
#endif /* #ifdef FORCE_EQUAL_TIMESTEPS */

  double IsoSoundSpeed;

  double CourantFac; /*!< Hydrodynamics-Courant factor */

#ifdef REGULARIZE_MESH_FACE_ANGLE
  double CellMaxAngleFactor;
#else  /* #ifdef REGULARIZE_MESH_FACE_ANGLE */
  double CellShapingFactor;
#endif /* #ifdef REGULARIZE_MESH_FACE_ANGLE #else */
  double CellShapingSpeed;

  int CPU_TimeBinCountMeasurements[TIMEBINS];
  double CPU_TimeBinMeasurements[TIMEBINS][NUMBER_OF_MEASUREMENTS_TO_RECORD];

  /* gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening
   * length)
   *
   */

  int SofteningTypeOfPartType[NTYPES];

  double SofteningComoving[NSOFTTYPES]; /*!< comoving gravitational softening lengths for each softeniung type */
  double SofteningMaxPhys[NSOFTTYPES];  /*!< maximum physical gravitational softening lengths for each softening type */

  double
      SofteningTable[NSOFTTYPES + NSOFTTYPES_HYDRO]; /*!< current (comoving) gravitational softening lengths for each softening type */
  double ForceSoftening[NSOFTTYPES + NSOFTTYPES_HYDRO + 1]; /*!<  current (comoving) gravitational softening lengths, multiplied by a
                                                               factor 2.8 - at that scale the force is Newtonian */

  /*! If particle masses are all equal for one type, the corresponding entry in MassTable is set to this
   *  value, * allowing the size of the snapshot files to be reduced
   */
  double MassTable[NTYPES];

#ifdef ADAPTIVE_HYDRO_SOFTENING
  double MinimumComovingHydroSoftening;
  double AdaptiveHydroSofteningSpacing;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */

  /* some filenames */
  char InitCondFile[MAXLEN_PATH], OutputDir[MAXLEN_PATH], SnapshotFileBase[MAXLEN_PATH], ResubmitCommand[MAXLEN_PATH],
      OutputListFilename[MAXLEN_PATH];

  /*! table with desired output times */
  double OutputListTimes[MAXLEN_OUTPUTLIST];
  char OutputListFlag[MAXLEN_OUTPUTLIST];
  int OutputListLength; /*!< number of times stored in table of desired output times */

#ifdef USE_SFR /* enable Springel & Hernquist model */
  double OverDensThresh;
  double CritOverDensity;
  double TemperatureThresh;
  double CritPhysDensity;
  double PhysDensThresh;
  double EgySpecSN;
  double EgySpecCold;
  double FactorEVP;
  double TempSupernova;
  double TempClouds;
  double MaxSfrTimescale;
  double FactorSN;
#endif /* #ifdef USE_SFR */

#ifdef MHD_POWELL
  double Powell_Momentum[3];
  double Powell_Angular_Momentum[3];
  double Powell_Energy;
#endif /* #ifdef MHD_POWELL */

#ifdef MHD_SEEDFIELD
  int B_dir;      /* flags for direction: x = 1, y = 2, z = 4 */
  double B_value; /* value for the chosen component(s) of the magnetic field */
#endif            /* #ifdef MHD_SEEDFIELD */

  MyIDType MaxID;

#ifdef REFINEMENT_VOLUME_LIMIT
  double MaxVolumeDiff;
  double MinVolume;
  double MaxVolume;
#endif /* #ifdef REFINEMENT_VOLUME_LIMIT */

#ifdef REDUCE_FLUSH
  double FlushCpuTimeDiff;
  double FlushLast;
#endif /* #ifdef REDUCE_FLUSH */

#ifdef TILE_ICS
  int TileICsFactor;
#endif /* #ifdef TILE_ICS */

#ifdef ADDBACKGROUNDGRID
  int GridSize;
#endif /* #ifdef ADDBACKGROUNDGRID */

#ifdef ONEDIMS_SPHERICAL
  double CoreMass;
  double CoreRadius;
#endif /* #ifdef ONEDIMS_SPHERICAL */

  double GlobalDisplacementVector[3];
} All;

/*****************************************************************************
 ** particle data ************************************************************
 ****************************************************************************/

/*! \brief This structure holds all the information that is
 *         stored for each particle of the simulation.
 */
extern struct particle_data
{
  MyDouble Pos[3];       /*!< particle position at its current time */
  MyDouble Mass;         /*!< particle mass */
  MyFloat Vel[3];        /*!< particle velocity at its current time */
  MySingle GravAccel[3]; /*!< particle acceleration due to gravity */

#ifdef EXTERNALGRAVITY
  MySingle dGravAccel; /*!< norm of spatial derivatives tensor of gravity accelerations due to external force */
#endif

#ifdef PMGRID
  MySingle GravPM[3]; /*!< particle acceleration due to long-range PM gravity force */
#endif                /* #ifdef PMGRID */

#ifdef FORCETEST
  MyFloat GravAccelDirect[3]; /*!< particle acceleration calculated by direct summation */
  MyFloat PotentialDirect;    /*!< potential computed with direct summation */
  MyFloat DistToID1;
#ifdef PMGRID
  MyFloat GravAccelShortRange[3]; /*!< short range component of gravitational acceleration */
  MyFloat GravAccelLongRange[3];  /*!< long range component of gravitational acceleration */
  MyFloat PotentialShortRange;    /*!< potential due to short-range forces */
  MyFloat PotentialLongRange;     /*!< potential due to long-range forces */
#endif                            /* #ifdef PMGRID */
#endif                            /* #ifdef FORCETEST  */

#if defined(EVALPOTENTIAL) || defined(OUTPUTPOTENTIAL)
  MySingle Potential; /*!< gravitational potential */
#if defined(PMGRID)
  MySingle PM_Potential; /*!< gravitational potential in Particle-Mesh */
#endif                   /* #if defined(PMGRID) */
#endif                   /* #if defined(EVALPOTENTIAL) || defined (OUTPUTPOTENTIAL) */

#ifdef OUTPUTGRAVINTERACTIONS
  int GravInteractions; /*!< number of gravitational ineractions calculated */
#endif                  /* #ifdef OUTPUTGRAVINTERACTIONS */

#ifdef EXTERNALGRAVITY
  MyFloat ExtPotential; /*!< value of external potential */
#endif                  /* #ifdef EXTERNALGRAVITY */

  MyIDType ID; /*!< unique ID of particle */

#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT)
  MyIDType FileOrder;
#endif /* #ifdefined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) */

  integertime Ti_Current; /*!< current time on integer timeline */

  float OldAcc; /*!< magnitude of old gravitational force. Used in relative opening criterion */

  float GravCost[GRAVCOSTLEVELS]; /*!< weight factors used for balancing the work-load */

  unsigned char Type; /*!< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
  unsigned char SofteningType;
  signed char TimeBinGrav;
  signed char TimeBinHydro;
} * P,              /*!< holds particle data on local processor */
    *DomainPartBuf; /*!< buffer for particle data used in domain decomposition */

/*****************************************************************************
 ** (sub)halo data ***********************************************************
 ****************************************************************************/

extern struct subfind_data
{
  int OriginIndex, OriginTask;
  int TargetIndex, TargetTask;
  int GrNr;

#ifdef SUBFIND
  int SubNr;
  int OldIndex;
  int submark;
  int originindex, origintask;
  MyFloat Utherm;
  MyFloat Density;
  MyFloat Potential;
  MyFloat Hsml;
  MyFloat BindingEnergy;

#ifdef CELL_CENTER_GRAVITY
  MyDouble Center[3];
#endif /* #ifdef CELL_CENTER_GRAVITY */

#ifdef SUBFIND_CALC_MORE
  MyFloat SubfindHsml;
  MyFloat SubfindDensity;   /* total matter density */
  MyFloat SubfindDMDensity; /* dark matter density */
  MyFloat SubfindVelDisp;   /* 3D DM velocity dispersion */
#endif                      /* #ifdef SUBFIND_CALC_MORE */

#endif /* #ifdef SUBFIND */
} * PS;

/*****************************************************************************
 ** cell data ****************************************************************
 ****************************************************************************/

/*! \brief Holds data that is stored for each hydro mesh cell in addition to
 *         the collisionless variables.
 */
extern struct sph_particle_data
{
  /* conserved variables */
  MyFloat Energy;
  MyFloat Momentum[3];
  MyFloat Volume;
  MyFloat OldMass;

  /* primitive variables */
  MyFloat Density;
  MyFloat Pressure; /*!< current pressure */
  MySingle Utherm;

#ifdef HIERARCHICAL_GRAVITY
  MySingle FullGravAccel[3];
#endif /* #ifdef HIERARCHICAL_GRAVITY */

  /* variables for mesh  */
  MyDouble Center[3];    /*!< center of mass of cell */
  MySingle VelVertex[3]; /*!< current vertex velocity (primitive variable) */

  MySingle MaxDelaunayRadius;
  MySingle Hsml; /* auxiliary search radius for points around a delaunay triangle */
  MySingle SurfaceArea;

#if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE)
  MySingle MaxFaceAngle;
#endif /* #if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE) */

  MySingle ActiveArea;

#if defined(OUTPUT_DIVVEL)
  MyFloat DivVel; /*!< divergence of the velocity field */
#endif            /* #if defined(OUTPUT_DIVVEL) */

#if defined(REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED) || defined(OUTPUT_CURLVEL)
  MySingle CurlVel; /*!< magnitude of the curl of the velocity field */
#endif              /* #if defined(REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED) || defined(OUTPUT_CURLVEL) */

#ifdef TREE_BASED_TIMESTEPS
  MySingle CurrentMaxTiStep;
  MySingle Csnd;
#endif /* #ifdef TREE_BASED_TIMESTEPS */

#if defined(REFINEMENT_HIGH_RES_GAS)
  MyFloat HighResMass;
  MyFloat HighResDensity;
#endif /* #if defined(REFINEMENT_HIGH_RES_GAS) */

#ifdef MHD
  MyFloat B[3];
  MyFloat BConserved[3];
  MyFloat DivB;
  MyFloat CurlB[3];
#endif /* #ifdef MHD */

#ifdef PASSIVE_SCALARS
  MyFloat PScalars[PASSIVE_SCALARS];
  MyFloat PConservedScalars[PASSIVE_SCALARS];
#endif /* #ifdef PASSIVE_SCALARS */

#ifdef OUTPUT_SURFACE_AREA
  int CountFaces;
#endif /* #ifdef OUTPUT_SURFACE_AREA */

#if defined(REFINEMENT_SPLIT_CELLS)
  MySingle MinimumEdgeDistance;
#endif /* #if defined(REFINEMENT_SPLIT_CELLS) */

#if defined(COOLING)
  MyFloat Ne; /* electron fraction, expressed as local electron number
                 density normalized to the hydrogen number density. Gives
                 indirectly ionization state and mean molecular weight. */
#endif        /* #if defined(COOLING) */

#ifdef USE_SFR
  MySingle Sfr;
#endif /* #ifdef USE_SFR */

#ifdef OUTPUT_COOLHEAT
  MyFloat CoolHeat;
#endif /* #ifdef OUTPUT_COOLHEAT */

  struct grad_data Grad;

  int first_connection;
  int last_connection;

#ifdef REFINEMENT_HIGH_RES_GAS
  int AllowRefinement;
#endif /* #ifdef REFINEMENT_HIGH_RES_GAS */

#ifdef REFINEMENT_SPLIT_CELLS
  MySingle SepVector[3];
#endif /* #ifdef REFINEMENT_SPLIT_CELLS */

#ifdef REFINEMENT_VOLUME_LIMIT
  MyFloat MinNgbVolume;
#endif /* #ifdef REFINEMENT_VOLUME_LIMIT */

  double TimeLastPrimUpdate;

#ifdef ADDBACKGROUNDGRID
  MyFloat Weight;
#endif /* #ifdef ADDBACKGROUNDGRID */

} * SphP,          /*!< holds SPH particle data on local processor */
    *DomainSphBuf; /*!< buffer for SPH particle data in domain decomposition */

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
extern struct special_particle_data
{
  MyIDType ID;
  double pos[3];
  double mass;
} * PartSpecialListGlobal;
#endif /* #ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE */

extern peanokey *DomainKeyBuf;

/*! global state of system
 */
extern struct state_of_system
{
  double Mass, EnergyKin, EnergyPot, EnergyInt, EnergyTot, Momentum[4], AngMomentum[4], CenterOfMass[4], MassComp[NTYPES],
      EnergyKinComp[NTYPES], EnergyPotComp[NTYPES], EnergyIntComp[NTYPES], EnergyTotComp[NTYPES], MomentumComp[NTYPES][4],
      AngMomentumComp[NTYPES][4], CenterOfMassComp[NTYPES][4];
} SysState, SysStateAtStart, SysStateAtEnd;

/*! \brief Struct used for passing the parameters during the mesh cell search.
 */
typedef struct
{
  MyDouble Pos[3];
  int Task;
  union
  {
    int Index;
    float hsmlguess;
  } u;

} mesh_search_data;

/*! \brief Struct used for sending positions to other tasks during the
 *         mesh cell search.
 */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Distance;
} mesh_search_request;

/*! \brief Struct used for receiving the results from other tasks during the
 *         mesh cell search.
 */
typedef struct
{
  MyDouble Distance;
  int Task;
  int Index;
} mesh_search_response;

extern struct data_partlist
{
  int Task;  /*!< The task the item was exported to. */
  int Index; /*!< The particle index of the item on the sending task. */
} * PartList;

extern struct datanodelist
{
  int Task;  /*!< target process */
  int Index; /*!< local index that wants to open this node */
  int Node;  /*!< node to be opened on foreign process */
} * NodeList;

#define FAC_AVG_NODES_PER_EXPORT 4.0 /*!< default choice for estimated average number of exported nodes per exported particle */

extern struct directdata
{
  MyDouble Pos[3];
  MyDouble Mass;
  unsigned char Type;
  unsigned char SofteningType;
} * DirectDataIn, *DirectDataAll;

extern struct accdata
{
  MyFloat Acc[3];
#ifdef EVALPOTENTIAL
  MyFloat Potential;
#endif /* #ifdef EVALPOTENTIAL */
} * DirectAccOut, *DirectAccIn;

#if defined(EVALPOTENTIAL) || defined(OUTPUTPOTENTIAL) || defined(SUBFIND)
extern struct potdata_out
{
  MyFloat Potential;
}
    /*! \brief Holds the partial results computed for imported particles. Note:
     *         We use GravDataResult = GravDataGet, such that the result replaces
     *         the imported data
     */
    * PotDataResult,
    /*! \brief Holds partial results received from other processors. This will
     *         overwrite the GravDataIn array
     */
    *PotDataOut;
#endif /* #if defined (EVALPOTENTIAL) || defined (OUTPUTPOTENTIAL) || defined(SUBFIND) */

/*! \brief Buffer of size NTask used for flagging whether a particle needs to
 *         be exported to the other tasks.
 */
extern int *Exportflag;
/*! \brief Buffer of size NTask used for counting how many nodes are to be
 *         exported to the other tasks?
 */
extern int *Exportnodecount;
/*! \brief Buffer of size NTask used for holding the index into the
 *         DataIndexTable.
 */
extern int *Exportindex;
/*! \brief Array of NTask size of the offset into the send array where the
 *         objects to be sent to the specified task starts.
 */
extern int *Send_offset,
    /*! \brief Array of NTask size of the number of objects to send to the
     *  tasks.
     */
    *Send_count,
    /*! \brief Array of NTask size of the number of objects to receive from the
     *         tasks.
     */
    *Recv_count,
    /*! \brief Array of NTask size of the offset into the receive array where the
     *         objects from the specified task starts.
     */
    *Recv_offset;

extern int *TasksThatSend, *TasksThatRecv, NSendTasks, NRecvTasks;

extern struct send_recv_counts
{
  int Count;
  int CountNodes;
} * Send, *Recv;

extern int *Send_offset_nodes, *Send_count_nodes, *Recv_count_nodes, *Recv_offset_nodes;

extern int Mesh_nimport, Mesh_nexport, *Mesh_Send_offset, *Mesh_Send_count, *Mesh_Recv_count, *Mesh_Recv_offset;

extern int Force_nimport, Force_nexport, *Force_Send_offset, *Force_Send_count, *Force_Recv_count, *Force_Recv_offset;

/*! \brief Header for the standard file format.
 */
#if(NTYPES == 7 || NTYPES == 8)
#define NTYPES_INT_HEADER 8
#else /* #if (NTYPES==7 || NTYPES==8) */
#define NTYPES_INT_HEADER NTYPES
#endif /* #if (NTYPES==7 || NTYPES==8) #else */
extern struct io_header
{
  int npart[NTYPES_INT_HEADER];                       /*!< number of particles of each type in this file */
  double mass[NTYPES];                                /*!< mass of particles of each type. If 0, then the masses are explicitly
                                                         stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                                        /*!< time of snapshot file */
  double redshift;                                    /*!< redshift of snapshot file */
  int flag_sfr;                                       /*!< flags whether the simulation was including star formation */
  int flag_feedback;                                  /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[NTYPES_INT_HEADER];         /*!< total number of particles of each type in this snapshot. This can be
                                         different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;                                   /*!< flags whether cooling was included  */
  int num_files;                                      /*!< number of files in multi-file snapshot */
  double BoxSize;                                     /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;                                      /*!< matter density in units of critical density */
  double OmegaLambda;                                 /*!< cosmological constant parameter */
  double HubbleParam;                                 /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;                                /*!< flags whether the file contains formation times of star particles */
  int flag_metals;                                    /*!< flags whether the file contains metallicity values for gas and star
                                                         particles */
  unsigned int npartTotalHighWord[NTYPES_INT_HEADER]; /*!< High word of the total number of particles of each type */
  int flag_entropy_instead_u;                         /*!< flags that IC-file contains entropy instead of u */
  int flag_doubleprecision;                           /*!< flags that snapshot contains double-precision instead of single precision */

  int flag_lpt_ics;        /*!< flag to signal that IC file contains 2lpt initial conditions */
  float lpt_scalingfactor; /*!< scaling factor for 2lpt initial conditions */

  int flag_tracer_field; /*!< flags presence of a tracer field */

  int composition_vector_length; /*!< specifies the length of the composition vector (0 if not present)  */

#if(NTYPES == 6)
  char fill[40];   /*!< fills to 256 Bytes */
#elif(NTYPES == 7) /* #if (NTYPES==6) */
  char fill[8]; /*!< fills to 256 Bytes */
#endif             /* #elif (NTYPES==7) */
} header;          /*!< holds header for snapshot files */

/*! \brief Header for the ICs file format, if NTYPES does not match.
 */
#ifdef NTYPES_ICS
extern struct io_header_ICs
{
  int npart[NTYPES_ICS];                       /*!< number of particles of each type in this file */
  double mass[NTYPES_ICS];                     /*!< mass of particles of each type. If 0, then the masses are explicitly
                                                  stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                                 /*!< time of snapshot file */
  double redshift;                             /*!< redshift of snapshot file */
  int flag_sfr;                                /*!< flags whether the simulation was including star formation */
  int flag_feedback;                           /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[NTYPES_ICS];         /*!< total number of particles of each type in this snapshot. This can be
                                          different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;                            /*!< flags whether cooling was included  */
  int num_files;                               /*!< number of files in multi-file snapshot */
  double BoxSize;                              /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;                               /*!< matter density in units of critical density */
  double OmegaLambda;                          /*!< cosmological constant parameter */
  double HubbleParam;                          /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;                         /*!< flags whether the file contains formation times of star particles */
  int flag_metals;                             /*!< flags whether the file contains metallicity values for gas and star
                                                  particles */
  unsigned int npartTotalHighWord[NTYPES_ICS]; /*!< High word of the total number of particles of each type */
  int flag_entropy_instead_u;                  /*!< flags that IC-file contains entropy instead of u */
  int flag_doubleprecision;                    /*!< flags that snapshot contains double-precision instead of single precision */

  int flag_lpt_ics;        /*!< flag to signal that IC file contains 2lpt initial conditions */
  float lpt_scalingfactor; /*!< scaling factor for 2lpt initial conditions */

  int flag_tracer_field; /*!< flags presence of a tracer field */

  int composition_vector_length; /*!< specifies the length of the composition vector (0 if not present)  */

#if(NTYPES_ICS == 6)
  char fill[40]; /*!< fills to 256 Bytes */
#else            /* #if (NTYPES_ICS==6) */
  terminate("NTYPES_ICS != 6")
#endif           /* #if (NTYPES_ICS==6) #else */
} header_ICs;    /*!< holds header for IC files */
#endif           /* #ifdef NTYPES_ICS */

enum iofields
{
  IO_POS,
  IO_VEL,
  IO_ID,
  IO_MASS,
  IO_U,
  IO_RHO,
  IO_VORT,
  IO_VOL,
  IO_CM,
  IO_VERTEXVEL,
  IO_FACEANGLE,
  IO_SAREA,
  IO_NFACES,

  IO_HIGHRESMASS,
  IO_PRESSURE,
  IO_CSND,
  IO_NE,
  IO_NH,
  IO_SFR,

  IO_POT,
  IO_ACCEL,
  IO_GRADP,
  IO_GRADR,
  IO_GRADV,
  IO_GRADB,

  IO_POT_MINI,
  IO_POS_MINI,

  IO_HI,
  IO_TSTP,
  IO_BFLD,
  IO_DIVB,
  IO_COOLRATE,
  IO_ALLOWREFINEMENT,

  IO_DIVVEL,
  IO_CURLVEL,
  IO_COOLHEAT,
  IO_PASS,

  IO_SUBFINDHSML,
  IO_SUBFINDDENSITY,
  IO_SUBFINDDMDENSITY,
  IO_SUBFINDVELDISP,
  IO_GROUPNR,

  IO_SOFTENING,
  IO_TASK,
  IO_TIMEBIN_HYDRO,

  IO_LASTENTRY /* This should be kept - it signals the end of the list */
};

enum arrays
{
  A_NONE,
  A_SPHP,
  A_P,
  A_PS
};

enum types_in_file
{
  FILE_NONE        = -1,
  FILE_INT         = 0,
  FILE_MY_ID_TYPE  = 2,
  FILE_MY_IO_FLOAT = 1,
  FILE_DOUBLE      = 3,
  FILE_FLOAT       = 4
};

enum types_in_memory
{
  MEM_INT,
  MEM_MY_ID_TYPE,
  MEM_FLOAT,
  MEM_DOUBLE,
  MEM_MY_SINGLE,
  MEM_MY_FLOAT,
  MEM_MY_DOUBLE,
  MEM_NONE
};

enum e_typelist
{
  GAS_ONLY                      = 1,
  STARS_ONLY                    = 16,
  GAS_AND_STARS                 = 17,
  BHS_ONLY                      = 32,
  ALL_TYPES                     = ((1 << NTYPES) - 1),
  SET_IN_GET_PARTICLES_IN_BLOCK = 0
};

enum sn_type
{
  SN_FULL      = 0,
  SN_MINI      = 1,
  SN_MINI_ONLY = 2,
  SN_NO_SUBBOX = 3
};

typedef struct
{
  enum iofields field;
  enum types_in_memory type_in_memory;
  enum types_in_file type_in_file_input;
  enum types_in_file type_in_file_output;
  int values_per_block;
  char label[4];
  char datasetname[256];
  void (*io_func)(int, int, void *, int);
  int typelist;
  enum arrays array;
  size_t offset;
  enum sn_type snap_type;

  char hasunit;
  double a;
  double h;
  double L;
  double M;
  double V;
  double c;
} IO_Field;

extern IO_Field *IO_Fields;
extern int N_IO_Fields;
extern int Max_IO_Fields;

extern char (*Parameters)[MAXLEN_PARAM_TAG];
extern char (*ParametersValue)[MAXLEN_PARAM_VALUE];
extern char *ParametersType;

/*! \brief The tree data structure.
 *
 *  Nodes points to the actual memory
 *  allocated for the internal nodes, but is shifted such that
 *  Nodes[All.MaxPart] gives the first allocated node. Note that node
 *  numbers less than All.MaxPart are the leaf nodes that contain a
 *  single particle, and node numbers >= MaxPart+MaxNodes are "pseudo
 *  particles" that hang off the toplevel leaf nodes belonging to
 *  other tasks. These are not represented by this structure. Instead,
 *  the tree traversal for these are saved in the Nextnode, Prevnode
 *  and Father arrays, indexed with the node number in the case of
 *  real particles and by nodenumber-MaxNodes for pseudo
 *  particles.
 */
extern struct NODE
{
  union
  {
    int suns[8]; /*!< temporary pointers to daughter nodes */
    struct
    {
      MyDouble s[3]; /*!< center of mass of node */
      MyDouble mass; /*!< mass of node */
      /*! The next node in the tree walk in case the current node does
       *  not need to be opened. This means that it traverses the 8
       *  subnodes of a node in a breadth-first fashion, and then goes
       *  to father->sibling.
       */
      int sibling;
      /*! The next node in case the current node needs to be
       *  opened. Applying nextnode repeatedly results in a pure
       *  depth-first traversal of the tree.
       */
      int nextnode;
      /*! The parent node of the node. (Is -1 for the root node.)
       */
      int father;
#if(NSOFTTYPES > 1)
      unsigned char maxsofttype; /**< hold the maximum gravitational softening of particles */
#if defined(MULTIPLE_NODE_SOFTENING) && defined(ADAPTIVE_HYDRO_SOFTENING)
      unsigned char maxhydrosofttype;
      unsigned char minhydrosofttype;
#endif /* #if defined(MULTIPLE_NODE_SOFTENING) && defined(ADAPTIVE_HYDRO_SOFTENING) */
#endif /* #if (NSOFTTYPES > 1) */
    } d;
  } u;

  MyDouble center[3]; /*!< geometrical center of node */
  MyFloat len;        /*!< sidelength of treenode */

} * Nodes;

#ifdef MULTIPLE_NODE_SOFTENING
extern struct ExtNODE
{
  MyDouble mass_per_type[NSOFTTYPES];
} * ExtNodes;
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

/*! Gives next node in tree walk for the "particle" nodes. Entries 0
 *  -- MaxPart-1 are the real particles, and the "pseudoparticles" are
 *  indexed by the node number-MaxNodes.
 */
extern int *Nextnode;

/*! Gives previous node in tree walk for the leaf (particle)
 *  nodes. Entries 0 -- MaxPart-1 are the real particles, and the
 *  "pseudoparticles" are indexed by the node number-MaxNodes.
 */
extern int *Father;

/*! Variables for neighbor tree */
extern int Ngb_MaxPart;
extern int Ngb_NumNodes;
extern int Ngb_MaxNodes;
extern int Ngb_FirstNonTopLevelNode;
extern int Ngb_NextFreeNode;
extern int *Ngb_Father;
extern int *Ngb_Marker;
extern int Ngb_MarkerValue;

extern int *Ngb_DomainNodeIndex;
extern int *DomainListOfLocalTopleaves;
extern int *DomainNLocalTopleave;
extern int *DomainFirstLocTopleave;
extern int *Ngb_Nextnode;

/*! The ngb-tree data structure
 */
extern struct NgbNODE
{
  union
  {
    int suns[8]; /*!< temporary pointers to daughter nodes */
    struct
    {
      int sibling;
      int nextnode;
      MyNgbTreeFloat range_min[3];
      MyNgbTreeFloat range_max[3];
    } d;
  } u;

  MyNgbTreeFloat vertex_vmin[3];
  MyNgbTreeFloat vertex_vmax[3];

  int father;

  integertime Ti_Current;

} * Ngb_Nodes;

extern struct ExtNgbNODE
{
  float vmin[3];
  float vmax[3];
  float MaxCsnd;
} * ExtNgb_Nodes;

#ifdef STATICNFW
extern double Rs, R200;
extern double Dc;
extern double RhoCrit, V200;
extern double fac;
#endif /* #ifdef STATICNFW */

extern int MaxThreads;

#endif /* #define ALLVARS_H */
