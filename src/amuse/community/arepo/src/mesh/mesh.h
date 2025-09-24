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
 * \file        src/mesh/mesh.h
 * \date        05/2018
 * \brief       Header for mesh structures.
 * \details
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 29.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifndef MESH_H
#define MESH_H

#define SCALAR_TYPE_PASSIVE 0   /*!< only advection */
#define SCALAR_TYPE_SPECIES 1   /*!< species are normalised to guarantee sum{species}=1 */
#define SCALAR_TYPE_NORMALIZE 2 /*!< the same normalisation factor as for species is applied, but no contribution to sum{species} */

#define REFL_X_FLAGS 115043766
#define REFL_Y_FLAGS 132379128
#define REFL_Z_FLAGS 134217216

#define OUTFLOW_X (1 << 27)
#define OUTFLOW_Y (1 << 28)
#define OUTFLOW_Z (1 << 29)

#if defined MAXSCALARS
extern struct scalar_elements
{
  int type;           /*!< scalar type, determines whether a normalization is applied */
  size_t offset;      /*!< offset of the primitive quantity in the SphP struct */
  size_t offset_mass; /*!< offset of the conserved quantity in the SphP struct */
} scalar_elements[MAXSCALARS];

extern struct scalar_index
{
#ifdef REFINEMENT_HIGH_RES_GAS
  int HighResMass;
#endif /* #ifdef REFINEMENT_HIGH_RES_GAS */
} ScalarIndex;

extern int N_Scalar; /*!< number of registered scalars */
#endif               /* #if defined MAXSCALARS */

#define GRADIENT_TYPE_NORMAL 0
#define GRADIENT_TYPE_VELX 1
#define GRADIENT_TYPE_VELY 2
#define GRADIENT_TYPE_VELZ 3
#define GRADIENT_TYPE_DENSITY 4
#define GRADIENT_TYPE_PRESSURE 5
#define GRADIENT_TYPE_UTHERM 6
#define GRADIENT_TYPE_AX 7
#define GRADIENT_TYPE_AY 8
#define GRADIENT_TYPE_AZ 9
#define GRADIENT_TYPE_FLD 10
#define GRADIENT_TYPE_RTF 11

extern struct grad_elements
{
  int type;           /*!< gradient type, ensures special treatment for velocities and speed of sound */
  size_t offset;      /*!< offset of the quantity in the SphP struct */
  size_t offset_exch; /*!< offset of the quantity in the PrimExch struct */
  size_t offset_grad; /*!< offset in the grad_data struct */
  double *min_value, *max_value;
  double value0, value1;
} grad_elements[MAXGRADIENTS], *GDensity, *GVelx, *GVely, *GVelz, *GPressure, *GUtherm;

extern int N_Grad; /*!< number of gradients to be calculated */

extern struct grad_data
{
  MySingle drho[3];

  MySingle dvel[3][3];
  MySingle dpress[3];

#ifdef MHD
  MySingle dB[3][3];
#endif /* #ifdef MHD */

#ifdef MAXSCALARS
  MySingle dscalars[MAXSCALARS][3];
#endif /* #ifdef MAXSCALARS */
} * GradExch;

extern struct primexch
{
  double Volume;
  MyFloat Density;

  MyFloat VelGas[3];
  MyFloat VelVertex[3];

#ifdef MHD
  MyFloat B[3];

#ifdef MHD_POWELL
  MyFloat DivB;
#endif /* #ifdef MHD_POWELL */

  MyFloat CurlB[3];
#endif /* #ifdef MHD */
  MyFloat Pressure;

#ifdef MAXSCALARS
  MyFloat Scalars[MAXSCALARS];
#endif /* #ifdef MAXSCALARS */

  double TimeLastPrimUpdate;

  MyDouble Center[3];
  MyFloat OldMass;
  MySingle Csnd;
  MySingle SurfaceArea;
  MySingle ActiveArea;
  /*  int task, index; */
  short int TimeBinHydro;
} * PrimExch;

#ifdef REFINEMENT
extern struct refdata
{
#ifdef REFINEMENT_VOLUME_LIMIT
  double Volume;
#endif /* #ifdef REFINEMENT_VOLUME_LIMIT */
  short int TimeBinHydro;
} * RefExch;
#endif /* #ifdef REFINEMENT */

typedef struct face_data
{
  int p1, p2;
#ifdef REFINEMENT_MERGE_CELLS
  int t, nr; /* delaunay tetra and edge number that generated this face */
#endif       /* #ifdef REFINEMENT_MERGE_CELLS */

#ifdef OPTIMIZE_MEMORY_USAGE
  MyFloat area;
  MyFloat cx, cy, cz; /* center-of-mass of face */
#else                 /* #ifdef OPTIMIZE_MEMORY_USAGE */
  double area;
  double cx, cy, cz; /* center-of-mass of face */
#endif                /* #ifdef OPTIMIZE_MEMORY_USAGE #else */

#ifdef VORONOI_BACKUP_RESTORE_FACE_AREAS
  double area_backup;
#endif /* #ifdef VORONOI_BACKUP_RESTORE_FACE_AREAS */
#ifdef TETRA_INDEX_IN_FACE
  int dt_index;
#endif /* #ifdef TETRA_INDEX_IN_FACE */
} face;

/*! left or right state of a face */
struct state
{
  double dx, dy, dz;
  double dt_half;
  short int timeBin;

  double rho;
  double velx, vely, velz;
  double press;
  double oldmass;
  double surfacearea;
  double activearea;
  double volume;

  MyFloat velGas[3];
  MyFloat velVertex[3];
  struct grad_data *grad;

  double csnd;
  double Energy;
#ifdef MHD
  double Bx, By, Bz;
#ifdef MHD_POWELL
  double divB;
#endif /* #ifdef MHD_POWELL */
  double CurlB[3];
#endif /* #ifdef MHD */

#if defined(GODUNOV_STATS)
  double mach;
#endif /* #if defined(GODUNOV_STATS) */

#ifdef MAXSCALARS
  double scalars[MAXSCALARS];
#endif /* #ifdef MAXSCALARS */
  MyIDType ID;

#ifdef ONEDIMS_SPHERICAL
  double radius;
#endif /* #ifdef ONEDIMS_SPHERICAL */

  double dtExtrapolation;
};

/*! state on a face determined by riemann solver */
extern struct state_face
{
  double rho;
  double velx, vely, velz;
  double press;
#ifdef MHD
  double Bx, By, Bz;
#endif /* #ifdef MHD */

#ifdef MAXSCALARS
  double *scalars;
#endif /* #ifdef MAXSCALARS */
} state_face;

/*! flux through a face */
extern struct fluxes
{
  double mass;
  double momentum[3];
  double energy;

#ifdef MHD
  double B[3];
#endif /* #ifdef MHD */

#ifdef MAXSCALARS
  double scalars[MAXSCALARS];
#endif /* #ifdef MAXSCALARS */
} fluxes, diffusionfluxes;

extern struct geometry
{
  double nn;
  double nx, ny, nz;
  double mx, my, mz;
  double px, py, pz;
  double cx, cy, cz;
} geom;

struct pv_update_data
{
  double atime;
  double hubble_a;
  double a3inv;
};
#endif /* MESH_H */

struct fvs_stat
{
  int count_disable_extrapolation;
};
