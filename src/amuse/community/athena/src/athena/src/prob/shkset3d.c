#include "copyright.h"
/*============================================================================*/
/*! \file shkset3d.c
 *  \brief Sets up shock at angle to grid to test 3D algorithms.
 *
 * PURPOSE: Sets up shock at angle to grid to test 3D algorithms.  Setting up
 *   the initial conditions to minimize grid noise is very complex, so this
 *   special function has been written to handle only 3D problems.
 *
 * This problem cannot be run with Static Mesh Refinement.
 *
 * This code is most easily understood in terms of a one dimensional
 * problem in the coordinate system (x1,x2,x3).  Two coordinate rotations are
 * applied to obtain a new wave vector in a 3D space in the (x,y,z)
 * coordinate system.
 *
 *   First rotate about the x2 axis:
 *  -  x1' = x1*cos(ang_2) - x3*sin(ang_2)
 *  -  x2' = x2
 *  -  x3' = x1*sin(ang_2) + x3*cos(ang_2)
 *
 *   Next rotate about the x3' axis:
 *   - x = x1'*cos(ang_3) - x2'*sin(ang_3)
 *   - y = x1'*sin(ang_3) + x2'*cos(ang_3)
 *   - z = x3'
 *
 *   Expanding this out we get:
 *   - x = x1*cos(ang_2)*cos(ang_3) - x2*sin(ang_3) - x3*sin(ang_2)*cos(ang_3)
 *   - y = x1*cos(ang_2)*sin(ang_3) + x2*cos(ang_3) - x3*sin(ang_2)*sin(ang_3)
 *   - z = x1*sin(ang_2)                            + x3*cos(ang_2)
 *
 *   This inverts to:
 *   - x1 =  x*cos(ang_2)*cos(ang_3) + y*cos(ang_2)*sin(ang_3) + z*sin(ang_2)
 *   - x2 = -x*sin(ang_3)            + y*cos(ang_3)
 *   - x3 = -x*sin(ang_2)*cos(ang_3) - y*sin(ang_2)*sin(ang_3) + z*cos(ang_2)
 * Note these transformations are the same used in linear_wave3d and cpaw3d
 *
 * The initial conditions must be translation invariant, i.e. 
 * q(x,y,z) = q(x+tx, y+ty, z+tz) for (tx,ty,tz) such that 
 * x1(x,y,z) = x1(x+tx, y+ty, z+ty).  Inserting this last expression for
 * x1(x,y,z) into the transformation above leads to
 *
 * - tx*cos(ang_2)*cos(ang_3) + ty*cos(ang_2)*sin(ang_3) + tz*sin(ang_2) = 0.
 *
 * Now, restrict the translation symmetry to be discrete by inserting
 * (tx, ty, tz) = (nx*dx, ny*dy, nz*dz) where (nx, ny, nz) are integers (not 
 * the number of grid cells in each direction) and (dx, dy, dz) are the grid
 * cell sizes.  With some simplification the translation symmetry becomes
 *
 * - nx + ny*tan(ang_3)*(dy/dx) + nz*(tan(ang_2)/cos(ang_3))*(dz/dx) = 0.
 *
 * Now, choose an integer unit cell size (rx, ry, rz) where
 *
 * - tan(ang_3)*(dy/dx) = (rx/ry)
 *
 * and
 *
 * - (tan(ang_2)/cos(ang_3))*(dz/dx) = (rx/rz).
 *
 * With this choice, or equation for discrete translation symmetry becomes
 *
 * - (nx/rx) + (ny/ry) + (nz/rz) = 0.
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - lx_bc() - use periodicity of "unit cell" to apply BCs to left-x edge
 * - rx_bc() - use periodicity of "unit cell" to apply BCs to right-x edge
 * - ly_bc() - use periodicity of "unit cell" to apply BCs to left-y edge
 * - ry_bc() - use periodicity of "unit cell" to apply BCs to right-y edge
 * - lz_bc() - use periodicity of "unit cell" to apply BCs to left-z edge
 * - rz_bc() - use periodicity of "unit cell" to apply BCs to right-z edge
 * - Ax() - x-component of vector potential for initial conditions
 * - Ay() - y-component of vector potential for initial conditions
 * - Az() - z-component of vector potential for initial conditions
 *
 * This is the equation for a set of discrete points which lie in a
 * plane and is a key relation for setting ghost cells in the boundary
 * condition routines below.  --  T. A. Gardiner  --  7/21/2006		      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef STATIC_MESH_REFINEMENT
#error The shkset3d test does not work with SMR.
#endif

/* Parameters which define initial solution -- made global so that they can be
 * shared with all private functions in this file */

static int rx, ry, rz;
static Real ang_2, ang_3; /* Rotation angles about the 2 and 3' axis */
static Real sin_a2, cos_a2, sin_a3, cos_a3;
/* The Riemann Problem Left and Right states */
static  Real dl, vxl, vyl, vzl;
static  Real dr, vxr, vyr, vzr;
#ifdef MHD
static  Real Bxl, Byl, Bzl, Bxr, Byr, Bzr;
#endif /* MHD */
#ifdef ADIABATIC
static  Real Pl, Pr;
#endif /* ADIABATIC */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * lx_bc() - use periodicity of "unit cell" to apply BCs to left-x edge
 * rx_bc() - use periodicity of "unit cell" to apply BCs to right-x edge
 * ly_bc() - use periodicity of "unit cell" to apply BCs to left-y edge
 * ry_bc() - use periodicity of "unit cell" to apply BCs to right-y edge
 * lz_bc() - use periodicity of "unit cell" to apply BCs to left-z edge
 * rz_bc() - use periodicity of "unit cell" to apply BCs to right-z edge
 * Ax() - x-component of vector potential for initial conditions
 * Ay() - y-component of vector potential for initial conditions
 * Az() - z-component of vector potential for initial conditions
 *============================================================================*/

static void lx_bc(GridS *pG);
static void rx_bc(GridS *pG);
static void ly_bc(GridS *pG);
static void ry_bc(GridS *pG);
static void lz_bc(GridS *pG);
static void rz_bc(GridS *pG);
static Real Ax(const Real x1, const Real x2, const Real x3);
static Real Ay(const Real x1, const Real x2, const Real x3);
static Real Az(const Real x1, const Real x2, const Real x3);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int n;
  int i, j, k, ix, jx, kx, mix, mjx, mkx;

  int Nx1, Nx2, Nx3; /* Complete Grid dimensions */
  int sx, sy, sz; /* Number of unit cells in each direction */
  div_t id;
  double d_ix;
  Real xl, x, xr, dt_usr;
  Real lx1,lx2,lx3,cx1,cx2,cx3,rx1,rx2,rx3;
  Real sp0_x1, sp0_x2, sp0_x3;

/* Variables assocaiated with a cubic array on which the Riemann problem
 * interface normal is aligned along the (1,1,1) direction. */
  int N, isqa, jsqa, ksqa;
  int scx, scy, scz; /* scx = N/rx, scy = N/ry, scz = N/rz */
  Real sdx, sdy, sdz;
  ConsS *pq, ***sqa=NULL; /* allocated to dimensions sqa[N][N][2N] */
#ifdef MHD
  Real ***sB1i=NULL, ***sB2i=NULL, ***sB3i=NULL;
#endif /* MHD */

/* The Riemann Problem Left and Right states in the Grid coordinate system */
  ConsS ql, qr;
/* The unit cell enclosing the Riemann problem interface. */
  ConsS ***qa=NULL; /* allocated to dimensions qa[rz][ry][rx] */
#ifdef MHD
  Real ***aB1i=NULL, ***aB2i=NULL, ***aB3i=NULL;
#endif /* MHD */
  int qa_min_ix, qa_min_jx, qa_min_kx;
  int qa_max_ix, qa_max_jx, qa_max_kx;

/*--- Step 1 -------------------------------------------------------------------
 * Read initial conditions, check that inputs are within valid range, allocate
 * memory for arrays.  We know this Domain must be the root Grid, since this
 * problem generator does not work with SMR.
 */

  Nx1 = pDomain->Nx[0];
  Nx2 = pDomain->Nx[1];
  Nx3 = pDomain->Nx[2];
  if(pDomain->Nx[2] <= 0)
    ath_error("[get_set_input]: This problem assumes a 3D domain\n");

/* Get dimensions of unit cell, check that total grid size is an integer number
 * of cells.  Reminder: the expression x % y in C returns the remainder when
 * x is divided by y */

  rx = par_geti("problem","rx");
  ry = par_geti("problem","ry");
  rz = par_geti("problem","rz");

  if(rx <= 0 || ry <= 0 || rz <= 0)
    ath_error("[get_set_input]: Only rx > 0, ry > 0 and rz > 0 is allowed\n");

  id = div(Nx1,rx);
  sx = id.quot;
  if(id.rem != 0)
    ath_error("[get_set_input]: Choose Nx1 % rx = 0 instead of %d\n",Nx1%rx);

  id = div(Nx2,ry);
  sy = id.quot;
  if(id.rem != 0)
    ath_error("[get_set_input]: Choose Nx2 % ry = 0 instead of %d\n",Nx2%ry);

  id = div(Nx3,rz);
  sz = id.quot;
  if(id.rem != 0)
    ath_error("[get_set_input]: Choose Nx3 % rz = 0 instead of %d\n",Nx3%rz);

  if(sy > sx || sz > sx)
    ath_error("Nx1/rx is assumed to be larger than the y- and z-dir. ratios\n");

/* Imposing translation symmetry over the domain (rx,ry,rz) gives */

  ang_3 = atan((rx*pGrid->dx1)/(ry*pGrid->dx2));
  sin_a3 = sin(ang_3);
  cos_a3 = cos(ang_3);

  ang_2 = atan((rx*pGrid->dx1*cos_a3)/(rz*pGrid->dx3));
  sin_a2 = sin(ang_2);
  cos_a2 = cos(ang_2);

/* Read in the left and right states in the Riemann problem */

  dl  = par_getd("problem","dl");
  vxl = par_getd("problem","vxl");
  vyl = par_getd("problem","vyl");
  vzl = par_getd("problem","vzl");
#ifdef MHD
  Bxl = par_getd("problem","Bxl");
  Byl = par_getd("problem","Byl");
  Bzl = par_getd("problem","Bzl");
#endif /* MHD */
#ifdef ADIABATIC
  Pl  = par_getd("problem","pl");
#endif /* ADIABATIC */

  dr  = par_getd("problem","dr");
  vxr = par_getd("problem","vxr");
  vyr = par_getd("problem","vyr");
  vzr = par_getd("problem","vzr");
#ifdef MHD
  Bxr = par_getd("problem","Bxr");
  Byr = par_getd("problem","Byr");
  Bzr = par_getd("problem","Bzr");
#endif /* MHD */
#ifdef ADIABATIC
  Pr  = par_getd("problem","pr");
#endif /* ADIABATIC */

/* Calculate the conservative (ConsS) L/R states in the Grid's coordinates */

/* This is the left state */
  ql.d  = dl;
  ql.M1 = dl*(vxl*cos_a2*cos_a3 - vyl*sin_a3 - vzl*sin_a2*cos_a3);
  ql.M2 = dl*(vxl*cos_a2*sin_a3 + vyl*cos_a3 - vzl*sin_a2*sin_a3);
  ql.M3 = dl*(vxl*sin_a2                     + vzl*cos_a2);

#ifdef MHD
  ql.B1c = (Bxl*cos_a2*cos_a3 - Byl*sin_a3 - Bzl*sin_a2*cos_a3);
  ql.B2c = (Bxl*cos_a2*sin_a3 + Byl*cos_a3 - Bzl*sin_a2*sin_a3);
  ql.B3c = (Bxl*sin_a2                     + Bzl*cos_a2);
#endif /* MHD */

#ifdef ADIABATIC
  ql.E  = Pl/Gamma_1
#ifdef MHD
    + 0.5*(Bxl*Bxl + Byl*Byl + Bzl*Bzl)
#endif /* MHD */
    + 0.5*dl*(vxl*vxl + vyl*vyl + vzl*vzl);
#endif

/* This is the right state */
  qr.d  = dr;
  qr.M1 = dr*(vxr*cos_a2*cos_a3 - vyr*sin_a3 - vzr*sin_a2*cos_a3);
  qr.M2 = dr*(vxr*cos_a2*sin_a3 + vyr*cos_a3 - vzr*sin_a2*sin_a3);
  qr.M3 = dr*(vxr*sin_a2                     + vzr*cos_a2);

#ifdef MHD
  qr.B1c = (Bxr*cos_a2*cos_a3 - Byr*sin_a3 - Bzr*sin_a2*cos_a3);
  qr.B2c = (Bxr*cos_a2*sin_a3 + Byr*cos_a3 - Bzr*sin_a2*sin_a3);
  qr.B3c = (Bxr*sin_a2                     + Bzr*cos_a2);
#endif /* MHD */

#ifdef ADIABATIC
  qr.E  = Pr/Gamma_1
#ifdef MHD
    + 0.5*(Bxr*Bxr + Byr*Byr + Bzr*Bzr)
#endif /* MHD */
    + 0.5*dr*(vxr*vxr + vyr*vyr + vzr*vzr);
#endif

/* Allocate 3D arrays with size of UNIT CELL */

  if((qa = (ConsS***)calloc_3d_array(rz, ry, rx+rx, sizeof(ConsS)))==NULL)
    ath_error("[shkset3d]: Error allocating qa[%d][%d][%d]\n",rz,ry,rx+rx);

#ifdef MHD
  if((aB1i = (Real***)calloc_3d_array(rz, ry, rx+rx+1, sizeof(Real))) == NULL)
    ath_error("[shkset3d]: Error allocating aB1i[%d][%d][%d]\n",rz,ry,rx+rx+1);

  if((aB2i = (Real***)calloc_3d_array(rz, ry+1, rx+rx, sizeof(Real))) == NULL)
    ath_error("[shkset3d]: Error allocating aB2i[%d][%d][%d]\n",rz,ry+1,rx+rx);

  if((aB3i = (Real***)calloc_3d_array(rz+1, ry, rx+rx, sizeof(Real))) == NULL)
    ath_error("[shkset3d]: Error allocating aB3i[%d][%d][%d]\n",rz+1,ry,rx+rx);
#endif /* MHD */

/* Initialize qa_min_co and qa_max_co */
  qa_min_jx = 0;  qa_max_jx = ry;
  qa_min_kx = 0;  qa_max_kx = rz;

/* qa_min_ix comes from calculating x=0 at jx=ry, kx=rz */
/* qa_max_ix comes from calculating x=0 at jx=0, kx=0 */
/* qa_max_ix - qa_min_ix = 2*rx! */

  d_ix = -pGrid->MinX[0]/pGrid->dx1 - rx*(pGrid->MinX[1]/(ry*pGrid->dx2) 
    + pGrid->MinX[2]/(rz*pGrid->dx3));
  qa_max_ix = ceil(d_ix);
  qa_min_ix = qa_max_ix - rx - rx;

/* Note that this code has the undesireable "feature" that when I wrote it I
 * assumed that d_ix would be equal to an integer.  This is not always the 
 * case..., so if it fails we should quit and give the user advice
 * T. A. Gardiner  --  7/21/2006  */

  if(qa_max_ix - d_ix > 1.0e-12){
    ath_perr(-1,"[shkset3d]: qa_max_ix - d_ix = %g\n",
	    qa_max_ix - d_ix);
    ath_error("[shkset3d]: Try setting the x2min and x3min = 0.0\n");
  }

/*--- Step 2 -------------------------------------------------------------------
 * Calculate the size of a square domain which can be restriced in an integral
 * fashion to the relevant grid qa[][][].  Note, one can multiply N, scx, scy,
 * and scz by some common integer to increase the resolution to the
 * approximation of the volume and interface integral averages.  Tests up to
 * now have shown no difference in the initial conditions when multiplying by
 * 3 or 4.*/

  N = rx*ry*rz;
  scx = ry*rz;
  scy = rx*rz;
  scz = rx*ry;

  sdx = pGrid->dx1/scx;
  sdy = pGrid->dx2/scy;
  sdz = pGrid->dx3/scz;

  if((sqa = (ConsS***)calloc_3d_array(N, N, N+N, sizeof(ConsS))) == NULL)
    ath_error("[shkset3d]: Error allocating sqa[%d][%d][%d]\n",N,N,N+N);

#ifdef MHD
  if((sB1i = (Real***)calloc_3d_array(N, N, N+N+1, sizeof(Real))) == NULL)
    ath_error("[shkset3d]: Error allocating sB1i[%d][%d][%d]\n",N,N,N+N+1);

  if((sB2i = (Real***)calloc_3d_array(N, N+1, N+N, sizeof(Real))) == NULL)
    ath_error("[shkset3d]: Error allocating sB2i[%d][%d][%d]\n",N,N+1,N+N);

  if((sB3i = (Real***)calloc_3d_array(N+1, N, N+N, sizeof(Real))) == NULL)
    ath_error("[shkset3d]: Error allocating sB3i[%d][%d][%d]\n",N+1,N,N+N);
#endif /* MHD */

/* Calculate the Position of the left-most corner of the sqa and qa grids. */
  sp0_x1 = pGrid->MinX[0] + qa_min_ix*pGrid->dx1;
  sp0_x2 = pGrid->MinX[1] + qa_min_jx*pGrid->dx2;
  sp0_x3 = pGrid->MinX[2] + qa_min_kx*pGrid->dx3;

/*--- Step 3 -------------------------------------------------------------------
 * First, initialize the interface magnetic fields in sqa
 */

#ifdef MHD
  for(k=0; k<=N; k++){
    for(j=0; j<=N; j++){
      for(i=0; i<=(N+N); i++){
/* Calculate the right corner position in sqa[k][j][i] */
	rx1 = sp0_x1 + (i+1)*sdx;
	rx2 = sp0_x2 + (j+1)*sdy;
	rx3 = sp0_x3 + (k+1)*sdz;

/* Calculate the left corner position in sqa[k][j][i] */
	lx1 = sp0_x1 + i*sdx;
	lx2 = sp0_x2 + j*sdy;
	lx3 = sp0_x3 + k*sdz;

/* Calculate the cell-center position and the x-position */
	cx1 = lx1 + 0.5*sdx;
	cx2 = lx2 + 0.5*sdy;
	cx3 = lx3 + 0.5*sdz;

/* Initialize B1i */
	if(j<N && k<N){ 
/* Calculate the left- and right-most x-position on this x1-face */
	  xl = lx1*cos_a2*cos_a3 + lx2*cos_a2*sin_a3 + lx3*sin_a2;
	  xr = lx1*cos_a2*cos_a3 + rx2*cos_a2*sin_a3 + rx3*sin_a2;

	  if(xl >= 0.0){ /* This cell is in the right state */
	    sB1i[k][j][i] = qr.B1c;
	  }
	  else if(xr <= 0.0){ /* This cell is in the left state */
	    sB1i[k][j][i] = ql.B1c;
	  }
	  else{
	    sB1i[k][j][i] = Bxl*cos_a2*cos_a3 +
	      (Az(lx1, rx2, cx3) - Az(lx1, lx2, cx3) )/sdy -
	      (Ay(lx1, cx2, rx3) - Ay(lx1, cx2, lx3) )/sdz;
	  }
	}

/* Initialize B2i */
	if(i<(N+N) && k<N){ 
/* Calculate the left- and right-most x-position on this x2-face */
	  xl = lx1*cos_a2*cos_a3 + lx2*cos_a2*sin_a3 + lx3*sin_a2;
	  xr = rx1*cos_a2*cos_a3 + lx2*cos_a2*sin_a3 + rx3*sin_a2;

	  if(xl >= 0.0){ /* This cell is in the right state */
	    sB2i[k][j][i] = qr.B2c;
	  }
	  else if(xr <= 0.0){ /* This cell is in the left state */
	    sB2i[k][j][i] = ql.B2c;
	  }
	  else{
	    sB2i[k][j][i] = Bxl*cos_a2*sin_a3 +
	      (Ax(cx1, lx2, rx3) - Ax(cx1, lx2, lx3) )/sdz -
	      (Az(rx1, lx2, cx3) - Az(lx1, lx2, cx3) )/sdx;
	  }
	}

/* Initialize B3i */
	if(i<(N+N) && j<N){ 
/* Calculate the left- and right-most x-position on this x3-face */
	  xl = lx1*cos_a2*cos_a3 + lx2*cos_a2*sin_a3 + lx3*sin_a2;
	  xr = rx1*cos_a2*cos_a3 + rx2*cos_a2*sin_a3 + lx3*sin_a2;

	  if(xl >= 0.0){ /* This cell is in the right state */
	    sB3i[k][j][i] = qr.B3c;
	  }
	  else if(xr <= 0.0){ /* This cell is in the left state */
	    sB3i[k][j][i] = ql.B3c;
	  }
	  else{
	    sB3i[k][j][i] = Bxl*sin_a2 +
	      (Ay(rx1, cx2, lx3) - Ay(lx1, cx2, lx3) )/sdx -
	      (Ax(cx1, rx2, lx3) - Ax(cx1, lx2, lx3) )/sdy;
	  }
	}
      }
    }
  }
#endif /* MHD */

/*--- Step 4 -------------------------------------------------------------------
 * Next, initialize the cell-centered quantities (density, momenta, total
 * energy, and cell-centered B) in sqa
 */

  for(k=0; k<N; k++){
    for(j=0; j<N; j++){
      for(i=0; i<(N+N); i++){
/* Calculate the right corner position in sqa[k][j][i] */
	rx1 = sp0_x1 + (i+1)*sdx;
	rx2 = sp0_x2 + (j+1)*sdy;
	rx3 = sp0_x3 + (k+1)*sdz;

	xr = rx1*cos_a2*cos_a3 + rx2*cos_a2*sin_a3 + rx3*sin_a2;

/* Calculate the left corner position in sqa[k][j][i] */
	lx1 = sp0_x1 + i*sdx;
	lx2 = sp0_x2 + j*sdy;
	lx3 = sp0_x3 + k*sdz;

	xl = lx1*cos_a2*cos_a3 + lx2*cos_a2*sin_a3 + lx3*sin_a2;

/* Calculate the cell-center position and the x-position */
	cx1 = lx1 + 0.5*sdx;
	cx2 = lx2 + 0.5*sdy;
	cx3 = lx3 + 0.5*sdz;
	x = cx1*cos_a2*cos_a3 + cx2*cos_a2*sin_a3 + cx3*sin_a2;

	if(xr <= 0.0){ /* This cell is in the left state */
	  sqa[k][j][i] = ql;
	  continue;
	}

	if(xl >= 0.0){ /* This cell is in the right state */
	  sqa[k][j][i] = qr;
	  continue;
	}

	pq = &(sqa[k][j][i]); /* get a pointer to the fluid element */

/* Initialize the density, the momenta and the total energy */
	if(x < 0.0){ /* Left state */
	  pq->d  = dl;
	  pq->M1 = dl*(vxl*cos_a2*cos_a3 - vyl*sin_a3 - vzl*sin_a2*cos_a3);
	  pq->M2 = dl*(vxl*cos_a2*sin_a3 + vyl*cos_a3 - vzl*sin_a2*sin_a3);
	  pq->M3 = dl*(vxl*sin_a2                     + vzl*cos_a2);
#ifdef ADIABATIC
	  pq->E  = Pl/Gamma_1 + 0.5*dl*(vxl*vxl + vyl*vyl + vzl*vzl);
#endif
	}
	else{ /* Right state */
	  pq->d  = dr;
	  pq->M1 = dr*(vxr*cos_a2*cos_a3 - vyr*sin_a3 - vzr*sin_a2*cos_a3);
	  pq->M2 = dr*(vxr*cos_a2*sin_a3 + vyr*cos_a3 - vzr*sin_a2*sin_a3);
	  pq->M3 = dr*(vxr*sin_a2                     + vzr*cos_a2);
#ifdef ADIABATIC
	  pq->E  = Pr/Gamma_1 + 0.5*dr*(vxr*vxr + vyr*vyr + vzr*vzr);
#endif
	}

/* Initialize the cell-center magnetic field components */
#ifdef MHD
        pq->B1c = 0.5*(sB1i[k][j][i] + sB1i[k  ][j  ][i+1]);
	pq->B2c = 0.5*(sB2i[k][j][i] + sB2i[k  ][j+1][i  ]);
	pq->B3c = 0.5*(sB3i[k][j][i] + sB3i[k+1][j  ][i  ]);

#ifdef ADIABATIC
	pq->E += 0.5*(pq->B1c*pq->B1c + pq->B2c*pq->B2c + pq->B3c*pq->B3c);
#endif /* ADIABATIC */
#endif /* MHD */
      }
    }
  }

/*--- Step 5 -------------------------------------------------------------------
 *  Now conservatively restrict sqa[][][] into qa[][][] 
 */

/* First, the interface magnetic fields */

#ifdef MHD
  for(k=0; k<=rz; k++){
    for(j=0; j<=ry; j++){
      for(i=0; i<=(rx+rx); i++){

/* Initialize B1i */
	if(j<ry && k<rz){ 
	  aB1i[k][j][i] = 0.0;
	  isqa = i*scx; /* The left-x1 interface field */
	  for(ksqa=k*scz; ksqa<(k+1)*scz; ksqa++){
	    for(jsqa=j*scy; jsqa<(j+1)*scy; jsqa++){
	      aB1i[k][j][i] += sB1i[ksqa][jsqa][isqa];
	    }
	  }
	  aB1i[k][j][i] /= (Real)(scy*scz);
	}

/* Initialize B2i */
	if(i<(rx+rx) && k<rz){ 
	  aB2i[k][j][i] = 0.0;
	  jsqa = j*scy; /* The left-x2 interface field */
	  for(ksqa=k*scz; ksqa<(k+1)*scz; ksqa++){
	    for(isqa=i*scx; isqa<(i+1)*scx; isqa++){
	      aB2i[k][j][i] += sB2i[ksqa][jsqa][isqa];
	    }
	  }
	  aB2i[k][j][i] /= (Real)(scx*scz);
	}

/* Initialize B3i */
	if(i<(rx+rx) && j<ry){ 
	  aB3i[k][j][i] = 0.0;
	  ksqa = k*scz; /* The left-x3 interface field */
	  for(jsqa=j*scy; jsqa<(j+1)*scy; jsqa++){
	    for(isqa=i*scx; isqa<(i+1)*scx; isqa++){
	      aB3i[k][j][i] += sB3i[ksqa][jsqa][isqa];
	    }
	  }
	  aB3i[k][j][i] /= (Real)(scx*scy);
	}
      }
    }
  }
#endif /* MHD */

/* Next, the volume averaged quantities */
  for(k=0; k<rz; k++){
    for(j=0; j<ry; j++){
      for(i=0; i<(rx+rx); i++){
	qa[k][j][i].d = 0.0;
	qa[k][j][i].M1 = qa[k][j][i].M2 = qa[k][j][i].M3 = 0.0;
#ifdef ADIABATIC
	qa[k][j][i].E = 0.0;
#endif

/* loop over the overlapping elements in the sqa array and sum up the
 * conservative variables. */

	for(ksqa=k*scz; ksqa<(k+1)*scz; ksqa++){
	  for(jsqa=j*scy; jsqa<(j+1)*scy; jsqa++){
	    for(isqa=i*scx; isqa<(i+1)*scx; isqa++){
	      qa[k][j][i].d  += sqa[ksqa][jsqa][isqa].d;
	      qa[k][j][i].M1 += sqa[ksqa][jsqa][isqa].M1;
	      qa[k][j][i].M2 += sqa[ksqa][jsqa][isqa].M2;
	      qa[k][j][i].M3 += sqa[ksqa][jsqa][isqa].M3;
#ifdef ADIABATIC
 	      qa[k][j][i].E  += sqa[ksqa][jsqa][isqa].E;
#endif
	    }
	  }
	}

/* Normalize them - Note, scx*scy*scz = N*N; */

	qa[k][j][i].d  /= (Real)(N*N);
	qa[k][j][i].M1 /= (Real)(N*N);
	qa[k][j][i].M2 /= (Real)(N*N);
	qa[k][j][i].M3 /= (Real)(N*N);
#ifdef ADIABATIC
	qa[k][j][i].E  /= (Real)(N*N);
#endif

/* Finally initialize the cell center field components. */
#ifdef MHD
	qa[k][j][i].B1c = 0.5*(aB1i[k][j][i] + aB1i[k  ][j  ][i+1]);
	qa[k][j][i].B2c = 0.5*(aB2i[k][j][i] + aB2i[k  ][j+1][i  ]);
	qa[k][j][i].B3c = 0.5*(aB3i[k][j][i] + aB3i[k+1][j  ][i  ]);
#endif /* MHD */
      }
    }
  }

/*--- Step 6 -------------------------------------------------------------------
 * Now initialize variables over the whole grid
 */

  for (k=0; k<=(pGrid->ke)+nghost; k++) {
    for (j=0; j<=(pGrid->je)+nghost; j++) {
      for (i=0; i<=(pGrid->ie)+nghost; i++) {

/* Calculate the coordinate of this fluid element and the mapped coordinate to
 * the equivalent fluid element on the interval (0 < iy <= ry) and
 * (0 < iz <= rz). */

	mix = ix = (i-(pGrid->is)) + pGrid->Disp[0];
	mjx = jx = (j-(pGrid->js)) + pGrid->Disp[1];
	mkx = kx = (k-(pGrid->ks)) + pGrid->Disp[2];

	n=0;
	while(mjx < 0){
	  mjx += ry;
	  n++;
	}
	if(n){
	  mix -= n*rx;
	}

	n=0;
	while(mjx >= ry){
	  mjx -= ry;
	  n++;
	}
	if(n){
	  mix += n*rx;
	}

	n=0;
	while(mkx < 0){
	  mkx += rz;
	  n++;
	}
	if(n){
	  mix -= n*rx;
	}

	n=0;
	while(mkx >= rz){
	  mkx -= rz;
	  n++;
	}
	if(n){
	  mix += n*rx;
	}

	if(mix < qa_min_ix){ /* This cell is in the left state */
	  pGrid->U[k][j][i] = ql;
#ifdef MHD
	  pGrid->B1i[k][j][i] = ql.B1c;
	  pGrid->B2i[k][j][i] = ql.B2c;
	  pGrid->B3i[k][j][i] = ql.B3c;
#endif
	}
	else if(mix >= qa_max_ix){ /* This cell is in the right state */
	  pGrid->U[k][j][i] = qr;
#ifdef MHD
	  pGrid->B1i[k][j][i] = qr.B1c;
	  pGrid->B2i[k][j][i] = qr.B2c;
	  pGrid->B3i[k][j][i] = qr.B3c;
#endif
	}
	else{ /* Otherwise it's in the qa[][][] array */
	  pGrid->U[k][j][i] = qa[mkx][mjx][mix - qa_min_ix];
#ifdef MHD
	  pGrid->B1i[k][j][i] = aB1i[mkx][mjx][mix - qa_min_ix];
	  pGrid->B2i[k][j][i] = aB2i[mkx][mjx][mix - qa_min_ix];
	  pGrid->B3i[k][j][i] = aB3i[mkx][mjx][mix - qa_min_ix];
#endif
	}
      }
    }
  }

/*--- Step 7 -------------------------------------------------------------------
 * set function pointers for BCs, and conclude
 */

  bvals_mhd_fun(pDomain, left_x1,  lx_bc);
  bvals_mhd_fun(pDomain, right_x1, rx_bc);
  bvals_mhd_fun(pDomain, left_x2,  ly_bc);
  bvals_mhd_fun(pDomain, right_x2, ry_bc);
  bvals_mhd_fun(pDomain, left_x3,  lz_bc);
  bvals_mhd_fun(pDomain, right_x3, rz_bc);

  free_3d_array((void***)sqa );  sqa  = NULL;
#ifdef MHD
  free_3d_array((void***)sB1i);  sB1i = NULL;
  free_3d_array((void***)sB2i);  sB2i = NULL;
  free_3d_array((void***)sB3i);  sB3i = NULL;
#endif /* MHD */

  return;
}


/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM,FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/*! \fn static void lx_bc(GridS *pG)
 *  \brief Apply boundary condition in left-x direction
 */

static void lx_bc(GridS *pG)
{
  int i, j, k, mi, mj, mk, is, js, je, ks, ke;
  is = pG->is;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;

  for(i=is-1; i>=0; i--){ /* Do NOT Change this loop ordering! */
    for(k = 0; k <= ke+nghost; k++){
      for(j = 0; j <= je+nghost; j++){
/* Test the nearest "2D" translation vectors */
	if(k - rz >= ks){
	  mi = i + rx;
	  mj = j;
	  mk = k - rz;
	}
	else if(j - ry >= js){
	  mi = i + rx;
	  mj = j - ry;
	  mk = k;
	}
        else continue;

	pG->U[k][j][i] = pG->U[mk][mj][mi];
#ifdef MHD
	pG->B1i[k][j][i] = pG->B1i[mk][mj][mi];
        pG->B2i[k][j][i] = pG->B2i[mk][mj][mi];
        pG->B3i[k][j][i] = pG->B3i[mk][mj][mi];
#endif
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void rx_bc(GridS *pG)
 *  \brief Apply boundary condition in right-x direction
 */

static void rx_bc(GridS *pG)
{
  int i, j, k, mi, mj, mk, ie, js, je, ks, ke;
  ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;

  for(i=ie+1; i<=ie+nghost; i++){ /* Do NOT Change this loop ordering! */
    for(k = 0; k <= ke+nghost; k++){
      for(j = 0; j <= je+nghost; j++){
/* Test the nearest "2D" translation vectors */
	if(k + rz <= ke){
	  mi = i - rx;
	  mj = j;
	  mk = k + rz;
	}
	else if(j + ry <= je){
	  mi = i;
	  mj = j;
	  mk = k;
	}
        else continue;

	pG->U[k][j][i] = pG->U[mk][mj][mi];
#ifdef MHD
	if(i > ie+1) pG->B1i[k][j][i] = pG->B1i[mk][mj][mi];
        pG->B2i[k][j][i] = pG->B2i[mk][mj][mi];
        pG->B3i[k][j][i] = pG->B3i[mk][mj][mi];
#endif
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void ly_bc(GridS *pG) 
 *  \brief Apply boundary condition in left-y direction
 */

static void ly_bc(GridS *pG)
{
  int i, j, k, mi, mj, mk, is, ie, js, ks, ke;
  is = pG->is; ie = pG->ie;
  js = pG->js;
  ks = pG->ks; ke = pG->ke;

  for(j=js-1; j>=0; j--){
    for(k = 0; k <= ke+nghost; k++){
      for(i = 0; i <= ie+nghost; i++){
/* Test the nearest "2D" translation vectors */
        if(i - rx >= is){
          mi = i - rx;
          mj = j + ry;
          mk = k;
        }
        else if(k - rz >= ks){
          mi = i;
          mj = j + ry;
          mk = k - rz;
        }
        else continue;

	pG->U[k][j][i] = pG->U[mk][mj][mi];
#ifdef MHD
	pG->B1i[k][j][i] = pG->B1i[mk][mj][mi];
        pG->B2i[k][j][i] = pG->B2i[mk][mj][mi];
        pG->B3i[k][j][i] = pG->B3i[mk][mj][mi];
#endif
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void ry_bc(GridS *pG)
 *  \brief  Apply boundary condition in right-y direction
 */

static void ry_bc(GridS *pG)
{
  int i, j, k, mi, mj, mk, is, ie, je, ks, ke;
  is = pG->is; ie = pG->ie;
  je = pG->je;
  ks = pG->ks; ke = pG->ke;

  for(j=je+1; j<=je+nghost; j++){
    for(k = 0; k <= ke+nghost; k++){
      for(i = 0; i <= ie+nghost; i++){
/* Test the nearest "2D" translation vectors */
        if(i + rx <= ie){
          mi = i + rx;
          mj = j - ry;
          mk = k;
        }
        else if(k + rz <= ke){
          mi = i;
          mj = j - ry;
          mk = k + rz;
        }
        else continue;

	pG->U[k][j][i] = pG->U[mk][mj][mi];
#ifdef MHD
	pG->B1i[k][j][i] = pG->B1i[mk][mj][mi];
        if(j > je+1) pG->B2i[k][j][i] = pG->B2i[mk][mj][mi];
        pG->B3i[k][j][i] = pG->B3i[mk][mj][mi];
#endif
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void lz_bc(GridS *pG)
 *  \brief Apply boundary condition in left-z direction
 */

static void lz_bc(GridS *pG)
{
  int i, j, k, mi, mj, mk, is, ie, js, je, ks;
  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks;

  for(k=ks-1; k>=0; k--){
    for(j = 0; j <= je+nghost; j++){
      for(i = 0; i <= ie+nghost; i++){
/* Test the nearest "2D" translation vectors */
	if(i - rx >= is){
	  mi = i - rx;
	  mj = j;
	  mk = k + rz;
	}
	else if(j - ry >= js){
	  mi = i;
	  mj = j - ry;
	  mk = k + rz;
	}
        else continue;

	pG->U[k][j][i] = pG->U[mk][mj][mi];
#ifdef MHD
	pG->B1i[k][j][i] = pG->B1i[mk][mj][mi];
        pG->B2i[k][j][i] = pG->B2i[mk][mj][mi];
        pG->B3i[k][j][i] = pG->B3i[mk][mj][mi];
#endif
      }
    }
  }

  return;
}

/*---------------------------------------------------------------------------*/
/*! \fn static void rz_bc(GridS *pG)
 *  \brief Apply boundary condition in right-z direction
 */

static void rz_bc(GridS *pG)
{
  int i, j, k, mi, mj, mk, is, ie, js, je, ke;
  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ke = pG->ke;

  for(k=ke+1; k<=ke+nghost; k++){
    for(j = 0; j <= je+nghost; j++){
      for(i = 0; i <= ie+nghost; i++){
/* Test the nearest "2D" translation vectors */
	if(i + rx <= ie){
	  mi = i + rx;
	  mj = j;
	  mk = k - rz;
	}
	else if(j + ry <= je){
	  mi = i;
	  mj = j + ry;
	  mk = k - rz;
	}
        else continue;

	pG->U[k][j][i] = pG->U[mk][mj][mi];
#ifdef MHD
	pG->B1i[k][j][i] = pG->B1i[mk][mj][mi];
        pG->B2i[k][j][i] = pG->B2i[mk][mj][mi];
        if(k > ke+1) pG->B3i[k][j][i] = pG->B3i[mk][mj][mi];
#endif
      }
    }
  }

  return;
}

/*-----------------------------------------------------------------------------
 * Ax,Ay,Az: x,y,z-components of vector potential
 */

#ifdef MHD
/*! \fn static Real Ax(const Real x, const Real y, const Real z)
 *  \brief s-component of vector potential */
static Real Ax(const Real x, const Real y, const Real z)
{
  Real x1 = x*cos_a2*cos_a3 + y*cos_a2*sin_a3 + z*sin_a2;
  Real A2 =  x1*(x1 < 0.0 ? Bzl : Bzr);
  Real A3 = -x1*(x1 < 0.0 ? Byl : Byr);

  return -A2*sin_a3 - A3*sin_a2*cos_a3;
}

/*! \fn static Real Ay(const Real x, const Real y, const Real z)
 *  \brief y-component of vector potential */
static Real Ay(const Real x, const Real y, const Real z)
{
  Real x1 = x*cos_a2*cos_a3 + y*cos_a2*sin_a3 + z*sin_a2;
  Real A2 =  x1*(x1 < 0.0 ? Bzl : Bzr);
  Real A3 = -x1*(x1 < 0.0 ? Byl : Byr);

  return A2*cos_a3 - A3*sin_a2*sin_a3;
}

/*! \fn static Real Az(const Real x, const Real y, const Real z) 
 *  \brief z-component of vector potential */
static Real Az(const Real x, const Real y, const Real z)
{
  Real x1 = x*cos_a2*cos_a3 + y*cos_a2*sin_a3 + z*sin_a2;
/* Real A2 =  x1*(x1 < 0.0 ? Bzl : Bzr); */
  Real A3 = -x1*(x1 < 0.0 ? Byl : Byr);

  return A3*cos_a2;
}
#endif /* MHD */
