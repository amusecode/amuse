#include "../copyright.h"
/*============================================================================*/
/*! \file selfg_multigrid.c
 *  \brief Contains functions to solve Poisson's equation for self-gravity in
 *   3D using multigrid.
 *
 *   These functions work for non-periodic domains.  A low-order multipole
 *   expansion is used to compute the potential on the boundaries.
 *
 * HISTORY:
 * - june-2007 - 2D and 3D solvers written by Irene Balmes
 * - july-2007 - routines incorporated into Athena by JMS and IB
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - selfg_by_multig_2d() - 2D Poisson solver using multigrid
 * - selfg_by_multig_3d() - 3D Poisson solver using multigrid
 * - selfg_by_multig_3d_init() - Initializes send/receive buffers for MPI */
/*============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef SELF_GRAVITY_USING_MULTIGRID

#ifdef STATIC_MESH_REFINEMENT
#error self gravity with multigrid not yet implemented to work with SMR
#endif

#ifdef MPI_PARALLEL
#error self gravity with multigrid not yet implemented to work with MPI
#endif

#ifdef MPI_PARALLEL
static double *send_buf=NULL, *recv_buf=NULL;
#endif

/*! \struct MGrid
 *  \brief Holds RHS, potential, and information about grid
 * size for a given level in the multi-grid hierarchy  */
typedef struct MGrid_s{
  Real ***rhs,***Phi;  /* RHS of elliptic equation, and solution */
  Real dx1,dx2,dx3;
  int Nx1,Nx2,Nx3;
  int is,ie;
  int js,je;
  int ks,ke; 
  int rx1_id, lx1_id;
  int rx2_id, lx2_id;
  int rx3_id, lx3_id;
}MGrid;

/* 3D temporary array needed for restriction of errors  */
Real ***error;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   multig_3d() -
 *============================================================================*/

void multig_3d(MGrid *pMG);
void Jacobi(MGrid *pMG);
void Restriction_3d(MGrid *pMG_fine, MGrid *pMG_coarse);
void Prolongation_3d(MGrid *pMG_coarse, MGrid *pMG_fine);

#ifdef MPI_PARALLEL
void set_mg_bvals(MGrid *pMG);
void swap_mg_ix1(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq);
void swap_mg_ox1(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq);
void swap_mg_ix2(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq);
void swap_mg_ox2(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq);
void swap_mg_ix3(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq);
void swap_mg_ox3(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq);
#endif


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void selfg_multig_1d(DomainS *pD)
 *  \brief 1-D multigrid gravity
 */

void selfg_multig_1d(DomainS *pD)
{
  ath_error("[selfg_multig_1d] 1D multigrid not implemented yet\n");
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void selfg_multig_2d(DomainS *pD)
 *  \brief 2-D multigrid gravity
 */

void selfg_multig_2d(DomainS *pD)
{
  ath_error("[selfg_multig_2d] 2D multigrid not implemented yet\n");
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void selfg_multig_3d(DomainS *pD)
 *  \brief Do not use with periodic BCs, uses multipole expansion
 *   to compute potential at boundary
 */

void selfg_multig_3d(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int Nx1z, Nx2z, Nx3z;
  MGrid Root_grid;
  Real mass = 0.0, tmass, dVol, rad, x1, x2, x3;
  Real Grav_const = four_pi_G/(4.0*PI);
#ifdef MPI_PARALLEL
  Real mpi_err;
#endif

/* Copy current potential into old */

  for (k=ks-nghost; k<=ke+nghost; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->Phi_old[k][j][i] = pG->Phi[k][j][i];
      }
    }
  }

/* Compute solution at boundaries using monopole expansion */

  dVol = pG->dx1*pG->dx2*pG->dx3;
  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        mass += pG->U[k][j][i].d*dVol;
      }
    }
  }

#ifdef MPI_PARALLEL
  mpi_err = MPI_Allreduce(&mass, &tmass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (mpi_err) ath_error("[selfg_multigrid]: MPI_Reduce returned err = %d\n",
    mpi_err);
#else
  tmass = mass;
#endif /* MPI_PARALLEL */

/*  Inner and outer x1 boundaries */

  if (pG->lx1_id < 0) {
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=1; i<=nghost; i++){
          cc_pos(pG,is-i,j,k,&x1,&x2,&x3);
          rad = sqrt(x1*x1 + x2*x2 + x3*x3);
          pG->Phi[k][j][is-i] = -Grav_const*tmass/rad;
        }
      }
    }
  }

  if (pG->rx1_id < 0) {
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=1; i<=nghost; i++){
          cc_pos(pG,ie+i,j,k,&x1,&x2,&x3);
          rad = sqrt(x1*x1 + x2*x2 + x3*x3);
          pG->Phi[k][j][ie+i] = -Grav_const*tmass/rad;
        }
      }
    }
  }

/*  Inner and outer x2 boundaries */

  if (pG->lx2_id < 0) {
    for (k=ks; k<=ke; k++){
      for (j=1; j<=nghost; j++){
        for (i=is-nghost; i<=ie+nghost; i++){
          cc_pos(pG,i,js-j,k,&x1,&x2,&x3);
          rad = sqrt(x1*x1 + x2*x2 + x3*x3);
          pG->Phi[k][js-j][i] = -Grav_const*tmass/rad;
        }
      }
    }
  }

  if (pG->rx2_id < 0) {
    for (k=ks; k<=ke; k++){
      for (j=1; j<=nghost; j++){
        for (i=is-nghost; i<=ie+nghost; i++){
          cc_pos(pG,i,je+j,k,&x1,&x2,&x3);
          rad = sqrt(x1*x1 + x2*x2 + x3*x3);
          pG->Phi[k][je+j][i] = -Grav_const*tmass/rad;
        }
      }
    }
  }

/*  Inner and outer x3 boundaries */

  if (pG->lx3_id < 0) {
    for (k=1; k<=nghost; k++){
      for (j=js-nghost; j<=je+nghost; j++){
        for (i=is-nghost; i<=ie+nghost; i++){
          cc_pos(pG,i,j,ks-k,&x1,&x2,&x3);
          rad = sqrt(x1*x1 + x2*x2 + x3*x3);
          pG->Phi[ks-k][j][i] = -Grav_const*tmass/rad;
        }
      }
    }
  }

  if (pG->rx3_id < 0) {
    for (k=1; k<=nghost; k++){
      for (j=js-nghost; j<=je+nghost; j++){
        for (i=is-nghost; i<=ie+nghost; i++){
          cc_pos(pG,i,j,ke+k,&x1,&x2,&x3);
          rad = sqrt(x1*x1 + x2*x2 + x3*x3);
          pG->Phi[ke+k][j][i] = -Grav_const*tmass/rad;
        }
      }
    }
  }

/* Initialize MGrid structure for top level (the root, or finest, level) */

  Root_grid.Nx1 = pG->Nx[0];
  Root_grid.Nx2 = pG->Nx[1];
  Root_grid.Nx3 = pG->Nx[2];
  Root_grid.is = 1;  Root_grid.ie = pG->Nx[0];
  Root_grid.js = 1;  Root_grid.je = pG->Nx[1];
  Root_grid.ks = 1;  Root_grid.ke = pG->Nx[2];
  Root_grid.dx1 = pG->dx1;
  Root_grid.dx2 = pG->dx2;
  Root_grid.dx3 = pG->dx3;
  Root_grid.rx1_id = pG->rx1_id; Root_grid.lx1_id = pG->lx1_id;
  Root_grid.rx2_id = pG->rx2_id; Root_grid.lx2_id = pG->lx2_id;
  Root_grid.rx3_id = pG->rx3_id; Root_grid.lx3_id = pG->lx3_id;

/* There is only one ghost zone needed at each level, not nghost */
  Nx1z = pG->Nx[0] + 2;
  Nx2z = pG->Nx[1] + 2;
  Nx3z = pG->Nx[2] + 2;
  Root_grid.rhs = (Real ***) calloc_3d_array(Nx3z,Nx2z,Nx1z,sizeof(Real));
  Root_grid.Phi = (Real ***) calloc_3d_array(Nx3z,Nx2z,Nx1z,sizeof(Real));
  if (Root_grid.rhs == NULL) {
    ath_error("[selfg_by_multig_3d]: Error allocating memory\n");
  }
  if (Root_grid.Phi == NULL) {
    ath_error("[selfg_by_multig_3d]: Error allocating memory\n");
  }

/* Initialize solution on root grid, including single ghost zone */
  for (k=ks-1; k<=ke+1; k++){
    for (j=js-1; j<=je+1; j++){
      for (i=is-1; i<=ie+1; i++){
        Root_grid.rhs[k-ks+1][j-js+1][i-is+1] = four_pi_G*pG->U[k][j][i].d;
        Root_grid.Phi[k-ks+1][j-js+1][i-is+1] = pG->Phi[k][j][i];
      }
    }
  }
#ifdef MPI_PARALLEL
  set_mg_bvals(&Root_grid);
#endif

/* Compute new potential.  Note multig_3d calls itself recursively. */

  multig_3d(&Root_grid);

/* copy solution for potential from MGrid into Grid structure.  Boundary
 * conditions for nghost ghost cells are set by set_bvals() call in main() */

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        pG->Phi[k][j][i] = Root_grid.Phi[k-ks+1][j-js+1][i-is+1];
      }
    }
  }

  free_3d_array(Root_grid.rhs);
  free_3d_array(Root_grid.Phi);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void multig_3d(MGrid *pMG)
 *  \brief Functions needed for the multigrid solver in 3D
 */

void multig_3d(MGrid *pMG)
{
  int i, is = pMG->is, ie = pMG->ie;
  int j, js = pMG->js, je = pMG->je;
  int k, ks = pMG->ks, ke = pMG->ke;
  int Nx1z, Nx2z, Nx3z;
  MGrid Coarse_grid;

/* If we are down to 4 cells in any dimension do 10 iterations and return */

  if (pMG->Nx1==4 || pMG->Nx2==4 || pMG->Nx3==4)
    Jacobi(pMG);

/* Else, do 10 iterations at this level, restrict to a coarser grid, and call
 * multig_3d again with this coarse grid */

  else { 
    Jacobi(pMG);

/* Allocate and initialize MGrid at next coarsest level */

    Coarse_grid.Nx1 = pMG->Nx1/2;
    Coarse_grid.Nx2 = pMG->Nx2/2;
    Coarse_grid.Nx3 = pMG->Nx3/2;
    Coarse_grid.is = 1;  Coarse_grid.ie = Coarse_grid.Nx1;
    Coarse_grid.js = 1;  Coarse_grid.je = Coarse_grid.Nx2;
    Coarse_grid.ks = 1;  Coarse_grid.ke = Coarse_grid.Nx3;
    Coarse_grid.dx1 = 2.0*pMG->dx1;
    Coarse_grid.dx2 = 2.0*pMG->dx2;
    Coarse_grid.dx3 = 2.0*pMG->dx3;
    Coarse_grid.rx1_id = pMG->rx1_id; Coarse_grid.lx1_id = pMG->lx1_id;
    Coarse_grid.rx2_id = pMG->rx2_id; Coarse_grid.lx2_id = pMG->lx2_id;
    Coarse_grid.rx3_id = pMG->rx3_id; Coarse_grid.lx3_id = pMG->lx3_id;

/* Again, only one ghost zone on each level */
    Nx1z = (pMG->Nx1)/2 + 2;
    Nx2z = (pMG->Nx2)/2 + 2;
    Nx3z = (pMG->Nx3)/2 + 2;
    Coarse_grid.rhs = (Real ***) calloc_3d_array(Nx3z,Nx2z,Nx1z,sizeof(Real));
    Coarse_grid.Phi = (Real ***) calloc_3d_array(Nx3z,Nx2z,Nx1z,sizeof(Real));
    if (Coarse_grid.rhs == NULL) {
      ath_error("[multig_3d]: Error allocating memory for some level\n");
    }
    if (Coarse_grid.Phi == NULL) {
      ath_error("[multig_3d]: Error allocating memory for some level\n");
    }

    Restriction_3d(pMG, &Coarse_grid);
#ifdef MPI_PARALLEL
    set_mg_bvals(&Coarse_grid);
#endif

    multig_3d(&Coarse_grid);

/* The following code is first reached after 10 iterations at the coarsest
 * level.  We then prolongate, do 10 iterations, and return.  This will return
 * execution to this same spot for the next coarsest level, so we will
 * prolongate, do 10 iterations, return, and so on.
 */

    Prolongation_3d(&Coarse_grid, pMG);
    free_3d_array(Coarse_grid.rhs);
    free_3d_array(Coarse_grid.Phi);
#ifdef MPI_PARALLEL
    set_mg_bvals(pMG);
#endif

/* End with 10 iterations at the finest level */

    Jacobi(pMG);
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void Jacobi(MGrid *pMG)
 *  \brief Jacobi iterations. 
 *
 *   Do not use with periodic BCs, uses multipole expansion
 *   to compute potential at boundary
 */
void Jacobi(MGrid *pMG)
{
  int i, is = pMG->is, ie = pMG->ie;
  int j, js = pMG->js, je = pMG->je;
  int k, ks = pMG->ks, ke = pMG->ke;
  int n;
  Real dx1sq = (pMG->dx1*pMG->dx1);
  Real dx2sq = (pMG->dx2*pMG->dx2);
  Real dx3sq = (pMG->dx3*pMG->dx3);

/* Jacobi iterations in 3D */

  if (pMG->Nx3 > 1) {
    for (n=0; n<=10; n++){  /* hardwired to do 10 iterations */
      for (k=ks; k<=ke; k++){
        for (j=js; j<=je; j++){
          for (i=is; i<=ie; i++){
            error[k][j][i] = -pMG->rhs[k][j][i];
            error[k][j][i] += (pMG->Phi[k][j][i+1] + pMG->Phi[k][j][i-1])/dx1sq;
            error[k][j][i] += (pMG->Phi[k][j+1][i] + pMG->Phi[k][j-1][i])/dx2sq;
            error[k][j][i] += (pMG->Phi[k+1][j][i] + pMG->Phi[k-1][j][i])/dx3sq;
          }
        }
      }

      for (k=ks; k<=ke; k++){
        for (j=js; j<=je; j++){
          for (i=is; i<=ie; i++){
            pMG->Phi[k][j][i] = error[k][j][i]/(2.0/dx1sq + 2.0/dx2sq + 2.0/dx3sq);
          }
        }
      }
#ifdef MPI_PARALLEL
      set_mg_bvals(pMG);
#endif
    }
  }

/* Jacobi iterations in 2D (x-y plane) */

  if (pMG->Nx3 == 1) {
    for (n=0; n<=10; n++){  /* hardwired to do 10 iterations */
      for (j=js; j<=je; j++){
        for (i=is; i<=ie; i++){
          error[ks][j][i] = -pMG->rhs[ks][j][i];
          error[ks][j][i] += (pMG->Phi[ks][j][i+1] + pMG->Phi[ks][j][i-1])/dx1sq;
          error[ks][j][i] += (pMG->Phi[ks][j+1][i] + pMG->Phi[ks][j-1][i])/dx2sq;
        }
      }

      for (j=js; j<=je; j++){
        for (i=is; i<=ie; i++){
          pMG->Phi[ks][j][i] = error[ks][j][i]/(2.0/dx1sq + 2.0/dx2sq);
        }
      }
#ifdef MPI_PARALLEL
      set_mg_bvals(pMG);
#endif
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void Restriction_3d(MGrid *pMG_fine, MGrid *pMG_coarse) 
 *  \brief Averages fine grid solution onto coarse
 */

void Restriction_3d(MGrid *pMG_fine, MGrid *pMG_coarse)
{
  int i, is = pMG_fine->is, ie = pMG_fine->ie;
  int j, js = pMG_fine->js, je = pMG_fine->je;
  int k, ks = pMG_fine->ks, ke = pMG_fine->ke;
  Real dx1sq = (pMG_fine->dx1*pMG_fine->dx1);
  Real dx2sq = (pMG_fine->dx2*pMG_fine->dx2);
  Real dx3sq = (pMG_fine->dx3*pMG_fine->dx3);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        error[k][j][i] = pMG_fine->rhs[k][j][i];
        error[k][j][i] -= (pMG_fine->Phi[k][j][i+1] + pMG_fine->Phi[k][j][i-1]
          - 2.0*pMG_fine->Phi[k][j][i]) / dx1sq;
        error[k][j][i] -= (pMG_fine->Phi[k][j+1][i] + pMG_fine->Phi[k][j-1][i]
          - 2.0*pMG_fine->Phi[k][j][i]) / dx2sq;
        error[k][j][i] -= (pMG_fine->Phi[k+1][j][i] + pMG_fine->Phi[k-1][j][i]
          - 2.0*pMG_fine->Phi[k][j][i]) / dx3sq;
      }
    }
  }

  for(k=ks; k<=pMG_coarse->ke; k++){
    for (j=js; j<=pMG_coarse->je; j++){
      for (i=is; i<=pMG_coarse->ie; i++){
        pMG_coarse->Phi[k][j][i] =
          (pMG_fine->Phi[2*k  ][2*j  ][2*i] + pMG_fine->Phi[2*k  ][2*j  ][2*i-1]
         + pMG_fine->Phi[2*k  ][2*j-1][2*i] + pMG_fine->Phi[2*k  ][2*j-1][2*i-1]
         + pMG_fine->Phi[2*k-1][2*j  ][2*i] + pMG_fine->Phi[2*k-1][2*j  ][2*i-1]
         + pMG_fine->Phi[2*k-1][2*j-1][2*i] + pMG_fine->Phi[2*k-1][2*j-1][2*i-1])/8.0 ;
        pMG_coarse->rhs[k][j][i] =
           (error[2*k  ][2*j  ][2*i] + error[2*k  ][2*j  ][2*i-1]
          + error[2*k  ][2*j-1][2*i] + error[2*k  ][2*j-1][2*i-1]
          + error[2*k-1][2*j  ][2*i] + error[2*k-1][2*j  ][2*i-1]
          + error[2*k-1][2*j-1][2*i] + error[2*k-1][2*j-1][2*i-1])/8.0 ;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void Prolongation_3d(MGrid *pMG_coarse, MGrid *pMG_fine)
 *  \brief Linear interpolation of coarse grid onto fine  
 */
void Prolongation_3d(MGrid *pMG_coarse, MGrid *pMG_fine)
{
  int i, is = pMG_coarse->is, ie = pMG_coarse->ie;
  int j, js = pMG_coarse->js, je = pMG_coarse->je;
  int k, ks = pMG_coarse->ks, ke = pMG_coarse->ke;

  for (k=ks; k<=ke; k++){
  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      pMG_fine->Phi[2*k  ][2*j  ][2*i  ] += 0.75*pMG_coarse->Phi[k  ][j  ][i  ];
      pMG_fine->Phi[2*k  ][2*j  ][2*i  ] += 0.25*pMG_coarse->Phi[k+1][j+1][i+1];

      pMG_fine->Phi[2*k  ][2*j  ][2*i-1] += 0.75*pMG_coarse->Phi[k  ][j  ][i  ];
      pMG_fine->Phi[2*k  ][2*j  ][2*i-1] += 0.25*pMG_coarse->Phi[k+1][j+1][i-1];

      pMG_fine->Phi[2*k  ][2*j-1][2*i  ] += 0.75*pMG_coarse->Phi[k  ][j  ][i  ];
      pMG_fine->Phi[2*k  ][2*j-1][2*i  ] += 0.25*pMG_coarse->Phi[k+1][j-1][i+1];

      pMG_fine->Phi[2*k-1][2*j  ][2*i  ] += 0.75*pMG_coarse->Phi[k  ][j  ][i  ];
      pMG_fine->Phi[2*k-1][2*j  ][2*i  ] += 0.25*pMG_coarse->Phi[k-1][j+1][i+1];

      pMG_fine->Phi[2*k  ][2*j-1][2*i-1] += 0.75*pMG_coarse->Phi[k  ][j  ][i  ];
      pMG_fine->Phi[2*k  ][2*j-1][2*i-1] += 0.25*pMG_coarse->Phi[k+1][j-1][i-1];

      pMG_fine->Phi[2*k-1][2*j-1][2*i  ] += 0.75*pMG_coarse->Phi[k  ][j  ][i  ];
      pMG_fine->Phi[2*k-1][2*j-1][2*i  ] += 0.25*pMG_coarse->Phi[k-1][j-1][i+1];

      pMG_fine->Phi[2*k-1][2*j]  [2*i-1] += 0.75*pMG_coarse->Phi[k  ][j  ][i  ];
      pMG_fine->Phi[2*k-1][2*j]  [2*i-1] += 0.25*pMG_coarse->Phi[k-1][j+1][i-1];

      pMG_fine->Phi[2*k-1][2*j-1][2*i-1] += 0.75*pMG_coarse->Phi[k  ][j  ][i  ];
      pMG_fine->Phi[2*k-1][2*j-1][2*i-1] += 0.25*pMG_coarse->Phi[k-1][j-1][i-1];
    }
  }}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void set_mg_bvals(MGrid *pMG)
 *  \brief Sets BC for Jacobi iterates for MPI parallel jobs.
 *
 *   With self-gravity using multigrid, the boundary conditions at the edge of
 *   the Domain are held fixed.  So only ghostzones associated with internal
 *   boundaries between MPI grids need to be passed.
 *
 * This routine is largely a copy of set_bvals().
 * Order for updating boundary conditions must always be x1-x2-x3 in order to
 * fill the corner cells properly
 */

#ifdef MPI_PARALLEL
void set_mg_bvals(MGrid *pMG)
{
  int cnt3, cnt, err;
  MPI_Status stat;
  MPI_Request rq;

/*--- Step 1. ------------------------------------------------------------------
 * Boundary Conditions in x1-direction */

  cnt3 = 1;
  if (pMG->Nx3 > 1) cnt3 = pMG->Nx3;
  cnt = pMG->Nx2*cnt3;

/* MPI blocks to both left and right */
  if (pMG->rx1_id >= 0 && pMG->lx1_id >= 0) {
    /* Post a non-blocking receive for the input data from the left grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->lx1_id,
      boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_mg_bvals]: MPI_Irecv error = %d\n",err);

    swap_mg_ox1(pMG,cnt,0,&rq);  /* send R */
    swap_mg_ix1(pMG,cnt,1,&rq);  /* listen L */

    /* Post a non-blocking receive for the input data from the right grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->rx1_id,
      boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_mg_bvals]: MPI_Irecv error = %d\n",err);

    swap_mg_ix1(pMG,cnt,0,&rq);  /* send L */
    swap_mg_ox1(pMG,cnt,1,&rq);  /* listen R */
  }

/* Physical boundary on left, MPI block on right */
  if (pMG->rx1_id >= 0 && pMG->lx1_id < 0) {
    /* Post a non-blocking receive for the input data from the right grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->rx1_id,
      boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_mg_bvals]: MPI_Irecv error = %d\n",err);

    swap_mg_ox1(pMG,cnt,0,&rq);  /* send R */
    swap_mg_ox1(pMG,cnt,1,&rq);  /* listen R */
  }

/* MPI block on left, Physical boundary on right */
  if (pMG->rx1_id < 0 && pMG->lx1_id >= 0) {
    /* Post a non-blocking receive for the input data from the left grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->lx1_id,
      boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_mg_bvals]: MPI_Irecv error = %d\n",err);

    swap_mg_ix1(pMG,cnt,0,&rq);  /* send L */
    swap_mg_ix1(pMG,cnt,1,&rq);  /* listen L */
  }

/*--- Step 2. ------------------------------------------------------------------
 * Boundary Conditions in x2-direction */

  cnt3 = 1;
  if (pMG->Nx3 > 1) cnt3 = pMG->Nx3;
  cnt = (pMG->Nx1 + 2)*cnt3;

/* MPI blocks to both left and right */
  if (pMG->rx2_id >= 0 && pMG->lx2_id >= 0) {
    /* Post a non-blocking receive for the input data from the left grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->lx2_id,
      boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_mg_bvals]: MPI_Irecv error = %d\n",err);

    swap_mg_ox2(pMG,cnt,0,&rq);  /* send R */
    swap_mg_ix2(pMG,cnt,1,&rq);  /* listen L */

    /* Post a non-blocking receive for the input data from the right grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->rx2_id,
      boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_mg_bvals]: MPI_Irecv error = %d\n",err);

    swap_mg_ix2(pMG,cnt,0,&rq);  /* send L */
    swap_mg_ox2(pMG,cnt,1,&rq);  /* listen R */
  }

/* Physical boundary on left, MPI block on right */
  if (pMG->rx2_id >= 0 && pMG->lx2_id < 0) {
    /* Post a non-blocking receive for the input data from the right grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->rx2_id,
      boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_mg_bvals]: MPI_Irecv error = %d\n",err);

    swap_mg_ox2(pMG,cnt,0,&rq);  /* send R */
    swap_mg_ox2(pMG,cnt,1,&rq);  /* listen R */
  }

/* MPI block on left, Physical boundary on right */
  if (pMG->rx2_id < 0 && pMG->lx2_id >= 0) {
    /* Post a non-blocking receive for the input data from the left grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->lx2_id,
      boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_mg_bvals]: MPI_Irecv error = %d\n",err);

    swap_mg_ix2(pMG,cnt,0,&rq);  /* send L */
    swap_mg_ix2(pMG,cnt,1,&rq);  /* listen L */
  }

/*--- Step 3. ------------------------------------------------------------------
 * Boundary Conditions in x3-direction */

  if (pMG->Nx3 > 1){

    cnt = (pMG->Nx1 + 2)*(pMG->Nx2 + 2);

/* MPI blocks to both left and right */
    if (pMG->rx3_id >= 0 && pMG->lx3_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->lx3_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_mg_bvals]: MPI_Irecv error = %d\n",err);

      swap_mg_ox3(pMG,cnt,0,&rq);  /* send R */
      swap_mg_ix3(pMG,cnt,1,&rq);  /* listen L */

      /* Post a non-blocking receive for the input data from the right grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->rx3_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_mg_bvals]: MPI_Irecv error = %d\n",err);

      swap_mg_ix3(pMG,cnt,0,&rq);  /* send L */
      swap_mg_ox3(pMG,cnt,1,&rq);  /* listen R */
    }

/* Physical boundary on left, MPI block on right */
    if (pMG->rx3_id >= 0 && pMG->lx3_id < 0) {
      /* Post a non-blocking receive for the input data from the right grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->rx3_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_mg_bvals]: MPI_Irecv error = %d\n",err);

      swap_mg_ox3(pMG,cnt,0,&rq);  /* send R */
      swap_mg_ox3(pMG,cnt,1,&rq);  /* listen R */
    }

/* MPI block on left, Physical boundary on right */
    if (pMG->rx3_id < 0 && pMG->lx3_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->lx3_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_mg_bvals]: MPI_Irecv error = %d\n",err);

      swap_mg_ix3(pMG,cnt,0,&rq);  /* send L */
      swap_mg_ix3(pMG,cnt,1,&rq);  /* listen L */
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void swap_mg_ix1(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
 *  \brief MPI_SWAP of boundary conditions, Inner x1 boundary
 *
 *   This function does either a send (swap_flag=0), or receive (swap_flag=1)
 *   Largely copied from set_bvals/send_ix1 and set_bvals/receive_ix1
 */

void swap_mg_ix1(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
{
  int j,jl,ju,k,kl,ku,err;
  MPI_Status stat;
  double *psb = send_buf;
  double *prb = recv_buf;

  jl = pMG->js;
  ju = pMG->je;

  if(pMG->Nx3 > 1){
    kl = pMG->ks;
    ku = pMG->ke;
  } else {
    kl = ku = pMG->ks;
  }

/* Pack single row i=is of iterates into send buffer */

  if (swap_flag == 0) {
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        *(psb++) = pMG->Phi[k][j][pMG->is];
      }
    }

    /* send contents of buffer to the neighboring grid on L-x1 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pMG->lx1_id,
                   boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[swap_mg_ix1]: MPI_Send error = %d\n",err);
  }

/* Receive message and unpack single row i=is-1 of iteratas into ghost zonee */

  if (swap_flag == 1) {

    /* Wait to receive the input data from the left grid */
    err = MPI_Wait(prq, &stat);
    if(err) ath_error("[swap_mg_ix1]: MPI_Wait error = %d\n",err);

    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        pMG->Phi[k][j][pMG->is-1] = *(prb++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void swap_mg_ox1(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
 *  \brief MPI_SWAP of boundary conditions, Outer x1 boundary
 *
 *   This function does either a send (swap_flag=0), or receive (swap_flag=1)
 *   Largely copied from set_bvals/send_ox1 and set_bvals/receive_ox1
 */

void swap_mg_ox1(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
{
  int j,jl,ju,k,kl,ku,err;
  MPI_Status stat;
  double *psb = send_buf;
  double *prb = recv_buf;

  jl = pMG->js;
  ju = pMG->je;

  if(pMG->Nx3 > 1){
    kl = pMG->ks;
    ku = pMG->ke;
  } else {
    kl = ku = pMG->ks;
  }

/* Pack single row i=ie of iterates into send buffer */

  if (swap_flag == 0) {
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        *(psb++) = pMG->Phi[k][j][pMG->ie];
      }
    }

    /* send contents of buffer to the neighboring grid on R-x1 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pMG->rx1_id,
      boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[swap_mg_ox1]: MPI_Send error = %d\n",err);
  }

/* Receive message and unpack single row i=ie+1 of iterates into ghost zone */

  if (swap_flag == 1) {

    /* Wait to receive the input data from the right grid */
    err = MPI_Wait(prq, &stat);
    if(err) ath_error("[swap_mg_ox1]: MPI_Wait error = %d\n",err);

    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        pMG->Phi[k][j][pMG->ie+1] = *(prb++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void swap_mg_ix2(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
 *  \brief MPI_SWAP of boundary conditions, Inner x2 boundary
 *
 *   This function does either a send (swap_flag=0), or receive (swap_flag=1)
 *   Largely copied from set_bvals/send_ix2 and set_bvals/receive_ix2
 */

void swap_mg_ix2(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
{
  int i,il,iu,k,kl,ku,err;
  MPI_Status stat;
  double *pd = send_buf;

  il = pMG->is - 1;
  iu = pMG->ie + 1;

  if(pMG->Nx3 > 1){
    kl = pMG->ks;
    ku = pMG->ke;
  } else {
    kl = ku = pMG->ks;
  }

/* Pack single row j=js of iterates into send buffer */

  if (swap_flag == 0) {
    for (k=kl; k<=ku; k++){
      for (i=il; i<=iu; i++){
        *(pd++) = pMG->Phi[k][pMG->js][i];
      }
    }

    /* send contents of buffer to the neighboring grid on L-x2 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pMG->lx2_id,
       boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[swap_mg_ix2]: MPI_Send error = %d\n",err);
  }

/* Receive message and unpack single row j=js-1 of iterates into ghost zone */

  if (swap_flag == 1) {

    /* Wait to receive the input data from the left grid */
    err = MPI_Wait(prq, &stat);
    if(err) ath_error("[swap_mg_ix2]: MPI_Wait error = %d\n",err);

    for (k=kl; k<=ku; k++){
      for (i=il; i<=iu; i++){
        pMG->Phi[k][pMG->js-1][i] = *(pd++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void swap_mg_ox2(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
 *  \brief MPI_SWAP of boundary conditions, Outer x2 boundary
 *
 *   This function does either a send (swap_flag=0), or receive (swap_flag=1)
 *   Largely copied from set_bvals/send_ix2 and set_bvals/receive_ix2
 */

void swap_mg_ox2(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
{
  int i,il,iu,k,kl,ku,err;
  MPI_Status stat;
  double *pd = send_buf;

  il = pMG->is - 1;
  iu = pMG->ie + 1;

  if(pMG->Nx3 > 1){
    kl = pMG->ks;
    ku = pMG->ke;
  } else {
    kl = ku = pMG->ks;
  }

/* Pack single row j=je of iterates into send buffer */

  if (swap_flag == 0) {
    for (k=kl; k<=ku; k++){
      for (i=il; i<=iu; i++){
        *(pd++) = pMG->Phi[k][pMG->je][i];
      }
    }

    /* send contents of buffer to the neighboring grid on R-x2 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pMG->rx2_id,
                   boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[swap_mg_ox2]: MPI_Send error = %d\n",err);
  }

/* Receive message and unpack single row j=je+1 of iterates into ghost zone */

  if (swap_flag == 1) {

    /* Wait to receive the input data from the right grid */
    err = MPI_Wait(prq, &stat);
    if(err) ath_error("[swap_mg_ox2]: MPI_Wait error = %d\n",err);

    for (k=kl; k<=ku; k++){
      for (i=il; i<=iu; i++){
        pMG->Phi[k][pMG->je+1][i] = *(pd++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void swap_mg_ix3(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
 *  \brief MPI_SWAP of boundary conditions, Inner x3 boundary
 *
 *   This function does either a send (swap_flag=0), or receive (swap_flag=1)
 *   Largely copied from set_bvals/send_ix3 and set_bvals/receive_ix3
 */

void swap_mg_ix3(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,err;
  MPI_Status stat;
  double *pd = send_buf;

  il = pMG->is - 1;
  iu = pMG->ie + 1;
  jl = pMG->js - 1;
  ju = pMG->je + 1;

/* Pack single row k=ks of iterates into send buffer */

  if (swap_flag == 0) {
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        *(pd++) = pMG->Phi[pMG->ks][j][i];
      }
    }

    /* send contents of buffer to the neighboring grid on L-x3 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pMG->lx3_id,
                   boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[swap_mg_ix3]: MPI_Send error = %d\n",err);
  }

/* Receive message and unpack single row k=ks-1 of iterates into ghost zone */

  if (swap_flag == 1) {

    /* Wait to receive the input data from the left grid */
    err = MPI_Wait(prq, &stat);
    if(err) ath_error("[swap_mg_ix3]: MPI_Wait error = %d\n",err);

    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        pMG->Phi[pMG->ks-1][j][i] = *(pd++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void swap_mg_ox3(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
 *  \brief MPI_SWAP of boundary conditions, Outer x3 boundary
 *
 *   This function does either a send (swap_flag=0), or receive (swap_flag=1)
 *   Largely copied from set_bvals/send_ix1 and set_bvals/receive_ix1
 */

void swap_mg_ox3(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,err;
  MPI_Status stat;
  double *pd = send_buf;

  il = pMG->is - 1;
  iu = pMG->ie + 1;
  jl = pMG->js - 1;
  ju = pMG->je + 1;

/* Pack single row k=ke of interates into send buffer */

  if (swap_flag == 0) {
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        *(pd++) = pMG->Phi[pMG->ke][j][i];
      }
    }

    /* send contents of buffer to the neighboring grid on R-x3 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pMG->rx3_id,
                   boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[swap_mg_ox3]: MPI_Send error = %d\n",err);
  }

/* Receive message and unpack single row k=ke+1 of iterates into ghost zone */

  if (swap_flag == 1) {

    /* Wait to receive the input data from the right grid */
    err = MPI_Wait(prq, &stat);
    if(err) ath_error("[swap_mg_ox3]: MPI_Wait error = %d\n",err);

    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        pMG->Phi[pMG->ke+1][j][i] = *(pd++);
      }
    }
  }

  return;
}
#endif /* MPI_PARALLEL */

/*----------------------------------------------------------------------------*/
/*! \fn void selfg_multig_3d_init(MeshS *pM)
 *  \brief Initialize send/receive buffers needed to swap
 *   iterates during Jacobi iterations.
 */

void selfg_multig_3d_init(MeshS *pM)
{
#ifdef MPI_PARALLEL
  int i,j,k,x1cnt,x2cnt,x3cnt;
  int nx1t,nx2t,nx3t, size;
  int NGrid_x1, NGrid_x2, NGrid_x3;
#endif
  int nmax,size1=0,size2=0,size3=0,nl,nd;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2, Nx3 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if (pM->Domain[nl][nd].Grid->Nx[0] > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
        if (pM->Domain[nl][nd].Grid->Nx[1] > size2){
          size2 = pM->Domain[nl][nd].Grid->Nx[1];
        }
        if (pM->Domain[nl][nd].Grid->Nx[2] > size3){
          size3 = pM->Domain[nl][nd].Grid->Nx[2];
        }
      }
    }
  }

  size1 += 2;
  size2 += 2;
  size3 += 2;
  error = (Real ***)calloc_3d_array(size3, size2, size1, sizeof(Real));
  if (error == NULL) {
    ath_error("[multig_3d]: Error allocating memory for error array\n");
  }

/* Allocate memory for send and receive buffers for Phi in MultiGrid
 * structure for MPI parallel.
 */
#ifdef MPI_PARALLEL
  NGrid_x1 = par_geti("parallel","NGrid_x1");
  NGrid_x2 = par_geti("parallel","NGrid_x2");
  NGrid_x3 = par_geti("parallel","NGrid_x3");

  x1cnt = x2cnt = x3cnt = 0;
        
  for (k=0; k<NGrid_x3; k++){
    for (j=0; j<NGrid_x2; j++){
      for (i=0; i<NGrid_x1; i++){
        if(NGrid_x1 > 1){
          nx2t = pD->grid_block[k][j][i].jde - pD->grid_block[k][j][i].jds + 1;
          nx3t = pD->grid_block[k][j][i].kde - pD->grid_block[k][j][i].kds + 1;
        
          x1cnt = nx2t*nx3t > x1cnt ? nx2t*nx3t : x1cnt;
        }
        
        if(NGrid_x2 > 1){
          nx1t = pD->grid_block[k][j][i].ide - pD->grid_block[k][j][i].ids + 1;
          if(nx1t > 1) nx1t += 2;
          nx3t = pD->grid_block[k][j][i].kde - pD->grid_block[k][j][i].kds + 1;

          x2cnt = nx1t*nx3t > x2cnt ? nx1t*nx3t : x2cnt;
        }


        if(NGrid_x3 > 1){
          nx1t = pD->grid_block[k][j][i].ide - pD->grid_block[k][j][i].ids + 1;
          if(nx1t > 1) nx1t += 2;
          nx2t = pD->grid_block[k][j][i].jde - pD->grid_block[k][j][i].jds + 1;
          if(nx2t > 1) nx2t += 2;

          x3cnt = nx1t*nx2t > x3cnt ? nx1t*nx2t : x3cnt;
        }
      }
    }
  }

  size = x1cnt > x2cnt ? x1cnt : x2cnt;
  size = x3cnt >  size ? x3cnt : size;

  if (size > 0) {
    if((send_buf = (double*)malloc(size*sizeof(double))) == NULL)
      ath_error("[selfg_by_multig_3d_init]: Failed to allocate send buffer\n");

    if((recv_buf = (double*)malloc(size*sizeof(double))) == NULL)
      ath_error("[selfg_by_multig_3d_init]: Failed to allocate recv buffer\n");
  }
#endif /* MPI_PARALLEL */
  return;
}

#endif /* SELF_GRAVITY_USING_MULTIGRID */
