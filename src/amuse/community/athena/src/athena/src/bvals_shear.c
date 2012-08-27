#include "copyright.h"
/*============================================================================*/
/*! \file bvals_shear.c
 *  \brief Shearing sheet boundary conditions at ix1 and ox1 for both 2D and 3D
 *
 * PURPOSE: Shearing sheet boundary conditions at ix1 and ox1 for both 2D and 3D
 *   Called by bvals_mhd.  Decomposition of the Domain into MPI grids in X,Y
 *   and/or Z is allowed.  The RemapEy() function (which applies the shearing
 *   sheet boundary conditions to the y-component of the EMF to keep <Bz>=const)
 *   is called directly by the 3D integrator.  Configure code with
 *   --enable-shearing-box to use.
 *
 * FARGO (orbital advection) algorithm is implemented in the Fargo() function
 *   called in the main loop directly after the integrator and before bvals_mhd.
 *   Code must be configured with --enable-shearing-box --enable-fargo to use
 *   this option.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - ShearingSheet_ix1() - shearing sheet BCs on ix1
 * - ShearingSheet_ox1() - shearing sheet BCs on ox1
 * - RemapEy_ix1()       - sets Ey at ix1 in integrator to keep <Bz>=const. 
 * - RemapEy_ox1()       - sets Ey at ox1 in integrator to keep <Bz>=const. 
 * - Fargo()             - implements FARGO algorithm for background flow
 * - bvals_shear_init() - allocates memory for arrays used here
 * - bvals_shear_destruct() - frees memory for arrays used here		      
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - RemapFlux() - 2nd or 3rd order reconstruction for remap in ghost zones   */
/*============================================================================*/

#include <float.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/* The functions in this file will only work with the shearing box */
#if defined(SHEARING_BOX) || (defined(CYLINDRICAL) && defined(FARGO))

/* Define number of variables to be remapped */
#ifdef BAROTROPIC /* BAROTROPIC EOS */
#ifdef HYDRO
 enum {NREMAP = 4, NFARGO = 4};
#endif
#ifdef MHD
 enum {NREMAP = 8, NFARGO = 6};
#endif
#else /* ADIABATIC or other EOS */
#ifdef HYDRO
 enum {NREMAP = 5, NFARGO = 5};
#endif
#ifdef MHD
 enum {NREMAP = 9, NFARGO = 7};
#endif
#endif /* EOS */

#ifdef MHD
#define NVAR_SHARE (NVAR + 3)
#else
#define NVAR_SHARE NVAR
#endif

/*! \struct Remap
 *  \brief Define structure which holds variables remapped 
 *  by shearing sheet BCs */
typedef struct Remap_s{
  Real U[NREMAP];
#if (NSCALARS > 0)
  Real s[NSCALARS];
#endif
}Remap;

/*! \struct FConsS
 *  \brief Define structure for variables used in FARGO algorithm */
typedef struct FCons_s{
  Real U[NFARGO];
#if (NSCALARS > 0)
  Real s[NSCALARS];
#endif
}FConsS;

/* The memory for all the arrays below is allocated in bvals_shear_init */
/* Arrays of ghost zones containing remapped conserved quantities */
static Remap ***GhstZns=NULL, ***GhstZnsBuf=NULL;
/* 1D vectors for reconstruction in conservative remap step */
static Real *U=NULL, *Flx=NULL;
/* Arrays of Ey remapped at ix1/ox1 edges of Domain */
#ifdef MHD
static Real **tEyBuf=NULL;
#endif
/* temporary vector needed for 3rd order reconstruction in ghost zones */
#if defined(THIRD_ORDER_CHAR) || defined(THIRD_ORDER_PRIM)
static Real *Uhalf=NULL;
#endif
/* MPI send and receive buffers */
#ifdef MPI_PARALLEL
static double *send_buf = NULL, *recv_buf = NULL;
#endif
#ifdef FARGO
int nfghost;
static FConsS ***FargoVars=NULL, ***FargoFlx=NULL;
#endif

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * RemapFlux() - 2nd or 3rd order reconstruction for remap in ghost zones
 *============================================================================*/

void RemapFlux(const Real *U,const Real eps,const int ji,const int jo, Real *F);

#endif

#ifdef SHEARING_BOX
/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void ShearingSheet_ix1(DomainS *pD)
 *  \brief 3D shearing-sheet BCs in x1.  
 *
 * It applies a remap
 * in Y after the ghost cells have been set by the usual periodic BCs in X and
 * Y implemented in bvals_mhd.c
 *
 * This is a public function which is called by bvals_mhd() inside a
 * SHEARING_BOX macro.							      */
/*----------------------------------------------------------------------------*/
void ShearingSheet_ix1(DomainS *pD)
{
  GridS *pG = pD->Grid;
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,ii,j,k,ku,n,joffset,jremap;
  Real xmin,xmax,Lx,Ly,qomL,yshear,deltay,epsi;
#ifdef MPI_PARALLEL
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids;
  int ierr,sendto_id,getfrom_id;
  double *pSnd,*pRcv;
  Remap *pRemap;
  ConsS *pCons;
  MPI_Request rq;
#endif

/*--- Step 1. ------------------------------------------------------------------
 * Compute the distance the computational domain has sheared in y */

  xmin = pD->RootMinX[0];
  xmax = pD->RootMaxX[0];
  Lx = xmax - xmin;

  xmin = pD->RootMinX[1];
  xmax = pD->RootMaxX[1];
  Ly = xmax - xmin;

  qomL = qshear*Omega_0*Lx;
  yshear = qomL*pG->time;

/* Split this into integer and fractional peices of the Domain in y.  Ignore
 * the integer piece because the Grid is periodic in y */

  deltay = fmod(yshear, Ly);

/* further decompose the fractional peice into integer and fractional pieces of
 * a grid cell.  Note 0.0 <= epsi < 1.0.  If Domain has MPI decomposition in Y,
 * then it is possible that:  pD->Nx2 > joffset > pG->Nx2   */

  joffset = (int)(deltay/pG->dx2);
  epsi = (fmod(deltay,pG->dx2))/pG->dx2;

/*--- Step 2. ------------------------------------------------------------------
 * Copy data into GhstZns array.  Note i and j indices are switched.
 * Steps 2-10 are for 3D or 2d xy runs.  Step 11 handles 2D xz separately */

  if (pG->Nx[2] > 1 || ShBoxCoord==xy) {  /* this if ends at end of step 10 */

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for(k=ks; k<=ku; k++) {
    for(j=js-nghost; j<=je+nghost; j++){
      for(i=0; i<nghost; i++){
        ii = is-nghost+i;
        GhstZns[k][i][j].U[0] = pG->U[k][j][ii].d;
        GhstZns[k][i][j].U[1] = pG->U[k][j][ii].M1;
        GhstZns[k][i][j].U[2] = pG->U[k][j][ii].M2;
#ifndef FARGO
        GhstZns[k][i][j].U[2] += qomL*pG->U[k][j][ii].d;
#endif
        GhstZns[k][i][j].U[3] = pG->U[k][j][ii].M3;
#ifdef ADIABATIC
/* No change in the internal energy */
        GhstZns[k][i][j].U[4] = pG->U[k][j][ii].E + (0.5/GhstZns[k][i][j].U[0])*
          (SQR(GhstZns[k][i][j].U[2]) - SQR(pG->U[k][j][ii].M2));
#endif /* ADIABATIC */
#ifdef MHD
        GhstZns[k][i][j].U[NREMAP-4] = pG->U[k][j][ii].B1c;
        GhstZns[k][i][j].U[NREMAP-3] = pG->B1i[k][j][ii];
        GhstZns[k][i][j].U[NREMAP-2] = pG->B2i[k][j][ii];
        GhstZns[k][i][j].U[NREMAP-1] = pG->B3i[k][j][ii];
#endif /* MHD */
#if (NSCALARS > 0)
        for(n=0; n<NSCALARS; n++) GhstZns[k][i][j].s[n] = pG->U[k][j][ii].s[n];
#endif
      }
    }
  }

/*--- Step 3. ------------------------------------------------------------------
 * Copy GhstZns into buffer, at the same time apply a conservative remap of
 * solution over the fractional part of grid cell */

  for(k=ks; k<=ku; k++) {
    for(i=0; i<nghost; i++){

      for (n=0; n<(NREMAP); n++) {
        for (j=js-nghost; j<=je+nghost; j++) U[j] = GhstZns[k][i][j].U[n];
        RemapFlux(U,epsi,js,je+1,Flx);
        for(j=js; j<=je; j++){
          GhstZnsBuf[k][i][j].U[n] = GhstZns[k][i][j].U[n] - (Flx[j+1]-Flx[j]);
        }
      }

#if (NSCALARS > 0)
      for (n=0; n<(NSCALARS); n++) {
        for (j=js-nghost; j<=je+nghost; j++) U[j] = GhstZns[k][i][j].s[n];
        RemapFlux(U,epsi,js,je+1,Flx);
        for(j=js; j<=je; j++){
          GhstZnsBuf[k][i][j].s[n] = GhstZns[k][i][j].s[n] - (Flx[j+1]-Flx[j]);
        }
      }
#endif

    }
  }

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from buffer back into GhstZns.  */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ku; k++) {
      for(j=js; j<=je; j++){
        jremap = j - joffset;
        if (jremap < (int)js) jremap += pG->Nx[1];

        for(i=0; i<nghost; i++){
          for (n=0; n<(NREMAP); n++) {
            GhstZns[k][i][j].U[n]  = GhstZnsBuf[k][i][jremap].U[n];
          }
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) { 
            GhstZns[k][i][j].s[n] = GhstZnsBuf[k][i][jremap].s[n];
          }
#endif
        }

      }
    }

#ifdef MPI_PARALLEL
  } else {

/*--- Step 5. ------------------------------------------------------------------
 * If Domain contains MPI decomposition in Y, then MPI calls are required for
 * the cyclic shift needed to apply shift over integer number of grid cells
 * during copy from buffer back into GhstZns.  */

    get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);

/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
    Ngrids = (int)(joffset/pG->Nx[1]);
    joverlap = joffset - Ngrids*pG->Nx[1];

/*--- Step 5a. -----------------------------------------------------------------
 * Find ids of processors that data in [je-(joverlap-1):je] is sent to, and
 * data in [js:js+(joverlap-1)] is received from.  Only execute if joverlap>0 */
/* This can result in send/receive to self -- we rely on MPI to handle this
 * properly */

    if (joverlap != 0) {

      jproc = my_jproc + (Ngrids + 1);
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1]; 
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc - (Ngrids + 1);
      if (jproc < 0) jproc += pD->NGrid[1]; 
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

/*--- Step 5b. -----------------------------------------------------------------
 * Pack send buffer and send data in [je-(joverlap-1):je] from GhstZnsBuf */

      cnt = nghost*joverlap*(ku-ks+1)*(NREMAP+NSCALARS);
/* Post a non-blocking receive for the input data */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_ix1_tag, pD->Comm_Domain, &rq);

      pSnd = send_buf;
      for (k=ks; k<=ku; k++) {
        for (j=je-(joverlap-1); j<=je; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZnsBuf[k][i][j]);

            for (n=0; n<NREMAP; n++) *(pSnd++) = pRemap->U[n];
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) *(pSnd++) = pRemap->s[n];
#endif
          }
        }
      }
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     shearing_sheet_ix1_tag, pD->Comm_Domain);

/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [je-(joverlap-1):je], and remap into cells in
 * [js:js+(joverlap-1)] in GhstZns */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for (k=ks; k<=ku; k++) {
        for (j=js; j<=js+(joverlap-1); j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZns[k][i][j]);
  
            for (n=0; n<NREMAP; n++) pRemap->U[n] = *(pRcv++);
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) pRemap->s[n] = *(pRcv++);
#endif
          }
        }
      }

    }

/*--- Step 5d. -----------------------------------------------------------------
 * If shear is less one full Grid, remap cells which remain on same processor
 * from GhstZnsBuf into GhstZns.  Cells in [js:je-joverlap] are shifted by
 * joverlap into [js+joverlap:je] */

    if (Ngrids == 0) {

      for(k=ks; k<=ku; k++) {
        for(j=js+joverlap; j<=je; j++){
          jremap = j-joverlap;
          for(i=0; i<nghost; i++){
            for (n=0; n<(NREMAP); n++) {
              GhstZns[k][i][j].U[n]  = GhstZnsBuf[k][i][jremap].U[n];
            }
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) { 
              GhstZns[k][i][j].s[n] = GhstZnsBuf[k][i][jremap].s[n];
            }
#endif
          }
        }
      }

/*--- Step 5e. -----------------------------------------------------------------
 * If shear is more than one Grid, pack and send data from [js:je-joverlap]
 * from GhstZnsBuf (this step replaces 5d) */

    } else {

/* index of sendto and getfrom processors in GData are -/+1 from Step 5a */

      jproc = my_jproc + Ngrids;
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1];
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc - Ngrids;
      if (jproc < 0) jproc += pD->NGrid[1];
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      cnt = nghost*(pG->Nx[1]-joverlap)*(ku-ks+1)*(NREMAP+NSCALARS);
/* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_ix1_tag, pD->Comm_Domain, &rq);

      pSnd = send_buf;
      for (k=ks; k<=ku; k++) {
        for (j=js; j<=je-joverlap; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZnsBuf[k][i][j]);
            for (n=0; n<NREMAP; n++) *(pSnd++) = pRemap->U[n];
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) *(pSnd++) = pRemap->s[n];
#endif
          }
        }
      }
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     shearing_sheet_ix1_tag, pD->Comm_Domain);

/* unpack data sent from [js:je-overlap], and remap into cells in
 * [js+joverlap:je] in GhstZns */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for (k=ks; k<=ku; k++) {
        for (j=js+joverlap; j<=je; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZns[k][i][j]);
            for (n=0; n<NREMAP; n++) pRemap->U[n] = *(pRcv++);
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) pRemap->s[n] = *(pRcv++);
#endif
          }
        }
      }
    } /* end of step 5e - shear is more than one Grid */

#endif /* MPI_PARALLEL */
  } /* end of step 5 - MPI decomposition in Y */

/*--- Step 6. ------------------------------------------------------------------
 * Now copy remapped variables back into ghost cells */

  for(k=ks; k<=ke; k++) {
    for(j=js; j<=je; j++){
      for(i=0; i<nghost; i++){
        pG->U[k][j][is-nghost+i].d  = GhstZns[k][i][j].U[0];
        pG->U[k][j][is-nghost+i].M1 = GhstZns[k][i][j].U[1];
        pG->U[k][j][is-nghost+i].M2 = GhstZns[k][i][j].U[2];
        pG->U[k][j][is-nghost+i].M3 = GhstZns[k][i][j].U[3];
#ifdef ADIABATIC
        pG->U[k][j][is-nghost+i].E  = GhstZns[k][i][j].U[4];
#endif /* ADIABATIC */
#ifdef MHD
        pG->U[k][j][is-nghost+i].B1c = GhstZns[k][i][j].U[NREMAP-4];
        pG->B1i[k][j][is-nghost+i] = GhstZns[k][i][j].U[NREMAP-3];
        pG->B2i[k][j][is-nghost+i] = GhstZns[k][i][j].U[NREMAP-2];
        pG->B3i[k][j][is-nghost+i] = GhstZns[k][i][j].U[NREMAP-1];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) {
          pG->U[k][j][is-nghost+i].s[n] = GhstZns[k][i][j].s[n];
        }
#endif
      }
    }
  }

/* Copy the face-centered B3 component of the field at k=ke+1 in 3D */
#ifdef MHD
  if (pG->Nx[2] > 1) {
    for(j=js; j<=je; j++){
      for(i=0; i<nghost; i++){
        pG->B3i[ke+1][j][is-nghost+i] = GhstZns[ke+1][i][j].U[NREMAP-1];
      }
    }
  }
#endif /* MHD */

/*--- Step 7. ------------------------------------------------------------------
 * compute cell-centered B as average of remapped face centered B, except B1.
 * The value of B2c at j=je is incorrect since B2i[je+1] not yet set -- fix in
 * step 10 below */

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for(i=is-nghost; i<is; i++){
        pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i]+pG->B2i[k][j+1][i]);
      }
    }
  }
  if (pG->Nx[2] > 1) {
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for(i=is-nghost; i<is; i++){
          pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i]+pG->B3i[k+1][j][i]);
        }
      }
    }
  }
#endif /* MHD */

/*--- Step 8. ------------------------------------------------------------------
 * With no MPI decomposition in Y, apply periodic BCs in Y (similar to
 * periodic_ix2() and periodic_ox2() in bvals_mhd.c) */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ke; k++) {
      for(j=1; j<=nghost; j++){
        for(i=is-nghost; i<is; i++){
          pG->U[k][js-j][i] = pG->U[k][je-(j-1)][i];
          pG->U[k][je+j][i] = pG->U[k][js+(j-1)][i];
#ifdef MHD
          pG->B1i[k][js-j][i] = pG->B1i[k][je-(j-1)][i];
          pG->B2i[k][js-j][i] = pG->B2i[k][je-(j-1)][i];
          pG->B3i[k][js-j][i] = pG->B3i[k][je-(j-1)][i];

          pG->B1i[k][je+j][i] = pG->B1i[k][js+(j-1)][i];
          pG->B2i[k][je+j][i] = pG->B2i[k][js+(j-1)][i];
          pG->B3i[k][je+j][i] = pG->B3i[k][js+(j-1)][i];
#endif /* MHD */
        }
      }
    }
#ifdef MHD
    if (pG->Nx[2] > 1) {
      for (j=1; j<=nghost; j++) {
        for (i=is-nghost; i<is; i++) {
          pG->B3i[ke+1][js-j][i] = pG->B3i[ke+1][je-(j-1)][i];
          pG->B3i[ke+1][je+j][i] = pG->B3i[ke+1][js+(j-1)][i];
        }
      }
    }
#endif /* MHD */

#ifdef MPI_PARALLEL
  } else {

/*--- Step 9. ------------------------------------------------------------------
 * With MPI decomposition in Y, use MPI calls to handle periodic BCs in Y (like
 * send_ox2/receive_ix1 and send_ix1/receive_ox2 pairs in bvals_mhd.c */


/* Post a non-blocking receive for the input data from the left grid */
    cnt = nghost*nghost*(ku-ks+1)*NVAR_SHARE;
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                    shearing_sheet_ix1_tag, pD->Comm_Domain, &rq);

    pSnd = send_buf;
    for (k=ks; k<=ku; k++){
      for (j=je-nghost+1; j<=je; j++){
        for (i=is-nghost; i<is; i++){
          /* Get a pointer to the ConsS cell */
          pCons = &(pG->U[k][j][i]);

          *(pSnd++) = pCons->d;
          *(pSnd++) = pCons->M1;
          *(pSnd++) = pCons->M2;
          *(pSnd++) = pCons->M3;
#ifndef BAROTROPIC
          *(pSnd++) = pCons->E;
#endif /* BAROTROPIC */
#ifdef MHD
          *(pSnd++) = pCons->B1c;
          *(pSnd++) = pCons->B2c;
          *(pSnd++) = pCons->B3c;
          *(pSnd++) = pG->B1i[k][j][i];
          *(pSnd++) = pG->B2i[k][j][i];
          *(pSnd++) = pG->B3i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pSnd++) = pCons->s[n];
#endif
        }
      }
    }

/* send contents of buffer to the neighboring grid on R-x2 */
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                   shearing_sheet_ix1_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ku; k++){
      for (j=js-nghost; j<=js-1; j++){
        for (i=is-nghost; i<is; i++){
          /* Get a pointer to the ConsS cell */
          pCons = &(pG->U[k][j][i]);

          pCons->d  = *(pRcv++);
          pCons->M1 = *(pRcv++);
          pCons->M2 = *(pRcv++);
          pCons->M3 = *(pRcv++);
#ifndef BAROTROPIC
          pCons->E  = *(pRcv++);
#endif /* BAROTROPIC */
#ifdef MHD
          pCons->B1c = *(pRcv++);
          pCons->B2c = *(pRcv++);
          pCons->B3c = *(pRcv++);
          pG->B1i[k][j][i] = *(pRcv++);
          pG->B2i[k][j][i] = *(pRcv++);
          pG->B3i[k][j][i] = *(pRcv++);
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) pCons->s[n] = *(pRcv++);
#endif
        }
      }
    }

/* Post a non-blocking receive for the input data from the right grid */
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                    shearing_sheet_ix1_tag, pD->Comm_Domain, &rq);

    pSnd = send_buf;
    for (k=ks; k<=ku; k++){
      for (j=js; j<=js+nghost-1; j++){
        for (i=is-nghost; i<is; i++){
          /* Get a pointer to the ConsS cell */
          pCons = &(pG->U[k][j][i]);

          *(pSnd++) = pCons->d;
          *(pSnd++) = pCons->M1;
          *(pSnd++) = pCons->M2;
          *(pSnd++) = pCons->M3;
#ifndef BAROTROPIC
          *(pSnd++) = pCons->E;
#endif /* BAROTROPIC */
#ifdef MHD
          *(pSnd++) = pCons->B1c;
          *(pSnd++) = pCons->B2c;
          *(pSnd++) = pCons->B3c;
          *(pSnd++) = pG->B1i[k][j][i];
          *(pSnd++) = pG->B2i[k][j][i];
          *(pSnd++) = pG->B3i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pSnd++) = pCons->s[n];
#endif
        }
      }
    }

/* send contents of buffer to the neighboring grid on L-x2 */
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                   shearing_sheet_ix1_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ku; k++){
      for (j=je+1; j<=je+nghost; j++){
        for (i=is-nghost; i<is; i++){
          /* Get a pointer to the ConsS cell */
          pCons = &(pG->U[k][j][i]);

          pCons->d  = *(pRcv++);
          pCons->M1 = *(pRcv++);
          pCons->M2 = *(pRcv++);
          pCons->M3 = *(pRcv++);
#ifndef BAROTROPIC
          pCons->E  = *(pRcv++);
#endif /* BAROTROPIC */
#ifdef MHD
          pCons->B1c = *(pRcv++);
          pCons->B2c = *(pRcv++);
          pCons->B3c = *(pRcv++);
          pG->B1i[k][j][i] = *(pRcv++);
          pG->B2i[k][j][i] = *(pRcv++);
          pG->B3i[k][j][i] = *(pRcv++);
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) pCons->s[n] = *(pRcv++);
#endif
        }
      }
    }
#endif /* MPI_PARALLEL */

  } /* end of step 9 - periodic BC in Y with MPI */

/*--- Step 10 ------------------------------------------------------------------
 * Fix B2c at j=je,js-1, now that B2i[je+1] has been set properly  */

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for(i=is-nghost; i<is; i++){
      pG->U[k][je  ][i].B2c = 0.5*(pG->B2i[k][je+1][i]+pG->B2i[k][je][i]);
      pG->U[k][js-1][i].B2c = 0.5*(pG->B2i[k][js-1][i]+pG->B2i[k][js][i]);
    }
  }
#endif /* MHD */

  } /* end of if */

/*--- Step 11 ------------------------------------------------------------------
 * Shearing sheet BC in 2D xz.  Periodic BC already applied in x1 and x2 in
 * bvals_mhd (including for MPI parallel jobs).  Now just have to add offset
 * to azimuthal velocity when FARGO not defined */

#ifndef FARGO
  if (pG->Nx[2] == 1 && ShBoxCoord==xz) {

    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=1; i<=nghost; i++) {
#ifdef ADIABATIC
/* No change in the internal energy */
        pG->U[ks][j][is-i].E += (0.5/pG->U[ks][j][is-i].d)*
         (SQR((pG->U[ks][j][is-i].M3 + qomL*pG->U[ks][j][is-i].d))
        - SQR(pG->U[ks][j][is-i].M3));
#endif
        pG->U[ks][j][is-i].M3 += qomL*pG->U[ks][j][is-i].d;
      }
    }

  }
#endif /* FARGO */

  return;
}


/*----------------------------------------------------------------------------*/
/*! \fn void ShearingSheet_ox1(DomainS *pD)
 *  \brief 3D shearing-sheet BCs in x1.  
 *
 * It applies a remap
 * in Y after the ghost cells have been set by the usual periodic BCs in X and
 * Y implemented in bvals_mhd.c
 *
 * This is a public function which is called by bvals_mhd() inside a
 * SHEARING_BOX macro.							      */
/*----------------------------------------------------------------------------*/

void ShearingSheet_ox1(DomainS *pD)
{
  GridS *pG = pD->Grid;
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,ii,j,k,ku,n,joffset,jremap;
  Real xmin,xmax,Lx,Ly,qomL,yshear,deltay,epso;
#ifdef MPI_PARALLEL
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids;
  int ierr,sendto_id,getfrom_id;
  double *pSnd,*pRcv;
  Remap *pRemap;
  ConsS *pCons;
  MPI_Request rq;
#endif

/*--- Step 1. ------------------------------------------------------------------
 * Compute the distance the computational domain has sheared in y */

  xmin = pD->RootMinX[0];
  xmax = pD->RootMaxX[0];
  Lx = xmax - xmin;

  xmin = pD->RootMinX[1];
  xmax = pD->RootMaxX[1];
  Ly = xmax - xmin;

  qomL = qshear*Omega_0*Lx;
  yshear = qomL*pG->time;

/* Split this into integer and fractional peices of the Domain in y.  Ignore
 * the integer piece because the Grid is periodic in y */

  deltay = fmod(yshear, Ly);

/* further decompose the fractional peice into integer and fractional pieces of
 * a grid cell.  Note 0.0 <= epsi < 1.0.  If Domain has MPI decomposition in Y,
 * then it is possible that:  pD->Nx2 > joffset > pG->Nx2   */

  joffset = (int)(deltay/pG->dx2);
  epso = -(fmod(deltay,pG->dx2))/pG->dx2;

/*--- Step 2. ------------------------------------------------------------------
 * Copy data into GhstZns array.  Note i and j indices are switched.
 * Steps 2-10 are for 3D or 2D xy runs.  Step 11 handles 2D xz separately */

  if (pG->Nx[2] > 1 || ShBoxCoord==xy) {  /* this if ends at end of step 10 */

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for(k=ks; k<=ku; k++) {
    for(j=js-nghost; j<=je+nghost; j++){
      for(i=0; i<nghost; i++){
        ii = ie+1+i;
        GhstZns[k][i][j].U[0] = pG->U[k][j][ii].d;
        GhstZns[k][i][j].U[1] = pG->U[k][j][ii].M1;
        GhstZns[k][i][j].U[2] = pG->U[k][j][ii].M2;
#ifndef FARGO
        GhstZns[k][i][j].U[2] -= qomL*pG->U[k][j][ii].d;
#endif
        GhstZns[k][i][j].U[3] = pG->U[k][j][ii].M3;
#ifdef ADIABATIC
/* No change in the internal energy */
        GhstZns[k][i][j].U[4] = pG->U[k][j][ii].E + (0.5/GhstZns[k][i][j].U[0])*
          (SQR(GhstZns[k][i][j].U[2]) - SQR(pG->U[k][j][ii].M2));
#endif /* ADIABATIC */
#ifdef MHD
        GhstZns[k][i][j].U[NREMAP-4] = pG->U[k][j][ii].B1c;
        GhstZns[k][i][j].U[NREMAP-3] = pG->B1i[k][j][ii];
        GhstZns[k][i][j].U[NREMAP-2] = pG->B2i[k][j][ii];
        GhstZns[k][i][j].U[NREMAP-1] = pG->B3i[k][j][ii];
#endif /* MHD */
#if (NSCALARS > 0)
        for(n=0; n<NSCALARS; n++) GhstZns[k][i][j].s[n] = pG->U[k][j][ii].s[n];
#endif
      }
    }
  }

/*--- Step 3. ------------------------------------------------------------------
 * Copy GhstZns into buffer, at the same time apply a conservative remap of
 * solution over the fractional part of grid cell */

  for(k=ks; k<=ku; k++) {
    for(i=0; i<nghost; i++){

      for (n=0; n<(NREMAP); n++) {
        for (j=js-nghost; j<=je+nghost; j++) U[j] = GhstZns[k][i][j].U[n];
        RemapFlux(U,epso,js,je+1,Flx);
        for(j=js; j<=je; j++){
          GhstZnsBuf[k][i][j].U[n] = GhstZns[k][i][j].U[n] - (Flx[j+1]-Flx[j]);
        }
      }

#if (NSCALARS > 0)
      for (n=0; n<(NSCALARS); n++) {
        for (j=js-nghost; j<=je+nghost; j++) U[j] = GhstZns[k][i][j].s[n];
        RemapFlux(U,epso,js,je+1,Flx);
        for(j=js; j<=je; j++){
          GhstZnsBuf[k][i][j].s[n] = GhstZns[k][i][j].s[n] - (Flx[j+1]-Flx[j]);
        }
      }
#endif

    }
  }

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from buffer back into GhstZns.  */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ku; k++) {
      for(j=js; j<=je; j++){
        jremap = j + joffset;
        if (jremap > (int)je) jremap -= pG->Nx[1];

        for(i=0; i<nghost; i++){
          for (n=0; n<(NREMAP); n++) {
            GhstZns[k][i][j].U[n]  = GhstZnsBuf[k][i][jremap].U[n];
          }
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) { 
            GhstZns[k][i][j].s[n] = GhstZnsBuf[k][i][jremap].s[n];
          }
#endif
        }

      }
    }

#ifdef MPI_PARALLEL
  } else {

/*--- Step 5. ------------------------------------------------------------------
 * If Domain contains MPI decomposition in Y, then MPI calls are required for
 * the cyclic shift needed to apply shift over integer number of grid cells
 * during copy from buffer back into GhstZns.  */

    get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);

/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
    Ngrids = (int)(joffset/pG->Nx[1]);
    joverlap = joffset - Ngrids*pG->Nx[1];

/*--- Step 5a. -----------------------------------------------------------------
 * Find ids of processors that data in [js:js+(joverlap-1)] is sent to, and
 * data in [je-(overlap-1):je] is received from.  Only execute if joverlap>0  */
/* This can result in send/receive to self -- we rely on MPI to handle this
 * properly */

    if (joverlap != 0) {

      jproc = my_jproc - (Ngrids + 1);
      if (jproc < 0) jproc += pD->NGrid[1]; 
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc + (Ngrids + 1);
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1]; 
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

/*--- Step 5b. -----------------------------------------------------------------
 * Pack send buffer and send data in [js:js+(joverlap-1)] from GhstZnsBuf */

      cnt = nghost*joverlap*(ku-ks+1)*(NREMAP+NSCALARS);
/* Post a non-blocking receive for the input data */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_ox1_tag, pD->Comm_Domain, &rq);

      pSnd = send_buf;
      for (k=ks; k<=ku; k++) {
        for (j=js; j<=js+(joverlap-1); j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZnsBuf[k][i][j]);

            for (n=0; n<NREMAP; n++) *(pSnd++) = pRemap->U[n];
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) *(pSnd++) = pRemap->s[n];
#endif
          }
        }
      }
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     shearing_sheet_ox1_tag, pD->Comm_Domain);


/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [js:js+(joverlap-1)], and remap into cells in
 * [je-(joverlap-1):je] in GhstZns
 */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for (k=ks; k<=ku; k++) {
        for (j=je-(joverlap-1); j<=je; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZns[k][i][j]);
  
            for (n=0; n<NREMAP; n++) pRemap->U[n] = *(pRcv++);
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) pRemap->s[n] = *(pRcv++);
#endif
          }
        }
      }

    }

/*--- Step 5d. -----------------------------------------------------------------
 * If shear is less one full Grid, remap cells which remain on same processor
 * from GhstZnsBuf into GhstZns.  Cells in [js+joverlap:je] are shifted by
 * joverlap into [js:je-joverlap] */

    if (Ngrids == 0) {

      for(k=ks; k<=ku; k++) {
        for(j=js; j<=je-joverlap; j++){
          jremap = j+joverlap;
          for(i=0; i<nghost; i++){
            for (n=0; n<(NREMAP); n++) {
              GhstZns[k][i][j].U[n]  = GhstZnsBuf[k][i][jremap].U[n];
            }
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) { 
              GhstZns[k][i][j].s[n] = GhstZnsBuf[k][i][jremap].s[n];
            }
#endif
          }
        }
      }

/*--- Step 5e. -----------------------------------------------------------------
 * If shear is more than one Grid, pack and send data from [js+joverlap:je]
 * from GhstZnsBuf (this step replaces 5d) */

    } else {

/* index of sendto and getfrom processors in GData are -/+1 from Step 5a */

      jproc = my_jproc - Ngrids;
      if (jproc < 0) jproc += pD->NGrid[1];
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc + Ngrids;
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1];
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      cnt = nghost*(pG->Nx[1]-joverlap)*(ku-ks+1)*(NREMAP+NSCALARS);
/* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_ox1_tag, pD->Comm_Domain, &rq);

      pSnd = send_buf;
      for (k=ks; k<=ku; k++) {
        for (j=js+joverlap; j<=je; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZnsBuf[k][i][j]);
            for (n=0; n<NREMAP; n++) *(pSnd++) = pRemap->U[n];
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) *(pSnd++) = pRemap->s[n];
#endif
          }
        }
      }
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     shearing_sheet_ox1_tag, pD->Comm_Domain);

/* unpack data sent from [js+joverlap:je], and remap into cells in
 * [js:je-joverlap] in GhstZns */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for (k=ks; k<=ku; k++) {
        for (j=js; j<=je-joverlap; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZns[k][i][j]);
            for (n=0; n<NREMAP; n++) pRemap->U[n] = *(pRcv++);
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) pRemap->s[n] = *(pRcv++);
#endif
          }
        }
      }
    } /* end of step 5e - shear is more than one Grid */

#endif /* MPI_PARALLEL */
  } /* end of step 5 - MPI decomposition in Y */

/*--- Step 6. ------------------------------------------------------------------
 * Now copy remapped variables back into ghost cells, except B1i[ie+1] */

  for(k=ks; k<=ke; k++) {
    for(j=js; j<=je; j++){
      for(i=0; i<nghost; i++){
        pG->U[k][j][ie+1+i].d  = GhstZns[k][i][j].U[0];
        pG->U[k][j][ie+1+i].M1 = GhstZns[k][i][j].U[1];
        pG->U[k][j][ie+1+i].M2 = GhstZns[k][i][j].U[2];
        pG->U[k][j][ie+1+i].M3 = GhstZns[k][i][j].U[3];
#ifdef ADIABATIC
        pG->U[k][j][ie+1+i].E  = GhstZns[k][i][j].U[4];
#endif /* ADIABATIC */
#ifdef MHD
        pG->U[k][j][ie+1+i].B1c = GhstZns[k][i][j].U[NREMAP-4];
        if(i>0) pG->B1i[k][j][ie+1+i] = GhstZns[k][i][j].U[NREMAP-3];
        pG->B2i[k][j][ie+1+i] = GhstZns[k][i][j].U[NREMAP-2];
        pG->B3i[k][j][ie+1+i] = GhstZns[k][i][j].U[NREMAP-1];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) {
          pG->U[k][j][ie+1+i].s[n] = GhstZns[k][i][j].s[n];
        }
#endif
      }
    }
  }

/* Copy the face-centered B3 component of the field at k=ke+1 in 3D */
#ifdef MHD
  if (pG->Nx[2] > 1) {
    for(j=js; j<=je; j++){
      for(i=0; i<nghost; i++){
        pG->B3i[ke+1][j][ie+1+i] = GhstZns[ke+1][i][j].U[NREMAP-1];
      }
    }
  }
#endif /* MHD */

/*--- Step 7. ------------------------------------------------------------------
 * compute cell-centered B as average of remapped face centered B, except B1.
 * The value of B2c at j=je is incorrect since B2i[je+1] not yet set -- fix in
 * step 10 below */

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for(i=ie+1; i<=ie+nghost; i++){
        pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i]+pG->B2i[k][j+1][i]);
      }
    }
  }
  if (pG->Nx[2] > 1) {
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for(i=ie+1; i<=ie+nghost; i++){
          pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i]+pG->B3i[k+1][j][i]);
        }
      }
    }
  }
#endif /* MHD */

/*--- Step 8. ------------------------------------------------------------------
 * With no MPI decomposition in Y, apply periodic BCs in Y (similar to
 * periodic_ix2() and periodic_ox2() in bvals_mhd.c) */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ke; k++) {
      for(j=1; j<=nghost; j++){
        for(i=ie+1; i<=ie+nghost; i++){
          pG->U[k][js-j][i] = pG->U[k][je-(j-1)][i];
          pG->U[k][je+j][i] = pG->U[k][js+(j-1)][i];
#ifdef MHD
          pG->B1i[k][js-j][i] = pG->B1i[k][je-(j-1)][i];
          pG->B2i[k][js-j][i] = pG->B2i[k][je-(j-1)][i];
          pG->B3i[k][js-j][i] = pG->B3i[k][je-(j-1)][i];

          pG->B1i[k][je+j][i] = pG->B1i[k][js+(j-1)][i];
          pG->B2i[k][je+j][i] = pG->B2i[k][js+(j-1)][i];
          pG->B3i[k][je+j][i] = pG->B3i[k][js+(j-1)][i];
#endif /* MHD */
        }
      }
    }
#ifdef MHD
    for (j=1; j<=nghost; j++) {
      for (i=ie+1; i<=ie+nghost; i++) {
        pG->B3i[ke+1][js-j][i] = pG->B3i[ke+1][je-(j-1)][i];
        pG->B3i[ke+1][je+j][i] = pG->B3i[ke+1][js+(j-1)][i];
      }
    }
#endif /* MHD */

#ifdef MPI_PARALLEL
  } else {

/*--- Step 9. ------------------------------------------------------------------
 * With MPI decomposition in Y, use MPI calls to handle periodic BCs in Y (like
 * send_ox2/receive_ix1 and send_ix1/receive_ox2 pairs in bvals_mhd.c */


/* Post a non-blocking receive for the input data from the left grid */
    cnt = nghost*nghost*(ku-ks+1)*NVAR_SHARE;
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                    shearing_sheet_ox1_tag, pD->Comm_Domain, &rq);

    pSnd = send_buf;
    for (k=ks; k<=ku; k++){
      for (j=je-nghost+1; j<=je; j++){
        for (i=ie+1; i<=ie+nghost; i++){
          /* Get a pointer to the ConsS cell */
          pCons = &(pG->U[k][j][i]);

          *(pSnd++) = pCons->d;
          *(pSnd++) = pCons->M1;
          *(pSnd++) = pCons->M2;
          *(pSnd++) = pCons->M3;
#ifndef BAROTROPIC
          *(pSnd++) = pCons->E;
#endif /* BAROTROPIC */
#ifdef MHD
          *(pSnd++) = pCons->B1c;
          *(pSnd++) = pCons->B2c;
          *(pSnd++) = pCons->B3c;
          *(pSnd++) = pG->B1i[k][j][i];
          *(pSnd++) = pG->B2i[k][j][i];
          *(pSnd++) = pG->B3i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pSnd++) = pCons->s[n];
#endif
        }
      }
    }

/* send contents of buffer to the neighboring grid on R-x2 */
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                   shearing_sheet_ox1_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ku; k++){
      for (j=js-nghost; j<=js-1; j++){
        for (i=ie+1; i<=ie+nghost; i++){
          /* Get a pointer to the ConsS cell */
          pCons = &(pG->U[k][j][i]);

          pCons->d  = *(pRcv++);
          pCons->M1 = *(pRcv++);
          pCons->M2 = *(pRcv++);
          pCons->M3 = *(pRcv++);
#ifndef BAROTROPIC
          pCons->E  = *(pRcv++);
#endif /* BAROTROPIC */
#ifdef MHD
          pCons->B1c = *(pRcv++);
          pCons->B2c = *(pRcv++);
          pCons->B3c = *(pRcv++);
          pG->B1i[k][j][i] = *(pRcv++);
          pG->B2i[k][j][i] = *(pRcv++);
          pG->B3i[k][j][i] = *(pRcv++);
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) pCons->s[n] = *(pRcv++);
#endif
        }
      }
    }

/* Post a non-blocking receive for the input data from the right grid */
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                    shearing_sheet_ox1_tag, pD->Comm_Domain, &rq);

    pSnd = send_buf;
    for (k=ks; k<=ku; k++){
      for (j=js; j<=js+nghost-1; j++){
        for (i=ie+1; i<=ie+nghost; i++){
          /* Get a pointer to the ConsS cell */
          pCons = &(pG->U[k][j][i]);

          *(pSnd++) = pCons->d;
          *(pSnd++) = pCons->M1;
          *(pSnd++) = pCons->M2;
          *(pSnd++) = pCons->M3;
#ifndef BAROTROPIC
          *(pSnd++) = pCons->E;
#endif /* BAROTROPIC */
#ifdef MHD
          *(pSnd++) = pCons->B1c;
          *(pSnd++) = pCons->B2c;
          *(pSnd++) = pCons->B3c;
          *(pSnd++) = pG->B1i[k][j][i];
          *(pSnd++) = pG->B2i[k][j][i];
          *(pSnd++) = pG->B3i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pSnd++) = pCons->s[n];
#endif
        }
      }
    }

/* send contents of buffer to the neighboring grid on L-x2 */
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                   shearing_sheet_ox1_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ku; k++){
      for (j=je+1; j<=je+nghost; j++){
        for (i=ie+1; i<=ie+nghost; i++){
          /* Get a pointer to the ConsS cell */
          pCons = &(pG->U[k][j][i]);

          pCons->d  = *(pRcv++);
          pCons->M1 = *(pRcv++);
          pCons->M2 = *(pRcv++);
          pCons->M3 = *(pRcv++);
#ifndef BAROTROPIC
          pCons->E  = *(pRcv++);
#endif /* BAROTROPIC */
#ifdef MHD
          pCons->B1c = *(pRcv++);
          pCons->B2c = *(pRcv++);
          pCons->B3c = *(pRcv++);
          pG->B1i[k][j][i] = *(pRcv++);
          pG->B2i[k][j][i] = *(pRcv++);
          pG->B3i[k][j][i] = *(pRcv++);
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) pCons->s[n] = *(pRcv++);
#endif
        }
      }
    }
#endif /* MPI_PARALLEL */

  } /* end of step 9 - periodic BC in Y with MPI */

/*--- Step 10 ------------------------------------------------------------------
 * Fix B2c at j=je,js-1, now that B2i[je+1] has been set properly  */

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (i=ie+1; i<=ie+nghost; i++){
      pG->U[k][je  ][i].B2c = 0.5*(pG->B2i[k][je+1][i]+pG->B2i[k][je][i]);
      pG->U[k][js-1][i].B2c = 0.5*(pG->B2i[k][js-1][i]+pG->B2i[k][js][i]);
    }
  }
#endif /* MHD */

  } /* end of if */

/*--- Step 11 ------------------------------------------------------------------
 * Shearing sheet BC in 2D.  Periodic BC already applied in x1 and x2 in
 * bvals_mhd (including for MPI parallel jobs).  Now just have to add offset
 * to azimuthal velocity when FARGO not defined */

#ifndef FARGO
  if (pG->Nx[2] == 1 && ShBoxCoord==xz) {

    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=1; i<=nghost; i++) {
#ifdef ADIABATIC
/* No change in the internal energy */
        pG->U[ks][j][ie+i].E += (0.5/pG->U[ks][j][ie+i].d)*
          (SQR((pG->U[ks][j][ie+i].M3 - qomL*pG->U[ks][j][ie+i].d))
         - SQR(pG->U[ks][j][ie+i].M3));
#endif
        pG->U[ks][j][ie+i].M3 -= qomL*pG->U[ks][j][ie+i].d;
      }
    }

  }
#endif /* FARGO */

  return;
}


/*----------------------------------------------------------------------------*/
/*! \fn void RemapEy_ix1(DomainS *pD, Real ***emfy, Real **tEy)
 *  \brief Remaps Ey at [is] due to background shear, and then
 *   averages remapped and original field.  This guarantees the sums of Ey
 * along the x1 boundaries at [is] and [ie+1] are identical -- thus net Bz is
 * conserved
 *
 * This is a public function which is called by integrator (inside a
 * SHEARING_BOX macro).							      */
/*----------------------------------------------------------------------------*/

#ifdef MHD
void RemapEy_ix1(DomainS *pD, Real ***emfy, Real **tEy)
{
  GridS *pG = pD->Grid;
  int ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int j,k,joffset,jremap;
  Real xmin,xmax,Lx,Ly,qomL,yshear,deltay,epsi;
#ifdef MPI_PARALLEL
  int is = pG->is;
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids;
  int ierr,sendto_id,getfrom_id;
  double *pSnd,*pRcv;
  Real *pEy;
  MPI_Request rq;
#endif

/* Compute the distance the computational domain has sheared in y in integer
 * and fractional pieces of a cell.  Same code as in ShearingSheet_ix1()  */

  xmin = pD->RootMinX[0];
  xmax = pD->RootMaxX[0];
  Lx = xmax - xmin;

  xmin = pD->RootMinX[1];
  xmax = pD->RootMaxX[1];
  Ly = xmax - xmin;

  qomL = qshear*Omega_0*Lx;
  yshear = qomL*pG->time;
  deltay = fmod(yshear, Ly);
  joffset = (int)(deltay/pG->dx2);
  epsi = (fmod(deltay,pG->dx2))/pG->dx2;

/*--- Step 1. ------------------------------------------------------------------
 * Copy Ey from [ie+1] into temporary array, using periodic BC in x1.
 * Requires MPI calls if NGrid_x1 > 1   */

  if (pD->NGrid[0] == 1) {
    for(k=ks; k<=ke+1; k++) {
      for(j=js; j<=je; j++){
        tEy[k][j] = emfy[k][j][ie+1];
      }
    }

#ifdef MPI_PARALLEL
  } else {

/* MPI calls to swap data */

    cnt = pG->Nx[1]*(pG->Nx[2]+1); /* Post a non-blocking receive for the input data from remapEy_ox1 (listen L) */
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx1_id,
                    remapEy_tag, pD->Comm_Domain, &rq);

/* send Ey at [is] to ox1 (send L) -- this data is needed by remapEy_ox1 */

    pSnd = send_buf;
    for (k=ks; k<=ke+1; k++) {
      for (j=js; j<=je; j++) {
        pEy = &(emfy[k][j][is]);
        *(pSnd++) = *pEy;
      }
    }
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx1_id,
                   remapEy_tag, pD->Comm_Domain);

/* Listen for data from ox1 (listen L), unpack and set temporary array */

    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ke+1; k++) {
      for (j=js; j<=je; j++) {
          pEy = &(tEy[k][j]);
          *pEy = *(pRcv++);
      }
    }
#endif /* MPI_PARALLEL */
  }

/*--- Step 2. ------------------------------------------------------------------
 * Apply periodic BC in x2 to temporary array.  Requires MPI calls if 
 * NGrid_x2 > 1 */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ke+1; k++) {
      for(j=1; j<=nghost; j++){
        tEy[k][js-j] = tEy[k][je-(j-1)];
        tEy[k][je+j] = tEy[k][js+(j-1)];
      }
    }

#ifdef MPI_PARALLEL
  } else {

/* MPI calls to swap data */

    cnt = nghost*(pG->Nx[2] + 1);
/* Post a non-blocking receive for the input data from the left grid */
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                    remapEy_tag, pD->Comm_Domain, &rq);

    pSnd = send_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=je-nghost+1; j<=je; j++){
        pEy = &(tEy[k][j]);
        *(pSnd++) = *pEy;
      }
    }

/* send contents of buffer to the neighboring grid on R-x2 */
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                   remapEy_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=js-nghost; j<=js-1; j++){
        pEy = &(tEy[k][j]);
        *pEy = *(pRcv++);
      }
    }

/* Post a non-blocking receive for the input data from the right grid */
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                    remapEy_tag, pD->Comm_Domain, &rq);

    pSnd = send_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=js; j<=js+nghost-1; j++){
        pEy = &(tEy[k][j]);
        *(pSnd++) = *pEy;
      }
    }

/* send contents of buffer to the neighboring grid on L-x2 */
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                   remapEy_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=je+1; j<=je+nghost; j++){
        pEy = &(tEy[k][j]);
        *pEy = *(pRcv++);
      }
    }
#endif /* MPI_PARALLEL */
  }

/*--- Step 3. ------------------------------------------------------------------
 * Copy tEy into buffer, at the same time apply a conservative remap of
 * solution over the fractional part of grid cell */

  for(k=ks; k<=ke+1; k++) {
    RemapFlux(tEy[k],epsi,js,je+1,Flx);
    for(j=js; j<=je; j++){
      tEyBuf[k][j] = tEy[k][j] - (Flx[j+1] - Flx[j]);
    }
  }

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from buffer back into tEy.  */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ke+1; k++) {
      for(j=js; j<=je; j++){
        jremap = j - joffset;
        if (jremap < (int)js) jremap += pG->Nx[1];
        tEy[k][j]  = tEyBuf[k][jremap];
      }
    }

#ifdef MPI_PARALLEL
  } else {

/*--- Step 5. ------------------------------------------------------------------
 * If Domain contains MPI decomposition in Y, then MPI calls are required for
 * the cyclic shift needed to apply shift over integer number of grid cells
 * during copy from buffer back into tEy.  */

    get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);

/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
    Ngrids = (int)(joffset/pG->Nx[1]);
    joverlap = joffset - Ngrids*pG->Nx[1];

/*--- Step 5a. -----------------------------------------------------------------
 * Find ids of processors that data in [je-(joverlap-1):je] is sent to, and
 * data in [js:js+(joverlap-1)] is received from.  Only execute if joverlap>0 */
/* This can result in send/receive to self -- we rely on MPI to handle this
 * properly */

    if (joverlap != 0) {

      jproc = my_jproc + (Ngrids + 1);
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1];
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc - (Ngrids + 1);
      if (jproc < 0) jproc += pD->NGrid[1];
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

/*--- Step 5b. -----------------------------------------------------------------
 * Pack send buffer and send data in [je-(joverlap-1):je] from tEyBuf */

      cnt = joverlap*(pG->Nx[2]+1);
/* Post a non-blocking receive for the input data */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      remapEy_tag, pD->Comm_Domain, &rq);

      pSnd = send_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=je-(joverlap-1); j<=je; j++) {
          pEy = &(tEyBuf[k][j]);
          *(pSnd++) = *pEy;
        }
      }
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     remapEy_tag, pD->Comm_Domain);


/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [je-(joverlap-1):je], and remap into cells in
 * [js:js+(joverlap-1)] in tEy
 */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js; j<=js+(joverlap-1); j++) {
            pEy = &(tEy[k][j]);
            *pEy = *(pRcv++);
        }
      }

    }

/*--- Step 5d. -----------------------------------------------------------------
 * If shear is less one full Grid, remap cells which remain on same processor
 * from tEyBuf into tEy.  Cells in [js:je-joverlap] are shifted by
 * joverlap into [js+joverlap:je] */

    if (Ngrids == 0) {

      for(k=ks; k<=ke+1; k++) {
        for(j=js+joverlap; j<=je; j++){
          jremap = j-joverlap;
          tEy[k][j]  = tEyBuf[k][jremap];
        }
      }

/*--- Step 5e. -----------------------------------------------------------------
 * If shear is more than one Grid, pack and send data from [js:je-joverlap]
 * from tEyBuf (this step replaces 5d) */

    } else {

/* index of sendto and getfrom processors in GData are -/+1 from Step 5a */

      jproc = my_jproc + Ngrids;
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1];
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc - Ngrids;
      if (jproc < 0) jproc += pD->NGrid[1];
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      cnt = (pG->Nx[1]-joverlap)*(pG->Nx[2]+1);
/* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      remapEy_tag, pD->Comm_Domain, &rq);

      pSnd = send_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js; j<=je-joverlap; j++) {
          pEy = &(tEyBuf[k][j]);
          *(pSnd++) = *pEy;
        }
      }
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     remapEy_tag, pD->Comm_Domain);

/* unpack data sent from [js:je-overlap], and remap into cells in
 * [js+joverlap:je] in tEy */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js+joverlap; j<=je; j++) {
          pEy = &(tEy[k][j]);
          *pEy = *(pRcv++);
        }
      }
    } /* end of step 5e - shear is more than one Grid */

#endif /* MPI_PARALLEL */
  } /* end of step 5 - MPI decomposition in Y */

/*--- Step 6. ------------------------------------------------------------------
 * Now return remapped Ey */

  return;
}
#endif /* MHD */


/*----------------------------------------------------------------------------*/
/*! \fn void RemapEy_ox1(DomainS *pD, Real ***emfy, Real **tEy)
 *  \brief Remaps Ey at [ie+1] due to background shear, and then
 * averages remapped and original field.  This guarantees the sums of Ey
 * along the x1 boundaries at [is] and [ie+1] are identical -- thus net Bz is
 * conserved
 *
 * This is a public function which is called by integrator (inside a
 * SHEARING_BOX macro).							      */
/*----------------------------------------------------------------------------*/

#ifdef MHD
void RemapEy_ox1(DomainS *pD, Real ***emfy, Real **tEy)
{
  GridS *pG = pD->Grid;
  int is = pG->is;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int j,k,joffset,jremap;
  Real xmin,xmax,Lx,Ly,qomL,yshear,deltay,epso;
#ifdef MPI_PARALLEL
  int ie = pG->ie;
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids;
  int ierr,sendto_id,getfrom_id;
  double *pSnd,*pRcv;
  Real *pEy;
  MPI_Request rq;
#endif

/* Compute the distance the computational domain has sheared in y in integer
 * and fractional pieces of a cell.  Same code as in ShearingSheet_ox1()  */

  xmin = pD->RootMinX[0];
  xmax = pD->RootMaxX[0];
  Lx = xmax - xmin;

  xmin = pD->RootMinX[1];
  xmax = pD->RootMaxX[1];
  Ly = xmax - xmin;

  qomL = qshear*Omega_0*Lx;
  yshear = qomL*pG->time;
  deltay = fmod(yshear, Ly);
  joffset = (int)(deltay/pG->dx2);
  epso = -(fmod(deltay,pG->dx2))/pG->dx2;

/*--- Step 1. ------------------------------------------------------------------
 * Copy Ey from [is] into temporary array, using periodic BC in x1.
 * Requires MPI calls if NGrid_x1 > 1   */

  if (pD->NGrid[0] == 1) {
    for(k=ks; k<=ke+1; k++) {
      for(j=js; j<=je; j++){
        tEy[k][j] = emfy[k][j][is];
      }
    }

#ifdef MPI_PARALLEL
  } else {

/* MPI calls to swap data */

    cnt = pG->Nx[1]*(pG->Nx[2]+1);
/* Post a non-blocking receive for the input data from remapEy_ix1 (listen R) */
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx1_id,
                    remapEy_tag, pD->Comm_Domain, &rq);

/* send Ey at [ie+1] to ix1 (send R) -- this data is needed by remapEy_ix1 */

    pSnd = send_buf;
    for (k=ks; k<=ke+1; k++) {
      for (j=js; j<=je; j++) {
        pEy = &(emfy[k][j][ie+1]);
        *(pSnd++) = *pEy;
      }
    }
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx1_id,
                   remapEy_tag, pD->Comm_Domain);

/* Listen for data from ix1 (listen R), unpack and set temporary array */

    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ke+1; k++) {
      for (j=js; j<=je; j++) {
          pEy = &(tEy[k][j]);
          *pEy = *(pRcv++);
      }
    }
#endif /* MPI_PARALLEL */
  }

/*--- Step 2. ------------------------------------------------------------------
 * Apply periodic BC in x2 to temporary array.  Requires MPI calls if 
 * NGrid_x2 > 1 */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ke+1; k++) {
      for(j=1; j<=nghost; j++){
        tEy[k][js-j] = tEy[k][je-(j-1)];
        tEy[k][je+j] = tEy[k][js+(j-1)];
      }
    }

#ifdef MPI_PARALLEL
  } else {

/* MPI calls to swap data */

    cnt = nghost*(pG->Nx[2] + 1);
/* Post a non-blocking receive for the input data from the left grid */
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                    remapEy_tag, pD->Comm_Domain, &rq);

    pSnd = send_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=je-nghost+1; j<=je; j++){
        pEy = &(tEy[k][j]);
        *(pSnd++) = *pEy;
      }
    }

/* send contents of buffer to the neighboring grid on R-x2 */
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                   remapEy_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=js-nghost; j<=js-1; j++){
        pEy = &(tEy[k][j]);
        *pEy = *(pRcv++);
      }
    }

/* Post a non-blocking receive for the input data from the right grid */
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                    remapEy_tag, pD->Comm_Domain, &rq);

    pSnd = send_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=js; j<=js+nghost-1; j++){
        pEy = &(tEy[k][j]);
        *(pSnd++) = *pEy;
      }
    }

/* send contents of buffer to the neighboring grid on L-x2 */
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                   remapEy_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=je+1; j<=je+nghost; j++){
        pEy = &(tEy[k][j]);
        *pEy = *(pRcv++);
      }
    }
#endif /* MPI_PARALLEL */
  }

/*--- Step 3. ------------------------------------------------------------------
 * Copy tEy into buffer, at the same time apply a conservative remap of
 * solution over the fractional part of grid cell */

  for(k=ks; k<=ke+1; k++) {
    RemapFlux(tEy[k],epso,js,je+1,Flx);
    for(j=js; j<=je; j++){
      tEyBuf[k][j] = tEy[k][j] - (Flx[j+1] - Flx[j]);
    }
  }

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from buffer back into tEy.  */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ke+1; k++) {
      for(j=js; j<=je; j++){
        jremap = j + joffset;
        if (jremap > (int)je) jremap -= pG->Nx[1];
        tEy[k][j]  = tEyBuf[k][jremap];
      }
    }

#ifdef MPI_PARALLEL
  } else {

/*--- Step 5. ------------------------------------------------------------------
 * If Domain contains MPI decomposition in Y, then MPI calls are required for
 * the cyclic shift needed to apply shift over integer number of grid cells
 * during copy from buffer back into tEy.  */

    get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);

/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
    Ngrids = (int)(joffset/pG->Nx[1]);
    joverlap = joffset - Ngrids*pG->Nx[1];

/*--- Step 5a. -----------------------------------------------------------------
 * Find ids of processors that data in [js:js+(joverlap-1)] is sent to, and
 * data in [je-(joverlap-1):je] is received from.  Only execute if joverlap>0 */
/* This can result in send/receive to self -- we rely on MPI to handle this
 * properly */

    if (joverlap != 0) {

      jproc = my_jproc - (Ngrids + 1);
      if (jproc < 0) jproc += pD->NGrid[1];
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc + (Ngrids + 1);
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1];
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

/*--- Step 5b. -----------------------------------------------------------------
 * Pack send buffer and send data in [js:js+(joverlap-1)] from tEyBuf */

      cnt = joverlap*(pG->Nx[2]+1);
/* Post a non-blocking receive for the input data */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      remapEy_tag, pD->Comm_Domain, &rq);

      pSnd = send_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js; j<=js+(joverlap-1); j++) {
          pEy = &(tEyBuf[k][j]);
          *(pSnd++) = *pEy;
        }
      }
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     remapEy_tag, pD->Comm_Domain);

/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [js:js+(joverlap-1)], and remap into cells in
 * [je-(joverlap-1):je] in tEy
 */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=je-(joverlap-1); j<=je; j++) {
            pEy = &(tEy[k][j]);
            *pEy = *(pRcv++);
        }
      }

    }

/*--- Step 5d. -----------------------------------------------------------------
 * If shear is less one full Grid, remap cells which remain on same processor
 * from tEyBuf into tEy.  Cells in [js+joverlap:je] are shifted by
 * joverlap into [js:je-overlap] */

    if (Ngrids == 0) {

      for(k=ks; k<=ke+1; k++) {
        for(j=js; j<=je-joverlap; j++){
          jremap = j+joverlap;
          tEy[k][j]  = tEyBuf[k][jremap];
        }
      }

/*--- Step 5e. -----------------------------------------------------------------
 * If shear is more than one Grid, pack and send data from [js+overlap:je]
 * from tEyBuf (this step replaces 5d) */

    } else {

/* index of sendto and getfrom processors in GData are -/+1 from Step 5a */

      jproc = my_jproc - Ngrids;
      if (jproc < 0) jproc += pD->NGrid[1];
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc + Ngrids;
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1];
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      cnt = (pG->Nx[1]-joverlap)*(pG->Nx[2]+1);
/* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      remapEy_tag, pD->Comm_Domain, &rq);

      pSnd = send_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js+joverlap; j<=je; j++) {
          pEy = &(tEyBuf[k][j]);
          *(pSnd++) = *pEy;
        }
      }
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     remapEy_tag, pD->Comm_Domain);

/* unpack data sent from [js+overlap:je], and remap into cells in
 * [js:je-joverlap] in tEy */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js; j<=je-joverlap; j++) {
          pEy = &(tEy[k][j]);
          *pEy = *(pRcv++);
        }
      }
    } /* end of step 5e - shear is more than one Grid */

#endif /* MPI_PARALLEL */
  } /* end of step 5 - MPI decomposition in Y */

/*--- Step 6. ------------------------------------------------------------------
 * Now return remapped Ey */

  return;
}
#endif /* MHD */

#endif /* Shearing Box */

/*----------------------------------------------------------------------------*/
/*! \fn void Fargo(DomainS *pD)
 *  \brief Implements FARGO algorithm.  Only works in 3D or 2D xy.  Called in
 * the main loop after the integrator (and before bvals_mhd).
 *  Written in Munich 24-27.4.2008					      */
/*----------------------------------------------------------------------------*/

#ifdef FARGO
void Fargo(DomainS *pD)
{
  GridS *pG = pD->Grid;
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int jfs = nfghost, jfe = pG->Nx[1] + nfghost - 1;
  int i,j,k,ku,jj,n,joffset;
  Real x1,x2,x3,yshear,eps;
#if defined(ADIABATIC) && defined(SHEARING_BOX)
  Real qom_dt = qshear*Omega_0*pG->dt;
#endif
#ifdef CYLINDRICAL
  const Real *r=pG->r, *ri=pG->ri;
#endif
#ifdef MPI_PARALLEL
  int ierr,cnt;
  double *pSnd,*pRcv;
  FConsS *pFCons;
  MPI_Request rq;
#endif

/*--- Step 1. ------------------------------------------------------------------
 * Copy data into FargoVars array.  Note i and j indices are switched.
 * Since the FargoVars array has extra ghost zones in y, it must be indexed
 * over [jfs:jfe] instead of [js:je]  */
  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for(k=ks; k<=ku; k++) {
    for(j=jfs; j<=jfe+1; j++){
      for(i=is; i<=ie+1; i++){
        jj = j-(jfs-js);
        FargoVars[k][i][j].U[0] = pG->U[k][jj][i].d;
        FargoVars[k][i][j].U[1] = pG->U[k][jj][i].M1;
        FargoVars[k][i][j].U[2] = pG->U[k][jj][i].M2;
        FargoVars[k][i][j].U[3] = pG->U[k][jj][i].M3;
#if defined(ADIABATIC) && defined(SHEARING_BOX)
#ifdef MHD
/* Add energy equation source term in MHD */
        pG->U[k][jj][i].E -= qom_dt*pG->U[k][jj][i].B1c*
         (pG->U[k][jj][i].B2c - (qom_dt/2.)*pG->U[k][jj][i].B1c);
#endif /* MHD */
        FargoVars[k][i][j].U[4] = pG->U[k][jj][i].E;
#endif /* ADIABATIC */
/* Only store Bz and Bx in that order.  This is to match order in FargoFlx:
 *  FargoFlx.U[NFARGO-2] = emfx = -Vy*Bz
 *  FargoFlx.U[NFARGO-1] = emfy = Vy*Bx  */
#ifdef MHD
        FargoVars[k][i][j].U[NFARGO-2] = pG->B3i[k][jj][i];
        FargoVars[k][i][j].U[NFARGO-1] = pG->B1i[k][jj][i];
#endif /* MHD */
#if (NSCALARS > 0)
        for(n=0;n<NSCALARS;n++) FargoVars[k][i][j].s[n] = pG->U[k][jj][i].s[n];
#endif
      }
    }
  }
/*--- Step 2. ------------------------------------------------------------------
 * With no MPI decomposition in Y, apply periodic BCs in Y to FargoVars array
 * (method is similar to periodic_ix2() and periodic_ox2() in bvals_mhd).
 * Note there are extra ghost zones in Y in the FargoVars array  */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ku; k++) {
      for(i=is; i<=ie+1; i++){
        for(j=1; j<=nfghost; j++){
          FargoVars[k][i][jfs-j] = FargoVars[k][i][jfe-(j-1)];
          FargoVars[k][i][jfe+j] = FargoVars[k][i][jfs+(j-1)];
        }
      }
    }

#ifdef MPI_PARALLEL
  } else {

/*--- Step 3. ------------------------------------------------------------------
 * With MPI decomposition in Y, use MPI calls to handle periodic BCs in Y (like
 * send_ox2/receive_ix1 and send_ix1/receive_ox2 pairs in bvals_mhd.c */

/* Post a non-blocking receive for the input data from the left grid */
    cnt = (pG->Nx[0]+1)*nfghost*(ku-ks+1)*(NFARGO + NSCALARS);
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                    fargo_tag, pD->Comm_Domain, &rq);

    pSnd = send_buf;
    for (k=ks; k<=ku; k++){
      for (i=is; i<=ie+1; i++){
        for (j=jfe-nfghost+1; j<=jfe; j++){
          /* Get a pointer to the FargoVars cell */
          pFCons = &(FargoVars[k][i][j]);
          for (n=0; n<NFARGO; n++) *(pSnd++) = pFCons->U[n];
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pSnd++) = pFCons->s[n];
#endif
        }
      }
    }

/* send contents of buffer to the neighboring grid on R-x2 */
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                   fargo_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ku; k++){
      for (i=is; i<=ie+1; i++){
        for (j=jfs-nfghost; j<=jfs-1; j++){
          /* Get a pointer to the FargoVars cell */
          pFCons = &(FargoVars[k][i][j]);
          for (n=0; n<NFARGO; n++) pFCons->U[n] = *(pRcv++);
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) pFCons->s[n] = *(pRcv++);
#endif
        }
      }
    }

/* Post a non-blocking receive for the input data from the right grid */
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                    fargo_tag, pD->Comm_Domain, &rq);

    pSnd = send_buf;
    for (k=ks; k<=ku; k++){
      for (i=is; i<=ie+1; i++){
        for (j=jfs; j<=jfs+nfghost-1; j++){
          /* Get a pointer to the FargoVars cell */
          pFCons = &(FargoVars[k][i][j]);
          for (n=0; n<NFARGO; n++) *(pSnd++) = pFCons->U[n];
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pSnd++) = pFCons->s[n];
#endif
        }
      }
    }

/* send contents of buffer to the neighboring grid on L-x2 */
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                   fargo_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ku; k++){
      for (i=is; i<=ie+1; i++){
        for (j=jfe+1; j<=jfe+nfghost; j++){
          /* Get a pointer to the FargoVars cell */
          pFCons = &(FargoVars[k][i][j]);
          for (n=0; n<NFARGO; n++) pFCons->U[n] = *(pRcv++);
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) pFCons->s[n] = *(pRcv++);
#endif
        }
      }
    }
#endif /* MPI_PARALLEL */
  }

/*--- Step 4. ------------------------------------------------------------------
 * Compute fluxes, including both (1) the fractional part of grid cell, and
 * (2) the sum over integer number of cells  */
  for(k=ks; k<=ku; k++) {
    for(i=is; i<=ie+1; i++){

/* Compute integer and fractional peices of a cell covered by shear */
      cc_pos(pG, i, js, ks, &x1,&x2,&x3);
#ifdef SHEARING_BOX
      yshear = -qshear*Omega_0*x1*pG->dt;
#endif
#ifdef CYLINDRICAL
			yshear = ((*OrbitalProfile)(x1))*pG->dt;
#endif
      joffset = (int)(yshear/pG->dx2);
      if (abs(joffset) > (jfs-js))
        ath_error("[bvals_shear]: FARGO offset exceeded # of gh zns\n");
      eps = (fmod(yshear,pG->dx2))/pG->dx2;

/* Compute fluxes of hydro variables  */
#ifdef MHD
      for (n=0; n<(NFARGO-2); n++) {
#else
      for (n=0; n<(NFARGO); n++) {
#endif
        for (j=jfs-nfghost; j<=jfe+nfghost; j++) U[j] = FargoVars[k][i][j].U[n];
        RemapFlux(U,eps,(jfs-joffset),(jfe+1-joffset),Flx);

        for(j=jfs; j<=jfe+1; j++){
          FargoFlx[k][i][j].U[n] = Flx[j-joffset];

/* Sum the flux from integer offset: +/- for +/- joffset */
          for (jj=1; jj<=joffset; jj++) {
            FargoFlx[k][i][j].U[n] += FargoVars[k][i][j-jj].U[n];
          }
          for (jj=(joffset+1); jj<=0; jj++) {
            FargoFlx[k][i][j].U[n] -= FargoVars[k][i][j-jj].U[n];
          }
        }
      }

/* Compute fluxes of passive scalars */
#if (NSCALARS > 0)
      for (n=0; n<(NSCALARS); n++) {

        for (j=jfs-nfghost; j<=jfe+nfghost; j++) U[j] = FargoVars[k][i][j].s[n];
        RemapFlux(U,eps,(jfs-joffset),(jfe+1-joffset),Flx);

        for(j=jfs; j<=jfe+1; j++){
          FargoFlx[k][i][j].s[n] = Flx[j-joffset];
          for (jj=1; jj<=joffset; jj++) {
            FargoFlx[k][i][j].s[n] += FargoVars[k][i][j-jj].s[n];
          }
          for (jj=(joffset+1); jj<=0; jj++) {
            FargoFlx[k][i][j].s[n] -= FargoVars[k][i][j-jj].s[n];
          }
        }

      }
#endif

#ifdef MHD
/* Compute emfx = -VyBz, which is at cell-center in x1-direction */

      for (j=jfs-nfghost; j<=jfe+nfghost; j++) {
        U[j] = FargoVars[k][i][j].U[NFARGO-2];
      }
      RemapFlux(U,eps,(jfs-joffset),(jfe+1-joffset),Flx);

      for(j=jfs; j<=jfe+1; j++){
        FargoFlx[k][i][j].U[NFARGO-2] = -Flx[j-joffset];
        for (jj=1; jj<=joffset; jj++) {
          FargoFlx[k][i][j].U[NFARGO-2] -= FargoVars[k][i][j-jj].U[NFARGO-2];
        }
        for (jj=(joffset+1); jj<=0; jj++) {
          FargoFlx[k][i][j].U[NFARGO-2] += FargoVars[k][i][j-jj].U[NFARGO-2];
        }
      }

/* Compute emfz =  VyBx, which is at cell-face in x1-direction  */
#ifdef SHEARING_BOX
      yshear = -qshear*Omega_0*(x1 - 0.5*pG->dx1)*pG->dt;
#endif
#ifdef CYLINDRICAL
			yshear = ((*OrbitalProfile)(ri[i]))*pG->dt;
#endif
      joffset = (int)(yshear/pG->dx2);
      if (abs(joffset) > (jfs-js))
        ath_error("[bvals_shear]: FARGO offset exceeded # of gh zns\n");
      eps = (fmod(yshear,pG->dx2))/pG->dx2;

      for (j=jfs-nfghost; j<=jfe+nfghost; j++) {
        U[j] = FargoVars[k][i][j].U[NFARGO-1];
      }
      RemapFlux(U,eps,(jfs-joffset),(jfe+1-joffset),Flx);

      for(j=jfs; j<=jfe+1; j++){
        FargoFlx[k][i][j].U[NFARGO-1] = Flx[j-joffset];
        for (jj=1; jj<=joffset; jj++) {
          FargoFlx[k][i][j].U[NFARGO-1] += FargoVars[k][i][j-jj].U[NFARGO-1];
        }
        for (jj=(joffset+1); jj<=0; jj++) {
          FargoFlx[k][i][j].U[NFARGO-1] -= FargoVars[k][i][j-jj].U[NFARGO-1];
        }
      }
#endif

    }
  }
/*--- Step 5. ------------------------------------------------------------------
 * Update cell centered variables with flux gradient.  Note i/j are swapped */

  for(k=ks; k<=ke; k++) {
    for(j=js; j<=je; j++){
      jj = j-js+jfs;
      for(i=is; i<=ie; i++){
        pG->U[k][j][i].d  -=(FargoFlx[k][i][jj+1].U[0]-FargoFlx[k][i][jj].U[0]);
        pG->U[k][j][i].M1 -=(FargoFlx[k][i][jj+1].U[1]-FargoFlx[k][i][jj].U[1]);
        pG->U[k][j][i].M2 -=(FargoFlx[k][i][jj+1].U[2]-FargoFlx[k][i][jj].U[2]);
        pG->U[k][j][i].M3 -=(FargoFlx[k][i][jj+1].U[3]-FargoFlx[k][i][jj].U[3]);
#ifdef ADIABATIC
        pG->U[k][j][i].E  -=(FargoFlx[k][i][jj+1].U[4]-FargoFlx[k][i][jj].U[4]);
#endif /* ADIABATIC */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) {
         pG->U[k][j][i].s[n]-=FargoFlx[k][i][jj+1].s[n]-FargoFlx[k][i][jj].s[n];
        }
#endif
      }
    }
  }

/*--- Step 6. ------------------------------------------------------------------
 * Update face centered field using CT.  Note i/j are swapped.
 *  FargoFlx.U[NFARGO-2] = emfx
 *  FargoFlx.U[NFARGO-1] = emfz  */

#ifdef MHD
  if (pG->Nx[2]==1) ath_error("[Fargo] only works in 3D with MHD\n");
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      jj = j-js+jfs;
      for (i=is; i<=ie; i++) {
        pG->B1i[k][j][i] -= 
          (FargoFlx[k][i][jj+1].U[NFARGO-1] - FargoFlx[k][i][jj].U[NFARGO-1]);
#ifdef SHEARING_BOX
        pG->B2i[k][j][i] += (pG->dx2/pG->dx1)* 
          (FargoFlx[k][i+1][jj].U[NFARGO-1] - FargoFlx[k][i][jj].U[NFARGO-1]);
        pG->B2i[k][j][i] -= (pG->dx2/pG->dx3)* 
          (FargoFlx[k+1][i][jj].U[NFARGO-2] - FargoFlx[k][i][jj].U[NFARGO-2]);
#endif 
#ifdef CYLINDRICAL
        pG->B2i[k][j][i] += (pG->dx2/pG->dx1)*
                  (ri[i+1]*FargoFlx[k][i+1][jj].U[NFARGO-1] - ri[i]*FargoFlx[k][i][jj].U[NFARGO-1]);
        pG->B2i[k][j][i] -= (r[i]*pG->dx2/pG->dx3)*
          (FargoFlx[k+1][i][jj].U[NFARGO-2] - FargoFlx[k][i][jj].U[NFARGO-2]);
#endif
        pG->B3i[k][j][i] += 
          (FargoFlx[k][i][jj+1].U[NFARGO-2] - FargoFlx[k][i][jj].U[NFARGO-2]);
      }
      pG->B1i[k][j][ie+1] -= 
        (FargoFlx[k][ie+1][jj+1].U[NFARGO-1]-FargoFlx[k][ie+1][jj].U[NFARGO-1]);
    }
    for (i=is; i<=ie; i++) {
#ifdef SHEARING_BOX
      pG->B2i[k][je+1][i] += (pG->dx2/pG->dx1)* 
        (FargoFlx[k][i+1][jfe+1].U[NFARGO-1]-FargoFlx[k][i][jfe+1].U[NFARGO-1]);
      pG->B2i[k][je+1][i] -= (pG->dx2/pG->dx3)* 
        (FargoFlx[k+1][i][jfe+1].U[NFARGO-2]-FargoFlx[k][i][jfe+1].U[NFARGO-2]);
#endif
#ifdef CYLINDRICAL
      pG->B2i[k][je+1][i] += (pG->dx2/pG->dx1)*
        (ri[i+1]*FargoFlx[k][i+1][jfe+1].U[NFARGO-1]-ri[i]*FargoFlx[k][i][jfe+1].U[NFARGO-1]);
      pG->B2i[k][je+1][i] -= (r[i]*pG->dx2/pG->dx3)*
        (FargoFlx[k+1][i][jfe+1].U[NFARGO-2]-FargoFlx[k][i][jfe+1].U[NFARGO-2]);
#endif
    }
  }
  for (j=js; j<=je; j++) {
    jj = j-js+jfs;
    for (i=is; i<=ie; i++) {
      pG->B3i[ke+1][j][i] += 
        (FargoFlx[ke+1][i][jj+1].U[NFARGO-2]-FargoFlx[ke+1][i][jj].U[NFARGO-2]);
    }
  }
#endif /* MHD */

/*--- Step 7. ------------------------------------------------------------------
 * compute cell-centered B  */

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for(i=is; i<=ie; i++){
#ifdef SHEARING_BOX
        pG->U[k][j][i].B1c = 0.5*(pG->B1i[k][j][i]+pG->B1i[k][j][i+1]);
#endif
#ifdef CYLINDRICAL
        pG->U[k][j][i].B1c = 0.5*(1.0/r[i])*(ri[i]*pG->B1i[k][j][i] + ri[i+1]*pG->B1i[k][j][i+1]);
#endif
        pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i]+pG->B2i[k][j+1][i]);
        pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i]+pG->B3i[k+1][j][i]);
      }
    }
  }
#endif /* MHD */

  return;
}
#endif /* FARGO */

#if defined(SHEARING_BOX) || (defined(CYLINDRICAL) && defined(FARGO))
/*----------------------------------------------------------------------------*/
/*! \fn void bvals_shear_init(MeshS *pM)
 *  \brief Allocates memory for temporary arrays/buffers
 */

void bvals_shear_init(MeshS *pM)
{
  GridS *pG;
  int nl,nd,nx1,nx2,nx3,max1=0,max2=0,max3=0;
#ifdef MPI_PARALLEL
  int size1=0,size2=0,size;
#endif
#ifdef FARGO
  Real xmin,xmax,x2,x3;
#endif
#ifdef CYLINDRICAL
	Real MachKep;
#endif

/* Loop over all Grids on this processor to find maximum size of arrays */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) { /* there is a Grid on this proc */
        pG=pM->Domain[nl][nd].Grid;          /* set pointer to Grid */

        nx1 = pG->Nx[0] + 2*nghost;
        nx2 = pG->Nx[1] + 2*nghost;
        nx3 = pG->Nx[2] + 2*nghost;
        max1 = MAX(max1,nx1);
        max2 = MAX(max2,nx2);
        max3 = MAX(max3,nx3);
      }
    }
  }

/* Allocate memory for temporary arrays and vectors */

  if((GhstZns=(Remap***)calloc_3d_array(max3,nghost,max2,sizeof(Remap)))==NULL)
    ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");

  if((GhstZnsBuf=(Remap***)calloc_3d_array(max3,nghost,max2,sizeof(Remap))) ==
    NULL) ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");

#ifdef MHD
  if ((tEyBuf=(Real**)calloc_2d_array(max3,max2,sizeof(Real))) == NULL)
    ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");
#endif /* MHD */

/* estimate extra ghost cells needed by FARGO, allocate arrays accordingly */
#ifdef FARGO
  xmin = pM->RootMinX[0];
  xmax = pM->RootMaxX[0];
#ifdef SHEARING_BOX
  nfghost = nghost + ((int)(1.5*CourNo*MAX(fabs(xmin),fabs(xmax))) + 1);
#endif
#ifdef CYLINDRICAL
  MachKep = MAX( xmin*(*OrbitalProfile)(xmin), xmax*(*OrbitalProfile)(xmax) )/Iso_csound;
  nfghost = nghost + 1 + ((int)(CourNo*MachKep));
#endif
  max2 = max2 + 2*nfghost;

  if((FargoVars=(FConsS***)calloc_3d_array(max3,max1,max2,sizeof(FConsS)))
    ==NULL) ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");

  if((FargoFlx=(FConsS***)calloc_3d_array(max3,max1,max2,sizeof(FConsS)))==NULL)
    ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");
#endif

  if((U = (Real*)malloc(max2*sizeof(Real))) == NULL)
    ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");

  if((Flx = (Real*)malloc(max2*sizeof(Real))) == NULL)
    ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");

#if defined(THIRD_ORDER_CHAR) || defined(THIRD_ORDER_PRIM)
  if ((Uhalf = (Real*)malloc(max2*sizeof(Real))) == NULL)
    ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");
#endif

/* allocate memory for send/receive buffers in MPI parallel calculations */

#ifdef MPI_PARALLEL
  size1 = nghost*pG->Nx[1]*(pG->Nx[2]+1)*(NREMAP+NSCALARS);
#ifdef FARGO
  size2 = nfghost*(pG->Nx[0]+1)*(pG->Nx[2]+1)*(NFARGO+NSCALARS);
#endif
  size = MAX(size1,size2);

  if((send_buf = (double*)malloc(size*sizeof(double))) == NULL)
    ath_error("[bvals_shear_init]: Failed to allocate send buffer\n");

  if((recv_buf = (double*)malloc(size*sizeof(double))) == NULL)
    ath_error("[bvals_shear_init]: Failed to allocate receive buffer\n");
#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void bvals_shear_destruct(void)
 *  \brief Free temporary arrays
 */

void bvals_shear_destruct(void)
{
  if (GhstZns    != NULL) free_3d_array(GhstZns);
  if (GhstZnsBuf != NULL) free_3d_array(GhstZnsBuf);
#ifdef MHD
  if (tEyBuf != NULL) free_2d_array(tEyBuf);
#endif
#if defined(THIRD_ORDER_CHAR) || defined(THIRD_ORDER_PRIM)
  if (Uhalf != NULL) free(Uhalf);
#endif
#ifdef FARGO
  if (FargoVars != NULL) free_3d_array(FargoVars);
  if (FargoFlx  != NULL) free_3d_array(FargoFlx);
#endif
  if (U   != NULL) free(U);
  if (Flx != NULL) free(Flx);
#ifdef MPI_PARALLEL
  if (send_buf != NULL) free(send_buf);
  if (recv_buf != NULL) free(recv_buf);
#endif /* MPI_PARALLEL */

  return;
}


/*=========================== PRIVATE FUNCTIONS ==============================*/

/*------------------------------------------------------------------------------
 * RemapFlux: computes "fluxes" of conserved variables for conservative remap
 * Input Arguments:
 *   U = 1D vector of conserved variable at cell centers along 1-D slice
 *   eps = fraction of a cell to be remapped
 * Output Arguments:
 *   Flux = fluxes of conserved variable at interfaces over [jinner:jouter]
 */

#if defined(SECOND_ORDER_CHAR) || defined (SECOND_ORDER_PRIM)
/*----------------------------------------------------------------------------*/
/*! \fn void RemapFlux(const Real *U, const Real eps,
 *             const int jinner, const int jouter, Real *Flux)
 *  \brief Second order reconstruction for conservative remap.
 *   using piecewise linear reconstruction and min/mod limiters
 */

void RemapFlux(const Real *U, const Real eps,
               const int jinner, const int jouter, Real *Flux)
{
  int j,jl,ju;
  Real dUc,dUl,dUr,dUm,lim_slope;

/* jinner,jouter are index range over which flux must be returned.  Set loop
 * limits depending on direction of upwind differences  */

  if (eps > 0.0) { /* eps always > 0 for inner i boundary */
    jl = jinner-1;
    ju = jouter-1;
  } else {         /* eps always < 0 for outer i boundary */
    jl = jinner;
    ju = jouter;
  }

  for (j=jl; j<=ju; j++) {
      dUc = U[j+1] - U[j-1];
      dUl = U[j  ] - U[j-1];
      dUr = U[j+1] - U[j  ];

      dUm = 0.0;
      if (dUl*dUr > 0.0) {
        lim_slope = MIN(fabs(dUl),fabs(dUr));
        dUm = SIGN(dUc)*MIN(0.5*fabs(dUc),2.0*lim_slope);
      }
 
    if (eps > 0.0) { /* eps always > 0 for inner i boundary */
      Flux[j+1] = eps*(U[j] + 0.5*(1.0 - eps)*dUm);
    } else {         /* eps always < 0 for outer i boundary */
      Flux[j  ] = eps*(U[j] - 0.5*(1.0 + eps)*dUm);
    }
  }

  return;
}

#endif /* SECOND_ORDER */

#if defined(THIRD_ORDER_CHAR) || defined(THIRD_ORDER_PRIM)
/*----------------------------------------------------------------------------*/
/*! \fn void RemapFlux(const Real *U, const Real eps,
 *             const int jinner, const int jouter, Real *Flux)
 *  \brief third order reconstruction for conservative remap 
 *   using Colella & Sekora extremum preserving algorithm (PPME)
 */

void RemapFlux(const Real *U, const Real eps,
               const int jinner, const int jouter, Real *Flux)
{
  int j,jl,ju;
  Real d2Uc,d2Ul,d2Ur,d2U,d2Ulim,lim_slope,Ulv,Urv,dU,U6,qa,qb,qc,qx;

/* jinner,jouter are index range over which flux must be returned.  Set loop
 * limits depending on direction of upwind differences  */

  if (eps > 0.0) { /* eps always > 0 for inner i boundary */
    jl = jinner-1;
    ju = jouter-1;
  } else {         /* eps always < 0 for outer i boundary */
    jl = jinner;
    ju = jouter;
  }

  for (j=jl; j<=ju+1; j++) {
    Uhalf[j]=(7.0*(U[j-1]+U[j]) - (U[j-2]+U[j+1]))/12.0;
    d2Uc = 3.0*(U[j-1] - 2.0*Uhalf[j] + U[j]);
    d2Ul = (U[j-2] - 2.0*U[j-1] + U[j  ]);
    d2Ur = (U[j-1] - 2.0*U[j  ] + U[j+1]);
    d2Ulim = 0.0;
    lim_slope = MIN(fabs(d2Ul),fabs(d2Ur));
    if (d2Uc > 0.0 && d2Ul > 0.0 && d2Ur > 0.0) {
      d2Ulim = SIGN(d2Uc)*MIN(1.25*lim_slope,fabs(d2Uc));
    }
    if (d2Uc < 0.0 && d2Ul < 0.0 && d2Ur < 0.0) {
      d2Ulim = SIGN(d2Uc)*MIN(1.25*lim_slope,fabs(d2Uc));
    }
    Uhalf[j] = 0.5*((U[j-1]+U[j]) - d2Ulim/3.0);
  }

  for (j=jl; j<=ju; j++) {
    Ulv = Uhalf[j  ];
    Urv = Uhalf[j+1];

    qa = (Urv-U[j])*(U[j]-Ulv);
    qb = (U[j-1]-U[j])*(U[j]-U[j+1]);
    if (qa <= 0.0 && qb <= 0.0) {
      qc = 6.0*(U[j] - 0.5*(Ulv+Urv));
      d2U  = -2.0*qc;
      d2Uc = (U[j-1] - 2.0*U[j  ] + U[j+1]);
      d2Ul = (U[j-2] - 2.0*U[j-1] + U[j  ]);
      d2Ur = (U[j  ] - 2.0*U[j+1] + U[j+2]);
      d2Ulim = 0.0;
      lim_slope = MIN(fabs(d2Ul),fabs(d2Ur));
      lim_slope = MIN(fabs(d2Uc),lim_slope);
      if (d2Uc > 0.0 && d2Ul > 0.0 && d2Ur > 0.0 && d2U > 0.0) {
        d2Ulim = SIGN(d2U)*MIN(1.25*lim_slope,fabs(d2U));
      }
      if (d2Uc < 0.0 && d2Ul < 0.0 && d2Ur < 0.0 && d2U < 0.0) {
        d2Ulim = SIGN(d2U)*MIN(1.25*lim_slope,fabs(d2U));
      }
      if (d2U == 0.0) {
        Ulv = U[j];
        Urv = U[j];
      } else {
        Ulv = U[j] + (Ulv - U[j])*d2Ulim/d2U;
        Urv = U[j] + (Urv - U[j])*d2Ulim/d2U;
      }
    }

    qa = (Urv-U[j])*(U[j]-Ulv);
    qb = Urv-Ulv;
    qc = 6.0*(U[j] - 0.5*(Ulv+Urv));
    if (qa <= 0.0) {
      Ulv = U[j];
      Urv = U[j];
    } else if ((qb*qc) > (qb*qb)) {
      Ulv = 3.0*U[j] - 2.0*Urv;
    } else if ((qb*qc) < -(qb*qb)) {
      Urv = 3.0*U[j] - 2.0*Ulv;
    }

    dU = Urv - Ulv;
    U6 = 6.0*(U[j] - 0.5*(Ulv + Urv));

    if (eps > 0.0) { /* eps always > 0 for inner i boundary */
      qx = TWO_3RDS*eps;
      Flux[j+1] = eps*(Urv - 0.75*qx*(dU - (1.0 - qx)*U6));

    } else {         /* eps always < 0 for outer i boundary */
      qx = -TWO_3RDS*eps;
      Flux[j  ] = eps*(Ulv + 0.75*qx*(dU + (1.0 - qx)*U6));
    }
  }

  return;
}

#endif /* THIRD_ORDER */

#endif /* SHEARING_BOX || Cylindrical + Fargo */
