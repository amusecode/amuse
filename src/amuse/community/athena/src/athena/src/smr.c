#include "copyright.h"
/*============================================================================*/
/*! \file smr.c
 *  \brief Functions to handle static mesh refinement (SMR).
 *
 * PURPOSE: Functions to handle static mesh refinement (SMR).
 *
 * REFERENCES:
 * - G. Toth and P.L. Roe, "Divergence and Curl-preserving prolongation and
 *   restriction formulas", JCP 180, 736 (2002)
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - RestrictCorrect(): restricts (averages) fine Grid solution to coarse, and 
 *    corrects cells at fine/coarse boundaries using restricted fine Grid fluxes
 * - Prolongate(): sets BC on fine Grid by prolongation (interpolation) of
 *     coarse Grid solution into fine grid ghost zones
 * - SMR_init(): allocates memory for send/receive buffers
 *
 * PRIVATE FUNCTION PROTOTYPES: 
 * - ProCon() - prolongates conserved variables
 * - ProFld() - prolongates face-centered B field using TR formulas
 * - mcd_slope() - returns monotonized central-difference slope		      */
/*============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef STATIC_MESH_REFINEMENT

static double **send_bufP= NULL; 
static double **send_bufRC=NULL; 
static double ***recv_bufP= NULL;
#ifdef MPI_PARALLEL
static double ***recv_bufRC=NULL;
static MPI_Request ***recv_rq=NULL;
static MPI_Request  **send_rq=NULL;
#endif
static int maxND, *start_addrP;

static ConsS ***GZ[3];
#ifdef MHD
Real **SMRemf1, **SMRemf2, **SMRemf3;
Real3Vect ***BFld[3];
#endif

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   ProCon - prolongates conserved variables
 *   ProFld - prolongates face-centered B field using TR formulas
 *   mcd_slope - returns monotonized central-difference slope
 *============================================================================*/

void ProCon(const ConsS Uim1,const ConsS Ui,  const ConsS Uip1,
            const ConsS Ujm1,const ConsS Ujp1,
            const ConsS Ukm1,const ConsS Ukp1, ConsS PCon[][2][2]);
#ifdef MHD
void ProFld(Real3Vect BGZ[][3][3], Real3Vect PFld[][3][3], 
  const Real dx1c, const Real dx2c, const Real dx3c);
#endif /* MHD */
#ifndef FIRST_ORDER
static Real mcd_slope(const Real vl, const Real vc, const Real vr);
#endif /* FIRST_ORDER */

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void RestrictCorrect(MeshS *pM)
 *  \brief Restricts (averages) fine Grid solution to coarse, and 
 *    corrects cells at fine/coarse boundaries using restricted fine Grid fluxes
 */

void RestrictCorrect(MeshS *pM)
{
  GridS *pG;
  int nl,nd,ncg,dim,nDim,npg,rbufN,start_addr,cnt,nCons,nFlx;
  int i,ii,ics,ice,ips,ipe;
  int j,jj,jcs,jce,jps,jpe;
  int k,kk,kcs,kce,kps,kpe;
  Real q1,q2,q3,fact;
  double *pRcv,*pSnd;
  GridOvrlpS *pCO, *pPO;
#if (NSCALARS > 0)
  int n;
#endif
#ifdef MHD
  int ib,jb,kb,nFld=0;
#endif
#ifdef MPI_PARALLEL
  int ierr,mAddress,mIndex,mCount;
#endif

/* number of dimensions in Grid. */
  nDim=1;
  for (i=1; i<3; i++) if (pM->Nx[i]>1) nDim++;

/* Loop over all Domains, starting at maxlevel */

  for (nl=(pM->NLevels)-1; nl>=0; nl--){

#ifdef MPI_PARALLEL
/* Post non-blocking receives at level nl-1 for data from child Grids at this
 * level (nl).  This data is sent in Step 3 below, and will be read in Step 1
 * at the next iteration of the loop. */ 

  if (nl>0) {
    for (nd=0; nd<(pM->DomainsPerLevel[nl-1]); nd++){
      if (pM->Domain[nl-1][nd].Grid != NULL) {
        pG=pM->Domain[nl-1][nd].Grid;

/* Recv buffer is addressed from 0 for first MPI message, even if NmyCGrid>0.
 * First index alternates between 0 and 1 for even/odd values of nl, since if
 * there are Grids on multiple levels there may be 2 receives posted at once */
        mAddress = 0;
        rbufN = ((nl-1) % 2);
        for (ncg=(pG->NmyCGrid); ncg<(pG->NCGrid); ncg++){
          mIndex = ncg - pG->NmyCGrid;
          ierr = MPI_Irecv(&(recv_bufRC[rbufN][nd][mAddress]),
            pG->CGrid[ncg].nWordsRC, MPI_DOUBLE, pG->CGrid[ncg].ID,
            pG->CGrid[ncg].DomN, pM->Domain[nl-1][nd].Comm_Children,
            &(recv_rq[nl-1][nd][mIndex]));
          mAddress += pG->CGrid[ncg].nWordsRC;
        }

      }
    }
  }
#endif /* MPI_PARALLEL */

/*=== Step 1. Get child solution, inject into parent Grid ====================*/
/* Loop over Domains and child Grids.  Maxlevel domains skip this step because
 * they have NCGrids=0 */

  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){

  if (pM->Domain[nl][nd].Grid != NULL) { /* there is a Grid on this processor */
    pG=pM->Domain[nl][nd].Grid;
    rbufN = (nl % 2);

    for (ncg=0; ncg<(pG->NCGrid); ncg++){

/*--- Step 1a. Get restricted solution and fluxes. ---------------------------*/

/* If child Grid is on this processor, set pointer to start of send buffer
 * loaded by this child Grid in Step 3 below during last iteration of loop. */

      if (ncg < pG->NmyCGrid) {
        pCO=(GridOvrlpS*)&(pG->CGrid[ncg]);
        pRcv = (double*)&(send_bufRC[pCO->DomN][0]);
      } else {

#ifdef MPI_PARALLEL
/* Check non-blocking receives posted above for restricted solution from child
 * Grids, sent in Step 3 during last iteration of loop over nl.  Accept messages
 * in any order. */

        mCount = pG->NCGrid - pG->NmyCGrid;
        ierr = MPI_Waitany(mCount,recv_rq[nl][nd],&mIndex,MPI_STATUS_IGNORE);
        if(mIndex == MPI_UNDEFINED){
          ath_error("[RestCorr]: Invalid request index nl=%i nd=%i\n",nl,nd);
        }
      
/* Recv buffer is addressed from 0 for first MPI message, even if NmyCGrid>0 */
        mAddress = 0;
        mIndex += pG->NmyCGrid;
        for (i=pG->NmyCGrid; i<mIndex; i++) mAddress += pG->CGrid[i].nWordsRC;
        pCO=(GridOvrlpS*)&(pG->CGrid[mIndex]);
        pRcv = (double*)&(recv_bufRC[rbufN][nd][mAddress]);
#else
/* If not MPI_PARALLEL, and child Grid not on this processor, then error */

        ath_error("[RestCorr]: no Child grid on Domain[%d][%d]\n",nl,nd);
#endif /* MPI_PARALLEL */
      }

/* Get coordinates ON THIS GRID of overlap region of child Grid */

      ics = pCO->ijks[0];
      ice = pCO->ijke[0];
      jcs = pCO->ijks[1];
      jce = pCO->ijke[1];
      kcs = pCO->ijks[2];
      kce = pCO->ijke[2];

/*--- Step 1b. Set conserved variables on parent Grid ------------------------*/

      for (k=kcs; k<=kce; k++) {
        for (j=jcs; j<=jce; j++) {
          for (i=ics; i<=ice; i++) {
            pG->U[k][j][i].d  = *(pRcv++);
            pG->U[k][j][i].M1 = *(pRcv++);
            pG->U[k][j][i].M2 = *(pRcv++);
            pG->U[k][j][i].M3 = *(pRcv++);
#ifndef BAROTROPIC
            pG->U[k][j][i].E  = *(pRcv++);
#endif /* BAROTROPIC */
#ifdef MHD
            pG->U[k][j][i].B1c = *(pRcv++);
            pG->U[k][j][i].B2c = *(pRcv++);
            pG->U[k][j][i].B3c = *(pRcv++);
#endif /* MHD */
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) pG->U[k][j][i].s[n] = *(pRcv++);
#endif
          }
        }
      }

#ifdef MHD
/*--- Step 1c. Set face-centered fields. -------------------------------------*/
/* Also recompute cell-centered fields based on new face-centered values.
 * Do not set face-centered fields on fine/coarse boundaries (e.g. ics and ice+1
 * for B1i), but use the flux-correction step below to do that */

      if (nDim == 1){ /* 1D problem - set face-centered to cell-centered */
        for (i=ics; i<=ice; i++) {
          pG->B1i[kcs][jcs][i] = pG->U[kcs][jcs][i].B1c;
          pG->B2i[kcs][jcs][i] = pG->U[kcs][jcs][i].B2c;
          pG->B3i[kcs][jcs][i] = pG->U[kcs][jcs][i].B3c;
        }
      } else {

        if (nDim == 2){  /* 2D problem */
/* Correct B1i */
          for (j=jcs  ; j<=jce; j++) {
          for (i=ics+1; i<=ice; i++) {
            /* Set B1i at ics if no flux correction will be made.  Increment
             * pointer even if value in Rcv pointer is ignored. */
            if (i==ics+1){
              if (pCO->myFlx[0] == NULL) {pG->B1i[kcs][j][ics] = *(pRcv++);}
              else {pRcv++;}
            }

            pG->B1i[kcs][j][i] = *(pRcv++);

            /* Set B1i at ice+1 if no flux correction will be made.  Increment
             * pointer even if value in Rcv pointer is ignored. */
            if (i==ice){
              if (pCO->myFlx[1] == NULL) {pG->B1i[kcs][j][ice+1] = *(pRcv++);}
              else {pRcv++;}
            }
          }}

/* Correct B2i */
          /* Set B2i at jcs if no flux correction will be made.  Increment
           * pointer even if value in Rcv pointer is ignored. */
          for (i=ics; i<=ice; i++) {
            if (pCO->myFlx[2] == NULL) {pG->B2i[kcs][jcs][i] = *(pRcv++);}
            else {pRcv++;}
          }

          for (j=jcs+1; j<=jce; j++) {
          for (i=ics  ; i<=ice; i++) {
            pG->B2i[kcs][j][i] = *(pRcv++);
          }}

          /* Set B2i at jce+1 if no flux correction will be made.  Increment
           * pointer even if value in Rcv pointer is ignored. */
          for (i=ics; i<=ice; i++) {
            if (pCO->myFlx[3] == NULL) {pG->B2i[kcs][jce+1][i] = *(pRcv++);}
            else {pRcv++;}
          }

/* Set cell-centered fields */
          for (j=jcs; j<=jce; j++) {
          for (i=ics; i<=ice; i++) {
            pG->U[kcs][j][i].B1c=0.5*(pG->B1i[kcs][j][i] +pG->B1i[kcs][j][i+1]);
            pG->U[kcs][j][i].B2c=0.5*(pG->B2i[kcs][j][i] +pG->B2i[kcs][j+1][i]);
            pG->B3i[kcs][j][i] = pG->U[kcs][j][i].B3c;
          }}

        } else { /* 3D problem */
/* Correct B1i */
          for (k=kcs  ; k<=kce; k++) {
          for (j=jcs  ; j<=jce; j++) {
          for (i=ics+1; i<=ice; i++) {
            /* Set B1i at ics if no flux correction will be made.  Increment
             * pointer even if value in Rcv pointer is ignored. */
            if (i==ics+1){
              if (pCO->myFlx[0] == NULL) {pG->B1i[k][j][ics] = *(pRcv++);}
              else {pRcv++;}
            }

            pG->B1i[k][j][i] = *(pRcv++);

            /* Set B1i at ice+1 if no flux correction will be made.  Increment
             * pointer even if value in Rcv pointer is ignored. */
            if (i==ice){
              if (pCO->myFlx[1] == NULL) {pG->B1i[k][j][ice+1] = *(pRcv++);}
              else {pRcv++;}
            }

          }}}

/* Correct B2i */
          for (k=kcs  ; k<=kce; k++) {
            /* Set B2i at jcs if no flux correction will be made.  Increment
             * pointer even if value in Rcv pointer is ignored. */
            for (i=ics; i<=ice; i++) {
              if (pCO->myFlx[2] == NULL) {pG->B2i[k][jcs][i] = *(pRcv++);}
              else {pRcv++;}
            }

            for (j=jcs+1; j<=jce; j++) {
            for (i=ics  ; i<=ice; i++) {
              pG->B2i[k][j][i] = *(pRcv++);
            }}

            /* Set B2i at jce+1 if no flux correction will be made.  Increment
             * pointer even if value in Rcv pointer is ignored. */
            for (i=ics; i<=ice; i++) {
              if (pCO->myFlx[3] == NULL) {pG->B2i[k][jce+1][i] = *(pRcv++);}
              else {pRcv++;}
            }
          }

/* Correct B3i */
          /* Set B3i at kcs if no flux correction will be made.  Increment
           * pointer even if value in Rcv pointer is ignored. */
          for (j=jcs; j<=jce; j++) {
          for (i=ics; i<=ice; i++) {
            if (pCO->myFlx[4] == NULL) {pG->B3i[kcs][j][i] = *(pRcv++);}
            else {pRcv++;}
          }}

          for (k=kcs+1; k<=kce; k++) {
          for (j=jcs  ; j<=jce; j++) {
          for (i=ics  ; i<=ice; i++) {
            pG->B3i[k][j][i] = *(pRcv++);
          }}}

          /* Set B3i at kce+1 if no flux correction will be made.  Increment
           * pointer even if value in Rcv pointer is ignored. */
          for (j=jcs; j<=jce; j++) {
          for (i=ics; i<=ice; i++) {
            if (pCO->myFlx[5] == NULL) {pG->B3i[kce+1][j][i] = *(pRcv++);}
            else {pRcv++;}
          }}

/* Set cell-centered fields */
          for (k=kcs; k<=kce; k++) {
          for (j=jcs; j<=jce; j++) {
          for (i=ics; i<=ice; i++) {
            pG->U[k][j][i].B1c = 0.5*(pG->B1i[k][j][i] + pG->B1i[k][j][i+1]);
            pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j+1][i]);
            pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i] + pG->B3i[k+1][j][i]);
          }}}

        }

      }
#endif /* MHD */

/*=== Step 2. Flux correction ================================================*/
/*--- Step 2a. Flux correction to conserved variables. -----------------------*/
/* Subtract off flux from this Grid, then add in restricted flux, both
 * multiplied by dt/dx with sign chosen depending on whether flux is L/R of cell
 * center. */

/* Correct solution at x1-faces */

      for (dim=0; dim<2; dim++){
        if (pCO->myFlx[dim] != NULL) {
          if (dim == 0) {i=ics-1; q1=-(pG->dt/pG->dx1);}
          if (dim == 1) {i=ice+1; q1= (pG->dt/pG->dx1);}
          for (k=kcs, kk=0; k<=kce; k++, kk++) {
          for (j=jcs, jj=0; j<=jce; j++, jj++) {
            pG->U[k][j][i].d  -= q1*(pCO->myFlx[dim][kk][jj].d  - *(pRcv++));
            pG->U[k][j][i].M1 -= q1*(pCO->myFlx[dim][kk][jj].M1 - *(pRcv++));
            pG->U[k][j][i].M2 -= q1*(pCO->myFlx[dim][kk][jj].M2 - *(pRcv++));
            pG->U[k][j][i].M3 -= q1*(pCO->myFlx[dim][kk][jj].M3 - *(pRcv++));
#ifndef BAROTROPIC
            pG->U[k][j][i].E  -= q1*(pCO->myFlx[dim][kk][jj].E  - *(pRcv++));
#endif /* BAROTROPIC */
#ifdef MHD
            pG->U[k][j][i].B1c -= q1*(pCO->myFlx[dim][kk][jj].B1c - *(pRcv++));
            pG->U[k][j][i].B2c -= q1*(pCO->myFlx[dim][kk][jj].B2c - *(pRcv++));
            pG->U[k][j][i].B3c -= q1*(pCO->myFlx[dim][kk][jj].B3c - *(pRcv++));
#endif /* MHD */
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) pG->U[k][j][i].s[n] -= 
              q1*(pCO->myFlx[dim][kk][jj].s[n] - *(pRcv++));
#endif
          }}
        }
      }

/* Correct solution at x2-faces */

      for (dim=2; dim<4; dim++){
        if (pCO->myFlx[dim] != NULL) {
          if (dim == 2) {j=jcs-1; q2=-(pG->dt/pG->dx2);}
          if (dim == 3) {j=jce+1; q2= (pG->dt/pG->dx2);}
          for (k=kcs, kk=0; k<=kce; k++, kk++) {
          for (i=ics, ii=0; i<=ice; i++, ii++) {
            pG->U[k][j][i].d  -= q2*(pCO->myFlx[dim][kk][ii].d  - *(pRcv++));
            pG->U[k][j][i].M1 -= q2*(pCO->myFlx[dim][kk][ii].M1 - *(pRcv++));
            pG->U[k][j][i].M2 -= q2*(pCO->myFlx[dim][kk][ii].M2 - *(pRcv++));
            pG->U[k][j][i].M3 -= q2*(pCO->myFlx[dim][kk][ii].M3 - *(pRcv++));
#ifndef BAROTROPIC
            pG->U[k][j][i].E  -= q2*(pCO->myFlx[dim][kk][ii].E  - *(pRcv++));
#endif /* BAROTROPIC */
#ifdef MHD
            pG->U[k][j][i].B1c -= q2*(pCO->myFlx[dim][kk][ii].B1c - *(pRcv++));
            pG->U[k][j][i].B2c -= q2*(pCO->myFlx[dim][kk][ii].B2c - *(pRcv++));
            pG->U[k][j][i].B3c -= q2*(pCO->myFlx[dim][kk][ii].B3c - *(pRcv++));
#endif /* MHD */
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) pG->U[k][j][i].s[n] -= 
              q2*(pCO->myFlx[dim][kk][ii].s[n] - *(pRcv++));
#endif
          }}
        }
      }

/* Correct solution at x3-faces */

      for (dim=4; dim<6; dim++){
        if (pCO->myFlx[dim] != NULL) {
          if (dim == 4) {k=kcs-1; q3=-(pG->dt/pG->dx3);}
          if (dim == 5) {k=kce+1; q3= (pG->dt/pG->dx3);}
          for (j=jcs, jj=0; j<=jce; j++, jj++) {
          for (i=ics, ii=0; i<=ice; i++, ii++) {
            pG->U[k][j][i].d  -= q3*(pCO->myFlx[dim][jj][ii].d  - *(pRcv++));
            pG->U[k][j][i].M1 -= q3*(pCO->myFlx[dim][jj][ii].M1 - *(pRcv++));
            pG->U[k][j][i].M2 -= q3*(pCO->myFlx[dim][jj][ii].M2 - *(pRcv++));
            pG->U[k][j][i].M3 -= q3*(pCO->myFlx[dim][jj][ii].M3 - *(pRcv++));
#ifndef BAROTROPIC
            pG->U[k][j][i].E  -= q3*(pCO->myFlx[dim][jj][ii].E  - *(pRcv++));
#endif /* BAROTROPIC */
#ifdef MHD
            pG->U[k][j][i].B1c -= q3*(pCO->myFlx[dim][jj][ii].B1c - *(pRcv++));
            pG->U[k][j][i].B2c -= q3*(pCO->myFlx[dim][jj][ii].B2c - *(pRcv++));
            pG->U[k][j][i].B3c -= q3*(pCO->myFlx[dim][jj][ii].B3c - *(pRcv++));
#endif /* MHD */
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) pG->U[k][j][i].s[n] -= 
              q3*(pCO->myFlx[dim][jj][ii].s[n] - *(pRcv++));
#endif
          }}
        }
      }

/*--- Step 2b. Flux correction to face-centered fields. ----------------------*/
/* Only needed in 2D/3D.  Add EMF on this Grid back in, then subtract new, both
 * multiplied by dt/dx with sign chosen depending on whether EMF is to L/R of
 * cell center. */

#ifdef MHD
      if (nDim > 1) {

/* On x1-faces; correct B1i, B2i, B3i */

      for (dim=0; dim<2; dim++){
        if (pCO->myEMF3[dim] != NULL) {
          if (dim == 0) {
            i=ics-1; ib = ics;
            q1 = -(pG->dt/pG->dx1);
            q2 = -(pG->dt/pG->dx2); 
            q3 = -(pG->dt/pG->dx3);
          } 
          if (dim == 1) {
            i=ice+1; ib = ice+1;
            q1 =  (pG->dt/pG->dx1);
            q2 = -(pG->dt/pG->dx2);
            q3 = -(pG->dt/pG->dx3);
          }

          for (k=kcs, kk=0; k<=kce  ; k++, kk++) {
            for (j=jcs, jj=0; j<=jce; j++, jj++) {
              SMRemf3[kk][jj] = *(pRcv++);
              pG->B2i[k][j][i]+= q1*(pCO->myEMF3[dim][kk][jj] -SMRemf3[kk][jj]);
            }
           /* only update B2i[jce+1] if ox2 [dim=3] boundary is at edge of child
            * Domain (so ox2 boundary is between fine/coarse Grids), or at edge
            * of this Grid. Otherwise ox2 boundary is between MPI Grids on child
            * Domain, and it will be updated as B2i[jcs] by other child Grid */
            jj = jce-jcs+1;
            SMRemf3[kk][jj] = *(pRcv++);
            if((pCO->myEMF3[3] != NULL) || (jce==pG->je)){
              pG->B2i[k][jce+1][i] +=
                q1*(pCO->myEMF3[dim][kk][jj] - SMRemf3[kk][jj]);
            }
          }

          for (k=kcs, kk=0; k<=kce; k++, kk++) {
          for (j=jcs, jj=0; j<=jce; j++, jj++) {
            pG->B1i[k][j][ib] -=
              q2*(pCO->myEMF3[dim][kk][jj+1] - SMRemf3[kk][jj+1]) -
              q2*(pCO->myEMF3[dim][kk][jj  ] - SMRemf3[kk][jj  ]);
          }}

          if (nDim == 3){ /* 3D problem */

            for (k=kcs, kk=0; k<=kce; k++, kk++) {
            for (j=jcs, jj=0; j<=jce; j++, jj++) {
              SMRemf2[kk][jj] = *(pRcv++);
              pG->B3i[k][j][i] -= q1*(pCO->myEMF2[dim][kk][jj]-SMRemf2[kk][jj]);
            }}
           /* only update B3i[kce+1] if ox3 [dim=5] boundary is at edge of child
            * Domain (so ox3 boundary is between fine/coarse Grids), or at edge
            * of this Grid. Otherwise ox3 boundary is between MPI Grids on child
            * Domain, and it will be updated as B3i[kcs] by other child Grid */
            kk = kce-kcs+1;
            for (j=jcs, jj=0; j<=jce  ; j++, jj++) {
              SMRemf2[kk][jj] = *(pRcv++); 
              if((pCO->myEMF2[5] != NULL) || (kce==pG->ke)){
                pG->B3i[kce+1][j][i] -= 
                  q1*(pCO->myEMF2[dim][kk][jj] - SMRemf2[kk][jj]);
              }
            }

            for (k=kcs, kk=0; k<=kce; k++, kk++) {
            for (j=jcs, jj=0; j<=jce; j++, jj++) {
              pG->B1i[k][j][ib] +=
                q3*(pCO->myEMF2[dim][kk+1][jj] - SMRemf2[kk+1][jj]) -
                q3*(pCO->myEMF2[dim][kk  ][jj] - SMRemf2[kk  ][jj]);
            }}

            for (k=kcs-1; k<=kce+1; k++) {
            for (j=jcs; j<=jce; j++) {
              pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i] + pG->B3i[k+1][j][i]);
            }}

          }

          for (k=kcs; k<=kce; k++) {
            for (j=jcs; j<=jce; j++) {
              pG->U[k][j][ib  ].B1c=0.5*(pG->B1i[k][j][ib]+pG->B1i[k][j][ib+1]);
              pG->U[k][j][ib-1].B1c=0.5*(pG->B1i[k][j][ib-1]+pG->B1i[k][j][ib]);
            }
            for (j=jcs-1; j<=jce+1; j++) {
              pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j+1][i]);
            }
          }
        }
      }

/* On x2-faces; correct B1i, B2i, B3i */

      for (dim=2; dim<4; dim++){
        if (pCO->myEMF3[dim] != NULL) {
          if (dim == 2) {
            j=jcs-1; jb=jcs; 
            q1 = -(pG->dt/pG->dx1);
            q2 = -(pG->dt/pG->dx2);
            q3 = -(pG->dt/pG->dx3);
          }
          if (dim == 3) {
            j=jce+1; jb=jce+1; 
            q1 = -(pG->dt/pG->dx1);
            q2 =  (pG->dt/pG->dx2);
            q3 = -(pG->dt/pG->dx3);
          }

          for (k=kcs, kk=0; k<=kce; k++, kk++) {
            for (i=ics, ii=0; i<=ice; i++, ii++) {
              SMRemf3[kk][ii] = *(pRcv++);
              pG->B1i[k][j][i]-= q2*(pCO->myEMF3[dim][kk][ii] -SMRemf3[kk][ii]);
            }
           /* only update B1i[ice+1] if ox1 [dim=1] boundary is at edge of child
            * Domain (so ox1 boundary is between fine/coarse Grids), or at edge
            * of this Grid. Otherwise ox1 boundary is between MPI Grids on child
            * Domain, and it will be updated as B1i[ics] by other child Grid */
            ii = ice-ics+1;
            SMRemf3[kk][ii] = *(pRcv++);
            if ((pCO->myEMF3[1] != NULL) || (ice==pG->ie)){
              pG->B1i[k][j][ice+1] -=
                q2*(pCO->myEMF3[dim][kk][ii] - SMRemf3[kk][ii]);
            }
          }
          for (k=kcs, kk=0; k<=kce; k++, kk++) {
          for (i=ics, ii=0; i<=ice; i++, ii++) {
            pG->B2i[k][jb][i] += 
              q1*(pCO->myEMF3[dim][kk][ii+1] - SMRemf3[kk][ii+1]) -
              q1*(pCO->myEMF3[dim][kk][ii  ] - SMRemf3[kk][ii  ]);
          }}

          if (nDim == 3){ /* 3D problem */

            for (k=kcs, kk=0; k<=kce; k++, kk++) {
            for (i=ics, ii=0; i<=ice; i++, ii++) {
              SMRemf1[kk][ii] = *(pRcv++);
              pG->B3i[k][j][i] += q2*(pCO->myEMF1[dim][kk][ii]-SMRemf1[kk][ii]);
            }}
           /* only update B3i[kce+1] if ox3 [dim=5] boundary is at edge of child
            * Domain (so ox3 boundary is between fine/coarse Grids), or at edge
            * of this Grid. Otherwise ox3 boundary is between MPI Grids on child
            * Domain, and it will be updated as B3i[kcs] by other child Grid */
            kk = kce-kcs+1;
            for (i=ics, ii=0; i<=ice  ; i++, ii++) {
              SMRemf1[kk][ii] = *(pRcv++);
              if((pCO->myEMF1[5] != NULL) || (kce==pG->ke)){
                pG->B3i[kce+1][j][i] +=
                  q2*(pCO->myEMF1[dim][kk][ii] - SMRemf1[kk][ii]);
              }
            }

            for (k=kcs, kk=0; k<=kce; k++, kk++) {
            for (i=ics, ii=0; i<=ice; i++, ii++) {
              pG->B2i[k][jb][i] -= 
                q3*(pCO->myEMF1[dim][kk+1][ii] - SMRemf1[kk+1][ii]) -
                q3*(pCO->myEMF1[dim][kk  ][ii] - SMRemf1[kk  ][ii]);
            }}

            for (k=kcs-1; k<=kce+1; k++) {
            for (i=ics; i<=ice; i++) {
              pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i] + pG->B3i[k+1][j][i]);
            }}

          }

          for (k=kcs; k<=kce; k++) {
            for (i=ics-1; i<=ice+1; i++) {
              pG->U[k][j][i].B1c = 0.5*(pG->B1i[k][j][i] + pG->B1i[k][j][i+1]);
            }
            for (i=ics; i<=ice; i++) {
              pG->U[k][jb  ][i].B2c=0.5*(pG->B2i[k][jb][i]+pG->B2i[k][jb+1][i]);
              pG->U[k][jb-1][i].B2c=0.5*(pG->B2i[k][jb-1][i]+pG->B2i[k][jb][i]);
            }
          }
        }
      }

/* On x3-faces; correct B1i, B2i, B3i.  Must be a 3D problem. */

      for (dim=4; dim<6; dim++){
        if (pCO->myEMF1[dim] != NULL) {
          if (dim == 4) {
            k=kcs-1; kb=kcs; 
            q1 = -(pG->dt/pG->dx1);
            q2 = -(pG->dt/pG->dx2);
            q3 = -(pG->dt/pG->dx3);
          }
          if (dim == 5) {
            k=kce+1; kb=kce+1; 
            q1 = -(pG->dt/pG->dx1);
            q2 = -(pG->dt/pG->dx2);
            q3 =  (pG->dt/pG->dx3);
          }

          for (j=jcs, jj=0; j<=jce; j++, jj++) {
          for (i=ics, ii=0; i<=ice; i++, ii++) {
            SMRemf1[jj][ii] = *(pRcv++);
            pG->B2i[k][j][i] -= q3*(pCO->myEMF1[dim][jj][ii] - SMRemf1[jj][ii]);
          }}
         /* only update B2i[jce+1] if ox2 [dim=3] boundary is at edge of child
          * Domain (so ox2 boundary is between fine/coarse Grids), or at edge
          * of this Grid. Otherwise ox2 boundary is between MPI Grids on child
          * Domain, and it will be updated as B2i[jcs] by other child Grid */
          jj = jce-jcs+1;
          for (i=ics, ii=0; i<=ice  ; i++, ii++) {
            SMRemf1[jj][ii] = *(pRcv++);
            if((pCO->myEMF1[3] != NULL) || (jce==pG->je)){
              pG->B2i[k][jce+1][i] -=
                q3*(pCO->myEMF1[dim][jj][ii] - SMRemf1[jj][ii]);
            }
          }

          for (j=jcs, jj=0; j<=jce; j++, jj++) {
            for (i=ics, ii=0; i<=ice; i++, ii++) {
              SMRemf2[jj][ii] = *(pRcv++);
              pG->B1i[k][j][i]+= q3*(pCO->myEMF2[dim][jj][ii] -SMRemf2[jj][ii]);
            }
           /* only update B1i[ice+1] if ox1 [dim=1] boundary is at edge of child
            * Domain (so ox1 boundary is between fine/coarse Grids), or at edge
            * of this Grid. Otherwise ox1 boundary is between MPI Grids on child
            * Domain, and it will be updated as B1i[ics] by other child Grid */
            ii = ice-ics+1;
            SMRemf2[jj][ii] = *(pRcv++);
            if ((pCO->myEMF2[1] != NULL) || (ice==pG->ie)){
              pG->B1i[k][j][ice+1] +=
                q3*(pCO->myEMF2[dim][jj][ii] -SMRemf2[jj][ii]);
            }
          }

          for (j=jcs, jj=0; j<=jce; j++, jj++) {
          for (i=ics, ii=0; i<=ice; i++, ii++) {
            pG->B3i[kb][j][i] +=
              q2*(pCO->myEMF1[dim][jj+1][ii  ] - SMRemf1[jj+1][ii  ]) -
              q2*(pCO->myEMF1[dim][jj  ][ii  ] - SMRemf1[jj  ][ii  ]) -
              q1*(pCO->myEMF2[dim][jj  ][ii+1] - SMRemf2[jj  ][ii+1]) +
              q1*(pCO->myEMF2[dim][jj  ][ii  ] - SMRemf2[jj  ][ii  ]);
          }}

          for (j=jcs; j<=jce; j++) {
          for (i=ics-1; i<=ice+1; i++) {
            pG->U[k][j][i].B1c = 0.5*(pG->B1i[k][j][i] + pG->B1i[k][j][i+1]);
          }}

          for (j=jcs-1; j<=jce+1; j++) {
          for (i=ics; i<=ice; i++) {
            pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j+1][i]);
          }}

          for (j=jcs; j<=jce; j++) {
          for (i=ics; i<=ice; i++) {
            pG->U[kb  ][j][i].B3c=0.5*(pG->B3i[kb][j][i] + pG->B3i[kb+1][j][i]);
            pG->U[kb-1][j][i].B3c=0.5*(pG->B3i[kb-1][j][i] + pG->B3i[kb][j][i]);
          }}
        }
      }

      } /* end test for 2D/3D */
#endif /* MHD */

    }  /* end loop over child grids */
  }} /* end loop over Domains */

/*=== Step 3. Restrict child solution and fluxes and send ====================*/
/* Loop over all Domains and parent Grids.  Maxlevel grids skip straight to this
 * step to start the chain of communication.  Root (level=0) skips this step
 * since it has NPGrid=0.  If there is a parent Grid on this processor, it will
 * be first in the PGrid array, so it will be at start of send_bufRC */

  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){

  if (pM->Domain[nl][nd].Grid != NULL) { /* there is a Grid on this processor */
    pG=pM->Domain[nl][nd].Grid;          /* set pointer to this Grid */
    start_addr=0;

    for (npg=0; npg<(pG->NPGrid); npg++){
      pPO=(GridOvrlpS*)&(pG->PGrid[npg]);    /* ptr to Grid overlap */
      cnt = 0;

/* Get coordinates ON THIS GRID of overlap region of parent Grid */

      ips = pPO->ijks[0];
      ipe = pPO->ijke[0];
      jps = pPO->ijks[1];
      jpe = pPO->ijke[1];
      kps = pPO->ijks[2];
      kpe = pPO->ijke[2];

/*--- Step 3a. Restrict conserved variables  ---------------------------------*/
/* 1D/2D/3D problem: Conservative average of conserved variables in x1. */

      pSnd = (double*)&(send_bufRC[nd][start_addr]);
      for (k=kps; k<=kpe; k+=2) {
      for (j=jps; j<=jpe; j+=2) {
      for (i=ips; i<=ipe; i+=2) {
        *(pSnd++) = pG->U[k][j][i].d  + pG->U[k][j][i+1].d;
        *(pSnd++) = pG->U[k][j][i].M1 + pG->U[k][j][i+1].M1;
        *(pSnd++) = pG->U[k][j][i].M2 + pG->U[k][j][i+1].M2;
        *(pSnd++) = pG->U[k][j][i].M3 + pG->U[k][j][i+1].M3;
#ifndef BAROTROPIC
        *(pSnd++) = pG->U[k][j][i].E  + pG->U[k][j][i+1].E;
#endif
#ifdef MHD
        *(pSnd++) = pG->U[k][j][i].B1c + pG->U[k][j][i+1].B1c;
        *(pSnd++) = pG->U[k][j][i].B2c + pG->U[k][j][i+1].B2c;
        *(pSnd++) = pG->U[k][j][i].B3c + pG->U[k][j][i+1].B3c;
#endif
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) 
          *(pSnd++) = pG->U[k][j][i].s[n] + pG->U[k][j][i+1].s[n];
#endif
      }}}
      fact = 0.5;
      nCons = (ipe-ips+1)*(NVAR)/2;

/* 2D/3D problem: Add conservative average in x2 */

      if (pG->Nx[1] > 1) {
        pSnd = (double*)&(send_bufRC[nd][start_addr]); /* restart pointer */
        for (k=kps; k<=kpe; k+=2) {
        for (j=jps; j<=jpe; j+=2) {
        for (i=ips; i<=ipe; i+=2) {
          *(pSnd++) += pG->U[k][j+1][i].d  + pG->U[k][j+1][i+1].d;
          *(pSnd++) += pG->U[k][j+1][i].M1 + pG->U[k][j+1][i+1].M1;
          *(pSnd++) += pG->U[k][j+1][i].M2 + pG->U[k][j+1][i+1].M2;
          *(pSnd++) += pG->U[k][j+1][i].M3 + pG->U[k][j+1][i+1].M3;
#ifndef BAROTROPIC
          *(pSnd++) += pG->U[k][j+1][i].E  + pG->U[k][j+1][i+1].E;
#endif
#ifdef MHD
          *(pSnd++) += pG->U[k][j+1][i].B1c + pG->U[k][j+1][i+1].B1c;
          *(pSnd++) += pG->U[k][j+1][i].B2c + pG->U[k][j+1][i+1].B2c;
          *(pSnd++) += pG->U[k][j+1][i].B3c + pG->U[k][j+1][i+1].B3c;
#endif
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) 
            *(pSnd++) += pG->U[k][j+1][i].s[n] + pG->U[k][j+1][i+1].s[n];
#endif
        }}}
        fact = 0.25;
        nCons = (jpe-jps+1)*(ipe-ips+1)*(NVAR)/4;
      }

/* 3D problem: Add conservative average in x3 */

      if (pG->Nx[2] > 1) {
        pSnd = (double*)&(send_bufRC[nd][start_addr]);  /* restart pointer */
        for (k=kps; k<=kpe; k+=2) {
        for (j=jps; j<=jpe; j+=2) {
        for (i=ips; i<=ipe; i+=2) {
          *(pSnd++) += pG->U[k+1][j  ][i].d  + pG->U[k+1][j  ][i+1].d +
                       pG->U[k+1][j+1][i].d  + pG->U[k+1][j+1][i+1].d;
          *(pSnd++) += pG->U[k+1][j  ][i].M1 + pG->U[k+1][j  ][i+1].M1 +
                       pG->U[k+1][j+1][i].M1 + pG->U[k+1][j+1][i+1].M1;
          *(pSnd++) += pG->U[k+1][j  ][i].M2 + pG->U[k+1][j  ][i+1].M2 +
                       pG->U[k+1][j+1][i].M2 + pG->U[k+1][j+1][i+1].M2;
          *(pSnd++) += pG->U[k+1][j  ][i].M3 + pG->U[k+1][j  ][i+1].M3 +
                       pG->U[k+1][j+1][i].M3 + pG->U[k+1][j+1][i+1].M3;
#ifndef BAROTROPIC
          *(pSnd++) += pG->U[k+1][j  ][i].E  + pG->U[k+1][j  ][i+1].E +
                       pG->U[k+1][j+1][i].E  + pG->U[k+1][j+1][i+1].E;
#endif
#ifdef MHD
          *(pSnd++) += pG->U[k+1][j  ][i].B1c + pG->U[k+1][j  ][i+1].B1c +
                       pG->U[k+1][j+1][i].B1c + pG->U[k+1][j+1][i+1].B1c;
          *(pSnd++) += pG->U[k+1][j  ][i].B2c + pG->U[k+1][j  ][i+1].B2c +
                       pG->U[k+1][j+1][i].B2c + pG->U[k+1][j+1][i+1].B2c;
          *(pSnd++) += pG->U[k+1][j  ][i].B3c + pG->U[k+1][j  ][i+1].B3c +
                       pG->U[k+1][j+1][i].B3c + pG->U[k+1][j+1][i+1].B3c;
#endif
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) 
            *(pSnd++) += pG->U[k+1][j  ][i].s[n] + pG->U[k+1][j  ][i+1].s[n] +
                         pG->U[k+1][j+1][i].s[n] + pG->U[k+1][j+1][i+1].s[n];
#endif
        }}}
        fact = 0.125;
        nCons = (kpe-kps+1)*(jpe-jps+1)*(ipe-ips+1)*(NVAR)/8;
      }

/* reset pointer to beginning and normalize averages */
      pSnd = (double*)&(send_bufRC[nd][start_addr]);  
      for (i=start_addr; i<(start_addr+nCons); i++) *(pSnd++) *= fact;
      cnt = nCons;

#ifdef MHD
/*--- Step 3b. Restrict face-centered fields  --------------------------------*/
/* Average face-centered magnetic fields.  Send fields at Grid boundaries
 * (e.g. ips/ipe+1 for B1i) in case they are needed at edges of MPI blocks on
 * same level.  Step 1c decides whether to use or ignore them. */

      if (nDim == 3) { /* 3D problem, restrict B1i,B2i,B3i */

        for (k=kps; k<=kpe  ; k+=2) {
        for (j=jps; j<=jpe  ; j+=2) {
        for (i=ips; i<=ipe+1; i+=2) {
          *(pSnd++) = 0.25*(pG->B1i[k  ][j][i] + pG->B1i[k  ][j+1][i]
                         +  pG->B1i[k+1][j][i] + pG->B1i[k+1][j+1][i]);
        }}}
        for (k=kps; k<=kpe  ; k+=2) {
        for (j=jps; j<=jpe+1; j+=2) {
        for (i=ips; i<=ipe  ; i+=2) {
          *(pSnd++) = 0.25*(pG->B2i[k  ][j][i] + pG->B2i[k  ][j][i+1]
                          + pG->B2i[k+1][j][i] + pG->B2i[k+1][j][i+1]);
        }}}
        for (k=kps; k<=kpe+1; k+=2) {
        for (j=jps; j<=jpe  ; j+=2) {
        for (i=ips; i<=ipe  ; i+=2) {
          *(pSnd++) = 0.25*(pG->B3i[k][j  ][i] + pG->B3i[k][j  ][i+1]
                          + pG->B3i[k][j+1][i] + pG->B3i[k][j+1][i+1]);
        }}}
        nFld = ((kpe-kps+1)/2    )*((jpe-jps+1)/2    )*((ipe-ips+1)/2 + 1) +
               ((kpe-kps+1)/2    )*((jpe-jps+1)/2 + 1)*((ipe-ips+1)/2    ) +
               ((kpe-kps+1)/2 + 1)*((jpe-jps+1)/2    )*((ipe-ips+1)/2    );

      } else {

        if (nDim == 2) { /* 2D problem, restrict B1i,B2i */
          for (j=jps; j<=jpe  ; j+=2) {
          for (i=ips; i<=ipe+1; i+=2) {
            *(pSnd++) = 0.5*(pG->B1i[kps][j][i] + pG->B1i[kps][j+1][i]);
          }}
          for (j=jps; j<=jpe+1; j+=2) {
          for (i=ips; i<=ipe  ; i+=2) {
            *(pSnd++) = 0.5*(pG->B2i[kps][j][i] + pG->B2i[kps][j][i+1]);
          }}
          nFld = ((jpe-jps+1)/2    )*((ipe-ips+1)/2 + 1) +
                 ((jpe-jps+1)/2 + 1)*((ipe-ips+1)/2    );
        }

      }
      cnt += nFld;
#endif /* MHD */

/*--- Step 3c. Restrict fluxes of conserved variables ------------------------*/
/*---------------- Restrict fluxes at x1-faces -------------------------------*/

      for (dim=0; dim<2; dim++){
      if (pPO->myFlx[dim] != NULL) {

      pSnd = (double*)&(send_bufRC[nd][(start_addr+cnt)]);  

      if (nDim == 1) {  /*----- 1D problem -----*/

        *(pSnd++) = pPO->myFlx[dim][kps][jps].d ;
        *(pSnd++) = pPO->myFlx[dim][kps][jps].M1;
        *(pSnd++) = pPO->myFlx[dim][kps][jps].M2;
        *(pSnd++) = pPO->myFlx[dim][kps][jps].M3;
#ifndef BAROTROPIC
        *(pSnd++) = pPO->myFlx[dim][kps][jps].E ;
#endif
#ifdef MHD
        *(pSnd++) = pPO->myFlx[dim][kps][jps].B1c;
        *(pSnd++) = pPO->myFlx[dim][kps][jps].B2c;
        *(pSnd++) = pPO->myFlx[dim][kps][jps].B3c;
#endif
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) *(pSnd++) = pPO->myFlx[dim][kps][jps].s[n];
#endif
        nFlx = NVAR;
      } else {  /*----- 2D or 3D problem -----*/

/* Conservative average in x2 of x1-fluxes */

        for (k=0; k<=(kpe-kps); k+=2) {
        for (j=0; j<=(jpe-jps); j+=2) {
          *(pSnd++) = pPO->myFlx[dim][k][j].d  + pPO->myFlx[dim][k][j+1].d;
          *(pSnd++) = pPO->myFlx[dim][k][j].M1 + pPO->myFlx[dim][k][j+1].M1;
          *(pSnd++) = pPO->myFlx[dim][k][j].M2 + pPO->myFlx[dim][k][j+1].M2;
          *(pSnd++) = pPO->myFlx[dim][k][j].M3 + pPO->myFlx[dim][k][j+1].M3;
#ifndef BAROTROPIC
          *(pSnd++) = pPO->myFlx[dim][k][j].E  + pPO->myFlx[dim][k][j+1].E;
#endif
#ifdef MHD
          *(pSnd++) = pPO->myFlx[dim][k][j].B1c + pPO->myFlx[dim][k][j+1].B1c;
          *(pSnd++) = pPO->myFlx[dim][k][j].B2c + pPO->myFlx[dim][k][j+1].B2c;
          *(pSnd++) = pPO->myFlx[dim][k][j].B3c + pPO->myFlx[dim][k][j+1].B3c;
#endif
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pSnd++) =
            pPO->myFlx[dim][k][j].s[n] + pPO->myFlx[dim][k][j+1].s[n];
#endif
        }}
        fact = 0.5;
        nFlx = ((jpe-jps+1)/2)*(NVAR);

/* Add conservative average in x3 of x1-fluxes */

        if (nDim == 3) {  /*----- 3D problem -----*/
          pSnd = (double*)&(send_bufRC[nd][start_addr+cnt]); /* restart ptr */
          for (k=0; k<=(kpe-kps); k+=2) {
          for (j=0; j<=(jpe-jps); j+=2) {
            *(pSnd++)+=pPO->myFlx[dim][k+1][j].d  +pPO->myFlx[dim][k+1][j+1].d;
            *(pSnd++)+=pPO->myFlx[dim][k+1][j].M1 +pPO->myFlx[dim][k+1][j+1].M1;
            *(pSnd++)+=pPO->myFlx[dim][k+1][j].M2 +pPO->myFlx[dim][k+1][j+1].M2;
            *(pSnd++)+=pPO->myFlx[dim][k+1][j].M3 +pPO->myFlx[dim][k+1][j+1].M3;
#ifndef BAROTROPIC
            *(pSnd++)+=pPO->myFlx[dim][k+1][j].E  +pPO->myFlx[dim][k+1][j+1].E;
#endif
#ifdef MHD
          *(pSnd++)+=pPO->myFlx[dim][k+1][j].B1c +pPO->myFlx[dim][k+1][j+1].B1c;
          *(pSnd++)+=pPO->myFlx[dim][k+1][j].B2c +pPO->myFlx[dim][k+1][j+1].B2c;
          *(pSnd++)+=pPO->myFlx[dim][k+1][j].B3c +pPO->myFlx[dim][k+1][j+1].B3c;
#endif
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) *(pSnd++) +=
              pPO->myFlx[dim][k+1][j].s[n] + pPO->myFlx[dim][k+1][j+1].s[n];
#endif
          }}
          fact = 0.25;
          nFlx = ((kpe-kps+1)*(jpe-jps+1)/4)*(NVAR);
        }

/* reset pointer to beginning of x1-fluxes and normalize averages */
        pSnd = (double*)&(send_bufRC[nd][(start_addr+cnt)]);  
        for (i=(start_addr+cnt); i<(start_addr+cnt+nFlx); i++) *(pSnd++) *=fact;
      }
      cnt += nFlx;

      }}

/*---------------- Restrict fluxes at x2-faces -------------------------------*/

      for (dim=2; dim<4; dim++){
      if (pPO->myFlx[dim] != NULL) {

      pSnd = (double*)&(send_bufRC[nd][(start_addr+cnt)]);

/* Conservative average in x1 of x2-fluxes */

      for (k=0; k<=(kpe-kps); k+=2) {
      for (i=0; i<=(ipe-ips); i+=2) {
        *(pSnd++) = pPO->myFlx[dim][k][i].d  + pPO->myFlx[dim][k][i+1].d;
        *(pSnd++) = pPO->myFlx[dim][k][i].M1 + pPO->myFlx[dim][k][i+1].M1;
        *(pSnd++) = pPO->myFlx[dim][k][i].M2 + pPO->myFlx[dim][k][i+1].M2;
        *(pSnd++) = pPO->myFlx[dim][k][i].M3 + pPO->myFlx[dim][k][i+1].M3;
#ifndef BAROTROPIC
        *(pSnd++) = pPO->myFlx[dim][k][i].E  + pPO->myFlx[dim][k][i+1].E;
#endif
#ifdef MHD
        *(pSnd++) = pPO->myFlx[dim][k][i].B1c + pPO->myFlx[dim][k][i+1].B1c;
        *(pSnd++) = pPO->myFlx[dim][k][i].B2c + pPO->myFlx[dim][k][i+1].B2c;
        *(pSnd++) = pPO->myFlx[dim][k][i].B3c + pPO->myFlx[dim][k][i+1].B3c;
#endif
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) *(pSnd++) =
          pPO->myFlx[dim][k][i].s[n] + pPO->myFlx[dim][k][i+1].s[n];
#endif
      }}
      fact = 0.5;
      nFlx = ((ipe-ips+1)/2)*(NVAR);

/* Add conservative average in x3 of x2-fluxes */

      if (nDim == 3) {  /*----- 3D problem -----*/
        pSnd = (double*)&(send_bufRC[nd][start_addr+cnt]); /* restart ptr */
        for (k=0; k<=(kpe-kps); k+=2) {
        for (i=0; i<=(ipe-ips); i+=2) {
          *(pSnd++) +=pPO->myFlx[dim][k+1][i].d  + pPO->myFlx[dim][k+1][i+1].d;
          *(pSnd++) +=pPO->myFlx[dim][k+1][i].M1 + pPO->myFlx[dim][k+1][i+1].M1;
          *(pSnd++) +=pPO->myFlx[dim][k+1][i].M2 + pPO->myFlx[dim][k+1][i+1].M2;
          *(pSnd++) +=pPO->myFlx[dim][k+1][i].M3 + pPO->myFlx[dim][k+1][i+1].M3;
#ifndef BAROTROPIC
          *(pSnd++) +=pPO->myFlx[dim][k+1][i].E  + pPO->myFlx[dim][k+1][i+1].E;
#endif
#ifdef MHD
          *(pSnd++)+= pPO->myFlx[dim][k+1][i].B1c+pPO->myFlx[dim][k+1][i+1].B1c;
          *(pSnd++)+= pPO->myFlx[dim][k+1][i].B2c+pPO->myFlx[dim][k+1][i+1].B2c;
          *(pSnd++)+= pPO->myFlx[dim][k+1][i].B3c+pPO->myFlx[dim][k+1][i+1].B3c;
#endif
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pSnd++) +=
            pPO->myFlx[dim][k+1][i].s[n] + pPO->myFlx[dim][k+1][i+1].s[n];
#endif
        }}
        fact = 0.25;
        nFlx = ((kpe-kps+1)*(ipe-ips+1)/4)*(NVAR);
      }

/* reset pointer to beginning of x2-fluxes and normalize averages */
      pSnd = (double*)&(send_bufRC[nd][(start_addr+cnt)]);  
      for (i=(start_addr+cnt); i<(start_addr+cnt+nFlx); i++) *(pSnd++) *= fact;
      cnt += nFlx;

      }}

/*---------------- Restrict fluxes at x3-faces -------------------------------*/

      for (dim=4; dim<6; dim++){
      if (pPO->myFlx[dim] != NULL) {

      pSnd = (double*)&(send_bufRC[nd][(start_addr+cnt)]);

/* Conservative average in x1 of x3-fluxes */

      for (j=0; j<=(jpe-jps); j+=2) {
      for (i=0; i<=(ipe-ips); i+=2) {
        *(pSnd++) = pPO->myFlx[dim][j][i].d  + pPO->myFlx[dim][j][i+1].d;
        *(pSnd++) = pPO->myFlx[dim][j][i].M1 + pPO->myFlx[dim][j][i+1].M1;
        *(pSnd++) = pPO->myFlx[dim][j][i].M2 + pPO->myFlx[dim][j][i+1].M2;
        *(pSnd++) = pPO->myFlx[dim][j][i].M3 + pPO->myFlx[dim][j][i+1].M3;
#ifndef BAROTROPIC
        *(pSnd++) = pPO->myFlx[dim][j][i].E  + pPO->myFlx[dim][j][i+1].E;
#endif
#ifdef MHD
        *(pSnd++) = pPO->myFlx[dim][j][i].B1c + pPO->myFlx[dim][j][i+1].B1c;
        *(pSnd++) = pPO->myFlx[dim][j][i].B2c + pPO->myFlx[dim][j][i+1].B2c;
        *(pSnd++) = pPO->myFlx[dim][j][i].B3c + pPO->myFlx[dim][j][i+1].B3c;
#endif
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) *(pSnd++) =
          pPO->myFlx[dim][j][i].s[n] + pPO->myFlx[dim][j][i+1].s[n];
#endif
      }}

/* Add conservative average in x2 of x3-fluxes */

      pSnd = (double*)&(send_bufRC[nd][(start_addr+cnt)]);
      for (j=0; j<=(jpe-jps); j+=2) {
      for (i=0; i<=(ipe-ips); i+=2) {
        *(pSnd++) +=pPO->myFlx[dim][j+1][i].d  + pPO->myFlx[dim][j+1][i+1].d;
        *(pSnd++) +=pPO->myFlx[dim][j+1][i].M1 + pPO->myFlx[dim][j+1][i+1].M1;
        *(pSnd++) +=pPO->myFlx[dim][j+1][i].M2 + pPO->myFlx[dim][j+1][i+1].M2;
        *(pSnd++) +=pPO->myFlx[dim][j+1][i].M3 + pPO->myFlx[dim][j+1][i+1].M3;
#ifndef BAROTROPIC
        *(pSnd++) +=pPO->myFlx[dim][j+1][i].E  + pPO->myFlx[dim][j+1][i+1].E;
#endif
#ifdef MHD
        *(pSnd++) +=pPO->myFlx[dim][j+1][i].B1c + pPO->myFlx[dim][j+1][i+1].B1c;
        *(pSnd++) +=pPO->myFlx[dim][j+1][i].B2c + pPO->myFlx[dim][j+1][i+1].B2c;
        *(pSnd++) +=pPO->myFlx[dim][j+1][i].B3c + pPO->myFlx[dim][j+1][i+1].B3c;
#endif
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) *(pSnd++) +=
          pPO->myFlx[dim][j+1][i].s[n] + pPO->myFlx[dim][j+1][i+1].s[n];
#endif
      }}
      fact = 0.25;
      nFlx = ((jpe-jps+1)*(ipe-ips+1)/4)*(NVAR);

/* reset pointer to beginning of x3-fluxes and normalize averages */
      pSnd = (double*)&(send_bufRC[nd][(start_addr+cnt)]);  
      for (i=(start_addr+cnt); i<(start_addr+cnt+nFlx); i++) *(pSnd++) *= fact;
      cnt += nFlx;

      }}

/*--- Step 3d. Restrict fluxes (EMFs) of face-centered fields ----------------*/

/*------------------ Restrict EMF at x1-faces --------------------------------*/
/* Only required for 2D or 3D problems.  Since EMF is a line integral, only
 * averaging along direction of EMF is required.  */

#ifdef MHD
      for (dim=0; dim<2; dim++){
        if (pPO->myEMF3[dim] != NULL) {

/* 2D problem -- Copy EMF3 */

          if (pG->Nx[2] == 1) {  
            for (k=0; k<=(kpe-kps)  ; k+=2) {
            for (j=0; j<=(jpe-jps)+1; j+=2) {
              *(pSnd++) = pPO->myEMF3[dim][k][j];
            }}

          } else {  

/* 3D problem -- Conservative averages of EMF3 and EMF2 */

            for (k=0; k<=(kpe-kps)  ; k+=2) {
            for (j=0; j<=(jpe-jps)+1; j+=2) {
              *(pSnd++) = 0.5*(pPO->myEMF3[dim][k][j]+pPO->myEMF3[dim][k+1][j]);
            }}

            for (k=0; k<=(kpe-kps)+1; k+=2) {
            for (j=0; j<=(jpe-jps)  ; j+=2) {
              *(pSnd++) = 0.5*(pPO->myEMF2[dim][k][j]+pPO->myEMF2[dim][k][j+1]);
            }}
          }
        }
      }

/*------------------- Restrict EMF at x2-faces -------------------------------*/

      for (dim=2; dim<4; dim++){
        if (pPO->myEMF3[dim] != NULL) {

/* 2D problem --  Copy EMF3 */

          if (pG->Nx[2] == 1) {
            for (k=0; k<=(kpe-kps)  ; k+=2) {
            for (i=0; i<=(ipe-ips)+1; i+=2) {
              *(pSnd++) = pPO->myEMF3[dim][k][i];
            }}

          } else {

/* 3D problem -- Conservative averages of EMF3 and EMF1 */

            for (k=0; k<=(kpe-kps)  ; k+=2) {
            for (i=0; i<=(ipe-ips)+1; i+=2) {
              *(pSnd++) = 0.5*(pPO->myEMF3[dim][k][i]+pPO->myEMF3[dim][k+1][i]);
            }}

            for (k=0; k<=(kpe-kps)+1; k+=2) {
            for (i=0; i<=(ipe-ips)  ; i+=2) {
              *(pSnd++) = 0.5*(pPO->myEMF1[dim][k][i]+pPO->myEMF1[dim][k][i+1]);
            }}
          }
        }
      }

/*------------------- Restrict EMF at x3-faces -------------------------------*/
/* Must be a 3D problem */

      for (dim=4; dim<6; dim++){
        if (pPO->myEMF1[dim] != NULL) {

/*----- 3D problem ----- Conservative averages of EMF1 and EMF2 */

          for (j=0; j<=(jpe-jps)+1; j+=2) {
          for (i=0; i<=(ipe-ips)  ; i+=2) {
            *(pSnd++) = 0.5*(pPO->myEMF1[dim][j][i] + pPO->myEMF1[dim][j][i+1]);
          }}

          for (j=0; j<=(jpe-jps)  ; j+=2) {
          for (i=0; i<=(ipe-ips)+1; i+=2) {
            *(pSnd++) = 0.5*(pPO->myEMF2[dim][j][i] + pPO->myEMF2[dim][j+1][i]);
          }}
        }
      }

#endif

#ifdef MPI_PARALLEL
/*--- Step 3e. Send rectricted soln and fluxes -------------------------------*/
/* non-blocking send with MPI, using Domain number as tag.  */

      if (npg >= pG->NmyPGrid){
        mIndex = npg - pG->NmyPGrid;
        ierr = MPI_Isend(&(send_bufRC[nd][start_addr]), pG->PGrid[npg].nWordsRC,
          MPI_DOUBLE, pG->PGrid[npg].ID, nd, pM->Domain[nl][nd].Comm_Parent,
          &(send_rq[nd][mIndex]));
      }
#endif /* MPI_PARALLEL */

      start_addr += pG->PGrid[npg].nWordsRC;

    }  /* end loop over parent grids */
  }} /* end loop over Domains per level */

#ifdef MPI_PARALLEL
/*--- Step 4. Check non-blocking sends completed. ----------------------------*/
/* For MPI jobs, wait for all non-blocking sends in Step 3e to complete.  This
 * is more efficient if there are multiple messages per Grid. */

  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {
      pG=pM->Domain[nl][nd].Grid;

      if (pG->NPGrid > pG->NmyPGrid) {
        mCount = pG->NPGrid - pG->NmyPGrid;
        ierr = MPI_Waitall(mCount, send_rq[nd], MPI_STATUS_IGNORE);
      }
    }
  }
#endif /* MPI_PARALLEL */

  } /* end loop over levels */
}

/*============================================================================*/
/*----------------------------------------------------------------------------*/
/*! \fn void Prolongate(MeshS *pM)
 *  \brief Sets BC on fine Grid by prolongation (interpolation) of
 *     coarse Grid solution into fine grid ghost zones */
void Prolongate(MeshS *pM)
{
  GridS *pG;
  int nDim,nl,nd,ncg,dim,npg,rbufN,id,l,m,n,mend,nend;
  int i,ii,ics,ice,ips,ipe,igzs,igze;
  int j,jj,jcs,jce,jps,jpe,jgzs,jgze;
  int k,kk,kcs,kce,kps,kpe,kgzs,kgze;
  int ngz1,ngz2,ngz3;
  double *pRcv,*pSnd;
  GridOvrlpS *pCO, *pPO;
  ConsS ProlongedC[2][2][2];
#if (NSCALARS > 0)
  int ns;
#endif
#ifdef MHD
  Real3Vect BGZ[3][3][3], ProlongedF[3][3][3];
#endif
#ifdef MPI_PARALLEL
  int ierr,mAddress,mIndex,mCount;
#endif

/* number of dimensions in Grid. */
  nDim=1;
  for (dim=1; dim<3; dim++) if (pM->Nx[dim]>1) nDim++;

/* Loop over all levels, starting at root level */

  for (nl=0; nl<(pM->NLevels); nl++){

#ifdef MPI_PARALLEL
/* Post non-blocking receives at level nl+1 for data from parent Grids at this
 * level (nl). This data is sent in Step 1 below,
 * and will be read in Step 2 during the next iteration of nl */

  if (nl<(pM->NLevels)-1) {
    for (nd=0; nd<(pM->DomainsPerLevel[nl+1]); nd++){
      if (pM->Domain[nl+1][nd].Grid != NULL) {
        pG=pM->Domain[nl+1][nd].Grid;

/* Recv buffer is addressed from PGrid[0].nWordsP if NmyPGrid>0 since data
 * from parent Grid on same processor begins at 0 (see Step 4 below).  First
 * index alternates between 0 and 1 for even/odd values of nl, since if
 * there are Grids on multiple levels there may be 2 receives posted at once */

        mAddress = 0;
        rbufN = ((nl+1) % 2);
        if (pG->NmyPGrid > 0) mAddress = pG->PGrid[0].nWordsP;

        for (npg=(pG->NmyPGrid); npg<(pG->NPGrid); npg++){
          mIndex = npg - pG->NmyPGrid;
          ierr = MPI_Irecv(&(recv_bufP[rbufN][nd][mAddress]),
            pG->PGrid[npg].nWordsP, MPI_DOUBLE, pG->PGrid[npg].ID,
            pG->PGrid[npg].DomN, pM->Domain[nl+1][nd].Comm_Parent,
            &(recv_rq[nl+1][nd][mIndex]));
          mAddress += pG->PGrid[npg].nWordsP;
        }

      }
    }
  }
#endif /* MPI_PARALLEL */

/*=== Step 1. Send step ======================================================*/
/* Loop over all Domains, and send ghost zones to all child Grids. */

  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){

  if (pM->Domain[nl][nd].Grid != NULL) { /* there is a Grid on this processor */
    pG=pM->Domain[nl][nd].Grid;
    for(i=0; i<maxND; i++) start_addrP[i] = 0;

    for (ncg=0; ncg<(pG->NCGrid); ncg++){
      pCO=(GridOvrlpS*)&(pG->CGrid[ncg]);    /* ptr to child Grid overlap */

/* index send_buf with DomN of child, since could be multiple child Domains on
 * same processor.  Start address must be different for each DomN */
      pSnd = (double*)&(send_bufP[pCO->DomN][start_addrP[pCO->DomN]]); 

      for (dim=0; dim<(2*nDim); dim++){
        if (pCO->myFlx[dim] != NULL) {

/* Get coordinates ON THIS GRID of zones that overlap child Grid ghost zones */

          ics = pCO->ijks[0] - (nghost/2) - 1;
          ice = pCO->ijke[0] + (nghost/2) + 1;
          if (pG->Nx[1] > 1) {
            jcs = pCO->ijks[1] - (nghost/2) - 1;
            jce = pCO->ijke[1] + (nghost/2) + 1;
          } else {
            jcs = pCO->ijks[1];
            jce = pCO->ijke[1];
          }
          if (pG->Nx[2] > 1) {
            kcs = pCO->ijks[2] - (nghost/2) - 1;
            kce = pCO->ijke[2] + (nghost/2) + 1;
          } else {
            kcs = pCO->ijks[2];
            kce = pCO->ijke[2];
          }
          if (dim == 0) ice = pCO->ijks[0];
          if (dim == 1) ics = pCO->ijke[0];
          if (dim == 2) jce = pCO->ijks[1];
          if (dim == 3) jcs = pCO->ijke[1];
          if (dim == 4) kce = pCO->ijks[2];
          if (dim == 5) kcs = pCO->ijke[2];

/*--- Step 1a. ---------------------------------------------------------------*/
/* Load send buffer with values in zones that overlap child ghost zones */

          for (k=kcs; k<=kce; k++) {
          for (j=jcs; j<=jce; j++) {
          for (i=ics; i<=ice; i++) {
            *(pSnd++) = pG->U[k][j][i].d;
            *(pSnd++) = pG->U[k][j][i].M1;
            *(pSnd++) = pG->U[k][j][i].M2;
            *(pSnd++) = pG->U[k][j][i].M3;
#ifndef BAROTROPIC
            *(pSnd++) = pG->U[k][j][i].E;
#endif
#ifdef MHD
            *(pSnd++) = pG->U[k][j][i].B1c;
            *(pSnd++) = pG->U[k][j][i].B2c;
            *(pSnd++) = pG->U[k][j][i].B3c;
            *(pSnd++) = pG->B1i[k][j][i];
            *(pSnd++) = pG->B2i[k][j][i];
            *(pSnd++) = pG->B3i[k][j][i];
#endif
#if (NSCALARS > 0)
            for (ns=0; ns<NSCALARS; ns++) {
               *(pSnd++) = pG->U[k][j][i].s[ns];
            }
#endif
          }}}
        }
      }

/*--- Step 1b. ---------------------------------------------------------------*/
/* non-blocking send of data to child, using Domain number as tag. */
#ifdef MPI_PARALLEL
      if (ncg >= pG->NmyCGrid) {
        mIndex = ncg - pG->NmyCGrid;
        ierr = MPI_Isend(&(send_bufP[pCO->DomN][start_addrP[pCO->DomN]]),
          pG->CGrid[ncg].nWordsP, MPI_DOUBLE, pG->CGrid[ncg].ID, nd,
          pM->Domain[nl][nd].Comm_Children, &(send_rq[nd][mIndex]));
      }
#endif /* MPI_PARALLEL */

      start_addrP[pCO->DomN] += pG->CGrid[ncg].nWordsP;

    } /* end loop over child grids */
  }} /* end loop over Domains */

/*=== Step 2. Get step =======================================================*/
/* Loop over all Domains, get data sent by parent Grids, and prolongate solution
 * into ghost zones. */


  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){

  if (pM->Domain[nl][nd].Grid != NULL) { /* there is a Grid on this processor */
    pG=pM->Domain[nl][nd].Grid;          /* set pointer to Grid */
    rbufN = (nl % 2);

    for (npg=0; npg<(pG->NPGrid); npg++){

/* If parent Grid is on this processor, data is at start of recv buffer */

      if (npg < pG->NmyPGrid) {
        pPO = (GridOvrlpS*)&(pG->PGrid[npg]);
        pRcv = (double*)&(recv_bufP[rbufN][nd][0]);
      } else {

#ifdef MPI_PARALLEL
/* Check non-blocking receives posted above for data in ghost zone from parent
 * Grids, sent in Step 1.  Accept messages in any order. */

        mCount = pG->NPGrid - pG->NmyPGrid;
        ierr = MPI_Waitany(mCount,recv_rq[nl][nd],&mIndex,MPI_STATUS_IGNORE);
        if(mIndex == MPI_UNDEFINED){
          ath_error("[Prolong]: Invalid request index nl=%i nd=%i\n",nl,nd);
        }

/* Recv buffer is addressed from PGrid[0].nWordsP for first MPI message
 * if NmyPGrid>0 */

        mAddress = 0;
        mIndex += pG->NmyPGrid;
        for (i=0; i<mIndex; i++) mAddress += pG->PGrid[i].nWordsP;
        pPO = (GridOvrlpS*)&(pG->PGrid[mIndex]); 
        pRcv = (double*)&(recv_bufP[rbufN][nd][mAddress]);

#else
/* If not MPI_PARALLEL, and parent Grid not on this processor, then error */

        ath_error("[Prolong]: no Parent Grid on Domain[%d][%d]\n",nl,nd);
#endif /* MPI_PARALLEL */
      }

/*=== Step 3. Set ghost zones ================================================*/
/* Loop over 6 boundaries, set ghost zones */

      for (dim=0; dim<(2*nDim); dim++){
        if (pPO->myFlx[dim] != NULL) {

/*--- Steps 3a.  Set GZ and BFld arrays --------------------------------------*/
/* Compute size of array containing ghost zones from parent Grid.  Set
 * starting and ending indices of GZ array. */

          if (dim == 0 || dim == 1) {
            ngz1 = (nghost/2) + 2;
            id = 0;
          } else {
            ngz1 = (pPO->ijke[0] - pPO->ijks[0] + 1)/2 + nghost + 2;
          }

          if (dim == 2 || dim == 3) {
            ngz2 = (nghost/2) + 2;
            id = 1;
          } else {
            ngz2 = (pPO->ijke[1] - pPO->ijks[1] + 1)/2 + nghost + 2;
          }

          if (dim == 4 || dim == 5) {
            ngz3 = (nghost/2) + 2;
            id = 2;
          } else {
            ngz3 = (pPO->ijke[2] - pPO->ijks[2] + 1)/2 + nghost + 2;
          }

          igzs = 0;
          igze = ngz1-1;
          if (pG->Nx[1] > 1) {
            jgzs = 0;
            jgze = ngz2-1;
            mend = 1;
          } else {
            ngz2 = 1;
            jgzs = 1;
            jgze = 1;
            mend = 0;
          }
          if (pG->Nx[2] > 1) {
            kgzs = 0;
            kgze = ngz3-1;
            nend = 1;
          } else {
            ngz3 = 1;
            kgzs = 1;
            kgze = 1;
            nend = 0;
          }

/* Load GZ array with values in receive buffer */

          for (k=kgzs; k<=kgze; k++) {
          for (j=jgzs; j<=jgze; j++) {
          for (i=igzs; i<=igze; i++) {
            GZ[id][k][j][i].d  = *(pRcv++);
            GZ[id][k][j][i].M1 = *(pRcv++);
            GZ[id][k][j][i].M2 = *(pRcv++);
            GZ[id][k][j][i].M3 = *(pRcv++);
#ifndef BAROTROPIC
            GZ[id][k][j][i].E = *(pRcv++);
#endif
#ifdef MHD
            GZ[id][k][j][i].B1c = *(pRcv++);
            GZ[id][k][j][i].B2c = *(pRcv++);
            GZ[id][k][j][i].B3c = *(pRcv++);
            BFld[id][k][j][i].x = *(pRcv++);
            BFld[id][k][j][i].y = *(pRcv++);
            BFld[id][k][j][i].z = *(pRcv++);
#endif
#if (NSCALARS > 0)
            for (ns=0; ns<NSCALARS; ns++) {
              GZ[id][k][j][i].s[ns] = *(pRcv++);
            }
#endif
          }}}

/* Set BC on GZ array in 1D; and on GZ and BFld arrays in 2D */

          if (nDim == 1) {
            for (i=igzs; i<=igze; i++) {
              GZ[id][1][0][i] = GZ[id][1][1][i];
              GZ[id][1][2][i] = GZ[id][1][1][i];
              GZ[id][0][1][i] = GZ[id][1][1][i];
              GZ[id][2][1][i] = GZ[id][1][1][i];
            }
          }

          if (nDim == 2) {
            for (j=jgzs; j<=jgze; j++) {
            for (i=igzs; i<=igze; i++) {
              GZ[id][0][j][i] = GZ[id][1][j][i];
              GZ[id][2][j][i] = GZ[id][1][j][i];
            }}
#ifdef MHD
            for (j=jgzs; j<=jgze; j++) {
            for (i=igzs; i<=igze; i++) {
              BFld[id][0][j][i] = BFld[id][1][j][i];
              BFld[id][2][j][i] = BFld[id][1][j][i];
            }}
#endif /* MHD */
          }

/*--- Steps 3b.  Prolongate cell-centered values -----------------------------*/
/* Get coordinates ON THIS GRID of ghost zones that overlap parent Grid */

          ips = pPO->ijks[0] - nghost;
          ipe = pPO->ijke[0] + nghost;
          if (pG->Nx[1] > 1) {
            jps = pPO->ijks[1] - nghost;
            jpe = pPO->ijke[1] + nghost;
          } else {
            jps = pPO->ijks[1];
            jpe = pPO->ijke[1];
          }
          if (pG->Nx[2] > 1) {
            kps = pPO->ijks[2] - nghost;
            kpe = pPO->ijke[2] + nghost;
          } else {
            kps = pPO->ijks[2];
            kpe = pPO->ijke[2];
          }
          if (dim == 0) {ipe = pPO->ijks[0] - 1;}
          if (dim == 1) {ips = pPO->ijke[0] + 1;}
          if (dim == 2) {jpe = pPO->ijks[1] - 1;}
          if (dim == 3) {jps = pPO->ijke[1] + 1;}
          if (dim == 4) {kpe = pPO->ijks[2] - 1;}
          if (dim == 5) {kps = pPO->ijke[2] + 1;}

/* Prolongate these values in ghost zones */

          for (k=kps, kk=1; k<=kpe; k+=2, kk++) {
          for (j=jps, jj=1; j<=jpe; j+=2, jj++) {
          for (i=ips, ii=1; i<=ipe; i+=2, ii++) {

            ProCon(GZ[id][kk][jj][ii-1],GZ[id][kk][jj][ii],GZ[id][kk][jj][ii+1],
                   GZ[id][kk][jj-1][ii],                   GZ[id][kk][jj+1][ii],
                   GZ[id][kk-1][jj][ii],                   GZ[id][kk+1][jj][ii],
                   ProlongedC);

/* 1D/2D/3D problem, set solution prolongated in x1 */

            for (n=0; n<=nend; n++) {
            for (m=0; m<=mend; m++) {
            for (l=0; l<=1; l++) {
              pG->U[k+n][j+m][i+l].d  = ProlongedC[n][m][l].d;
              pG->U[k+n][j+m][i+l].M1 = ProlongedC[n][m][l].M1;
              pG->U[k+n][j+m][i+l].M2 = ProlongedC[n][m][l].M2;
              pG->U[k+n][j+m][i+l].M3 = ProlongedC[n][m][l].M3;
#ifndef BAROTROPIC
              pG->U[k+n][j+m][i+l].E  = ProlongedC[n][m][l].E;
#endif
#ifdef MHD
              pG->U[k+n][j+m][i+l].B1c = ProlongedC[n][m][l].B1c;
              pG->U[k+n][j+m][i+l].B2c = ProlongedC[n][m][l].B2c;
              pG->U[k+n][j+m][i+l].B3c = ProlongedC[n][m][l].B3c;
#endif
#if (NSCALARS > 0)
              for (ns=0; ns<NSCALARS; ns++) 
                pG->U[k+n][j+m][i+l].s[ns] = ProlongedC[n][m][l].s[ns];
#endif
            }}}

#ifdef MHD
/*--- Steps 3c.  Prolongate face-centered B ----------------------------------*/
/* Set prolonged face-centered B fields for 1D (trivial case)  */

            if (nDim == 1) {
              for (l=0; l<=1; l++) {
                pG->B1i[k][j][i+l] = pG->U[k][j][i+l].B1c;
                pG->B2i[k][j][i+l] = pG->U[k][j][i+l].B2c;
                pG->B3i[k][j][i+l] = pG->U[k][j][i+l].B3c;
              }
            } else {
              for (n=0; n<3; n++) {
              for (m=0; m<3; m++) {
              for (l=0; l<3; l++) {
                ProlongedF[n][m][l].x = 0.0;
                ProlongedF[n][m][l].y = 0.0;
                ProlongedF[n][m][l].z = 0.0;
              }}}
            }

/* Load B-field ghost zone array with values read from Rcv buffer in 2D/3D */

            if (nDim == 2 || nDim ==3) {
              for (n=0; n<3; n++) {
              for (m=0; m<3; m++) {
              for (l=0; l<3; l++) {
                BGZ[n][m][l].x = BFld[id][kk+(n-1)][jj+(m-1)][ii+(l-1)].x;
                BGZ[n][m][l].y = BFld[id][kk+(n-1)][jj+(m-1)][ii+(l-1)].y;
                BGZ[n][m][l].z = BFld[id][kk+(n-1)][jj+(m-1)][ii+(l-1)].z;
              }}}

/* If edge of cell touches fine/coarse boundary, use fine grid fields for the
 * normal component at interface. ProFld will not overwrite these values.  If
 * the start/end of boundary is between MPI Grids (pPO->myFlx[]==NULL), then use
 * fine grid fields in corners as well. */

/* inner x1 boundary */
              if ((dim == 0) &&
                  (i == (ipe-1)) &&
                  ((j >= (jps+nghost)) || (pPO->myFlx[2]==NULL)) &&
                  ((j <  (jpe-nghost)) || (pPO->myFlx[3]==NULL)) ){
                ProlongedF[0][0][2].x = pG->B1i[k][j  ][i+2];
                ProlongedF[0][1][2].x = pG->B1i[k][j+1][i+2];
                ProlongedF[1][0][2].x = pG->B1i[k][j  ][i+2];
                ProlongedF[1][1][2].x = pG->B1i[k][j+1][i+2];
                if ((nDim == 3) &&
                    ((k >= (kps+nghost)) || (pPO->myFlx[4]==NULL)) &&
                    ((k <  (kpe-nghost)) || (pPO->myFlx[5]==NULL)) ){
                  ProlongedF[1][0][2].x = pG->B1i[k+1][j  ][i+2];
                  ProlongedF[1][1][2].x = pG->B1i[k+1][j+1][i+2];
                }
              }

/* outer x1 boundary */
              if ((dim == 1) &&
                  (i == ips) &&
                  ((j >= (jps+nghost)) || (pPO->myFlx[2]==NULL)) &&
                  ((j <  (jpe-nghost)) || (pPO->myFlx[3]==NULL)) ){
                ProlongedF[0][0][0].x = pG->B1i[k][j  ][i];
                ProlongedF[0][1][0].x = pG->B1i[k][j+1][i];
                ProlongedF[1][0][0].x = pG->B1i[k][j  ][i];
                ProlongedF[1][1][0].x = pG->B1i[k][j+1][i];
                if ((nDim == 3) &&
                    ((k >= (kps+nghost)) || (pPO->myFlx[4]==NULL)) &&
                    ((k <  (kpe-nghost)) || (pPO->myFlx[5]==NULL)) ){
                  ProlongedF[1][0][0].x = pG->B1i[k+1][j  ][i];
                  ProlongedF[1][1][0].x = pG->B1i[k+1][j+1][i];
                }
              }

/* inner x2 boundary */
              if ((dim == 2) &&
                  (j == (jpe-1)) &&
                  ((i >= (ips+nghost)) || (pPO->myFlx[0]==NULL)) &&
                  ((i <  (ipe-nghost)) || (pPO->myFlx[1]==NULL)) ){
                ProlongedF[0][2][0].y = pG->B2i[k][j+2][i  ];
                ProlongedF[0][2][1].y = pG->B2i[k][j+2][i+1];
                ProlongedF[1][2][0].y = pG->B2i[k][j+2][i  ];
                ProlongedF[1][2][1].y = pG->B2i[k][j+2][i+1];
                if ((nDim == 3) &&
                    ((k >= (kps+nghost)) || (pPO->myFlx[4]==NULL)) &&
                    ((k <  (kpe-nghost)) || (pPO->myFlx[5]==NULL)) ){
                  ProlongedF[1][2][0].y = pG->B2i[k+1][j+2][i  ];
                  ProlongedF[1][2][1].y = pG->B2i[k+1][j+2][i+1];
                }
              }

/* outer x2 boundary */
              if ((dim == 3) &&
                  (j == jps) &&
                  ((i >= (ips+nghost)) || (pPO->myFlx[0]==NULL)) &&
                  ((i <  (ipe-nghost)) || (pPO->myFlx[1]==NULL)) ){
                ProlongedF[0][0][0].y = pG->B2i[k][j][i  ];
                ProlongedF[0][0][1].y = pG->B2i[k][j][i+1];
                ProlongedF[1][0][0].y = pG->B2i[k][j][i  ];
                ProlongedF[1][0][1].y = pG->B2i[k][j][i+1];
                if ((nDim == 3) &&
                    ((k >= (kps+nghost)) || (pPO->myFlx[4]==NULL)) &&
                    ((k <  (kpe-nghost)) || (pPO->myFlx[5]==NULL)) ){
                  ProlongedF[1][0][0].y = pG->B2i[k+1][j][i  ];
                  ProlongedF[1][0][1].y = pG->B2i[k+1][j][i+1];
                }
              }

/* inner x3 boundary */
              if ((dim == 4) &&
                  (k == (kpe-1)) &&
                  ((i >= (ips+nghost)) || (pPO->myFlx[0]==NULL)) &&
                  ((i <  (ipe-nghost)) || (pPO->myFlx[1]==NULL)) &&
                  ((j >= (jps+nghost)) || (pPO->myFlx[2]==NULL)) &&
                  ((j <  (jpe-nghost)) || (pPO->myFlx[3]==NULL)) ){
                ProlongedF[2][0][0].z = pG->B3i[k+2][j  ][i  ];
                ProlongedF[2][0][1].z = pG->B3i[k+2][j  ][i+1];
                ProlongedF[2][1][0].z = pG->B3i[k+2][j+1][i  ];
                ProlongedF[2][1][1].z = pG->B3i[k+2][j+1][i+1];
              }

/* outer x3 boundary */
              if ((dim == 5) && 
                  (k == kps) &&
                  ((i >= (ips+nghost)) || (pPO->myFlx[0]==NULL)) &&
                  ((i <  (ipe-nghost)) || (pPO->myFlx[1]==NULL)) &&
                  ((j >= (jps+nghost)) || (pPO->myFlx[2]==NULL)) &&
                  ((j <  (jpe-nghost)) || (pPO->myFlx[3]==NULL)) ){
                ProlongedF[0][0][0].z = pG->B3i[k][j  ][i  ];
                ProlongedF[0][0][1].z = pG->B3i[k][j  ][i+1];
                ProlongedF[0][1][0].z = pG->B3i[k][j+1][i  ];
                ProlongedF[0][1][1].z = pG->B3i[k][j+1][i+1];
              }

              ProFld(BGZ, ProlongedF, pG->dx1, pG->dx2, pG->dx3);

              for (n=0; n<=nend; n++) {
              for (m=0; m<=mend; m++) {
              for (l=0; l<=1; l++) {
                if (dim != 1 || (i+l) != ips)
                  pG->B1i[k+n][j+m][i+l] = ProlongedF[n][m][l].x;
                if (dim != 3 || (j+m) != jps)
                  pG->B2i[k+n][j+m][i+l] = ProlongedF[n][m][l].y;
                if (dim != 5 || (k+n) != kps)
                  pG->B3i[k+n][j+m][i+l] = ProlongedF[n][m][l].z;

                pG->U[k+n][j+m][i+l].B1c = 
                  0.5*(ProlongedF[n][m][l].x + ProlongedF[n][m][l+1].x);
                pG->U[k+n][j+m][i+l].B2c = 
                  0.5*(ProlongedF[n][m][l].y + ProlongedF[n][m+1][l].y);
                pG->U[k+n][j+m][i+l].B3c = 
                  0.5*(ProlongedF[n][m][l].z + ProlongedF[n+1][m][l].z);
              }}}
            }

#endif /* MHD */
          }}}

        }
      } /* end loop over dims */

    } /* end loop over parent grids */
  }} /* end loop over Domains */

/*=== Step 4. Clear send_bufP ================================================*/
/* For child/parent Grids on same processor, copy send_bufP into recv_bufP to
 * prevent "send" in Step 1 above from over-writing data in buffer on the next
 * iteration of the loop over levels (for nl=nl+1). */

  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
    if (pM->Domain[nl][nd].Grid != NULL) { 
      pG=pM->Domain[nl][nd].Grid; 
      rbufN = ((nl+1) % 2);

/* For each Domain nd, the data for child grid on same processor must be in
 * first element of CGrid array */ 

      for (ncg=0; ncg<(pG->NmyCGrid); ncg++){
        pCO=(GridOvrlpS*)&(pG->CGrid[ncg]);    /* ptr to child Grid overlap */

        for (i=0; i<pCO->nWordsP; i++) {
          recv_bufP[rbufN][pCO->DomN][i]=send_bufP[pCO->DomN][i];
        }
      }
    }
  }

#ifdef MPI_PARALLEL
/* For MPI jobs, wait for all non-blocking sends in Step 1 to complete */

  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {
      pG=pM->Domain[nl][nd].Grid;

      if (pG->NCGrid > pG->NmyCGrid) {
        mCount = pG->NCGrid - pG->NmyCGrid;
        ierr = MPI_Waitall(mCount, send_rq[nd], MPI_STATUS_IGNORE);
      }
    }
  }
#endif /* MPI_PARALLEL */

  } /* end loop over levels */
}

/*============================================================================*/
/*----------------------------------------------------------------------------*/
/*! \fn void SMR_init(MeshS *pM)
 *  \brief Allocates memory for send/receive buffers
 */

void SMR_init(MeshS *pM)
{
  int nl,nd,sendRC,recvRC,sendP,recvP,npg,ncg;
  int max_sendRC=1,max_recvRC=1,max_sendP=1,max_recvP=1;
  int max1=0,max2=0,max3=0,maxCG=1;
#ifdef MHD
  int ngh1;
#endif
  GridS *pG;
  
  maxND=1;
  for (nl=0; nl<(pM->NLevels); nl++) maxND=MAX(maxND,pM->DomainsPerLevel[nl]);
  if((start_addrP = (int*)calloc_1d_array(maxND,sizeof(int))) == NULL)
    ath_error("[SMR_init]:Failed to allocate start_addrP\n");

/* Loop over all parent Grids of Grids on this processor to find maximum total
 * number of words communicated */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      sendRC=0;
      recvRC=0;
      sendP =0;
      recvP =0;

      if (pM->Domain[nl][nd].Grid != NULL) { /* there is a Grid on this proc */
        pG=pM->Domain[nl][nd].Grid;          /* set pointer to Grid */

        for (npg=0; npg<pG->NPGrid; npg++){
          sendRC += pG->PGrid[npg].nWordsRC;
          recvP  += pG->PGrid[npg].nWordsP;
        }
        for (ncg=0; ncg<pG->NCGrid; ncg++){
          recvRC += pG->CGrid[ncg].nWordsRC;
          sendP  += pG->CGrid[ncg].nWordsP;
        }

        max_sendRC = MAX(max_sendRC,sendRC);
        max_recvRC = MAX(max_recvRC,recvRC);
        max_sendP  = MAX(max_sendP ,sendP );
        max_recvP  = MAX(max_recvP ,recvP );
        max1 = MAX(max1,(pG->Nx[0]+1));
        max2 = MAX(max2,(pG->Nx[1]+1));
        max3 = MAX(max3,(pG->Nx[2]+1));
        maxCG = MAX(maxCG,pG->NCGrid);
      }
    }
  }

/* Allocate memory for send/receive buffers and EMFs used in RestrictCorrect */

  if((send_bufRC =
    (double**)calloc_2d_array(maxND,max_sendRC,sizeof(double))) == NULL)
    ath_error("[SMR_init]:Failed to allocate send_bufRC\n");

#ifdef MPI_PARALLEL
  if((recv_bufRC =
    (double***)calloc_3d_array(2,maxND,max_recvRC,sizeof(double))) == NULL)
    ath_error("[SMR_init]: Failed to allocate recv_bufRC\n");
  if((recv_rq = (MPI_Request***)
    calloc_3d_array(pM->NLevels,maxND,maxCG,sizeof(MPI_Request))) == NULL)
    ath_error("[SMR_init]: Failed to allocate recv MPI_Request array\n");
  if((send_rq = (MPI_Request**)
    calloc_2d_array(maxND,maxCG,sizeof(MPI_Request))) == NULL)
    ath_error("[SMR_init]: Failed to allocate send MPI_Request array\n");
#endif /* MPI_PARALLEL */

#ifdef MHD
  if((SMRemf1 =
    (Real**)calloc_2d_array(MAX(max2,max3),max1,sizeof(Real))) == NULL)
    ath_error("[smr_init]Failed to calloc_2d_array for SMRemf1\n");;
  if((SMRemf2 =
    (Real**)calloc_2d_array(MAX(max2,max3),MAX(max1,max2),sizeof(Real))) ==NULL)
    ath_error("[smr_init]Failed to calloc_2d_array for SMRemf2\n");;
  if((SMRemf3 =
    (Real**)calloc_2d_array(max3,MAX(max1,max2),sizeof(Real))) == NULL)
    ath_error("[smr_init]Failed to calloc_2d_array for SMRemf3\n");;
#endif /* MHD */

/* Allocate memory for send/receive buffers and GZ arrays used in Prolongate */

  if((send_bufP =
    (double**)calloc_2d_array(maxND,max_sendP,sizeof(double))) == NULL)
    ath_error("[SMR_init]:Failed to allocate send_bufP\n");

  if((recv_bufP =
    (double***)calloc_3d_array(2,maxND,max_recvP,sizeof(double))) == NULL)
    ath_error("[SMR_init]: Failed to allocate recv_bufP\n");

  max1 += 2*nghost;
  max2 += 2*nghost;
  max3 += 2*nghost;

  if((GZ[0]=(ConsS***)calloc_3d_array(max3,max2,nghost,sizeof(ConsS)))
    ==NULL) ath_error("[SMR_init]:Failed to allocate GZ[0]C\n");
  if((GZ[1]=(ConsS***)calloc_3d_array(max3,nghost,max1,sizeof(ConsS)))
    ==NULL) ath_error("[SMR_init]:Failed to allocate GZ[1]C\n");
  if((GZ[2]=(ConsS***)calloc_3d_array(nghost,max2,max1,sizeof(ConsS)))
    ==NULL) ath_error("[SMR_init]:Failed to allocate GZ[2]C\n");
#ifdef MHD
  ngh1 = nghost + 1;
  if((BFld[0]=(Real3Vect***)calloc_3d_array(max3,max2,ngh1,sizeof(Real3Vect)))
    ==NULL) ath_error("[SMR_init]:Failed to allocate BFld[0]C\n");
  if((BFld[1]=(Real3Vect***)calloc_3d_array(max3,ngh1,max1,sizeof(Real3Vect)))
    ==NULL) ath_error("[SMR_init]:Failed to allocate BFld[1]C\n");
  if((BFld[2]=(Real3Vect***)calloc_3d_array(ngh1,max2,max1,sizeof(Real3Vect)))
    ==NULL) ath_error("[SMR_init]:Failed to allocate BFld[2]C\n");
#endif /* MHD */

  return;
}
/*=========================== PRIVATE FUNCTIONS ==============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void ProCon(const ConsS Uim1,const ConsS Ui,  const ConsS Uip1,
 *            const ConsS Ujm1,const ConsS Ujp1,
 *            const ConsS Ukm1,const ConsS Ukp1, ConsS PCon[][2][2])
 *  \brief Prolongates conserved variables in a 2x2x2 cube.
 */

void ProCon(const ConsS Uim1,const ConsS Ui,  const ConsS Uip1,
            const ConsS Ujm1,const ConsS Ujp1,
            const ConsS Ukm1,const ConsS Ukp1, ConsS PCon[][2][2])
{
  int i,j,k;
  Real dq1,dq2,dq3,Pim1,Pi,Pip1,Pjm1,Pjp1,Pkm1,Pkp1;
#ifdef SPECIAL_RELATIVITY
  PrimS W;
  Real Vsq;
  int fail,dfail,Pfail,Vfail;
#endif
#if (NSCALARS > 0)
  int n;
#endif

/* First order prolongation -- just copy values */
#ifdef FIRST_ORDER

  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].d  = Ui.d;
    PCon[k][j][i].M1 = Ui.M1;
    PCon[k][j][i].M2 = Ui.M2;
    PCon[k][j][i].M3 = Ui.M3;
#ifndef BAROTROPIC
    PCon[k][j][i].E  = Ui.E;
#endif /* BAROTROPIC */
#ifdef MHD
    PCon[k][j][i].B1c = Ui.B1c;
    PCon[k][j][i].B2c = Ui.B2c;
    PCon[k][j][i].B3c = Ui.B3c;
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) PCon[k][j][i].s[n] = Ui.s[n];
#endif /* NSCALARS */
  }}}

/* second order prolongation -- apply limited slope reconstruction */
#else /* SECOND_ORDER or THIRD_ORDER */

/* density */
  dq1 = mcd_slope(Uim1.d, Ui.d, Uip1.d);
  dq2 = mcd_slope(Ujm1.d, Ui.d, Ujp1.d);
  dq3 = mcd_slope(Ukm1.d, Ui.d, Ukp1.d);
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].d  = Ui.d 
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;;
  }}}

/* 1-momentum */
  dq1 = mcd_slope(Uim1.M1, Ui.M1, Uip1.M1);
  dq2 = mcd_slope(Ujm1.M1, Ui.M1, Ujp1.M1);
  dq3 = mcd_slope(Ukm1.M1, Ui.M1, Ukp1.M1);
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].M1 = Ui.M1 
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;;
  }}}

/* 2-momentum */
  dq1 = mcd_slope(Uim1.M2, Ui.M2, Uip1.M2);
  dq2 = mcd_slope(Ujm1.M2, Ui.M2, Ujp1.M2);
  dq3 = mcd_slope(Ukm1.M2, Ui.M2, Ukp1.M2);
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].M2 = Ui.M2 
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;;
  }}}

/* 3-momentum */
  dq1 = mcd_slope(Uim1.M3, Ui.M3, Uip1.M3);
  dq2 = mcd_slope(Ujm1.M3, Ui.M3, Ujp1.M3);
  dq3 = mcd_slope(Ukm1.M3, Ui.M3, Ukp1.M3);
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].M3 = Ui.M3 
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;;
  }}}

#ifdef MHD
/* 1-cell-centered magnetic field */
  dq1 = mcd_slope(Uim1.B1c, Ui.B1c, Uip1.B1c);
  dq2 = mcd_slope(Ujm1.B1c, Ui.B1c, Ujp1.B1c);
  dq3 = mcd_slope(Ukm1.B1c, Ui.B1c, Ukp1.B1c);
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].B1c = Ui.B1c 
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;;
  }}}

/* 2-cell-centered magnetic field */
  dq1 = mcd_slope(Uim1.B2c, Ui.B2c, Uip1.B2c);
  dq2 = mcd_slope(Ujm1.B2c, Ui.B2c, Ujp1.B2c);
  dq3 = mcd_slope(Ukm1.B2c, Ui.B2c, Ukp1.B2c);
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].B2c = Ui.B2c 
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;;
  }}}

/* 3-cell-centered magnetic field */
  dq1 = mcd_slope(Uim1.B3c, Ui.B3c, Uip1.B3c);
  dq2 = mcd_slope(Ujm1.B3c, Ui.B3c, Ujp1.B3c);
  dq3 = mcd_slope(Ukm1.B3c, Ui.B3c, Ukp1.B3c);
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].B3c = Ui.B3c 
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;;
  }}}
#endif /* MHD */

#ifndef BAROTROPIC
#ifdef SPECIAL_RELATIVITY

/* Prolongate E not P. Otherwise we'd need lots & lots of calls to
 * Con_to_Prim, or a complete rewrite of the code here, so we'll just
 * do this for now */

  dq1 = mcd_slope(Uim1.E, Ui.E, Uip1.E);
  dq2 = mcd_slope(Ujm1.E, Ui.E, Ujp1.E);
  dq3 = mcd_slope(Ukm1.E, Ui.E, Ukp1.E);
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].E = Ui.E 
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;;
  }}}

#else

/* Prolongate P not E.   This is intentionally non-conservative. */

  Pi   = Ui.E   - 0.5*(SQR(Ui.M1  ) + SQR(Ui.M2  ) + SQR(Ui.M3  ))/Ui.d;
  Pim1 = Uim1.E - 0.5*(SQR(Uim1.M1) + SQR(Uim1.M2) + SQR(Uim1.M3))/Uim1.d;
  Pip1 = Uip1.E - 0.5*(SQR(Uip1.M1) + SQR(Uip1.M2) + SQR(Uip1.M3))/Uip1.d;
#ifdef MHD
  Pi   -= 0.5*(SQR(Ui.B1c  ) + SQR(Ui.B2c  ) + SQR(Ui.B3c  ));
  Pim1 -= 0.5*(SQR(Uim1.B1c) + SQR(Uim1.B2c) + SQR(Uim1.B3c));
  Pip1 -= 0.5*(SQR(Uip1.B1c) + SQR(Uip1.B2c) + SQR(Uip1.B3c));
#endif /* MHD */
  dq1 = mcd_slope(Pim1, Pi, Pip1);

  Pjm1 = Ujm1.E - 0.5*(SQR(Ujm1.M1) + SQR(Ujm1.M2) + SQR(Ujm1.M3))/Ujm1.d;
  Pjp1 = Ujp1.E - 0.5*(SQR(Ujp1.M1) + SQR(Ujp1.M2) + SQR(Ujp1.M3))/Ujp1.d;
#ifdef MHD
  Pjm1 -= 0.5*(SQR(Ujm1.B1c) + SQR(Ujm1.B2c) + SQR(Ujm1.B3c));
  Pjp1 -= 0.5*(SQR(Ujp1.B1c) + SQR(Ujp1.B2c) + SQR(Ujp1.B3c));
#endif /* MHD */
  dq2 = mcd_slope(Pjm1, Pi, Pjp1);

  Pkm1 = Ukm1.E - 0.5*(SQR(Ukm1.M1) + SQR(Ukm1.M2) + SQR(Ukm1.M3))/Ukm1.d;
  Pkp1 = Ukp1.E - 0.5*(SQR(Ukp1.M1) + SQR(Ukp1.M2) + SQR(Ukp1.M3))/Ukp1.d;
#ifdef MHD
  Pkm1 -= 0.5*(SQR(Ukm1.B1c) + SQR(Ukm1.B2c) + SQR(Ukm1.B3c));
  Pkp1 -= 0.5*(SQR(Ukp1.B1c) + SQR(Ukp1.B2c) + SQR(Ukp1.B3c));
#endif /* MHD */
  dq3 = mcd_slope(Pkm1, Pi, Pkp1);

  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].E = Pi
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;
    PCon[k][j][i].E += 0.5*(SQR(PCon[k][j][i].M1) + SQR(PCon[k][j][i].M2) +
      SQR(PCon[k][j][i].M3))/PCon[k][j][i].d;
#ifdef MHD
    PCon[k][j][i].E += 0.5*(SQR(PCon[k][j][i].B1c) + SQR(PCon[k][j][i].B2c) +
      SQR(PCon[k][j][i].B3c));
#endif /* MHD */
  }}}

#endif /* SPECIAL_RELATIVITY */
#endif /* BAROTROPIC */

#if (NSCALARS > 0)
/* passive scalars */
  for (n=0; n<NSCALARS; n++) {
    dq1 = mcd_slope(Uim1.s[n], Ui.s[n], Uip1.s[n]);
    dq2 = mcd_slope(Ujm1.s[n], Ui.s[n], Ujp1.s[n]);
    dq3 = mcd_slope(Ukm1.s[n], Ui.s[n], Ukp1.s[n]);
    for (k=0; k<2; k++){
    for (j=0; j<2; j++){
    for (i=0; i<2; i++){
      PCon[k][j][i].s[n] = Ui.s[n] 
        + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;;
    }}}
  }
#endif /* NSCALARS */

#ifdef SPECIAL_RELATIVITY
/* With SR, we need to ensure that the new state is physical, otherwise
 * everything will fall apart at the next time step */
  dfail = 0;
  Pfail = 0;
  Vfail = 0;
  fail = 0;
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    W = check_Prim(&(PCon[k][j][i]));
    Vsq = SQR(W.V1) + SQR(W.V2) + SQR(W.V3);
    if (W.d < 0.0){
      dfail++;
      fail = 1;
    }
    if (W.P < 0.0){
      Pfail++;
      fail = 1;
    }
    if (Vsq > 1.0){
      Vfail++;
      fail = 1;
    }
  }}}

/* If the state is unphysical, revert to first order prologongation */

  if (fail) {
    
    for (k=0; k<2; k++){
      for (j=0; j<2; j++){
        for (i=0; i<2; i++){
          PCon[k][j][i].d  = Ui.d;
          PCon[k][j][i].M1 = Ui.M1;
          PCon[k][j][i].M2 = Ui.M2;
          PCon[k][j][i].M3 = Ui.M3;
#ifndef BAROTROPIC
          PCon[k][j][i].E  = Ui.E;
#endif /* BAROTROPIC */
#ifdef MHD
          PCon[k][j][i].B1c = Ui.B1c;
          PCon[k][j][i].B2c = Ui.B2c;
          PCon[k][j][i].B3c = Ui.B3c;
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) PCon[k][j][i].s[n] = Ui.s[n];
#endif /* NSCALARS */
        }}}
  }

  dfail = 0;
  Pfail = 0;
  Vfail = 0;
  fail = 0;

#endif SPECIAL_RELATIVITY

#endif /* FIRST_ORDER */
}

/*----------------------------------------------------------------------------*/
/*! \fn void ProFld(Real3Vect BGZ[][3][3], Real3Vect PFld[][3][3], 
 *            const Real dx1c, const Real dx2c, const Real dx3c)
 *  \brief Uses the divergence-preserving prolongation operators of
 * Toth & Roe (JCP, 180, 736, 2002) to interpolate the face centered fields
 * in a 3x2x2 block for Bx, 2x3x2 block for By, and 2x2x3 block for Bz.
 */

#ifdef MHD
void ProFld(Real3Vect BGZ[][3][3], Real3Vect PFld[][3][3], 
            const Real dx1c, const Real dx2c, const Real dx3c)
{
  int i,j,k;
  Real dBdx,dBdy,dBdz,Uxx,Vyy,Wzz,Uxyz,Vxyz,Wxyz;

/* initialize Bx on left-x1 boundry, if not set already */

  if (PFld[0][0][0].x == 0.0) {
    dBdy = mcd_slope(BGZ[1][0][1].x, BGZ[1][1][1].x, BGZ[1][2][1].x);
    dBdz = mcd_slope(BGZ[0][1][1].x, BGZ[1][1][1].x, BGZ[2][1][1].x);

    PFld[0][0][0].x = BGZ[1][1][1].x - 0.25*dBdy - 0.25*dBdz;
    PFld[0][1][0].x = BGZ[1][1][1].x + 0.25*dBdy - 0.25*dBdz;
    PFld[1][0][0].x = BGZ[1][1][1].x - 0.25*dBdy + 0.25*dBdz;
    PFld[1][1][0].x = BGZ[1][1][1].x + 0.25*dBdy + 0.25*dBdz;
  }

/* initialize Bx on right-x1 boundry, if not set already */

  if (PFld[0][0][2].x == 0.0) {
    dBdy = mcd_slope(BGZ[1][0][2].x, BGZ[1][1][2].x, BGZ[1][2][2].x);
    dBdz = mcd_slope(BGZ[0][1][2].x, BGZ[1][1][2].x, BGZ[2][1][2].x);

    PFld[0][0][2].x = BGZ[1][1][2].x - 0.25*dBdy - 0.25*dBdz;
    PFld[0][1][2].x = BGZ[1][1][2].x + 0.25*dBdy - 0.25*dBdz;
    PFld[1][0][2].x = BGZ[1][1][2].x - 0.25*dBdy + 0.25*dBdz;
    PFld[1][1][2].x = BGZ[1][1][2].x + 0.25*dBdy + 0.25*dBdz;
  }

/* initialize By on left-x2 boundry, if not set already */

  if (PFld[0][0][0].y == 0.0) {
    dBdx = mcd_slope(BGZ[1][1][0].y, BGZ[1][1][1].y, BGZ[1][1][2].y);
    dBdz = mcd_slope(BGZ[0][1][1].y, BGZ[1][1][1].y, BGZ[2][1][1].y);

    PFld[0][0][0].y = BGZ[1][1][1].y - 0.25*dBdx - 0.25*dBdz;
    PFld[0][0][1].y = BGZ[1][1][1].y + 0.25*dBdx - 0.25*dBdz;
    PFld[1][0][0].y = BGZ[1][1][1].y - 0.25*dBdx + 0.25*dBdz;
    PFld[1][0][1].y = BGZ[1][1][1].y + 0.25*dBdx + 0.25*dBdz;
  }

/* initialize By on right-x2 boundry, if not set already */

  if (PFld[0][2][0].y == 0.0) {
    dBdx = mcd_slope(BGZ[1][2][0].y, BGZ[1][2][1].y, BGZ[1][2][2].y);
    dBdz = mcd_slope(BGZ[0][2][1].y, BGZ[1][2][1].y, BGZ[2][2][1].y);

    PFld[0][2][0].y = BGZ[1][2][1].y - 0.25*dBdx - 0.25*dBdz;
    PFld[0][2][1].y = BGZ[1][2][1].y + 0.25*dBdx - 0.25*dBdz;
    PFld[1][2][0].y = BGZ[1][2][1].y - 0.25*dBdx + 0.25*dBdz;
    PFld[1][2][1].y = BGZ[1][2][1].y + 0.25*dBdx + 0.25*dBdz;
  }

/* initialize Bz on left-x3 boundry, if not set already */

  if (PFld[0][0][0].z == 0.0) {
    dBdx = mcd_slope(BGZ[1][1][0].z, BGZ[1][1][1].z, BGZ[1][1][2].z);
    dBdy = mcd_slope(BGZ[1][0][1].z, BGZ[1][1][1].z, BGZ[1][2][1].z);

    PFld[0][0][0].z = BGZ[1][1][1].z - 0.25*dBdx - 0.25*dBdy;
    PFld[0][0][1].z = BGZ[1][1][1].z + 0.25*dBdx - 0.25*dBdy;
    PFld[0][1][0].z = BGZ[1][1][1].z - 0.25*dBdx + 0.25*dBdy;
    PFld[0][1][1].z = BGZ[1][1][1].z + 0.25*dBdx + 0.25*dBdy;
  }

/* initialize Bz on right-x3 boundry, if not set already */

  if (PFld[2][0][0].z == 0.0) {
    dBdx = mcd_slope(BGZ[2][1][0].z, BGZ[2][1][1].z, BGZ[2][1][2].z);
    dBdy = mcd_slope(BGZ[2][0][1].z, BGZ[2][1][1].z, BGZ[2][2][1].z);

    PFld[2][0][0].z = BGZ[2][1][1].z - 0.25*dBdx - 0.25*dBdy;
    PFld[2][0][1].z = BGZ[2][1][1].z + 0.25*dBdx - 0.25*dBdy;
    PFld[2][1][0].z = BGZ[2][1][1].z - 0.25*dBdx + 0.25*dBdy;
    PFld[2][1][1].z = BGZ[2][1][1].z + 0.25*dBdx + 0.25*dBdy;
  }

/* Fill in the face-centered fields in the interior of the cell using the
 * interpolation formulae of T&R, eqs. 8-12.  The k=0,1 terms have been written
 * out explicetely so they can be grouped to reduce round-off error  */

  Uxx = Vyy = Wzz = 0.0;
  Uxyz = Vxyz = Wxyz = 0.0;
  for(j=0; j<2; j++){
  for(i=0; i<2; i++){
    Uxx += (2*i-1)*((2*j-1)*dx3c*(PFld[0][2*j][i].y + PFld[1][2*j][i].y) +
                            dx2c*(PFld[2][j  ][i].z - PFld[0][j  ][i].z) );

    Vyy += (2*j-1)*(        dx1c*(PFld[2][j][i  ].z - PFld[0][j][i  ].z) +
                    (2*i-1)*dx3c*(PFld[0][j][2*i].x + PFld[1][j][2*i].x) );

    Wzz += ((2*i-1)*dx2c*(PFld[1][j][2*i].x - PFld[0][j][2*i].x) +
            (2*j-1)*dx1c*(PFld[1][2*j][i].y - PFld[0][2*j][i].y) );

    Uxyz += (2*i-1)*(2*j-1)*(PFld[1][j][2*i].x - PFld[0][j][2*i].x);
    Vxyz += (2*i-1)*(2*j-1)*(PFld[1][2*j][i].y - PFld[0][2*j][i].y);
    Wxyz += (2*i-1)*(2*j-1)*(PFld[2][j][i].z - PFld[0][j][i].z);
  }}

/* Multiply through by some common factors */

  Uxx *= 0.125*dx1c;
  Vyy *= 0.125*dx2c;
  Wzz *= 0.125*dx3c;
  Uxyz *= 0.125*dx2c*dx3c/(dx2c*dx2c + dx3c*dx3c);
  Vxyz *= 0.125*dx1c*dx3c/(dx1c*dx1c + dx3c*dx3c);
  Wxyz *= 0.125*dx1c*dx2c/(dx1c*dx1c + dx2c*dx2c);

/* Initialize B1i on interior faces */

  for(k=0; k<2; k++){
  for(j=0; j<2; j++){
    PFld[k][j][1].x = 0.5*(PFld[k][j][0].x + PFld[k][j][2].x) + Uxx/(dx2c*dx3c)
       + (2*k-1)*(dx3c/dx2c)*Vxyz + (2*j-1)*(dx2c/dx3c)*Wxyz;
  }}

/* Initialize B2i on interior faces */

  for(k=0; k<2; k++){
  for(i=0; i<2; i++){
    PFld[k][1][i].y = 0.5*(PFld[k][0][i].y + PFld[k][2][i].y) + Vyy/(dx3c*dx1c)
      + (2*i-1)*(dx1c/dx3c)*Wxyz + (2*k-1)*(dx3c/dx1c)*Uxyz;
  }}

/* Initialize B3i on interior faces */

  for(j=0; j<2; j++){
  for(i=0; i<2; i++){
    PFld[1][j][i].z = 0.5*(PFld[0][j][i].z + PFld[2][j][i].z) + Wzz/(dx1c*dx2c)
      + (2*j-1)*(dx2c/dx1c)*Uxyz + (2*i-1)*(dx1c/dx2c)*Vxyz;
  }}

}
#endif /* MHD */

/*----------------------------------------------------------------------------*/
/*! \fn static Real mcd_slope(const Real vl, const Real vc, const Real vr)
 *  \brief Computes monotonized linear slope.
 */

#ifndef FIRST_ORDER
static Real mcd_slope(const Real vl, const Real vc, const Real vr){

  Real dvl = (vc - vl), dvr = (vr - vc);
  Real dv, dvm;

  if(dvl > 0.0 && dvr > 0.0){
    dv = 2.0*(dvl < dvr ? dvl : dvr);
    dvm = 0.5*(dvl + dvr);
    return (dvm < dv ? dvm : dv);
  }
  else if(dvl < 0.0 && dvr < 0.0){
    dv = 2.0*(dvl > dvr ? dvl : dvr);
    dvm = 0.5*(dvl + dvr);
    return (dvm > dv ? dvm : dv);
  }

  return 0.0;
}
#endif /* FIRST_ORDER */

#endif /* STATIC_MESH_REFINEMENT */
