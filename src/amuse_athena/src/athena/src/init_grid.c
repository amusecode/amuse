#include "copyright.h"
/*============================================================================*/
/*! \file init_grid.c 
 *  \brief Initializes most variables in the Grid structure.
 *
 * PURPOSE: Initializes most variables in the Grid structure.  Allocates memory
 *   for 3D arrays of Cons, interface B, etc.  With SMR, finds all overlaps
 *   between child and parent Grids, and initializes data needed for restriction
 *   flux-correction, and prolongation steps.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - init_grid()
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - checkOverlap() - checks for overlap of cubes, and returns overlap coords
 * - checkOverlapTouch() - same as above, but checks for overlap and/or touch */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *  checkOverlap() - checks for overlap of cubes, and returns overlap coords
 *  checkOverlapTouch() - same as above, but checks for overlap and/or touch
 *============================================================================*/
#ifdef STATIC_MESH_REFINEMENT
/*! \fn int checkOverlap(SideS *pC1, SideS *pC2, SideS *pC3);
 *  \brief Checks for overlap of cubes, and returns overlap coords */
int checkOverlap(SideS *pC1, SideS *pC2, SideS *pC3);
/*! \fn int checkOverlapTouch(SideS *pC1, SideS *pC2, SideS *pC3);
 *  \brief Same as above, but checks for overlap and/or touch */
int checkOverlapTouch(SideS *pC1, SideS *pC2, SideS *pC3);
#endif

/*----------------------------------------------------------------------------*/
/*! \fn void init_grid(MeshS *pM)
 *  \brief Initializes most variables in the Grid structure.
 *
 * PURPOSE: Initializes most variables in the Grid structure.  Allocates memory
 *   for 3D arrays of Cons, interface B, etc.  With SMR, finds all overlaps
 *   between child and parent Grids, and initializes data needed for restriction
 *   flux-correction, and prolongation steps.				      */

void init_grid(MeshS *pM)
{
  DomainS *pD;
  GridS *pG;
  int nDim,nl,nd,myL,myM,myN;
  int i,l,m,n,n1z,n2z,n3z,n1p,n2p,n3p;
#ifdef STATIC_MESH_REFINEMENT
  DomainS *pCD,*pPD;
  SideS D1,D2,D3,G1,G2,G3;
  int isDOverlap,isGOverlap,irefine,ncd,npd,dim,iGrid;
  int ncg,nCG,nMyCG,nCB[6],nMyCB[6],nb;
  int npg,nPG,nMyPG,nPB[6],nMyPB[6];
#endif

/* number of dimensions in Grid. */
  nDim=1;
  for (i=1; i<3; i++) if (pM->Nx[i]>1) nDim++;

/* Loop over all levels and domains per level */

  for (nl=0; nl<pM->NLevels; nl++){
  for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {
      pD = (DomainS*)&(pM->Domain[nl][nd]);  /* set ptr to Domain */
      pG = pM->Domain[nl][nd].Grid;          /* set ptr to Grid */

      pG->time = pM->time;

/* get (l,m,n) coordinates of Grid being updated on this processor */

      get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);

/* ---------------------  Intialize grid in 1-direction --------------------- */
/* Initialize is,ie,dx1
 * Compute Disp, MinX[0], and MaxX[0] using displacement of Domain and Grid
 * location within Domain */

      pG->Nx[0] = pD->GData[myN][myM][myL].Nx[0];

      if(pG->Nx[0] > 1) {
        pG->is = nghost;
        pG->ie = pG->Nx[0] + nghost - 1;
      }
      else
        pG->is = pG->ie = 0;
    
      pG->dx1 = pD->dx[0];
    
      pG->Disp[0] = pD->Disp[0];
      pG->MinX[0] = pD->MinX[0];
      for (l=1; l<=myL; l++) {
        pG->Disp[0] +=        pD->GData[myN][myM][l-1].Nx[0];
        pG->MinX[0] += (Real)(pD->GData[myN][myM][l-1].Nx[0])*pG->dx1;
      }
      pG->MaxX[0] = pG->MinX[0] + (Real)(pG->Nx[0])*pG->dx1;
    
/* ---------------------  Intialize grid in 2-direction --------------------- */
/* Initialize js,je,dx2
 * Compute Disp, MinX[1], and MaxX[1] using displacement of Domain and Grid
 * location within Domain */

      pG->Nx[1] = pD->GData[myN][myM][myL].Nx[1];
    
      if(pG->Nx[1] > 1) {
        pG->js = nghost;
        pG->je = pG->Nx[1] + nghost - 1;
      }
      else
        pG->js = pG->je = 0;
    
      pG->dx2 = pD->dx[1];

      pG->Disp[1] = pD->Disp[1];
      pG->MinX[1] = pD->MinX[1];
      for (m=1; m<=myM; m++) {
        pG->Disp[1] +=        pD->GData[myN][m-1][myL].Nx[1];
        pG->MinX[1] += (Real)(pD->GData[myN][m-1][myL].Nx[1])*pG->dx2;
      }
      pG->MaxX[1] = pG->MinX[1] + (Real)(pG->Nx[1])*pG->dx2;

/* ---------------------  Intialize grid in 3-direction --------------------- */
/* Initialize ks,ke,dx3
 * Compute Disp, MinX[2], and MaxX[2] using displacement of Domain and Grid
 * location within Domain */

      pG->Nx[2] = pD->GData[myN][myM][myL].Nx[2];

      if(pG->Nx[2] > 1) {
        pG->ks = nghost;
        pG->ke = pG->Nx[2] + nghost - 1;
      }
      else
        pG->ks = pG->ke = 0;

      pG->dx3 = pD->dx[2];

      pG->Disp[2] = pD->Disp[2];
      pG->MinX[2] = pD->MinX[2];
      for (n=1; n<=myN; n++) {
        pG->Disp[2] +=        pD->GData[n-1][myM][myL].Nx[2];
        pG->MinX[2] += (Real)(pD->GData[n-1][myM][myL].Nx[2])*pG->dx3;
      }
      pG->MaxX[2] = pG->MinX[2] + (Real)(pG->Nx[2])*pG->dx3;

/* ---------  Allocate 3D arrays to hold Cons based on size of grid --------- */

      if (pG->Nx[0] > 1)
        n1z = pG->Nx[0] + 2*nghost;
      else
        n1z = 1;

      if (pG->Nx[1] > 1)
        n2z = pG->Nx[1] + 2*nghost;
      else
        n2z = 1;

      if (pG->Nx[2] > 1)
        n3z = pG->Nx[2] + 2*nghost;
      else
        n3z = 1;

/* Build a 3D array of type ConsS */

      pG->U = (ConsS***)calloc_3d_array(n3z, n2z, n1z, sizeof(ConsS));
      if (pG->U == NULL) goto on_error1;
    
/* Build 3D arrays to hold interface field */

#ifdef MHD
      pG->B1i = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->B1i == NULL) goto on_error2;

      pG->B2i = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->B2i == NULL) goto on_error3;

      pG->B3i = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->B3i == NULL) goto on_error4;
#endif /* MHD */

/* Build 3D arrays to magnetic diffusivities */

#ifdef RESISTIVITY
      pG->eta_Ohm = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->eta_Ohm == NULL) goto on_error5;

      pG->eta_Hall = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->eta_Hall == NULL) goto on_error6;

      pG->eta_AD = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->eta_AD == NULL) goto on_error7;
#endif /* RESISTIVITY */

/* Build 3D arrays to gravitational potential and mass fluxes */

#ifdef SELF_GRAVITY
      pG->Phi = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->Phi == NULL) goto on_error9;

      pG->Phi_old = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->Phi_old == NULL) goto on_error10;

      pG->x1MassFlux = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->x1MassFlux == NULL) goto on_error11;

      pG->x2MassFlux = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->x2MassFlux == NULL) goto on_error12;

      pG->x3MassFlux = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
      if (pG->x3MassFlux == NULL) goto on_error13;
#endif /* SELF_GRAVITY */

/* Allocate and initialize cylindrical scaling factors */
#ifdef CYLINDRICAL
      pG->r = (Real*)calloc_1d_array(n1z, sizeof(Real));
      if (pG->r == NULL) goto on_error14;

      pG->ri = (Real*)calloc_1d_array(n1z, sizeof(Real));
      if (pG->ri == NULL) goto on_error15;
      for (i=pG->is-nghost; i<=pG->ie+nghost; i++) {
        pG->ri[i] = pG->MinX[0] + ((Real)(i - pG->is))*pG->dx1;
        pG->r[i]  = pG->ri[i] + 0.5*pG->dx1;
      }
#endif /* CYLINDRICAL */


/*-- Get IDs of neighboring Grids in Domain communicator ---------------------*/
/* If Grid is at the edge of the Domain (so it is either a physical boundary,
 * or an internal boundary between fine/coarse grids), then ID is set to -1
 */

/* Left-x1 */
      if(myL > 0) pG->lx1_id = pD->GData[myN][myM][myL-1].ID_Comm_Domain;
      else pG->lx1_id = -1;

/* Right-x1 */
      if(myL <(pD->NGrid[0])-1)
        pG->rx1_id = pD->GData[myN][myM][myL+1].ID_Comm_Domain;
      else pG->rx1_id = -1;

/* Left-x2 */
      if(myM > 0) pG->lx2_id = pD->GData[myN][myM-1][myL].ID_Comm_Domain;
      else pG->lx2_id = -1;

/* Right-x2 */
      if(myM <(pD->NGrid[1])-1)
        pG->rx2_id = pD->GData[myN][myM+1][myL].ID_Comm_Domain;
      else pG->rx2_id = -1;

/* Left-x3 */
      if(myN > 0) pG->lx3_id = pD->GData[myN-1][myM][myL].ID_Comm_Domain;
      else pG->lx3_id = -1;

/* Right-x3 */
      if(myN <(pD->NGrid[2])-1)
        pG->rx3_id = pD->GData[myN+1][myM][myL].ID_Comm_Domain;
      else pG->rx3_id = -1;
   
#ifdef STATIC_MESH_REFINEMENT
/*---------------------- Initialize variables for SMR ------------------------*/
/* Number of child/parent grids, and data about overlap regions. */

      pG->NCGrid = 0;
      pG->NPGrid = 0;
      pG->NmyCGrid = 0;  /* can be as large as # of child Domains */
      pG->NmyPGrid = 0;  /* must be 0 or 1 */
      pG->CGrid = NULL;
      pG->PGrid = NULL;
#endif

    }
  }}

#ifdef STATIC_MESH_REFINEMENT
/*-------------- Count number of child Grids, and boundaries -----------------*/
/* First we have to count the number of child Grids, and fine/coarse boundaries
 * with child Grids, including child Grids on this and other processors, before
 * we can allocate the CGrid array.  This is a shortcoming of making these
 * structures 1D arrays. */ 
/* Loop over levels (up to next to last level: maxlevel-1), and domains/level */

  for (nl=0; nl<(pM->NLevels)-1; nl++){
  for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {
      pD = (DomainS*)&(pM->Domain[nl][nd]);  /* set ptr to this Domain */
      pG = pM->Domain[nl][nd].Grid;          /* set ptr to this Grid */

/* edges of this Domain */
      for (i=0; i<3; i++) {
        D1.ijkl[i] = pD->Disp[i];
        D1.ijkr[i] = pD->Disp[i] + pD->Nx[i];
      }

/* edges of this Grid */
      for (i=0; i<3; i++) {
        G1.ijkl[i] = pG->Disp[i];
        G1.ijkr[i] = pG->Disp[i] + pG->Nx[i];
      }

/* For this Domain, check if there is a Domain at next level that overlaps.
 * Divide by two to check in units of this (not the child) domain coordinates */

      for (ncd=0; ncd<pM->DomainsPerLevel[nl+1]; ncd++){
        pCD = (DomainS*)&(pM->Domain[nl+1][ncd]);  /* ptr to potential child  */

/* edges of potential child Domain */
        for (i=0; i<3; i++) {
          D2.ijkl[i] = pCD->Disp[i]/2;
          D2.ijkr[i] = 1;
          if (pCD->Nx[i] > 1) D2.ijkr[i] = (pCD->Disp[i] + pCD->Nx[i])/2;
        }

        isDOverlap = checkOverlap(&D1, &D2, &D3);
        if (isDOverlap == 1){

/*----------------------------------------------------------------------------*/
/* There is a child Domain that overlaps. So on the child Domain, find all the
 * Grids that overlap this Grid. */

          for (n=0; n<pCD->NGrid[2]; n++){
          for (m=0; m<pCD->NGrid[1]; m++){
          for (l=0; l<pCD->NGrid[0]; l++){

/* edges of child Grid */
/* Divide by two to check in units of this (not the child) grid coordinates */
            for (i=0; i<3; i++) {
              G2.ijkl[i] = pCD->GData[n][m][l].Disp[i]/2;
              G2.ijkr[i] = 1;
              if (pCD->Nx[i] > 1) 
                G2.ijkr[i] = (pCD->GData[n][m][l].Disp[i]
                            + pCD->GData[n][m][l].Nx[i])/2;
            }

            isGOverlap = checkOverlap(&G1, &G2, &G3);

/* If Grid overlaps, increment Child counters */

            if (isGOverlap == 1){
              pG->NCGrid++;
              if (pCD->GData[n][m][l].ID_Comm_world == myID_Comm_world)
                pG->NmyCGrid++;
            }

          }}} /* end loops over [n,m,l]: all Grids in Domain[nl+1][ncd] */
        }
      } /* end loop over all child Domains at level [nl+1] */

/*------------------------- Allocate CGrid array -----------------------------*/
/* Now we know how many child Grids there are for the Grid in Domain[nd] at
 * level nl on both this and other processors.    */

      if (pG->NCGrid > 0) {
        pG->CGrid =(GridOvrlpS*)calloc_1d_array(pG->NCGrid, sizeof(GridOvrlpS));
        if(pG->CGrid==NULL) ath_error("[init_grid]:failed to allocate CGrid\n");

        for (ncg=0; ncg<pG->NCGrid; ncg++){
          for (dim=0; dim<6; dim++) {
            pG->CGrid[ncg].myFlx[dim] = NULL;
#ifdef MHD
            pG->CGrid[ncg].myEMF1[dim] = NULL;
            pG->CGrid[ncg].myEMF2[dim] = NULL;
            pG->CGrid[ncg].myEMF3[dim] = NULL;
#endif /* MHD */
          }
        }
      }

/*-------------------------- Fill in CGrid array -----------------------------*/
/* Repeat loop over all domains at next level, and all the logic to find
 * overlapping Grids, to fill in data about overlap regions in CGrid array.
 * This isn't a particularly pretty way to do it. */

      nMyCG = 0;
      nCG = pG->NmyCGrid;

      for (ncd=0; ncd<pM->DomainsPerLevel[nl+1]; ncd++){
        pCD = (DomainS*)&(pM->Domain[nl+1][ncd]);   /* ptr to potential child */

/* edges of potential child Domain */
        for (i=0; i<3; i++) {
          D2.ijkl[i] = pCD->Disp[i]/2;
          D2.ijkr[i] = 1;
          if (pCD->Nx[i] > 1) D2.ijkr[i] = (pCD->Disp[i] + pCD->Nx[i])/2;
        }

        isDOverlap = checkOverlap(&D1, &D2, &D3);
        if (isDOverlap == 1){

/*----------------------------------------------------------------------------*/
/* Found the Domain that overlaps, so on the child Domain check if there
 * is a Grid that overlaps */

          for (n=0; n<pCD->NGrid[2]; n++){
          for (m=0; m<pCD->NGrid[1]; m++){
          for (l=0; l<pCD->NGrid[0]; l++){

/* edges of child Grid */
/* Divide by two to check in units of this (not the child) grid coordinates */
            for (i=0; i<3; i++) {
              G2.ijkl[i] = pCD->GData[n][m][l].Disp[i]/2;
              G2.ijkr[i] = 1;
              if (pCD->Nx[i] > 1)
                G2.ijkr[i] = (pCD->GData[n][m][l].Disp[i]
                            + pCD->GData[n][m][l].Nx[i])/2;
            }

            isGOverlap = checkOverlap(&G1, &G2, &G3);
            if (isGOverlap == 1){

/* If Grid OVERLAPS, then:
 * (1) fill-in data in CGrid array */
/* Index CGrid array so that child Grids on this processor come first */

              if (pCD->GData[n][m][l].ID_Comm_world == myID_Comm_world) {
                ncg=nMyCG;
                nMyCG++;
              } else {
                ncg=nCG;
                nCG++;
              }

              pG->CGrid[ncg].ijks[0] = G3.ijkl[0] - pG->Disp[0] + pG->is;
              pG->CGrid[ncg].ijke[0] = G3.ijkr[0] - pG->Disp[0] + pG->is - 1;
              pG->CGrid[ncg].ijks[1] = G3.ijkl[1] - pG->Disp[1] + pG->js;
              pG->CGrid[ncg].ijke[1] = G3.ijkr[1] - pG->Disp[1] + pG->js - 1;
              pG->CGrid[ncg].ijks[2] = G3.ijkl[2] - pG->Disp[2] + pG->ks;
              pG->CGrid[ncg].ijke[2] = G3.ijkr[2] - pG->Disp[2] + pG->ks - 1;
  
              pG->CGrid[ncg].DomN = ncd;
              pG->CGrid[ncg].ID = pCD->GData[n][m][l].ID_Comm_Parent;

              n1z = pG->CGrid[ncg].ijke[0] - pG->CGrid[ncg].ijks[0] + 1;
              n2z = pG->CGrid[ncg].ijke[1] - pG->CGrid[ncg].ijks[1] + 1;
              n3z = pG->CGrid[ncg].ijke[2] - pG->CGrid[ncg].ijks[2] + 1;
              pG->CGrid[ncg].nWordsRC = n1z*n2z*n3z*(NVAR);
              pG->CGrid[ncg].nWordsP  = 0;
if(myID_Comm_world==0){
printf("\nCGrid overlap is %d x %d x %d\n",n1z,n2z,n3z);
}
#ifdef MHD
              if (nDim==3) {
                pG->CGrid[ncg].nWordsRC += 
                  (n1z+1)*n2z*n3z + n1z*(n2z+1)*n3z + n1z*n2z*(n3z+1);
              } else {
                if (nDim==2) {
                  pG->CGrid[ncg].nWordsRC += (n1z+1)*n2z + n1z*(n2z+1);
                }
              }
#endif /* MHD */

/* (2) if edge of the child Grid is at edge of the child Domain, then allocate
 * memory for fluxes and EMFs for Correction, and count GZ for Prolongation */

              for (dim=0; dim<nDim; dim++){    /* only checks nDim directions */
                if (dim == 0) iGrid=l;
                if (dim == 1) iGrid=m;
                if (dim == 2) iGrid=n;

/* inner x1/x2/x3 boundary */
/* First check that edge of child Grid is at edge of overlap (otherwise boundary
 * is between MPI blocks in parent, and is internal to child Grid).
 * Then check that edge of child Grid is not at l-edge of root (so that physical
 * BCs are applied), but is at l-edge of own Domain (so it is not an internal
 * MPI boundary on the child Domain). */

                if ((G2.ijkl[dim] == G3.ijkl[dim]) &&
                    (pCD->Disp[dim] != 0) &&
                    (iGrid == 0)) {

                  if (dim == 0) {
                    n1z = G3.ijkr[1] - G3.ijkl[1];
                    n2z = G3.ijkr[2] - G3.ijkl[2];
                    n1p = n1z;
                    n2p = n2z;
                    if (pG->Nx[1] > 1) n1p += nghost + 2;
                    if (pG->Nx[2] > 1) n2p += nghost + 2;
                  }
                  if (dim == 1) {
                    n1z = G3.ijkr[0] - G3.ijkl[0];
                    n2z = G3.ijkr[2] - G3.ijkl[2];
                    n1p = n1z + nghost + 2;
                    n2p = n2z;
                    if (pG->Nx[2] > 1) n2p += nghost + 2;
                  }
                  if (dim == 2) {
                    n1z = G3.ijkr[0] - G3.ijkl[0];
                    n2z = G3.ijkr[1] - G3.ijkl[1];
                    n1p = n1z + nghost + 2;
                    n2p = n2z + nghost + 2;
                  }

                  pG->CGrid[ncg].nWordsRC += n1z*n2z*(NVAR); 
                  pG->CGrid[ncg].nWordsP  += ((nghost/2)+2)*n1p*n2p*(NVAR); 

/* Allocate memory for myFlx and my EMFs */

                  pG->CGrid[ncg].myFlx[2*dim] = (ConsS**)calloc_2d_array(
                    n2z,n1z, sizeof(ConsS));
                  if(pG->CGrid[ncg].myFlx[2*dim] == NULL) ath_error(
                   "[init_grid]:failed to allocate CGrid ixb myFlx\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for ixb CGrid.myFlx\n",n2z,n1z);
}
#ifdef MHD
                  pG->CGrid[ncg].nWordsP += 6*((nghost/2)+2)*n1p*n2p;

                  if (pG->Nx[1] > 1 && dim != 2) {
                    pG->CGrid[ncg].nWordsRC += (n1z+1)*n2z; 
                    pG->CGrid[ncg].myEMF3[2*dim] = (Real**)calloc_2d_array(
                      n2z,n1z+1, sizeof(Real));
                    if(pG->CGrid[ncg].myEMF3[2*dim] == NULL) ath_error(
                      "[init_grid]:failed to allocate CGrid ixb myEMF3\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for ixb CGrid.myEMF3\n",n2z,n1z+1);
}
                  }

                  if (pG->Nx[2] > 1  && dim == 0) {
                    pG->CGrid[ncg].nWordsRC += n1z*(n2z+1); 
                    pG->CGrid[ncg].myEMF2[2*dim] = (Real**)calloc_2d_array(
                      n2z+1,n1z, sizeof(Real));
                    if(pG->CGrid[ncg].myEMF2[2*dim] == NULL) ath_error(
                      "[init_grid]:failed to allocate CGrid ixb myEMF2\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for ixb CGrid.myEMF2\n",n2z+1,n1z);
}
                  }

                  if (pG->Nx[2] > 1  && dim == 1) {
                    pG->CGrid[ncg].nWordsRC += n1z*(n2z+1); 
                    pG->CGrid[ncg].myEMF1[2*dim] = (Real**)calloc_2d_array(
                      n2z+1,n1z, sizeof(Real));
                    if(pG->CGrid[ncg].myEMF1[2*dim] == NULL) ath_error(
                      "[init_grid]:failed to allocate CGrid ixb myEMF1\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for ixb CGrid.myEMF1\n",n2z+1,n1z);
}
                  }

                  if (pG->Nx[2] > 1  && dim == 2) {
                    pG->CGrid[ncg].nWordsRC += n1z*(n2z+1) + (n1z+1)*n2z; 
                    pG->CGrid[ncg].myEMF1[2*dim] = (Real**)calloc_2d_array(
                      n2z+1,n1z, sizeof(Real));
                    if(pG->CGrid[ncg].myEMF1[2*dim] == NULL) ath_error(
                      "[init_grid]:failed to allocate CGrid ixb myEMF1\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for ixb CGrid.myEMF1\n",n2z+1,n1z);
}
                    pG->CGrid[ncg].myEMF2[2*dim] = (Real**)calloc_2d_array(
                      n2z,n1z+1, sizeof(Real));
                    if(pG->CGrid[ncg].myEMF2[2*dim] == NULL) ath_error(
                      "[init_grid]:failed to allocate CGrid ixb myEMF2\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for ixb CGrid.myEMF2\n",n2z,n1z+1);
}
                  }
#endif /* MHD */

                }

/* outer x1/x2/x3 boundary */
/* First check that edge of child Grid is at edge of overlap (otherwise boundary
 * is between MPI blocks in parent, and is internal to child Grid).
 * Then check that edge of child Grid is not at r-edge of root (so that physical
 * BCs are applied), but is at r-edge of own Domain (so it is not an internal
 * MPI boundary on the child Domain). */

                irefine = 1;
                for (i=1;i<=(nl+1);i++) irefine *= 2; /* child refinement lev */
                if ( (G2.ijkr[dim] == G3.ijkr[dim]) &&
                    ((pCD->Disp[dim] + pCD->Nx[dim])/irefine != pM->Nx[dim]) &&
                     (iGrid == (pCD->NGrid[dim]-1)) ) {

                  if (dim == 0) {
                    n1z = G3.ijkr[1] - G3.ijkl[1];
                    n2z = G3.ijkr[2] - G3.ijkl[2];
                    n1p = n1z;
                    n2p = n2z;
                    if (pG->Nx[1] > 1) n1p += nghost + 2;
                    if (pG->Nx[2] > 1) n2p += nghost + 2;
                  }
                  if (dim == 1) {
                    n1z = G3.ijkr[0] - G3.ijkl[0];
                    n2z = G3.ijkr[2] - G3.ijkl[2];
                    n1p = n1z + nghost + 2;
                    n2p = n2z;
                    if (pG->Nx[2] > 1) n2p += nghost + 2;
                  }
                  if (dim == 2) {
                    n1z = G3.ijkr[0] - G3.ijkl[0];
                    n2z = G3.ijkr[1] - G3.ijkl[1];
                    n1p = n1z + nghost + 2;
                    n2p = n2z + nghost + 2;
                  }

                  pG->CGrid[ncg].nWordsRC += n1z*n2z*(NVAR); 
                  pG->CGrid[ncg].nWordsP  += ((nghost/2)+2)*n1p*n2p*(NVAR); 

/* Allocate memory for myFlx and myEMFs*/

                  pG->CGrid[ncg].myFlx[(2*dim)+1] = (ConsS**)calloc_2d_array(
                    n2z,n1z, sizeof(ConsS));
                  if(pG->CGrid[ncg].myFlx[(2*dim)+1] == NULL) ath_error(
                    "[init_grid]:failed to allocate CGrid oxb myFlx\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for oxb CGrid.myFlx\n",n2z,n1z);
}
#ifdef MHD
                  pG->CGrid[ncg].nWordsP += 6*((nghost/2)+2)*n1p*n2p;

                  if (pG->Nx[1] > 1 && dim != 2) {
                    pG->CGrid[ncg].nWordsRC += (n1z+1)*n2z;
                    pG->CGrid[ncg].myEMF3[(2*dim)+1] =(Real**)calloc_2d_array(
                      n2z,n1z+1, sizeof(Real));
                    if(pG->CGrid[ncg].myEMF3[(2*dim)+1] == NULL) ath_error(
                      "[init_grid]:failed to allocate CGrid oxb myEMF3\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for oxb CGrid.myEMF3\n",n2z,n1z+1);
}
                  }

                  if (pG->Nx[2] > 1  && dim == 0) {
                    pG->CGrid[ncg].nWordsRC += n1z*(n2z+1);
                    pG->CGrid[ncg].myEMF2[(2*dim)+1] =(Real**)calloc_2d_array(
                      n2z+1,n1z, sizeof(Real));
                    if(pG->CGrid[ncg].myEMF2[(2*dim)+1] == NULL) ath_error(
                      "[init_grid]:failed to allocate CGrid oxb myEMF2\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for oxb CGrid.myEMF2\n",n2z+1,n1z);
}
                  }

                  if (pG->Nx[2] > 1  && dim == 1) {
                    pG->CGrid[ncg].nWordsRC += n1z*(n2z+1);
                    pG->CGrid[ncg].myEMF1[(2*dim)+1] =(Real**)calloc_2d_array(
                      n2z+1,n1z, sizeof(Real));
                    if(pG->CGrid[ncg].myEMF1[(2*dim)+1] == NULL) ath_error(
                      "[init_grid]:failed to allocate CGrid oxb myEMF1\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for oxb CGrid.myEMF1\n",n2z+1,n1z);
}
                  }

                  if (pG->Nx[2] > 1  && dim == 2) {
                    pG->CGrid[ncg].nWordsRC += n1z*(n2z+1) + (n1z+1)*n2z;
                    pG->CGrid[ncg].myEMF1[(2*dim)+1] =(Real**)calloc_2d_array(
                      n2z+1,n1z, sizeof(Real));
                    if(pG->CGrid[ncg].myEMF1[(2*dim)+1] == NULL) ath_error(
                      "[init_grid]:failed to allocate CGrid oxb myEMF1\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for oxb CGrid.myEMF1\n",n2z+1,n1z);
}
                    pG->CGrid[ncg].myEMF2[(2*dim)+1] =(Real**)calloc_2d_array(
                      n2z,n1z+1, sizeof(Real));
                    if(pG->CGrid[ncg].myEMF2[(2*dim)+1] == NULL) ath_error(
                      "[init_grid]:failed to allocate CGrid oxb myEMF2\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for oxb CGrid.myEMF2\n",n2z,n1z+1);
}
                  }
#endif /* MHD */

                }

              } /* end loop over boundary directions */
            } /* end if child Grid overlaps or touches */
          }}} /* end loops over all Grids in Domain[nl+1][ncd] */
        } /* end if child Domain */
      } /* end loop over all Domains at level [nl+1] */
    } /* end if Grid on this processor */
  }} /* end loops over all levels and domains per level */


/************************************/
for (nl=0; nl<(pM->NLevels); nl++){
for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
  if (pM->Domain[nl][nd].Grid != NULL) {
  pG = pM->Domain[nl][nd].Grid;          /* set ptr to this Grid */

printf("\nProcID=%d level=%d Domain=%d NCgrid=%d NmyCGrid=%d\n",
  myID_Comm_world,nl,nd,pG->NCGrid,pG->NmyCGrid);
for (i=0;i<pG->NCGrid; i++){
printf("CGrid=%d, [is,ie,js,je,ks,ke]=%d %d %d %d %d %d\n",i,
  pG->CGrid[i].ijks[0],pG->CGrid[i].ijke[0],
  pG->CGrid[i].ijks[1],pG->CGrid[i].ijke[1],
  pG->CGrid[i].ijks[2],pG->CGrid[i].ijke[2]);
printf("Child_ID=%d DomN=%d nWordsRC=%d nWordsP=%d\n",
  pG->CGrid[i].ID,pG->CGrid[i].DomN,pG->CGrid[i].nWordsRC,
  pG->CGrid[i].nWordsP);
}

}}}
/************************************/


/*--------------- Count number of parent Grids, and boundaries ---------------*/
/* Now we have to count the number of parent Grids, and fine/coarse boundaries
 * with parent Grids, including parent Grids on this and other processors,
 * before we can allocate the PGrid array.  This is a shortcoming of making
 * these structures 1D arrays. */
/* Loop over levels (except root level=0), and domains per level */

  for (nl=1; nl<(pM->NLevels); nl++){
  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {
      pD = (DomainS*)&(pM->Domain[nl][nd]);  /* set ptr to Domain */
      pG = pM->Domain[nl][nd].Grid;          /* set ptr to Grid */
      get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);

/* edges of this Domain */
      for (i=0; i<3; i++) {
        D1.ijkl[i] = pD->Disp[i];
        D1.ijkr[i] = pD->Disp[i] + pD->Nx[i];
      }

/* edges of this Grid */
      for (i=0; i<3; i++) {
        G1.ijkl[i] = pG->Disp[i];
        G1.ijkr[i] = pG->Disp[i] + pG->Nx[i];
      }

/* For this domain, find domain at last level that overlaps (parent Domain).
 * Multiply by two to check in units of this (not the parent) domain coordinates
 */

      for (npd=0; npd<(pM->DomainsPerLevel[nl-1]); npd++){
        pPD = (DomainS*)&(pM->Domain[nl-1][npd]); /* ptr to potential parent */

/* edges of potential parent Domain */
        for (i=0; i<3; i++) {
          D2.ijkl[i] = 2*(pPD->Disp[i]);
          D2.ijkr[i] = 2*(pPD->Disp[i] + pPD->Nx[i]);
        }

        isDOverlap = checkOverlap(&D1, &D2, &D3);
        if (isDOverlap == 1){

/*----------------------------------------------------------------------------*/
/* Found parent Domain that overlaps.  So on the parent Domain, find all the
 * Grids that overlap this Grid. */

          for (n=0; n<pPD->NGrid[2]; n++){
          for (m=0; m<pPD->NGrid[1]; m++){
          for (l=0; l<pPD->NGrid[0]; l++){

/* edges of parent Grid */
/* Multiply by two to check in units of this (not the parent) grid coord */
            for (i=0; i<3; i++) {
              G2.ijkl[i] = 2*(pPD->GData[n][m][l].Disp[i]);
              G2.ijkr[i] = 2*(pPD->GData[n][m][l].Disp[i] + 
                              pPD->GData[n][m][l].Nx[i]);
            }

            isGOverlap = checkOverlap(&G1, &G2, &G3);

/* If Grid overlaps, increment Parent counters */

            if (isGOverlap == 1){
              pG->NPGrid++;
              if (pPD->GData[n][m][l].ID_Comm_world == myID_Comm_world)
                pG->NmyPGrid++;
            }

          }}} /* end loops over [n,m,l]: all Grids in Domain[nl-1][npd] */
        }
      } /* end loop over all parent Domains at level [nl-1] */

/*-------------------------- Allocate PGrid array ----------------------------*/
/* Now we know how many parent Grids there are for the Grid in Domain[nd] at
 * level nl on both this and other processors. */

      if (pG->NPGrid > 0) {
        pG->PGrid =(GridOvrlpS*)calloc_1d_array(pG->NPGrid, sizeof(GridOvrlpS));
        if(pG->PGrid==NULL) ath_error("[init_grid]:failed to allocate PGrid\n");

        for (npg=0; npg<pG->NPGrid; npg++){
          for (dim=0; dim<6; dim++) {
            pG->PGrid[npg].myFlx[dim] = NULL;
#ifdef MHD
            pG->PGrid[npg].myEMF1[dim] = NULL;
            pG->PGrid[npg].myEMF2[dim] = NULL;
            pG->PGrid[npg].myEMF3[dim] = NULL;
#endif /* MHD */
          }
        }
      }

/*--------------------------- Fill in PGrid array ----------------------------*/
/* Repeat loop over all domains at last level, and all the logic to find
 * overlapping Grids, to fill in data about overlap regions in PGrid array */

      nMyPG = 0;
      nPG = pG->NmyPGrid;

      for (npd=0; npd<pM->DomainsPerLevel[nl-1]; npd++){
        pPD = (DomainS*)&(pM->Domain[nl-1][npd]);   /* ptr to potential parent*/

/* edges of potential parent Domain */
        for (i=0; i<3; i++) {
          D2.ijkl[i] = 2*(pPD->Disp[i]);
          D2.ijkr[i] = 2*(pPD->Disp[i] + pPD->Nx[i]);
        }

        isDOverlap = checkOverlap(&D1, &D2, &D3);
        if (isDOverlap == 1){

/*----------------------------------------------------------------------------*/
/* Found the Domain that overlaps, so on the parent Domain check if there is a
 * Grid that overlaps */

          for (n=0; n<pPD->NGrid[2]; n++){
          for (m=0; m<pPD->NGrid[1]; m++){
          for (l=0; l<pPD->NGrid[0]; l++){

/* edges of parent Grid */
/* Multiply by two to check in units of this (not the parent) grid coord */
            for (i=0; i<3; i++) {
              G2.ijkl[i] = 2*(pPD->GData[n][m][l].Disp[i]);
              G2.ijkr[i] = 2*(pPD->GData[n][m][l].Disp[i] +
                              pPD->GData[n][m][l].Nx[i]);
            }

            isGOverlap = checkOverlap(&G1, &G2, &G3);

            if (isGOverlap == 1){

/* If Grid OVERLAPS, then:
 * (1) fill-in data in PGrid array */
/* Index PGrid array so that parent Grids on this processor come first */

              if (pPD->GData[n][m][l].ID_Comm_world == myID_Comm_world) {
                npg=nMyPG;
                nMyPG++;
              } else {
                npg=nPG;
                nPG++;
              }

              pG->PGrid[npg].ijks[0] = G3.ijkl[0] - pG->Disp[0] + pG->is;
              pG->PGrid[npg].ijke[0] = G3.ijkr[0] - pG->Disp[0] + pG->is - 1;
              pG->PGrid[npg].ijks[1] = G3.ijkl[1] - pG->Disp[1] + pG->js;
              pG->PGrid[npg].ijke[1] = G3.ijkr[1] - pG->Disp[1] + pG->js - 1;
              pG->PGrid[npg].ijks[2] = G3.ijkl[2] - pG->Disp[2] + pG->ks;
              pG->PGrid[npg].ijke[2] = G3.ijkr[2] - pG->Disp[2] + pG->ks - 1;
  
              pG->PGrid[npg].DomN = npd;
              pG->PGrid[npg].ID = pPD->GData[n][m][l].ID_Comm_Children;

              n1z =    (pG->PGrid[npg].ijke[0]-pG->PGrid[npg].ijks[0] + 1)/2;
              n2z =MAX((pG->PGrid[npg].ijke[1]-pG->PGrid[npg].ijks[1] + 1)/2,1);
              n3z =MAX((pG->PGrid[npg].ijke[2]-pG->PGrid[npg].ijks[2] + 1)/2,1);
/*
              if (pG->Nx[1]>1) n2z /= 2;
              if (pG->Nx[2]>1) n3z /= 2;
*/
              pG->PGrid[npg].nWordsRC = n1z*n2z*n3z*(NVAR);
              pG->PGrid[npg].nWordsP  = 0;
#ifdef MHD
              if (pG->Nx[2]>1) {
                pG->PGrid[npg].nWordsRC += 
                  (n1z+1)*n2z*n3z + n1z*(n2z+1)*n3z + n1z*n2z*(n3z+1);
              } else {
                if (pG->Nx[1]>1) {
                  pG->PGrid[npg].nWordsRC += (n1z+1)*n2z + n1z*(n2z+1);
                }
              }
#endif /* MHD */

/* (2) If any edge of this Grid is at the edge of this Domain, then allocate
 * memory for fluxes and EMFs for Correction step, and count GZ data passed in
 * Prolongation step. */

              for (dim=0; dim<nDim; dim++){    /* only checks nDim directions */
                if (dim == 0) iGrid=myL;
                if (dim == 1) iGrid=myM;
                if (dim == 2) iGrid=myN;

/* inner x1/x2/x3 boundary */
/* First check that edge of this Grid is at edge of overlap (otherwise boundary
 * is between MPI blocks in parent, and is internal to child Grid).
 * Then check that edge of this Grid is not at l-edge of root (so that physical
 * BCs are applied), but is at l-edge of own Domain (so it is not an internal
 * MPI boundary on this Domain). */

                  if ((G1.ijkl[dim] == G3.ijkl[dim]) &&
                      (pD->Disp[dim] != 0) &&
                      (iGrid == 0)) {

                    if (dim == 0) {
                      n1z = G3.ijkr[1] - G3.ijkl[1];
                      n2z = G3.ijkr[2] - G3.ijkl[2];
                      n1p = MAX((n1z/2),1);
                      n2p = MAX((n2z/2),1);
                      if (pG->Nx[1] > 1) n1p += nghost + 2;
                      if (pG->Nx[2] > 1) n2p += nghost + 2;
                    }
                    if (dim == 1) {
                      n1z = G3.ijkr[0] - G3.ijkl[0];
                      n2z = G3.ijkr[2] - G3.ijkl[2];
                      n1p = (n1z/2) + nghost + 2;
                      n2p = MAX((n2z/2),1);
                      if (pG->Nx[2] > 1) n2p += nghost + 2;
                    }
                    if (dim == 2) {
                      n1z = G3.ijkr[0] - G3.ijkl[0];
                      n2z = G3.ijkr[1] - G3.ijkl[1];
                      n1p = (n1z/2) + nghost + 2;
                      n2p = (n2z/2) + nghost + 2;
                    }

                    pG->PGrid[npg].nWordsRC += 
                      MAX((n1z/2),1)*MAX((n2z/2),1)*(NVAR);
                    pG->PGrid[npg].nWordsP  += ((nghost/2)+2)*n1p*n2p*(NVAR);

/* Allocate memory for myFlx and my EMFS.  Note they have dimension of the
 * parent overlap on THIS Grid, which is 2x the transverse dimension of the
 * overlap on the parent Grid (the actual number of words sent). */

                    pG->PGrid[npg].myFlx[2*dim] = (ConsS**)calloc_2d_array(
                      n2z,n1z, sizeof(ConsS));
                    if(pG->PGrid[npg].myFlx[2*dim] == NULL) ath_error(
                      "[init_grid]:failed to allocate PGrid ixb myFlx\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for ixb PGrid.myFlx\n",n2z,n1z);
}
#ifdef MHD
                    pG->PGrid[npg].nWordsP += 6*((nghost/2)+2)*n1p*n2p;

                    if (pG->Nx[1] > 1 && dim != 2) {
                      pG->PGrid[npg].nWordsRC +=(n1z/2+1)*MAX((n2z/2),1);
                      pG->PGrid[npg].myEMF3[2*dim] = (Real**)calloc_2d_array(
                        n2z,n1z+1, sizeof(Real));
                      if(pG->PGrid[npg].myEMF3[2*dim]==NULL) ath_error(
                        "[init_grid]:failed to allocate PGrid ixb myEMF3\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for ixb PGrid.myEMF3\n",n2z,n1z+1);
}
                    }

                    if (pG->Nx[2] > 1  && dim == 0) {
                      pG->PGrid[npg].nWordsRC += (n1z/2)*(n2z/2+1);
                      pG->PGrid[npg].myEMF2[2*dim] = (Real**)calloc_2d_array(
                        n2z+1,n1z, sizeof(Real));
                      if(pG->PGrid[npg].myEMF2[2*dim]==NULL) ath_error(
                        "[init_grid]:failed to allocate PGrid ixb myEMF2\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for ixb PGrid.myEMF2\n",n2z+1,n1z);
}
                    }

                    if (pG->Nx[2] > 1  && dim == 1) {
                      pG->PGrid[npg].nWordsRC += (n1z/2)*(n2z/2+1);
                      pG->PGrid[npg].myEMF1[2*dim] = (Real**)calloc_2d_array(
                        n2z+1,n1z, sizeof(Real));
                      if(pG->PGrid[npg].myEMF1[2*dim]==NULL) ath_error(
                        "[init_grid]:failed to allocate PGrid ixb myEMF1\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for ixb PGrid.myEMF1\n",n2z+1,n1z);
}
                    }

                    if (pG->Nx[2] > 1  && dim == 2) {
                      pG->PGrid[npg].nWordsRC += (n1z/2)*(n2z/2+1) 
                                             + (n1z/2+1)*(n2z/2);
                      pG->PGrid[npg].myEMF1[2*dim] = (Real**)calloc_2d_array(
                        n2z+1,n1z, sizeof(Real));
                      if(pG->PGrid[npg].myEMF1[2*dim]==NULL) ath_error(
                        "[init_grid]:failed to allocate PGrid ixb myEMF1\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for ixb PGrid.myEMF1\n",n2z+1,n1z);
}
                      pG->PGrid[npg].myEMF2[2*dim] = (Real**)calloc_2d_array(
                        n2z,n1z+1, sizeof(Real));
                      if(pG->PGrid[npg].myEMF2[2*dim]==NULL) ath_error(
                        "[init_grid]:failed to allocate PGrid ixb myEMF2\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for ixb PGrid.myEMF2\n",n2z,n1z+1);
}
                    }
#endif /* MHD */

                }

/* outer x1/x2/x3 boundary */
/* First check that edge of this Grid is at edge of overlap (otherwise boundary
 * is between MPI blocks in parent, and is internal to child Grid).
 * Then check that edge of this Grid is not at r-edge of root (so that physical
 * BCs are applied), but is at r-edge of own Domain (so it is not an internal
 * MPI boundary on this Domain). */

                irefine = 1;
                for (i=1;i<=nl;i++) irefine *= 2; /* this level refinement */
                if ( (G1.ijkr[dim] == G3.ijkr[dim]) &&
                    ((pD->Disp[dim] + pD->Nx[dim])/irefine != pM->Nx[dim]) &&
                     (iGrid == (pD->NGrid[dim]-1)) ) {
  
                    if (dim == 0) {
                      n1z = G3.ijkr[1] - G3.ijkl[1];
                      n2z = G3.ijkr[2] - G3.ijkl[2];
                      n1p = MAX((n1z/2),1);
                      n2p = MAX((n2z/2),1);
                      if (pG->Nx[1] > 1) n1p += nghost + 2;
                      if (pG->Nx[2] > 1) n2p += nghost + 2;
                    }
                    if (dim == 1) {
                      n1z = G3.ijkr[0] - G3.ijkl[0];
                      n2z = G3.ijkr[2] - G3.ijkl[2];
                      n1p = (n1z/2) + nghost + 2;
                      n2p = MAX((n2z/2),1);
                      if (pG->Nx[2] > 1) n2p += nghost + 2;
                    }
                    if (dim == 2) {
                      n1z = G3.ijkr[0] - G3.ijkl[0];
                      n2z = G3.ijkr[1] - G3.ijkl[1];
                      n1p = (n1z/2) + nghost + 2;
                      n2p = (n2z/2) + nghost + 2;
                    }

                    pG->PGrid[npg].nWordsRC += 
                      MAX((n1z/2),1)*MAX((n2z/2),1)*(NVAR);
                    pG->PGrid[npg].nWordsP  += ((nghost/2)+2)*n1p*n2p*(NVAR);

/* Allocate memory for myFlx and my EMFS.  Note they have dimension of the
 * parent overlap on THIS Grid, which is 2x the transverse dimension of the
 * overlap on the parent Grid (the actual number of words sent). */

                    pG->PGrid[npg].myFlx[(2*dim)+1] =
                      (ConsS**)calloc_2d_array(n2z,n1z, sizeof(ConsS));
                    if(pG->PGrid[npg].myFlx[(2*dim)+1] == NULL) ath_error(
                      "[init_grid]:failed to allocate PGrid oxb myFlx\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for oxb PGrid.myFlx\n",n2z,n1z);
}
#ifdef MHD
                    pG->PGrid[npg].nWordsP += 6*((nghost/2)+2)*n1p*n2p;

                    if (pG->Nx[1] > 1 && dim != 2) {
                      pG->PGrid[npg].nWordsRC += (n1z/2+1)*MAX(n2z/2,1);
                      pG->PGrid[npg].myEMF3[(2*dim)+1] =(Real**)calloc_2d_array(
                        n2z,n1z+1, sizeof(Real));
                      if(pG->PGrid[npg].myEMF3[(2*dim)+1]==NULL) ath_error(
                        "[init_grid]:failed to allocate PGrid oxb myEMF3\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for oxb PGrid.myEMF3\n",n2z,n1z+1);
}
                    }

                    if (pG->Nx[2] > 1  && dim == 0) {
                      pG->PGrid[npg].nWordsRC += (n1z/2)*(n2z/2+1);
                      pG->PGrid[npg].myEMF2[(2*dim)+1] =(Real**)calloc_2d_array(
                        n2z+1,n1z, sizeof(Real));
                      if(pG->PGrid[npg].myEMF2[(2*dim)+1]==NULL) ath_error(
                        "[init_grid]:failed to allocate PGrid oxb myEMF2\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for oxb PGrid.myEMF2\n",n2z+1,n1z);
}
                    }

                    if (pG->Nx[2] > 1  && dim == 1) {
                      pG->PGrid[npg].nWordsRC += (n1z/2)*(n2z/2+1);
                      pG->PGrid[npg].myEMF1[(2*dim)+1] =(Real**)calloc_2d_array(
                        n2z+1,n1z, sizeof(Real));
                      if(pG->PGrid[npg].myEMF1[(2*dim)+1]==NULL) ath_error(
                        "[init_grid]:failed to allocate PGrid oxb myEMF1\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for oxb PGrid.myEMF1\n",n2z+1,n1z);
}
                    }

                    if (pG->Nx[2] > 1  && dim == 2) {
                      pG->PGrid[npg].nWordsRC += (n1z/2)*(n2z/2+1) 
                                               + (n1z/2+1)*(n2z/2);
                      pG->PGrid[npg].myEMF1[(2*dim)+1] =(Real**)calloc_2d_array(
                        n2z+1,n1z, sizeof(Real));
                      if(pG->PGrid[npg].myEMF1[(2*dim)+1]==NULL) ath_error(
                        "[init_grid]:failed to allocate PGrid oxb myEMF1\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for oxb PGrid.myEMF1\n",n2z+1,n1z);
}
                      pG->PGrid[npg].myEMF2[(2*dim)+1] =(Real**)calloc_2d_array(
                        n2z,n1z+1, sizeof(Real));
                      if(pG->PGrid[npg].myEMF2[(2*dim)+1]==NULL) ath_error(
                        "[init_grid]:failed to allocate PGrid oxb myEMF2\n");
if(myID_Comm_world==0){
printf("Allocated %d x %d array for oxb PGrid.myEMF2\n",n2z,n1z+1);
}
                    }
#endif /* MHD */

                }

              } /* end loop over boundary directions */
            } /* end if parent Grid overlaps or touches */
          }}} /* end loops over all Grids in Domain[nl-1][npd] */
        } /* end if parent Domain */
      } /* end loop over all Domains at level [nl-1] */
    } /* end if Grid on this processor */
  }} /* end loops over all levels and domains per level */

/************************************/
for (nl=0; nl<(pM->NLevels); nl++){
for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
  if (pM->Domain[nl][nd].Grid != NULL) {
  pG = pM->Domain[nl][nd].Grid;          /* set ptr to this Grid */

printf("\nProcID=%d level=%d Domain=%d NPgrid=%d NmyPGrid=%d\n",
  myID_Comm_world,nl,nd,pG->NPGrid,pG->NmyPGrid);
for (i=0;i<pG->NPGrid; i++){
printf("PGrid=%d, [is,ie,js,je,ks,ke]=%d %d %d %d %d %d\n",i,
  pG->PGrid[i].ijks[0],pG->PGrid[i].ijke[0],
  pG->PGrid[i].ijks[1],pG->PGrid[i].ijke[1],
  pG->PGrid[i].ijks[2],pG->PGrid[i].ijke[2]);
printf("Parent_ID=%d DomN=%d nWordsRC=%d nWordsP=%d\n",
  pG->PGrid[i].ID,pG->PGrid[i].DomN,pG->PGrid[i].nWordsRC,
  pG->PGrid[i].nWordsP);
}

}}}
/************************************/


#endif /* STATIC_MESH_REFINEMENT */

  return;

/*--- Error messages ---------------------------------------------------------*/

#ifdef CYLINDRICAL
  on_error15:
    free_1d_array(pG->ri);
  on_error14:
    free_1d_array(pG->r);
#endif
#ifdef SELF_GRAVITY
  on_error13:
    free_3d_array(pG->x3MassFlux);
  on_error12:
    free_3d_array(pG->x2MassFlux);
  on_error11:
    free_3d_array(pG->x1MassFlux);
  on_error10:
    free_3d_array(pG->Phi_old);
  on_error9:
    free_3d_array(pG->Phi);
#endif
#ifdef RESISTIVITY
  on_error7:
    free_3d_array(pG->eta_AD);
  on_error6:
    free_3d_array(pG->eta_Hall);
  on_error5:
    free_3d_array(pG->eta_Ohm);
#endif
#ifdef MHD
  on_error4:
    free_3d_array(pG->B3i);
  on_error3:
    free_3d_array(pG->B2i);
  on_error2:
    free_3d_array(pG->B1i);
#endif
  on_error1:
    free_3d_array(pG->U);
    ath_error("[init_grid]: Error allocating memory\n");
}

#ifdef STATIC_MESH_REFINEMENT
/*=========================== PRIVATE FUNCTIONS ==============================*/
/*----------------------------------------------------------------------------*/
/*! \fn int checkOverlap(SideS *pC1, SideS *pC2, SideS *pC3)
 *  \brief Checks if two cubes are overlapping.
 *
 *  - If yes returns true and sides of overlap region in Cube3
 *  - If no  returns false and -1 in Cube3
 *
 * Arguments are Side structures, containing indices of the 6 sides of cube */

int checkOverlap(SideS *pC1, SideS *pC2, SideS *pC3)
{
  int isOverlap=0;

  if (pC1->ijkl[0] < pC2->ijkr[0] && pC1->ijkr[0] > pC2->ijkl[0] &&
      pC1->ijkl[1] < pC2->ijkr[1] && pC1->ijkr[1] > pC2->ijkl[1] &&
      pC1->ijkl[2] < pC2->ijkr[2] && pC1->ijkr[2] > pC2->ijkl[2]) isOverlap=1;

  if (isOverlap==1) {
    pC3->ijkl[0] = MAX(pC1->ijkl[0], pC2->ijkl[0]);
    pC3->ijkr[0] = MIN(pC1->ijkr[0], pC2->ijkr[0]);
    pC3->ijkl[1] = MAX(pC1->ijkl[1], pC2->ijkl[1]);
    pC3->ijkr[1] = MIN(pC1->ijkr[1], pC2->ijkr[1]);
    pC3->ijkl[2] = MAX(pC1->ijkl[2], pC2->ijkl[2]);
    pC3->ijkr[2] = MIN(pC1->ijkr[2], pC2->ijkr[2]);
  } else {
    pC3->ijkl[0] = -1;
    pC3->ijkr[0] = -1;
    pC3->ijkl[1] = -1;
    pC3->ijkr[1] = -1;
    pC3->ijkl[2] = -1;
    pC3->ijkr[2] = -1;
  }

  return isOverlap;
}

/*----------------------------------------------------------------------------*/
/*! \fn int checkOverlapTouch(SideS *pC1, SideS *pC2, SideS *pC3)
 *  \brief Checks if two cubes are overlapping or touching.
 *
 *  - If yes returns true and sides of overlap region in Cube3
 *  - If no  returns false and -1 in Cube3
 *
 * Arguments are Side structures, containing indices of the 6 sides of cube */

int checkOverlapTouch(SideS *pC1, SideS *pC2, SideS *pC3)
{
  int isOverlap=0;

  if (pC1->ijkl[0] <= pC2->ijkr[0] && pC1->ijkr[0] >= pC2->ijkl[0] &&
      pC1->ijkl[1] <= pC2->ijkr[1] && pC1->ijkr[1] >= pC2->ijkl[1] &&
      pC1->ijkl[2] <= pC2->ijkr[2] && pC1->ijkr[2] >= pC2->ijkl[2]) isOverlap=1;

  if (isOverlap==1) {
    pC3->ijkl[0] = MAX(pC1->ijkl[0], pC2->ijkl[0]);
    pC3->ijkr[0] = MIN(pC1->ijkr[0], pC2->ijkr[0]);
    pC3->ijkl[1] = MAX(pC1->ijkl[1], pC2->ijkl[1]);
    pC3->ijkr[1] = MIN(pC1->ijkr[1], pC2->ijkr[1]);
    pC3->ijkl[2] = MAX(pC1->ijkl[2], pC2->ijkl[2]);
    pC3->ijkr[2] = MIN(pC1->ijkr[2], pC2->ijkr[2]);
  } else {
    pC3->ijkl[0] = -1;
    pC3->ijkr[0] = -1;
    pC3->ijkl[1] = -1;
    pC3->ijkr[1] = -1;
    pC3->ijkl[2] = -1;
    pC3->ijkr[2] = -1;
  }

  return isOverlap;
}
#endif /* STATIC_MESH_REFINEMENT */
