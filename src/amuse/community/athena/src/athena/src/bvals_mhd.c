#include "copyright.h"
/*============================================================================*/
/*! \file bvals_mhd.c
 *  \brief Sets boundary conditions (quantities in ghost zones) on each edge
 *   of a Grid for the MHD variables.
 *
 * PURPOSE: Sets boundary conditions (quantities in ghost zones) on each edge
 *   of a Grid for the MHD variables.  Each edge of a Grid represents either:
 *  - (1) a physical boundary at the edge of the Mesh; in which case BCs are
 *        specified by an integer flag input by user (or by user-defined BC
 *        function in the problem file)
 *  - (2) the boundary between Grids resulting from decomposition of a larger
 *        Domain using MPI; in which case BCs are handled by MPI calls
 *  - (3) an internal boundary between fine/coarse grid levels in a nested Mesh;
 *        in which case the BCs require a prolongation operator from the parent
 *
 *   This file contains functions that can handle the first two cases.  Case (3)
 *   is handled in the Prolongate function called from main loop.
 *   The naming convention of the integer flags for BCs is:
 *    -  bc_ix1 = Boundary Condition for Inner x1 (at i=is)
 *    -  bc_ox1 = Boundary Condition for Outer x1 (at i=ie)
 *   similarly for bc_ix2; bc_ox2; bc_ix3; bc_ox3
 *
 * For case (1) -- PHYSICAL BOUNDARIES
 *   The values of the integer flags (bc_ix1, etc.) are:
 *   - 1 = reflecting; 2 = outflow; 4 = periodic; 5 = conductor
 *   For flow-in bondaries (ghost zones set to pre-determined values), pointers
 *   to user-defined functions in the problem file are used. 
 *
 * For case (2) -- MPI BOUNDARIES
 *   We do the parallel synchronization by having every grid:
 *   - 1) Post non-blocking receives for data from both L and R Grids
 *   - 2) Pack and send data to the Grids on both L and R
 *   - 3) Check for receives and unpack data in order of first to finish
 *   If the Grid is at the edge of the Domain, we set BCs as in case (1) or (3).
 *
 * For case (3) -- INTERNAL GRID LEVEL BOUNDARIES
 *   This step is complicated and must be handled separately, in the function
 *   Prolongate() called from the main loop.  In the algorithm below, nothing is
 *   done for this case; the BCs are left to be set later.
 *
 * With SELF-GRAVITY: BCs for Phi are set independently of the MHD variables
 *   in a separate function bvals_grav(). 
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - bvals_mhd()      - calls appropriate functions to set ghost cells
 * - bvals_mhd_init() - sets function pointers used by bvals_mhd()
 * - bvals_mhd_fun()  - enrolls a pointer to a user-defined BC function
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - reflect_ix1()  - reflecting BCs at boundary ix1
 * - reflect_ox1()  - reflecting BCs at boundary ox1
 * - reflect_ix2()  - reflecting BCs at boundary ix2
 * - reflect_ox2()  - reflecting BCs at boundary ox2
 * - reflect_ix3()  - reflecting BCs at boundary ix3
 * - reflect_ox3()  - reflecting BCs at boundary ox3
 * - outflow_ix1()  - outflow BCs at boundary ix1
 * - outflow_ox1()  - outflow BCs at boundary ox1
 * - outflow_ix2()  - outflow BCs at boundary ix2
 * - outflow_ox2()  - outflow BCs at boundary ox2
 * - outflow_ix3()  - outflow BCs at boundary ix3
 * - outflow_ox3()  - outflow BCs at boundary ox3
 * - periodic_ix1() - periodic BCs at boundary ix1
 * - periodic_ox1() - periodic BCs at boundary ox1
 * - periodic_ix2() - periodic BCs at boundary ix2
 * - periodic_ox2() - periodic BCs at boundary ox2
 * - periodic_ix3() - periodic BCs at boundary ix3
 * - periodic_ox3() - periodic BCs at boundary ox3
 * - conduct_ix1()  - conducting BCs at boundary ix1
 * - conduct_ox1()  - conducting BCs at boundary ox1
 * - conduct_ix2()  - conducting BCs at boundary ix2
 * - conduct_ox2()  - conducting BCs at boundary ox2
 * - conduct_ix3()  - conducting BCs at boundary ix3
 * - conduct_ox3()  - conducting BCs at boundary ox3
 * - pack_ix1()     - pack data for MPI non-blocking send at ix1 boundary
 * - pack_ox1()     - pack data for MPI non-blocking send at ox1 boundary
 * - pack_ix2()     - pack data for MPI non-blocking send at ix2 boundary
 * - pack_ox2()     - pack data for MPI non-blocking send at ox2 boundary
 * - pack_ix3()     - pack data for MPI non-blocking send at ix3 boundary
 * - pack_ox3()     - pack data for MPI non-blocking send at ox3 boundary
 * - unpack_ix1()   - unpack data for MPI non-blocking receive at ix1 boundary
 * - unpack_ox1()   - unpack data for MPI non-blocking receive at ox1 boundary
 * - unpack_ix2()   - unpack data for MPI non-blocking receive at ix2 boundary
 * - unpack_ox2()   - unpack data for MPI non-blocking receive at ox2 boundary
 * - unpack_ix3()   - unpack data for MPI non-blocking receive at ix3 boundary
 * - unpack_ox3()   - unpack data for MPI non-blocking receive at ox3 boundary*/
/*============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef MPI_PARALLEL
/* MPI send and receive buffers */
static double **send_buf = NULL, **recv_buf = NULL;
static MPI_Request *recv_rq, *send_rq;
#endif /* MPI_PARALLEL */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   reflect_???()  - reflecting BCs at boundary ???
 *   outflow_???()  - outflow BCs at boundary ???
 *   periodic_???() - periodic BCs at boundary ???
 *   conduct_???()  - conducting BCs at boundary ???
 *   pack_???()     - pack data for MPI non-blocking send at ??? boundary
 *   unpack_???()   - unpack data for MPI non-blocking receive at ??? boundary
 *============================================================================*/

static void reflect_ix1(GridS *pG);
static void reflect_ox1(GridS *pG);
static void reflect_ix2(GridS *pG);
static void reflect_ox2(GridS *pG);
static void reflect_ix3(GridS *pG);
static void reflect_ox3(GridS *pG);

static void outflow_ix1(GridS *pG);
static void outflow_ox1(GridS *pG);
static void outflow_ix2(GridS *pG);
static void outflow_ox2(GridS *pG);
static void outflow_ix3(GridS *pG);
static void outflow_ox3(GridS *pG);

static void periodic_ix1(GridS *pG);
static void periodic_ox1(GridS *pG);
static void periodic_ix2(GridS *pG);
static void periodic_ox2(GridS *pG);
static void periodic_ix3(GridS *pG);
static void periodic_ox3(GridS *pG);

static void conduct_ix1(GridS *pG);
static void conduct_ox1(GridS *pG);
static void conduct_ix2(GridS *pG);
static void conduct_ox2(GridS *pG);
static void conduct_ix3(GridS *pG);
static void conduct_ox3(GridS *pG);

static void ProlongateLater(GridS *pG);

#ifdef MPI_PARALLEL
static void pack_ix1(GridS *pG);
static void pack_ox1(GridS *pG);
static void pack_ix2(GridS *pG);
static void pack_ox2(GridS *pG);
static void pack_ix3(GridS *pG);
static void pack_ox3(GridS *pG);

static void unpack_ix1(GridS *pG);
static void unpack_ox1(GridS *pG);
static void unpack_ix2(GridS *pG);
static void unpack_ox2(GridS *pG);
static void unpack_ix3(GridS *pG);
static void unpack_ox3(GridS *pG);
#endif /* MPI_PARALLEL */

/*=========================== PUBLIC FUNCTIONS ===============================*/

/*----------------------------------------------------------------------------*/
/*! \fn void bvals_mhd(DomainS *pD)
 *  \brief Calls appropriate functions to set ghost zones.  
 *
 *   The function
 *   pointers (*(pD->???_BCFun)) are set by bvals_init() to be either a
 *   user-defined function, or one of the functions corresponding to reflecting,
 *   periodic, or outflow.  If the left- or right-Grid ID numbers are >= 1
 *   (neighboring grids exist), then MPI calls are used.
 *
 * Order for updating boundary conditions must always be x1-x2-x3 in order to
 * fill the corner cells properly
 */

void bvals_mhd(DomainS *pD)
{
  GridS *pGrid = (pD->Grid);
#ifdef SHEARING_BOX
  int myL,myM,myN,BCFlag;
#endif
#ifdef MPI_PARALLEL
  int cnt, cnt2, cnt3, ierr, mIndex;
#endif /* MPI_PARALLEL */

/*--- Step 1. ------------------------------------------------------------------
 * Boundary Conditions in x1-direction */

  if (pGrid->Nx[0] > 1){

#ifdef MPI_PARALLEL
    cnt = nghost*(pGrid->Nx[1])*(pGrid->Nx[2])*(NVAR);
#ifdef MHD
    cnt2 = (pGrid->Nx[1] > 1) ? (pGrid->Nx[1] + 1) : 1;
    cnt3 = (pGrid->Nx[2] > 1) ? (pGrid->Nx[2] + 1) : 1;
    cnt += (nghost-1)*(pGrid->Nx[1])*(pGrid->Nx[2]);
    cnt += nghost*cnt2*(pGrid->Nx[2]);
    cnt += nghost*(pGrid->Nx[1])*cnt3;
#endif

/* MPI blocks to both left and right */
    if (pGrid->rx1_id >= 0 && pGrid->lx1_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix1(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox1(pGrid); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1(pGrid);
      if (mIndex == 1) unpack_ox1(pGrid);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1(pGrid);
      if (mIndex == 1) unpack_ox1(pGrid);

    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx1_id >= 0 && pGrid->lx1_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox1(pGrid); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix1_BCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox1(pGrid);

    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx1_id < 0 && pGrid->lx1_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix1(pGrid); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox1_BCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix1(pGrid);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx1_id < 0 && pGrid->lx1_id < 0) {
      (*(pD->ix1_BCFun))(pGrid);
      (*(pD->ox1_BCFun))(pGrid);
    } 

  }

/*--- Step 2. ------------------------------------------------------------------
 * Boundary Conditions in x2-direction */

  if (pGrid->Nx[1] > 1){

#ifdef MPI_PARALLEL
    cnt = (pGrid->Nx[0] + 2*nghost)*nghost*(pGrid->Nx[2])*(NVAR);
#ifdef MHD
    cnt3 = (pGrid->Nx[2] > 1) ? (pGrid->Nx[2] + 1) : 1;
    cnt += (pGrid->Nx[0] + 2*nghost - 1)*nghost*(pGrid->Nx[2]);
    cnt += (pGrid->Nx[0] + 2*nghost)*(nghost-1)*(pGrid->Nx[2]);
    cnt += (pGrid->Nx[0] + 2*nghost)*nghost*cnt3;
#endif

/* MPI blocks to both left and right */
    if (pGrid->rx2_id >= 0 && pGrid->lx2_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix2(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox2(pGrid); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2(pGrid);
      if (mIndex == 1) unpack_ox2(pGrid);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2(pGrid);
      if (mIndex == 1) unpack_ox2(pGrid);

    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx2_id >= 0 && pGrid->lx2_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox2(pGrid); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix2_BCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox2(pGrid);

    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx2_id < 0 && pGrid->lx2_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix2(pGrid); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox2_BCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix2(pGrid);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx2_id < 0 && pGrid->lx2_id < 0) {
      (*(pD->ix2_BCFun))(pGrid);
      (*(pD->ox2_BCFun))(pGrid);
    } 

/* shearing sheet BCs; function defined in problem generator.
 * Enroll outflow BCs if perdiodic BCs NOT selected.  This assumes the root
 * level grid is specified by the <domain1> block in the input file */
#ifdef SHEARING_BOX
    BCFlag = par_geti_def("domain1","bc_ix1",0);
    get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);
    if (myL == 0 && BCFlag == 4) {
      ShearingSheet_ix1(pD);
    }
    BCFlag = par_geti_def("domain1","bc_ox1",0);
    if (myL == ((pD->NGrid[0])-1) && BCFlag == 4) {
      ShearingSheet_ox1(pD);
    }
#endif

  }

/*--- Step 3. ------------------------------------------------------------------
 * Boundary Conditions in x3-direction */

  if (pGrid->Nx[2] > 1){

#ifdef MPI_PARALLEL
    cnt = (pGrid->Nx[0] + 2*nghost)*(pGrid->Nx[1] + 2*nghost)*nghost*(NVAR);
#ifdef MHD
    cnt += (pGrid->Nx[0] + 2*nghost - 1)*(pGrid->Nx[1] + 2*nghost)*nghost;
    cnt += (pGrid->Nx[0] + 2*nghost)*(pGrid->Nx[1] + 2*nghost - 1)*nghost;
    cnt += (pGrid->Nx[0] + 2*nghost)*(pGrid->Nx[1] + 2*nghost)*(nghost-1);
#endif

/* MPI blocks to both left and right */
    if (pGrid->rx3_id >= 0 && pGrid->lx3_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix3(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox3(pGrid); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3(pGrid);
      if (mIndex == 1) unpack_ox3(pGrid);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3(pGrid);
      if (mIndex == 1) unpack_ox3(pGrid);

    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx3_id >= 0 && pGrid->lx3_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox3(pGrid); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix3_BCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox3(pGrid);

    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx3_id < 0 && pGrid->lx3_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix3(pGrid); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox3_BCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix3(pGrid);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx3_id < 0 && pGrid->lx3_id < 0) {
      (*(pD->ix3_BCFun))(pGrid);
      (*(pD->ox3_BCFun))(pGrid);
    } 

  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void bvals_mhd_init(MeshS *pM)
 *  \brief Sets function pointers for physical boundaries during
 *   initialization, allocates memory for send/receive buffers with MPI.
 */

void bvals_mhd_init(MeshS *pM)
{
  GridS *pG;
  DomainS *pD;
  int i,nl,nd,irefine;
#ifdef MPI_PARALLEL
  int myL,myM,myN,l,m,n,nx1t,nx2t,nx3t,size;
  int x1cnt=0, x2cnt=0, x3cnt=0; /* Number of words passed in x1/x2/x3-dir. */
#endif /* MPI_PARALLEL */

/* Cycle through all the Domains that have active Grids on this proc */

  for (nl=0; nl<(pM->NLevels); nl++){
  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
  if (pM->Domain[nl][nd].Grid != NULL) {
    pD = (DomainS*)&(pM->Domain[nl][nd]);  /* ptr to Domain */
    pG = pM->Domain[nl][nd].Grid;          /* ptr to Grid */
    irefine = 1;
    for (i=1;i<=nl;i++) irefine *= 2;   /* C pow fn only takes doubles !! */
#ifdef MPI_PARALLEL
/* get (l,m,n) coordinates of Grid being updated on this processor */
    get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);
#endif /* MPI_PARALLEL */

/* Set function pointers for physical boundaries in x1-direction -------------*/

    if(pG->Nx[0] > 1) {

/*---- ix1 boundary ----------------------------------------------------------*/

      if(pD->ix1_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if(pD->Disp[0] != 0) {      
          pD->ix1_BCFun = ProlongateLater;

/* Domain is at L-edge of root Domain, but not R-edge and periodic BC  */
        } else {
          if(((pD->Disp[0] + pD->Nx[0])/irefine != pM->Nx[0]) && 
               pM->BCFlag_ix1 == 4) {
            ath_error("[bvals_init]:level=%d Domain touching ix1b but not ox1b and periodic BC not allowed\n",nl); 

/* Domain is at L-edge of root Domain */
          } else {                    
            switch(pM->BCFlag_ix1){

            case 1: /* Reflecting, B_normal=0 */
              pD->ix1_BCFun = reflect_ix1;
            break;

            case 2: /* Outflow */
              pD->ix1_BCFun = outflow_ix1;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->ix1_BCFun = periodic_ix1;
#ifdef MPI_PARALLEL
              if(pG->lx1_id < 0 && pD->NGrid[0] > 1){
                pG->lx1_id = pD->GData[myN][myM][pD->NGrid[0]-1].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;

            case 5: /* Reflecting, B_normal!=0 */
              pD->ix1_BCFun = reflect_ix1;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ix1=%d unknown\n",pM->BCFlag_ix1);
              exit(EXIT_FAILURE);
            }
          }
        }
      }

/*---- ox1 boundary ----------------------------------------------------------*/

      if(pD->ox1_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if((pD->Disp[0] + pD->Nx[0])/irefine != pM->Nx[0]) {
          pD->ox1_BCFun = ProlongateLater;

/* Domain is at R-edge of root Domain, but not L-edge and periodic BC */
        } else {
          if((pD->Disp[0] != 0) && (pM->BCFlag_ox1 == 4)) {      
            ath_error("[bvals_init]:level=%d Domain touching ox1b but not ix1b and periodic BC not allowed\n",nl); 


/* Domain is at R-edge of root Domain */
          } else {
            switch(pM->BCFlag_ox1){

            case 1: /* Reflecting, B_normal=0 */
              pD->ox1_BCFun = reflect_ox1;
            break;

            case 2: /* Outflow */
              pD->ox1_BCFun = outflow_ox1;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->ox1_BCFun = periodic_ox1;
#ifdef MPI_PARALLEL
              if(pG->rx1_id < 0 && pD->NGrid[0] > 1){
                pG->rx1_id = pD->GData[myN][myM][0].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;

            case 5: /* Reflecting, B_normal!=0 */
              pD->ox1_BCFun = reflect_ox1;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ox1=%d unknown\n",pM->BCFlag_ox1);
              exit(EXIT_FAILURE);
            }
          }
        }
      }
    }

/* Set function pointers for physical boundaries in x2-direction -------------*/

    if(pG->Nx[1] > 1) {

/*---- ix2 boundary ----------------------------------------------------------*/

      if(pD->ix2_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if(pD->Disp[1] != 0) {
          pD->ix2_BCFun = ProlongateLater;

/* Domain is at L-edge of root Domain, but not R-edge and periodic BC  */
        } else {
          if(((pD->Disp[1] + pD->Nx[1])/irefine != pM->Nx[1]) &&
               pM->BCFlag_ix2 == 4) {
            ath_error("[bvals_init]:level=%d Domain touching ix2b but not ox2b and periodic BC not allowed\n",nl); 


/* Domain is at L-edge of root Domain */
          } else {
            switch(pM->BCFlag_ix2){

            case 1: /* Reflecting, B_normal=0 */
              pD->ix2_BCFun = reflect_ix2;
            break;

            case 2: /* Outflow */
              pD->ix2_BCFun = outflow_ix2;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->ix2_BCFun = periodic_ix2;
#ifdef MPI_PARALLEL
              if(pG->lx2_id < 0 && pD->NGrid[1] > 1){
                pG->lx2_id = pD->GData[myN][pD->NGrid[1]-1][myL].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;
  
            case 5: /* Reflecting, B_normal!=0 */
              pD->ix2_BCFun = reflect_ix2;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ix2=%d unknown\n",pM->BCFlag_ix2);
              exit(EXIT_FAILURE);
            }
          }
        }
      }

/*---- ox2 boundary ----------------------------------------------------------*/

      if(pD->ox2_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if((pD->Disp[1] + pD->Nx[1])/irefine != pM->Nx[1]) {
          pD->ox2_BCFun = ProlongateLater;

/* Domain is at R-edge of root Domain, but not L-edge and periodic BC */
        } else {
          if((pD->Disp[1] != 0) && (pM->BCFlag_ox2 == 4)) {
            ath_error("[bvals_init]:level=%d Domain touching ox2b but not ix2b and periodic BC not allowed\n",nl); 

/* Domain is at R-edge of root Domain */
          } else {
            switch(pM->BCFlag_ox2){

            case 1: /* Reflecting, B_normal=0 */
              pD->ox2_BCFun = reflect_ox2;
            break;

            case 2: /* Outflow */
              pD->ox2_BCFun = outflow_ox2;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->ox2_BCFun = periodic_ox2;
#ifdef MPI_PARALLEL
              if(pG->rx2_id < 0 && pD->NGrid[1] > 1){
                pG->rx2_id = pD->GData[myN][0][myL].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;

            case 5: /* Reflecting, B_normal!=0 */
              pD->ox2_BCFun = reflect_ox2;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ox2=%d unknown\n",pM->BCFlag_ox2);
              exit(EXIT_FAILURE);
            }
          }
        }
      }
    }

/* Set function pointers for physical boundaries in x3-direction -------------*/

    if(pG->Nx[2] > 1) {

/*---- ix3 boundary ----------------------------------------------------------*/

      if(pD->ix3_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if(pD->Disp[2] != 0) {
          pD->ix3_BCFun = ProlongateLater;

/* Domain is at L-edge of root Domain, but not R-edge and periodic BC  */
        } else {
          if(((pD->Disp[2] + pD->Nx[2])/irefine != pM->Nx[2]) &&
               pM->BCFlag_ix3 == 4) {
            ath_error("[bvals_init]:level=%d Domain touching ix3b but not ox3b and periodic BC not allowed\n",nl); 

/* Domain is at L-edge of root Domain */
          } else {
            switch(pM->BCFlag_ix3){

            case 1: /* Reflecting, B_normal=0 */
              pD->ix3_BCFun = reflect_ix3;
            break;

            case 2: /* Outflow */
              pD->ix3_BCFun = outflow_ix3;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->ix3_BCFun = periodic_ix3;
#ifdef MPI_PARALLEL
              if(pG->lx3_id < 0 && pD->NGrid[2] > 1){
                pG->lx3_id = pD->GData[pD->NGrid[2]-1][myM][myL].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;

            case 5: /* Reflecting, B_normal!=0 */
              pD->ix3_BCFun = reflect_ix3;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ix3=%d unknown\n",pM->BCFlag_ix3);
              exit(EXIT_FAILURE);
            }
          }
        }
      }

/*---- ox3 boundary ----------------------------------------------------------*/

      if(pD->ox3_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if((pD->Disp[2] + pD->Nx[2])/irefine != pM->Nx[2]) {
          pD->ox3_BCFun = ProlongateLater;

/* Domain is at R-edge of root Domain, but not L-edge and periodic BC */
        } else {
          if((pD->Disp[2] != 0) && (pM->BCFlag_ox3 == 4)) {
            ath_error("[bvals_init]:level=%d Domain touching ox3b but not ix3b and periodic BC not allowed\n",nl); 

/* Domain is at R-edge of root Domain */
          } else {
            switch(pM->BCFlag_ox3){

            case 1: /* Reflecting, B_normal=0 */
              pD->ox3_BCFun = reflect_ox3;
            break;

            case 2: /* Outflow */
              pD->ox3_BCFun = outflow_ox3;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->ox3_BCFun = periodic_ox3;
#ifdef MPI_PARALLEL
              if(pG->rx3_id < 0 && pD->NGrid[2] > 1){
                pG->rx3_id = pD->GData[0][myM][myL].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;

            case 5: /* Reflecting, B_normal!=0 */
              pD->ox3_BCFun = reflect_ox3;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ox3=%d unknown\n",pM->BCFlag_ox3);
              exit(EXIT_FAILURE);
            }
          }
        }
      }
    }

/* Figure out largest size needed for send/receive buffers with MPI ----------*/

#ifdef MPI_PARALLEL

    for (n=0; n<(pD->NGrid[2]); n++){
    for (m=0; m<(pD->NGrid[1]); m++){
      for (l=0; l<(pD->NGrid[0]); l++){

/* x1cnt is surface area of x1 faces */
	if(pD->NGrid[0] > 1){
	  nx2t = pD->GData[n][m][l].Nx[1];
	  if(nx2t > 1) nx2t += 1;

	  nx3t = pD->GData[n][m][l].Nx[2];
	  if(nx3t > 1) nx3t += 1;

          if(nx2t*nx3t > x1cnt) x1cnt = nx2t*nx3t;
	}

/* x2cnt is surface area of x2 faces */
	if(pD->NGrid[1] > 1){
	  nx1t = pD->GData[n][m][l].Nx[0];
	  if(nx1t > 1) nx1t += 2*nghost;

	  nx3t = pD->GData[n][m][l].Nx[2];
	  if(nx3t > 1) nx3t += 1;

          if(nx1t*nx3t > x2cnt) x2cnt = nx1t*nx3t;
	}

/* x3cnt is surface area of x3 faces */
	if(pD->NGrid[2] > 1){
	  nx1t = pD->GData[n][m][l].Nx[0];
	  if(nx1t > 1) nx1t += 2*nghost;

	  nx2t = pD->GData[n][m][l].Nx[1];
	  if(nx2t > 1) nx2t += 2*nghost;

          if(nx1t*nx2t > x3cnt) x3cnt = nx1t*nx2t;
	}
      }
    }}
#endif /* MPI_PARALLEL */

  }}}  /* End loop over all Domains with active Grids -----------------------*/

#ifdef MPI_PARALLEL
/* Allocate memory for send/receive buffers and MPI_Requests */

  size = x1cnt > x2cnt ? x1cnt : x2cnt;
  size = x3cnt >  size ? x3cnt : size;

#ifdef MHD
  size *= nghost*((NVAR)+3);
#else
  size *= nghost*(NVAR);
#endif

  if (size > 0) {
    if((send_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      ath_error("[bvals_init]: Failed to allocate send buffer\n");

    if((recv_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      ath_error("[bvals_init]: Failed to allocate recv buffer\n");
  }

  if((recv_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    ath_error("[bvals_init]: Failed to allocate recv MPI_Request array\n");
  if((send_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    ath_error("[bvals_init]: Failed to allocate send MPI_Request array\n");

#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void bvals_mhd_fun(DomainS *pD, enum BCDirection dir, VGFun_t prob_bc)
 *  \brief Sets function ptrs for user-defined BCs.
 */

void bvals_mhd_fun(DomainS *pD, enum BCDirection dir, VGFun_t prob_bc)
{
  switch(dir){
  case left_x1:
    pD->ix1_BCFun = prob_bc;
    break;
  case right_x1:
    pD->ox1_BCFun = prob_bc;
    break;
  case left_x2:
    pD->ix2_BCFun = prob_bc;
    break;
  case right_x2:
    pD->ox2_BCFun = prob_bc;
    break;
  case left_x3:
    pD->ix3_BCFun = prob_bc;
    break;
  case right_x3:
    pD->ox3_BCFun = prob_bc;
    break;
  default:
    ath_perr(-1,"[bvals_fun]: Unknown direction = %d\n",dir);
    exit(EXIT_FAILURE);
  }
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/* Following are the functions:
 *   reflecting_???:   where ???=[ix1,ox1,ix2,ox2,ix3,ox3]
 *   outflow_???
 *   periodic_???
 *   conduct_???
 *   pack_???
 *   unpack_???
 */

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix1(GridS *pGrid)
 *  \brief  REFLECTING boundary conditions, Inner x1 boundary (bc_ix1=1) */

static void reflect_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju,ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i]    =  pGrid->U[k][j][is+(i-1)];
        pGrid->U[k][j][is-i].M1 = -pGrid->U[k][j][is-i].M1; /* reflect 1-mom. */
#ifdef MHD
        pGrid->U[k][j][is-i].B1c= -pGrid->U[k][j][is-i].B1c;/* reflect 1-fld. */
#endif
      }
    }
  }

#ifdef MHD
/* reflect normal component of B field, B1i not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      pGrid->B1i[k][j][is] = 0.0;
      for (i=1; i<=nghost-1; i++) {
        pGrid->B1i[k][j][is-i]   = -pGrid->B1i[k][j][is+i];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][is+(i-1)];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][is+(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox1(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Outer x1 boundary (bc_ox1=1). */

static void reflect_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju,ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i]    =  pGrid->U[k][j][ie-(i-1)];
        pGrid->U[k][j][ie+i].M1 = -pGrid->U[k][j][ie+i].M1; /* reflect 1-mom. */
#ifdef MHD
        pGrid->U[k][j][ie+i].B1c= -pGrid->U[k][j][ie+i].B1c;/* reflect 1-fld. */
#endif
      }
    }
  }

#ifdef MHD
/* reflect normal component of B field */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      pGrid->B1i[k][j][ie+1] = 0.0;
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i] = -pGrid->B1i[k][j][ie-(i-2)];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][ie-(i-1)];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][ie-(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix2(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Inner x2 boundary (bc_ix2=1) */

static void reflect_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][js-j][i]    =  pGrid->U[k][js+(j-1)][i];
        pGrid->U[k][js-j][i].M2 = -pGrid->U[k][js-j][i].M2; /* reflect 2-mom. */
#ifdef MHD
        pGrid->U[k][js-j][i].B2c= -pGrid->U[k][js-j][i].B2c;/* reflect 2-fld. */
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[k][js-j][i] = pGrid->B1i[k][js+(j-1)][i];
      }
    }
  }

/* reflect normal component of B field, B2i not set at j=js-nghost */
  for (k=ks; k<=ke; k++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      pGrid->B2i[k][js][i] = 0.0;
    }
    for (j=1; j<=nghost-1; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][js-j][i] = -pGrid->B2i[k][js+j][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][js-j][i] = pGrid->B3i[k][js+(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox2(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Outer x2 boundary (bc_ox2=1) */

static void reflect_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][je+j][i]    =  pGrid->U[k][je-(j-1)][i];
        pGrid->U[k][je+j][i].M2 = -pGrid->U[k][je+j][i].M2; /* reflect 2-mom. */
#ifdef MHD
        pGrid->U[k][je+j][i].B2c= -pGrid->U[k][je+j][i].B2c;/* reflect 2-fld. */
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[k][je+j][i] = pGrid->B1i[k][je-(j-1)][i];
      }
    }
  }

/* reflect normal component of B field */
  for (k=ks; k<=ke; k++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      pGrid->B2i[k][je+1][i] = 0.0;
    }
    for (j=2; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][je+j][i] = -pGrid->B2i[k][je-(j-2)][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][je+j][i] = pGrid->B3i[k][je-(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix3(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Inner x3 boundary (bc_ix3=1) */

static void reflect_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ks-k][j][i]    =  pGrid->U[ks+(k-1)][j][i];
        pGrid->U[ks-k][j][i].M3 = -pGrid->U[ks-k][j][i].M3; /* reflect 3-mom. */
#ifdef MHD
        pGrid->U[ks-k][j][i].B3c= -pGrid->U[ks-k][j][i].B3c;/* reflect 3-fld.*/
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ks+(k-1)][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ks-k][j][i] = pGrid->B2i[ks+(k-1)][j][i];
      }
    }
  }

/* reflect normal component of B field, B3i not set at k=ks-nghost */
  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      pGrid->B3i[ks][j][i] = 0.0;
    }
  }
  for (k=1; k<=nghost-1; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ks-k][j][i]   = -pGrid->B3i[ks+k][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox3(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Outer x3 boundary (bc_ox3=1) */

static void reflect_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ke+k][j][i]    =  pGrid->U[ke-(k-1)][j][i];
        pGrid->U[ke+k][j][i].M3 = -pGrid->U[ke+k][j][i].M3; /* reflect 3-mom. */
#ifdef MHD
        pGrid->U[ke+k][j][i].B3c= -pGrid->U[ke+k][j][i].B3c;/* reflect 3-fld. */
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke-(k-1)][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke-(k-1)][j][i];
      }
    }
  }

/* reflect normal component of B field */
  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      pGrid->B3i[ke+1][j][i] = 0.0;
    }
  }
  for (k=2; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ke+k][j][i]   = -pGrid->B3i[ke-(k-2)][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ix1(GridS *pGrid)
 *  \brief OUTFLOW boundary condition, Inner x1 boundary (bc_ix1=2) */

static void outflow_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i] = pGrid->U[k][j][is];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost-1; i++) {
        pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j][is];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][is];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][is];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ox1(GridS *pGrid)
 *  \brief OUTFLOW boundary conditions, Outer x1 boundary (bc_ox1=2) */

static void outflow_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][ie];
      }
    }
  }

#ifdef MHD
/* i=ie+1 is not a boundary condition for the interface field B1i */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][ie];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][ie];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][ie];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ix2(GridS *pGrid)
 *  \brief OUTFLOW boundary conditions, Inner x2 boundary (bc_ix2=2) */

static void outflow_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][js-j][i] = pGrid->U[k][js][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[k][js-j][i] = pGrid->B1i[k][js][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost-1; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][js-j][i] = pGrid->B2i[k][js][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][js-j][i] = pGrid->B3i[k][js][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ox2(GridS *pGrid)
 *  \brief OUTFLOW boundary conditions, Outer x2 boundary (bc_ox2=2) */

static void outflow_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][je+j][i] = pGrid->U[k][je][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[k][je+j][i] = pGrid->B1i[k][je][i];
      }
    }
  }

/* j=je+1 is not a boundary condition for the interface field B2i */
  for (k=ks; k<=ke; k++) {
    for (j=2; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][je+j][i] = pGrid->B2i[k][je][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][je+j][i] = pGrid->B3i[k][je][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ix3(GridS *pGrid)
 *  \brief OUTFLOW boundary conditions, Inner x3 boundary (bc_ix3=2) */

static void outflow_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ks-k][j][i] = pGrid->U[ks][j][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ks][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ks-k][j][i] = pGrid->B2i[ks][j][i];
      }
    }
  }

/* B3i is not set at k=ks-nghost */
  for (k=1; k<=nghost-1; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ks-k][j][i] = pGrid->B3i[ks][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ox3(GridS *pGrid)
 *  \brief OUTFLOW boundary conditions, Outer x3 boundary (bc_ox3=2) */

static void outflow_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ke+k][j][i] = pGrid->U[ke][j][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke][j][i];
      }
    }
  }

/* k=ke+1 is not a boundary condition for the interface field B3i */
  for (k=2; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ke+k][j][i] = pGrid->B3i[ke][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix1(GridS *pGrid)
 *  \brief PERIODIC boundary conditions, Inner x1 boundary (bc_ix1=4) */

static void periodic_ix1(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i] = pGrid->U[k][j][ie-(i-1)];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost-1; i++) {
        pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j][ie-(i-1)];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][ie-(i-1)];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][ie-(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox1(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Outer x1 boundary (bc_ox1=4) */

static void periodic_ox1(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][is+(i-1)];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=ie+1 */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][is+(i-1)];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][is+(i-1)];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][is+(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix2(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Inner x2 boundary (bc_ix2=4) */

static void periodic_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][js-j][i] = pGrid->U[k][je-(j-1)][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[k][js-j][i] = pGrid->B1i[k][je-(j-1)][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost-1; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][js-j][i] = pGrid->B2i[k][je-(j-1)][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][js-j][i] = pGrid->B3i[k][je-(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox2(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Outer x2 boundary (bc_ox2=4) */

static void periodic_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][je+j][i] = pGrid->U[k][js+(j-1)][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[k][je+j][i] = pGrid->B1i[k][js+(j-1)][i];
      }
    }
  }

/* B2i is not set at j=je+1 */
  for (k=ks; k<=ke; k++) {
    for (j=2; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][je+j][i] = pGrid->B2i[k][js+(j-1)][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][je+j][i] = pGrid->B3i[k][js+(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix3(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Inner x3 boundary (bc_ix3=4) */

static void periodic_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ks-k][j][i] = pGrid->U[ke-(k-1)][j][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ke-(k-1)][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ks-k][j][i] = pGrid->B2i[ke-(k-1)][j][i];
      }
    }
  }

/* B3i is not set at k=ks-nghost */
  for (k=1; k<=nghost-1; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ks-k][j][i] = pGrid->B3i[ke-(k-1)][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox3(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Outer x3 boundary (bc_ox3=4) */

static void periodic_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ke+k][j][i] = pGrid->U[ks+(k-1)][j][i];
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ks+(k-1)][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ke+k][j][i] = pGrid->B2i[ks+(k-1)][j][i];
      }
    }
  }

/* B3i is not set at k=ke+1 */
  for (k=2; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ke+k][j][i] = pGrid->B3i[ks+(k-1)][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ix1(GridS *pGrid)
 *  \brief CONDUCTOR boundary conditions, Inner x1 boundary (bc_ix1=5) */

static void conduct_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju,ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i]    =  pGrid->U[k][j][is+(i-1)];
        pGrid->U[k][j][is-i].M1 = -pGrid->U[k][j][is-i].M1; /* reflect 1-mom. */
#ifdef MHD
        pGrid->U[k][j][is-i].B2c= -pGrid->U[k][j][is-i].B2c;/* reflect fld */
        pGrid->U[k][j][is-i].B3c= -pGrid->U[k][j][is-i].B3c;
#endif
      }
    }
  }

#ifdef MHD
/* B1i not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost-1; i++) {
        pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j][is+i];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i]   = -pGrid->B2i[k][j][is+(i-1)];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i]   = -pGrid->B3i[k][j][is+(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ox1(GridS *pGrid)
 *  \brief CONDUCTOR boundary conditions, Outer x1 boundary (bc_ox1=5) */

static void conduct_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju,ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i]    =  pGrid->U[k][j][ie-(i-1)];
        pGrid->U[k][j][ie+i].M1 = -pGrid->U[k][j][ie+i].M1; /* reflect 1-mom. */
#ifdef MHD
        pGrid->U[k][j][ie+i].B2c= -pGrid->U[k][j][ie+i].B2c;/* reflect fld */
        pGrid->U[k][j][ie+i].B3c= -pGrid->U[k][j][ie+i].B3c;
#endif
      }
    }
  }

#ifdef MHD
/* i=ie+1 is not set for the interface field B1i */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][ie-(i-2)];
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i]   = -pGrid->B2i[k][j][ie-(i-1)];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i]   = -pGrid->B3i[k][j][ie-(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ix2(GridS *pGrid)
 *  \brief CONDUCTOR boundary conditions, Inner x2 boundary (bc_ix2=5) */

static void conduct_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][js-j][i]    =  pGrid->U[k][js+(j-1)][i];
        pGrid->U[k][js-j][i].M2 = -pGrid->U[k][js-j][i].M2; /* reflect 2-mom. */
#ifdef MHD
        pGrid->U[k][js-j][i].B1c= -pGrid->U[k][js-j][i].B1c;/* reflect fld */
        pGrid->U[k][js-j][i].B3c= -pGrid->U[k][js-j][i].B3c;
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B1i[k][js-j][i] = -pGrid->B1i[k][js+(j-1)][i];
      }
    }
  }

/* B2i not set at j=js-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost-1; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][js-j][i] = pGrid->B2i[k][js+j][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][js-j][i] = -pGrid->B3i[k][js+(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ox2(GridS *pGrid)
 *  \brief CONDUCTOR boundary conditions, Outer x2 boundary (bc_ox2=5) */

static void conduct_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[k][je+j][i]    =  pGrid->U[k][je-(j-1)][i];
        pGrid->U[k][je+j][i].M2 = -pGrid->U[k][je+j][i].M2; /* reflect 2-mom. */
#ifdef MHD
        pGrid->U[k][je+j][i].B1c= -pGrid->U[k][je+j][i].B1c;/* reflect fld */
        pGrid->U[k][je+j][i].B3c= -pGrid->U[k][je+j][i].B3c;
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[k][je+j][i]   = -pGrid->B1i[k][je-(j-1)][i];
      }
    }
  }

/* j=je+1 is not set for the interface field B2i */
  for (k=ks; k<=ke; k++) {
    for (j=2; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[k][je+j][i]   = pGrid->B2i[k][je-(j-2)][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][je+j][i]   = -pGrid->B3i[k][je-(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ix3(GridS *pGrid)
 *  \brief CONDUCTOR boundary conditions, Inner x3 boundary (bc_ix3=5) */

static void conduct_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ks-k][j][i]    =  pGrid->U[ks+(k-1)][j][i];
        pGrid->U[ks-k][j][i].M3 = -pGrid->U[ks-k][j][i].M3; /* reflect 3-mom. */
#ifdef MHD
        pGrid->U[ks-k][j][i].B1c= -pGrid->U[ks-k][j][i].B1c;/* reflect fld */
        pGrid->U[ks-k][j][i].B2c= -pGrid->U[ks-k][j][i].B2c;
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ks-k][j][i]   = -pGrid->B1i[ks+(k-1)][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ks-k][j][i]   = -pGrid->B2i[ks+(k-1)][j][i];
      }
    }
  }

/* B3i is not set at k=ks-nghost */
  for (k=1; k<=nghost-1; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ks-k][j][i] = pGrid->B3i[ks+k][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void conduct_ox3(GridS *pGrid)
 *  \brief CONDUCTOR boundary conditions, Outer x3 boundary (bc_ox3=5) */

static void conduct_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ke+k][j][i]    =  pGrid->U[ke-(k-1)][j][i];
        pGrid->U[ke+k][j][i].M3 = -pGrid->U[ke+k][j][i].M3; /* reflect 3-mom. */
#ifdef MHD
        pGrid->U[ke+k][j][i].B1c= -pGrid->U[ke+k][j][i].B1c;/* reflect fld */
        pGrid->U[ke+k][j][i].B2c= -pGrid->U[ke+k][j][i].B2c;
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ke+k][j][i]   = -pGrid->B1i[ke-(k-1)][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ke+k][j][i]   = -pGrid->B2i[ke-(k-1)][j][i];
      }
    }
  }

/* k=ke+1 is not set for the interface field B3i */
  for (k=2; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ke+k][j][i]   = pGrid->B3i[ke-(k-2)][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void ProlongateLater(GridS *pGrid)
 *  \brief PROLONGATION boundary conditions.  
 *
 *  Nothing is actually done here, the
 * prolongation is actually handled in ProlongateGhostZones in main loop, so
 * this is just a NoOp Grid function.  */

static void ProlongateLater(GridS *pGrid)
{
  return;
}

#ifdef MPI_PARALLEL  /* This ifdef wraps the next 12 funs; ~800 lines */
/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ix1(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Inner x1 boundary */

static void pack_ix1(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif
#if (NSCALARS > 0)
  int n;
#endif
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=is+(nghost-1); i++){
        *(pSnd++) = pG->U[k][j][i].d;
        *(pSnd++) = pG->U[k][j][i].M1;
        *(pSnd++) = pG->U[k][j][i].M2;
        *(pSnd++) = pG->U[k][j][i].M3;
#ifndef BAROTROPIC
        *(pSnd++) = pG->U[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        *(pSnd++) = pG->U[k][j][i].B1c;
        *(pSnd++) = pG->U[k][j][i].B2c;
        *(pSnd++) = pG->U[k][j][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) *(pSnd++) = pG->U[k][j][i].s[n];
#endif
      }
    }
  }

#ifdef MHD
/* B1i at i=is maps to B1i at i=ie+1 and is not passed */
  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is+1; i<=is+(nghost-1); i++){
        *(pSnd++) = pG->B1i[k][j][i];
      }
    }
  }

  if (pG->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++){
    for (j=js; j<=ju; j++){
      for (i=is; i<=is+(nghost-1); i++){
        *(pSnd++) = pG->B2i[k][j][i];
      }
    }
  }

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=is+(nghost-1); i++){
        *(pSnd++) = pG->B3i[k][j][i];
      }
    }
  }
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ox1(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Outer x1 boundary */

static void pack_ox1(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif
#if (NSCALARS > 0)
  int n;
#endif
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=ie-(nghost-1); i<=ie; i++){
        *(pSnd++) = pG->U[k][j][i].d;
        *(pSnd++) = pG->U[k][j][i].M1;
        *(pSnd++) = pG->U[k][j][i].M2;
        *(pSnd++) = pG->U[k][j][i].M3;
#ifndef BAROTROPIC
        *(pSnd++) = pG->U[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        *(pSnd++) = pG->U[k][j][i].B1c;
        *(pSnd++) = pG->U[k][j][i].B2c;
        *(pSnd++) = pG->U[k][j][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) *(pSnd++) = pG->U[k][j][i].s[n];
#endif
      }
    }
  }

#ifdef MHD
/* B1i at i=ie-(nghost-1) maps to B1i at i=is-nghost and is not passed */
  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=ie-(nghost-2); i<=ie; i++){
        *(pSnd++) = pG->B1i[k][j][i];
      }
    }
  }

  if (pG->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++){
    for (j=js; j<=ju; j++){
      for (i=ie-(nghost-1); i<=ie; i++){
        *(pSnd++) = pG->B2i[k][j][i];
      }
    }
  }

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++){
    for (j=js; j<=je; j++){
      for (i=ie-(nghost-1); i<=ie; i++){
        *(pSnd++) = pG->B3i[k][j][i];
      }
    }
  }
#endif /* MHD */
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ix2(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Inner x2 boundary */

static void pack_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif
#if (NSCALARS > 0)
  int n;
#endif
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=js+(nghost-1); j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        *(pSnd++) = pG->U[k][j][i].d;
        *(pSnd++) = pG->U[k][j][i].M1;
        *(pSnd++) = pG->U[k][j][i].M2;
        *(pSnd++) = pG->U[k][j][i].M3;
#ifndef BAROTROPIC
        *(pSnd++) = pG->U[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        *(pSnd++) = pG->U[k][j][i].B1c;
        *(pSnd++) = pG->U[k][j][i].B2c;
        *(pSnd++) = pG->U[k][j][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) *(pSnd++) = pG->U[k][j][i].s[n];
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=js+(nghost-1); j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        *(pSnd++) = pG->B1i[k][j][i];
      }
    }
  }

/* B2i at j=js maps to B2i at j=je+1 and is not passed */
  for (k=ks; k<=ke; k++) {
    for (j=js+1; j<=js+(nghost-1); j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        *(pSnd++) = pG->B2i[k][j][i];
      }
    }
  }

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=js+(nghost-1); j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        *(pSnd++) = pG->B3i[k][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ox2(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Outer x2 boundary */

static void pack_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif
#if (NSCALARS > 0)
  int n;
#endif
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

  for (k=ks; k<=ke; k++){
    for (j=je-(nghost-1); j<=je; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        *(pSnd++) = pG->U[k][j][i].d;
        *(pSnd++) = pG->U[k][j][i].M1;
        *(pSnd++) = pG->U[k][j][i].M2;
        *(pSnd++) = pG->U[k][j][i].M3;
#ifndef BAROTROPIC
        *(pSnd++) = pG->U[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        *(pSnd++) = pG->U[k][j][i].B1c;
        *(pSnd++) = pG->U[k][j][i].B2c;
        *(pSnd++) = pG->U[k][j][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) *(pSnd++) = pG->U[k][j][i].s[n];
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=je-(nghost-1); j<=je; j++){
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        *(pSnd++) = pG->B1i[k][j][i];
      }
    }
  }

/* B2i at j=je-(nghost-1) maps to B2i at j=js-nghost and is not passed */
  for (k=ks; k<=ke; k++) {
    for (j=je-(nghost-2); j<=je; j++){
      for (i=is-nghost; i<=ie+nghost; i++) {
        *(pSnd++) = pG->B2i[k][j][i];
      }
    }
  }

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=je-(nghost-1); j<=je; j++){
      for (i=is-nghost; i<=ie+nghost; i++) {
        *(pSnd++) = pG->B3i[k][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ix3(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Inner x3 boundary */

static void pack_ix3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
#if (NSCALARS > 0)
  int n;
#endif
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

  for (k=ks; k<=ks+(nghost-1); k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        *(pSnd++) = pG->U[k][j][i].d;
        *(pSnd++) = pG->U[k][j][i].M1;
        *(pSnd++) = pG->U[k][j][i].M2;
        *(pSnd++) = pG->U[k][j][i].M3;
#ifndef BAROTROPIC
        *(pSnd++) = pG->U[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        *(pSnd++) = pG->U[k][j][i].B1c;
        *(pSnd++) = pG->U[k][j][i].B2c;
        *(pSnd++) = pG->U[k][j][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) *(pSnd++) = pG->U[k][j][i].s[n];
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ks+(nghost-1); k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        *(pSnd++) = pG->B1i[k][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=ks; k<=ks+(nghost-1); k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        *(pSnd++) = pG->B2i[k][j][i];
      }
    }
  }

/* B3i at k=ks maps to B3i at k=ke+1 and is not passed */
  for (k=ks+1; k<=ks+(nghost-1); k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        *(pSnd++) = pG->B3i[k][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ox3(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Outer x3 boundary */

static void pack_ox3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
#if (NSCALARS > 0)
  int n;
#endif
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

  for (k=ke-(nghost-1); k<=ke; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        *(pSnd++) = pG->U[k][j][i].d;
        *(pSnd++) = pG->U[k][j][i].M1;
        *(pSnd++) = pG->U[k][j][i].M2;
        *(pSnd++) = pG->U[k][j][i].M3;
#ifndef BAROTROPIC
        *(pSnd++) = pG->U[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        *(pSnd++) = pG->U[k][j][i].B1c;
        *(pSnd++) = pG->U[k][j][i].B2c;
        *(pSnd++) = pG->U[k][j][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) *(pSnd++) = pG->U[k][j][i].s[n];
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ke-(nghost-1); k<=ke; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        *(pSnd++) = pG->B1i[k][j][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=ke-(nghost-1); k<=ke; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        *(pSnd++) = pG->B2i[k][j][i];
      }
    }
  }

/* B3i at k=ke-(nghost-1) maps to B3i at k=ks-nghost and is not passed */
  for (k=ke-(nghost-2); k<=ke; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        *(pSnd++) = pG->B3i[k][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ix1(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Inner x1 boundary */

static void unpack_ix1(GridS *pG)
{
  int is = pG->is;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif
#if (NSCALARS > 0)
  int n;
#endif
  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is-nghost; i<=is-1; i++){
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
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is-(nghost-1); i<=is-1; i++){
        pG->B1i[k][j][i] = *(pRcv++);
      }
    }
  }

  if (pG->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=is-nghost; i<=is-1; i++){
        pG->B2i[k][j][i] = *(pRcv++);
      }
    }
  }

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=is-nghost; i<=is-1; i++){
        pG->B3i[k][j][i] = *(pRcv++);
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ox1(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Outer x1 boundary */

static void unpack_ox1(GridS *pG)
{
  int ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif
#if (NSCALARS > 0)
  int n;
#endif
  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=ie+1; i<=ie+nghost; i++) {
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
/* B1i is not set at i=ie+1 */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=ie+2; i<=ie+nghost; i++) {
        pG->B1i[k][j][i] = *(pRcv++);
      }
    }
  }

  if (pG->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=ie+1; i<=ie+nghost; i++) {
        pG->B2i[k][j][i] = *(pRcv++);
      }
    }
  }

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=ie+1; i<=ie+nghost; i++) {
        pG->B3i[k][j][i] = *(pRcv++);
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ix2(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Inner x2 boundary */

static void unpack_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif
#if (NSCALARS > 0)
  int n;
#endif
  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

  for (k=ks; k<=ke; k++) {
    for (j=js-nghost; j<=js-1; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
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
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js-nghost; j<=js-1; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pG->B1i[k][j][i] = *(pRcv++);
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js-(nghost-1); j<=js-1; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->B2i[k][j][i] = *(pRcv++);
      }
    }
  }

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js-nghost; j<=js-1; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->B3i[k][j][i] = *(pRcv++);
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ox2(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Outer x2 boundary */

static void unpack_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif
#if (NSCALARS > 0)
  int n;
#endif
  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

  for (k=ks; k<=ke; k++) {
    for (j=je+1; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
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
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=je+1; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pG->B1i[k][j][i] = *(pRcv++);
      }
    }
  }

/* B2i is not set at j=je+1 */
  for (k=ks; k<=ke; k++) {
    for (j=je+2; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->B2i[k][j][i] = *(pRcv++);
      }
    }
  }

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=je+1; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->B3i[k][j][i] = *(pRcv++);
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ix3(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Inner x3 boundary */

static void unpack_ix3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks;
  int i,j,k;
#if (NSCALARS > 0)
  int n;
#endif
  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

  for (k=ks-nghost; k<=ks-1; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
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
/* B1i is not set at i=is-nghost */
  for (k=ks-nghost; k<=ks-1; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pG->B1i[k][j][i] = *(pRcv++);
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=ks-nghost; k<=ks-1; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->B2i[k][j][i] = *(pRcv++);
      }
    }
  }

/* B3i is not set at k=ks-nghost */
  for (k=ks-(nghost-1); k<=ks-1; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->B3i[k][j][i] = *(pRcv++);
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ox3(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Outer x3 boundary */

static void unpack_ox3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ke = pG->ke;
  int i,j,k;
#if (NSCALARS > 0)
  int n;
#endif
  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

  for (k=ke+1; k<=ke+nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
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
/* B1i is not set at i=is-nghost */
  for (k=ke+1; k<=ke+nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pG->B1i[k][j][i] = *(pRcv++);
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=ke+1; k<=ke+nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->B2i[k][j][i] = *(pRcv++);
      }
    }
  }

/* B3i is not set at k=ke+1 */
  for (k=ke+2; k<=ke+nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->B3i[k][j][i] = *(pRcv++);
      }
    }
  }
#endif /* MHD */

  return;
}
#endif /* MPI_PARALLEL */
