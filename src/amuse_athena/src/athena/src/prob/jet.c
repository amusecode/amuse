#include "copyright.h"
/*============================================================================*/
/*! \file jet.c
 *  \brief Sets up a jet introduced through L-x1 boundary (left edge) */
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
 * jet_iib() - sets BCs on L-x1 (left edge) of grid.
 *============================================================================*/

void jet_iib(GridS *pGrid);

/* Make radius of jet and jet variables static global so they can be accessed
   by BC functions */
static Real rjet,Bxjet;
static Prim1DS Wjet;
static Cons1DS Ujet;
static Real x1_mid,x2_mid,x3_mid;

/*----------------------------------------------------------------------------*/
/* problem */

void problem(DomainS *pDomain){
   GridS *pGrid=(pDomain->Grid);
   int i, is = pGrid->is, ie = pGrid->ie;
   int j, js = pGrid->js, je = pGrid->je;
   int k, ks = pGrid->ks, ke = pGrid->ke;
   int il,iu,jl,ju,kl,ku;
   Prim1DS W;
   Cons1DS U;
   Real x1_min, x1_max;
   Real x2_min, x2_max;
   Real x3_min, x3_max;
   Real Bx=0.0;
   Bxjet = 0.0;

/* read parameters from input file */

   W.d  = par_getd("problem", "d");
   W.P  = par_getd("problem", "p");
   W.Vx = par_getd("problem", "vx");
   W.Vy = par_getd("problem", "vy");
   W.Vz = par_getd("problem", "vz");
#ifdef MHD
   Bx   = par_getd("problem", "bx");
   W.By = par_getd("problem", "by");
   W.Bz = par_getd("problem", "bz");
#endif

   Wjet.d  = par_getd("problem", "djet");
   Wjet.P  = par_getd("problem", "pjet");
   Wjet.Vx = par_getd("problem", "vxjet");
   Wjet.Vy = par_getd("problem", "vyjet");
   Wjet.Vz = par_getd("problem", "vzjet");
#ifdef MHD
   Bxjet   = par_getd("problem", "bxjet");
   Wjet.By = par_getd("problem", "byjet");
   Wjet.Bz = par_getd("problem", "bzjet");
#endif

   rjet = par_getd("problem", "rjet");
   
   U = Prim1D_to_Cons1D(&W, &Bx);
   Ujet = Prim1D_to_Cons1D(&Wjet, &Bxjet);

   x1_min = pDomain->RootMinX[0];
   x1_max = pDomain->RootMaxX[0];
   x2_min = pDomain->RootMinX[1];
   x2_max = pDomain->RootMaxX[1];
   x3_min = pDomain->RootMinX[2];
   x3_max = pDomain->RootMaxX[2];

   x1_mid = 0.5 * (x1_max + x1_min);
   x2_mid = 0.5 * (x2_max + x2_min);
   x3_mid = 0.5 * (x3_max + x3_min);

/* initialize index bounds assuming problem in xy plane */
   iu = ie + nghost;
   il = is - nghost;
   ju = je + nghost;
   jl = js - nghost;
   if(pGrid->Nx[2] > 1){
      ku = pGrid->ke + nghost;
      kl = pGrid->ks - nghost;
   }
   else{
      ku = pGrid->ke;
      kl = pGrid->ks;
   }

/* initialize conserved variables */
   
   for(k=kl; k<=ku; k++){
      for(j=jl; j<=ju; j++){
         for(i=il; i<=iu; i++){
            pGrid->U[k][j][i].d  = U.d;
            pGrid->U[k][j][i].M1 = U.Mx;
            pGrid->U[k][j][i].M2 = U.My;
            pGrid->U[k][j][i].M3 = U.Mz;
            pGrid->U[k][j][i].E  = U.E;
#ifdef MHD
            pGrid->U[k][j][i].B1c = Bx;
            pGrid->U[k][j][i].B2c = U.By;
            pGrid->U[k][j][i].B3c = U.Bz;
            pGrid->B1i[k][j][i] = Bx;
            pGrid->B2i[k][j][i] = U.By;
            pGrid->B3i[k][j][i] = U.Bz;
#endif
         }
      }
   }

/* Set boundary value function pointers */

   if (pDomain->Disp[0] == 0) bvals_mhd_fun(pDomain,left_x1,jet_iib);

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

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  int nl,nd;

  Wjet.d  = par_getd("problem", "djet");
  Wjet.P  = par_getd("problem", "pjet");
  Wjet.Vx = par_getd("problem", "vxjet");
  Wjet.Vy = par_getd("problem", "vyjet");
  Wjet.Vz = par_getd("problem", "vzjet");
#ifdef MHD
  Bxjet   = par_getd("problem", "bxjet");
  Wjet.By = par_getd("problem", "byjet");
  Wjet.Bz = par_getd("problem", "bzjet");
#endif

  rjet = par_getd("problem", "rjet");
  Ujet = Prim1D_to_Cons1D(&Wjet, &Bxjet);

  for (nl=0; nl<pM->NLevels; nl++){
  for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {
      if (pM->Domain[nl][nd].Disp[0] == 0) 
        bvals_mhd_fun(&(pM->Domain[nl][nd]),left_x1,jet_iib);
    }
  }}

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


/*===================== PRIVATE FUNCTIONS ====================================*/

/******************************************************************************/
/*! \fn void jet_iib(GridS *pGrid) 
 *  \brief Sets ghost zones to either outflow BC or
 * if cell is within jet radius, to jet values */
/******************************************************************************/

void jet_iib(GridS *pGrid){
  int i, is = pGrid->is;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real rad,x1,x2,x3;

  for(k=ks; k<=ke; k++){
    for(j=js; j<=je; j++){
      for(i=1; i<=nghost; i++){
        cc_pos(pGrid,(is-i),j,k,&x1,&x2,&x3);
        rad = sqrt(SQR(x2 - x2_mid) + SQR(x3 - x3_mid));
            
        if(rad <= rjet){
          pGrid->U[k][j][is-i].d  = Ujet.d;
          pGrid->U[k][j][is-i].M1 = Ujet.Mx;
          pGrid->U[k][j][is-i].M2 = Ujet.My;
          pGrid->U[k][j][is-i].M3 = Ujet.Mz;
          pGrid->U[k][j][is-i].E  = Ujet.E;
#ifdef MHD
          pGrid->U[k][j][is-i].B1c = Bxjet;
          pGrid->U[k][j][is-i].B2c = Ujet.By;
          pGrid->U[k][j][is-i].B3c = Ujet.Bz;
          pGrid->B1i[k][j][is-i] = Bxjet;
          pGrid->B2i[k][j][is-i] = Ujet.By;
          pGrid->B3i[k][j][is-i] = Ujet.Bz;
#endif
        } else{
          pGrid->U[k][j][is-i] = pGrid->U[k][j][is+(i-1)];
          pGrid->U[k][j][is-i].M1 = -pGrid->U[k][j][is-i].M1;
#ifdef MHD
          pGrid->U[k][j][is-i].B2c = -pGrid->U[k][j][is-i].B2c;
          pGrid->U[k][j][is-i].B3c = -pGrid->U[k][j][is-i].B3c;
          pGrid->B1i[k][j][is-i] =  pGrid->B1i[k][j][is+i];
          pGrid->B2i[k][j][is-i] = -pGrid->B2i[k][j][is+(i-1)];
          pGrid->B3i[k][j][is-i] = -pGrid->B3i[k][j][is+(i-1)];
#endif
        }
      }
    }
  }
}
