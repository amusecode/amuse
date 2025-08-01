#include "copyright.h"
/*============================================================================*/
/*! \file dump_vtk.c 
 *  \brief Function to write a dump in VTK "legacy" format.
 *
 * PURPOSE: Function to write a dump in VTK "legacy" format.  With SMR,
 *   dumps are made for all levels and domains, unless nlevel and ndomain are
 *   specified in <output> block.  Works for BOTH conserved and primitives.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - dump_vtk() - writes VTK dump (all variables).			      */
/*============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"
#ifdef PARTICLES
#include "particles/particle.h"
#endif

/*----------------------------------------------------------------------------*/
/*! \fn void dump_vtk(MeshS *pM, OutputS *pOut)
 *  \brief Writes VTK dump (all variables).				      */

void dump_vtk(MeshS *pM, OutputS *pOut)
{
  GridS *pGrid;
  PrimS ***W;
  FILE *pfile;
  char *fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
/* Upper and Lower bounds on i,j,k for data dump */
  int i,j,k,il,iu,jl,ju,kl,ku,nl,nd;
  int big_end = ath_big_endian();
  int ndata0,ndata1,ndata2;
  float *data;   /* points to 3*ndata0 allocated floats */
  double x1, x2, x3;
#if (NSCALARS > 0)
  int n;
#endif

/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){

/* write files if domain and level match input, or are not specified (-1) */
      if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
          (pOut->ndomain == -1 || pOut->ndomain == nd)){
        pGrid = pM->Domain[nl][nd].Grid;

        il = pGrid->is, iu = pGrid->ie;
        jl = pGrid->js, ju = pGrid->je;
        kl = pGrid->ks, ku = pGrid->ke;

#ifdef WRITE_GHOST_CELLS
        iu = pGrid->ie + nghost;
        il = pGrid->is - nghost;

        if(pGrid->Nx[1] > 1) {
          ju = pGrid->je + nghost;
          jl = pGrid->js - nghost;
        }

        if(pGrid->Nx[2] > 1) {
          ku = pGrid->ke + nghost;
          kl = pGrid->ks - nghost;
        }
#endif /* WRITE_GHOST_CELLS */

        ndata0 = iu-il+1;
        ndata1 = ju-jl+1;
        ndata2 = ku-kl+1;

/* calculate primitive variables, if needed */

        if(strcmp(pOut->out,"prim") == 0) {
          if((W = (PrimS***)calloc_3d_array(ndata2,ndata1,ndata0,sizeof(PrimS)))
             == NULL) ath_error("[dump_vtk]: failed to allocate Prim array\n");

          for (k=kl; k<=ku; k++) {
          for (j=jl; j<=ju; j++) {
          for (i=il; i<=iu; i++) {
            W[k-kl][j-jl][i-il] = Cons_to_Prim(&(pGrid->U[k][j][i]));
          }}}
        }

/* construct filename, open file */
        if (nl>0) {
          plev = &levstr[0];
          sprintf(plev,"lev%d",nl);
        }
        if (nd>0) {
          pdom = &domstr[0];
          sprintf(pdom,"dom%d",nd);
        }
        if((fname = ath_fname(plev,pM->outfilename,plev,pdom,num_digit,
            pOut->num,NULL,"vtk")) == NULL){
          ath_error("[dump_vtk]: Error constructing filename\n");
        }

        if((pfile = fopen(fname,"w")) == NULL){
          ath_error("[dump_vtk]: Unable to open vtk dump file\n");
          return;
        }

/* Allocate memory for temporary array of floats */

        if((data = (float *)malloc(3*ndata0*sizeof(float))) == NULL){
          ath_error("[dump_vtk]: malloc failed for temporary array\n");
          return;
        }

/* There are five basic parts to the VTK "legacy" file format.  */
/*  1. Write file version and identifier */

        fprintf(pfile,"# vtk DataFile Version 2.0\n");

/*  2. Header */

        if (strcmp(pOut->out,"cons") == 0){
          fprintf(pfile,"CONSERVED vars at time= %e, level= %i, domain= %i\n",
            pGrid->time,nl,nd);
        } else if(strcmp(pOut->out,"prim") == 0) {
          fprintf(pfile,"PRIMITIVE vars at time= %e, level= %i, domain= %i\n",
            pGrid->time,nl,nd);
        }

/*  3. File format */

        fprintf(pfile,"BINARY\n");

/*  4. Dataset structure */

/* Set the Grid origin */

        x1 = pGrid->MinX[0];
        x2 = pGrid->MinX[1];
        x3 = pGrid->MinX[2];

        fprintf(pfile,"DATASET STRUCTURED_POINTS\n");
        if (pGrid->Nx[1] == 1) {
          fprintf(pfile,"DIMENSIONS %d %d %d\n",iu-il+2,1,1);
        } else {
          if (pGrid->Nx[2] == 1) {
            fprintf(pfile,"DIMENSIONS %d %d %d\n",iu-il+2,ju-jl+2,1);
          } else {
            fprintf(pfile,"DIMENSIONS %d %d %d\n",iu-il+2,ju-jl+2,ku-kl+2);
          }
        }
        fprintf(pfile,"ORIGIN %e %e %e \n",x1,x2,x3);
        fprintf(pfile,"SPACING %e %e %e \n",pGrid->dx1,pGrid->dx2,pGrid->dx3);

/*  5. Data  */

        fprintf(pfile,"CELL_DATA %d \n", (iu-il+1)*(ju-jl+1)*(ku-kl+1));

/* Write density */

        fprintf(pfile,"SCALARS density float\n");
        fprintf(pfile,"LOOKUP_TABLE default\n");
        for (k=kl; k<=ku; k++) {
          for (j=jl; j<=ju; j++) {
            for (i=il; i<=iu; i++) {
              if (strcmp(pOut->out,"cons") == 0){
                data[i-il] = (float)pGrid->U[k][j][i].d;
              } else if(strcmp(pOut->out,"prim") == 0) {
                data[i-il] = (float)W[k-kl][j-jl][i-il].d;
              }
            }
            if(!big_end) ath_bswap(data,sizeof(float),iu-il+1);
            fwrite(data,sizeof(float),(size_t)ndata0,pfile);
          }
        }

/* Write momentum or velocity */

        if (strcmp(pOut->out,"cons") == 0){
          fprintf(pfile,"\nVECTORS momentum float\n");
        } else if(strcmp(pOut->out,"prim") == 0) {
          fprintf(pfile,"\nVECTORS velocity float\n");
        }
        for (k=kl; k<=ku; k++) {
          for (j=jl; j<=ju; j++) {
            for (i=il; i<=iu; i++) {
              if (strcmp(pOut->out,"cons") == 0){
                data[3*(i-il)  ] = (float)pGrid->U[k][j][i].M1;
                data[3*(i-il)+1] = (float)pGrid->U[k][j][i].M2;
                data[3*(i-il)+2] = (float)pGrid->U[k][j][i].M3;
              } else if(strcmp(pOut->out,"prim") == 0) {
                data[3*(i-il)  ] = (float)W[k-kl][j-jl][i-il].V1;
                data[3*(i-il)+1] = (float)W[k-kl][j-jl][i-il].V2;
                data[3*(i-il)+2] = (float)W[k-kl][j-jl][i-il].V3;
              }
            }
            if(!big_end) ath_bswap(data,sizeof(float),3*(iu-il+1));
            fwrite(data,sizeof(float),(size_t)(3*ndata0),pfile);
          }
        }

/* Write total energy or pressure */

#ifndef BAROTROPIC
        if (strcmp(pOut->out,"cons") == 0){
          fprintf(pfile,"\nSCALARS total_energy float\n");
        } else if(strcmp(pOut->out,"prim") == 0) {
          fprintf(pfile,"\nSCALARS pressure float\n");
        }
        fprintf(pfile,"LOOKUP_TABLE default\n");
        for (k=kl; k<=ku; k++) {
          for (j=jl; j<=ju; j++) {
            for (i=il; i<=iu; i++) {
              if (strcmp(pOut->out,"cons") == 0){
                data[i-il] = (float)pGrid->U[k][j][i].E;
              } else if(strcmp(pOut->out,"prim") == 0) {
                data[i-il] = (float)W[k-kl][j-jl][i-il].P;
              }
            }
            if(!big_end) ath_bswap(data,sizeof(float),iu-il+1);
            fwrite(data,sizeof(float),(size_t)ndata0,pfile);
          }
        }
#endif

/* Write cell centered B */

#ifdef MHD
        fprintf(pfile,"\nVECTORS cell_centered_B float\n");
        for (k=kl; k<=ku; k++) {
          for (j=jl; j<=ju; j++) {
            for (i=il; i<=iu; i++) {
              data[3*(i-il)] = (float)pGrid->U[k][j][i].B1c;
              data[3*(i-il)+1] = (float)pGrid->U[k][j][i].B2c;
              data[3*(i-il)+2] = (float)pGrid->U[k][j][i].B3c;
            }
            if(!big_end) ath_bswap(data,sizeof(float),3*(iu-il+1));
            fwrite(data,sizeof(float),(size_t)(3*ndata0),pfile);
          }
        }
#endif

/* Write gravitational potential */

#ifdef SELF_GRAVITY
        fprintf(pfile,"\nSCALARS gravitational_potential float\n");
        fprintf(pfile,"LOOKUP_TABLE default\n");
        for (k=kl; k<=ku; k++) {
          for (j=jl; j<=ju; j++) {
            for (i=il; i<=iu; i++) {
              data[i-il] = (float)pGrid->Phi[k][j][i];
            }
            if(!big_end) ath_bswap(data,sizeof(float),iu-il+1);
            fwrite(data,sizeof(float),(size_t)ndata0,pfile);
          }
        }
#endif

/* Write binned particle grid */

#ifdef PARTICLES
        if (pOut->out_pargrid) {
          fprintf(pfile,"\nSCALARS particle_density float\n");
          fprintf(pfile,"LOOKUP_TABLE default\n");
          for (k=kl; k<=ku; k++) {
            for (j=jl; j<=ju; j++) {
              for (i=il; i<=iu; i++) {
                data[i-il] = pGrid->Coup[k][j][i].grid_d;
              }
              if(!big_end) ath_bswap(data,sizeof(float),iu-il+1);
              fwrite(data,sizeof(float),(size_t)ndata0,pfile);
            }
          }
          fprintf(pfile,"\nVECTORS particle_momentum float\n");
          for (k=kl; k<=ku; k++) {
            for (j=jl; j<=ju; j++) {
              for (i=il; i<=iu; i++) {
                data[3*(i-il)] = pGrid->Coup[k][j][i].grid_v1;
                data[3*(i-il)+1] = pGrid->Coup[k][j][i].grid_v2;
                data[3*(i-il)+2] = pGrid->Coup[k][j][i].grid_v3;
              }
              if(!big_end) ath_bswap(data,sizeof(float),3*(iu-il+1));
              fwrite(data,sizeof(float),(size_t)(3*ndata0),pfile);
            }
          }
        }
#endif

/* Write passive scalars */

#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++){
          if (strcmp(pOut->out,"cons") == 0){
            fprintf(pfile,"\nSCALARS scalar[%d] float\n",n);
          } else if(strcmp(pOut->out,"prim") == 0) {
            fprintf(pfile,"\nSCALARS specific_scalar[%d] float\n",n);
          }
          fprintf(pfile,"LOOKUP_TABLE default\n");
          for (k=kl; k<=ku; k++) {
            for (j=jl; j<=ju; j++) {
              for (i=il; i<=iu; i++) {
                if (strcmp(pOut->out,"cons") == 0){
                  data[i-il] = (float)pGrid->U[k][j][i].s[n];
                } else if(strcmp(pOut->out,"prim") == 0) {
                  data[i-il] = (float)W[k-kl][j-jl][i-il].r[n];
                }
              }
              if(!big_end) ath_bswap(data,sizeof(float),iu-il+1);
              fwrite(data,sizeof(float),(size_t)ndata0,pfile);
            }
          }
        }
#endif

/* close file and free memory */

        fclose(pfile);
        free(data);
        if(strcmp(pOut->out,"prim") == 0) free_3d_array(W);
      }}
    }
  }
  return;
}
