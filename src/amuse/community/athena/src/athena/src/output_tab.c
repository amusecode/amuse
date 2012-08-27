#include "copyright.h"
/*============================================================================*/
/*! \file output_tab.c
 *  \brief Functions for writing output in tabular format.
 *
 * PURPOSE: Functions for writing output in tabular format.  With SMR,
 *   dumps are made for all levels and domains, unless nlevel and ndomain are
 *   specified in <output> blocks.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - output_tab() - opens file and calls appropriate 1D/2D/3D output function
 *     Uses OutData1,2,3() to extract appropriate section to be output.
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - output_tab_1d() - write tab file for 1D slice of data
 * - output_tab_2d() - write tab file for 2D plane of data
 * - output_tab_3d() - write tab file for 3D section of data		      */
/*============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   output_tab_1d() - write tab file for 1D slice of data
 *   output_tab_2d() - write tab file for 2D plane of data
 *   output_tab_3d() - write tab file for 3D section of data
 *============================================================================*/

void output_tab_1d(MeshS *pM, OutputS *pOut, int nl, int nd);
void output_tab_2d(MeshS *pM, OutputS *pOut, int nl, int nd);
void output_tab_3d(MeshS *pM, OutputS *pOut, int nl, int nd);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void output_tab(MeshS *pM, OutputS *pOut)
 *  \brief Open file, call 1D/2D/3D writer; called by data_ouput  */

void output_tab(MeshS *pM, OutputS *pOut)
{
  int nl,nd;

/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){

/* write files if domain and level match input, or are not specified (-1) */
      if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
          (pOut->ndomain == -1 || pOut->ndomain == nd)){

        if (pOut->ndim == 3) {
          output_tab_3d(pM,pOut,nl,nd);
        } else if (pOut->ndim == 2) {
          output_tab_2d(pM,pOut,nl,nd);
        } else if (pOut->ndim == 1) {
          output_tab_1d(pM,pOut,nl,nd);
        }
      }}
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void output_tab_1d(MeshS *pM, OutputS *pOut, int nl, int nd) 
 *  \brief Writes 1D data.  Note x-coordinate is just i-index.  */

void output_tab_1d(MeshS *pM, OutputS *pOut, int nl, int nd)
{
  GridS *pGrid=pM->Domain[nl][nd].Grid;
  int i,nx1;
  FILE *pFile;
  char fmt[80],*fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
  Real *data=NULL;
  Real dmin, dmax, xworld;

/* Add a white space to the format, setup format for integer zone columns */
  if(pOut->dat_fmt == NULL){
     sprintf(fmt," %%12.8e"); /* Use a default format */
  }
  else{
    sprintf(fmt," %s",pOut->dat_fmt);
  }

/* compute 1D array of data */
  data = OutData1(pGrid,pOut,&nx1);
  if (data == NULL) return;  /* slice not in range of Grid */

  minmax1(data,nx1,&dmin,&dmax);

/* construct output filename */
  if (nl>0) {
    plev = &levstr[0];
    sprintf(plev,"lev%d",nl);
  }
  if (nd>0) {
    pdom = &domstr[0];
    sprintf(pdom,"dom%d",nd);
  }

  if((fname = ath_fname(plev,pM->outfilename,plev,pdom,num_digit,pOut->num,
      pOut->id,"tab")) == NULL){
    ath_error("[output_tab]: Error constructing filename\n");
  }

/* open filename */
  pFile = fopen(fname,"w");
  if (pFile == NULL) {
    ath_error("[output_tab]: Unable to open tab file %s\n",fname);
  }

/* write data */
  for (i=0; i<nx1; i++) {
    xworld = (float)(i);  /* just i index for now */
    fprintf(pFile,fmt,xworld);
    fprintf(pFile,fmt,data[i]);
    fprintf(pFile,"\n");
  }
  
/* Compute and store global min/max, for output at end of run */
  pOut->gmin = MIN(dmin,pOut->gmin);
  pOut->gmax = MAX(dmax,pOut->gmax);

  fclose(pFile);
  free_1d_array(data); /* Free the memory we malloc'd */
}

/*----------------------------------------------------------------------------*/
/*! \fn void output_tab_2d(MeshS *pM, OutputS *pOut, int nl, int nd)
 *  \brief Writes 2D data.  Note x/y-coordinate is just i/j-index.  */

void output_tab_2d(MeshS *pM, OutputS *pOut, int nl, int nd)
{
  GridS *pGrid=pM->Domain[nl][nd].Grid;
  int i,j,nx1,nx2;
  FILE *pFile;
  char fmt[80],*fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
  Real **data=NULL;
  Real dmin, dmax, xworld, yworld;

/* Add a white space to the format, setup format for integer zone columns */
  if(pOut->dat_fmt == NULL){
     sprintf(fmt," %%12.8e"); /* Use a default format */
  }
  else{
    sprintf(fmt," %s",pOut->dat_fmt);
  }

/* compute 2D array of data */
  data = OutData2(pGrid,pOut,&nx1,&nx2);
  if (data == NULL) return;  /* slice not in range of Grid */

  minmax2(data,nx2,nx1,&dmin,&dmax);

/* construct output filename */
  if (nl>0) {
    plev = &levstr[0];
    sprintf(plev,"lev%d",nl);
  }
  if (nd>0) {
    pdom = &domstr[0];
    sprintf(pdom,"dom%d",nd);
  }

  if((fname = ath_fname(plev,pM->outfilename,plev,pdom,num_digit,pOut->num,
      pOut->id,"tab")) == NULL){
    ath_error("[output_tab]: Error constructing filename\n");
  }

/* open filename */
  pFile = fopen(fname,"w");
  if (pFile == NULL) {
    ath_error("[output_tab]: Unable to open tab file %s\n",fname);
  }

/* write data */
  for (j=0; j<nx2; j++) {
    for (i=0; i<nx1; i++) {
      xworld = (float)(i); /* just i index for now */
      yworld = (float)(j); /* just j index for now */
      fprintf(pFile,fmt,xworld);
      fprintf(pFile,fmt,yworld);
      fprintf(pFile,fmt,data[j][i]);
      fprintf(pFile,"\n");
    }
  }
  
/* Compute and store global min/max, for output at end of run */
  if (pOut->num == 0) {
    pOut->gmin = dmin;
    pOut->gmax = dmax;
  } else {
    pOut->gmin = MIN(dmin,pOut->gmin);
    pOut->gmax = MAX(dmax,pOut->gmax);
  }

  fclose(pFile);
  free_2d_array(data); /* Free the memory we malloc'd */
}

/*----------------------------------------------------------------------------*/
/*! \fn void output_tab_3d(MeshS *pM, OutputS *pOut, int nl, int nd)
 *  \brief Writes 3D data.  Note x/y/z-coordinate is just i/j/k-index  */

void output_tab_3d(MeshS *pM, OutputS *pOut, int nl, int nd)
{
  GridS *pGrid=pM->Domain[nl][nd].Grid;
  int i,j,k,nx1,nx2,nx3;
  FILE *pFile;
  char fmt[80],*fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
  Real ***data, dmin, dmax, xworld, yworld, zworld;

/* Add a white space to the format, setup format for integer zone columns */
  if(pOut->dat_fmt == NULL){
     sprintf(fmt," %%12.8e"); /* Use a default format */
  }
  else{
    sprintf(fmt," %s",pOut->dat_fmt);
  }

/* compute 3D array of data */
  data = OutData3(pGrid,pOut,&nx1,&nx2,&nx3);
  minmax3(data,nx3,nx2,nx1,&dmin,&dmax);

/* construct output filename */
  if (nl>0) {
    plev = &levstr[0];
    sprintf(plev,"lev%d",nl);
  }
  if (nd>0) {
    pdom = &domstr[0];
    sprintf(pdom,"dom%d",nd);
  }

  if((fname = ath_fname(plev,pM->outfilename,plev,pdom,num_digit,pOut->num,
      pOut->id,"tab")) == NULL){
    ath_error("[output_tab]: Error constructing filename\n");
  }

/* open filename */
  pFile = fopen(fname,"w");
  if (pFile == NULL) {
    ath_error("[output_tab]: Unable to open tab file %s\n",fname);
  }

/* write data */
  for (k=0; k<nx3; k++) {
    for (j=0; j<nx2; j++) {
      for (i=0; i<nx1; i++) {
        xworld = (float)(i);  /* just i-index for now */
        yworld = (float)(j);  /* just j-index for now */
        zworld = (float)(k);  /* just k-index for now */
        fprintf(pFile,fmt,xworld);
        fprintf(pFile,fmt,yworld);
        fprintf(pFile,fmt,zworld);
        fprintf(pFile,fmt,data[k][j][i]);
        fprintf(pFile,"\n");
      }
    }
  }
  
/* Compute and store global min/max, for output at end of run */
  if (pOut->num == 0) {
    pOut->gmin = dmin;
    pOut->gmax = dmax;
  } else {
    pOut->gmin = MIN(dmin,pOut->gmin);
    pOut->gmax = MAX(dmax,pOut->gmax);
  }

  fclose(pFile);
  free_3d_array(data); /* Free the memory we malloc'd */
}
