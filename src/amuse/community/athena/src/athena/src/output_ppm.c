#include "copyright.h"
/*============================================================================*/
/*! \file output_ppm.c
 *  \brief Writes single variable as a PPM image with color table.
 *
 * PURPOSE: Writes single variable as a PPM image with color table.  With SMR,
 *   dumps are made for all levels and domains, unless nlevel and ndomain are
 *   specified in <output> block.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - output_ppm()
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - compute_rgb()							      */
/*============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

static Real **data=NULL;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   compute_rgb()  
 *============================================================================*/

static void compute_rgb(double data, double min, double max, int *pR, int *pG,
                    int *pB, OutputS *p);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void output_ppm(MeshS *pM, OutputS *pOut) 
 *  \brief Output PPM image */
void output_ppm(MeshS *pM, OutputS *pOut)
{
  GridS *pGrid;
  FILE *pfile;
  char *fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
  int nl,nd,nx1,nx2,i,j;
  Real dmin, dmax;
  int red,green,blue;

/* check output data is 2D (output must be a 2D slice for 3D runs) */
  if (pOut->ndim != 2) {
    ath_error("[output_ppm:] Data must be 2D\n");
  }

/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){
    
/* write files if domain and level match input, or are not specified (-1) */
      if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
          (pOut->ndomain == -1 || pOut->ndomain == nd)){
        pGrid = pM->Domain[nl][nd].Grid;

/* Extract 2D data from 3D data,  Can either be slice or average along axis,
 * depending on range of ix1,ix2,ix3 in <ouput> block.  If OutData2 returns
 * a NULL pointer, then slice is outside of range of data in pGrid, so skip */

        data = OutData2(pGrid,pOut,&nx1,&nx2);
        if (data != NULL){

/* construct output filename.  pOut->id will either be name of variable,
 * if 'id=...' was included in <ouput> block, or 'outN' where N is number of
 * <output> block.  */
          if (nl>0) {
            plev = &levstr[0];
            sprintf(plev,"lev%d",nl);
          }
          if (nd>0) {
            pdom = &domstr[0];
            sprintf(pdom,"dom%d",nd);
          }

          if((fname = ath_fname(plev,pM->outfilename,plev,pdom,num_digit,
              pOut->num,pOut->id,"ppm")) == NULL){
            ath_error("[output_ppm]: Error constructing filename\n");
          }

/* open output file */
          if((pfile = fopen(fname,"w")) == NULL){
            ath_error("[output_ppm]: Unable to open ppm file %s\n",fname);
          }

/* Store the global min / max, for output at end of run */
          minmax2(data,nx2,nx1,&dmin,&dmax);
          pOut->gmin = MIN(dmin,pOut->gmin);
          pOut->gmax = MAX(dmax,pOut->gmax);

          fprintf(pfile,"P6\n");
          fprintf(pfile,"# dmin= %.7e, dmax= %.7e, gmin= %.7e, gmax= %.7e\n",
  	  dmin,dmax,pOut->gmin,pOut->gmax);
          fprintf(pfile,"%d %d\n255\n",nx1,nx2);

/* Override auto-scale? */
          if (pOut->sdmin != 0) dmin = pOut->dmin;
          if (pOut->sdmax != 0) dmax = pOut->dmax;

          for (j=nx2-1; j>=0; j--) {
            for (i=0; i<nx1; i++) {
              compute_rgb(data[j][i],dmin,dmax,&red,&green,&blue,pOut);
              fputc(red,pfile);
              fputc(green,pfile);
              fputc(blue,pfile);
            }
          }

/* Close the file, free memory */
          fclose(pfile);
          free_2d_array(data);
          free(fname);
          data = NULL;
        }
      }}
    }
  }
}

/*----------------------------------------------------------------------------*/
/* compute_rgb: converts data into RGB values using palette in Output struct  */

/*! \fn static void compute_rgb(double data, double min, double max, int *pR, 
 *			        int *pG, int *pB, OutputS *p);
 *  \brief Converts data into RGB values using palette in Output struct  */
static void compute_rgb(double data, double min, double max,
  int *R, int *G, int *B, OutputS *pOut)
{
  int i;
  float x, *rgb = pOut->rgb, *der = pOut->der;

  if (rgb == 0) {
    *R = *G = *B = (data > max ? 255 : 0);    
    return;
  }

  if (min==max) {
    *R = *G = *B = (data > max ? 255 : 0);    
    return;
  }
#if 1
  x = (data-min)*255.0/(max-min);
  if (x<=0.0 || x>=255.0) {         /* out of bounds */
    i=  (x <= 0.0 ? 0 : 255);
    *R = (int) (rgb[i*3+0] * 255.0);
    *G = (int) (rgb[i*3+1] * 255.0);
    *B = (int) (rgb[i*3+2] * 255.0);
    return;
  }
  i = (int) x;
  *R = (int)  ((rgb[i*3+0] + (x-i)*der[i*3+0])*255.0);
  *G = (int)  ((rgb[i*3+1] + (x-i)*der[i*3+1])*255.0);
  *B = (int)  ((rgb[i*3+2] + (x-i)*der[i*3+2])*255.0);
#else
  i = (int) ((data-min)*255.0/(max-min));
  if (i<0) i=0;
  if (i>255) i=255;
  *R = (int) (rgb[i*3+0] * 255.0);
  *G = (int) (rgb[i*3+1] * 255.0);
  *B = (int) (rgb[i*3+2] * 255.0);
#endif
}
