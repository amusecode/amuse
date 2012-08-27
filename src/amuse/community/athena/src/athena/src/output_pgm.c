#include "copyright.h"
/*============================================================================*/
/*! \file output_pgm.c
 *  \brief Writes Portable Gray Map (PGM) outputs.
 *
 * PURPOSE: Writes Portable Gray Map (PGM) outputs.  These are extremely simple
 *   grayscale 2D images, see e.g. sourceforge for documentation. With SMR,
 *   dumps are made for all levels and domains, unless nlevel and ndomain are
 *   specified in <output> block.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - output_pgm() -  outputs 2D PGM images
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/*! \fn void output_pgm(MeshS *pM, OutputS *pOut)
 *  \brief  Output 2D PGM image   */

void output_pgm(MeshS *pM, OutputS *pOut)
{
  GridS *pGrid;
  FILE *pfile;
  char *fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
  int nl,nd,nx1,nx2,gray,i,j;
  Real **data, dmin, dmax, max_min, sfact;

/* check output data is 2D (output must be a 2D slice for 3D runs) */
  if (pOut->ndim != 2) {
    ath_error("[output_pgm]: Output must be a 2D slice\n");
    return;
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
 * NULL pointer, then slice is outside range of data in pGrid, so skip */

        data = OutData2(pGrid,pOut,&nx1,&nx2);
        if (data != NULL) {

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
            pOut->num,pOut->id,"pgm")) == NULL){
            ath_error("[output_pgm]: Error constructing filename\n");
          }

/* open output file */
          if((pfile = fopen(fname,"w")) == NULL){
            ath_error("[output_pgm]: Unable to open pgm file %s\n",fname);
            return;
          }

          fprintf(pfile,"P5\n%d %d\n255\n",nx1,nx2);

/* Store the global min / max, for output at end of run */
          minmax2(data,nx2,nx1,&dmin,&dmax);
          pOut->gmin = MIN(dmin,pOut->gmin);
          pOut->gmax = MAX(dmax,pOut->gmax);

/* Override auto-scale? */
          if (pOut->sdmin != 0) dmin = pOut->dmin;
          if (pOut->sdmax != 0) dmax = pOut->dmax;
  
          max_min = (dmax - dmin)*(1.0 + FLT_EPSILON);

/* map the data which satisfies [min <= data <= max] to the range 
 * [0.0 , 256.0] -- Not inclusive of 256 */

          if(max_min > 0.0) {
            sfact = 256.0/max_min;
            for (j=nx2-1; j>=0; j--) {
              for (i=0; i<nx1; i++) {
/* Map the data to an 8 bit int, i.e. 0 -> 255 */
                gray = (int)(sfact*(data[j][i] - dmin));
/* Out of bounds data is mapped to the min or max integer value */
                gray = gray >   0 ? gray :   0;
                gray = gray < 255 ? gray : 255;

                fputc(gray,pfile);
              }
            }

/* else, if max=min set image to constant */

          } else {
            gray = 0;
            for (j=0; j<nx2; j++) {
              for (i=0; i<nx1; i++) {
                fputc(gray,pfile);
              }
            }
          }

/* Close the file, free memory */

          fclose(pfile); 
          free_2d_array(data);
          free(fname);
        }
      }}
    }
  }
}
