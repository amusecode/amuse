/*============================================================================*/
/*! \file ath_fft.c
 *  \brief Simple wrappers for 2D and 3D FFT functions. 
 *
 * PURPOSE:  Simple wrappers for 2D and 3D FFT functions.  These exist to
 *   hide the differences in function calls needed for single processor vs
 *   MPI FFT calls.  If you're concerned about performance, or want
 *   additional functionality, use these functions as examples to either
 *   write your own wrappers or to use the FFTW and/or block decomposition
 *   libraries directly.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - ath_3d_fft_quick_plan()   - create a plan for global 3D grid
 * - ath_3d_fft_create_plan()  - create a more flexible plan for 3D FFT
 * - ath_3d_fft_malloc()       - allocate memory for 3D FFT data
 * - ath_3d_fft()              - perform a 3D FFT
 * - ath_3d_fft_free()         - free memory for 3D FFT data
 * - ath_3d_fft_destroy_plan() - free up memory
 * - ath_2d_fft_quick_plan()   - create a plan for global 2D grid (Nx3=1)
 * - ath_2d_fft_create_plan()  - create a more flexible plan for 2D FFT
 * - ath_2d_fft_malloc()       - allocate memory for 2D FFT data
 * - ath_2d_fft()              - perform a 2D FFT
 * - ath_2d_fft_free()         - free memory for 2D FFT data
 * - ath_2d_fft_destroy_plan() - free up memory				      */
/*============================================================================*/

#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "../globals.h"

#ifdef FFT_ENABLED

#ifdef FFT_BLOCK_DECOMP
/* Include Steve Plimpton's FFTW interface code */
#include "fft_3d.h"
#include "fft_2d.h"
#else /* FFT_BLOCK_DECOMP */
/* For a single processor, use FFTW directly */
#include "fftw3.h"
#endif /* FFT_BLOCK_DECOMP */

/**************************************************************************
 *
 *  Athena 3D FFT functions
 *
 **************************************************************************/

/*! \fn struct ath_3d_fft_plan *ath_3d_fft_quick_plan(DomainS *pD,
 *	ath_fft_data *data, ath_fft_direction dir)
 *  \brief Sets up an FFT plan for the entire 3D grid, using
 *      ath_3d_fft_create_plan()
 */

struct ath_3d_fft_plan *ath_3d_fft_quick_plan(DomainS *pD,
				ath_fft_data *data, ath_fft_direction dir)
{
  GridS *pGrid = (pD->Grid);
  /* Get size of global FFT grid */
  int gnx1 = pD->Nx[0];
  int gnx2 = pD->Nx[1];
  int gnx3 = pD->Nx[2];

  /* Get extents of local FFT grid in global coordinates */
  int gis = pGrid->Disp[0];
  int gie = pGrid->Disp[0] + pGrid->Nx[0] - 1;
  int gjs = pGrid->Disp[1];
  int gje = pGrid->Disp[1] + pGrid->Nx[1] - 1;
  int gks = pGrid->Disp[2];
  int gke = pGrid->Disp[2] + pGrid->Nx[2] - 1;

  /* Create the plan using a more generic function.
   * If the data hasn't already been allocated, it will now */
  return ath_3d_fft_create_plan(gnx3, gnx2, gnx1, gks, gke, gjs, gje,
				   gis, gie, data, 0, dir);
}

/*! \fn struct ath_3d_fft_plan *ath_3d_fft_create_plan(int gnx3, int gnx2,
 *				int gnx1, int gks, int gke, int gjs, int gje,
 *				int gis, int gie, ath_fft_data *data, int al,
 *				ath_fft_direction dir)
 *  \brief Sets up a 3D FFT plan
 *
 *  - gnx3, gnx2, gnx1 are the dimensions of the GLOBAL data
 *  - gks, gke, gjs, gje, gis, gie are the starting and ending indices of
 *      the LOCAL data in GLOBAL coordinates
 *  - data is any array of type ath_fft_data big enough to hold entire
 *      transform, for use in planning (contents will be trashed)
 *  - al != 0 means allocate data if it doesn't exist (otherwise temporary)
 *  - dir is either ATH_FFT_FOWARD or ATH_FFT_BACKWARD
 *  FFTs will be done in place (overwrite data)
 */

struct ath_3d_fft_plan *ath_3d_fft_create_plan(int gnx3, int gnx2,
				int gnx1, int gks, int gke, int gjs, int gje,
				int gis, int gie, ath_fft_data *data, int al,
				ath_fft_direction dir)
{
  int nbuf, tmp;
  struct ath_3d_fft_plan *ath_plan;

  if ((dir != ATH_FFT_FORWARD) && (dir != ATH_FFT_BACKWARD)) {
    ath_error("Invalid Athena FFT direction.\n");
  }

  /* Allocate memory for the plan */
  ath_plan = (struct ath_3d_fft_plan *)malloc(sizeof(struct ath_3d_fft_plan));
  if (ath_plan == NULL) {
    ath_error("Couldn't malloc for FFT plan.");
  }
  /* Set forward/backward FFT */
  ath_plan->dir = dir;
  /* Set element count (for easy malloc and memset) */
  ath_plan->cnt = (gke-gks+1)*(gje-gjs+1)*(gie-gis+1);
  ath_plan->gcnt = gnx3*gnx2*gnx1;

  tmp = (al==0 ? 1 : 0);
  if (data != NULL) tmp = 0;

  /* If data == NULL, then allocate something (temporarily if tmp=1) */
  if (data == NULL)
    data = (ath_fft_data *)ath_3d_fft_malloc(ath_plan);
  if (data == NULL)
    ath_error("Couln't malloc for FFT plan data.");

  /* Create the plan */
#ifdef FFT_BLOCK_DECOMP
  /* Block decomp library plans don't care if forward or backward */
  ath_plan->plan = fft_3d_create_plan(MPI_COMM_WORLD, gnx3, gnx2, gnx1, 
					gks, gke, gjs, gje, gis, gie, 
			    		gks, gke, gjs, gje, gis, gie, 
                            		0, 0, &nbuf);
#else /* FFT_BLOCK_DECOMP */
  if (dir == ATH_FFT_FORWARD) {
    ath_plan->plan = fftw_plan_dft_3d(gnx1, gnx2, gnx3, data, data,
					FFTW_FORWARD, FFTW_MEASURE);
  } else {
    ath_plan->plan = fftw_plan_dft_3d(gnx1, gnx2, gnx3, data, data,
					FFTW_BACKWARD, FFTW_MEASURE);
  }
#endif /* FFT_BLOCK_DECOMP */

  if (tmp) ath_3d_fft_free(data);

  return ath_plan;
}

/*! \fn ath_fft_data *ath_3d_fft_malloc(struct ath_3d_fft_plan *ath_plan)
 *  \brief Easy allocation of data array needed for particular 3D plan 
 */

ath_fft_data *ath_3d_fft_malloc(struct ath_3d_fft_plan *ath_plan)
{
  return (ath_fft_data *)fftw_malloc(sizeof(ath_fft_data) * ath_plan->cnt);
}

/*! \fn void ath_3d_fft(struct ath_3d_fft_plan *ath_plan, ath_fft_data *data)
 *  \brief Performs a 3D FFT in place
 */

void ath_3d_fft(struct ath_3d_fft_plan *ath_plan, ath_fft_data *data)
{
#ifdef FFT_BLOCK_DECOMP
  fft_3d(data, data, ath_plan->dir, ath_plan->plan);
#else /* FFT_BLOCK_DECOMP */
  /* Plan already includes forward/backward */
  fftw_execute_dft(ath_plan->plan, data, data);
#endif /* FFT_BLOCK_DECOMP */

  return;
}

/*! \fn void ath_3d_fft_free(ath_fft_data *data)
 *  \brief Frees memory used to hold data for 3D FFT
 */

void ath_3d_fft_free(ath_fft_data *data)
{
  if (data != NULL) fftw_free((void*)data);

  return;
}

/*! \fn void ath_3d_fft_destroy_plan(struct ath_3d_fft_plan *ath_plan)
 *  \brief Frees a 3D FFT plan
 */

void ath_3d_fft_destroy_plan(struct ath_3d_fft_plan *ath_plan)
{
  if (ath_plan != NULL) {
#ifdef FFT_BLOCK_DECOMP
    fft_3d_destroy_plan(ath_plan->plan);
#else /* FFT_BLOCK_DECOMP */
    fftw_destroy_plan(ath_plan->plan);
#endif /* FFT_BLOCK_DECOMP */
    free(ath_plan);
  }

  return;
}

/**************************************************************************
 *
 *  Athena 2D FFT functions
 *
 **************************************************************************/

/*! \fn struct ath_2d_fft_plan *ath_2d_fft_quick_plan(DomainS *pD,
 *				ath_fft_data *data, ath_fft_direction dir)
 *  \brief Sets up an FFT plan for the entire 2D grid, assuming NX3=1, using
 *      ath_2d_fft_create_plan()
 */

struct ath_2d_fft_plan *ath_2d_fft_quick_plan(DomainS *pD,
				ath_fft_data *data, ath_fft_direction dir)
{
  GridS *pGrid=(pD->Grid);
  if (pGrid->Nx[2] != 1) {
    ath_error("ath_2d_fft_quick_plan only works for Nx3=1.\n");
  }

  /* Get size of global FFT grid */
  int gnx1 = pD->Nx[0];
  int gnx2 = pD->Nx[1];

  /* Get extents of local FFT grid in global coordinates */
  int gis = pGrid->Disp[0];
  int gie = pGrid->Disp[0] + pGrid->Nx[0] - 1;
  int gjs = pGrid->Disp[1];
  int gje = pGrid->Disp[1] + pGrid->Nx[1] - 1;

  /* Create the plan using a more generic function
   * If the data hasn't already been allocated, it will now */
  return ath_2d_fft_create_plan(gnx2, gnx1, gjs, gje, gis, gie, data, 0, dir);
}

/*! \fn struct ath_2d_fft_plan *ath_2d_fft_create_plan(int gnx2, int gnx1,
 *				int gjs, int gje, int gis, int gie,
 *				ath_fft_data *data, int al,
 *				ath_fft_direction dir)
 *  \brief Sets up a 2D FFT plan
 *
 *  - gnx2, gnx1 are the dimensions of the GLOBAL data
 *  - gjs, gje, gis, gie are the starting and ending indices of the
 *      LOCAL data in GLOBAL coordinates
 *  - dir is either ATH_FFT_FOWARD or ATH_FFT_BACKWARD
 *  FFTs will be done in place (overwrite data)
 */

struct ath_2d_fft_plan *ath_2d_fft_create_plan(int gnx2, int gnx1,
				int gjs, int gje, int gis, int gie,
				ath_fft_data *data, int al,
				ath_fft_direction dir)
{
  int nbuf, tmp;
  struct ath_2d_fft_plan *ath_plan;

  if ((dir != ATH_FFT_FORWARD) && (dir != ATH_FFT_BACKWARD)) {
    ath_error("Invalid Athena FFT direction.\n");
  }

  /* Allocate memory for plan */
  ath_plan = (struct ath_2d_fft_plan *)malloc(sizeof(struct ath_2d_fft_plan));
  if (ath_plan == NULL) {
    ath_error("Couldn't malloc for FFT plan.");
  }
  /* Set forward/backward FFT */
  ath_plan->dir = dir;
  /* Set element count (for easy malloc and memset) */
  ath_plan->cnt = (gje-gjs+1)*(gie-gis+1);
  ath_plan->gcnt = gnx2*gnx1;

  tmp = (al==0 ? 1 : 0);
  if (data != NULL) tmp = 0;

  /* If data == NULL, then allocate something (temporarily if tmp=1) */
  if (data == NULL)
    data = (ath_fft_data *)ath_2d_fft_malloc(ath_plan);
  if (data == NULL)
    ath_error("Couln't malloc for FFT plan data.");

  /* Create the plan */
#ifdef FFT_BLOCK_DECOMP
  /* Block decomp plans don't care if forward/backward */
  ath_plan->plan = fft_2d_create_plan(MPI_COMM_WORLD, gnx2, gnx1, gjs, gje,
					gis, gie, gjs, gje, gis, gie, 
                            		0, 0, &nbuf);
#else /* FFT_BLOCK_DECOMP */
  if (dir == ATH_FFT_FORWARD) {
    ath_plan->plan = fftw_plan_dft_2d(gnx1, gnx2, data, data, FFTW_FORWARD,
					FFTW_MEASURE);
  } else {
    ath_plan->plan = fftw_plan_dft_2d(gnx1, gnx2, data, data, FFTW_BACKWARD,
					FFTW_MEASURE);
  }
#endif /* FFT_BLOCK_DECOMP */

  if (tmp) ath_2d_fft_free(data);

  return ath_plan;
}

/*! \fn ath_fft_data *ath_2d_fft_malloc(struct ath_2d_fft_plan *ath_plan)
 *  \brief Easy allocation of data array needed for particular 2D plan
 */

ath_fft_data *ath_2d_fft_malloc(struct ath_2d_fft_plan *ath_plan)
{
  return (ath_fft_data *)fftw_malloc(sizeof(ath_fft_data) * ath_plan->cnt);
}


/*! \fn void ath_2d_fft(struct ath_2d_fft_plan *ath_plan, fftw_complex *data)
 *  \brief Performs a 2D FFT in place
 */

void ath_2d_fft(struct ath_2d_fft_plan *ath_plan, fftw_complex *data)
{
#ifdef FFT_BLOCK_DECOMP
  fft_2d(data, data, ath_plan->dir, ath_plan->plan);
#else /* FFT_BLOCK_DECOMP */
  /* Plan already includes forward/backward */
  fftw_execute_dft(ath_plan->plan, data, data);
#endif /* FFT_BLOCK_DECOMP */

  return;
}

/*! \fn void ath_2d_fft_free(ath_fft_data *data)
 *  \brief Frees memory used to hold data for 2D FFT
 */

void ath_2d_fft_free(ath_fft_data *data)
{
  if (data != NULL) fftw_free((void*)data);

  return;
}

/*! \fn void ath_2d_fft_destroy_plan(struct ath_2d_fft_plan *ath_plan)
 *  \brief Frees a 2D FFT plan
 */

void ath_2d_fft_destroy_plan(struct ath_2d_fft_plan *ath_plan)
{
  if (ath_plan != NULL) {
#ifdef FFT_BLOCK_DECOMP
    fft_2d_destroy_plan(ath_plan->plan);
#else /* FFT_BLOCK_DECOMP */
    fftw_destroy_plan(ath_plan->plan);
#endif /* FFT_BLOCK_DECOMP */
    free(ath_plan);
  }

  return;
}

#endif /* FFT_ENABLED */
