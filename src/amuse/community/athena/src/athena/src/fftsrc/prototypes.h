#ifndef ATHENA_FFT_H
#define ATHENA_FFT_H
/**************************************************************************
 *
 *  Athena 2D and 3D FFT functions
 *
 *  These provide extremely limited functionality, and exist only to hide
 *  the differences in function calls needed for single processor vs MPI.
 *
 *  If you're concerned about performance, or want additional functionality,
 *  use these functions as examples to either write your own wrappers or to
 *  use the FFTW and/or block decomposition libraries directly.
 *
 *  The indexing convention of the FFT data is DIFFERENT than Athena's.
 *  See F3DI() and F2DI() macros below for how to access the data.
 *
 *  Written by Nicole Lemaster on February 25, 2007
 *
 *  Last updated June 14, 2007
 *
 **************************************************************************/

#include "../athena.h"
#include "../defs.h"

#ifdef FFT_ENABLED

#ifndef MPI_PARALLEL
/* For a single processor, use FFTW 3.1.2 directly */
#define FFT_SLAB_DECOMP
#else /* MPI_PARALLEL */
/* Do multi-processor FFTs using block decomposition */
#define FFT_BLOCK_DECOMP
#endif /* MPI_PARALLEL */

#ifdef FFT_BLOCK_DECOMP
/* Include Steve Plimpton's FFTW interface code */
#include "fft_3d.h"
#else /* FFT_BLOCK_DECOMP */
/* For a single processor, use FFTW 3.1.2 directly */
#include "fftw3.h"
#endif /* FFT_BLOCK_DECOMP */

#define ath_fft_data fftw_complex

/* Indexing convention of FFT data
 * FFT Nfast=k, Nmid=j, Nslow=i (opposite to Athena) */
#define F3DI(i, j, k, nx1, nx2, nx3) ((k) + (nx3)*((j) + (nx2)*(i)))
#define F2DI(i, j, nx1, nx2) ((j) + (nx2)*(i))

/* Any component of wavenumber k
 * e.g. KCOMP(i-is, is+idisp, nx1)
 * where nx1 is the size of the global grid */
#define KCOMP(a, gas, gnxa) ((double)(((a)+(gas))-(int)(2*((a)+(gas))/(gnxa))*(gnxa)))

typedef enum {
  ATH_FFT_FORWARD=-1, ATH_FFT_BACKWARD=1
} ath_fft_direction;

struct ath_3d_fft_plan {
#ifdef FFT_BLOCK_DECOMP
  struct fft_plan_3d *plan;
#else /* FFT_BLOCK_DECOMP */
  fftw_plan plan;
#endif /* FFT_BLOCK_DECOMP */
  ath_fft_direction dir;
  long int cnt;
  long int gcnt;
};

struct ath_2d_fft_plan {
#ifdef FFT_BLOCK_DECOMP
  struct fft_plan_2d *plan;
#else /* FFT_BLOCK_DECOMP */
  fftw_plan plan;
#endif /* FFT_BLOCK_DECOMP */
  ath_fft_direction dir;
  long int cnt;
  long int gcnt;
};

/**************************************************************************
 *
 *  Athena 3D FFT functions
 *
 **************************************************************************/

struct ath_3d_fft_plan *ath_3d_fft_quick_plan(DomainS *pD,
				ath_fft_data *data, ath_fft_direction dir);
struct ath_3d_fft_plan *ath_3d_fft_create_plan(int gnx3, int gnx2,
				int gnx1, int gks, int gke, int gjs, int gje,
				int gis, int gie, ath_fft_data *data, int al,
				ath_fft_direction dir);
ath_fft_data *ath_3d_fft_malloc(struct ath_3d_fft_plan *ath_plan);
void ath_3d_fft(struct ath_3d_fft_plan *ath_plan, ath_fft_data *data);
void ath_3d_fft_free(ath_fft_data *data);
void ath_3d_fft_destroy_plan(struct ath_3d_fft_plan *ath_plan);

/**************************************************************************
 *
 *  Athena 2D FFT functions
 *
 **************************************************************************/

struct ath_2d_fft_plan *ath_2d_fft_quick_plan(DomainS *pD,
				ath_fft_data *data, ath_fft_direction dir);
struct ath_2d_fft_plan *ath_2d_fft_create_plan(int gnx2, int gnx1,
				int gjs, int gje, int gis, int gie,
				ath_fft_data *data, int al,
				ath_fft_direction dir);
ath_fft_data *ath_2d_fft_malloc(struct ath_2d_fft_plan *ath_plan);
void ath_2d_fft(struct ath_2d_fft_plan *ath_plan, ath_fft_data *data);
void ath_2d_fft_free(ath_fft_data *data);
void ath_2d_fft_destroy_plan(struct ath_2d_fft_plan *ath_plan);

#endif /* FFT_ENABLED */

#endif /* ATHENA_FFT_H */
