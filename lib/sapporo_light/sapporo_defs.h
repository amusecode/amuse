#ifndef _SAPPORODEFS_H_
#define _SAPPORODEFS_H_

#define MAXCUDADEVICES 4
#define NBODIES_MAX 524288
#define NBLOCKS 16        /* number of block which can be run simultaneously */

#ifdef NGB
#define NTHREADS 256   /* max number of threads which can run per block */
#else
#define NTHREADS 256   /* max number of threads which can run per block */
#endif

#define NGB_PP      256     /* max number of neighbours per particle */
#define NGB_PB      NGB_PP  /* max number of neighbours per particle */

 
#endif
