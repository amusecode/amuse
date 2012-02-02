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





//Defines taken from the cutil header files
//

#if CUDART_VERSION >= 4000
#define CUT_DEVICE_SYNCHRONIZE( )   cudaDeviceSynchronize();
#else
#define CUT_DEVICE_SYNCHRONIZE( )   cudaThreadSynchronize();
#endif




#  define CUDA_SAFE_CALL_NO_SYNC( call) {                                    \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } }

#  define CUDA_SAFE_CALL( call)     CUDA_SAFE_CALL_NO_SYNC(call);                        


    //! Check for CUDA error
#ifdef _DEBUG
#  define CUT_CHECK_ERROR(errorMessage) {                                    \
    cudaError_t err = cudaGetLastError();                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    err = CUT_DEVICE_SYNCHRONIZE();                                           \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    }
#else
#  define CUT_CHECK_ERROR(errorMessage) {                                    \
    cudaError_t err = cudaGetLastError();                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    }
#endif






 
#endif





