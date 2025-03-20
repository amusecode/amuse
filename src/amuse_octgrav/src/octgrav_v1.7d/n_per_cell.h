#ifndef N_PER_CELL

#define N_PER_CELL 8
#define NCRIT   128
#define NTHREADS 256
#define NBLOCKS  128

extern "C" {
  extern int  cuda_interaction_list_len;
  extern int  cuda_interaction_node_len;
  extern int  cuda_interaction_leaf_len;
  extern int  cuda_interaction_node_list;
  extern int  cuda_interaction_leaf_list;
  extern int  cuda_n_node;
  extern int  cuda_n_leaf;
  extern int  n_alloc;

  extern int3 *dev_interaction_list_len;
  extern int2 *dev_interaction_node_len;
  extern int2 *dev_interaction_leaf_len;
  extern int  *dev_interaction_node_list;
  extern int  *dev_interaction_leaf_list;
  extern int  *dev_n_node;
  extern int  *dev_n_leaf;
}


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
