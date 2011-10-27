#!/bin/sh

#flags=-DNGB
flags=""

CUDAINC="-I/usr/local/cuda/include -I/usr/local/cuda-sdk/common/inc"
CUDALIB="-L/usr/local/cuda/lib64 -lcuda"
# #  -L/usr/local/cuda-sdk/lib -lcutil"

# # nvcc -O0 -g  --device-emulation $flags -maxrregcount=32  -o host_evaluate_gravity.cu_o -c host_evaluate_gravity.cu -I/home/nvidia/NVIDIA_CUDA_SDK/common/inc/

# nvcc -O0 -g -D_DEBUG $flags -maxrregcount=32  -o host_evaluate_gravity.cu_o -c host_evaluate_gravity.cu $CUDAINC

# g++ -O3 $flags -g -c GPUWorker.cc $CUDAINC
# g++ -O3 $flags -g -c sapporo.cpp $CUDAINC
# g++ -O3 $flags -g -c send_fetch_data.cpp $CUDAINC
# g++ -O3 $flags -g -c sapporoG6lib.cpp $CUDAINC
# /bin/rm -rf libsapporo.a
# ar qv ./libsapporo.a sapporo.o send_fetch_data.o sapporoG6lib.o host_evaluate_gravity.cu_o GPUWorker.o
# ranlib ./libsapporo.a 

g++ -O3 $flags -g -o test_gravity_block test_gravity_block.cpp -lsapporo -L.  $CUDAINC $CUDALIB -fopenmp  -lOpenCL

g++ -O3 $flags -g -o test_gravity_block_6th test_gravity_block_6th.cpp -lsapporo -L.  $CUDAINC $CUDALIB -fopenmp  -lOpenCL


#g++ -O3 $flags -g -o test_jeroen test_jeroen.cpp -L. -lsapporo $CUDAINC $CUDALIB -lboost_thread-mt

#g++ -O3 $flags -g -o test_gravity_N2ngb test_gravity_N2ngb.cpp -L. -lsapporo $CUDAINC $CUDALIB -lboost_thread

# g++ -O3 $flags -g -o test_gravity_N2g6 test_gravity_N2g6.cpp -L. -lsapporo $CUDAINC $CUDALIB -lboost_thread -lpthread


