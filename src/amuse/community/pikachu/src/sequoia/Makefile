# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX   ?= mpicxx
CXXFLAGS ?= -Wall -g -O2


CXX = $(MPICXX)
LD  = $(MPICXX) 

.SUFFIXES: .o .cpp .ptx .cu

CUDA_TK ?= /usr/local/cuda

# OFLAGS = -O0 -g -Wall 
OFLAGS = -O3 -g -Wall -fopenmp 
CXXFLAGS +=  -fPIC $(OFLAGS) -I$(CUDA_TK)/include 



# NVCC      = $(CUDA_TK)/bin/nvcc  --device-emulation
# NVCCFLAGS = -D_DEBUG -O0 -g -I$(CUDA_SDK)/common/inc -arch=sm_12 --maxrregcount=64  --opencc-options -OPT:Olimit=0 -I$(CUDPP)/cudpp/include
NVCC ?= $(CUDA_TK)/bin/nvcc  
#NVCCFLAGS = -arch sm_20 -ftz=true -prec-div=false -prec-sqrt=false 
#NVCCFLAGS = -arch sm_20 -g -G 
NVCCFLAGS ?= -arch sm_20 


# Use with Mac OS X
# NVCCFLAGS = -arch sm_12 -Xcompiler="-Duint=unsigned\ int"

#LDFLAGS = -lcuda -lOpenCL  
LDFLAGS += -lcuda -lOpenCL -fopenmp 

#For tipsy support uncomment the following line
#CXXFLAGS += -DTIPSYOUTPUT


INCLUDEPATH = ./include
CXXFLAGS  += -I$(INCLUDEPATH)
NVCCFLAGS += -I$(INCLUDEPATH)

CUDAKERNELSPATH = CUDAkernels
#CUDAKERNELSPATH = CUDA

#CUDAKERNELS = build_tree.cu \
#	compute_properties.cu \
#	compute_propertiesD.cu \
#	dev_approximate_gravity.cu \
#	dev_approximate_gravity_let.cu \
#	dev_get_ngb.cu \
#	scanKernels.cu \
#	sortKernels.cu \
#	support_kernels.cu \
#	timestep.cu \
#	parallel.cu

CUDAKERNELS = build_tree.cu \
	compute_properties.cu \
	compute_propertiesD.cu \
	dev_approximate_gravity.cu \
	dev_get_ngb.cu \
	scanKernels.cu \
	sortKernels.cu \
	support_kernels.cu \
	timestep.cu \
	parallel.cu


CUDAPTX = $(CUDAKERNELS:%.cu=%.ptx)

SRCPATH = src
#SRC = main.cpp octree.cpp load_kernels.cpp scanFunctions.cpp build.cpp compute_properties.cpp sort_bodies_gpu.cpp gpu_iterate.cpp parallel.cpp libraryInterface.cpp
SRC = octree.cpp load_kernels.cpp scanFunctions.cpp sort_bodies_gpu.cpp sequoiaInterface.cpp build.cpp compute_properties.cpp useTreeFunctions.cpp
OBJ = $(SRC:%.cpp=%.o)

PROG = main

AR = ar ruv
RANLIB = ranlib
RM = rm

CODELIB = libsequoia.a

#all:	  $(OBJ) $(CUDAPTX) $(PROG) $(CODELIB)
all:	  $(OBJ) $(CUDAPTX) $(CODELIB)
kernels:  $(CUDAPTX)

$(CODELIB): $(OBJ) $(CUDAPTX)
	$(RM) -f $@
	$(AR) $@ $(OBJ)
	$(RANLIB) $@


$(PROG): $(OBJ)
	$(LD) $(LDFLAGS) $^ -o $@ 

%.o: $(SRCPATH)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.ptx: $(CUDAKERNELSPATH)/%.cu
	$(NVCC) $(NVCCFLAGS) -ptx $< -o $@

clean:
	$(RM) -rf *.o *.ptx main *.a

$(OBJ): $(INCLUDEPATH)/*.h

support_kernels.ptx: $(INCLUDEPATH)/node_specs.h
build_tree.ptx: $(CUDAKERNELSPATH)/support_kernels.cu $(INCLUDEPATH)/node_specs.h
compute_properties.ptx: $(CUDAKERNELSPATH)/support_kernels.cu $(INCLUDEPATH)/node_specs.h
compute_propertiesD.ptx: $(CUDAKERNELSPATH)/support_kernels.cu $(INCLUDEPATH)/node_specs.h
dev_approximate_gravity.ptx: $(CUDAKERNELSPATH)/support_kernels.cu $(CUDAKERNELSPATH)/dev_shared_traverse_functions.cu $(INCLUDEPATH)/node_specs.h
#dev_approximate_gravity_let.ptx: $(CUDAKERNELSPATH)/support_kernels.cu $(INCLUDEPATH)/node_specs.h
timestep.ptx: $(CUDAKERNELSPATH)/support_kernels.cu $(INCLUDEPATH)/node_specs.h
sortKernels.ptx: $(CUDAKERNELSPATH)/scanKernels.cu  $(CUDAKERNELSPATH)/support_kernels.cu $(INCLUDEPATH)/node_specs.h
parallel.ptx: $(CUDAKERNELSPATH)/support_kernels.cu $(INCLUDEPATH)/node_specs.h








