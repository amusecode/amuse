#======== modify the following values ==================

# The main files to be compiled
PROGRAMS = genOperations
#========= the following commands should work on most Unix systems ========

# The main CUDA root path
CUDA_PATH = /usr/local/cuda-7.5

# The Name of the C++/CUDA compiler
CXX=g++
NVCC=$(CUDA_PATH)/bin/nvcc

# Optional flags to give to the C++ compiler
CXXOPTS=

# Extra flags to give to the C++ compiler
CXXFLAGS=$(CXXINC) $(CXXOPTS)

# Additional include paths
CUDASDK = $(CUDA_PATH)/samples/common/inc
CUDAINC = -I$(CUDA_PATH)/include -I$(CUDASDK)

# Additional library paths
CUDALIB = -L$(CUDA_PATH)/lib64 -lcudart
MPFRLIB = -lmpfr
MPFILIB = -lmpfi
PTHREADLIB = -L/lib64 -lpthread
THREADPOOL = ./util/cthpool/thpool.c
#========== you shouldn't modify anything below ===========================

.SUFFIXES:
.SUFFIXES: .cpp .cc .hpp .o

default:
	@echo
	@echo
	@echo "Usage: make all | mw_henon_multi_gpu_0D_dbl ..."
	@echo
	@echo "    tidy / clean                 (tidy / clean up)"
	@echo "";

# -------------------------------------------------------------------------

all: $(PROGRAMS)

genOperations: genOperations.cpp
	$(CXX) -lglut -O2 -g -o $@ $@.cpp 
	
.PHONY: default all

# -----------------------------------------------------------------------
# TIDY UP
# -----------------------------------------------------------------------

tidy:
	@ - rm -f *.o *~

#-----------------------------------------------------------------------
# CLEAN UP
# -----------------------------------------------------------------------

clean:
	@ - rm -f *.o *~ $(PROGRAMS)

# -----------------------------------------------------------------------

