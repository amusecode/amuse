# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90
FC      = $(MPIFC)

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIBDIR = src/phantom/bin
CODELIB = src/phantom/bin/libphantom-amuse.a

#FCFLAGS += -L$(CODELIBDIR) -lphantom-amuse

CODEDIR = src/phantom
BUILDDIR = $(CODEDIR)/build

SETUP = amuse
SYSTEM = gfortran

DOWNLOAD_FROM_WEB = $(PYTHON) ./download.py


all: phantom_worker 

$(CODEDIR)/Makefile: download
	#make -C . download

download:
	echo "not downloading"
	#$(RM) -Rf src
	#mkdir src
	#$(DOWNLOAD_FROM_WEB)

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code worker_code.f90
	$(RM) $(CODELIB)
	make -C $(CODEDIR) clean

$(CODELIB): $(CODEDIR)/Makefile
	make -C $(CODEDIR) libphantom-amuse SETUP=$(SETUP) SYSTEM=$(SYSTEM) OPENMP=yes FC="$(FC)"

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py PhantomInterface -o $@

interface.o: interface.f90
	$(MPIFC) $(FCFLAGS) $(SC_FLAGS) -I$(BUILDDIR) -c -o $@ $<

phantom_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(FCFLAGS) $(SC_FLAGS) $(FS_FLAGS) $(LDFLAGS) -I$(BUILDDIR) $< interface.o -o $@ $(CODELIB) $(SC_FCLIBS) $(FS_LIBS) -fopenmp

%.o: %.f90
	$(FC) $(SC_FLAGS) $(FCFLAGS) -c -o $@ $<
