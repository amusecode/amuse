# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90
FC      = $(MPIFC)
AR = ar ruv
RANLIB = ranlib
RM = rm

FCFLAGS+= -Isrc/SYCLIST
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODEDIR = src/SYCLIST
CODELIB = $(CODEDIR)/libsyclist.a


all: syclist_worker 

clean:
	$(RM) -rf __pycache__
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code worker_code.f90
	make -C $(CODEDIR) clean

distclean: clean
	make -C $(CODEDIR) distclean

$(CODELIB):
	make -C $(CODEDIR) libsyclist.a

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py SyclistInterface -o $@

syclist_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(FCFLAGS) $(FS_FLAGS) $< $(OBJS) $(CODELIB) $(FS_LIBS) -o $@

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<
