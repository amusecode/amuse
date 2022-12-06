# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90
FC = $(MPIFC)

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

FCFLAGS += -Isrc/GENEC/code
# FCFLAGS += -O0

OBJS = interface.o

CODEDIR = src/GENEC/code
INI_CODEDIR = src/GENEC/tools/makeini

CODELIB = $(CODEDIR)/libgenec.a

all: genec_worker 

clean:
	$(RM) -rf __pycache__
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ genec_worker worker_code.f90
	$(RM) $(CODELIB)
	$(RM) $(CODEDIR)/*.o
	$(RM) $(INI_CODEDIR)/*.o
	$(RM) *.mod
	make -C $(CODEDIR) clean

distclean: clean
	make -C src distclean

$(CODELIB): #$(OBJECTS)
	make -C $(CODEDIR) libgenec.a
	# @$(COMPILER) -shared $(FFLAGS) $(OBJECTS) $(LDFLAGS) -o $(CODELIB)

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py GenecInterface -o $@

genec_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(FCFLAGS) $(FS_FLAGS) $< $(OBJS) -o $@ $(CODELIB) $(FS_LIBS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<
