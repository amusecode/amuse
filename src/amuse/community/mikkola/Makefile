MPIF90 ?= mpif90
FC      = $(MPIF90)

FFLAGS   += -Wall -g
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libMikkola.a

AMUSE_DIR?=../../amuse

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: mikkola_worker 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code worker_code.f90
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py MikkolaInterface -o $@

mikkola_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIF90) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<
