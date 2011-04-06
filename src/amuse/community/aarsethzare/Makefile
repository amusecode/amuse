MPIF90 ?= mpif90
FC      = $(MPIF90)

FFLAGS   += -Wall -g
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libaarsethzare.a

AMUSE_DIR ?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: aarsethzare_worker 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ aarsethzare_worker worker_code.f90
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py AarsethZareInterface -o $@

aarsethzare_worker: worker_code.f90 $(CODELIB) $(OBJS)
	mpif90 $(CXXFLAGS) worker_code.f90 $(OBJS) $(CODELIB) -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<
