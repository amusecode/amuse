MPIF90 ?= mpif90
FC      = $(MPIF90)

FFLAGS   += -Wall -g
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libmmc.a

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: worker_code $(CODELIB)

clean:
	$(RM) -f *.o *.pyc worker_code.f90  
	$(RM) *~ worker_code worker_code.f90
	make -C src clean

$(CODELIB):
	make -C src

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py mmcInterface -o $@

worker_code: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIF90) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<
