MPIF90 ?= mpif90
FC      = $(MPIF90)

FFLAGS   += -Wall -g
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libmmc.a

AMUSE_DIR?=../../../..
MMC_DATA_DIR = ../../../../data/mmc
CODE_GENERATOR = $(AMUSE_DIR)/build.py

SRCINC = -I./src

all: compile worker_code 

clean:
	$(RM) -f *.o *.pyc worker_code.f90  
	$(RM) *~ worker_code worker_code.f90
	$(RM) $(MMC_DATA_DIR) *.dat
	make -C src clean

compile:
	make -C src

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py mmcInterface -o $@

worker_code: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIF90) $(CXXFLAGS) $(SRCINC) $< $(OBJS) $(CODELIB) -o $@

%.o: %.f90
	$(FC) $(FFLAGS) $(SRCINC) -c -o $@ $<
