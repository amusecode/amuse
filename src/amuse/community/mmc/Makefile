MPIF90 ?= mpif90
FC      = $(MPIF90)

FFLAGS   += -Wall -g -fno-automatic
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libmmc.a

AMUSE_DIR?=../../../..
MMC_DATA_DIR = ../../../../data/mmc
CODE_GENERATOR = $(AMUSE_DIR)/build.py

SRCINC = -I./src

all: compile interface.o worker_code 

clean:
	$(RM) -f *.o *.pyc worker_code.f90  
	$(RM) *~ worker_code worker_code.f90 interface.f90
	$(RM) $(MMC_DATA_DIR)/*.dat
	make -C src clean

compile:
	make -C src

interface.f90 /src/amuse_interface2.f: makeinterface.py
	python makeinterface.py

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py mmcInterface -o $@

worker_code: worker_code.f90 $(CODELIB) 
	$(MPIF90) $(CXXFLAGS) $(OBJS) $(SRCINC) $< $(CODELIB) -o $@

interface.o: interface.f90 
	$(FC) $(FFLAGS) $(SRCINC) -c -o $@ $<

