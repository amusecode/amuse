# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90
FC      = $(MPIFC)
CODEDIR = src

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o 

CODELIB = src/libflash.a

CODE_GENERATOR ?= $(AMUSE_DIR)/build.py

drive_loc = ./src/FLASH4.2.2/object
#drive_loc = ./src/FLASH4.2-mod/object
#drive_loc = ./src/flash4-ita-amnh-amuse/object

#flash_obj =  $(drive_loc)/Driver_interface.o $(drive_loc)/Driver_evolveFlash.o $(drive_loc)/Driver_initParallel.o $(drive_loc)/Driver_initFlash.o
allflash_obj =  $(wildcard $(drive_loc)/*.o)
ignore_obj = $(drive_loc)/Flash.o $(drive_loc)/UNIT_TEST.o
flash_obj = $(filter-out  $(ignore_obj), $(allflash_obj))

all: flash_worker 

ifdef DOWNLOAD_CODES
$(CODEDIR)/makefile:
	make -C . download
else
$(CODEDIR)/makefile:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set. FLASH will not be downloaded and build."
	@echo "If you do want FLASH, set DOWNLOAD_CODES to 1."
	@echo "bash> export DOWNLOAD_CODES=1"
	@echo "csh> setenv DOWNLOAD_CODES 1"
	@echo "*** Unfortunately due to flash download restrictions;"
	@echo "*** FLASH CANNOT be download automatically by this makefile"
	@echo "*** Please register and download flash yourself"
	@echo "*** You can do so at http://flash.uchicago.edu/site/flashcode/coderequest/"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error
endif

download:
	@echo ""
	@echo "*** Unfortunately due to flash download restrictions;"
	@echo "*** FLASH CANNOT be download automatically by this makefile"
	@echo "*** Please register and download flash yourself"
	@echo "*** You can do so at http://flash.uchicago.edu/site/flashcode/coderequest/"
	@echo ""
	@make -s --no-print-directory -C . raise_error

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h #flash.par
	$(RM) *~ flash_worker worker_code.f90
	make -C src clean

$(CODELIB): $(CODEDIR)/makefile
	make -C src all

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py FlashInterface -o $@

flash_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(FCFLAGS) $(FS_FLAGS) $(LDFLAGS) $< $(OBJS)  $(FS_LIBS) -o $@ $(flash_obj) $(HDF5_LIBS)
	cp $(drive_loc)/flash.par .

%.o: %.F90
	$(FC) $(FCFLAGS) -I$(drive_loc) -c -o $@ $<
