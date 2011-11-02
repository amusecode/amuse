AMUSE_DIR?=../../../..
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90
PYTHON ?= python
FC      = $(MPIFC)

FFLAGS   +=  -g -ffree-form
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o
OBJS2D = interface_2d.o

LIBNAME = libmpiamrvac.a
LIBNAME2D = libmpiamrvac2d.a


CODEDIR = src/mpiamrvac/src
CODELIB = src/mpiamrvac/src/$(LIBNAME)
CODELIB2D = src/mpiamrvac/src/$(LIBNAME2D)

CODE_GENERATOR = $(AMUSE_DIR)/build.py

DOWNLOAD_FROM_WEB = $(PYTHON) ./download_svn.py
PATCH_FILES = $(PYTHON) ./patch_files.py
all: mpiamrvac_worker 

ifdef DOWNLOAD_CODES
$(CODEDIR)/makefile:
	make -C . download
else
$(CODEDIR)/makefile:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set. MpiAmrVac will not be downloaded and build."
	@echo "If you do want MpiAmrVac, set DOWNLOAD_CODES to 1."
	@echo "bash> export DOWNLOAD_CODES=1"
	@echo "csh> setenv DOWNLOAD_CODES 1"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error
endif

download:
	$(RM) -Rf .pc
	$(RM) -Rf src
	mkdir src
	$(DOWNLOAD_FROM_WEB)
	$(PATCH_FILES)

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h *.mod
	$(RM) *~ mpiamrvac_worker worker_code.f90
	-make -C $(CODEDIR) clean
	$(RM) -f $(CODELIB) $(CODELIB2D)
	
distclean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h *.mod
	$(RM) *~ mpiamrvac_worker worker_code.f90
	$(RM) -Rf src
	$(RM) -Rf .pc

$(CODELIB): $(CODEDIR)/makefile
	cd src/mpiamrvac/src; ./setamrvac -d=33 -g=14,14,14 -phi=3 -z=2 -p=hd -u=amuse -cp=openmpi
	make -C $(CODEDIR) clean amrvac $(LIBNAME) F90FLAGS="$(FFLAGS)"

$(CODELIB2D): $(CODEDIR)/makefile
	cd src/mpiamrvac/src; ./setamrvac -d=22 -g=14,14 -phi=3 -z=2 -p=hd -u=amuse -cp=openmpi
	make -C $(CODEDIR) clean amrvac $(LIBNAME2D) F90FLAGS="$(FFLAGS)"
	
worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py MpiAmrVacInterface -o $@

mpiamrvac_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@
	
mpiamrvac_worker_2d: worker_code.f90 $(CODELIB2D) $(OBJS2D)
	$(MPIFC) $(CXXFLAGS) $< $(OBJS2D) $(CODELIB2D) -o $@
	


%.o: %.f90
	$(FC) $(FFLAGS) -I$(CODEDIR) -c -o $@ $<
