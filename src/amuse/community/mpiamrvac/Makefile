AMUSE_DIR?=../../../..
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90
PYTHON ?= python
FC      = $(MPIFC)

FFLAGS   +=  -g -ffree-form
LDFLAGS  += -lm $(MUSE_LD_FLAGS)


BUILDDIR = build
TWOD_BUILDDIR = build_2d

OBJS = $(BUILDDIR)/interface.o
OBJS2D = $(TWOD_BUILDDIR)/interface.o

LIBNAME = libmpiamrvac.a

CODEDIR = src/mpiamrvac/src

CODELIB = $(BUILDDIR)/$(LIBNAME)
CODELIB2D = $(TWOD_BUILDDIR)/$(LIBNAME)

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
	$(RM) -f mpiamrvac_worker_2d
	-make -C $(CODEDIR) clean
	$(RM) -Rf $(BUILDDIR)
	$(RM) -Rf $(TWOD_BUILDDIR)
	$(RM) -Rf *.dSYM
	
distclean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h *.mod
	$(RM) *~ mpiamrvac_worker worker_code.f90
	$(RM) -f mpiamrvac_worker_2d
	$(RM) -Rf src
	$(RM) -Rf src.*
	$(RM) -Rf .pc

$(BUILDDIR)/makefile: $(CODEDIR)/makefile
	cp -r $(CODEDIR) $(BUILDDIR)

$(TWOD_BUILDDIR)/makefile: $(CODEDIR)/makefile
	cp -r $(CODEDIR) $(TWOD_BUILDDIR)
    
$(CODELIB): $(BUILDDIR)/makefile
	cd $(BUILDDIR); ./setamrvac -d=33 -g=14,14,14 -phi=3 -z=2 -p=hd -u=amuse -cp=openmpi
	make -C $(BUILDDIR) clean amrvac $(LIBNAME) F90FLAGS="$(FFLAGS)" F90="$(MPIFC)"

$(CODELIB2D): $(TWOD_BUILDDIR)/makefile
	cd $(TWOD_BUILDDIR); ./setamrvac -d=22 -g=14,14 -phi=3 -z=2 -p=hd -u=amuse -cp=openmpi
	make -C $(TWOD_BUILDDIR) clean amrvac $(LIBNAME) F90FLAGS="$(FFLAGS)" F90="$(MPIFC)"
	
worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py MpiAmrVacInterface -o $@

mpiamrvac_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@
	
mpiamrvac_worker_2d: worker_code.f90 $(CODELIB2D) $(OBJS2D)
	$(MPIFC) $(CXXFLAGS) $< $(OBJS2D) $(CODELIB2D) -o $@
	

$(BUILDDIR)/interface.f90: interface.t
	$(BUILDDIR)/vacpp.pl $< > $@
    
$(BUILDDIR)/interface.o: $(BUILDDIR)/interface.f90
	$(FC) $(FFLAGS) -I$(BUILDDIR) -c -o $@ $<
    

$(TWOD_BUILDDIR)/interface.f90: interface.t
	$(TWOD_BUILDDIR)/vacpp.pl $< > $@
    
$(TWOD_BUILDDIR)/interface.o: $(TWOD_BUILDDIR)/interface.f90
	$(FC) $(FFLAGS) -I$(TWOD_BUILDDIR) -c -o $@ $<
    

%.o: %.f90
	$(FC) $(FFLAGS) -I$(CODEDIR) -c -o $@ $<
