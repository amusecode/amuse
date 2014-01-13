# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx
MPICC ?= mpicc
PYTHON ?= python

VERSION = 4.1

ATHENA_DIR = src/athena

CC = $(MPICC)
CFLAGS   += -DAMUSE -g  
CXXFLAGS += $(CFLAGS) 
UNAME = $(shell uname)

ifneq (,$(findstring MINGW,$(UNAME)))
  LDFLAGS += -lm
else
  LDFLAGS += -ldl -lm
endif


BUILDDIR = build
SELFGRAV_BUILDDIR = build_selfgrav
MHD_BUILDDIR = build_mhd
SCALAR_BUILDDIR = build_scalar

CODELIB = $(BUILDDIR)/bin/libathena.a
SELFGRAV_CODELIB = $(SELFGRAV_BUILDDIR)/bin/libathena.a
MHD_CODELIB = $(MHD_BUILDDIR)/bin/libathena.a
SCALAR_CODELIB = $(SCALAR_BUILDDIR)/bin/libathena.a

AMUSE_OBJS = $(BUILDDIR)/interface_$(VERSION).o
SELFGRAV_AMUSE_OBJS = $(SELFGRAV_BUILDDIR)/interface_$(VERSION).o
MHD_AMUSE_OBJS = $(MHD_BUILDDIR)/interface_$(VERSION).o
SCALAR_AMUSE_OBJS = $(SCALAR_BUILDDIR)/interface_$(VERSION).o

CODE_GENERATOR = $(AMUSE_DIR)/build.py
DOWNLOAD_FROM_WEB = $(PYTHON) ./downloadhttp.py
PATCH_FILES = $(PYTHON) ./patch_files.py

ATHENA_OBJS = $(subst $(ATHENA_DIR)/src/main.o,,$(wildcard $(ATHENA_DIR)/src/*.o $(ATHENA_DIR)/src/*/*.o))

FFTW_LIBS ?= -L/usr/lib -lfftw3
FFTW_FLAGS ?= -I/usr/include
FFTWLIB = $(FFTW_LIBS)
FFTWINC = $(FFTW_FLAGS)
 
.PHONY : download

all: athena_worker 

clean:
	$(RM) -f *.so *.o *.selfgrav_o *.pyc worker_code.cc worker_code.h 
	$(RM) -f *~ athena_worker  athena_worker_selfgrav athena_worker_mhd
	$(RM) -f athena_worker*.exe
	$(RM) -f athena_worker_scalar
	$(RM) -f worker_code-sockets.cc athena_worker_sockets
	$(RM) -f $(CODELIB) $(AMUSE_OBJS_SELFGRAV)
	$(RM) -f *.mhd_o
	$(RM) -rf *.dSYM
	$(RM) -rf $(BUILDDIR) $(SELFGRAV_BUILDDIR) $(MHD_BUILDDIR)
	$(RM) -rf $(SCALAR_BUILDDIR)
	
codeclean: clean
	make -C $(ATHENA_DIR) clean

distclean: clean
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) -Rf .pc
	$(RM) *~ athena_worker

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py AthenaInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py -i amuse.support.codes.stopping_conditions.StoppingConditionInterface AthenaInterface -o $@

$(BUILDDIR)/Makefile: $(ATHENA_DIR)/Makefile.in
	cp -r $(ATHENA_DIR) $(BUILDDIR)

$(CODELIB): $(BUILDDIR)/Makefile
ifeq ($(MPI_ENABLED),no)
	cd $(BUILDDIR); ./configure --with-gas=hydro --disable-mpi --enable-smr 
else
	cd $(BUILDDIR); ./configure --with-gas=hydro --enable-mpi --enable-smr 
endif
	cd $(BUILDDIR); ln -fs ../../amuse_problem_$(VERSION).c src/problem.c
	exit
	cd $(BUILDDIR); make all CC=$(MPICC) LDR=$(MPICC)
	rm -f $(CODELIB)
	cd $(BUILDDIR); make all CC=$(MPICC) LDR=$(MPICC)

athena_worker: worker_code.cc worker_code.h $(CODELIB) $(AMUSE_OBJS)
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) worker_code.cc $(AMUSE_OBJS) $(CODELIB) $(LDFLAGS) -o $@ $(SC_CLIBS)  $(LIBS)

$(SELFGRAV_BUILDDIR)/Makefile: $(ATHENA_DIR)/Makefile.in
	cp -r $(ATHENA_DIR) $(SELFGRAV_BUILDDIR)

$(SELFGRAV_CODELIB): $(SELFGRAV_BUILDDIR)/Makefile
ifeq ($(MPI_ENABLED),no)
	cd $(SELFGRAV_BUILDDIR); ./configure --with-gas=hydro --disable-mpi --with-gravity=fft --enable-fft --disable-cooling
else
	cd $(SELFGRAV_BUILDDIR); ./configure --with-gas=hydro --with-gravity=fft --enable-fft --disable-cooling
endif
	cd $(SELFGRAV_BUILDDIR); ln -fs ../../amuse_problem_$(VERSION).c src/problem.c
	exit
	cd $(SELFGRAV_BUILDDIR); make all CC=$(MPICC) LDR=$(MPICC)
	rm -f $(CODELIB)
	cd $(SELFGRAV_BUILDDIR); make all CC=$(MPICC) LDR=$(MPICC)
	
athena_worker_selfgrav: worker_code.cc worker_code.h $(SELFGRAV_CODELIB) $(SELFGRAV_AMUSE_OBJS)
	$(MPICXX) $(CXXFLAGS) $(FFTW_FLAGS) $(SC_FLAGS) worker_code.cc $(SELFGRAV_AMUSE_OBJS) $(SELFGRAV_CODELIB) $(LDFLAGS) -o $@ $(SC_CLIBS) $(FFTWLIB)  $(LIBS)

$(MHD_BUILDDIR)/Makefile: $(ATHENA_DIR)/Makefile.in
	cp -r $(ATHENA_DIR) $(MHD_BUILDDIR)

$(MHD_CODELIB): $(MHD_BUILDDIR)/Makefile
ifeq ($(MPI_ENABLED),no)
	cd $(MHD_BUILDDIR); ./configure --with-gas=mhd --with-eos=adiabatic --disable-mpi --enable-smr --with-order=3 --with-flux=hlld
else
	cd $(MHD_BUILDDIR); ./configure --with-gas=mhd --with-eos=adiabatic --enable-mpi --enable-smr --with-order=3 --with-flux=hlld
endif
	cd $(MHD_BUILDDIR); ln -fs ../../amuse_problem_$(VERSION).c src/problem.c
	exit
	cd $(MHD_BUILDDIR); make all CC=$(MPICC) LDR=$(MPICC)
	rm -f $(CODELIB)
	cd $(MHD_BUILDDIR); make all CC=$(MPICC) LDR=$(MPICC)

athena_worker_mhd: worker_code.cc worker_code.h $(MHD_CODELIB) $(MHD_AMUSE_OBJS)
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) worker_code.cc $(MHD_AMUSE_OBJS) $(MHD_CODELIB) $(LDFLAGS) -o $@ $(SC_CLIBS)  $(LIBS)


$(SCALAR_BUILDDIR)/Makefile: $(ATHENA_DIR)/Makefile.in
	cp -r $(ATHENA_DIR) $(SCALAR_BUILDDIR)

$(SCALAR_CODELIB): $(SCALAR_BUILDDIR)/Makefile
ifeq ($(MPI_ENABLED),no)
	cd $(SCALAR_BUILDDIR); ./configure --with-gas=hydro --disable-mpi --enable-smr --with-nscalars=1
else
	cd $(SCALAR_BUILDDIR); ./configure --with-gas=hydro --enable-mpi --enable-smr --with-nscalars=1
endif
	cd $(SCALAR_BUILDDIR); ln -fs ../../amuse_problem_$(VERSION).c src/problem.c
	exit
	cd $(SCALAR_BUILDDIR); make all CC=$(MPICC) LDR=$(MPICC)
	rm -f $(CODELIB)
	cd $(SCALAR_BUILDDIR); make all CC=$(MPICC) LDR=$(MPICC)

athena_worker_scalar: worker_code.cc worker_code.h $(SCALAR_CODELIB) $(SCALAR_AMUSE_OBJS)
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) worker_code.cc $(SCALAR_AMUSE_OBJS) $(SCALAR_CODELIB) $(LDFLAGS) -o $@ $(SC_CLIBS)  $(LIBS)

$(BUILDDIR)/%.o: %.c
	$(CC) $(CXXFLAGS) -I$(BUILDDIR)/src $(SC_FLAGS) -c -o $@ $< 
	
$(SELFGRAV_BUILDDIR)/%.o: %.c
	$(CC) $(CXXFLAGS) $(FFTW_FLAGS) -I$(SELFGRAV_BUILDDIR)/src $(SC_FLAGS) -c -o $@ $< 
	
$(MHD_BUILDDIR)/%.o: %.c
	$(CC) $(CXXFLAGS) -I$(MHD_BUILDDIR)/src $(SC_FLAGS) -c -o $@ $< 
	
$(SCALAR_BUILDDIR)/%.o: %.c
	$(CC) $(CXXFLAGS) -I$(SCALAR_BUILDDIR)/src $(SC_FLAGS) -c -o $@ $< 
	


