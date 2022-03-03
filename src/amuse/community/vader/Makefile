ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include ${AMUSE_DIR}/config.mk

MPICXX   ?= mpicxx
MPICC    ?= mpicc

CFLAGS   += -Wall -g 
CXXFLAGS += $(CFLAGS) $(GSL_FLAGS)
LDFLAGS  += $(GSL_LIBS) $(MUSE_LD_FLAGS) 

OBJS = interface.o

DOWNLOAD_FROM_WEB = $(PYTHON) ./download.py

PROB = none

BUILDDIR = build

.PHONY: all code download

ifdef MAKE_CODE

TARGETS := $(patsubst src/prob/userFunc_%.c, vader_worker_%, ${wildcard src/prob/userFunc_*.c})
all: $(TARGETS)

${TARGETS}: vader_worker_%: worker_code.cc worker_code.h interface_%.o $(BUILDDIR)_%/libvader.a
	$(MPICXX) $(CXXFLAGS) $< interface_$*.o $(BUILDDIR)_$*/libvader.a -o $@ $(LDFLAGS)

else

all: src code

endif

code:
	make -C . MAKE_CODE=True


download:
	$(DOWNLOAD_FROM_WEB)

src: vader_src
	cp -r vader_src/vader_csrc src
	cp Makefile_vader src/Makefile


ifdef DOWNLOAD_CODES
vader_src:
	make -C . download
else
vader_src:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set. VADER will not be downloaded and build."
	@echo "If you do want VADER, set DOWNLOAD_CODES to 1."
	@echo "bash> export DOWNLOAD_CODES=1"
	@echo "csh> setenv DOWNLOAD_CODES 1"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@exit -1
endif 


clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ vader_worker_* worker_code.cc
	$(RM) -Rf build_*
#~ make -C src clean

vader_worker: download src code

$(BUILDDIR)_%: src
	mkdir $@
	
$(BUILDDIR)_%/Makefile: src/Makefile | $(BUILDDIR)_%
	cp src/Makefile $@
	
$(BUILDDIR)_%/userFunc.h: src/userFunc.h | $(BUILDDIR)_%
	cp src/userFunc.h $@ 

$(BUILDDIR)_%/vader_common.h: src/vader_common.h | $(BUILDDIR)_%
	cp src/vader_common.h $@ 
	
$(BUILDDIR)_%/libvader.a: $(BUILDDIR)_%/Makefile $(BUILDDIR)_%/userFunc.h $(BUILDDIR)_%/vader_common.h
	make -C $(BUILDDIR)_$* all PROB=$* VPATH=../src

interface_%.o: interface.cc worker_code.h
	$(MPICC) $(CXXFLAGS) -c -o $@ $<

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py VaderInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py VaderInterface -o $@

.cc.o: $<
	$(MPICC) $(CXXFLAGS) -c -o $@ $< 

distclean: clean
	$(RM) -Rf src vader_src *.tar.gz

.PRECIOUS: $(BUILDDIR)_% $(BUILDDIR)_%/Makefile $(BUILDDIR)_%/userFunc.h $(BUILDDIR)_%/vader_common.h
