# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx

CFLAGS   += -Isrc/include -Isrc/include/star
CXXFLAGS += $(CFLAGS)
LDFLAGS  += -Lsrc/node -lnode -Lsrc/std -lstd -lm

AMUSE_OBJS = interface.o

CODELIB = src/sstar/libsstar.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py
DOWNLOAD_FROM_SVN = ./downloadsvn.py

all: 
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "The SEBA build is currently broken"
	@echo "The SEBA svn is not available"
	@echo "To fix it we need a way to download the seba code first"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

src:
ifdef DOWNLOAD_CODES
		$(DOWNLOAD_FROM_SVN)
else
		@echo ""
		@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@echo ""
		@echo "DOWNLOAD_CODES is not set. SeBa will not be downloaded and build."
		@echo "If you do want SeBa, set DOWNLOAD_CODES to 1."
		@echo "bash> export DOWNLOAD_CODES=1"
		@echo "csh> setenv DOWNLOAD_CODES 1"
		@echo ""
		@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@echo ""
		@make -s --no-print-directory -C . raise_error
endif

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h
	$(RM) *~ worker_code

SEBAclean: clean
	@-make -C src clean

distclean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h
	$(RM) *~ worker_code
	$(RM) -Rf src

$(CODELIB): src
	cd src; make Makefile.inc
	cd src; make all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py SebaInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py SebaInterface -o $@

seba_worker: worker_code.cc worker_code.h $(CODELIB) $(AMUSE_OBJS)
	$(MPICXX) $(CXXFLAGS) $< $(AMUSE_OBJS) $(CODELIB) $(LDFLAGS) -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<
