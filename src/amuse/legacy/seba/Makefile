CFLAGS   += -Isrc/include -Isrc/include/star
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -Lsrc/node -lnode -Lsrc/std -lstd -lm

AMUSE_OBJS = interface.o

CODELIB = src/sstar/libsstar.a

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py
DOWNLOAD_FROM_SVN = ./downloadsvn.py

all: worker_code 

src:
ifdef DO_INSTALL_SEBA
		$(DOWNLOAD_FROM_SVN)
else
		@echo ""
		@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@echo ""
		@echo "DO_INSTALL_SEBA is not set. SeBa will not be downloaded and build."
		@echo "If you do want SeBa, set DO_INSTALL_SEBA to 1."
		@echo "bash> export DO_INSTALL_SEBA=1"
		@echo "csh> setenv DO_INSTALL_SEBA 1"
		@echo ""
		@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@echo ""
		@make -s --no-print-directory -C . raise_error
endif

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code
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

worker_code: worker_code.cc worker_code.h $(CODELIB) $(AMUSE_OBJS)
	mpicxx $(CXXFLAGS) $@.cc $(AMUSE_OBJS) $(CODELIB) $(LDFLAGS) -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
	
	
