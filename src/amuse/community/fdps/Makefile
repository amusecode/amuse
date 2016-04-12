# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../../../../amuse
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libfdps.a
CODEDIR = src/FDPS
LIB_OBJ = $(CODEDIR)/sample/nbody/nbody.o

DOWNLOAD_FROM_WEB = $(PYTHON) ./download_http.py
CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: fdps_worker 

ifdef DOWNLOAD_CODES
$(CODEDIR)/sample/nbody/Makefile:
	make -C . download
else
$(CODEDIR)/sample/nbody/Makefile:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set. FDPS will not be downloaded and built."
	@echo "If you do want FDPS, set DOWNLOAD_CODES to 1."
	@echo "bash> export DOWNLOAD_CODES=1"
	@echo "csh> setenv DOWNLOAD_CODES 1"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error
endif

download:
	$(DOWNLOAD_FROM_WEB)
	tr -d '\r' < src/FDPS/sample/nbody/Makefile > src/FDPS/sample/nbody/Makefile.tmp
	mv src/FDPS/sample/nbody/Makefile.tmp src/FDPS/sample/nbody/Makefile
	$(PATCH_FILES)

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ fdps_worker worker_code.cc
	$(RM) -Rf $(CODEDIR)/sample/nbody/*.o
	$(RM) -Rf src/libfdps.a
	#make -C $(CODEDIR) clean

$(CODELIB): $(CODEDIR)/sample/nbody/Makefile
	make -C $(CODEDIR)/sample/nbody nbody.o CC=$(MPICXX) 
	ar rv $(CODELIB) $(LIB_OBJ)
	ranlib $(CODELIB)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py FDPSInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py FDPSInterface -o $@

fdps_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) -I$(CODEDIR)/src -I$(CODEDIR)/sample/nbody $< $(OBJS) $(CODELIB) -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -I$(CODEDIR)/src -I$(CODEDIR)/sample/nbody -c -o $@ $< 
