# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx

CFLAGS   += -Isrc/include -Isrc/include/star
CXXFLAGS += $(CFLAGS)
LDFLAGS  += -Lsrc/node -lnode -Lsrc/std -lstd -lm

unexport CFLAGS CXXFLAGS LDFLAGS

AMUSE_OBJS = interface.o

CODELIB = src/sstar/libsstar.a src/dstar/libdstar.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py
DOWNLOAD_FROM_SVN = ./downloadsvn.py

all: seba_worker

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h
	$(RM) *~ worker_code seba_worker
	@-make -C src clean

distclean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h
	$(RM) *~ worker_code

$(CODELIB): src
	cd src; make -C. CXX="$(CXX) $(CXXFLAGS)"

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py SeBaInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py SeBaInterface -o $@

seba_worker: worker_code.cc worker_code.h $(CODELIB) $(AMUSE_OBJS)
	$(MPICXX) $(CXXFLAGS) $< $(AMUSE_OBJS) $(CODELIB) $(LDFLAGS) -o $@  $(LIBS)

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<
