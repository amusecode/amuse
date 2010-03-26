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
	$(DOWNLOAD_FROM_SVN)

clean: src
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code
	cd src; make clean
	
distclean: src
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code
	$(RM) -Rf src

$(CODELIB): src
	cd src; make Makefile.inc
	cd src; make all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py sebaInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py sebaInterface -o $@

worker_code: worker_code.cc worker_code.h $(CODELIB) $(AMUSE_OBJS)
	mpicxx $(CXXFLAGS) $@.cc $(AMUSE_OBJS) $(CODELIB) $(LDFLAGS) -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
