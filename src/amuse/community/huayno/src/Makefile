# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../../..
-include ${AMUSE_DIR}/config.mk


MPICXX ?= mpicxx
MPICC ?= mpicc

OPENMP_CFLAGS?= -fopenmp

ifneq (,$(findstring xlc, $(notdir $(CC))))
CFLAGS += -qlanglvl=extc99
else
CFLAGS += -std=gnu99
endif
LIBS += -lm
INCLUDE =

AR = ar ruv
RANLIB = ranlib
RM = rm

#---------------------------------------------

OBJS = evolve.o evolve_shared.o evolve_sf.o evolve_cc.o \
  evolve_ok.o evolve_kepler.o universal_variable_kepler.o evolve_bs.o \
  evolve_shared_collisions.o 

all: libhuayno.a

clean:
	rm -f *.o *.bck *.pyc *.clh worker_code.cc worker_code.h worker_code-sockets.cc

distclean: clean
	rm -f huayno_worker huayno_worker_cl huayno_worker_mp huayno_worker_sockets

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py HuaynoInterface -o $@
	
worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py HuaynoInterface -o $@


libhuayno.a: $(OBJS)
	$(RM) -f $@
	$(AR) $@ $(OBJS)
	$(RANLIB) $@
	
libhuayno_cl.a: evolve_kern.clh $(OBJS) evolve_cl.o 
	$(RM) -f $@
	$(AR) $@ $(OBJS) evolve_cl.o 
	$(RANLIB) $@
	
.cc.o: $<
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(INCLUDE) -c -o $@ $< 

.c.o: $<
	$(MPICC) $(CFLAGS) $(SC_FLAGS) -I../build_cl $(INCLUDE) -c -o $@ $< 

%.clh: %.cl
	awk 'BEGIN{print "const char srcstr[]=" } {if(substr($$0,length($0))=="\\"){$$0=$$0"\\"};print "\""$$0"\\n\""} END{print ";"}' $< > $@
