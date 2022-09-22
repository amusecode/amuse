# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
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
#~ INCLUDE =

AR = ar ruv
RANLIB = ranlib
RM = rm

#---------------------------------------------

OBJS = evolve.o evolve_shared.o evolve_sf.o evolve_cc.o \
  evolve_ok.o evolve_kepler.o universal_variable_kepler.o evolve_bs.o \
  evolve_shared_collisions.o evolve_error_control.o simple_map.o simple_hash.o

all: libhuayno.a

clean:
	rm -f *.o *.bck *.pyc *.clh 

distclean: clean
	
libhuayno.a: $(OBJS)
	$(RM) -f $@
	$(AR) $@ $(OBJS)
	$(RANLIB) $@
	
libhuayno_cl.a: evolve_kern.clh $(OBJS) evolve_cl.o 
	$(RM) -f $@
	$(AR) $@ $(OBJS) evolve_cl.o 
	$(RANLIB) $@
	
.cc.o: $<
	$(CXX) $(CXXFLAGS) $(SC_FLAGS) $(INCLUDE) -c -o $@ $< 

.c.o: $<
	$(CC) $(CFLAGS) $(SC_FLAGS) $(INCLUDE) -c -o $@ $< 

%.clh: %.cl
	awk 'BEGIN{print "const char *srcstr=" } {if(substr($$0,length($0))=="\\"){$$0=$$0"\\"};print "\""$$0"\\n\""} END{print ";"}' $< > $@
