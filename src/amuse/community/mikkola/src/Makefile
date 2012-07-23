# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../../..
-include ${AMUSE_DIR}/config.mk

MPIFC ?= mpif90
FC    = $(MPIFC)

FFLAGS   += -O3
FFLAGS   += $(F77FLAGS)
LDFLAGS  += -lm 

CODELIB = libMikkola.a

CODEOBJS = Mikkola.o

AR = ar ruv
RANLIB = ranlib
RM = rm

all: $(CODELIB) 

clean:
	$(RM) -f *.o *.a

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

