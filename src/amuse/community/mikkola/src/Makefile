FC = mpif90

FFLAGS   += -Wall -g -O4
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

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

