FC ?= mpif90

FCFLAGS   += -Wall -g
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CODELIB = libaarsethzare.a

CODEOBJS = AarsethZare.o transform.o

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
	$(FC) $(FCFLAGS) -c -o $@ $<

