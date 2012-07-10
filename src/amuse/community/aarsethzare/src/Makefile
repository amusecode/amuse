FC ?= mpif90

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CODELIB = libaarsethzare.a

CODEOBJS = AarsethZare.o transform.o

AR = ar ruv
RANLIB = ranlib
RM = rm

ifneq (,$(findstring xlf, $(notdir $(FORTRAN))))
FCFLAGS += -qfixed
endif


all: $(CODELIB) 

clean:
	$(RM) -f *.o *.a

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<
%.o: %.f
	$(FC) $(FCFLAGS) -c -o $@ $<

