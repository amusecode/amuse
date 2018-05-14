IGNORE:
CMPLR ?= gfortran
FFLAGS += -fPIC
LFLAGS = const_mobse.h zdata.h 

.f.o:
	$(CMPLR) -c $(FFLAGS) $<
 
SRCE1 = \
mosse.f deltat.f evolv1.f  \
hrdiag.f kick.f mlwind.f mrenv.f ran3.f \
star.f zcnsts.f zfuncs.f pisn.f eddington.f fallback.f

OBJT1 = $(SRCE1:.f=.o)

mosse:	$(OBJT1) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT1) -o mosse  
 
