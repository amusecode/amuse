IGNORE:
CMPLR ?= gfortran
FFLAGS ?= -fPIC
LFLAGS = const_mobse.h zdata.h 

.f.o:
	$(CMPLR) -c $(FFLAGS) $<
 
SRCE1 = \
mobse.f comenv.f corerd.f deltat.f dgcore.f evolv3.f gntage.f \
hrdiag.f instar.f kick.f mix.f mlwind.f mrenv.f ran3.f rl.f \
star.f zcnsts.f zfuncs.f pisn.f eddington.f fallback.f

OBJT1 = $(SRCE1:.f=.o)

mobse:	$(OBJT1) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT1) -o mobse  
 