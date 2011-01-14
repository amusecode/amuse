IGNORE:
CMPLR ?= gfortran
FFLAGS ?= -fPIC
LFLAGS = const_bse.h zdata.h 

.f.o:
	$(CMPLR) -c $(FFLAGS) $<

SRCE1 = \
sse.f deltat.f evolv1.f hrdiag.f kick.f mlwind.f mrenv.f \
ran3.f star.f zcnsts.f zfuncs.f

OBJT1 = $(SRCE1:.f=.o)

sse:    $(OBJT1) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT1) -o sse 

SRCE2 = \
bse.f comenv.f corerd.f deltat.f dgcore.f evolv2.f gntage.f \
hrdiag.f instar.f kick.f mix.f mlwind.f mrenv.f ran3.f rl.f \
star.f zcnsts.f zfuncs.f
 
OBJT2 = $(SRCE2:.f=.o)

bse:    $(OBJT2) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT2) -o bse 

SRCE3 = \
popbin.f comenv.f corerd.f deltat.f dgcore.f evolv2.f gntage.f \
hrdiag.f instar.f kick.f mix.f mlwind.f mrenv.f ran3.f rl.f \
star.f zcnsts.f zfuncs.f
 
OBJT3 = $(SRCE3:.f=.o)

popbin: $(OBJT3) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT3) -o popbin

