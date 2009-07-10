IGNORE:
CMPLR = gfortran
FFLAGS = -fPIC
LFLAGS = const_bse.h zdata.h 

.f.o:
	$(CMPLR) -c $(FFLAGS) $<

SRCE1 = \
sse.f deltat.f evolv1.f hrdiag.f kick.f mlwind.f mrenv.f \
ran3.f star.f zcnsts.f zfuncs.f

OBJT1 = $(SRCE1:.f=.o)

sse:    $(OBJT1) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT1) -o sse 

