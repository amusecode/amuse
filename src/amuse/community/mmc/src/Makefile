.KEEP_STATE:
#FFLAGS = -O3 -fbounds-check
#FFLAGS = -O3
#FFLAGS = -O3 -ipo
#FFLAGS = -O2 -static
#FFLAGS = -g -static
#FFLAGS = -O3 -Bstatic
#FFLAGS = -fast -Bstatic -Mbounds
#FFLAGS = -fast -Bstatic
#FFLAGS = -O3 -Mbounds -Bstatic
#FFLAGS = -O3 -Bstatic
#FFLAGS = -fast
#FFLAGS = -O2
#FFLAGS = -O2 -pg
##ceFFLAGS = -O3 -fno-automatic
FFLAGS = -fno-automatic -g
#FFLAGS = -CB
#FFLAGS = -C
#FFLAGS = -ftrapuv -check all -check noarg_temp_created -debug all -Winline -vec-report0 -WB -g -traceback -warn all
#FFLAGS = -Mbounds
#FFLAGS = -BC -check all
#FFLAGS = -fbounds-check [not ifort]
#FC = g77
#FC = ifort
#FC = gfortran

# AMUSE STUFF======
MPIFC ?= mpif90
FC      = $(MPIFC)
#FFLAGS += -Wall -g
# =================
#FC = pgf77
#FC = pgf90
#FC = g95

#==========================
#c-stuff
CC=/usr/bin/gcc
#CC=icc
HEADERS=.
LIBS=rand.o
#PROGS=mufu plumix
#CFLAGS=-I$(HEADERS) -O3 -Wall -ffast-math
CFLAGS=-I$(HEADERS) -O3 
LDFLAGS=-lm $(MUSE_LD_FLAGS)
#===========================

# AMUSE
#AR = ar crs
AR = ar ruv
RANLIB = ranlib
RM = rm

#===============

SOURCE = checke.f coefpot.f collabr.f  core.f data.f \
denbin.f density.f energy.f escape.f formb3.f input.f intb3b3.f \
intb3b3a.f intb3b3b.f intb3f.f intcol.f interface_bse.f jarrod.f \
kick.f king.f lagrad.f mloss.f  mloss_single.f montcarl.f mydump.f \
newpos.f newvel.f output.f potent.f ran2.f relaxall.f relax.f \
relvel.f runtime.f scale0.f sort2.f sort3.f start.f timepot.f \
velocity.f zero.f zone.f jh_routines.f poly.f 

SOURCE_AMUSE = $(SOURCE:montcarl.f=amuse_montecarl.f)
SOURCE_AMUSE += amuse_output.f
SOURCE_AMUSE += amuse_interface.f
SOURCE_AMUSE += amuse_interface2.f

OBJECTS = $(SOURCE:.f=.o) plumix.o rand.o
OBJECTS_AMUSE = $(SOURCE_AMUSE:.f=.o) plumix.o rand.o

amuse: $(OBJECTS_AMUSE) params.h
	$(RM) -f libmmc.a
	$(AR) libmmc.a $(OBJECTS_AMUSE)
	$(RANLIB) libmmc.a

remark: 
	@echo $(SOURCE_AMUSE)

mont-car: $(OBJECTS) params.h
	$(FC) $(FFLAGS) $(OBJECTS) -o mont-car

standalone: standalone.o interface_bse.o jarrod.o
	$(FC) $(FFLAGS) standalone.o interface_bse.o jarrod.o -o standalone

.f.o:
	$(FC) -c $(FFLAGS)  $< -o $@

clean:
	rm -f *.o *.a

tar:checke.f coefpot.f collabr.f  core.f data.f \
	denbin.f density.f energy.f escape.f formb3.f input.f intb3b3.f \
	intb3b3a.f intb3b3b.f intb3f.f intcol.f interface_bse.f jarrod.f \
	kick.f king.f lagrad.f mloss.f  mloss_single.f montcarl.f mydump.f \
	newpos.f newvel.f output.f potent.f ran2.f relaxall.f relax.f \
	relvel.f runtime.f scale0.f sort2.f sort3.f start.f timepot.f \
	velocity.f zero.f zone.f  mont.run common.h  const_bse.h  interface_bse.h \
  params.h  zdata.h Makefile standalone.f jh_routines.f Kurucz.dat wdhyd.dat
	tar cvf mc.tar checke.f coefpot.f collabr.f  core.f data.f \
denbin.f density.f energy.f escape.f formb3.f input.f intb3b3.f \
intb3b3a.f intb3b3b.f intb3f.f intcol.f interface_bse.f* jarrod.f kick.f \
king.f lagrad.f mloss.f  mloss_single.f montcarl.f mydump.f \
newpos.f newvel.f output.f* potent.f ran2.f relaxall.f relax.f \
relvel.f runtime.f scale0.f sort2.f sort3.f start.f timepot.f \
velocity.f zero.f zone.f  mont.run common.h  const_bse.h  interface_bse.h \
params.h  zdata.h Makefile standalone.f jh_routines.f Kurucz.dat wdhyd.dat poly.f

