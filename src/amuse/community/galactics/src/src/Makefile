#Makefile for particle generators of the combined df - June 8/94 - JJD
#
.f.o:
	$(F77) $(FFLAGS) -c $*.f

F77?=gfortran
CC?=gcc
#F77=ifort
#CC=icc
#CFLAGS = -O3
#CFLAGS = -O -DRINGASCII -DASCII
LIBS = -lm
#These flags for SunOS:
#FFLAGS =  -e -C -O 
#FLIBS = -L/usr/X11/lib/ -lX11 -lpgplot 

#Replace FFLAGS and FLIBS above with the following for linux:
#FFLAGS = -ffixed-line-length-0 -O -fno-backslash
FFLAGS += -O3 -132 -nofor-main
#FLIBS  = -L/usr/X11R6/lib64/ -lX11
#FLIBS  = -L/usr/X11R6/lib/ -lX11

ICPROGS= dbh genhalo genbulge getfreqs gendisk diskdf toascii

#FFLAGS = -O -extend_source
#CFLAGS = -O2 -DRINGASCII -DASCII -DTP
#LIBS = -lm -lots -lfor

all: ics parallel

ics: $(ICPROGS)

dbh: dbh.o polardens.o bulgepotential.o totdens.o halopotential.o pot.o diskdens.o dens.o appdiskpot.o plgndr1.o bulgedenspsi.o halodenspsi.o gendenspsihalo.o gendenspsibulge.o polarbulgedens.o polarhalodens.o appdiskdens.o halodens.o dfhalo.o dfbulge.o erfcen.o modstamp.o force.o appdiskforce.o gendf.o getpsi.o dpolardens.o diskpotentialestimate.o halopotentialestimate.o nfwprofiles.o sersicprofiles.o
	$(F77) $(FFLAGS) $(FLIBS) dbh.o polardens.o bulgepotential.o totdens.o halopotential.o pot.o diskdens.o dens.o appdiskpot.o plgndr1.o bulgedenspsi.o halodenspsi.o gendenspsihalo.o gendenspsibulge.o polarbulgedens.o polarhalodens.o appdiskdens.o halodens.o dfhalo.o dfbulge.o erfcen.o modstamp.o force.o appdiskforce.o gendf.o getpsi.o dpolardens.o diskpotentialestimate.o halopotentialestimate.o nfwprofiles.o sersicprofiles.o -o dbh

omegakappa: omegakappa.o polardens.o bulgepotential.o totdens.o halopotential.o pot.o diskdens.o dens.o appdiskpot.o plgndr1.o bulgedenspsi.o halodenspsi.o gendenspsihalo.o gendenspsibulge.o polarbulgedens.o polarhalodens.o appdiskdens.o halodens.o dfhalo.o dfbulge.o  erfcen.o modstamp.o force.o appdiskforce.o splined.o splintd.o
	$(F77) $(FFLAGS) $(FLIBS) omegakappa.o polardens.o bulgepotential.o totdens.o halopotential.o pot.o diskdens.o dens.o appdiskpot.o plgndr1.o bulgedenspsi.o halodenspsi.o gendenspsihalo.o gendenspsibulge.o polarbulgedens.o polarhalodens.o appdiskdens.o halodens.o dfhalo.o dfbulge.o  erfcen.o modstamp.o force.o appdiskforce.o splined.o splintd.o -o omegakappa 

gendisk: gendisk.o ran1.o query.o \
	appdiskforce.o force.o velocityfactors.o readdiskdf.o golden.o simpson.o diskdf5ez.o diskdensf.o readharmfile.o sigr2.o sigz2.o omekap.o diskdens.o splined.o splintd.o diskdf3ez.o diskdenspsi.o pot.o fnamidden.o appdiskpot.o plgndr1.o rcirc.o getpsi.o sersicprofiles.o
	$(F77) $(FFLAGS) gendisk.o ran1.o query.o appdiskforce.o force.o velocityfactors.o diskdf5ez.o diskdensf.o readharmfile.o sigr2.o sigz2.o omekap.o diskdens.o splined.o splintd.o diskdf3ez.o diskdenspsi.o pot.o fnamidden.o appdiskpot.o plgndr1.o rcirc.o getpsi.o readdiskdf.o golden.o simpson.o sersicprofiles.o -o gendisk

genbulge: genbulge.o readmassrad.o query.o ran1.o readharmfile.o pot.o dfbulge.o dfhalo.o  appdiskpot.o plgndr1.o bulgedenspsi.o readdenspsi.o erfcen.o getpsi.o sersicprofiles.o force.o appdiskforce.o 
	$(F77) $(FFLAGS) genbulge.o readmassrad.o query.o ran1.o readharmfile.o pot.o dfbulge.o dfhalo.o  appdiskpot.o plgndr1.o bulgedenspsi.o readdenspsi.o erfcen.o getpsi.o sersicprofiles.o force.o appdiskforce.o -o genbulge

genhalo: genhalo.o readmassrad.o dfcorrection.o query.o ran1.o readharmfile.o pot.o halodens.o dfhalo.o dfbulge.o  appdiskpot.o plgndr1.o halodenspsi.o readdenspsi.o erfcen.o getpsi.o sersicprofiles.o force.o appdiskforce.o 
	$(F77) $(FFLAGS) genhalo.o readmassrad.o dfcorrection.o query.o ran1.o readharmfile.o appdiskpot.o plgndr1.o halodenspsi.o readdenspsi.o pot.o halodens.o dfhalo.o dfbulge.o  erfcen.o getpsi.o sersicprofiles.o force.o appdiskforce.o -o genhalo


diskdf: diskdf.o diskdf5intez.o diskdensf.o splintd.o splined.o readharmfile.o sigr2.o sigz2.o omekap.o pot.o diskdf3intez.o appdiskpot.o plgndr1.o fnamidden.o rcirc.o diskdens.o modstamp.o
	$(F77) $(FFLAGS)  $(FLIBS) diskdf.o diskdf5intez.o diskdf3intez.o appdiskpot.o plgndr1.o diskdensf.o splintd.o splined.o readharmfile.o sigr2.o sigz2.o fnamidden.o rcirc.o diskdens.o omekap.o pot.o modstamp.o  -o diskdf

getfreqs: getfreqs.o pot.o appdiskpot.o plgndr1.o erfcen.o
	$(F77) $(FFLAGS)  $(FLIBS) getfreqs.o pot.o appdiskpot.o plgndr1.o erfcen.o -o getfreqs

toascii: toascii.o
	$(CC) $(CFLAGS) toascii.o -o toascii $(LIBS)

#Parallel versions

MPICC?=mpicc
MPIFC?=mpif90
MPIF77=$(MPIFC)

PARPROGS=pgendisk pgenhalo pgenbulge
parallel: $(PARPROGS)

pgendisk.o: pgendisk.c
	$(MPICC) $(CFLAGS) -c pgendisk.c

pgenbulge.o: pgenbulge.c
	$(MPICC) $(CFLAGS) -c pgenbulge.c

pgenhalo.o: pgenhalo.c
	$(MPICC) $(CFLAGS) -c pgenhalo.c


pgendisk: pgendisk.o ran1.o query.o appdiskforce.o force.o velocityfactors.o \
readdiskdf.o golden.o simpson.o diskdf5ez.o diskdensf.o readharmfile.o \
sigr2.o sigz2.o omekap.o diskdens.o splined.o splintd.o diskdf3ez.o \
diskdenspsi.o pot.o fnamidden.o appdiskpot.o plgndr1.o rcirc.o getpsi.o \
sersicprofiles.o
	$(MPIF77) $(FFLAGS) pgendisk.o ran1.o query.o appdiskforce.o force.o \
velocityfactors.o diskdf5ez.o diskdensf.o readharmfile.o sigr2.o sigz2.o \
omekap.o diskdens.o splined.o splintd.o diskdf3ez.o diskdenspsi.o pot.o \
fnamidden.o appdiskpot.o plgndr1.o rcirc.o getpsi.o readdiskdf.o golden.o \
simpson.o sersicprofiles.o -o pgendisk

pgenbulge: pgenbulge.o readmassrad.o query.o ran1.o readharmfile.o pot.o dfbulge.o dfhalo.o  appdiskpot.o plgndr1.o bulgedenspsi.o readdenspsi.o erfcen.o getpsi.o sersicprofiles.o force.o appdiskforce.o 
	$(MPIF77) $(FFLAGS) pgenbulge.o readmassrad.o query.o ran1.o \
readharmfile.o pot.o dfbulge.o dfhalo.o  appdiskpot.o plgndr1.o \
bulgedenspsi.o readdenspsi.o erfcen.o getpsi.o sersicprofiles.o force.o \
appdiskforce.o -o pgenbulge

pgenhalo: pgenhalo.o readmassrad.o dfcorrection.o query.o ran1.o readharmfile.o pot.o halodens.o dfhalo.o dfbulge.o  appdiskpot.o plgndr1.o halodenspsi.o readdenspsi.o erfcen.o getpsi.o sersicprofiles.o force.o appdiskforce.o 
	$(MPIF77) $(FFLAGS) pgenhalo.o readmassrad.o dfcorrection.o query.o \
ran1.o readharmfile.o appdiskpot.o plgndr1.o halodenspsi.o readdenspsi.o \
pot.o halodens.o dfhalo.o dfbulge.o  erfcen.o getpsi.o sersicprofiles.o \
force.o appdiskforce.o -o pgenhalo


clean:
	rm -f *.o $(ICPROGS) $(PARPROGS)
	(cd ../bin; rm -f *.o $(ICPROGS) $(PARPROGS))

install:
	cp -f $(ICPROGS) $(PARPROGS) ../bin


readmassrad.o: readmassrad.c
simpson.o: simpson.c
ran1.o: ran1.c
query.o: query.c
golden.o: golden.c
tobinary.o: tobinary.c
centre1.o: centre1.c
