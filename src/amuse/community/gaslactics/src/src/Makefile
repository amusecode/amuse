#Makefile for particle generators of the combined df - June 8/94 - JJD
#
.f.o:
	$(FC) $(FFLAGS) -c $*.f

FLIBS = -L/usr/X11/lib/ -lX11 -L${PGPLOT_DIR} -lpgplot 


ALLFILES= dbh genhalo genbulge gendisk diskdf  \
     getfreqs vcirc plothalo plotforce gengas toascii

#FFLAGS = -O -extend_source
#CFLAGS = -O2 -DRINGASCII -DASCII -DTP
#LIBS = -lm -lots -lfor


all: $(ALLFILES)

%.o: %.f90
	$(FC) $(FFLAGS2) -g -c -o $@ $< 

%.o: %.f
	$(FC) $(FFLAGS) -g -c -o $@ $< 

dbh: dbh.o dbhplot.o polardens.o bulgepotential.o totdens.o halopotential.o \
	pot.o thickdisk2.o diskdens.o dens.o appdiskpot.o plgndr1.o bulgedenspsi.o densrpsi.o \
	polarbulgedens.o polarhalodens.o appdiskdens.o halodens.o erfcen.o \
	modstamp.o fixedhalopot.o sphericaldf_nfw.o bessj01.o  \
	erfcc.o bulgedens.o fixedbulgepot.o sphericaldf_hernquist.o alldens.o	
	$(FC) $(FFLAGS)  -o $@ $^ $(FLIBS) $(LIBS)

genhalo: genhalo.o readmassrad.o query.o ran1.o readharmfile.o pot.o \
	halodens.o appdiskpot.o plgndr1.o densrpsi.o erfcen.o fixedhalopot.o thickdisk2.o \
	bessj01.o erfcc.o sphericaldf_nfw.o diskdens.o \
	fixedbulgepot.o sphericaldf_hernquist.o alldens.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

genbulge: genbulge.o readmassrad.o query.o ran1.o readharmfile.o pot.o \
	bulgedens.o appdiskpot.o plgndr1.o bulgedenspsi.o erfcen.o \
	fixedhalopot.o thickdisk2.o bessj01.o erfcc.o sphericaldf_nfw.o \
	bulgedens.o fixedbulgepot.o sphericaldf_hernquist.o diskdens.o alldens.o
	$(CC) $(CFLAGS)  -o $@ $^ $(LIBS)

gendisk: gendisk.o ran1.o query.o readdiskdf.o golden.o simpson.o diskdf5ez.o \
	 diskdensf.o readharmfile.o sigr2.o sigz2.o omekap.o diskdens.o \
	 splined.o splintd.o diskdf3ez.o diskdenspsi.o pot.o fnamidden.o \
	 appdiskpot.o plgndr1.o rcirc.o fixedhalopot.o  \
	 thickdisk2.o bessj01.o erfcc.o sphericaldf_nfw.o \
	 fixedbulgepot.o sphericaldf_hernquist.o alldens.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

gengas:  readharmfile.o diskinfo.o thickdisk2.o bessj01.o erfcc.o gasdisk3.o \
	diskdens.o sphericaldf_nfw.o bessj01.o erfcc.o pot.o force.o \
	fixedhalopot.o appdiskpot.o plgndr1.o appdiskforce.o \
	fixedbulgepot.o sphericaldf_hernquist.o alldens.o
	$(FC) $(FFLAGS)  -o $@ $^ $(LIBS)

diskdf: diskdf.o diskdf5intez.o diskdensf.o splintd.o splined.o readharmfile.o \
	sigr2.o sigz2.o omekap.o pot.o diskdf3intez.o appdiskpot.o plgndr1.o \
	fnamidden.o rcirc.o diskdens.o modstamp.o fixedhalopot.o \
	thickdisk2.o bessj01.o erfcc.o sphericaldf_nfw.o \
	fixedbulgepot.o sphericaldf_hernquist.o alldens.o
	$(FC) $(FFLAGS)  $(FLIBS) $(LIBS) -lstdc++ -o $@ $^
	
getfreqs: newfreqs.o pot.o appdiskpot.o plgndr1.o erfcen.o thickdisk2.o bessj01.o \
	erfcc.o fixedhalopot.o\
	fixedbulgepot.o
	$(FC) $(FFLAGS)  $(FLIBS) $(LIBS) -o $@ $^

#tobinary: tobinary.o
#	$(FC) $(CFLAGS) $(LIBS) -o $@ $^

vcirc: vcirc.o query.o readharmfile.o force.o appdiskforce.o bessj01.o erfcc.o\
	appdiskpot.o plgndr1.o sphericaldf_nfw.o diskdens.o fixedhalopot.o \
	thickdisk2.o pot.o fixedbulgepot.o sphericaldf_hernquist.o alldens.o
	$(CC) $(CFLAGS) $(LIBS) -o $@ $^
	
plothalo: plothalo.o contourden.o readdiskdf.o diskdens.o halodens.o bulgedens.o \
	bulgedenspsi.o diskdensf.o readharmfile.o sigr2.o densrpsi.o omekap.o pot.o \
	splined.o appdiskpot.o splintd.o diskdenspsi.o erfcen.o plgndr1.o  modstamp.o \
	thickdisk2.o sphericaldf_nfw.o fixedhalopot.o bessj01.o erfcc.o \
	bulgedens.o fixedbulgepot.o sphericaldf_hernquist.o alldens.o
	$(FC) $(FFLAGS2) -o $@ $^ $(FLIBS) $(LIBS)

plotforce: plotforce.o readharmfile.o readdiskdf.o force.o appdiskforce.o appdiskpot.o \
	diskdensf.o sigr2.o omekap.o plgndr1.o diskdens.o splined.o splintd.o diskdenspsi.o \
	pot.o modstamp.o fixedhalopot.o sphericaldf_nfw.o thickdisk2.o bessj01.o \
	fixedbulgepot.o sphericaldf_hernquist.o erfcc.o alldens.o
	$(FC) $(FFLAGS) $(FLIBS) $(LIBS)  -o $@ $^

testgas: testgas.o ran1.o readharmfile.o pot.o \
	halodens.o appdiskpot.o plgndr1.o densrpsi.o erfcen.o fixedhalopot.o thickdisk2.o \
	bessj01.o erfcc.o sphericaldf_nfw.o diskdens.o \
	fixedbulgepot.o sphericaldf_hernquist.o alldens.o
	$(FC) $(CFLAGS) -o $@ $^ $(LIBS)

toascii: toascii.o
	$(CC) $(CFLAGS) toascii.o -o toascii $(LIBS)


clean:
	rm -f *.o *.mod *.bck

purge: clean
	rm -f $(ALLFILES)

install: all
	cp -f $(ALLFILES) ../bin

#flatdiskpot.o: flatdiskpot.cc
#	g++ -I/software/local/include
