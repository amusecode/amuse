.f.o:
	$(F77) $(FFLAGS) -c $*.f

#F77=ifort
#CC=icc
F77=gfortran
CC=gcc
CFLAGS = -O3
LIBS = -lm
RIGIDOBJS=rigidpot.o appdiskforce.o appdiskpot.o force.o plgndr1.o readharmfile.o
RECIPEOBJS=bsstep.o mmid.o rzextr.o nrutil.o spline.o splint.o

FFLAGS = -ffixed-line-length-0 -O -fno-backslash 
#FFLAGS = -extend_source -nofor_main -O

ALLFILES= testforce testp checkvirial vcirc energysort

all: $(ALLFILES) librigidpot.a

checkvirial: checkvirial.o librigidpot.a
	$(F77) $(FFLAGS) $(FLIBS) checkvirial.o librigidpot.a -o checkvirial

energysort: energysort.o librigidpot.a
	$(F77) $(FFLAGS) $(FLIBS) energysort.o librigidpot.a -o energysort

testforce: testforce.o librigidpot.a
	$(F77) $(FFLAGS) $(FLIBS) testforce.o librigidpot.a -o testforce

vcirc: vcirc.o librigidpot.a
	$(F77) $(FFLAGS) $(FLIBS) vcirc.o librigidpot.a -o vcirc

librigidpot.a: $(RIGIDOBJS)
	ar r librigidpot.a $(RIGIDOBJS)
	ranlib librigidpot.a

testp: testp.o
	$(CC) $(CFLAGS) testp.o -o testp -lm

clean:
	rm -f *.o $(ALLFILES) librigidpot.a
	(cd ../bin; rm -f $(ALLFILES))

install:
	cp energysort testforce checkvirial vcirc ../bin
