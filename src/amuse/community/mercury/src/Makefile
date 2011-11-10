# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../../..
-include ${AMUSE_DIR}/config.mk

#IGNORE:
CMPLR ?= $(FC)
FFLAGS += -O 

SRCE = mercury6_2.for mercury_main.for 
DFLE = big.in  message.in  small.in element.in param.in swift.inc \
	close.in files.in mercury.inc mercury6.man  
TARFILES = $(SRCE) $(DFLE) README.txt Makefile

EXECUTABLES = mercury
OBJT = $(SRCE:.for=.o)

SCINC ?= -I$(AMUSE_DIR)/lib/stopcond
SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi

%.o: %.for Makefile
	$(CMPLR) $(FFLAGS) $(SCINC) -c -o $@ $< 

%.o: %.f90 Makefile
	$(CMPLR) $(FFLAGS) $(SCINC) -c -o $@ $< 

all: $(EXECUTABLES)

amuse_interface: libmercury.a

libmercury.a: $(OBJT) hash.o amuse_helpers.o
	ar crs $@ $^

mercury: $(OBJT) 
	$(CMPLR) $(FFLAGS) $(OBJT) -o $@ $(SCLIBS)

#mercury6_2.o:
#	$(CMPLR) -c $(FFLAGS) $(SRCE) -c $< -DTOOLBOX

amuse: mercury6_2.o

clean:
	rm -rf *.o *~ *.mod *.bck $(OUTPUTFILES) $(EXECUTABLES)

tar:
	make clean; tar cvf mercury.tar $(TARFILES); gzip mercury.tar

