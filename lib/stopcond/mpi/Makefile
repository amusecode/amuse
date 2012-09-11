AMUSE_DIR?=../../..
include ${AMUSE_DIR}/config.mk


all:stopcondf.o
	
%.o: ../%.F90
	$(FC) $(FCFLAGS) -c -o $@ $<


clean:
	rm -Rf *.o *.lo *.a *.mo *.mod
