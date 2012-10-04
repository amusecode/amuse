AMUSE_DIR?=../../..
include ${AMUSE_DIR}/config.mk


all:stopcondf.o
	
%.o: ../%.F90
	$(MPIFC) $(FCFLAGS) -DMPILIB -c -o $@ $<


clean:
	rm -Rf *.o *.lo *.a *.mo *.mod
