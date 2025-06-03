all: f_worker

# TestInterfaceFortran
f_code.o: f_implementation_test.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

f_interface.f90: interface.py
	amusifier --mode=mpi --type=f90 -o $@ $< ForTestingInterface

f_interface.o: f_interface.f90
	$(MPIFC) -c -o $@ $(FORSOCKETS_CFLAGS) $(FCFLAGS) $<

f_worker: f_code.o f_interface.o
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(FORSOCKETS_LIBS) $(LIBS)


.PHONY: clean
clean:
	rm -f f_interface.f90 *.mod *.o *_worker

.PHONY: distclean
distclean: clean

