all: c_worker

# TestInterface
c_code.o: c_implementation_test.c
	$(MPICC) -c -o $@ $(CFLAGS) -I. $<

c_interface.h: interface.py
	amusifier --mode=mpi --type=h -o $@ $< ForTestingInterface

c_interface.c: interface.py
	amusifier --mode=mpi --type=c -o $@ $< ForTestingInterface

c_interface.o: c_interface.c c_interface.h
	$(MPICXX) -c -o $@ $(CFLAGS) -I. $<

c_worker: c_code.o c_interface.o
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LIBS)


.PHONY: clean
clean:
	rm -f c_interface.h c_interface.c *.o *_worker

.PHONY: distclean
distclean: clean

