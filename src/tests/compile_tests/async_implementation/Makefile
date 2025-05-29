all: async_worker

# TestInterface
c_code.o: ../c_implementation/c_implementation_test.c
	$(MPICC) -c -o $@ $(CFLAGS) -I. $<

async_code.o: async_test.c ../c_implementation/c_implementation_test.c
	$(MPICC) -c -o $@ $(CFLAGS) -I../c_implementation $<

c_interface.h: interface.py
	amusifier --mode=mpi --type=h -o $@ $< ForTestingInterface

c_interface.c: interface.py
	amusifier --mode=mpi --type=c -o $@ $< ForTestingInterface

c_interface.o: c_interface.c c_interface.h
	$(MPICXX) -c -o $@ $(CFLAGS) -I. $<

async_worker: async_code.o c_code.o c_interface.o
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LIBS)


.PHONY: clean
clean:
	rm -f c_interface.h c_interface.c *.o *_worker

.PHONY: distclean
distclean: clean

