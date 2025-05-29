all: c_worker

# TestInterface
cpp_code.o: speed_test.cpp
	$(MPICXX) -c -o $@ $(CFLAGS) -I. $<

c_interface.h: interface.py
	amusifier --mode=mpi --type=h -o $@ $< TestCode

c_interface.c: interface.py
	amusifier --mode=mpi --type=c -o $@ $< TestCode

c_interface.o: c_interface.c c_interface.h
	$(MPICXX) -c -o $@ $(CFLAGS) -I. $<

c_worker: cpp_code.o c_interface.o
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LIBS)


.PHONY: clean
clean:
	rm -f c_interface.h c_interface.c *.o *_worker

.PHONY: distclean
distclean: clean

