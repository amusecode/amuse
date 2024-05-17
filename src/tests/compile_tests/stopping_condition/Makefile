all: c_worker c_mpi_ext_worker f_worker f_ext_worker f_mpi_ext_worker

# TestInterface
c_code.o: stopcond_test.c
	$(MPICC) -c -o $@ $(STOPCOND_CFLAGS) $(CFLAGS) -I. $<

c_interface.h: interface.py
	amusifier --mode=mpi --type=h -o $@ $< ForTestingInterface

c_interface.c: interface.py
	amusifier --mode=mpi --type=c -o $@ $< ForTestingInterface

c_interface.o: c_interface.c c_interface.h
	$(MPICXX) -c -o $@ $(STOPCOND_CFLAGS) $(CFLAGS) -I. $<

c_worker: c_code.o c_interface.o
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(STOPCOND_LIBS) $(LIBS)


# TestInterfaceMP
c_mpi_ext_code.o: stopcond_test.c
	$(MPICC) -c -o $@ $(STOPCONDMPI_CFLAGS) $(CFLAGS) -I. $<

c_ext_interface.h: interface.py
	amusifier --mode=mpi --type=h -o $@ $< ForTestingInterfaceFortranModule

c_ext_interface.c: interface.py
	amusifier --mode=mpi --type=c -o $@ $< ForTestingInterfaceFortranModule

c_mpi_ext_interface.o: c_ext_interface.c c_ext_interface.h
	$(MPICXX) -c -o $@ $(STOPCONDMPI_CFLAGS) $(CFLAGS) -I. $<

c_mpi_ext_worker: c_mpi_ext_code.o c_mpi_ext_interface.o
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(STOPCONDMPI_LIBS) $(LIBS)


# TestInterfaceFortran
f_code.o: stopcond_test.f
	$(MPIFC) -c -o $@ $(STOPCOND_CFLAGS) $(FCFLAGS) $<

f_interface.f90: interface.py
	amusifier --mode=mpi --type=f90 -o $@ $< ForTestingInterface

f_interface.o: f_interface.f90
	$(MPIFC) -c -o $@ $(STOPCOND_CFLAGS) $(FORSOCKETS_CFLAGS) $(FCFLAGS) $<

f_worker: f_code.o f_interface.o
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(STOPCOND_LIBS) $(FORSOCKETS_LIBS) $(LIBS)


# TestInterfaceFortranModule
f_ext_interface.f90: interface.py
	amusifier --mode=mpi --type=f90 -o $@ $< ForTestingInterfaceFortranModule

amuseinterface.o amuseinterface.mod &: amuseinterface.f90
	$(MPIFC) -c -o amuseinterface.o $(STOPCOND_CFLAGS) $(FCFLAGS) $<

f_ext_interface.o: f_ext_interface.f90 amuseinterface.mod
	$(MPIFC) -c -o $@ $(STOPCOND_CFLAGS) $(FORSOCKETS_CFLAGS) $(FCFLAGS) $<

f_ext_worker: f_code.o f_ext_interface.o amuseinterface.o
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(STOPCOND_LIBS) $(FORSOCKETS_LIBS) $(LIBS)


# TestInterfaceFortranModuleMultiprocess
f_mpi_ext_interface.o: f_ext_interface.f90 amuseinterface.mod
	$(MPIFC) -c -o $@ $(STOPCONDMPI_CFLAGS) $(FORSOCKETS_CFLAGS) $(FCFLAGS) $<

f_mpi_ext_worker: f_code.o f_mpi_ext_interface.o amuseinterface.o
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(STOPCONDMPI_LIBS) $(FORSOCKETS_LIBS) $(LIBS)


.PHONY: clean
clean:
	rm -f c_*.h c_*.c f_*.f90 *.mod *.o *_worker

.PHONY: distclean
distclean: clean

