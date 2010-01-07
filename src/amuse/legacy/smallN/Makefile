ifdef STARLAB_INSTALL_PATH
    export STARLAB_INSTALL_PATH
    TARGS = .compile
else
    TARGS = .warning
endif

OBJS	=	amuse_interface.o integrate.o analyze.o
EXE	=	smallN2 triple

all:
	@echo "Nothing build by default"


CXXFLAGS +=	-I$(STARLAB_INSTALL_PATH)/include/starlab -DHAVE_CONFIG_H -O -g
LDLIBS	 +=	-L$(STARLAB_INSTALL_PATH)/lib/starlab \
			-lhdyn -l_dyn_ -ldstar -lhdyn -lstubs \
			-ldyn -lsstar -lnode -lpthread -lstd -lm

LIBS	  =	$(STARLAB_INSTALL_PATH)/lib/starlab/libhdyn.a \
		$(STARLAB_INSTALL_PATH)/lib/starlab/lib_dyn_.a \
		$(STARLAB_INSTALL_PATH)/lib/starlab/libdyn.a \
		$(STARLAB_INSTALL_PATH)/lib/starlab/libnode.a \
		$(STARLAB_INSTALL_PATH)/lib/starlab/libstd.a \
		$(STARLAB_INSTALL_PATH)/lib/starlab/libstubs.a \
		$(STARLAB_INSTALL_PATH)/lib/starlab/libdstar.a \
		$(STARLAB_INSTALL_PATH)/lib/starlab/libhdyn.a \
		$(STARLAB_INSTALL_PATH)/lib/starlab/libsstar.a

compile:	$(TARGS)

.compile: muse_worker 

.warning:
	@echo ""
	@echo STARLAB_INSTALL_PATH = $(STARLAB_INSTALL_PATH)
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "Starlab not found; skipping smallN build for now.  Please"
	@echo "install starlab and set STARLAB_INSTALL_PATH to the install path."
	@echo "To download starlab, see"
	@echo ""
	@echo "    http://www.ids.ias.edu/~starlab/install"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo ""

# Standalone test program:

smallN2:	smallN2.o integrate.o analyze.o
	$(CXX) $(CXXFLAGS) smallN2.o integrate.o analyze.o \
	       $(LDLIBS) -o $@

mytest0:	$(EXE)
	$(STARLAB_INSTALL_PATH)/bin/makeplummer -n 5 -i \
		| ./smallN2 -s 10 -t 2000

mytest1:	$(EXE)
	$(STARLAB_INSTALL_PATH)/bin/makeplummer -n 5 -i -s 1191011580 \
		| ./smallN2 -s 10 -t 2000

mytest2:	$(EXE)
	$(STARLAB_INSTALL_PATH)/bin/makeplummer -n 5 -i -s 1191033642 \
		| ./smallN2 -s 10 -t 2000

triple:	triple.o amuse_interface.o integrate.o analyze.o
	$(CXX) $(CXXFLAGS) $@.o amuse_interface.o integrate.o analyze.o \
	       $(LDLIBS) -o $@

cleanall: clean
	$(RM) muse_worker $(EXE) *~
	rm -Rf muse_worker.dSYM/

clean:
	rm -f *.so *.o *.pyc muse_worker.cc

muse_worker.cc: muse_dynamics_mpi.py
	../../../../bin/create_c_worker.py muse_dynamics_mpi.py SmallN > $@

muse_worker_starlab:        muse_worker.cc $(OBJS)
	mpicxx $(CXXFLAGS) $(LDLIBS) $@.cc $(OBJS) -o $@

.cc.o: $<
	g++ $(CXXFLAGS) -c -o $@ $<

.c.o: $<
	g++ $(CXXFLAGS) -c -o $@ $<

