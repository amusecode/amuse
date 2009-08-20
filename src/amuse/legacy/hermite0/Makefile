
CFLAGS = -Wall -g -DTOOLBOX  $(MUSE_INCLUDE_DIR) -fPIC
CXXFLAGS = $(CFLAGS) 
LDFLAGS = -lm $(MUSE_LD_FLAGS)

PYTHONINC := $(shell python -c 'import distutils.sysconfig; \
                                print distutils.sysconfig.get_python_inc()')

UNAME = $(shell uname)

ifeq ($(UNAME),Darwin)
	SHAREDFLAGS = -bundle -undefined suppress -flat_namespace
else
	SHAREDFLAGS = -shared
endif

CXXFLAGS += -fPIC



OBJS = muse_dynamics.o

all:	 _muse_dynamics_interface.so muse_driver muse_worker

cleanall:	clean
	$(RM) muse_driver *~

# Standalone test program:

muse_driver:	muse_driver.cc $(OBJS)
	g++ $@.cc $(OBJS) -o $@

muse_worker.cc: muse_dynamics_mpi.py
	../../../../bin/create_c_worker.py muse_dynamics_mpi.py Hermite > $@

muse_worker:	muse_worker.cc $(OBJS)
	mpicxx $@.cc $(OBJS) -o $@

_%_interface.so: %_wrap.o $(OBJS)
	g++ $(SHAREDFLAGS) $^ $(LIBS) -o $@

compare:	all
	@echo comparing C++ and python drivers
	muse_driver < ../../muse/in.20 > C++.out
	python test.py > python.out
	@echo output files should be identical
	diff -bwi C++.out python.out
	@rm -i C++.out python.out



.cc.o: $<
	g++ $(CXXFLAGS) -c -o $@ $< 
.C.o:
	$(CC) $(CFLAGS) -c -o $@ $?

clean: standalones_clean swigclean

standalones_clean:
	rm -f *.so *.o *.pyc

%_wrap.cxx %.py:: %.i %.h
	swig -c++ -python $<

%_wrap.o: %_wrap.cxx
	g++ -c $(CXXFLAGS) $< -o $@ -I $(PYTHONINC)



.PHONY = swigclean

SWIGI=$(shell ls *.i)
SWIGPY=$(SWIGI:.i=_interface.py)
SWIGPYC=$(SWIGI:.i=_interface.pyc)
SWIGO=$(SWIGI:.i=.o)
SWIGSO=$(patsubst %.i,_%_interface.so,$(SWIGI))
SWIGCXX=$(SWIGI:.i=_wrap.cxx)

swigclean:
	rm -f $(SWIGPY) $(SWIGPYC) $(SWIGO) $(SWIGSO) $(SWIGCXX)



