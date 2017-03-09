ifndef AMUSE_DIR
    AMUSE_DIR=../../../..
endif

CODE_GENERATOR = $(AMUSE_DIR)/build.py

CXXFLAGS = -Wall -g -DTOOLBOX  -O3 $(MUSE_INCLUDE_DIR)
LDFLAGS = -lm $(MUSE_LD_FLAGS)

OBJS = interface.o src/evolve.o src/ODE_system.o src/types.o src/cvode/cvode.o src/cvode/cvode_dense.o src/cvode/cvode_direct.o src/cvode/cvode_io.o src/cvode/nvector_serial.o src/cvode/sundials_dense.o src/cvode/sundials_direct.o src/cvode/sundials_math.o src/cvode/sundials_nvector.o 


all: worker_code

distclean: clean
	$(RM) worker_code 

clean:
	rm -f *.so *.o *.pyc worker_code.cc src/*.o*

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py SecularMultipleInterface -o $@

worker_code: worker_code.cc $(OBJS)
	mpicxx $@.cc $(OBJS) -o $@

.cc.o: $<
	g++ $(CXXFLAGS) -c -o $@ $<

.c.o: $<
	g++ $(CXXFLAGS) -c -o $@ $<
