PS_PATH = ../../../src/
INC = -I$(PS_PATH)

#CXX = g++
CXX = mpicxx
CXXFLAGS = -std=c++17 -O3 -ffast-math -funroll-loops -Wall
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
# CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

# Simulation control macros 
CXXFLAGS += -DENABLE_VARIABLE_SMOOTHING_LENGTH
CXXFLAGS += -DUSE_ENTROPY
CXXFLAGS += -DUSE_BALSARA_SWITCH
#CXXFLAGS += -DUSE_PRESCR_OF_THOMAS_COUCHMAN_1992
CXXFLAGS += -DISOTHERMAL_EOS

CXXFLAGS += -DINITIAL_CONDITION=0
#CXXFLAGS += -DINITIAL_CONDITION=1
#CXXFLAGS += -DINITIAL_CONDITION=2
#CXXFLAGS += -DINITIAL_CONDITION=3

use_phantom_grape_x86 = yes

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

ifeq ($(use_phantom_grape_x86),yes)
PG_ROOT = $(PS_PATH)/phantom_grape_x86/G5/newton/libpg5
INC += -I$(PG_ROOT)
CXXFLAGS += -DENABLE_PHANTOM_GRAPE_X86
LIBS = -L$(PG_ROOT) -lpg5
PG_BUILD = cd $(PG_ROOT) && $(MAKE) distclean libpg5.a
PG_CLEAN = cd $(PG_ROOT) && $(MAKE) distclean
else
PG_BUILD =
PG_CLEAN = 
endif

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CPPHDRS = $(wildcard *.h)
PROGRAM = mizuki.out

.PHONY:	clean all

all: $(CPPOBJS) $(CPPHDRS)
	$(PG_BUILD)
	$(CXX) $(CXXFLAGS) $(CPPOBJS) -o $(PROGRAM) $(LIBS) $(INC)

%.o: %.cpp $(CPPHDRS)
	$(CXX) -c $< $(CXXFLAGS) $(INC)

clean:
	rm -f $(CPPOBJS)

distclean: clean
	$(PG_CLEAN)
	rm -f $(PROGRAM)
	rm -rf result

# fdps-autotest-run (DO NOT CHANGE THIS LINE)

