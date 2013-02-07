AMUSE_DIR?=../../../../..
-include ${AMUSE_DIR}/config.mk

CODELIB = libbhtree.a

CXXFLAGS += -DTREE -DNOGRAPHICS 

OBJ = BHTC.o nbody.o BHtree.o gravity.o pgetopt.o second.o

CODEOBJS = BHtree.o gravity.o pgetopt.o second.o nbody.o

AR = ar ruv
RANLIB = ranlib
RM = rm

SC_LIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcond
SC_FLAGS ?= -I$(AMUSE_DIR)/lib/stopcond

all: BHTC $(CODELIB)

clean:
	rm -f *.o BHTC *~ $(CODELIB)

BHTC: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -lm -o BHTC $(SC_LIBS)

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SC_FLAGS) -c $?

.C.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(SC_FLAGS) -c $?

test: compile

compile: BHTC

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

