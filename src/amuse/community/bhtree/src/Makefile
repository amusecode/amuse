CODELIB = libbhtree.a

CXXFLAGS += -Wall -DTREE -DNOGRAPHICS 

OBJ = BHTC.o nbody.o BHtree.o gravity.o pgetopt.o second.o

CODEOBJS = BHtree.o gravity.o pgetopt.o second.o nbody.o

AR = ar ruv
RANLIB = ranlib
RM = rm

SCLIBS ?= -L../../../../../lib/stopcond -lstopcond
SCINC ?= -I../../../../../lib/stopcond

all: BHTC $(CODELIB)

clean:
	rm -f *.o BHTC *~ $(CODELIB)

BHTC: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -lm -o BHTC $(SCLIBS)

.c.o:
	$(CC) $(CFLAGS) $(SCINC) -c $?

.C.o:
	$(CXX) $(CXXFLAGS) $(SCINC) -c $?

test: compile

compile: BHTC

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

