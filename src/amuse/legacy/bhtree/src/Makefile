
CODELIB = libnearestneighbor.a

CXXFLAGS += -Wall -DTREE -DNOGRAPHICS 

OBJ = BHTC.o nbody.o BHtree.o gravity.o pgetopt.o second.o

CODEOBJS = BHtree.o gravity.o pgetopt.o second.o nbody.o

AR = ar ruv
RANLIB = ranlib
RM = rm

all: BHTC $(CODELIB)

clean:
	rm -f *.o BHTC *~ $(CODELIB)

BHTC: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -lm -o BHTC

.c.o:
	$(CC) $(CFLAGS) -c $?

.C.o:
	$(CXX) $(CXXFLAGS) -c $?

test: compile

compile: BHTC

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

