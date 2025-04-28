CODELIB = libbhtree.a

CXXFLAGS += -DTREE -DNOGRAPHICS 

OBJ = BHTC.o nbody.o BHtree.o gravity.o pgetopt.o second.o

CODEOBJS = BHtree.o gravity.o pgetopt.o second.o nbody.o

AR = ar ruv
RANLIB = ranlib
RM = rm

all: BHTC $(CODELIB)

clean:
	rm -f *.o BHTC *~ $(CODELIB)

BHTC: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -lm -o BHTC $(LDLIBS)

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $?

.C.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $?

test: compile

compile: BHTC

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

