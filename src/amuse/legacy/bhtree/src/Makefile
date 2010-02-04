

CXXFLAGS += -Wall -DTREE -DNOGRAPHICS 

OBJ = BHTC.o nbody.o BHtree.o gravity.o pgetopt.o second.o

all: BHTC

clean:
	rm -f *.o BHTC *~

BHTC: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -lm -o BHTC

.c.o:
	$(CC) $(CFLAGS) -c $?

.C.o:
	$(CXX) $(CXXFLAGS) -c $?

test: compile

compile: BHTC
