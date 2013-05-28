CODELIB = libhop.a
CODEOBJS = hop.o kd.o smooth.o

all: $(CODELIB)

clean:
	rm -f *.o $(CODELIB)

$(CODELIB): $(CODEOBJS)
	ar ruv $@ $(CODEOBJS)
	ranlib $@

.c.o: $<
	$(CC) $(CFLAGS) -c -o $@ $< 
