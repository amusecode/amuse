CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS)

CODELIB = libtsunami.a

CODEOBJS = tsunami.o


all: $(CODELIB)


clean:
	rm -f *.o *.a

distclean: clean

$(CODELIB): $(CODEOBJS)
	rm -f $@
	$(AR) -ruv $@ $(CODEOBJS)
	$(RANLIB) $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<

