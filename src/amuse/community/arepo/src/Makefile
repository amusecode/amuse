CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CODELIB = libarepo.a

CODEOBJS = test.o

AR = ar ruv
RANLIB = ranlib
RM = rm

all: $(CODELIB) 


clean:
	$(RM) -f *.o *.a

distclean: clean

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
