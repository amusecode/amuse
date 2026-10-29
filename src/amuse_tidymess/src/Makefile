CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS)

RANLIB = ranlib

CODELIB := libtidymess.a

SRCDIR := tidymess/integrator/src
INCDIR := tidymess/integrator/include

SRC := $(wildcard $(SRCDIR)/*.cpp)
CODEOBJS := $(SRC:$(SRCDIR)/%.cpp=%.o)


all: $(CODELIB)


clean:
	rm -f *.o *.a

distclean: clean

$(CODELIB): $(CODEOBJS)
	rm -f $@
	$(AR) -ruv $@ $(CODEOBJS)
	$(RANLIB) $@

%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ -I$(INCDIR) $<
