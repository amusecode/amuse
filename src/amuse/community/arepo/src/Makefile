ifeq ($(origin AMUSE_DIR), undefined)
	AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

GSL_FLAGS  = -I/opt/Homebrew/include  # TODO: Need to make GSL_FLAGS generalisable.

CFLAGS   += -Wall -g $(GSL_FLAGS)
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CODELIB = libarepo.a

SUBDIR = add_backgroundgrid cooling debug_md5 domain fof gitversion gravity \
	hydro init io main mesh mpi_utils ngbtree star_formation subfind \
	time_integration utils
SRCS = $(foreach fd, $(SUBDIR), $(wildcard $(fd)/*.c))

CODEOBJS = test.o $(SRCS:c=o)

AR = ar ruv
RANLIB = ranlib
RM = rm

all: $(CODELIB)

clean:
	$(RM) -f *.o *.a
	$(RM) $(SRCS:c=o)

distclean: clean

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
