ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

all:libg6.a

CFLAGS += -O2

libg6.a: g6lib.o
	ar -r  libg6.a g6lib.o
	ranlib libg6.a


clean:
	rm -Rf *.o *.lo *.a
