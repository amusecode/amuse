CXXFLAGS  += -I$(STARLAB_INSTALL_PATH)/include \
	     -I$(STARLAB_INSTALL_PATH)/include/star
LDLIBS	  += -L$(STARLAB_INSTALL_PATH)/usr/lib/starlab \
	     -lsstar -ldsatr -ldyn -lnode -lstd -lm

CFLAGS    +=	-O

CXXFLAGS += -D_SRC_='"no_source_available"'