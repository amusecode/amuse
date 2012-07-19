CXXFLAGS  += -I../include -I../include/star \
	     -D_SRC_='"no_source_available"' \
	     -DHAVE_CONFIG_H  -DTOOLBOX 
LDLIBS	  += \
	     -L. -ldstar  \
	     -L../sstar -lsstar \
	     -L../node -lnode \
	     -L../node/dyn -ldyn \
	     -L../std -lstd  \
	     -lm

CFLAGS    += -O

EXE	= SeBa
DIRS	= init io starclass stardyn util

all:	lib $(EXE)

Makefile.inc:	Makefile.inc.conf
	@sed s#__BASE__#$(PWD)#g Makefile.inc.conf >Makefile.inc

lib:	Makefile.inc
	@for d in $(DIRS) ; do echo '\nmake' $@ in $$d; make -C $$d $@ ; done

clean:	Makefile.inc
	@for d in $(DIRS) ; do echo '\nmake' $@ in $$d; make -C $$d $@ ; done
	/bin/rm -f *~ $(EXE) Makefile.inc

