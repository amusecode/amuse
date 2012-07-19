CXXFLAGS  += -I../include -I../include/star \
	     -D_SRC_='"no_source_available"' \
	     -DHAVE_CONFIG_H  -DTOOLBOX 
LDLIBS	  += \
	     -L. -lrdc  \
	     -L../dstar -ldstar  \
	     -L../sstar -lsstar \
	     -L../node -lnode \
	     -L../node/dyn -ldyn \
	     -L../std -lstd  \
	     -lm

CFLAGS    += -O

EXE	= rdc_SeBa rs_birthrate rs_countsn rs_findtype rs_population rs_reorder rs_scenarios rs_snapshot
DIRS	= starclass 

all:	lib $(EXE)

Makefile.inc:	Makefile.inc.conf
	@sed s#__BASE__#$(PWD)#g Makefile.inc.conf >Makefile.inc

lib:	Makefile.inc
	@for d in $(DIRS) ; do echo '\nmake' $@ in $$d; make -C $$d $@ ; done

clean:	Makefile.inc
	@for d in $(DIRS) ; do echo '\nmake' $@ in $$d; make -C $$d $@ ; done
	/bin/rm -f *~ *.a $(EXE) Makefile.inc

