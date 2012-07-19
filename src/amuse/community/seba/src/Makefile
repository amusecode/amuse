CXXFLAGS  += -I./include -I./include/star -DHAVE_CONFIG_H -DTOOLBOX
LDLIBS	  += -Lsstar -lsstar -Lnode -lnode -Lstd -lstd -lm
CFLAGS    += -O

EXE	= starev
#DIRS	= node std sstar dstar
DIRS	= node std sstar dstar rdc

all:	lib $(EXE)

Makefile.inc:	Makefile.inc.conf
	@sed s#__BASE__#$(PWD)#g Makefile.inc.conf >Makefile.inc

lib:	Makefile.inc
	@for d in $(DIRS) ; do echo '\nmake' $@ in $$d; make -C $$d $@ ; done

clean:	Makefile.inc
	@for d in $(DIRS) ; do echo '\nmake' $@ in $$d; make -C $$d $@ ; done
	/bin/rm -f *~ $(EXE) Makefile.inc
