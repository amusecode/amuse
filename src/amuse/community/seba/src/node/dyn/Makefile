DIRS	= evolve init kepler util

lib:
	@ for d in $(DIRS) ; do echo star: make $@ in $$d ; make -C $$d $@ ; done

clean:
	@for d in $(DIRS) ; do echo star: make $@ in $$d ; make -C $$d $@; done
	/bin/rm -f *~ *.a
