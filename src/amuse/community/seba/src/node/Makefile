DIRS	= util dyn

lib:
	@for d in $(DIRS) ; do echo node: make $@ in $$d ; make -C $$d $@ ; done

clean:
	@for d in $(DIRS) ; do echo node: make $@ in $$d ; make -C $$d $@; done
	/bin/rm -f *~ *.a
