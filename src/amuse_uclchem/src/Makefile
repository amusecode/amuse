CODELIB = libuclchem.a

CODEOBJS = uclchem.o


all: $(CODELIB) 

clean:
	rm -f *.o *.a *.mod

distclean: clean

$(CODELIB): $(CODEOBJS)
	rm -f $@
	$(AR) -ruv $@ $(CODEOBJS)
	$(RANLIB) $@

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

