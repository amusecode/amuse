.PHONY: default test clean install

ifndef evpath
	evpath=${PWD}/../../
endif

ifndef Z
	Z=02
endif

ifndef EV
	EV=../../code/ev
endif

default:zahb.mod

test: zahb.mod test.dat test.run
	${EV} test ${Z}

clean:
	rm -rf *.mod *.mdl1 *.plt1 *.io12 *.last* fort.* *.mas *.out *.list *.log *.out1

m225.mod:m225.run m225.dat
	${EV} m225 ${Z}

m225.last1:m225.mod
	tail -n200 $< > $@

strip.mod:strip.dat strip.run m225.last1
	${EV} strip ${Z}

relax.start: strip.mod
	tail -n200 $< > $@

relax.mod: relax.start relax.run relax.dat
	${EV} relax ${Z}

zahb.mod: relax.mod
	tail -n200 $< > $@

install: zahb.mod
	install zahb.mod ${evpath}/input/zahb${Z}.mod
