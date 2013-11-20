# Makefile for RAMSES utils
F90=ftn
BINDIR=.
FLAGS=-traceback

# Make all targets
all: $(patsubst %.f90, %, $(wildcard *.f90))

# Make a specific f90 file
%: %.f90
	$(F90) $(FLAGS) $< -o $(BINDIR)/$@
