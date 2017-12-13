ETICS
=====

This is ETICS (Expansion Techniques In Collisionless Systems), a GPU (currently
**CUDA only**) N-body code which uses series expansion to calculate the
gravitational field. See more details in this publication:

Meiron, Y., Li, B., Holley-Bockelmann, K., & Spurzem, R. 2014, ApJ, 792, 98

http://adsabs.harvard.edu/abs/2014ApJ...792...98M


What's inside
-------------

- ETICS standalone program
- ETICS static library
- ETICS module for AMUSE


Prerequisites
-------------

- CUDA (>= 6; mandatory)
- HDF5 (optional)
- Boost (optional)
- AMUSE (mandatory only for AMUSE module)

To disable HDF5 [insert explanation here]
To disable Boost [insert explanation here]


Compilation
-----------

### Standalone program

    make standalone

builds in `src/`.


### Static library

    make library

builds in `src/`.


### AMUSE module

    make

builds in top level directory. The whole `etics` directory has to be placed (or
linked) inside:

    $AMUSE_DIR/src/amuse/community/


How to use
----------

The `file.ini` is a self-explanatory input file; if compiled without Boost, fill
in the relevant variables in file `noboost.inc` which is compiled into the
executable (any change of parameters requires re-compilation). Start simulation
with:

    ./etics file.ini

Any input file name is acceptable.


Known issues
------------

* No MEX

The MEX (Multipole Expansion) method is not available in this version; thus, the
SCF (Self-Consistent Field) method is the only expansion technique available.
The ETICS program has been heavily restructured and the MEX routines are no
longer compatible. Hopefully this will be fixed.

* Hardcoded launch configuration

For at least one of the CUDA kernels, for various reasons, it seems that "brute
force" search is needed to find the optimal launch configuration. Currently it
is hardcoded, and a primitive search routine is provided.

* Problem for particles with |θ| << 1

Due to using an unstable recursion relation to calculate the Associated Legendre
polynomials, particles in a narrow cone around the z-axis cannot be considered
accurately. This means that they give an erroneous contribution to the
gravitational field and also are assigned erroneous force and potential. The
size of this cone increases with the angular term of the expansion. To partly
solve this, the current (ugly) fix is to only consider particles with cos(θ) >
0.999 at the monopole level. This is not so bad because the monopole is always
the most dominant term (and is error free) and the number of particles in this
cone is small and they are transient (i.e. they come out of it usually in a
small number of steps). A somewhat better solution is to make this arbitrary
cutoff of 0.999 l-dependent, and an even better solution would be to use an
asymptotic expression for the Associated Legendre polynomials.
