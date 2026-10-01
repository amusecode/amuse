.. _adding_example_code_tidymess:

Adding a C++ N-Body Code to AMUSE
=================================

In this tutorial, we will create an interface from scratch for the
TIdal DYnamics of Multi-body ExtraSolar Systems code, or TIDYMESS,
written by Dr. Tjarda Boekholt and Dr. Alexandre Correia
(https://doi.org/10.48550/arXiv.2209.03955). This code implements detailed
tidal forces into an N-body code to track the deformation of bodies.
This community code has already been implemented into AMUSE so you
can follow along this tutorial.

.. NOTE::

    In this guide, TIDYMESS refers to the standalone simulation package,
    while ``Tidymess`` refers to the TIDYMESS package inside of AMUSE.


Getting Started
===============

This tutorial assumes you have a working amuse or amuse development build,
preferably in seperated environment (virtualenv, venv or conda etc).
Please ensure that amuse is set up correctly, this can be verified by running the
``amusifier`` .

.. code-block:: bash

    > amusifier --help

Naming our project
~~~~~~~~~~~~~~~~~~
PEP-8 naming conventions for classes follows PascalCase, so we will name our project ``Tidymess``.

Creating the initial directory structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To start, we need to create the directory structure for ``Tidymess``, along with the
necessary files to build our interface. The fastest method to setup the directory is
by using the ``amusifier`` script with ``--mode=dir``.

Since TIDYMESS is a native C++ code with no other dependencies, we will specify
``--type=c``, but the ``amusifier`` can also build the interface directory for
``f90`` and ``python`` codes. Make sure you run the ``amusifier`` in the ``amuse/src/`` directory.

.. code-block:: bash

    > cd amuse/src/
    > amusifier --type=c --mode=dir Tidymess

Having run the ``amusifier``, we now have our new directory in ``amuse/src/amuse_tidymess``.
There should be all the required folders for building our interface, as well as a few code
stubs to expand upon.

Building the code
=================
Before we start working on the interface, we should try to install and compile TIDYMESS
inside of ``AMUSE``.

Defining dependencies
~~~~~~~~~~~~~~~~~~~~~
The AMUSE build system needs to know what packages and libraries our project depends on.
Navigate to ``amuse_tidymess/packages/amuse_tidymess.amuse_deps``, which is where we define every
dependency we will need. By default it will look like:

.. code-block:: text

   c c++ fortran java python cmake install download mpi openmp cuda opencl x11 opengl blas lapack gsl gmp mpfr fftw hdf5 netcdf4

Since ``TIDYMESS`` is a standalone C++ code, we can delete most of those and simplify our dependencies to:

.. code-block:: text

   c++

To test that we did everything properly, we can run ``./setup`` from the top-level directory of AMUSE

.. code-block:: text

    > cd amuse/
    > ./setup

    Checking for dependencies, one moment please...

    *** Configuration complete ***
    Detected features: c c++ fortran python python-dev gmake cmake install
    download patch tar unzip gunzip bunzip2 unxz perl bison mpi openmp blas
    lapack gsl gmp mpfr fftw libz hdf5 netcdf4 qhull healpix-cxx

    ** Enabled packages **

    Packages marked i) are currently installed.

    i) amuse-framework           amuse-hop                 amuse-petar
       amuse-aarsethzare         amuse-huayno-openmp       amuse-ph4
       amuse-adaptb              amuse-huayno              amuse-phantom
       amuse-athena              amuse-kepler              amuse-phigrape
       amuse-bhtree              amuse-kepler-orbiters     amuse-rebound
       amuse-brutus              amuse-krome               amuse-sakura
       amuse-bse                 amuse-mameclot            amuse-seba
       amuse-capreole            amuse-mercury             amuse-secularmultiple
       amuse-evtwin              amuse-mesa-r15140         amuse-sei
       amuse-fastkick            amuse-mesa-r2208          amuse-simplex
       amuse-fi                  amuse-mi6                 amuse-smalln
       amuse-fractalcluster      amuse-mikkola             amuse-sphray
       amuse-gadget2             amuse-mmams               amuse-sse
       amuse-galactics           amuse-mobse               amuse-symple
       amuse-galaxia             amuse-mocassin            amuse-tidymess
       amuse-halogen             amuse-mosse               amuse-tupan
       amuse-hermite             amuse-mpiamrvac           amuse-twobody
       amuse-hermite-grx         amuse-nbody6xx            amuse-vader

We can see that the AMUSE framework is correctly installed, indicated by the ``i)`` symbol. We can
also see our package ``amuse-tidymess`` is listed meaning the AMUSE build system now knows of our project
and we are ready to move on!

Setting up Autoconf
~~~~~~~~~~~~~~~~~~~
We now need to determine what compilers and libraries are on the system and how to use them. For this we will edit
the ``configure.ac`` file in ``amuse_tidymess/support/``. This file contains a set of macros which will
detect the tools and libraries needed to build our package. The template should contain all the macros
needed for our package, so its just a matter of deleting what we don't need. Delete any comment prefaced
with ``#####``, but only after following the direction of the comment. Since TIDYMESS is a native C++ code,
most of the optional macros can be deleted, especially the ones related to external libraries like ``CUDA``,
``MPI``, ``FFTW``, etc... The only optional macros we will keep are ``AMUSE_LIB_STOPCOND()`` for enabling stopping
conditions in our project, as well as ``AMUSE_DOWNLOAD()`` and ``AC_CHECK_TOOL(TAR, tar)`` for dynamically
downloading TIDYMESS from github into our project when the user tries to install it (so we don't need to package
the source code directly into AMUSE).

Once ``configure.ac`` is setup correctly, we can edit ``config.mk.in``, and remove any unneeded variables
as well as any ``#####`` comment. This file is a template for ``config.mk``, which will contain a description
of all the compiler and library variables needed for our package. Just like ``configure.ac``, the only optional variables
we will keep are again related to the ``STOPCOND`` library and downloading our package.

Once these files are cleaned up, run ``autoreconf`` to (re)create the ``configure`` script, then run
``./configure``. This will test the detection and check for errors. As a sanity check, run ``cat config.mk``
and ensure that there are no ``@VARIABLE@`` symbols left! If there are, check that ``configure.ac`` and
``config.mk.in`` were setup correctly.

.. code-block:: text

    > cd amuse/src/amuse_tidymess/support/
    > autoreconf
    > ./configure
    > cat config.mk

    # Compilers
    CXX = arm64-apple-darwin20.0.0-clang++

    MPICXX = mpicxx

    CPU_COUNT = 8

    # Tools
    AR = arm64-apple-darwin20.0.0-ar
    RANLIB = arm64-apple-darwin20.0.0-ranlib
    DOWNLOAD = curl -L
    TAR = tar

    # AMUSE framework libraries
    STOPCOND_CFLAGS =  -I/Users/astro/miniconda3/envs/tidymess/include
    STOPCOND_LIBS = -lstopcond

.. WARNING::

    Make sure the the ``amuse_tidymess/support/shared/`` folder is a symlink to ``amuse/support/shared/``
    to ensure that there is no code duplication in the codebase, and that bug fixes are propagated to each
    package automatically. This should be done automatically by the ``amusifier`` but can be a source of bugs
    if not setup correctly. You can verify it with:

    .. code-block:: console

        > cd amuse/src/amuse_tidymess/support/
        > ls -l

    You should see:

    .. code-block:: text

        lrwxr-xr-x  1 user user  24 May 13 16:16 shared -> ../../../support/shared/


Setting up the Makefile
=======================
The ``amusifier`` already created our ``Makefile``
for us in ``src/amuse_tidymess/``, which has most of the code we will need for this step. The ``Makefile``
has several responsibilities: downloading and patching the community code source, compiling it into a static library,
and linking that library with the auto-generated worker stub to produce the ``tidymess_worker``.

Downloading and patching the source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AMUSE does not bundle community code sources. The ``Makefile`` downloads a
pinned release at install time, which decouples AMUSE from upstream changes.
Three edits are needed:

1. Set ``VERSION`` to the upstream commit hash or tag.
2. Set the URL in the ``tidymess.tar.gz`` rule to match the upstream archive.
3. Set the extraction target. This rule unpacks the archive and moves the
   resulting directory into ``src/``. The extracted directory is named
   ``tidymess-$(VERSION)``, so the generated target
   ``src/tidymess-$(VERSION)`` changes with every version. Renaming it with
   ``mv`` to the fixed name ``src/tidymess`` keeps all later paths
   independent of ``VERSION``.

.. code-block:: make

    VERSION = 4f97bfe11e8c638fdda744ca288e57565efe718a

    tidymess.tar.gz:
        $(DOWNLOAD) https://github.com/tidymess-code/tidymess/archive/$(VERSION).tar.gz >$@

    src/tidymess: tidymess.tar.gz
        tar xf $<
        mv tidymess-$(VERSION) src/tidymess

If for some reason the source code needs to be modified, the ``Makefile`` allows us to
include ``.patch`` files to override the source code. Tidymess needs no patches, so we can
safely skip this step.

Selecting compiler flags
~~~~~~~~~~~~~~~~~~~~~~~~
The generated ``Makefile`` adds the flags and libraries of every dependency
AMUSE supports (``DEPFLAGS`` for compilation, ``LDLIBS`` for linking). The
values come from ``config.mk``, which ``./configure`` generates. Edits:

1. Delete every dependency the code does not use. Each retained one adds a
   link-time requirement and can fail the build on systems lacking that library.
   Since Tidymess has no external libraries, delete any flag related to CUDA,
   GSL, GMP, FFTW, etc...
2. Keep only the variants matching the code's language. TIDYMESS is a C++ code
   so any Fortran related flags should be deleted.
3. Append ``DEPFLAGS`` to the variable matching the compiler: ``CXXFLAGS``
   for C++ (the template uses ``CFLAGS``).
4. Remove linker-related lines the code does not need (``LDLIBS``).
   We keep ``STOPCOND`` to support AMUSE stopping conditions.

.. code-block:: make

    DEPFLAGS += $(STOPCOND_CFLAGS)
    CXXFLAGS += $(DEPFLAGS)

    LDLIBS += -lm $(STOPCOND_LIBS)

Building the static library
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: make

    CODELIB = src/libtidymess.a

    .PHONY: $(CODELIB)
    $(CODELIB): | src/tidymess
        $(MAKE) -C src -j $(CPU_COUNT) all

- ``CODELIB``: path of the static library, relative to this ``Makefile``.
- ``.PHONY``: the outer ``Makefile`` cannot see which files the inner build
  depends on. Marking the library phony forces ``make -C src`` to run on every
  invocation. The inner ``Makefile`` then decides what to recompile.
- ``| src/tidymess``: an order-only prerequisite. Normally, ``make`` rebuilds a target whenever one
  of its prerequisites is newer than the target. An order-only prerequisite only
  has to exist before the recipe runs, and its timestamp is ignored. In our case,
  this means we will only redownload the source code if it is missing. This is very important
  as otherwise any local edits would be overwritten when compiling!
- ``-j $(CPU_COUNT)``: compile in parallel.

The inner Makefile
~~~~~~~~~~~~~~~~~~
The outer ``Makefile`` delegates compilation to a second ``Makefile`` in
``src/amuse_tidymess/src/``. The ``amusifier`` also generates this one. It
must provide:

- ``all``: builds the static library ``libtidymess.a`` in ``src/``.
- ``clean``: removes the build products.

.. code-block:: make

    CFLAGS   += -Wall -g
    CXXFLAGS += $(CFLAGS)

    CODELIB := libtidymess.a

    SRCDIR := tidymess/integrator/src
    INCDIR := tidymess/integrator/include

    SRC := $(wildcard $(SRCDIR)/*.cpp)
    CODEOBJS := $(SRC:$(SRCDIR)/%.cpp=%.o)

    all: $(CODELIB)

    clean:
    	rm -f *.o *.a

    distclean: clean

    $(CODELIB): $(CODEOBJS)
    	rm -f $@
    	$(AR) -ruv $@ $(CODEOBJS)
    	$(RANLIB) $@

    %.o: $(SRCDIR)/%.cpp
        $(CXX) $(CXXFLAGS) -c -o $@ -I$(INCDIR) $<

TIDYMESS keeps its sources and headers in separate directories inside the
downloaded tree. ``SRCDIR`` tells ``make`` where to find the ``.cpp`` files,
and ``INCDIR`` is passed to the compiler with ``-I`` so that ``#include``
directives resolve. Both are relative to the inner ``Makefile`` and must be
adapted to the layout of the community code.

Building the worker
~~~~~~~~~~~~~~~~~~~
The worker is the executable AMUSE spawns to communicate with the community
code. Building it takes three steps.

1. Generate the stub. ``amusifier`` reads ``interface.py`` and writes
   ``tidymess_worker.h`` and ``tidymess_worker.cc``. Both are rebuilt when
   ``interface.py`` changes.
2. ``tidymess_worker.o`` compiles the stub. ``interface.o``
   compiles ``interface.cc``, which implements the functions the stub calls.
   Its ``-I`` flag must point to the code's headers. It depends on
   ``| src/tidymess`` so that the headers exist before compilation.
3. Link. ``tidymess_worker`` combines both objects with ``$(CODELIB)``.

Package targets
~~~~~~~~~~~~~~~

.. code-block:: make

    # Building the workers
    tidymess_worker.h: interface.py
    	amusifier --type=h interface.py TidymessInterface -o $@

    tidymess_worker.cc: interface.py
    	amusifier --type=c interface.py TidymessInterface -o $@

    tidymess_worker.o: tidymess_worker.cc tidymess_worker.h
    	$(MPICXX) -c -o $@ $(CXXFLAGS) $<

    tidymess_worker: interface.o tidymess_worker.o $(CODELIB)
    	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

    interface.o: interface.cc tidymess_worker.h | src/tidymess
    	$(MPICXX) -o $@ -c -I src/tidymess/integrator/include $(CXXFLAGS) $<

- ``amuse-tidymess_contains``: lists the workers the package ships. Add one worker per
    line for multi-worker packages. It guarantees the workers are built before ``pip`` installs the package.

The final bits of code at the end of the Makefile are for building and installing the package, as well
as defining how to uninstall and cleanup the package.

.. code-block:: make

    # Which packages contain which workers?
    amuse-tidymess_contains: tidymess_worker

    # Building and installing packages
    develop-%: %_contains
    	support/shared/uninstall.sh $*
    	python -m pip install -e packages/$*

    install-%: %_contains
    	support/shared/uninstall.sh $*
    	python -m pip install packages/$*

    package-%: %_contains
    	python3 -m pip install -vv --no-cache-dir --no-deps --no-build-isolation --prefix ${PREFIX} packages/$*

    test-%:
    	cd packages/$* && pytest

    # Cleaning up
    .PHONY: clean
    clean:
	$(MAKE) -C src clean
	rm -rf *.o *worker*

    .PHONY: distclean
    distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache


Creating the Interfaces
=======================
With ``TIDYMESS`` compiled into ``AMUSE``, we can now begin the process of creating our interface!
The interface system allows community codes, which are all unique and depend on diverse libraries and
programming languages, to communicate with the ``AMUSE`` framework, which is native Python. ``AMUSE``
interfaces define a number of interface functions, which provide a standardized way for ``AMUSE`` to
communicate with each community code. This way, the experience of using any ``AMUSE`` code is identical:
all codes can be evolved with the ``evolve_model`` method, particles are represented as
``amuse.datamodel.particles``, etc... The strength of ``AMUSE`` lies in its ability to prototype quickly:
If the user wants to see what solution a different code would give for the same calculation, all they have
to do is switch which code they are using, and the script most likely does not need to change.

Therefore, our job when creating an interface is to map the community code functions to the ``AMUSE`` interface
functions. The amusifier already created all the files we need: the ``interface.py`` and the ``interface.cc``.


Defining the Python interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``interface.py`` actually defines two classes: the high-level and low-level Python interfaces. The high-level
interface is what the user interacts with when using the community code. This defines the methods, parameters,
and properties of each community code, and defines how the ``amuse.datamodel.particles`` work within that code.
The ``amusifer`` generates a minimal code stub for the ``interface.py``, which is enough to get us started.
The process of defining our interface starts with figuring out what type of code we are adding to ``AMUSE``.
``AMUSE`` has a set of predefined Python interfaces we can use to build our interface from.


+-----------------------------------+-------------------------------------------+
| Interface:                        | Example codes:                            |
+-----------------------------------+-------------------------------------------+
| ``GravitationDynamicsInterface``  | N-body: Tidymess, Ph4, Huayno             |
+-----------------------------------+-------------------------------------------+
| ``HydrodynamicsInterface``        | Hydrodynamical: Capreole                  |
+-----------------------------------+-------------------------------------------+
| ``MagnetoHydrodynamicsInterface`` | MHD: Athena                               |
+-----------------------------------+-------------------------------------------+
| ``StellarEvolutionInterface``     | Stellar: MESA, EVtwin, SeBa               |
+-----------------------------------+-------------------------------------------+

These interfaces define many
of the required interface functions so that we don't have to.
In this case, ``TIDYMESS`` is a N-Body code, and therefore falls under the
