==================================
Installation of the AMUSE software
==================================

Before installing AMUSE the prerequisite software must be downloaded and
installed, see :ref:`prerequisite-label`.

AMUSE can be installed from online packages using the python ``pip`` utility
or from a checkout of the repository. For the latter option, you can either
install it in a pip managed environment, or a self managed python environment.

Package install
===============

In this case, install the non-python prerequisites first then:

.. code-block:: sh
    
    > pip install amuse

this will attempt to fetch and install the amuse-framework package,
and the code packages (amuse-<code>) and python prerequisites. If one of the code fails to build,
you can install packages by hand.

Pip managed install from repository checkout
============================================

If you have a python environment where you can install packages with pip (it is 
recommended to compartementalize environments of different projects with virtualenv),
again install non-python prerequisites and checkout the AMUSE repository:

.. code-block:: sh

    > git clone https://github.com/amusecode/amuse
    > cd amuse
    > pip install -e .
    > python setup.py develop_build 

In this case, a link to the repository is created and codes are compiled in place. This is 
the recommended install for development.

Self-managed environment
========================

Lastly, you can install without pip. Again all code is build 
in the AMUSE source directories. This setup also allows you to easily edit
the code and run it, without the need for an extra installation step. There are bootstrap 
scripts in ``doc/install/`` to download and build all prerequisites if c++/fortran compilers are available.

Environment variables
---------------------

you need to tell python where to find the amuse package and the command line interpreter where to find
the ``amusifier`` executable.

.. code-block:: sh

    > export PYTHONPATH=${PYTHONPATH}:/path/to/amuse/src
    > export PATH=${PATH}:/path/to/amuse/bin
 
where you replace ``/path/to/amuse`` with the correct path to AMUSE.
 
Configuring the code
--------------------
The code is configured using the ``configure`` command. 
Before building the code, run 'configure' in the AMUSE
root directory.

.. code-block:: sh
    
    > ./configure
    
The 'configure' script will check for all prerequisite software
and report if any are missing.

Building the code
-----------------

The code is build using a  ``Makefile``. To build the code run 'make'
in the AMUSE root directory.

.. code-block:: sh
    
    > make clean
    > make
    ...
    community codes build
    ==================
    * sse
    * hermite0
    * bhtree
    * phiGRAPE
    running generate_main

If everything goes well all community codes will be build (e.g. sse, hermite0, 
bhtree, phiGRAPE and many others).

In order to use codes not stored in the AMUSE repository (e.g. MESA, ATHENA, Rebound and some others), the codes must be downloaded additionally.
This is done automatically after setting the environment variable DOWNLOAD_CODES to 1.
Alternatively, instead of a plain 'make' like in the example above you could do:

.. code-block:: sh

    > make DOWNLOAD_CODES=1

or:

.. code-block:: sh

    > make mesa.code DOWNLOAD_CODES=1
    > make athena.code DOWNLOAD_CODES=1

Running the code
----------------

You can quickly test your installation by importing some AMUSE module and running a code

.. code-block:: sh

    > python
    Python 3.7.3 (default, Dec  2 2019, 17:46:02) 
    [GCC 9.2.1 20190903 [gcc-9-branch revision 275330]] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>> from amuse.units import units
    >>> units.m
    unit<m>
    >>> from amuse.community.bhtree.interface  import Bhtree
    >>> b=Bhtree()
    >>> print(b)
    <amuse.community.bhtree.interface.BHTree object at 0x7f20f02e6dd8>


Testing the build
-----------------

The tests are run using the nosetests program.

.. code-block:: sh
    
    > nosetests
    ............................................
    Ran 91 tests in 12.013s

    OK


.. warning::

    If you have an MPICH2 installation but no mpd program your MPICH2
    installation has been configured for the Hydra process manager. 
    To run amuse scripts with the hydra process manager you must start
    every command with ``mpiexec``:
    
    .. code-block:: sh
        
        > mpiexec nosetests -v
    
    
    If you do not run under mpiexec you get an error with a usage statement.    

    

