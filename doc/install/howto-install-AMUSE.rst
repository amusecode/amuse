==================================
Installation of the AMUSE software
==================================

Before installing AMUSE the prerequisite software must be downloaded and
installed, see :ref:`prerequisite-label`.

In the current stage of development AMUSE will not be installed in 
the python ``site-packages`` library. Instead, all code is build 
in the AMUSE source directories. With this setup we can easily edit
the code and run it, without the need for an extra installation step.

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


Testing the build
-----------------

.. warning::
    
    For MPICH2 installations, the `mpd` process daemon must be 
    started befor testing the code. The `mpd` application manages 
    the creation of MPI processes. If this is the first time the 
    MPICH2 daemon is run it will complain about a missing 
    ``.mpd.conf`` file. Please follow the instructions printed by 
    the mpd daemon.

    .. code-block:: sh
        
        > mpd &

    If the mpd deamon only complains with 'no mpd.conf', these
    are the steps to take, to create a mpd.conf file:
    
    .. code-block:: sh
        
        > echo 'MPD_SECRETWORD=secret' > ~/.mpd.conf
        > chmod 600 ~/.mpd.conf
        
    Please make sure to replace '''secret'''.
    
    After starting `mpd` we can start the tests.
    

    
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
    The error starts like this:
    
    .. code-block:: sh
    
        unable to parse user arguments

        Usage: ./mpiexec [global opts] [exec1 local opts] : [exec2 local opts] : ...
    
    
.. warning::
    
    On some laptops the hostname will not point
    to the correct internet address. For these laptops 
    you can start the mpd daemon on the localhost ip. To do so,
    you need to set the ``--ifhn`` option:
    
    .. code-block:: sh
    
        > mpd --ifhn=localhost &
        
    
    
.. warning::
    
    On OS X, when you install the prerequisites with macports, 
    ``nosetests`` will not have a standard name. It will be named
    ``nosetests-<python-version>``. So for python2.7 you'll need to 
    use *nosetests-2.7*
    
    .. code-block:: sh
    
        > nosetests-2.7
        ............................................
        
        OK
        



Real-time testing
~~~~~~~~~~~~~~~~~
The code includes support for real-time testing. The real-time testing 
application monitors the files in the source directories ('src' 
and 'test'). Every time a file is changed it will run most of the tests.
After each test a report is created, this report can be viewed with
a web browser.

.. code-block:: sh

    # go to the AMUSE root directory
    # display help information of the realtime_test script
    >  python -m support.realtime_test --help
    Usage: realtime_test.py [options]

    Options:
      -h, --help            show this help message and exit
      -p PORT, --port=PORT  start serving on PORT
      
    # start the python realtime_test script on port 9080
    > python -m support.realtime_test -p 9080
    starting server on port:  9080
    start test run
    ...
    # open a browser to view the results
    > firefox http://localhost:9080/
    

Running the code
----------------
A python script will not find the AMUSE code as the code is not 
installed into the python 'site-packages' directory or any other 
directory that can be found by python automatically. 

During a build a shell script is created to run the AMUSE code. To 
use this script you first have to copy it to a directory in your PATH.
The script is called ''amuse.sh''. After copying this script you can run
amuse code from anywhere on your disk by starting 'amuse.sh'. This
script has exactly the same command line parameters as the normal python
application.

.. code-block:: sh

    > amuse.sh
    Python 2.6.2 (r262:71600, Sep  1 2009, 16:14:27) 
    [GCC 4.3.2 20081105 (Red Hat 4.3.2-7)] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> from amuse.units import units
    >>> units.m
    unit<m>
    
    

    

