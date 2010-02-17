This directory contains the AMUSE software. With AMUSE you
can write scripts to simulate astrophysical problems in
different domains.

AMUSE documentation and software can be found at:

    http://castle.strw.leidenuniv.nl/trac/amuse


Getting Started
===============

To build amuse you need a working build environment and
install some prerequisites. This documet contains 
the quick install instructions, if these fail please look
at the detailed descriptions of the installation 
procedure in the documents in the 'doc/install' directory. 

Compilers
=========

To build AMUSE from source you need to have a working 
build environment. The AMUSE build system needs
a C++ and fortan 90 compiler. Please check first if you
have a working build environment on your system.

In Ubuntu you can setup the environment with (as root):

    apt-get install build-essential curl g++ gfortran gettext zlib1g-dev


In Fedora you can setup the environment with (as root)::

    yum groupinstall "Development Tools" "Development Libraries"
    
In OS X you can install fink or the Apple Developer Tools. 
The Apple Developer Tools do not include a gfortran compiler, 
you can find one at:

    http://r.research.att.com/tools/

Installing Prerequisites
========================

This document describes installation of the pre-requisite 
software packages to a user directory. If you have an
recent Ubuntu or Fedora distribution you can follow the installation
instructions in `doc/install/howto-install-prerequisites.txt` 
to install the packages as part of the system.

1. Make a prerequisite software directory (can be set to any directory)

    mkdir ~/amuse/prerequsites
    
2. Set the PREFIX, PATH and LD_LIBRARY_PATH environement variables

    export PREFIX=~/amuse/prerequisites
    export PATH=${PREFIX}/bin:${PATH}
    export LD_LIBRARY_PATH=${PREFIX}/lib:${LD_LIBRARY_PATH}

2b. If you have ifort and/or icc, you need to also set:  

    export LIBRARY_PATH=${PREFIX}/lib:${LIBRARY_PATH}

3. Download and install python

    cd doc/install
    ./install-python.sh

4. Download and install the other pre-requisites
   (script is also in the `doc/install` directory)

    ./install.py install
    
    
Set Environment
===============
You can set the the PREFIX, PATH and LD_LIBRARY_PATH
environment variables in you bashrc file. Please make sure
the ${PREFIX}/bin directory is first in the path.

In bash, you can extend your `.bashrc` file with:

export PREFIX=~/amuse/prerequisites
export PATH=${PREFIX}/bin:${PATH}
export LD_LIBRARY_PATH=${PREFIX}/lib:${LD_LIBRARY_PATH}

Building AMUSE
==============

All modules can be build with a configure and make command. Start make
from the main amuse directory (directory this README file lives in).
The main task of the configure script is to check if the
prerequisite packages have been installed.

1. Configure the source code 

    ./configure 
    
1. Build the code with make

    make 
    

Starting the MPI daemon process
===============================
For MPICH2, you need to start the mpd daemon process. If run for
the firt time, the damon will complain about a .mpd.conf file. 
If so, please create the .mpd.conf file as instructed by the 
mpd command.

1. Start the mpd daemon

    mpd &
    

Running the tests
=================
AMUSE comes with a large set of tests, most can be run automatically. 
To run these tests start the nosetests command from the main 
amuse directory (directory this README file lives in).

To run these tests do:

1. Run the automatic tests

    nosetests -v
    
