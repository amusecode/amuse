This directory contains the AMUSE software. With AMUSE you
can write scripts to simulate astrophysical problems in
different domains.

Getting Started
===============

A detailed description of the installation procedure can be
found in the documents in the 'doc/install' directory. 

Compilers
=========

To build AMUSE from source you need to have a working 
build environment. The AMUSE build system needs
a C++ and fortan 90 compiler. Please check first if you
have a working build environment on your system.

In Ubuntu you can setup the environment with (as root):

    apt-get install build-essential curl g++ gfortran gettext 


Installing Prerequisites
========================

This document describes installation of all 
prerequisite software in a user directory. If you have
Ubuntu or Fedora distribution you can follow the installation
instuctrions in `doc/install/howto-install-prerequisites.txt` 
to install the packages as part of the system.

1. Make prerequisite software directory (can be set to any directory)

    mkdir ~/amuse/prerequsites
    
2. Set the PREFIX, PATH and LD_LIBRARY_PATH environement variables

    export PREFIX=~/amuse/prerequisites
    export PATH=${PREFIX}/bin:${PATH}
    export LD_LIBRARY_PATH=${PREFIX}/lib:${LD_LIBRARY_PATH}

3. Download and install python

    cd doc/install
    ./install-python.sh

4. Download and install other prerequisites
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

All modules can be build with a configure and make commands. Start make
from the main amuse directory (directory this README file lives in).
The main task of the configure script is to check if the
prerequisite packages have been installed.

1. Configure the source code 

    ./configure 
    
1. Build the code with make

    make 
    

Starting the MPI daemon process
===============================
For MPICH2, you need to start the mpd daemon process. If this
is the firt time the damon is run, it will complain about a
.mpd.conf file. If so, please create the .mpd.conf file as
instructed by the mpd command.

1. Start the mpd daemon

    mpd &
    

Running the tests
=================
AMUSE comes with a large set of tests, most can be run automatically. 
Start the nosetests command from the main 
amuse directory (directory this README file lives in).

To run these tests do:

1. Run the unit tests

    nosetests -v
    
