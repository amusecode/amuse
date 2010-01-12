This directory contains the AMUSE software. With AMUSE you
can write scripts to simulate astrophysical problems in
different domains.


Getting Started
===============

A detailed description of the installation procedure can be
found in the documents in the 'doc/install' directory. 

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

All modules can be build with one make command.

1. Build the code with make

    make 
