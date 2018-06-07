Installing on Arch Linux
========================

In this section we asume a default Arch Linux installation.

All
---
The prerequisites can be installed with a couple of commands
on Arch Linux. 

To install the prerequisites do (for base-devel select *all* members)::

    > sudo pacman -Syu base-devel curl gcc-fortran gettext zlib

Install python and dependencies::

    > sudo pacman -Syu python2 python2-numpy \
      hdf5 docutils openmpi
      python2-mpi4py python2-nose\
      fftw gsl cmake gmp mpfr

To install h5py, first install distribute and then run easy_install::

    > sudo pacman -Syu python2-distribute
    
    > sudo easy_install-2.7 h5py
