PYTHON=/usr/bin/python


#
# PYTHON embedding
#
PYTHON_DEV=no
CYTHON=PYTHON_DEV='no'
PYTHONCONFIG=/usr/bin/python-config
PYTHONDEV_CFLAGS=
PYTHONDEV_LDFLAGS=


#
# Compilers, for compiling without MPI
# these compilers must be the same or compatible with
# the mpi compilers
#
CXX=g++
CC=gcc
FC=gfortran
GFORTRAN_VERSION=4.7.2
IFORT_VERSION=
FC_ISO_C_AVAILABLE=yes

#
# Default flags, append to these in the makefiles
#
CXXFLAGS=-g -O2
CFLAGS=-g -O2
FCFLAGS=-g -O2

#
# MPI Compilers
#
MPI_ENABLED=yes
MPICXX=mpic++
MPICC=mpicc
MPIFC=mpif90
MPIEXEC=/usr/bin/mpiexec

MPI_CFLAGS=-I/usr/lib/openmpi/include/openmpi -pthread
MPI_CXXFLAGS=-I/usr/lib/openmpi/include/openmpi -pthread
MPI_FCFLAGS=-pthread -I/usr/lib/openmpi/lib
MPI_CLIBS=-L/usr/lib/openmpi/lib -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl
MPI_CXXLIBS=-L/usr/lib/openmpi/lib -lmpi_cxx -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl
MPI_FCLIBS=-I/usr/lib/openmpi/lib -L/usr/lib/openmpi/lib -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl

#
# Java
#
JAVA_ENABLED=yes
JAVA=/usr/bin/java
JAVAC=/usr/bin/javac
JAR=/usr/bin/jar
JAVA_VERSION=1.7

#
# OpenMP
#

OPENMP_ENABLED=yes
OPENMP_FCFLAGS=-fopenmp
OPENMP_CFLAGS=-fopenmp
    
#
# Libraries: FFTW, GSL, HDF5, CUDA, SAPPORO
#

FOUND_FFTW=yes
FFTW_FLAGS= 
FFTW_LIBS=-lfftw3 -lm   -lfftw3_threads
    
FOUND_GSL=yes
GSL_FLAGS= 
GSL_LIBS=-lgsl -lgslcblas -lm  

FOUND_GMP=yes
GMP_FLAGS=
GMP_LIBS=-lgmp

FOUND_MPFR=yes
MPFR_FLAGS=
MPFR_LIBS=-lmpfr

HDF5_FLAGS=
HDF5_LIBS=-lhdf5

CUDA_ENABLED=no
NVCC=
NVCC_FLAGS=
CUDA_TK=/NOCUDACONFIGURED
CUDA_LIBS=-L/NOCUDACONFIGURED cuda cudart

FOUND_CL=no
CL_LIBS=
CL_FLAGS=

SAPPORO_VERSION=light
SAPPORO_LIBS=-L${AMUSE_DIR}/lib/sapporo_light -lsapporo 

FS_FLAGS=-I$(AMUSE_DIR)/lib/forsockets
FS_LIBS=-L$(AMUSE_DIR)/lib/forsockets -lforsockets -lforsocketsf

SC_FLAGS=-I$(AMUSE_DIR)/lib/stopcond
SC_CLIBS=-L$(AMUSE_DIR)/lib/stopcond -lstopcond
SC_FCLIBS=-L$(AMUSE_DIR)/lib/stopcond -lstopcondf
SC_MPI_CLIBS=-L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi
SC_MPI_FCLIBS=-L$(AMUSE_DIR)/lib/stopcond -lstopcondfmpi

export PYTHON CXX CC FC CXXFLAGS CFLAGS FCFLAGS MPICXX MPICC MPIFC
export JAVA JAVAC JAVAH JAR JAVA_FLAGS
export FS_FLAGS FS_LIBS SC_FLAGS SC_CLIBS SC_FCLIBS SC_MPI_CLIBS SC_MPI_FCLIBS
export FFTW_FLAGS FFTW_LIBS
export GSL_FLAGS GSL_LIBS
export HDF5_FLAGS HDF5_LIBS
export OPENMP_FCFLAGS OPENMP_CFLAGS
export NVCC CUDA_TK CUDA_LIBS CUDA_ENABLED
export SAPPORO_VERSION SAPPORO_LIBS
export MPICFLAGS MPICCFLAGS MPIFCFLAGS
export IFORT_VERSION
export NVCC_FLAGS
export PYTHONDEV_LIBS PYTHONDEV_CFLAGS PYTHONDEV_LIBS
export FC_ISO_C_ENABLED
