#
# support/config.py.  Generated from configpy.in by configure.
#

class interpreters(object):
    python = r'/disks/strw4/vanelteren/env/python3/bin/python3' 

class compilers(object):
    cxx = 'g++' 
    cc  = 'gcc'
    fc = 'gfortran'
    
    cxx_flags = '-g -O2'
    cc_flags  = '-g -O2'
    fc_flags = '-g -O2'
    
    found_fftw = 'yes'
    fftw_flags = '-I/disks/strw4/vanelteren/env/python3/include'
    fftw_libs = '-L/disks/strw4/vanelteren/env/python3/lib -lfftw3  -lfftw3_threads'
    
    found_gsl = 'yes'
    gsl_flags = '-I/disks/strw4/vanelteren/env/python3/include  '
    gsl_libs = '-L/disks/strw4/vanelteren/env/python3/lib -lgsl -lgslcblas -lm  '
    
    gfortran_version = '4.6.3'
    ifort_version = ''
       

class mpi(object):
    is_enabled = 'yes'=='yes'
    mpicxx = 'mpic++' 
    mpicc  = 'mpicc'
    mpif95 = 'mpif90'

class java(object):
    jni_includes = 'none'
    jdk = '/usr'

class cuda(object):
    is_enabled   = 'no'=='yes'
    compiler     = ''
    toolkit_path = '/NOCUDACONFIGURED'
    sdk_path     = '@CUDA_SDK@'
    cuda_libs = '-L/NOCUDACONFIGURED cuda cudart'
    sapporo_version = 'light'
    
class openmp(object):
    is_enabled   = 'yes'=='yes'
    fcflags = '-fopenmp'
    cflags = '-fopenmp' 
    
class modules(object):
    have_matplotlib = 0==1
    have_h5py       = 1==1
