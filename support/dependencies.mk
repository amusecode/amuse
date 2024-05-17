# Dependencies
#
# These should list all the packages required to make the following available,
# including development headers.
#
# - GCC C compiler
# - GCC C++ compiler
# - GCC Fortran compiler
# - Python 3
# - PkgConfig
# - CMake
# - OpenMPI (or another implementation)
# - GNU Scientific Library
# - FFTW3
# - GNU MP
# - GNU MPFR
# - HDF5
# - NetCDF4
# - BLAS
# - LAPACK

DEPS_macports := gcc12 python39 cmake openmpi-gcc12 gsl fftw-3 gmp mpfr hdf5 netcdf
DEPS_macports += netcdf-fortran openblas lapack

define DEPS_POST_macports :=

sudo port select --set gcc mp-gcc12
sudo port select --set python3 python39
sudo port select --set mpi openmpi-gcc12-fortran

endef

DEPS_homebrew := gcc@12 gfortran@12 pkg-config cmake open-mpi gsl fftw gmp mpfr hdf5
DEPS_homebrew += netcdf netcdf-cxx netcdf-fortran openblas lapack

DEPS_apt := gcc g++ gfortran python3 pkg-config cmake libopenmpi-dev openmpi-bin
DEPS_apt += libgsl-dev libfftw3-3 libfftw3-dev libgmp3-dev libmpfr6 libmpfr-dev
DEPS_apt += libhdf5-serial-dev hdf5-tools libnetcdf-dev liblapack-dev libblas-dev

DEPS_dnf := gcc gcc-c++ gcc-gfortran python3 pkgconf-pkg-config cmake openmpi-devel
DEPS_dnf += gsl-devel fftw-devel gmp-devel mpfr-devel hdf5-devel netcdf-devel
DEPS_dnf += netcdf-cxx-devel netcdf-fortran-devel blas-devel lapack-devel

DEPS_conda := gcc gxx gfortran openjdk python pkgconfig coreutils cmake openmpi
DEPS_conda += openmpi-mpicc openmpi-mpicxx openmpi-mpifort openmp gsl fftw gmp mpfr hdf5
DEPS_conda += netcdf4 libopenblas liblapack


# Help messages for the user, used by help.mk

define CONDA_CMDS :=

    conda install -c conda-forge $(DEPS_conda)

endef


define MACPORTS_CMDS :=

    sudo port install $(DEPS_macports)
$(DEPS_POST_macports)

endef


define HOMEBREW_CMDS :=

    sudo brew install $(DEPS_homebrew)

endef


define APT_CMDS :=

    sudo apt install $(DEPS_apt)

endef


define DNF_CMDS :=

    sudo dnf install $(DEPS_dnf)

endef

