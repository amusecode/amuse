# Dependencies
#
# These should list all the packages required to make the following available,
# including development headers.
#
# - GCC C compiler
# - GCC C++ compiler
# - GCC Fortran compiler
# - TODO: coreutils
# - Python 3
# - PkgConfig
# - wget or curl
# - TODO: patch
# - GNU tar
# - unzip
# - gunzip
# - bunzip2
# - unxz
# - Perl
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

DEPS_conda="gcc gxx gfortran openjdk python pkgconfig coreutils curl tar unzip"
DEPS_conda="${DEPS_conda} gzip bzip2 xz perl cmake openmpi"
DEPS_conda="${DEPS_conda} openmpi-mpicc openmpi-mpicxx openmpi-mpifort openmp gsl fftw"
DEPS_conda="${DEPS_conda} gmp mpfr hdf5 netcdf4 libopenblas liblapack"

DEPS_macports="gcc12 python312 pkgconfig curl gnutar unzip gzip bzip2 xz perl5 cmake"
DEPS_macports="${DEPS_macports} openmpi-gcc12 gsl fftw-3 gmp mpfr hdf5 netcdf"
DEPS_macports="${DEPS_macports} netcdf-fortran openblas lapack"

DEPS_POST_macports="

sudo port select --set gcc mp-gcc12
sudo port select --set python3 python312
sudo port select --set mpi openmpi-gcc12-fortran
"

DEPS_homebrew="gcc@12 gfortran@12 pkg-config curl gnu-tar unzip gzip bzip2 xz perl"
DEPS_homebrew="${DEPS_homebrew} cmake open-mpi gsl fftw gmp mpfr hdf5"
DEPS_homebrew="${DEPS_homebrew} netcdf netcdf-cxx netcdf-fortran openblas lapack"

DEPS_apt="gcc g++ gfortran python3 pkg-config curl tar unzip gzip bzip2 xz-utils perl"
DEPS_apt="${DEPS_apt} cmake libopenmpi-dev openmpi-bin"
DEPS_apt="${DEPS_apt} libgsl-dev libfftw3-3 libfftw3-dev libgmp3-dev libmpfr6"
DEPS_apt="${DEPS_apt} libmpfr-dev libhdf5-serial-dev hdf5-tools libnetcdf-dev"
DEPS_apt="${DEPS_apt} liblapack-dev libblas-dev"

DEPS_dnf="gcc gcc-c++ gcc-gfortran python3 pkgconf-pkg-config curl tar unzip gzip"
DEPS_dnf="${DEPS_dnf} bzip2 xz perl-core cmake openmpi-devel"
DEPS_dnf="${DEPS_dnf} gsl-devel fftw-devel gmp-devel mpfr-devel hdf5-devel netcdf-devel"
DEPS_dnf="${DEPS_dnf} netcdf-cxx-devel netcdf-fortran-devel blas-devel lapack-devel"


# Help messages for the user, used by help.sh

CONDA_CMDS="${BOLD}conda install -c conda-forge ${DEPS_conda}${END_BOLD}\n"
MACPORTS_CMDS="${BOLD}sudo port install ${DEPS_macports}\n${DEPS_POST_macports}${END_BOLD}\n"
HOMEBREW_CMDS="${BOLD}sudo brew install ${DEPS_homebrew}${END_BOLD}\n"
APT_CMDS="${BOLD}sudo apt install ${DEPS_apt}${END_BOLD}\n"
DNF_CMDS="${BOLD}sudo dnf install ${DEPS_dnf}${END_BOLD}\n"

