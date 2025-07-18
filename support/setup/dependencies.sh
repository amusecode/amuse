# Dependencies
#
# These should list all the packages required to make the following available,
# including development headers.
#
# - GCC C compiler
# - GCC C++ compiler
# - GCC Fortran compiler
# - coreutils
# - Python 3
# - PkgConfig
# - wget or curl
# - patch
# - GNU tar
# - unzip
# - gunzip
# - bunzip2
# - unxz
# - Perl
# - bison
# - CMake
# - OpenMPI (or another implementation)
# - GNU Scientific Library
# - FFTW3
# - GNU MP
# - GNU MPFR
# - HDF5
# - NetCDF4
# - QHull
# - HEALPix C++
# - BLAS
# - LAPACK

# MESA r15140 fails its tests when compiled with gfortran 14

DEPS_conda="c-compiler cxx-compiler fortran-compiler 'gfortran<14' python pkgconfig"
DEPS_conda="${DEPS_conda} coreutils patch"
DEPS_conda="${DEPS_conda} curl tar unzip gzip bzip2 xz perl bison make cmake openmpi"
DEPS_conda="${DEPS_conda} gsl fftw gmp mpfr hdf5 netcdf4 qhull healpix_cxx libopenblas"
DEPS_conda="${DEPS_conda} liblapack zlib 'docutils>=0.6' 'mpi4py>=1.1.0' 'numpy>=1.2.2'"
DEPS_conda="${DEPS_conda} 'h5py>=1.1.0'"

DEPS_macports="gcc12 python312 pkgconfig curl gpatch gnutar unzip gzip bzip2 xz perl5"
DEPS_macports="${DEPS_macports} gmake cmake openmpi-gcc12 gsl fftw-3 gmp mpfr hdf5"
DEPS_macports="${DEPS_macports} netcdf netcdf-fortran qhull healpix-cxx openblas lapack"

DEPS_POST_macports="

sudo port select --set gcc mp-gcc12
sudo port select --set python3 python312
sudo port select --set mpi openmpi-gcc12-fortran
"

DEPS_homebrew="gcc@13 python pkg-config curl gpatch gnu-tar unzip gzip"
DEPS_homebrew="${DEPS_homebrew} bzip2 xz"
DEPS_homebrew="${DEPS_homebrew} perl bison make cmake open-mpi gsl fftw gmp mpfr hdf5"
DEPS_homebrew="${DEPS_homebrew} netcdf netcdf-cxx netcdf-fortran qhull healpix openblas"
DEPS_homebrew="${DEPS_homebrew} lapack"

DEPS_apt="gcc g++ gfortran python3 python3-dev pkg-config curl patch tar unzip gzip"
DEPS_apt="${DEPS_apt} bzip2 xz-utils"
DEPS_apt="${DEPS_apt} perl bison make cmake libopenmpi-dev openmpi-bin"
DEPS_apt="${DEPS_apt} libgsl-dev libfftw3-dev libgmp3-dev libmpfr6"
DEPS_apt="${DEPS_apt} libmpfr-dev libhdf5-dev hdf5-tools libnetcdf-dev libqhull-dev"
DEPS_apt="${DEPS_apt} libhealpix-cxx-dev liblapack-dev libblas-dev"

DEPS_dnf="gcc gcc-c++ gcc-gfortran python3 python3-devel pkgconf-pkg-config curl patch"
DEPS_dnf="${DEPS_dnf} tar unzip gzip"
DEPS_dnf="${DEPS_dnf} bzip2 xz perl-core bison make cmake openmpi-devel"
DEPS_dnf="${DEPS_dnf} gsl-devel fftw-devel gmp-devel mpfr-devel hdf5-devel netcdf-devel"
DEPS_dnf="${DEPS_dnf} netcdf-cxx-devel netcdf-fortran-devel healpix-devel"
DEPS_dnf="${DEPS_dnf} libqhull-devel blas-devel lapack-devel"


# Help messages for the user, used by help.sh

CONDA_CMDS="${BOLD}conda install -c conda-forge ${DEPS_conda}${END_BOLD}\n"
MACPORTS_CMDS="${BOLD}sudo port install ${DEPS_macports}\n${DEPS_POST_macports}${END_BOLD}\n"
HOMEBREW_CMDS="${BOLD}brew install ${DEPS_homebrew}${END_BOLD}\n"
APT_CMDS="${BOLD}sudo apt install ${DEPS_apt}${END_BOLD}\n"
DNF_CMDS="${BOLD}sudo dnf install ${DEPS_dnf}${END_BOLD}\n"

