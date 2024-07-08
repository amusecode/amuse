AC_CONFIG_MACRO_DIRS([shared/m4])
AC_CONFIG_AUX_DIR([shared])
AC_INIT([amuse], [1.0])

AC_CANONICAL_TARGET

# If we have a Conda environment active, then we're going to do everything via Conda,
# including the dependencies. So in this case, we'll detect installed Conda packages,
# rather than searching the system for libraries in the usual way.
#
# If we're not in a Conda environment, then we'll detect features from the system, and
# hope that the user hasn't made too big a mess of things.
#
# Conda build environments for conda-forge have a different setup than normal Conda
# environments, with different packages, so that the normal Conda mode doesn't work.
# We use the normal detection as well in that case.
if test "x" != "x$CONDA_BUILD"
then
    # Conda list doesn't work in conda build environments, and normal detection
    # can be problematic because it picks up e.g. compilers from the system. So
    # we assume that the conda metadata is correct, and that if it isn't there's
    # an expert packages who can fix things.
    AS_ECHO(["Conda package build detected, skipping checks for dependencies."])
    ASSUME_FEATURES="gcc g++ gfortran java python install download mpi openmp cuda opencl blas lapack gsl gmp mpfr fftw hdf5 netcdf4"

elif test "x" != "x$CONDA_DEFAULT_ENV"
then
    AS_ECHO([])
    AS_ECHO(['    Conda environment active, dependencies will be searched in Conda only,'])
    AS_ECHO(['    and must be installed via Conda to be detected.'])
    AS_ECHO([])

    AX_CONDA_LIST()

    if test "x$target_os" == "xlinux-gnu"
    then
        AX_CONDA_PACKAGE(CC, gcc)
        AX_CONDA_PACKAGE(CXX, gxx)
    else
        AX_CONDA_PACKAGE(CC, clang)
        AX_CONDA_PACKAGE(CXX, clangxx)
    fi
    AX_CONDA_PACKAGE(FC, gfortran)
    AX_CONDA_PACKAGE(JAVAC, openjdk)
    AX_CONDA_PACKAGE(PYTHON, python)
    AX_CONDA_PACKAGE(INSTALL, coreutils)

    AMUSE_DOWNLOAD()
    AX_CONDA_PACKAGE(CMAKE, cmake)

    AX_CONDA_PACKAGE(MPICC, mpi)
    AX_CONDA_PACKAGE(HAVE_OPENMP, openmp)
    AX_CONDA_PACKAGE(CUDA_TK, cuda-toolkit)
    AX_CONDA_PACKAGE(FOUND_CL, khronos-opencl-icd-loader)
    AX_CONDA_PACKAGE(LAPACK_LIBS, liblapack)
    AX_CONDA_PACKAGE(BLAS, libblas)
    AX_CONDA_PACKAGE(FOUND_GSL, gsl)
    AX_CONDA_PACKAGE(FOUND_GMP, gmp)
    AX_CONDA_PACKAGE(FOUND_MPFR, mpfr)
    AX_CONDA_PACKAGE(FOUND_FFTW, fftw)
    AX_CONDA_PACKAGE(with_hdf5, hdf5)
    AX_CONDA_PACKAGE(FOUND_NETCDF, netcdf4)
    AX_CONDA_PACKAGE(SAPPORO_LIGHT_LIBS, sapporo_light)

    # Export conda list output so make doesn't have to rerun it to detect pip and wheel
    CONDA_LIST="$ax_conda_list"
    AC_SUBST(CONDA_LIST)

else
    AS_ECHO(['No Conda environment active, detecting dependencies from the system'])

    # The standard compiler checks AC_PROG_* crash if the compiler is not found, which
    # is not what we want here. Also, we only support specific compilers. So we use
    # AC_CHECK_TOOL instead. If the user sets the corresponding variable explicitly,
    # then this will just assume that it's a working compiler.
    AC_CHECK_TOOL(CC, gcc)
    AC_CHECK_TOOL(CXX, g++)
    AC_CHECK_TOOL(FC, gfortran)
    AC_CHECK_TOOL(JAVAC, javac)
    AC_CHECK_TOOL(PYTHON, python)
    AC_CHECK_TOOL(INSTALL, install)

    # Not all compilers have a f77 symlink, so if you have more than one installed,
    # you risk getting gcc and gfortran from one, and f77 from another, and then we
    # crash on a link error when runnig AX_LAPACK() below. So we set f77 to be the same
    # as FC, and avoid that problem, hopefully not causing any others...
    F77=${FC}

    AC_CHECK_TOOL(APT, apt)
    AC_CHECK_TOOL(DNF, dnf)
    AC_CHECK_TOOL(MACPORTS, port)
    AC_CHECK_TOOL(HOMEBREW, brew)

    AMUSE_DOWNLOAD()
    AC_CHECK_TOOL(CMAKE, cmake)

    AX_MPI()
    AX_OPENMP([HAVE_OPENMP="yes"])
    AX_CUDA()
    AX_CHECK_CL()

    # Avoid AX_BLAS crashing out due to lack of a Fortran compiler
    if test "x" != "x$FC"
    then
        AX_LAPACK()     # calls AX_BLAS automatically
    else
        AS_ECHO(["Skipping LAPACK/BLAS detection because we don't have Fortran"])
    fi

    AX_PATH_GSL(1.0, [FOUND_GSL=yes])
    AX_GMP()
    AX_MPFR()
    AX_FFTW()
    AX_LIB_HDF5()
    AX_LIB_NETCDF4()

    AMUSE_LIB_SAPPORO_LIGHT()
fi

# Create the FEATURES variable used by amuse.mk to decide what to build

# This macro adds its first argument to FEATURES if the variable
# in the second argument is defined.
m4_define([ENABLE_FEATURE_IF_DEFINED],
          [AS_IF([test "x" != "x${$2}"], [FEATURES="${FEATURES} $1"])])

# Add the first argument to FEATURES if the variable in the second argument has the
# value specified in the third argument.
m4_define([ENABLE_FEATURE_IF_EQUALS],
          [AS_IF([test "x${$2}" = "x$3"], [FEATURES="${FEATURES} $1"])])

FEATURES="$ASSUME_FEATURES"

ENABLE_FEATURE_IF_DEFINED([gcc], [CC])
ENABLE_FEATURE_IF_DEFINED([g++], [CXX])
ENABLE_FEATURE_IF_DEFINED([gfortran], [FC])
ENABLE_FEATURE_IF_DEFINED([java], [JAVAC])
ENABLE_FEATURE_IF_DEFINED([python], [PYTHON])
ENABLE_FEATURE_IF_DEFINED([install], [INSTALL])
ENABLE_FEATURE_IF_DEFINED([download], [DOWNLOAD])
ENABLE_FEATURE_IF_DEFINED([mpi], [MPICC])
ENABLE_FEATURE_IF_DEFINED([openmp], [HAVE_OPENMP])
ENABLE_FEATURE_IF_DEFINED([cuda], [CUDA_TK])
ENABLE_FEATURE_IF_EQUALS([opencl], [FOUND_CL], [yes])
ENABLE_FEATURE_IF_DEFINED([blas], [BLAS_LIBS])
ENABLE_FEATURE_IF_DEFINED([lapack], [LAPACK_LIBS])
ENABLE_FEATURE_IF_DEFINED([gsl], [FOUND_GSL])
ENABLE_FEATURE_IF_EQUALS([gmp], [FOUND_GMP], [yes])
ENABLE_FEATURE_IF_EQUALS([mpfr], [FOUND_MPFR], [yes])
ENABLE_FEATURE_IF_EQUALS([fftw], [FOUND_FFTW], [yes])
ENABLE_FEATURE_IF_DEFINED([hdf5], [HDF5_VERSION])
ENABLE_FEATURE_IF_EQUALS([netcdf4], [FOUND_NETCDF], [yes])
ENABLE_FEATURE_IF_DEFINED([sapporo_light], [SAPPORO_LIGHT_LIBS])

AC_SUBST(FEATURES)

AC_CONFIG_FILES([features.mk])
AC_OUTPUT
