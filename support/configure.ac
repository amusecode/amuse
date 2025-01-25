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
    # an expert packager who can fix things.
    AS_ECHO(["Conda package build detected, skipping checks for dependencies."])
    ASSUME_FEATURES="c c++ fortran java python tar unzip gunzip bunzip2 unxz perl cmake install download mpi openmp cuda opencl blas lapack gsl gmp mpfr fftw hdf5 netcdf4 sapporo_light"

else
    if test "x" != "x$CONDA_DEFAULT_ENV"
    then
        AS_ECHO([])
        AS_ECHO(['    Conda environment active, dependencies will be searched in Conda only,'])
        AS_ECHO(['    and must be installed via Conda to be detected.'])
        AS_ECHO([])

        # Doesn't have a conda package yet, checked directly below.
        # AMUSE_CONDA_PACKAGE(SAPPORO_LIGHT_LIBS, sapporo_light)

        AMUSE_CONDA_LIST()

        if test "x$target_os" == "xlinux-gnu"
        then
            # CC needs to be a valid compiler for the sapporo check at the end
            AMUSE_CONDA_PACKAGE(CC, gcc, gcc)
            AMUSE_CONDA_PACKAGE(CXX, gxx)
        else
            AMUSE_CONDA_PACKAGE(CC, clang)
            AMUSE_CONDA_PACKAGE(CXX, clangxx)
        fi
        AMUSE_CONDA_PACKAGE(FC, gfortran)
        AMUSE_CONDA_PACKAGE(JAVAC, openjdk)
        AMUSE_CONDA_PACKAGE(PYTHON, python)
        AMUSE_CONDA_PACKAGE(PYTHON_DEV, python)
        AMUSE_CONDA_PACKAGE(INSTALL, coreutils)
        AMUSE_CONDA_PACKAGE(PYTEST, pytest)

        AMUSE_DOWNLOAD()
        AC_CHECK_TOOL(PATCH, patch)
        AC_CHECK_TOOL(TAR, tar)
        AC_CHECK_TOOL(UNZIP, unzip)
        AC_CHECK_TOOL(GUNZIP, gunzip)
        AC_CHECK_TOOL(BUNZIP2, bunzip2)
        AC_CHECK_TOOL(UNXZ, unxz)
        AC_CHECK_TOOL(PERL, perl)

        AMUSE_CONDA_PACKAGE(CMAKE, cmake)

        AMUSE_CONDA_PACKAGE(HAVE_MPI, mpi)
        AMUSE_CONDA_PACKAGE(HAVE_OPENMP, openmp)
        AMUSE_CONDA_PACKAGE(CUDA_TK, cuda-toolkit)
        AMUSE_CONDA_PACKAGE(FOUND_OPENCL, ocl-icd)
        if test "$FOUND_OPENCL" == "no"
        then
            # ocl-icd is only for Linux, macOS and Windows use this one
            AMUSE_CONDA_PACKAGE(FOUND_OPENCL, khronos-opencl-icd-loader)
        fi
        AMUSE_CONDA_PACKAGE(HAVE_LAPACK, liblapack)
        AMUSE_CONDA_PACKAGE(BLAS_LIBS, libblas)
        AMUSE_CONDA_PACKAGE(FOUND_GSL, gsl)
        AMUSE_CONDA_PACKAGE(FOUND_GMP, gmp)
        AMUSE_CONDA_PACKAGE(FOUND_MPFR, mpfr)
        AMUSE_CONDA_PACKAGE(FOUND_FFTW, fftw)
        AMUSE_CONDA_PACKAGE(with_hdf5, hdf5)
        AMUSE_CONDA_PACKAGE(FOUND_NETCDF, netcdf4)

        # Export conda list output so setup doesn't have to rerun it to detect pip and wheel
        CONDA_LIST="$amuse_conda_list"
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
        AC_CHECK_TOOL(PYTHON, python3)
        AC_CHECK_TOOL(PYTHON_DEV, python3-config)
        AC_CHECK_TOOL(INSTALL, install)
        AC_CHECK_TOOL(PYTEST, pytest)

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
        AC_CHECK_TOOL(PATCH, patch)
        AC_CHECK_TOOL(TAR, tar)
        AC_CHECK_TOOL(UNZIP, unzip)
        AC_CHECK_TOOL(GUNZIP, gunzip)
        AC_CHECK_TOOL(BUNZIP2, bunzip2)
        AC_CHECK_TOOL(UNXZ, unxz)
        AC_CHECK_TOOL(PERL5, perl)

        AC_CHECK_TOOL(CMAKE, cmake)

        if test "x" != "x$CC"
        then
            # These tests require compiling things, so we need the compiler
            AC_PROG_CC()
            AC_LANG_PUSH([C])

            AX_MPI([HAVE_MPI="yes"])
            AX_OPENMP([HAVE_OPENMP="yes"])
            AMUSE_CUDA()
            AMUSE_OPENCL()

            # Avoid AX_BLAS crashing out due to lack of compilers
            if test "x" != "x$FC"
            then
                AX_LAPACK([HAVE_LAPACK=yes])     # calls AX_BLAS automatically
            else
                AS_ECHO(["Skipping LAPACK/BLAS detection because we don't have Fortran"])
            fi

            AX_PATH_GSL(1.0, [FOUND_GSL=yes])
            AX_GMP()
            AX_MPFR()
            AX_FFTW()
            AX_LIB_HDF5()
            AX_LIB_NETCDF4()

            AC_LANG_POP([C])
        else
            AS_ECHO(["Skipping library detection because we don't have C"])
        fi
    fi

    # This needs to be factored out, because autoconf won't instantiate the compiler
    # checks this calls twice, causing the compiler to appear broken even though it
    # works.
    AMUSE_LIB_SAPPORO_LIGHT()
fi

# Create the FEATURES variable used by amuse.mk to decide what we can build

# This macro adds its first argument to FEATURES if the variable
# in the second argument is defined.
m4_define([ENABLE_FEATURE_IF_DEFINED],
          [AS_IF([test "x" != "x${$2}"], [FEATURES="${FEATURES} $1"])])

# Add the first argument to FEATURES if the variable in the second argument has the
# value specified in the third argument.
m4_define([ENABLE_FEATURE_IF_EQUALS],
          [AS_IF([test "x${$2}" = "x$3"], [FEATURES="${FEATURES} $1"])])

FEATURES="$ASSUME_FEATURES"

ENABLE_FEATURE_IF_DEFINED([c], [CC])
ENABLE_FEATURE_IF_DEFINED([c++], [CXX])
ENABLE_FEATURE_IF_DEFINED([fortran], [FC])
ENABLE_FEATURE_IF_DEFINED([java], [JAVAC])
ENABLE_FEATURE_IF_DEFINED([python], [PYTHON])
ENABLE_FEATURE_IF_DEFINED([python-dev], [PYTHON_DEV])
ENABLE_FEATURE_IF_DEFINED([cmake], [CMAKE])
ENABLE_FEATURE_IF_DEFINED([install], [INSTALL])
ENABLE_FEATURE_IF_DEFINED([download], [DOWNLOAD])
ENABLE_FEATURE_IF_DEFINED([patch], [PATCH])
ENABLE_FEATURE_IF_DEFINED([tar], [TAR])
ENABLE_FEATURE_IF_DEFINED([unzip], [UNZIP])
ENABLE_FEATURE_IF_DEFINED([gunzip], [GUNZIP])
ENABLE_FEATURE_IF_DEFINED([bunzip2], [BUNZIP2])
ENABLE_FEATURE_IF_DEFINED([unxz], [UNXZ])
ENABLE_FEATURE_IF_DEFINED([perl], [PERL])
ENABLE_FEATURE_IF_DEFINED([mpi], [HAVE_MPI])
ENABLE_FEATURE_IF_DEFINED([openmp], [HAVE_OPENMP])
ENABLE_FEATURE_IF_DEFINED([cuda], [CUDA_TK])
ENABLE_FEATURE_IF_EQUALS([opencl], [FOUND_OPENCL], [yes])
ENABLE_FEATURE_IF_DEFINED([blas], [BLAS_LIBS])
ENABLE_FEATURE_IF_DEFINED([lapack], [HAVE_LAPACK])
ENABLE_FEATURE_IF_DEFINED([gsl], [FOUND_GSL])
ENABLE_FEATURE_IF_EQUALS([gmp], [FOUND_GMP], [yes])
ENABLE_FEATURE_IF_EQUALS([mpfr], [FOUND_MPFR], [yes])
ENABLE_FEATURE_IF_EQUALS([fftw], [FOUND_FFTW], [yes])
ENABLE_FEATURE_IF_DEFINED([hdf5], [with_hdf5])
ENABLE_FEATURE_IF_EQUALS([netcdf4], [FOUND_NETCDF], [yes])
ENABLE_FEATURE_IF_EQUALS([sapporo_light], [FOUND_SAPPORO_LIGHT], [yes])

AC_SUBST(PYTEST)

AC_SUBST(APT)
AC_SUBST(DNF)
AC_SUBST(MACPORTS)
AC_SUBST(HOMEBREW)

AC_SUBST(FEATURES)

AC_SUBST(CONDA_LIST)

AC_CONFIG_FILES([configuration.sh])
AC_OUTPUT
