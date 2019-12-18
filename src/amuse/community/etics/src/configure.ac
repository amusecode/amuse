#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

### TODO handle --without-mpi and --without-hdf5

AC_PREREQ([2.69])
AC_INIT([etics], [1.0.0], [sahmes@gmail.com])

AC_MSG_WARN([This is a rudimentary and experimental configure script!])

AC_CONFIG_SRCDIR([ic.hpp])
AC_CONFIG_HEADERS([config.h])


AS_IF([test -z $CXXFLAGS], [CXXFLAGS="-O3"])
AC_PROG_CXX
AC_PROG_RANLIB

# We are going to politely ignore the GNU standards and explicitly ask the user
# to specify unusual installation locations using the --with commands.

AC_ARG_WITH([a00000],
[
 For the following --with-xxx options, you can specify the root directory of
 these special packages that may be installed in unusual locations.

    --with-xxx=DIR     - Specify only the directory which contains the
                         include/ and lib/ subdirectories
])

AC_ARG_WITH(
    [cuda],
    [AS_HELP_STRING([--with-cuda=DIR], [CUDA installation directory])],
    [CUDA_ROOT="$withval"]
)
AS_IF([test "$CUDA_ROOT" == "no"], [AC_MSG_ERROR([This package needs CUDA])])

AC_ARG_WITH(
    [mpi],
    [AS_HELP_STRING([--with-mpi=DIR], [MPI installation directory])],
    [MPI_ROOT="$withval"]
)

AC_ARG_WITH(
    [boost],
    [AS_HELP_STRING([--with-boost=DIR], [BOOST installation directory])],
    [BOOST_ROOT="$withval"]
)

AC_ARG_WITH(
    [gpu-arch],
    [AS_HELP_STRING([--with-gpu-arch=ARCH], [Specify device architechture [default=35]])],
    [GPU_ARCH="$withval"]
)
AS_IF([test -z $GPU_ARCH], [GPU_ARCH=35])

# Some macros to help the script to detect libraries in unusual locations

AC_DEFUN([AX_LDFLAGS_PUSH], [
    LDFLAGS_temporary_storage=$LDFLAGS
    LDFLAGS="$1 $LDFLAGS"
])

AC_DEFUN([AX_LDFLAGS_POP], [
    LDFLAGS=$LDFLAGS_temporary_storage
])

AC_DEFUN([AX_CPPFLAGS_PUSH], [
    CPPFLAGS_temporary_storage=$CPPFLAGS
    CPPFLAGS="$1 $CPPFLAGS"
])

AC_DEFUN([AX_CPPFLAGS_POP], [
    CPPFLAGS=$CPPFLAGS_temporary_storage
])


# Check CUDA
AS_IF([test -z $CUDA_ROOT],
[
    AC_CHECK_PROG([NVCC], [nvcc], [nvcc], [no])
    AS_IF([test "$NVCC" == "no"], [AC_MSG_ERROR([This package needs CUDA])],
    [
        AC_PROG_SED
        CUDA_ROOT=`which nvcc | sed 's:/bin/nvcc$::'`
    ])
],
[
    AC_CHECK_PROG([NVCC], [nvcc], [$CUDA_ROOT/bin/nvcc], [no], [$CUDA_ROOT/bin])
    AS_IF([test "$NVCC" == "no"], [AC_MSG_ERROR([No nvcc executable was found in $CUDA_ROOT/bin])])
])

AX_LDFLAGS_PUSH([-L$CUDA_ROOT/lib64])
AC_CHECK_LIB([cuda], [main], [], [AC_MSG_ERROR([Error while checking for libcuda])])
AC_CHECK_LIB([cudart], [main], [], [AC_MSG_ERROR([Error while checking for libcudart])])
AC_CHECK_LIB([cudadevrt], [main], [], [AC_MSG_ERROR([Error while checking for cudadevrt])])
AX_LDFLAGS_POP

# Check MPI
# AS_IF([test -z $MPI_ROOT],
# [
#     AC_CHECK_LIB([mpi], [main], [],
#     [
#         AS_UNSET([ac_cv_lib_mpi_main])
#         AC_CHECK_PROG([MPICC], [mpicc], [mpicc], [:])
#         AS_IF([test "$MPICC" == :], [AC_MSG_ERROR([This package needs MPI])],
#             [
#             AC_PROG_SED
#             AC_MSG_NOTICE([MPICC executable found, trying to find library nearby])
#             MPI_ROOT=`which mpicc | sed 's:/bin/mpicc$::'`
#             AX_LDFLAGS_PUSH([-L$MPI_ROOT/lib])
#             AC_CHECK_LIB([mpi], [main], [], [AC_MSG_ERROR([This package needs MPI])])
#             AX_LDFLAGS_POP
#         ])
#     ])
# ],
# [
#     AX_LDFLAGS_PUSH([-L$MPI_ROOT/lib])
#     AC_CHECK_LIB([mpi], [main], [], [AC_MSG_ERROR([No compatible MPI library was found in $MPI_ROOT/lib])])
#     AX_LDFLAGS_POP
# ])


AS_IF([test -z $MPI_ROOT],
[
    AC_CHECK_PROG([MPICXX], [mpic++], [mpic++], [no])
    AS_IF([test "$MPICXX" == "no"], [AC_MSG_ERROR([This package needs MPI (not really but this configure script does)])])
],
[
    AC_CHECK_PROG([MPICXX], [mpic++], [$MPI_ROOT/bin/mpic++], [no], [$MPI_ROOT/bin])
    AS_IF([test "$MPICXX" == "no"], [AC_MSG_ERROR([Cound not find MPI compiler in $MPI_ROOT/bin])])
])

AS_IF([test -z $HDF5_ROOT],
[
    AC_CHECK_PROG([H5CXX], [h5c++], [h5c++], [no])
    AS_IF([test "$H5CXX" == "no"], [AC_MSG_ERROR([This package needs HDF5 (not really but this configure script does)])])
    HDF5_ROOT=`which h5c++ | sed 's:/bin/h5c++$::'`
],
[
    AC_CHECK_PROG([H5CXX], [h5c++], [$HDF5_ROOT/bin/h5c++], [no], [$HDF5_ROOT/bin])
    AS_IF([test "$H5CXX" == "no"], [AC_MSG_ERROR([Cound not find HDF5 compiler in $HDF5_ROOT/bin])])
])


# # Check HDF5
# AS_IF([test -z $HDF5_ROOT],
# [
#     AC_CHECK_LIB([hdf5], [main], [],
#     [
#         AS_UNSET([ac_cv_lib_hdf5_main])
#         AC_CHECK_PROG([HDF5CC], [h5cc], [h5cc], [:])
#         AS_IF([test "$HDF5CC" == :], [AC_MSG_ERROR([This package needs HDF5])],
#             [
#             AC_PROG_SED
#             AC_MSG_NOTICE([H5CC executable found, trying to find library nearby])
#             HDF5_ROOT=`which h5cc | sed 's:/bin/h5cc$::'`
#             AC_MSG_NOTICE([XXXXX $HDF5_ROOT])
# #             AX_LDFLAGS_PUSH([dgsfgdfgd /home/ymeiron/local/lib/libhdf5.a])
# #             AC_CHECK_LIB([hdf5], [main], [], [AC_MSG_ERROR([This package needs HDF5])])
# #             AX_LDFLAGS_POP
#         ])
#     ])
# ],
# [
#     AX_LDFLAGS_PUSH([-L$HDF5_ROOT/lib])
#     AC_CHECK_LIB([hdf5], [main], [], [AC_MSG_ERROR([No compatible HDF5 library was found in $HDF5_ROOT/lib])])
#     AX_LDFLAGS_POP
# ])



H5CXX_SHOW=`h5c++ -show`
MPICXX_SHOW=`mpic++ -show`

MPICXX_COMPILER=`echo $MPICXX_SHOW | awk '{print $1;}'`
H5CXX_COMPILER=`echo $H5CXX_SHOW | awk '{print $1;}'`

H5CXX_SHOW=`echo "$H5CXX_SHOW" | cut -d " " -f2-`
MPICXX_SHOW=`echo "$MPICXX_SHOW" | cut -d " " -f2-`

AS_IF([test $MPICXX_COMPILER == $H5CXX_COMPILER], [], [AC_MSG_ERROR([mpic++ and h5c++ seem to wrap different compilers!])])
AS_IF([test $MPICXX_COMPILER == $CXX], [], [AC_MSG_ERROR([mpic++ and h5c++ seem to wrap a different compiler than the requested $CXX])])

# MPICXX_SHOW_UNIQUE=""
# for OPTION in $MPICXX_SHOW
# do
#     if test "${H5CXX_SHOW/$OPTION}" == "$H5CXX_SHOW"; then
#         MPICXX_SHOW_UNIQUE="$MPICXX_SHOW_UNIQUE $OPTION"
#     fi
# done
MPICXX_SHOW_UNIQUE=$MPICXX_SHOW

CPPFLAGS=`echo $H5CXX_SHOW $MPICXX_SHOW_UNIQUE | grep -o '\-[[DI]]\S\+' | tr '\n' ' ' | sed -e 's/\ \+/\ /g'`" "$CPPFLAGS
LDFLAGS=`echo $H5CXX_SHOW $MPICXX_SHOW_UNIQUE | grep -o '\-L\S\+\|\-Wl\S\+' | tr '\n' ' ' | sed -e 's/\ \+/\ /g'`" "$LDFLAGS
LIBS=`echo $H5CXX_SHOW $MPICXX_SHOW_UNIQUE | grep -o '\-l\S\+' | tr '\n' ' ' | sed -e 's/\ \+/\ /g'`" "$LIBS
STATICLIBS=`echo $H5CXX_SHOW $MPICXX_SHOW_UNIQUE | grep -o '\s[[^-]].\+' | tr '\n' ' ' | sed -e 's/\ \+/\ /g'`" "$STATICLIBS

KNOWNFLAGS="$CPPFLAGS $LDFLAGS $LIBS $STATICLIBS"
ALLOPTIONS="$H5CXX_SHOW $MPICXX_SHOW_UNIQUE"
for OPTION in $ALLOPTIONS
do
    if test "${KNOWNFLAGS/$OPTION}" == "$KNOWNFLAGS"; then
        CXXFLAGS="$CXXFLAGS $OPTION"
    fi
done


# echo -e "$CXX \n $CPPFLAGS $LDFLAGS $STATICLIBS $LIBS $CXXFLAGS" | tr '\n' ' ' | sed -e 's/\ \+/\ /g'


# Check BOOST
CONFIGSCRIPT_HAVE_BOOST="no"
AC_LANG_PUSH([C++])
AS_IF([test -z $BOOST_ROOT],
[
    AC_CHECK_HEADER([boost/version.hpp],
    [
        AC_DEFINE([HAVE_BOOST], [1], [Define to 1 if you have the boost libraries.])
        CONFIGSCRIPT_HAVE_BOOST="yes"
    ])
],
[
    AS_IF([test "$BOOST_ROOT" != "no"],
    [
        AX_CPPFLAGS_PUSH([-I$BOOST_ROOT])
        AC_CHECK_HEADER([boost/version.hpp],
        [
            AC_DEFINE([HAVE_BOOST], [1], [Define to 1 if you have the boost libraries.])
            CONFIGSCRIPT_HAVE_BOOST="yes"
        ],
        [
            AC_MSG_ERROR([Cound not find BOOST headers in $BOOST_ROOT])
        ])
        AX_CPPFLAGS_POP
    ])
])
AC_LANG_POP([C++])

AS_IF([test "$CONFIGSCRIPT_HAVE_BOOST" == "yes"], [BOOSTINC=$BOOST_ROOT], [BOOSTINC="NOBOOST"])


AC_SUBST([CUDA_ROOT])
AC_SUBST([NVCC])
AC_SUBST([GPU_ARCH])
AC_SUBST([STATICLIBS])
AC_SUBST([HDF5_ROOT])
AC_SUBST([BOOSTINC])
AC_CONFIG_FILES([Makefile])


echo "CPPFLAGS=$CPPFLAGS"
echo "CXXFLAGS=$CXXFLAGS"
echo "LDFLAGS=$LDFLAGS"
echo "LIBS=$LIBS"
echo "STATICLIBS=$STATICLIBS"



AC_MSG_NOTICE([printing summary

Compilers:
----------
                   C++ compiler: $CXX
                  CUDA compiler: $NVCC

Libraries:
----------
                          boost: $CONFIGSCRIPT_HAVE_BOOST
                           cuda: yes
                           hdf5: yes
                            mpi: yes

Parameters:
-----------
              GPU architechture: $GPU_ARCH
])

AC_OUTPUT
