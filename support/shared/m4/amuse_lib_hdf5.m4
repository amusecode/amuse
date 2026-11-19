# AMUSE_LIB_HDF5_CHECK_ROOT()
#
# Basic check to see if the given path contains an HDF5 installation.
#
# Args:
#     source  Source of this path, for printing
#     path    Root of the potential installation tree
#
# Sets amuse_h5_root to the given directory if an installation was found.
#
AC_DEFUN([AMUSE_LIB_HDF5_CHECK_ROOT], [
    AC_MSG_CHECKING([whether $1 points to an HDF5 installation])
    if test -z "$2"
    then
        AC_MSG_RESULT([no])
    elif test -f "$2/bin/h5pcc"
    then
        AC_MSG_RESULT([no, because this seems to be a parallel build])
    else
        if test -f "$2/include/hdf5.h" && test -f "$2/lib/libhdf5.so"
        then
            AC_MSG_RESULT([yes])
            amuse_h5_result="yes"
            amuse_h5_root="$2"
        else
            AC_MSG_RESULT([no])
            amuse_h5_result="no"
        fi
    fi
])


# AMUSE_LIB_HDF5_VERIFY_FLAGS()
#
# Checks that we can compile and link a simple application.
#
# Set a language first using AC_LANG_PUSH, then call this to try compiling and linking.
#
# Args:
#   name    Name of the language/scenario we're testing
#   prefix  Prefix for compiler flags variable, e.g. C or CXX or FC
#   target  Name of target variable to set to yes on success
#   program Program to try to compile and link, use AC_LANG_PROGRAM
#   libs    Libraries to link
#
# Sets amuse_h5_detect_<target> to "yes" or "no" depending on the result.
#
AC_DEFUN([AMUSE_LIB_HDF5_VERIFY_FLAGS], [
    amuse_h5_name="$1"
    amuse_h5_prefix="$2"
    amuse_h5_target="$3"
    amuse_h5_program="$4"
    amuse_h5_libs="$5"

    amuse_h5_flags_var="${amuse_h5_prefix}FLAGS"
    amuse_h5_flags_source_var="amuse_h5_${amuse_h5_flags_var}"
    amuse_h5_flags_target_var="amuse_h5_detect_${amuse_h5_target}"

    eval "amuse_h5_save_FLAGS=\"\$$amuse_h5_flags_var\""
    amuse_h5_save_LDFLAGS="$LDFLAGS"
    amuse_h5_save_LIBS="$LIBS"

    AC_MSG_CHECKING([whether we can build a program with HDF5 in $amuse_h5_name])

    eval "$amuse_h5_flags_var=\"\$$amuse_h5_flags_source_var \$$amuse_h5_flags_var\""
    LDFLAGS="$amuse_h5_LDFLAGS $LDFLAGS"
    LIBS="$amuse_h5_libs"

    dnl that AC_LANG_DEFINES_PROVIDED call silences a spurious warning
    AC_LINK_IFELSE([AC_LANG_DEFINES_PROVIDED()][$amuse_h5_program], [
        eval "$amuse_h5_flags_target_var=\"yes\""
        AC_MSG_RESULT([yes])
    ], [
        eval "$amuse_h5_flags_target_var=\"no\""
        AC_MSG_RESULT([no])
    ])

    LIBS="$amuse_h5_save_LIBS"
    LDFLAGS="$amuse_h5_save_LDFLAGS"
    eval "$amuse_h5_flags_var=\"\$amuse_h5_save_FLAGS\""
])


# AMUSE_LIB_HDF5_VERIFY_FORTRAN()
#
# Checks that we can compile and link with HDF5 in Fortran
#
# See AMUSE_LIB_HDF5_VERIFY_FEATURES below, it's split off from that because we're
# calling it directly in AMUSE_LIB_HDF5_DETECT_ENV.
#
# Args:
#     variant: Text to print describing which variant we're trying
#
# Sets amuse_h5_detect_FC, amuse_h5_detect_HL_FC, and, if amuse_h5_detect_HL_FC is
# yes, sets amuse_h5_HL_FC_LIB to the name of the high-level Fortran library to use.
#
AC_DEFUN([AMUSE_LIB_HDF5_VERIFY_FORTRAN], [
    AC_LANG_PUSH(Fortran)
    AMUSE_LIB_HDF5_VERIFY_FLAGS(
        [Fortran $1], [FC], [FC], [AC_LANG_PROGRAM([], [
            use h5lib

            integer :: error

            call h5open_f(error)
            call h5close_f(error)])],
        [-lhdf5_fortran])

    AMUSE_LIB_HDF5_VERIFY_FLAGS(
        [Fortran with the high-level API in hdf5_hl_fortran $1], [FC], [HL_FC], AC_LANG_PROGRAM([], [
            use h5lt

            integer :: error

            error = h5ltfind_dataset_f(0_hid_t, 'test')
        ]),
        [-lhdf5_hl_fortran])

    if test $amuse_h5_detect_HL_FC = "yes"
    then
        amuse_h5_HL_FC_LIB="-lhdf5_hl_fortran"
    else
        AMUSE_LIB_HDF5_VERIFY_FLAGS(
            [Fortran with the high-level API in hdf5hl_fortran $1], [FC], [HL_FC], AC_LANG_PROGRAM([], [
                use h5lt

                integer :: error

                error = h5ltfind_dataset_f(0_hid_t, 'test')
            ]),
            [-lhdf5hl_fortran])

        if test $amuse_h5_detect_HL_FC = "yes"
        then
            amuse_h5_HL_FC_LIB="-lhdf5hl_fortran"
        fi
    fi

    AC_LANG_POP(Fortran)
])


# AMUSE_LIB_HDF5_DETECT_FEATURES_AT()
#
# Checks that we can use (compile and link) various features of HDF5
#
# Args:
#     path    Root path of the installation to detect features of
#
# Sets amuse_h5_detect_C, amuse_h5_detect_HL_C, amuse_h5_detect_CXX,
# amuse_h5_detect_FC, and amuse_h5_detect_HL_FC to "yes" or "no" indicating whether
# the corresponding features are found to work.
#
# Sets amuse_h5_CFLAGS, amuse_h5_CXXFLAGS, amuse_h5_FCFLAGS, amuse_h5_LDFLAGS to any
# additional (on top of CFLAGS/CXXFLAGS/FCFLAGS/LDFLAGS) flags needed to compile and
# link applications.
#
AC_DEFUN([AMUSE_LIB_HDF5_DETECT_FEATURES_AT], [
    amuse_h5_LDFLAGS="-L$1/lib"

    AC_LANG_PUSH(C)
    amuse_h5_CFLAGS="-I$1/include"
    AMUSE_LIB_HDF5_VERIFY_FLAGS(
        [C], [C], [C], AC_LANG_PROGRAM([#include <hdf5.h>], [H5open(); H5close();]),
        [-lhdf5])

    AMUSE_LIB_HDF5_VERIFY_FLAGS(
        [C with the high-level API], [C], [HL_C],
        AC_LANG_PROGRAM([#include <hdf5_hl.h>], [H5LTfind_dataset(0, \"test\")]),
        [-lhdf5_hl])
    AC_LANG_POP(C)

    AC_LANG_PUSH(C++)
    amuse_h5_CXXFLAGS="-I$1/include"
    AMUSE_LIB_HDF5_VERIFY_FLAGS(
        [C++], [CXX], [CXX],
        AC_LANG_PROGRAM([#include <H5Cpp.h>], [H5::H5Library::open();]),
        [-lhdf5_cpp])
    AC_LANG_POP(C++)

    amuse_h5_FCFLAGS="-I$1/include"
    AMUSE_LIB_HDF5_VERIFY_FORTRAN([with mods in /include])

    if test $amuse_h5_detect_FC = "no"
    then
        amuse_h5_FCFLAGS="-I$1/mod/shared"
        AMUSE_LIB_HDF5_VERIFY_FORTRAN([with mods in /mod/shared])
    fi

    if test $amuse_h5_detect_FC = "no"
    then
        amuse_h5_FCFLAGS=""
    fi
])


# AMUSE_LIB_HDF5_FIND_ROOT()
#
# Finds an HDF5 installation and detects its features.
#
# Tries HDF5_ROOT, HDF5_DIR, and CMAKE_PREFIX_PATH, in that order.
#
# Sets amuse_h5_root to a non-empty value if a working root is found. Sets
# amuse_h5_detect_<feature> to "yes" or "no" indicating whether the feature is present
# and working.
#
AC_DEFUN([AMUSE_LIB_HDF5_FIND_ROOT], [
    amuse_h5_result="no"

    if test -n "$CONDA_PREFIX" 
    then
        AC_MSG_NOTICE([Conda environment detected, ignoring HDF5_ROOT, HDF5_DIR and CMAKE_PREFIX_PATH])
    else
        AMUSE_LIB_HDF5_CHECK_ROOT([HDF5_ROOT], [$HDF5_ROOT])

        if test "$amuse_h5_result" = "no"
        then
            AMUSE_LIB_HDF5_CHECK_ROOT([HDF5_DIR], [$HDF5_DIR])
        fi

        if test "$amuse_h5_result" = "no" && test -n "$CMAKE_PREFIX_PATH"
        then
            amuse_h5_save_IFS="$IFS"
            IFS=:
            for amuse_h5_path in $CMAKE_PREFIX_PATH
            do
                # restore IFS immediately as it may mess up the below
                IFS="$amuse_h5_save_IFS"
                if test -z "$amuse_h5_path"
                then
                    continue
                fi

                AMUSE_LIB_HDF5_CHECK_ROOT(
                        [$amuse_h5_path from CMAKE_PREFIX_PATH], [$amuse_h5_path])

                if test "$amuse_h5_result" = "yes"
                then
                    break
                fi

            done
            IFS="$amuse_h5_save_IFS"
        fi

        if test "$amuse_h5_result" = "yes"
        then
            AMUSE_LIB_HDF5_DETECT_FEATURES_AT([$amuse_h5_root])
        fi
    fi
])


# AMUSE_LIB_HDF5_INCLUDE_PATHS_FROM_FLAGS()
#
# Args:
#   flags: a string containing compiler flags with possible embedded includes
#
# Parses the given string for -isystem <path> and -I<path> flags and sets
# amuse_h5_result to a colon-separated list of those paths.
#
AC_DEFUN([AMUSE_LIB_HDF5_INCLUDE_PATHS_FROM_FLAGS], [
    amuse_h5_result=""
    amuse_h5_in_isystem="no"
    for flag in $1
    do
        if test "$amuse_h5_in_isystem" = "yes"
        then
            if ! echo "$amuse_h5_result" | grep -e "$flag" >/dev/null 2>&1
            then
                amuse_h5_result="$amuse_h5_result:$flag"
            fi
            amuse_h5_in_isystem="no"
            continue
        fi

        if test "$flag" = "-isystem"
        then
            amuse_h5_in_isystem="yes"
        elif test "${flag#-I}" != "$flag"
        then
            amuse_h5_path="${flag#-I}"
            if ! echo "$amuse_h5_result" | grep -e "$amuse_h5_path" >/dev/null 2>&1
            then
                amuse_h5_result="$amuse_h5_result:${flag#-I}"
            fi
        fi
    done
])


# AMUSE_LIB_HDF5_FIX_FORTRAN_USING()
#
# Args:
#   name: description of the source of the paths for printing
#   paths: colon-separated list of include paths to try
#
# Will try to get Fortran to work by adding the given paths to the module search
# directory. Will also try alternative directories (i.e. mod/shared/) based on the given
# paths.
#
AC_DEFUN([AMUSE_LIB_HDF5_FIX_FORTRAN_USING], [
    amuse_h5_save_IFS="$IFS"
    IFS=:
    for amuse_h5_path in $2
    do
        # need to restore IFS immediately as it messes up the below
        IFS="$amuse_h5_save_IFS"
        if test -z "$amuse_h5_path"
        then
            continue
        fi

        amuse_h5_FCFLAGS="$FCFLAGS -I$amuse_h5_path"
        AMUSE_LIB_HDF5_VERIFY_FORTRAN([with $amuse_h5_path from $1])
        if test "$amuse_h5_detect_FC" = "yes"
        then
            amuse_h5_FCFLAGS="-I$amuse_h5_path"
            break
        fi

        # The new Fortran CMake-based build system uses the non-standard and
        # not-searched-unless-explicitly-specified /mod/shared subdirectory by
        # default
        amuse_h5_FCFLAGS="$FCFLAGS -I$amuse_h5_path/../mod/shared"
        AMUSE_LIB_HDF5_VERIFY_FORTRAN([with $amuse_h5_path from $1 and /mod/shared])
        if test "$amuse_h5_detect_FC" = "yes"
        then
            amuse_h5_FCFLAGS="-I$amuse_h5_path/../mod/shared"
            break
        fi

        if test "$amuse_h5_detect_FC" = "no"
        then
            amuse_h5_FCFLAGS=""
        fi
    done
    IFS="$amuse_h5_save_IFS"
])


# AMUSE_LIB_HDF5_DETECT_ENV()
#
# Checks that we can use (compile and link) various features without additional flags
#
# Sets amuse_h5_detect_C, amuse_h5_detect_HL_C, amuse_h5_detect_CXX,
# amuse_h5_detect_FC, and amuse_h5_detect_HL_FC to "yes" or "no" indicating whether
# the corresponding features are found to work.
#
# Sets amuse_h5_FCFLAGS in case additional Fortran flags are needed to work around
# broken installations.
#
AC_DEFUN([AMUSE_LIB_HDF5_DETECT_ENV], [
    amuse_h5_LDFLAGS=""

    AC_LANG_PUSH(C)
    amuse_h5_CFLAGS=""
    AMUSE_LIB_HDF5_VERIFY_FLAGS(
        [C with no additional flags], [C], [C],
        AC_LANG_PROGRAM([#include <hdf5.h>], [H5open(); H5close();]),
        [-lhdf5])
    AC_LANG_POP(C)

    if test $amuse_h5_detect_C = "yes"
    then
        AC_LANG_PUSH(C)
        AMUSE_LIB_HDF5_VERIFY_FLAGS(
            [C with the high-level API], [C], [HL_C],
            AC_LANG_PROGRAM([#include <hdf5_hl.h>], [H5LTfind_dataset(0, \"test\")]),
            [-lhdf5_hl])
        AC_LANG_POP(C)

        AC_LANG_PUSH(C++)
        AMUSE_LIB_HDF5_VERIFY_FLAGS(
            [C++], [CXX], [CXX],
            AC_LANG_PROGRAM([#include <H5Cpp.h>], [H5::H5Library::open();]),
            [-lhdf5_cpp])
        AC_LANG_POP(C++)

        AMUSE_LIB_HDF5_VERIFY_FORTRAN([])

        if test $amuse_h5_detect_FC = "no"
        then
            AMUSE_LIB_HDF5_INCLUDE_PATHS_FROM_FLAGS([$FFLAGS])
            AMUSE_LIB_HDF5_FIX_FORTRAN_USING([FFLAGS], [$amuse_h5_result])
        fi

        if test $amuse_h5_detect_FC = "no"
        then
            AMUSE_LIB_HDF5_INCLUDE_PATHS_FROM_FLAGS([$CFLAGS])
            AMUSE_LIB_HDF5_FIX_FORTRAN_USING([CFLAGS], [$amuse_h5_result])
        fi

        if test $amuse_h5_detect_FC = "no"
        then
            AMUSE_LIB_HDF5_FIX_FORTRAN_USING([CPATH], [$CPATH])
        fi

        if test $amuse_h5_detect_FC = "no"
        then
            AMUSE_LIB_HDF5_FIX_FORTRAN_USING([C_INCLUDE_PATH], [$C_INCLUDE_PATH])
        fi
    fi
])


# AMUSE_LIB_HDF5_IL_FLAGS()
#
# Extracts the -I and -L flags from the given input
#
# Sets amuse_h5_include_flags to a string containing only the -I flags,
# space-separated, and amuse_h5_link_flags to a string containing only the -L flags,
# space-separated.
#
AC_DEFUN([AMUSE_LIB_HDF5_IL_FLAGS], [
    amuse_h5_include_flags=""
    amuse_h5_link_flags=""

    for amuse_h5_flag in $1
    do
        if test "${amuse_h5_flag#-I}" != "$amuse_h5_flag"
        then
            if ! echo "$amuse_h5_include_flags" | grep -e "$amuse_h5_flag" >/dev/null 2>&1
            then
                amuse_h5_include_flags="$amuse_h5_include_flags $amuse_h5_flag"
            fi
        elif test "${amuse_h5_flag#-L}" != "$amuse_h5_flag"
        then
            if ! echo "$amuse_h5_result" | grep -e "$amuse_h5_flag" >/dev/null 2>&1
            then
                amuse_h5_link_flags="$amuse_h5_link_flags $amuse_h5_flag"
            fi
        fi
    done
])


# AMUSE_LIB_HDF5_FROM_WRAPPERS()
#
# Tries to use the h5cc etc. compiler wrappers 
AC_DEFUN([AMUSE_LIB_HDF5_DETECT_FROM_WRAPPERS], [
    AC_PATH_PROG([H5CC], [h5cc])
    if test -n "$H5CC"
    then
        amuse_h5_h5cc_flags=$(HDF5_USE_SHLIB=yes $H5CC -show)
        dnl work around bug in h5cc 1.14.6 where it only shows linking flags unless you
        dnl add -c
        amuse_h5_h5cc_extra_flags=$(HDF5_USE_SHLIB=yes $H5CC -show -c)
        AMUSE_LIB_HDF5_IL_FLAGS([$amuse_h5_h5cc_flags $amuse_h5_h5cc_extra_flags])

        AC_LANG_PUSH(C)
        amuse_h5_CFLAGS="$amuse_h5_include_flags"
        amuse_h5_LDFLAGS="$amuse_h5_link_flags"
        AMUSE_LIB_HDF5_VERIFY_FLAGS(
            [C with flags from h5cc], [C], [C],
            AC_LANG_PROGRAM([#include <hdf5.h>], [H5open(); H5close();]),
            [-lhdf5])
        AMUSE_LIB_HDF5_VERIFY_FLAGS(
            [C with flags from h5cc and the high-level API], [C], [HL_C],
            AC_LANG_PROGRAM([#include <hdf5_hl.h>], [H5LTfind_dataset(0, \"test\")]),
            [-lhdf5_hl])
        AC_LANG_POP(C)
    fi

    AC_PATH_PROG([H5CXX], [h5c++])
    if test -n "$H5CXX"
    then
        amuse_h5_h5cxx_flags=$(HDF5_USE_SHLIB=yes $H5CXX -show)
        dnl Work around a bug(?) in h5cc 1.14.6 where it only shows flags unless you
        dnl add -c _and_ pass a C++ file that exists and has a nonzero length
        dnl Note that conftest.* gets cleaned up automatically by autoconf in case of a
        dnl crash.
        echo "broken" >conftest.cpp
        amuse_h5_h5cxx_extra_flags=$(HDF5_USE_SHLIB=yes $H5CXX -show -c conftest.cpp)
        rm conftest.cpp
        AMUSE_LIB_HDF5_IL_FLAGS([$amuse_h5_h5cxx_flags $amuse_h5_h5cxx_extra_flags])

        AC_LANG_PUSH(C++)
        amuse_h5_CXXFLAGS="$amuse_h5_include_flags"
        amuse_h5_LDFLAGS="$amuse_h5_link_flags"
        AMUSE_LIB_HDF5_VERIFY_FLAGS(
            [C++ with flags from h5c++], [CXX], [CXX],
            AC_LANG_PROGRAM([#include <H5Cpp.h>], [H5::H5Library::open();]),
            [-lhdf5_cpp])
        AC_LANG_POP(C++)
    fi

    AC_PATH_PROG([H5FC], [h5fc])
    if test -n "$H5FC"
    then
        amuse_h5_h5fc_flags=$(HDF5_USE_SHLIB=yes h5fc -show)
        dnl work around bug in h5cc 1.14.6 where it only shows linking flags unless you
        dnl add -c
        amuse_h5_h5fc_extra_flags=$(HDF5_USE_SHLIB=yes $H5CC -show -c)
        AMUSE_LIB_HDF5_IL_FLAGS([$amuse_h5_h5fc_flags $amuse_h5_h5fc_extra_flags])

        amuse_h5_FCFLAGS="$amuse_h5_include_flags"
        amuse_h5_LDFLAGS="$amuse_h5_link_flags"
        AMUSE_LIB_HDF5_VERIFY_FORTRAN([with flags from h5fc])

        if test "$amuse_h5_detect_FC" = "no"
        then
            dnl The extra brackets in the regex get removed by m4
            amuse_h5_extra_flags=$(echo "$amuse_h5_include_flags" | sed -e 's^\([[[:graph:]]]*\)\( \|$\)^\1/../mod/shared ^g')
            amuse_h5_FCFLAGS="$amuse_h5_include_flags $amuse_h5_extra_flags"
            AMUSE_LIB_HDF5_VERIFY_FORTRAN([with flags from h5fc and mod/shared])
        fi

        if test "$amuse_h5_detect_FC" = "no"
        then
            amuse_h5_FCFLAGS=""
        fi
    fi
])



# AMUSE_LIB_HDF5()
#
# Detect HDF5 installation
#
# Sets HDF5_C_FOUND, HDF5_HL_C_FOUND, HDF5_CXX_FOUND, HDF5_FC_FOUND, HDF5_HL_FC_FOUND
# to "yes" or "no" depending on whether the corresponding features are available.
#
# Sets HDF5_TYPE, HDF5_CFLAGS, HDF5_CXXFLAGS, HDF5_FCFLAGS, and HDF5_LDFLAGS to any
# additional flags needed on top of the existing CFLAGS, CXXFLAGS, FCFLAGS and LDFLAGS
# to compile and link with HDF5.
#
# Limitations:
#
# - This only detects serial HDF5 installations, not parallel ones
# - Spaces in paths probably don't work
#
AC_DEFUN([AMUSE_LIB_HDF5], [
    amuse_h5_detect_C="no"

    dnl Detect using HDF5_ROOT, HDF5_DIR, CMAKE_PREFIX_PATH
    AMUSE_LIB_HDF5_FIND_ROOT()

    dnl Check whether the environment is already configured correctly
    if test "$amuse_h5_detect_C" = "no"
    then
        AMUSE_LIB_HDF5_DETECT_ENV()
    fi

    dnl Try to use the h5cc and co compiler wrappers
    if test "$amuse_h5_detect_C" = "no"
    then
        AMUSE_LIB_HDF5_DETECT_FROM_WRAPPERS()
    fi

    if test "$amuse_h5_detect_C" = "yes"
    then
        AC_MSG_NOTICE([HDF5 installation detected])

        HDF5_C_FOUND="$amuse_h5_detect_C"
        HDF5_HL_C_FOUND="$amuse_h5_detect_HL_C"
        HDF5_CXX_FOUND="$amuse_h5_detect_CXX"
        HDF5_FC_FOUND="$amuse_h5_detect_FC"
        HDF5_HL_FC_FOUND="$amuse_h5_detect_HL_FC"

        HDF5_TYPE="serial"
        HDF5_CFLAGS="$amuse_h5_CFLAGS"
        HDF5_CXXFLAGS="$amuse_h5_CXXFLAGS"
        HDF5_FCFLAGS="$amuse_h5_FCFLAGS"
        HDF5_LDFLAGS="$amuse_h5_LDFLAGS"
        HDF5_HLFC_LIB="$amuse_h5_HL_FC_LIB"
    else
        AC_MSG_WARN([no HDF5 installation detected])
    fi

    AC_SUBST(HDF5_TYPE)
    AC_SUBST(HDF5_CFLAGS)
    AC_SUBST(HDF5_CXXFLAGS)
    AC_SUBST(HDF5_FCFLAGS)
    AC_SUBST(HDF5_LDFLAGS)
    AC_SUBST(HDF5_HLFC_LIB)
])

