# generated automatically by aclocal 1.16.3 -*- Autoconf -*-

# Copyright (C) 1996-2020 Free Software Foundation, Inc.

# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

m4_ifndef([AC_CONFIG_MACRO_DIRS], [m4_defun([_AM_CONFIG_MACRO_DIRS], [])m4_defun([AC_CONFIG_MACRO_DIRS], [_AM_CONFIG_MACRO_DIRS($@)])])
# ===========================================================================
#       https://www.gnu.org/software/autoconf-archive/ax_lib_hdf5.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_HDF5([serial/parallel])
#
# DESCRIPTION
#
#   This macro provides tests of the availability of HDF5 library.
#
#   The optional macro argument should be either 'serial' or 'parallel'. The
#   former only looks for serial HDF5 installations via h5cc. The latter
#   only looks for parallel HDF5 installations via h5pcc. If the optional
#   argument is omitted, serial installations will be preferred over
#   parallel ones.
#
#   The macro adds a --with-hdf5 option accepting one of three values:
#
#     no   - do not check for the HDF5 library.
#     yes  - do check for HDF5 library in standard locations.
#     path - complete path to the HDF5 helper script h5cc or h5pcc.
#
#   If HDF5 is successfully found, this macro calls
#
#     AC_SUBST(HDF5_VERSION)
#     AC_SUBST(HDF5_CC)
#     AC_SUBST(HDF5_CFLAGS)
#     AC_SUBST(HDF5_CPPFLAGS)
#     AC_SUBST(HDF5_LDFLAGS)
#     AC_SUBST(HDF5_LIBS)
#     AC_SUBST(HDF5_FC)
#     AC_SUBST(HDF5_FFLAGS)
#     AC_SUBST(HDF5_FLIBS)
#     AC_SUBST(HDF5_TYPE)
#     AC_DEFINE(HAVE_HDF5)
#
#   and sets with_hdf5="yes".  Additionally, the macro sets
#   with_hdf5_fortran="yes" if a matching Fortran wrapper script is found.
#   Note that Autoconf's Fortran support is not used to perform this check.
#   H5CC and H5FC will contain the appropriate serial or parallel HDF5
#   wrapper script locations.
#
#   If HDF5 is disabled or not found, this macros sets with_hdf5="no" and
#   with_hdf5_fortran="no".
#
#   Your configuration script can test $with_hdf to take any further
#   actions. HDF5_{C,CPP,LD}FLAGS may be used when building with C or C++.
#   HDF5_F{FLAGS,LIBS} should be used when building Fortran applications.
#
#   To use the macro, one would code one of the following in "configure.ac"
#   before AC_OUTPUT:
#
#     1) dnl Check for HDF5 support
#        AX_LIB_HDF5()
#
#     2) dnl Check for serial HDF5 support
#        AX_LIB_HDF5([serial])
#
#     3) dnl Check for parallel HDF5 support
#        AX_LIB_HDF5([parallel])
#
#   One could test $with_hdf5 for the outcome or display it as follows
#
#     echo "HDF5 support:  $with_hdf5"
#
#   You could also for example, override the default CC in "configure.ac" to
#   enforce compilation with the compiler that HDF5 uses:
#
#     AX_LIB_HDF5([parallel])
#     if test "$with_hdf5" = "yes"; then
#             CC="$HDF5_CC"
#     else
#             AC_MSG_ERROR([Unable to find HDF5, we need parallel HDF5.])
#     fi
#
#   The HDF5_TYPE environment variable returns "parallel" or "serial",
#   depending on which type of library is found.
#
# LICENSE
#
#   Copyright (c) 2009 Timothy Brown <tbrown@freeshell.org>
#   Copyright (c) 2010 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 20

AC_DEFUN([AX_LIB_HDF5], [

AC_REQUIRE([AC_PROG_SED])
AC_REQUIRE([AC_PROG_AWK])
AC_REQUIRE([AC_PROG_GREP])

dnl Check first argument is one of the recognized values.
dnl Fail eagerly if is incorrect as this simplifies case statements below.
if   test "m4_normalize(m4_default([$1],[]))" = ""        ; then
    : # Recognized value
elif test "m4_normalize(m4_default([$1],[]))" = "serial"  ; then
    : # Recognized value
elif test "m4_normalize(m4_default([$1],[]))" = "parallel"; then
    : # Recognized value
else
    AC_MSG_ERROR([
Unrecognized value for AX[]_LIB_HDF5 within configure.ac.
If supplied, argument 1 must be either 'serial' or 'parallel'.
])
fi

dnl Add a default --with-hdf5 configuration option.
AC_ARG_WITH([hdf5],
  AS_HELP_STRING(
    [--with-hdf5=[yes/no/PATH]],
    m4_case(m4_normalize([$1]),
            [serial],   [location of h5cc for serial HDF5 configuration],
            [parallel], [location of h5pcc for parallel HDF5 configuration],
            [location of h5cc or h5pcc for HDF5 configuration])
  ),
  [if test "$withval" = "no"; then
     with_hdf5="no"
   elif test "$withval" = "yes"; then
     with_hdf5="yes"
   else
     with_hdf5="yes"
     H5CC="$withval"
   fi],
   [with_hdf5="yes"]
)

dnl Set defaults to blank
HDF5_CC=""
HDF5_VERSION=""
HDF5_CFLAGS=""
HDF5_CPPFLAGS=""
HDF5_LDFLAGS=""
HDF5_LIBS=""
HDF5_FC=""
HDF5_FFLAGS=""
HDF5_FLIBS=""
HDF5_TYPE=""

dnl Try and find hdf5 compiler tools and options.
if test "$with_hdf5" = "yes"; then
    if test -z "$H5CC"; then
        dnl Check to see if H5CC is in the path.
        AC_PATH_PROGS(
            [H5CC],
            m4_case(m4_normalize([$1]),
                [serial],   [h5cc],
                [parallel], [h5pcc],
                [h5cc h5pcc]),
            [])
    else
        AC_MSG_CHECKING([Using provided HDF5 C wrapper])
        AC_MSG_RESULT([$H5CC])
    fi
    AC_MSG_CHECKING([for HDF5 type])
    AS_CASE([$H5CC],
        [*h5pcc], [HDF5_TYPE=parallel],
        [*h5cc], [HDF5_TYPE=serial],
        [HDF5_TYPE=neither])
    AC_MSG_RESULT([$HDF5_TYPE])
    AC_MSG_CHECKING([for HDF5 libraries])
    if test ! -f "$H5CC" || test ! -x "$H5CC"; then
        AC_MSG_RESULT([no])
        AC_MSG_WARN(m4_case(m4_normalize([$1]),
            [serial],  [
Unable to locate serial HDF5 compilation helper script 'h5cc'.
Please specify --with-hdf5=<LOCATION> as the full path to h5cc.
HDF5 support is being disabled (equivalent to --with-hdf5=no).
],            [parallel],[
Unable to locate parallel HDF5 compilation helper script 'h5pcc'.
Please specify --with-hdf5=<LOCATION> as the full path to h5pcc.
HDF5 support is being disabled (equivalent to --with-hdf5=no).
],            [
Unable to locate HDF5 compilation helper scripts 'h5cc' or 'h5pcc'.
Please specify --with-hdf5=<LOCATION> as the full path to h5cc or h5pcc.
HDF5 support is being disabled (equivalent to --with-hdf5=no).
]))
        with_hdf5="no"
        with_hdf5_fortran="no"
    else
        dnl Get the h5cc output
        HDF5_SHOW=$(eval $H5CC -show)

        dnl Get the actual compiler used
        HDF5_CC=$(eval $H5CC -show | head -n 1 | $AWK '{print $[]1}')
        if test "$HDF5_CC" = "ccache"; then
            HDF5_CC=$(eval $H5CC -show | head -n 1 | $AWK '{print $[]2}')
        fi

        dnl h5cc provides both AM_ and non-AM_ options
        dnl depending on how it was compiled either one of
        dnl these are empty. Lets roll them both into one.

        dnl Look for "HDF5 Version: X.Y.Z"
        HDF5_VERSION=$(eval $H5CC -showconfig | $GREP 'HDF5 Version:' \
            | $AWK '{print $[]3}')

        dnl A ideal situation would be where everything we needed was
        dnl in the AM_* variables. However most systems are not like this
        dnl and seem to have the values in the non-AM variables.
        dnl
        dnl We try the following to find the flags:
        dnl (1) Look for "NAME:" tags
        dnl (2) Look for "H5_NAME:" tags
        dnl (3) Look for "AM_NAME:" tags
        dnl
        HDF5_tmp_flags=$(eval $H5CC -showconfig \
            | $GREP 'FLAGS\|Extra libraries:' \
            | $AWK -F: '{printf("%s "), $[]2}' )

        dnl Find the installation directory and append include/
        HDF5_tmp_inst=$(eval $H5CC -showconfig \
            | $GREP 'Installation point:' \
            | $AWK '{print $[]NF}' )

        dnl Add this to the CPPFLAGS
        HDF5_CPPFLAGS="-I${HDF5_tmp_inst}/include"

        dnl Now sort the flags out based upon their prefixes
        for arg in $HDF5_SHOW $HDF5_tmp_flags ; do
          case "$arg" in
            -I*) echo $HDF5_CPPFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || HDF5_CPPFLAGS="$HDF5_CPPFLAGS $arg"
              ;;
            -L*) echo $HDF5_LDFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || HDF5_LDFLAGS="$HDF5_LDFLAGS $arg"
              ;;
            -l*) echo $HDF5_LIBS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || HDF5_LIBS="$HDF5_LIBS $arg"
              ;;
          esac
        done

        HDF5_LIBS="-lhdf5 $HDF5_LIBS"
        AC_MSG_RESULT([yes (version $[HDF5_VERSION])])

        dnl See if we can compile
        AC_LANG_PUSH([C])
        ax_lib_hdf5_save_CC=$CC
        ax_lib_hdf5_save_CPPFLAGS=$CPPFLAGS
        ax_lib_hdf5_save_LIBS=$LIBS
        ax_lib_hdf5_save_LDFLAGS=$LDFLAGS
        CC=$HDF5_CC
        CPPFLAGS=$HDF5_CPPFLAGS
        LIBS=$HDF5_LIBS
        LDFLAGS=$HDF5_LDFLAGS
        AC_CHECK_HEADER([hdf5.h], [ac_cv_hadf5_h=yes], [ac_cv_hadf5_h=no])
        AC_CHECK_LIB([hdf5], [H5Fcreate], [ac_cv_libhdf5=yes],
                     [ac_cv_libhdf5=no])
        if test "$ac_cv_hadf5_h" = "no" && test "$ac_cv_libhdf5" = "no" ; then
          AC_MSG_WARN([Unable to compile HDF5 test program])
        fi
        dnl Look for HDF5's high level library
        AC_HAVE_LIBRARY([hdf5_hl], [HDF5_LIBS="-lhdf5_hl $HDF5_LIBS"], [], [])

        CC=$ax_lib_hdf5_save_CC
        CPPFLAGS=$ax_lib_hdf5_save_CPPFLAGS
        LIBS=$ax_lib_hdf5_save_LIBS
        LDFLAGS=$ax_lib_hdf5_save_LDFLAGS
        AC_LANG_POP([C])

        AC_MSG_CHECKING([for matching HDF5 Fortran wrapper])
        dnl Presume HDF5 Fortran wrapper is just a name variant from H5CC
        H5FC=$(eval echo -n $H5CC | $SED -n 's/cc$/fc/p')
        if test -x "$H5FC"; then
            AC_MSG_RESULT([$H5FC])
            with_hdf5_fortran="yes"
            AC_SUBST([H5FC])

            dnl Again, pry any remaining -Idir/-Ldir from compiler wrapper
            for arg in `$H5FC -show`
            do
              case "$arg" in #(
                -I*) echo $HDF5_FFLAGS | $GREP -e "$arg" >/dev/null \
                      || HDF5_FFLAGS="$HDF5_FFLAGS $arg"
                  ;;#(
                -L*) echo $HDF5_FFLAGS | $GREP -e "$arg" >/dev/null \
                      || HDF5_FFLAGS="$HDF5_FFLAGS $arg"
                     dnl HDF5 installs .mod files in with libraries,
                     dnl but some compilers need to find them with -I
                     echo $HDF5_FFLAGS | $GREP -e "-I${arg#-L}" >/dev/null \
                      || HDF5_FFLAGS="$HDF5_FFLAGS -I${arg#-L}"
                  ;;
              esac
            done

            dnl Make Fortran link line by inserting Fortran libraries
            for arg in $HDF5_LIBS
            do
              case "$arg" in #(
                -lhdf5_hl) HDF5_FLIBS="$HDF5_FLIBS -lhdf5hl_fortran $arg"
                  ;; #(
                -lhdf5)    HDF5_FLIBS="$HDF5_FLIBS -lhdf5_fortran $arg"
                  ;; #(
                *) HDF5_FLIBS="$HDF5_FLIBS $arg"
                  ;;
              esac
            done
        else
            AC_MSG_RESULT([no])
            with_hdf5_fortran="no"
        fi

	AC_SUBST([HDF5_VERSION])
	AC_SUBST([HDF5_CC])
	AC_SUBST([HDF5_CFLAGS])
	AC_SUBST([HDF5_CPPFLAGS])
	AC_SUBST([HDF5_LDFLAGS])
	AC_SUBST([HDF5_LIBS])
	AC_SUBST([HDF5_FC])
	AC_SUBST([HDF5_FFLAGS])
	AC_SUBST([HDF5_FLIBS])
	AC_SUBST([HDF5_TYPE])
	AC_DEFINE([HAVE_HDF5], [1], [Defined if you have HDF5 support])
    fi
fi
])

# ===========================================================================
#      https://www.gnu.org/software/autoconf-archive/ax_lib_netcdf4.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_NETCDF4([serial/parallel])
#
# DESCRIPTION
#
#   This macro provides tests of the availability of the NetCDF v4 library.
#
#   The optional macro argument should be either 'serial' or 'parallel'. The
#   macro will call nc-config to check the output of the '--has-pnetcdf'
#   option and error out if the requested parallel isn't supported.
#
#   If the optional argument is omitted, no check is made to see if NetCDF
#   has parallel support.
#
#   The macro adds a --with-netcdf4 option accepting one of three values:
#
#     no   - do not check for the NetCDF4 library.
#     yes  - do check for NetCDF4 library in standard locations.
#     path - installation prefix for NetCDF version 4.
#
#   If NetCDF4 is successfully found, this macro calls
#
#     AC_SUBST(NETCDF4_VERSION)
#     AC_SUBST(NETCDF4_CC)
#     AC_SUBST(NETCDF4_CFLAGS)
#     AC_SUBST(NETCDF4_CPPFLAGS)
#     AC_SUBST(NETCDF4_LDFLAGS)
#     AC_SUBST(NETCDF4_LIBS)
#     AC_SUBST(NETCDF4_FC)
#     AC_SUBST(NETCDF4_FFLAGS)
#     AC_SUBST(NETCDF4_FLIBS)
#     AC_DEFINE(HAVE_NETCDF4)
#
#   It also sets
#
#     with_netcdf4="yes"
#     with_netcdf4_fortran="yes"    (if NetCDF has Fortran support)
#     with_netcdf4_parallel="yes"   (if NetCDF has MPI support)
#
#   If NetCDF4 is disabled or not found, this macros sets
#
#     with_netcdf4="no"
#     with_netcdf4_fortran="no"
#
#   Note it does not set with_netcdf4_parallel in this case.
#
#   Your configuration script can test $with_netcdf4 to take any further
#   actions. NETCDF4_{C,CPP,LD}FLAGS may be used when building with C or
#   C++. NETCDF4_F{FLAGS,LIBS} and NETCDF4_LDFLAGS should be used when
#   building Fortran applications.
#
#   To use the macro, one would code one of the following in "configure.ac"
#   before AC_OUTPUT:
#
#     1) dnl Check for NetCDF4 support
#        AX_LIB_NETCDF4()
#
#     2) dnl Check for serial NetCDF4 support
#        AX_LIB_NETCDF4([serial])
#
#     3) dnl Check for parallel NetCDF4 support
#        AX_LIB_NETCDF4([parallel])
#
#   One could test $with_netcdf4 for the outcome or display it as follows
#
#     echo "NetCDF v4 support:  $with_netcdf4"
#
#   One could also for example, override the default CC in "configure.ac" to
#   enforce compilation with the compiler that NetCDF v4 was built with:
#
#     AX_LIB_NETCDF4([parallel])
#     if test "$with_netcdf4" = "yes"; then
#             CC="$NETCDF4_CC"
#     else
#             AC_MSG_ERROR([Unable to find NetCDF4, we need parallel NetCDF4.])
#     fi
#
# LICENSE
#
#   Copyright (c) 2016 Timothy Brown <tbrown@freeshell.org>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 2

AC_DEFUN([AX_LIB_NETCDF4], [

AC_REQUIRE([AC_PROG_SED])
AC_REQUIRE([AC_PROG_AWK])
AC_REQUIRE([AC_PROG_GREP])

dnl Check first argument is one of the recognized values.
dnl Fail eagerly if is incorrect as this simplifies case statements below.
if   test "m4_normalize(m4_default([$1],[]))" = ""        ; then
    netcdf4_requested_mode="serial"
elif test "m4_normalize(m4_default([$1],[]))" = "serial"  ; then
    netcdf4_requested_mode="serial"
elif test "m4_normalize(m4_default([$1],[]))" = "parallel"; then
    netcdf4_requested_mode="parallel"
else
    AC_MSG_ERROR([
Unrecognized value for AX[]_LIB_NETCDF4 within configure.ac.
If supplied, argument 1 must be either 'serial' or 'parallel'.
])
fi

dnl Add a default --with-netcdf4 configuration option.
AC_ARG_WITH([netcdf4],
  AS_HELP_STRING(
    [--with-netcdf4=[yes/no/PATH]],
    m4_case(m4_normalize([$1]),
            [serial],   [base directory of serial NetCDF4 installation],
            [parallel], [base directory of parallel NetCDF4 installation],
            [base directory of NetCDF4 installation])
  ),
  [if test "$withval" = "no"; then
     with_netcdf4="no"
   elif test "$withval" = "yes"; then
     with_netcdf4="yes"
   else
     with_netcdf4="yes"
     NETCDF4_PREFIX="${withval}"
     NC_CONFIG="${withval}/bin/nc-config"
   fi],
   [with_netcdf4="yes"]
)

dnl Set defaults to blank
NETCDF4_CC=""
NETCDF4_VERSION=""
NETCDF4_CFLAGS=""
NETCDF4_CPPFLAGS=""
NETCDF4_LDFLAGS=""
NETCDF4_LIBS=""
NETCDF4_FC=""
NETCDF4_FFLAGS=""
NETCDF4_FLIBS=""

dnl Try and find NetCDF4 tools and options.
if test "$with_netcdf4" = "yes"; then
    if test -z "$NC_CONFIG"; then
        dnl Check to see if NC_CONFIG is in the path.
        AC_PATH_PROGS([NC_CONFIG], [nc-config], [])
        NETCDF4_PREFIX=$(AS_DIRNAME([$(AS_DIRNAME(["$NC_CONFIG"]))]))
    else
        AC_MSG_CHECKING([Using provided NetCDF4 prefix])
        AC_MSG_RESULT([$NC_CONFIG])
    fi

    AC_MSG_CHECKING([for NetCDF4 libraries])

    if test ! -f "$NC_CONFIG" || test ! -x "$NC_CONFIG"; then
        AC_MSG_RESULT([no])
        AC_MSG_WARN([

Unable to locate NetCDF4 compilation helper script 'nc-config'.
Please specify --with-netcdf4=<LOCATION> as the full path prefix
where NetCDF4 has been installed.
NetCDF4 support is being disabled (equivalent to --with-netcdf4=no).
])
        with_netcdf4="no"
        with_netcdf4_fortran="no"
    else
        dnl Get the actual compiler used
        NETCDF4_CC=$(eval $NC_CONFIG --cc | $AWK '{print $[]1}')
        if test "$NETCDF4_CC" = "ccache"; then
            NETCDF4_CC=$(eval $NC_CONFIG --cc | $AWK '{print $[]2}')
        fi

        dnl Look for version
        NETCDF4_VERSION=$(eval $NC_CONFIG --version | $AWK '{print $[]2}')

        dnl Look for the CFLAGS
        NETCDF4_CFLAGS=$(eval $NC_CONFIG --cflags)

        dnl Look for the LIBS and LDFLAGS
        NETCDF4_tmp_clibs=$(eval $NC_CONFIG --libs)

        dnl Sort out the tmp libs based on their prefixes
        for arg in $NETCDF4_tmp_clibs ; do
          case "$arg" in
            -L*) echo $NETCDF4_LDFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                   || NETCDF4_LDFLAGS="$arg $NETCDF4_LDFLAGS"
              ;;
            -l*) echo $NETCDF4_LIBS | $GREP -e "$arg" 2>&1 >/dev/null \
                   || NETCDF4_LIBS="$arg $NETCDF4_LIBS"
              ;;
          esac
        done

        AC_MSG_RESULT([yes (version $[NETCDF4_VERSION])])

        dnl See if we need (and have) parallel support
        if test "$netcdf4_requested_mode" = "parallel" ; then
            with_netcdf4_parallel=$(eval $NC_CONFIG --has-pnetcdf)
            if test "$with_netcdf4_parallel" = "no" ; then
                AC_MSG_ERROR([
parallel NetCDF4 is not supported (while it was requested)
])
            fi
        fi

        dnl See if we can compile
        ax_lib_netcdf4_save_CC=$CC
        ax_lib_netcdf4_save_CPPFLAGS=$CPPFLAGS
        ax_lib_netcdf4_save_LIBS=$LIBS
        ax_lib_netcdf4_save_LDFLAGS=$LDFLAGS
        CC=$NETCDF4_CC
        CFLAGS=$NETCDF4_CFLAGS
        LIBS=$NETCDF4_LIBS
        LDFLAGS=$NETCDF4_LDFLAGS
        AC_CHECK_HEADER([netcdf.h], [ac_cv_netcdf4_h=yes], [ac_cv_netcdf4_h=no])
        AC_CHECK_LIB([netcdf], [nc_create], [ac_cv_libnetcdf4=yes],
                     [ac_cv_libnetcdf4=no])
        if test "$ac_cv_netcdf4_h" = "no" && \
           test "$ac_cv_libnetcdf4" = "no" ; then
            AC_MSG_WARN([Unable to compile NetCDF4 test program])
        fi

        CC=$ax_lib_netcdf4_save_CC
        CFLAGS=$ax_lib_netcdf4_save_CFLAGS
        LIBS=$ax_lib_netcdf4_save_LIBS
        LDFLAGS=$ax_lib_hdf5_save_LDFLAGS


        AC_MSG_CHECKING([for matching NetCDF4 Fortran libraries])
        NF_CONFIG="${NETCDF4_PREFIX}/bin/nf-config"
        if test ! -f "$NF_CONFIG" || test ! -x "$NF_CONFIG"; then
            AC_MSG_RESULT([no])
            with_netcdf4_fortran="no"
        else
            NETCDF_FVERSION=$(eval $NF_CONFIG --version | $AWK '{print $[]2}')
            AC_MSG_RESULT([yes (version $[NETCDF_FVERSION])])
            NETCDF4_FC=$(eval $NF_CONFIG --fc | $AWK '{print $[]1}')
            if test "$NETCDF4_FC" = "ccache"; then
                NETCDF4_FC=$(eval $NF_CONFIG --fc | $AWK '{print $[]2}')
            fi
            dnl Look for the FFLAGS
            NETCDF4_FFLAGS=$(eval $NC_CONFIG --fflags)

            dnl Look for the FLIBS and LDFLAGS
            NETCDF4_tmp_flibs=$(eval $NC_CONFIG --flibs)

            dnl Sort out the tmp libs based on their prefixes
            for arg in $NETCDF4_tmp_flibs ; do
              case "$arg" in
                -L*) echo $NETCDF4_LDFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                       || NETCDF4_LDFLAGS="$arg $NETCDF4_LDFLAGS"
                  ;;
                -l*) echo $NETCDF4_FLIBS | $GREP -e "$arg" 2>&1 >/dev/null \
                       || NETCDF4_FLIBS="$arg $NETCDF4_FLIBS"
                  ;;
               esac
            done
            with_netcdf4_fortran="yes"
        fi

        AC_SUBST([NETCDF4_VERSION])
        AC_SUBST([NETCDF4_CC])
        AC_SUBST([NETCDF4_CFLAGS])
        AC_SUBST([NETCDF4_LDFLAGS])
        AC_SUBST([NETCDF4_LIBS])
        AC_SUBST([NETCDF4_FC])
        AC_SUBST([NETCDF4_FFLAGS])
        AC_SUBST([NETCDF4_FLIBS])
        AC_DEFINE([HAVE_NETCDF4], [1], [Defined if you have NETCDF4 support])
    fi
fi
])

# ===========================================================================
#          https://www.gnu.org/software/autoconf-archive/ax_mpi.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_MPI([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro tries to find out how to compile programs that use MPI
#   (Message Passing Interface), a standard API for parallel process
#   communication (see http://www-unix.mcs.anl.gov/mpi/)
#
#   On success, it sets the MPICC, MPICXX, MPIF77, or MPIFC output variable
#   to the name of the MPI compiler, depending upon the current language.
#   (This may just be $CC/$CXX/$F77/$FC, but is more often something like
#   mpicc/mpiCC/mpif77/mpif90.) It also sets MPILIBS to any libraries that
#   are needed for linking MPI (e.g. -lmpi or -lfmpi, if a special
#   MPICC/MPICXX/MPIF77/MPIFC was not found).
#
#   Note that this macro should be used only if you just have a few source
#   files that need to be compiled using MPI. In particular, you should
#   neither overwrite CC/CXX/F77/FC with the values of
#   MPICC/MPICXX/MPIF77/MPIFC, nor assume that you can use the same flags
#   etc. as the standard compilers. If you want to compile a whole program
#   using the MPI compiler commands, use one of the macros
#   AX_PROG_{CC,CXX,FC}_MPI.
#
#   ACTION-IF-FOUND is a list of shell commands to run if an MPI library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run if it is not
#   found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_MPI.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Julian C. Cummings <cummings@cacr.caltech.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 9

AU_ALIAS([ACX_MPI], [AX_MPI])
AC_DEFUN([AX_MPI], [
AC_PREREQ(2.50) dnl for AC_LANG_CASE

AC_LANG_CASE([C], [
	AC_REQUIRE([AC_PROG_CC])
	AC_ARG_VAR(MPICC,[MPI C compiler command])
	AC_CHECK_PROGS(MPICC, mpicc hcc mpxlc_r mpxlc mpcc cmpicc, $CC)
	ax_mpi_save_CC="$CC"
	CC="$MPICC"
	AC_SUBST(MPICC)
],
[C++], [
	AC_REQUIRE([AC_PROG_CXX])
	AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
	AC_CHECK_PROGS(MPICXX, mpic++ mpicxx mpiCC hcp mpxlC_r mpxlC mpCC cmpic++, $CXX)
	ax_mpi_save_CXX="$CXX"
	CXX="$MPICXX"
	AC_SUBST(MPICXX)
],
[Fortran 77], [
	AC_REQUIRE([AC_PROG_F77])
	AC_ARG_VAR(MPIF77,[MPI Fortran 77 compiler command])
	AC_CHECK_PROGS(MPIF77, mpif77 hf77 mpxlf_r mpxlf mpf77 cmpifc, $F77)
	ax_mpi_save_F77="$F77"
	F77="$MPIF77"
	AC_SUBST(MPIF77)
],
[Fortran], [
	AC_REQUIRE([AC_PROG_FC])
	AC_ARG_VAR(MPIFC,[MPI Fortran compiler command])
	AC_CHECK_PROGS(MPIFC, mpif90 mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpf90 cmpif90c, $FC)
	ax_mpi_save_FC="$FC"
	FC="$MPIFC"
	AC_SUBST(MPIFC)
])

if test x = x"$MPILIBS"; then
	AC_LANG_CASE([C], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
		[C++], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
		[Fortran 77], [AC_MSG_CHECKING([for MPI_Init])
			AC_LINK_IFELSE([AC_LANG_PROGRAM([],[      call MPI_Init])],[MPILIBS=" "
				AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])],
		[Fortran], [AC_MSG_CHECKING([for MPI_Init])
			AC_LINK_IFELSE([AC_LANG_PROGRAM([],[      call MPI_Init])],[MPILIBS=" "
				AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])])
fi
AC_LANG_CASE([Fortran 77], [
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(fmpi, MPI_Init, [MPILIBS="-lfmpi"])
	fi
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(fmpich, MPI_Init, [MPILIBS="-lfmpich"])
	fi
],
[Fortran], [
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(fmpi, MPI_Init, [MPILIBS="-lfmpi"])
	fi
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(mpichf90, MPI_Init, [MPILIBS="-lmpichf90"])
	fi
])
if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpi, MPI_Init, [MPILIBS="-lmpi"])
fi
if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpich, MPI_Init, [MPILIBS="-lmpich"])
fi

dnl We have to use AC_TRY_COMPILE and not AC_CHECK_HEADER because the
dnl latter uses $CPP, not $CC (which may be mpicc).
AC_LANG_CASE([C], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[C++], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[Fortran 77], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpif.h])
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[      include 'mpif.h'])],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[Fortran], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpif.h])
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[      include 'mpif.h'])],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi])

AC_LANG_CASE([C], [CC="$ax_mpi_save_CC"],
	[C++], [CXX="$ax_mpi_save_CXX"],
	[Fortran 77], [F77="$ax_mpi_save_F77"],
	[Fortran], [FC="$ax_mpi_save_FC"])

AC_SUBST(MPILIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x = x"$MPILIBS"; then
        $2
        :
else
        ifelse([$1],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$1])
        :
fi
])dnl AX_MPI

# ===========================================================================
#        https://www.gnu.org/software/autoconf-archive/ax_openmp.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_OPENMP([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro tries to find out how to compile programs that use OpenMP a
#   standard API and set of compiler directives for parallel programming
#   (see http://www-unix.mcs/)
#
#   On success, it sets the OPENMP_CFLAGS/OPENMP_CXXFLAGS/OPENMP_F77FLAGS
#   output variable to the flag (e.g. -omp) used both to compile *and* link
#   OpenMP programs in the current language.
#
#   NOTE: You are assumed to not only compile your program with these flags,
#   but also link it with them as well.
#
#   If you want to compile everything with OpenMP, you should set:
#
#     CFLAGS="$CFLAGS $OPENMP_CFLAGS"
#     #OR#  CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
#     #OR#  FFLAGS="$FFLAGS $OPENMP_FFLAGS"
#
#   (depending on the selected language).
#
#   The user can override the default choice by setting the corresponding
#   environment variable (e.g. OPENMP_CFLAGS).
#
#   ACTION-IF-FOUND is a list of shell commands to run if an OpenMP flag is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_OPENMP.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2015 John W. Peterson <jwpeterson@gmail.com>
#   Copyright (c) 2016 Nick R. Papior <nickpapior@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 13

AC_DEFUN([AX_OPENMP], [
AC_PREREQ([2.69]) dnl for _AC_LANG_PREFIX

AC_CACHE_CHECK([for OpenMP flag of _AC_LANG compiler], ax_cv_[]_AC_LANG_ABBREV[]_openmp, [save[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
ax_cv_[]_AC_LANG_ABBREV[]_openmp=unknown
# Flags to try:  -fopenmp (gcc), -mp (SGI & PGI),
#                -qopenmp (icc>=15), -openmp (icc),
#                -xopenmp (Sun), -omp (Tru64),
#                -qsmp=omp (AIX),
#                none
ax_openmp_flags="-fopenmp -openmp -qopenmp -mp -xopenmp -omp -qsmp=omp none"
if test "x$OPENMP_[]_AC_LANG_PREFIX[]FLAGS" != x; then
  ax_openmp_flags="$OPENMP_[]_AC_LANG_PREFIX[]FLAGS $ax_openmp_flags"
fi
for ax_openmp_flag in $ax_openmp_flags; do
  case $ax_openmp_flag in
    none) []_AC_LANG_PREFIX[]FLAGS=$save[]_AC_LANG_PREFIX[] ;;
    *) []_AC_LANG_PREFIX[]FLAGS="$save[]_AC_LANG_PREFIX[]FLAGS $ax_openmp_flag" ;;
  esac
  AC_LINK_IFELSE([AC_LANG_SOURCE([[
@%:@include <omp.h>

static void
parallel_fill(int * data, int n)
{
  int i;
@%:@pragma omp parallel for
  for (i = 0; i < n; ++i)
    data[i] = i;
}

int
main()
{
  int arr[100000];
  omp_set_num_threads(2);
  parallel_fill(arr, 100000);
  return 0;
}
]])],[ax_cv_[]_AC_LANG_ABBREV[]_openmp=$ax_openmp_flag; break],[])
done
[]_AC_LANG_PREFIX[]FLAGS=$save[]_AC_LANG_PREFIX[]FLAGS
])
if test "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp" = "xunknown"; then
  m4_default([$2],:)
else
  if test "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp" != "xnone"; then
    OPENMP_[]_AC_LANG_PREFIX[]FLAGS=$ax_cv_[]_AC_LANG_ABBREV[]_openmp
  fi
  m4_default([$1], [AC_DEFINE(HAVE_OPENMP,1,[Define if OpenMP is enabled])])
fi
])dnl AX_OPENMP

# Configure path for the GNU Scientific Library
# Christopher R. Gabriel <cgabriel@linux.it>, April 2000


AC_DEFUN([AX_PATH_GSL],
[
AC_ARG_WITH(gsl-prefix,[  --with-gsl-prefix=PFX   Prefix where GSL is installed (optional)],
            gsl_prefix="$withval", gsl_prefix="")
AC_ARG_WITH(gsl-exec-prefix,[  --with-gsl-exec-prefix=PFX Exec prefix where GSL is installed (optional)],
            gsl_exec_prefix="$withval", gsl_exec_prefix="")
AC_ARG_ENABLE(gsltest, [  --disable-gsltest       Do not try to compile and run a test GSL program],
		    , enable_gsltest=yes)

  if test "x${GSL_CONFIG+set}" != xset ; then
     if test "x$gsl_prefix" != x ; then
         GSL_CONFIG="$gsl_prefix/bin/gsl-config"
     fi
     if test "x$gsl_exec_prefix" != x ; then
        GSL_CONFIG="$gsl_exec_prefix/bin/gsl-config"
     fi
  fi

  AC_PATH_PROG(GSL_CONFIG, gsl-config, no)
  min_gsl_version=ifelse([$1], ,0.2.5,$1)
  AC_MSG_CHECKING(for GSL - version >= $min_gsl_version)
  no_gsl=""
  if test "$GSL_CONFIG" = "no" ; then
    no_gsl=yes
  else
    GSL_CFLAGS=`$GSL_CONFIG --cflags`
    GSL_LIBS=`$GSL_CONFIG --libs`

    gsl_major_version=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${gsl_major_version}" = "x" ; then
       gsl_major_version=0
    fi

    gsl_minor_version=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${gsl_minor_version}" = "x" ; then
       gsl_minor_version=0
    fi

    gsl_micro_version=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${gsl_micro_version}" = "x" ; then
       gsl_micro_version=0
    fi

    if test "x$enable_gsltest" = "xyes" ; then
      ac_save_CFLAGS="$CFLAGS"
      ac_save_LIBS="$LIBS"
      CFLAGS="$CFLAGS $GSL_CFLAGS"
      LIBS="$LIBS $GSL_LIBS"

      rm -f conf.gsltest
      AC_TRY_RUN([
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* my_strdup (const char *str);

char*
my_strdup (const char *str)
{
  char *new_str;
  
  if (str)
    {
      new_str = (char *)malloc ((strlen (str) + 1) * sizeof(char));
      strcpy (new_str, str);
    }
  else
    new_str = NULL;
  
  return new_str;
}

int main (void)
{
  int major = 0, minor = 0, micro = 0;
  int n;
  char *tmp_version;

  system ("touch conf.gsltest");

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  tmp_version = my_strdup("$min_gsl_version");

  n = sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) ;

  if (n != 2 && n != 3) {
     printf("%s, bad version string\n", "$min_gsl_version");
     exit(1);
   }

   if (($gsl_major_version > major) ||
      (($gsl_major_version == major) && ($gsl_minor_version > minor)) ||
      (($gsl_major_version == major) && ($gsl_minor_version == minor) && ($gsl_micro_version >= micro)))
     { 
       exit(0);
     }   
   else
     {
       exit(1);
     }
}

],, no_gsl=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
       CFLAGS="$ac_save_CFLAGS"
       LIBS="$ac_save_LIBS"
     fi
  fi
  if test "x$no_gsl" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test "$GSL_CONFIG" = "no" ; then
       echo "*** The gsl-config script installed by GSL could not be found"
       echo "*** If GSL was installed in PREFIX, make sure PREFIX/bin is in"
       echo "*** your path, or set the GSL_CONFIG environment variable to the"
       echo "*** full path to gsl-config."
     else
       if test -f conf.gsltest ; then
        :
       else
          echo "*** Could not run GSL test program, checking why..."
          CFLAGS="$CFLAGS $GSL_CFLAGS"
          LIBS="$LIBS $GSL_LIBS"
          AC_TRY_LINK([
#include <stdio.h>
],      [ return 0; ],
        [ echo "*** The test program compiled, but did not run. This usually means"
          echo "*** that the run-time linker is not finding GSL or finding the wrong"
          echo "*** version of GSL. If it is not finding GSL, you'll need to set your"
          echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
          echo "*** to the installed location  Also, make sure you have run ldconfig if that"
          echo "*** is required on your system"
	  echo "***"
          echo "*** If you have an old version installed, it is best to remove it, although"
          echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH"],
        [ echo "*** The test program failed to compile or link. See the file config.log for the"
          echo "*** exact error that occured. This usually means GSL was incorrectly installed"
          echo "*** or that you have moved GSL since it was installed. In the latter case, you"
          echo "*** may want to edit the gsl-config script: $GSL_CONFIG" ])
          CFLAGS="$ac_save_CFLAGS"
          LIBS="$ac_save_LIBS"
       fi
     fi
#     GSL_CFLAGS=""
#     GSL_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  AC_SUBST(GSL_CFLAGS)
  AC_SUBST(GSL_LIBS)
  rm -f conf.gsltest
])

AU_ALIAS([AM_PATH_GSL], [AX_PATH_GSL])

# Copyright (C) 1999-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.


# AM_PATH_PYTHON([MINIMUM-VERSION], [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ---------------------------------------------------------------------------
# Adds support for distributing Python modules and packages.  To
# install modules, copy them to $(pythondir), using the python_PYTHON
# automake variable.  To install a package with the same name as the
# automake package, install to $(pkgpythondir), or use the
# pkgpython_PYTHON automake variable.
#
# The variables $(pyexecdir) and $(pkgpyexecdir) are provided as
# locations to install python extension modules (shared libraries).
# Another macro is required to find the appropriate flags to compile
# extension modules.
#
# If your package is configured with a different prefix to python,
# users will have to add the install directory to the PYTHONPATH
# environment variable, or create a .pth file (see the python
# documentation for details).
#
# If the MINIMUM-VERSION argument is passed, AM_PATH_PYTHON will
# cause an error if the version of python installed on the system
# doesn't meet the requirement.  MINIMUM-VERSION should consist of
# numbers and dots only.
AC_DEFUN([AM_PATH_PYTHON],
 [
  dnl Find a Python interpreter.  Python versions prior to 2.0 are not
  dnl supported. (2.0 was released on October 16, 2000).
  m4_define_default([_AM_PYTHON_INTERPRETER_LIST],
[python python2 python3 dnl
 python3.9 python3.8 python3.7 python3.6 python3.5 python3.4 python3.3 dnl
 python3.2 python3.1 python3.0 dnl
 python2.7 python2.6 python2.5 python2.4 python2.3 python2.2 python2.1 dnl
 python2.0])

  AC_ARG_VAR([PYTHON], [the Python interpreter])

  m4_if([$1],[],[
    dnl No version check is needed.
    # Find any Python interpreter.
    if test -z "$PYTHON"; then
      AC_PATH_PROGS([PYTHON], _AM_PYTHON_INTERPRETER_LIST, :)
    fi
    am_display_PYTHON=python
  ], [
    dnl A version check is needed.
    if test -n "$PYTHON"; then
      # If the user set $PYTHON, use it and don't search something else.
      AC_MSG_CHECKING([whether $PYTHON version is >= $1])
      AM_PYTHON_CHECK_VERSION([$PYTHON], [$1],
			      [AC_MSG_RESULT([yes])],
			      [AC_MSG_RESULT([no])
			       AC_MSG_ERROR([Python interpreter is too old])])
      am_display_PYTHON=$PYTHON
    else
      # Otherwise, try each interpreter until we find one that satisfies
      # VERSION.
      AC_CACHE_CHECK([for a Python interpreter with version >= $1],
	[am_cv_pathless_PYTHON],[
	for am_cv_pathless_PYTHON in _AM_PYTHON_INTERPRETER_LIST none; do
	  test "$am_cv_pathless_PYTHON" = none && break
	  AM_PYTHON_CHECK_VERSION([$am_cv_pathless_PYTHON], [$1], [break])
	done])
      # Set $PYTHON to the absolute path of $am_cv_pathless_PYTHON.
      if test "$am_cv_pathless_PYTHON" = none; then
	PYTHON=:
      else
        AC_PATH_PROG([PYTHON], [$am_cv_pathless_PYTHON])
      fi
      am_display_PYTHON=$am_cv_pathless_PYTHON
    fi
  ])

  if test "$PYTHON" = :; then
  dnl Run any user-specified action, or abort.
    m4_default([$3], [AC_MSG_ERROR([no suitable Python interpreter found])])
  else

  dnl Query Python for its version number.  Although site.py simply uses
  dnl sys.version[:3], printing that failed with Python 3.10, since the
  dnl trailing zero was eliminated. So now we output just the major
  dnl and minor version numbers, as numbers. Apparently the tertiary
  dnl version is not of interest.

  AC_CACHE_CHECK([for $am_display_PYTHON version], [am_cv_python_version],
    [am_cv_python_version=`$PYTHON -c "import sys; print('%u.%u' % sys.version_info[[:2]])"`])
  AC_SUBST([PYTHON_VERSION], [$am_cv_python_version])

  dnl Use the values of $prefix and $exec_prefix for the corresponding
  dnl values of PYTHON_PREFIX and PYTHON_EXEC_PREFIX.  These are made
  dnl distinct variables so they can be overridden if need be.  However,
  dnl general consensus is that you shouldn't need this ability.

  AC_SUBST([PYTHON_PREFIX], ['${prefix}'])
  AC_SUBST([PYTHON_EXEC_PREFIX], ['${exec_prefix}'])

  dnl At times (like when building shared libraries) you may want
  dnl to know which OS platform Python thinks this is.

  AC_CACHE_CHECK([for $am_display_PYTHON platform], [am_cv_python_platform],
    [am_cv_python_platform=`$PYTHON -c "import sys; sys.stdout.write(sys.platform)"`])
  AC_SUBST([PYTHON_PLATFORM], [$am_cv_python_platform])

  # Just factor out some code duplication.
  am_python_setup_sysconfig="\
import sys
# Prefer sysconfig over distutils.sysconfig, for better compatibility
# with python 3.x.  See automake bug#10227.
try:
    import sysconfig
except ImportError:
    can_use_sysconfig = 0
else:
    can_use_sysconfig = 1
# Can't use sysconfig in CPython 2.7, since it's broken in virtualenvs:
# <https://github.com/pypa/virtualenv/issues/118>
try:
    from platform import python_implementation
    if python_implementation() == 'CPython' and sys.version[[:3]] == '2.7':
        can_use_sysconfig = 0
except ImportError:
    pass"

  dnl Set up 4 directories:

  dnl pythondir -- where to install python scripts.  This is the
  dnl   site-packages directory, not the python standard library
  dnl   directory like in previous automake betas.  This behavior
  dnl   is more consistent with lispdir.m4 for example.
  dnl Query distutils for this directory.
  AC_CACHE_CHECK([for $am_display_PYTHON script directory],
    [am_cv_python_pythondir],
    [if test "x$prefix" = xNONE
     then
       am_py_prefix=$ac_default_prefix
     else
       am_py_prefix=$prefix
     fi
     am_cv_python_pythondir=`$PYTHON -c "
$am_python_setup_sysconfig
if can_use_sysconfig:
    sitedir = sysconfig.get_path('purelib')
else:
    from distutils import sysconfig
    sitedir = sysconfig.get_python_lib(0, 0)
sys.stdout.write(sitedir)"`
     case $am_cv_python_pythondir in
     $am_py_prefix*)
       am__strip_prefix=`echo "$am_py_prefix" | sed 's|.|.|g'`
       am_cv_python_pythondir=`echo "$am_cv_python_pythondir" | sed "s,^$am__strip_prefix,$PYTHON_PREFIX,"`
       ;;
     *)
       case $am_py_prefix in
         /usr|/System*) ;;
         *)
	  am_cv_python_pythondir=$PYTHON_PREFIX/lib/python$PYTHON_VERSION/site-packages
	  ;;
       esac
       ;;
     esac
    ])
  AC_SUBST([pythondir], [$am_cv_python_pythondir])

  dnl pkgpythondir -- $PACKAGE directory under pythondir.  Was
  dnl   PYTHON_SITE_PACKAGE in previous betas, but this naming is
  dnl   more consistent with the rest of automake.

  AC_SUBST([pkgpythondir], [\${pythondir}/$PACKAGE])

  dnl pyexecdir -- directory for installing python extension modules
  dnl   (shared libraries)
  dnl Query distutils for this directory.
  AC_CACHE_CHECK([for $am_display_PYTHON extension module directory],
    [am_cv_python_pyexecdir],
    [if test "x$exec_prefix" = xNONE
     then
       am_py_exec_prefix=$am_py_prefix
     else
       am_py_exec_prefix=$exec_prefix
     fi
     am_cv_python_pyexecdir=`$PYTHON -c "
$am_python_setup_sysconfig
if can_use_sysconfig:
    sitedir = sysconfig.get_path('platlib')
else:
    from distutils import sysconfig
    sitedir = sysconfig.get_python_lib(1, 0)
sys.stdout.write(sitedir)"`
     case $am_cv_python_pyexecdir in
     $am_py_exec_prefix*)
       am__strip_prefix=`echo "$am_py_exec_prefix" | sed 's|.|.|g'`
       am_cv_python_pyexecdir=`echo "$am_cv_python_pyexecdir" | sed "s,^$am__strip_prefix,$PYTHON_EXEC_PREFIX,"`
       ;;
     *)
       case $am_py_exec_prefix in
         /usr|/System*) ;;
         *)
	   am_cv_python_pyexecdir=$PYTHON_EXEC_PREFIX/lib/python$PYTHON_VERSION/site-packages
	   ;;
       esac
       ;;
     esac
    ])
  AC_SUBST([pyexecdir], [$am_cv_python_pyexecdir])

  dnl pkgpyexecdir -- $(pyexecdir)/$(PACKAGE)

  AC_SUBST([pkgpyexecdir], [\${pyexecdir}/$PACKAGE])

  dnl Run any user-specified action.
  $2
  fi

])


# AM_PYTHON_CHECK_VERSION(PROG, VERSION, [ACTION-IF-TRUE], [ACTION-IF-FALSE])
# ---------------------------------------------------------------------------
# Run ACTION-IF-TRUE if the Python interpreter PROG has version >= VERSION.
# Run ACTION-IF-FALSE otherwise.
# This test uses sys.hexversion instead of the string equivalent (first
# word of sys.version), in order to cope with versions such as 2.2c1.
# This supports Python 2.0 or higher. (2.0 was released on October 16, 2000).
AC_DEFUN([AM_PYTHON_CHECK_VERSION],
 [prog="import sys
# split strings by '.' and convert to numeric.  Append some zeros
# because we need at least 4 digits for the hex conversion.
# map returns an iterator in Python 3.0 and a list in 2.x
minver = list(map(int, '$2'.split('.'))) + [[0, 0, 0]]
minverhex = 0
# xrange is not present in Python 3.0 and range returns an iterator
for i in list(range(0, 4)): minverhex = (minverhex << 8) + minver[[i]]
sys.exit(sys.hexversion < minverhex)"
  AS_IF([AM_RUN_LOG([$1 -c "$prog"])], [$3], [$4])])

# Copyright (C) 2001-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_RUN_LOG(COMMAND)
# -------------------
# Run COMMAND, save the exit status in ac_status, and log it.
# (This has been adapted from Autoconf's _AC_RUN_LOG macro.)
AC_DEFUN([AM_RUN_LOG],
[{ echo "$as_me:$LINENO: $1" >&AS_MESSAGE_LOG_FD
   ($1) >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD
   ac_status=$?
   echo "$as_me:$LINENO: \$? = $ac_status" >&AS_MESSAGE_LOG_FD
   (exit $ac_status); }])

m4_include([m4/ax_mpi.m4])
m4_include([m4/ax_openmp.m4])
m4_include([m4/cuda.m4])
m4_include([m4/fftw.m4])
m4_include([m4/fortran.m4])
m4_include([m4/gmp.m4])
m4_include([m4/gsl.m4])
m4_include([m4/java.m4])
m4_include([m4/m4_ax_lib_hdf5.m4])
m4_include([m4/m4_ax_lib_netcdf4.m4])
m4_include([m4/mpfr.m4])
m4_include([m4/opencl.m4])
m4_include([m4/pic_flags.m4])
m4_include([m4/pkg.m4])
m4_include([m4/python_dev.m4])
m4_include([m4/python_module.m4])
m4_include([m4/utils.m4])
