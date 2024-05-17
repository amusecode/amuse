# generated automatically by aclocal 1.16.5 -*- Autoconf -*-

# Copyright (C) 1996-2021 Free Software Foundation, Inc.

# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

m4_ifndef([AC_CONFIG_MACRO_DIRS], [m4_defun([_AM_CONFIG_MACRO_DIRS], [])m4_defun([AC_CONFIG_MACRO_DIRS], [_AM_CONFIG_MACRO_DIRS($@)])])
# ===========================================================================
#         https://www.gnu.org/software/autoconf-archive/ax_blas.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for a library that implements the BLAS linear-algebra
#   interface (see http://www.netlib.org/blas/). On success, it sets the
#   BLAS_LIBS output variable to hold the requisite library linkages.
#
#   To link with BLAS, you should link with:
#
#     $BLAS_LIBS $LIBS $FLIBS
#
#   in that order. FLIBS is the output variable of the
#   AC_F77_LIBRARY_LDFLAGS macro (called if necessary by AX_BLAS), and is
#   sometimes necessary in order to link with F77 libraries. Users will also
#   need to use AC_F77_DUMMY_MAIN (see the autoconf manual), for the same
#   reason.
#
#   Many libraries are searched for, from ATLAS to CXML to ESSL. The user
#   may also use --with-blas=<lib> in order to use some specific BLAS
#   library <lib>. In order to link successfully, however, be aware that you
#   will probably need to use the same Fortran compiler (which can be set
#   via the F77 env. var.) as was used to compile the BLAS library.
#
#   ACTION-IF-FOUND is a list of shell commands to run if a BLAS library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_BLAS.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2019 Geoffrey M. Oxberry <goxberry@gmail.com>
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

#serial 17

AU_ALIAS([ACX_BLAS], [AX_BLAS])
AC_DEFUN([AX_BLAS], [
AC_PREREQ([2.55])
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
AC_REQUIRE([AC_CANONICAL_HOST])
ax_blas_ok=no

AC_ARG_WITH(blas,
	[AS_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) ax_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.dylib | *.dylib.* | *.o)
		BLAS_LIBS="$with_blas"
	;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac

# Get fortran linker names of BLAS functions to check for.
AC_F77_FUNC(sgemm)
AC_F77_FUNC(dgemm)

ax_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $ax_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for $sgemm in $BLAS_LIBS])
	AC_LINK_IFELSE([AC_LANG_CALL([], [$sgemm])], [ax_blas_ok=yes], [BLAS_LIBS=""])
	AC_MSG_RESULT($ax_blas_ok)
	LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $ax_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_MSG_CHECKING([if $sgemm is being linked in already])
	AC_LINK_IFELSE([AC_LANG_CALL([], [$sgemm])], [ax_blas_ok=yes])
	AC_MSG_RESULT($ax_blas_ok)
	LIBS="$save_LIBS"
fi

# BLAS in OpenBLAS library? (http://xianyi.github.com/OpenBLAS/)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(openblas, $sgemm, [ax_blas_ok=yes
			                BLAS_LIBS="-lopenblas"])
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(f77blas, $sgemm,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[ax_blas_ok=yes
			 BLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(dgemm, $dgemm,
		[AC_CHECK_LIB(sgemm, $sgemm,
			[ax_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
fi

# BLAS in Intel MKL library?
if test $ax_blas_ok = no; then
	# MKL for gfortran
	if test x"$ac_cv_fc_compiler_gnu" = xyes; then
		# 64 bit
		if test $host_cpu = x86_64; then
			AC_CHECK_LIB(mkl_gf_lp64, $sgemm,
			[ax_blas_ok=yes;BLAS_LIBS="-lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread"],,
			[-lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread])
		# 32 bit
		elif test $host_cpu = i686; then
			AC_CHECK_LIB(mkl_gf, $sgemm,
				[ax_blas_ok=yes;BLAS_LIBS="-lmkl_gf -lmkl_sequential -lmkl_core -lpthread"],,
				[-lmkl_gf -lmkl_sequential -lmkl_core -lpthread])
		fi
	# MKL for other compilers (Intel, PGI, ...?)
	else
		# 64-bit
		if test $host_cpu = x86_64; then
			AC_CHECK_LIB(mkl_intel_lp64, $sgemm,
				[ax_blas_ok=yes;BLAS_LIBS="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread"],,
				[-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread])
		# 32-bit
		elif test $host_cpu = i686; then
			AC_CHECK_LIB(mkl_intel, $sgemm,
				[ax_blas_ok=yes;BLAS_LIBS="-lmkl_intel -lmkl_sequential -lmkl_core -lpthread"],,
				[-lmkl_intel -lmkl_sequential -lmkl_core -lpthread])
		fi
	fi
fi
# Old versions of MKL
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(mkl, $sgemm, [ax_blas_ok=yes;BLAS_LIBS="-lmkl -lguide -lpthread"],,[-lguide -lpthread])
fi

# BLAS in Apple vecLib library?
if test $ax_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="-framework vecLib $LIBS"
	AC_MSG_CHECKING([for $sgemm in -framework vecLib])
	AC_LINK_IFELSE([AC_LANG_CALL([], [$sgemm])], [ax_blas_ok=yes;BLAS_LIBS="-framework vecLib"])
	AC_MSG_RESULT($ax_blas_ok)
	LIBS="$save_LIBS"
fi

# BLAS in Alpha CXML library?
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(cxml, $sgemm, [ax_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(dxml, $sgemm, [ax_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $ax_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, $sgemm,
				[BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 ax_blas_ok=yes],[],[-lsunmath])])
	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(scs, $sgemm, [ax_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $sgemm,
		     [ax_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(essl, $sgemm,
			[ax_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
fi

# Generic BLAS library?
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm, [ax_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

AC_SUBST(BLAS_LIBS)

LIBS="$ax_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_blas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        ax_blas_ok=no
        $2
fi
])dnl AX_BLAS

# ===========================================================================
#        https://www.gnu.org/software/autoconf-archive/ax_lapack.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for a library that implements the LAPACK linear-algebra
#   interface (see http://www.netlib.org/lapack/). On success, it sets the
#   LAPACK_LIBS output variable to hold the requisite library linkages.
#
#   To link with LAPACK, you should link with:
#
#     $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
#
#   in that order. BLAS_LIBS is the output variable of the AX_BLAS macro,
#   called automatically. FLIBS is the output variable of the
#   AC_F77_LIBRARY_LDFLAGS macro (called if necessary by AX_BLAS), and is
#   sometimes necessary in order to link with F77 libraries. Users will also
#   need to use AC_F77_DUMMY_MAIN (see the autoconf manual), for the same
#   reason.
#
#   The user may also use --with-lapack=<lib> in order to use some specific
#   LAPACK library <lib>. In order to link successfully, however, be aware
#   that you will probably need to use the same Fortran compiler (which can
#   be set via the F77 env. var.) as was used to compile the LAPACK and BLAS
#   libraries.
#
#   ACTION-IF-FOUND is a list of shell commands to run if a LAPACK library
#   is found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it
#   is not found. If ACTION-IF-FOUND is not specified, the default action
#   will define HAVE_LAPACK.
#
# LICENSE
#
#   Copyright (c) 2009 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2019 Geoffrey M. Oxberry <goxberry@gmail.com>
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

#serial 10

AU_ALIAS([ACX_LAPACK], [AX_LAPACK])
AC_DEFUN([AX_LAPACK], [
AC_REQUIRE([AX_BLAS])
ax_lapack_ok=no

AC_ARG_WITH(lapack,
        [AS_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])
case $with_lapack in
        yes | "") ;;
        no) ax_lapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.dylib | *.dylib.* | *.o)
                 LAPACK_LIBS="$with_lapack"
        ;;
        *) LAPACK_LIBS="-l$with_lapack" ;;
esac

# Get fortran linker name of LAPACK function to check for.
AC_F77_FUNC(cheev)

# We cannot use LAPACK if BLAS is not found
if test "x$ax_blas_ok" != xyes; then
        ax_lapack_ok=noblas
        LAPACK_LIBS=""
fi

# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
        AC_MSG_CHECKING([for $cheev in $LAPACK_LIBS])
        AC_LINK_IFELSE([AC_LANG_CALL([], [$cheev])], [ax_lapack_ok=yes], [LAPACK_LIBS=""])
        AC_MSG_RESULT($ax_lapack_ok)
        LIBS="$save_LIBS"
        if test $ax_lapack_ok = no; then
                LAPACK_LIBS=""
        fi
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $ax_lapack_ok = no; then
        save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
        AC_CHECK_FUNC($cheev, [ax_lapack_ok=yes])
        LIBS="$save_LIBS"
fi

# Generic LAPACK library?
for lapack in lapack lapack_rs6k; do
        if test $ax_lapack_ok = no; then
                save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
                AC_CHECK_LIB($lapack, $cheev,
                    [ax_lapack_ok=yes; LAPACK_LIBS="-l$lapack"], [], [$FLIBS])
                LIBS="$save_LIBS"
        fi
done

AC_SUBST(LAPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_lapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
        :
else
        ax_lapack_ok=no
        $2
fi
])dnl AX_LAPACK

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

m4_include([shared/m4/ax_mpi.m4])
m4_include([shared/m4/conda_packages.m4])
m4_include([shared/m4/cuda.m4])
m4_include([shared/m4/fftw.m4])
m4_include([shared/m4/gmp.m4])
m4_include([shared/m4/gsl.m4])
m4_include([shared/m4/mpfr.m4])
m4_include([shared/m4/opencl.m4])
m4_include([shared/m4/pkg.m4])
