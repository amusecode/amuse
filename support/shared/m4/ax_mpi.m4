# ===========================================================================
#          http://www.gnu.org/software/autoconf-archive/ax_mpi.html
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
#   If you want to compile everything with MPI, you should use something
#   like this for C:
#
#     if test -z "$CC" && test -n "$MPICC"; then
#       CC="$MPICC"
#     fi
#     AC_PROG_CC
#     AX_MPI
#     CC="$MPICC"
#     LIBS="$MPILIBS $LIBS"
#
#   and similar for C++ (change all instances of CC to CXX), Fortran 77
#   (with F77 instead of CC) or Fortran (with FC instead of CC).
#
#   NOTE: The above assumes that you will use $CC (or whatever) for linking
#   as well as for compiling. (This is the default for automake and most
#   Makefiles.)
#
#   The user can force a particular library/compiler by setting the
#   MPICC/MPICXX/MPIF77/MPIFC and/or MPILIBS environment variables.
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
#   with this program. If not, see <http://www.gnu.org/licenses/>.
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

#serial 7

AU_ALIAS([ACX_MPI], [AX_MPI])
AC_DEFUN([AX_MPI], [
AC_PREREQ([2.71]) dnl for AC_LANG_CASE

AC_LANG_CASE([C], [
	AC_REQUIRE([AC_PROG_CC])
	AC_ARG_VAR(MPICC,[MPI C compiler command])
	AC_CHECK_PROGS(MPICC, mpicc hcc mpxlc_r mpxlc mpcc cmpicc, $CC)
	ax_mpi_save_CC="$CC"
	CC="$MPICC"
	AC_SUBST(MPICC)
 	AC_MSG_CHECKING([checking MPI C flags])
 	ax_mpi_c_flags="`$MPICC -showme:compile 2>/dev/null| cut -d\  -f2-`"
 	ax_mpi_c_libs="`$MPICC -showme:link 2>/dev/null| cut -d\  -f2-`"
        AS_IF([test "x$ax_mpi_c_flags" = "x"],[
          ax_mpi_c_flags="`$MPICC -show -c 2>/dev/null| cut -d\  -f2-|sed s/-c\ //`"
          ax_mpi_c_libs="`$MPICC -show 2>/dev/null| cut -d\  -f2-`"
	   AS_IF([test "x$ax_mpi_c_flags" = "x"],[AC_MSG_RESULT([could not determine c flags from show functions])],[AC_MSG_RESULT([flags found])])
          
        ], [
          AC_MSG_RESULT([flags found])])
	MPI_CFLAGS="$ax_mpi_c_flags"
	MPI_CLIBS="$ax_mpi_c_libs"
	AC_SUBST(MPI_CFLAGS)
	AC_SUBST(MPI_CLIBS)

],
[C++], [
	AC_REQUIRE([AC_PROG_CXX])
	AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
	AC_CHECK_PROGS(MPICXX, mpicxx mpiCC mpic++ hcp mpxlC_r mpxlC mpCC cmpic++, $CXX)
	ax_mpi_save_CXX="$CXX"
	CXX="$MPICXX"
	AC_SUBST(MPICXX)

 	AC_MSG_CHECKING([checking MPI C++ flags])
 	ax_mpi_cc_flags="`$MPICXX -showme:compile 2>/dev/null| cut -d\  -f2-`"
 	ax_mpi_cc_libs="`$MPICXX -showme:link 2>/dev/null| cut -d\  -f2-`"
        AS_IF([test "x$ax_mpi_cc_flags" = "x"],[
          ax_mpi_cc_flags="`$MPICXX -show -c 2>/dev/null| cut -d\  -f2-|sed s/-c\ //`"
          ax_mpi_cc_libs="`$MPICXX -show 2>/dev/null| cut -d\  -f2-`"
	   AS_IF([test "x$ax_mpi_cc_flags" = "x"],[AC_MSG_RESULT([could not determine C++ flags from show functions])],[AC_MSG_RESULT([flags found])])
          
        ], [
          AC_MSG_RESULT([flags found])])
	MPI_CXXFLAGS="$ax_mpi_cc_flags"
	MPI_CXXLIBS="$ax_mpi_cc_libs"
	AC_SUBST(MPI_CXXFLAGS)
	AC_SUBST(MPI_CXXLIBS)

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
	AC_MSG_CHECKING([checking MPI Fortran flags])
 	ax_mpi_fc_flags="`$MPIFC -showme:compile 2>/dev/null| cut -d\  -f2-`"
 	ax_mpi_fc_libs="`$MPIFC -showme:link 2>/dev/null| cut -d\  -f2-`"
        AS_IF([test "x$ax_mpi_fc_flags" = "x"],[
          ax_mpi_fc_flags="`$MPIFC -show -c 2>/dev/null| cut -d\  -f2-|sed s/-c\ //`"
          ax_mpi_fc_libs="`$MPIFC -show 2>/dev/null| cut -d\  -f2-`"
	   AS_IF([test "x$ax_mpi_fc_flags" = "x"],[AC_MSG_RESULT([could not determine c flags from show functions])],[AC_MSG_RESULT([flags found])])
          
        ], [
          AC_MSG_RESULT([flags found])])
	MPI_FCFLAGS="$ax_mpi_fc_flags"
	MPI_FCLIBS="$ax_mpi_fc_libs"
	AC_SUBST(MPI_FCFLAGS)
	AC_SUBST(MPI_FCLIBS)
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

dnl We have to use AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[]])],[],[]) and not AC_CHECK_HEADER because the
dnl latter uses $CPP, not $CC (which may be mpicc).
AC_LANG_CASE([C], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]], [[]])],[AC_MSG_RESULT(yes)],[MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[C++], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]], [[]])],[AC_MSG_RESULT(yes)],[MPILIBS=""
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
