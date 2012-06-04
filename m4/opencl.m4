AC_DEFUN([AX_LANG_COMPILER_MS],
[AC_CACHE_CHECK([whether we are using the Microsoft _AC_LANG compiler],
                [ax_cv_[]_AC_LANG_ABBREV[]_compiler_ms],
[AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[#ifndef _MSC_VER
       choke me
#endif
]])],
                   [ax_compiler_ms=yes],
                   [ax_compiler_ms=no])
ax_cv_[]_AC_LANG_ABBREV[]_compiler_ms=$ax_compiler_ms
])])

# -*- mode: autoconf -*-
#
# AX_CHECK_CL
#
# Check for an OpenCL implementation.  If CL is found, the required compiler
# and linker flags are included in the output variables "CL_CFLAGS" and
# "CL_LIBS", respectively.  If no usable CL implementation is found, "no_cl"
# is set to "yes".
#
# If the header "CL/cl.h" is found, "HAVE_CL_CL_H" is defined.  If the header
# "OpenCL/cl.h" is found, HAVE_OPENCL_CL_H is defined.  These preprocessor
# definitions may not be mutually exclusive.
#
# Based on AX_CHECK_GL, version: 2.4 author: Braden McDaniel
# <braden@endoframe.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
#
# As a special exception, the you may copy, distribute and modify the
# configure scripts that are the output of Autoconf when processing
# the Macro.  You need not follow the terms of the GNU General Public
# License when using or distributing such scripts.
#
AC_DEFUN([AX_CHECK_CL],
[
AC_ARG_WITH(opencl,
             AC_HELP_STRING([--with-opencl=PFX], [Prefix where OpenCl has been installed] ),
             [
                test "$withval" = no && AC_MSG_WARN([opencl is an optional package for sapporo2 ])
                
                if test "x$withval" != xyes; then
                   opencl_prefix="$withval" 
                    with_opencl=yes
                fi
                if test "x$withval" == xno; then
                   with_opencl=no
                fi
                ],
             [ with_opencl=yes ] 
        )
    
        AS_IF([test "x$with_opencl" != xno ],
        [
        
if test x$opencl_prefix == x; then
   if test x$PREFIX != x; then
        opencl_prefix=$PREFIX 
   fi  
fi

if test x$opencl_prefix != x; then
    CL_CFLAGS="-I$opencl_prefix/include"
    CL_LIBS="-L$opencl_prefix/lib -L$opencl_prefix/lib64"
fi

AC_LANG_PUSH([C])
AX_LANG_COMPILER_MS

ax_save_CPPFLAGS=$CPPFLAGS
CPPFLAGS="$CL_CFLAGS $CPPFLAGS"
AC_CHECK_HEADERS([CL/cl.h OpenCL/cl.h])
CPPFLAGS=$ax_save_CPPFLAGS

AC_CHECK_HEADERS([windows.h])

m4_define([AX_CHECK_CL_PROGRAM],
          [AC_LANG_PROGRAM([[
# if defined(HAVE_WINDOWS_H) && defined(_WIN32)
#   include <windows.h>
# endif
# ifdef HAVE_CL_CL_H
#   include <CL/cl.h>
# elif defined(HAVE_OPENCL_CL_H)
#   include <OpenCL/cl.h>
# else
#   error no cl.h
# endif]],
                           [[clFinish(0)]])])

AC_CACHE_CHECK([for OpenCL library], [ax_cv_check_cl_libcl],
[ax_cv_check_cl_libcl=no
case $host_cpu in
  x86_64) ax_check_cl_libdir=lib64 ;;
  *)      ax_check_cl_libdir=lib ;;
esac
ax_save_CPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS $CL_CFLAGS"
ax_save_LIBS=$LIBS
LIBS=""
ax_check_libs="-lOpenCL -lCL"
for ax_lib in $ax_check_libs; do
  AS_IF([test X$ax_compiler_ms = Xyes],
        [ax_try_lib=`echo $ax_lib | $SED -e 's/^-l//' -e 's/$/.lib/'`],
        [ax_try_lib=$ax_lib])
  LIBS="$ax_try_lib $CL_LIBS $ax_save_LIBS"
AC_LINK_IFELSE([AX_CHECK_CL_PROGRAM],
               [ax_cv_check_cl_libcl=$ax_try_lib; break],
               [ax_check_cl_nvidia_flags="-L/usr/$ax_check_cl_libdir/nvidia" LIBS="$ax_try_lib $ax_check_cl_nvidia_flags $CL_LIBS $ax_save_LIBS"
               AC_LINK_IFELSE([AX_CHECK_CL_PROGRAM],
                              [ax_cv_check_cl_libcl="$ax_try_lib $ax_check_cl_nvidia_flags"; break],
                              [ax_check_cl_dylib_flag='-dylib_file /System/Library/Frameworks/OpenCL.framework/Versions/A/Libraries/libCL.dylib:/System/Library/Frameworks/OpenCL.framework/Versions/A/Libraries/libCL.dylib' LIBS="$ax_try_lib $ax_check_cl_dylib_flag $CL_LIBS $ax_save_LIBS"
                              AC_LINK_IFELSE([AX_CHECK_CL_PROGRAM],
                                             [ax_cv_check_cl_libcl="$ax_try_lib $ax_check_cl_dylib_flag"; break])])])
done

AS_IF([test "X$ax_cv_check_cl_libcl" = Xno -a X$no_x = Xyes],
      [LIBS='-framework OpenCL'
      AC_LINK_IFELSE([AX_CHECK_CL_PROGRAM],
                     [ax_cv_check_cl_libcl=$LIBS])])

LIBS=$ax_save_LIBS
CPPFLAGS=$ax_save_CPPFLAGS])

AS_IF([test "X$ax_cv_check_cl_libcl" = Xno],
      [have_cl=no; CL_CFLAGS=""; CL_LIBS=""],
      [have_cl=yes; CL_LIBS="$CL_LIBS $ax_cv_check_cl_libcl"])
AC_LANG_POP([C])
], [have_cl=no
CL_CFLAGS=""
CL_LIBS=""
])

FOUND_CL=$have_cl

AC_SUBST([FOUND_CL])
AC_SUBST([CL_CFLAGS])
AC_SUBST([CL_LIBS])

])dnl
