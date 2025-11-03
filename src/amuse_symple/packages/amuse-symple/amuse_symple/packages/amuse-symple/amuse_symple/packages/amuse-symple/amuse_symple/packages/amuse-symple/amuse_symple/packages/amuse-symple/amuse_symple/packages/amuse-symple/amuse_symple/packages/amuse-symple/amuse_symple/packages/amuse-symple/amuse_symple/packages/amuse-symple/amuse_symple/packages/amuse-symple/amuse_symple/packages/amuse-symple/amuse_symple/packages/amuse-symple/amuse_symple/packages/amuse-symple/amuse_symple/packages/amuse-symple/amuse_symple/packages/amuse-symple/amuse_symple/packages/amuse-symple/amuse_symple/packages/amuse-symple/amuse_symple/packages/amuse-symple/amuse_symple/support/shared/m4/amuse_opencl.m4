# Helper macro for detecting OpenCL for C.

# AMUSE_OPENCL()
#
# Sets FOUND_OPENCL to "yes" if OpenCL is found, to "no" if it is not.
# Sets OPENCL_CFLAGS and OPENCL_LIBS as needed.
# Calls AC_SUBST on all of the above.
#
AC_DEFUN([AMUSE_OPENCL], [
    amuse_save_CFLAGS="$CFLAGS"
    amuse_save_LIBS="$LIBS"

    AC_MSG_CHECKING([for OpenCL])

    # If the user overrode the variables, then we'll use what they set and verify that
    # it works.
    AS_IF([test "x$OPENCL_CFLAGS" != "x"], [
        amuse_user_OPENCL_CFLAGS="$OPENCL_CFLAGS"
        CFLAGS="$OPENCL_CFLAGS $CFLAGS"
    ])
    AS_IF([test "x$OPENCL_LIBS" != "x"], [
        amuse_user_OPENCL_LIBS="$OPENCL_LIBS"
        LIBS="$OPENCL_LIBS $LIBS"
    ])

    # We only use OpenCL from C, in huayno
    AC_LANG_PUSH([C])

    # Search for the headers
    m4_define([amuse_opencl_test_program], [AC_LANG_PROGRAM([
        #if defined(__APPLE__) || defined(__MACOSX)
        #include <OpenCL/cl.h>
        #else
        #include <CL/cl.h>
        #endif
    ], [clFinish(0)])])

    AC_COMPILE_IFELSE([amuse_opencl_test_program], [
        # Search for the library as well
        LIBS="-lOpenCL $LIBS"
        AC_LINK_IFELSE([amuse_opencl_test_program], [
            FOUND_OPENCL="yes"
            OPENCL_LIBS="-lOpenCL"
            OPENCL_CFLAGS=""
        ], [
            FOUND_OPENCL="no"
        ])
    ], [
        FOUND_OPENCL="no"
    ])

    AC_LANG_POP([C])

    LIBS="$amuse_save_LIBS"
    CFLAGS="$amuse_save_CFLAGS"

    AC_MSG_RESULT([$FOUND_OPENCL])

    AC_SUBST([FOUND_OPENCL])
    AC_SUBST([OPENCL_CFLAGS])
    AC_SUBST([OPENCL_LIBS])
])

