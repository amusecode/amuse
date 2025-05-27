# Helper macro for detecting libhealpix_cxx
#
# AMUSE_LIB_HEALPIX_CXX()
#
# Searches for libhealpix_cxx and sets FOUND_HEALPIX_CXX to "yes" if found.
#
# Also sets HEALPIX_CXX_LIBS to any needed link commands, and HEALPIX_CXX_CFLAGS to any
# needed compiler flags.
#
AC_DEFUN([AMUSE_LIB_HEALPIX_CXX], [
    amuse_lib_healpix_cxx_save_libs="$LIBS"

    AC_MSG_CHECKING([for HEALPix])

    AC_LANG_PUSH([C++])

    LIBS="-lhealpix_cxx"

    AC_LINK_IFELSE([
        AC_LANG_PROGRAM([
            #include <healpix_cxx/healpix_base.h>
        ], [
            Healpix_Base hb;
        ])
    ], [
        FOUND_HEALPIX_CXX="yes"
        HEALPIX_CXX_CFLAGS=""
        HEALPIX_CXX_LIBS="$LIBS"

        AC_MSG_RESULT([yes])
    ],
    [
        # Not found, try pkg-config instead
        PKG_CHECK_MODULES([HEALPIX_CXX], [healpix_cxx], [
            FOUND_HEALPIX_CXX="yes"
            AC_MSG_RESULT([yes])
        ], [
            FOUND_HEALPIX_CXX="no"
            AC_MSG_RESULT([no])
        ])
    ])

    AC_LANG_POP([C++])

    AC_SUBST([HEALPIX_CXX_CFLAGS])
    AC_SUBST([HEALPIX_CXX_LIBS])
])

