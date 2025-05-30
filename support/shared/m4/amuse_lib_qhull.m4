# Helper macro for detecting libqhull
#
# AMUSE_LIB_QHULL()
#
# Searches for libqhull and sets FOUND_QHULL to "yes" if found.
#
# Also sets QHULL_LIBS to any needed link commands, and QHULL_FLAGS to any needed
# compiler flags.
#
# Recent versions of libqhull have a reentrant libqhull_r and also a libqhull_p, but we
# test only for libqhull.
#
AC_DEFUN([AMUSE_LIB_QHULL], [
    amuse_lib_qhull_save_libs="$LIBS"

    LIBS=""
    AC_SEARCH_LIBS([qh_freeqhull], [qhull_r], [
        FOUND_QHULL="yes"
    ], [])
    QHULL_LIBS="$LIBS"
    QHULL_FLAGS=""

    LIBS="$amuse_lib_qhull_save_libs"

    AC_CHECK_HEADERS([libqhull_r/libqhull_r.h], [], [FOUND_QHULL="no"])

    AC_SUBST([QHULL_LIBS])
    AC_SUBST([QHULL_FLAGS])
])

