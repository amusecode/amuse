# Macro for detecting the platform we're on
#
# Sets AMUSE_ON_MACOS if we're on macOS, or AMUSE_ON_LINUX if we're on GNU/Linux.
AC_DEFUN([AMUSE_DETECT_OS], [
    AC_CANONICAL_TARGET
    if test "x$target_os" == "xlinux-gnu"
    then
        AMUSE_ON_LINUX=1
    fi

    if test "x$target_vendor" == "xapple"
    then
        AMUSE_ON_MACOS=1
    fi

    AC_SUBST([AMUSE_ON_LINUX])
    AC_SUBST([AMUSE_ON_MACOS])
])

