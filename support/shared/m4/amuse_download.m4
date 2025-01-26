# Helper macros for detecting download tools
#
# Some of the community codes aren't included with AMUSE, but are downloaded at build
# time. We need a tool for that, and here is where we find one.
#
# AMUSE_DOWNLOAD()
#
# Searches for a download tool. Two popular and commonly available command line
# download tools are curl and wget. This macro tries to find either one of them
# and sets DOWNLOAD to a command that will take a URL and download its contents to
# standard output. This makes it easier to write that output to a file with a known
# name, which you usually want in a Makefile.
#
# To download a file, use $(DOWNLOAD) https://example.com >example.html
#
AC_DEFUN([AMUSE_DOWNLOAD], [
    AC_CHECK_TOOL(WGET, wget)
    AC_CHECK_TOOL(CURL, curl)

    AC_MSG_CHECKING([for a wget or curl to download files with])
    if test "x$WGET" != "x"
    then
        # The MESA SDK server rejects wget, this is the official work-around
        DOWNLOAD="$WGET --progress=bar:force:noscroll --user-agent='' -O -"
        AC_MSG_RESULT([yes])
    else
        if test "x$CURL" != "x"
        then
            DOWNLOAD="$CURL -L"
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no])
        fi
    fi

    AC_SUBST([DOWNLOAD])
])

