# AX_CUDA_VERIFY_HEADERS(PATH)
#
# Checks that the given path contains a CUDA installation with include/cuda.h.
#
# Sets ax_cuda_verify_FLAGS to -I/path/to/cuda/include if successful.
AC_DEFUN([AX_CUDA_VERIFY_HEADERS], [
    AS_IF([test "x$1" = x], [
        ax_cuda_verify_msg="for cuda.h"
        ax_cuda_flags=
    ], [
        ax_cuda_verify_msg="for cuda.h in $1/include"
        ax_cuda_flags="-I$1/include"
    ])

    AC_MSG_CHECKING([$ax_cuda_verify_msg])
    ax_cuda_verify_FLAGS=

    ax_save_CFLAGS="$CFLAGS"

    CFLAGS="$ax_cuda_flags"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([#include <cuda.h>])], [
        ax_cuda_verify_FLAGS="$CFLAGS"
        ax_cuda_verify_headers_found="yes"
        AC_MSG_RESULT([yes])
    ], [
        AC_MSG_RESULT([no])
    ])

    CFLAGS="$ax_save_CFLAGS"
])

# AX_CUDA_VERIFY_LIBS(PATH)
#
# Checks that the given path contains a CUDA installation with lib/libcudart.so or
# lib64/libcudart.so.
#
# Sets ax_cuda_verify_LDFLAGS to -L/path/to/cuda/lib{64} if successful.
AC_DEFUN([AX_CUDA_VERIFY_LIBS], [
    AS_IF([test "x$1" = x], [
        ax_cuda_verify_msg="for libcudart"
        ax_cuda_ldflags=
    ], [
        ax_cuda_verify_msg="for libcudart in $1/lib"
        ax_cuda_verify_msg64="for libcudart in $1/lib64"
        ax_cuda_ldflags="-L$1/lib"
        ax_cuda_ldflags64="-L$1/lib64"
    ])

    AC_MSG_CHECKING([$ax_cuda_verify_msg])
    ax_cuda_verify_LDFLAGS=

    ax_save_LDFLAGS="$LDFLAGS"
    ax_save_LIBS="$LIBS"

    LDFLAGS="$ax_cuda_ldflags"
    LIBS="-lcudart"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <cuda.h>], [cudaFree(0);])], [
        ax_cuda_verify_LDFLAGS="$LDFLAGS"
        ax_cuda_verify_libs_found="yes"
        AC_MSG_RESULT([yes])
    ], [
        AC_MSG_RESULT([no])
    ])

    AS_IF([test "x$ax_cuda_verify_libs_found" = x], [
        AC_MSG_CHECKING([$ax_cuda_verify_msg64])
        LDFLAGS="$ax_cuda_ldflags64"
        AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <cuda.h>], [cudaFree(0);])], [
            ax_cuda_verify_LDFLAGS="$LDFLAGS"
            ax_cuda_verify_libs_found="yes"
            AC_MSG_RESULT([yes])
        ], [
            AC_MSG_RESULT([no])
        ])
    ])

    LIBS="$ax_save_LIBS"
    LDFLAGS="$ax_save_LDFLAGS"
])

# AX_CUDA_VERIFY_SET_VARS(CUDA_TK)
#
# Sets CUDA_TK, NVCC, CUDA_FLAGS and CUDA_LDFLAGS from the internal values, if all
# of them have been found.
AC_DEFUN([AX_CUDA_VERIFY_SET_VARS], [
    AS_IF(
        [test \(x"$ax_cuda_verify_NVCC" != x\) -a \(x"$ax_cuda_verify_headers_found" != x\) -a \(x"$ax_cuda_verify_libs_found" != x\)],
        [
            CUDA_TK="$1"
            NVCC="$ax_cuda_verify_NVCC"
            CUDA_FLAGS="$ax_cuda_verify_FLAGS"
            CUDA_LDFLAGS="$ax_cuda_verify_LDFLAGS"
        ]
    )
])

# AX_CUDA_VERIFY(PATH)
#
# Checks that the given path contains a CUDA installation with bin/nvcc,
# lib/libcudart.so, and include/cuda.h.
#
# Sets CUDA_TK=<path> if successful, and if successful also sets NVCC to the location
# of nvcc, CUDA_LDFLAGS to -L/path/to/cuda/libs, and CUDA_FLAGS to
# -I/path/to/cuda/include.
AC_DEFUN([AX_CUDA_VERIFY], [
    ax_cuda_verify_headers_found=
    ax_cuda_verify_libs_found=

    AS_IF([test -d $1], [
        AC_PATH_PROG([ax_cuda_verify_NVCC], [nvcc], [], [$1/bin])

        AS_IF([test x"$ax_cuda_verify_NVCC" != x], [
            AX_CUDA_VERIFY_HEADERS([$1])
        ])

        AS_IF([test x"$ax_cuda_verify_headers_found" != x], [
            AX_CUDA_VERIFY_LIBS([$1])
        ])

        AX_CUDA_VERIFY_SET_VARS([$1])
    ], [
        AC_MSG_NOTICE([$1 does not exist or is not a directory, CUDA not found there])
    ])
])

# AX_CUDA_VERIFY_DEFAULT()
#
# Checks for NVCC on the PATH, then finds the CUDA directory from there.
#
# Sets CUDA_TK=<path> if successful, and if successful also sets NVCC to the location
# of nvcc, CUDA_LDFLAGS to -L/path/to/cuda/libs, and CUDA_FLAGS to
# -I/path/to/cuda/include.
AC_DEFUN([AX_CUDA_VERIFY_DEFAULT], [
    AC_PATH_PROG([ax_cuda_verify_NVCC], [nvcc])

    AS_IF([test x"$ax_cuda_verify_NVCC" != x], [
        # Got nvcc, verify that we have the rest too
        ax_cuda_verify_ctk_rel=$(dirname -- "$ax_cuda_verify_NVCC")/..
        # Canonicalise path portably so that it looks nicer
        ax_cuda_verify_CUDA_TK=$(test -d "$ax_cuda_verify_ctk_rel" && CDPATH= cd -P -- "$ax_cuda_verify_ctk_rel" && pwd -P)

        AS_IF([test x"$ax_cuda_verify_NVCC" != x], [
            AX_CUDA_VERIFY_HEADERS()
        ])

        AS_IF([test x"$ax_cuda_verify_headers_found" != x], [
            AX_CUDA_VERIFY_LIBS()
        ])

        AX_CUDA_VERIFY_SET_VARS([$ax_cuda_verify_CUDA_TK])
    ], [
        AC_MSG_NOTICE([Could not find CUDA via PATH])
    ])
])

# AX_CUDA
# ----------------------------------------------------------
# set up for CUDA
#
# will set:
#
# CUDA_TK
# NVCC
#
# CUDA_LDFLAGS  (with -L/path/to/cuda/lib)
# CUDA_FLAGS    (with -I/path/to/cuda/include)
#
# and call AC_SUBST on them.
#
AC_DEFUN([AX_CUDA], [
    AC_LANG_PUSH([C])
    AS_IF([test x"$CUDA_TK" != x], [
        # User set CUDA_TK, verify and use it if it works
        AX_CUDA_VERIFY($CUDA_TK)

        AS_IF([test x"$NVCC" = x], [
            AC_MSG_ERROR([CUDA_TK is set, but there is no nvcc in $CUDA_TK/bin. Please set CUDA_TK to the CUDA installation directory, or unset it to autodetect.])
        ])
    ], [
        # CUDA_TK not set, try to discover CUDA via PATH
        AX_CUDA_VERIFY_DEFAULT()

        # Not in PATH, try default directory
        AS_IF([test x"$CUDA_TK" = x], [
            AX_CUDA_VERIFY([/usr/local/cuda])
        ])

        # Try /usr/local/cuda-#.#, but only if there's exactly one match
        AS_IF([test x"$CUDA_TK" = x], [
            # List CUDA installations
            ax_cuda_installs=$(ls -1 -d /usr/local/cuda-*.* 2>/dev/null)
            # If there's more than one, the above will have newlines, and change when we do this
            ax_cuda_installs2=$(echo "x${ax_cuda_installs}" | tr -d '[:space:]')
            AS_IF([test "x${ax_cuda_installs}" != x], [
                AS_IF([test "x${ax_cuda_installs}" = "${ax_cuda_installs2}"], [
                    # Here, there's exactly one match
                    AX_CUDA_VERIFY([$ax_cuda_installs])
                ], [
                    AC_MSG_NOTICE([Multiple CUDA installations found in /usr/local, without a /usr/local/cuda symlink. Please set CUDA_TK to specify which one you want to use.])
                ])
            ], [
                AC_MSG_NOTICE([No CUDA installations found in /usr/local])
            ])
        ])

        # Try some other locations
        for dir in /opt/cuda /usr/local/cuda/cuda /opt/cuda/cuda ; do
            AS_IF([test x"$CUDA_TK" = x], [
                AX_CUDA_VERIFY([$dir])
            ])
        done

        # These directories were checked by the old macro, but they're highly obsolete
        # so we're not trying them anymore:
        # /usr/lib/nvidia, /usr/include/nvidia, /usr/lib/nvidia-current, /usr/include/nvidia-current

        AS_IF([test x"$CUDA_TK" = x], [
            AC_MSG_NOTICE([CUDA not found. Set CUDA_TK if you have it in an odd location.])
        ])
    ])

    AS_IF([test x"$CUDA_TK" != x], [
        AC_MSG_NOTICE([CUDA found at $CUDA_TK])
        AC_SUBST(CUDA_TK)
        AC_SUBST(NVCC)
        AC_SUBST(CUDA_FLAGS)
        AC_SUBST(CUDA_LDFLAGS)
    ])

    AC_LANG_POP([C])
])

