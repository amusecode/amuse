AC_INIT([amuse], [3.1])
AC_CONFIG_SRCDIR([build.py])

AC_CONFIG_FILES([src/amuse/config.py:support/configpy.in])
AC_CONFIG_FILES([support/config.py:support/configpy.in])

AC_PROG_CC()
AC_PROG_CXX()
AC_PROG_FC([],95)
AC_PATH_PROG([PYTHON], [python], [])
AC_ARG_VAR([NVCC], [CUDA compiler command])

#
#  Cuda libraries and sdk
#
AC_ARG_ENABLE([cuda],
    [AS_HELP_STRING([--enable-cuda],[Enable CUDA for codes])],
    [WITH_CUDA=yes],
    [WITH_CUDA=no]
)
AC_ARG_VAR([CUDA_SDK], [CUDA sdk directory])
AC_ARG_VAR([CUDA_TK], [CUDA toolkit directory])

AS_IF([test x"$WITH_CUDA" != xno],[
    WITH_CUDA=yes
    AC_PATH_PROG([NVCC], [nvcc],
            [NVCC=no]
            [$PATH$ /usr/local/cuda/cuda/bin /usr/local/cuda/bin]
            [CUDA_TK]
            [/opt/cuda/cuda/bin /opt/cuda/bin])
    
    AS_IF([test x"$CUDA_TK" = x],
    [
        CUDA_TK_BIN=`AS_DIRNAME("$NVCC")`
        CUDA_TK=`AS_DIRNAME("$CUDA_TK_BIN")`
        AS_IF([test x"$NVCC" != xno],
        [],
        [AC_MSG_ERROR([CUDA_TK path is not set, and could not find nvcc in path, please set CUDA_TK variable])]
        )
    ]
    )
    
    AS_IF([test x"$CUDA_SDK" = xno],
        [AC_MSG_ERROR([CUDA_SDK path is not set, please set the CUDA_SDK variable first or disable CUDA])]
    )
    AC_CHECK_FILE([$CUDA_TK/lib], [],
        [AC_MSG_ERROR([cuda toolkit path is incorrect, must have lib directory])], 
        [])
        
    AC_CHECK_FILE([$CUDA_SDK/common/inc/cutil.h], [],
        [AC_MSG_ERROR([cuda sdk path is incorrect, must have common/inc/cutil.h in the CUDA_SDK path])], 
        [])
    
    PREV_LDFLAGS=${LDFLAGS}
    LDFLAGS="-L$CUDA_TK/lib -L$CUDA_TK/lib64"
    AC_CHECK_LIB(
        [cudart],
        [main],
        [],
        [AC_MSG_ERROR([CUDA_TK, cuda toolkit path is incorrect])]
    )
    LDFLAGS=${PREV_LDFLAGS}
])

AC_SUBST(WITH_CUDA)
AC_SUBST(CUDA_SDK)
AC_SUBST(CUDA_TK)

#
# Sapporo 
#

AC_ARG_WITH([sapporo],
    [AS_HELP_STRING([--with-sapporo],[Use sapporo library in given directory])],
    [WITH_SAPPORO=yes],
    [WITH_SAPPORO=no]
)

AC_SUBST(WITH_SAPPORO)


AC_LANG(C++)
AX_MPI()
AC_LANG(C)
AX_MPI()
AC_LANG(Fortran)
AX_MPI()

AC_CHECK_PYTHON_MODULE(numpy, __version__)
AS_IF([test x"$PYTHON_numpy" != x1],
    [AC_MSG_ERROR([numpy not found, please install numpy first])]
)
AS_VERSION_COMPARE([$PYTHON_numpy_VERSION], [1.3.0],
        [AC_MSG_ERROR([numpy version $PYTHON_numpy_VERSION found, but at least 1.3.0 needed])]
)
AC_CHECK_PYTHON_MODULE(mpi4py, __version__)
AS_IF([test x"$PYTHON_mpi4py" != x1],
    [AC_MSG_ERROR([mpi4py not found, please install mpi4py first])]
)
AS_VERSION_COMPARE([$PYTHON_mpi4py_VERSION], [1.1.0],
        [AC_MSG_ERROR([numpy version $PYTHON_mpi4py_VERSION found, but at least 1.1.0 needed])]
)
AC_CHECK_PYTHON_MODULE(h5py, version.version)
AS_VERSION_COMPARE([$PYTHON_h5py_VERSION], [1.1.0],
        [AC_MSG_WARN([h5py version $PYTHON_h5py_VERSION found, but at least 1.1.0 needed])]
)
AC_CHECK_PYTHON_MODULE(nose, __version__)
AS_VERSION_COMPARE([$PYTHON_nose_VERSION], [0.11],
        [AC_MSG_ERROR([nose version $PYTHON_nose_VERSION found, but at least 0.11 needed])]
)
AC_CHECK_PYTHON_MODULE(docutils, __version__)
AS_VERSION_COMPARE([$PYTHON_docutils_VERSION], [0.6],
        [AC_MSG_ERROR([docutils version $PYTHON_docutils_VERSION found, but at least 0.6 needed])]
)
AC_CHECK_PYTHON_MODULE(zlib, __version__)
AS_IF([test x"$PYTHON_zlib" != x1],
    [AC_MSG_ERROR([zlib not found, please recompile python check that development libraries for zlib ar installed first])]
)
$PYTHON -c "import zlib; exit(zlib.crc32('amuse')&0xffffffff == 0xc0cc9367)" 2>/dev/null
AS_IF([test x"$?" != x1],
    [AC_MSG_ERROR([crc32 function of zlib not found, please recompile python])]
)
AS_VERSION_COMPARE([$PYTHON_zlib_VERSION], [1.0],
        [AC_MSG_ERROR([numpy version $PYTHON_zlib_VERSION found, but at least 1.0 needed])]
)
AC_CHECK_PYTHON_MODULE(matplotlib)
AC_OUTPUT()
