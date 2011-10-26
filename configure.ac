AC_INIT([amuse], [3.1])
AC_CONFIG_SRCDIR([src/amuse/__init__.py])

AC_CONFIG_FILES([src/amuse/config.py:support/configpy.in])
AC_CONFIG_FILES([support/config.py:support/configpy.in])
AC_CONFIG_FILES([config.mk])
AC_CONFIG_FILES([build.py],[chmod +x build.py])

AC_PROG_CC()
AC_PROG_CXX()
AC_PROG_FC([],95)
AC_ARG_VAR([PYTHON], [Python interpreter])


AS_IF([test "x$PYTHON" = "x"], [PYTHON=python])

AC_PATH_PROG([PYTHON], [$PYTHON], [])
AC_ARG_VAR([NVCC], [CUDA compiler command])


AC_ARG_VAR([JAVA_HOME], [Java sdk directory])

AMUSE_CHECK_JDK()

AC_SUBST(JDK)
AC_SUBST(JNI_INCLUDES)


AC_ARG_ENABLE(mpi,
    [AS_HELP_STRING([--enable-mpi],[Enable MPI support for communicating with codes, by default enabled])],
    [WITH_MPI=no
    AS_IF([test "x$enable_mpi" != "xno"], [
	    WITH_MPI=yes
    ])
    ],
    [WITH_MPI=yes]
)

AC_SUBST(WITH_MPI)

AS_IF([test x"$WITH_MPI" != xno],[
    AC_LANG(C++)
    AX_MPI()
    AC_LANG(C)
    AX_MPI()
    AC_LANG(Fortran)
    AX_MPI()
], [
    MPICC="$CC"
    MPICXX="$CXX"
    MPIFC="$FC"
    
    CFLAGS="$CFLAGS -DNOMPI"
    CXXFLAGS="$CXXFLAGS -DNOMPI"
    FCFLAGS="$FCFLAGS -DNOMPI"
    
    AC_SUBST(MPICC)
    AC_SUBST(MPICXX)
    AC_SUBST(MPIFC)
]    
)
AC_LANG(C)
AX_OPENMP(
     [HAVE_OPENMP="yes"]
    ,[HAVE_OPENMP="no"]
)
AC_LANG(Fortran)
AX_OPENMP()
    
AC_SUBST(HAVE_OPENMP)
AC_SUBST(OPENMP_CFLAGS)
AC_SUBST(OPENMP_FCFLAGS)

AC_LANG(C)

AX_CUDA()

#
# Sapporo 
#

AC_ARG_WITH([sapporo],
    [AS_HELP_STRING([--with-sapporo],[Use sapporo library in given directory])],
    [WITH_SAPPORO=yes],
    [WITH_SAPPORO=no]
)

AC_ARG_ENABLE(sapporo2,
    [AS_HELP_STRING([--enable-sapporo2],[Use the provided sapporo library version 2, instead of sapporo_light])],
    [SAPPORO_VERSION="light"
    AS_IF([test "x$enable_sapporo2" != "xno"], [
	SAPPORO_VERSION="2"
    ])
    ],
    [SAPPORO_VERSION="light"]
)

AC_SUBST(WITH_SAPPORO)
AC_SUBST(SAPPORO_VERSION)

PKG_CHECK_MODULES([FFTW],[fftw3 >= 3.2],
    [
    FFTW_FLAGS="$FFTW_CFLAGS"
    FFTW_LIBS="$FFTW_LIBS -lfftw3_threads"
    FOUND_FFTW=yes
    AC_SUBST(FFTW_FLAGS)
    AC_SUBST(FFTW_LIBS)
    ],
    [AX_FFTW()]
)
AC_SUBST(FOUND_FFTW)

PKG_CHECK_MODULES([GSL],[gsl >= 1.0],
    [
    GSL_FLAGS="$GSL_CFLAGS"
    FOUND_GSL=yes
    AC_SUBST(GSL_FLAGS)
    ],
    [
    AX_PATH_GSL(1.0, [FOUND_GSL=yes], [
	FOUND_GSL=no
	AC_MSG_WARN([gsl not found, please specify --with-gsl])
    ])
])

AC_SUBST(FOUND_GSL)

AC_CHECK_PYTHON_MODULE(numpy, __version__)
AS_IF([test x"$PYTHON_numpy" != x1],
    [AC_MSG_ERROR([numpy not found, please install numpy first])]
)
AS_VERSION_COMPARE([$PYTHON_numpy_VERSION], [1.3.0],
        [AC_MSG_ERROR([numpy version $PYTHON_numpy_VERSION found, but at least 1.3.0 needed])]
)


AS_IF([test x"$WITH_MPI" != xno],[
    AC_CHECK_PYTHON_MODULE(mpi4py, __version__)
    AS_IF([test x"$PYTHON_mpi4py" != x1],
        [AC_MSG_ERROR([mpi4py not found, please install mpi4py first])]
    )

    AS_VERSION_COMPARE([$PYTHON_mpi4py_VERSION], [1.1],
            [AC_MSG_ERROR([numpy version $PYTHON_mpi4py_VERSION found, but at least 1.1.0 needed])]
    )
])


AC_CHECK_PYTHON_MODULE(h5py, version.version)
AS_VERSION_COMPARE([$PYTHON_h5py_VERSION], [1.1.0],
        [AC_MSG_WARN([h5py version $PYTHON_h5py_VERSION found, but at least 1.1.0 needed])]
)
AC_CHECK_PYTHON_MODULE(nose, __version__)
AS_VERSION_COMPARE([$PYTHON_nose_VERSION], [0.11],
        [AC_MSG_ERROR([nose version $PYTHON_nose_VERSION found, but at least 0.11 needed])]
)
AC_CHECK_PYTHON_MODULE(docutils, __version__)
AS_VERSION_COMPARE([$PYTHON_docutils_VERSION], [0.5],
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
