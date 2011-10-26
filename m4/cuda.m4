# AX_CUDA
# ----------------------------------------------------------
# set up for CUDA
#
# will set:
#
# CUDA_SDK
# CUDA_TK
#
# CUDA_LIBS
# CUDA_FLAGS
#
# CUDART_LIBS
# CUDART_FLAGS
#
AC_DEFUN([AX_CUDA],[
    AC_ARG_ENABLE(cuda,
        [AS_HELP_STRING([--enable-cuda],[Enable CUDA for codes])],
        [WITH_CUDA=no
        AS_IF([test "x$enable_cuda" != "xno"], [
        WITH_CUDA=yes
        ])
        ],
        [WITH_CUDA=no]
    )
    AC_ARG_VAR([CUDA_SDK], [CUDA sdk directory])
    AC_ARG_VAR([CUDA_TK], [CUDA toolkit directory])

    AS_IF([test x"$WITH_CUDA" != xno],[
        WITH_CUDA=yes
        AC_PATH_PROG(
        [NVCC], 
        nvcc,
        [no],
        $PATH:/usr/local/cuda/cuda/bin:/usr/local/cuda/bin:$CUDA_TK/bin:/opt/cuda/cuda/bin:/opt/cuda/bin
        )
        
        
        AS_IF([test x"$CUDA_TK" = x],
        [
        
        AS_IF([test x"$NVCC" = xno],
            [AC_MSG_ERROR([CUDA_TK path is not set, and could not find nvcc in path, please set CUDA_TK variable])]
        )
        
        CUDA_TK_BIN=`AS_DIRNAME("$NVCC")`
        CUDA_TK=`AS_DIRNAME("$CUDA_TK_BIN")`
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
        
        
        save_LIBS="$LIBS"
        CUDART_LIBS=''
        AC_FIND_LIB(
            [cudart],
            [main],
            [$CUDA_TK/lib $CUDA_TK/lib64],
            [CUDART_LIBS=$LIBS],
            [AC_MSG_ERROR([cannot find cuda runtime libraries in $CUDA_TK/lib $CUDA_TK/lib64])]
        )
        LIBS="$save_LIBS"
        
        save_LIBS="$LIBS"
        CUDA_PATH=[]
        CUDA_LIBS=''
        AC_FIND_LIB(
            [cuda],
            [main],
            [/usr/lib /usr/lib64 /usr/lib/nvidia /usr/lib64/nvidia  $CUDA_TK/lib $CUDA_TK/lib64],
            [CUDA_LIBS="$LIBS $CUDART_LIBS"],
            [AC_MSG_ERROR([cannot find cuda library])]
        )
        LIBS="$save_LIBS"
       
    ])

    AC_SUBST(WITH_CUDA)
    
    AC_SUBST(CUDA_SDK)
    AC_SUBST(CUDA_TK)
    
    AC_SUBST(CUDA_LIBS)

])
