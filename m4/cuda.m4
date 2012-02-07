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
        [WITH_CUDA=no
        CUDA_TK=/NOCUDACONFIGURED
        CUDA_LIBS="-L$CUDA_TK cuda cudart"	
	]
    )
    AC_ARG_WITH(
        cuda-libdir, [  --with-cuda-libdir=PFX   Directory where libcuda.so is installed (optional)],
            cuda_libdir="$withval", cuda_libdir=".")
            
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
        
        
        AC_CHECK_FILE([$CUDA_TK/lib], [],
        [AC_MSG_ERROR([cuda toolkit path is incorrect, must have lib directory])], 
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
            [$cuda_libdir /usr/lib /usr/lib64 /usr/lib/nvidia /usr/lib64/nvidia /usr/lib/nvidia-current /usr/lib64/nvidia-current $CUDA_TK/lib $CUDA_TK/lib64],
            [CUDA_LIBS="$LIBS $CUDART_LIBS"],
            [AC_MSG_ERROR([cannot find cuda library])]
        )
        LIBS="$save_LIBS"
       
    ])

    AC_SUBST(WITH_CUDA)
    
    AC_SUBST(CUDA_TK)
    
    AC_SUBST(CUDA_LIBS)

])
