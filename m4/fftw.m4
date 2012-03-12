# AX_FFTW
# ----------------------------------------------------------
# set up for FFTW
#
AC_DEFUN([AX_FFTW],[
        AC_ARG_WITH(fftw,
             AC_HELP_STRING([--with-fftw=PFX], [Prefix where FFTW has been installed] ),
             [
                test "$withval" = no && AC_MSG_WARN([fftw is a required package for some modules])
                test "$withval" = yes || fftw_prefix="$withval" 
                with_fftw=yes ],
             [ with_fftw=yes ] 
        )
    
        AS_IF([test "x$with_fftw" != xno ],
        [
            #user override
            AS_IF([test "x$FFTW_LIBS" != x && test "x$FFTW_FLAGS" != x ],
            [
                have_fftw=yes
                FOUND_FFTW="yes"
            ],
            [
                saved_LIBS="$LIBS"
                saved_CXXFLAGS="$CXXFLAGS"
                FFTW_LIBS=""
                FFTW_FLAGS=""
                FOUND_FFTW="no"
                if test x$fftw_prefix == x; then
                   if test x$PREFIX != x; then
		      fftw_prefix=$PREFIX 
                   fi  
                fi
                if test x$fftw_prefix != x; then
                    ac_FFTW_CFLAGS="-I$fftw_prefix/include"
                    ac_FFTW_LDOPTS="-L$fftw_prefix/lib"
                    
                    
                    save_CFLAGS="$CFLAGS"
                    save_CPPFLAGS="$CPPFLAGS"
                    CFLAGS="$ac_FFTW_CFLAGS $save_CFLAGS"
                    CPPFLAGS="$ac_FFTW_CFLAGS $save_CPPFLAGS"
                    AC_CHECK_HEADER(
                        [fftw3.h],
                        [FFTW_FLAGS="$ac_FFTW_CFLAGS"
                        FOUND_FFTW="yes"],
                        [AC_MSG_WARN([Cannot find headers (fftw3.h) of the library FFTW in $fftw_prefix/include.])]
                    )
                    CFLAGS="$save_CFLAGS"
                    CPPFLAGS="$save_CPPFLAGS"
                    
                    save_LIBS="$LIBS"
                    LIBS="$ac_FFTW_LDOPTS -lfftw3  -lfftw3_threads $save_LIBS"
                    
                    cache_var=AS_TR_SH([ac_cv_lib_fftw3_fftw_plan_dft_r2c])
                    AC_CHECK_LIB([fftw3], [fftw_plan_dft_r2c],
                            [FFTW_LIBS="$ac_FFTW_LDOPTS -lfftw3  -lfftw3_threads"],
                            [FOUND_FFTW="no"]
                    )
                    $as_unset $cache_var
                    if test x$FOUND_FFTW != xyes; then
                        LIBS="$ac_FFTW_LDOPTS -lfftw3 -lfftw3_threads -lm $save_LIBS"
                        AC_CHECK_LIB([fftw3], [fftw_plan_dft_r2c],
                            [FFTW_LIBS="$ac_FFTW_LDOPTS -lfftw3 -lfftw3_threads -lm"],
                            [AC_MSG_WARN([libfftw3 : library missing. (Cannot find symbol fftw_plan_dft_r2c) in $fftw_prefix. Check if libfftw3 is installed and if the version is correct])]
                        )
                        $as_unset $cache_var
                    fi
                    LIBS="$save_LIBS"
	fi
    
	if test x$FOUND_FFTW != xyes; then
        PKG_CHECK_MODULES([FFTW],[fftw3 >= 3.2],
            [
            FFTW_FLAGS="$FFTW_CFLAGS"
            FFTW_LIBS="$FFTW_LIBS -lfftw3_threads"
            FOUND_FFTW=yes
            ],
            []
        )
    fi
    
	if test x$FOUND_FFTW != xyes; then
                    AC_CHECK_HEADER(
                        [fftw3.h],
                        [FFTW_FLAGS=""
                        FOUND_FFTW="yes"],
                        [AC_MSG_WARN([Cannot find headers (fftw3.h) of the library FFTW in $fftw_prefix/include.])]
                    )
                    cache_var=AS_TR_SH([ac_cv_lib_fftw3_fftw_plan_dft_r2c])
                    AC_CHECK_LIB(
                        [fftw3],
                        [fftw_plan_dft_r2c],
                        [FFTW_LIBS="-lfftw3 -lfftw3_threads"],
                        [FOUND_FFTW="no"]
                    )
                    $as_unset $cache_var
                    
                    if test x$FOUND_FFTW != xyes; then
                    
                        save_LIBS="$LIBS"
                        LIBS="-lfftw3_threads -lm $save_LIBS"
                        AC_CHECK_LIB([fftw3], [fftw_plan_dft_r2c],
                                [FFTW_LIBS="$ac_FFTW_LDOPTS -lfftw3  -lfftw3_threads -lm"],
                                [FOUND_FFTW="no"
                                AC_MSG_WARN([libfftw3 : library missing. (Cannot find symbol fftw_plan_dft_r2c) in $fftw_prefix. Check if libfftw3 is installed and if the version is correct])]
                        )
                        LIBS="$save_LIBS"
                    fi
    fi
                
               
            ])
        ])

        AC_SUBST(FOUND_FFTW)
        AC_SUBST(FFTW_FLAGS)
        AC_SUBST(FFTW_LIBS)
    ]
)
