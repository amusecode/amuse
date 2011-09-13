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
    
        AS_IF([test "$with_fftw" = yes ],
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

                if test x$fftw_prefix != x; then
                    ac_FFTW_CFLAGS="-I$fftw_prefix/include"
                    ac_FFTW_LDOPTS="-L$fftw_prefix/lib"
                    
                    
                    save_CFLAGS="$CFLAGS"
                    CFLAGS="$ac_FFTW_CFLAGS $save_CFLAGS"
                    AC_CHECK_HEADER(
                        [fftw3.h],
                        [FFTW_FLAGS="$ac_FFTW_CFLAGS"
                        FOUND_FFTW="yes"],
                        [AC_MSG_WARN([Cannot find headers (fftw3.h) of the library FFTW in $fftw_prefix/include.])]
                    )
                    CFLAGS="$save_CFLAGS"
                    
                     save_LIBS="$LIBS"
                    LIBS="-L$ac_FFTW_LDOPTS -lfftw3 $save_LIBS"
                    AC_CHECK_LIB([fftw3], [fftw_plan_dft_r2c],
                            [FFTW_LIBS="$ac_FFTW_LDOPTS -lfftw3"],
                            [AC_MSG_WARN([libfftw3 : library missing. (Cannot find symbol fftw_plan_dft_r2c) in $with_fftw_library. Check if libfftw3 is installed and if the version is correct])]
                            )
                    LIBS="$save_LIBS"
                    
                else
                    AC_CHECK_HEADER(
                        [fftw3.h],
                        [FFTW_FLAGS=""
                        FOUND_FFTW="yes"],
                        [AC_MSG_WARN([Cannot find headers (fftw3.h) of the library FFTW in $fftw_prefix/include.])]
                    )
                    AC_CHECK_LIB(
                        [fftw3],
                        [fftw_plan_dft_r2c],
                        [FFTW_LIBS="-lfftw3"],
                        [AC_MSG_WARN([libfftw3 : library missing. (Cannot find symbol fftw_plan_dft_r2c) in $with_fftw_library. Check if libfftw3 is installed and if the version is correct])]
                    )
                fi

                
               
            ])
        ])

        AC_SUBST(FOUND_FFTW)
        AC_SUBST(FFTW_FLAGS)
        AC_SUBST(FFTW_LIBS)
    ]
)
