# AX_GSL
# ----------------------------------------------------------
# set up for GSL
#
AC_DEFUN([AX_GSL],[
        AC_ARG_WITH(gsl,
             AC_HELP_STRING([--with-gsl=PFX], [Prefix where GSL has been installed] ),
             [
                test "$withval" = no && AC_MSG_WARN([gsl is a required package for some modules])
                test "$withval" = yes || gsl_prefix="$withval" 
                with_gsl=yes ],
             [ with_gsl=yes ] 
        )
    
        AS_IF([test "$with_gsl" = yes ],
        [
            #user override
            AS_IF([test "x$GSL_LIBS" != x && test "x$GSL_FLAGS" != x ],
            [
                have_gsl=yes
                FOUND_GSL=yes
            ],
            [
                saved_LIBS="$LIBS"
                saved_CXXFLAGS="$CXXFLAGS"
                GSL_LIBS=""
                GSL_FLAGS=""
                FOUND_GSL=no
                

                if test x$gsl_prefix != x; then
                    ac_GSL_CFLAGS="-I$gsl_prefix/include"
                    ac_GSL_LDOPTS="-L$gsl_prefix/lib"
                    
                    
                    save_CFLAGS="$CFLAGS"
                    CFLAGS="$ac_GSL_CFLAGS $save_CFLAGS"
                    AC_CHECK_HEADER(
                        [gsl/gsl_sys.h],
                        [GSL_FLAGS="$ac_GSL_CFLAGS",
                        FOUND_GSL="yes"
                        ],
                        [AC_MSG_WARN([Cannot find headers of the library GSL in $gsl_prefix/include.])]
                    )
                    CFLAGS="$save_CFLAGS"
                    
                     save_LIBS="$LIBS"
                    LIBS="$ac_GSL_LDOPTS -lgsl $save_LIBS"
                    AC_CHECK_LIB([gsl], [gsl_blas_dgemm],
                            [GSL_LIBS="$ac_GSL_LDOPTS -lgsl"],
                            [AC_MSG_WARN([libgsl : library missing. (Cannot find symbol gsl_blas_dgemm) in $with_gsl_library. Check if libgsl is installed and if the version is correct])]
                            )
                    LIBS="$save_LIBS"
                    
                else
                    AC_CHECK_HEADER(
                        [gsl/gsl_sys.h],
                        [GSL_FLAGS=""
                        FOUND_GSL="yes"],
                        [AC_MSG_WARN([Cannot find headers of the library GSL, please provide the --with-gsl argument.])]
                    )
                    AC_CHECK_LIB(
                        [gsl],
                        [gsl_blas_dgemm],
                        [GSL_LIBS="-lgsl"],
                        [AC_MSG_WARN([libgsl : library missing. (Cannot find symbol gsl_blas_dgemm) in $with_gsl_library. Check if libgsl is installed and if the version is correct])]
                    )
                fi

                
               
            ])
        ])

        AC_SUBST(FOUND_GSL)
        AC_SUBST(GSL_FLAGS)
        AC_SUBST(GSL_LIBS)
    ]
)
