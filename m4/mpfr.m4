# AX_MPFR
# ----------------------------------------------------------
# set up for MPFR
#
AC_DEFUN([AX_MPFR],[
        AC_ARG_WITH(mpfr,
             AC_HELP_STRING([--with-mpfr=PFX], [Prefix where MPFR has been installed] ),
             [
                test "$withval" = no && AC_MSG_WARN([mpfr is a required package for some modules])
                test "$withval" = yes || mpfr_prefix="$withval" 
                with_mpfr=yes ],
             [ with_mpfr=yes ] 
        )
    
        AS_IF([test "x$with_mpfr" != xno ],
        [
            #user override
            AS_IF([test "x$MPFR_LIBS" != x ],
            [
                have_mpfr=yes
                FOUND_MPFR="yes"
            ],
            [
                
                
                
                saved_LIBS="$LIBS"
                saved_CXXFLAGS="$CXXFLAGS"
                MPFR_LIBS=""
                MPFR_FLAGS=""
                FOUND_MPFR="no"
                if test x$mpfr_prefix == x; then
                   if test x$PREFIX != x; then
		                mpfr_prefix=$PREFIX 
                   fi  
                fi
                if test x$mpfr_prefix != x; then
                    ac_mpfr_CFLAGS="-I$mpfr_prefix/include"
                    ac_mpfr_LDOPTS="-L$mpfr_prefix/lib"
                    
                    
                    save_CFLAGS="$CFLAGS"
                    save_CPPFLAGS="$CPPFLAGS"
                    CFLAGS="$ac_mpfr_CFLAGS $save_CFLAGS"
                    CPPFLAGS="$ac_mpfr_CFLAGS $save_CPPFLAGS"
                    
                    cache_var=AS_TR_SH([ac_cv_header_mpfr.h])
                    AC_CHECK_HEADER(
                        [mpfr.h],
                        [MPFR_FLAGS="$ac_mpfr_CFLAGS"
                        FOUND_MPFR="yes"],
                        [AC_MSG_WARN([Cannot find headers (mpfr.h) of the library MPFR in $mpfr_prefix/include.])]
                    )
                    $as_unset $cache_var
                    CFLAGS="$save_CFLAGS"
                    CPPFLAGS="$save_CPPFLAGS"
                    
	                if test x$FOUND_MPFR == xyes; then
                        save_LIBS="$LIBS"
                        LIBS="$ac_mpfr_LDOPTS -lmpfr $save_LIBS"
                        
                        cache_var=AS_TR_SH([ac_cv_lib_mpfr_mpfr_init])
                        AC_CHECK_LIB([mpfr], [mpfr_init],
                            [MPFR_LIBS="$ac_mpfr_LDOPTS -lmpfr"],
                            [FOUND_MPFR="no"
                            AC_MSG_WARN([libmpfr : library missing. (Cannot find symbol mpfr_init) in $mpfr_prefix. Check if libmpfr is installed and if the version is correct])]
                        )
                        LIBS="$save_LIBS"
                        $as_unset $cache_var
                    fi
	fi
	if test x$FOUND_MPFR != xyes; then
	    PKG_CHECK_MODULES([MPFR],[mpfr >= 3],
                    [
                    MPFR_FLAGS="$MPFR_FLAGS"
                    MPFR_LIBS="$MPFR_LIBS"
                    FOUND_MPFR=yes
                    ],
                    []
                )
            
    fi
	
	if test x$FOUND_MPFR != xyes; then
                    AC_CHECK_HEADER(
                        [mpfr.h],
                        [MPFR_FLAGS=""
                        FOUND_MPFR="yes"],
                        [AC_MSG_WARN([Cannot find headers (mpfr.h) of the library MPFR.])]
                    )
                    AC_CHECK_LIB(
                        [mpfr],
                        [mpfr_init],
                        [MPFR_LIBS="-lmpfr"],
                        [
                        FOUND_MPFR="no"
                        AC_MSG_WARN([libmpfr : library missing. (Cannot find symbol mpfr_init). Check if libmpfr is installed and if the version is correct])]
                    )
		fi
                
               
            ])
        ])

        AC_SUBST(FOUND_MPFR)
        AC_SUBST(MPFR_FLAGS)
        AC_SUBST(MPFR_LIBS)
    ]
)
