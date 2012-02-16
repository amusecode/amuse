# AX_GMP
# ----------------------------------------------------------
# set up for GMP
#
AC_DEFUN([AX_GMP],[
        AC_ARG_WITH(gmp,
             AC_HELP_STRING([--with-gmp=PFX], [Prefix where GMP has been installed] ),
             [
                test "$withval" = no && AC_MSG_WARN([gmp is a required package for some modules])
                test "$withval" = yes || gmp_prefix="$withval" 
                with_gmp=yes ],
             [ with_gmp=yes ] 
        )
    
        AS_IF([test "x$with_gmp" != xno ],
        [
            #user override
            AS_IF([test "x$GMP_LIBS" != x ],
            [
                have_gmp=yes
                FOUND_GMP="yes"
            ],
            [
                
                
                
                saved_LIBS="$LIBS"
                saved_CXXFLAGS="$CXXFLAGS"
                GMP_LIBS=""
                GMP_FLAGS=""
                FOUND_GMP="no"
                if test x$gmp_prefix == x; then
                   if test x$PREFIX != x; then
		                gmp_prefix=$PREFIX 
                   fi  
                fi
                if test x$gmp_prefix != x; then
                    ac_gmp_CFLAGS="-I$gmp_prefix/include"
                    ac_gmp_LDOPTS="-L$gmp_prefix/lib"
                    
                    
                    save_CFLAGS="$CFLAGS"
                    save_CPPFLAGS="$CPPFLAGS"
                    CFLAGS="$ac_gmp_CFLAGS $save_CFLAGS"
                    CPPFLAGS="$ac_gmp_CFLAGS $save_CPPFLAGS"
                    AC_CHECK_HEADER(
                        [gmp.h],
                        [GMP_FLAGS="$ac_gmp_CFLAGS"
                        FOUND_GMP="yes"],
                        [AC_MSG_WARN([Cannot find headers (gmp.h) of the library GMP in $gmp_prefix/include.])]
                    )
                    CFLAGS="$save_CFLAGS"
                    CPPFLAGS="$save_CPPFLAGS"
                    
	                if test x$FOUND_GMP == xyes; then
                        save_LIBS="$LIBS"
                        LIBS="$ac_gmp_LDOPTS -lgmp $save_LIBS"
                        AC_CHECK_LIB([gmp], [__gmpz_init],
                            [GMP_LIBS="$ac_gmp_LDOPTS -lgmp"],
                            [FOUND_GMP="no"
                            AC_MSG_WARN([libgmp : library missing. (Cannot find symbol gmpz_init) in $gmp_prefix. Check if libgmp is installed and if the version is correct])]
                        )
                        LIBS="$save_LIBS"
                    fi
	fi
	if test x$FOUND_GMP != xyes; then
	    PKG_CHECK_MODULES([GMP],[gmp >= 3],
                    [
                    GMP_FLAGS="$GMP_FLAGS"
                    GMP_LIBS="$GMP_LIBS"
                    FOUND_GMP=yes
                    ],
                    []
                )
            
    fi
	
	if test x$FOUND_GMP != xyes; then
                    AC_CHECK_HEADER(
                        [gmp.h],
                        [GMP_FLAGS=""
                        FOUND_GMP="yes"],
                        [AC_MSG_WARN([Cannot find headers (gmp.h) of the library GMP.])]
                    )
                    AC_CHECK_LIB(
                        [gmp],
                        [__gmpz_init],
                        [GMP_LIBS="-lgmp"],
                        [
                        FOUND_GMP="no"
                        AC_MSG_WARN([libgmp : library missing. (Cannot find symbol gmpz_init). Check if libgmp is installed and if the version is correct])]
                    )
		fi
                
               
            ])
        ])

        AC_SUBST(FOUND_GMP)
        AC_SUBST(GMP_FLAGS)
        AC_SUBST(GMP_LIBS)
    ]
)
