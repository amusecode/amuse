# Configure path for the GNU Scientific Library
# Christopher R. Gabriel <cgabriel@linux.it>, April 2000


AC_DEFUN([AX_PATH_GSL],
[
AC_ARG_WITH(gsl-prefix,[  --with-gsl-prefix=PFX   Prefix where GSL is installed (optional)],
            gsl_prefix="$withval", gsl_prefix="")
AC_ARG_WITH(gsl-exec-prefix,[  --with-gsl-exec-prefix=PFX Exec prefix where GSL is installed (optional)],
            gsl_exec_prefix="$withval", gsl_exec_prefix="")
AC_ARG_ENABLE(gsltest, [  --disable-gsltest       Do not try to compile and run a test GSL program],
		    , enable_gsltest=yes)

  if test "x${GSL_CONFIG+set}" != xset ; then
     if test "x$gsl_prefix" != x ; then
         GSL_CONFIG="$gsl_prefix/bin/gsl-config"
     fi
     if test "x$gsl_exec_prefix" != x ; then
        GSL_CONFIG="$gsl_exec_prefix/bin/gsl-config"
     fi
  fi
  
  
  if test "x${GSL_CONFIG+set}" == xset ; then
    if test -x "${GSL_CONFIG}" ; then
        gsl_config_exists=yes
    else
        gsl_config_exists=no
    fi
  fi
  FOUND_GSL=no
  if test "x$GSL_CFLAGS" != x ; then
    if test "x$GSL_LIBS" != x ; then
        GSL_FLAGS="$GSL_CFLAGS"
        GSL_LIBS="$GSL_LIBS"
        FOUND_GSL=yes
        AC_SUBST(GSL_FLAGS)
    fi
  fi
  if test "$FOUND_GSL" = "no"; then
      if test "x$gsl_config_exists" != xyes ; then
        PKG_CHECK_MODULES([GSL],[gsl >= 1.0],
            [
            GSL_FLAGS="$GSL_CFLAGS"
            GSL_LIBS="$GSL_LIBS"
            FOUND_GSL=yes
            AC_SUBST(GSL_FLAGS)
            
            ifelse([$2], , :, [$2]) 
            ],
            [
            ]
        )
      fi
  fi
  if test "$FOUND_GSL" = "no"; then
      AC_PATH_PROG(GSL_CONFIG, gsl-config, no)
      min_gsl_version=ifelse([$1], ,0.2.5,$1)
      AC_MSG_CHECKING(for GSL - version >= $min_gsl_version)
      no_gsl=""
      if test "$GSL_CONFIG" = "no" ; then
        no_gsl=yes
      else
        GSL_FLAGS=`$GSL_CONFIG --cflags`
        GSL_LIBS=`$GSL_CONFIG --libs`

        gsl_version=`$GSL_CONFIG --version`
        AS_VERSION_COMPARE(
            [$gsl_version], 
            [$min_gsl_version],
            [no_gsl=yes],
            
        )
      fi
      if test "x$no_gsl" = x ; then
         AC_MSG_RESULT(yes)
         ifelse([$2], , :, [$2])     
      else
         AC_MSG_RESULT(no)
         if test "$GSL_CONFIG" = "no" ; then
           echo "*** The gsl-config script installed by GSL could not be found"
           echo "*** If GSL was installed in PREFIX, make sure PREFIX/bin is in"
           echo "*** your path, or set the GSL_CONFIG environment variable to the"
           echo "*** full path to gsl-config."
         else
           if test -f conf.gsltest ; then
            :
           else
              echo "*** Could not run GSL test program, checking why..."
              CFLAGS="$CFLAGS $GSL_FLAGS"
              LIBS="$LIBS $GSL_LIBS"
              AC_LINK_IFELSE([AC_LANG_PROGRAM([[
    #include <stdio.h>
    ]], [[ return 0; ]])],[ echo "*** The test program compiled, but did not run. This usually means"
              echo "*** that the run-time linker is not finding GSL or finding the wrong"
              echo "*** version of GSL. If it is not finding GSL, you'll need to set your"
              echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
              echo "*** to the installed location  Also, make sure you have run ldconfig if that"
              echo "*** is required on your system"
          echo "***"
              echo "*** If you have an old version installed, it is best to remove it, although"
              echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH"],[ echo "*** The test program failed to compile or link. See the file config.log for the"
              echo "*** exact error that occured. This usually means GSL was incorrectly installed"
              echo "*** or that you have moved GSL since it was installed. In the latter case, you"
              echo "*** may want to edit the gsl-config script: $GSL_CONFIG" ])
              CFLAGS="$ac_save_CFLAGS"
              LIBS="$ac_save_LIBS"
           fi
         fi
    #     GSL_FLAGS=""
    #     GSL_LIBS=""
         ifelse([$3], , :, [$3])
      fi
    
  fi
  AC_SUBST(GSL_FLAGS)
  AC_SUBST(GSL_LIBS)
  rm -f conf.gsltest

])

AU_ALIAS([AM_PATH_GSL], [AX_PATH_GSL])
