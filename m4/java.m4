### java.m4 -- macros for Java environment detection    -*- Autoconf -*-
###
### Copyright (C) 2005-7 R Core Team
###
### This file is part of R.
###
### R is free software; you can redistribute it and/or modify it under
### the terms of the GNU General Public License as published by the Free
### Software Foundation; either version 2 of the License, or (at your
### option) any later version.
###
### R is distributed in the hope that it will be useful, but WITHOUT ANY
### WARRANTY; without even the implied warranty of MERCHANTABILITY or
### FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
### License for more details.
###
### You should have received a copy of the GNU General Public License
### along with R; if not, a copy is available at
### http://www.r-project.org/Licenses/

### This java.m4 file adjusted for use in the AMUSE project by Niels Drost

## R_RUN_JAVA(variable for the result, parameters)
## ----------
## runs the java interpreter ${JAVA} with specified parameters and
## saves the output to the supplied variable. The exit value is ignored.
AC_DEFUN([R_RUN_JAVA],
[
  acx_java_result=
  if test -z "${JAVA}"; then
    echo "$as_me:$LINENO: JAVA is not set, cannot run java $2" >&AS_MESSAGE_LOG_FD
  else
    echo "$as_me:$LINENO: running ${JAVA} $2" >&AS_MESSAGE_LOG_FD
    acx_java_result=`${JAVA} $2 2>&AS_MESSAGE_LOG_FD`
    echo "$as_me:$LINENO: output: '$acx_java_result'" >&AS_MESSAGE_LOG_FD
  fi
  $1=$acx_java_result
])


## R_JAVA
## -----------
## Looks for Java JRE/JDK and sets:
## have_java to yes/no; if it is yes then also sets:
## JAVA to Java interpreter path
## JAVA_HOME to the home directory of the Java runtime/jdk
## JAVA_LD_LIBRARY_PATH to the path necessary for Java runtime
## JAVAC to Java compiler path (optional)
## JAR to Java archiver (optional)
##
## (*) - those variables are modified for use in make files
##       to rely on $(JAVA_HOME) and substituted with 0 suffix
##
## JAVA_HOME env var is honored during the search and the search
## will fail if it is set incorrectly.
AC_DEFUN([R_JAVA],
[
have_java=no

## find java compiler binaries
if test -z "${JAVA_HOME}" ; then
  JAVA_PATH=${PATH}
else
  ## try jre/bin first just in case we don't have full JDK
  JAVA_PATH=${JAVA_HOME}:${JAVA_HOME}/jre/bin:${JAVA_HOME}/bin:${JAVA_HOME}/../bin:${PATH}
fi
## if 'java' is not on the PATH or JAVA_HOME, add some guesses as of
## where java could live
JAVA_PATH=${JAVA_PATH}:/usr/java/bin:/usr/jdk/bin:/usr/lib/java/bin:/usr/lib/jdk/bin:/usr/local/java/bin:/usr/local/jdk/bin:/usr/local/lib/java/bin:/usr/local/lib/jdk/bin
AC_PATH_PROGS(JAVA,java,,${JAVA_PATH})
## FIXME: we may want to check for jikes, kaffe and others...
## (however, most of them have compatibility wrappers by now)
AC_PATH_PROGS(JAVAC,javac,,${JAVA_PATH})
AC_PATH_PROGS(JAR,jar,,${JAVA_PATH})

## we don't require a compiler, but it would be useful
AC_CACHE_CHECK([whether Java compiler works], [r_cv_javac_works],
[r_cv_javac_works=no
if test -n "${JAVAC}"; then
  rm -f A.java A.class
  echo "public class A { }" > A.java
  if "${JAVAC}" A.java 2>&AS_MESSAGE_LOG_FD; then
    if test -f A.class; then
      r_cv_javac_works=yes
    fi
  fi
  rm -rf A.java A.class
fi])

## this is where our test-class lives (in the AMUSE support dir)
getsp_cp=support

AC_CACHE_CHECK([whether Java interpreter works], [r_cv_java_works],
[r_cv_java_works=no
if test -n "${JAVA}" ; then
  R_RUN_JAVA(acx_jc_result,[-classpath ${getsp_cp} getsp -test])
  if test "${acx_jc_result}" = "Test1234OK"; then
    r_cv_java_works=yes
  fi
  acx_jc_result=
fi])


if test ${r_cv_java_works} = yes; then
  AC_CACHE_CHECK([Java environment], [r_cv_java_home], [
    ## find JAVA_HOME from Java itself unless specified
    if test -z "${JAVA_HOME}" ; then
      R_RUN_JAVA(JAVA_HOME,[-classpath ${getsp_cp} getsp java.home])
    fi
    r_cv_java_home="${JAVA_HOME}"
  ])
  JAVA_HOME="${r_cv_java_home}"

  AC_CACHE_CHECK([Java version], [r_cv_java_version], [
  ## Detect the java version
  R_RUN_JAVA(r_cv_java_version,[-classpath ${getsp_cp} getsp java.specification.version])

  ])
  JAVA_VERSION="${r_cv_java_version}"

  have_java="yes"
else  ## not r_cv_java_works
    have_java="no"
    AC_MSG_WARN([Java not found, Java codes disabled])
fi

## AC_SUBST(JAVA_HOME) # not needed? is precious now
AC_SUBST(have_java)
AC_SUBST(JAVA)
AC_SUBST(JAVAC)
AC_SUBST(JAR)
AC_SUBST(JAVA_VERSION)
])
