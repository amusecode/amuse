##  FLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)


##  The Fortran_FLAGS_* are used in addition to, after and thus overruling the Fortran_FLAGS
##  Overruling is especially likely for the different -O flags
##  To see exactly which flags are passed to the compiler, select CMAKE_VERBOSE_MAKEFILE, e.g. with cmake-gui
set (OFLAGS "-O2")
if (WANT_O3)
   set (OFLAGS "-O3")
endif (WANT_O3)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   EXEC_PROGRAM(${CMAKE_Fortran_COMPILER}
        ARGS --version
        OUTPUT_VARIABLE _evtwin_COMPILER_VERSION
   )
   set (FPP_FLAGS "-cpp")
    STRING(REGEX REPLACE ".* ([0-9])\\.([0-9])\\.[0-9].*" "\\1\\2" _evtwin_COMPILER_VERSION ${_evtwin_COMPILER_VERSION})
   message ("-- Fortran compiler version: " ${_evtwin_COMPILER_VERSION})
   EXEC_PROGRAM(${CMAKE_Fortran_COMPILER}
        ARGS -dumpmachine
        OUTPUT_VARIABLE _evtwin_MACHINE
   )
   string(REPLACE "-" ";" _evtwin_MACHINE_LIST ${_evtwin_MACHINE})
   message ("-- Fortran arch: " ${_evtwin_MACHINE_LIST})
   list(GET _evtwin_MACHINE_LIST 0 _evtwin_MACHINE_ARCH) 
   message ("-- Fortran arch: " ${_evtwin_MACHINE_ARCH})
   
   if(NOT _evtwin_COMPILER_VERSION LESS 43 )
    set (CMAKE_Fortran_FLAGS "-finit-local-zero")
   endif(NOT _evtwin_COMPILER_VERSION LESS 43 )
   
   set (CMAKE_Fortran_FLAGS_RELEASE "${OFLAGS} -pipe -funroll-all-loops")
   set (CMAKE_Fortran_FLAGS_DEBUG "-g -ffpe-trap=zero,invalid -fsignaling-nans")
   set (CMAKE_Fortran_FLAGS_PROFILE "-g -gp")

   
   if(NOT _evtwin_COMPILER_VERSION LESS 43 )
    if(_evtwin_MACHINE_ARCH STREQUAL "x86_64")
     if(WANT_SSE42)
         set (SSE_FLAGS "-msse4.2")
     endif(WANT_SSE42)
    endif(_evtwin_MACHINE_ARCH STREQUAL "x86_64")
    
   endif(NOT _evtwin_COMPILER_VERSION LESS 43 )
   
   set (OPENMP_FLAGS "-fopenmp")
   set (STATIC_FLAGS "-static")
   set (FFLAGS "-std=f2008")
   set (FFLAGS "-std=gnu") # Needed for isnan work-around!

   set (CHECK_FLAGS "-O0 -ffpe-trap=zero,invalid,overflow -fsignaling-nans -g -fbacktrace")
   if( COMPILER_VERSION VERSION_GREATER "4.4.99" )
      set( CHECK_FLAGS "-fcheck=all ${CHECK_FLAGS}" )    # >= v.4.5
   else( COMPILER_VERSION VERSION_GREATER "4.4.99" )
      set( CHECK_FLAGS "-fbounds-check ${CHECK_FLAGS}" ) # <= v.4.4
   endif( COMPILER_VERSION VERSION_GREATER "4.4.99" )

   set( FPP_FLAGS "-Dieee_arithmetic=ieee_arithmetic_stub ${FPP_FLAGS}" )

   set (WARN_FLAGS "-Wall")

   if (WANT_STRICT_WARNINGS)
      set (WARN_FLAGS "${WARN_FLAGS} -Wextra -fall-intrinsics -pedantic -Werror")
   endif (WANT_STRICT_WARNINGS)

   set (PROF_FLAGS "-O0 -pg")

   set (LIB_FLAGS "-fPIC")

elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")

   
   set (FPP_FLAGS "-fpp")
   set (CMAKE_Fortran_FLAGS "-nogen-interfaces")
   set (CMAKE_Fortran_FLAGS_RELEASE "${OFLAGS} -vec-guard-write -fpconstant -extend_source -funroll-loops -align all -ip")
   set (CMAKE_Fortran_FLAGS_DEBUG "-g -debug -traceback")
   set (CMAKE_Fortran_FLAGS_PROFILE "-g -gp")
   set (FFLAGS "-stand f08")

   # HOST_OPT overrules SSE42, as ifort would
   set (SSE_FLAGS "-axSSE4.2,SSSE3")
   if(WANT_HOST_OPT)
     set (SSE_FLAGS "-xHost")
   endif(WANT_HOST_OPT)

   if(WANT_IPO)
      set (IPO_FLAGS "-ipo")
   endif(WANT_IPO)

   if (WANT_STRICT_FLOATS)
      set (FP_FLAGS "-fp-model strict")
   endif (WANT_STRICT_FLOATS)

   set (OPENMP_FLAGS "-openmp -openmp-report2")
   set (STATIC_FLAGS "-static")

   if (WANT_AUTOPARALLEL)
      set (OPENMP_FLAGS "-parallel")
      set (WANT_OPENMP true)
   endif (WANT_AUTOPARALLEL)

   set (CHECK_FLAGS "-O0 -ftrapuv -fpe0 -check all -check noarg_temp_created -traceback")

   if (WANT_WARNINGS)
      set (WARN_FLAGS "-warn -stand f08")
      message( STATUS "Compiling with warnings" )
   endif (WANT_WARNINGS)

   set (PROF_FLAGS "-O0 -p -g")

  set (LIB_FLAGS "-fPIC")

elseif( CMAKE_Fortran_COMPILER_ID MATCHES "G95" )
  
  set( FPP_FLAGS "-cpp")
  set( CMAKE_Fortran_FLAGS "" )
  set( CMAKE_Fortran_FLAGS_RELEASE ${OFLAGS})
  set( CMAKE_Fortran_FLAGS_DEBUG "-O0 -g" )
  set( CHECK_FLAGS "-O0 -fbounds-check -ftrace=full" )
  set( FFLAGS "-std=f2008" )
  
  if( WANT_WARNINGS )
    set( WARN_FLAGS "-Wall -Wno=102,112,136,140,161,163,165" )
    # Ignore the following warnings:
    #   102: unused module procedures
    #   112: variables that are set, but never used
    #   136: unused module variables
    #   140: type cast may lead to precision loss
    #   161: unused module types
    #   163: arguments without intent
    #   165: procedure with implicit interface called
    
    if (WANT_STRICT_WARNINGS)
      set (WARN_FLAGS "${WARN_FLAGS} -Wextra -Werror")
      message( STATUS "Compiling with strict warnings" )
    else (WANT_STRICT_WARNINGS)
      message( STATUS "Compiling with warnings" )
    endif (WANT_STRICT_WARNINGS)
  endif( WANT_WARNINGS )
  
   set( LIB_FLAGS "-fPIC -g" )


elseif (Fortran_COMPILER_NAME MATCHES "xlf.*")
   ##set (CMAKE_Fortran_FLAGS "-brename:flush,flush_")
   set (CMAKE_Fortran_FLAGS "-qextname")

else (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")


   message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
   message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
   message ("Fortran compilervendor : " ${CMAKE_Fortran_COMPILER_ID})
   message ("No optimized Fortran compiler flags are known, we just try ${OFLAGS}...")
   set (CMAKE_Fortran_FLAGS "-O2")
   set (CMAKE_Fortran_FLAGS_RELEASE ${OFLAGS})
   set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")

endif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")


if(WANT_ASSERT)
   set(FPP_FLAGS "${FPP_FLAGS} -DDEBUG")
endif(WANT_ASSERT)

##  For all compilers, combine the flags.
##  The compiler uses (e.g. for release) CMAKE_Fortran_FLAGS CMAKE_Fortran_FLAGS_RELEASE, 
##  so Fortran_FLAGS and USER_FLAGS are overruled by Fortran_FLAGS_*

set (USER_FLAGS "${FPP_FLAGS} ${FFLAGS} ${IPO_FLAGS} ${FP_FLAGS} ")

if(WANT_SSE42)
   set (USER_FLAGS "${USER_FLAGS} ${SSE_FLAGS} ")
endif(WANT_SSE42)

if (WANT_OPENMP)
   set (USER_FLAGS "${USER_FLAGS} ${OPENMP_FLAGS} ")
endif (WANT_OPENMP)

if (WANT_STATIC)
   set (USER_FLAGS "${USER_FLAGS} ${STATIC_FLAGS} ")
endif (WANT_STATIC)

if (WANT_CHECKS)
   set (USER_FLAGS "${USER_FLAGS} ${CHECK_FLAGS} ")
   message( STATUS "Compiling with run-time checks" )
endif (WANT_CHECKS)

if (WANT_WARNINGS)
   set (USER_FLAGS "${USER_FLAGS} ${WARN_FLAGS} ")
   message( STATUS "Compiling with warnings" )
endif (WANT_WARNINGS)

if( WANT_LIBRARY )
   set (USER_FLAGS "${USER_FLAGS} ${LIB_FLAGS} ")
   message( STATUS "Compiling with library options" )
endif( WANT_LIBRARY )

if (WANT_PROFILING)
   set (USER_FLAGS "${USER_FLAGS} ${PROF_FLAGS} ")
endif (WANT_PROFILING)

set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  ${USER_FLAGS}")
#set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}  ${USER_FLAGS}")  #Fortran_FLAGS are used in addtion to Fortran_FLAGS_RELEASE

set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE}  -g")


message( STATUS "Using Fortran compiler:  ${Fortran_COMPILER_NAME}  (${CMAKE_Fortran_COMPILER})" )
message( STATUS "Compiler vendor: ${CMAKE_Fortran_COMPILER_ID}")
message( STATUS "Compiler flags used:  ${CMAKE_Fortran_FLAGS}" )

