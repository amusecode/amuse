##  FLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)


##  The Fortran_FLAGS_* are used in addition to, after and thus overruling the Fortran_FLAGS
##  Overruling is especially likely for the different -O flags
##  To see exactly which flags are passed to the compiler, select CMAKE_VERBOSE_MAKEFILE, e.g. with cmake-gui

if (Fortran_COMPILER_NAME STREQUAL "gfortran")
   EXEC_PROGRAM(${CMAKE_Fortran_COMPILER}
    ARGS --version
    OUTPUT_VARIABLE _evtwin_COMPILER_VERSION
   )
   STRING(REGEX REPLACE ".* ([0-9])\\.([0-9])\\.[0-9].*" "\\1\\2" _evtwin_COMPILER_VERSION ${_evtwin_COMPILER_VERSION})
   if(_evtwin_COMPILER_VERSION LESS 43 )
      set (CMAKE_Fortran_FLAGS "-W -O2")
   else (_evtwin_COMPILER_VERSION LESS 43 )
      set (CMAKE_Fortran_FLAGS "-O2 -finit-local-zero")
   endif(_evtwin_COMPILER_VERSION LESS 43 )
   set (CMAKE_Fortran_FLAGS_RELEASE "-pipe -funroll-all-loops")
   set (CMAKE_Fortran_FLAGS_DEBUG "-g -ffpe-trap=zero,invalid -fsignaling-nans")
   set (CMAKE_Fortran_FLAGS_PROFILE "-g -gp")

   if(NOT _evtwin_COMPILER_VERSION LESS 43 )
    if(WANT_SSE42)
         set (SSE_FLAGS "-msse4.2")
    endif(WANT_SSE42)
   endif(NOT _evtwin_COMPILER_VERSION LESS 43 )
   
   if (WANT_OPENMP)
      set (OPENMP_FLAGS "-fopenmp")
   endif (WANT_OPENMP)

   if (WANT_STATIC)
      set (STATIC_FLAGS "-static")
   endif (WANT_STATIC)

   if (WANT_CHECKS)
      set (CHECK_FLAGS "-O0 -fbounds-check -ffpe-trap=zero,invalid -fsignaling-nans -fbacktrace") # v.4.4
      #set (CHECK_FLAGS "-O0 -fcheck=all -ffpe-trap=zero,invalid -fsignaling-nans -fbacktrace")  # From v.4.5
   endif (WANT_CHECKS)

   if (WANT_WARNINGS)
      set (WARN_FLAGS "-Wall")
   endif (WANT_WARNINGS)

   if (WANT_PROFILING)
      set (PROF_FLAGS "-O0 -pg")
   endif (WANT_PROFILING)

   if (WANT_LIBRARY)
      set (LIB_FLAGS "-fPIC")
   endif (WANT_LIBRARY)

elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
   
   set (CMAKE_Fortran_FLAGS "-O2")
   set (CMAKE_Fortran_FLAGS_RELEASE "-vec-guard-write -fpconstant -extend_source -funroll-loops -align all -ip")
   set (CMAKE_Fortran_FLAGS_DEBUG "-g -debug -traceback")
   set (CMAKE_Fortran_FLAGS_PROFILE "-g -gp")

   # HOST_OPT overrules SSE42, as ifort would
   if(WANT_SSE42)
     set (SSE_FLAGS "-axSSE4.2,SSSE3")
   endif(WANT_SSE42)
   if(WANT_HOST_OPT)
     set (SSE_FLAGS "-xHost")
   endif(WANT_HOST_OPT)

   if(WANT_IPO)
      set (IPO_FLAGS "-ipo")
   endif(WANT_IPO)

   if (WANT_STRICT_FLOATS)
      set (FP_FLAGS "-fp-model strict")
   endif (WANT_STRICT_FLOATS)

   if (WANT_OPENMP)
      set (OPENMP_FLAGS "-openmp -openmp-report2")
   endif (WANT_OPENMP)

   if (WANT_STATIC)
      set (STATIC_FLAGS "-static")
   endif (WANT_STATIC)

   if (WANT_CHECKS)
      set (CHECK_FLAGS "-O0 -ftrapuv -stand f03 -check all -check noarg_temp_created -traceback")
   endif (WANT_CHECKS)

   if (WANT_WARNINGS)
      set (WARN_FLAGS "-warn -stand f03")
   endif (WANT_WARNINGS)

   if (WANT_PROFILING)
      set (PROF_FLAGS "-O0 -p -g")
   endif (WANT_PROFILING)

   if (WANT_LIBRARY)
      set (LIB_FLAGS "-fPIC")
    endif (WANT_LIBRARY)

else (Fortran_COMPILER_NAME STREQUAL "gfortran")

   message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
   message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
   message ("No optimized Fortran compiler flags are known, we just try -O2...")
   set (CMAKE_Fortran_FLAGS "-O2")
   set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
   set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")

endif (Fortran_COMPILER_NAME STREQUAL "gfortran")


##  For all compilers, combine the flags.
##  The compiler uses (e.g. for release) CMAKE_Fortran_FLAGS CMAKE_Fortran_FLAGS_RELEASE, 
##  so Fortran_FLAGS and USER_FLAGS are overruled by Fortran_FLAGS_*

set (USER_FLAGS "${LIB_FLAGS}  ${WARN_FLAGS}  ${SSE_FLAGS}  ${IPO_FLAGS}  ${FP_FLAGS}  ${OPENMP_FLAGS}  ${STATIC_FLAGS}  ${CHECK_FLAGS} ${PROF_FLAGS} ")

set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  ${USER_FLAGS}")
#set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}  ${USER_FLAGS}")  #Fortran_FLAGS are used in addtion to Fortran_FLAGS_RELEASE

set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE}  -g")

