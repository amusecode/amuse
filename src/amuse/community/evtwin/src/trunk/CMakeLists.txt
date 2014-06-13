cmake_minimum_required(VERSION 2.4 FATAL_ERROR)

# Set build type. Do this *before* we set the project name
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo Profile."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)
set(CMAKE_CONFIGURATION_TYPES "${CMAKE_BUILD_TYPE}" CACHE INTERNAL "internal")

project(TWIN Fortran)

include(CMakeDependentOption)

# Various options that we may want to enable
option(WANT_ASSERT            "Use assertions (run-time sanity checks; requires preprocessor)" off)
option(WANT_OPENMP            "Use OpenMP parallelisation (experimental)" off)
option(WANT_AUTOPARALLEL      "Use the auto-paralleliser to parallelise the code (experimental)" off)
option(WANT_SSE42             "Enable generation of SSE4.2 code" on)
option(WANT_HOST_OPT          "Enable host-specific optimisation. Choose only when compiling and running on the same machine! Overrides WANT_SSE42" off)
option(WANT_IPO               "Inter-procedural optimisation" off)
option(WANT_O3                "Use O3 optimisations rather than O2. May be unsafe" off)
option(WANT_STATIC            "Generate statically linked executable" off)
option(WANT_CHECKS            "Activate runtime checks (array bounds, NaNs)" off)
option(WANT_WARNINGS          "Issue compile-time warnings" on)
option(WANT_STRICT_WARNINGS   "Issue strict compile-time warnings (for development)" off)
option(WANT_LIBRARY           "Compile the code as a library" on)
option(WANT_STRICT_FLOATS     "Enable strict floating point calculations" off)
option(WANT_PROFILING         "Compile the code with profiling symbols" off)
#cmake_dependent_option(WANT_MKMERGERMOD "Create merger import utility" on "WANT_LIBRARY" off)
option(WANT_MKMERGERMOD "Create merger import utility" off)
cmake_dependent_option(WANT_DYNAMIC_LIBRARY "Create dynamic (.so) library" on "WANT_LIBRARY" off)

# Search in the `cmake' directory for additional CMake modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Main evolution code source files
include(FileList)

# Set FORTRAN compiler flags
include(CompilerFlags)

# Tell CMake the source won't be available until build time.
SET_SOURCE_FILES_PROPERTIES(${CMAKE_SOURCE_DIR}/code/svn_version.f90 PROPERTIES GENERATED 1)
SET_SOURCE_FILES_PROPERTIES(${CMAKE_SOURCE_DIR}/code/installation_path.f90 PROPERTIES GENERATED 1)

# SVN version generator
ADD_CUSTOM_COMMAND(
   OUTPUT ${CMAKE_SOURCE_DIR}/code/svn_version.f90
   COMMAND ${CMAKE_SOURCE_DIR}/detect_svn_version.pl ${CMAKE_SOURCE_DIR}
   DEPENDS ${TWIN_SRC_FILES}
)
ADD_CUSTOM_COMMAND(
   OUTPUT ${CMAKE_SOURCE_DIR}/code/installation_path.f90
   COMMAND ${CMAKE_SOURCE_DIR}/set_standard_path.pl ${CMAKE_SOURCE_DIR} ${CMAKE_INSTALL_PREFIX}
)

# Code library
add_library ("twincore" STATIC ${TWIN_SRC_FILES})
set_target_properties(twincore PROPERTIES COMPILE_FLAGS "-fPIC" )

# We need to changed the flags passed to ranlib on OS X, or we get
# undefined references for all symbols.
if(APPLE)
   set(CMAKE_Fortran_ARCHIVE_FINISH "<CMAKE_RANLIB> -c <TARGET>")
endif(APPLE)

if (WANT_LIBRARY)
   if (WANT_DYNAMIC_LIBRARY)
      add_library ("dynamic_twin" SHARED ${TWIN_SRC_FILES} ${TWIN_LIBRARY_FILES})
      set_target_properties(dynamic_twin PROPERTIES OUTPUT_NAME "twin")
   endif(WANT_DYNAMIC_LIBRARY)
   add_library ("twin" STATIC ${TWIN_LIBRARY_FILES})
   target_link_libraries("twin" twincore)

   add_executable ("twintest" code/twin.f90)
   target_link_libraries("twintest" twin)
   set_target_properties(twintest PROPERTIES OUTPUT_NAME "twin")
endif (WANT_LIBRARY)

# Main TWIN executable
#add_executable ("ev" ${TWIN_SRC_FILES} code/svn_version.f90 code/ev_main.f90)
add_executable ("ev" code/svn_version.f90 code/ev_main.f90)
target_link_libraries("ev" twincore)

# Merger import utility
if(WANT_MKMERGERMOD)
   add_executable ("mkmergermod" code/mkmergermod.f90)
   target_link_libraries("mkmergermod" twin)
endif(WANT_MKMERGERMOD)

# Installation targets
install (TARGETS "ev" RUNTIME DESTINATION "bin")
install (FILES "${CMAKE_SOURCE_DIR}/input/nrates.dat"                      DESTINATION "stars/input/")
install (FILES "${CMAKE_SOURCE_DIR}/input/nucdata.dat"                     DESTINATION "stars/input/")
install (FILES "${CMAKE_SOURCE_DIR}/input/rates.dat"                       DESTINATION "stars/input/")
install (FILES "${CMAKE_SOURCE_DIR}/input/mutate.dat"                      DESTINATION "stars/input/")
install (FILES "${CMAKE_SOURCE_DIR}/input/physinfo.dat"                    DESTINATION "stars/input/")
install (FILES "${CMAKE_SOURCE_DIR}/input/zahb.dat"                        DESTINATION "stars/input/")
install (FILES "${CMAKE_SOURCE_DIR}/input/zahb02.mod"                      DESTINATION "stars/input/")
install (FILES "${CMAKE_SOURCE_DIR}/input/nrates.dat"                      DESTINATION "stars/input/")
install (FILES "${CMAKE_SOURCE_DIR}/input/lt2ubv.dat"                      DESTINATION "stars/input/")
install (FILES "${CMAKE_SOURCE_DIR}/input/poly32.dat"                      DESTINATION "stars/input/")

install (FILES "${CMAKE_SOURCE_DIR}/input/COtables/COtables_z00001"        DESTINATION "stars/input/COtables/")
install (FILES "${CMAKE_SOURCE_DIR}/input/COtables/COtables_z00003"        DESTINATION "stars/input/COtables/")
install (FILES "${CMAKE_SOURCE_DIR}/input/COtables/COtables_z0001"         DESTINATION "stars/input/COtables/")
install (FILES "${CMAKE_SOURCE_DIR}/input/COtables/COtables_z0003"         DESTINATION "stars/input/COtables/")
install (FILES "${CMAKE_SOURCE_DIR}/input/COtables/COtables_z001"          DESTINATION "stars/input/COtables/")
install (FILES "${CMAKE_SOURCE_DIR}/input/COtables/COtables_z002"          DESTINATION "stars/input/COtables/")
install (FILES "${CMAKE_SOURCE_DIR}/input/COtables/COtables_z003"          DESTINATION "stars/input/COtables/")
install (FILES "${CMAKE_SOURCE_DIR}/input/COtables/COtables_z004"          DESTINATION "stars/input/COtables/")
install (FILES "${CMAKE_SOURCE_DIR}/input/COtables/COtables_z006"          DESTINATION "stars/input/COtables/")
install (FILES "${CMAKE_SOURCE_DIR}/input/COtables/COtables_z008"          DESTINATION "stars/input/COtables/")
install (FILES "${CMAKE_SOURCE_DIR}/input/COtables/COtables_z01"           DESTINATION "stars/input/COtables/")
install (FILES "${CMAKE_SOURCE_DIR}/input/COtables/COtables_z02"           DESTINATION "stars/input/COtables/")
install (FILES "${CMAKE_SOURCE_DIR}/input/COtables/COtables_z03"           DESTINATION "stars/input/COtables/")
install (FILES "${CMAKE_SOURCE_DIR}/input/COtables/COtables_z04"           DESTINATION "stars/input/COtables/")
install (FILES "${CMAKE_SOURCE_DIR}/input/COtables/COtables_z05"           DESTINATION "stars/input/COtables/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z0/phys.z0"               DESTINATION "stars/input/metals/z0/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z0/phys_a04.z0"           DESTINATION "stars/input/metals/z0/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z0001/phys.z0001"         DESTINATION "stars/input/metals/z0001/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z0001/phys_a04.z0001"     DESTINATION "stars/input/metals/z0001/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z0001/zams0001_start.mod" DESTINATION "stars/input/metals/z0001/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z0003/phys.z0003"         DESTINATION "stars/input/metals/z0003/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z0003/phys_a04.z0003"     DESTINATION "stars/input/metals/z0003/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z0003/zams0003_start.mod" DESTINATION "stars/input/metals/z0003/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z001/phys.z001"           DESTINATION "stars/input/metals/z001/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z001/phys_a04.z001"       DESTINATION "stars/input/metals/z001/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z001/zams001_start.mod"   DESTINATION "stars/input/metals/z001/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z002/phys.z002"           DESTINATION "stars/input/metals/z002/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z002/phys_a04.z002"       DESTINATION "stars/input/metals/z002/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z004/phys.z004"           DESTINATION "stars/input/metals/z004/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z004/phys_a04.z004"       DESTINATION "stars/input/metals/z004/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z004/zams004_start.mod"   DESTINATION "stars/input/metals/z004/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z01/phys.z01"             DESTINATION "stars/input/metals/z01/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z01/phys_a04.z01"         DESTINATION "stars/input/metals/z01/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z01/zams01_start.mod"     DESTINATION "stars/input/metals/z01/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z02/phys.z02"             DESTINATION "stars/input/metals/z02/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z02/phys_a04.z02"         DESTINATION "stars/input/metals/z02/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z02/zams02_start.mod"     DESTINATION "stars/input/metals/z02/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z03/phys.z03"             DESTINATION "stars/input/metals/z03/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z03/phys_a04.z03"         DESTINATION "stars/input/metals/z03/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z03/zams03_start.mod"     DESTINATION "stars/input/metals/z03/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z04/phys.z04"             DESTINATION "stars/input/metals/z04/")
install (FILES "${CMAKE_SOURCE_DIR}/input/metals/z04/phys_a04.z04"         DESTINATION "stars/input/metals/z04/")
install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams02.mas"                 DESTINATION "stars/input/zams/")
install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams02.mod"                 DESTINATION "stars/input/zams/")
install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams02.out"                 DESTINATION "stars/input/zams/")

install (FILES "${CMAKE_SOURCE_DIR}/stars_standards/init.dat"              DESTINATION "stars/stars_standards/")
install (FILES "${CMAKE_SOURCE_DIR}/stars_standards/init.run"              DESTINATION "stars/stars_standards/")
install (FILES "${CMAKE_SOURCE_DIR}/stars_standards/init_twin.dat"         DESTINATION "stars/stars_standards/")
install (FILES "${CMAKE_SOURCE_DIR}/stars_standards/init_binary.dat"       DESTINATION "stars/stars_standards/")

#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams0001.mas" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams0001.mod" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams0001.out" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams0003.mas" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams0003.mod" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams0003.out" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams001.mas" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams001.mod" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams001.out" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams004.mas" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams004.mod" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams004.out" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams01.mas" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams01.mod" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams01.out" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams016.mas" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams016.mod" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams016.out" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams03.mas" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams03.mod" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zams/zams03.out" DESTINATION "stars/input/")

#install (FILES "${CMAKE_SOURCE_DIR}/input/zahb001.mod" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zahb004.mod" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zahb0001.mod" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zahb002.mod" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zahb006.mod" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zahb0003.mod" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zahb003.mod" DESTINATION "stars/input/")
#install (FILES "${CMAKE_SOURCE_DIR}/input/zahb008.mod" DESTINATION "stars/input/")

# Place the executable in code/ in the top-level STARS directory, for
# historical reasons.
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/code)
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)
