cmake_minimum_required(VERSION 2.4 FATAL_ERROR)

# Set build type. Do this *before* we set the project name
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo Profile."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)
set(CMAKE_CONFIGURATION_TYPES "${CMAKE_BUILD_TYPE}" CACHE INTERNAL "internal")

project(TWIN Fortran)

# Various options that we may want to enable
option(WANT_OPENMP "Use OpenMP parallelisation (experimental)" off)
option(WANT_SSE42  "Enable generation of SSE4.2 code" on)
option(WANT_HOST_OPT  "Enable host-specific optimisation. Choose only when compiling and running on the same machine! Overrides WANT_SSE42" off)
option(WANT_IPO    "Inter-procedural optimisation" off)
option(WANT_STATIC "Generate statically linked executable" off)
option(WANT_CHECKS "Activate runtime checks (array bounds, NaNs)" off)
option(WANT_WARNINGS "Activate warnings" off)
option(WANT_LIBRARY "Compile the code as a library" on)
option(WANT_STRICT_FLOATS "Enable strict floating point calculations" off)
option(WANT_PROFILING "Compile the code with profiling symbols" off)

# Search in the `cmake' directory for additional CMake modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Main evolution code source files
include(FileList)

# Set FORTRAN compiler flags
include(CompilerFlags)

# Tell CMake the source won't be available until build time.
SET_SOURCE_FILES_PROPERTIES(${CMAKE_SOURCE_DIR}/code/svn_version.f90 PROPERTIES GENERATED 1)

# SVN version generator
ADD_CUSTOM_COMMAND(
   OUTPUT ${CMAKE_SOURCE_DIR}/code/svn_version.f90
   COMMAND ${CMAKE_SOURCE_DIR}/detect_svn_version.pl ${CMAKE_SOURCE_DIR}
   DEPENDS ${TWIN_SRC_FILES}
)


# Main TWIN executable
add_executable ("ev" ${TWIN_SRC_FILES} code/svn_version.f90 code/ev_main.f90)

# MUSE library
if(WANT_LIBRARY)
   add_library ("twin" STATIC ${TWIN_SRC_FILES} ${TWIN_LIBRARY_FILES})
endif(WANT_LIBRARY)

# Installation targets
install (TARGETS "ev" RUNTIME DESTINATION "bin")

# Place the executable in code/ in the top-level STARS directory, for
# historical reasons.
set(EXECUTABLE_OUTPUT_PATH ${BUILDDIR})
set(LIBRARY_OUTPUT_PATH ${BUILDDIR})
