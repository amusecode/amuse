CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

PROJECT(MMaS2)

IF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.3)

#------------------------------
# MMaS2 vesion number.
SET(MMaS2_MAJOR_VERSION 2)
SET(MMaS2_MINOR_VERSION 0)
SET(MMaS2_VERSION "${MMaS2_MAJOR_VERSION}.${MMaS2_MINOR_VERSION}")

# Platform configuration tests.
INCLUDE(${CMAKE_ROOT}/Modules/CMakeBackwardCompatibilityC.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/TestForANSIStreamHeaders.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/TestForSTDNamespace.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/TestForANSIForScope.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckTypeSize.cmake)

# Simulate old CMakeBackwardCompatibilityCXX test.
# In CMake 2.4 and up this could be just
# INCLUDE(${CMAKE_ROOT}/Modules/TestForSSTREAM.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckIncludeFileCXX.cmake)
CHECK_INCLUDE_FILE_CXX("sstream" CMAKE_HAS_ANSI_STRING_STREAM)
IF(NOT CMAKE_HAS_ANSI_STRING_STREAM)
  SET(CMAKE_NO_ANSI_STRING_STREAM 1 CACHE INTERNAL
    "Does the compiler support sstream or stringstream.")
ENDIF(NOT CMAKE_HAS_ANSI_STRING_STREAM)

SET(CMAKE_C_FLAGS "-Wall -O3")
#SET(CMAKE_CXX_FLAGS "-Wall -O0 -g")
SET(CMAKE_INSTALL_PREFIX ./usr)

SUBDIRS(src)

ELSE("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.3)
  MESSAGE(SEND_ERROR
          "This version of CMake is too old to build MMaS2.  "
          "Please upgrade to CMake 2.3.")


ENDIF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.3)

