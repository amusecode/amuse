###############
# EXECUTABLES #
###############

#############
# LIBRARIES #
#############


#######
# GSL #
#######

IF (NOT GSL_LIBS)
  FIND_LIBRARY(GSL_LIB gsl PATHS /usr/lib /usr/local/lib /opt/local/lib /sw/lib)
  IF (NOT GSL_LIB)
    MESSAGE(FATAL_ERROR "Cannot find libgsl. Perhaps GSL developer packages are not installed.")
  ENDIF (NOT GSL_LIB)
  FIND_LIBRARY(GSLCBLAS_LIB gslcblas PATHS /usr/lib /usr/local/lib /opt/local/lib /sw/lib)
  SET(GSL_LIBS  ${GSL_LIB} ${GSLCBLAS_LIB})
ELSE (NOT GSL_LIBS)
ENDIF (NOT GSL_LIBS)

IF (NOT GSL_FLAGS)
FIND_PATH(GSL_INC gsl/gsl_sys.h PATHS /usr/include /usr/local/include /opt/local/include /sw/include)
SET(INC ${INC} ${GSL_INC})
ELSE (NOT GSL_FLAGS)
set_source_files_properties(
        merge_stars_gsl.cpp
        merge_stars_gsl_adaptive.cpp
        merge_stars_consistently.cpp
        shock_heating.cpp
        compute_energy_gsl.cpp
        mixing.cpp
	PROPERTIES COMPILE_FLAGS "${GSL_FLAGS}")
ENDIF (NOT GSL_FLAGS)

SET(LIBS ${LIBS} ${GSL_LIBS})


INCLUDE_DIRECTORIES(./ ../ ${INC})
ADD_LIBRARY(_mmas2
#  merge_stars_gsl.cpp
  merge_stars_gsl_adaptive.cpp
  merge_stars_consistently.cpp
  compute_extra.cpp
  sort_entropy.cpp
  shock_heating.cpp
  compute_energy_gsl.cpp
  mixing.cpp
  mass_loss.cpp
  )
TARGET_LINK_LIBRARIES(_mmas2 ${LIBS})

ADD_EXECUTABLE(mmas2 mmas2.cpp)
TARGET_LINK_LIBRARIES(mmas2 _mmas2 eos usm std)

# SET_SOURCE_FILES_PROPERTIES(merge_stars_gsl.cpp mmas2.cpp mixing.cpp
#  PROPERTIES COMPILE_FLAGS "-O0 -gstabs+")

################
# INSTALLATION #
################

INSTALL(TARGETS
  mmas2
  RUNTIME DESTINATION bin
)
