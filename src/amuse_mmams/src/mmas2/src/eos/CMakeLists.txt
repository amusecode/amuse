###############
# EXECUTABLES #
###############

#############
# LIBRARIES #
#############

INCLUDE_DIRECTORIES(./ ../)
ADD_LIBRARY(eos eos.cpp)

SET_SOURCE_FILES_PROPERTIES(eos.cpp
  PROPERTIES COMPILE_FLAGS "-O0 -g")

################
# INSTALLATION #
################

