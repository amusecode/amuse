###############
# EXECUTABLES #
###############

#############
# LIBRARIES #
#############

INCLUDE_DIRECTORIES(./ ../)
ADD_LIBRARY(usm usm.cpp)

ADD_EXECUTABLE(ez2usm ez2usm.cpp)
TARGET_LINK_LIBRARIES(ez2usm usm eos)

ADD_EXECUTABLE(usm2col usm2col.cpp)
TARGET_LINK_LIBRARIES(usm2col usm eos)

ADD_EXECUTABLE(usm2quad usm2quad.cpp)
TARGET_LINK_LIBRARIES(usm2quad usm eos)

################
# INSTALLATION #
################

# INSTALL(TARGETS ez2usm usm2col usm2quad
#  RUNTIME DESTINATION bin
# )
