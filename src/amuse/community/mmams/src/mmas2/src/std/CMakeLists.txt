###############
# EXECUTABLES #
###############

#############
# LIBRARIES #
#############

INCLUDE_DIRECTORIES(./ ../)
ADD_LIBRARY(std pgetopt.cpp trapfpe.c)

################
# INSTALLATION #
################
