      SUBROUTINE initgrape
C===  init the local GRAPE
      IMPLICIT NONE
#include "general.inc"

C---  local vars
      INTEGER new_xunit,new_tunit,aflag,jflag,pflag

C---  functions
      INTEGER g6_npipes

C---  set x- and tunit
      new_xunit = 48
      new_tunit = 48
      aflag     = 1
      jflag     = 1
      pflag     = 1
      clusterid = 0

C---  open GRAPE and set some parameters
      CALL  g6_open(clusterid)
C           -------
      npipe = g6_npipes()

      CALL g6_set_tunit(new_tunit)
      CALL g6_set_xunit(new_xunit)


#ifdef DEBUG
      WRITE(STDERR,600) myRank
 600  FORMAT('PE',I2.2,': local GRAPE is initialized')
      CALL flush(STDERR)
#endif

      RETURN
      END
