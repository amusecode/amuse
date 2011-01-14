       subroutine start_viewer
#ifndef NOMPI
#include "general.inc"
       if(myRank.EQ.rootRank) then
#endif
        print*,myRank,rootRank
        call viewer
#ifndef NOMPI
       endif
#endif
       end subroutine

