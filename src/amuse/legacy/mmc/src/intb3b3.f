*
      subroutine inb3b3(nup)
*
*       calculate interaction between two three-body binaries according to
*       ------------------------------------------------------------------
*       Sdodolkiewicz's prescription (Acta Astronomica, 1986, 36, 19)
*       -------------------------------------------------------------
*       ib3f = 3 - Pmax used to calculate encounter probability
*       -------------------------------------------------------
*       ib3f = 4 - Pmax used to calculate encounter probability and
*       -----------------------------------------------------------
*                  numerical integration of four-body interaction
*                  ----------------------------------------------
*
*
      include 'common.h'
*
      integer nup
*
      if((nbin3-nescb3-ndist3-ndist4-ndiste-nmerge).le.1) return
*
*     if number of binaries presented in the system is less then 100
*     compute probability for binary-binary interaction using
*     Stodolkiewicz procedure (Acta Astronomica 1986), otherwise
*     compute probability for binary-binary interaction using Pmax
*     (Hut and Bahcal ApJ 268, 1983 and Bacon and Sigurdsson 
*     astro-ph96/03036)
*
      if((nbin3-nescb3-ndist3-ndist4-ndiste-nmerge).le.100) then
*
        call intb3b3a(nup)
*
      else
*
        call intb3b3b(nup)
*
      endif
*
      return
*
      end
*
*
*
