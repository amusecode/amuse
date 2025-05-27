      subroutine resort(length, nxtlst, tm)
c
c     ------------------------------------------
c     Maintain particle list sorted by next time
c     ------------------------------------------
c     length (o) length of array with minimum timenw
c     nxtlist (o) index array
c     tm (o) minimum time
      include 'common6.h'
      integer length,j,nxtlst(*)
c
*     Include new times equal to TSMALL (R.Sp. scheme modified 10/98).
      inxt = 0
      do 10 j = IFIRST, NTOT
         if (timenw(j) .eq. TSMALL) then
            inxt = inxt + 1
            nxtlst(inxt) = j
         end if
 10   continue
      length = inxt
*     Copy smallest time from common6.
      tm = TSMALL
      return
      end
