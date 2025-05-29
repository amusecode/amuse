      subroutine inext(NQ,LISTQ,TMIN,length,nxtlst)
*
*       Next time-step block.
*       ---------------------
*
      include 'common6.h'
      integer length,j,NQ,nxtlst(*),LISTQ(*)
*
*       Find all new times equal to TMIN.
      length = 0
      do 10 L = 1,NQ
         j = LISTQ(L)
         if (TNEW(j) .eq. TMIN) then
            length = length + 1
            nxtlst(length) = j
         end if
 10   continue
*
      return
      end
