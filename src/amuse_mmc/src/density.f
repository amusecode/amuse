      subroutine densi(k,den)
*
*
*       determine density in the vicinity of the 'k' star
*       --------------------------------------------------
*
*
      include 'common.h'
*
      real*8 den,co
*
      integer k,nspan,n,i1,i2
*
*
      n = nt
c      nspan = nminzo/2
c      nspan = 7
      nspan = 4
 10   co = 3.0d0*float(nspan)

      if(k.gt.nspan) then
        i1 = k - nspan
        i2 = k + nspan
*
        if(i2.gt.n) then
          i2 = n
 20       i1 = n - 2*nspan -1 
          if(r(i1).gt.1.d8) then
            nspan = nspan + 1
            co = 3.0d0*float(nspan)
            print*,'dens  nspan,r1,r2,i1,i2=',nspan,r(i1),r(i2),i1,i2
            go to 20
          else  
            den = co/(r(i2)**3 - r(i1)**3)
            return
          endif
        endif  
*           
        if(r(i1).gt.1.d8.or.r(i2).gt.1.d8) then
          nspan = nspan + 1
          print*,'density  nspan,r1,r2,i1,i2=',nspan,r(i1),r(i2),i1,i2
          go to 10
        else
          den = co/(r(i2)**3 - r(i1)**3)
          return
        endif  
*        
      else
*
 30     i2 = 2*nspan + 1
*
        if(r(i2).gt.1.d8) then
          nspan = nspan + 1
          print*,'density  nspan,r2,i2=',nspan,r(i2),i2
          co = 3.0d0*float(nspan)
          go to 30
        else
          den = co/r(i2)**3
          return
        endif
      endif  
*
*
      return
*
      end
*
*
*
*

