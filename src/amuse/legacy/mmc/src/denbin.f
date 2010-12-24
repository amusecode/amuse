*
*
      subroutine denbin(nup,nb,denb)
*
*     calculate binary number density in the case of primordial binares
*     --------------------------------------------------------------------
*     density is calculated for binaries which position is inside r(lmax).
*     --------------------------------------------------------------------
*     lmax = nzon(nup)
*     ----------------
*
*
      include 'common.h'
*
cChanged following ftncheck
c      real*8 rb(nmax),denb(nmax),co
      real*8 rb(nmax),denb(nmax),co
*
      integer nup,lmax,i,im1,ib,nspan,nspan2,nmin1,nmax1,nb(nmax)
*
      nspan = 2
      nspan2 = 2*nspan
      co = 0.75d0*(nspan2 + 1)/pi
      lmax = nzst(nup)
      ib = 0
*
      do 10 i = 1,lmax
*
         nb(i) = 0
         im1 = iname(i)
         if(ikind(im1).eq.2) then
           ib = ib + 1
           rb(ib) = r(i)
           nb(i) = ib
           denb(ib) = 0.d0
         endif
*
 10   continue
*
      do 20 i = 1,ib
*
         if(i.le.nspan) then
           nmin1 = 1
           nmax1 = nspan2 + 1
           denb(i) = co/(rb(nmax1)**3 - rb(nmin1)**3)
*
         else
*
           if(i+nspan.gt.ib) then
c             nmin1 = 2*i - ib - nspan
c             denb(i) = co/(rb(ib)**3 - rb(nmin1)**3)
             nmin1 = ib - 2*nspan
             if(nmin1.le.0) then
               print*,' wrong-nmin0 i,ib,nmin,nspan =',i,ib,nmin,nspan
               call flush(6)
               nmin1 = 1
             endif
             nmax1 = ib
             denb(i) = co/(rb(nmax1)**3 - rb(nmin1)**3)
*
           else
*
             nmin1 = i - nspan
             if(nmin1.le.0) then 
               print*,' wrong-nmin1 i,ib,nmin,nspan =',i,ib,nmin,nspan
               call flush(6)
               nmin1 = 1
             endif
             nmax1 = i + nspan
             denb(i) = co/(rb(nmax1)**3 - rb(nmin1)**3)
*
           endif
*
         endif
*
 20      continue
*
*       maybe is good to smooth the binary density
* 
         return
*
         end
*
*