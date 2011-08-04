      subroutine splintd(xa,ya,y2a,n,x,y)
      integer n
      real*4 x,y,xa(n),y2a(n),ya(n)
      integer k,khi,klo
      real*4 a,b,h
      data klo,khi /1,1/
      save khi,klo
c     check first of previous klo and khi are still ok.
      if (khi.le.n) then
         if (xa(khi).gt.x .and. xa(klo).lt.x .and. khi-klo.eq.1) goto 2
      endif
c     if not, search by bisection
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif
 2    h=xa(khi)-xa(klo)
      if (h.le.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*
     *     (h**2)/6.
      return
      END
