      subroutine coepot
*
*
*         compute smooth potential for all particles
*         ------------------------------------------
*
*
      include 'common.h'
*
      real*8 ux,uy,smx,rx
*
      integer i,n,im
*
*
      n = nt
*
      smx = smt
      uy = 0.0d0
      rx = r(n)
      ux = -smx/rx
      u(n) = ux
*
      do 10 i=n-1,1,-1
*    
         im = iname(i+1)
         smx = smx - body(im)
         ux = -smx/r(i)
         uy = uy - body(im)/rx
         rx = r(i)
         u(i) = ux + uy
*
 10   continue
*
*     
      return
*
      end
*
*
*
*    
