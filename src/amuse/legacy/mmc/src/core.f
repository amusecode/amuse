       subroutine core
*
*
*     calculate core parameters: rc,vc,nc,smc,roc
*     -------------------------------------------
*
*
      include 'common.h'
*
      real*8 sx,v2x
*
      integer i,n,nx,im,ncut
*
*
      n = nt
*     
*       determination number of stars to calculate the central parameters
*
      ncut = 0.01*n
      if(ncor.le.ncut) then
        nx = ncor
      else
        if(ncut.gt.nmin) then
          nx = ncut
        else
          nx = nmin
        endif
      endif    
*
*
*       determination of the central parameters
*
      sx = 0.0d0
      v2x = 0.0d0
*
      do 10 i=1,nx
         im = iname(i)
         sx = sx + body(im)
   10    v2x = v2x + body(im)*(vr(im)*vr(im) + vt(im)*vt(im))
*
      v2x = v2x/sx
      roc = 3.0d0*sx/(4.0d0*pi*r(nx)**3)
      rc = 3.0d0*v2x/(4.0d0*pi*roc)
      vc = sqrt(v2x)
      rc = sqrt(rc)
*
      i = 0
      sx = 0.0d0
*
   20 continue
      i = i + 1
      im = iname(i)
      sx = sx + body(im)
      if(rc.gt.r(i)) go to 20
      smc = sx - body(im)
      nc = i-1
*
*
      return
*
      end
*
*
*
*
