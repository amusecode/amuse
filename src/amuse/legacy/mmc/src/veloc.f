*
*
      subroutine veloc(k,velave,smave)
*
*       determine average velocity and mass of stars in the vicinity
*       ------------------------------------------------------------
*       of 'k' star
*       -----------
*
*
      include 'common.h'
*
      real*8 velave,smave,v2
*
      integer k,i,i1,i2,im,n,nspan,inum
*
*
      velave = 0.0d0
      smave = 0.0d0
      n = nt
c      nspan = int(1.5*nminzo/2.0)
c      nspan = 10 
       nspan = 4
*
      if(k.gt.nspan) then
        i1 = k-nspan
        i2 = k+nspan
*
        if(i2.ge.n) then
          i1 = n - 2*nspan - 1
          i2 = n
        endif
*
      else
        i1 = 1
        i2 = 2*nspan + 1
      endif
*
*      compute average velocity and mass of stars from 2*nspan + 1 stars
*
 20   inum = 0
      do 10 i = i1,i2
         im = iname(i)
         inum = inum + 1
         print*,'i,im,inum,vr,vt,sms=',i,im,inum,vr(im),vt(im),body(im)
         v2 = vr(im)**2 + vt(im)**2
         velave = velave + v2*body(im)
         smave = smave + body(im)
 10   continue
*
      velave = velave/smave
      smave = smave/float(inum)
*
      return
      end
*
*
