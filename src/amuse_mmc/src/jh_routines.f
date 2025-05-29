C ***
C       PROGRAM colour
C c-------------------------------------------------------------c
C c
C c     Example program to show how to compute colours 
C c     using Kurucz and Bergeron (WDs) model atmosphere data. 
C c     Requires data files: 
C c         Kurucz.dat
C c         (Kurucz, 1992, Proc. IAU Symp. 149, p. 225)
C c         wdhyd.dat
C c         (Bergeron, Wesemael & Beauchamp, 1995, PASP, 107, 1047)
C c     Written by Jarrod Hurley 10/08/06. 
C c
C c-------------------------------------------------------------c
C       implicit none
C *
C       INTEGER kw
C *
C       REAL*8 z,massi,logl,logr,logte,logz,bolc
C       REAL*8 mv,mb,bminv,vmini,uminb
C       REAL*8 mv1,mv2,bminv1,bminv2
C *
C * Set metallicity. 
C *
C       z = 0.001d0
C       logz = log10(z/0.02)
C *
C * Read in the bolometric correction tables. 
C *
C       CALL ubvtab
C *
C * Mass, luminosity and radius of a star (as well as metallicity) 
C * are required to compute colours. Here are some test values: 
C * 2 Msun turn-off star, Z = 0.001
C       massi = 2.0
C       logl = 1.8355
C       logr = 0.3955
C * 1.3 Msun turn-off star, Z = 0.02
C *     massi = 1.3 
C *     logl = 0.7574    
C *     logr = 0.3369
C * 0.55 Msun cool WD 
C *     massi = 0.55
C *     logl = -3.1694   
C *     logr = -1.8712
C *
C       logte = 3.762d0 + 0.25d0*logl - 0.5d0*logr
C       if(kw.ge.10.and.kw.le.12)then
C          CALL wd2ubv(logl,logte,massi,bolc,mv,
C      &               uminb,bminv,vmini)
C       else
C          CALL lt2ubv(logl,logte,massi,logz,bolc,mv,
C      &               uminb,bminv,vmini)
C       endif
C       mb = mv + bminv
C *
C       WRITE(*,*)
C       WRITE(*,*)' Single Star:'
C       WRITE(*,*)'              Mv ',mv
C       WRITE(*,*)'              Mb ',mb
C       WRITE(*,*)'             U-B ',uminb
C       WRITE(*,*)'             B-V ',bminv
C       WRITE(*,*)'             V-I ',vmini
C *
C * Here is an example for an unresolved binary. 
C *
C       mv1 = mv
C       bminv1 = bminv
C *
C * Take a 0.55 Msun WD companion. 
C       massi = 0.55
C       logl = -3.1694   
C       logr = -1.8712
C *
C       logte = 3.762d0 + 0.25d0*logl - 0.5d0*logr
C       if(kw.ge.10.and.kw.le.12)then
C          CALL wd2ubv(logl,logte,massi,bolc,mv2,
C      &               uminb,bminv2,vmini)
C       else
C          CALL lt2ubv(logl,logte,massi,logz,bolc,mv2,
C      &               uminb,bminv2,vmini)
C       endif
C *
C * Add fluxes to get combined magnitudes. 
C       mv = -2.5d0*log10(10.d0**(-0.4d0*mv1) + 10.d0**(-0.4d0*mv2))
C       mb = -2.5d0*log10(10.d0**(-0.4d0*(mv1+bminv1)) +
C      &                  10.d0**(-0.4d0*(mv2+bminv2)))
C       bminv = mb - mv
C *
C       WRITE(*,*)
C       WRITE(*,*)' Binary Star:'
C       WRITE(*,*)'              Mv ',mv
C       WRITE(*,*)'              Mb ',mb
C       WRITE(*,*)'             B-V ',bminv
C *
C       STOP
C       END
***
      SUBROUTINE ubvtab
      implicit none
*
      INTEGER i,j,k,n
      INTEGER nzgr,ntgr,nggr
      PARAMETER(nzgr=8,ntgr=61,nggr= 11)
      INTEGER ntgr2,nggr2
      PARAMETER(ntgr2=91,nggr2=5)
*
      REAL feh
      REAL*8 zgr(nzgr),tgr(ntgr),ggr(nggr),ubv(nzgr,ntgr,nggr,5)
      COMMON /ubvdata/ zgr,tgr,ggr,ubv
      REAL*8 wtgr(ntgr2),wggr(nggr2),wubv(ntgr2,nggr2,5)
      COMMON /wubvdata/ wtgr,wggr,wubv
      CHARACTER*200 datadir
      COMMON /AMUSE/ datadir

*
*cello, dir mod
      OPEN(23,file=trim(datadir)//'/static/Kurucz.dat',
     &        form='formatted',status='old')
      do k = 1, nzgr
         do i = 1, ntgr
            do j = 1, nggr
               read(23,*)feh,tgr(i),ggr(j),(ubv(k,i,j,n),n=1,5)
            end do
            tgr(i) = log10(tgr(i))
         end do
c....... zgr=log(Z/0.02), assuming X=0.76-3*Z and Z(sun)=0.02
         zgr(k) = -log10((3.d0 + 37.425d0*10.d0**(-feh))/38.d0)
*        zgr(k) = -log10(0.07895 + 0.92105*10.0**(-feh))
      end do
      CLOSE(23)
*
      OPEN(24,file=trim(datadir)//'/static/wdhyd.dat',
     &        form='formatted',status='old')
      do j = 1,nggr2
         do i = 1,ntgr2
            read(24,*)wtgr(i),wggr(j),(wubv(i,j,k),k=1,5)
         enddo
      enddo
      do i = 1,ntgr2
         wtgr(i) = log10(wtgr(i))
      enddo
      CLOSE(24)
*
      RETURN
      END
***
      subroutine lt2ubv(logl,logt,mass,logz,bolc,mv,uminb,bminv,vmini)
      implicit none
c.... computes values of Mv, U-B, B-V and V-I for given log L, log T,
c.... mass and log(Z/0.02)
*
      integer k,ig,it,iz
      integer nzgr,ntgr,nggr
      parameter(nzgr=8,ntgr=61,nggr=11)
      integer indx
      external indx
      real*8 zgr(nzgr),tgr(ntgr),ggr(nggr),ubv(nzgr,ntgr,nggr,5)
      common /ubvdata/ zgr,tgr,ggr,ubv
      real*8 cm(5),gconst
      parameter(gconst=-10.6071d0)
      real*8 logl,logt,mass,logz,mv,uminb,bminv,vmini
      real*8 logm,logg,dg1,dg2,dt1,dt2,dz1,dz2,mbol,bolc
*
      logm = dlog10(mass)
      logg = logm + 4.d0*logt - logl + gconst
c.... find indices of log Z, log g and log T to interpolate between.
c.... don't allow extrapolation outside log Z and log g grid.
      ig = indx(logg,ggr,nggr)
      it = indx(logt,tgr,ntgr)
      iz = indx(logz,zgr,nzgr)
      dg1 = (logg - ggr(ig-1))/(ggr(ig) - ggr(ig-1))
      dg1 = max(0.d0, min(1.d0, dg1))
      dg2 = 1.d0 - dg1
      dt1 = (logt - tgr(it-1))/(tgr(it) - tgr(it-1))
      dt2 = 1.d0 - dt1
      dz1 = (logz - zgr(iz-1))/(zgr(iz) - zgr(iz-1))
      dz1 = max(0.d0, min(1.d0, dz1))
      dz2 = 1.d0 - dz1
      do k = 1, 5
         cm(k) = ((ubv(iz,it,ig,k)*dg1 + ubv(iz,it,ig-1,k)*dg2)*dt1
     &            + (ubv(iz,it-1,ig,k)*dg1 +
     &            ubv(iz,it-1,ig-1,k)*dg2)*dt2)*dz1 +
     &           ((ubv(iz-1,it,ig,k)*dg1 +
     &            ubv(iz-1,it,ig-1,k)*dg2)*dt1 +
     &           (ubv(iz-1,it-1,ig,k)*dg1 +
     &            ubv(iz-1,it-1,ig-1,k)*dg2)*dt2)*dz2
      enddo
      mbol = 4.75d0 - 2.5d0*logl
      bolc = cm(1)
      mv = mbol - bolc
      uminb = cm(2)
      bminv = cm(3)
      vmini = cm(4) + cm(5)
*
      return
      end
***
      integer function indx(ax,xx,nx)
c.....finds index of ax in monotonously increasing or decreasing array xx
      implicit none
      integer nx,j,jl,jh
      real*8 ax,xx(nx),sx
*
      sx = xx(nx) - xx(1)
      jl = 1
      jh = nx
 1    if (jh-jl.gt.1) then
         j = (jh + jl)/2
         if ((ax-xx(j))*sx.gt.0.d0) then
            jl = j
         else
            jh = j
         end if
         goto 1
      end if
      indx = jh
*
      return
      end
***
      subroutine wd2ubv(logl,logt,mass,bolc,mv,uminb,bminv,vmini)
c.... computes values of Mv, U-B, B-V and V-I for given log L, log T,
c.... mass and log(Z/0.02)
      implicit none
      integer k,ig,it
      integer ntgr,nggr
      parameter(ntgr=91,nggr=5)
      integer indx
      external indx
      real*8 tgr(ntgr),ggr(nggr),ubv(ntgr,nggr,5)
      common /wubvdata/ tgr,ggr,ubv
      real*8 cm(5),gconst
      parameter(gconst = -10.6071d0)
      real*8 mbol,bolc,mv,uminb,bminv,vmini
      real*8 logm,logg,dg1,dg2,dt1,dt2,mass,logt,logl
*
      logm = dlog10(mass)
      logg = logm + 4.d0*logt - logl + gconst
c.... find indices of log g and log T to interpolate between.
c.... don't allow extrapolation outside log g grid.
      ig = indx(logg,ggr,nggr)
      it = indx(logt,tgr,ntgr)
      dg1 = (logg - ggr(ig-1))/(ggr(ig) - ggr(ig-1))
      dg1 = MAX(0.d0,MIN(1.d0,dg1))
      dg2 = 1.d0 - dg1
      dt1 = (logt - tgr(it-1))/(tgr(it) - tgr(it-1))
      dt2 = 1.d0 - dt1
      do k = 1,5
         cm(k) = (ubv(it,ig,k)*dg1 + ubv(it,ig-1,k)*dg2)*dt1
     &         + (ubv(it-1,ig,k)*dg1 + ubv(it-1,ig-1,k)*dg2)*dt2
      enddo
      mbol = 4.75d0 - 2.5d0*logl
      bolc = cm(1)
      mv = mbol - bolc
      uminb = cm(2)
      bminv = cm(3)
      vmini = cm(4) + cm(5)
      return
      end
***
