
      subroutine initialize(z_in,
     &			    neta_in, bwind_in, hewind_in, sigma_in,
     &			    ifflag_in, wdflag_in, bhflag_in,
     &			    nsflag_in, mxns_in,
     &			    pts1_in, pts2_in, pts3_in,
     &     		    status)
cf2py intent(in) z_in
cf2py intent(in) neta_in, bwind_in, hewind_in, sigma_in
cf2py intent(in) ifflag_in, wdflag_in, bhflag_in, nsflag_in, mxns_in
cf2py intent(in) pts1_in, pts2_in, pts3_in
cf2py intent(out) status
      implicit none
      real*8 z_in, z, zpars(20)
      real*8 neta_in, bwind_in, hewind_in, sigma_in
      integer ifflag_in, wdflag_in, bhflag_in, nsflag_in, mxns_in
      real*8 pts1_in, pts2_in, pts3_in
      integer status
      include 'src/const_bse.h'
      common /SSE_init/ z, zpars

c     Input parameters are passed from MUSE, rather than being read here.

      z = z_in
      neta = neta_in
      bwind = bwind_in
      hewind = hewind_in
      sigma = sigma_in
      ifflag = ifflag_in
      wdflag = wdflag_in
      bhflag = bhflag_in
      nsflag = nsflag_in
      mxns = mxns_in
      pts1 = pts1_in
      pts2 = pts2_in
      pts3 = pts3_in
      
      call zcnsts(z, zpars)
      if(idum.gt.0) idum = -idum

      status = 0
      return
      end

      subroutine evolve0(kw, mass, mt, r, lum, mc, rc, menv, renv,
     &                   ospin, epoch, tm, tphys, tphysf)

c     Should work, accordig to the f2py documentation, but doesn't...

      implicit none
      integer,intent(inout) :: kw
      real*8,intent(inout) :: mass,mt,r,lum,mc,rc,menv,renv,ospin
      real*8,intent(inout) :: epoch,tm,tphys,tphysf
      real*8 z, zpars(20)
      real*8 dtp
      common /SSE_init/ z, zpars
cf2py intent(inout) mass
cf2py intent(inout) kw, mt, lum, mc, rc, menv, renv, ospin
cf2py intent(inout) epoch, tm, tphys, tphysf
cf2py intent(inout) r

      dtp = tphys+1
      print *, kw, mass, mt, r, lum
      call evolv1(kw, mass, mt, r, lum, mc, rc, menv, renv,
     &            ospin, epoch, tm, tphys, tphysf, dtp, z, zpars)
      print *, kw, mass, mt, r, lum
      kw = 2

      return
      end
      
      subroutine evolve(kw, mass, mt, r, lum, mc, rc, menv, renv,
     &                  ospin, epoch, tm, tphys, tphysf,
     &              kw1, mass1, mt1, r1, lum1, mc1, rc1, menv1, renv1,
     &              ospin1, epoch1, tm1, tphys1, tphysf1)

c     Ugly, but it works!

      implicit none
      integer kw
      real*8 mass,mt,r,lum,mc,rc,menv,renv,ospin
      real*8 epoch,tm,tphys,tphysf
cf2py intent(in) kw, mass, mt, r, lum, mc, rc, menv, renv, ospin
cf2py intent(in) epoch, tm, tphys, tphysf
      integer kw1
      real*8 mass1,mt1,r1,lum1,mc1,rc1,menv1,renv1,ospin1
      real*8 epoch1,tm1,tphys1,tphysf1
cf2py intent(out) kw1, mass1, mt1, r1, lum1, mc1, rc1, menv1, renv1, ospin1
cf2py intent(out) epoch1, tm1, tphys1, tphysf1

      real*8 z, zpars(20)
      real*8 dtp
      common /SSE_init/ z, zpars

c     Copying everything is probably overkill, but not wrong.

      kw1 = kw
      mass1 = mass
      mt1 = mt
      r1 = r
      lum1 = lum
      mc1 = mc
      rc1 = rc
      menv1 = menv
      renv1 = renv
      ospin1 = ospin
      epoch1 = epoch
      tm1 = tm
      tphys1 = tphys
      tphysf1 = tphysf

      dtp = tphys1+1
      call evolv1(kw1, mass1, mt1, r1, lum1, mc1, rc1, menv1, renv1,
     &            ospin1, epoch1, tm1, tphys1, tphysf1, dtp, z, zpars)

      return
      end
      
      subroutine get_time_step(kw, mass, age, mt, tm, epoch, dt)
cf2py intent(out) dt
cf2py intent(in) kw, mass, age, mt, tm, epoch
      implicit none
      integer kw
      real*8 mass, age, mt, tm, epoch
      real*8 dt
      real*8 tscls(20), lums(10), GB(10), tn, dtm, dtr
      real*8 z, zpars(20)
      common /SSE_init/ z, zpars

!     Call star fuction to get stellar parameters
      call star(kw, mass, mt, tm, tn, tscls, lums, GB, zpars)

!     Call deltat function to get next timestep
      call deltat(kw, age-epoch, tm, tn, tscls, dtm, dtr)
      dt = min(dtr, dtm)

      return
      end
