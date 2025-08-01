      subroutine initialize(z_in,
     &			    neta_in, bwind_in, hewind_in, sigma1_in,
     &			    sigma2_in, ifflag_in, wdflag_in, bhflag_in,
     &			    nsflag_in, piflag_in, mxns_in, idum_in,
     &			    pts1_in, pts2_in, pts3_in,
     &     		    status)
cf2py intent(in) z_in
cf2py intent(in) neta_in, bwind_in, hewind_in, sigma1_in, sigma2_in
cf2py intent(in) ifflag_in, wdflag_in, bhflag_in, nsflag_in, piflag_in
cf2py intent(in) mxns_in, idum_in
cf2py intent(in) pts1_in, pts2_in, pts3_in
cf2py intent(out) status
      implicit none
      real*8 z_in, z, zpars(20)
      real*8 neta_in, bwind_in, hewind_in, sigma1_in, sigma2_in
      integer ifflag_in, wdflag_in, bhflag_in, nsflag_in, piflag_in
      integer idum_in 
      real*8 mxns_in, pts1_in, pts2_in, pts3_in
      integer status
      include 'src/const_mobse.h'
      common /SSE_init/ z, zpars

c     Input parameters are passed from MUSE, rather than being read here.

      z = z_in
      neta = neta_in
      bwind = bwind_in
      hewind = hewind_in
      sigma1 = sigma1_in
      sigma2 = sigma2_in
      ifflag = ifflag_in
      wdflag = wdflag_in
      bhflag = bhflag_in
      nsflag = nsflag_in
      mxns = mxns_in
      idum = idum_in
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

      dtp = tphysf+1
      if (mass.eq.0.0) mass=mt
      call evolv1(kw, mass, mt, r, lum, mc, rc, menv, renv,
     &            ospin, epoch, tm, tphys, tphysf, dtp, z, zpars)

      return
      end
      
      subroutine evolve_star(kw, mass, mt, r, lum, mc, rc, menv, renv,
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

! fix for the very small timestep problem (where dt ~ 1.e-16 * age )
! this line mirrors the corresponding lne in evolve1.f      
      dt = MAX(dt,1.0d-07*age)
      return
      end
      
      
      subroutine get_mass_loss_wind(kw, lum, r, mt, mc, mlout)
      implicit none
      integer kw
      real*8 lum,r,mt,mc,rl
      real*8 mlout
      real*8 mlwind
      real*8 z, zpars(20)
      common /SSE_init/ z, zpars
      rl = 0.0
      mlout = mlwind(kw, lum, r, mt, mc, rl, z)
      
      return
      end


      subroutine get_gyration_radius(kw, mass, mt, r, lum, epoch, tm, 
     &                              tphys, rg)
      implicit none
      integer kw
      real*8 mass, aj, mt, tm, epoch, r, lum, tphys, mc, rc, menv, renv
      real*8 tscls(20), lums(10), GB(10), tn
      real*8 z, zpars(20)
      real*8 rg2,rg
      common /SSE_init/ z, zpars

!     Call star fuction to get stellar parameters
      call star(kw, mass, mt, tm, tn, tscls, lums, GB, zpars)
      
      aj = tphys - epoch

!     Call hrdiag fuction to get additional stellar parameters
      CALL hrdiag(mass, aj, mt, tm, tn, tscls, lums, GB, zpars,
     &               r, lum, kw, mc, rc, menv, renv, rg2)

      rg = SQRT(rg2)

      return
      end
